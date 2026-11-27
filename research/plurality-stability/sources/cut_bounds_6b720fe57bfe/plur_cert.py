#!/usr/bin/env python3
"""
Uniform certificate for the radial operators in the three-set simplex proof.

The program verifies the two weighted Hilbert--Schmidt estimates uniformly for

    1/43 <= sigma <= 1/2.

It uses Arb balls (through python-flint) for every transcendental table and
outward-rounded binary64 interval arithmetic for the remaining elementary
operations.  The finite square is enclosed by interval Riemann sums.  The
complement is controlled by the analytic estimates printed at the end.

Install the sole nonstandard dependency with

    python3 -m pip install "python-flint>=0.8"

The default parameters are intended for the manuscript certificate.  The
``--quick`` option is useful only as a smoke test; it deliberately uses a
coarser subdivision and need not certify the desired inequalities.
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from dataclasses import dataclass

import numpy as np

try:
    from flint import arb, ctx
except ImportError as exc:  # pragma: no cover - exercised only without flint
    raise SystemExit(
        "plur_cert.py requires python-flint>=0.8; install it with\n"
        '  python3 -m pip install "python-flint>=0.8"'
    ) from exc


SIGMA_MIN = 1.0 / 43.0
SIGMA_MAX = 0.5
UNIT_ROUNDOFF = 2.0**-53


def down(x):
    """Round a scalar or ndarray one binary64 step toward -infinity."""
    return np.nextafter(x, -np.inf)


def up(x):
    """Round a scalar or ndarray one binary64 step toward +infinity."""
    return np.nextafter(x, np.inf)


def add_lo(a, b):
    return down(a + b)


def add_hi(a, b):
    return up(a + b)


def sub_lo(a, b):
    return down(a - b)


def sub_hi(a, b):
    return up(a - b)


def mul_pos_lo(a, b):
    # Every caller uses quantities proved nonnegative.  Clamping the directed
    # result at zero prevents a one-ulp negative artifact at an exact zero
    # from being reused as though it were a genuinely nonnegative endpoint.
    return np.maximum(0.0, down(a * b))


def mul_pos_hi(a, b):
    return up(a * b)


def div_pos_lo(a, b):
    return np.maximum(0.0, down(a / b))


def div_pos_hi(a, b):
    return up(a / b)


def sqrt_pos_lo(x):
    lo, _ = bounds(ball(float(x), float(x)).sqrt())
    return lo


def sqrt_pos_hi(x):
    _, hi = bounds(ball(float(x), float(x)).sqrt())
    return hi


def ball(lo: float, hi: float) -> arb:
    """Return an Arb ball containing the closed binary64 interval [lo, hi]."""
    if not (math.isfinite(lo) and math.isfinite(hi) and lo <= hi):
        raise ValueError(f"invalid interval [{lo}, {hi}]")
    # arb(float) contains the exact binary64 input.  Taking the interval hull
    # avoids a potentially unsafe binary64 midpoint/radius construction when
    # the interval is narrow but far from the origin.
    return arb(lo).union(arb(hi))


def bounds(x: arb) -> tuple[float, float]:
    """Extract outward-rounded binary64 endpoints from an Arb ball."""
    lo = math.nextafter(float(x.lower()), -math.inf)
    hi = math.nextafter(float(x.upper()), math.inf)
    return lo, hi


@dataclass(frozen=True)
class Constants:
    """Directed binary64 enclosures of the transcendental constants."""

    j1_lo: float
    j1_hi: float
    c0_lo: float
    e_lo: float


def build_constants() -> Constants:
    pi = arb.pi()
    j1 = (9 - pi * arb(3).sqrt()) / (4 * pi * pi)
    c0 = 5 * pi * (3 * pi).sqrt() / 1728
    e = arb(1).exp()
    j1_lo, j1_hi = bounds(j1)
    c0_lo, _ = bounds(c0)
    e_lo, _ = bounds(e)
    if min(j1_lo, c0_lo, e_lo) <= 0:
        raise ArithmeticError("failed to enclose the positive constants")
    return Constants(j1_lo, j1_hi, c0_lo, e_lo)


def sum_enclosure(values: np.ndarray) -> tuple[float, float]:
    """Enclose a binary64 sum, including all round-to-nearest accumulation error.

    For n terms, the standard gamma_n bound gives
    |fl(sum)-sum| <= gamma_n sum |x_i|.  When n*u < 1/8, the deliberately
    looser radius 4*n*u*max(1, fl(sum |x_i|)) dominates that error as well as
    the error in the auxiliary absolute-value sum.
    """

    flat = np.asarray(values, dtype=np.float64).ravel()
    n = int(flat.size)
    if n == 0:
        return 0.0, 0.0
    if n * UNIT_ROUNDOFF >= 0.125:
        raise ArithmeticError("summation is too long for the certified error bound")
    total = float(np.sum(flat, dtype=np.float64))
    abs_total = float(np.sum(np.abs(flat), dtype=np.float64))
    scale = max(1.0, abs_total)
    error_ball = (
        arb(4)
        * n
        * (arb(1) / (2**53))
        * ball(scale, scale)
    )
    _, error_hi = bounds(error_ball)
    return sub_lo(total, error_hi), add_hi(total, error_hi)


def arb_exp_neg_square_upper(x: float) -> float:
    """Upper bound for exp(-x^2), with x interpreted as a binary64 endpoint."""
    xb = ball(x, x)
    _, value_hi = bounds((-(xb * xb)).exp())
    return value_hi


@dataclass
class SparseRange:
    """Static range-minimum/range-maximum tables."""

    mins: list[np.ndarray]
    maxs: list[np.ndarray]

    @classmethod
    def build(cls, lo: np.ndarray, hi: np.ndarray) -> "SparseRange":
        mins = [np.asarray(lo, dtype=np.float64)]
        maxs = [np.asarray(hi, dtype=np.float64)]
        step = 1
        while 2 * step <= len(lo):
            prev_lo = mins[-1]
            prev_hi = maxs[-1]
            mins.append(np.minimum(prev_lo[:-step], prev_lo[step:]))
            maxs.append(np.maximum(prev_hi[:-step], prev_hi[step:]))
            step *= 2
        return cls(mins, maxs)

    def query(self, first: np.ndarray, last: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Enclose every inclusive index range [first, last]."""
        length = last - first + 1
        k = np.floor(np.log2(length)).astype(np.int64)
        lo = np.empty(first.shape, dtype=np.float64)
        hi = np.empty(first.shape, dtype=np.float64)
        for level in range(len(self.mins)):
            mask = k == level
            if not np.any(mask):
                continue
            span = 1 << level
            i = first[mask]
            j = last[mask] - span + 1
            lo[mask] = np.minimum(self.mins[level][i], self.mins[level][j])
            hi[mask] = np.maximum(self.maxs[level][i], self.maxs[level][j])
        return down(lo), up(hi)


@dataclass
class UniformTable:
    """Interval table on adjacent closed cells, with range queries."""

    edges: np.ndarray
    ranges: dict[str, SparseRange]

    @property
    def bins(self) -> int:
        return len(self.edges) - 1

    @property
    def xmax(self) -> float:
        return float(self.edges[-1])

    def index_range(self, lo, hi) -> tuple[np.ndarray, np.ndarray]:
        lo = np.asarray(lo, dtype=np.float64)
        hi = np.asarray(hi, dtype=np.float64)
        if np.any(hi > self.xmax):
            excess = float(np.max(hi - self.xmax))
            raise ValueError(
                f"table query exceeds its certified domain by {excess:.3e}"
            )
        # All table arguments are proved nonnegative.  A negative lower
        # endpoint can only be an outward-rounding artifact at zero.
        lo = np.maximum(0.0, lo)
        if np.any(lo > hi):
            raise ValueError("table query lies outside its domain")
        first = np.searchsorted(self.edges, down(lo), side="right") - 1
        last = np.searchsorted(self.edges, up(hi), side="right") - 1
        return np.clip(first, 0, self.bins - 1), np.clip(last, 0, self.bins - 1)

    def query(self, name: str, lo, hi) -> tuple[np.ndarray, np.ndarray]:
        first, last = self.index_range(lo, hi)
        return self.ranges[name].query(first, last)


@dataclass
class Jet2:
    """Value and first two derivatives with respect to one real variable."""

    value: arb
    d1: arb
    d2: arb

    @classmethod
    def constant(cls, value) -> "Jet2":
        return cls(arb(value), arb(0), arb(0))

    @classmethod
    def variable(cls, value) -> "Jet2":
        return cls(value if isinstance(value, arb) else arb(value), arb(1), arb(0))

    @staticmethod
    def promote(other) -> "Jet2":
        return other if isinstance(other, Jet2) else Jet2.constant(other)

    def __add__(self, other):
        other = self.promote(other)
        return Jet2(
            self.value + other.value, self.d1 + other.d1, self.d2 + other.d2
        )

    __radd__ = __add__

    def __neg__(self):
        return Jet2(-self.value, -self.d1, -self.d2)

    def __sub__(self, other):
        return self + (-self.promote(other))

    def __rsub__(self, other):
        return self.promote(other) - self

    def __mul__(self, other):
        other = self.promote(other)
        return Jet2(
            self.value * other.value,
            self.d1 * other.value + self.value * other.d1,
            self.d2 * other.value
            + 2 * self.d1 * other.d1
            + self.value * other.d2,
        )

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        if isinstance(scalar, Jet2):
            raise TypeError("Jet2 division is implemented only for constants")
        return Jet2(self.value / scalar, self.d1 / scalar, self.d2 / scalar)

    def exp(self):
        value = self.value.exp()
        return Jet2(
            value,
            value * self.d1,
            value * (self.d2 + self.d1 * self.d1),
        )

    def sin(self):
        s = self.value.sin()
        c = self.value.cos()
        return Jet2(
            s,
            c * self.d1,
            c * self.d2 - s * self.d1 * self.d1,
        )

    def cos(self):
        s = self.value.sin()
        c = self.value.cos()
        return Jet2(
            c,
            -s * self.d1,
            -s * self.d2 - c * self.d1 * self.d1,
        )

    def erfc(self):
        sqrtpi = arb.pi().sqrt()
        gaussian = (-(self.value * self.value)).exp()
        first = -2 * gaussian / sqrtpi
        second = 4 * self.value * gaussian / sqrtpi
        return Jet2(
            self.value.erfc(),
            first * self.d1,
            second * self.d1 * self.d1 + first * self.d2,
        )


def max_abs_upper(value: arb) -> float:
    lo, hi = bounds(value)
    return math.nextafter(max(abs(lo), abs(hi)), math.inf)


def midpoint_integral_piece(
    midpoint_value: arb, second_derivative: arb, width: arb
) -> arb:
    """Rigorous midpoint-rule enclosure on one interval."""
    derivative_bound = max_abs_upper(second_derivative)
    error_ball = (
        width
        * width
        * width
        * ball(derivative_bound, derivative_bound)
        / 24
    )
    _, error_hi = bounds(error_ball)
    return width * midpoint_value + arb(0, error_hi)


def angular_integrands(parameter: arb, v: Jet2, kind: str) -> tuple[Jet2, ...]:
    """Angular integrands and their first two v-derivatives."""
    pi = arb.pi()
    sqrtpi = pi.sqrt()
    c = (pi * v).cos()
    if kind == "a":
        return (((2 + 3 * v) * (-parameter * c).exp()) / 2,)

    b = parameter * c
    erfcx = (b * b).exp() * b.erfc()
    psi = 1 - sqrtpi * b * erfcx
    neg_dpsi = sqrtpi * (1 + 2 * b * b) * erfcx - 2 * b
    g = v * (1 - 3 * v / 2)
    return (
        (2 + 3 * v) * psi / 2,
        pi * g * (pi * v).sin() * neg_dpsi / 2,
    )


def arb_v_integrals(parameter: arb, v_panels: int, kind: str) -> tuple[arb, ...]:
    """Second-order interval midpoint enclosure on [0,2/3]."""
    denominator = 3 * v_panels
    width = arb(2) / denominator
    values = [arb(0), arb(0)] if kind == "q" else [arb(0)]
    for j in range(v_panels):
        vlo = arb(2 * j) / denominator
        vhi = arb(2 * (j + 1)) / denominator
        midpoint = arb(2 * j + 1) / denominator
        center_jets = angular_integrands(
            parameter, Jet2.variable(midpoint), kind
        )
        range_jets = angular_integrands(
            parameter, Jet2.variable(vlo.union(vhi)), kind
        )
        for index in range(len(values)):
            values[index] += midpoint_integral_piece(
                center_jets[index].value, range_jets[index].d2, width
            )
    return tuple(values)


def graded_edges(xmax: float, bins: int, power: float) -> np.ndarray:
    mesh = np.linspace(0.0, 1.0, bins + 1)
    edges = xmax * mesh**power
    edges[0] = 0.0
    edges[-1] = xmax
    if np.any(np.diff(edges) <= 0):
        raise ValueError("the graded table mesh is not strictly increasing")
    return edges


def build_a_table(
    amax: float, bins: int, v_panels: int, table_power: float
) -> UniformTable:
    t_lo = np.empty(bins)
    t_hi = np.empty(bins)
    p_lo = np.empty(bins)
    p_hi = np.empty(bins)
    edges = graded_edges(amax, bins, table_power)
    for i in range(bins):
        lo, hi = float(edges[i]), float(edges[i + 1])
        a = ball(lo, hi)
        (t_value,) = arb_v_integrals(a, v_panels, "a")
        p_value = a * (-a).exp()
        t_lo[i], t_hi[i] = bounds(t_value)
        p_lo[i], p_hi[i] = bounds(p_value)
        # The exact functions are nonnegative for a >= 0.  Intersecting with
        # that proved range also makes subsequent positive interval products
        # legitimate when Arb dependency gives a slightly negative endpoint.
        t_lo[i] = max(0.0, t_lo[i])
        p_lo[i] = max(0.0, p_lo[i])
    return UniformTable(
        edges,
        {"T": SparseRange.build(t_lo, t_hi), "P": SparseRange.build(p_lo, p_hi)},
    )


def build_q_table(
    qmax: float, bins: int, v_panels: int, table_power: float
) -> UniformTable:
    u_lo = np.empty(bins)
    u_hi = np.empty(bins)
    h_lo = np.empty(bins)
    h_hi = np.empty(bins)
    s_lo = np.empty(bins)
    s_hi = np.empty(bins)
    edges = graded_edges(qmax, bins, table_power)
    sqrtpi = arb.pi().sqrt()
    for i in range(bins):
        lo, hi = float(edges[i]), float(edges[i + 1])
        q = ball(lo, hi)
        u_value, h_value = arb_v_integrals(q, v_panels, "q")
        erfcx = (q * q).exp() * q.erfc()
        s_value = sqrtpi * (1 + 2 * q * q) * erfcx - 2 * q
        u_lo[i], u_hi[i] = bounds(u_value)
        h_lo[i], h_hi[i] = bounds(h_value)
        s_lo[i], s_hi[i] = bounds(s_value)
        # U, H, and S have nonnegative integral representations on q >= 0.
        u_lo[i] = max(0.0, u_lo[i])
        h_lo[i] = max(0.0, h_lo[i])
        s_lo[i] = max(0.0, s_lo[i])
    return UniformTable(
        edges,
        {
            "U": SparseRange.build(u_lo, u_hi),
            "H": SparseRange.build(h_lo, h_hi),
            "S": SparseRange.build(s_lo, s_hi),
        },
    )


def build_exp_table(zmax: float, bins: int, table_power: float) -> UniformTable:
    lo_values = np.empty(bins)
    hi_values = np.empty(bins)
    edges = graded_edges(zmax, bins, table_power)
    for i in range(bins):
        z = ball(float(edges[i]), float(edges[i + 1]))
        lo_values[i], hi_values[i] = bounds((-z).exp())
        lo_values[i] = max(0.0, lo_values[i])
    return UniformTable(
        edges, {"E": SparseRange.build(lo_values, hi_values)}
    )


def abs_square_upper(lo, hi):
    return up(np.maximum(np.abs(lo), np.abs(hi)) ** 2)


def mul_interval(lo1, hi1, lo2, hi2):
    """Outward product of two real binary64 intervals."""
    candidates = np.stack(
        (lo1 * lo2, lo1 * hi2, hi1 * lo2, hi1 * hi2), axis=0
    )
    return down(np.min(candidates, axis=0)), up(np.max(candidates, axis=0))


def gaussian_outside_square(lam: float, cutoff: float) -> float:
    """Union/Mills upper bound for int outside [0,X]^2 exp(-lam|x|^2)."""
    if lam <= 0 or cutoff <= 0:
        raise ValueError("the Gaussian tail parameters must be positive")
    lb = ball(lam, lam)
    xb = ball(cutoff, cutoff)
    value = (
        arb.pi().sqrt()
        * (-(lb * xb * xb)).exp()
        / (2 * lb * lb.sqrt() * xb)
    )
    _, value_hi = bounds(value)
    return value_hi


@dataclass
class SlabResult:
    sigma_lo: float
    sigma_hi: float
    c_lo: float
    c_hi: float
    upper_finite: float
    upper_tail: float
    lower_finite: float
    lower_tail: float

    @property
    def upper_total(self):
        return up(self.upper_finite + self.upper_tail)

    @property
    def lower_total(self):
        return up(self.lower_finite + self.lower_tail)


def certify_slab(
    sigma_lo: float,
    sigma_hi: float,
    cutoff: float,
    x_panels: int,
    chunk: int,
    a_table: UniformTable,
    q_table: UniformTable,
    exp_table: UniformTable,
    constants: Constants,
) -> SlabResult:
    x_edges = np.linspace(0.0, cutoff, x_panels + 1, dtype=np.float64)
    x_edges[0] = 0.0
    x_edges[-1] = cutoff
    if np.any(np.diff(x_edges) <= 0):
        raise ValueError("the radial mesh is not strictly increasing")
    x_lo = x_edges[:-1]
    x_hi = x_edges[1:]
    width_lo = sub_lo(x_hi, x_lo)
    width_hi = sub_hi(x_hi, x_lo)

    sigma_lo_sq = mul_pos_lo(sigma_lo, sigma_lo)
    sigma_hi_sq = mul_pos_hi(sigma_hi, sigma_hi)
    d_lo = sub_lo(1.0, sigma_hi_sq)
    d_hi = sub_hi(1.0, sigma_lo_sq)
    invd_lo = div_pos_lo(1.0, d_hi)
    invd_hi = div_pos_hi(1.0, d_lo)
    kap_lo = div_pos_lo(sigma_lo_sq, d_hi)
    kap_hi = div_pos_hi(sigma_hi_sq, d_lo)
    alpha_lo = div_pos_lo(sigma_lo, sqrt_pos_hi(d_hi))
    alpha_hi = div_pos_hi(sigma_hi, sqrt_pos_lo(d_lo))
    beta_lo = mul_pos_lo(mul_pos_lo(2.0, sigma_lo), invd_lo)
    beta_hi = mul_pos_hi(mul_pos_hi(2.0, sigma_hi), invd_hi)

    q_lo = mul_pos_lo(alpha_lo, x_lo)
    q_hi = mul_pos_hi(alpha_hi, x_hi)
    ubase_lo, ubase_hi = q_table.query("U", q_lo, q_hi)
    hbase_lo, hbase_hi = q_table.query("H", q_lo, q_hi)
    sbase_lo, sbase_hi = q_table.query("S", q_lo, q_hi)

    zx_lo = mul_pos_lo(kap_lo, mul_pos_lo(x_lo, x_lo))
    zx_hi = mul_pos_hi(kap_hi, mul_pos_hi(x_hi, x_hi))
    ex_lo, ex_hi = exp_table.query("E", zx_lo, zx_hi)

    # u = exp(-kappa x^2) [U(q) + 2 j_(1) q S(q)].
    two_j1_lo = mul_pos_lo(2.0, constants.j1_lo)
    two_j1_hi = mul_pos_hi(2.0, constants.j1_hi)
    corr_u_lo = mul_pos_lo(two_j1_lo, mul_pos_lo(q_lo, sbase_lo))
    corr_u_hi = mul_pos_hi(two_j1_hi, mul_pos_hi(q_hi, sbase_hi))
    u_lo = mul_pos_lo(ex_lo, add_lo(ubase_lo, corr_u_lo))
    u_hi = mul_pos_hi(ex_hi, add_hi(ubase_hi, corr_u_hi))

    # J_*(x^2)/x = exp(-kappa x^2) sigma/sqrt(d)
    #                  [H(q) + j_(1) S(q)/2].
    half_j1_lo = div_pos_lo(constants.j1_lo, 2.0)
    half_j1_hi = div_pos_hi(constants.j1_hi, 2.0)
    jt_inner_lo = add_lo(hbase_lo, mul_pos_lo(half_j1_lo, sbase_lo))
    jt_inner_hi = add_hi(hbase_hi, mul_pos_hi(half_j1_hi, sbase_hi))
    jt_lo = mul_pos_lo(ex_lo, mul_pos_lo(alpha_lo, jt_inner_lo))
    jt_hi = mul_pos_hi(ex_hi, mul_pos_hi(alpha_hi, jt_inner_hi))
    if np.any(jt_lo <= 0):
        raise ArithmeticError("the interval lower bound for J_*/x is not positive")

    # Compute C_*-1 by integrating u_*-1.  This algebraically preserves the
    # cancellation at sigma=0; integrating u_* directly would introduce a
    # sigma-independent rectangle error and make a uniform certificate near
    # sigma=1/43 unnecessarily expensive.
    wx_z_lo = mul_pos_lo(x_lo, x_lo)
    wx_z_hi = mul_pos_hi(x_hi, x_hi)
    wexp_lo, wexp_hi = exp_table.query("E", wx_z_lo, wx_z_hi)
    weight_c_lo = mul_pos_lo(mul_pos_lo(2.0, x_lo), wexp_lo)
    weight_c_hi = mul_pos_hi(mul_pos_hi(2.0, x_hi), wexp_hi)
    centered_u_lo = down(u_lo - 1.0)
    centered_u_hi = up(u_hi - 1.0)
    c_integrand_lo, c_integrand_hi = mul_interval(
        weight_c_lo, weight_c_hi, centered_u_lo, centered_u_hi
    )
    c_piece_lo, c_piece_hi = mul_interval(
        width_lo, width_hi, c_integrand_lo, c_integrand_hi
    )
    c_sum_lo, _ = sum_enclosure(c_piece_lo)
    _, c_sum_hi = sum_enclosure(c_piece_hi)
    # Since 0 <= u_* <= 27/25, one has |u_*-1| <= 1, so the
    # omitted centered C_* integral is at most exp(-cutoff^2).
    c_tail = arb_exp_neg_square_upper(cutoff)
    c_lo = sub_lo(add_lo(1.0, c_sum_lo), c_tail)
    c_hi = add_hi(add_hi(1.0, c_sum_hi), c_tail)
    if c_lo <= 0:
        raise ArithmeticError("the interval lower bound for C_* is not positive")

    # The Gaussian factor in the transformed Hilbert--Schmidt integrals.
    weight_lo, weight_hi = exp_table.query(
        "E", mul_pos_lo(x_lo, x_lo), mul_pos_hi(x_hi, x_hi)
    )

    upper_sum = 0.0
    lower_sum = 0.0
    ylo = x_lo[None, :]
    yhi = x_hi[None, :]
    for start in range(0, x_panels, chunk):
        stop = min(start + chunk, x_panels)
        xlo = x_lo[start:stop, None]
        xhi = x_hi[start:stop, None]

        prod_lo = mul_pos_lo(xlo, ylo)
        prod_hi = mul_pos_hi(xhi, yhi)
        a_lo = mul_pos_lo(beta_lo, prod_lo)
        a_hi = mul_pos_hi(beta_hi, prod_hi)
        t_lo, t_hi = a_table.query("T", a_lo, a_hi)
        p_lo, p_hi = a_table.query("P", a_lo, a_hi)

        sumsq_lo = add_lo(mul_pos_lo(xlo, xlo), mul_pos_lo(ylo, ylo))
        sumsq_hi = add_hi(mul_pos_hi(xhi, xhi), mul_pos_hi(yhi, yhi))
        z_lo = mul_pos_lo(kap_lo, sumsq_lo)
        z_hi = mul_pos_hi(kap_hi, sumsq_hi)
        ek_lo, ek_hi = exp_table.query("E", z_lo, z_hi)

        k0_lo = mul_pos_lo(mul_pos_lo(invd_lo, ek_lo), t_lo)
        k0_hi = mul_pos_hi(mul_pos_hi(invd_hi, ek_hi), t_hi)
        # 4 B = (2 j_(1)/d) exp(-kappa(x^2+y^2)) a exp(-a).
        bcoef_lo = mul_pos_lo(two_j1_lo, invd_lo)
        bcoef_hi = mul_pos_hi(two_j1_hi, invd_hi)
        b_lo = mul_pos_lo(
            mul_pos_lo(bcoef_lo, ek_lo), p_lo
        )
        b_hi = mul_pos_hi(
            mul_pos_hi(bcoef_hi, ek_hi), p_hi
        )
        k_lo = add_lo(k0_lo, b_lo)
        k_hi = add_hi(k0_hi, b_hi)

        ui_lo = u_lo[start:stop, None]
        ui_hi = u_hi[start:stop, None]
        uj_lo = u_lo[None, :]
        uj_hi = u_hi[None, :]

        # E = K - u(x) - u(y) + C.
        e_lo = add_lo(sub_lo(sub_lo(k_lo, ui_hi), uj_hi), c_lo)
        e_hi = add_hi(sub_hi(sub_hi(k_hi, ui_lo), uj_lo), c_hi)
        e2_hi = abs_square_upper(e_lo, e_hi)

        # R = K - u(x)u(y)/C.
        rank_lo = div_pos_lo(mul_pos_lo(ui_lo, uj_lo), c_hi)
        rank_hi = div_pos_hi(mul_pos_hi(ui_hi, uj_hi), c_lo)
        r_lo = sub_lo(k_lo, rank_hi)
        r_hi = sub_hi(k_hi, rank_lo)
        r2_hi = abs_square_upper(r_lo, r_hi)

        denom_lo = mul_pos_lo(jt_lo[start:stop, None], jt_lo[None, :])
        gauss_hi = mul_pos_hi(weight_hi[start:stop, None], weight_hi[None, :])
        common_hi = div_pos_hi(mul_pos_hi(4.0, gauss_hi), denom_lo)
        area_hi = mul_pos_hi(
            width_hi[start:stop, None], width_hi[None, :]
        )
        upper_terms = mul_pos_hi(
            area_hi, mul_pos_hi(common_hi, e2_hi)
        )
        lower_terms = mul_pos_hi(
            area_hi, mul_pos_hi(common_hi, r2_hi)
        )
        _, upper_chunk_hi = sum_enclosure(upper_terms)
        _, lower_chunk_hi = sum_enclosure(lower_terms)
        upper_sum = add_hi(upper_sum, upper_chunk_hi)
        lower_sum = add_hi(lower_sum, lower_chunk_hi)

    # Uniform analytic tails.  Put alpha=sigma/sqrt(1-sigma^2).
    jtail_lo = mul_pos_lo(constants.c0_lo, alpha_lo)
    lam_kernel = div_pos_lo(1.0, add_hi(1.0, sigma_hi))
    four_j1_hi = mul_pos_hi(4.0, constants.j1_hi)
    three_e_lo = mul_pos_lo(3.0, constants.e_lo)
    correction_hi = div_pos_hi(four_j1_hi, three_e_lo)
    k_envelope_hi = mul_pos_hi(add_hi(1.0, correction_hi), invd_hi)
    kernel_ratio_hi = div_pos_hi(k_envelope_hi, jtail_lo)
    kernel_tail = mul_pos_hi(
        mul_pos_hi(kernel_ratio_hi, kernel_ratio_hi),
        gaussian_outside_square(lam_kernel, cutoff),
    )

    lam_add = sub_lo(1.0, kap_hi)
    add_mass = gaussian_outside_square(lam_add, cutoff)
    e_add = div_pos_hi(81.0, 25.0)
    u_global_hi = div_pos_hi(27.0, 25.0)
    r_add = div_pos_hi(mul_pos_hi(u_global_hi, u_global_hi), c_lo)
    # 4 * (K+A)^2 <= 8 K^2 + 8 A^2.
    e_ratio_hi = div_pos_hi(e_add, jtail_lo)
    r_ratio_hi = div_pos_hi(r_add, jtail_lo)
    e_add_tail = mul_pos_hi(
        mul_pos_hi(e_ratio_hi, e_ratio_hi), add_mass
    )
    r_add_tail = mul_pos_hi(
        mul_pos_hi(r_ratio_hi, r_ratio_hi), add_mass
    )
    upper_tail = add_hi(
        mul_pos_hi(8.0, kernel_tail), mul_pos_hi(8.0, e_add_tail)
    )
    lower_tail = add_hi(
        mul_pos_hi(8.0, kernel_tail), mul_pos_hi(8.0, r_add_tail)
    )

    return SlabResult(
        sigma_lo,
        sigma_hi,
        c_lo,
        c_hi,
        upper_sum,
        upper_tail,
        lower_sum,
        lower_tail,
    )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cutoff", type=float, default=5.0)
    parser.add_argument("--x-panels", type=int, default=1000)
    parser.add_argument("--sigma-panels", type=int, default=640)
    parser.add_argument(
        "--sigma-power",
        type=float,
        default=1.3,
        help="power-law grading of sigma slabs toward sigma=1/43",
    )
    parser.add_argument("--v-panels", type=int, default=64)
    parser.add_argument("--q-bins", type=int, default=4096)
    parser.add_argument("--a-bins", type=int, default=8192)
    parser.add_argument("--exp-bins", type=int, default=65536)
    parser.add_argument(
        "--table-power",
        type=float,
        default=2.0,
        help="power-law grading of q, a, and exponential tables toward zero",
    )
    parser.add_argument("--chunk", type=int, default=48)
    parser.add_argument("--precision", type=int, default=100)
    parser.add_argument(
        "--quick",
        action="store_true",
        help="use a coarse smoke-test grid (not the manuscript certificate)",
    )
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.quick:
        args.x_panels = 120
        args.sigma_panels = 12
        args.v_panels = 64
        args.q_bins = 512
        args.a_bins = 1024
        args.exp_bins = 4096
        args.chunk = 40

    if args.cutoff <= 0 or min(
        args.x_panels,
        args.sigma_panels,
        args.v_panels,
        args.q_bins,
        args.a_bins,
        args.exp_bins,
        args.chunk,
        args.precision,
    ) <= 0:
        raise SystemExit("all cutoffs and subdivision counts must be positive")
    if args.sigma_power <= 0 or args.table_power <= 0:
        raise SystemExit("both mesh-grading powers must be positive")

    ctx.prec = args.precision
    constants = build_constants()
    started = time.perf_counter()
    sigma_table_hi = math.nextafter(SIGMA_MAX, math.inf)
    sigma_sq_hi = mul_pos_hi(sigma_table_hi, sigma_table_hi)
    d_table_lo = sub_lo(1.0, sigma_sq_hi)
    invd_table_hi = div_pos_hi(1.0, d_table_lo)
    alpha_table_hi = div_pos_hi(
        sigma_table_hi, sqrt_pos_lo(d_table_lo)
    )
    beta_table_hi = mul_pos_hi(
        mul_pos_hi(2.0, sigma_table_hi), invd_table_hi
    )
    cutoff_sq_hi = mul_pos_hi(args.cutoff, args.cutoff)
    qmax = mul_pos_hi(alpha_table_hi, args.cutoff)
    amax = mul_pos_hi(beta_table_hi, cutoff_sq_hi)
    kap_table_hi = mul_pos_hi(sigma_sq_hi, invd_table_hi)
    radial_sum_hi = add_hi(cutoff_sq_hi, cutoff_sq_hi)
    zmax = max(
        cutoff_sq_hi,
        mul_pos_hi(kap_table_hi, radial_sum_hi),
    )
    zmax = up(zmax)

    print(
        f"building Arb tables: q bins={args.q_bins}, a bins={args.a_bins}, "
        f"exp bins={args.exp_bins}, angular panels={args.v_panels}",
        flush=True,
    )
    q_table = build_q_table(qmax, args.q_bins, args.v_panels, args.table_power)
    a_table = build_a_table(amax, args.a_bins, args.v_panels, args.table_power)
    exp_table = build_exp_table(zmax, args.exp_bins, args.table_power)

    mesh = np.linspace(0.0, 1.0, args.sigma_panels + 1)
    edges = SIGMA_MIN + (SIGMA_MAX - SIGMA_MIN) * mesh**args.sigma_power
    edges[0] = SIGMA_MIN
    edges[-1] = SIGMA_MAX
    if np.any(np.diff(edges) <= 0):
        raise SystemExit("the sigma mesh is not strictly increasing")
    results: list[SlabResult] = []
    for i in range(args.sigma_panels):
        result = certify_slab(
            math.nextafter(float(edges[i]), -math.inf),
            math.nextafter(float(edges[i + 1]), math.inf),
            args.cutoff,
            args.x_panels,
            args.chunk,
            a_table,
            q_table,
            exp_table,
            constants,
        )
        results.append(result)
        if (
            args.quick
            or i < 3
            or i + 1 == args.sigma_panels
            or (i + 1) % max(1, args.sigma_panels // 20) == 0
        ):
            print(
                f"slab {i + 1:3d}/{args.sigma_panels}: "
                f"sigma=[{result.sigma_lo:.10f},{result.sigma_hi:.10f}] "
                f"upper H^2<={result.upper_total:.10f} "
                f"lower H^2<={result.lower_total:.10f}",
                flush=True,
            )

    upper = max(results, key=lambda item: item.upper_total)
    lower = max(results, key=lambda item: item.lower_total)
    c_lo = min(item.c_lo for item in results)
    c_hi = max(item.c_hi for item in results)
    elapsed = time.perf_counter() - started

    print()
    print(
        f"sigma interval                  [{SIGMA_MIN:.12f}, {SIGMA_MAX:.12f}]"
    )
    print(
        f"x cutoff / panels              {args.cutoff:g} / {args.x_panels}"
    )
    print(f"sigma panels                   {args.sigma_panels}")
    print(f"uniform C_* enclosure          [{c_lo:.12g}, {c_hi:.12g}]")
    print(
        f"upper H^2 finite worst          {upper.upper_finite:.12g} "
        f"on [{upper.sigma_lo:.10f},{upper.sigma_hi:.10f}]"
    )
    print(f"upper H^2 analytic tail         {upper.upper_tail:.12g}")
    print(f"upper H^2 certified             {upper.upper_total:.12g}")
    print(f"upper H certified               {sqrt_pos_hi(upper.upper_total):.12g}")
    print(
        f"lower H^2 finite worst          {lower.lower_finite:.12g} "
        f"on [{lower.sigma_lo:.10f},{lower.sigma_hi:.10f}]"
    )
    print(f"lower H^2 analytic tail         {lower.lower_tail:.12g}")
    print(f"lower H^2 certified             {lower.lower_total:.12g}")
    print(f"lower H certified               {sqrt_pos_hi(lower.lower_total):.12g}")
    print(f"runtime                         {elapsed:.2f} seconds")

    certified = upper.upper_total < 4.0 and lower.lower_total < 4.0
    if certified:
        print(
            "CERTIFIED: uniform upper mean-zero H<2 and unrestricted lower H<2."
        )
        return 0
    qualifier = "QUICK GRID" if args.quick else "FULL GRID"
    print(f"NOT CERTIFIED on the {qualifier}; refine the subdivisions.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
