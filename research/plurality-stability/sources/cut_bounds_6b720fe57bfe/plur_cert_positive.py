#!/usr/bin/env python3
"""
Uniform directed certificate for the positive-correlation three-set appendix.

Together with the previously known range 0 < rho < 1/10, this program
certifies the weighted Hilbert--Schmidt bound for

    1/10 <= rho <= 2/5.

The script imports the Arb and outward-rounded binary64 interval primitives
from plur_cert.py.  It uses the sharpened first angular mode, the
rho-dependent additive center displayed in the appendix, interval Riemann
sums on [0,7]^2, and analytic Gaussian tails.
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from dataclasses import dataclass

import numpy as np
from flint import arb, ctx

import plur_cert as base


RHO_MIN = 0.1
RHO_MAX = 0.4

# If ell_rho(x)=sum_{n=0}^4 c_n(rho)x^n, each row below lists the
# coefficients of 1,rho,...,rho^4 in c_n.  The terminating decimals are
# interpreted as exact binary64 numbers; no approximation property of the
# center is assumed by the proof.
CENTER = (
    (
        0.5001276609531384,
        -0.04839192730072833,
        0.3420752274381793,
        -0.02840693452661011,
        0.3929215251696178,
    ),
    (
        0.0020711604332492072,
        0.11278766504365427,
        0.42921652164107776,
        -1.4633497642462785,
        2.392753763954644,
    ),
    (
        -0.00657914712574693,
        0.1533613436474771,
        -1.6052209990260096,
        4.917034579755538,
        -7.5025494501524435,
    ),
    (
        0.00509284321106939,
        -0.11817713264541818,
        0.9939294369107601,
        -4.02653695935449,
        4.798761639440846,
    ),
    (
        -0.00111590812093279,
        0.02578645463758624,
        -0.21531279149872498,
        0.776584990984738,
        -0.6977289068083498,
    ),
)


@dataclass(frozen=True)
class PositiveConstants:
    j1_lo: float
    j1_hi: float
    e_lo: float
    cw_lo: float


def build_constants() -> PositiveConstants:
    pi = arb.pi()
    j1 = (9 - pi * arb(3).sqrt()) / (4 * pi * pi)
    cw = (6 - pi) / (4 * pi * pi.sqrt())
    e = arb(1).exp()
    j1_lo, j1_hi = base.bounds(j1)
    cw_lo, _ = base.bounds(cw)
    e_lo, _ = base.bounds(e)
    if min(j1_lo, cw_lo, e_lo) <= 0:
        raise ArithmeticError("failed to enclose the positive constants")
    return PositiveConstants(j1_lo, j1_hi, e_lo, cw_lo)


def positive_a_integrand(parameter: arb, v: base.Jet2) -> base.Jet2:
    pi = arb.pi()
    return ((2 + 3 * v) * (parameter * (pi * v).cos()).exp()) / 2


def positive_q_integrand(parameter: arb, v: base.Jet2) -> base.Jet2:
    """Integrand for W_0(x^2)/(x exp(-kappa x^2) rho/sqrt(d))."""
    pi = arb.pi()
    sqrtpi = pi.sqrt()
    b = -parameter * (pi * v).cos()
    erfcx = (b * b).exp() * b.erfc()
    neg_dpsi = sqrtpi * (1 + 2 * b * b) * erfcx - 2 * b
    g = v * (1 - 3 * v / 2)
    return pi * g * (pi * v).sin() * neg_dpsi / 2


def midpoint_integral(
    parameter: arb, v_panels: int, integrand
) -> arb:
    denominator = 3 * v_panels
    width = arb(2) / denominator
    value = arb(0)
    for j in range(v_panels):
        vlo = arb(2 * j) / denominator
        vhi = arb(2 * (j + 1)) / denominator
        midpoint = arb(2 * j + 1) / denominator
        center = integrand(parameter, base.Jet2.variable(midpoint))
        enclosure = integrand(
            parameter, base.Jet2.variable(vlo.union(vhi))
        )
        value += base.midpoint_integral_piece(
            center.value, enclosure.d2, width
        )
    return value


def build_positive_a_table(
    amax: float, bins: int, v_panels: int, table_power: float
) -> base.UniformTable:
    t_lo = np.empty(bins)
    t_hi = np.empty(bins)
    p_lo = np.empty(bins)
    p_hi = np.empty(bins)
    edges = base.graded_edges(amax, bins, table_power)
    for i in range(bins):
        lo, hi = float(edges[i]), float(edges[i + 1])
        a = base.ball(lo, hi)
        value = midpoint_integral(a, v_panels, positive_a_integrand)
        damping = a * (-a).exp()
        t_lo[i], t_hi[i] = base.bounds(value)
        p_lo[i], p_hi[i] = base.bounds(damping)
        t_lo[i] = max(0.0, t_lo[i])
        p_lo[i] = max(0.0, p_lo[i])
    return base.UniformTable(
        edges,
        {
            "TPLUS": base.SparseRange.build(t_lo, t_hi),
            "P": base.SparseRange.build(p_lo, p_hi),
        },
    )


def build_positive_q_table(
    qmax: float, bins: int, v_panels: int, table_power: float
) -> base.UniformTable:
    h_lo = np.empty(bins)
    h_hi = np.empty(bins)
    s_lo = np.empty(bins)
    s_hi = np.empty(bins)
    edges = base.graded_edges(qmax, bins, table_power)
    sqrtpi = arb.pi().sqrt()
    for i in range(bins):
        lo, hi = float(edges[i]), float(edges[i + 1])
        q = base.ball(lo, hi)
        h_value = midpoint_integral(q, v_panels, positive_q_integrand)
        erfcx = (q * q).exp() * q.erfc()
        s_value = sqrtpi * (1 + 2 * q * q) * erfcx - 2 * q
        h_lo[i], h_hi[i] = base.bounds(h_value)
        s_lo[i], s_hi[i] = base.bounds(s_value)
        h_lo[i] = max(0.0, h_lo[i])
        s_lo[i] = max(0.0, s_lo[i])
    return base.UniformTable(
        edges,
        {
            "HPLUS": base.SparseRange.build(h_lo, h_hi),
            "S": base.SparseRange.build(s_lo, s_hi),
        },
    )


def interval_horner(coefficients, x_lo, x_hi):
    value_lo = base.down(float(coefficients[-1]))
    value_hi = base.up(float(coefficients[-1]))
    for coefficient in reversed(coefficients[:-1]):
        value_lo, value_hi = base.mul_interval(
            value_lo, value_hi, x_lo, x_hi
        )
        value_lo = base.add_lo(value_lo, float(coefficient))
        value_hi = base.add_hi(value_hi, float(coefficient))
    return value_lo, value_hi


def center_coefficients(rho_lo: float, rho_hi: float):
    return tuple(
        interval_horner(row, rho_lo, rho_hi) for row in CENTER
    )


def center_interval(rho_lo: float, rho_hi: float, x_lo, x_hi):
    coefficient_intervals = center_coefficients(rho_lo, rho_hi)
    value_lo, value_hi = coefficient_intervals[-1]
    for coefficient_lo, coefficient_hi in reversed(
        coefficient_intervals[:-1]
    ):
        value_lo, value_hi = base.mul_interval(
            value_lo, value_hi, x_lo, x_hi
        )
        value_lo = base.add_lo(value_lo, coefficient_lo)
        value_hi = base.add_hi(value_hi, coefficient_hi)
    return value_lo, value_hi


def gaussian_even_moments(lam: float, cutoff: float, degree: int = 4):
    """Upper bounds for full and tail moments of x^(2n) exp(-lam*x^2)."""
    lb = base.ball(lam, lam)
    xb = base.ball(cutoff, cutoff)
    exponential = (-(lb * xb * xb)).exp()
    full = [arb.pi().sqrt() / (2 * lb.sqrt())]
    tail = [exponential / (2 * lb * xb)]
    for n in range(1, degree + 1):
        full.append((2 * n - 1) * full[-1] / (2 * lb))
        tail.append(
            xb ** (2 * n - 1) * exponential / (2 * lb)
            + (2 * n - 1) * tail[-1] / (2 * lb)
        )
    return (
        [base.bounds(value)[1] for value in full],
        [base.bounds(value)[1] for value in tail],
    )


def correlated_gaussian_outside(
    rho_hi: float, d_hi: float, cutoff: float
) -> float:
    """Bound the complement using completion of the square.

    For Q=(x^2+y^2-4*rho*x*y)/(1-rho^2), extending the inner half-line
    to the whole real line and then using the union bound gives

      int_outside exp(-Q)
      <= sqrt(pi*d) exp(-lambda*X^2)/(lambda*X),

    where lambda=(1-4*rho^2)/d.
    """
    rho_ball = base.ball(rho_hi, rho_hi)
    d_ball = base.ball(d_hi, d_hi)
    x_ball = base.ball(cutoff, cutoff)
    numerator = 1 - 4 * rho_ball * rho_ball
    lambda_ball = numerator / d_ball
    lambda_lo, _ = base.bounds(lambda_ball)
    if lambda_lo <= 0:
        raise ArithmeticError("positive tail completion requires rho<1/2")
    lb = base.ball(lambda_lo, lambda_lo)
    value = (
        arb.pi().sqrt()
        * d_ball.sqrt()
        * (-(lb * x_ball * x_ball)).exp()
        / (lb * x_ball)
    )
    return base.bounds(value)[1]


def analytic_tail(
    rho_lo: float,
    rho_hi: float,
    d_hi: float,
    invd_hi: float,
    kappa_hi: float,
    alpha_lo: float,
    cutoff: float,
    constants: PositiveConstants,
) -> tuple[float, float, float]:
    jtail_lo = base.mul_pos_lo(constants.cw_lo, alpha_lo)

    correction_hi = base.div_pos_hi(
        base.mul_pos_hi(2.0, constants.j1_hi), constants.e_lo
    )
    envelope_hi = base.mul_pos_hi(
        base.add_hi(1.0, correction_hi), invd_hi
    )
    ratio_hi = base.div_pos_hi(envelope_hi, jtail_lo)
    kernel_tail = base.mul_pos_hi(
        base.mul_pos_hi(ratio_hi, ratio_hi),
        correlated_gaussian_outside(rho_hi, d_hi, cutoff),
    )

    lam_add = base.sub_lo(1.0, kappa_hi)
    full_moments, tail_moments = gaussian_even_moments(
        lam_add, cutoff
    )
    coefficient_intervals = center_coefficients(rho_lo, rho_hi)
    coefficient_squares = []
    for lo, hi in coefficient_intervals:
        magnitude = base.up(max(abs(float(lo)), abs(float(hi))))
        coefficient_squares.append(base.mul_pos_hi(magnitude, magnitude))

    ell_full = 0.0
    ell_tail = 0.0
    for coefficient_square, full, tail in zip(
        coefficient_squares, full_moments, tail_moments
    ):
        ell_full = base.add_hi(
            ell_full, base.mul_pos_hi(coefficient_square, full)
        )
        ell_tail = base.add_hi(
            ell_tail, base.mul_pos_hi(coefficient_square, tail)
        )
    # (sum_{n=0}^4 c_n x^n)^2 <= 5 sum c_n^2 x^(2n).
    ell_full = base.mul_pos_hi(5.0, ell_full)
    ell_tail = base.mul_pos_hi(5.0, ell_tail)
    mass_full = full_moments[0]
    mass_tail = tail_moments[0]
    pair_numerator = base.add_hi(
        base.mul_pos_hi(ell_tail, mass_full),
        base.mul_pos_hi(ell_full, mass_tail),
    )
    jtail_square_lo = base.mul_pos_lo(jtail_lo, jtail_lo)
    additive_pair = base.div_pos_hi(pair_numerator, jtail_square_lo)

    # 4 E^2 <= 12 K^2 + 12 ell(x)^2 + 12 ell(y)^2.
    kernel_contribution = base.mul_pos_hi(12.0, kernel_tail)
    additive_contribution = base.mul_pos_hi(24.0, additive_pair)
    return (
        base.add_hi(kernel_contribution, additive_contribution),
        kernel_contribution,
        additive_contribution,
    )


@dataclass
class SlabResult:
    rho_lo: float
    rho_hi: float
    finite: float
    tail: float
    kernel_tail: float
    additive_tail: float

    @property
    def total(self):
        return base.add_hi(self.finite, self.tail)


def certify_slab(
    rho_lo: float,
    rho_hi: float,
    cutoff: float,
    x_panels: int,
    chunk: int,
    a_table: base.UniformTable,
    q_table: base.UniformTable,
    exp_table: base.UniformTable,
    constants: PositiveConstants,
) -> SlabResult:
    x_edges = np.linspace(0.0, cutoff, x_panels + 1, dtype=np.float64)
    x_edges[0] = 0.0
    x_edges[-1] = cutoff
    x_lo = x_edges[:-1]
    x_hi = x_edges[1:]
    width_hi = base.sub_hi(x_hi, x_lo)

    rho_lo_sq = base.mul_pos_lo(rho_lo, rho_lo)
    rho_hi_sq = base.mul_pos_hi(rho_hi, rho_hi)
    d_lo = base.sub_lo(1.0, rho_hi_sq)
    d_hi = base.sub_hi(1.0, rho_lo_sq)
    invd_lo = base.div_pos_lo(1.0, d_hi)
    invd_hi = base.div_pos_hi(1.0, d_lo)
    kappa_lo = base.mul_pos_lo(rho_lo_sq, invd_lo)
    kappa_hi = base.mul_pos_hi(rho_hi_sq, invd_hi)
    alpha_lo = base.div_pos_lo(rho_lo, base.sqrt_pos_hi(d_hi))
    alpha_hi = base.div_pos_hi(rho_hi, base.sqrt_pos_lo(d_lo))
    beta_lo = base.mul_pos_lo(base.mul_pos_lo(2.0, rho_lo), invd_lo)
    beta_hi = base.mul_pos_hi(base.mul_pos_hi(2.0, rho_hi), invd_hi)

    q_lo = base.mul_pos_lo(alpha_lo, x_lo)
    q_hi = base.mul_pos_hi(alpha_hi, x_hi)
    h_lo, h_hi = q_table.query("HPLUS", q_lo, q_hi)
    s_lo, s_hi = q_table.query("S", q_lo, q_hi)

    half_j_lo = base.div_pos_lo(constants.j1_lo, 2.0)
    half_j_hi = base.div_pos_hi(constants.j1_hi, 2.0)
    inner_lo = base.add_lo(h_lo, base.mul_pos_lo(half_j_lo, s_lo))
    inner_hi = base.add_hi(h_hi, base.mul_pos_hi(half_j_hi, s_hi))
    jbar_lo = base.mul_pos_lo(alpha_lo, inner_lo)
    jbar_hi = base.mul_pos_hi(alpha_hi, inner_hi)
    if np.any(jbar_lo <= 0):
        raise ArithmeticError("the lower interval bound for scaled J_+/x is not positive")

    center_lo, center_hi = center_interval(
        rho_lo, rho_hi, x_lo, x_hi
    )
    weight_lo, weight_hi = exp_table.query(
        "E", base.mul_pos_lo(x_lo, x_lo), base.mul_pos_hi(x_hi, x_hi)
    )

    total = 0.0
    ylo = x_lo[None, :]
    yhi = x_hi[None, :]
    for start in range(0, x_panels, chunk):
        stop = min(start + chunk, x_panels)
        xlo = x_lo[start:stop, None]
        xhi = x_hi[start:stop, None]

        product_lo = base.mul_pos_lo(xlo, ylo)
        product_hi = base.mul_pos_hi(xhi, yhi)
        a_lo = base.mul_pos_lo(beta_lo, product_lo)
        a_hi = base.mul_pos_hi(beta_hi, product_hi)
        t_lo, t_hi = a_table.query("TPLUS", a_lo, a_hi)
        p_lo, p_hi = a_table.query("P", a_lo, a_hi)

        sumsq_lo = base.add_lo(
            base.mul_pos_lo(xlo, xlo), base.mul_pos_lo(ylo, ylo)
        )
        sumsq_hi = base.add_hi(
            base.mul_pos_hi(xhi, xhi), base.mul_pos_hi(yhi, yhi)
        )
        two_j_lo = base.mul_pos_lo(2.0, constants.j1_lo)
        two_j_hi = base.mul_pos_hi(2.0, constants.j1_hi)
        bracket_lo = base.sub_lo(t_lo, base.mul_pos_hi(two_j_hi, p_hi))
        bracket_hi = base.sub_hi(t_hi, base.mul_pos_lo(two_j_lo, p_lo))
        kbar_lo, kbar_hi = base.mul_interval(
            invd_lo, invd_hi, bracket_lo, bracket_hi
        )

        half_z_lo = base.mul_pos_lo(
            0.5, base.mul_pos_lo(kappa_lo, sumsq_lo)
        )
        half_z_hi = base.mul_pos_hi(
            0.5, base.mul_pos_hi(kappa_hi, sumsq_hi)
        )
        eminus_lo, eminus_hi = exp_table.query(
            "E", half_z_lo, half_z_hi
        )
        eplus_lo = base.div_pos_lo(1.0, eminus_hi)
        eplus_hi = base.div_pos_hi(1.0, eminus_lo)
        scaled_kernel_lo, scaled_kernel_hi = base.mul_interval(
            eminus_lo, eminus_hi, kbar_lo, kbar_hi
        )

        center_pair_lo = base.add_lo(
            center_lo[start:stop, None], center_lo[None, :]
        )
        center_pair_hi = base.add_hi(
            center_hi[start:stop, None], center_hi[None, :]
        )
        scaled_center_lo, scaled_center_hi = base.mul_interval(
            eplus_lo, eplus_hi, center_pair_lo, center_pair_hi
        )
        residual_lo = base.sub_lo(scaled_kernel_lo, scaled_center_hi)
        residual_hi = base.sub_hi(scaled_kernel_hi, scaled_center_lo)
        residual_square_hi = base.abs_square_upper(
            residual_lo, residual_hi
        )

        denominator_lo = base.mul_pos_lo(
            jbar_lo[start:stop, None], jbar_lo[None, :]
        )
        gaussian_hi = base.mul_pos_hi(
            weight_hi[start:stop, None], weight_hi[None, :]
        )
        common_hi = base.div_pos_hi(
            base.mul_pos_hi(4.0, gaussian_hi), denominator_lo
        )
        area_hi = base.mul_pos_hi(
            width_hi[start:stop, None], width_hi[None, :]
        )
        terms = base.mul_pos_hi(
            area_hi, base.mul_pos_hi(common_hi, residual_square_hi)
        )
        row_indices = np.arange(start, stop, dtype=np.int64)[:, None]
        column_indices = np.arange(x_panels, dtype=np.int64)[None, :]
        upper_triangle = column_indices >= row_indices
        strict_upper = column_indices > row_indices
        symmetric_terms = np.where(
            upper_triangle,
            np.where(strict_upper, base.mul_pos_hi(2.0, terms), terms),
            0.0,
        )
        _, chunk_hi = base.sum_enclosure(symmetric_terms)
        total = base.add_hi(total, chunk_hi)

    tail, kernel_tail, additive_tail = analytic_tail(
        rho_lo,
        rho_hi,
        d_hi,
        invd_hi,
        kappa_hi,
        alpha_lo,
        cutoff,
        constants,
    )
    return SlabResult(
        rho_lo,
        rho_hi,
        total,
        tail,
        kernel_tail,
        additive_tail,
    )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho-min", type=float, default=RHO_MIN)
    parser.add_argument("--rho-max", type=float, default=RHO_MAX)
    parser.add_argument("--cutoff", type=float, default=7.0)
    parser.add_argument("--x-panels", type=int, default=1600)
    parser.add_argument(
        "--rho-panels",
        type=int,
        default=30,
        help="initial rho panels; failing slabs are bisected adaptively",
    )
    parser.add_argument(
        "--min-rho-width",
        type=float,
        default=1.0e-5,
        help="stop with failure if a noncertifying slab reaches this width",
    )
    parser.add_argument("--v-panels", type=int, default=64)
    parser.add_argument("--q-bins", type=int, default=4096)
    parser.add_argument("--a-bins", type=int, default=16384)
    parser.add_argument("--exp-bins", type=int, default=65536)
    parser.add_argument("--table-power", type=float, default=2.0)
    parser.add_argument("--chunk", type=int, default=48)
    parser.add_argument("--precision", type=int, default=100)
    parser.add_argument("--quick", action="store_true")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    if args.quick:
        args.x_panels = 120
        args.rho_panels = 12
        args.v_panels = 48
        args.q_bins = 512
        args.a_bins = 1024
        args.exp_bins = 4096
        args.chunk = 40

    ctx.prec = args.precision
    if not (RHO_MIN <= args.rho_min < args.rho_max <= RHO_MAX):
        raise SystemExit(
            f"require {RHO_MIN} <= rho-min < rho-max <= {RHO_MAX}"
        )
    constants = build_constants()
    started = time.perf_counter()

    rho_hi = math.nextafter(args.rho_max, math.inf)
    rho_hi_sq = base.mul_pos_hi(rho_hi, rho_hi)
    d_lo = base.sub_lo(1.0, rho_hi_sq)
    invd_hi = base.div_pos_hi(1.0, d_lo)
    alpha_hi = base.div_pos_hi(rho_hi, base.sqrt_pos_lo(d_lo))
    beta_hi = base.mul_pos_hi(base.mul_pos_hi(2.0, rho_hi), invd_hi)
    cutoff_sq = base.mul_pos_hi(args.cutoff, args.cutoff)
    qmax = base.mul_pos_hi(alpha_hi, args.cutoff)
    amax = base.mul_pos_hi(beta_hi, cutoff_sq)
    kappa_hi = base.mul_pos_hi(rho_hi_sq, invd_hi)
    zmax = base.up(
        max(
            cutoff_sq,
            base.mul_pos_hi(
                kappa_hi, base.add_hi(cutoff_sq, cutoff_sq)
            ),
        )
    )

    print(
        f"building positive Arb tables: q bins={args.q_bins}, "
        f"a bins={args.a_bins}, exp bins={args.exp_bins}, "
        f"angular panels={args.v_panels}",
        flush=True,
    )
    q_table = build_positive_q_table(
        qmax, args.q_bins, args.v_panels, args.table_power
    )
    a_table = build_positive_a_table(
        amax, args.a_bins, args.v_panels, args.table_power
    )
    exp_table = base.build_exp_table(
        zmax, args.exp_bins, args.table_power
    )

    edges = np.linspace(args.rho_min, args.rho_max, args.rho_panels + 1)
    edges[0] = args.rho_min
    edges[-1] = args.rho_max
    results = []
    pending = [
        (float(edges[i]), float(edges[i + 1]))
        for i in reversed(range(args.rho_panels))
    ]
    attempts = 0
    failed_at_minimum = False
    while pending:
        rho_lo_exact, rho_hi_exact = pending.pop()
        result = certify_slab(
            math.nextafter(rho_lo_exact, -math.inf),
            math.nextafter(rho_hi_exact, math.inf),
            args.cutoff,
            args.x_panels,
            args.chunk,
            a_table,
            q_table,
            exp_table,
            constants,
        )
        attempts += 1
        if result.total < 4.0 or args.quick:
            results.append(result)
        elif rho_hi_exact - rho_lo_exact > args.min_rho_width:
            midpoint = (rho_lo_exact + rho_hi_exact) / 2.0
            # LIFO order: append the right child first so the left is checked next.
            pending.append((midpoint, rho_hi_exact))
            pending.append((rho_lo_exact, midpoint))
        else:
            results.append(result)
            failed_at_minimum = True
        if args.quick or attempts <= 3 or attempts % 20 == 0 or not pending:
            print(
                f"attempt {attempts:4d}: "
                f"rho=[{result.rho_lo:.10f},{result.rho_hi:.10f}] "
                f"H_+^2<={result.total:.10f} "
                f"accepted={result.total < 4.0 or args.quick}",
                flush=True,
            )

    worst = max(results, key=lambda item: item.total)
    elapsed = time.perf_counter() - started
    print()
    print(
        f"rho interval                    "
        f"[{args.rho_min:.12f}, {args.rho_max:.12f}]"
    )
    print(f"x cutoff / panels              {args.cutoff:g} / {args.x_panels}")
    print(f"initial rho panels             {args.rho_panels}")
    print(f"accepted rho slabs             {len(results)}")
    print(f"total slab attempts            {attempts}")
    print(
        f"H_+^2 finite worst              {worst.finite:.12g} "
        f"on [{worst.rho_lo:.10f},{worst.rho_hi:.10f}]"
    )
    print(f"H_+^2 kernel tail               {worst.kernel_tail:.12g}")
    print(f"H_+^2 additive tail             {worst.additive_tail:.12g}")
    print(f"H_+^2 certified                 {worst.total:.12g}")
    print(f"H_+ certified                   {base.sqrt_pos_hi(worst.total):.12g}")
    print(f"runtime                         {elapsed:.2f} seconds")

    if worst.total < 4.0 and not failed_at_minimum:
        print("CERTIFIED: uniform positive mean-zero H_+<2.")
        return 0
    qualifier = "QUICK GRID" if args.quick else "FULL GRID"
    print(f"NOT CERTIFIED on the {qualifier}; refine the subdivisions.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
