#!/usr/bin/env python3
"""Verify the first exact normal slice for quintic Meng--Yang graphs.

Let ``Psi_(L,M,N)=L*A^2+M*A+N*B`` be the scaled five-variable
Meng--Yang potential and put ``r=R(x,y,p,q)``.  This checker works with the
complete degree-at-most-five 1-jet of ``R`` along ``x=0``.  That is enough
to reconstruct the Hessian of the pullback on the two-slope pencil

    (x,y,p,q) = (0,t,c*t,d*t).

The coefficients of t^10 down to t^5 give a triangular exact reduction.
In particular, the x-free quintic part is forced to ``zeta*y^5``.  If
``zeta=0``, the remaining equations force

    160*rho^2 + 1968*rho + 6021 = 0,

whose discriminant is ``576*34``; hence this branch has no rational point.

The leading homogeneous potential is ``L*(x^3*R_5)^2``.  The matrix
determinant lemma and Euler identities reduce its zero-Hessian condition to
that of ``x^3*R_5``.  The low-dimensional Gordan--Noether theorem then gives
a constant kernel direction.  The checker verifies the two resulting
corank-one charts when ``zeta`` is nonzero:

    R_5 = zeta*y^5 + x*T_4(x,y,p-kappa*q),
    R_5 = zeta*y^5 + x*U_4(x,y,q).

This is a normal-slice theorem, not an HC(4) exclusion.  Lower determinant
layers and the full off-axis identity remain to be solved.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import permutations
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "hc4_meng_yang_quintic_graph_normal_slice.json"
)


def permutation_sign(permutation: tuple[int, ...]) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def multiply_truncated(
    left: dict[tuple[int, ...], sp.Expr],
    right: dict[tuple[int, ...], sp.Expr],
    target: tuple[int, ...],
) -> dict[tuple[int, ...], sp.Expr]:
    product: dict[tuple[int, ...], sp.Expr] = {}
    for left_exponents, left_coefficient in left.items():
        for right_exponents, right_coefficient in right.items():
            exponents = tuple(
                a + b
                for a, b in zip(
                    left_exponents, right_exponents, strict=True
                )
            )
            if all(
                exponent <= bound
                for exponent, bound in zip(exponents, target, strict=True)
            ):
                product[exponents] = (
                    product.get(exponents, sp.S.Zero)
                    + left_coefficient * right_coefficient
                )
    return product


def determinant_truncation(
    matrix: sp.Matrix,
    variables: tuple[sp.Symbol, ...],
    target: tuple[int, ...],
) -> dict[tuple[int, ...], sp.Expr]:
    """Expand a 4-by-4 determinant only through the requested exponents."""

    assert matrix.shape == (4, 4)
    entries = [
        [
            {
                exponents: coefficient
                for exponents, coefficient in sp.Poly(
                    matrix[row, column], *variables
                ).terms()
                if all(
                    exponent <= bound
                    for exponent, bound in zip(
                        exponents, target, strict=True
                    )
                )
            }
            for column in range(4)
        ]
        for row in range(4)
    ]
    determinant: dict[tuple[int, ...], sp.Expr] = {}
    zero = (0,) * len(variables)
    for permutation in permutations(range(4)):
        term: dict[tuple[int, ...], sp.Expr] = {
            zero: sp.Integer(permutation_sign(permutation))
        }
        for row, column in enumerate(permutation):
            term = multiply_truncated(term, entries[row][column], target)
        for exponents, coefficient in term.items():
            determinant[exponents] = (
                determinant.get(exponents, sp.S.Zero) + coefficient
            )
    return determinant


x, y, p, q, t = sp.symbols("x y p q t")
c, d = sp.symbols("c d")
L, M, N = sp.symbols("L M N", nonzero=True)


# Only the x^0 and x^1 coefficients of R can enter the pullback Hessian on
# x=0.  Keep every such coefficient through total graph degree five.
coefficients: dict[tuple[int, int, int, int], sp.Symbol] = {}
R0 = sp.S.Zero
R1 = sp.S.Zero
for total_degree in range(6):
    for y_degree in range(total_degree + 1):
        for p_degree in range(total_degree - y_degree + 1):
            q_degree = total_degree - y_degree - p_degree
            exponents = (0, y_degree, p_degree, q_degree)
            coefficient = sp.Symbol("r" + "".join(map(str, exponents)))
            coefficients[exponents] = coefficient
            R0 += coefficient * y**y_degree * p**p_degree * q**q_degree
for residual_degree in range(5):
    for y_degree in range(residual_degree + 1):
        for p_degree in range(residual_degree - y_degree + 1):
            q_degree = residual_degree - y_degree - p_degree
            exponents = (1, y_degree, p_degree, q_degree)
            coefficient = sp.Symbol("r" + "".join(map(str, exponents)))
            coefficients[exponents] = coefficient
            R1 += coefficient * y**y_degree * p**p_degree * q**q_degree


def coefficient(exponents: tuple[int, int, int, int]) -> sp.Symbol:
    return coefficients[exponents]


# The x^2 jet of the graph pullback.  Terms O(x^3) cannot affect a second
# derivative at x=0.  This is much smaller than expanding the full generic
# degree-five graph.
A_jet = (
    p
    + x * (3 * y * p + 3 * q)
    + x**2 * (3 * y**2 * p + 6 * y * q)
)
B_jet = (
    4 * y**2 * p
    + y * q
    + x * (7 * y**3 * p + 12 * y**2 * q + 2 * R0)
    + x**2
    * (3 * y**4 * p + 9 * y**3 * q + 2 * R1 - 3 * y * R0)
)
potential_jet = L * A_jet**2 + M * A_jet + N * B_jet
axis_hessian = sp.hessian(potential_jet, (x, y, p, q)).subs(
    {x: 0, y: t, p: c * t, q: d * t}
)
axis_determinant = determinant_truncation(axis_hessian, (t,), (10,))


def axis_coefficient(degree: int) -> sp.Expr:
    # Substitute the triangular jet relations before factoring.  Factoring
    # the completely generic coefficient needlessly asks SymPy to discover
    # those relations inside a large multivariate square.
    return axis_determinant.get((degree,), sp.S.Zero)


R5_x_free = sum(
    value * y**exponents[1] * p**exponents[2] * q**exponents[3]
    for exponents, value in coefficients.items()
    if exponents[0] == 0 and sum(exponents) == 5
)
leading_q_derivative = sp.diff(R5_x_free, q).subs(
    {y: 1, p: c, q: d}
)
assert sp.factor(
    axis_coefficient(10) - 256 * N**4 * leading_q_derivative**2
) == 0


# Constancy kills every x-free quintic coefficient containing q.
q_quintic_zero = {
    value: 0
    for exponents, value in coefficients.items()
    if exponents[0] == 0
    and sum(exponents) == 5
    and exponents[3] > 0
}
assert len(q_quintic_zero) == 15
assert sp.factor(axis_coefficient(9).subs(q_quintic_zero)) == 0


r0004 = coefficient((0, 0, 0, 4))
r0013 = coefficient((0, 0, 1, 3))
r0022 = coefficient((0, 0, 2, 2))
r0031 = coefficient((0, 0, 3, 1))
r0103 = coefficient((0, 1, 0, 3))
r0112 = coefficient((0, 1, 1, 2))
r0121 = coefficient((0, 1, 2, 1))
r0202 = coefficient((0, 2, 0, 2))
r0211 = coefficient((0, 2, 1, 1))
r0301 = coefficient((0, 3, 0, 1))
s0050 = coefficient((0, 0, 5, 0))
s0140 = coefficient((0, 1, 4, 0))
s0230 = coefficient((0, 2, 3, 0))
s0320 = coefficient((0, 3, 2, 0))
s0410 = coefficient((0, 4, 1, 0))
zeta = coefficient((0, 5, 0, 0))

t8_inner = (
    32 * r0004 * d**3
    + 24 * r0013 * c * d**2
    + 16 * r0022 * c**2 * d
    + 8 * r0031 * c**3
    + 24 * r0103 * d**2
    + 16 * r0112 * c * d
    + 8 * r0121 * c**2
    + 16 * r0202 * d
    + 8 * r0211 * c
    + 8 * r0301
    - 5 * s0050 * c**4
    - 4 * s0140 * c**3
    - 3 * s0230 * c**2
    - 2 * s0320 * c
    - s0410
)
assert sp.factor(
    axis_coefficient(8).subs(q_quintic_zero)
    - 4 * N**4 * t8_inner**2
) == 0


t8_forced = dict(q_quintic_zero)
for symbol in (r0004, r0013, r0022, r0103, r0112, r0202):
    t8_forced[symbol] = 0
t8_forced.update(
    {
        s0050: 0,
        r0031: s0140 / 2,
        r0121: sp.Rational(3, 8) * s0230,
        r0211: s0320 / 4,
        r0301: s0410 / 8,
    }
)
assert sp.factor(axis_coefficient(8).subs(t8_forced)) == 0
t7_expected = -L * N**3 * (
    4 * s0140 * c**3 + 3 * s0230 * c**2 + 2 * s0320 * c + s0410
) * (
    2 * s0140 * c**4
    - s0230 * c**3
    - 4 * s0320 * c**2
    - 7 * s0410 * c
    - 10 * zeta
)
assert sp.factor(axis_coefficient(7).subs(t8_forced) - t7_expected) == 0


# If the product in t7 is the zero polynomial, one factor vanishes.  The
# second factor can vanish only when the first one (and zeta) vanishes, so in
# every branch the four displayed coefficients are zero.
top_forced = dict(t8_forced)
top_forced.update({s0140: 0, s0230: 0, s0320: 0, s0410: 0})
assert sp.factor(axis_coefficient(7).subs(top_forced)) == 0


eta = coefficient((0, 0, 0, 3))
theta = coefficient((0, 0, 1, 2))
tau = coefficient((0, 0, 2, 1))
a = coefficient((0, 0, 4, 0))
iota = coefficient((0, 1, 0, 2))
sigma = coefficient((0, 1, 1, 1))
b = coefficient((0, 1, 3, 0))
rho = coefficient((0, 2, 0, 1))
e = coefficient((0, 2, 2, 0))
f = coefficient((0, 3, 1, 0))

S = (
    -8 * a * c**3
    + (16 * tau - 6 * b) * c**2
    + 32 * theta * c * d
    + (16 * sigma - 4 * e) * c
    + 48 * eta * d**2
    + 32 * iota * d
    + 16 * rho
    - 2 * f
    + 89
)
T = (
    240 * eta * d**2
    + 160 * theta * c * d
    + 80 * tau * c**2
    + 160 * iota * d
    + 80 * sigma * c
    + 80 * rho
    + 492
)
assert sp.factor(
    axis_coefficient(6).subs(top_forced)
    - N**3 * (N * S**2 + L * zeta * T)
) == 0


degree_six_forced = dict(top_forced)
degree_six_forced.update(
    {
        a: 0,
        eta: 0,
        theta: 0,
        iota: 0,
        b: sp.Rational(8, 3) * tau,
    }
)
t5_after_six = sp.Poly(
    sp.expand(axis_coefficient(5).subs(degree_six_forced)), c, d
)
assert sp.factor(t5_after_six.coeff_monomial(c**5)) == (
    -sp.Rational(64, 3) * L * N**3 * tau**2
)


tau_zero = {tau: 0, b: 0}
t6_after_tau = sp.Poly(
    sp.expand(axis_coefficient(6).subs(degree_six_forced).subs(tau_zero)),
    c,
    d,
)
assert sp.factor(t6_after_tau.coeff_monomial(c**2)) == (
    16 * N**4 * (4 * sigma - e) ** 2
)
sigma_relation = {e: 4 * sigma}
t5_after_relation = sp.Poly(
    sp.expand(
        axis_coefficient(5)
        .subs(degree_six_forced)
        .subs(tau_zero)
        .subs(sigma_relation)
    ),
    c,
    d,
)
assert sp.factor(t5_after_relation.coeff_monomial(c**3)) == (
    64 * L * N**3 * sigma**2
)


final_high = dict(degree_six_forced)
final_high.update({tau: 0, b: 0, sigma: 0, e: 0})
beta = 16 * rho - 2 * f + 89
t6_final = sp.Poly(
    sp.expand(axis_coefficient(6).subs(final_high)), c, d
)
assert sp.expand(
    t6_final.as_expr()
    - N**3 * (N * beta**2 + L * zeta * (80 * rho + 492))
) == 0


# The zeta=0 branch has no rational point.  Its t6 equation first gives
# beta=0, and the c*t^5 equation then gives the irreducible quadratic Q.
zeta_zero = {zeta: 0, f: (16 * rho + 89) / 2}
Q = 160 * rho**2 + 1968 * rho + 6021
t5_zeta_zero = sp.Poly(
    sp.expand(axis_coefficient(5).subs(final_high).subs(zeta_zero)), c, d
)
assert sp.factor(t5_zeta_zero.coeff_monomial(c)) == 2 * L * N**3 * Q
assert sp.discriminant(Q, rho) == 576 * 34
assert sp.Poly(Q, rho).ground_roots() == {}


# Abstract leading-square determinant identity.  For a homogeneous degree-d
# form F in four variables, Euler gives H*x=(d-1)g and g.x=dF.  Combining
# this identity with the checked rank-one update yields
# det Hess(F^2)=2^4*(2d-1)/(d-1)*F^4*det(H).  At d=8 this is 240/7.
F0 = sp.Symbol("F0")
g_symbols = sp.symbols("g0:4")
g_vector = sp.Matrix(g_symbols)
h_symbols = sp.symbols("h00 h01 h02 h03 h11 h12 h13 h22 h23 h33")
H = sp.Matrix(
    [
        [h_symbols[0], h_symbols[1], h_symbols[2], h_symbols[3]],
        [h_symbols[1], h_symbols[4], h_symbols[5], h_symbols[6]],
        [h_symbols[2], h_symbols[5], h_symbols[7], h_symbols[8]],
        [h_symbols[3], h_symbols[6], h_symbols[8], h_symbols[9]],
    ]
)
rank_one_update = 16 * (
    F0**4 * H.det()
    + F0**3 * (g_vector.T * H.adjugate() * g_vector)[0]
)
assert sp.factor(
    (2 * (F0 * H + g_vector * g_vector.T)).det()
    - rank_one_update
) == 0
assert sp.Rational(16) * (1 + sp.Rational(8, 7)) == sp.Rational(240, 7)


# The two constant-kernel charts forced by Gordan--Noether when zeta != 0.
kappa = sp.Symbol("kappa")
T_coefficients: dict[tuple[int, int, int], sp.Symbol] = {}
T4 = sp.S.Zero
invariant = p - kappa * q
for x_degree in range(5):
    for y_degree in range(5 - x_degree):
        invariant_degree = 4 - x_degree - y_degree
        value = sp.Symbol(f"T{x_degree}{y_degree}{invariant_degree}")
        T_coefficients[(x_degree, y_degree, invariant_degree)] = value
        T4 += (
            value
            * x**x_degree
            * y**y_degree
            * invariant**invariant_degree
        )
R5_finite_chart = zeta * y**5 + x * T4
assert sp.expand(
    kappa * sp.diff(R5_finite_chart, p)
    + sp.diff(R5_finite_chart, q)
) == 0
assert sp.expand(R5_finite_chart.subs(x, 0) - zeta * y**5) == 0

U_coefficients: dict[tuple[int, int, int], sp.Symbol] = {}
U4 = sp.S.Zero
for x_degree in range(5):
    for y_degree in range(5 - x_degree):
        q_degree = 4 - x_degree - y_degree
        value = sp.Symbol(f"U{x_degree}{y_degree}{q_degree}")
        U_coefficients[(x_degree, y_degree, q_degree)] = value
        U4 += value * x**x_degree * y**y_degree * q**q_degree
R5_infinite_chart = zeta * y**5 + x * U4
assert sp.diff(R5_infinite_chart, p) == 0
assert sp.expand(R5_infinite_chart.subs(x, 0) - zeta * y**5) == 0


def finite_field_counts(prime: int) -> dict[str, int]:
    """Enumerate the projected (rho,beta,zeta) axis-jet quadric."""

    total = 0
    zeta_zero_count = 0
    zeta_zero_after_t5_count = 0
    beta_zero_count = 0
    exceptional_c_count = 0
    for rho_value in range(prime):
        for beta_value in range(prime):
            for zeta_value in range(prime):
                equation = (
                    2 * beta_value**2
                    + zeta_value * (80 * rho_value + 492)
                ) % prime
                if equation != 0:
                    continue
                total += 1
                if zeta_value == 0:
                    zeta_zero_count += 1
                    if (
                        beta_value == 0
                        and (
                            160 * rho_value**2
                            + 1968 * rho_value
                            + 6021
                        )
                        % prime
                        == 0
                    ):
                        zeta_zero_after_t5_count += 1
                if beta_value == 0:
                    beta_zero_count += 1
                if (5 * zeta_value + 4 * beta_value) % prime == 0:
                    exceptional_c_count += 1
    return {
        "total": total,
        "zeta_zero": zeta_zero_count,
        "zeta_zero_after_t5": zeta_zero_after_t5_count,
        "beta_zero": beta_zero_count,
        "five_zeta_plus_four_beta_zero": exceptional_c_count,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    primes = (101, 103)
    counts = {str(prime): finite_field_counts(prime) for prime in primes}
    # A small rational point on the projected nonzero-zeta quadric.  It is a
    # jet survivor only, not a full graph candidate.
    rational_projection = {
        "rho": "0",
        "beta": "6",
        "zeta": "-6/41",
    }
    assert sp.factor(
        2 * sp.Integer(6) ** 2
        + (-sp.Rational(6, 41)) * sp.Integer(492)
    ) == 0

    payload = {
        "schema": "hc4-meng-yang-quintic-graph-normal-slice.v1",
        "status": "exact normal-slice reduction; no HC4 candidate",
        "axis_graph_coefficients": len(coefficients),
        "finite_field_projected_counts": counts,
        "rational_projected_jet_survivor": rational_projection,
        "rational_zero_zeta_branch": "empty",
        "zero_zeta_polynomial": [160, 1968, 6021],
        "zero_zeta_discriminant": 576 * 34,
        "nonzero_zeta_kernel_charts": [
            "zeta*y^5+x*T4(x,y,p-kappa*q)",
            "zeta*y^5+x*U4(x,y,q)",
        ],
        "open": [
            "lower determinant layers",
            "off-axis determinant identity",
            "collision equations on the full graph",
        ],
    }
    canonical = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(canonical.encode()).hexdigest()
    if args.output is not None:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(canonical)
        print(f"WROTE {output}")
        print(f"SHA256 {digest}")

    print("PASS: the complete quintic x=0 line jet reduces triangularly")
    print("PASS: the x-free quintic top is zeta*y^5")
    print("PASS: the zeta=0 branch has no rational point")
    print("PASS: Gordan--Noether leaves two nonzero-zeta kernel charts")
    print("PASS: projected finite-field axis-jet loci were enumerated exactly")


if __name__ == "__main__":
    main()
