#!/usr/bin/env python3
"""Verify the exact all-degree HC4 frontend identities.

This checker supports three written arguments.

* HC4FSD1: for m>=2, the diagonal ternary reverse-Schur equation with
  Hess(f)=diag(x^m,y^m,z^m) has only
  a=alpha*x^m+beta*y^m+gamma*z^m.
* HC4FSD2: for m>=3, the normalized minimal diagonal three-layer tower has
  no constant rank-two quadratic determinant-preserving direction on any
  nonzero channel support.
* HC4FSD3: for every nonzero channel vector, after arbitrary compatible lower
  homogeneous layers are restored, no constant determinant-preserving
  direction can have rank at least two.  In degree five, HC4MR4 then closes
  rank one as well.
* HC4MYGJ2: every normal order k>=1 of a Meng--Yang graph has multiplier
  -4*L*N^3*k*(k+1).

The universal quantifiers are discharged in the accompanying written proofs.
This script checks their polynomial identities, parity inequalities, support
injections, and representative coefficient ideals exactly over QQ.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BOUNDED_ARTIFACT = (
    ROOT / "artifacts" / "generated-results" / "hc4_all_degree_frontend_experiments.json"
)
BOUNDED_SOURCE = ROOT / "scripts" / "research_hc4_all_degree_frontends.py"
EXPECTED_BOUNDED_ARTIFACT_SHA256 = (
    "7176670a76f152f23c6f9b56264c39e2345d91641fdfe275130b9b32e99ffedb"
)
EXPECTED_BOUNDED_SOURCE_SHA256 = (
    "0b6e0a1272d1201645dce96d8a039fde3c9d4c883476f40d1fed57dac446333f"
)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_existing() -> None:
    assert file_sha256(BOUNDED_ARTIFACT) == EXPECTED_BOUNDED_ARTIFACT_SHA256
    assert file_sha256(BOUNDED_SOURCE) == EXPECTED_BOUNDED_SOURCE_SHA256
    payload = json.loads(BOUNDED_ARTIFACT.read_text(encoding="utf-8"))
    assert payload["schema_version"] == 1
    assert payload["status"] == "exact_bounded_regression"
    assert payload["source_sha256"] == EXPECTED_BOUNDED_SOURCE_SHA256

    schur_rows = payload["fermat_ternary_schur"]["results"]
    assert [row["potential_degree"] for row in schur_rows] == list(range(4, 9))
    assert all(row["radical_equals_pure_power_locus"] for row in schur_rows)

    rank_two_rows = payload["diagonal_rank_two_recognition"]["results"]
    assert [
        (row["potential_degree"], row["active_channel_count"])
        for row in rank_two_rows
    ] == [(degree, channels) for degree in range(5, 9) for channels in range(1, 4)]
    assert all(row["rank_two_scheme_empty"] for row in rank_two_rows)

    normal_rows = payload["meng_yang_normal_symbol"]["verified_orders"]
    assert [row["graph_jet_order"] for row in normal_rows] == list(range(1, 13))
    warning = payload["scope_warning"]
    assert "finite tables are regressions" in warning
    assert "minimal tower" in warning
    assert "formal solvability rather than polynomial termination" in warning
    print(
        "PASS: committed bounded HC4 all-degree frontend regression is intact "
        "and correctly scoped; no symbolic or Singular replay"
    )


parser = argparse.ArgumentParser()
parser.add_argument(
    "--audit-existing-only",
    action="store_true",
    help="validate the committed bounded regression without mathematical replay",
)
arguments = parser.parse_args()
if arguments.audit_existing_only:
    audit_existing()
    raise SystemExit(0)

import sympy as sp


def homogeneous_exponents(degree: int) -> tuple[tuple[int, int, int], ...]:
    return tuple(
        (first, second, degree - first - second)
        for first in range(degree + 1)
        for second in range(degree - first + 1)
    )


def x_valuation(expression: sp.Expr, variable: sp.Symbol) -> int:
    polynomial = sp.Poly(sp.expand(expression), variable)
    if polynomial.is_zero:
        return 10**9
    return min(monomial[0] for monomial, _ in polynomial.terms())


# ---------------------------------------------------------------------------
# HC4FSD1: uniform diagonal ternary-Schur rigidity.
# ---------------------------------------------------------------------------

# If x^m divides a_x^2, then x^ceil(m/2) divides a_x.  Hence a monomial of a
# with positive x exponent has x exponent at least ceil(m/2)+1.  The same
# holds for y and z.  Two positive exponents would exceed total degree m.
parity_parameter = sp.symbols("parity_parameter", integer=True, nonnegative=True)
assert sp.expand(2 * (parity_parameter + 1) - 2 * parity_parameter) == 2
assert (
    sp.expand(2 * (parity_parameter + 2) - (2 * parity_parameter + 1))
    == 3
)

# Exact support regression over a deliberately wider range than the discovery
# sweep.  The written parity calculation above is the degree-free step.
for m in range(2, 65):
    threshold = (m + 1) // 2
    admissible = []
    for exponent in homogeneous_exponents(m):
        derivative_divisibility = all(
            coordinate == 0 or coordinate - 1 >= threshold
            for coordinate in exponent
        )
        if derivative_divisibility:
            admissible.append(exponent)
    assert set(admissible) == {(m, 0, 0), (0, m, 0), (0, 0, m)}

# Check the displayed converse and quotient formula exactly in several degrees.
x, y, z = sp.symbols("x y z")
alpha, beta, gamma = sp.symbols("alpha beta gamma")
for m in range(2, 13):
    a_form = alpha * x**m + beta * y**m + gamma * z**m
    determinant = x**m * y**m * z**m
    numerator = sp.expand(
        sp.diff(a_form, x) ** 2 * y**m * z**m
        + sp.diff(a_form, y) ** 2 * x**m * z**m
        + sp.diff(a_form, z) ** 2 * x**m * y**m
    )
    quotient = m**2 * (
        alpha**2 * x ** (m - 2)
        + beta**2 * y ** (m - 2)
        + gamma**2 * z ** (m - 2)
    )
    assert sp.expand(numerator - determinant * quotient) == 0


# ---------------------------------------------------------------------------
# HC4FSD2: rank-two obstruction on the minimal diagonal tower.
# ---------------------------------------------------------------------------

q11, q12, q13, q22, q23, q33 = q_variables = sp.symbols(
    "q11 q12 q13 q22 q23 q33"
)
q_matrix = sp.Matrix(
    ((q11, q12, q13), (q12, q22, q23), (q13, q23, q33))
)
pencil = sp.symbols("pencil")

# A general rank-two constant symmetric 4-by-4 direction first loses its
# fourth row and column.  The coefficient of pencil in the top determinant
# is q44*x^m*y^m*z^m.  After q44=0, the three degree-2m coefficients of
# pencil^2 are -q14^2, -q24^2, and -q34^2 times distinct monomials.
q14, q24, q34, q44 = sp.symbols("q14 q24 q34 q44")
q_matrix4 = sp.Matrix(
    (
        (q11, q12, q13, q14),
        (q12, q22, q23, q24),
        (q13, q23, q33, q34),
        (q14, q24, q34, q44),
    )
)
top_hessian4 = sp.diag(x**3, y**3, z**3, 0)
top_pencil = sp.Poly(
    sp.expand((top_hessian4 + pencil * q_matrix4).det()),
    pencil,
    x,
    y,
    z,
)
assert top_pencil.coeff_monomial(pencil * x**3 * y**3 * z**3) == q44
assert sp.expand(
    top_pencil.coeff_monomial(pencil**2 * x**3 * y**3).subs(q44, 0)
    + q34**2
) == 0
assert sp.expand(
    top_pencil.coeff_monomial(pencil**2 * x**3 * z**3).subs(q44, 0)
    + q24**2
) == 0
assert sp.expand(
    top_pencil.coeff_monomial(pencil**2 * y**3 * z**3).subs(q44, 0)
    + q14**2
) == 0


def canonical_polynomial(expression: sp.Expr) -> sp.Expr:
    polynomial = sp.Poly(sp.expand(expression), *q_variables)
    return sp.factor(polynomial.monic().as_expr())


def top_face_equations(channels: tuple[int, int, int]) -> set[sp.Expr]:
    """Return the m=3 representative of the degree-free top support."""

    m = 3
    diagonal = sp.diag(x**m, y**m, z**m)
    direction = sp.Matrix(
        tuple(
            channel * variable ** (m - 1)
            for channel, variable in zip(channels, (x, y, z), strict=True)
        )
    )
    corner = sum(
        channel**2 * variable ** (m - 2)
        for channel, variable in zip(channels, (x, y, z), strict=True)
    )
    shifted = diagonal + pencil * q_matrix
    difference = sp.expand(
        shifted.det() * corner
        - (direction.T * shifted.adjugate() * direction)[0]
    )
    equations = {
        canonical_polynomial(coefficient)
        for monomial, coefficient in sp.Poly(
            difference, pencil, x, y, z
        ).terms()
        if monomial[0] > 0 and coefficient != 0
    }
    return equations


det_q = canonical_polynomial(q_matrix.det())
one_channel = top_face_equations((1, 0, 0))
assert {
    canonical_polynomial(q11),
    canonical_polynomial(q11 * q22 - q12**2),
    canonical_polynomial(q11 * q33 - q13**2),
    det_q,
}.issubset(one_channel)

two_channels = top_face_equations((1, 1, 0))
assert {
    canonical_polynomial(q11),
    canonical_polynomial(q22),
    canonical_polynomial(q12),
    canonical_polynomial(q11 * q33 - q13**2),
    canonical_polynomial(q22 * q33 - q23**2),
    det_q,
}.issubset(two_channels)

three_channels = top_face_equations((1, 1, 1))
assert {
    canonical_polynomial(variable) for variable in q_variables
}.issubset(three_channels)

# On the one-channel stratum the preceding equations force the first row and
# column of Q to vanish.  Rank two then means that the inactive 2-by-2 block
# has nonzero determinant.  Its complementary Hessian minor on (x,t) is the
# following degree-free expression.  Its t^2 coefficient is nonzero for
# every m>=3 in characteristic zero.
m_symbol = sp.symbols("m", integer=True, positive=True)
linear_coefficient = sp.expand((m_symbol - 1) - 2 * (m_symbol - 2))
quadratic_coefficient = sp.expand(
    (m_symbol - 2) * (m_symbol - 3) / 2 - (m_symbol - 2) ** 2
)
assert linear_coefficient == 3 - m_symbol
assert sp.simplify(
    quadratic_coefficient + (m_symbol - 1) * (m_symbol - 2) / 2
) == 0
for m in range(3, 65):
    assert quadratic_coefficient.subs(m_symbol, m) != 0
    t = sp.symbols("t")
    h_xx = (
        x**m
        + (m - 1) * t * x ** (m - 2)
        + sp.Rational(1, 2)
        * (m - 2)
        * (m - 3)
        * t**2
        * x ** (m - 4)
    )
    h_xt = x ** (m - 1) + (m - 2) * t * x ** (m - 3)
    h_tt = x ** (m - 2)
    complementary_minor = sp.expand(h_xx * h_tt - h_xt**2)
    expected_minor = sp.expand(
        (3 - m) * t * x ** (2 * m - 4)
        - sp.Rational(1, 2)
        * (m - 1)
        * (m - 2)
        * t**2
        * x ** (2 * m - 6)
    )
    assert complementary_minor == expected_minor

# The scope limitation is real.  At m=4 the permitted next homogeneous term
# t^3/6 adds t to H_tt.  It cancels both displayed minimal-tower faces and
# leaves t^3, so those faces are not invariant under arbitrary lower layers.
t = sp.symbols("t")
augmented_m4_minor = sp.expand(
    (x**4 + 3 * t * x**2 + t**2) * (x**2 + t)
    - (x**3 + 2 * t * x) ** 2
)
assert augmented_m4_minor == t**3


# ---------------------------------------------------------------------------
# HC4FSD3: arbitrary-lower-layer propagation.
# ---------------------------------------------------------------------------

# Once the leading faces reduce a rank-two direction to a nondegenerate
# inactive 2-by-2 block, its exact s^2 coefficient is the inactive determinant
# times the complementary binary Hessian determinant.  This identity does not
# see any of the lower layers or any of the cross entries of the full Hessian.
hxx, hxy, hxz, hxt = sp.symbols("hxx hxy hxz hxt")
hyy, hyz, hyt = sp.symbols("hyy hyz hyt")
hzz, hzt, htt = sp.symbols("hzz hzt htt")
generic_hessian = sp.Matrix(
    (
        (hxx, hxy, hxz, hxt),
        (hxy, hyy, hyz, hyt),
        (hxz, hyz, hzz, hzt),
        (hxt, hyt, hzt, htt),
    )
)
inactive_block = sp.Matrix(((q22, q23), (q23, q33)))
inactive_direction = sp.zeros(4)
inactive_direction[1:3, 1:3] = inactive_block
inactive_pencil = sp.Poly(
    sp.expand((generic_hessian + pencil * inactive_direction).det()), pencil
)
assert sp.factor(
    inactive_pencil.coeff_monomial(pencil**2)
    - inactive_block.det() * (hxx * htt - hxt**2)
) == 0

# A rank-three direction has already lost its t row and column.  Its exact
# s^3 coefficient is det(Q_3)*H_tt, so the nonzero leading Schur coefficient
# of H_tt excludes rank three independently of all later layers.
ternary_direction = sp.zeros(4)
ternary_direction[:3, :3] = q_matrix
ternary_pencil = sp.Poly(
    sp.expand((generic_hessian + pencil * ternary_direction).det()), pencil
)
assert sp.factor(
    ternary_pencil.coeff_monomial(pencil**3) - q_matrix.det() * htt
) == 0

# If the complementary binary Hessian vanished, p=phi_x and q=phi_t would
# have zero Jacobian and hence a common polynomial generator.  Their degrees
# are consecutive, D-1 and D-2, so degree multiplicativity forces that
# generator to have degree one.  The written proof then uses the pure x^(D-1)
# leading form of p and equality p_t=q_x to obtain the contradiction.
for potential_degree in range(5, 65):
    assert sp.gcd(potential_degree - 1, potential_degree - 2) == 1

# Exact bounded regressions retain every binary term in the unrestricted
# homogeneous layers.  These computations are checks of the degree-free
# common-generator proof, not its logical basis.
binary_x, binary_t = sp.symbols("binary_x binary_t")
for potential_degree in range(5, 10):
    binary_form = (
        sp.Rational(1, potential_degree * (potential_degree - 1))
        * binary_x**potential_degree
        + sp.Rational(1, potential_degree - 2)
        * binary_t
        * binary_x ** (potential_degree - 2)
        + sp.Rational(1, 2)
        * binary_t**2
        * binary_x ** (potential_degree - 4)
    )
    lower_coefficients = []
    # The unrestricted ternary pieces in degrees D-1 and D-2 restrict to
    # these three additional binary monomials.
    for label, monomial in (
        ("b", binary_x ** (potential_degree - 1)),
        ("c", binary_t * binary_x ** (potential_degree - 3)),
        ("e", binary_x ** (potential_degree - 2)),
    ):
        coefficient = sp.symbols(f"{label}_{potential_degree}")
        lower_coefficients.append(coefficient)
        binary_form += coefficient * monomial
    for total_degree in range(2, potential_degree - 2):
        for t_exponent in range(total_degree + 1):
            coefficient = sp.symbols(
                f"lower_{potential_degree}_{total_degree}_{t_exponent}"
            )
            lower_coefficients.append(coefficient)
            binary_form += (
                coefficient
                * binary_x ** (total_degree - t_exponent)
                * binary_t**t_exponent
            )
    binary_hessian = sp.expand(
        sp.diff(binary_form, binary_x, 2)
        * sp.diff(binary_form, binary_t, 2)
        - sp.diff(binary_form, binary_x, binary_t) ** 2
    )
    binary_equations = sp.Poly(binary_hessian, binary_x, binary_t).coeffs()
    binary_ideal = sp.groebner(
        binary_equations, *lower_coefficients, order="grevlex"
    )
    assert binary_ideal.contains(sp.Integer(1))


# ---------------------------------------------------------------------------
# HC4MYGJ2: exact all-normal multiplier for Meng--Yang graphs.
# ---------------------------------------------------------------------------

x, y, p, q, r = sp.symbols("x y p q r")
L, M, N = sp.symbols("L M N")
u = 1 + x * y
q0 = y**2 * (1 + 3 * u)
meng_a = u**3 * p + 3 * x * u**2 * q - x**3 * r
meng_b = (
    u * q0 * p
    + (y + 3 * x * q0) * q
    + x * (5 - 3 * u) * r
)
potential = sp.expand(L * meng_a**2 + M * meng_a + N * meng_b)

# The potential is exactly quadratic in r with second derivative 2*L*x^6.
# Its omitted gradient has universal first x-jet 2*N*x, independently of R.
potential_r = sp.diff(potential, r)
assert sp.factor(sp.diff(potential, r, 2) - 2 * L * x**6) == 0
assert sp.factor(potential_r.subs(x, 0)) == 0
assert sp.factor(sp.diff(potential_r, x).subs(x, 0) - 2 * N) == 0
assert x_valuation(sp.expand(potential_r - 2 * N * x), x) >= 2

# Exact Taylor identity in the omitted coordinate; this avoids any analytic
# assumption in the x-adic argument.
delta = sp.symbols("delta")
assert sp.expand(
    potential.subs(r, r + delta)
    - potential
    - potential_r * delta
    - L * x**6 * delta**2
) == 0

# The graph Hessian's complementary cofactor at x=0 is independent of every
# trace and normal jet.  This is the determinant factor multiplying the first
# changed (x,x) entry.
plane_trace, plane_normal = sp.symbols("plane_trace plane_normal")
trace_y, trace_p, trace_q = sp.symbols("trace_y trace_p trace_q")
ambient_hessian = sp.hessian(potential, (x, y, p, q, r))
ambient_plane = ambient_hessian.subs({x: 0, r: plane_trace})
source_unit = sp.Matrix((1, 0, 0, 0))
trace_gradient = sp.Matrix((plane_normal, trace_y, trace_p, trace_q))
graph_plane_hessian = (
    ambient_plane[:4, :4]
    + 2
    * N
    * (source_unit * trace_gradient.T + trace_gradient * source_unit.T)
)
complementary_cofactor = sp.factor(graph_plane_hessian[1:, 1:].det())
assert complementary_cofactor == -2 * L * N**2

# For symbolic k, twice differentiating the leading potential difference
# 2*N*x^(k+1)*U gives 2*N*k*(k+1)*x^(k-1)*U.  All other Hessian entries have
# x-order at least k.  Multiplication by the cofactor gives HC4MYGJ2.
k = sp.symbols("k", integer=True, positive=True)
U = sp.symbols("U")
leading_xx = sp.factor(sp.diff(2 * N * x ** (k + 1) * U, x, 2))
assert sp.factor(leading_xx - 2 * N * k * (k + 1) * x ** (k - 1) * U) == 0
normal_multiplier = sp.factor(
    complementary_cofactor * 2 * N * k * (k + 1)
)
assert normal_multiplier == -4 * L * N**3 * k * (k + 1)

print("PASS HC4FSD1: diagonal ternary Schur rigidity holds for every m>=2")
print("PASS HC4FSD2: minimal diagonal rank-two obstruction holds for every m>=3")
print("PASS HC4FSD3: arbitrary lower layers do not restore rank >= 2")
print("PASS HC4MYGJ2: Meng--Yang normal recursion has a unit at every order")
print("COROLLARY: with HC4MR4, the diagonal quintic has no nonzero constant Q")
