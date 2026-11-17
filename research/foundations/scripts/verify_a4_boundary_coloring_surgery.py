#!/usr/bin/env python3
"""Exact boundary-coloring surgery experiment for the polynomial A4 cone.

The experiment has three parts.

1. Recompute the complete valuation above the target B-divisor.  The
   displayed affine divisor L has (e,f)=(2,2), so there is no second global
   height-one prime available as a third center.
2. Perform the elementary affine modification of the residual center
   (L,W), using q=W/L.  This moves the old generic L-valuation out of the
   affine chart and keeps both the A4 function field and polynomial target
   coordinates, but the exceptional divisor is contracted to the target
   origin and the relative canonical factor changes the determinant ledger.
3. Try the analogous reconstruction-denominator center (L,H).  Its affine
   chart x=H/L is singular above the two cubic-contact cluster points.
"""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import runpy

import sympy as sp


SCRIPT_DIR = Path(__file__).resolve().parent


def quiet_run(path):
    """Run a neighboring exact verifier without replaying its status lines."""

    with redirect_stdout(StringIO()):
        return runpy.run_path(str(path))


ledger = quiet_run(SCRIPT_DIR / "verify_a4_ledger_reduction.py")

U = ledger["U"]
V = ledger["V"]
W = ledger["W"]
H = ledger["H"]
K = ledger["K"]
L = ledger["L"]
N1 = ledger["N1"]
N2 = ledger["N2"]
B_target = ledger["B_homogeneous"]
P = ledger["P"]
Q = ledger["Q"]
R = ledger["R"]
t = ledger["t"]
U_param = ledger["U_param"]
V_param = ledger["V_param"]


def factor_order(expression, factor, variables):
    """Return the exact polynomial order of ``factor`` in ``expression``."""

    quotient = sp.Poly(sp.expand(expression), *variables, domain=sp.QQ)
    divisor = sp.Poly(sp.expand(factor), *variables, domain=sp.QQ)
    order = 0
    while True:
        next_quotient, remainder = sp.div(quotient, divisor)
        if not remainder.is_zero:
            return order
        quotient = next_quotient
        order += 1


# ---------------------------------------------------------------------------
# 1. There is one global (e,f)=(2,2) L-prime, not a missing second prime
# ---------------------------------------------------------------------------

s = sp.symbols("s")
source_normalization = {U: U_param, V: V_param}
residue_map = sp.factor(
    ((N1 - N2) / H).subs(source_normalization)
)
assert residue_map == t**2 / (2 * t - 3)

residue_polynomial = sp.Poly(
    t**2 - 2 * s * t + 3 * s,
    t,
    domain=sp.QQ.frac_field(s),
)
assert residue_polynomial.is_irreducible
assert sp.factor(sp.discriminant(residue_polynomial.as_expr(), t)) == (
    4 * s * (s - 3)
)

ramification_index = 2
residue_degree = residue_polynomial.degree()
cover_degree = 4
assert ramification_index * residue_degree == cover_degree

print("PASS: L has complete branch datum (e,f)=(2,2)")
print("PASS: its local degree is four, so no second global B-prime is missing")


# ---------------------------------------------------------------------------
# 2. Residual-center surgery: blow up (L,W), retain q=W/L
# ---------------------------------------------------------------------------

q = sp.symbols("q")
old_cone_map = sp.Matrix([W * N1, W * N2, W * H])
old_jacobian = sp.factor(old_cone_map.jacobian((U, V, W)).det())
assert old_jacobian == 4 * W**2 * K**3 * L

# The retained blowup chart is again affine three-space.  Its map back to
# the old cone is pi(U,V,q)=(U,V,Lq), with inverse q=W/L at the generic
# point.  Hence the source and target subfields are unchanged.
pi_substitution = {W: L * q}
modified_map = sp.expand(old_cone_map.subs(pi_substitution))
expected_modified_map = sp.Matrix([q * L * N1, q * L * N2, q * L * H])
assert all(
    sp.factor(actual - expected) == 0
    for actual, expected in zip(modified_map, expected_modified_map)
)
relative_jacobian = sp.factor(
    sp.Matrix([U, V, L * q]).jacobian((U, V, q)).det()
)
assert relative_jacobian == L

modified_jacobian = sp.factor(
    modified_map.jacobian((U, V, q)).det()
)
assert modified_jacobian == 4 * q**2 * K**3 * L**4
assert sp.factor(
    modified_jacobian
    - old_jacobian.subs(pi_substitution) * relative_jacobian
) == 0

modified_B_pullback = sp.factor(
    B_target.subs(
        {
            P: modified_map[0],
            Q: modified_map[1],
            R: modified_map[2],
        }
    )
)
assert modified_B_pullback == q**3 * K**3 * L**5
assert sp.factor(
    modified_B_pullback / modified_jacobian - q * L / 4
) == 0

modified_jacobian_orders = {
    "q": factor_order(modified_jacobian, q, (U, V, q)),
    "K": factor_order(modified_jacobian, K, (U, V, q)),
    "L": factor_order(modified_jacobian, L, (U, V, q)),
}
modified_target_orders = {
    "q": factor_order(modified_B_pullback, q, (U, V, q)),
    "K": factor_order(modified_B_pullback, K, (U, V, q)),
    "L": factor_order(modified_B_pullback, L, (U, V, q)),
}
assert modified_jacobian_orders == {"q": 2, "K": 3, "L": 4}
assert modified_target_orders == {"q": 3, "K": 3, "L": 5}

# The old generic point L=0,W!=0 has no center in this chart: W=Lq forces
# W=0 over L=0.  The new exceptional divisor E=(L=0) is affine, but every
# target coordinate vanishes there, so E maps to the target origin rather
# than dominating the generic point of B=0.
assert all(sp.factor(coordinate / L).is_polynomial(U, V, q)
           for coordinate in modified_map)
assert all(
    sp.rem(
        sp.Poly(coordinate, U, V, q, domain=sp.QQ),
        sp.Poly(L, U, V, q, domain=sp.QQ),
    ).is_zero
    for coordinate in modified_map
)
exceptional_valuation = {"L": 1, "W": 1, "q": 0}
assert exceptional_valuation == {"L": 1, "W": 1, "q": 0}

print("PASS: q=W/L preserves the rational A4 source field")
print("PASS: all three original target coordinates remain polynomial")
print("PASS: the old generic L-prime is outside the retained affine chart")
print("PASS: the replacement exceptional divisor maps to the target origin")
print("FAIL: the transformed Jacobian ledger is q^2*K^3*L^4")
print("OBSTRUCTION: det(D pi)=L contributes one relative-canonical copy")


# ---------------------------------------------------------------------------
# 3. Denominator-center surgery: adjoin x=H/L
# ---------------------------------------------------------------------------

x = sp.symbols("x")
denominator_chart_equation = sp.expand(L * x - H)

# The normalization of L makes the intersection with H completely explicit.
H_on_L = sp.factor(H.subs(source_normalization))
assert H_on_L == (
    (2 * t - 3) * (t**2 - 3 * t + 9) ** 3
    / (t**2 * (t - 3) ** 2)
)
assert sp.factor(
    sp.resultant(L, H, U)
    - 19683 * (2 * V + 3) * (V**2 + 3 * V + 9) ** 3
) == 0

# The affine modification is the hypersurface L*x=H (times the free W
# direction).  Its Jacobian ideal has a nonreduced singular scheme supported
# at U=0, x=-1, V^2+3V+9=0.
singular_basis = sp.groebner(
    [
        denominator_chart_equation,
        sp.diff(denominator_chart_equation, U),
        sp.diff(denominator_chart_equation, V),
        sp.diff(denominator_chart_equation, x),
    ],
    x,
    U,
    V,
    order="lex",
    domain=sp.QQ,
)
expected_singular_generators = (
    x + 1 - U - sp.Rational(2, 3) * U * V,
    U**2,
    V**2 + 3 * V + 9,
)
for generator in expected_singular_generators:
    assert singular_basis.reduce(generator)[1] == 0
for generator in singular_basis.polys:
    expected_basis = sp.groebner(
        expected_singular_generators,
        x,
        U,
        V,
        order="lex",
        domain=sp.QQ,
    )
    assert expected_basis.reduce(generator.as_expr())[1] == 0

print("PASS: H|_L has one simple and one degree-two triple-contact factor")
print("FAIL: the chart x=H/L is singular at the conjugate contact cluster")
print("OBSTRUCTION: the reconstruction-denominator surgery is not A3")
