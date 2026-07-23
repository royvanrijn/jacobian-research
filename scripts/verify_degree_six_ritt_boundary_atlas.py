#!/usr/bin/env python3
"""Exact witnesses and finite cuts for the degree-six Ritt boundary atlas."""

from __future__ import annotations

import sympy as sp


w = sp.symbols("w")
h2, h3, h4, h5, h6 = sp.symbols("h2 h3 h4 h5 h6")
x, y, z = sp.symbols("x y z")

# Normalized coefficient chart h4=x,h5=y,h6=z.
h3_normal = -1 - 2 * x - 3 * y - 4 * z
h2_normal = 1 + x + 2 * y + 3 * z
H_generic = sp.expand(
    h2_normal * w**2
    + h3_normal * w**3
    + x * w**4
    + y * w**5
    + z * w**6
)
phi23 = sp.expand(
    32 * h3 * h5 * h6**2
    + 64 * h3 * h6**3
    + 16 * h4**2 * h6**2
    - 24 * h4 * h5**2 * h6
    + 64 * h4 * h6**3
    + 5 * h5**4
    + 64 * h5 * h6**3
    + 64 * h6**4
).subs({h3: h3_normal, h4: x, h5: y, h6: z})
phi32 = (
    27 * h3 * h6**2 - 18 * h4 * h5 * h6 + 5 * h5**3
).subs({h3: h3_normal, h4: x, h5: y, h6: z})
psi33 = 9 * h2_normal * z**2 - 3 * x**2 * z + x * y**2

P_generic = sp.cancel(H_generic / w**2)
affine_boundary = sp.factor(
    sp.resultant(P_generic, 2 * P_generic + w * sp.diff(P_generic, w), w)
)
extra_root = sp.factor(affine_boundary / (z * h2_normal))
assert sp.Poly(extra_root, x, y, z).total_degree() == 4


def coefficients(H: sp.Expr) -> dict[sp.Symbol, sp.Expr]:
    """Return normalized chart coefficients for a sextic seed."""
    polynomial = sp.Poly(sp.expand(H), w)
    return {
        x: polynomial.coeff_monomial(w**4),
        y: polynomial.coeff_monomial(w**5),
        z: polynomial.coeff_monomial(w**6),
    }


def audit_seed(H: sp.Expr, *, zero_order: int, gcd_degree: int) -> None:
    """Check the common weighted, Hessian, and primitive-root conditions."""
    H = sp.factor(H)
    assert sp.Poly(H, w).degree() == 6
    assert H.subs(w, 0) == sp.diff(H, w).subs(w, 0) == 0
    assert H.subs(w, 1) == 0
    assert sp.diff(H, w).subs(w, 1) == -1
    assert sp.diff(H, w, 2).subs(w, 1) != -2
    assert sp.discriminant(sp.diff(H, w, 2), w) != 0
    polynomial = sp.Poly(H, w)
    actual_zero_order = min(
        exponent[0] for exponent, coefficient in polynomial.terms() if coefficient != 0
    )
    assert actual_zero_order == zero_order
    primitive_gcd = sp.gcd(
        sp.Poly(H / w**2, w), sp.Poly(sp.diff(H, w) / w, w)
    )
    assert primitive_gcd.degree() == gcd_degree


# -------------------------------------------------------------------------
# Dense decomposition charts and their factored affine boundaries.
a, b, c, p, q = sp.symbols("a b c p q")
D23 = 2 * p**2 + 5 * p - q + 3
a23 = sp.cancel(-(p + 1) / ((p + q + 1) * D23))
b23 = sp.cancel((p + q + 1) / D23)
B23 = w**3 + p * w**2 + q * w
H23 = sp.factor(a23 * B23**2 + b23 * B23 - b23 * q * w)
U23 = p**3 + 2 * p**2 * q + 2 * p**2 + 2 * p * q + p - q**2
V23 = 4 * p**2 * q + p**2 - 10 * p * q - 2 * p - 27 * q**2 - 14 * q - 3

b32 = sp.cancel(
    -(a * p**3 + 6 * a * p**2 + 9 * a * p + 4 * a + 1)
    / (2 * (p + 1))
)
c32 = sp.cancel(
    (
        a * p**4
        + 5 * a * p**3
        + 9 * a * p**2
        + 7 * a * p
        + 2 * a
        + p
        + 1
    )
    / 2
)
B32 = w**2 + p * w
H32 = sp.factor(a * B32**3 + b32 * B32**2 + c32 * B32 - c32 * p * w)
X32 = a * p**3 + 4 * a * p**2 + 5 * a * p + 2 * a + 1
Y32 = 5 * a * p**3 + 12 * a * p**2 + 9 * a * p + 2 * a + 2 * p + 1
Z32 = sp.expand(
    25 * a**2 * p**6
    + 155 * a**2 * p**5
    + 379 * a**2 * p**4
    + 457 * a**2 * p**3
    + 272 * a**2 * p**2
    + 64 * a**2 * p
    + 23 * a * p**3
    + 47 * a * p**2
    + 24 * a * p
    - 2
)

# Five rational witnesses, one on every irreducible affine-boundary curve.
boundary_witnesses = {
    "D23-zero": sp.factor(H23.subs({p: 0, q: 0})),
    "D23-extra": sp.factor(H23.subs({p: -2, q: -sp.Rational(5, 27)})),
    "D32-extra-X": sp.factor(H32.subs({p: -6, a: sp.Rational(1, 100)})),
    "D32-zero": sp.factor(H32.subs({p: -6, a: -sp.Rational(11, 700)})),
    "D32-extra-Z": sp.factor(H32.subs({p: sp.Rational(1, 2), a: sp.Rational(8, 147)})),
}
expected_witnesses = {
    "D23-zero": -w**3 * (w - 1) * (w**2 + w + 1) / 3,
    "D23-extra": -w**2 * (w - 1) * (9 * w - 17) ** 2 * (9 * w + 7) / 1024,
    "D32-extra-X": w**2 * (w - 6) ** 2 * (w - 5) * (w - 1) / 100,
    "D32-zero": -w**3 * (w - 1) * (11 * w**2 - 187 * w + 876) / 700,
    "D32-extra-Z": w**2 * (w - 1) * (2 * w - 5) * (2 * w + 5) ** 2 / 147,
}
for name, witness in boundary_witnesses.items():
    assert sp.expand(witness - expected_witnesses[name]) == 0
    audit_seed(
        witness,
        zero_order=3 if name.endswith("zero") else 2,
        gcd_degree=1,
    )

substitutions = {name: coefficients(H) for name, H in boundary_witnesses.items()}

# Membership and exclusion assertions separate all five components.
assert phi23.subs(substitutions["D23-zero"]) == 0
assert h2_normal.subs(substitutions["D23-zero"]) == 0
assert extra_root.subs(substitutions["D23-zero"]) != 0
assert phi32.subs(substitutions["D23-zero"]) != 0

assert phi23.subs(substitutions["D23-extra"]) == 0
assert h2_normal.subs(substitutions["D23-extra"]) != 0
assert extra_root.subs(substitutions["D23-extra"]) == 0
assert phi32.subs(substitutions["D23-extra"]) != 0

for name in ("D32-extra-X", "D32-extra-Z"):
    assert phi32.subs(substitutions[name]) == 0
    assert h2_normal.subs(substitutions[name]) != 0
    assert extra_root.subs(substitutions[name]) == 0
    assert phi23.subs(substitutions[name]) != 0
    assert psi33.subs(substitutions[name]) != 0

assert phi32.subs(substitutions["D32-zero"]) == 0
assert h2_normal.subs(substitutions["D32-zero"]) == 0
assert extra_root.subs(substitutions["D32-zero"]) != 0
assert phi23.subs(substitutions["D32-zero"]) != 0
assert psi33.subs(substitutions["D32-zero"]) != 0

# Verify the chart factors at the same points.
assert U23.subs({p: 0, q: 0}) == 0
assert V23.subs({p: 0, q: 0}) != 0
assert U23.subs({p: -2, q: -sp.Rational(5, 27)}) != 0
assert V23.subs({p: -2, q: -sp.Rational(5, 27)}) == 0
assert X32.subs({p: -6, a: sp.Rational(1, 100)}) == 0
assert Y32.subs({p: -6, a: -sp.Rational(11, 700)}) == 0
assert Z32.subs({p: sp.Rational(1, 2), a: sp.Rational(8, 147)}) == 0

# The D23 boundary witnesses have explicit all-double omitted values.
assert sp.factor(
    boundary_witnesses["D23-zero"]
    - sp.Rational(1, 3) / 4
    + sp.Rational(1, 3) * (w**3 - sp.Rational(1, 2)) ** 2
) == 0
assert sp.factor(
    boundary_witnesses["D23-extra"]
    + sp.Rational(5, 27) * w
    - sp.Rational(256, 729)
    + sp.Rational(729, 1024) * (B23.subs({p: -2, q: -sp.Rational(5, 27)}) + sp.Rational(512, 729)) ** 2
) == 0

# -------------------------------------------------------------------------
# A centered global chart for the common Ritt curve.
center = sp.symbols("center")
g = ((w - center) ** 3 + q * (w - center)) ** 2
g0 = g.subs(w, 0)
gp0 = sp.diff(g, w).subs(w, 0)
endpoint = sp.factor(g.subs(w, 1) - g0 - gp0)
D = sp.factor(sp.diff(g, w).subs(w, 1) - gp0)
Qcommon = sp.cancel((g - g0 - gp0 * w) / w**2)

expected_endpoint = (
    15 * center**4
    - 20 * center**3
    + 12 * center**2 * q
    + 15 * center**2
    - 8 * center * q
    - 6 * center
    + q**2
    + 2 * q
    + 1
)
assert sp.expand(endpoint - expected_endpoint) == 0
zero_common = sp.factor(Qcommon.subs(w, 0))
assert zero_common == 15 * center**4 + 12 * center**2 * q + q**2
residual_discriminant = sp.factor(sp.discriminant(Qcommon, w))
assert residual_discriminant == (
    16
    * center**2
    * (center**2 + q) ** 2
    * (3 * center**2 + q)
    * (225 * center**4 + 285 * center**2 * q + 64 * q**2)
)
assert sp.factor(sp.discriminant(sp.diff(g, w, 2), w)) == 108380160 * q**6

# Six reduced zero-cluster points, forming one irreducible rational orbit.
zero_cut = sp.factor(sp.resultant(endpoint, zero_common, q))
expected_zero_cut = -(
    560 * center**6
    - 840 * center**5
    + 411 * center**4
    - 20 * center**3
    - 42 * center**2
    + 12 * center
    - 1
)
assert zero_cut == expected_zero_cut
assert sp.factor(zero_cut) == zero_cut
assert sp.gcd(zero_cut, sp.diff(zero_cut, center)) == 1
assert sp.discriminant(zero_cut, center) == 83901850583040
zero_q = (
    1120 * center**5
    - 1400 * center**4
    + 472 * center**3
    + 73 * center**2
    - 62 * center
    + 7
) / 2

# The six extra-root points split as two rational points and a quartic orbit.
assert sp.expand(
    endpoint.subs(q, -3 * center**2)
    + (center + 1) * (2 * center - 1) ** 2 * (3 * center - 1)
) == 0
for rational_point in (
    {center: -1, q: -3},
    {center: sp.Rational(1, 3), q: -sp.Rational(1, 3)},
):
    assert endpoint.subs(rational_point) == 0
    assert D.subs(rational_point) != 0
    assert zero_common.subs(rational_point) != 0
    assert residual_discriminant.subs(rational_point) == 0

quartic_cut = 180 * center**4 - 300 * center**3 - 49 * center**2 + 208 * center - 64
quartic_q = (60 * center**3 - 80 * center**2 - 43 * center + 28) / 12
assert sp.factor(quartic_cut) == quartic_cut
assert sp.gcd(quartic_cut, sp.diff(quartic_cut, center)) == 1
assert sp.discriminant(quartic_cut, center) == 13168189440000
for equation in (endpoint, 225 * center**4 + 285 * center**2 * q + 64 * q**2):
    numerator = sp.together(equation.subs(q, quartic_q)).as_numer_denom()[0]
    assert sp.rem(sp.Poly(numerator, center), sp.Poly(quartic_cut, center)) == 0

# The apparent remaining factors in the residual discriminant lie on the
# normalization-pole divisor D=0 and hence do not define normalized seeds.
assert sp.factor(endpoint.subs(center, 0)) == (q + 1) ** 2
assert D.subs({center: 0, q: -1}) == 0
assert sp.factor(endpoint.subs(q, -center**2)) == (
    center - 1
) ** 2 * (2 * center - 1) ** 2
assert D.subs({center: 1, q: -1}) == 0
assert D.subs({center: sp.Rational(1, 2), q: -sp.Rational(1, 4)}) == 0

# Four reduced Hessian collisions on the common curve.  They are exactly the
# type-(6) omitted-component support, and they remain affine-boundary clean.
admissibility_numerator = sp.factor(sp.diff(g, w, 2).subs(w, 1) - 2 * D)
type6_cut = sp.factor(endpoint.subs(q, 0))
assert type6_cut == 15 * center**4 - 20 * center**3 + 15 * center**2 - 6 * center + 1
assert sp.gcd(type6_cut, sp.diff(type6_cut, center)) == 1
assert sp.discriminant(type6_cut, center) == 10800
Hcommon = sp.cancel(-(g - g0 - gp0 * w) / D)
Pcommon = sp.cancel(Hcommon / w**2)
common_affine_boundary = sp.factor(
    sp.resultant(Pcommon, 2 * Pcommon + w * sp.diff(Pcommon, w), w).subs(q, 0)
)
for nonvanishing in (D.subs(q, 0), admissibility_numerator.subs(q, 0), common_affine_boundary):
    numerator = sp.together(nonvanishing).as_numer_denom()[0]
    assert sp.gcd(sp.Poly(type6_cut, center), sp.Poly(numerator, center)).degree() == 0

# Every algebraic orbit avoids normalization poles and weighted
# inadmissibility.  The zero and extra cuts also avoid the Hessian divisor.
for orbit, q_expression in ((zero_cut, zero_q), (quartic_cut, quartic_q)):
    for forbidden in (q_expression, D.subs(q, q_expression), admissibility_numerator.subs(q, q_expression)):
        numerator = sp.together(forbidden).as_numer_denom()[0]
        assert sp.gcd(sp.Poly(orbit, center), sp.Poly(numerator, center)).degree() == 0

# The common clean witness lies outside all three finite deletion divisors.
clean_point = {center: sp.Rational(2, 5), q: -sp.Rational(1, 5)}
assert endpoint.subs(clean_point) == 0
assert D.subs(clean_point) != 0
assert zero_common.subs(clean_point) != 0
assert residual_discriminant.subs(clean_point) != 0
assert q.subs(clean_point) != 0

print("PASS degree-six Ritt boundary: five rational component witnesses")
print("PASS degree-six Ritt boundary: every witness is Hessian-clean and weighted-admissible")
print("PASS degree-six Ritt boundary: common zero-cluster cut is one reduced sextic orbit")
print("PASS degree-six Ritt boundary: common extra-root cut splits as two rational points plus a quartic orbit")
print("PASS degree-six Ritt boundary: four Hessian collisions are exactly the type-(6) support")
