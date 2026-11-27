#!/usr/bin/env sage-python
"""Exact HG10 certificate for the repeated marked elliptic-factor HLP chart.

Run from the repository root with SageMath 10.9:
    sage -python research/arithmetic/verify_genus_three_hlp_repeated_factor.py

This classifies one parameter curve. It does not classify the full gluing
surface, prove geometric simplicity, or produce an infinite family.
"""

from sage.all import (QQ, Curve, EllipticCurve, EllipticCurve_from_cubic,
                      PolynomialRing, PowerSeriesRing, matrix, prod, vector)


P = PolynomialRing(QQ, names=("x", "y", "z"))
x, y, z = P.gens()
A = z**2 - 6*z*x - 3*x**2
F0 = (2*A*(z**2*y**4 - 6*z*y**2*x**3 - 3*x**6)
      - (A**2 + 64*z*x**3)*y*x**3)
assert F0.is_homogeneous() and F0.degree() == 8
factors = F0.factor()
assert len(factors) == 1 and factors[0][1] == 1

# Homogeneous coordinates are (x:y:z)=(1:s:b). Each substitution is a
# quadratic Cremona transformation conjugated by an invertible matrix.
M1 = matrix(QQ, [[0, 0, 1], [0, 1, 1], [1, 0, 1]])
M2 = matrix(QQ, [[-1, -QQ(1)/2, QQ(1)/8],
                 [0, 1, -QQ(3)/4], [1, 0, 1]])
assert M1.det() and M2.det()


def cremona(q):
    a, b, c = q
    return vector([b*c, a*c, a*b])


def strict_transform(f, mat, multiplicities):
    raw = f(*(mat*cremona(P.gens())))
    monomial = prod(v**m for v, m in zip(P.gens(), multiplicities))
    quotient, remainder = raw.quo_rem(monomial)
    assert not remainder
    return quotient / quotient.content()


F1 = strict_transform(F0, M1, [4, 4, 3])
F2 = strict_transform(F1, M2, [2, 3, 2])
assert F1.degree() == 5 and F2.degree() == 3
cubic = (32*x**3 + 144*x**2*y + 162*x*y**2 - 96*x**2*z
         - 216*x*y*z - 81*y**2*z + 72*y*z**2 - 16*z**3)
assert F2 == cubic or F2 == -cubic
assert not Curve(F2).is_singular()

# Every rational singular point on the original plane curve lies outside
# the open chart b not in {0,1,9}. Thus a rational point in that chart has
# a rational lift to the smooth normalization.
singular = Curve(F0).singular_points()
assert len(singular) == 4
for q in singular:
    assert not q[0] or q[2]/q[0] in (1, 9)

phi = EllipticCurve_from_cubic(F2, [0, 1, 0], morphism=True)
E = EllipticCurve(QQ, [0, 1, 0, -24, 36])
assert E.is_isomorphic(phi.codomain())
iso = E.isomorphism_to(phi.codomain())
assert E.rank(proof=True, use_database=False) == 0
assert list(E.torsion_subgroup().invariants()) == [2, 4]
points = E.torsion_points()
assert len(points) == 8

original_map = M1*cremona(M2*cremona(P.gens()))


def normalize(q):
    return q / next(a for a in q if a)


def resolve_original_point(cp):
    """Extend the rational map at its base points on the smooth cubic.

    The three sections have degree four. A nonzero section on a cubic
    has at most twelve zeros, counted with multiplicity, so precision 20
    suffices to recover its first nonzero coefficient. Hensel lifting and
    all coefficient operations take place over QQ, without approximation.
    """
    q = normalize(vector(QQ, list(cp)))
    direct = vector(QQ, [h(*q) for h in original_map])
    if direct:
        return normalize(direct)
    pivot = next(i for i, a in enumerate(q) if a)
    k = next(i for i in range(3)
             if i != pivot and F2.derivative(P.gen(i))(*q))
    j = next(i for i in range(3) if i not in (pivot, k))
    ps = PowerSeriesRing(QQ, "T", default_prec=20)
    qq = [ps(a) for a in q]
    qq[j] += ps.gen()
    for _ in range(5):
        qq[k] -= F2(*qq) / F2.derivative(P.gen(k))(*qq)
    assert not F2(*qq) or F2(*qq).valuation() >= 20
    sections = [h(*qq) for h in original_map]
    m = min(h.valuation() for h in sections if h)
    assert m <= 12
    return normalize(vector(QQ, [h[m] for h in sections]))


images = {tuple(resolve_original_point(phi.inverse()(iso(ep))))
          for ep in points}
expected = {(QQ(1), -3, -QQ(3)/5), (0, 1, 0), (1, -3, 1),
            (1, 1, 9), (1, -QQ(3)/2, 6), (1, -QQ(1)/3, 9),
            (0, 0, 1), (1, 2, 0)}
assert images == expected
assert all(F0(*q) == 0 for q in images)
admissible = {(q[2], q[1]) for q in images
              if q[0] and q[2] not in (0, 1, 9) and q[1]}
assert admissible == {(QQ(6), -QQ(3)/2), (-QQ(3)/5, QQ(-3))}

for b, s in sorted(admissible):
    a1 = (b**2 - 6*b - 3)/4
    a3 = (b**2*s**4 - 6*b*s**2 - 3)/4
    delta1 = a1**2 - 4*b
    delta3 = a3**2 - 4*b*s**2
    d3 = s*delta1/(2*a1)
    product_root = delta1*d3
    assert d3**2 == delta3 and delta1*delta3 != 0
    assert (product_root*(2*a1**2/delta1 + a3**2/delta3 - 1)
            - 2*a1**2*a3) == 0
    # HLP Proposition 14's quartic coefficients and conic signs.
    qa = product_root*b/2*(b*s**2/delta3)
    qb = qa
    qc = product_root*b*s**2/2*(2*b/delta1 - b*s**2/delta3)
    d = e = 1/(b*s)
    f = -1/b
    assert (b*b*s**2*d**2, b*b*s**2*e**2, b*b*f**2) == (1, 1, 1)
    assert (-qa*e*f, -qb*d*f, -qc*d*e) == (a1, a1, a3)
    assert e/d == 1 and f/d == -s and s < 0
    h1, h3 = (b - 3)/2, (b*s**2 - 3)/2
    print(f"b={b}, s={s}: h1=h2={h1}, h3={h3}; "
          f"conic X^2+Y^2+{-s}*Z^2=0 has no real point")

print("HG10 PASS: irreducible parameter curve; exact Cremona normalization; "
      "rank zero; all eight rational points; two admissible definite conics.")
