# Cubic contact interpolation and the historical trielliptic gate

**Subsequent geometric closure:** the [F1 adjoint obstruction](F1_TRIELLIPTIC_ADJOINT_OBSTRUCTION_2026-09-15.md)
now excludes the degree-three elliptic map discussed below and proves
finiteness of this image-genus layer modulo inherited translations. The
explicit interpolation criterion remains applicable to individual exceptional
cubic contacts. The necessary-condition discussion below is retained as
the derivation of the now-closed moving-supply route.

**Written structural criterion with symbolic checks; no new MW17 pair.**
For bisections of arithmetic genus four and normalization genus one on a
24-I1 elliptic K3, the relevant contact divisor has degree three on a
genus-nine halving curve. In a finite, unramified quotient chart its
interpolation reduces to the rank of an explicit three-by-two matrix.
A successful member has anti-trace height 28.

Unlike the preceding two image-genus layers, finiteness of the arithmetic
contact set does not follow from the genus and the degree-four ruling alone.
An infinite collection in this layer would force a degree-three map from
one of the halving curves to an elliptic curve of positive rank over the
ground number field. Such a map is necessary, not sufficient: its fibres
must still pass the contact interpolation and residual-cover checks.
No such map on an actual MW17 parent's halving curve is supplied here.

## Geometry of the layer

Use the [bisection height identity](SINGULAR_BISECTION_TRACE_GEOMETRY_2026-09-15.md)
and [contact rigidity](QUADRATIC_CONTACT_INTERPOLATION_2026-09-15.md).
For arithmetic genus a=4 and z=B.O, the trace satisfies

```
h(T)=4z+2,       h(Q-sigma(Q))=28.
```

The first identity makes T geometrically primitive modulo twice MW: the
height of twice a section in this even lattice is divisible by eight.
The quotient by P -> T-P is a split Hirzebruch surface F_n, n<=2, with
smooth branch curve R=D_T of genus nine. The image C of B satisfies

```
C²=3,        R.C=10,        deg Z=3,
Z=sum floor(local contact multiplicity/2)*P.
```

Writing C=C0+kF gives -n+2k=3, hence n=1 and k=2. Thus the candidate
linear system is |C0+2F| on F1, of projective dimension four. In the plane
obtained by contracting C0 it consists of conics through the contracted
point. Removing twice the degree-three contact divisor leaves branch
degree four. For smooth parent branch fibres the anti-trace height is 28.

If its literal quadratic extension matches a verified image of arithmetic
genus 0, 1, 2 or 3, the height ratios are respectively 7/3, 7/4, 7/5 and 7/6.
All are nonsquares, so the new direction is independent of that earlier
direction. Matching the actual extension, rather than only its geometric
branch support, is essential. Infinite rational points on the covering
base remain a separate requirement.

## Exact cubic-contact constructor

Work in an F1 chart in which the candidate section is

```
A+B*t+C*t²+(D+E*t)*u=0.                         (1)
```

Let f(t) be monic separable cubic, specifying three distinct original base
values. At the actual contact divisor on R, assume the original base
coordinate generates the residue algebra and dt is nonzero. Write u and
v=du/dt as polynomials of degree less than three modulo f. This covers an
irreducible cubic contact point and split separable cubic contact divisors.
It does not create a point on R from arbitrary field elements u,v.

There is a unique polynomial U of degree less than six with

```
U mod f=u,        U' mod f=v.
U=u+f*((v-u')/f' mod f).                        (2)
```

Write f²=t⁶+c5*t⁵+...+c0 and U=z5*t⁵+...+z0. Equation (1) has the
required value and tangent at every contact precisely when

```
A+B*t+C*t²+(D+E*t)*U = 0 mod f².
```

Eliminating A,B,C leaves

```
       [ z3   z2-z5*c3 ]
M  =   [ z4   z3-z5*c4 ],      M * [D,E]^T = 0.  (3)
       [ z5   z4-z5*c5 ]
```

Thus all three two-by-two minors must vanish. Rank one gives the unique
projective denominator (D,E); take the negative of the remainder of
(D+E*t)*U modulo f² for the numerator A+B*t+C*t². Require D+E*t to be a
unit modulo f and require the homogeneous numerator of degree two and
denominator of degree one to have no common zero. This excludes vertical
components, including a common factor at infinity. With these conditions,
the reconstructed divisor is an integral section in the specified class.

Rank two gives no candidate. Rank zero cannot be accepted as a family of
integral candidates: contact rigidity would force any two integral members
to be identical, while the kernel supplies two independent equations.
In this chart it occurs when U has degree at most one; every equation then
factors and fails homogeneous coprimality.

Finally restrict the actual double-cover branch equation to C, divide by
the exact square contact factor, and check that the residual binary quartic
is squarefree with the correct literal scalar. Squarefreeness permits
contacts of order three as well as two, but rejects additional delta. It
also proves that the inverse image is geometrically integral. Transport
omitted points and poles with the ruled surface. Repeated contact divisors
require higher jets and are covered by rigidity, but not by the separable
cubic formula (2).

## Why a moving arithmetic supply would have to be trielliptic

Fix a trace representative modulo twice E(k(t)). Contact rigidity gives at
most one candidate per effective degree-three divisor on R. The preceding
[quadratic-contact theorem](QUADRATIC_CONTACT_INTERPOLATION_2026-09-15.md)
proves finiteness of degree-at-most-two points on R. Therefore infinitely
many candidates for this trace would require infinitely many closed cubic
points, not merely repetitions or combinations of lower-degree points.

The [Kadets–Vogt genus bound, Theorem 1.3](https://arxiv.org/pdf/2208.01067)
then forces a degree-three map R -> Y with infinitely many k-rational
points: its exceptional genus bound for minimal density degree three is
four, below genus nine. Thus Y has genus zero or one. A map of degree
three to genus zero is impossible by Castelnuovo–Severi against the ruling
of degree four: the two maps generate the function field since 3 and 4
are coprime, and their genus bound is six.

For genus one the same bound is exactly nine, so it gives no contradiction.
The only remaining possibility is a degree-three map to an elliptic curve
of positive rank over k. There are finitely many trace representatives;
hence the implication also holds for infinitely many such bisections
modulo inherited translations. Conversely, absence of these elliptic maps
for all relevant traces proves finiteness of this entire image-genus layer.
This is not an effective list or an emptiness result for individual covers.

A degree-three elliptic map by itself is insufficient. Its fibres are only
potential contact divisors. One must check (3) on actual fibres, integral
reconstruction, the quartic, and any match to a second independent direction.
The elliptic curve parametrizing contacts and the elliptic covering base
of the original fibration are different curves; positive rank of the former
does not prove positive rank of the latter.

## Evidence and boundary

The [checker](scripts/verify_cubic_contact_interpolation.py) verifies (2),
(3), the three symbolic minors, equality with the six original jet
conditions, and the integer intersection/genus premises. Its controls use
Q[t]/(t³-t-1): one valid rational section, one incompatible derivative,
one reducible rank-zero case, and one valid polynomial section. These are
linear-algebra controls, not points on a K3 halving curve.

The [certificate](../artifacts/generated-results/elkies-k3-cubic-contact-interpolation-v1/result.json)
records these checks. The quotient geometry, contact rigidity and arithmetic
classification are written or inherited, without formal verification or
external review. No Jacobian computation, map search, or arithmetic point
campaign is run or scheduled by this note.

```
.venv/bin/python research/elkies-k3/scripts/verify_cubic_contact_interpolation.py
```

The [original MW17 two-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains open. This criterion identifies a possible source of a moving contact
supply and an exact way to test it; it does not assert that this source exists.
