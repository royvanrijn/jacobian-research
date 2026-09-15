# Determinant 852: exactly eight rational CM points

The subsequent [full rational NS theorem](DET852_RATIONAL_MARKING_SOURCE_2026-09-15.md)
now completes primitive lattice and actual divisor descent. Remaining
admission gaps stated below belong to this earlier gate alone.

The subsequent [intrinsic CM-orbit doubling theorem](DET852_INTRINSIC_NONCM_POINT_2026-09-15.md)
now specifies a rational non-CM period without numerical coordinates.
The point-identification gap below records this earlier gate alone;
primitive NS and rational divisor descent remain open.

The full marked curve C=X_0^6(71)/<w426> has exactly eight rational CM
points: four of discriminant-67 and four of discriminant-163. Their
coordinates on the [426b1 elliptic model](DET852_FULL_MARKED_CURVE_2026-09-15.md)
remain unidentified. This result does not label(7,-17), or any other
particular elliptic-model point, as non-CM.

## Complete finite reduction

A rational point of C lifts to a point P of X_0^6(71) of degree at most two.
If P is CM by an imaginary quadratic order R in K, the ring class field
satisfies H_R=K Q(P), by
[Gonzalez–Rotger, Theorem5.8](https://arxiv.org/pdf/math/0612732).
Thus h(R)<=2. The complete class-number-one and class-number-two lists
contain13 and29 orders, respectively. The checker retains these lists
from the existing determinant1236 certificate and verifies every order's
class number using reduced primitive binary quadratic forms. Completeness
of the lists is an imported theorem, not inferred from a bounded search.

For each order, the optimal-embedding count upstairs is

    h(R) (1-(R/2)) (1-(R/3)) (1+(R/71)),

where the Eichler symbol is1 at primes dividing the conductor and otherwise
the quadratic symbol of K. The checker then applies the residue-field
alternatives of Corollary5.14 of the same paper, with m=426. Write
A=D(R)N*(R) and q=m/m_r, using that paper's definitions.

For h=1, a nonempty row has rational quotient image exactly when A=1 or
q=A. For h=2, rationality can occur only in the two-involution fixed-field
case A=q=1. No nonempty h=2 row reaches this case here, so no unresolved
Artin ideal-class choice remains. The exact local factors and residue-field
parameters for all42 orders are in the certificate.

Only two rows survive:

| Order discriminant | h | Local factors at2,3,71 | Points upstairs | A=q | Rational images |
|---|---:|---|---:|---:|---:|
| -67 | 1 | 2,2,2 | 8 | 426 | 4 |
| -163 | 1 | 2,2,2 | 8 | 426 | 4 |

The w426 fixed locus has discriminant-1704, with class number24, as checked
in the full-curve certificate. Neither surviving CM locus is fixed, so
each eight-point locus descends in pairs to four distinct rational points.
Their orders differ, so the two image loci are disjoint. The count eight
is complete.

## Finite non-CM test set and exact remaining gap

Let E be426b1 and P=(7,-17), already proved non-torsion. The nine points

    O, P, 2P, ..., 8P

are distinct. The checker writes all nine exact projective coordinates.
Under any Q-isomorphism C->E, at least one of these points is non-CM,
because C has only eight rational CM points. This gives a finite certified
set containing a non-CM point, not an identified member or a normalized
moduli isomorphism. In particular, infinite order on E does not itself
exclude CM on C.

The next task is to normalize the elliptic model at a specified CM point
and locate the two four-point CM orbits, or otherwise certify a particular
rational period as non-CM. Primitive NS and actual rational divisor descent
then complete arithmetic admission. No852 rootless-frame calculation has
been authorized by this partial gate, and the601 unadmitted count is unchanged.

```sh
sage -python research/elkies-k3/scripts/certify_det852_rational_cm_locus.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det852-cm-gate-v1/certificate.json)
pins all42 orders, local counts, residue-field tests, eight-point total and
the nine-point test set. The preflight is retained in the same packet.
The class-number classification, CM field formula and optimal-embedding
formula are named proof inputs. Independent implementation, formal
verification, external review and novelty are unclaimed.
