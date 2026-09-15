# Determinant 852: the full marked curve is a rank-one elliptic curve

The subsequent [full rational NS theorem](DET852_RATIONAL_MARKING_SOURCE_2026-09-15.md)
now completes primitive lattice and actual divisor descent. Remaining
admission gaps stated below belong to this earlier gate alone.

The subsequent [intrinsic CM-orbit doubling theorem](DET852_INTRINSIC_NONCM_POINT_2026-09-15.md)
now specifies a rational non-CM period without numerical coordinates.
The point-identification gap below records this earlier gate alone;
primitive NS and rational divisor descent remain open.

For K3-4ff75fec54d01662 the literal transcendental lattice is

```text
-2 0   1
 0 4   0
 1 0 106
```

Its full stable marked curve over Q is X_0^6(71)/<w_426>. This curve has a
rational point and Jacobian 426b1, hence is Q-isomorphic to

    y² + xy = x³ + x² - 286x + 1780.

The point (7,-17) has infinite order. Thus the full marked curve has
infinitely many rational points. A particular point has not yet been
identified as non-CM on the moduli curve, and primitive NS/divisor descent
has not been completed. **No new MW17 surface or arithmetic admission is
claimed by this gate.** The displayed point is not being called non-CM.

## Literal order and full stable group

The even Clifford basis (1,e0e1,e0e2,e1e2) has reduced trace pairing

```text
2 0   1    0
0 4   0    2
1 0 107    0
0 2   0 -212
```

Its reduced discriminant is426. The rational algebra (2,213/4) has
ramified primes2 and3. Squarefree reduced discriminant therefore identifies
this order as maximal at2 and3 and Eichler of level71 at the split prime71.
There is no extra integral level. This uses the local classification of
quaternion orders of squarefree reduced discriminant.

The discriminant group of T is cyclic of order852. Its orthogonal units
are1,143,283,425,427,569,709,851. Exact positive-norm normalizers are:

| Label | Even-basis coordinates | Norm | Discriminant multiplier |
|---|---|---:|---:|
| w2 | (28,-11,-10,-7) | 2 | 427 |
| w6 | (20,-11,-10,7) | 6 | 143 |
| w426 | (0,1,0,-2) | 426 | 851 |

They preserve both the order and T integrally. Their generated actions
exhaust all eight orthogonal units. The checker also proves a polynomial
identity: every coefficient of the numerator of (Ad(q)-I)G^-1 is integral.
Consequently all local order units act trivially on the local discriminant
group. Their reduced norms exhaust the local units, including at the
Eichler prime71. This checks the arithmetic level, not just a few global
norm-one elements.

The normalizer theorem leaves exactly the Atkin–Lehner cosets. Only1 and
w426 act by a simultaneous global sign on the discriminant group. Moreover
-Ad(w426) is the stable reflection in e0 of square-2. Thus projectivizing
the entire stable orthogonal group adds precisely w426 to the norm-one
curve. As in the [determinant388 arithmetic argument](DET388_RATIONAL_MARKING_SOURCE_2026-09-14.md#1-full-stable-group-with-the-literal-arithmetic-level),
the full local unit groups and global orientation identify the canonical
Q-curve, not just its complex uniformization. The arithmetic period-map
and stable-group theorem inputs are the same ones cited there.

Ogg's formula gives genus13 upstairs and24 fixed points of w426,
from the order of discriminant-1704 and class number24. Riemann–Hurwitz
gives quotient genus1. The checker reuses only the retained exact
fixed-point helper from the [determinant1236 calculation](DET1236_MARKED_SHIMURA_CURVE_2026-09-04.md).

## Published rational-point and elliptic-model inputs

[Padurariu–Saia](https://arxiv.org/abs/2509.25368) provide the complete
low-genus quotient tables. Their
[rational positive-rank table](https://github.com/fsaia/GenusAtMost2/blob/main/genus_1_AL_quotients_rat_pts_pos_rank.m)
contains (6,71,{426}); their
[Jacobian isomorphism table](https://github.com/fsaia/GenusAtMost2/blob/main/genus_1_AL_quotient_jacobian_isomorphism_classes.m)
assigns426b1. Both exact entries are checked against retained downloaded
bytes with source hashes. Rational-point existence is a separate input
from Jacobian identification: a genus-one torsor is not identified with its
Jacobian merely because the latter has positive rank.

The checker verifies the displayed equation, point, conductor and a
non-torsion proof using good-reduction group orders at5,7,11. The resulting
torsion bound does not kill(7,-17). No analytic rank assumption is needed
for the infinite rational locus. The Q-isomorphism to the moduli curve is
not normalized at a specified CM origin; therefore these coordinates alone
do not identify a particular non-CM period.

## Next admission gate and retained discovery

Identify a rational non-CM period with a certified moduli interpretation,
then verify the primitive NS complement and actual rational divisor descent.
Only after that admission should the NS enter a rootless-frame gate.
The601 historical unadmitted rows are not reduced by this Phase-1 result.
Determinants388 and622 remain globally excluded from MW17.

Discovery matched only literal quaternion/order data to the published
single full-Fricke quotient tables; it did not use frame rootlessness.
Among the retained surviving rows the positive-rank genus-one full-Fricke
match is this852 row. Quotient matches alone were not promoted to full
marking statements. The normalizer discovery bounds were b,c,d in[-12,12],
with the remaining coordinate solved from the quadratic norm equation.
The final checker pins its witnesses and does not replay that search.

```sh
sage -python research/elkies-k3/scripts/certify_det852_marked_curve.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det852-marking-gate-v1/certificate.json)
contains the literal order, normalizers, all local-unit polynomial
coefficients, full discriminant action, genus calculation and infinite-order
point check. The [source receipts](../artifacts/generated-results/elkies-k3-shimura-positive-source-preflight-v1/sources.json)
pin the imported tables. Published moduli/classification results remain
named theorem inputs; independent implementation, formal verification,
external review and novelty are unclaimed.
