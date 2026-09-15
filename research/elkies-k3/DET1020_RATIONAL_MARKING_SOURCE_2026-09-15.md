# Determinant 1020: a new full rational NS marking source

There exists a projective K3 surface over Q with geometric Picard rank19
and full saturated rational Neron–Severi marking

    S = (-T) + E8(-1) + E8(-1),
    T = [[-2,1,0],[1,2,0],[0,0,204]].

Its determinant is1020. This source lies outside the retained827-row
catalogue. Its [rootless frame and arithmetic MW17 existence](DET1020_ARITHMETIC_MW17_EXISTENCE_2026-09-15.md)
are now proved. Neither a K3 equation nor numerical period coordinates
have been computed.

## Full canonical marked curve

The even Clifford order in basis(1,e0e1,e0e2,e1e2) has reduced trace Gram

```text
2 1    0    0
1 3    0    0
0 0  204 -102
0 0 -102 -204
```

Its reduced discriminant is510. Its rational quaternion algebra(5,102)
is ramified precisely at2,3,5,17, hence this literal order is maximal.
The discriminant group of T is cyclic of order1020 and has16 orthogonal
units. The following exact integral normalizers generate their action:

| Label | Even-basis coordinates | Norm | Discriminant multiplier |
|---|---|---:|---:|
| w2 | (46,-16,-12,-7) | 2 | 511 |
| w3 | (174,-15,-15,4) | 3 | 341 |
| w5 | (54,-13,-14,-8) | 5 | 409 |
| w510 | (0,0,1,2) | 510 | 1019 |

All order and T actions are integral. A polynomial numerator identity
proves that every local order unit is stable on the discriminant module.
Maximal-order local reduced norms exhaust local units. Only the identity
and w510 have simultaneous global signs on the discriminant group;
-Ad(w510) is the stable reflection in the first vector of square-2.
The full canonical marking curve is therefore

    C = X^510 / <w510> over Q.

The maximal-order normalizer and canonical arithmetic period inputs are
the same ones used in the [388 full-group proof](DET388_RATIONAL_MARKING_SOURCE_2026-09-14.md#1-full-stable-group-with-the-literal-arithmetic-level).
This checks the entire marked group, not merely a norm-one or coarse
quotient. Ogg's formula gives genus9 upstairs and16 Fricke fixed points,
hence quotient genus1.

[Padurariu–Saia](https://arxiv.org/abs/2509.25368) separately prove rational
point existence for(510,1,{510}) and identify its Jacobian as510d2 in their
[published tables](https://github.com/fsaia/GenusAtMost2/blob/main/genus_1_AL_quotient_jacobian_isomorphism_classes.m).
Thus C is Q-isomorphic to

    y² + xy + y = x³ + x² - 421x - 3157.

Retained source bytes are checked for both entries. Exact PARI descent
gives rank bounds[1,1]. Rational torsion is(Z/2)^2, checked using its four
points and good-reduction group orders12,16,16 at7,11,13.

## Complete CM locus and an intrinsic non-CM point

A rational point of C lifts upstairs to degree at most two.
[Gonzalez–Rotger, Theorem5.8 and Corollary5.14](https://arxiv.org/pdf/math/0612732)
reduce rational CM orders to the complete class-number-one and two lists.
The checker tests all42 orders with the local optimal-embedding factors
at2,3,5,17 and the exact residue-field alternatives. No nonempty
class-number-two row reaches the possible rational two-involution case.
The only surviving orders are:

| Discriminant | Points upstairs | Rational points on C |
|---|---:|---:|
| -3 | 8 | 4 |
| -163 | 16 | 8 |

The Fricke fixed order has discriminant-2040 and class number16, so it
fixes neither of these loci. Dividing by two gives the complete counts.
The w3 fixed-point calculation has precisely the eight CM(-3) points
upstairs. Thus residual w3 fixes all four rational CM(-3) points on C.

Choose any rational CM(-3) point O as elliptic origin, and any rational
CM(-163) point R. Then

    O + 2(R-O)

is a particular rational non-CM point. Every such choice works, without
first computing numerical CM coordinates.

To prove this, the residual Atkin–Lehner group is faithfully(Z/2)^3:
it is the quotient of the full sixteen-element group by<w510>. Its
involution w3 fixes O and hence is inversion. Its four-element subgroup
acting trivially on differentials consists of translations by the entire
rational2-torsion group. Consequently the orbit of O is exactly the four
rational2-torsion points, and exactly the CM(-3) locus.

Since all rational torsion is2-torsion, R is non-torsion. Its orbit under
these translations and inversion is the eight distinct points
+/-R+T, T in C(Q)[2], hence the entire CM(-163) locus. The double2R is not
torsion. If it equalled R+T or-R+T, then R or3R would be torsion,
a contradiction. Thus2R belongs to neither CM locus and is non-CM.
This specifies an actual period intrinsically, rather than merely proving
that an unidentified member of a finite sample is non-CM.

## Primitive K3 embedding and rational divisor descent

There is a direct integral embedding certificate. On U^3 use the Gram
H=[[0,I],[I,0]], and let B be upper triangular with B+B^t=T. The rows

    [I | B]       and       [I | -B^t]

have Grams T and-T, are orthogonal, and are both primitive because each
has an identity leading minor. Adding two negative E8 summands to H gives
the K3 lattice. Thus S is the actual primitive orthogonal complement of T,
not a finite-index approximation. The checker stores these bases and the
full rank19 NS Gram.

At the intrinsic non-CM period, geometric NS is exactly S. Choose its ample
chamber. As in the [852 descent proof](DET852_RATIONAL_MARKING_SOURCE_2026-09-15.md),
a finite marked automorphism on the rank-three non-CM Hodge structure can
only act as+I or-I on T; -I is incompatible with the identity on S across
the exponent1020 discriminant glue. Torelli gives trivial marked inertia,
so the rational marked period yields effective descent of the projective
surface with all NS classes invariant.

Two adjacent negative E8 roots r,s have squares-2 and intersection1.
Their Euler characteristics are1, killing their line-bundle Brauer
obstructions. Their descended intersection gives a zero-cycle of degree
one, making Br(Q)->Br(X) injective. The Picard-Brauer exact sequence then
descends every invariant geometric line bundle. This proves full actual
saturated rational divisor marking.

S also contains an explicit U: if e0,e1 are its first two coordinates and
r is the first root of the first negative E8, put v=e0+r and w=v-e1.
Their squares are zero and v.w=1. This gives a primitive abstract elliptic
frame to inspect after arithmetic admission; no rootlessness is asserted.

## Replay and scope

```sh
sage -python research/elkies-k3/scripts/certify_det1020_rational_marking.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det1020-marking-gate-v1/certificate.json)
contains the full Clifford action, local-unit polynomial check, all CM
orders, involution fixed-point data, primitive embedding and descent pair.
Discovery normalizers used b,c,d in[-16,16] and solved the remaining
coordinate from the norm equation; the final checker pins the witnesses.
Moduli descent, the group-orbit non-CM argument and the Brauer argument are
written proof inputs. Independent implementation, formal verification,
external review and literature novelty are unclaimed.

The rootless-U gate is now closed positively by the linked existence proof.
The next task is an explicit equation and17 saturated rational sections.
The600 unadmitted historical rows remain a separate queue. The closed
388,622,852 branches stay excluded; this source does not reopen them.
