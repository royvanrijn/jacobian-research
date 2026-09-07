# Six additional parents have geometric NS determinant 468

For each fixed Mestre outer parameter `u=11,13,17,19,23,29`, the constructed
K3 surface has geometric Picard rank **19**, rational Picard rank **18**,
geometric Néron–Severi determinant **+468**, and rational fixed-lattice
determinant **−468**. Its displayed elliptic fibration has Mordell–Weil rank
exactly **11 over Q(T)** and **12 over Qbar(T)**. The supplied eleven-point
basis is the full arithmetic Mordell–Weil group; torsion is zero.

These are additional parent surfaces, not additional names for the production
surface X948. Their different geometric NS determinant excludes geometric
isomorphism with X948. The [earlier parent audit](PARENT_PORTFOLIO_AND_SECTION_LABEL_AUDIT_2026-09-07.md)
already separated the six from each other over Q. We do **not** claim that
their geometric lattices differ from one another, or that these surfaces
are absent from the literature.

Every elliptic fibration with section over Q on one of these parents has
generic MW rank at most **16**. A useful rank-16 fibration remains
**unconstructed**. The subsequent [full NS matrix proof](MESTRE_FULL_NS_GRAMS_2026-09-07.md)
now gives identical rational and geometric Gram matrices on all six.
The marked transcendental lattice remains **UNKNOWN**. This is different lattice
territory, but not admission to the arithmetic-MW17 foundry, a new
near-record elliptic curve, or an upper bound on specialized elliptic ranks.

## Models and evidence

The six fixed models are the exact Jacobians produced by
[`mestre_parent_adapter.py`](../cas/mestre_parent_adapter.py), with
`v=(u^2+u+2)/u` and coherent root labels and ordinate branches. The
[self-contained replay input](../../artifacts/generated-results/elliptic-curves/mestre_468_replay_bundle_v1.json)
contains all six rational-function equations, section coordinates, integral
relation words, height matrices, finite reduction witnesses, and all twelve
finite-field count records. Thus the proof does not depend on an unavailable
coefficient table or on regenerating a point search.

The [twelve-fibre pilot](MESTRE_PARENT_CALIBRATION_2026-09-07.md) is unchanged:
all 588 boxes completed and only `(u,T)=(11,1),(13,1)` gained a direction.
With the present generic upper bound, their certified ranks at least 12
are now proved arithmetic specialization jumps from generic rank 11.
Neither is a near-record result. No new parameter sweep or point exposure
is part of this proof.

## Exact height lattice and base involution

The minimal geometric fibre configuration is `16 I1 + 2 I2 + I4`.
The two I2 base values are rational. Their nodes are nonsplit, but each
unique nonidentity component class is Galois fixed. The I4 fibre at infinity
is split. Hence the rational trivial lattice is
`U + A1(-1)^2 + A3(-1)`, of rank 7 and absolute determinant 16.

[`mestre_generic_height_audit.sage`](../cas/mestre_generic_height_audit.sage)
computes the exact heights of the fourteen coherently labelled covariant
images. On a K3 elliptic surface,

```
h(P) = 4 + 2(P.O) - sum_v contr_v(P),
contr_v(P) = i(n-i)/n  for component i of an In fibre.
```

Finite node tests and the first blowup at infinity determine the components.
The computation checks doubling, component addition, and all 105 pairwise
sum/difference parallelograms. The fourteen-image height matrix has rank 11.
Three exact rational-function group relations verify its kernel. After
passing to the divisor-section cloud `C0=P0`, `Cj=(Pj-P0)/2`, the fixed
indices `[0,1,2,3,4,5,6,7,8,9,12]` give the same positive definite
eleven-dimensional height matrix on all six parents, with determinant
`117/4`. These sections and the trivial lattice generate an integral
rank-18 NS sublattice of absolute determinant `16*(117/4)=468`.

The standalone V2 replay directly solves all thirteen halving problems
over Q(T) on each parent, verifies every doubling identity, and checks
that their T=1 images are precisely the specialized points used for
saturation. Thus these are rational generic sections, not merely formal
halves in a rational vector space or halves found on one specialization.

The [independent height replay](../cas/verify_mestre_generic_heights.sage)
avoids the producer's component selector. Multiplication by four puts
every section in the narrow subgroup, so it computes
`h(P)=(4+2(4P.O))/16` directly. All **714** checks pass, 119 per parent.

All six equations descend under `s=T^2` to a geometrically rational elliptic
surface with `8 I1 + 2 I2`. Its geometric MW rank is 6. The
[base-involution computation](../cas/mestre_base_involution.sage)
verifies every column of the `T -> -T` action by exact Q(T) group law,
as well as `M^2=I` and height preservation. The arithmetic invariant and
anti-invariant ranks are respectively **5 and 6**. Pulling back the full
geometric quotient MW group supplies six invariant directions, independent
of the six known anti-invariant directions. Thus geometric MW rank is at
least 12 and geometric Picard rank at least 19.

The height, base-change and Shioda–Tate formulas used here are standard;
see [Schütt–Shioda, Elliptic Surfaces, especially section 11](https://arxiv.org/abs/0907.0298).
The model-specific identities and independence checks are in the replay.

## Two-prime geometric upper bounds and the rational upper bound

Before counting, the protocol fixes the first two good primes in `5..131`
for each parent. Good reduction is checked using the semistable
discriminant multiplicities and c4 coprimality. The count worker uses at
most 120 seconds and 2 GiB per pair, one worker at a time. All twelve
pairs complete in 14.850189106 summed supervised seconds.

For every base point over Fp and Fp², the worker saves its fibre cardinality.
Independent Sage/PARI elliptic cardinalities agree on every smooth fibre;
singular counts are checked directly. Resolving the two I2 and one I4
fibres adds `q,q,3q` to the corresponding Weierstrass counts.

| u | p | #X(Fp) | #X(Fp²) | Residual quadratic trace t | NS squareclass representative |
|---|---:|---:|---:|---:|---:|
| 11 | 47 | 2964 | 4927908 | −92 | −372 |
| 11 | 53 | 3647 | 7941165 | −11 | −11115 |
| 13 | 41 | 2303 | 2857245 | −35 | −5499 |
| 13 | 47 | 2964 | 4927908 | −92 | −372 |
| 17 | 53 | 3647 | 7941165 | −11 | −11115 |
| 17 | 59 | 4427 | 12180021 | 1 | −13923 |
| 19 | 41 | 2303 | 2857245 | −35 | −5499 |
| 19 | 53 | 3647 | 7941165 | −11 | −11115 |
| 23 | 47 | 2939 | 4919973 | −23 | −8307 |
| 23 | 67 | 5596 | 20241924 | −100 | −7956 |
| 29 | 41 | 2303 | 2857245 | −35 | −5499 |
| 29 | 47 | 2939 | 4919973 | −23 | −8307 |

Subtract the known eighteen eigenvalues p. If `a` is the remaining first
trace and `b=(a^2-trace_2)/2`, Poincaré duality makes the residual polynomial

```
z^4 - a*z^3 + b*z^2 - epsilon*p^2*a*z + epsilon*p^4.
```

When `b != 0`, reciprocity forces `epsilon=+1`. In every retained `b=0`
case, `p < |a| <= 2p`; the plus case violates the Weil bound (it would
require `|a|<=p`), so `epsilon=-1`. Factoring the resulting polynomial
leaves two eigenvalues in `{p,-p}` and a quadratic
`z^2-t*z+p^2`. Its roots divided by p are not roots of unity: for a
rational quadratic of this form the cyclotomic possibilities are exactly
`t in {-2p,-p,0,p,2p}`, all excluded.

At least one prime per parent has both extra algebraic eigenvalues `-p`.
Its p-eigenspace has dimension exactly 18. Specialization of rational
divisor classes gives the **unconditional rational Picard upper bound 18**;
this step does not require the Tate conjecture.

For elliptic K3 surfaces in these finite characteristics, the Tate theorem
gives geometric Picard rank 20. Over Fp² all twenty algebraic eigenvalues
are p². Artin–Tate then gives NS discriminant squareclass `t^2-4p^2`.
For each parent the ratio of its two displayed representatives is a
nonsquare. A characteristic-zero rank-20 NS lattice would specialize with
finite index into both rank-20 lattices and force equal squareclasses.
Therefore geometric Picard rank is at most 19. Together with the previous
lower bounds, the ranks are exactly **19 geometrically and 18 rationally**.
Shioda–Tate gives exact generic MW ranks **12 and 11** respectively.

For the finite-characteristic theorem and Artin–Tate context see
[Milne, On the conjecture of Artin and Tate](https://www.jmilne.org/math/articles/1975a.html);
the two-prime specialization method is also used in
[van Luijk, K3 surfaces with Picard number one and infinitely many rational points](https://arxiv.org/abs/math/0506416).
The saved counts, rather than a heuristic rank estimate, supply the
model-specific upper bounds here.

## Saturation and the full arithmetic lattice

For any nonzero geometric section, the correction sum is at most
`1/2+1/2+1=2`; thus its height is at least 2. This excludes geometric
torsion. The rank-18 integral sublattice above can have only indices
`1,2,3,6` in the full rational fixed NS lattice, because the square of
the index must divide 468.

At T=1, complete finite elliptic groups and their cosets modulo 2 and 3
are reconstructed for the actual selected eleven sections. Both reduction
maps on their eleven-dimensional coefficient spaces are injective.
If a rational generic section divided a nontrivial word by 2 or 3,
specialization and these reduction maps would force all its coefficients
to be divisible by that prime, contradicting a nontrivial saturation
class (torsion is zero). Thus the basis is saturated at both primes,
excluding every nontrivial possible index. The full arithmetic MW height
determinant is `117/4`, and the full rational fixed NS determinant is
**−468**.

The integral base involution on this full arithmetic MW basis has fixed
rank-5 height lattice of determinant **16**, and anti-invariant rank-6
height lattice of determinant **1872**. These integral kernels, rather
than rational eigenspace bases, are used in the following glue calculation.

## The missing geometric divisor and determinant 468

Let L be the full geometric MW lattice of the rational elliptic quotient.
Its NS lattice is unimodular, its reducible-fibre root determinant is 4,
and its torsion is zero because nonzero heights are at least 1. Hence
`rank(L)=6` and `det(L)=1/4`. Its rational fixed sublattice has rank 5.
Degree-two base change doubles heights, so the fixed determinant computed
above gives `det(L+)=16/2^5=1/2`.

The constant-field Galois action fixes a codimension-one subspace. Since
it is a finite group of isometries, its nontrivial image is the reflection
on the remaining line. Write P for the primitive anti-invariant generator.
The index `e=[L:L+ + ZP]` is 1 or 2: projection to the one-dimensional
anti-invariant factor injects the glue group into its reduction modulo 2.
Orthogonality gives

```
(1/2)*h(P)/e^2 = 1/4,  hence h(P)=e^2/2.
```

The lower bound `h(P)>=1` forces `e=2`, `h(P)=2`. On the quotient,
`h(P)=2+2(P.O)-contr(P)` with `0<=contr(P)<=1`; therefore `P.O=0`
and `contr(P)=0`. In particular P is narrow.

Its pullback to the K3 is narrow of height 4. Its Shioda class
`phi(P)=P-O-2F` is an integral anti-invariant NS class of square −4.
It is primitive in the even integral K3 lattice: division by 2 would
give square −1, and a larger division cannot have even integral square.
Since the geometric and rational NS ranks differ by one, it generates
the full anti-invariant NS lattice.

The quotient index-two glue supplies a geometric section Q with
`2Q=R+P`, R rational. Pullback is injective on sections and preserves this
relation. Q cannot lie in the rational MW group plus ZP, by comparing
anti-invariant coefficients. Consequently its NS class does not lie in
the sum of the rational fixed lattice and Z phi(P). The full K3 NS
eigenspace index is therefore 2, since a rank-one anti-invariant space
also bounds it by 2. We obtain

```
det NS(X_Qbar) = (-468)*(-4)/2^2 = +468.
```

This determines the discriminant without constructing a full geometric
Gram matrix. The subsequent full NS matrix proof linked above now supplies
that matrix; a rank16 fibration remains unconstructed. Rational Picard rank 18 bounds
every Q-fibration with section by `18-2=16`, irrespective of its reducible
fibres. It does not bound rank jumps on rational specializations.

The preserved `mestre_geometric_discriminant_v1.json` correctly computed
absolute determinant 468 but printed −468 for geometric NS. That signed
record is rejected: rank 19 has signature `(1,18)` and positive determinant.
The [V2 certificate](../../artifacts/generated-results/elliptic-curves/mestre_geometric_discriminant_v2.json)
corrects the sign and preserves V1 as evidence. Arithmetic NS has signature
`(1,17)` and negative determinant. No status entry relies on the V1 sign.

## Reproduction and search decision

Copy the replay input and
[`verify_mestre_468_bundle_v2.sage`](../cas/verify_mestre_468_bundle_v2.sage)
into an empty directory and run, using Sage's Python:

```sh
sage -python verify_mestre_468_bundle_v2.sage --input mestre_468_replay_bundle_v1.json
```

The fresh-directory replay passed in **30.780057396 seconds**, with a
180-second and 2-GiB cap. It directly constructs and checks all 78 rational generic halves, and
independently checks 714 heights, generic
group relations, involution action, specialized images and finite
saturation witnesses, every saved finite fibre count, and the exact
Picard and determinant arithmetic. The mathematical implications above
use the cited theorems; this is not formal proof-assistant verification.
The source/input hashes and terminal transcript are retained in the
[portable replay record](../../artifacts/generated-results/elliptic-curves/mestre_468_portable_replay_v2.json).
No search process remains active for this cohort.

The productive next construction question is whether one of these exact
rank-18 rational lattices admits a useful rank-16 fibration over Q.
The result rules out an arithmetic rank-17 fibration on all six. Until
explicit fibrations and independent rational sections are available,
enlarging point sweeps on their current rank-11 presentations does not
resolve that construction gap. The already completed score comparison
and its validation primes remain untouched.
