# Curve 302: primary construction-recovery project

Curve 302 is the main construction-recovery target. The endpoint is an
explicit elliptic surface over `Q(t)`, its **full generic MW basis**, and a
rational parameter with an exact `Q`-isomorphism and section transport to
`E302`. Generic rank 17 is not a requirement. The additive rank-17 candidate
is evidence to test, not a prescribed generic subgroup.

**Current priority: mathematically reconstruct a structured alternative
parent using calibrated inverse-family and MW recognition.** The user has
reopened this construction route; original-provenance recovery is useful but
is not a prerequisite for an alternative parent. The MW9 result below is a
**completed baseline**, not a production-search priority merely because it
contains 302.
**No outreach:** the user has instructed that no messages be sent. The local
request draft is inactive. Mathematical inverse construction proceeds locally.

The boxed surface/basis/specialization endpoint is now complete for a
**constructed generic-rank-nine K3 parent**: the full arithmetic MW basis
and `t0=1` specialization to `E302` are explicit and certified. The original
parent remains **UNKNOWN**. This construction works for general nondegenerate
point configurations. Neither that construction nor the additive core
establishes the discoverers' method. A parent of generic rank 17–20 would
be valuable, with the surface type and its possible rank checked first;
the existing K3 has arithmetic MW rank exactly nine. Every elliptic
fibration over `Q` on this same K3 has arithmetic MW rank at most nine.

## From a recovered MW subspace to a family, 2026-09-07

The geometric control now works **without supplying the known section-image
set or parity class**. Start from the rank-twelve primitive subspace selected
by the existing graph-consensus calibration (source index 65), enumerate
vectors in that subspace, and partition them by parity. The
[packet selector](../cas/select_mw_quartic_packet.sage) retains thirteen
vectors per eligible class in numerical height order. It uses the recovered
subspace, curve and displayed points; it does not load the known generic
image list or construction parameter during selection.

The development height bound is 84, with a 200,000-representative cap and
300-second time limit. PARI at 280-bit working precision returns 171,874
signed vectors, or **85,937 representatives up to sign**, occupying all
4,096 parity classes. There are 4,092 classes with at least thirteen
representatives. This is numerical candidate selection, not an
interval-certified exhaustive canonical-height ball. The bound was chosen
during control development, so this is not a held-out statistical test.

For each thirteen-point packet, construct the degree-two coordinates below
and test for six disjoint equal-gap pairs over finite fields. At each prime
test every projective Mobius pole, including infinity. Distinct coordinate
reductions are required for exclusion; collisions and exceptional chart
values pass as unresolved to the next prime. A used coordinate cannot be
the reduced pole when all thirteen reductions are distinct: comparing its
pair difference with a disjoint pair gives a nonzero numerator in the
cross-multiplied equation. Consequently the finite pole test remains
necessary even if an unnormalized rational Mobius matrix has bad reduction.

Primes 1019, 1031 and 1033 leave respectively 360, 5 and 2 packets. Exact
rational reconstruction finds one six-pair configuration in class 1061
and none in class 3728. The positive result is precisely the 245 roots
`(0,106,344,475,594,731)` and `T=5801/10`. The replay verifies the generic
Mestre identity and matches all twelve visible covariant images, up to
sign, to the selected MW points on the target fibre. Its
[certificate](../../artifacts/generated-results/elliptic-curves/curve302_mw_packet_245_control_v1.json)
records the packet selection, finite exclusions and exact reconstruction.
The first twelve height-ordered points of the positive packet happen to be
the twelve known visible images; the thirteenth is not one of those images.
Neither that fact nor the known class was used to choose the packet.

This recovers a family from the previously recovered MW subspace. It does
not yet infer the full generic rank-twelve basis solely from that packet:
the thirteenth Fermigier section and its earlier exact subgroup calibration
remain the separate control results below.

Before switching to this subspace route, a full exact XOR-convolution pilot
counted triple decompositions of all rank-twenty parity classes. Ranking
previously absent classes by decreasing count put the known control class
199,356th: 199,355 scored higher and 474 tied its count of 1,459. Popularity
of triple decompositions alone therefore did not repair selection. Local
pilot inputs and outputs remain under
`artifacts/local/elliptic-curves/mw-parity-selection/`.

Replay the packet control with:

```
sage -python elliptic-curves/cas/select_mw_quartic_packet.sage --check
```

The [target driver](../cas/search_curve302_mw_packets.sage) also checks that
its optimized array extraction reproduces every control packet exactly:

```
sage -python elliptic-curves/cas/search_curve302_mw_packets.sage --target 245
```

The subsequent bounded 302 protocol takes the minimum-determinant retained
rank-seventeen finalist, the existing additive core, and uses height bound
190, at most four million vector representatives, thirteen points per
eligible parity class, primes 1019/1031/1033/1039, at most 32 exact survivor
packets, one worker and 900 seconds. The height bound is a computational
coverage choice, not a prediction of generic rank or a completeness bound.
Checkpoints are under `artifacts/local/elliptic-curves/curve302-mw-packets/`.
An initial decimal-input conversion failure is retained separately; it
produced no exclusion.

The [302 screen](../../artifacts/generated-results/elliptic-curves/curve302_mw_packet_rank17_h190_v1.json)
completed in approximately 202 seconds. It enumerated **2,724,221 vector
representatives up to sign**, occupying 131,071 of the 131,072 parity
classes. Of these, **130,706 classes** supplied thirteen points. The four
modular stages excluded respectively 113,741, 15,437, 1,378 and 129 packets.
All 21 remaining packets were unresolved because of modular chart issues
or coordinate collisions, rather than positive finite-field configurations;
exact rational reconstruction excluded all 21. There were **no six-pair
hits and no unresolved packets** in the completed screen.

This excludes the pair configuration in those specific packets. In
particular it does not exclude other point selections in the same classes,
the 366 classes that did not supply a packet, other subspaces, or every
Mestre parent. The rank-seventeen candidate is a tested container, not a
generic-rank requirement. The earlier forced-rank-seventeen 245 candidate
intersects its actual rank-twelve family space in only rank nine; that is
a further reason not to infer global provenance from this 302 core.

The exact [packet inputs](../../artifacts/generated-results/elliptic-curves/curve302_mw_packet_rank17_h190_inputs_v1.npz)
store every selected integer vector. The
[certificate replay](../cas/verify_curve302_mw_packets.sage) checks their
hashes, all parity classes, the primitive rank-seventeen embedding, every
finite exclusion and all 21 rational exclusions. It avoids repeating the
numerical height enumeration and shares the original chart/detector code.

```
sage -python elliptic-curves/cas/verify_curve302_mw_packets.sage --check
```

The construction objective remains open. The new gate is now passed for
the recovered control subspace. Further construction work must change the
tested packet or subspace, or use a different family architecture; rerunning
these frozen packets cannot add information.

## Recovering quartic pairs from MW images, 2026-09-07

The next geometric recognizer now passes an exact supplied-point control.
Authority: `EC-CURVE302-MW-QUARTIC-PAIR-CONTROL`. The
[checker](../cas/recover_mw_quartic_pairs.sage) takes the thirteen known 245
section images in the displayed MW coordinates, scrambles their order, and
withholds their quartic coordinates and pair labels. Its
[certificate](../../artifacts/generated-results/elliptic-curves/curve302_mw_quartic_pairs_245_control_v1.json)
recovers the unique six-pair configuration, the roots
`(0,106,344,475,594,731)`, and `T=5801/10`. It then checks the generic
Mestre polynomial identity, a rational isomorphism of the recovered fibre,
and equality of all twelve visible covariant images with the selected input
images up to sign. This calibrates geometric recognition **given the correct
MW images**, not blind selection of those images or a full generic basis.

Here is the coordinate issue resolved by that control. Covariant images of
the quartic points have the form `Q_i=2R_i-C`. They lie in one class modulo
twice the displayed MW group. Choose `R_0=O`, so `C=-Q_0`, and compute
`R_i=(Q_i-Q_0)/2` using the integer MW coordinates. On a short Weierstrass
model the function

```
z(R) = (y(R)+y(C))/(x(R)-x(C))
```

has poles at `O,C` and involution `R -> C-R`. It is therefore a suitable
degree-two coordinate, up to a rational Mobius transformation. The checker
verifies this recovery on the control; it fails closed at exceptional
coordinates not implemented in this control. Directly matching the
Weierstrass x-coordinates of `Q_i` would omit both the halving and the Mobius
freedom.

After choosing a finite chart, write four abscissas as `a,b,c,d`. Equality
of two transformed pair differences determines the Mobius pole `k` through

```
(a-b)(c-k)(d-k) = +/-(c-d)(a-k)(b-k).
```

This is quadratic in `k`; also test the affine case `k=infinity`. Every
six-pair solution contains two disjoint pairs, so this enumerates all
possible poles for the supplied thirteen coordinates. Here 4,290 equations
leave 17 candidate poles including infinity, and exactly one pole has six
disjoint pairs with a common separation. All matchings are enumerated.
The pair centers recover the roots; half the separation recovers `T`, up
to the usual affine scaling and reflection. The subsequent polynomial and
point checks are essential: equal gaps alone do not certify a Mestre parent.

The same replay diagnoses the old point selector exactly. The thirteen
control image rays occur **zero times** in either the 1,928-vector short
ball or the 6,798-vector sparse pool. More decisively, **none of the 1,928
ball seeds has their parity class**. For every such seed `c`, no vector
`c+2z`, for any integral displayed-basis vector `z`, can equal a known
control image, even up to sign. Enlarging only the old shift bound cannot
repair that failure. This does not exclude alternative presentations of
the family or other parent families.

A cheap retrospective parity audit points to a specific change of selector:
the missing class is the sum of no two distinct ball-ray classes, but is
the sum of **1,459 triples** of distinct ball-ray classes. The replay checks
all 1,857,628 unordered pairs against a parity lookup and stores witness
triples. This is a diagnostic using the known control class, not a blind
method for prioritizing unknown classes on 302.

The packet control above now reaches the missing class by enumerating inside
the previously recovered MW subspace. Triple sums remain a possible source
of additional classes when such a subspace is unavailable; count popularity
alone failed the development pilot. Successful configurations must still
give a surface, a certified generic MW basis and exact specialization. This
supplied-image control alone asserts no alternative 302 parent or exclusion.

```
sage -python elliptic-curves/cas/recover_mw_quartic_pairs.sage --check
```

The declared control limit is thirteen points, 4,290 quadratic equations,
1,857,628 parity comparisons, one worker and 180 seconds. Replay takes
approximately one second in the local Sage runtime.

## Calibrated inverse Fermigier recognition, 2026-09-07

The first executable inverse run is complete. It recovers the known 245
family and its actual transported rank-twelve section subgroup, while
excluding 302 from the declared height-eight parameter box. The authority
is `EC-CURVE302-INVERSE-FERMIGIER-H8`; this narrows the parent objective
without completing it.

Fermigier's two rational parameters `u,v` give six labelled roots. The
constructor builds `q(X)=product(X-a_i)` and the unique monic degree-six
`g(X,T)` with

```
q(X-T)q(X+T) = g(X,T)^2 - T^2 R(X,T),   degree_X R <= 4.
```

The identity and quartic condition are checked symbolically for each root
configuration. If `I,J` are the binary-quartic invariants, its Jacobian is
`y^2=x^3-27I(T)x-27J(T)` and

```
j(T) = 6912 I(T)^3 / (4 I(T)^3-J(T)^2).
```

The even parameter dependence permits exact factorization in `Z=T^2`.
For each root configuration, cancel the common factors of the j-map and
compare it with the target's exact j-invariant. Modular rejection checks
all square `Z` residues and the point at infinity; degree loss or a zero
reduction cannot be counted as an exclusion. Surviving comparisons are
factored over `Q`. Only rational square linear roots can give rational
finite `T`; both signs are checked, followed by an exact `Q`-isomorphism
to the target. Singular raw models and matching infinity limits remain
unresolved unless local minimalization settles them. None occurs in the
completed 302 run.

The frozen box is `max(|numerator|,denominator)<=8` for each of `u,v`, in
lowest terms: 87 rationals and **7,569 ordered pairs**. Of these, 383 have
repeated roots and lie outside the declared nondegenerate six-root locus.
Affine normalization and reflection leave **2,578 root configurations**.
This is not asserted to be a complete fibration-equivalence classification.
There is **no height cutoff on rational T**. The fixed prime list is
`101,103,107,109,127,131,137,139,149,151,157,163`, with one worker,
900 seconds per run and per-configuration checkpoints.

The [245 recognition certificate](../../artifacts/generated-results/elliptic-curves/curve302_inverse_fermigier_245_h8_v1.json)
finds the roots `(0,106,344,475,594,731)` and `T=+-5801/10`. Its equivalent
parameter presentations include `(u,v,T_native)=(3/2,2,5801/160)`.
The known parameter was used only in the post-search control evaluator.
The [MW control replay](../../artifacts/generated-results/elliptic-curves/curve302_inverse_fermigier_245_mw_control_v1.json)
then reconstructs the twelve visible quartic points and the extra affine
section from the recovered parameters, verifies the generic identities,
and transports their covariant images to the public 245 model. Eleven
visible images and the extra image give exactly the retained integer
rank-twelve subgroup. Its Smith factors are `1,2,2,2,2,2,2,2,2,2,2,2` and
its index in the displayed primitive closure is **2,048**. Fresh complete
finite-group quotients modulo doubles independently verify the public
twenty-point independence needed for this comparison. This checks the
actual subgroup, not just its rank or primitive closure; it is not a new
generic-saturation theorem for 245.

The [302 result](../../artifacts/generated-results/elliptic-curves/curve302_inverse_fermigier_302_h8_v1.json)
has **2,577 modular exclusions and one exact exclusion**. The sole modular
survivor has roots `(0,24662,468768,473773,502957,581750)` and an irreducible
degree-twelve comparison in `Z`, so it has no rational `T`. There are no
unresolved configurations and no 302 match. Comparing canonical root tuples
with the earlier 2,333-six-root screen gives 30 overlaps and **2,548 newly
tested root configurations**. The old census and this parameter box are
different finite spaces.

This is a calibrated inverse recognizer with a bounded negative result,
not a solution of the full rational moduli problem. The next construction
step must reach parameters outside this box or impose genuine section
incidence constraints on the moduli equations. Another search for an
independent section on the completed MW9 K3 cannot achieve that step.

Replay:

```
sage -python elliptic-curves/cas/inverse_fermigier_parent.sage --target 245 --height 8 --check
sage -python elliptic-curves/cas/verify_inverse_fermigier_245_control.sage --check
sage -python elliptic-curves/cas/inverse_fermigier_parent.sage --target 302 --height 8 --check
```

Local logs and checkpoints are in
`artifacts/local/elliptic-curves/inverse-fermigier-parent/`. An initial
Fraction/Sage coercion error and a Smith-form return-type error are retained;
neither produced a mathematical exclusion. The generic-recognition script
also preserves unresolved local-minimalization cases instead of rejecting
a potentially removable singular model.

## Full arithmetic MW9 certificate

The [full-basis certificate](../../artifacts/generated-results/elkies-k3-curve302-full-mw-basis-v1.json)
stores `A(t),B(t)`, all nine pairs of rational-function Weierstrass
coordinates, the height matrix, and the exact specialization map. Its
authority is `EC-K3-CURVE302-FULL-MW9` in
[`MATH_STATUS.json`](../../MATH_STATUS.json). It completes the earlier
saturation result below with an unconditional generic rank upper bound.

The surface `y^2=x^3+A(t)x+B(t)` has degrees eight and twelve and 24
geometric `I1` fibres. It has good reduction at 47. The Newton tetrahedron
has vertices `(0,0,0),(12,0,0),(0,3,0),(0,0,2)`, unique interior point
`(1,1,1)`, and primitive facets at distance one from that point. Thus it
gives a toric K3 compactification. Nondegeneracy on all faces follows from
the checked squarefreeness of `B` and the degree-24 discriminant, nonzero
endpoint coefficients, and smooth fibres at zero and infinity. The smooth
minimal model is the same K3 as the elliptic presentation.

[ToricControlledReduction](https://github.com/edgarcosta/ToricControlledReduction),
at commit `74cda9e8148cd8e9a3928fc15a558c9a70b67cc1`, computes its
20-dimensional primitive factor, with Hodge vector `[1,18,1]`. The method
is [Costa–Harvey–Kedlaya controlled reduction](https://arxiv.org/abs/1806.00368).
At 47 its characteristic polynomial is

```
P(T) = (T-47)^9 (T+47) F(T),
F(T) = T^10 - 32 T^9 + 3149 T^8 - 53016 T^7
       + 8098194 T^6 - 78074896 T^5 + 17888910546 T^4
       - 258701167896 T^3 + 33943749071021 T^2
       - 761961173176352 T + 52599132235830049.
```

Exact factor division, the weight-two functional equation and a certified
real-algebraic Weil-circle test pass. Exhaustive cyclotomic division of
`47^-20 P(47 Z)` finds only `(Z-1)^9 (Z+1)`: orders through 800 suffice
by `phi(m)>=sqrt(m/2)`. The residual factor has no root of unity. Every
input and output field is bound to the exact reduced model; raw backend
input, output, source commit and original binary hash are retained.

Independent fibrewise PARI counts over `F47` and `F47^2`, including infinity,
give primitive Frobenius power sums `408,16816`, agreeing with the first
two polynomial moments. They give surface point counts `2712,4900916`.
These independent moments supplement the complete controlled-reduction
calculation; they do not by themselves determine its degree-20 polynomial.
The complete polynomial has also been recomputed with the pinned backend.

Here is the rank argument, without using the Tate conjecture in the converse
direction. Good-reduction specialization and the cycle class map inject
the Néron–Severi group into the reduction's cohomology. Each divisor defined
over `Q` contributes a `47`-eigenvector. The primitive factor has exactly
nine such eigenvalues; its complement in the K3's 22-dimensional `H^2`
has dimension two. Consequently

```
rank NS(X_Qbar)^Gal(Qbar/Q) <= 9+2 = 11.
```

The fibre, zero section and nine independent rational sections give the
matching lower bound eleven. Shioda–Tate, with no reducible fibres, therefore
gives **rank MW(Q(t))=9**. The previously certified geometric saturation and
trivial torsion turn the nine sections into a **full arithmetic basis**.
Neither BSD nor GRH enters. The cyclotomic calculation also gives geometric
Picard rank at most twelve, hence geometric MW rank between nine and ten;
the full geometric rank remains **UNKNOWN**.

The quadratic source twist now has arithmetic rank exactly one, since the
base-change rank is `8+rank(twist)`. Moreover every elliptic fibration over
`Q` on this K3 consumes at least the two rational fibre/zero classes, so
its arithmetic MW rank is at most nine. This excludes recovering generic
rank 17–20 by fibration changes on this particular surface. It does not
exclude a different parent surface through the same fibre.

The expanded coordinate numerator/denominator degrees are `(4/0,6/0)` for
the first seven sections, `(8/4,12/6)` for the eighth, and `(6/2,9/3)` for
the ninth. At `t=1`, set `X=x/D(1)^2`, `Y=y/D(1)^3`; the map to the public
302 model is `xp=(X-15)/36`, `yp=(Y/108-xp-1)/2`. All nine exact generic
point equations and their specialized public group-law identities replay.
Thus the explicit family is calibrated on a fibre of certified rank at
least 31, with 22 independent directions beyond its generic image. No new
rank claim for a rational fibre is made.

## Certified improvement: full source basis and primitive K3 span

The [new certificate](../../artifacts/generated-results/elkies-k3-curve302-section-saturation-v1.json)
completes the rational pencil's generic rank-eight MW basis and enlarges
the constructed K3's nine-section group by index three. The enlarged
rank-nine group is saturated even in the geometric MW group. The subsequent
Frobenius certificate above proves it is the full arithmetic MW group.

The status authority is `EC-K3-CURVE302-SECTION-SATURATION` in
[`MATH_STATUS.json`](../../MATH_STATUS.json). The original
[nine-direction certificate and geometric proof](ICARM_CURVE302_INVERSE_FIBRATIONS_2026-09-06.md)
remain dependencies; their coordinates and historical assertions are preserved.

Let the source pencil be `F0+u F1=0`, with zero `P1`, and write

```
C1,...,C7 = P2-P1,...,P8-P1,
C8 = R-P1,       R=-(P1+...+P8) on E302's short model.
```

Let `Q12(u)` be the third intersection with the line through `P1,P2`.
If `lambda` parametrizes that line as `P1+lambda(P2-P1)` in affine
coordinates, its restriction is a cubic with roots 0 and 1. Its third root
is exactly `-coefficient(lambda^2)/coefficient(lambda^3)-1`. This gives
explicit rational coordinates for `Q12` in the certificate. Hereafter `Q`
denotes its pointed section `Q12-P1`.

Write `h` for the hyperplane divisor class minus `3P1` on the generic cubic.
The nine basepoints cut out the second cubic, hence their sum as divisors
is linearly equivalent to three hyperplanes. The secant cuts out one
hyperplane. Consequently

```
C1+...+C8 = 3h,          Q = h-C1.
```

Thus replace `C8` by `Q`, retaining `C1,...,C7` and the ninth moving section
`S`. This changes basis by a rational matrix `T` which is the identity
except for row eight, which is

```
(-2/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 0).
```

Conversely `C8=3Q+2C1-C2-...-C7`, so the old group has index three in
the new one. Applying `T` to the old height matrix gives

```
H = [4 2 2 2 2 2 2 2 3]
    [2 4 2 2 2 2 2 4 3]
    [2 2 4 2 2 2 2 4 3]
    [2 2 2 4 2 2 2 4 3]
    [2 2 2 2 4 2 2 4 3]
    [2 2 2 2 2 4 2 4 3]
    [2 2 2 2 2 2 4 4 3]
    [2 4 4 4 4 4 4 8 5]
    [3 3 3 3 3 3 3 5 6],       det(H)=512.
```

The source height matrix is half the upper-left eight-by-eight block. It
is positive definite, even, integral, and has determinant one. The rational
elliptic surface has twelve `I1` fibres, Picard rank ten, and geometric MW
rank eight by Shioda–Tate. Its MW lattice is integral and has no torsion.
The determinant-one subgroup therefore has index one: these eight rational
sections form the **full geometric and arithmetic generic MW basis**.
This supplies a complete rank-eight baseline surface, basis, and `u=0`
specialization to 302. It is not provenance recovery.

For the K3, every nonzero geometric section has height
`4+2(P.O)>=4`, because there are no reducible fibres. The MW lattice is
even integral and torsion-free. Any proper overlattice of the rank-nine
group inside its rational span would have index a power of two, since
`det(H)=2^9`, and would contain an order-two first step. Such a step is
represented by a nonzero `v` in `{0,1}^9` satisfying

```
H v = 0 mod 2,             v^T H v = 0 mod 8.
```

There are exactly 71 such cosets. For **every** one, the certificate supplies
an integer vector `w=v mod 2` with `w^T H w=8`. Adjoining `v/2` would
therefore also adjoin `w/2`, a section of height two, a contradiction.
All 512 parity vectors and all 71 witnesses are checked exactly. This is
a complete finite overlattice proof, independent of whether the short-vector
routine found every short vector. Odd-prime saturation is already forced
by the determinant. The resulting rank-nine group is primitive in the full
geometric MW group, whose rank is now bounded between nine and ten.

The height and rational-surface facts used here are the standard results in
[Schütt–Shioda, sections 8 and 11](https://arxiv.org/html/0907.0298v3#S11).

## Explicit transport and specialization

The certificate contains a Weierstrass model over `Q(u)`, all eight generic
basis points on that model, and the rational cubic-to-Weierstrass map.
The map uses a tangent frame at `P1`, rather than finding rational flexes.
Its substitutions are verified as polynomial identities; each of the eight
Weierstrass point equations is checked over `Q(u)`. At `u=0`, all eight
points select and verify the same isomorphism to the short 302 model,
resolving the possible sign. The fixed short-to-public change is retained.

The source equation is the pinned compact Jacobian. Seven basis sections
have polynomial coordinates of degrees `(2,3)`; the eighth has a small
rational denominator. The K3 equation and sections are explicitly obtained
by substituting `u=N(t)/D(t)` and multiplying point coordinates by
`D(t)^2,D(t)^3`. Its ninth section is the same cubic map applied to the
pinned `L(t)`, followed by these scalings. These formulas are stored as
rational-function coefficients and a composition rule. At `t=1`, first
undo the `D(1)` scalings and then apply the source fibre isomorphism.
All nine specialized section transports are checked.

At this fibre, `Q=-2P1-P2`, so the new specialization rows are

```
P2-P1, ..., P8-P1, -2P1-P2, P9-P1.
```

Their Smith factors are `1,1,1,1,1,1,1,1,3`. With the independently
certified public points, specialization is injective on this group and
the displayed public quotient is `Z^22 + Z/3`. Generic saturation therefore
does not imply saturated specialization. No additional independent direction
on the K3 or on the rational fibre was discovered.

## Original construction record: separate open route

1. Inspect available original records for the family or generation procedure,
   exact parameter/input values, and generators before basis reduction.
   No messages may be sent. Do not assume a one-parameter elliptic surface
   or generic rank seventeen. The catalogue source audit below distinguishes
   the submitted witness from the still-unknown pre-upload generators.
2. Preserve any supplied record in its original form with attribution and
   hash. Replay its model and parameter to the exact public curve, transport
   its original points, and compare their span with the published basis.
   Recover a generic basis and rank only when the supplied construction
   actually defines a family to which those questions apply.
3. Use that provenance to investigate the source of the exceptional rank
   and possible improvements. Retain the MW9 theorem and replay as a
   completed baseline. Its successful 302 specialization alone does not
   justify promoting it to a major production search. Additional imposed
   parents are not a substitute for provenance. The calibrated mathematical
   inverse-construction route above proceeds independently.

The arithmetic rank gate for the baseline is closed. Every elliptic
fibration over `Q` on that K3 has MW rank at most nine. This remains useful
evidence, independently of the original construction's provenance.

### Outreach preparation, 2026-09-07

The [catalogue](https://elliptic-rank.icarm.cloud/curve/302) credits Claude,
Levent Alpöge and Ava Howell. Levent's [linked CV](https://alpo.ge/cv.pdf)
publishes `alpoge@fas.harvard.edu`; his [homepage](https://alpo.ge/) retains
the matching abbreviated address. Ava's [website](https://avahowell.me/pages/about/)
has a contact link, but its protected address was not recoverable in this
session. No address has been guessed or taken from a namesake's profile.

The [prepared request](CURVE302_PROVENANCE_REQUEST.md) addresses Levent and
asks him to include Ava if she holds the relevant record. **Not sent and
inactive:** the user subsequently instructed that no messages be sent.
Connecting a mail service would not authorize delivery. No reply or
recovered provenance is claimed.

### Catalogue source audit, 2026-09-07

The [source-review evidence](../../artifacts/generated-results/elliptic-curves/curve302_catalogue_provenance_audit_v1.json)
pins the public JSON and eight source blobs. The live model and all 31
ordered points exactly match the retained certified inputs. This is a
provenance audit, not a new mathematical certificate.

Two code revisions were inspected: current commit
`15957f7d7d2f539ff9a457264544caac9e7567bf` and the latest repository commit
preceding the recorded submission,
`4c9ef36265c3cbc10f5da84fcade3f4bbe6854cc`. The actual deployed revision at
submission time has not been established.

- **The inspected application does not LLL-reduce its stored witnesses.**
  Its successful [verifier output](https://github.com/icarm/elliptic-rank/blob/4c9ef36265c3cbc10f5da84fcade3f4bbe6854cc/src/verify.ts#L778)
  maps each input point through the minimal-model change, preserving order.
  [Storage](https://github.com/icarm/elliptic-rank/blob/4c9ef36265c3cbc10f5da84fcade3f4bbe6854cc/src/store.ts#L234)
  serializes those points. The independence algorithm's internal point
  replacements do not become the stored witness. Thus reduction performed
  before upload remains a distinct unknown; the catalogue is not evidence
  that it happened, nor that it did not.
- **Empty history is inconclusive.** The
  [contribution-log migration](https://github.com/icarm/elliptic-rank/blob/18bf9b3ce79ea1bf1e58c983dadcbcd2b65cfc48/migrations/0012_curve_events.sql)
  was introduced on September 2, after 302's August 23 submission, without
  backfilling earlier events. The public `history: []` cannot certify that
  the present witness was the first submitted list. The recorded update
  at `20:16:24` could reflect a rank improvement or a prime/conductor update;
  the public export does not identify which.
- **The public store cannot reconstruct the original coordinate model.**
  Its inspected schema and export retain the minimal model and current
  witness, without the original input model, its transformation tuple,
  family parameter, or generation procedure. This says nothing about
  private logs or backups. No such records were accessed.

A bounded commit-message search found only a coordinate-length limit change
mentioning a doubled point on a rank-31 curve, not a generation recipe.
The public issues about subfamilies and literature also supplied no 302
construction. Public Zulip messages and the gist landing page were not
readable in the available view and remain unaudited.

This narrows the useful missing evidence to **the pre-upload construction
record**, rather than an assumed catalogue basis-reduction transform. It
does not identify a parent. Preserve point ordering for comparison with any
future source record, while retaining the earlier negative calibration of
the first-seventeen boundary.

For replay, use the recorded git commits and `git show COMMIT:PATH`, compare
each blob's SHA-256 with the evidence file, and compare the embedded public
snapshot with `icarm_curve302.py`. The local source checkout is under
`artifacts/local/tools/icarm-elliptic-rank-provenance/`; no code from that
checkout was executed. Only public read operations were used.

Large neighbour, descent, or specialization campaigns require a separate
finite protocol with a mathematical gate, limits, checkpoints, and replay
inputs. This step runs no such campaign and does not feed record data into
the separate target-free MW16 experiment.

## Parent search, 2026-09-07

A fresh public-source search found no parent formula or specialization
parameter in the inspected records. The [302 entry](https://elliptic-rank.icarm.cloud/curve/302)
still supplies the curve, points and attribution; its
[commentary history](https://elliptic-rank.icarm.cloud/curve/302/commentary-history)
has one edit, containing the attribution and conditional exact-rank claim.
Neither [Alpöge's publication page](https://alpo.ge/) nor
[Howell's website](https://avahowell.me/pages/about/) or the inspected
[public repository listing](https://github.com/avahowell?tab=repositories)
provided a construction. Exact-coefficient and author/topic searches also
produced no usable family. The public Zulip landing page did not expose its
messages in the available web view, so those discussions were not audited.
This is a dated search result, not a claim that no public account exists.

A new finite mathematical search tested conics beyond the earlier line
screen. Fix the existing nine-basepoint cubic pencil. For each of the
`binomial(9,4)=126` four-basepoint subsets and each of `P9,...,P31`, take
the unique conic through those five points when it exists. A smooth conic
has total cubic intersection degree six; its four prescribed basepoints
leave a degree-two map to the pencil parameter. The frozen space is
**2,898 conics**, with prime budget `101,103,107,109,127`; no other pencils,
point combinations, or parameter sweeps are included.

At each usable prime, solve the five conic conditions exactly, require rank
five and a nonsingular conic, and parametrize it by lines through its first
basepoint. Substitution into `F0,F1` and removal of their common factor must
give a degree-two map. For residual polynomials `g0(s),g1(s)`, take the
discriminant in `s` of `g0(s)+u*g1(s)`. Compare its full `F_p(u)` squareclass
with the discriminant of the pinned `N(t)-uD(t)`, including the relative
constant. A differing branch divisor or nonsquare constant excludes equality
of the two covers over `Q(u)`. Unusable reductions are deferred to the next
prime, never counted as exclusions. The conic parametrization identities
are checked, as is the required square discriminant at `u=0` supplied by
the selected rational point.

All **2,898** rows are excluded: `1832,916,125,12,13` respectively at the five
primes. No row remains unresolved. The
[certificate](../../artifacts/generated-results/elkies-k3-curve302-conic-parent-overlap-v1.json)
retains each row's subset, point, prime, conic and branch coefficients; exact
replay reconstructs every witness from the rational inputs.

This excludes only these conics as sections after the **existing** quadratic
base change. It does not exclude their other quadratic covers as different
K3s, common covers among conics on another pencil, higher multisections,
combinations of public points, or an original parent. It proves no generic
rank upper bound by itself; the later Frobenius certificate supplies that
bound for this K3. Original-parent recovery remains open. The authority for
this bounded result is `EC-K3-CURVE302-CONIC-COVER-SCREEN`.

## Replay

```
sage -python elkies-k3/scripts/construct_curve302_nine_direction_k3.sage --check
sage -python elkies-k3/scripts/certify_curve302_section_saturation.sage --check
python3 elkies-k3/scripts/verify_curve302_saturation_witnesses.py
sage -python elkies-k3/scripts/search_curve302_conic_parent_overlap.sage --check
sage -python elkies-k3/scripts/certify_curve302_full_mw_basis.sage --check
# Repeat the complete external Frobenius computation as well:
sage -python elkies-k3/scripts/certify_curve302_full_mw_basis.sage --check --recompute-frobenius
```

The new computation is one exact section construction, eight symbolic point
checks, one fibre transport, and a fixed 512-element parity enumeration.
The Python command independently replays the entire finite lattice proof using
only Python's standard library, without Sage or a short-vector enumerator.
It does not replace the geometric and symbolic section checks.
Local construction/replay logs are retained under
`artifacts/local/elkies-k3/curve302-construction-recovery/`. The initial
unsupported multivariate division attempt and two Sage interface/coercion
failures are retained; the final run uses exact monomial division and direct
rational coordinate formulas. No failed run is a mathematical exclusion.
Conic search and replay logs are under
`artifacts/local/elkies-k3/curve302-parent-search/`.

The full-basis replay expands all nine coordinates and recomputes both
finite-field moments. The optional external replay uses the pinned backend
checkout under `artifacts/local/tools/` and a 300-second timeout. A fresh
clone at the recorded commit must be built there if the local executable
is absent. The exact input/output is embedded, so ordinary `--check` does
not require that checkout. It still relies on the retained full Frobenius
computation as evidence and is not an independent implementation of it.
The frozen Frobenius budget was three primes `47,53,127`, 300 seconds each,
one worker, stopping when rank nine closed. Only the first prime needed a
complete Frobenius computation; its initial run took about seven seconds.
Logs, unused inputs, and pilot moment checks are retained under
`artifacts/local/elkies-k3/curve302-full-mw/`.
