# Exact norm-twelve bisection-character exhaustion (2026-09-03)

<!-- status-consumer: EC-K3-R17-NORM12-11952-COMPLETE-BISECTION-CHARACTER-EXHAUSTION 6a83ca559fed9c2b -->
<!-- status-consumer: EC-K3-R17-NORM12-103B2-DIRECT-BISECTION-CHARACTER-EXHAUSTION c79f496f7a669ae4 -->

## Status

This note records an exact negative result for the proposed alternate-Q80
equal-cover and product-character constructions.  It does **not** prove a new
rank-19 or rank-20 family.

For the direct `norm12-orbit-11952` alternate-Q80 fibration:

- the canonical equation has `24 I1` fibres and a saturated rank-17 section
  basis;
- the 121 inherited height-four curves give 121 distinct quadratic
  squareclasses;
- all 39,147 section-translation classes in the complete rootless rational
  bisection frame have exact equations and verified lifts, and their 39,147
  squareclasses are pairwise distinct;
- no product of two complete-frame characters is a third complete-frame
  character, and none matches the committed old-base character catalogue.

The prescribed fallback was also completed.  The hidden
`norm12-orbit-103b2` marking now has an explicit `24 I1` equation and a
saturated rank-17 basis in the published-R17 frame.  Its 82 inherited
height-four bisections and all 39,120 complete-frame rational bisection classes
again give pairwise distinct extensions and no three-character closure.

Thus the requested two-by-two anti-invariant height Gram never arises: there
is no equal-extension bucket on either complete frame.  Likewise there is no
three-character input from these frames with which to certify the requested
`V4` rank-20 base.

The exact boundary is important.  This exhausts the smooth rootless
`(-2)`-bisection mechanism represented by the two finite frame tables.  It
does not exclude higher-arithmetic-genus degree-two curves, singular rational
bisections outside those tables, or unrelated quadratic characters carrying
new twist sections.

## Literature control

Elkies's current paper proves rank 18 from a rational quadratic section and
rank 19 from the compositum of two distinct quadratic base changes over a
positive-rank elliptic curve; see [Elkies, *An elliptic K3 surface over
Q(t) with Mordell--Weil rank 17, I*](https://arxiv.org/abs/2608.25406),
Theorem 3 and Section 3.  That positive construction is the published-R17
two-character genus-one route.  It neither supplies two independent new
directions on one rational quadratic cover nor supplies the third nonzero
character needed for rank 20.  The computations below test exactly those two
stronger requirements after changing to the alternate-Q80 and hidden
`0x103b2` fibrations.

## 1. Canonical alternate-Q80 model

The equation source is
[`elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json),
SHA-256
`76c54483c93c7090def42a8dad256838eb9510cd8479d07c5e3123eefa5cfe66`.
The compiler verifies:

- `(deg A,deg B,deg Delta)=(8,12,24)`;
- squarefree discriminant, coprime to `A`, and a smooth fibre at infinity,
  hence `24 I1`;
- determinant-948 rootless alternate-Q80 frame;
- seventeen exact rational sections whose frame-coordinate determinant is
  `-1`, so the displayed basis is saturated.

```bash
sage -python elkies-k3/scripts/compile_r17_norm12_orbit11952_qq.sage --check
```

## 2. Inherited cheap set

The exact inherited compiler transports all old height-four sections and
selects precisely those whose divisor class has degree two over the new
fibre.  On the alternate frame this gives 121 classes.  It recovers the two
requested seeds

\[
B_1=-P_3+P_8,\qquad B_2=-P_1-P_2-P_7.
\]

For each class it constructs the quadratic equation over `QQ(u)`, completes
the square to `s^2=q(u)`, retains the rational constant squareclass, and
verifies the lifted point and its conjugate coefficientwise on the direct
Weierstrass equation.  The resulting anti-invariant direction has height 12.

The cover artifact has SHA-256
`6ba1db53ed02de18910456e75dcdad55e87205c7355986080dab1955d83050b9`.
Exact hashing returns `121` extensions and `0` collisions.  The collision
artifact has SHA-256
`b3b5befad90367fb7cabd01fcb96c3e972ce0d3b71aa592c0607419a8f48e35e`.

The 7,260 unordered products are also all distinct.  None equals another
inherited character, and none formally equals one of the eleven committed
rank-28 quartics or `q_103b2`.  The regenerated product artifact has SHA-256
`fe97c92d5ea8609dd31f2e68a9cfe157ee93530d6e0b1e35f4714e59b0d6409d`.

```bash
sage -python \
  elkies-k3/scripts/construct_r17_norm12_11952_inherited_bisections.sage
.venv/bin/python elkies-k3/scripts/hash_bisection_extensions.py \
  --compact \
  --input artifacts/generated-results/elkies-k3-r17-norm12-11952-inherited-bisection-covers-v1.json \
  --output artifacts/generated-results/elkies-k3-r17-norm12-11952-inherited-bisection-collisions-v1.json
.venv/bin/python \
  elkies-k3/scripts/analyze_r17_norm12_11952_inherited_products.py \
  --source-label norm12-orbit-11952
```

## 3. Complete alternate frame

Exact norm-ten enumeration gives 39,147 section-nonnegative translation
classes in the alternate-Q80 frame.  The priority compiler selects a cheap
equation representative in each class.  The equation compiler then constructs
the cover and verifies the lifted section for every class, using a reciprocal
chart when required.  Merging is accepted only after matching the complete
priority table by orbit mask and frame vector.

The merged equation artifact has SHA-256
`fbf979bfe7d92528405c62330a80dbfd7742dc27c41a0426dcd4014f6865c8ce`.
The independent squareclass pass returns

```text
bisections=39147  extensions=39147  collisions=0
```

and its compact collision artifact has SHA-256
`e5cc1207b6d60b3ad9aefef44afc1634f3f808e964c2569e54b8bed8b7bef178`.
This is the exact obstruction to Step 4: no two distinct frame classes define
the same quadratic extension, so there is no two-dimensional anti-invariant
height space to certify.

## 4. Complete product-character closure

The complete closure checker does not enumerate 766 million pair products.
It factors every canonical branch exactly and uses the following exhaustive
support argument.

Every one of the 39,147 alternate branch polynomials is irreducible quadratic
over `QQ`.  Hence its squareclass contains exactly one finite irreducible
polynomial atom, and these atoms are pairwise distinct.  The product of two
distinct characters contains two such atoms, whereas every character in the
complete frame contains one.  Therefore

\[
q_iq_j\ne q_k\quad\text{in }\mathbf Q(u)^*/\mathbf Q(u)^{*2}
\]

for all distinct frame characters.  The checker still compares all rational
constant atoms and performs exact lookup against all twelve committed old-base
characters.  It finds no formal match.  The closure artifact has SHA-256
`e72e191e1f8c9df0579dfd258deae9ff733b0a5705336f1f0f808a50b98e59b1`.

The Neron--Severi compatibility gate is separate from the coefficient lookup.
The rank-28 source curve and the `q_103b2` source curve have degrees 8 and 9
over the alternate `u`-line.  Neither is a degree-two character curve for the
alternate fibration, so a formal replacement `t -> u` could not transfer its
known twist section even if the coefficient strings matched.

```bash
.venv/bin/python \
  elkies-k3/scripts/analyze_r17_norm12_complete_character_closure.py \
  --source-label norm12-orbit-11952
```

## 5. Exact `0x103b2` fallback

The generalized direct compiler uses

\[
w=(0,-1,1,1,-1,1,-1,0,1,1,0,1,1,-1,-1,0,-1),
\qquad D=(3,2,w).
\]

It constructs an exact `24 I1` Weierstrass model.  Fifteen transported old
sections together with the two rational old bisections `orbit-1d5f2` and
`orbit-0abc2` have frame-coordinate determinant `+/-1`; they are a saturated
rank-17 basis of the hidden fibration's published-R17 frame.  The direct
artifact has SHA-256
`3f676dd0ce76da7f3092b073519b11af964c6982bfbc7262057d0f5c66234b9f`.

The inherited height-four cheap set has 82 members.  It recovers the requested
seeds `P4` and `P4-P3`, constructs all covers and lifts exactly, and finds 82
distinct extensions with no collision.  The cover and collision artifacts
have SHA-256 hashes
`43b47129d078522c3a68812115247fd0d0c4d67d42af0886b0d70b115f814830`
and
`f84c00fa4fcc9c861c8011ca8e19440eaa82f804cc1dce60cee0caaf8b1dac0f`.
All 3,321 cheap-set products are distinct, with no third-character or formal
catalogue match; that artifact has SHA-256
`52d6f5a346601290bdb307c084e3445a7b40fda00e4892ecd6e12baa08042c2c`.

The full published-R17 frame has 39,120 classes.  Its merged equation artifact
has SHA-256
`faa98745f9bc8fb304493533ca93ecbb4b36606981ed2899282de381305fdeed`.
Exact hashing again returns

```text
bisections=39120  extensions=39120  collisions=0
```

with collision-artifact SHA-256
`03b57f945ef563bff7f5657317e9859c6463d5a583ce32c751145fb3a0634bcc`.
Every hidden-frame branch is again an irreducible quadratic with a unique
polynomial atom, so complete three-character closure is empty.  The closure
artifact has SHA-256
`7a905e23cb441b028795a2dd264d8063cf9f2d192bdcea4301b63d41ac585928`.
Here the old rank-28 and `q_103b2` source curves have target degrees 12 and 0,
respectively, so they fail the hidden-base quadratic-character gate as well.

```bash
sage -python elkies-k3/scripts/compile_r17_norm12_orbit11952_qq.sage \
  --source-label norm12-orbit-103b2 --check
sage -python \
  elkies-k3/scripts/construct_r17_norm12_11952_inherited_bisections.sage \
  --source-label norm12-orbit-103b2
.venv/bin/python \
  elkies-k3/scripts/analyze_r17_norm12_11952_inherited_products.py \
  --source-label norm12-orbit-103b2
.venv/bin/python \
  elkies-k3/scripts/analyze_r17_norm12_complete_character_closure.py \
  --source-label norm12-orbit-103b2
```

## Consequence and next exact frontier

The requested conclusions

\[
\operatorname{rank}E/\mathbf Q(r)\ge 19,
\qquad
\operatorname{rank}E/\mathbf Q(C)\ge 20
\]

are **not proved**.  The exact computations instead show that neither can be
obtained from an equal cover or a three-character relation inside either
complete smooth rational-bisection frame.  No `V4` curve, genus-one Jacobian,
or rank computation is emitted because the required character triple does not
exist.

The next controlled search would have to enlarge the geometry, not repeat the
same hashes: construct degree-two curves of positive arithmetic genus on the
alternate equation, or produce an independently certified same-base twist
character.  Only after such a character exists does a `V4` Jacobian and exact
three-eigenspace height calculation become meaningful.
