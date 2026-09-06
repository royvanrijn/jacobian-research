# Curve 302: primary construction-recovery project

Curve 302 is the main construction-recovery target. The endpoint is an
explicit elliptic surface over `Q(t)`, its **full generic MW basis**, and a
rational parameter with an exact `Q`-isomorphism and section transport to
`E302`. Generic rank 17 is not a requirement. The additive rank-17 candidate
is evidence to test, not a prescribed generic subgroup.

The boxed surface/basis/specialization endpoint is now complete for a
**constructed generic-rank-nine K3 parent**: the full arithmetic MW basis
and `t0=1` specialization to `E302` are explicit and certified. The original
parent remains **UNKNOWN**. This construction works for general nondegenerate
point configurations. Neither that construction nor the additive core
establishes the discoverers' method. A parent of generic rank 17–20 would
be valuable, with the surface type and its possible rank checked first;
the existing K3 has arithmetic MW rank exactly nine. Every elliptic
fibration over `Q` on this same K3 has arithmetic MW rank at most nine.

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

## Next mathematical gates

1. The baseline arithmetic rank gate is closed. Additional sections or
   fibration changes on this K3 cannot produce a higher arithmetic generic
   rank. Recovering a higher-rank parent requires a different underlying
   surface, such as a different pencil or base change, and an exact reason
   its section configuration reflects the complete 31-point group.
2. For construction recovery beyond this baseline, test candidate parent
   configurations with equations,
   rational sections, a generic rank upper bound, saturation, and the exact
   302 specialization map. Use the additive rank-17 core only as a
   retrospective overlap test, allowing different generic ranks.
3. The constructed MW9 family now meets the equation/basis/specialization
   gate for a separate 302-calibrated rank-search lane. Its known exceptional
   points calibrate recovery at 302; fresh-fibre selection and independent
   point proofs must be distinguished from that retrospective calibration.
   The original construction's provenance remains a separate open problem.

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
