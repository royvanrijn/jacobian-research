# Reproducing the results

The Makefile is the public verification interface.  Run commands from the
repository root after creating the Python environment described in the main
[README](README.md).

## Fast structural check

```bash
make check
```

This compiles the active Python code, checks local Markdown links, and audits
the single status ledger.

## Characteristic-two plane Keller counterexample

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_characteristic_two.py
```

The first block independently replays the Huq--Kuruvilla threefold Jacobian,
collision, inverse cubic, projective normalization charts, wild-radicial
boundary, reconstruction, and zero-pole determinant ledger.  The second block
internally replays Mondello's external theorem
[arXiv:2608.02634v1](https://arxiv.org/abs/2608.02634): the source and target
coordinate changes, skew product and preserved plane fibre, plane Jacobian,
three-point collision, hidden cubic, rational recovery identities,
irreducibility coprimality certificate, and separability witness.  The
degree-one-in-the-actual-target-parameter irreducibility argument in the
canonical note is not inferred from a bounded search.  The same argument is
written separately for arbitrary characteristic-two base fields as
`HKM2-ALLFIELDS`; no perfectness is assumed.  The command also checks that
neither displayed integer formula is a characteristic-zero Keller map.  This
internal replay is distinct from the Lean kernel-check and computationally
separate Harmonic Aristotle replay recorded by Mondello; none is independent
human peer review.

## Characteristic-two plane normalization and wild boundary

Requires Singular:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_plane_boundary.py
```

The dependency-light block verifies the hidden-cubic discriminant, the
normalized source presentation, three compatible reconstruction charts, the
global normalized formula for the second source coordinate, the retained and
missing primes over `Q=0`, their three reduced intersections, and the local
equation giving generic different exponent one.  The Singular block certifies
the integral closure, the primitive-order conductor `(P,T)`, the two upstairs
conductor branches, and the exact reconstruction-boundary ideal.  This is a
finite exact calculation, not a bounded search.  The checker assumes a
working Singular installation in addition to the repository Python
environment.  The SymPy identities and Singular certificate do not constitute
a second independent implementation of the normalization algorithm; the
separate audit requested in the canonical note remains open.

## Characteristic-two plane modulo-four lift obstruction

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_plane_w2_obstruction.py
```

This exact symbolic check isolates the `xy` coefficient of the first
Jacobian variation \(\mathcal D_F\) for an arbitrary polynomial correction
modulo two.  It
proves an all-degree obstruction to a constant-Jacobian lift over `Z/4`; it
is not a bounded correction search.  It also verifies the explicit
determinant-one lift after adjoining one identity coordinate and the
geometric-series identity producing a compatible tower over all finite Witt
levels.

## Six-variable quartic HN Waring rigidity

Requires SymPy in the repository Python environment; the loop-closure
command also requires Singular:

```bash
.venv/bin/python scripts/verify_quartic_hn_waring_rigidity.py
.venv/bin/python scripts/verify_quartic_hn_rank9_one_zero.py
.venv/bin/python scripts/verify_quartic_hn_rank9_top_determinant.py
.venv/bin/python scripts/verify_quartic_hn_rank10_parallel_obstructions.py
.venv/bin/python scripts/verify_quartic_hn_rank10_matroid_survivors.py
.venv/bin/python scripts/verify_quartic_hn_rank10_loop_survivors.py
.venv/bin/python scripts/verify_quartic_hn_rank10_loop_closure.py
.venv/bin/python scripts/verify_quartic_hn_rank10_simple_survivors.py
```

These commands exactly replay the finite codimension-two Gale calculation,
the rank-nine one-zero trace coefficients and complementary-minor gate, the
rank-ten parallel-class obstructions, and the rank-ten nonsimple-survivor
audit.  The last check gives literal and loopless characteristic-zero
counterexamples to the proposed cyclic-complement lemma, constructs the
normalized realization ideals for the frozen catalogue slice, and excludes
all 35 characteristic-zero types from the full Gram branch using the first
two HN traces.  The
catalogue completeness is relative to `matroid-database==0.3`.  The
penultimate command replays the frozen complete Gale-loop census: 115
abstract coloured types, 111 characteristic-zero types, 103 universally
closed by splitting, the self-square support obstruction, or exact SymPy
saturation, with eight special types left.  The final Singular-backed command
closes those eight by four exact characteristic-zero saturations, one
disconnected Waring splitting, and three loop--triple Witt obstructions.
The final simple-survivor command freezes the exact 23-extension/five-type
census, constructs a rational realization of every type, and closes all five
rank-six Gram branches using the universal six-point-flat obstruction.  Thus
the rank-ten branch is empty and every essential six-variable quartic HN
counterexample has Waring rank at least eleven.  See
[`QUARTIC_HN_RANK10_MATROID_SURVIVORS.md`](extended-geometry/QUARTIC_HN_RANK10_MATROID_SURVIVORS.md).
Pass `--details` to either survivor verifier to print normalized matrices,
ideal generators, basis-minor data, and rational witnesses as exact JSON.

To regenerate the external catalogue extraction rather than replay its
pinned artifact, run

```bash
.venv/bin/python scripts/enumerate_quartic_hn_rank10_loop_survivors.py \
  --database-root /path/to/matroid_database/_all
.venv/bin/python scripts/enumerate_quartic_hn_rank10_simple_survivors.py \
  --database-root /path/to/matroid_database/_all
```

The source must be the `matroid-database==0.3` rank-four files.  The
loop enumerator uses the files through nine elements; the simple enumerator
uses the complete nine-element file.  The canonical note records the
source-wheel, catalogue-file, and generated-artifact hashes.

## LND-image Mathieu finite-fiber replay

```bash
.venv/bin/python scripts/verify_lnd_radical_slice_fibers.py
```

For the slice LND `D=d/ds`, this constructs a radical complete intersection
of six points in three vertical fibers.  It checks reducedness, verifies on
a generic degree window that primitive membership is exactly the three
vertical interval conditions, and replays a safe seed and multiplier through
exponent twelve.  It also replays two nonreduced length-two residual schemes
for the carrier `q=s^2`, one on and one off the carrier.  The powers are
bounded regression checks; the all-order arbitrary finite-residual and
monic-carrier finite-residual theorems are proved in
[`LND_MATHIEU_SLICE_CONDUCTOR_FRONTIER.md`](extended-geometry/LND_MATHIEU_SLICE_CONDUCTOR_FRONTIER.md).

## LND-image nonmonic degree-drop search

```bash
.venv/bin/python scripts/search_lnd_nonmonic_degree_drop.py
```

This uses exact rational arithmetic for the carrier `p=x*s-1`, powers
`x^c*p^d` with `(c,d)=(0,1),(0,2),(1,1)`, three primary residual schemes
at the degree-drop point, 256 sparse seeds, and six fixed multipliers.  The
primitive-carrier assertions replay the all-order support-weight exclusion
theorem.  All reported finite prefixes are bounded regressions or candidate
generation only.

## LND-image plinth-divisor search

```bash
.venv/bin/python scripts/search_lnd_plinth_ideal_images.py
```

This searches the linear LND `D=x*d/dy+y*d/dz` on `Q[x,y,z]`, whose local
slice `y/x` has plinth element `x`.  In each required homogeneous degree it
constructs `D(I_n)` exactly by rational linear algebra.  Five homogeneous
ideals, 45 sparse seeds, pure powers through six, and four multipliers are
tested.  The output is candidate generation only.

## LND-image reducible-plinth search

Requires Singular:

```bash
.venv/bin/python scripts/search_lnd_reducible_plinth.py
```

This searches
`D_r=x*(x-1)*d/dy+y^r*d/dz` for `r=1,2` on `Q[x,y,z]`.  With weights
`wt(x)=0`, `wt(y)=1`, and `wt(z)=r+1`, image membership in each weight is
an exact finite module calculation over `Q[x]`; Singular computes the
module standard bases.  Five ideals couple the fibers `x=0,1`.  The
`r=1` profile tests 205 seeds through weight three, and the `r=2` profile
tests 210 through weight four; both use pure powers through six and five
multipliers.  A third profile tests
`D=x*(x-1)*d/dy+(y^2+x)*d/dz` on eight ideals and 350 mixed-weight seeds.
It uses bounded normalized primitive lifts and finite quotient/kernel
images, so individual membership decisions remain exact even though the
grading is broken.  Every exponent range is bounded.

## LND-image crossing-plinth search

Requires Singular:

```bash
.venv/bin/python scripts/search_lnd_crossing_plinth.py
```

This searches the four-variable LND
`D=u*v*d/dy+(y^2+u)*d/dz`, whose plinth divisor has intersecting
components `u=0` and `v=0`.  The compiler normalizes primitives modulo
`Q[u,v,3*u*v*z-y^3-3*u*y]`, computes an exact bounded lift by a Singular
module standard basis over `Q[u,v]`, and decides the remaining kernel
correction in each finite quotient exactly.  Five zero-dimensional
crossing ideals and 956 sparse seeds are tested through six pure powers;
mixed powers four through six use the multipliers `1,u,v,y,z`.  Individual
membership decisions are exact, but both exponent windows are bounded.

## LND-image nonprincipal-plinth search

Requires Singular:

```bash
.venv/bin/python scripts/search_lnd_nonprincipal_plinth.py
```

This searches `D=u*d/dx+v*d/dy` on `Q[u,v,x,y]`, with invariant
`w=u*y-v*x` and nonprincipal plinth ideal `(u,v)` in the kernel.
Homogeneous primitive membership is an exact Singular module calculation
over `Q[u,v]`; kernel corrections are decided exactly in five finite
quotients by closure under `u,v,w`.  The search tests 1,055 sparse seeds,
pure powers through six, and mixed powers four through six for
`1,u,v,x,y,w`.  Individual membership decisions are exact; the exponent
windows are bounded.

## LND-image positive-dimensional plinth search

Requires Singular:

```bash
.venv/bin/python scripts/search_lnd_positive_dimensional_plinth.py
```

This retains a free `y`-direction for `D=u*d/dx+v*d/dy`.  The five ideals
begin with `(u,v,x)` and include four nilpotent or tilted plinth jets.
Although their quotients are positive-dimensional, the images of
`u,v,w=u*y-v*x` are nilpotent, so the kernel-image span and every
membership decision are exact without a `y`-degree cutoff.  The checker
also verifies through total degree eight the coefficient-functional
identity used in the free-line corollary of the all-order plinth-power
saturation theorem, and the filtration identity through total degree
seven.  Pure powers through six and mixed powers four through six are
bounded regressions; the written theorem independently proves all five
displayed ideals safe.

## LND-image principal-conductor search

Requires Singular:

```bash
.venv/bin/python scripts/search_lnd_principal_conductor.py
```

For `D=u*d/dx+v*d/dy` and `I=(x)`, the image of the invariant ring modulo
`I` is exactly `Q[u,v,u*y]`.  Consequently a primitive residue monomial
`u^a*v^b*y^c` is correctable exactly when `a>=c`.  The script combines
this exact valuation-face test with Singular primitive lifts for 1,055
seeds, pure powers through six, and mixed powers four through six.  The
membership decisions are exact and untruncated; the exponent windows are
bounded.  The support census also checks 43 `y`-free survivors, thirteen
invariant survivors, four forms touching the slope-zero face, and exact
membership of all 48 survivors in `u*ker(D)[x]`.  It additionally verifies
the exact algebraic square-gate failure: a homogenized shifted-Legendre
face gives `f,f^2 in D((x))` but `f^3 not in D((x))`.  This is not an
LNED counterexample.  It also checks the local-slice identity
`T(D(x*y*a))=(u*y)*a` underlying the exact two-branch criterion, together
with `D(u*x*a)=u*D(x*a)` and `D(u*y*a)=u*D(y*a)` for the aligned and
crossed invariant-content carriers.  For the rational-root carrier
`q_1=u*x+w`, it normalizes `A/(q_1)` by
`x=u*t, y=-(u-v)*t`; the invariant image is
`Q[u,v,u^2*t]`, so the exact residue test is `a>=2*c`.  Seventeen seeds
survive the genuinely eventual window `m=4,5,6`; all seventeen lie in
`u*ker(D)[q_1]`, and no bounded mixed-tail obstruction occurs.  The
checker also verifies the normalization identities through the ladder
`q_n=u^n*x+w`, `1<=n<=4`.  For the first tied carrier
`q=u*v*x+w`, its exact quotient cone is `a>=2*c, b>=c`.  Eight seeds
survive powers `4,5,6`; all eight lie in `u*v*ker(D)[q]`, with no bounded
mixed-tail obstruction.  The checker verifies the two-prime
normalization identities for `1<=r,s<=3`.  Finally, for a sample
invariant-affine coordinate `h=b_0+b_1*x+b_2*y`, it verifies
`D(h)=b_1*u+b_2*v`, `D^2(h)=0`, and both inverse-chart identities
expressing `x,y` in `ker(D)_(D(h))[h]`.

The all-order proof in the canonical note uses the
full eventual-power hypothesis and the one-variable polynomial moment
lemma to prove that `D((x))`, and hence `D((ell(x,y)))` for every nonzero
linear form `ell`, is Mathieu--Zhao.  Divisibility bootstrapping closes
all powers `(ell^d)`, while primitive evaluation at two generic roots
proves zero Mathieu radical for carriers such as `(x*y)` and
`(x*(x-1))`.  A fixed-denominator extension of the local-slice proof also
closes `(a*ell)` for `a` in `{u,v}` and `ell` in `{x,y}`.  A lowest-face
moment argument closes every rational-root ladder carrier
`q_n=u^n*x+w`, `n>=1`, and the paired lowest-face argument closes the
two-prime grid `q_(r,s)=u^r*v^s*x+w`, `r,s>=1`.  None of these all-order
results is inferred from the bounded census.  The same lowest-face proof
closes `u^r*v^s*x+b` for every invariant intercept `b in ker(D)`, with
`r>=1` and `s>=0`; this arbitrary-intercept extension is a written
theorem rather than a separate bounded search.  Prime-by-prime lowest
faces close every coprime `q=a*x+b`, `a,b in ker(D)`.  On the aligned
condition `v*(b mod u)=(a mod u)*w`, one has `q in u*A`, so pure
membership itself forces the missing `u`-content.  A nontrivial common
invariant factor is removed by the Mathieu scaling lemma: if `M` is
Mathieu--Zhao and `c*M` is contained in `M`, then `c*M` is
Mathieu--Zhao.  Consequently every `q=a*x+b` with nonzero
`a in ker(D)` is closed, without a coprimality assumption.  The intrinsic
condition for generic orbit degree one is `D^2(h)=0`, `D(h)!=0`.
Prime-local inverse charts and the same moment-face argument close every
irreducible such `h`; divisibility bootstrapping closes all powers, and
the scaling lemma restores invariant content.  Factoring an arbitrary
principal carrier over the generic orbit now gives either at least two
distinct roots (zero eventual-power radical) or one invariant-affine
irreducible factor with multiplicity.  Thus the canonical note proves
`D(q*A)` Mathieu--Zhao for every nonzero `q in A` for this model
derivation.  The checker only replays the coordinate identities and
bounded searches; the all-order conclusion is deductive.

## Checked-in Lean projects

All four local Lean packages use the pinned Lean/Mathlib `v4.32.1` release.
Build their default targets and audit their source policies with:

```bash
make verify-lean-local
```

This rejects `sorry` and `admit` throughout `formal/`, rejects unexpected
explicit axioms, checks the finite-étale publication-certificate import
boundary, and builds `discriminant-pencils`, `finite-etale-keller`, `gmc2`,
and `gvc`.  The GMC(2) package deliberately exposes exactly two mathematical
inputs as axioms; their names and roles are documented in
[`formal/gmc2/README.md`](formal/gmc2/README.md).  The GVC package contains no
explicit axioms, but is currently a partial audit.  It constructs and proves
both the minimal concrete and arbitrary cusp-profile quadric phase bridges,
giving the full profile failure family over every characteristic-zero field
and an unconditional counterexample in every finite dimension at least
three.  Its sole remaining bridge structure isolates the unformalized binary
envelope obligations: global envelope closure and the `delta = 0` equal-face
ordering from shifted-ray separation.  The complete finite-support
common-threshold cutoff is now proved in Lean: finite maximization constructs
the integral coordinate cut, which is refined to a strict positive weight
whose unit gap is amplified under powers.  Empty equality faces are handled
directly.  The finite core of Lemma 3.1 is also checked: Mathlib's Hall
theorem extracts a deficient set, two-dimensional linear algebra localizes it
in one direction class, and exact counting yields the sharp `d - e + 1`
annihilator bound.  Lean also proves the coordinate-free power divisibilities
for both displayed normal forms.  The translated
Duistermaat--van der Kallen/polarization
step that supplies the absence of a matching remains explicit.  Its
checked core now includes the algebraic beta and full endpoint-profile
coefficient calculations, the literal multivariate profile family's
degree/order formula, the rational `p`-adic factorial-valuation lemma,
the coefficientwise Reynolds/Laurent phase identity and concrete endpoint
extraction, the apolar contraction and operator-composition laws,
characteristic-zero base change, unused-variable padding, and the final
negative-slope envelope crossing.  The full-profile proof enforces the
manuscript's declared-degree condition `S.natDegree <= e`, proves arbitrary
even-phase extraction and the shifted primitive identity, and constructs
the previously missing bridge.  Its exact coverage is
documented in [`formal/gvc/README.md`](formal/gvc/README.md).

Build only the GVC audit with:

```bash
make verify-gvc-lean
```

## Stable core

```bash
make verify-minimal
make verify-core
make verify-foundations
```

`verify-minimal` uses only the Python standard library for the foundational
map.  `verify-core` adds the cubic marked-root and exact-image implementations.
Its normalized-factorization certificate checks both polynomial compositions
across `a=0`, residual-torus equivariance, determinant `-1` for normalized
multiplication, and the two explicit linear changes recovering the announced
map.  The same target runs the scoped ordinary-degree-six boundary audit;
it can be replayed separately with:

```bash
.venv/bin/python scripts/verify_ordinary_degree_six_boundary_audit.py
```

This verifies the exact asymmetric `(1,2)` determinant and degree floor,
both balanced `2+2` resolution charts and their transition in addition to
the cone/Cox identities, the affine-linear Wronskian reduction, the
impossibility of the `(1,2)` cubic profile, the removable-jet residue
obstruction for `(0,3)`, the rational-map-degree obstruction for a
nonconstant boundary jet, and the double-pole residue obstruction for a
nonzero constant boundary jet.  It also checks the pure-`C` weighted `D^3`
divisor ledger, the degree-eleven floor for every `z`-linear standard
reciprocal clearing, and the nodal conductor character rank.  It does not
enumerate arbitrary affine modifications, nonmonomial Wronskian profiles,
or all maps with two boundary relations.

It also checks the normalized `(2,3)` factorization slice: the unimodular
boundary lattice, class `L^5-L^3`, direct counts `q^5-q^3` for four small
prime fields, and generic degree ten.  A separate two-chart certificate
checks the Euclidean quadratic norm, its affine-modification presentation
over `A^2 x SL_2`, the complementary-chart transition, and the integral
residue coefficients `1`, `2`, and the nonzero mod-two boundary used to prove
that integral cohomology is `Z` in degrees zero and three only.  The scripts
check the algebraic, arithmetic, and Gysin inputs; the written audit supplies
the localization sequences and homotopy argument.
The same slice now has an exact invariant-kernel certificate.  It constructs
the primitive saturated LND `D7`, verifies
`ker(D10) intersect ker(D7)=k[K,H,V]`,
`ker(D7)=k[K,H,V,s]`, and `ker(D22)=k[K,H,V,W]`, and checks the boundary
identities used by the minimal-pole proof.  The finite/non-finite control
target also replays Maubach's cusp-base ladder:

```bash
make verify-hilbert14-invariants
```

The bounded ladder replay is regression evidence.  The uniform
modulo-`T^4` degree argument in
[`HILBERT14_INVARIANT_KERNEL_PROGRAM.md`](extended-geometry/HILBERT14_INVARIANT_KERNEL_PROGRAM.md)
is the non-finite-generation proof.  The same target verifies the next
normalized `(2,4)` experiment: three triangular gauge LNDs, the exact
generic quotient `M^2+4*U*N^2=256*a^4`, the regular boundary classes
`C,Q,S`, and the boundary-linear relation that terminates the saturation
ladder.  It also checks the induced third LND and the exact finitely
generated triple intersection `k[a,U,N,M,C]`, whose boundary is a cusp.
The written minimal-pole proofs are in
[`QUADRATIC_QUARTIC_HILBERT14_SLICE.md`](extended-geometry/QUADRATIC_QUARTIC_HILBERT14_SLICE.md).
Finally, the target runs the genuine multiboundary control.  It checks the
two commuting cusp LNDs, the invariant grid
`s^2*t^2*(X+sY)^m*(U+tV)^n`, and the conductor-square replay modulo
`(s^4,t^4)`.  The arbitrary-bidegree rectangle escape in
[`MULTIBOUNDARY_HILBERT14_CONTROL.md`](extended-geometry/MULTIBOUNDARY_HILBERT14_CONTROL.md)
is the non-finite-generation proof.  The same note computes the exact
finite-generation ideal as the conductor
`s^2*t^2*k[s,t,P,Q]` to the normalized-ambient invariant algebra and gives
its four infinite return ladders and infinite monomial SAGBI basis.  The
checker replays the monomial conductor criterion in a configurable box; the
written localization/specialization argument proves the arbitrary-degree
statement.  The same argument is proved for every
`tensor_i(k+t_i^2*k[t_i,P_i])`: its finite-generation ideal is the product
conductor and has `2^r` infinite return ladders.  The checker replays this
formula through `r=4` by default.  The note also proves that the two leading
divisors in a tangent-normalized factorization slice are disjoint and
explains why the coupled three-boundary Cox fill is a different branch.
The same target verifies the general weight-`(1,-1,-k)` invariant-coordinate
Jacobian reduction for `k=1,2,3,4`, including the foundational
`(-2,-1,1)` output weights.  It then reconstructs the complete
sixteen-monomial coefficient ideal, proves the gauge-fixed dual-number
presentation, extracts the infinitesimal deformation and its quadratic
obstruction, and separates it from the affine left--right orbit.  The same
target independently rewrites that normalized ideal as three univariate
weighted-Wronskian layers and checks its exact Poisson-square and tangent-pencil
identities.  The leading layer exposes the quadratic obstruction directly,
and two further unit-ideal checks eliminate both one-sided nonconstant-`C`
boundary charts.

## Global low-degree support census below `(7,6,4)`

Generate the eight support-first stage ledgers and their manifest with:

```bash
.venv/bin/python scripts/compile_global_low_degree_census.py
```

This enumerates the 74 invariant degree flags, the complete raw-degree-seven
exact-support strata through six nonlinear monomial occurrences, every
determinant bucket, and every integer infinity weight modulo exposed Newton
faces and coordinate strata.  It then runs the sign SMT gate, exact Singular
coefficient-torus algebra over `F_11`, `F_13`, `F_17`, and `QQ`, plus an
independent SymPy rational Gröbner replay.  The pinned result has `30`, `85`,
and `1694` determinant-balanced labelled supports in sizes four, five, and
six; their `913` residual-symmetry representatives all have unit exact ideals.
The dense quadratic collision ideal is also `(1)`.

Replay every pinned JSON decision with:

```bash
make verify-global-low-degree-census
```

The result is complete only through nonlinear support six and for the dense
degree-at-most-two row.  It proves a support lower bound of seven below
`(7,6,4)`, without asserting attainment at seven;
it does not claim the cardinality-unbounded census is complete.

`verify-foundations` adds the weighted construction and its clean-room checker.
It also runs the all-degree rational-fiber checker, whose symbolic odd/even
identities prove uniform admissibility and whose exact degrees `3,...,100`
remain as a regression:

```bash
.venv/bin/python scripts/verify_padic_inverse_branches.py
.venv/bin/python scripts/verify_foundational_arithmetic_dynamics.py
.venv/bin/python scripts/verify_composite_degree_twelve.py
.venv/bin/python scripts/verify_degree_twelve_wreath_elimination.py
.venv/bin/python scripts/verify_all_degree_rational_fibers.py
.venv/bin/python scripts/verify_finite_etale_keller_fibers.py
.venv/bin/python scripts/verify_common_arithmetic_fibers.py
.venv/bin/python scripts/verify_locally_prescribed_common_fibers.py
.venv/bin/python scripts/search_cross_family_collision.py
.venv/bin/python scripts/verify_universal_quartic_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_quartic_gauge_multiplicity.py
.venv/bin/python scripts/verify_universal_cubic_gauge_multiplicity.py
.venv/bin/python scripts/verify_universal_power_shifted_gauge_multiplicity.py
.venv/bin/python scripts/verify_whole_plane_stable_multiplicity.py
.venv/bin/python scripts/verify_universal_quintic_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_higher_degree_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_multiplicity_witness_cards.py
.venv/bin/python scripts/verify_universal_relative_keller_map.py
.venv/bin/python scripts/verify_generic_tschirnhaus_non_descent.py
.venv/bin/python scripts/verify_rank_five_tschirnhaus_transition_locus.py
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --module-resolution
.venv/bin/python scripts/verify_keller_tschirnhaus_descent_567.py
.venv/bin/python scripts/verify_rank_three_collision_descent.py
.venv/bin/python scripts/verify_rank_four_collision_cross_ratio.py
.venv/bin/python scripts/verify_low_rank_multiplicity_boundaries.py
.venv/bin/python scripts/verify_real_fiber_spectrum.py
.venv/bin/python scripts/verify_adelic_fiber_engineering.py
.venv/bin/python scripts/verify_local_global_keller_fibers.py
.venv/bin/python scripts/verify_a5_grunwald_keller_fiber.py
.venv/bin/python scripts/verify_hasse_keller_fiber.py
.venv/bin/python scripts/verify_infinite_hasse_keller_fibers.py
.venv/bin/python scripts/verify_multiplicative_hasse_artifact.py
```

The finite-étale Keller-fiber checker includes the exact ordered collision
fiber. For the degree-three, degree-four, and degree-five presentations it
verifies the diagonal/off-diagonal Chinese-remainder decomposition, its
explicit separability idempotent, and ranks `N^2`, `N`, and `N*(N-1)`. It
also verifies three cubic `S_3` normal-closure sheets and the optimal Hasse
fiber decomposition `A5 tensor A5 = A5 times (N6^3 times L2)`. The
presentation-independent collision algebra, diagonal kernel, ordered-pair
universal property, and obstruction rank are checked in Lean.

The universal-relative checker applies this interface to the Osada
`T^N-T-1` root covers in ranks three through eight.  It verifies the
divided-difference idempotent and the exact collision, diagonal, and
off-diagonal standard-monomial ranks `N^2`, `N`, and `N*(N-1)`.  It also
enumerates every `S_N` orbit of ordered distinct `m`-tuples for
`1<=m<=N<=8`, checking the rank `N!/(N-m)!` and stabilizer `(N-m)!`.
The rank-three descent checker then identifies the cubic ordered-pair sheet
with the full `S_3` frame torsor and verifies the exact projective
interpolation cocycle, quadratic Tschirnhaus boundary ledger,
target-localized factorization transport, saturated global stabilizer, and
fixed-map scaling equivariance.  It does not classify nonlinear polynomial
self-equivalences outside the canonical factorization transport.
The rank-four continuation checks that ordered triples, rather than ordered
pairs, give the full `S_4` frame.  It factors the fourth-root projective
interpolation residual and labeled cross-ratio difference by the same exact
defect, separates that defect from the primitive-element boundary, and
clears it into the universal quartic Keller target coordinates.  It does not
assume that every Keller-incidence equivalence is projective on the root
line.

The all-rank continuation is:

```bash
.venv/bin/python scripts/verify_all_rank_collision_projective_descent.py
```

It completes `Conf_(N-1)` to the full `S_N` frame, verifies the intrinsic
rank-at-most-three criterion for the columns `1,r,u,r*u`, constructs its
normalized polynomial coefficient matrix, recovers the automatic cubic and
quartic cross-ratio cases, and checks the `N-3` independent framed residuals.
It also supplies exact projective and primitive-nonprojective witnesses in
every tested rank.  The bounded replay supports the written all-rank
linear-algebra proof; it does not claim that every Keller equivalence acts
projectively on the root line.

The generic stable non-descent continuation is:

```bash
.venv/bin/python scripts/verify_generic_tschirnhaus_non_descent.py
```

It verifies that the split change `r -> r+r^2` has a nonzero projective
minor in every rank, checks the quintic `I_5` base case, and proves
symbolically for every `N>=6` that

```text
J_N(P_(r+r^2))-J_N(P_r)
=-(N-1)(7N+11)/(30N(N+1)(N+2)).
```

Direct root-polynomial calculations through rank twenty are regressions.
Dominance of the presentation-to-boundary map and the resulting generic
codimension statement are written geometric proofs, not bounded searches.

The rank-five transition-locus continuation is:

```bash
.venv/bin/python scripts/verify_rank_five_tschirnhaus_transition_locus.py
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py
```

It computes the ambient stable-equivalence hypersurface and the two
labelled projective residuals, verifies local dimensions `4`, `3`, and `2`
for the ambient, projective, and intersection loci, and checks the explicit
coefficient-torus equivalence.  It also proves that this canonical
equivalence carries the selected complete fibre exactly on the root-scaling
locus.  The canonical note then determines the fixed map's standard marked
stable target orbit completely.  It does not classify vertical
automorphisms of the added identity factors.

The second command supplies the exact fixed-map calculations.  It
factors the prime quintic ramified discriminant, computes the exact
logarithmic-vector-field spaces through quotient degree twelve using FLINT
integer nullspaces, and applies a three-point characteristic-zero Jacobian
Groebner test.  An exact triple-root point forces the boundary multiplier to
be one in every degree.  Exact recursive Newton-face pruning proves that the
stable marked-target orbit is a point through total target degree twenty-eight;
every unstabilized target self-equivalence in every degree is the identity.
The all-degree unstabilized conclusion uses the coordinate-polynomial
intruder theorem at `P^2*B^5*C`.  Kuroda's stable-invariant theorem, applied
to conjugates of every stable translation by the target automorphism and its
inverse, makes the standard marked orbit a point for arbitrary stabilization.
All bounded stable
branches expose one of `P^12*C^4` and `P^2*B^5*C`.
The checker also proves in all degrees that these are the only positive
upper Newton vertices.  They tie on
`10*w_P-5*w_B+3*w_C=0`, and it verifies an explicit logarithmic field whose
two leading contributions cancel there.  More sharply, its `P`-zero Koszul
ladder first ties at target degree fifty, where the UFD cube condition
fails, and first admits leading cancellation at degree fifty-five.  Thus
unrestricted monomial avoidance is false; exact boundary preservation on
that one binomial wall remains open.  This is a comparatively expensive
exact regression.  The optional third command requires Macaulay2.  It
proves that the homogenized logarithmic module has two generators in
quotient degree seven, thirteen in degree eight, eighteen first relations
in degree nine, and six second relations in degree ten.  Equivalently, its
filtered Hilbert numerator is `2*t^7+13*t^8-18*t^9+6*t^10`.

The optional Newton-topology certificate requires Singular:

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --generic-fibre-newton
```

It checks all forty-six nontrivial Newton faces of `H-h` and its coordinate
restrictions, obtains normalized-volume contributions
`8,2,0;-38,-52,-2;328`, and certifies `chi(H=h)=246`.  Thus the
vanishing-`H^2` stable-rigidity shortcut does not apply; this calculation is
independent of the Kuroda descent proof.

The optional all-degree wall research calculations are:

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-koszul-homology
make verify-rank-five-singular-support
```

The singular-support target is the short form of these three commands:

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-triple-root-prime
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-two-double-root-prime
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-singular-boundary
```

These use exact characteristic zero and require Macaulay2.  Add
`--research-characteristic=1000003` for the much faster good-prime
discovery pass.  The homology command proves that the non-Koszul quotient
has dimension two, degree 296, and the Betti table recorded in the canonical
note.  The next two commands construct the prime triple-root and
two-double-root curves by contraction; their projective degrees are
seventeen and nineteen.  The final command proves that the affine
`P=0` chart is empty and that the radical at infinity is `(Z,P*C)`.
Together with the root-partition argument in the canonical note, these
three targeted commands give the four minimal supports without asking
Macaulay2 for a blind primary decomposition.  The older
`--research-singular-primes` mode remains available as an expensive
independent comparison, but is not part of the proof chain.

The first exact continuation of the cancellable target-degree-55 Koszul
wall is:

```bash
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-koszul-hensel --research-depth=4
.venv/bin/python scripts/verify_rank_five_stable_target_stabilizer.py \
  --research-koszul-hensel --research-depth=9 \
  --research-zero-l=l_3_0 --research-continue-constraints
```

The first command projects to `B`-degree zero and forces `l_3_0^4=0`.
After imposing `l_3_0=0`, the second uniquely lifts all divisible residuals
with lower homogeneous pieces of `G`, records the depth-eight plane
constraint, and proves that it and the five depth-nine equations generate
the unit ideal.  This excludes the normalized two-generator family
`V_B=L*H_C`, `V_C=-L*H_B-eta*G*H`.  It does not include the independent
`H`-multiple in `V_B` or any non-Koszul singular-support class.

The first nonprojective rank-four continuation is:

```bash
.venv/bin/python scripts/verify_rank_four_nonprojective_keller_lift.py
```

It verifies the exact fifth-power ground-field orbit class, shows why the
`(1,2,3,4)` quadratic witness carries a separate rational Kummer twist, and
constructs an arithmetic-neutral witness with seed ratio `4^5`.  It checks
the explicit source--target scaling to the single fixed map
`F_(-124416)`, the two endpoint fibers, and their residual target
translation.  The straight target line has inverse polynomial
`-(S-12)(S+12)(S^2+24S+108*lambda)/3456`; its discriminant and wrong framed
sheet partition are checked exactly.  A label-preserving rational path has
a divergence-free polynomial first-order source lift of degrees
`(55,53,55)`, while the fixed-map target translation has first-order degrees
`(31,29,31)`.  All-finite-order liftability is supplied by the separate
formal-orbit theorem.  At `lambda=-4`, the checker reconstructs the exact
two-point affine fiber and uses integer-orbit fiber invariance to rule out a
polynomial lift of the straight target translation.  It does not claim a
high-degree endpoint self-equivalence: the prime ramified discriminant has
ordinary degree thirteen, so divisibility and the exact `mu_5` orbit test
exclude every target candidate through degree twelve.

The next exact frontier requires Singular:

```bash
.venv/bin/python scripts/verify_rank_four_degree_eighteen_target_obstruction.py
```

It computes the logarithmic-derivation nullities
`(0,0,0,0,1,7)` through multiplier degree five.  For target degree at most
eighteen, the endpoint condition leaves four parameters.  Ten exact
constant-Jacobian evaluations in those parameters generate the unit ideal
over `QQ`; Singular returns the reduced basis `[1]`.  This excludes all
endpoint target symmetries through degree eighteen and leaves degree
nineteen as the first unresolved case.  It is an exact Gröbner
inconsistency, not a bounded coefficient search.

Build the Lean interface with:

```bash
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.CollisionFiber
lake build FiniteEtaleKeller.PaperCertificate
```

The relative whole-plane statement and the exact stable-separation
certificates are checked independently in Lean:

```bash
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.WholePlaneStableMultiplicity
lake build FiniteEtaleKeller.StableSeparationCertificates
```

The first build treats all `(B,C)` simultaneously and restricts naturally to
the universal inverse-discriminant open.  The second proves the Fitting
shoelace formula, Laurent translation and unimodular invariance, strict
power-shift separation, and cubic boundary-count separation.  The geometric
normalization/Fitting and boundary-exhaustion inputs remain explicit
interfaces.

The first reconstruction step beyond the common `P=1` plane is checked by:

```bash
.venv/bin/python scripts/verify_two_marked_fiber_gauge_reconstruction.py
```

For a root-marked fiber over `P=p`, the checker extracts the linear
coefficient `g_1/(g_N*p^(N+m))` of the monic root annihilator.  It verifies
that the planes `P=1,2` recover the normalized seed, the gauge exponent, and
therefore the stable Fitting area.  It also checks the sharp periodic
counterexamples on finite collections of torsion planes and recovers `m`
from the pole order `N+m` on a transverse affine line.  The result retains
the inverse-root generator and the base character `P`; it is not an
unmarked-cover Torelli theorem.

The unrestricted finite-sampling question is answered negatively by:

```bash
.venv/bin/python scripts/verify_finite_marked_plane_nonreconstruction.py
```

The checker replaces the common monomial power shift by an arbitrary
polynomial multiplier `R(P)`.  Exact interpolation makes `R=1` on every
prescribed sample plane, while simple roots of `R` each create a Newton
block `(0,0)->(3,0)->(N,1)` and one boundary prime with
`(e,f)=(N-3,1)`.  Squarefree interpolants of increasing degree therefore
give maps agreeing on any finite collection of complete marked inverse
planes but having strictly increasing stable boundary counts
`deg(R)+2`.  The direct quartic calculation independently checks the full
determinant, inverse, and reconstruction identities.

The final full-boundary reconstruction layer is checked by:

```bash
.venv/bin/python scripts/verify_polynomial_gauge_decorated_torelli.py
```

On the clean polynomial-multiplier locus, the boundary ledger selects
`P=0` intrinsically and hence recovers `P` up to scalar.  The checker then
verifies that the unmarked ramified-stratum Fitting divisor excludes root
inversion, kills every punctured-base unit twist, reconstructs the seed and
`R(P)` coefficientwise, and yields exactly the ordinary source--target
scaling action.  The written theorem also records the general converse:
the full finite-normalization morphism carrying its reconstruction boundary
restricts to a left--right equivalence, so that top layer is a complete
stable invariant without a chosen inverse root.

To intentionally refresh the pinned count artifact after changing its
generator, run:

```bash
.venv/bin/python scripts/count_multiplicative_hasse_parameters.py \
  --bound 1000000 \
  --output artifacts/generated-results/multiplicative_hasse_parameters_1000000.json
```

Then record the changed file hash in the canonical note and the Hasse paper
before committing. Continue the broader verifier catalogue with:

```bash
.venv/bin/python scripts/verify_fixed_quintic_arithmetic_zoo.py
.venv/bin/python scripts/verify_stratified_adelic_engineering.py
```

The first command generates the explicit decomposable degree-twelve map
`F_4 o F_3`, checks both determinant-one factors, and records the expanded
coordinate fingerprints and the `4*3` intermediate-field tower.  Pass
`--print-map` to print all three expanded coordinates.  The second command
reduces the pulled-back cubic discriminant modulo the quartic inverse
equation, factors the saturated resultant as `C^8 Q` with `Q` irreducible
of exponent one, separates the other boundary image, and certifies
`Mon(F_4 o F_3)=S_3 wr S_4`.

The cubic and power-shifted multiplicity checkers also exercise the public
`compile_polynomial_to_keller_fiber(..., stable_parameter=k)` path.  They
compare its maps with the symbolic constructions, preserve the selected
inverse polynomial, and audit the returned boundary-count or Newton-area
record.  Stable functoriality is supplied by the corresponding written
theorems rather than inferred from these finite regressions.

The local-to-global checker audits the ramified quintic coefficient CRT:
the prescribed algebras at `2` and `3`, signature `(1,2)`, cycle types `(5)`
at `5` and `(2,2,1)` at `7`, and the determinant-one quadratic-gauge
compilation with complete target `(1,0,-98/809)`.  It first reconstructs the
polynomial through the generic prime-power coefficient synthesizer in
`jcsearch.local_global`, including its common denominator `1261`, and then
feeds it through the shared end-to-end compiler in `jcsearch.keller_fiber`.
It also derives the universal radii `2^5` and `3^3` from the two local
discriminants and reconstructs the fully automatic witness with common
denominator `30241`.

The fixed-map Hasse checker verifies the determinant, target-line
factorization, modulo-`9` Hensel reduction, and the first prime and composite
parameters.  The multiplicative enumerator then lists both the full sufficient
family and the clean prime-support subfamily through `a=10^6`, checks
primitive coordinates and height `32*a`, counts one- and two-prime members,
and records stable SHA-256 digests in the generated JSON certificate.

The locally prescribed common-fiber checker keeps both maps fixed.  It
derives parameter radii `2^9`, `3^3`, and `5`, constructs
`u=95231/69121`, verifies the ramified completions at `2` and `3`, proves
inertness at `5` and signature `(2,2)`, and checks both transported common
targets.

The common-fiber checker synthesizes the arithmetic transfer and stable
boundary results.  It verifies the fixed all-degree pair over `Q`, the
fixed quartic triple over `Q(sqrt(-2))`, the small rational quartic, and the
mod-`17` irreducibility certificates for the connected triple fibers.  The
following search command enumerates the declared rational tangent-chord,
scale, and constant-term boxes, checks the weighted and quadratic
presentation gates, filters at a split prime of `Q(sqrt(-2))`, and recovers
the coefficient-minimal shared polynomial
`9W^4-19W^3+10W^2-8W-4`.

The universal-quartic checker verifies the trace-zero quartic
tangent-chord factorization, its diagonal rank-five trace quadric, the three
possible indefinite real signatures, the normalized weighted parameter
`alpha=u/e-1/2`, and the finite clean-locus exclusions.  Over a number
field, the uniform existence and infinitude proof uses local isotropy,
Hasse--Minkowski, and rational-point density on the resulting smooth quadric,
as recorded in
[`verified/UNIVERSAL_QUARTIC_FIBER_MULTIPLICITY.md`](verified/UNIVERSAL_QUARTIC_FIBER_MULTIPLICITY.md);
it is not a bounded-search conclusion.

The universal quartic gauge checker verifies the second, unconditional
rank-four mechanism.  It raises the quartic lift from `P^4*S^4` to
`P^(m+4)*S^4`, checks the denominator-free determinant and inverse identities,
keeps the selected inverse polynomial fixed at `P=1`, and computes the
ramified-stratum Fitting-support indices `2*m+5`.  The written stable
normalization argument proves infinite multiplicity for every rank-four
finite etale algebra over every characteristic-zero field, including the
anisotropic trace-chord example.

The universal cubic gauge checker verifies the final low-rank step.  It adds
`g_3(P^n-P^3)S^3` to the cubic lift, checks the paired polynomial corrections
and determinant, and confirms that the selected inverse cubic at `P=1` is
unchanged.  The degree-drop polynomial `1+P^(n-1)-P^2` has `n-1` simple
nonzero geometric roots; the written Newton and boundary-exhaustion argument
turns them into `n-1` intrinsic unramified boundary target components.  Their
count separates the maps stably and proves universal infinite multiplicity
in rank three.

The all-degree power-shift checker verifies the uniform replacement
`g_j*P^j*S^j -> g_j*P^(j+m)*S^j` for every `j>=4`.  Representative
three-variable expansions check the determinant, inverse, and reconstruction
identities.  Exact convex-hull calculations verify that the normalized
ramified Fitting Newton polygon has area
`2*N-3+(N-2)*m`; the written stable-normalization argument makes this a
strict stable invariant for every `N>=4`.  This unifies the quartic and
higher-rank multiplicity mechanisms without trace or translation input.

The universal-quintic checker verifies the translated quintic derivative
jets, the primitive relation `(-1,-6,5)` among the three quadratic-gauge
stable-moduli weights, the invariant
`a_5^5/(a_3*a_4^6)=g_5^5*g_1^2/(g_3*g_4^6)`, and its forced pole after
choosing a trace-zero primitive generator with nonzero second trace moment.
The written argument then gives infinitely many stable classes for every
rank-five finite etale algebra over every characteristic-zero field.  The
higher-degree checker verifies the universal top-weight relation and

```text
J_N=a_(N-2)*a_N/a_(N-1)^2
   =(N-1)/(2N)+c_(N-2)/(N^2*s^2),
```

which proves the same conclusion in every rank `N>=6`.  Together with the
power-shifted quartic argument, this gives infinite universal multiplicity
over every characteristic-zero field in every rank at least three, as recorded in
[`verified/UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md`](verified/UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md).
The next checker supplies connected degree-four, degree-five, and degree-six
three-map witness cards, including modular irreducibility, exact targets,
complete inverse identities, and distinct stable invariant values.  The
universal-relative checker verifies the block-triangular determinant through
the compact reciprocal chart, the unchanged-coordinate promotion to one
absolute map of `A^N`, the exact `U_N=V_N x A1` normalization, and the sharp
`N-3`-parameter inverse-polynomial specialization through degree twelve.  The
written theorem proves all-rank finite-etale universality and imports `S_N`
monodromy and primitive-monodromy atomicity; those last two steps are not
inferred from a bounded symbolic computation.  It also distinguishes
presentation dominance from stack descent and essential dimension.  The
adversarial extension gives a closed-form all-rank target for the Osada
`S_N` family `T^N-T-1`, checks it through rank twelve, pins additional
connected, split, and disconnected targets in ranks three through six, and
verifies the genuine degree-drop, bad-translation, and repeated-root
boundaries.

The formalized promoted map, coefficient compiler, and witness cards are:

```bash
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.UniversalPromotedBlock
lake build FiniteEtaleKeller.UniversalParameterCompiler
lake build FiniteEtaleKeller.UniversalParameterQuotient
lake build FiniteEtaleKeller.UniversalParameterWitnesses
lake build FiniteEtaleKeller.UniversalPromotedMap
lake build FiniteEtaleKeller.UniversalPromotedGauge
```

These modules prove the abstract unchanged-coordinate block determinant, the
promoted inverse-polynomial identity, its selected degree, nonvanishing of
compiled top parameters, automatic admissible translation, invariance of the
quotient algebra under the compiler's nonzero normalization and translation,
three exact quartic targets, and the literal promoted map on an `N`-element
coordinate type with its actual full Jacobian block and determinant-one
identity.  They do not formalize the literal promoted full-fiber/compiler
bridge, its geometric degree, `S_N` monodromy, or stable atomicity.  The
low-rank checker verifies the collapse of all three present cubic mechanisms
and the exact biquadratic trace form used in the written two-step Springer
anisotropy proof.

The Hasse-fiber command expands an explicit degree-eight weighted map, checks
its determinant `-38`, proves that its complete target fiber has no rational
point, and audits roots over `R` and every `Q_p` through the elementary
quadratic-residue covering and the two exceptional Hensel lifts.

The normal-covering front end and its first two exact certificates are
replayed by

```bash
.venv/bin/python scripts/verify_normal_covering_certificates.py
python3 scripts/verify_banks_degree_5_10_candidates.py
.venv/bin/python scripts/verify_degree_six_normal_cover_keller.py
```

The first command independently enumerates the groups and all subgroups in
the `S_3` quintic and `C_2^2` sextic actions, proves conjugate coverage and
trivial common core, computes `gamma(S_3)=2` and `gamma(C_2^2)=3`, and checks
the exact ramified-prime Hensel witnesses.  The second validates the pinned
necessary-candidate transcription of Banks' Table C.1; it does not assert
that every candidate row is arithmetically realized.  The third compiles
`(T^2-2)(T^2-17)(T^2-34)` with the shared quadratic-gauge compiler, verifies
the target `(1,0,528/577)`, and expands its determinant-one Keller map.
For larger finite groups, the GAP front end is loaded with
`Read("scripts/normal_covering_certificate.g");`; the checked-in small
certificates deliberately have a dependency-free Python replay.

The fixed-quintic commands check one determinant-`-2` map and its finite
certificate ledger: all three real quintic signatures, all five transitive
groups `C_5`, `D_5`, `F_{20}`, `A_5`, and `S_5`, split and
quadratic-times-cubic fibers, and all seven unramified partitions modulo
`7`.  The group certificates use witness-prime factor patterns, an explicit
order-five automorphism, a pair-sum resolvent, and Cayley's sextic resolvent.
The clean Hasse row has normalized polynomial
`(T^2-8T+47)(T^3+8T^2+12T+8)` and common quadratic resolvent
`Q(sqrt(-31))`; only `2` and `31` need special local witnesses.  The
original `Q(sqrt(-3))` row remains an independent regression and supplies
the `Q_5` trace obstruction to the standard pure-cubic infinitude route.
Infinitude inside this particular split-seed pencil remains open.

The local `Q_2` action-certificate branch has a separate three-layer replay:

```bash
python3 scripts/verify_gq2_permutation_action.py \
  arithmetic/certificates/gq2_s3_x3_minus_2.json --json
python3 scripts/verify_gq2_permutation_action.py \
  arithmetic/certificates/gq2_s4_mixed_action.json --json
python3 scripts/verify_gq2_permutation_action.py \
  arithmetic/certificates/gq2_common_quintic_stable_pair.json --json
.venv/bin/python scripts/verify_gq2_action_first_keller.py
.venv/bin/python scripts/verify_gq2_s4_quartic_keller.py
.venv/bin/python scripts/verify_marked_q2_stable_separation.py
gp -q scripts/verify_gq2_s4_local_models.gp
gp -q scripts/verify_gq2_local_decompositions.gp
```

The first command is dependency-free and evaluates the exact Roe--Turturean
word ledger, including the finite `omega_2` powers and the normal 2-core
condition.  Its named comparison proves that the marked `S_3` action is the
splitting action of the tame Eisenstein cubic `T^3-2` over `Q_2`.  The second
command translates that polynomial into the degree-three quadratic-gauge
formula, verifies determinant one after output scaling, and checks the complete
inverse polynomial `(S+1)^3-2`.  The mixed `S_4` checker enumerates the three
marked `x_0` orbits over the fixed tame frame, evaluates the candidate
quadratic obstruction as `(1,0,0)`, and compiles
`T^4+4T^2-4T+2` into a determinant-one complete quartic fiber.  The first
PARI/GP command proves that the three classified local quartics have closure
group `S_4`, inertia `A_4`, wild inertia `V_4`, the displayed exact
ramification groups, and normalized relative Stiefel--Whitney bits `(0,0,1)`.
It also verifies the three resolvent Kummer square classes and their sole
product relation over the tame `S_3` closure.  The unique nonzero obstruction
therefore matches the worked `x_0=(12)(34)` orbit to
`T^4+4T^2-4T+2`; the other two orbits remain unordered.
The common-quintic checker verifies the exact unramified marking
`sigma=(1234)(5), tau=x_0=x_1=1`, global irreducibility modulo `17`, both
determinant-one inverse equations, and the stable unit-rank separation
`1 != 2`.  The final PARI/GP
command recomputes exact
ramification-index/residue-degree decompositions at `2` for the selected
quartics, all ten fixed-quintic zoo rows, and the separate quintic witness
card.  It uses maximal-order prime-ideal decomposition rather than bounded
`2`-adic precision.  The checked-in table was last recomputed with PARI/GP
2.17.4 on arm64 Darwin (GMP 6.3.0).  The combined target is:

```bash
make verify-gq2-local-fibers
```

The height-`21` five-row witness card and its separate bounded discovery
audit are reproduced by

```bash
.venv/bin/python scripts/verify_universal_quintic_calculator.py
.venv/bin/python scripts/search_universal_quintic_calculator.py --bound 21
```

The first command uses only exact rational arithmetic and finite-field
factorization.  The second requires PARI/GP, enumerates primitive
projective targets through height `21` modulo the sign involution, and uses
`polgalois` only after exact discriminant and Frobenius-pattern prefilters.
It is bounded computational minimality evidence, separate from the
oracle-free certificates for the five displayed rows.

The mechanically generated finite ledger is checked by

```bash
.venv/bin/python scripts/verify_fixed_quintic_certificate_ledger.py
```

It recomputes all ten rows, their real-root counts and witness-prime
patterns, the seven modulo-`7` partitions, and the `-48*Pi^8` coefficient
Jacobian.  It also runs the three canonical exact checkers and compares both
the Markdown table and
`artifacts/generated-results/fixed_quintic_certificate_ledger.json` with
the generated data.  Pass `--write` only when intentionally refreshing both
generated forms.

Its bounded height search requires PARI/GP:

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_targets.py
.venv/bin/python scripts/search_fixed_quintic_hasse_curves.py
.venv/bin/python scripts/search_fixed_quintic_hasse_discriminants.py
.venv/bin/python scripts/verify_fixed_quintic_hasse_minus_seven.py
.venv/bin/python scripts/verify_fixed_quintic_hasse_minus_thirty_one.py
.venv/bin/python scripts/analyze_fixed_quintic_hasse_minus_thirty_one.py
.venv/bin/python scripts/search_fixed_quintic_hasse_rational_curves.py
.venv/bin/python scripts/search_fixed_quintic_hasse_elliptic_slice.py
.venv/bin/python scripts/search_fixed_quintic_hasse_rank_one_slice.py
```

The first command's default box is stated in its output.  It reports the two
sign-related Hasse targets of projective height `257280` and no other target
below the previous height `458080` in that box.  This is search evidence,
not a global height-minimality claim.  The second command verifies an exact
rational parametrization of the common-quadratic-resolvent incidence, checks
rank-two and rank-one elliptic slices, and searches a bounded proportional
family for irreducible candidates having cubic roots over `Q_2`, `Q_3`, and
`Q_5`.  Its only small-prime survivors are four presentations of the known
Hasse target; it does not test every completion and is not an infinitude
proof.
The rank-two elliptic-slice command closes the
`kappa/A=-1, R=1` route exactly:
the Mordell equation forces `v_2(Pi)=-(2m+1)` and
`v_2(A)=-(3m+1)`, while the cubic factor has a single Newton-polygon
slope `3m+7/3`.  Neither it nor the discriminant-`-3` quadratic has a
`Q_2` root.  Its default 624-point Mordell--Weil enumeration is only a
regression for that all-points proof.
The final command closes the rank-one `kappa/A=5/4, R=4` route at `Q_5`:
the elliptic equation forces `v_5(Pi)=-2m` and `v_5(A)=-3m`; after
translating the cubic by `T=2A+Y`, its unique Newton slope is
`3m-1/3`.  The quadratic has discriminant `-48`, so neither factor has a
`Q_5` root.  Its 24-point multiple regression is again secondary to the
uniform valuation proof.
The third command varies squarefree shared quadratic resolvents.  In its
default integral box it finds a new `Q(sqrt(-7))` Hasse target
`(-7,387/14,400/2401)` of projective height `132741`.  The fourth command
independently checks its factorization, irreducibility, common resolvent, and
exact local witnesses at `2`, `5`, `7`, and `79`.
The wider command

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_discriminants.py \
  --d-bound 60 --r-bound 20 --a-bound 50 --pi-bound 40
```

also finds a `Q(sqrt(-31))` target
`(5,-144/5,-188/3125)` of projective height `90000`.  The fifth command
above audits it independently; only the exceptional primes `2` and `31`
require local witnesses.
The sixth command verifies two exact continuation reductions for the
`Q(sqrt(-31))` row: a genus-two cube curve on the fixed normalized-factor
slice, and a rational trace quadric with a quartic cube condition for affine
variation of the two field generators.  It also enumerates rational
coordinates of height at most `600` on the genus-two slice; only the known
coordinate `Pi=5` occurs.  This last statement is bounded search evidence.
The seventh command proves an exact obstruction to every base line through
`(A,R,Pi)=(-8,2,5)` on the fixed-`-31` common-resolvent double cover: after
recursive square-root reconstruction, the residual ideals have Groebner
basis `[1]` on all three projective direction charts.  It also excludes
every degree-at-most-two curve on each coordinate-fixed slice
`A=-8`, `R=2`, and `Pi=5`: all twelve weighted-projective charts have empty
fiber modulo the good prime `32003`, hence empty characteristic-zero
generic fiber by properness.  Finally, it tests 15024 genuine general
quadratic parametrizations with six integral coefficients in `[-2,2]` and
finds no square pullback.  The line and coordinate-slice results are exact;
the general quadratic result is only bounded evidence.  New primes ramifying
in cubic specializations remain an additional all-prime local obstacle.

The larger fixed-discriminant integral search is

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_discriminants.py \
  --d-value -31 --r-bound 100 --a-bound 200 --pi-bound 200
```

It finds only the certified point and its sign mate.  The low-denominator
rational search is

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_discriminants.py \
  --d-value -31 \
  --r-bound 20 --r-denominator 4 \
  --a-bound 30 --a-denominator 4 \
  --pi-bound 30 --pi-denominator 4 \
  --show-failures
```

It finds two further common-resolvent presentations, both failing at `17`.
Both commands are bounded search evidence, not finiteness theorems.

The last command audits the constructive CRT/weak-approximation lift and an
explicit nonsurjective type-`(3,2)` quintic seed with trivial Hessian
stabilizer and complete fibers of all three quintic signatures, each with
cycle types `(5)` at `7` and `(2,2,1)` at `11`. The preceding adelic command
audits an explicit totally imaginary quartic complete fiber that is inert at
`7` and has unramified splitting type `(2,1,1)` at `11`.

The linear-torus-free quadratic-gauge specialization has a separate exact
certificate:

```bash
make verify-linear-torus-free
```

It checks determinant one, a four-point rational collision, all 734
coefficient equations in `B F = JF A x`, and a displayed `18 x 18` primitive
integer minor of determinant `-5`.  Thus every infinitesimal linear
source-target symmetry vanishes; conjugation makes the result invariant
under independent linear coordinate changes.  A dependency-free clean-room
replay rebuilds the sparse rational polynomial calculation and verifies the
matrix by Bareiss elimination.  A separate parameter calculation proves
within the SymPy checker that the same rational minor has determinant
`(10935/4) g_4^6/g_1^6`, so every admissible quartic quadratic-gauge map has
the same linear-symmetry exclusion.  For the displayed small-coefficient
map, both implementations also verify that the complete `785 x 24` system
allowing constant terms in both vector fields has full column rank.  Hence
the example remains free of affine-linear torus equivariance after
independent affine coordinate changes.

The intrinsic algebraic-torus strengthening is checked by:

```bash
make verify-algebraic-torus-free
```

On the canonical ramified normalization stratum it reconstructs
\(J(P,r)=-1-3Pr^2+4P^4r^3\), computes the scheme-theoretic stabilizer
\(\beta=\alpha^{-2}\), \(\alpha^5=1\), and verifies that its tangent matrix
has determinant \(5\).  The Newton-support pass checks all six permutations:
only the identity and the involution
\(\left(\begin{smallmatrix}-2&-1\\3&2\end{smallmatrix}\right)\) are integral
unimodular, and the ordered second-boundary image rejects the involution
because it does not preserve the intrinsic base character \(P\).  The
checker also verifies the resulting explicit
\(\mu _5\) source--target symmetries of the displayed map.  The canonical
boundary argument then upgrades the conclusion from affine-linear actions:
a connected torus acts trivially on the decorated stratum, hence fixes the
prime nonnormal discriminant hypersurface pointwise; the weight-space lemma
in the canonical note forces the target action to be trivial, and
\(S_4\) deck rigidity forces the source action to be trivial.  Thus no
polynomial left--right representative is algebraic-torus-equivariant.  This
does not classify all discrete or unipotent polynomial self-equivalences,
and literal symmetry-freeness after identity stabilization is impossible:
the added identity coordinates carry tautological torus actions.  The exact
stable conclusion is that every connected action on the pulled-back
decoration is vertical over its intrinsic two-torus; no splitting of such a
vertical action is claimed.

The same checker expands the rational-root sparse representative from
`G(S)=S(S-1)(S-2)(3S+2)`, verifies determinant one and its displayed
four-point rational fiber, and obtains component support counts `(7,51,38)`
and ordinary degrees `(7,26,24)`.  A symbolic coefficient audit shows that
exactly seven generic support coefficients contain `g_2/g_1`, while every
other coefficient is a nonzero Laurent monomial in the admissible
`g_3/g_1,g_4/g_1`.  Hence `g_2=0` is support-minimal in this fixed normal
form.  No absolute sparsity claim under arbitrary polynomial left--right
changes is made.

The bounded exact polynomial left--right sparsity search is:

```bash
.venv/bin/python scripts/search_quartic_lr_sparsity.py \
  --source-degree 2 \
  --target-degree 2 \
  --tadic-max-exponent 12 \
  --scaling-bound 16 \
  --output artifacts/generated-results/quartic_lr_sparsity_search.json
```

It tests every rational exceptional parameter in 15 one-monomial source
shears and 15 one-monomial target shears, all 76 two-term representatives of
the essential source jet through exponent 12, and 25281 rational diagonal
scalings.  No searched nonidentity shear improves support `(7,51,38)`.
This is bounded computational evidence, not an absolute minimum.  The
scaling search finds the balanced exact height improvement
`(alpha,beta)=(1/4,12/5)`.  The generated record has SHA-256
`3fea20be042106fb5fe452ebe241dc5c3316eed6a893b00bf7ca2bcc0bef1b70`.

The continued two-move circuit search is:

```bash
.venv/bin/python scripts/search_quartic_lr_two_move_circuits.py \
  --source-degree 2 \
  --source-parameter-bound 4 \
  --jet-max-exponent 12 \
  --workers 4 \
  --output artifacts/generated-results/quartic_lr_two_move_circuits.json
```

It checks 330 rational source-shear/optimal-linear-target circuits and all
286 three-term essential jets through exponent 12, exactly at every rational
exceptional jet parameter.  The best nonidentity jet is
`z -> z+(16/7)y^2`, with support `(7,51,39)`.  Target monomial cleanup
through degree three does not improve it; among all exceptional elementary
second source shears through degree two, only the literal inverse returns to
96.  Adding a second structured monomial through exponent 12 also fails to
improve the 97-term near miss.  This remains bounded evidence.  The generated
record has SHA-256
`28c9e1c2ed9c765fef7c51d7e8ace3262c6fac337c2b11d723f3dccdb3781826`.

The remaining constant-`C` boundary has a separate exact Singular
certificate:

```bash
make verify-weighted-boundary
```

It computes exactly two primary components, checks their declared radicals,
and verifies that the reduced affine-three-space components meet in an
affine plane.

The reduced global attachment of the open torus orbit is checked by

```bash
Singular -q scripts/verify_foundational_reduced_gluing.sing
```

This verifies the degree-ten toric closure and its two boundary lines.

The heavier regression target also checks the explicit degree-five family and
its rank-two symplectic descent:

```bash
.venv/bin/python scripts/verify_degree_five_rank_two_descent.py
```

This exact calculation constructs the relative Hamiltonian over
`Q(lambda)`, extracts all four negative-`X` residue coefficients, proves the
unique parameter-dependent shear cancels them, and verifies the normalized
base brackets and polynomial source automorphism.  It normally takes roughly
half a minute in the pinned symbolic environment.

The smaller classical degree-drop viability test is:

```bash
.venv/bin/python scripts/verify_quartic_weighted_map.py
```

Besides the quartic inverse and collision, it identifies the seed with
\((\kappa,\tau)=(-5,0)\), specializes the rank-two completion, checks all six
Poisson brackets and the canonical coordinate change, and transports the
generic degree-four cover and an explicit two-point collision.  It reports
fiber orders \((4,3)\).  This is an exact classical certificate, not an
\(A_2\) quantization.

The rebuilt restricted quantization test for those exact \((4,3)\) symbols
is:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_quartic_degree_drop_quantization.py \
  --certificate artifacts/generated-results/quartic_degree_drop_quantization.json
```

It derives the specialized Bernstein bounds, solves the complete
parity-preserving \(\hbar^3\) affine equation, and writes an exact six-term
dual cocycle proving the rank jump \(143\to144\) at \(\hbar^5\).  It also
rebuilds the unrestricted first-order kernel, removes the complete
target-Hamiltonian gauge, projects the next Maurer--Cartan quadrics, and
tests all five surviving coordinate axes through the coupled
\(\hbar^2/\hbar^3\) equations.  Its bounded low-support pass classifies one
coordinate \(\mathbf P^4\), nine isolated rational directions, and no
algebraic support-two directions.  All nine isolated directions fail, while
a uniform third-order relaxation reduces the \(\mathbf P^4\) to an explicit
residual \(\mathbf P^2\).  Its genuine compatibility obstruction factors as
\((21a+28b+64c)^3/21^3\), so one rational \(\mathbf P^1\) reaches
\(\hbar^3\).  Parameterizing that line by
\((4,-3,0)+t(0,16,-7)\), the complete 38-dimensional lower-lift calculation
over \(\mathbb Q(t)\) gives the exact fourth-order rank jump
\(143\to144\).  Its six-term cocycle has sole denominator factor \(t\);
exact audits at \(t=0\) and projective infinity give the same jump.  Thus
the complete projective resonance line is eliminated at \(\hbar^4\).
The discrete base ranks are computed over \(\mathbb Q\) and independently
repeated over \(\mathbf F_{32003}\).  This certifies obstructions only for
the displayed normal ordering and inherited filtration; it is not an
\(A_2\) nonexistence theorem.  The recorded JSON certificate has SHA-256
`04646808c526697e7538a268605da6df1b5e3a66c51a6a5e3c1d68c80ab20ab9`.

The degree-five filtered contact problem has a separate two-invariant audit:

```bash
.venv/bin/python scripts/verify_degree_five_torus_module.py
```

It verifies the torus-gauge root recurrence over `Q[u,gamma]`, proves the
all-order profile `24m+1`, and checks survival of the candidate class in the
invariant-ring-saturated equivariant target quotient.

The minimal opposite-weight quadratic Rees witness is reproduced by

```bash
.venv/bin/python scripts/search_rees_torsion_witnesses.py --max-target-degree 0
```

It finds the unique constant opposite-weight pair
`(partial_B,partial_C)` of weights `(1,-1)`, computes its exact second
fundamental form by two independent formulas, and returns the nonzero leading
normal symbol `-146880u^5/7` in the third saturated summand `R/(gamma)`.

The finite LR Rees/SAGBI module calculation is reproduced by

```bash
.venv/bin/python scripts/compute_lr_rees_sagbi_modules.py
# or, including the dependency-free separator replay:
make verify-lr-rees-sagbi
```

It constructs the three-generator target-invariant SAGBI basis, the target
modules and initial lifts in weights `p=+-1,+-2`, and the saturated normal
quotient.  It certifies a linear weight-one degree drop `39 -> 34`, computes
the further subduction to a new degree-`29` initial-module generator, computes
the complete `3 x 24` quadratic matrices for `p=1,2`, proves the structural
cutoff `|p|>=3`, and performs exact Singular module membership.  The sole
new `p=2` column modulo the full `p=1` image is
`II_(F,2,-2)(partial_A,A^2 partial_A)`, with remainder
`-987/395*e_C`.  The generated JSON certificate is
`artifacts/generated-results/lr_rees_sagbi_module_computation.json`.

The decisive independence statement has a dependency-free replay:

```bash
python3 scripts/audit_lr_rees_sagbi_module_certificate.py
```

The all-order constant-direction rooted-tree normal classes are reproduced by

```bash
.venv/bin/python scripts/compile_lr_rooted_tree_classes.py --max-order 12
python3 scripts/audit_lr_rooted_tree_normal_classes.py
```

The compiler works in exact torus semi-invariant coordinates, reproduces the
known `II_F(partial_B,partial_C)` residue at order two, and constructs the
weight-zero ladders `tau_2=B(C)`, `tau_3=A(C(C))`,
`tau_(n+2)=B(C(tau_n))`.  A fixed `3 x 3` transfer matrix at
`(u,gamma)=(1/6,0)` and a positive-coefficient Cayley--Hamilton recurrence
prove that the third saturated normal residue, hence its associated-graded
symbol, is nonzero for every order `n>=2`.  This is an all-order theorem for
the individual tree classes, not a proof that the same class survives the
sum and lower-jet variation in a mixed BCH/LR forcing coefficient.

The balanced linear-in-`X` mixed BCH sector is reproduced by

```bash
.venv/bin/python scripts/compile_lr_mixed_bch_classes.py --max-k 3
python3 scripts/audit_lr_mixed_bch_classes.py
```

Here `X=N*(x,0,-3z)`, `D_B=ell_F(partial_B)`, and
`D_C=ell_F(partial_C)`.  The checker proves `[D_B,D_C]=0`, collapses the
balanced order-`2k+1` BCH sum to
`binomial(2k,k)*(ad(D_B)*ad(D_C))^k*X`, and derives the exact third-normal
recurrence

```text
c_(k+1) = -73440*(k+3)*(2k+7)*c_k,
c_1 = 14438891520/2401.
```

Thus the actual multihomogeneous BCH coefficient is nonzero in every odd
order.  It survives the saturated linear target quotient, but with target
amplitudes `s,t` it is multiplied by `s^k*t^k`; consequently this sector
alone is not universal over the full lower-jet scheme.

At `(u,gamma)=(1/6,0)`, the covector `(0,-144/79,1)` descends through the
saturated normal relations, kills all 24 `p=1` columns, and takes the value
`-987/395` on the new `p=2` column.  This matches the Singular normal
remainder without using SymPy or Singular in the replay.  The main checker
also computes the exact annihilator `(gamma,6*u-1)`, proving that the
`p=2` image modulo `p=1` is one reduced residue-field copy of `Q`.

The full normalized degree-five seed surface is checked by

```bash
.venv/bin/python scripts/verify_degree_five_flux_surface.py
```

This exact two-parameter calculation works over `Q(a,tau,s_2)`, verifies the
uniform adapted coordinate and quotient brackets, extracts the complete four
term Laurent obstruction, and proves that its unique quadratic shear makes
the Hamiltonian polynomial.  It takes several minutes in the pinned symbolic
environment.

The exceptional `kappa=-1` chart and its pole-filtered monomial shear
responses through cubic degree are replayed by

```bash
for degree in 0 1 2 3; do
  .venv/bin/python scripts/explore_kappa_minus_one_flux.py --shear-degree "$degree"
done
```

Each run verifies the replacement determinant and quotient brackets, all
three Hamiltonian components, and every negative-`X` residue coefficient.
The full exceptional-divisor completion is checked by

```bash
.venv/bin/python scripts/explore_kappa_minus_one_flux.py \
  --x-degree 1 --shear-degree 1
```

It proves that the complete principal part is canceled by
`2(2*tau^2-15*tau-18)*X*Q/105`.

The full degree-six generic chart, exceptional divisor, and fixed-`gamma`
specialization are checked by

```bash
.venv/bin/python scripts/verify_degree_six_flux_surface.py
.venv/bin/python scripts/verify_degree_six_kappa_minus_one_descent.py
.venv/bin/python scripts/verify_degree_six_fixed_gamma_descent.py
```

These verify the three-parameter generic seed chart and the full exceptional
divisor, componentwise Hamiltonian identities, complete residues, and unique
completing shears.  The generic symbolic replay is a heavy calculation.

The all-degree Laurent recurrence and exact fixed-`kappa=-9` probes in degrees
seven and eight are checked by

```bash
.venv/bin/python scripts/verify_four_residue_recurrence.py
.venv/bin/python scripts/explore_all_degree_fixed_gamma.py 7
.venv/bin/python scripts/explore_all_degree_fixed_gamma.py 8
```

The direct second-Weyl-algebra parity test is replayed by

```bash
.venv/bin/python scripts/explore_degree_five_a2_subprincipal.py
```

It solves the `hbar^3` equation exactly, retains its full 42-dimensional
solution space, and proves that the `hbar^5` cokernel contains `1=0`.  This is
an obstruction only to the declared parity-preserving filtered ansatz.

The parameter-uniform third-order lift, four bounded fifth-order periods,
their common cubic locus, and the genuine nonlinear fifth-order equations on
that locus are replayed by the commands in
[`extended-geometry/QUANTUM_RESIDUE_OBSTRUCTION.md`](extended-geometry/QUANTUM_RESIDUE_OBSTRUCTION.md#10-reproduction).
The small relative-family package, including the two exact coprime Fitting
charts, the interior Kuranishi shadow, and the root-at-infinity valuation
filtration, is checked independently by

```bash
make verify-degree-five-relative-quantization-family
```

It verifies the valuation weights
`(X,Q,W,R,gamma)=(1,-1,-2,1,0)` and the induced pure correction weights
`(S_2,T_2)=(4,5)` and `(S_4,T_4)=(10,11)`.  It does not promote the modular
length-218 Fitting computation to characteristic zero.
The decisive exact cubic-field check is

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --exact-cubic --print-radical-basis --seventh-line
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --seventh-component-elimination
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --exact-cubic --seventh-component-elimination \
  --seventh-component-program-output \
    artifacts/generated-results/degree_five_cubic_h7_unit_certificate.sing
.venv/bin/python \
  scripts/verify_degree_five_cubic_h7_unit_certificate.py
```

It verifies all 680 projected quadratic equations at an explicit
two-coordinate lower lift and solves the unreduced fifth-order correction
equation with a 14-term particular solution.  The exact radical is one
affine 27-space with a six-linear-form nonlinear core, and the entire
explicit one-parameter line is obstructed at order seven.  The second
command proves over `GF(32003)` that the full 20-column order-seven matrix
has constant rank six and that its 401-polynomial consistency ideal, in
only ten effective parameters, is the unit ideal.  Repeat it with
`--prime 31991 --a 109 --tau 28672` for the second good-prime certificate.
The third command performs the characteristic-zero lift over the cubic
field.  Batched constant-field elimination replaces the former
27-variable function-field solve.  The resulting 401 equations contain 27
nonzero constants; the selected \(X^{18}\) residual has an explicit Bézout
inverse, and Singular verifies both a direct one-generator identity and a
one-term degree-zero lift of \(1\).  The final command replays the pinned
1.3 MB Singular certificate in under a second.  Hence the complete reduced
affine 27-space of fifth-order lifts is obstructed at order seven.

The bounded audit of the standard-support parameter Fitting scheme is

```bash
.venv/bin/python scripts/compute_degree_five_qper_fitting.py \
  --timeout 120
```

It finishes at the three fixed good primes, checks the common
21-generator leading-monomial staircase, saturation exponent 12, dimension
zero, and length 218.  This is a stable modular certificate, not a
characteristic-zero proof.  The opt-in rational reconstruction experiment is
checkpointed by

```bash
.venv/bin/python scripts/compute_degree_five_qper_fitting.py \
  --prime 0 --method modular-rebuild --timeout 900 \
  --basis-output \
    artifacts/generated-results/degree_five_qper_fitting_basis_Q.sing

.venv/bin/python scripts/verify_degree_five_qper_fitting_basis.py \
  --check shape
.venv/bin/python scripts/verify_degree_five_qper_fitting_basis.py \
  --check groebner --jobs 8 --timeout 600
.venv/bin/python scripts/verify_degree_five_qper_fitting_basis.py \
  --check boundary-unit --timeout 600
```

The rebuilt 20,840,615-byte rational candidate has SHA-256
`25788668021f563e17373b55703a08ef5693576077ebdbe53c4c3f2c659d98e6`.
The 20 adjacent staircase \(S\)-pairs give an exact Gröbner certificate, and
the exact boundary-unit check proves that the candidate itself is saturated.
These checks do not prove that it equals the saturated maximal-minor ideal.
That last identification requires fraction-free quotient identities in both
containment directions.

The bounded modular audit of the 16 input-containment quotients is built
incrementally by:

```bash
.venv/bin/python \
  scripts/reconstruct_degree_five_qper_input_quotients.py \
  --no-resume --skip-diagnostic

.venv/bin/python \
  scripts/reconstruct_degree_five_qper_input_quotients.py \
  --primes 70001 70003 70009 70019 70039 70051 \
    70061 70067 70079 70099 70111 70117 \
  --jobs 8 --checkpoint-every 64 --skip-diagnostic --timeout 600

.venv/bin/python \
  scripts/reconstruct_degree_five_qper_input_quotients.py \
  --prime-start 1000000000 --prime-count 100 \
  --jobs 8 --checkpoint-every 64 --skip-diagnostic --timeout 600

.venv/bin/python \
  scripts/reconstruct_degree_five_qper_input_quotients.py \
  --prime-start 1000002043 --prime-count 500 \
  --jobs 8 --checkpoint-every 64 --timeout 600
```

It resumes
`artifacts/generated-results/degree_five_qper_input_quotients_modular.json`,
uses deterministic prime order even with parallel Singular workers, and
keeps one good image out of the CRT pool.  The recorded run has 613
support-stable good primes, two support-unlucky primes, 18,116 CRT bits, and
30 of 11,701 balanced reconstructions confirmed at the held-out prime.
These data are modular evidence and a coefficient-height diagnostic, not an
exact containment proof.  Use `--skip-diagnostic` when only extending the
checkpoint; the held-out reconstruction pass uses GMP-backed FLINT
arithmetic.  The recorded compact checkpoint has SHA-256
`70a690fd53b4b3a15d4eebf5116acf57b7d0079a8f96a1aadfb2826da86d0481`.

The low-support unrestricted odd audit is replayed by

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_mixed_quantization.py
```

It reconstructs the 38-dimensional gauge quotient and its 41 quadratic
obstruction equations and classifies every exact-support-two branch,
including the nine quadratic closed points.  The only three mixed
support-two directions reaching the simultaneous third-order equation retain
63 lower-lift parameters; after adding all 2079 enlarged obstruction
coefficients and every bounded next correction, the constant raises the span
rank from 626 to 627 over each of the good primes 31991, 32003, and 65521.

The generic residual-line and exact support-three audits are replayed by

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_rank_two_odd_mixed_function_field.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_support_three.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_rank_two_odd_support_three_points.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_support_three_curves.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_rank_two_odd_support_three_curves.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_residual_five_space.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_rank_two_odd_residual_support_three.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_rank_two_odd_residual_fourth_identity.py
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_normal_cones.py L1 --exact
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_normal_cones.py L2 --exact
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_normal_branch_lifts.py L1 --order 8
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_normal_branch_lifts.py L2 --order 8
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_rank_two_odd_l1_high_support.py
```

The first constructs a fixed three-monomial residue over `QQ(r)`, finds its
sole lower-basis pole at `r=-3/4`, and eliminates that exceptional
specialization exactly; hence the complete residual projective line is
closed.  The next two commands classify all 8436 support-three coordinate
charts and eliminate all 66 isolated closed-point classes.  The curve
commands compress the 149 positive-dimensional line/conic charts to 23
closed points, force all 19 quadratic classes to zero scale, and eliminate
the four rational survivors by exact `646->647` rank jumps.  Thus every
exact-support-three branch is closed.

The last three commands treat the residual projective four-space inside `L_2`.
The first proves exactly that its nonzero-scale locus is the union of the hyperplane
`z2+2*z3-9838*z4/105=0` and one explicit primitive quadric.  Uniform
fourth-order obstruction on those two threefolds is supplied by the last
command: the quadric is a binary form in the two coupling coordinates and
splits into conjugate hyperplanes over `QQ(sqrt(-2))`; a fixed three-term
residue handles every nonzero-coupling chart, and an exact 16-term residue
handles their rank-zero intersection plane.  The middle command
intersects that locus with every residual exact-support-three chart: the
already-closed residual line and twelve closed points result, and all twelve
points have exact `626->627` next-order rank.  Thus exact support three is
closed even inside the coordinate planes contained in `L_1 union L_2`, and
the final command eliminates every nonzero-scale branch in the residual
projective four-space, in every support.  The final five commands attack
higher support: they compute the exact `L1` and `L2` normal cones; show that
a generic `L1` normal branch is an exact 26-support classical solution while
the analogous `L2` branch is obstructed at its next Kuranishi equation; and
isolate a projective high-support `P6` that reaches nonzero filtered scale.
Uniform fourth-order continuation on that `P6` is the next open calculation.
These are statements about the displayed classical symbol only; they do not
prove `(DC_2)`.

The separately authored foundational Lean certificate is optional because it
downloads an external pinned checkout:

```bash
make verify-lean-foundational
```

GitHub Actions runs this target in the required `formal-lean` job using the
pinned upstream commit and Lean action.  The `formal-local-lean` job
audits both publication certificates and runs `make verify-lean-local`, so
every checked-in Lean package is built.  The `papers` job compiles the
finalized and active manuscripts listed in `papers/README.md`, while parked
manuscripts remain available for direct local builds.  The
`macaulay2-independent-check` runs the pinned Macaulay2 comparison.  Together
with the four Python matrix jobs, these are the complete CI verification
pipeline.  The final `verification-complete` job is the single aggregation
check intended for GitHub branch protection.

## Cancellation programme

```bash
make verify-master
```

The direct three-puncture reciprocal ledger and its all-degree
polynomiality obstruction are replayed separately by:

```bash
.venv/bin/python scripts/verify_puncture_rank_frontier.py
```

The checker derives the universal determinant ledger, enumerates all
two-character coefficient matrices in `[-2,2]` into 129 rank failures,
392 nonsaturated class-lattice failures, and 104 saturated bases, and lists
the 44 primitive positive `(r,a,b)` ledgers with coefficients at most four.
It then verifies the boundary-moment eigenvalue recurrence used by the
written all-degree proof and retains the eleven degree-four through
degree-seven coefficient factorizations as regressions.  No candidate
reaches polynomiality, so the calculation does not claim a Keller map or a
complete collision.  Singular is not required.

The surviving two-reconstruction-variable `A^6` core and its nonlinear
screens are checked by:

```bash
.venv/bin/python scripts/verify_three_puncture_nonlinear_frontier.py
```

In addition to the two dimension-free rank-drop gates and the 80 fixed
degree-at-most-three skeletons, this computes the zero-modification slice of
the proposed coupled ansatz.  For arbitrary affine `P,Q` over `Q(c,v)`,
arbitrary polynomial `H,S`, and two arbitrary affine transverse outputs, it
forms eleven determinant coefficient equations and the Plücker quadric.
Their exact Gröbner basis is `(1)`.  This proves that the next search must
give at least one transverse output degree at least two; it does not exclude
the remaining degree-two through degree-four systems.  The checker then
keeps `P,Q` arbitrary affine, makes the fourth output a completely general
degree-at-most-two polynomial, and proves unit coefficient ideals for the
eight transverse skeletons
`u,z,w,u+z,u+w,z+w,D0+z,D1+w`.  The two-general-quadratic-output system
remains open.  Uniformly for every nonconstant affine direction `C`, four
projective pivot trees avoid a slow monolithic Gröbner calculation.  The
nonzero-`r` chart `C=r+g*u+a*z+b*w` follows the exceptional divisors
`p1,q1,p2,q2,g,p0-a`; the `r`-free charts
`C=u+a*z+b*w,z+b*w,w` have chains `(p1,q1)`, `(q1-b*p1)`, and `(p1)`.
The coefficient and augmented ranks differ on every open and terminal
branch.  Thus no affine `C` can be paired with a general
degree-at-most-two fourth output.  Arbitrary quadratic `C`, two
simultaneously general quadratic outputs, and degree-at-least-three fourth
outputs remain open.

The checker also closes the exposed-`r` simultaneous-quadratic boundary
`P=Q=0`: for two general degree-at-most-two outputs, an eight-step
coefficient pivot tree ends with coefficient/augmented ranks `6/7`.
Degree three is the first zero-slice escape.  With `q=1-c*v` and `C=w`, it
verifies
`D3=(u*(-q^2-v*r*q+2*v^2*r^2)-6*v^2*z)/q^3` has slice determinant one.
The polynomial numerator instead has full determinant
`q^3+6*r*v^2*w`; its cofactor derivation has an explicit common zero.
Thus this exact rational survivor neither satisfies polynomiality nor lifts
to a full Keller map, and a next quadratic/cubic search must use nonzero
`P,Q` or `H`.

This target includes the exact quadratic-gauge/cancellation intersection
regression.  To run its symbolic `N=4,5,6,7` discriminant and all-factorization
checks directly:

```bash
.venv/bin/python scripts/verify_quadratic_cancellation_intersection.py
```

The all-rank clean quadratic-gauge stable-moduli and marked-stabilizer
certificate is:

```bash
.venv/bin/python scripts/verify_quadratic_gauge_stable_moduli.py
```

Besides the two-torus quotient and its saturated compiler-slice invariants,
it verifies the weight-one global receiver slice `lambda=u_5/u_4`, the
exact finite-etale descent identity, the universal discriminant inequalities
through rank 128, and exact discriminant supports in ranks four through
eight.  These checks certify the all-rank written proof that
`D_N=(2,N,1)`, corresponding to `P^2*B^N*C`, is uniquely exposed by
`(1,N+1,N)`.  Kuroda's and
Derksen--Hadas--Makar-Limanov's theorems are external mathematical inputs,
not re-proved by the script.

The minimal-boundary gateway and classification program has a separate fast
cubic certificate:

```bash
make verify-minimal-boundary
```

The eight-predicate invariant pipeline on finite canonical-normalization
exports has its own dependency-light target:

```bash
make verify-minimal-boundary-pipeline
```

It checks weighted degrees `3,...,8`, six cancellation parameter pairs,
quadratic-gauge degrees `3,...,8`, five single-defect perturbation/spectator
records, and a relabeling/order blindness regression.  Regenerate
`artifacts/generated-results/minimal_boundary_pipeline.json` with

```bash
.venv/bin/python scripts/verify_minimal_boundary_pipeline.py \
  --write-artifact
```

This target starts from exact finite exports; it does not compute a canonical
normalization or its intrinsic marking from a bare polynomial map.

It proves that the weighted geometric-degree-three seed has no modulus, that
the cancellation degree equation forces `(m,r)=(1,1)` and `h=3+9A`, and that
both maps are carried to the foundational polynomial by explicit diagonal
source and target automorphisms.  It also verifies the cubic two-place toric
defect atlas and the diagonal reciprocal-lift obstruction.  The accompanying
proof uses Abhyankar--Moh to make the one-place plane-core marking automatic.
It also checks the positive quotient tower and the target-polynomiality jet
that forces `gamma=1-3xy/2 mod x^2`; the written LND/Stein argument supplies
the slice under explicit intrinsic saturation labels.  On the reciprocal
side it checks the coefficient valuations `(n-1,2n-1)` and the extraction of
`Y=Q-Ps` from a primitive quadratic conormal coefficient.  The eight
minimal-boundary predicates are formalized in the accompanying note, but
this checker does not construct their finite-normalization witness, verify
`PC`, `NC`, or `CS` for an arbitrary boundary-minimal map, or extract a
suspension from the unmarked canonical normalization.

The same target also checks the finite-normalization frontend: the
Deligne--Faddeev cubic-algebra table and discriminant, the codimension-three
reflexive-module warning and its minimal excess-length-four special fiber,
whose exact module-theoretic defect is `Fitt_3=(x,y,z)`, the unique
critical-divisor DVR budget `(2,1)+(1,1)`, and the
tangent-hyperplane quotient coordinates.  The written local argument proves
that cubic point-flatness is equivalent to every canonical scheme fiber
having length three.  The cited nonflat triple-cover correspondence also
shows that normal cubic algebra structure alone cannot remove this defect.
The local structure theorem further identifies every defect with an
`(s+2)`-by-`s` determinantal presentation, where the excess fiber length is
exactly `s`; the checker includes the origin-primary `s=2`, length-five
rung in addition to the minimal Koszul rung.
For a reduced minimal defect it also verifies the linear-algebra inputs
forcing the unique square-zero fiber `k plus k^3`.  The written corollary
then identifies such a defect with a closed-point collision of the
ramified boundary sheet and the affine sheet over the critical divisor.
It distinguishes this from the allowed foundational collision, whose
triple-root fiber is curvilinear of length three.
The maximal-minor order argument proves that every reduced defect is
automatically this minimal Koszul rung; only nonreduced Fitting defects
remain outside the square-zero classification.
The local monogenicity theorem then closes all of those cases
simultaneously under intrinsic curvilinearity of the collision fibers:
Nakayama lifts a fiber generator and the resulting monic cubic algebra is
free.
Before either saturation test, the written tame-local proposition removes
every simple-normal-crossing point of the critical discriminant: after
strict henselization the cubic normalization is the finite-free sum
`R[s]/(s^2-t_1...t_r) plus R`.  Consequently the point-defect computation
first reduces to closed non-SNC points of the discriminant.  The written
ordinary-cusp proposition then classifies the two possible three-sheet
braid representations.  Equal meridian transpositions give the finite-free
`2+1` Kummer algebra, while distinct transpositions give the finite-free
monic cubic root cover.  Thus only worse-than-ordinary-cusp points remain.
For a reduced Koszul defect, the next written proposition identifies the
projectivized branch tangent cone with the discriminant of line sections of
the ternary cubic `h`.  The frontend checker verifies the complete
degree-six factor table for smooth, nodal, cuspidal, conic-plus-line,
triangle, and concurrent-line symbols, together with vanishing for double
and triple components.  Thus every reduced defect forces branch
multiplicity six, or at least seven in the non-squarefree case.
The frontend checker verifies that the foundational discriminant's singular
ideal is the expected triple-root locus.  The local cusp models and all nine
three-letter braid pairs are checked by:

```bash
.venv/bin/python plane-jc/cas/test_cubic_cusp_local_model.py
```

The remaining nonzero ternary-cubic symbol strata and both canonical
saturation modules are audited by:

```bash
.venv/bin/python scripts/verify_cubic_symbol_double_saturation.py
.venv/bin/python scripts/verify_cubic_symbol_deformation_saturation.py
.venv/bin/python scripts/verify_cubic_symbol_quartic_tangent_saturation.py
.venv/bin/python scripts/verify_smooth_cubic_quartic_plane_saturation.py
.venv/bin/python scripts/verify_singular_cubic_quartic_plane_saturation.py
.venv/bin/python scripts/verify_smooth_cubic_quartic_three_space_saturation.py
.venv/bin/python scripts/research_universal_cubic_quartic_kernel_saturation.py
.venv/bin/python scripts/verify_universal_cubic_filtered_syzygy_frontier.py
.venv/bin/python scripts/verify_cubic_quartic_ext_tail_absorption.py
.venv/bin/python scripts/verify_universal_cubic_quartic_different_complex.py
.venv/bin/python scripts/verify_universal_cubic_kahler_annihilator.py
.venv/bin/python scripts/verify_cubic_symbol_dense_quartic_plane_saturation.py
.venv/bin/python scripts/verify_cubic_symbol_affine_dense_quartic_plane_saturation.py
.venv/bin/python scripts/verify_universal_cubic_cotangent_saturation.py
.venv/bin/python scripts/verify_cubic_formal_gauge_cokernel_atlas.py
.venv/bin/python scripts/verify_nodal_cubic_formal_slice.py
```

The last checker also quotients the degree-five curvature by the
five-dimensional kernel of the quartic gauge lift.  It verifies a
rank-four action on the six slice--gauge coefficients, extracts two
intrinsic cross forms and three intrinsic pure-gauge quadrics, and checks
that the pure-curvature zero scheme is two reduced rational planes with
one embedded quadratic socle class.  Relative to the stored quartic lift,
it continues both reduced planes through degree six, constructs the exact
quadratic corrections, and obtains the two cubic Veronese classes.
Quartic-lift independence and the embedded-socle continuation remain open.

For the homogeneous tensor, all seven squarefree strata have saturated
cotangent presentation and a length-six `Ext_A^2(T,A)` support defect with
Hilbert function `3+3t` (three-dimensional top and zero `m^2` action);
double and triple lines instead have a one-dimensional support defect.
The zero homogeneous tensor passes both module tests but is nowhere
generically étale.  One explicit order-four kernel tensor makes the support
defect finite of length six in all ten strata, while cotangent saturation
still passes.  This is an exact leading-model computation.  It neither
proves lift-independence nor constructs a normal lift with a Keller open.
The support presentation is computed directly as the module preimage
`modulo(H1,H2)`, namely the kernel of the four action columns in the
threefold direct sum of the cotangent quotient.  This is exactly equivalent
to extracting four coordinates from the combined syzygy module, but avoids
the unnecessarily large 97-column elimination.
The second command works over `Q[t,x,y,z]`.  On each of the seven
squarefree lines `phi_h+t*psi_4`, it verifies uniform cotangent saturation,
no parameter torsion in relative `Ext_A^2(T,A)`, radical support equal to
the collision axis, multiplicity six, and equality of the relative
presentation with the scalar extension of its central specialization.
This proves constancy on those lines, not on the full 24-parameter
order-four space.
The filtered-syzygy frontier command resolves the smooth central
unit-pruned cotangent presentation with cokernel ranks `7 -> 13 -> 6`,
checks that `x+y+z` is regular there, and then applies the 24 unchanged
central input syzygies to the universal 6-by-25 matrix.  Twelve exact
remainders survive modulo the central image.  This certifies that
entrywise collision-order growth and two-jet agreement do not themselves
give a coefficient-independent Rees-strict resolution.  It does not prove
boundary torsion or disprove the smooth universal saturation theorem,
which the formal-gauge command below now proves.
The quartic-tangent command is the longer tangent-direction audit.  It
tests every
one of the 24 exact nullspace-basis axes for every squarefree symbol.  All
168 families have uniform cotangent saturation, no parameter torsion,
collision-axis radical support, and relative multiplicity six.  Literal
presentation equality changes in four rows, without changing any of those
invariants.  The basis axes span the order-four kernel, but this computation
does not test all their linear combinations.
The fourth command is the four-worker mixed-direction audit for the smooth
symbol.  Over `Q[u,v,x,y,z]` it tests all 276 full coordinate planes.  On
every plane the cotangent presentation is saturated and the relative
length-six `Ext^2` presentation is pulled back from `u=v=0`; this includes
every specialization on the plane.  Directions supported on three or more
basis tensors are not tested.
The fifth command is the four-worker coordinate-plane audit for the six
singular squarefree symbols.  It verifies that 1,652 pruned presentations
are pulled back from the origin.  Four ambient presentations jump, so the
checker forms their exact finite `Q[p0,p1]` presentations using the
verified zero `m^2` action and proves `Fitt_6=(1), Fitt_5=(0)`.  Thus all
1,656 singular-squarefree planes are flat of relative rank six; the four
exceptional modules are in fact free by Quillen--Suslin.
The sixth command is the longer four-worker three-space audit.  It tests
all 2,024 smooth coordinate three-spaces over
`Q[p0,p1,p2,x,y,z]`.  After pruning contractible free summands, every
relative rank-three `Ext^2` presentation is pulled back from the parameter
origin with multiplicity six.  Directions supported on four or more basis
tensors are not tested.
The seventh command is an exact frontier calculation, not a universal
theorem.  It checks four full-support dense lines for every squarefree
symbol and the full first-ten-coordinate subspace for the smooth symbol.
All 28 lines have uniform cotangent saturation, no parameter torsion,
collision-axis support, multiplicity six, and central Ext presentation.
The smooth parameter ten-space has the central pruned rank-three
presentation.  It also constructs the full universal cotangent matrix,
checks that its parameter-dependent terms have bidegrees `(1,3)`, `(1,5)`,
and `(2,6)` in parameters/collision variables, and removes six
parameter-independent unit pivots to obtain a cokernel-equivalent
6-by-25 presentation.  This command is an exact input reduction, not a
saturation or Ext calculation over all 24 parameters; the final command
below supplies the universal theorem by a different method.
The eighth command computes the last nonzero differential in the minimal
support resolution on the seven full-support squarefree planes.  Six rows
are parameter-independent and linear; their transposes present a
length-six quotient killed by `(x,y,z)^2`.  The seventh row has a central
quadratic part and a parameter-linear cubic part, so it reduces to zero
modulo those six rows.  The resulting 12-generator parameter module has
six independent constant relations and therefore
`Fitt_6=(1), Fitt_5=(0)`.  It also verifies that the seven canonical
different generators equal the complete annihilator on each plane.
The ninth command constructs the canonical different matrix
`[(0,z,-y,x),(s_ij,2*mu_ij)]` over all 24 quartic-kernel parameters and
all seven squarefree symbols.  Its explicit universal syzygy matrix
satisfies the Buchsbaum--Eisenbud grade conditions.  The resulting
canonical-different support has constant length-six `Ext^2` and
`Fitt_6=(1), Fitt_5=(0)` over the full parameter ring.  Identifying these
seven generators with the complete annihilator `Ann(Omega)` is conditional
at this stage; the eighth command verifies that equality on the seven
full-support planes, and the final command closes it universally.
The tenth command computes the universal Deligne--Faddeev locally free
cubic algebra and proves that the Kähler different `Fitt_0(Omega)` equals
the full annihilator `Ann(Omega)`.  On the punctured Koszul base this
identifies the canonical different with the actual support ideal.  Together
with the depth of the ninth command's exact complex, relative cotangent
saturation extends the equality across the collision axis and closes the
actual Ext Fittings.  The command does not prove that remaining universal
cotangent saturation; the final command does.
The eleventh command tests one low-height full-support plane for all seven
squarefree symbols.  If `psi_plus` is the sum of the 24 fixed kernel-basis
tensors and `psi_minus` their alternating sum, it computes the complete
family `phi_h+u*psi_plus+v*psi_minus` over `Q[u,v,x,y,z]`.  In every row
the cotangent presentation is saturated and the pruned rank-three relative
Ext presentation is pulled back from the origin with multiplicity six.
On `u^2-v^2!=0`, all 24 basis coordinates are nonzero.  This proves a
two-parameter full-support result, not a Zariski-open theorem in the
24-dimensional kernel.  Low-height dense parameter-three-space and
higher-height parameter-four-space runs reached their declared 600-second
timeouts and provide no mathematical evidence.
The twelfth command translates the same sum/alternating-sum plane by the
deterministic generic quartic lift and tests all nine nonzero cubic-symbol
orbits plus the zero symbol.  Over `Q[u,v,x,y,z]`, all ten rows have
saturated cotangent presentation, support exactly equal to the parameter
plane, relative Ext multiplicity six, and pruned rank-three presentation
pulled back from `u=v=0`.  The checker also reconstructs the generic lift
in the fixed primitive 24-element basis and verifies that every coordinate
is nonzero.  This proves uniform purity restoration for the double-line,
triple-line, and zero symbols on one affine plane; it does not prove
normality or Keller-open compatibility.

The final command proves the smooth-symbol 24-parameter theorem without
computing the universal saturation.  For the graded module `K` of all
compatible tensor corrections and the exact `10`-by-`9` simultaneous
coordinate/coefficient gauge matrix `G`, it verifies
`K=im(G)+A*eta` and `(x,y,z)*eta subset im(G)`.  An explicit matrix `L`
satisfies `G*L=[x*eta,y*eta,z*eta]`, so every compatible tensor term of
collision degree at least four is gauge.  It independently derives all
nine columns of `G` by expanding the determinant-twisted finite action over
the dual numbers.  It also stores an explicit linear-polynomial
`9`-by-`24` matrix `Q` with `G*Q=[psi_1,...,psi_24]`; both the quartic
compatible space and this gauge image have rank `24`, with gauge kernel
dimension three.  Successive homogeneous changes formally identify the
universal quartic family with its saturated central fiber.  Since
completion detects `(x,y,z)`-power torsion, this proves
`H^0_(x,y,z)(Omega)=0`; the canonical-different argument then gives the
universal annihilator equality and actual-support
`Fitt_6=(1), Fitt_5=(0)`.
The atlas command derives the determinant-twisted gauge differential over
the dual numbers for all ten ternary-cubic symbols and computes the exact
graded modules `ker(C)/im(G_h)`.  Their Hilbert series prove that smooth is
the unique symbol formally rigid above collision degree three.  The exact
quartic nongauge dimensions are `0`; `2,4,4,6,6,8` on the six singular
squarefree symbols; and `11,16,24` on the double-line, triple-line, and
zero symbols.  It also proves the exact singular-squarefree annihilator
sequence `(x),(x^2),(yz),(y^3),(xyz),(x^3)`; the three non-squarefree
quotients have zero annihilator and generic ranks `1,2,4`.  These data
delimit the formal-triviality method; they do not assert failure of
cotangent saturation.
The final nodal command proves the cyclic refinement
`ker(C)/im(G_nodal)=Q[y,z](-3)` with generator given by the tensor of
`Z^3`.  In quartic degree the 24-dimensional compatible space splits as a
22-dimensional gauge image plus the slice generated by `y*eta,z*eta`.
Only the first two fixed quartic basis directions survive in the quotient.
The sum/alternating-sum plane is a second transverse slice, with
change-of-slice determinant two.  Both slices replay the saturated
cotangent presentation and constant length-six Ext block.  This is a
first-stage slice theorem.  For the stored row-reduced quartic gauge lift,
the command also computes the complete degree-five normal curvature in the
basis `y^2*eta,y*z*eta,z^2*eta`: its components have `14,16,13` quadratic
terms and 30 nonzero cross-parameter pairs.  It vanishes on the coordinate
slice and has the explicitly recorded nonzero restriction on the dense
slice.  The command then quotients changes in the five-dimensional
gauge-lift kernel: its rank-four action on the six slice--gauge
coefficients leaves two intrinsic cross forms, while the three
pure-gauge quadrics are already invariant.  Their reduced zero scheme is
two rational planes, and their unreduced ideal has one embedded quadratic
socle class.  On both reduced planes the checker constructs an exact
quadratic correction of the degree-five term.  Its 15-dimensional
ambiguity acts trivially on the degree-six quotient, and the resulting
classes are `27/8*(q*y+p*z)^3*eta` and
`27/8*(q*y-p*z)^3*eta`.  Independence from the earlier quartic-lift
ambiguity and continuation of the embedded socle remain open.

The equivalent coordinate-free test is that each collision cotangent module
has unit first Fitting ideal (or vanishing second exterior power); the
checker separates the cyclic triple-root cotangent from the three-generator
square-zero cotangent.
The equivalent nilradical test has one generator and nilpotency index three
for the foundational collision, versus three generators and index two for
the reduced defect.
The written Hartogs extension theorem proves that a primitive cotangent
generator in codimension one extends through closed collisions whenever the
pure two-dimensional ramification support is `S_2` and its rank-one
cotangent module is `S_1`.  The companion two-`Ext` theorem identifies the
only closed-point obstruction modules as `Ext_A^2(T,A)` and
`Ext_A^3(Omega_{B/A},A)`.  Its double-saturation refinement forms the
canonical `S_2` hull `C=Ext_A^1(Ext_A^1(T,A),A)` and identifies those
obstructions successively with the canonical duals of `C/T` and
`Omega_{B/A}/T tau`.  The coupled local-cohomology sequence shows that
after `C=T` the latter is exactly the closed-point torsion of
`Omega_{B/A}`.  If `N` is the image of a free presentation and
`I=Fitt_3(B)`, the exact test is `N:I^infinity=N`; the Singular regression
checks this module saturation directly.
The phantom-boundary theorem identifies the quotient between the reduced
nonproperness and branch equations as the exact extra-divisor detector.
The checker calibrates it on the foundational map: boundary elimination and
the cubic discriminant give the same irreducible equation, so the quotient
is one.  The written boundary-minimality corollary then closes this
certificate for every boundary-minimal cubic: the foundational competitor
gives upper bound one for the number of target boundary components, while
nonproperness gives the matching lower bound.  Thus no second unramified
target divisor remains in the minimality problem.
This is not a global arbitrary-cubic closure.  Proposition 1.4 shows that a
second boundary sheet cannot lie over the critical divisor, because the
ramified `(2,1)` and affine `(1,1)` sheets exhaust degree three.  An
arbitrary cubic can still have a distinct unramified nonproperness divisor;
excluding that factor, or reducing it to the minimal stratum while
preserving genuine ungradedness, is the separate `OP-UG3` obligation.

The universal flat ungraded coefficient cell is checked separately by:

```bash
.venv/bin/python scripts/verify_universal_cubic_ungraded_testbed.py
```

It verifies the seven-parameter degree-at-most-four cell, and the written
argument extends the same identities to arbitrary
`A(P),gamma(P) in k[P]`.  The checker proves the determinant-minus-two
identity, inverse cubic and derivative reconstruction, reciprocal chart,
finite free Deligne--Faddeev multiplication table, discriminant, smooth
Laurent ramification parametrization, universal polynomial `GL_2`
discriminant transformation, and the exact equivalence
`G automorphic <=> phantom factor unit`.  It also separates this flat cell
from the 24-dimensional Koszul order-four tensor kernel by their intrinsic
third Fitting ideals.  It does not construct a Keller open for an arbitrary
Koszul-kernel combination.

The written no-global-monogenicity proposition then shows why these local
generators cannot be patched into one root coordinate: the derivative would
be a constant unit on `A^3` and would contradict cubic degree.
The written theorem proves uniqueness without a supplied suspension when
the intrinsic flatness defect is empty, the binary-cubic coefficient map is
affine-linear of full rank, and no extra simple boundary is omitted.
It also checks the nonlinear gauge-straightening theorem: every slice
`C_1=q-3C_0h` with `q!=0` and translation-invariant `h` is carried to
`C_1=q` by explicit polynomial source and target automorphisms.  It checks
the symmetric lower-unipotent family, the discriminant invariant, and the
variable-time Jacobian formula `1+D(h)`, which makes invariance necessary
for a single shear to be an automorphism.  The stress-test family
`C_1+tC_0^2=1` is verified directly to have source `A^3` and Jacobian `-1`
before being reduced to the foundational class.
The Borel corollary now exhausts every polynomial upper- or
lower-triangular `GL_2` gauge as well: its diagonal entries must be
constants, leaving exactly one classified invariant shear.
On the invariant coefficient hyperplane `C_0=0`, the same checker restricts
the time `h=4*C_1*C_3-C_2^2` to an explicit `A^3` automorphism, verifies its
inverse, determinant, and multidegree `(1,3,5)`, and checks the exact linear
conjugacy to the Nagata automorphism with parameter `-4`.  Wildness uses the
external Shestakov--Umirbaev theorem.  The first three swapped iterates are
expanded exactly and have multidegrees `(4s-3,4s-1,4s+1)`; this is a
known-family calibration, not a resolution of the open `(7,8,12)` case.
The alternating regression verifies the exact two-shear rank-two Jacobian
formula.  When the first time is invariant, conjugation reduces the second
factor to the single-shear theorem and gives an if-and-only-if transported
kernel criterion.  A Gröbner coefficient audit excludes every normalized
linear-time cancellation between two individually noninvertible factors.
The all-degree support theorem excludes every pair of nonzero monomial
times; the checker exhausts all 1,156 pairs through degree three.  Exact
graded ranks through degree eight give cokernel dimensions
`0,0,0,1,0,0,0,1`, confirming that the general recursive gauge equation has
only one discriminant obstruction in every fourth degree.  A second exact
checker parametrizes the ten-dimensional quadratic cancellation kernel and
proves that its degree-four discriminant projection vanishes identically;
all coupled basis directions admit recursive corrections through degree
eight.  The written `sl_2` divergence identity proves the bilinear
vanishing in every degree.  The ranked next attacks are recorded in
[`cancellation/CUBIC_CLOSURE_ATTACKS.md`](cancellation/CUBIC_CLOSURE_ATTACKS.md).

This runs the construction, parameter arithmetic, boundary, monodromy, and
current-ansatz rigidity regressions.  It includes the endpoint-moment
reduction of the cancellation contact resultant: the general triangular
identity is checked exactly on a bounded grid, while the complete
`r=1,2,3,4` columns are proved uniformly in `m`.  It also checks the
irreducibility transfer proving every `1<=m<=1000` column uniformly in `r` and
an explicit effective `r`-tail for each fixed `m`:

```bash
.venv/bin/python scripts/verify_parameter_irreducibility.py
.venv/bin/python scripts/verify_parameter_irreducibility_dusart_frontier.py
.venv/bin/python scripts/verify_parameter_irreducibility_sharp_dusart_frontier.py
.venv/bin/python scripts/verify_parameter_irreducibility_adaptive_dusart_frontier.py
.venv/bin/python scripts/verify_contact_resultant_irreducible_ranges.py
```

The second command is a slow four-process exact replay of the 2192 residual
pairs in `301<=m<=499`; it is kept out of the ordinary `verify-master`
target.  The third is a slower six-process exact replay of the 2899 residual
pairs in `500<=m<=741`; the fourth replays the 3335 adaptive residual pairs
in `742<=m<=1000`.  Both are likewise kept out of `verify-master`.

The `r=3` certificate checks
coefficientwise positivity of all six principal minors of the reciprocal
eliminant's Schur--Cohn matrix.  The heavier `r=4` certificate computes the
degree-eleven eliminant's `(9,2)` Schur--Cohn inertia, runs a 228-cell rational
Rouche localization, and proves the remaining argument separation by exact
angle and Bernstein-sign certificates.

The complete `r=5` column is a separate, substantially heavier exact replay.
It requires Singular for its boundary resultants:

```bash
.venv/bin/python scripts/verify_contact_resultant_r5.py
```

The formerly first open fixed-`r` column has an exact bounded-degree
reduction: the following Singular-backed checker constructs the
quintic--sextic endpoint equations and verifies that their residual eliminant
has degree 29 in `y` and degree 90 in `m`.

```bash
.venv/bin/python scripts/verify_contact_resultant_r6_reduction.py
```

The branch-at-infinity replay then proves eventual nonvanishing in that
column.  It checks the complete Newton edge after `y=1+c/m`, the squarefree
degree-29 edge polynomial, and the linear reconstruction of the limiting
`z`.  Lindemann--Weierstrass separates algebraic `z` from `exp(c)`.  This
intermediate argument does not by itself provide an explicit threshold in
`m`.

```bash
.venv/bin/python scripts/verify_contact_resultant_r6_asymptotic.py
```

For an effective certificate on the limiting edge itself, the following
atlas computes the first two `y` terms and first two `z` terms for all 29
branches (compressed to 15 complex-conjugation classes).  It uses 29 disjoint rational
Rouche disks and rational exponential enclosures to prove the strict modulus
gap `|z_0|!=|exp(c)|` branch by branch.  This limiting-edge atlas does not by
itself extract a positive tail threshold in `t=1/m`.

```bash
.venv/bin/python scripts/explore_contact_resultant_r6_branch_atlas.py
```

The effective continuation requires the pinned `python-flint` dependency in
`requirements.txt`.  It certifies 29 disjoint Rouche tubes on each of 256
rational cells covering `0<=t<=1/41`, separates the sixth-power identity by
modulus or phase on every tube, and checks the finite range `1<=m<=40`
modulo `1,000,003`:

```bash
.venv/bin/python scripts/verify_contact_resultant_r6_effective.py
```

A separate bounded structural audit constructs the limiting endpoint systems
for `r=5,6,7,8`.  At `r=7,8` it verifies squarefree branch polynomials of
degrees 42 and 55, excludes `c=0` and `z=infinity`, and reconstructs a unique
finite `z`.  The first command is the cross-column limiting-system audit.
The next two commands construct the full bidegree-`(42,126)` and
`(55,200)` endpoint eliminants, identify their complete top Newton edges,
and prove eventual nonvanishing.  These three commands do not claim an
effective threshold or a continuation uniform in `r`:

```bash
.venv/bin/python scripts/verify_contact_resultant_fixed_r_branch_schema.py
.venv/bin/python scripts/verify_contact_resultant_r7_asymptotic.py
.venv/bin/python scripts/verify_contact_resultant_r8_asymptotic.py
```

The effective fixed-column template is replayed for `r=7,8` by

```bash
.venv/bin/python scripts/verify_contact_resultant_r7_effective.py
.venv/bin/python scripts/verify_contact_resultant_r8_effective.py
```

For `r=8`, the second command certifies 55 disjoint roots on each of 1,024
rational cells covering `0<=t<=1/1001`.  Its 56,320 Arb tubes split into
52,224 modulus and 4,096 phase separations.  It then uses FLINT over
`GF(1,000,003)` for 1,000 degree-preserving endpoint gcd certificates,
closing every integer `m>=1`.  The shared checker is configuration-driven
for `r=6,7,8`; subsequent columns reuse the endpoint-chart construction,
Rouche tubes, logarithmic separation, and finite-field completion after
supplying their exact degree/edge data and a certified partition.

The additional finite `5<=r<=12` endpoint grid is quick to replay.  It checks
203 monic gcd certificates modulo `1,000,003`, including denominator and
leading-coefficient unit conditions:

```bash
.venv/bin/python scripts/verify_contact_resultant_modular_grid.py
```

It also runs the log-geometric bridge regression, including the reciprocal
determinant, canonical Jacobian-LND exponent, the degree-two plinth/Stein
countermodel, spectral squarefreeness, and Laurent-tail descent.  The reusable
classifier additionally checks exact prime valuations, both localized chart
compositions, boundary elimination, the displayed residue degree, the full
Stein field via local-slice invariantization, hidden covers, and the spectral
gcd obstruction.  It also checks the unsliced divided-difference Hensel
multiplier which upgrades the boundary value to the complete cancellation
jet and global slice.  Its built-in examples can be inspected directly:

```bash
.venv/bin/python scripts/classify_reciprocal_link.py cancellation
.venv/bin/python scripts/classify_reciprocal_link.py masuda
.venv/bin/python scripts/classify_reciprocal_link.py masuda-hidden
```

The arithmetic portion also checks the fixed-row Newton-ramification
extraction:

```bash
.venv/bin/python scripts/verify_fixed_r_newton_ramification.py
```

It verifies the reciprocal numerator and prime-power congruence on a bounded
grid and exact cyclotomic-cluster Newton edges for derivative orders one
through eight.  The analytic density estimate is the cited external theorem
input, not a finite computation.

## External quartic islands

Juntang Zhuang's pinned `F4a`, `F4b`, and `F4c` examples have an independent
compact reconstruction and canonical-boundary audit:

```bash
.venv/bin/python scripts/verify_external_quartic_islands.py
```

This command is also part of `make verify-regressions`.  It requires no
network access and does not copy or execute the upstream checker.

## Decorated normalization, affine-mark faithfulness, and Hurwitz--LL calculations

The LL critical-value incidence, low-pole filtration, contravariant
triangular target convention, affine pencil transport, higher-zero Newton
polygons, nonzero multiple-root collisions, and normalized rerooting
identities are checked exactly by

```bash
.venv/bin/python scripts/verify_stable_generator_rigidity.py
.venv/bin/python scripts/verify_generic_affine_mark_faithfulness.py
.venv/bin/python scripts/verify_intrinsic_selector_attack.py
.venv/bin/python scripts/verify_hasse_typical_seed_recovery.py
python3 scripts/verify_positive_characteristic_deformation_landscape.py
.venv/bin/python scripts/verify_multicluster_ll_comparison.py
.venv/bin/python scripts/verify_labelled_node_saturation.py
.venv/bin/python scripts/verify_branch_wonderful_pullback.py
.venv/bin/python scripts/verify_source_vertex_rigidity.py
.venv/bin/python scripts/verify_general_radial_source_atlas.py
.venv/bin/python scripts/verify_polynomial_monodromy_forests.py
.venv/bin/python scripts/verify_monodromy_inertia_characters.py
.venv/bin/python scripts/verify_recursive_resonance_atlas.py
.venv/bin/python scripts/verify_h1_h2_comparison_obstruction.py
.venv/bin/python scripts/verify_branch_scale_fan.py
.venv/bin/python scripts/verify_degree_six_branch_target_graph.py
.venv/bin/python scripts/verify_degree_six_admissible_equal_scale.py
.venv/bin/python scripts/verify_degree_six_admissible_radial_atlas.py
.venv/bin/python scripts/verify_degree_six_admissible_maxwell_atlas.py
.venv/bin/python scripts/verify_degree_six_central_hurwitz_selection.py
.venv/bin/python scripts/verify_degree_six_stack_inertia.py
.venv/bin/python scripts/verify_degree_six_stacky_fan_descent.py
.venv/bin/python scripts/verify_rerooting_groupoid_boundary.py
.venv/bin/python scripts/verify_coarse_affine_mark_descent.py
.venv/bin/python scripts/verify_restricted_ll_degree.py
.venv/bin/python scripts/verify_caustic_maxwell_boundary.py
```

For at least eight labelled finite branch values, the target receiver in the
wonderful-pullback construction is `Mbar_0,b+2` with `b+2>=10`.  Its Cox ring
is not finitely generated by the externally cited non-polyhedral
pseudo-effective-cone theorem.  This corollary has no local checker and does
not assert the same statement for the normalized pullback graph.

These checks support the generic affine-mark faithfulness theorem: the coarse
fiber is the exact rerooting orbit and every nontrivial rerooting moves the
unique unramified affine sheet into the reconstruction boundary.  The
Hasse-typical checker separately proves the sharp
`floor(log_p(N))+1`-channel coefficient repair in positive characteristic
and replays the degree-eight `F_5` collision.  It also verifies a clean
five-member degree-twelve `F_5` family whose distinct seeds have identical ordinary
derivatives and therefore define the identical weighted polynomial map,
proving that the channels cannot be made map-intrinsic without enriching
the construction.  The same checker then constructs the dimension-preserving
correction `A -> A-K(W)/(cC^2)`, verifies polynomiality, determinant one, and
all five intended inverse pencils, and certifies the marked transverse node
in the `c=2` member.  Five distinct reduced equal-image Groebner bases then
prove that the enriched maps are pairwise stably left--right inequivalent.
The written theorem upgrades this example to every odd-characteristic tame
clean degree by recovering the complete primitive-root factor from the
intrinsic second-boundary edge data.  It also records the characteristic-two
parity reconstruction, the identically singular old suspension parameter,
and the scalar-ansatz no-go theorem.  The checker then verifies the
weight-redistributed replacement
`u=1+x^2y`, `gamma=1+xz`: its coordinates are polynomial, its Jacobian is
one, and its inverse pencil is the prescribed normalized seed.  The final
characteristic-two block verifies the radicial discriminant factor
`W^2-T` and explicit squarefree, compressed-birational wild-clean witnesses
in every degree from five through sixteen; the formulas prove the resulting
stable-faithfulness theorem uniformly for all `N>=5`.  It also verifies the
complete symbolic quartic slice
`(1+lambda)W^2+W^3+lambda W^4`; the two affine marks on the radicial edge
remove the former low-support ambiguity, while the normalized cubic is
unique.  The same weight redistribution is checked in characteristics
three, five, and seven on the universal singular-parameter quartic
`2W^2-3W^3+W^4`, confirming that it complements the original chart exactly
on `2+H''(1)=0`.  The full-edge theorem then removes the old Hessian
degree/support restriction on every boundary-clean generically birational
locus.  Finally `d(WH'-H)=W dH'` proves that critical birationality is
automatic for every odd-characteristic exact-double seed; the
characteristic-two checker exhaustively regresses the corresponding
clean-implies-compressed-birational lemma through degree twelve.
Repeated-root examples in characteristics two, three, and five then verify
that the normalized second-boundary prime retains the complete primitive
root divisor with multiplicities, even when its critical image collides
with the zero cluster.  This is the executable collision regression for the
full theorem: the smaller marked-edge quotient
`(A^1_W;(W),(W-1),div(H/(W^2(W-1))))` reconstructs the normalized seed
exactly on every declared stratum.
The
selected root extends on the marked corrected graph, and the
normalized-Stein, completed-chart, and conductor comparisons are complete at
arbitrary simultaneous collisions.  Coarse affine-mark descent is also
complete over that graph: the marked invariant ring is the universal
monic-root incidence,
and the total-collision fiber `k[T]/(T^mu)` has one geometric point.  The
specialized restricted-LL
degree and caustic/Maxwell boundary-class calculations have no recorded
external review.

The companion affine-stratum audit verifies that the root-one component is
regular and that a nontrivial rerooting sends it to an extra-root boundary
component.  The multicluster audit checks distinct tangent lines, all pairwise
intersection numbers, the conductor exponent
`e_i(sum_j e_j-1)`, and regularity of the full marked-root incidence at
collisions.  The H1/H2 obstruction checker recovers the degree-five
`(x^3,y^2)` normalized blowup.  The branch-scale checker then computes the
degree-six `(2,2,2)` moving critical values, all six weighted braid-fan
chambers, and a triple-resonance cross-ratio proving that the radial fan is
only the first layer of the full logarithmic graph.  The wonderful-pullback
checker enumerates the complete `Mbar_0,n` boundary building set and maximal
nested sets for four through seven target marks, verifies permutation
equivariance, and recovers both the degree-five weighted blowup and the
degree-six six-line/four-center target from that one construction.  The
source-vertex checker exhausts 2,024 zero/pole divisor profiles in degrees
one through seven and proves computationally that two fibers reconstruct a
rational component map up to scale while one third-flag point fixes the
scale.  The general radial-source
checker then verifies the connector/local-polynomial-tail/identity-strand
rule for 780 multiplicity profiles and 48,580 ordered scale types, including
all component degrees, Riemann--Hurwitz identities, node partitions, lcm
saturations, label permutations, and independent dynamic verification of
the full-chain inertia formula.  It finds 42,158 nontrivial
unequal-multiplicity types in this range; equal multiplicities remain
trivial.  The monodromy-forest checker then
exhausts all 1,441 reduced polynomial transposition factorizations through
degree six and proves that every nested resonance source tree and node
partition is the corresponding edge subforest; it recovers pairwise
Maxwell, triple Maxwell, and caustic nodes from one rule.  The recursive
resonance-atlas checker then verifies framed residue coordinates on all 534
nested families with two through five branch labels, all 534 affine gauge
changes, 1,453 one-step and
2,926 two-step contractions, normalized flag equations through degree seven,
84 source/target frame transitions, 63 nonfactorized smoothing families,
automatic tame character extraction, 76 bounded full-centralizer radial
charts, all 89 degree-six
interval-nested families, and the order-four pair--triple inertia.  This
closes the former
explicit-stack gap.  The finite-normalization theorem
uses finiteness of the fully marked admissible-cover branch morphism to prove
that the normal wonderful graph is already the complete coarse source graph;
no additional source-side coarse blowup is possible, and corrected H2/H3
are unconditional.  The monodromy-centralizer checker computes all
polynomial tree deck groups through degree six, all cyclic connector groups
through degree eight, and anchored/unanchored inertia on every collision
node.  The recursive checker corrects the full-chain radial calculation:
equal multiplicities have trivial inertia, while an ordered partition
\(B_0|\cdots|B_k\) of arbitrary multiplicities has order
\(\prod_jL_j/M_j\); it checks 76 bounded equal and unequal charts.  The
centralizer checker gives one generic formula covering Maxwell and caustic
resonance.  The
complete-target
checker identifies the radial target with the three-coordinate-point blowup
of `P^2`, the stable target with the additional diagonal-point blowup
`Mbar_0,5`, and its source pullback with four reduced triple-Maxwell
branches.  The equal-scale admissible checker constructs the central
degree-six component and three quadratic tails, verifies all
Riemann--Hurwitz counts, and proves that the three index-two source nodes
normalize into exactly the same four Kummer branches.  The radial-atlas
checker then enumerates all thirteen ordered scale types, verifies degree six
and Riemann--Hurwitz on every target-bubble preimage, and checks every
node-index partition and Kummer saturation count.  The Maxwell-atlas checker
handles all three pairwise collision divisors and the triple collision,
matches their two- and four-branch source-node normalizations, and proves
that their residual radial intersections are transverse while their
coordinate intersections are already radial equality faces.  The central
Hurwitz-selection checker finds two ambient degree-six cover classes with
the required profiles, then proves by an exact square-cubic branch invariant
that the labelled source-root cross-ratio selects the polynomial class as a
reduced local branch.  The stack-inertia checker separates normalization
branches from genuine label-preserving cover inertia: every radial lift in
the equal-multiplicity degree-six chart has trivial inertia, while pairwise
and triple Maxwell lifts each retain one diagonal `mu_2`.  The stacky-fan
checker constructs the four-divisor Maxwell
root complex, proves all pair--triple face inclusions and `S_3` equivariance,
computes the four radial quotient orbit types, verifies the pair--triple and
radial--Maxwell inertia ranks needed for smooth tame-stack reconstruction,
and keeps the local
`(S_2)^3 semidirect S_3` pair-root stabilizer separate.  The
general labelled-node checker exhausts 1,554 index profiles, proves the
phase-quotient and label-preserving inertia formulas, checks permutation
equivariance, and verifies that the corrected marked/unmarked quotient over
any labelled normalized graph has degree `N-2`.  This makes label gluing and
the finite H2 factor formal over the `H1-COARSE` graph, independently of the
substantially stronger `H1-STACK` theorem.  The
rerooting-groupoid
audit separately checks the quotient degree `N-2`,
the selected-in/selected-out boundary pullbacks, generic transposition
ramification after coefficient contraction, and the distinction between a
cyclic total-collision slice and generic divisor inertia.  These three audits
and the companion affine-stratum audit are part of `make verify-regressions`.
The restricted-LL audit checks the Cayley/marking count and independently
computes degrees `8` and `75` from the quartic and quintic critical-value
eliminants.  The caustic--Maxwell audit checks the unique invariant Keel
relation, every collision and infinity valuation, both boundary
presentations, and the exact factorization `LL-discriminant=C^3 M^2` in
degrees four and five.  All displayed commands are part of
`make verify-regressions`.

## External consequence identities

Christopher D. Long's direct Gaussian-moment, `(xz)`, `SU(2)`, and `SO(3)`
identities, together with the exact normalization of the foundational map
used in his BCW discussion, have a dedicated target:

```bash
make verify-external-consequences
```

The Gaussian, `(xz)`, spherical `SO(3)`, and algebraic Haar scripts use only
the Python standard library.  Their bounded exact regressions are
distinguished from the all-exponent proofs in the canonical notes.  The
`SO(3)` replay checks the displayed moments through order fifteen and the
endpoint-jet identity through order one hundred.  The two algebraic checkers
verify the unique normalized functional on `UV+T^2=1`, its three
infinitesimal `so3` identities, the factorial functional on
`k[SL2]`, all six left/right `sl2` identities, and the explicit
`SL2/T` pullback.  The proof and the quotient/transfer theorem are in
[`ALGEBRAIC_HAAR_QUADRIC_AND_SL2.md`](extended-geometry/ALGEBRAIC_HAAR_QUADRIC_AND_SL2.md).
A separate symbolic checker proves the `SU(2)=S^3` Haar density in Hopf
coordinates, retaining an independent compact integration proof.  The same target
also performs all 18 balanced BCW steps and checks the resulting 79-variable
cubic-homogeneous collision, writes its sparse artifact, and replays it with a
separate standard-library implementation.  It then runs the shared-factor
optimization, which introduces 13 variables, reaches degree three in
dimension 16, and writes and replays a 33-variable baseline artifact.  It then
computes the exact rational rank 7 of the cubic component vector, constructs
the rank-compressed 24-variable cubic collision, and independently replays
the factorization, sparse map, and collision using only the standard library.
It then removes the two-dimensional constant Jacobian kernel, constructs the
22-variable quotient, and independently replays `BK=0`, `BC=I`, `H=HCB`,
cubic homogeneity, the descended collision, and the triangular determinant
factorization using only the standard library.
Finally, the essential-dimension search freezes a different 17-dimensional
trace of cubic-output rank six, homogenizes it in 24 variables, removes its
three-dimensional constant kernel, and independently replays the resulting
21-variable collision from the original map using only the standard library.
The backward-cubic continuation keeps the nonhomogeneous and homogeneous
dimension objectives separate.  It audits MacFarlane's displayed `F13` and
`G20`, restricts the sole fixed covector `tau` at the collision level, and
verifies the exact stable factorization
`M19=A_B o (F13 x I_6) o S_gamma`:

```bash
make verify-backward-cubic-reduction
```

The generated records are
[`macfarlane_g20_dimension_reduction_audit.json`](artifacts/generated-results/macfarlane_g20_dimension_reduction_audit.json),
[`macfarlane_f13_low_degree_invariants.json`](artifacts/generated-results/macfarlane_f13_low_degree_invariants.json),
and
[`backward_cubic_reduction_calibration.json`](artifacts/generated-results/backward_cubic_reduction_calibration.json).
The same target applies the two backward objectives and the pair-aware
collision policy to the retained restricted-minima archives, reconstructs
two current representatives exactly, and writes
[`backward_cubic_current_applications.json`](artifacts/generated-results/backward_cubic_current_applications.json).
The generic calibration also proves that the parent is isotrivial over
`t!=0`, its `t=0` fiber is triangular and injective, and every parent
collision can therefore be normalized to `t=1`; an exact `t=2` MacFarlane
collision is replayed as a regression.
It also runs the established `16 -> 24` rank-compressed BCW route with the
new reverse-companion regression enabled.
With the pinned external determinant certificate, the same audit updates the
external-certificate frontiers to `n_cub<=20` and, by homogeneous cotangent
lift, `n_HN,4<=40`; the internal dependency-free replay endpoints remain 21
and 42.
These commands calibrate the backward compiler and close stated direct
linear/degree-at-most-three routes; they do not construct a twelve-variable
map.

The next coordinate-pair reduction goes beyond pullback-fixed invariants.
It uses `s=F13_13=x13+x2^2` as a source coordinate and the target square
completion `y4-y8^2`.  The resulting exact relative form restricts at
`s=0` to a 12-variable degree-three Keller collision.  A direct sparse
determinant expansion and a separate standard-library implementation replay
the theorem.  Its cubic-output rank is six, giving a 19-variable
cubic-homogeneous parent and the updated bounds `n_cub<=19` and
`n_HN,4<=38`:

```bash
make verify-macfarlane-f12
```

The generated record is
[`macfarlane_f12_coordinate_pair_reduction.json`](artifacts/generated-results/macfarlane_f12_coordinate_pair_reduction.json).

The first exact continuation toward eleven variables classifies every linear
target coordinate whose pullback is a polynomial graph coordinate.  The raw
degree-three coefficient ideal is the unit ideal in all nine possible pivot
families.  At the literal triangular coordinates, every graph deletion has
degree four or five and at least one high-degree defect lies outside the
complete degree-at-most-three target-shear span in the other raw retained
outputs.  The two closest literal cases remain outside that span through
target degree four:

```bash
make verify-k12-coordinate-pair-frontier
```

The generated record is
[`k12_coordinate_pair_frontier.json`](artifacts/generated-results/k12_coordinate_pair_frontier.json).
This is a bounded obstruction, not a dimension-eleven lower bound; nonlinear
source coordinates and ordered multi-stage target automorphisms remain open.

The parameterized continuation then retains every linear target coordinate
whose pullback has a quadratic graph. Fixed full-column and augmented
minors, together with unit-ideal covers in the graph parameters, exclude
quadratic target completion for all six pivot families and cubic target
completion for all five single-defect families:

```bash
make verify-k12-parameterized-completion
```

The generated record is
[`k12_parameterized_completion_frontier.json`](artifacts/generated-results/k12_parameterized_completion_frontier.json).
This remains a bounded theorem. The multi-defect `z8` cubic completion is
handled by the next command; cubic graph corrections and ordered target
stages remain outside the combined scope.

The remaining multi-defect `z8` cubic system is assembled without a full
fraction-field expansion. Sparse modular elimination selects three minors,
which are then reconstructed exactly over the rational parameter ring.
Their determinant opens generate the unit ideal and every augmented
determinant is `9/7` times its column determinant:

```bash
make verify-k12-z8-cubic-completion
```

The generated record is
[`k12_z8_cubic_completion_frontier.json`](artifacts/generated-results/k12_z8_cubic_completion_frontier.json).
Together with the preceding command, this excludes one-stage cubic target
completion for all six quadratic graph-coordinate families.

The same sparse compiler extends through target degree four on every
single-defect family. Each parameter family has a nonzero constant
`990 x 990` column minor and a nonzero constant augmented minor:

```bash
make verify-k12-single-defect-quartic-completion
```

The generated record is
[`k12_single_defect_quartic_completion_frontier.json`](artifacts/generated-results/k12_single_defect_quartic_completion_frontier.json).
Only the much larger multi-defect `z8` quartic family remains in this graph
class.

The cross-construction audit compares the public dimension-38 route with the
independent `K12` route.  It checks the shared compressed cost `n+r=18`,
obstructs all seven source-affine linear/quadratic pivot completions of the
public eleven-variable lift.  On `K12`, it also obstructs quadratic target
completion of the nonlinear `z8` pivot and finds a fourteen-parameter
family of coordinated degree-preserving quadratic source shears.  A fixed
minor and an inconsistent exact Schur system of ranks `(5,6)` prove that
no member of that family lowers the cubic-output rank from six to five:

```bash
make verify-hvc38-cross-frontier
```

The generated record is
[`hvc38_cross_construction_frontier.json`](artifacts/generated-results/hvc38_cross_construction_frontier.json).
This is a bounded frontier computation, not a lower bound for quartic HVC.

The next gap-closure audit uses the square identities at the public `d`
pivot and local `z8` pivot to reduce nonlinear completion to a filtered
pullback calculation.  Good-prime ranks and matching exact kernels exclude
both pivots through target degree eight.  It then combines 140 quadratic
source columns with 792 elementary quadratic target columns.  The
36-dimensional exact high-degree kernel contains seventeen directions that
integrate to genuine triangular one-parameter source-target families.  On
their combined seventeen-parameter degree-three locus, the ideal obtained by
adjoining a selected cubic rank-six minor has Gröbner basis `[1]`:

```bash
make verify-hvc38-gap-closure
```

The generated record is
[`hvc38_gap_closure.json`](artifacts/generated-results/hvc38_gap_closure.json).
This excludes only the stated bounded pivot algebras and quadratic
source-target family; it does not prove minimality at dimension 38.

The maximal-block continuation enumerates all six maximal jointly affine
source blocks of `K12`.  For each block it combines every complementary
quadratic source shear with all 792 elementary quadratic target directions,
lifts the complete good-prime high-degree kernel over `QQ`, and verifies a
linearized rank-six Schur witness.  It then integrates every kernel
direction into one full triangular source-target family—including
source-only directions and directions that fail to preserve degree three
individually—and asks whether the exact degree-three locus can have
cubic-output rank at most five.  Pinned packets of at most 32 cubic minors
give unit ideals in Singular:

```bash
make verify-hvc38-maximal-block-closure
```

The generated record is
[`hvc38_maximal_block_closure.json`](artifacts/generated-results/hvc38_maximal_block_closure.json).
This closes the full quadratic left-right kernel class on all maximal
jointly affine blocks.  It remains a bounded theorem, not a dimension-38
minimality result.

The tensor continuation computes both natural coefficient flattenings of
the `K12` quadratic and cubic tensors and of the cubic-homogeneous `G19`
tensor.  Exact rational row reduction gives input-directional ranks `12`,
`12`, and `19`; hence all three common right kernels are zero.  The `G19`
output rank is `18`, with sole left annihilator the fixed `tau` output.  The
same checker replays both collisions and the companion scaling identity:

```bash
make verify-k12-tensor-module-frontier
```

The exact generated record is
[`k12_tensor_module_frontier.json`](artifacts/generated-results/k12_tensor_module_frontier.json).
This excludes constant linear tensor quotients of the displayed maps and
gives pure-cube decomposition lower bounds `12` and `19`; it does not
exclude nonlinear graphs, nonconstant modules, Schur elimination, or a
different tensor.

A separate finite-field scout enters the larger linear-coordinate families
whose graph corrections are genuinely cubic.  It searches all parameter
supports of size at most two with values in `{-2,-1,1,2}`, plus 250
deterministic random points per parameter count, over both `GF(101)` and
`GF(103)`.  Every bad retained output is tested against all 10 linear and 55
bilinear target monomials in the other raw outputs:

```bash
.venv/bin/python scripts/search_k12_cubic_graph_bilinear_completions.py \
  --support-max 2 --values=-2,-1,1,2 --random-samples 250 \
  --random-seed 20260804 --primes 101,103 --keep-closest 5 \
  --output artifacts/generated-results/k12_cubic_graph_bilinear_modular_search.json
```

The 15,688 modular evaluations contain no survivor at either prime.  The
generated
[`k12_cubic_graph_bilinear_modular_search.json`](artifacts/generated-results/k12_cubic_graph_bilinear_modular_search.json)
is a bounded discovery experiment, not a rational obstruction or a
dimension-eleven lower bound.

The exact continuation lifts the modular echelon rows across all nine
complete cubic graph-coordinate parameter spaces.  It combines constant
full-column/augmented minors with rank-stratified determinant covers over
`QQ`:

```bash
make verify-k12-cubic-graph-bilinear-obstruction
```

The generated record is
[`k12_cubic_graph_bilinear_obstruction.json`](artifacts/generated-results/k12_cubic_graph_bilinear_obstruction.json).
Five families have constant determinant ratios `-3,-9,-1,-1/2,-1/2`.
For `z5,z6,z7,z8`, the checker verifies exact determinant-open covers and
the residual closed strata, using explicit column relations for `z5,z6,z7`
and the complete BCR5 quadratic-family obstruction for the closed `z8`
stratum.  Together these certificates close all nine normalized linear graph
families through one-stage target degree two.

The cubic-target continuation first quotients by the certified bilinear
column space and adds the 220 cubic target monomials.  Its bounded two-prime
discovery pass uses all parameter supports of size at most one with values
in `{-2,-1,1,2}` plus 250 deterministic random points per parameter count:

```bash
make search-k12-schur-cubic-completions
```

All 2,902 modular evaluations have augmented-rank increment one and there is
no survivor.  The exact continuation is:

```bash
make verify-k12-graph-cubic-completion-obstruction
```

It reconstructs the complete 285-element target basis over the rational
parameter rings.  The pivots `z4,z8,z9,z10,z11,z12` have constant
full-column and augmented minors.  Determinant-open covers plus exact
closed-stratum relations handle `z5,z6,z7`.  Hence all nine normalized
linear graph families are obstructed through one-stage target degree three.
The generated records are
[`k12_schur_cubic_completion_modular_search.json`](artifacts/generated-results/k12_schur_cubic_completion_modular_search.json)
and
[`k12_graph_cubic_completion_obstruction.json`](artifacts/generated-results/k12_graph_cubic_completion_obstruction.json).

The restricted-minima continuation then changes the BCW circuit before
homogenization: it exposes two polynomial gates, cancels one complete
multi-term circuit block, and scores every partial trace by its Jacobian
power-rank profile.  The two frozen winners are a 22-variable
cubic-homogeneous collision of exact index 18 and a 24-variable collision of
exact generic rank 17 and index 18.  Singular certifies their generic kernel
dimensions, while independent standard-library audits multiply the full
polynomial Jacobians and verify `(JH)^17!=0`, `(JH)^18=0`.  A fifth grouped
atom cancels the first-coordinate circuit `x^2(3y+xz)`.  The expanded
32-family Pareto search finds a separate 22-variable cubic source whose
44-variable HN lift has exact generic Hessian rank 37:

These computations form the now-frozen upper-bound track.  The commands below
reproduce the recorded search; they are not an active broad-search queue.
If the program is reopened, it should begin with a theorem-directed
five-dimensional cubic classification or the invertibility-only question for
arbitrary cubic-homogeneous Keller maps with `(JH)^3=0`.

```bash
.venv/bin/python scripts/search_restricted_bcw_circuits.py \
  --width 64 --max-steps 24 --prebeam-factor 2 --partial-power-depth 8 \
  --skip-terminal-hessian-power \
  --enable-atom x2s --enable-atom v2r --enable-atom qb \
  --enable-atom v2h --enable-atom y2vb \
  --output artifacts/generated-results/restricted_bcw_circuit_search_v2_w64.json
.venv/bin/python scripts/search_rank37_gate_perturbations.py \
  --width 16 --max-steps 10 --prebeam-factor 3 \
  --partial-power-depth 8 \
  --output artifacts/generated-results/rank37_gate_perturbation_search.json
.venv/bin/python scripts/verify_index_reduced_bcw_22_route.py
python3 scripts/audit_index_reduced_bcw_22_independent.py
.venv/bin/python scripts/verify_rank_reduced_bcw_24_route.py
python3 scripts/audit_rank_reduced_bcw_24_independent.py
.venv/bin/python scripts/verify_hessian_rank_reduced_bcw_22_route.py
python3 scripts/audit_hessian_rank_reduced_bcw_22_independent.py
.venv/bin/python scripts/verify_hessian_rank_35_identity_slice.py
python3 scripts/audit_hessian_rank_35_identity_slice_independent.py
.venv/bin/python scripts/search_identity_slice_hessian_rank.py
.venv/bin/python scripts/search_identity_slice_local_perturbations.py
.venv/bin/python scripts/verify_hessian_rank_34_double_identity_slice.py
python3 scripts/audit_hessian_rank_34_double_identity_slice_independent.py
.venv/bin/python scripts/verify_index_three_inverse_model.py
.venv/bin/python scripts/verify_index_three_degree_bound_counterexample.py
.venv/bin/python scripts/derive_index_three_tree_obstruction.py
.venv/bin/python scripts/verify_restricted_minima_frontier.py
```

The second index-three command replays van den Essen's dimension-five
generic-rank-three automorphism, proves `(JH)^3=0`, verifies both inverse
compositions, and extracts the nonzero degree-eleven and degree-thirteen
terms.  The tree command independently evaluates the degree-eleven normal
form on the same tensor.  Together they disprove the proposed uniform
inverse-degree-nine bound while leaving the full-class invertibility-only
question open.

The resolved two-real theorem, the first classified three-real minimality
island, and the cross-conjecture minimum ledger have their own fast exact
target:

```bash
make verify-counterexample-scoreboard
```

The generated
`artifacts/generated-results/minimal_counterexample_scoreboard.json` records
the exact unrestricted-GVC failure dimension as three.  Its current
whole-file SHA-256 is
`09df0e398def5df799243c906066f0b469b17ccf63f7d9261e8944a96fe8f8b1`.

This proves GMC for every quadratic Gaussian polynomial in every dimension,
checks the two-weight and affine-circular-source obstructions in two real
variables, exactly excludes all 27 mixed-sign cubic three-weight supports on
their 72 nonvanishing charts using moments through order eight, excludes 29
of the 33 mixed-sign cubic four-weight supports on 97 charts using moments
through order six, and then excludes all four charts of the symmetric
exceptional support by three good-prime quotient-algebra certificates: in
each representative the tenth moment acts with rank 84 on the
84-dimensional order-eight quotient, and circular-coordinate reflection
supplies the fourth chart.  Seven further exact rational unit-ideal
calculations exclude the last three supports and 20 charts through moment
six.  Thus all 121 mixed-sign four-weight cubic charts are closed and a
cubic GMC(2) counterexample needs at least five rotational weights.  The
target also proves the Bessel--factorial moment formula for the three-level
family with support `{-1,0,1}` and computes 31 exact rational unit ideals:
6 charts in degree four through moment six, 10 charts in degree five through
moment eight, and 15 charts in degree six through moment nine.  A
prime-endpoint theorem now closes that family in every degree: at odd
prime \(p\), the orders \(p\) and \(2p\) isolate the \(C^p\) and \(D^p\)
endpoints according to the two possible \(U\)-adic order inequalities.
The bounded charts remain finite-cutoff regressions.  A companion exact
arithmetic check verifies the prime coefficient and factorial congruences;
an independent pure-Python audit reconstructs both polynomial endpoint
identities and both normalized factorial cases.  The unit-star regression
then checks the primitive invariants, the \(p,2p,3p\) endpoints, and all
three normalized-order cases for the smallest star; the theorem covers
every support `{0,1,-d_1,...,-d_q}` and its reflection.  Another regression derives
the finite radial-moment recurrence, constructs the four-dimensional
resolvent differential system for a centered degree-\((2,3)\) pair, and
checks it against the factorial series.  The same target now verifies the
all-degree first-cycle theorem for support `{-2,-1,1,2}`: it enumerates the
toric invariant moments, checks
`CT(P^(kp)) = CT(P^k)^p (mod p)`, and eliminates every unique, adjacent-tie,
four-way-tie, and boundary valuation face using invariant moments through
degree twelve.  The target also recomputes three
unit Groebner bases for the direct Long-style
collapse, checks the Dvorsky--Long five-variable GVC and five-pair SIC
identities by dependency-free exact sparse arithmetic, and writes the
dimension/rank/index/degree scoreboard.  It does not
use bounded support enumeration to settle GMC(2): the accompanying
lower-face theorem handles arbitrary rotational support.  Its audit checks
supporting-line minima, Frobenius constant-term dilation, and normalized
factorial isolation on representative stars, mixed semigroups, and cycles.
The same target now verifies the normalized rank-one three-real ansatz:
the first three moments cut out Long's family scheme-theoretically, a formal
square identity proves all-order vanishing, and deletion saturation proves
five-term and degree-four minimality inside that ansatz.  Global minimality
outside the ansatz remains open.

The credited factorially weighted multitorus regression is:

```bash
python3 scripts/audit_factorially_weighted_multitorus.py
```

It checks a rank-two exposed coefficient, the normalized congruence at two
prime dilations, strict torus separation with a mixed cutoff, and the
circular Gaussian embedding.  These finite exact identities are regressions,
not the all-order proof.  Long's theorem, the local proof audit, the
prime-separating arbitrary-torus synthesis, and the one-radial search sieve
are in
[`FACTORIALLY_WEIGHTED_MULTITORUS_THEOREM.md`](extended-geometry/FACTORIALLY_WEIGHTED_MULTITORUS_THEOREM.md).
The checked-in Lean development still formalizes only the rank-one Gaussian
specialization.

The companion one-profile Hopf classification is checked directly by

```bash
.venv/bin/python scripts/verify_hopf_lift_classification.py
.venv/bin/python scripts/verify_hopf_lift_classification.py --require-singular
```

This exact regression expands the phase-integrated polynomial for
several windings and endpoint multiplicities, checks the full lower-jet
binomial ladder through order twenty, and includes a non-power endpoint
profile.  It also checks the quadratic profile
`R(z)=(1-z)(z-1/5)`, whose first adjacent detector vanishes, showing that
the mixed-moment nonvanishing hypothesis is essential.  The all-order result
and exact polynomiality criterion are proved in
`extended-geometry/HOPF_LIFT_CLASSIFICATION.md`; the command performs no
search in a general `V_d`.  Finally it verifies the three triangular pure
jets in the complete class
`p=x^(-1)(C(x)+D(x)t^2)`, `deg(C)<=1`, `deg(D)<=3`; they force
`D=-C^3`, after normalization, and hence reduce the all-order statement to
the endpoint theorem.  In the next complete rectangle
`deg(C)<=2`, `deg(D)<=4`, it verifies the four triangular coefficient
solutions, the fifth-jet exceptional branch `e=-a^2/12`, and the nonzero
sixth residual `1280*a^6/81081` that removes that branch.
For `deg(C)<=3`, `deg(D)<=5`, it verifies the five triangular solutions,
the residual jets `P_6,P_7,P_8`, their two exact quadratic resultants in
the final parameter, and the Euclidean certificate `gcd(Q_5,Q_6)=1`.
This forces the two remaining weighted parameters to vanish and proves
eight-jet uniqueness in that rectangle.
The same checker verifies through order twenty the uniform triangular
coefficient
`2^(m-1)*m!/(2m+1)!!` of the new `b_m` term in the `m`-th pure jet.
The written proof uses this to reconstruct `D` successively from `C` in
every fixed numerator degree and identifies the surviving eventual tangent
directions that obstruct a first-order uniform proof.
For `deg(C)<=4`, `deg(D)<=6`, it derives the first six triangular solutions
and uses exact SymPy arithmetic over `QQ` to verify both containments in the
displayed residual Groebner ideal.  That ideal is supported only at Long's
point but has the predicted two-dimensional tangent space.
For `deg(C)<=5`, `deg(D)<=7`, it derives seven triangular solutions and five
residual jets, verifies zero-dimensionality, performs exact FGLM conversion,
and checks both containments in the displayed lexicographic ideal.  The
support is again Long's point and the tangent dimension is three.
For `deg(C)<=6`, `deg(D)<=8`, the default command reconstructs all eight
coefficients of `D`, constructs the six exact residuals through jet
fourteen, checks their term counts `[19,28,37,51,64,83]`, and verifies
residual Jacobian rank one.  The second command additionally requires
Singular and computes the exact rational modular standard basis and its
FGLM lexicographic conversion.  It checks quotient length 32 and the
triangular support relations `E^8`, `F^5`, then nonzero pure powers of
`G`, `H`, and `I` after successive substitution.  Singular's modular
routine does not give a deterministic ideal-equality guarantee for this
nonhomogeneous input, so this remains computational evidence for
fourteen-jet uniqueness, not a proof.  The same command separately performs
two deterministic rational checks: the boundary
`84*E+54*F+5=0` gives the unit ideal, and exact lift matrices put
`E^6`, `48*F^3+7*E^5`, and
`972*G^2+864*F^2*E+29*E^5-108*E^4` in the residual specialization
`H=I=0`.  The remaining compact certificate is `H^3,I^3` in the full
residual ideal; equivalently, one may eliminate `I` on the certified
principal open and exclude the resulting four-variable `H!=0` chart.  The
command also computes the exact original `H=0` slice: its quotient has
length 17 and an eight-element lex basis successively forcing
`E=F=G=I=0`.  Therefore the sole possible extra locus is the principal
chart `H!=0`, equivalently the unresolved saturation `J_6:H^infinity`.
Finite-field unit lifts on that chart must first be normalized modulo the
syzygy module before CRT: the raw coefficients are prime-dependent.  At
primes 32003 and 32009 the normalized seven-multiplier support is identical,
with term counts `[778,834,814,732,646,487,692]`; ten normalized small
primes give a 150-bit modulus but do not yet suffice for balanced rational
reconstruction.  This longer reconstruction is not part of the command
above and should use a resumable checkpoint before being promoted to a
reproducer.

A final group of checks uses the first collision-coordinate values `0,1,-1`
to fix the multiplier, expands the homogeneous 42-variable quartic, descends
the contraction to `SIC(20)`, independently reconstructs the 628-term
40-variable Laplacian witness, and verifies an all-order inverse recurrence:

```bash
.venv/bin/python scripts/generate_image_vanishing_counterexamples.py
.venv/bin/python scripts/generate_identity_slice_counterexamples.py
python3 scripts/audit_identity_slice_counterexamples_independent.py
.venv/bin/python scripts/verify_inverse_coordinate_recurrence.py
```

The provenance-preserving compression audit is separate from the absolute
two-pair construction:

```bash
python3 scripts/audit_bcw_21_low_degree_invariants.py
python3 scripts/audit_bcw_21_sextic_defect_sectors.py
Singular -q scripts/audit_bcw_21_vertical_ideal.sing
python3 scripts/audit_bcw_21_septic_component_screen.py
.venv/bin/python scripts/audit_keller_near_invariant_backtrace.py
python3 scripts/audit_keller_observable_quotients.py
python3 scripts/audit_keller_provenance_compression.py
```

The first command gives characteristic-zero rank and Lie-image certificates
through degree five for the stored 21-variable map. It also records the near-invariant
`Q=X_18*X_20-X_6*X_8` and its one-term pullback defect. The second uses two
exact torus gradings to exclude both sextic correction channels and classify
all 220 sextic Lie sectors. Unique-row peeling completely certifies 25 of
the 28 dense sectors and reduces the remaining three to small exact cores.
The full sextic Lie kernel is generated by `X_20^6`, `X_20^4*Q`,
`X_20^2*Q^2`, and `Q^3`; their pullback defects extend the fixed-space
identity to `Q[X_20]_{<=6}`. It also verifies that the two degree-seven
correction sectors are exactly `X_20` times the already excluded sextic
sectors. For the remaining pure septic problem, reduction modulo `X_20`
has 657800 columns in 204 sectors; exact unique-row peeling removes 451891
and leaves 205909 columns in 79 sectors. This last statement is a support
reduction, not a classification of the residual kernel or its lifts. The
checker also verifies that the reduced derivation is constant vertical over
the fourteen-variable base and that `X_9^7` has a nonzero first lifting
obstruction because its required sextic correction sector is empty.
The Singular command proves that the six vertical coefficients generate a
height-two ideal with five displayed minimal primes, verifies their
intersection as the radical, computes
`dim_Q(B/(A_14,...,A_19))_8=158412`, and independently confirms that
`X_0^2*X_9^6` survives in the obstruction quotient. It requires Singular
with `primdec.lib`.
The component-screen command evaluates all 77520 degree-seven base
monomials on the five minimal components. It proves that 71588 support-one
base septics have nonzero first lifting obstruction, leaving 5932 monomials
for sectorwise cancellation and higher-order analysis. Exactly eight
monomials, `X_3^(7-j)*X_5^j`, have zero first obstruction.
Stacking all five restrictions in the 29 bidegree sectors gives modular rank
61060 on the full 77520-dimensional base-septic space. Hence over `Q` the
radical-level survivor space has dimension at most 16460. Embedded torsion
and higher `X_20`-adic lifts remain open.
The backtrace command reconstructs the frozen
17-step circuit and identifies `Q=c_4*s-v_3*v_5` as a determinantal
shared-factor gate residual whose stable-source restriction is `x^2*y*z`.
The observable command proves that any rational semiconjugate quotient carrying either
`X_0` or the restricted
quadratic observable has dimension at least 13; its longer rank plateau is
printed as experiment only, while a stacked rank-20 certificate excludes a
common constant translation direction behind that plateau. The final provenance command verifies that the normalized
three-variable canonical contraction fails at its first pure moment, checks
the full twenty-coordinate inverse-recurrence dependency closure after the
known identity slice, and reports the finite stored-circuit census.  Only
the degree-at-most-six invariant statement is a nonlinear quotient-class
obstruction; the circuit census is computation, not minimality.

The independent small-witness audit can also be run directly:

```bash
python3 scripts/audit_dvorsky_gvc5_counterexample.py
```

It verifies the two pre-\(\partial_t\) identities and the resulting GVC(5)
and SIC(5) failures through order eight.  The all-order binomial proof and
the separation from ordinary-Laplacian GVC are documented in
[`DVORSKY_GVC5_COUNTEREXAMPLE.md`](extended-geometry/DVORSKY_GVC5_COUNTEREXAMPLE.md).

The smaller three-pair Image-Mathieu witness is checked by

```bash
python3 scripts/verify_three_pair_image_mathieu_counterexample.py
```

For
\(f=\tau(t-y)(wz+vt)\) and \(g=y\), the dependency-free
checker verifies exact sparse contractions
\(\mathcal E(f^m)=0\) and
\([t]\mathcal E(gf^m)=(-1)^{m-1}(m+1)!m!\) through order ten, records the
four-term bidegree-\((2,2)\) artifact, and replays the two binomial identities used in
the all-order proof.  The same script independently reads the dehomogenized
seed as the four-term cubic Gaussian polynomial
\[
P=(1-Z_2)(W_1Z_1+W_2),\qquad Q=Z_2,
\]
and verifies by exact Wick contraction through order ten that
\(\mathbb E(P^m)=0\) and
\(\mathbb E(QP^m)=(-1)^{m-1}m!\).  The all-order proof and the comparison
with Long's displayed six-term four-real cubic are in
[`THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](extended-geometry/THREE_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).

The sharp two-pair witness is checked independently by

```bash
python3 scripts/verify_two_pair_image_mathieu_counterexample.py
python3 scripts/audit_two_pair_image_mathieu_coefficient_extraction.py
cd formal/finite-etale-keller
lake build FiniteEtaleKeller.SIC2C4FiniteSum
```

For
\[
\begin{gathered}
R=\xi _1z_1+\xi _2z_2,\quad Z=\xi _1z_2,\quad
W=2\xi _2z_1,\quad T=\xi _1z_1-\xi _2z_2,\\
F=(R+Z)\left(R^2W-\frac12(2R+Z)T^2\right),\qquad Q=Z,
\end{gathered}
\]
the dependency-free checker verifies \(T^2=R^2-2ZW\), the sixteen-term
bidegree-\((4,4)\) expansion, the full-rank coefficient-matrix determinant
\(48\), and exact sparse contractions through order eight.  It also
replays the phase-extracted finite sums through order \(99\).  The written
Hopf-coordinate and beta-integral argument proves for every \(m\geq1\)
\[
\mathcal E_2(F^m)=0,\qquad
\mathcal E_2(QF^m)=\frac{(4m+2)!\,m!}{(2m+1)!!}.
\]
The second command audits a non-Gaussian all-order proof.  A formal
constant-term formula for contraction gives the two sums directly; the
pure sum is an \(m\)-th finite difference of a degree-\((m-1)\) polynomial,
while polynomial division reduces the mixed sum to
\[
B_m=\sum_{k=0}^m\frac{(-1)^k\binom mk}{2k+1},\qquad
(2m+1)B_m=2mB_{m-1}.
\]
Thus \(B_m=2^m m!/(2m+1)!!\).  The audit separately checks the chart
identity, both chart constant-term expansions, the divisibility certificate,
finite differences, the general denominator-remainder invariance, and the
recurrence with exact rational arithmetic.  Its cutoff is a regression; the
displayed degree bounds and termwise identities in the written proof are
all-order.  The Lean command
formalizes the general finite-difference cancellation, the
denominator-remainder theorem for alternating quotient sums and its rank-one
endpoint-residue specialization, the specialized normalized products, and
the generalized repeated-pole beta recurrence, the finite
remainder-to-jet identity, and the order-one
factorial/double-factorial evaluation; it also formalizes the scalar chart
identity, the coefficient functional and algebraic beta identity, the
monomial and balanced-array contraction/coefficient-extraction equalities,
the selected chart coefficients, the formal integrals of the resulting
chart polynomials, the product-polynomial evaluation at natural numbers,
and both normalized displayed binomial-sum identities.  It now also
represents (4.3) literally in \(\mathbb Q[v][x,x^{-1}]\), proves the
all-order binomial expansion of its powers, identifies the pure and mixed
constant terms, and proves their final formal-integral values.  It also
defines the original four-variable \(F,Q\) as `MvPolynomial` objects and
proves that their displayed substitution gives the Laurent witness and its
mixed powers.  The remaining Lean integration seam is only a wrapper between
the generic balanced coefficient array in the contraction theorem and
`MvPolynomial.coeff`; the algebraic chart proof itself is formalized on both
sides of that representation boundary.

The positive-characteristic phase diagram is replayed separately by

```bash
.venv/bin/python scripts/verify_two_pair_sic_characteristic_p.py
```

It clears denominators with \(\widetilde F=2F\), checks that the quadric
chart is nondegenerate for every odd prime and that the coefficient tensor
has full rank away from \(2,3\), and verifies
\[
\mathcal E_{2,p}(\widetilde F^m)=0,\qquad
\mathcal E_{2,p}(Z\widetilde F^m)
=\overline{2^m(4m+2)!m!/(2m+1)!!}.
\]
Legendre floor sums, base-\(p\) digit sums, and the Lucas--Kummer carry
criterion all give the exact nonvanishing condition \(4m+2<p\).  The
same audit treats every \(R^k(2F)^r\) of degree \(d=4r+k\) through
degree twenty and verifies the uniform criterion \(dm+2<p\).
At prime-power level it checks the integral radial quotient
\[
A_{s+1}/A_s=16(4s+3)(4s+5)(s+1)^2
\]
and the resulting valuation monotonicity.  It also verifies the
non-radial re-entry for \(R(2F)\): the order-four moment is zero modulo
\(11^2\), but the order-five moment is \(22\) modulo \(11^2\).
The general signed consecutive-order valuation recurrence (4.14e) is
audited through degree sixteen, prime \(31\), and order forty.
The checker also computes the exact coefficient determinants and every
exceptional modular rank for \(R^k(2F)\), \(0\leq k\leq4\), through
degree eight; these are tabulated in (7.6) of the written proof.  The
four binomial-convolution diagonal symbols and their lower-Hessenberg
determinant recurrence are independently compared with expanded
polynomials for every \(0\leq k\leq20\).  The closed characteristic-two
rank formula
\[
\operatorname {rank}C_{R^k(2F)}
=2^{1+s_2(\lfloor(k+2)/2\rfloor)}
\]
is checked through \(k=128\).
For the non-power profiles it checks the universal necessary cutoff
\(4hm+2<p\), proves the closed height-two formula
\[
C_{2,m}=\frac{4^m m!}{\prod_{j=0}^m(4j+1)},
\]
and records the first higher-profile numerator-prime holes at
\((h,m,p)=(6,1,47)\) and \((4,5,89)\).
The checker also verifies the characteristic-two one-sided degeneration, the
characteristic-three Hilbert--Mumford unit ideal, and the naive Hasse
formulas
\(\mathcal H_2(\widetilde F^m)=16^m\) and
\(\mathcal H_2(Z\widetilde F^m)=2m16^m\).
It additionally checks the binomial intertwining and Lucas no-carry units
behind the complete Hasse Image-kernel theorem and its reduction to the
\(p\)-typical operator orders \(1,p,p^2,\ldots\).
For that modified Image it verifies the one-pair counterexample
\[
f=\xi z^p,\qquad g=z,
\]
whose pure moments vanish because \(\binom{pm}{m}=0\), while the mixed
moments are nonzero at every \(m=(p^e-1)/(p-1)\).
The written proof in
[`TWO_PAIR_SIC_CHARACTERISTIC_P.md`](extended-geometry/TWO_PAIR_SIC_CHARACTERISTIC_P.md)
shows division-freely that the Image-kernel identity survives in every
characteristic.  Frobenius gives
\(\mathcal E_{r,p}(f^p)=f(0,z)^p\), so the single \(p\)-th pure moment
forces the dual-degree-zero part to vanish; every fixed mixed contraction
then vanishes for all \(m\geq p\).  This proves ordinary
\(\operatorname{SIC}(r)\) for every \(r\) and \(p>0\), with sharp cutoff
witness \(f=\xi _1,\ g=z_1^{p-1}\).  The finite replay is not being used
as an all-order or periodicity argument.

The Frobenius/\(p\)-curvature bridge is tested by

```bash
.venv/bin/python scripts/research_two_pair_sic_frobenius_curvature.py
```

For
\[
M_{d,r}(m)=4^{rm}(dm+2)!((rm)!)^2/(2rm+1)!,
\]
it derives the coprime minimal order-one recurrence
\(A_{d,r}(m)M(m+1)=B_{d,r}(m)M(m)\) at nine radial rows.
At every good prime its recurrence-operator \(p\)-curvature is proved and
directly replayed as
\[
\prod_{i=0}^{p-1}\frac{B_{d,r}(m+i)}{A_{d,r}(m+i)}
=d^d(m^p-m)^d.
\]
Separately, the normalized angular beta period has Picard--Fuchs operator
\(\theta(2\theta+1)-x(\theta+1)^2\).  Its differential \(p\)-curvature is
computed at every odd prime through \(101\); it is always nonzero,
square-zero of rank one, with poles only at \(0,2\).  This bounded
differential calculation is not an all-prime proof.  The correlation audit
shows why neither curvature recovers the exact phase diagram: first radial
lifts have one common Picard--Fuchs operator, while the recurrence shift
norm cancels the separate zero/pole factors responsible for prime-power
re-entry.  The exact reusable mechanism is instead the local rule
\[
v_p(M(m+1))-v_p(M(m))=v_p(B_{d,r}(m))-v_p(A_{d,r}(m)).
\]
The status and the degree-eight same-curvature/different-phase control are
in
[`TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md`](extended-geometry/TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md).
The resulting integral-lattice postprocessing stage is incorporated in
[`HOLONOMIC_HYPERGEOMETRIC_ALGORITHMS.md`](extended-geometry/HOLONOMIC_HYPERGEOMETRIC_ALGORITHMS.md)
and in the next-step protocols for the bidegree-\((3,3)\) and rank-two
bidegree-\((4,4)\) recurrence programmes.

The exact local geometry of this displayed \(F\) is checked by

```bash
.venv/bin/python scripts/verify_two_pair_counterexample_local_moduli.py
```

The all-order fourth-order continuation is checked separately by

```bash
.venv/bin/python scripts/verify_two_pair_counterexample_fourth_order.py
```

It derives the combined shifted beta-tail for numerator degrees
\(4,8,12,16\), rather than imposing an independent cutoff on each degree.
The exact \(3220\)-by-\(455\) universal polynomial-section system has rank
jump \(90\) to \(91\).  Restoring the eleven free cubic-lift parameters at
the reduced direction \((1,2,3,4,5)\) leaves one affine-linear equation
and one rank-one quadric, of dimension nine and degree two.  Its
discriminant square class is \(41\), so the fiber is two conjugate affine
\(9\)-planes over \(\mathbb Q(\sqrt {41})\) and has no rational point.
This certifies a nonradial geometric fourth-order lift, not a formal arc.

One explicit conjugate pair of fourth-order lifts is continued by

```bash
.venv/bin/python scripts/verify_two_pair_counterexample_fifth_order.py
```

The checker reconstructs the complete
\(\mathbb Q(\sqrt {41})\)-valued jet, verifies its four coefficients
against the all-order fourth tail, and derives the combined
degree-\(4,8,12,16,20\) fifth tail of rank \(56\).  After all thirteen
fourth-tangent corrections, the fifth rank jumps \(2\) to \(3\).  A
primitive obstruction supported only on residual rows \(12,13,14\) has
coefficients
\((2727113757934325760,-407042494824,17047)\).
This obstructs the selected conjugate jets, not the entire
nine-dimensional fourth-lift components.

The algebraization samples and their component-wide fifth obstructions are
computed by

```bash
.venv/bin/python scripts/research_two_pair_counterexample_algebraization.py
.venv/bin/python scripts/research_two_pair_counterexample_fifth_component.py
.venv/bin/python scripts/analyze_two_pair_counterexample_fifth_factor.py
```

For the generic rational direction \((2,-1,3,1,-2)\), the documented
direction \((1,2,3,4,5)\), and the pure apolar-odd direction
\((0,1,0,0,0)\), the exact fourth fiber again has dimension nine, degree
two, and discriminant square class \(41\).  At one point on each conjugate
pair, restoring the previously omitted eleven-dimensional cubic-tangent
kernel changes the fifth coefficient/augmented ranks from \(2/3\) to
\(4/5\), so the selected points remain obstructed.  The component-wide
commands then parameterize all nine coordinates on one component and
restore the full eleven-dimensional kernel.  After the rank-two
fourth-tangent image is eliminated, the remaining coefficient and
augmented ranks over the component function field are \(2/3\).  Every
coefficient \(3\)-by-\(3\) minor vanishes identically, while an augmented
\(3\)-by-\(3\) minor is a nonzero constant with nonzero quadratic norm.
Thus both conjugate affine \(9\)-planes are uniformly obstructed at fifth
order for all three directions.  The exact \(F_{1+s,1}\) control is
polynomial of parameter degree three, so its coefficients at orders
\(4,\ldots,12\) vanish.

The final command is a fast exact replay from nine stored component samples
on \(h_3=h_4=0\).  Within a quadratic projective ansatz, eight samples
reconstruct a selected constant augmented minor, up to a rational chart
unit, as
\(h_0A_1+\sqrt{41}B_2\), with \(A_1\) linear and \(B_2\) a nonsingular
quadratic form.  The ninth sample checks the reconstruction exactly.  Its
vanishing locus is a smooth conic over \(\mathbb Q(\sqrt{41})\); this
candidate exceptional locus has no rational projective point, by exact binary
discriminants.  This finite reconstruction does not prove the ansatz or a
universal identity, and it does not yet prove that the other augmented minors
have no common zero on the conic over \(\mathbb Q(\sqrt{41})\).
See
[`TWO_PAIR_COUNTEREXAMPLE_ALGEBRAIZATION_RESEARCH.md`](extended-geometry/TWO_PAIR_COUNTEREXAMPLE_ALGEBRAIZATION_RESEARCH.md).

The minimum-degree separating invariant and the low-degree invariant-ring
calculation are checked independently by

```bash
.venv/bin/python scripts/verify_two_pair_counterexample_missing_invariant.py
```

The degree-four moment-field continuation is replayed by

```bash
.venv/bin/python scripts/verify_degree_four_tau_even_parameters.py
.venv/bin/python scripts/research_degree_four_moment_field.py \
  --max-weight 16 --targets odd-square
.venv/bin/python scripts/verify_degree_four_diagonal_moment_field.py
.venv/bin/python scripts/verify_degree_four_single_phase_moment_fields.py
.venv/bin/python scripts/research_degree_four_phase_one_chart.py \
  --prime 101 --orders 1 2 3 4 5 6 7 8 9 10 \
  --threads 6 --timeout 300 --groebner-basis 2 \
  --compare-even-parameters 2 3 5 7 11 69 3 6 \
  --test-apolar-orbit 2 3 5 7 11 69 3 6 --certify-example
```

The first command constructs twenty-two algebraically independent
apolar-even trace invariants of degrees \(1,2^4,3^9,4^8\).  Their exact
modular Jacobian rank is \(22\), and their combined cotangent matrix with
\(\mu_1,\ldots,\mu_{22}\) still has rank \(22\).  The second command is a
bounded exact search: it proves that no relation
\(Q(\mu)+c_{234}^2P(\mu)=0\) of invariant weight at most sixteen exists
with that support.  It does not exclude a higher-weight denominator or a
higher even minimal polynomial.  The third command requires Singular.  On
the diagonal quartic slice it proves that the first five moments give a
finite parameter ring of quotient length \(120\), and that the complete
first-six-moment fiber through \((2,3,5,7,11)\) consists exactly of that
point and its reversal.  Finiteness then proves that the full diagonal
moment field is the reversal-fixed field, of exact generic degree two.
The reversal is also the \(\operatorname{SL}_2\) Weyl action on the
diagonal space, so these are two raw parameter points but one invariant
quotient point; this is a fixed-locus control, not a degree-two quotient
test.
The fourth command also requires Singular.  It repeats the finite-parameter
and exact-fiber calculation on all ten coordinate choices of a
\(\tau\)-even positive/negative direction pair in phases \(1,2,3,4\).
Every resulting six-dimensional parameter space has parameter quotient length
\(360\), and its first-seven-moment fiber through
\((2,3,5,7,11,221)\) is exactly the reduced reversal pair.  Openness
then gives exact degree two and fixed-field equality for a nonempty
Zariski-open family of raw direction-pair parameter spaces in every
phase.  The odd cubic is nonzero on exactly four coordinate cross-pairs
in phases one and two, so only those four are genuinely apolar-moving
quotient tests; the other six lie in the fixed locus.
The fifth command is an exact \(F_{101}\) experiment using `msolve`.  On
the eight-dimensional chart containing both positive and both negative
phase-one directions, the first ten moments have a reduced four-point
fiber.  All moments through order eleven and all twenty-two known even
parameters agree on the extra branch, while \(c_{234}\) changes from
\(11\) to \(-11\).  A four-equation orbit basis proves that this branch
is \(\operatorname{SL}_2\)-conjugate to the apolar reversal.  Thus the
four raw points form two candidate quotient points.  The same checker
verifies the exact rational branch \(u=5/3,w=6\): its first eleven moments
agree over \(\mathbb Q\), its odd cubic is \(-1728\) versus \(1728\), and
the matrix with rows \((0,-1/\sqrt3)\), \((\sqrt3,0)\) conjugates
\(\tau(p)\) to it.  It also reconstructs the four rational points
\(p,q,\tau(p),\tau(q)\) and verifies that the first-eight-moment Jacobian
is nonzero at each, so all four are reduced and isolated in
characteristic zero.  Fiber completeness remains proved only modulo
\(101\): additional characteristic-zero components have not been
excluded, and no characteristic-zero generic-degree conclusion is
claimed.
The full \(22\)-dimensional degree-four moment-field equality remains
open; see
[`DEGREE_FOUR_MOMENT_FIELD.md`](extended-geometry/DEGREE_FOUR_MOMENT_FIELD.md).

The completed-coordinate comparison in degrees three through five is
replayed by

```bash
.venv/bin/python scripts/research_completed_moment_algebra.py \
  --degrees 3 4 5 --max-weight 10
```

This exploratory checker constructs the quadratic Casimir decompositions,
verifies the coefficient rows
\(\binom{2d+1}{d-r}\) in \(\mu_2\), and gives exact modular moment-Jacobian
ranks \(13,22,33\).  It runs linear-denominator nonrelation searches for
the missing quadratics over the moments and over the moments with \(q_2\).
For the square of a first apolar-odd invariant it additionally tests the
proposed \((q_2,q_6)\) base and the full quadratic completion.  It also
records bounded Hilbert-series necessary tests for natural and minimally
corrected parameter-degree sequences, and checks the propagated
all-moment-zero witnesses in degrees four and five.  A zero relation
intersection is an exact bounded nonexistence certificate; Hilbert
compatibility is not a proof of a nullcone zero fiber.  See
[`COMPLETED_MOMENT_ALGEBRA_RESEARCH.md`](extended-geometry/COMPLETED_MOMENT_ALGEBRA_RESEARCH.md).

The automatic missing-invariant and \(d=6\) extension is replayed by

```bash
.venv/bin/python scripts/research_completed_moment_algebra.py \
  --degrees 3 4 5 6 --invariant-cutoff 6 \
  --skip-relation-tests --power-witness-cutoff 12 \
  --ladder-beta-check 32 \
  --output artifacts/generated-results/automatic_missing_invariants_d3_d6.json
```

The refined weight-zero-minus-weight-two calculation splits the invariant
spaces by the apolar involution through polynomial degree six and subtracts
the moment-monomial subspace.  It proves that the first missing degree is
two, with even multiplicity \(d-1\), and enumerates the first odd cubic
triples in degrees four through six.  It also certifies full modular
moment-Jacobian ranks \(13,22,33,46\), runs the Hilbert necessary tests for
candidate augmented parameter systems, and evaluates \(q_2\) on the
propagated all-order witnesses through \(d=6\).  The conclusion that \(q_2\)
removes the recorded witnesses is not a classification of every semistable
moment-zero component.  The same run verifies the radial Casimir recurrence
\[
q_{2r}^{(d+1)}(Rf)=(d-r+1)(d+r+2)q_{2r}^{(d)}(f),
\]
the resulting all-degree closed formulas for \(q_2,q_4\) on
\(R^{d-4}F_4\), and the exact power-witness pattern through \(F_4^{12}\).
A finite-difference proof using the chart expansion of \(F_4^m\) shows
for every \(m\geq1\) that all earlier quadratics vanish and
\(q_{2\lceil m/2\rceil}(F_4^m)\ne0\).  The run regresses the exact beta
sums used in that proof through \(m=32\).  The stronger formula listing
every surviving torus phase is recorded through \(m=12\) and remains
bounded evidence.

The completed-invariant comparison specialized to bidegree \((3,3)\) is
replayed by

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_casimir_fiber.py
```

The default run does not launch the older direct boundary standard-basis
calculations.  It compares
\((\mu_1,\ldots,\mu_{12},\mu_{14})\),
\((\mu_1,\ldots,\mu_{12},q_2)\), and full-rank mixed
moment/Casimir systems of the same total invariant degree \(92\).  It
also evaluates the complete weight-\(14\) monomial spaces generated by
the lower moments with \(q_2\), and with \(q_2,q_4\).  The modular ranks
prove over characteristic zero that \(\mu_{14}\) is independent from
the pure Casimir span modulo the lower-moment span.  Hilbert compatibility
does not prove that any displayed zero fiber is the nullcone.

The new null-quadratic synchronization experiment is included explicitly by

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_casimir_fiber.py \
  --run-normal-rank-locus \
  --run-residual-normal-probe \
  --run-complete-normal-fiber \
  --run-exceptional-normal-fibers \
  --timeout 300 \
  --power-bound 50
```

On the normalized chart \(F_2=X^2\), the synchronization ideal has the
seven normal coordinates \(s_3,s_4,s_5,s_6,t_2,t_3,t_4\).  The checker
computes the exact linear normal-symbol matrix of
\(\mu_2,\mu_3,\mu_4\), proves its generic rank is three, and proves that
the common divisor of its nonzero maximal minors is the displayed
irreducible cubic \(P\).  It then records two finite-field calculations
at \(p=32003\): the residual rank-drop locus after dividing by \(P\),
and the complete normal fiber of \(\mu_2,\ldots,\mu_{12}\) above the
fixed allowed-coordinate point \((20,27,36,47,60)\).  The latter has
dimension zero and quotient length \(195\).  Good reduction therefore
proves characteristic-zero transverse isolation at this point and on a
nonempty open subset of the synchronized chart.  The recorded coordinate
powers are finite-field certificates only.  The same run decomposes the
residual rank support exactly into two quadratic-field components and one
lower rational locus, proves that all three are disjoint from \(P=0\),
and tests full normal fibers at good reductions of exact algebraic points
on all four exceptional strata.  The three top-dimensional exceptional
fibers have quotient length \(195\); the lower locus has length \(197\).
Consequently transverse isolation holds on a nonempty open subset of
every linear-rank stratum.  Prefix computations through orders nine and
ten timed out and make no minimal-cutoff claim; proper closed subsets
inside the exceptional strata, the \(F_2=0\) chart, and global nullcone
equality remain open.

The exact diagonal fixed-field theorem in all three degrees is replayed
by

```bash
.venv/bin/python scripts/verify_completed_moment_diagonal_fields.py
```

For \(d=3,4,5\), the checker proves that the first \(d+1\) diagonal
moments form a parameter system of quotient length
\((d+1)!=24,120,720\).  Homogeneous finite-field standard bases at
\(32003\), followed by projective properness and the regular-sequence
Hilbert series, give the characteristic-zero finiteness statement.
An invertible midpoint/direction change then proves over \(\mathbb Q\)
that the first \(d+2\) moment fiber through the selected integral point is
exactly
\((y_0,\ldots,y_{d-1},s^2-1)\).  Hence the full diagonal moment field is
the reversal-fixed field and has exact generic degree two in each degree.
These are slice theorems, not statements about the full invariant
quotients.

The exact single-phase extension in degrees three and five is replayed
by

```bash
.venv/bin/python scripts/verify_completed_moment_single_phase_fields.py
```

For one matching apolar-eigendirection pair in every nonzero phase, the
checker proves that the first \(d+2\) moments have full Jacobian rank and
that adding \(\mu_{d+3}\) makes the moment-origin fiber finite.  Exact
standard bases over \(\mathbb F_{32003}\), weighted-projective
properness, and Nakayama lift the displayed two-point reversal fiber to
characteristic zero.  The quotient lengths at the special moment origin
are \(54\) for \(d=3\) and \(1934\) for \(d=5\).  Quintic cross-direction
slices in phases one and two have
\(c_{234}=-273686400/7\), so those reversal pairs are genuinely distinct
invariant-quotient points.  The remaining slices are raw parameter-space
fixed-field controls.  None of these slice certificates determines the
generic degree on the full invariant quotient.

The first branchwise global-\(d=4\) \(q_2\)-augmented nullcone attack is
replayed by

```bash
.venv/bin/python scripts/research_degree_four_q2_augmented_nullcone.py \
  --prime 32003 --max-jet 4 --composition native \
  --ordering dp --timeout 300
```

On the branch \(q_2=0,F_2\ne0\), normalize \(F_2\) to a highest-weight
square.  At one deterministic synchronized point, the checker expands
the first twenty-one moments in the twelve forbidden weight coordinates.
Moments two through five provide four formal pivots.  In the remaining
eight variables, the exact quadratic and cubic jet ideals over
\(\mathbb F_{32003}\) have dimensions six and four.  This is a bounded
normal-jet frontier, not a formal-isolation or global-nullcone theorem.
The native eight-variable quartic basis reaches its declared
\(300\)-second timeout, which supplies no mathematical evidence.
The cubic support decomposition and the cheaper dominant-sheet quartic
restrictions are replayed by

```bash
.venv/bin/python \
  scripts/research_degree_four_q2_cubic_decomposition.py \
  --prime 32003 --timeout 180
```

After the forced equation \(x_7=0\), the checker proves the radical of
the \(x_6=0\) cubic support, factors its binary cubic into a rational
sheet and an irreducible quadratic sheet over \(\mathbb F_{32003}\),
and records a distinct degree-nine off-\(x_6\) saturation with an
irreducible generic degree-nine fiber.  Quartic restriction collapses
both dominant sheets to the same three-plane.  The generic quartic
restriction on the off-axis saturation still times out, so formal
isolation, quintic normal terms, and all \(F_2=0\) boundary branches
remain open.  See
[`DEGREE_FOUR_Q2_AUGMENTED_NULLCONE.md`](extended-geometry/DEGREE_FOUR_Q2_AUGMENTED_NULLCONE.md).

The explicit first \(d=2\) moment relation is reconstructed and verified
by

```bash
.venv/bin/python scripts/verify_two_pair_d2_moment_field_relation.py
```

With factorial-normalized centered moments \(x_m\), this produces the
monic quintic relation for \(x_7\) over
\(\mathbb Q[x_2,\ldots,x_6]\).  Five good primes reconstruct its \(241\)
nonzero rational coefficients.  The checker then verifies the relation
exactly at \(418\) fixed integral coefficient matrices whose full ansatz
evaluation matrix has rank \(418\) modulo \(1000003\).  The generated
artifact stores the five sparse coefficient blocks.  The same checker
also proves, by five further modular rank jumps, that the centered
moments \(x_8,\ldots,x_{12}\) are each outside the algebra generated by
their predecessors.  Thus the first seven moments generate the moment
field but not the polynomial moment algebra.  Finally, it certifies that
the affine operator slice
\[
 \begin{pmatrix}0&1&0\\a&b&c\\d&e&f\end{pmatrix}
\]
has dense \(\mathrm{SL}_2\)-saturation and expands reconstructed formulas
for \(x_{13},\ldots,x_{18}\) identically to zero in
\(\mathbb Q[a,b,c,d,e,f]\).  Hence every moment through order \(18\)
belongs to the algebra generated by the first twelve.

The missing-invariant checker constructs the five Casimir projectors for
\(\operatorname{End}(\operatorname{Sym}^4)\), verifies the five quadratic
values and the exact expression of \(\mu_2\) in that basis, evaluates the
ten primitive cubic contractions, and certifies modulo \(1000003\) that
\(\mu_1,\ldots,\mu_{22}\) have Jacobian rank \(22\).  The nonzero modular
minor is an exact characteristic-zero algebraic-independence certificate,
not a numerical rank test.  It also constructs the equivariant apolar
adjoint, verifies its alternating signs on the five irreducible summands,
and checks that the first odd cubic is nonzero.  The written coefficient
reindexing proves that every moment is adjoint-even in all orders; hence
the moment and invariant fields differ and the algebra conductor is zero.
The same weight calculation finds the first odd invariant in bidegrees
two and three, and gives degree ten for the proved first-six-moment
parameter map on \(V_2\).  A nonzero fifteen-by-fifteen modular evaluation
determinant shows that \(\mu_7\) leaves the first-six moment field; the
prime degree-five tower then proves that the complete \(d=2\) moment field
is exactly the apolar-adjoint fixed field, of generic degree two in the
full invariant field.

The local-moduli checker computes the contraction-preserving stabilizer
and orbit, uses the
all-order Hopf coefficient identity to obtain the thirteen-dimensional
all-moment tangent space, records the seven independent quadratic lifting
obstructions, proves that their quotient radical is a five-plane, and
splits the quotient moment tangent as \(2+8\) and that five-plane as
\(1+4\) under the apolar involution.  It constructs a polynomial cubic
lift for every direction on that plane and verifies the defect-preserving
family
\[
F_{a,b}=\frac{aR+bZ}{2}
\left(2W(aR+bZ)^2-2abR^3-b^2R^2Z\right).
\]
The parity split reduces the certified cubic correction system from
\(490\times195\) to a consistent \(154\times73\) equivariant system and
reduces the proposed fourth-order continuation from \(980\times455\) to
at most \(644\times273\).
The generated artifact contains the exact tangent matrix, kernel basis,
obstruction quadrics, radical equations, and cubic correction.  The family
proves a positive-dimensional local quotient but is not claimed to exhaust
every reduced local branch; the first unresolved local equations are fourth
order.

The same checker verifies the propagated forms
\(F_d=R^{d-4}F\) for \(4\leq d\leq10\) through moment four.  The written
radial-shift identity proves, for every \(d\geq4\) and \(m\geq1\),
\[
\mathcal E_2(F_d^m)=0,\qquad
\mathcal E_2(QF_d^m)=\frac{(dm+2)!\,m!}{(2m+1)!!}.
\]
It also verifies through moment three and degrees \(4\leq d\leq15\) the
bounded-radial-order family
\[
G_{r,k}=R^kF^r,\qquad d=4r+k,\quad 0\leq k\leq3,
\]
whose all-order formula is
\[
\mathcal E_2(G_{r,k}^m)=0,\qquad
\mathcal E_2(QG_{r,k}^m)
=\frac{(dm+2)!\,(rm)!}{(2rm+1)!!}.
\]
The written proof shows that \(G_{r,k}\) has exact \(R\)-adic order \(k\).
Thus the family is \(R\)-primitive in degrees divisible by four and has
radial order at most three in every degree.
Finally, the checker verifies for \(1\leq h\leq4\) and \(m\leq3\) the
explicit non-power Hopf-profile family \(\Phi_h\) of degree \(4h\), whose
all-order mixed formula is
\[
\mathcal E_2(Q\Phi_h^m)
=(4hm+2)!\int_0^1
(1-v^2)^m(1+v^2)^{(h-1)m}\,dv.
\]
The written endpoint-contact proof shows that every pure moment vanishes,
the displayed detector is positive, and \(\Phi_h\) is \(R\)-primitive and
not a proper power.

The degree-five bilinear-multiplier obstruction is checked by

```bash
.venv/bin/python scripts/verify_two_pair_degree_five_multiplier_obstruction.py
```

For \(L=aR+bZ+cW+eT\), the checker derives the first four pure moments of
\(LF\).  On the nonradial branch, moments one and two eliminate \(c,a\);
moment three gives
\[
q(u)=8019u^4-623736u^2+3219760,
\]
and moment four gives
\[
p(u)=136323u^6-5359284u^4-174020976u^2-802761152.
\]
Their exact rational gcd is one.  Hence the first four moments force
\(L=aR\), excluding every nonradial bilinear lift of the quartic seed in
degree five.  This does not classify the full \(V_5\).

The sharp primitive finite-prefix family is replayed by

```bash
python3 scripts/verify_two_pair_primitive_prefix_obstruction.py
```

For \(G_{d,\lambda}=R^{d-4}F+\lambda Z^d\), the checker directly verifies
in degrees \(4\leq d\leq8\) that moments one through \(d\) vanish and
\[
\mathcal E_2(G_{d,\lambda}^{d+1})
=(d+1)\lambda(d(d+1)+1)!\frac{d!}{(2d+1)!!}.
\]
The written phase-support proof establishes the formula for every
\(d\geq4\).  For \(\lambda\ne0\), these are \(R\)-primitive sharp-prefix
points, not all-order counterexamples.
The same checker tests the stronger triangular statement for
\(4\leq d\leq7\): if
\(H=\sum c_jR^{d-j}Z^j\) has least nonzero phase \(s\), then
\(R^{d-4}F+H\) is detected exactly at moment \(s+1\).  The written proof is
all-order and excludes every nonzero positive-phase triangular correction.
It also checks the two-sided degree-five extreme ansatz
\(RF+aZ^5+bW^5\): moment two is \(921600ab\), and the remaining two
branches are excluded by explicit nonzero moments four and six.
Finally, for \(4\leq d\leq7\), it verifies the odd-height family
\[
J_{d,\lambda}=R^{d-4}F+\lambda Z^{d-1}T,
\]
whose all-order phase-parity proof gives zero moments below \(2d\) and
\[
\mathcal E_2(J_{d,\lambda}^{2d})
=\binom{2d}{2}\lambda^2(2d^2+1)!
\frac{(2d-2)!}{(4d-1)!!}.
\]
The same replay checks all opposite odd-height monomial pairs, all
opposite even-height pairs of phase \(s\geq3\), and the phase-two
exceptional branch in the displayed degree range against the written
all-degree formulas.  It also verifies the degree-five phase-one
elimination for
\(RF+aZT^4+bWT^4\): the normalized moments \(2,3,4\), their exact
lexicographic remainder, and the nonzero resultant
\(-418538718730248905250\).

The remaining uniform phase-one elimination is reproduced by

```bash
.venv/bin/python scripts/verify_two_pair_phase_one_uniform_obstruction.py
```

It derives the symbolic height-dependent moments, eliminates to a
quadratic and cubic in the remaining correction parameter, and factors
their resultant.  The only nonlinear height factor has degree \(31\)
and no root modulo \(29\).  Together with the preceding checker, this
excludes every correction
\(aZ^sT^{d-s}+bW^sT^{d-s}\) for every \(d>4\).

The first multi-pair local obstruction is checked by

```bash
.venv/bin/python scripts/verify_two_pair_degree_five_odd_height_quadratic_rigidity.py
```

For the degree-five odd-height correction
\(a_2Z^2T^3+a_4Z^4T+b_2W^2T^3+b_4W^4T\), it constructs the quadratic
coefficients of moments \(2,\ldots,10\).  Their exact Gröbner basis
contains the cube of every parameter, so their radical is the homogeneous
maximal ideal.  It also computes the Hilbert vector \((1,4,1)\), length
six, and a nondegenerate socle pairing, identifying the quotient as a
compressed quadratic Artin--Gorenstein algebra.  This is a local formal
obstruction.  The checker then treats the full ten-dimensional monomial
correction space: moments \(2,\ldots,7\) eliminate the six even-height
second-order variables, and projected moments \(8,\ldots,11\) give four
quadrics forming a length-sixteen complete intersection with Hilbert
vector \((1,4,6,4,1)\).  Thus \(RF\) is formally isolated in this full
space; this is not a global classification of finite corrections away
from \(RF\).

The higher-degree continuation is checked by

```bash
.venv/bin/python scripts/verify_two_pair_higher_degree_monomial_formal_rigidity.py
```

Using exact Hopf-angular integration, it treats the full monomial
correction spaces in degrees six and seven.  The consecutive linear
blocks have sizes six and eight, respectively.  In both degrees, the
projected odd-height obstruction consists of six quadrics forming a
complete intersection with Hilbert vector
\((1,6,15,20,15,6,1)\) and length \(64\).

The all-degree linear pivot theorem has the dependency-free regression

```bash
python3 scripts/verify_two_pair_linear_pivot.py
```

The written proof uses the endpoint factorization
`2p=-(1+x)(2+x)+(1-u)x^(-1)(1+x)^3` over
\(\mathbb Z_{(2)}\).  With
\(\delta_s=s+\nu_2(s!)\), the scaled local Smith exponents are two copies
of \(\delta_2,\delta_4,\ldots,\delta_d\) in even degree.  In degree
\(d=2h-1\), they are one \(\delta_1\), two copies of
\(\delta_3,\ldots,\delta_{2h-1}\), and one \(\delta_{2h+1}\).  The
checker independently reconstructs the terminating binomial entries and
verifies this pattern exactly in degrees 5 through 25.  This finite range
is a regression for the filtered proof, not the proof itself.

Degrees eight through eleven use the scalable good-prime replay

```bash
.venv/bin/python scripts/verify_two_pair_degree_eight_monomial_formal_rigidity.py
```

The checker constructs the exact rational pivot and projected systems,
reduces them at \(1000003\), and invokes Singular.  Degrees eight and nine
give eight-quadric quotients of dimension \(256\).  Degree ten gives a
ten-quadric quotient of dimension \(1024\), and degree eleven gives the
same ten-quadric dimension.  The corresponding ninth and eleventh
variable powers reduce to zero, proving the characteristic-zero complete
intersections.
See
[`TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](extended-geometry/TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).
This proves `not MN_d` for every `d>=4`, `not SIC(2)`, and, with the known
one-pair theorem, the exact minimum failing pair dimension two.  The finite
replay is not being used as the all-order proof.

The coefficient-rank frontier inside bidegree \((4,4)\) is replayed by

```bash
python3 scripts/verify_two_pair_sic_bidegree44_rank_frontier.py
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_invariants.py
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_all_order_audit.py
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_direct_chart.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree44_two_row_channel.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree44_two_row_boundaries.py
.venv/bin/python scripts/verify_two_pair_sic_bidegree44_rank_two_swap_slice.py
.venv/bin/python scripts/verify_two_variable_quartic_squarefree_pivot.py
.venv/bin/python scripts/verify_two_variable_quartic_two_root_finite.py
```

The dependency-free checker verifies the determinantal dimensions
\(r(10-r)\) for \(1\leq r\leq4\), exact representatives of ranks one
through four, and the coefficient formulas for pure contractions and all
four bilinear mixed multipliers through order four. It also checks the
nilpotent endomorphism trace screen on a fixed-flag one-sided chart and
replays the known rank-five determinant and mixed formula. The written
split-symbol argument excludes rank one for arbitrary SIC multipliers,
giving the rigorous interval \(2\leq r_{\min}\leq5\). Ranks two, three,
and four remain open; finite rank-two residuals are not treated as an
exact counterexample or a lower bound. See
[`TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md`](extended-geometry/TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md).

The swap-slice checker treats
\(F_P=\xi_1^4P(z_1,z_2)-\xi_2^4P(z_2,z_1)\).  It proves that odd moments
vanish by involution and that the even moments of orders
\(2,4,6,8,10\) have a length-twelve zero scheme whose radical consists
of \(P=(z_2\pm z_1)^4\).  Seven good primes reconstruct the triangular
lex ideal; exact rational reductions, the matching special-fiber length,
and four unit charts on the complete \(a_4=0\) projective boundary certify
the characteristic-zero result.  Both reduced points have coefficient
rank one.  The parity-even three-parameter core is already excluded by
orders \(2,4,6\).  This is an exact structured-slice theorem, not a global
rank-two exclusion.

The direct quotient checker fixes \(U=[I_2;A]\), verifies the pure and
mixed relative-period identities, and certifies the one-sided
function-field recurrence and mixed cutoff.  The two-row checker then
computes the exact degree-\(604\) seven-moment scheme and the
eighth-moment unit ideal on the dense support.  Its boundary companion
classifies all 135 proper support orbits: every mixed orbit has a sharp
cutoff at most seven, including separate rank-one and exact-rank-two
localizations on the three mixed-rank tori.

The cross-degree rank-stratified finite-prefix census is replayed by

```bash
python3 scripts/verify_rank_stratified_moment_census.py
```

For every determinantal rank stratum through rank four in balanced
degrees two through four, the checker computes the exact
diagonal-\(\mathrm{SL}_2\) invariant Hilbert coefficients through degree
\(85\) and a good-prime Jacobian of the dimension-sized consecutive
moment system. It recovers the certified degree-four rank-two coefficient
\(-5266\) and ambient degree-three coefficient \(-2186\). The new
rank-two cubic calculation has full nine-moment Jacobian and
\[
[t^{29}]H_{3,2}(t)\prod_{m=1}^{9}(1-t^m)=-58.
\]
Since the first four moments already cut the cubic rank-one Segre cone
down to the nullcone, the resulting semistable nine-moment point has
exact coefficient rank two. It is existential and finite-prefix only.
The first tested single replacement with full Jacobian and a nonnegative
candidate numerator through degree \(85\) is
\(\mu_1,\ldots,\mu_8,\mu_{12}\); orders ten and eleven do not pass. See
[`RANK_STRATIFIED_MOMENT_PROGRAM.md`](extended-geometry/RANK_STRATIFIED_MOMENT_PROGRAM.md).

The first holonomic probe on the degree-three rank-two row is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_holonomic_probe.py
```

It verifies the exact beta/constant-term period formula at two integral
exact-rank-two points, compiles the adjacent C++ moment engine in a
temporary directory, and computes 501 normalized moments at both points
over three primes.  In all six cases an order-\(27\), \(m\)-degree-\(11\)
recurrence ansatz uses 335 equations and passes 139 unused equations.  The
monic forward coefficient is common to every probe, reconstructs exactly,
and has no nonnegative integer root.  An order-\(27\), degree-\(10\)
ansatz fails at both points modulo \(1000003\).  These are exact bounded
computations and modular evidence for a universal recurrence shape, not a
creative-telescoping certificate.  The universal parameter denominator,
exceptional-locus stratification, and corrected-system bridge moments
remain open.  See
[`TWO_PAIR_SIC_BIDEGREE33_RANK_TWO_HOLONOMIC_PROBE.md`](extended-geometry/TWO_PAIR_SIC_BIDEGREE33_RANK_TWO_HOLONOMIC_PROBE.md).

The relative-cohomology refinement is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_relative_jacobian.py
```

It computes the exact logarithmic Jacobian length
\(18=2_{t=0}+2_{t=1}+14_{\rm interior}\), including the saturation
exponent six and an explicit eighteen-monomial basis.  At the same two
points and three primes it finds an order-\(18\), \(m\)-degree-\(18\)
recurrence with 83 unused equations; degree 17 fails at both points
modulo \(1000003\).  Its forward coefficient has eight common linear
factors and a point-dependent decic.  The leading \(m^{18}\) coefficient
has nonzero remainder in all eighteen relative-Jacobian coordinates, so
a naive degree-17 polynomial divergence certificate cannot prove the
recurrence.  The result remains exact quotient data plus modular
recurrence evidence, not a universal telescoping certificate.

The characteristic-zero cyclic splitting at the first integral rank-two
point is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_interior_cyclic_split.py
```

The checker proves the pairwise-comaximal critical-algebra decomposition
\(18=14_{\rm interior}+2_{t=0}+2_{t=1}\).  Eliminating
\(P=Q/u^3\) gives pairwise-coprime squarefree polynomials of degrees
\(14,2,2\), whose product is the degree-\(18\) relative eliminant up to
scalar.  Hence \(1,P,\ldots,P^{13}\) is an exact cyclic basis of the
interior algebra.  It also proves
\[
 ((uQ_u-3Q,tQ_t):(ut)^\infty)=I_{\rm interior},
\]
so the exact toric logarithmic critical rank at this fiber is \(14\).
Exact evaluation of the Chinese-remainder idempotents
against \(\nu_0,\ldots,\nu_{18}\) gives nonzero values on both endpoint
pairs.  Therefore the raw period does not descend through the ordinary
Jacobian quotient.  The same run computes the exact first divergence
seed
\[
u^{47}p_{\rm int}(P)=X(uQ_u-3Q)+YQ_t.
\]
The lift has 6750 and 6791 terms, its divergence has 6749 terms, and the
two nonzero endpoint restrictions have 45 and 85 terms.  An
\(m\)-dependent reduction of these retained terms is still required.
Their first ordinary normal forms occupy all \(14+2+2\) coordinates;
the discarded gradient parts must be lifted recursively rather than
silently set to zero.
The artifact is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_interior_cyclic_split.json`.

The exact fixed-fiber rational \(D\)-module seed is replayed by

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_dmodule_picard_fuchs.py \
  --annihilator-only
```

This command requires Macaulay2, the `BernsteinSato` package, and
Normaliz.  It rewrites the generating integrand as
\(u^2/(u^3-zQ(u,t))\), computes 34 first-order annihilators of
\((u^3-zQ)^{-1}\), and then computes the 76-generator annihilator of the
specific numerator \(u^2\).  Both left ideals are checked exactly to be
holonomic of rank one.  The artifact is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_dmodule_picard_fuchs_research.json`.
This is an all-order integrand certificate, not yet the integrated
Picard--Fuchs operator.  Omitting `--annihilator-only` requests the long
sequential \(t,u\) pushforward; its output must still pass the relative
endpoint audit.

The exact shift-Ore comparison of the two sampled recurrence shapes is
replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_ore_gcd.py
```

At the same two points and three primes, the checker recomputes the
order-\(18\) and order-\(27\) operators through moment 500 and performs
left Euclidean division in
\(\mathbb F_p(m)[S;\,Sf(m)=f(m+1)S]\).  The order-\(27\) operator is not
a left multiple of the order-\(18\) operator.  Their greatest common
right divisor has order \(14\), primitive coefficient degree \(58\), and
left quotient orders \(4\) and \(13\), respectively.  The primitive
order-\(14\) operator is checked directly on all 487 available moment
rows at every sample.  Its order matches the exact interior length in
the \(2+2+14\) relative-Jacobian decomposition.  The artifact is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_gcd.json`.
This is an exact bounded modular factor calculation, not a universal
Picard--Fuchs certificate.

Three bounded research probes test shortcuts from the sampled factor to
a relative telescoping certificate:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_ore_reconstruct.py \
  --primes 1000003 1000033
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --steps 1
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_leading_syzygy.py
```

The first exposes all 885 modular coefficients of the fixed-fiber
order-\(14\), degree-\(58\) factor and tests balanced reconstruction at a
held-out prime.  The latter two work in the exact length-eighteen
relative quotient modulo \(1000003\).  They prove that the leading
\(m^{58}\) class is nonzero and that the divergence classes of all
leading Koszul corrections \(R(C,-A)\) have rank zero.  Therefore the
direct zero-boundary polynomial certificate cannot start; the full
\(14+2+2\) endpoint-extended connection is required.  These artifacts
record exact modular no-go calculations for those ansätze, not a
characteristic-zero Picard--Fuchs certificate.

The completed all-order certificate at the first point modulo \(1000003\)
is produced in three restartable Laurent-reduction chunks:

The chunks are intentionally local cache material (large checkpoints and raw
Singular certificates), not Git artifacts.  Run the following target; override
`LOCAL_CERTIFICATE_CACHE` to place the cache elsewhere:

```bash
make generate-rank-two-divergence-local
```

Its expanded commands are:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --mode interior --mapped-quotient --steps 20 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m38.poly \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_certificate_m58_m38.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_research.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --mode interior --mapped-quotient \
  --checkpoint-input \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m38.poly \
  --steps 20 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m18.poly \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_certificate_m38_m18.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_research_m38_m19.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --mode interior --mapped-quotient \
  --checkpoint-input \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m18.poly \
  --steps 18 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m0.poly \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_certificate_m18_m0.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_interior_divergence_research_m18_m1.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_terminal_syzygy_block.py \
  --output artifacts/local/two_pair_sic_bidegree33_rank_two_terminal_syzygy_block_research.json \
  --R-output artifacts/local/two_pair_sic_bidegree33_rank_two_terminal_syzygy_R.sing
```

The uncorrected descent has a 307276-term terminal residual.  The final
command constructs the 298606-term Koszul correction \(R\) and the exact
identity \(T=QH(R)\).  The complete independent replay is

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_all_order_certificate.py
```

It replays all 58 coefficient identities in Singular, checks every
restart residual, independently verifies the terminal correction, and
proves that both corrected endpoint exponential-polynomials vanish.
The resulting theorem is the order-\(14\) recurrence for every
\(m\geq0\) at this fixed fiber over \(\mathbb F_{1000003}\).  The summary
artifact is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_all_order_certificate.json`.
This does not reconstruct characteristic zero or a generic rank-two
parameter identity.

The finite characteristic-zero lift at the same fixed fiber is rebuilt in
three stages.  The first command is resumable but expensive: from an empty
cache, computing the 205 exact modular images takes roughly 45 minutes on
the reference machine.

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_ore_reconstruct.py \
  --point 0 --prime-count 205 --prime-start 1000000 \
  --holdout-count 5 \
  --image-cache \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_reconstruct_images.json \
  --output \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_reconstruct_research.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_simultaneous_reconstruct.py \
  --kind common --prime-count 205 --holdout-count 5 \
  --cache \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_reconstruct_images.json \
  --output \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_characteristic_zero_lift.json
```

The simultaneous LLL reconstruction uses the first 200 images and reserves
five as fresh holdouts.  It returns a primitive integer order-\(14\),
degree-\(58\) operator whose largest coefficient has 2397 bits.  Generate
the first exact characteristic-zero divergence level with

```bash
make generate-rank-two-char0-leading-local
```

Equivalently:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --operator \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_ore_characteristic_zero_lift.json \
  --mode interior --mapped-quotient --steps 1 --timeout 900 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_char0_interior_divergence_checkpoint_m57.poly \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_char0_interior_divergence_certificate_m58_m57.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_char0_interior_divergence_research_m58_m57.json
```

On the reference machine this producer takes about 255 seconds and peaks
near 1.6 GB.  The combined independent verification is

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_characteristic_zero_lift.py
```

It replays all 205 modular images, 27 exact rational moment identities, the
known forward factor, and the \(m^{58}\) divergence identity over
\(\mathbb Q\).  Its Singular replay takes about 93 seconds and peaks near
1.2 GB.  The manifest is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_characteristic_zero_lift.json`.
This is exact finite verification, not an all-order characteristic-zero
certificate: levels \(m^{57}\) through \(m^1\), the terminal syzygy, and
both endpoint identities remain open.

An experimental reduction-based Picard--Fuchs route is retained in
`scripts/research_two_pair_sic_bidegree33_rank_two_picard_fuchs.jl`.  With
`MultivariateCreativeTelescoping.jl` 0.1.3 installed in an isolated Julia
environment, the generic projective CRT run is

```bash
timeout 900 julia --project=/tmp/sic33-mct-env \
  scripts/research_two_pair_sic_bidegree33_rank_two_picard_fuchs.jl \
  crt 1 original
```

The original-coordinate reference run reached the 900-second cap without an
operator.  The exact beta compression \(x=u, y=ut/(1-t)\) instead gives the
sixteen-term form

\[
 \frac{x+y}{(x+y)^3-z\Phi(x,y)}.
\]

The compact closed-cycle calculation finishes in roughly eight minutes and
peaks near 2 GB:

```bash
timeout 900 julia --project=/tmp/sic33-mct-env \
  scripts/research_two_pair_sic_bidegree33_rank_two_picard_fuchs.jl \
  crt 1 compact \
  artifacts/local/two_pair_sic_bidegree33_rank_two_compact_picard_fuchs.ore
```

It returns a differential order-eight closed-cycle operator.  The interval
period has a different inhomogeneous order-eight relation.  Generate 100 exact
modular images of that relation, the reference-prime structural comparison,
and its resumable local cache with

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_compact_relative_picard_fuchs.py \
  --prime-count 100 --prime-start 1000000 --jobs 4
```

From an empty cache this takes about six minutes on the reference machine.
It proves at \(p=1000003\) that the differential residual has degree 55,
converts its tail to an order-64, \(m\)-degree-eight shift operator \(R\), and
checks the exact modular factorization \(R=Q_{50}G_{14}\).  Reconstruct the
characteristic-zero differential operator from 95 images with five holdouts:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_simultaneous_reconstruct.py \
  --kind compact --prime-count 100 --holdout-count 5 \
  --dimensions 8 12 16 24 32 --offsets 8 \
  --output \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_compact_relative_pf_characteristic_zero_lift.json
```

The exact bridge verifier replays all 100 images, exact rational forcing and
moment rows, 50 exact initial \(G_{14}\)-identities, positivity of the quotient
forward denominator, and the rational shift-Ore identity
\(R=Q_{50}G_{14}\):

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_compact_relative_pf_lift.py
```

The exact Ore division takes about eight minutes but stays below 150 MB.  Use
`--skip-exact-ore-division` for the two-second reconstruction, moment, and
forward-coefficient checks only.

The compact modular all-order certificate is generated in two divergence
chunks and one terminal solve:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --operator \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_compact_relative_pf_research.json \
  --mode interior --mapped-quotient --steps 1 --timeout 600 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_checkpoint_m7.json \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_certificate_m8_m7.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_divergence_m8_m7.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_relative_divergence.py \
  --operator \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_compact_relative_pf_research.json \
  --mode interior --mapped-quotient \
  --checkpoint-input \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_checkpoint_m7.json \
  --steps 7 --timeout 900 \
  --checkpoint-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_checkpoint_m0.json \
  --certificate-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_certificate_m7_m0.sing \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_divergence_m7_m1.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_terminal_syzygy_block.py \
  --operator \
    artifacts/generated-results/two_pair_sic_bidegree33_rank_two_compact_relative_pf_research.json \
  --certificate \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_certificate_m7_m0.sing \
  --terminal \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_checkpoint_m0.poly \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_terminal_syzygy_research.json \
  --R-output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_shift_mod1000003_terminal_syzygy_R.sing
```

The independent combined replay is

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_compact_relative_modular_all_order.py
```

It verifies all eight coefficient identities, the 132615-term terminal
Koszul correction, and the zero endpoint trace, proving \(R\nu=0\) for every
\(m\geq0\) over \(\mathbb F_{1000003}\).  The compact characteristic-zero
operator and its exact factorization are proved, but its divergence and
endpoint identities remain open.  A first expanded rational descent level
reached 900 seconds and 8.6 GB; ordinary top-pole Griffiths reduction leaves
an 18-term remainder, so the next engine must use extended relative reduction
or reconstruct the eight modular certificate levels.

For the two-prime support scout, generate a second operator artifact with

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_compact_relative_picard_fuchs.py \
  --prime-count 1 --prime-start 1000003 --jobs 1 \
  --output \
    artifacts/local/two_pair_sic_bidegree33_rank_two_compact_relative_pf_research_mod1000033.json
```

Repeat the two divergence commands and the terminal block solve above with
that operator, replacing `mod1000003` by `mod1000033` in every local output.
Then compare the exact Laurent supports with

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_compact_certificate_support.py
```

All \(Y_r\) supports, five of eight \(X_r\) supports, and the terminal support
agree; the other three \(X_r\) supports differ by one monomial each.  This is
an exact two-prime feasibility scout, not a rational reconstruction.

The exact border-basis calculation on the generic factor pencil is
replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_pencil_border.py
```

Over each of three finite rational-function fields it verifies
saturation exponent six, quotient length eighteen, and six reduced
nineteen-term border relations.  The three distinct monic coefficient
denominators have degrees \(74,88,94\), common gcd degree \(74\), and
coprime quotient degrees \(14,20\); their lcm is squarefree of degree
\(108\).  The artifact is
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_pencil_border.json`.
The same run specializes all \(4,6,5\) base-field roots of that polynomial.
Every specialization remains coefficient-rank two.  Four roots of the
degree-\(74\) component preserve the full \(2+2+14\) relative profile,
while the ten accessible degree-\(14\) roots and one accessible
degree-\(20\) root lower one endpoint or interior length and have total
length seventeen.  This is an exact modular one-pencil classification, not
a universal parameter-space determinant or a classification of the
non-linear exceptional closed points.

The generic-pencil interpolation stress test is

```bash
.venv/bin/python \
  scripts/explore_two_pair_sic_bidegree33_rank_two_recurrence_line.py \
  --samples 256 --holdout 12 --jobs 4 --maximum-moment 390
```

The scaling-family control reconstructs the predicted homogeneous
degrees.  On the generic quadratic factor pencil, eight representative
coefficients have no rational interpolant in the sampled degree window
(combined degree at most 243).  This redirects the universal computation
to an eighteen-dimensional relative connection and determinant
representation instead of an expanded parameter interpolation.

The second checker works on the sixteen-dimensional rank-at-most-two
determinantal variety. At one exact rank-two factor chart it proves modulo
the good prime \(1000003\) that both
\(\mu_1,\ldots,\mu_{13}\) and
\(\mu_1,\ldots,\mu_{12},\mu_{14}\) have Jacobian rank thirteen. It then
computes the diagonal-\(\mathrm{SL}_2\) invariant Hilbert coefficients
from the two-row Cauchy decomposition. Degrees \(1,\ldots,13\) fail the
homogeneous-parameter test exactly:
\[
[t^{69}]H(t)\prod_{m=1}^{13}(1-t^m)=-5266.
\]
Thus their common rank-at-most-two zero fiber necessarily contains a
semistable point. The corrected degrees \(1,\ldots,12,14\) pass the
necessary Hilbert test through degree \(100\). Their numerator is
palindromic through degree \(82\), the top degree predicted by the
determinantal invariant ring's \(a\)-invariant \(-10\) and the candidate
parameter-degree sum \(92\), but their zero fiber remains open. The
rank-one finite-cutoff checks below close all collided-root strata, but
one squarefree uniform-specialization chart remains open. Thus exact rank
two is not yet forced for the semistable thirteen-moment point. No
all-order counterexample is claimed.

The third dependency-free command audits why recurrence work is currently
parked. It verifies
\[
\frac{\mu_m}{(4m+1)!}
=\operatorname{CT}_u\int_0^1
\Phi_C(1,u,t,(1-t)/u)^m\,dt
\]
on the rank-two factor chart and records the corresponding rational
generating function. It also proves that the displayed exact rank-two
Jacobian point has \(\mu_1=7414\), so it is not the existential
thirteen-moment survivor, and computes the generic rank-two Newton polygon
with normalized volume \(48\). The semistable fiber has no recorded closed
point or residue field, so no coefficient-specialized scalar recurrence
is claimed and \(\mu_{14}\) is not evaluated. See
[`TWO_PAIR_SIC_BIDEGREE44_RANK_TWO_ALL_ORDER_AUDIT.md`](extended-geometry/TWO_PAIR_SIC_BIDEGREE44_RANK_TWO_ALL_ORDER_AUDIT.md).

On the rank-one boundary, the same checker proves that moments one through
six have exact Jacobian rank six and computes their nonnegative Hilbert
numerator, of coefficient sum \(50\). Rank-one annihilator tensors are
Hilbert--Mumford unstable: for the one-parameter subgroup adapted to the
annihilating root, the missing extremal coefficient makes every surviving
tensor weight strictly positive.

The third checker refines that rank-one boundary on the squarefree
operator chart. After \(\mu_1\) and the normalization \(e=1\), it proves
that the coefficient pivot of the \(a\)-linear second moment has empty
common boundary with \(\mu_2,\ldots,\mu_6\) over
\(\mathbb Q(\lambda)\). On the pivot-open chart the three expected
annihilator sections have ideal
\[
(8c-3d^2,\ d(d-4)(d-4\lambda)),
\]
and the last four moments contain the cubic factor with multiplicities
\(1,1,2,2\) after imposing the first relation. The checker then computes
the gcd of all pairwise \(c\)-resultants. Its only extra quadratic factor
is supported on the pivot divisor because \(p^3\) reduces to zero there,
and the other three branches force
\((c,d)=(0,0),(6,4),(6\lambda^2,4\lambda)\).
Thus the generic squarefree fiber is exactly the annihilator sections.
The same checker closes the pivot-annihilator orbit
\(\lambda^2+4\lambda+1=0\) and the equianharmonic orbit
\(\lambda^2-\lambda+1=0\): each projective fiber has degree four and
the expected four-point radical, with eighth-power certificates.  The
harmonic orbit is closed by the separate exact \(\lambda=2\) anchor.
It also extracts the complete pivot-boundary exceptional gcd and the
three expected-branch gcds. Their sole new \(S_3\)-orbit is represented by
\[
22\lambda^4-54\lambda^3+\lambda^2-54\lambda+22=0;
\]
the checker closes that quartic-field fiber with the same degree and
radical and with eighth-power certificates. On the chart
\(8c-3d^2=0,\ d(d-4)(d-4\lambda)\ne0\), direct substitution and exact
division by the invertible cubic powers reduce the problem to a
three-variable unit ideal over \(\mathbb Q\). One Rabinowitsch membership
remains: on \(p\ne0\), the \(8c-3d^2\ne0\) chart should be supported only
at \(\lambda=0,1\).  Writing \(q=\lambda^4(\lambda-1)^4\) and
\(M=p(8c-3d^2)\), the checker reduces this to the target-only membership
\(qM^5\in(f_3,f_4,f_5,f_6)\).  Modulo \(101,103,107\), the least
saturation exponent is consistently \(5\) and the degree-order basis has
size \(87\).  At \(101\), the four lifted multipliers have degree/term
profiles `(34,5356)`, `(29,3679)`, `(27,3037)`, and `(22,1853)`.
These are exact finite-field identities, but the checker does not promote
them to a characteristic-zero certificate.

The target-only lift experiments can be run separately:

```bash
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py --prime 101
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py \
  --prime-lift 101 5
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py --direct 5
```

The first two commands reproduce the exponent, basis size, and finite-field
lift profile.  The third asks for the exact rational target-only lift; the
recorded run exceeded its 1,200-second bound.  A timeout is not evidence
against membership.

The resumable large-prime CRT experiment is recorded in
`artifacts/generated-results/two_variable_quartic_squarefree_crt.json`.
To rebuild it independently, use a new checkpoint path:

```bash
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py \
  --crt-lift 5 /tmp/two-variable-quartic-squarefree-crt.json \
  1000003 1000033 1000037 1000039 1000081 1000099 1000117 1000121 \
  1000133 1000151 1000159 1000171 1000183 1000187 1000193
```

The first thirteen images have a common 14,508-term support.  Twelve
build primes give a 240-bit modulus and `1000183` is the holdout; only
three balanced rational reconstructions agree at the holdout.  The final
two primes have different supports.  This rejects coefficientwise CRT of
these arbitrary lifts at the recorded bound, not membership over
\(\mathbb Q\).  The two attempted normalizations are:

```bash
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py \
  --compare-tracked-lifts 5 101 103 107
.venv/bin/python \
  scripts/research_two_variable_quartic_squarefree_membership.py \
  --compare-syzygy-lifts 5 std component 101
```

The tracked transformations finish but have three different support
hashes.  The component-order syzygy normalization reached its 600-second
per-prime timeout at \(101\).

The fourth checker handles the remaining at-most-two-root normal forms
\(u^rv^{4-r}\). For \(r=0,4\), the first moment gives the one-sided
hyperplane. For \(r=1,2,3\), the first four moments have exactly the
expected two one-sided linear components, with eighth-power radical
certificates. Combined with the existing five-moment three-root theorem,
this closes every collided-root rank-one stratum. The single squarefree
uniform-specialization gate described above remains.

The dual-linear two-pair theorem is replayed by

```bash
python3 scripts/verify_dual_linear_sic2.py
```

This dependency-free audit accompanies the
[dual-linear `SIC(2)` theorem](extended-geometry/DUAL_LINEAR_SIC2.md).
For every \(p=w\mathbin{\cdot}H\), the first two contractions force
\(\operatorname{tr}JH=\det JH=0\), hence
\(H=c+(b,-a)f(ax+by)\).  If \(d=\deg f\) and
\(G=\deg_{x,y}g\), the proof gives
\(\mathcal E_2(gp^m)=0\) for \(m>(d+2)G\).  The normalized Keller case
needs only the first contraction and retains the sharper cutoff \(m>G\).
The checker verifies the second-contraction identity and replays both
cutoffs on exact integer examples.

The unrestricted bidegree-\((2,2)\) theorem is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_bidegree22_frontier.py
```

The checker proves that the natural four-parameter linear compression of
the three-pair witness is forced onto a strict one-sided weight branch by
its first two pure contractions.  Globally, it eliminates a four-parameter
finite-flag chart of the pair-linear one-sided nullcone and obtains twelve
generators for its prime ideal.  If \(I\) is generated by the first six pure
contractions and \(J\) is this nullcone ideal, exact reductions prove
\(I\subseteq J\), \(j_1\in I\), and \(j_r^7\in I\) for
\(2\le r\le12\).  Thus \(\sqrt I=J\), proving SIC(2) for every
bidegree-\((2,2)\) form, including dense eight- and nine-term forms.

As an independent sparse regression, the checker also enumerates all exact
supports of size at most seven, certifies the eight six-term and twelve
seven-term Laurent-curve charts, and verifies their hidden one-sided
factorizations.  The proof and exact claim boundary are in
[`TWO_PAIR_SIC_BIDEGREE22_FRONTIER.md`](extended-geometry/TWO_PAIR_SIC_BIDEGREE22_FRONTIER.md).

The first mixed-bidegree ordinary-degree search is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_mixed_22_13_31.py
```

For
\[
F=F_{2,2}+F_{1,3}+F_{3,1},
\]
the checker enumerates all \(35\) nonconstant bidegrees of ordinary degree
at most seven and all \(256\) primitive positive/negative central circuits,
then treats this complete three-block stratum exactly.  On
the nonzero \(F_{1,3}\) branch, the dual-linear normal form gives
\(F_{1,3}=\xi _2z_1^3\).  Polynomial-valued contractions through order
three force the aligned triangular flag.  The remaining scalar
contractions through order four reduce to three polynomials whose exact
lexicographic basis contains \(h^3\), and whose \(h=0\) boundary contains
fourth powers of the other two central variables.  The surviving support
has a strict two-step exponent cone, giving an explicit all-multiplier
cutoff rather than a bounded moment prefix.  The \(F_{1,3}=0\) branch
reduces to the bidegree-\((2,2)\) theorem.  The proof and the still-open
ordinary-degree-\(<8\) collection classes are in
[`TWO_PAIR_SIC_ORDINARY_DEGREE_LT8_MIXED_FRONTIER.md`](extended-geometry/TWO_PAIR_SIC_ORDINARY_DEGREE_LT8_MIXED_FRONTIER.md).
The same replay checks the central factorial formula for
\(V_{1,d}\oplus V_{d,1}\), \(2\leq d\leq6\); the written exponent-cone
proof is uniform for every \(d\geq2\).

The next dual-linear mixed branch is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_mixed_33_14_41.py
```

For
\[
F=F_{3,3}+F_{1,4}+F_{4,1},\qquad F_{1,4}\ne0,
\]
the checker normalizes \(F_{1,4}=\xi _2z_1^4\), verifies the maximal-weight
linear ladder, and uses two full-polynomial output coefficients from
moments three and four to force the \(F_{3,3}\) coefficient matrix upper
triangular.  It constructs the four-variable diagonal/opposite-corner
core and its moments through order five.  Their exact 28-element rational
Gröbner basis contains tenth powers of all three diagonal parameters and
the fifth power of the corner parameter.  The residual exponent cone
then gives the explicit multiplier cutoff.  The checker does not close
the \(F_{1,4}=0\) boundary: the written proof identifies that boundary
exactly with the existing balanced bidegree-\((3,3)\) problem.

The first degree-eight diagonal-core probe is

```bash
.venv/bin/python scripts/explore_two_pair_sic_mixed_diagonal_core.py
```

The script uses a direct factorial composition formula for the core of
\(V_{4,4}\oplus V_{1,5}\oplus V_{5,1}\), avoiding expanded Wick
polynomials.  Moments two through six have exact rational Jacobian rank
five.  Singular gives basis size \(132\), quotient dimension \(360\), and
origin-support power reductions over both \(\mathbb F_{101}\) and
\(\mathbb F_{1009}\).  These are exact finite-field computations but only
evidence for the rational radical.  The script neither reconstructs a
characteristic-zero certificate nor checks the preceding \(d=4\)
triangularization equations.

The first non-dual-linear positive block is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_mixed_23_32_pure_summands.py
```

The checker constructs
\(V_{2,3}=\operatorname{Sym}^5\oplus\operatorname{Sym}^3
\oplus\operatorname{Sym}^1\) exactly.  It identifies the first contraction
with the linear projection, the pure cubic second-moment ideal with the
rational normal cubic, and the pure quintic ideal with the prime
tangential variety \(L^4M\), including an exact Singular primary
decomposition.  It then inserts all three nonzero orbit normal forms into
\(V_{2,3}\oplus V_{3,2}\).  Moments two through four eliminate every
negative-block coefficient capable of an unbounded central contraction,
and explicit exponent bounds prove eventual mixed vanishing.  The branch
with both positive irreducible summands nonzero remains open.  At one
explicit mixed-positive point, the same checker constructs the full
moment-one-through-four Jacobian.  Its two transverse excess directions
are obstructed at deformation orders two and four, respectively.  For the
second direction it retains all six parameters in both higher correction
spaces and verifies that the two order-four compatibility quadratics have
nonzero resultant
\(2283980165392458318151680000\).  This is an exact local jet exclusion,
not a global classification of the remaining branch.

The one-parameter generic strengthening is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_mixed_23_32_generic_local_gate.py
```

The raising operator removes the third incidence coefficient on
\(b(a-b)\ne0\), leaving \(B(u)=VZ^2(uVY+WZ)\).  Over
\(\mathbb Q(u)\), the checker retains all second- and third-order
corrections and factors the fourth-order resultant as
\[
648u^8(u+6)^4(4u+3)^2
(2u^3+10u^2+21u+9)^2S(u),
\]
where the displayed artifact records the irreducible sextic \(S\).
The linear subresultant also gives the unique fourth-order correction
\(\tau=-B(u)/A(u)\) at every root of \(S\), since the checker verifies
\(\operatorname{Res}(A,S)\ne0\).
It also recomputes \(u=0,-6,-3/4,\infty\) and one \(a=b,c\ne0\)
representative without the rational-function chart; every resulting
fourth-order obstruction ideal is the unit ideal.  This is a generic
local first-four-moment theorem, not a global SIC classification.

The balanced cubic two-variable GVC theorem is replayed by

```bash
.venv/bin/python scripts/verify_two_variable_cubic_gvc.py
```

For the three nonzero binary-cubic symbol orbits \(u^3\), \(u^2v\), and
\(uv(u+v)\), the checker derives the apolar moments directly from
coefficient expansion.  It verifies that moments through orders one,
three, and four, respectively, leave only the displayed one-sided normal
forms.  In the squarefree orbit it checks all three exact branch
factorizations and the constant annihilating direction of every surviving
pure cube.  The written degree argument gives
\(\Lambda^m(QP^m)=0\) for \(m>\deg Q\); this all-order conclusion, not a
bounded replay, proves the theorem.  See
[`TWO_VARIABLE_CUBIC_GVC_THEOREM.md`](extended-geometry/TWO_VARIABLE_CUBIC_GVC_THEOREM.md).

The quartic three-root continuation is replayed by

```bash
.venv/bin/python scripts/verify_two_variable_quartic_three_root_gvc.py
.venv/bin/python scripts/verify_two_variable_quartic_squarefree_generic.py
```

After normalizing the \((2,1,1)\) symbol to \(u^2v(u+v)\), the first
moment eliminates one polynomial coefficient.  The checker derives moments
two through five, proves containment in the three-component one-sided
ideal, and verifies that the fourth powers of all five generators of that
radical lie in the moment ideal.  The same
[`low-root note`](extended-geometry/TWO_VARIABLE_LOW_ROOT_GVC_THEOREMS.md)
proves the all-degree theorem for symbols with at most two roots via the
one-variable Duistermaat--van der Kallen constant-term theorem; that part
is a written proof rather than a bounded computation.  The squarefree
checker verifies the four symbolic annihilator sections and proves that at
cross-ratio \(2\), moments one through six have exactly the four expected
reduced projective zeros.  Proper-family upper semicontinuity then proves
the same equality on a nonempty Zariski-open set of cross-ratios.  This is
a generic theorem; finitely many exceptional squarefree orbits remain
possible.

The later
[`SPLIT_SYMBOL_GVC_THEOREM.md`](extended-geometry/SPLIT_SYMBOL_GVC_THEOREM.md)
requires no computer algebra.  It factors a homogeneous operator symbol
as directional derivatives and proves
\[
 \Lambda^m(P^m)
 =(m!)^d\operatorname{CT}
 \left(\frac{P(z+\sum_i t_i v_i)}{\prod_i t_i}\right)^m.
\]
Choose one generic \(z\) exposing the full finite \(t\)-support.
Duistermaat--van der Kallen gives an integral weight separating that
support from the origin, and the same weight works at every specialization
of \(z\).  A fixed translated multiplier cannot cross the linearly growing
gap.  Thus every homogeneous binary operator satisfies GVC for arbitrary
\(P\), including \(\deg P>\operatorname{ord}\Lambda\).  There is no
bounded replay to confuse with the proof.

The nonhomogeneous lowest-order extension and the rank obstruction to
natural conversions of the two-pair witness are checked by

```bash
python3 scripts/verify_separable_gvc_escape_obstructions.py
.venv/bin/python scripts/verify_binary_heat_quadratic_gvc.py
```

The first dependency-free command also rewrites
`artifacts/generated-results/separable_gvc_escape_obstructions.json`; its
current whole-file SHA-256 is
`3343e46cca1b9459f0a3f113278d1db610379e1c7083370290f30f32e420f226`.

The finite checker verifies that the witness matrix has determinant \(48\)
and rank five, that four separated rank-one channels cannot reach it, and
that coefficient extraction is not multiplicative.  The written proof in
[`SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md`](extended-geometry/SEPARABLE_GVC_ESCAPE_OBSTRUCTIONS.md)
shows that if \(r\) is the lowest positive order of a nonhomogeneous binary
operator, then \(\deg P\leq r\) implies GVC.  Higher operator pieces can
enter a mixed value only a bounded number of times, and the split-symbol
Newton separator absorbs those bounded defects.  This is an all-order
argument, not a bounded search.  The same note records the arbitrary-degree
factor-unit extension \(\Lambda=\Lambda_0\Gamma\), with \(\Lambda_0\)
homogeneous split and \(\Gamma(0)\ne0\).  For a completely general
nonhomogeneous binary operator, every fixed number of leading homogeneous
mixed layers still vanishes eventually; any defect must move to unbounded
filtration depth.  The note also proves that fixed linear translation plus
one diagonal coefficient is possible exactly for split homogeneous
symbols.  Finally, it closes the binary linear-plus-quadratic heat class
for \(\deg P\le2\): after normalizing the linear part, the only irreducible
coefficient is \(C\partial_y^2\), the first equation gives
\(P=c(y^2-2Cx)+ey+f\), and the second is \(16C^2c^2\).  The surviving
affine transverse form has cutoff
\(m>2\deg_xQ+\deg_yQ\).  The second command derives the first-moment
normal form and checks the universal second-moment identity symbolically.
It also constructs the generic degree-six heat-harmonic polynomial and
checks
\(\Lambda^2(P^2)=4C^2(P_{yy})^2\).  The written product-rule proof gives
the corresponding highest-\(y\)-degree square in every degree and closes
every binary operator with nonzero linear part and no terms above order
two.  The checker also verifies the generic cubic drift--diffusion normal
form.  The cutoff is the written derivative count.
The same checker tests the degree-six family for
\(\partial_x+C\partial_y^2+E\partial_y^3\).  The written iterated
product-defect proof closes every
\(\partial_x+h(\partial_y)\), with cutoff
\(m>r\deg_xQ+\deg_yQ\) when \(r\) is the lowest order of \(h\).
The proof permits formal \(h\), since its differential action is locally
finite.  Formal Weierstrass division writes every binary symbol with
nonzero linear part as \(U(\xi,\eta)(\xi+q(\eta))\); the unit
\(U(\partial)\) and its inverse are locally finite differential
automorphisms.  This reduces the entire lowest-order-one class to the
separated theorem for arbitrary \(P\).
The symbolic checker also adds a completely general cubic operator piece
and verifies that the quadratic-\(P\) second moment remains
\(16C^2c^2\).  The written order argument extends this to arbitrary
higher pieces, proving safety for every \(\deg P\le2\) when the linear
operator part is nonzero.
For a general operator and cubic polynomial, the checker solves the first
equation and retains every order-four and order-five term that can occur
in the second moment.  It verifies the decisive coefficients
\[
144C^2p_9^2,\qquad 16C^2p_5^2,\qquad 648G^2p_9^2.
\]
Higher operator orders kill \(P^2\), and the retained terms cannot change
these branches.  This gives the written theorem for every binary operator
with nonzero linear part and every \(P\) of degree at most three.
The checker additionally retains the complete operator \(7\)-jet for a
quartic \(P\) and verifies the successive branch coefficients
\(2304C^4p_4^2\), \(15552G^2p_4^2\), and \(39168L^2p_4^2\).
These finite-jet calculations are regressions for the stronger formal
straightening theorem, rather than the source of its all-degree proof.
For lowest positive order two and cubic \(P\), it separately checks both
quadratic-symbol orbits.  In the double-line orbit it verifies the unique
second-moment cancellation and the decisive third moment
\(-4608C^3p_{xy^2}^3\); in the distinct-line orbit it verifies the
triangular second-moment branches.  The written strict weighted-degree
cutoffs turn the surviving branches into an all-order theorem.
For quartic \(P\), the checker retains every operator jet that can enter
the displayed moments.  It verifies the triangular distinct-root
second-moment closure, the double-line \(xy^3\) third-moment coefficient
\(-3604176H^3\), and the two residual polynomials \(S,T\) on the
double-line \(y^4\) branch.  Their exact resultant is a monomial times the
homogeneous sextic printed in Proposition 3.8, leaving only finitely many
nonweighted parameter ratios for the next moment.
It then evaluates moment four, eliminates the remaining coefficient, and
verifies that the resulting octavic has gcd one with the sextic.  The
written \((2,1)\)-weighted-face identity shows that higher operator jets
and lower polynomial terms cannot alter any pure moment, upgrading this
calculation to closure of the full \(r=2,\deg P\le4\) cell.
For \(r=3,\deg P=4\) with triple-root leading symbol, the checker verifies
the three stabilizer types.  In the \(x^2y^2\) type it checks the forced
weighted correction chain and terminal coefficient
\(3361505280U^4\); the \(xy^3\) and \(y^4\) types terminate at moments
three and two.  The written weighted-face argument supplies the all-order
mixed cutoffs.  It also checks all double-root branches and the squarefree
orbit.  For the latter it computes the leading-moment Gröbner basis,
reduces the three fourth-power tips by root permutation, retains the full
order-four and order-five jets at the \(x^4\) tip, and verifies the
terminal weighted-face value \(129392640T^3\).  Together with the written
orbit and weight arguments this proves GVC for every binary operator and
every \(P\) of degree at most four.  Finally, the checker verifies the
leading reduction in the first open \(r=2,\deg P=5\) row: the only three
top pairs are
\((\partial_x\partial_y,x^5)\),
\((\partial_x^2,xy^4)\), and
\((\partial_x^2,y^5)\).  It then retains the complete jets that can enter
the first two moments.  The distinct-root branch reduces to
\(P=f(x)+ay\) and
\(\Lambda=\partial_y\Gamma+H(\partial_x)\), with
\(\operatorname{ord}_{\min}H\ge6>\deg f\).  The written binomial
derivative count gives its all-order mixed cutoff.  For the \(xy^4\)
branch the checker verifies two nested third-moment faces, terminating in
\(-553153536H^3\) and \(-5430509568J^3\).  For the \(y^5\) branch it
checks the six weighted second-moment ratios and the four nonzero
third-moment residuals; the other two ratios are one-sided.  This closes
the full \(r=2,\deg P\le5\) cell.  The checker finally derives the next
\(r=3,\deg P=5\) leading reduction: four triple-root, three double-root,
and one squarefree top-form normal forms.  This checker supplies the
leading-face classification only; the separate degree-five frontier
checker below closes all eight nonhomogeneous correction systems.
The accompanying written no-go for formal umbral straightening is
proof-theoretic: conjugation by an algebra automorphism preserves the
Leibniz rule, while a locally finite formal constant-coefficient operator
is a derivation only when its symbol is linear.  No bounded computation is
used for that statement.

The eight cubic-leading quintic normal forms and the squarefree
quartic-leading cross-ratio row are closed by
[`BINARY_DEGREE_FIVE_GVC_FRONTIER.md`](extended-geometry/BINARY_DEGREE_FIVE_GVC_FRONTIER.md).
The default exact checker replays the triangular moment eliminations,
terminal one-sided faces, and explicit mixed-multiplier tails:

```bash
make verify-binary-degree-five-gvc
```

For an exploratory dump of all eight complete second-moment jets, with the
optional product-splitting and unit-pivot heuristic, run:

```bash
.venv/bin/python scripts/explore_binary_degree_five_gvc_frontier.py \
  --triangular-components
```

With Singular, replay both residual radicals and the uniform
squarefree-quartic projective saturation:

```bash
.venv/bin/python scripts/verify_binary_degree_five_gvc_frontier.py \
  --singular --singular-top
```

The separate modular screen enumerates \(6{,}696{,}142\) residual triples
over \(p=101,103,107\) and \(2{,}082{,}612\) squarefree-quartic projective
top forms over seven smaller primes:

```bash
make search-binary-degree-five-gvc
```

It regenerates
[`binary_degree_five_gvc_face_search.json`](artifacts/generated-results/binary_degree_five_gvc_face_search.json).
The modular record is an exhaustive bounded experiment on the displayed
faces; the characteristic-zero radical computations and weighted
face-separation argument are the proof.  The non-squarefree
quartic-leading nonhomogeneous rows remain outside this particular theorem,
so it alone does not give a universal binary degree-five corollary.

The quadruple-root partition \((4)\) is closed separately, with arbitrary
lower polynomial pieces and arbitrary higher Weierstrass operator jets, by:

```bash
.venv/bin/python scripts/verify_binary_quartic_quadruple_root_gvc.py
```

This exact checker replays the defect-one radical, all three minimal
branches, the terminal weight-eight equality chain through pure moment
five, and the final strict or one-sided weight separators.  The proof and
normalizations are documented in
[`BINARY_QUARTIC_QUADRUPLE_ROOT_GVC.md`](extended-geometry/BINARY_QUARTIC_QUADRUPLE_ROOT_GVC.md).

The other repeated-root quartic partitions are closed by:

```bash
.venv/bin/python scripts/verify_binary_quartic_triple_simple_root_gvc.py
.venv/bin/python scripts/verify_binary_quartic_double_root_gvc.py
```

The first command verifies the \((3+1)\) defect-one radical, every
projective branch through defect three, and its final separators.  The
second closes \((2+2)\) and \((2+1+1)\), including both two-parameter
weight-six equality systems and the isolated fifth-power components.
Together with the earlier degree-five checker, these exact calculations
prove binary GVC for every polynomial of degree at most five.  See
[`BINARY_QUARTIC_TRIPLE_SIMPLE_ROOT_GVC.md`](extended-geometry/BINARY_QUARTIC_TRIPLE_SIMPLE_ROOT_GVC.md)
and
[`BINARY_QUARTIC_DOUBLE_ROOT_GVC.md`](extended-geometry/BINARY_QUARTIC_DOUBLE_ROOT_GVC.md).

The first sextic frontier cell, with lowest symbol of order five and root
partition \((5)\), is closed by:

```bash
.venv/bin/python scripts/verify_binary_quintic_quintuple_root_gvc.py
```

This exact checker verifies the defect-one radical, all five projective
top-form branches, and the complete terminal equality chain through
operator order ten and pure moment six.  The proof includes arbitrary
lower pieces of \(P\) and arbitrary higher operator jets.  See
[`BINARY_QUINTIC_QUINTUPLE_ROOT_GVC.md`](extended-geometry/BINARY_QUINTIC_QUINTUPLE_ROOT_GVC.md).
The complete replay below incorporates this longest component and closes
the other quintic root partitions.

All six remaining quintic root partitions, and hence the complete
\((r,\deg P)=(5,6)\) row, are closed by:

```bash
.venv/bin/python scripts/verify_binary_quintic_all_root_partitions_gvc.py
```

This checker verifies the Hall-matching classification of the leading
pure-zero locus, the local correction systems for root multiplicities one
through five, generic strict quintic cofactors, and every final weighted
separator.  It invokes the quintuple-root replay for the multiplicity-five
cell.  See
[`BINARY_QUINTIC_ALL_ROOT_PARTITIONS_GVC.md`](extended-geometry/BINARY_QUINTIC_ALL_ROOT_PARTITIONS_GVC.md).
The complete \((r,\deg P)=(4,6)\) row is closed by:

```bash
.venv/bin/python scripts/verify_binary_quartic_all_root_partitions_gvc.py
```

This exact checker verifies the Hall leading-locus classification, every
repeated-root terminal face, both coupled pure-sixth-power endpoint
radicals, the triple-root terminal tail coefficients, and the simple-root
defect layers through defect four.  The finite-tail inequalities covering
arbitrary later operator jets are documented in
[`BINARY_QUARTIC_ALL_ROOT_PARTITIONS_GVC.md`](extended-geometry/BINARY_QUARTIC_ALL_ROOT_PARTITIONS_GVC.md).
The complete \((r,\deg P)=(3,6)\) row is closed by:

```bash
.venv/bin/python scripts/verify_binary_cubic_all_root_partitions_gvc.py
```

This exact checker verifies the cubic Hall locus, the triple-, double-,
and simple-root Newton ladders, and the two coupled weighted-face chart
covers.  It requires Singular and `msolve`; the latter is used over
characteristic zero to prove explicit affine saturations empty.  See
[`BINARY_CUBIC_ALL_ROOT_PARTITIONS_GVC.md`](extended-geometry/BINARY_CUBIC_ALL_ROOT_PARTITIONS_GVC.md).
The complete \((r,\deg P)=(2,6)\) row, and hence universal binary GVC
through polynomial degree six, is closed by:

```bash
.venv/bin/python scripts/verify_binary_quadratic_all_root_partitions_gvc.py
```

This exact checker verifies the quadratic Hall locus, the full
distinct-root first-equation reduction and second-moment ladder, all
half-integral and integral double-line Newton faces, and every primary and
secondary radical at the pure-sixth-power endpoint.  It requires Singular
over characteristic zero.  The arbitrary-jet weight-defect argument is in
[`BINARY_QUADRATIC_ALL_ROOT_PARTITIONS_GVC.md`](extended-geometry/BINARY_QUADRATIC_ALL_ROOT_PARTITIONS_GVC.md).
The quadratic-leading row on the next polynomial-degree-seven frontier is
closed by:

```bash
.venv/bin/python scripts/verify_binary_quadratic_septic_gvc.py
```

This exact checker verifies the seven-factor Hall locus, the complete
distinct-root first-equation reduction and second-moment ladder, every
half-integral and integral double-line face, and the full branch tree over
the pure seventh-power endpoint.  In particular it checks the extra
slope-three axis migration which does not occur in degree six.  The written
proof in
[`BINARY_QUADRATIC_SEPTIC_GVC.md`](extended-geometry/BINARY_QUADRATIC_SEPTIC_GVC.md)
uses the final common-threshold coordinate deficits to cover arbitrary
higher jets and a fixed multiplier.  A degree-seven counterexample must
therefore have lowest positive operator order three through six.

The cubic-leading row on the same frontier is closed by:

```bash
.venv/bin/python scripts/verify_binary_cubic_septic_gvc.py
```

This exact checker derives every crossing from the cubic two-wing normal
form, verifies the four exceptional monomial radicals and every origin
radical in their child intervals over \(\mathbb Q\), and audits all terminal
common thresholds.  It requires Singular.  The written proof in
[`BINARY_CUBIC_SEPTIC_GVC.md`](extended-geometry/BINARY_CUBIC_SEPTIC_GVC.md)
uses the uniform terminal-face theorem to absorb arbitrary strict jets and
the fixed multiplier.  A degree-seven counterexample must therefore have
lowest positive operator order four through six.

The remaining three septic rows, and hence binary GVC through polynomial
degree seven, are closed by:

```bash
.venv/bin/python scripts/verify_binary_high_order_septic_gvc.py
```

This exact checker constructs all 46 Hall charts for lowest orders four
through six, verifies 287 initial face radicals and 98 child face radicals
over (mathbb Q), checks the fifteen squarefree axis exceptions, and
audits strict marked-gap descent to every final common threshold.  It
requires Singular.  The proof and complete census are in
[`BINARY_DEGREE_SEVEN_GVC.md`](extended-geometry/BINARY_DEGREE_SEVEN_GVC.md).

Unrestricted binary GVC is proved without a degree census in
[`BINARY_GVC_ENVELOPE_CLOSURE.md`](extended-geometry/BINARY_GVC_ENVELOPE_CLOSURE.md).
The proof combines Hall localization, the shifted-ray endpoint theorem, and
the unequal common-threshold theorem.  It has no new computational premise:
the finite global lower and upper Newton envelopes cannot exchange horizontal
order while their weight gap is positive, and finite support forces that gap
to reach zero.  The first octic Ferrers face which exposed this argument has
the exact regression:

```bash
.venv/bin/python scripts/verify_binary_gvc_ferrers_regression.py
```

It uses `msolve` over characteristic zero to verify the radical
((A,S,T,BP,BQ,CQ)).  The first degree-nine gap-four staircase is a longer
optional replay:

```bash
.venv/bin/python scripts/verify_binary_gvc_ferrers_regression.py --gap-four
```

That optional command checks eight affine saturations and can take several
minutes.  These Ferrers calculations are regressions, not dependencies of
the unrestricted envelope proof.

The all-degree Hall localization and unequal-weight terminal-face theorem
have a dependency-free exact regression:

```bash
.venv/bin/python scripts/verify_binary_gvc_uniform_face_termination.py
```

It exhausts the Hall inequality through order twelve, checks the
coefficient-independent prime-valuation inequalities on small weighted
lattice segments, and verifies an exact weight-\((3,2)\) example whose
first moment cancels but whose fifth moment has the uniquely predicted
5-adic endpoint, together with the generic-translation multifactorial
identity underlying the Newton-intersection criterion.  It also checks the
homogeneous beta-integral identity, the failure of constant-term
extraction to commute with powers in the first two-channel example, and
the exponent-dependent factorial distortion under a toric blow-up.  It
also verifies the exact minimal Bernstein-circuit formula
\(\Phi(F^m)=(ac+bd)^m/(m+1)\), which identifies Long's rank-one
beta--torus circuit as a linear Hall annihilator rather than a GVC
counterexample.  Finally it replays the all-degree primitive cusp
parallelogram obstruction
\[
 E_{r,s}=T_r-T_s+\frac92(C_s-C_r)(C_r+C_s-2)
\]
and the moment-two closure of every sparse four-channel dilation.  The
same checker verifies the five-channel warning
\[
 A=X^2-\frac13Y^3,\qquad
 P=x^2+y^3+\frac{13}{30}+\frac{11}{2}x+xy^3:
\]
its first three scalar moments vanish, no four-channel restriction
inherits those three vanishings, and its fourth moment is \(1\,205\,760\).
This is a finite-prefix obstruction, not a pure-moment-zero pair.  The
checker also replays the all-even unit-line half-bridge theorem.  After
moments two and three determine \(u=de\) and \(v=ce^2\), it checks
\[
 M_4=\frac12\left(
 2Q_n-40C_nT_n+48T_n+81C_n^3-180C_n^2+132C_n-48
 \right),
\]
the exceptional value \(H_2=-480\), and the increasing multinomial ratio
\(Q_n/(C_nT_n)>20\) from \(n=4\) onward.  Finally it constructs the six
quadratic--cubic and two double-quadratic obstruction formulas from the
factorial weights \(W_{m,k}\), checks the double-quadratic determinants,
and finds all eight obstructions nonzero for every unequal
\(1\leq r,s\leq30\).  This last window is an exact regression, not the
all-\((r,s)\) inequality proof for the five expressions left open below.

The return classification and the first unbounded arithmetic certificate
are replayed by

```bash
.venv/bin/python scripts/research_binary_gvc_eight_obstructions.py --prove-h00
.venv/bin/python scripts/research_binary_gvc_eight_obstructions.py --prove-negative-corners
.venv/bin/python scripts/research_binary_gvc_eight_obstructions.py --prove-three-more
.venv/bin/python scripts/research_binary_gvc_eight_obstructions.py --verify-opposite-packet --limit 40
.venv/bin/python scripts/research_binary_gvc_frobenius_carry.py
.venv/bin/python scripts/research_binary_gvc_frobenius_carry.py \
  --radial-limit 3 --order-limit 2 --bridge-limit 29 --residue-limit 2
.venv/bin/python scripts/research_binary_gvc_ghost_shell.py
.venv/bin/python scripts/research_binary_gvc_quotient_graver.py
.venv/bin/python scripts/research_binary_gvc_witt_rees.py
.venv/bin/python scripts/research_binary_gvc_equal_radial_union.py
.venv/bin/python scripts/research_binary_gvc_reversal_union.py --prove-width-two
.venv/bin/python scripts/research_binary_gvc_torsion_torus_trace.py
.venv/bin/python scripts/verify_binary_gvc_torsion_torus_digit_separation.py
.venv/bin/python scripts/verify_binary_gvc_weighted_trace_obstruction.py
.venv/bin/python scripts/verify_binary_gvc_first_ghost_source_collapse_and_ray_rigidity.py
.venv/bin/python scripts/research_binary_gvc_eight_obstructions.py --limit 100
```

The second Frobenius-carry command is the focused replay for the affine
singleton-localization theorem.  Its final blocks verify endpoint-fibre
signed-digit stabilization and factorial units at two successive primes,
replay the sharp interior family
\[
 ((p-1)/2,1,(p-1)/2)-(0,p,0)=((p-1)/2)(1,-2,1)
\]
through every odd prime at most \(29\), and verify the radial-carry Hasse
compression formula.  The finite checks are regressions for the written
all-prime Graver, exposed-vertex, and Wilson-unit proofs.  The same note's
global Hilbert-module theorem is a direct Dickson--Gordan argument: all
\((p,y)\) fibres are a finite module over the Hilbert basis of
\(B(p,y)=0\).  It justifies one fixed toric/Graver presentation per support;
the computation is not the source of that finiteness statement.

The ghost-shell command also constructs the complete primitive three-level
affine ghost
\[
 G_{p;u,v}(X)=\sum_{t=1}^{\lfloor p/(u+v)\rfloor}
 \frac{(p-1)!X^t}{(ut)!(vt)!(p-(u+v)t)!}.
\]
For the 12 coprime types \(1\leq u\leq v\leq6\), it tests the primes through
43, the reduced rational window with numerator of absolute value at most 40
and denominator at most 20, every cyclotomic order through 80, and all 2,139
primitive irreducible polynomials of degrees two and three with coefficient
height at most four.  The only non-support survivor is the centered value
\((u,v,X)=(1,1,1)\).  The note proves this for everywhere-good,
Galois-stable root-of-unity candidates at every width and torsion order by
the first-admissible-prime and Euler-\(\varphi\)-degree argument.  It does not
remove an arbitrary extra finite exceptional-prime set.  The finite search is still not an
all-prime classification of arbitrary non-torsion algebraic roots; the
constant-term, Wilson-coefficient, and cyclotomic-separation arguments printed
in the note are the exact theorems.

The first command checks the two closed determinant formulas and the exact
coefficient-positive proof that
\(\mathcal H_{0,0}(r,s)>0\) for every unequal endpoint pair.  Its forward
difference certificate has 932 nonnegative terms and its \(s=1,\ r\ge2\)
base certificate has 19.  This is an all-order symbolic certificate.  The
second command proves \(\mathcal H_{0,3}<0\) and
\(\mathcal D_{0,2}<0\) for every unequal endpoint pair.  It verifies the
three exact step-ratio bounds defining the coupled cone and four
coefficient-positive numerator expansions, of sizes 266, 361, 2236, and
2236.  The third command proves the uniform nonvanishing of
\(\mathcal H_{0,1},\mathcal H_{0,2},\mathcal H_{1,1}\), using monotone
ordered-tail cones, fixed-ray cones, and exact finite complements.  It
then closes \(\mathcal H_{1,0},\mathcal D_{0,1}\) on the last wedge
\(r>s\ge4\), using the increasing product \(L_nM_n/C_n^3\) and
coefficient-positive expansions of 408 and 1692 terms.  The fourth
command verifies the exact second-coefficient identity for the opposite
three-by-three packet and checks its central-binomial ratio dichotomy
through degree and endpoint order 40.  This is a bounded regression; the
unbounded proof is the strict Vandermonde supermultiplicativity argument
in Theorem 7.5 of the canonical note.  The fifth command checks the
binary Frobenius carry gap and its normalized unit formulas on 28,858
exact homogeneous return types, checks the nonhomogeneous jet--carry
score on 477 mixed-degree types, and constructs 111,930 one-sided triples and
5,787,067 two-sided pair-pairs through radial degree 40.  The unbounded
proof in Lemma 7.4 ter and Corollary 7.4 sexies is the corresponding
Legendre--Kummer and finite-field kernel calculation.  The sixth command
constructs the centered-triple and two-by-two ghost diagonal blocks.  It
verifies their universal factors \(X-1\) and \(X(X+1)\), respectively,
and displays their prime-dependent residual factors through prime 43.
It constructs the exact characteristic-zero beta diagonals and verifies
their common gcd \(X(X+1)(X^2+X+1)\) across the tested primes.  Its rational
cross-prime search of height 20 leaves only \(1\) in the centered block
and \(0,-1\) in the beta block, but does not see the algebraic cube-root
branch; that window is only a regression.  Proposition 7.4 quater proves
the universal and persistent factors for all primes at least five.  The
same command verifies
the terminal augmented blocks: the beta ordinary row \(1+X\) leaves only
the Hall value, and the centered Bessel endpoint rows \(U,U^2+2V\) have
Jacobian determinant \(2\) and force support loss.  Corollary 7.4
quinquies closes the isolated atom arithmetic; compatibility with the
common high-digit quotient remains unproved inside the parked Hall/carry
route.  The Hall-envelope proof does not need it.  The seventh command verifies
the first obstruction to circuit-only quotient peeling,
\[
 R_3B_1B_2=R_0B_3^2.
\]
On the projected support \(R=\{0,3\}\), \(B=\{1,2,3\}\), these are the
only two states of their color-count/level fiber, and their support-five
difference has no circuit move.  The checker verifies primitivity and
finds this as the first such two-color identity, up to reversal.  It also
checks the two terminal completions through \(R_2B_2^2\) and
\(R_3B_0B_3\).  Proposition 7.4 septies records the projected-scroll
obstruction and the finite Gröbner bound; Corollary 7.4 octies closes
this first block by the circuit-completion/radial-reversal dichotomy.
The checker also replays the explicit \(S(6)\) and \(S(5,4)\)
non-universal-Gröbner witnesses from
Bogart--Hemmecke--Petrović.  The projected \(S(5,4)\) fiber has exactly
four states, circuit-component sizes \(1+3\), and repeated-ray
factorial signatures \((2,2),(2,1,1),(1,1,1,1)\) with Stirling bases
\(16,4,1\).  Proposition 7.4 nonies proves that every whole exposed
scroll profile is terminal by a one-variable constant-term reduction.
The eighth command tests the first two normalized \(p\)-typical Witt
coordinates at \(p=13\) on two unit coefficient specializations of the
support-five, \(S(6)\), \(S(5,4)\), and first larger
reversal-symmetric packets.  In every case the first residual has
valuation exactly one, while both \(\mathcal G_2-\mathcal G_1\) and the
second Witt residual have valuation exactly two.  Proposition 7.4
decies proves the all-height result: the Laurent constant-term factor
and the signed normalized radial factorial are Gauss sequences, so
every exposed packet has an integral \(p\)-typical Witt recursion.
Ghost injectivity does not split a vanishing sum of several profile
Witt vectors; the remaining open step is a profile-separating Rees
initial idempotent or a separator/support-loss consequence of
least-profile cancellation.  The ninth command verifies the complete
equal-radial union identity directly through scale four.  It also finds
the smallest persistent failure of color-count saturation:
\(R=\{0,2\}\), \(B=\{1\}\), whose achievable red counts at scale \(N\)
are \(0,2,\ldots,2N\).  Despite those holes, the whole union is one
Laurent constant-term sequence; its first two symbolic rows reduce to
\(-14a^4\).  Proposition 7.4 undecies proves that color-count
saturation is unnecessary once all states at one oriented radial
vector are exposed.  It does not prove Hall/jet exposure of that
complete union or separate a coordinate-reversed pair.  The tenth
command saturates the four characteristic-zero endpoint charts for the
first coordinate-reversed Laurent width.  For support in \([-2,2]\)
and target slopes \(\pm1\), the charts
\((-1,1),(-1,2),(-2,1),(-2,2)\) close at rows \(2,4,4,8\).
It also exhausts the projective five-coefficient space modulo
\(5,7,11\), with no survivor through row eight.  The rational
saturations prove the width-two statement; the modular census is only
a regression.  Arbitrary reversal width is closed by the finite-trace
digit theorem below.  The eleventh
command verifies the exact regular-representation identity
\[
 \operatorname{CT}_{\mathbb Z\times C_q}(u^N)
 =q^{-1}\operatorname{CT}_{\mathbb Z}
   \operatorname{Tr}(\operatorname{Reg}_{C_q}(u)^N),
\]
its log-determinant generating series, and the compatible Frobenius
identity for \(p\equiv1\pmod q\), on exact \(C_2\) and \(C_3\)
examples.  It also diagonalizes the width-two reversal packet into its
two \(C_2\)-character components.  The twelfth command verifies signed
base-\(p\) digit uniqueness directly in one and two free Laurent
variables:
\[
 \operatorname {CT}(f^{n_0+n_1p+\cdots+n_sp^s})
 \equiv
 \prod_j\operatorname {CT}(f^{n_j})^{p^j}\pmod p.
\]
It also replays the Newton-identity endpoint for two, three, and five
trace components.  The characteristic-zero proof uses arbitrarily large
completely split primes: repeated equal digits recover every power sum
of the component moments, forcing componentwise vanishing.  Together
with character orthogonality and Duistermaat--van der Kallen, this closes
the identity-coefficient torsion--torus trace lemma and every
scale-compatible carry packet.  It does not prove that a Hall--jet shell
has that form.  The thirteenth command verifies the weighted mixed-digit
extension and two exact obstructions to that promotion.  The dilation
pair \(z+z^{-1}\), \(z^2+z^{-2}\) has equal constant-term power
sequences, so affine \(C_2\) character weights cancel at every pure row
while a fixed multiplier detects every odd row.  It also checks
\[
 v_{11}\!\left(\mathcal L((y^2+4xy+2x^2)^{12})\right)=3,
\]
where naive Laurent repeated-digit factorization predicts a
valuation-two nonzero residue.  The fourteenth command is
retained as an exact sign regression through 100; all eight all-order
obstruction conclusions are supplied independently by the first three
proof modes.  The accompanying
transverse-lattice proof shows that every order-four return in the 14
reduced types is generated by the primitive order-two and order-three
rows, so there is no separate primitive order-four branch to search.

The bounded counterexample probes which led to the repeated-digit theorem
can be replayed by

```bash
.venv/bin/python scripts/search_binary_gvc_torsion_torus_counterexample.py \
  pair --width 3 --height 1 --depth 14
.venv/bin/python scripts/search_binary_gvc_torsion_torus_counterexample.py \
  shared --width 4 --height 1 --depth 16
.venv/bin/python scripts/search_binary_gvc_translation_isoperiodic_twists.py \
  --max-degree 12 --extra-depth 4
.venv/bin/python scripts/search_binary_gvc_translation_isoperiodic_twists.py \
  --max-degree 3 --extra-depth 4 --rectangles
python3 scripts/verify_binary_gvc_translation_tangent_rigidity.py
python3 scripts/verify_positive_return_semigroup_jet_rigidity.py \
  --require-singular
.venv/bin/python scripts/search_binary_gvc_translation_isoperiodic_twists.py \
  --character-order 3 --max-degree 8 --extra-depth 4
.venv/bin/python scripts/search_binary_gvc_translation_isoperiodic_twists.py \
  --character-order 4 --max-degree 7 --extra-depth 4
```

They find no survivor.  Their longest zero prefixes are sparse binomial
lattice-delay examples, first obstructed at rows 10 and 14 respectively.
The third command tests character twists of one translated binomial row.
After quotienting scalar sign, the torus action \(t\mapsto-t\), and
coefficient reversal, it checks 8,188 twists and 81,924 exact
twist--slope rows through depth \(2d+4\), for every \(d\leq12\).  Every
collision is explained by one of those scale-compatible symmetries.
The fourth command adds all \(C_2\) twists on the \((2,2)\) and \((3,2)\)
binomial Taylor rectangles.  Its 2,304 twists and 4,352 exact moving rows
have no collision outside scalar, two-torus, reversal, and
coordinate-exchange symmetry through depths 12 and 14.

The fifth command replays the proved primitive translation-tangent theorem.
For every primitive \((d,r)\) through degree 12 it verifies that the matrix
\[
 \binom dj\binom{d(N-1)}{rN-j}
\]
has the one-dimensional kernel \(j-r\); it also checks the nonprimitive
power/subsequence rank jump, finite-field ranks away from the chosen integer
minors through prime 97, all 2,550 displayed spanning generators of the
universal blind tangent module on 225 rectangular slope cases through
bidegree \((6,6)\), the exact two-dimensional quadratic-Hessian kernel on
100 slope cases through bidegree \((5,5)\), and the factorially weighted
two-free-translate counterexample to module-only inheritance.  The theorems
and the two-direction no-go are unbounded and are proved in
[`BINARY_GVC_TRANSLATION_TANGENT_RIGIDITY.md`](extended-geometry/BINARY_GVC_TRANSLATION_TANGENT_RIGIDITY.md).
It implies flatness of every \(q^a\)-order character collision at primitive
one-direction slope once the underlying prime \(q\) is sufficiently large.
The bilinear ghost proves that two-dimensional separation cannot follow from
the first cyclotomic neighbourhood; its first potentially effective row is
quadratic or higher.  The positive return-word proof then shows that this
quadratic row separates every fixed nonflat integer label and does so modulo
all but finitely many primes.

The sixth command replays the broader
[positive return-semigroup jet theorem](extended-geometry/POSITIVE_RETURN_SEMIGROUP_JET_RIGIDITY.md).
Its dependency-free part checks the positive-return group-completion
mechanism on 27 bounded Cartesian configurations.  Singular then constructs
the centered \((2,2)\) derivative ideals from their return words: their
affine dimensions through jet orders one, two, three, and four are
\(4,2,1,0\), and the final quotient has dimension \(40\) with an
eleven-element standard basis.  This is an exact replay of the finite
four-jet certificate and was replayed with Singular 4.3.2.  The
all-configuration full-jet theorem and its
all-torsion-order corollary are proved in the note, not inferred from
the bounded computation.  The same note proves that independently marked
return polynomials at finitely many total degrees separate arbitrary paired
coefficient points modulo the coefficient torus.

The final two commands extend the earlier sign search using exact arithmetic
in \(\mathbb Z[\zeta_3]\) and \(\mathbb Z[i]\).  They test 63,972 and
123,792 moving rows and find no unexplained collision.  These \(C_3,C_4\)
searches, like the rectangular search, are bounded evidence for whether the
signed Hall shell inherits enough common marks before promotion.

The exact primitive translation-orbit census is

```bash
python3 scripts/research_binary_gvc_translation_observability.py \
  --radial-degree 4 \
  --output artifacts/generated-results/binary_gvc_translation_observability_span4.json
python3 scripts/research_binary_gvc_translation_observability.py \
  --radial-degree 5 \
  --modes operator,polynomial \
  --output artifacts/generated-results/binary_gvc_translation_observability_span5_one_colour.json
```

Normaliz 3.10.2 computes the projected Graver bases and Singular 4.3.2
computes the characteristic-zero derivative-orbit ideals after torus
saturation.  At span four, 48 of the 65 normalized mixed packets are already
factorially obstructed and all 17 survivors have empty one-colour and
independent-translation torus ideals.  Two survive only in the weaker
common-diagonal mode; they are the quartic Veronese identities displayed in
[`BINARY_GVC_PRIMITIVE_TRANSLATION_OBSERVABILITY.md`](extended-geometry/BINARY_GVC_PRIMITIVE_TRANSLATION_OBSERVABILITY.md).
The span-five one-colour run has 404 mixed packets, 125 factorial survivors,
and no torus survivor.  The written translation-degree proof gives the
zero-survivor conclusion for every projected span; these commands are exact
bounded regressions, not an extrapolated GVC theorem.  Both outputs also
replay the three-state span-two discriminant orbit identity and certify its
nonzero scale-two factorial obstruction under both fixed-external and
character-power phase laws, which shows why pairwise primitive separation
alone cannot promote a signed linear shell.

The exact small-shell all-scale-prefix eliminations are

```bash
.venv/bin/python scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --structural-certificate --maximum-wronskian-rank 7 \
  --maximum-affine-slope 5 --maximum-affine-offset 3 \
  --output artifacts/generated-results/binary_gvc_all_scale_structural_certificate.json
python3 scripts/verify_binary_gvc_cobham_carry_obstruction.py \
  --primes 3,5,7,11,13 --maximum-period 512 \
  --output artifacts/generated-results/binary_gvc_cobham_carry_obstruction.json
python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 2 --state-count 3 \
  --maximum-operator-count 3 --maximum-polynomial-count 3 \
  --maximum-scale 3 \
  --output artifacts/generated-results/binary_gvc_all_scale_orbit_circuits_span2_counts3.json
python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 3 --state-count 3 \
  --maximum-operator-count 2 --maximum-polynomial-count 2 \
  --maximum-scale 3 \
  --output artifacts/generated-results/binary_gvc_all_scale_orbit_circuits_span3_counts2.json
python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 2 --state-count 4 \
  --maximum-operator-count 2 --maximum-polynomial-count 2 \
  --maximum-scale 4 \
  --output artifacts/generated-results/binary_gvc_all_scale_orbit_circuits_span2_four_state_counts2.json
python3 scripts/research_binary_gvc_all_scale_orbit_circuits.py \
  --span 2 --state-count 4 \
  --maximum-operator-count 3 --maximum-polynomial-count 3 \
  --maximum-scale 3 --require-factorial-pairing \
  --output artifacts/generated-results/binary_gvc_all_scale_factorial_pair_circuits_span2_counts3.json
```

The first command is the structural replay.  It checks the exact confluent
Vandermonde leading symbol for all 44 block partitions through rank seven,
uses Singular to prove that the coefficient-torus smooth-conic divisor branch
of the arbitrary four-state scale-1/2/3 ideal is empty, and verifies the
rank-two pair-block converse.  It also checks exact coefficient-span ranks in
the full-span and projectively constant common-base calibrations.  Normaliz
computes the 106 Graver moves of the 35-atom integer-affine factorial
universe with slopes at most five and offsets from -3 through 3; all are
certified sums of exact same-rational-boundary transfers.  The canonical note
also replays the residue-wise rational-minor criterion for an eventually
periodic additive factorial law and detects a one-entry perturbation.  The
canonical note proves the Wronskian identity at every rank, the common-base
differential-ideal collapse, the singular-conic classification, the
boundary-transfer presentation for every integer-affine factorial profile,
and the periodic additive reduction.  Cobham's theorem makes the same
reduction for one finite state sequence automatic in two multiplicatively
independent bases.

The second command proves why complete Hall carry states do not supply that
common sequence.  It checks the Kummer and two-state digit formulas for the
central-binomial carry indicator through \(N=10{,}000\), constructs
equal-residue/opposite-output witnesses for all proposed periods through
\(512\), and checks that every sparse ray \(q p^e\) is stationary.  The
arbitrary-period witness in the canonical note proves nonperiodicity for
every odd prime.  Hence a common two-base automatic refinement retaining
the carry state is impossible by Cobham; this is a proof-route obstruction,
not a GVC counterexample.
Consequently every fixed finite affine-ray family
\(h_j(t)f_j(t)^N\) splits by proportional bases for all large \(N\), with
arbitrary scale-dependent scalar coefficients and changing active support.
The residual identity inside one proportional-base correction space is not
split by this theorem.  The computation is an exact regression for the stated
cross-base proof, not the unrestricted GVC(2) certificate; that certificate
is the written Hall-envelope theorem.

For three states, the note proves at every span that scales one and two force
all orbit-function ratios to be constant; every transferring shell is
therefore torus-empty.  The commands audit the residual zero-transfer block.
The span-two/count-three run tests 11,988 signed triples and kills its 416
scale-one survivors by scale three.  The span-three/count-two run kills all
240 survivors among 8,408 candidates at scale two.  The four-state pilot
tests 928 signed quartets and kills all 40 scale-one survivors at scale two.
These are characteristic-zero Singular saturations.  The counts audit the
proved three-state theorem and the arbitrary-coefficient four-state
scale-1/2/3 theorem; they are not those proofs.  Unrestricted GVC(2) is
proved separately by Hall-envelope separation.

The final command uses the proved all-scale factorial-ray splitting to retain
the 2,882 fixed-sign quartet rows which admit an opposite-sign equal-factorial
pairing.  It finds 142 scale-one survivors.  Character-power scale two kills
all of them.  Fixed signs leave 60 exact all-scale pair cancellations; 12
survive the other one-colour tower and six survive both independent towers.
For those six the script checks both ideal containments and proves that the
saturated pair ideal equals
`(R1^2-4*R0*R2, B1^2-4*B0*B2)`.  This is the already-safe product Veronese
block.  To replay the larger unfiltered finite-prefix census of all 52,416
signed quartets (about ten minutes on the recorded machine), omit
`--require-factorial-pairing` and `--output`; it has 928 scale-one survivors,
868 fixed-sign scale-two obstructions, 60 fixed-sign scale-three survivors,
and no character-power scale-two survivor.

The fast final regression suite is

```bash
.venv/bin/python scripts/verify_binary_gvc_all.py
```

It runs the uniform Hall/weighted-face checker, the regular-trace checker,
the repeated-digit/Newton checker, the weighted affine/factorial obstruction
checker, the default degree-ten translation-twist search, and the structural
Wronskian/four-state certificate.  It is a historical-route regression, not
the written proof of unrestricted GVC(2).  The module-only version of the
former final lemma is
disproved, while distinct bases in every fixed finite affine-ray template are
now separated, and every common-base all-scale ideal has collapsed to a
finite coefficient span.  The exact hypotheses left inside the parked route
are projective carry-rank/safe-rank-drop classification of the bounded
correction circuits and uniform extraction of such a circuit from the growing
positive-density Cartesian face.

## Factorial trace independence

```bash
make verify-factorial-trace-independence
```

To classify two rays directly, encode gamma factors as
`slope:offset[:multiplicity]`:

```bash
python3 scripts/verify_factorial_trace_independence.py \
  --compare '2:1' '1:1/2,1:1'
python3 scripts/verify_factorial_trace_independence.py \
  --compare '2:1' '1:1:2'
```

The first comparison certifies Gauss duplication and reports exponential
base `4`; the second reports the separating signature `[1/2]-[0]`.

The optional independent symbolic replay is

```bash
make verify-factorial-trace-independence-sympy
```

It uses the repository SymPy environment to simplify 322 Gauss/shift and
1,384 census certificates directly as rational-function identities.  The
default target remains dependency-free.

The dependency-free checker reconstructs all 78,124 nonzero signed slope
vectors on slopes `1,...,7` with coefficients in `[-2,2]`.  It integrates
8,134 exact zero-sum translation-orbit divisors and classifies 67,524 products
of at most three rational-offset gamma atoms into 66,140 canonical
signatures, certifying all 1,384 collisions.  It also verifies 161 Gauss
refinements before and after integer shifts, 1,000 seeded signed
transformation cases, and 24 integer-affine reductions.  Its exact
successor-divisor census classifies 82,250 products of at most four
integer-offset atoms into 72,383 classes and decomposes all 9,867 collisions
into boundary transfers.  It also replays the one-scale and entropy
collisions, `m`-fold periodic/rational-slope symmetries, and 2,187 signed
Frobenius-dilation valuation profiles at `p=2,3,5`, and separates all 276 SIC
radial-moment families `(d,r)` with `d<=48`.  These are finite regressions for
the formulas; the complete characteristic-zero gamma-affine classifier,
integer-affine boundary presentation, and exact characteristic-`p` valuation
obstruction are proved in
`extended-geometry/FACTORIAL_TRACE_INDEPENDENCE.md`.
The same proof shows that a factorial ratio with finite value set is
constant, so finite nonzero carry, automaton, sign, or torsion alphabets
introduce no additional projective class.

## Binary GVC prime-power tomography

Requires Normaliz; the pinned run used Normaliz 3.10.2.

```bash
.venv/bin/python scripts/research_binary_gvc_prime_power_tomography.py \
  --output artifacts/generated-results/binary_gvc_prime_power_tomography.json.gz \
  --summary-output artifacts/generated-results/binary_gvc_prime_power_tomography_summary.json
```

The Lawrence lifting gives the exact two-colour Graver basis on levels
`0,...,6`: 8,559 raw relations, 1,584 after translation/dilation and finite
symmetries, and 1,490 mixed support-at-least-five packet candidates on 868
projected supports.  The probes use primes `2,3,5,7,11`, exponents
`m=q*p^e` with `1<=e<=2` and `1<=q<=3`, factorial units modulo `p^2`, and
the `C2,C3` marked-character traces.  Exactly two symmetry classes survive
the configured signature.  They represent

```text
R6*B2*B3 = R0*B5*B6
R6*B1*B2 = R0*B4*B5
```

and together are the reversal orbits of
`R6*B_a*B_(a+1)=R0*B_(a+3)*B_(a+4)`, `a=0,1,2`.  Each exact projected
fibre has two states and no support-at-most-four path.  Equality of the
partitions proves that all factorial, digit, and carry data agree at every
scale, not just in the finite prime window; `C2,C3` are also identically
blind.  A `C4` character separates both classes.  This is a projected
collision census, not a GVC(2) counterexample; the fixed-character Franel
theorem below subsequently closes the family after packet exposure.

The logical JSON result SHA-256 is
`685f60b5843bca33d32034a16a8b599dcfbf46c35e80e62f747f6d9715e285eb`.
The whole-file SHA-256 hashes of the compressed full artifact and compact
summary are, respectively,
`b1cf105d161b83a5e0c23dab0ed3b2cbc0a2726970e20db760ca5ae78eb5c09b`
and
`c9380d66bdad08cb30896c0c1b31d9c5211397d85e962b48d808d90baa832522`.
The precise model and caveats are in
[`BINARY_GVC_PRIME_POWER_TOMOGRAPHY.md`](extended-geometry/BINARY_GVC_PRIME_POWER_TOMOGRAPHY.md).

The strengthened primitive-relation census is:

```bash
.venv/bin/python scripts/research_binary_gvc_prime_power_tomography.py \
  --radial-degree 7 --primitive-only \
  --primes 5,7,11,13 --max-exponent 3 --max-quotient 3 \
  --unit-power 3 --torsion-orders 2,3 \
  --output artifacts/generated-results/binary_gvc_adelic_tomography_span7.json.gz \
  --summary-output artifacts/generated-results/binary_gvc_adelic_tomography_span7_summary.json
```

It computes 34,890 raw Graver relations, 6,601 normalized relations, 6,401
mixed primitive candidates, and the exact bases of 3,107 represented support
semigroups.  Of the candidates, 4,750 are separated by a configured total
valuation.  The remaining 1,651 are exactly the all-scale scalar factorial
collisions, so there are no accidental finite-window collisions and no row
first separated only by a unit residue.  Marked digit/carry data separate
221, `C2,C3` separate another 1,427, and three normalized relations remain.
The first valuation separators occur at `p=5` for 4,557 relations, at `p=7`
for 148 which survive every configured `p=5` probe, and at `p=11` for 45
which survive both earlier primes.  No first separator needs `p=13` or
`e>1` in this span.
They are precisely the span-seven orbits of

```text
R_(s+6)*B_a*B_(a+1) = R_s*B_(a+3)*B_(a+4).
```

Every member is an exact two-state primitive collision, and `C4` separates
it.  Among all 1,430 equal marked-partition relations, the first character
separator distribution is `C2:1244`, `C3:183`, `C4:3`.

The span-seven logical result SHA-256 is
`a8128639805f9e6e0047dc39e70b20f8e939b6a76213930ab44bd0b26a35dde3`.
The whole-file SHA-256 hashes of its compressed full artifact and summary
are
`dc64e57cac395b4f98cfcb8cc0ac1cdf03c598c5fba67253f67c18e379f7f035`
and
`b268bf8ed6cb7c564e154efdcd12b7e60659533fc8cd85e134b6c82870ccabd9`.
This is an exact bounded result in the projected two-colour model, not a Hall
promotion theorem or a counterexample.  It is not the unrestricted proof;
that is the Hall-envelope theorem.

The surviving family has an exact fixed-character termination theorem.
Replay its fibre and coefficient identities with:

```bash
python3 scripts/verify_binary_gvc_six_step_packet_termination.py
```

After row normalization the support is `R6,R0,B0,B1,B3,B4`.  At scale `N`,
the complete `C2,C3`-blind fibre is
`(N-t,t,N-t,N-t,t,t)`, `0<=t<=N`, and its normalized coefficient is

```text
binom(2N,N) * sum_t binom(N,t)^3 * U^(N-t) * V^t.
```

If a further fixed finite character has relative order `h`, the
endpoint-containing rows at scales `h` and `2h` are
`U^h+V^h` and
`U^(2h)+binom(2h,h)^3*U^h*V^h+V^(2h)`.  Their only common zero in
characteristic zero is `U=V=0`, because substitution from the first row
leaves `(2-binom(2h,h)^3)*U^(2h)=0`.  Thus every fixed finite-character
promotion of the six-step family is separated, terminal, or loses support.
The proof is general; the dependency-free script enumerates blind fibres
through scale 12 and checks character orders through 32 as regressions.
Prime-dependent affine-carry promotion to one fixed packet remains unproved
inside the parked route, but is no longer required for binary GVC.

## Binary GVC nonfree-factorization tomography

The all-span consecutive-residue theorem has a quick exact regression:

```bash
.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py \
  --verify-consecutive-residues
```

The written incidence-forest proof is uniform in the cyclic order.  The
regression checks orders \(2,\ldots,16\): the \(C_q,C_{q+1}\) histograms are
injective through span \(2q-1\), while at span \(2q\) their kernel is

```text
(1,...,1,0,-1,...,-1)
```

and decomposes into the \(q\) safe beta swaps
\(R_iB_{i+q+1}=R_{i+q+1}B_i\).  This proves that the fixed marked nonfree
factorization quotient is injective at every span; it does not prove that a
prime-dependent Hall shell inherits the required fixed markings.

Run the complete span-four Hilbert/factorial/Graver census with:

```bash
.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py
```

It computes 426 raw and 90 normalized projected Graver relations, giving 65
nonfree profiles.  The exact factorial map is injective on 52; one collision
lattice is reversal-only and 12 contain nonreversal same-vector relations.
For the 11 collision profiles with at most 20 atoms, the complete
factorial-compatible Graver basis has 308 primitive moves.  Packet partitions
separate 15, \(C_2\) separates 207, and \(C_3\) separates 86.  Atom labels
separate every profile, including the two larger deferred Graver lattices.
The first factorial-only square has
\(R=\{0,4\}\), \(B=\{0,1,2,4\}\), counts \((1,3)\), and radial vector
\((8,8)\); \(C_3\) separates it.

Run the larger atom-signature censuses with:

```bash
.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py \
  --radial-degree 5 --signature-only --torsion-orders 2,3,4 \
  --output artifacts/generated-results/binary_gvc_nonfree_factorization_span5_signature.json

.venv/bin/python scripts/research_binary_gvc_nonfree_factorization.py \
  --radial-degree 6 --signature-only --torsion-orders 2,3,4 \
  --output artifacts/generated-results/binary_gvc_nonfree_factorization_span6_signature.json
```

The span-five run has 400 profiles and first-injective distribution
\(169,80,143,8\) at partitions, \(C_2,C_3,C_4\); the span-six run has
1,469 profiles and distribution \(382,358,599,130\).  Neither has an
unresolved atom signature.  The logical result hashes for spans four, five,
and six are, respectively,
`97bcdb8049b34ef0fed2bd0c9f70102e7b06d828e57a27cb09d6036395402627`,
`567712c213029dc01ca749888d21fc838b49bf21604e0a92b73997575d2dc8fd`,
and
`b8f334669719ac829b70182c4c648b8a66a92f409dd65b8945c6309ea2a6ecde`.
Their whole-file hashes are
`888f1f465c67045a2a39a157d6d6cf4872f1cb1a76a79b248e31aa658fd21d2d`,
`991f531db27e707a2155080f84751f2b3c12f9f88c957775cf6e62f26776a18d`,
and
`052928f81694395d1369a4e2c1e9973ebd1fc4d28d10e9c2216f549d4c1b1e99`.
These are exact bounded projected-semigroup computations; the all-span
claim comes from the incidence proof, not extrapolation from the census.

The accompanying exact finite-moment search is

```bash
.venv/bin/python scripts/search_binary_gvc_five_channel_descent.py
.venv/bin/python scripts/search_binary_gvc_five_channel_descent.py --frontier-suite
.venv/bin/python scripts/search_binary_gvc_five_channel_descent.py \
  --frontier-suite \
  --json-output artifacts/generated-results/binary_gvc_five_channel_pivot_clusters.json
```

The first command runs the quicker \((1,2)\) and \((1,3)\) cases.  The
frontier suite adds \((1,4)\) and \((2,3)\).  It enumerates respectively
13,288, 41,728, 95,368, and 253,576 supports with zero through three added
channels in the order-four balance boxes.  Successive saturation finds
the first unit at moments two, three, and four with distributions
\[
\begin{array}{c|rrr}
(r,s)&M_2&M_3&M_4\\ \hline
(1,2)&13082&173&33\\
(1,3)&41299&394&35\\
(1,4)&94636&679&53\\
(2,3)&252442&1074&60
\end{array}
\]
and no torus survivor.  Points outside those boxes cannot occur in a
balanced return through order four, so the computations cover arbitrary
nonnegative support with at most five channels for all four endpoint
pairs.  They are exact rational finite-moment computations, not a proof
for arbitrary endpoint orders or three operator endpoints, and not the
Hall-envelope proof of unrestricted GVC(2).  The third command additionally records every fourth-pivot
support and canonicalizes its balanced-selection rows under endpoint
exchange and permutation of the three added channels.  The 181 supports
collapse to 14 return-matrix types, all already realized at `(1,2)`.  The
generated JSON has SHA-256
`59436a3617671c4ca47cd354b45cb74abc7b9787352e725c97ebce1a304ffa16`.

The
proof and the precise remaining restricted beta--torus
coupled-convolution obstruction are in
[`BINARY_GVC_UNIFORM_FACE_TERMINATION.md`](extended-geometry/BINARY_GVC_UNIFORM_FACE_TERMINATION.md).
The checker is a regression, not the proof.

The explicit homogeneous GVC(3) counterexample and its exact consequences
are replayed by

```bash
python3 scripts/verify_gvc3_homogeneous_counterexample.py
python3 scripts/verify_gvc3_homogeneous_spillovers.py
.venv/bin/python scripts/verify_gvc3_power_tail_and_minimum.py
.venv/bin/python scripts/verify_gvc3_independent_parity_quartic.py
.venv/bin/python scripts/verify_gvc3_isotropic_harmonic_channels.py
.venv/bin/python scripts/verify_gvc3_four_coherent_channels.py
.venv/bin/python scripts/verify_gvc3_cusp_profile_suspension.py
.venv/bin/python scripts/research_gvc3_harmonic_cubic_profile.py \
  --cases alpha1_n3:8 alpha1_d_k:7 alpha1_dk:7 \
          alpha0:7 alpha0_radical:7 \
  --primes 101 103 107 --timeout 900 \
  --exact-all --msolve-threads 4
```

The first checker verifies polynomiality, homogeneity, primitivity, the
closed-form detector, and two exact finite replays through moment six.  The
all-order counterexample is proved in
[`THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md`](extended-geometry/THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md),
not inferred from the bounded replay.  The second checker replays the
homogeneous dimension and quadratic-rank spillovers.  The third checks the
maximal shifted-power formulas and the scoped one-profile minimum.  The
fourth verifies the exact characteristic-zero elimination for the
independent-linear quartic repair.  The fifth checks the full
two/three-channel isotropic coherent-state obstruction in balanced degrees
four through ten: it compiles invariant Reynolds moments, uses three modular
Gröbner runs to discover every chart cutoff, repeats all forty saturated
chart eliminations exactly over \(\mathbb Q\), and rejects the
\(\mathcal H_2\oplus\mathcal H_4\oplus\mathcal H_6\) near survivor at pure
moment five despite its surviving multiplier channel.  The sixth checks the
unique four-channel degree-eight coherent family.  It records the finite
proper-hypergeometric occupation/edge sum, proves all thirteen direction
collisions and the four-distinct `B=0` boundary exactly over \(\mathbb Q\),
and retains the `B!=0` unit bases at \(101,103,107\) as unpromoted modular
evidence.  The optional `--exact-open` flag attempts that remaining
characteristic-zero chart and fails rather than promoting a timeout.  The
seventh checks the full winding--profile--radial suspension, its cusp identity,
complete phase ladder, top Reynolds--apolar contractions, exact trace depths,
and direct
shifted-power detectors for the non-power profile \(S=1+z\).  The final
command is an exact search with modular discovery replay in the complete
harmonic-cubic repair:
it compiles seven invariant weight channels, covers the nonzero-even-part
chart by three pivot strata, and audits the radical of the zero-even-part
boundary.  It requires Singular and msolve.  All nine projective and
boundary saturations are exact over \(\mathbb Q\); the three-prime runs are
retained as discovery replay.
The first seven commands reproduce seven `gvc3_*.json` artifacts.  The new
four-channel artifact deliberately has mixed exact/modular status; the other
six are exact.  The final command writes the separate exact harmonic-cubic
artifact in `artifacts/generated-results/`.  These results
disprove GVC in every dimension at least three but do not disprove the
ordinary-Laplacian/Hessian-nilpotent conjecture.

The integration replay on 2026-08-02 used system Python 3.12.3 for the first
two dependency-free commands and the locked `.venv` Python 3.13.5 for the
SymPy commands.  Regenerating the first two artifacts normalized JSON list
formatting without changing their parsed content.  The counterexample
artifact changed from
`sha256:f05e5bee5c9b9aab5e245026f99af30a4f379bd88b9f3163fbdd7859f56aba06`
to
`sha256:ef44b4d7390ca261c432d23bcfc7b262062d3027b4f46ca9ccef1b9c556ec04d`;
the spillover artifact changed from
`sha256:ef27d18337a34a527140a799e15a3a242b008bc1b6c22e403ea775d86df31b50`
to
`sha256:6810131f43f822c39d4abed682e97570254faed6b1e19080cba1559821f2a666`.
`jq` comparison of each old/new pair is exactly `true`.

The earlier exact three-variable tagged-lift reduction and its bounded
extensions are replayed by

```bash
.venv/bin/python scripts/research_three_variable_gvc_tagged_lift.py
```

The target first replays the coordinate-only detector extracted from the
two-pair witness through order six and uses the proved equal-channel
degree-tag identity.  Over
\(\mathbb Q\), it then computes the complete binary-cubic operator-jet
moment ideals for the literal Long tag: moments one through four remain
nonunit and moment five gives the unit ideal.  The same run performs three
explicitly experimental calculations over \(\mathbf F_{101}\): the
canonical rank-five auxiliary chart, the normalized factor-compatible
cubic-profile chart, and 200 deterministic general cubic-profile fibers.
Their output is written to
[`three_variable_gvc_tagged_lift.json`](artifacts/generated-results/three_variable_gvc_tagged_lift.json).
The exact formulas, scope distinctions, and remaining mixed-order target
are documented in
[`THREE_VARIABLE_GVC_TAGGED_LIFT.md`](extended-geometry/THREE_VARIABLE_GVC_TAGGED_LIFT.md).
These are retained architecture exclusions; the homogeneous counterexample
above now settles GVC(3) negatively.

The coupled order/degree-\((2,3,4)\) continuation is replayed by

```bash
.venv/bin/python scripts/research_three_channel_gvc_lift.py
```

It enumerates all \(56\) oriented rank-three parallelograms on the positive
weighted-quartic plane, computes their exact rational moment ideals, and
then checks the persistent five-term radical, the complete quartic repair
on the polynomial side, the sparse activated operator endpoints, and the
complete odd-quartic operator/polynomial jet.  The remaining radicals have
written all-order mixed cutoffs.  The generated
[`three_channel_gvc_lift.json`](artifacts/generated-results/three_channel_gvc_lift.json)
records these characteristic-zero results.  The odd chart closes at
moment six with radical \((A,S,RU)\); the simultaneous complete even-and-odd
quartic total space remains open.

The first repeated-root continuation is the migrating defect-one ansatz
\((\Lambda_4+\Lambda_5,P_5+P_4)\).  Run its faithful-characteristic samples,
followed by the conditioned defect-two \((\Lambda_6,P_3)\) search, with:

```bash
.venv/bin/python \
  scripts/search_binary_repeated_quartic_gvc_jets_mod_p.py
```

The pinned run takes roughly two minutes.  Use `--quick` for a small
non-pinning regression.  The generated
[`binary_repeated_quartic_gvc_jet_search.json`](artifacts/generated-results/binary_repeated_quartic_gvc_jet_search.json)
is bounded experimental evidence only.  Its scope, the support-separator
proof for the two fifth-moment survivors, and the unique defect-two
fourth-moment failure are documented in
[`BINARY_REPEATED_QUARTIC_GVC_JET_SEARCH.md`](extended-geometry/BINARY_REPEATED_QUARTIC_GVC_JET_SEARCH.md).

The first rank-efficient ordinary-Laplacian lift is excluded by

```bash
.venv/bin/python scripts/verify_dvorsky_one_pair_schur_obstruction.py
```

Pairing \(t\) with one new variable \(s\) makes
\[
 \partial_a\partial_d-\partial_b\partial_c+\partial_t\partial_s
\]
a nondegenerate quadratic operator in six variables.  The checker first
retains the homogeneous cubic regression, then parametrizes the
unrestricted transverse two-jet of an arbitrary polynomial or formal
harmonic lift.  It proves the exact axis identity
\[
 \widetilde\Delta^2(F^2)=12t^2-8\rho t
 \quad\text{modulo }(a,b,c,d,s).
\]
Thus no degree mixture can repair the canonical six-variable hyperplane
lift.  Different quadratic completions, additional blocks, and nonlinear
specializations remain open.  See
[`DVORSKY_ONE_PAIR_SCHUR_OBSTRUCTION.md`](extended-geometry/DVORSKY_ONE_PAIR_SCHUR_OBSTRUCTION.md).

The still-open bidegree-\((3,3)\) classification is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_bidegree33_sextic_slice.py
.venv/bin/python scripts/verify_two_pair_sic_bidegree33_anchor_jacobians.py
.venv/bin/python scripts/verify_two_pair_sic_bidegree33_boundary_family.py
```

The checker eliminates the full seven-parameter one-sided locus in the
sixteen-dimensional bidegree-\((3,3)\) coefficient space, giving its exact
dimension seven, projective degree twenty, and Gröbner-basis size \(148\).
It also evaluates the full moment Jacobian without expanding the
sixteen-variable moments and verifies a displayed nonzero
\(13\times13\) integer minor for \(\mu_1,\ldots,\mu_{13}\). Thus these
moments are algebraically independent and attain the invariant-quotient
dimension bound. It then computes the invariant Hilbert series from the
sixteen exact \(\mathrm{SL}_2\)-weights: the proposed numerator for
degrees \(1,\ldots,13\) has coefficient \(-2186\) in degree \(63\), so
those moments cannot be a homogeneous system of parameters and their zero
fiber has an extra semistable component. Replacing \(\mu_{13}\) by
\(\mu_{14}\) gives the least-total-degree corrected candidate; the checker
verifies another exact rank-thirteen minor and a nonnegative proposed
Hilbert numerator through degree \(100\), while making no zero-fiber
claim. It also inverts the full Clebsch--Gordan basis to recover the global
quadratic discriminant and proves that moments \(1,\ldots,4\) have only
the origin on the maximal-torus fixed diagonal slice, with four
seventh-power certificates.
It then constructs the highest \(\operatorname{Sym}^6\) summand by the
\(\mathfrak{sl}_2\) lowering chain.  On this binary-sextic slice, exact
elimination and power reductions prove that moments \(2,4,6,10\) have the
same radical as the \(L^4Q\) nullcone; the ten nullcone generators have
power certificate \((1,5,5,5,5,5,5,5,5,5)\).  Thus this slice contains no
SIC(2) counterexample.  The same checker proves direct ideal equality with
the \(L^3R\) quartic nullcone from moments \(2,3\), and with the \(L^2\)
quadratic nullcone from moment \(2\).  Thus all three pure irreducible
summands are closed.  On the normalized non-null quadratic branch, moments
through order six prove \(c^6\) lies in the
\(\operatorname{Sym}^4\oplus\operatorname{Sym}^2\) moment ideal over
\(\mathbb Q\), excluding that branch.  For
\(\operatorname{Sym}^6\oplus\operatorname{Sym}^2\), the checker records
only the finite-field result over \(\mathbb F_{32003}\): even moments
through order fourteen give a basis of size \(7576\) and contain \(c^{25}\)
but not \(c^{24}\).  The exact characteristic-zero lift and the full
mixed-summand problem remain open; see
[`TWO_PAIR_SIC_BIDEGREE33_FRONTIER.md`](extended-geometry/TWO_PAIR_SIC_BIDEGREE33_FRONTIER.md).
The second checker normalizes a non-null quadratic to \(2XT\), covers the
full higher \(\operatorname{Sym}^6\oplus\operatorname{Sym}^4\) locus by
five residual-torus chart orbits, and proves by displayed nonzero exact
eleven-by-eleven Jacobian determinants that moments \(2,\ldots,12\) are
algebraically independent on every chart orbit.  This is a
dimension-sized coordinate theorem, not a zero-fiber exclusion. It also
verifies the exact normalized formula for \(\mu_2\), whose derivatives in
the five opposite-weight chart variables are respectively
\(-72,432,-1080,336,-1344\). Thus one variable is eliminated without
saturation on every chart. On \(s_0=1\), the checker then eliminates
\(s_6\) and proves that \(\mu_3\) is affine in \(s_5,t_4\), recording its
two explicit pivot coefficients and their common boundary. On the natural
two-parameter plane
\(t_0=a,t_3=3a,t_4=b,s_6=(14ab+70)/3\) in that boundary, it verifies
\(\mu_3=1866240a^3\), computes the fourth moment, and checks an explicit
unit certificate modulo \(a^3\). Thus this sparse boundary plane contains
no moment-zero point. Finally it evaluates the full chart Jacobian at
exact rational points in \(A\ne0\), \(A=0,B\ne0\), and \(A=B=0\).
Together with the independent gradients of \(A,B\), the nonzero
determinants prove maximal restricted differential ranks \(11,10,9\).
The checker also verifies the constant triangular pivots
\(\partial A/\partial t_3=1\) and
\(\partial B/\partial s_4=-3\), with both cross derivatives zero.
Thus \(A=0\) eliminates \(t_3\) globally.  The fully substituted
\(A=0,B\ne0\) export has nine effective variables after the \(\mu_3\)
pivot, and the \(A=B=0\) export eliminates \(t_3,s_4\) and has eight
effective variables.
The third checker enlarges the common-boundary plane to the exact
four-parameter family
\[
 s_4=-4q^2,\quad s_5=h,\quad
 (t_0,t_1,t_3,t_4)=(a,q,3a,b),\quad
 s_6=(14ab-168aq+70)/3.
\]
Here \(A=B=\mu_2=0\) identically.  Moments \(3,\ldots,6\) leave a
zero-dimensional quotient of length \(372\), and adjoining moment \(7\)
gives the unit ideal over \(\mathbb Q\).

The reduced common-boundary fiber calculation is replayed by

```bash
.venv/bin/python scripts/verify_two_pair_sic_bidegree33_boundary_fiber.py
```

After the constant \(t_3,s_4,s_6\) substitutions, it records a rational
\(\mu_3=0\) base point at which
\((\mu_4,\mu_5)\) cuts out a length-six quotient in the two fiber
variables \(s_5,t_4\).  Its six standard monomials are
\(1,t_4,t_4^2,t_4^3,s_5,s_5t_4\).  Openness promotes this to an exact
characteristic-zero rank-six theorem on a nonempty open of the
\(\mu_3=0\) base; it is not a zero-fiber exclusion.

The generic boundary quotient and its denominator strata are replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_boundary_generic_quotient.py
```

The checker constructs the sparse moments integrally and works over the
characteristic-zero rational-function field of the six base variables.
It obtains a three-element basis and quotient length six for
\((\mu_4,\mu_5)\).  Both \(\mu_6,\mu_7\) reduce to six fiber coordinates.
The \(\mu_6\) denominators are supported on
\(L=s_1t_0-t_1\) and
\(Q=s_1^2-s_2-(13/3)t_0^2\); \(\mu_7\) introduces the additional
quartic divisor \(J\) recorded in the artifact.  The checker also proves
that the \(t_4^3\)-coefficient of the \(\mu_6\) normal form is
affine-linear in \(t_2\), with constant derivative
\(-100078239744000\), so it globally eliminates \(t_2\) on the
generic \(LQJ\)-open.  It is also affine-linear in \(s_3\), with the
explicit alternate eight-term pivot \(H\).
In the adapted coordinates it verifies
\(J=(99Q+155t_0^2)^2+30420L^2\) and
\(H=32J+1179Q(99Q+155t_0^2)\), giving the exact localized exclusion
\((J,H):(LQ)^\infty=(1)\).
The checker also constructs an irreducible quartic number field and an
explicit point on \(J=\mu_3=0\).  At that point
\((\mu_4,\mu_5)\) has length five, with initial ideal
\((s_5^2,t_4^3,s_5t_4^2)\) and basis
\(1,t_4,t_4^2,s_5,s_5t_4\).  It then proves the same statement at the
generic point of \(J=0\).  Over
\(\mathbb Q(\alpha)(s_1,s_3,t_0,L,t_2)\), \(\alpha^2=-30420\), a
quadratic-pair fraction-free calculation constructs a three-element
Gröbner basis with supports \(6,7,6\) and leading monomials
\(s_5^2,s_5t_4^2,t_4^3\).  One pair is removed by the product criterion
and the final pair reduces exactly to zero in five steps.  Generic
rational-function-field calculations at both split roots modulo 47 and
101 independently reproduce the basis.  The same exact quadratic-pair
reducer sends \(\mu_6,\mu_7\) to normal forms supported on all five
standard monomials in respectively five and ten pseudo-reduction steps.
After solving the constant \(t_2\)-pivot, the checker forms the
642-term cubic \(P(s_3)\) coming from \(\mu_3\) and the cubic
\(s_5t_4\)-coefficient \(C(s_3)\) of the remaining \(\mu_6\) normal
form.  Their exact resultant factors as
\[
 \operatorname{Res}_{s_3}(P,C)=L^6Q^6\mathcal R_{63}.
\]
The residual factor has degree 63 and 6702 rational terms.  Its
degree-preserving reduction modulo 47 is irreducible, proving
\(\mathcal R_{63}\) irreducible over \(\mathbb Q\); reduction modulo
101 independently reproduces the degree and irreducibility.
For the next \(t_4^2\)-coefficient, the checker obtains
\(\operatorname{Res}_{s_3}(P,C_2)=L^9Q^6\mathcal T_{66}\).
Both reductions of \(\mathcal T_{66}\) are irreducible, and the checker
verifies modulo both primes that
\(\gcd(\mathcal R_{63},\mathcal T_{66})=1\).  Thus the principal
\(\mu_6\)-zero base has codimension at least two rather than a surviving
degree-63 divisorial component.
For the two cubics \(P=as_3^3+bs_3^2+cs_3+d\) and
\(C=es_3^3+fs_3^2+gs_3+h\), it also constructs the direct linear
pseudo-remainder \(V_1s_3+V_0\).  Over \(\mathbb Q\), \(V_1,V_0\)
have degrees \(60,63\) and respectively 2105 and 5170 terms.
Their reductions modulo both primes are coprime to
\(\mathcal R_{63}\).  Consequently the degree-63 incidence branch has
the dense rational pivot \(s_3=-V_0/V_1\); the checker does not assert
that this pivot covers every component of the residual codimension-two
intersection.
On the divisor strata \(L=0\), \(Q=0\), and \(L=Q=0\), the exact
characteristic-zero quotient lengths are respectively six, five, and
five.  Their changed standard monomial bases are recorded, and the
normal forms of both \(\mu_6,\mu_7\) occupy every basis coordinate.
Independent calculations over \(\mathbb F_{47}\) and
\(\mathbb F_{101}\) replay all quotient shapes.  This is a finite
quotient certificate, not a full boundary unit certificate.

The content-preserving corrected-boundary continuation is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_corrected_boundary.py \
  --branch LQ --prime 0 --through 10 --export-only \
  --include-branch-table --deepest-ffnf --deepest-ffnf-through 7 \
  --t0-zero-branch-table --t0-open-rank-six --timeout 300 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_corrected_boundary_deepest.json
```

It reconstructs the generic, \(L=0\), \(Q=0\), and \(L=Q=0\)
\((\mu_4,\mu_5)\) algebras over \(\mathbb Q\), imports the separate exact
rank-five \(J=0\) result, and proves that on \(L=Q=0\) the algebra after
\((\mu_3,\mu_4,\mu_5)\) has length \(15\) over
\(\mathbb Q(s_1,s_3,t_0)\).  The checker deliberately does not call
`cleardenom` on these exports: that command removes base content and can
turn a nonconstant base equation into \(1\).  In the rank-fifteen algebra,
it pseudo-reduces \(\mu_6,\mu_7\) exactly without division by a base
polynomial; both primitive normal forms occupy all fifteen standard
monomials.  It also proves
\((\mu_3,\ldots,\mu_{10})=(1)\) over \(\mathbb Q\) on
\(L=Q=t_0=0\).  More generally it saturates the generic, \(L\), \(Q\),
\(J\), and \(L=Q\) strata by their specialized principal opens and returns
the unit ideal through \(\mu_{10}\) on every \(t_0=0\) stratum.  These five
opens partition the adapted \((L,Q)\)-plane, so this closes the entire
branchwise \(t_0=0\) divisor.  On the remaining open it normalizes
\(t_0=1\), sets \(u=s_0^{-1}\), eliminates \(t_4,t_3,s_4\), and verifies
that \(\mu_3\) is fiber-independent while \((\mu_4,\mu_5)\) cuts out an
exact rank-six algebra in \(s_6,s_5\).  The exact \(\mu_6\) normal form
occupies all six standard monomials; the three leading coefficients have
the explicit \(K,H,Q_*KJ_*H\) factorization recorded in the artifact,
with \(K=4A_*-15Q_*\) and \(H=4J_*-15A_*Q_*\).  On \(K=0\), a changed
basis retains length six away from \(\ell J_*=0\); on the reduced
\(K=H=0\) linear slice, another changed basis has leading ideal
\((s_5^2,s_6^3)\), length six, and a six-coordinate \(\mu_6\) normal
form.  The separate exact \(J_*=0\) calculation covers the remaining
\(K=J_*=0\) intersection.  Finally it parametrizes the rational conic
\(H=0\) by (5.12t) and obtains the exact generic leading ideal
\((s_6^2,s_5^3)\) and length six; the omitted parametrization point lies
on \(J_*=0\).
Exact and mod-\(47\) pseudo-reductions of
\(\mu_8\), and the reordered mod-\(47\) full deepest solve, reach the
recorded \(600\)-second bounds.  The \(t_0\)-open common-root equations,
lower-dimensional coefficient specializations, inherited \(Q_*,J_*\)
branch radicals, orders \(7,8,9,10,11,12,14\), and the rational radical
remain open.

The first \(t_0\)-open common-root step has a separate subsecond exact
certificate:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_corrected_boundary.py \
  --branch generic --prime 0 --t0-open-fixed-fiber --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_open_fixed_fiber.json
```

At the rational base (5.12u), it verifies \(\mu_3=0\), nonvanishing of
\(Q_*J_*KH\), length six for
\(\mathbb Q[s_6,s_5]/(\mu_4,\mu_5)\), and the exact unit ideal
\((\mu_4,\mu_5,\mu_6)=(1)\).  This proves that the first norm is not
identically zero on the local base component, but does not expand or
classify its exceptional divisor.

The same fixed point extends to the exact rational \(\mu_3=0\) curve
(5.12v).  Its first norm and the next Fitting coefficient are replayed
by:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_corrected_boundary.py \
  --branch generic --prime 0 --t0-open-curve-norm --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_open_curve_norm.json
```

The checker obtains an irreducible degree-\(198\) numerator for
\(\det M_{\mu_6}\), with denominator degree \(144\).  The coefficient of
\(z\) in \(\det(M_{\mu_6}+zM_{\mu_7})\) has numerator and denominator
degrees \(209,153\), and its numerator is coprime to the degree-\(198\)
norm.  Thus \((\mu_6,\mu_7)\) has no common root on the norm divisor
where the curve and border chart are defined.  The norm denominator
factors into the specialized \(Q_*=0\), curve-pole, and \(J_*=0\)
factors with degrees \(2,3,4\).  Exact degree-two and degree-four
number-field calculations give length five for \((\mu_4,\mu_5)\) and
the unit ideal after adjoining \(\mu_6,\mu_7\) on both \(Q_*=0\) and
\(J_*=0\).  The cubic factor is a genuine pole with coprime numerator,
not an affine point of the parametrized curve.  Thus every defined point
of this rational curve is excluded.

Directional modular interpolation shards for the first two Fitting
coefficients use paired roots of the quadratic \(\mu_3(s_3)\):

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 1019 --variable s1 --sample-count 450 \
  --training-count 400 --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_fitting_s1_mod1019.json
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_t0_fitting_degree_scout.py
```

The generated scout aggregates all five directions through two base
points, at primes \(1019,2039\): 6750 paired samples fit fifteen rational
line reconstructions and 750 unused pairs verify them.  The common
observed denominator models are (5.12x).  This supplies stable degree
bounds only; it does not reconstruct the dense five-variable
numerators.

The same engine can evaluate the complete degree-six determinant pencil
on deterministic random paired-root shards and replay only its common-zero
candidates:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 47 --random-seed 102 --sample-count 450 \
  --max-attempts 2000 --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_pencil_random_p47_seed102.json
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_t0_pencil_random_scout.py
```

The expanded aggregate also includes the smallest admissible prime and
the direct divisor mode:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 43 --random-seed 404 --sample-count 450 \
  --max-attempts 2000 --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_pencil_random_p43_seed404.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 43 --random-seed 501 --stratum Q --sample-count 450 \
  --max-attempts 2000 --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_random_p43_seed501.json
```

The aggregate contains forty-four generic shards at primes
\(43,47,59,71\): 19800 accepted paired bases, 39600 evaluated roots of
\(\mu_3\), and twenty direct length-one common roots through \(\mu_7\).
All twenty have block rank five through \(\mu_7\) and rank six after
adding \(\mu_8\); nineteen use the leading \(M_6\) pivot and one uses a
second pivot.  Four direct divisor scouts evaluate another 3600 roots.
The generic quotient lengths are five on \(Q,J\) and six on \(K,H\);
one reduced through-\(\mu_7\) point occurs on each of \(Q,J\), and both
are excluded by \(\mu_8\).  This is a bounded modular scout, not a
reconstruction or a characteristic-zero exclusion.

The specialization, rank-complement, leading-border, and projected-border
continuation uses only short guarded jobs.  Representative producers and
the aggregate verifier are:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 43 --random-seed 602 --stratum QJH \
  --sample-count 450 --max-attempts 10000 --retain-pairs --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_QJH_random_p43_seed602.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_fitting_interpolation.py \
  --prime 43 --random-seed 901 --sample-count 225 \
  --max-attempts 6000 --pivot-scout --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_rank_complement_random_p43_seed901.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_stratum_leading.py \
  --prime 43 --stratum Q --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_leading_mod43.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_stratum_border_resultant.py \
  --prime 43 --stratum Q --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_border_resultant_mod43.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_stratum_leading.py \
  --prime 0 --stratum Q --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_leading_exact.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_stratum_border_resultant.py \
  --prime 0 --stratum Q --timeout 20 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_border_resultant_exact.json
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_t0_strata_rank_continuation.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_Q_residual_slice.py
```

The verifier pins 31 direct-stratum artifacts.  Across \(Q,J,K,H,KH\),
\(QJH,JH,JK\), \(a_2=0\), and the repeated-root stratum, seventeen
sampled reduced common roots through \(\mu_7\) are all excluded by direct
\(\mu_8\) evaluation.  Twelve rank-complement shards test another 6300
\(\mu_3\)-roots: four miss both selected pivots, but all four have full
rank in \([M_6\ M_7]\), and no joint-rank-at-most-four point occurs.

The `stratum_leading` command uses `liftstd`; the LCM of its leading
coefficients is the specialization border.  Modulo \(43\) it gives the
nine irreducible border profiles in (5.12z) and contains every retained
length-drop point.  On \(Q,J,JK\) its sampled zeros are exactly the
length drops.  The degree-at-most-four \(Q,J\) point-cloud kernels contain
only the ambient cubic equation and its five linear multiples.

The border-resultant command projects the quadratic \(\mu_3(s_3)\)
against this border and factors the result.  The residual factors have
the degree, term-count, and multiplicity profiles in (5.12aa).  The
linear pseudo-remainder \(A s_3+B\) is coprime to every residual factor,
so all have the dense pivot \(s_3=-B/A\).  These are exact finite-field
calculations, not characteristic-zero component certificates.

The two `--prime 0 --stratum Q` commands promote the \(Q\)-row.  Over
\(\mathbb Q\), its leading border is irreducible of degree \(36\) with
588 terms.  The exact degree-\(76\) resultant factors as
\(c\,u^{20}J_Q^4R_{20}^2\), where \(J_Q\) is the inherited four-term
quartic and \(R_{20}\) is irreducible of degree \(20\) with 200 terms.
The exact linear pivot is coprime to \(R_{20}\), and reduction modulo
43 matches the modular artifacts up to units.  The next step is custom
arithmetic in this degree-five extension and the remaining modular
degree-five/degree-six extensions; direct Gröbner recomputation is
deliberately not part of this command sequence.

The final verifier closes the exact \(s_1=\ell=0\) one-parameter slice
of that \(Q\)-residual component.  The degree-\(100\) factor of the
\(\mu_6\) norm is removed by the first
\(\det(M_{\mu_6}+zM_{\mu_7})\) coefficient.  On the only remaining
cubic factor, the Kummer modulus specializes to \(u^4\), so the whole
scheme-theoretic fiber is supported at \(u=0\) and is empty after the
chart localization \(u=s_0^{-1}\).  It writes
`artifacts/generated-results/two_pair_sic_bidegree33_Q_cubic_exceptional_factor.json`;
the whole-file SHA-256 is
`69caf2d4b83fc2d70e1ce46b945471b6d0666b0ebbe12de6d20ec660a6e7114a`.
This is an exact slice exclusion, not an exclusion of the full
\(Q\)-residual component.

Three corrected exact closed-fibre calibrations away from that slice use
\(\ell=s_1u-t_1\):

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual.py \
  --through 7 --timeout 20 \
  --specialize s1=5 --specialize ell=7 --specialize u=2 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_fiber_s1_5_ell_7_u_2_exact.json
```

The analogous artifacts at \((s_1,\ell,u)=(7,4,3)\) and
\((11,9,5)\) give the same exact degree-five extension, length-four
fibre, and unit ideal through \(\mu_7\).  They prove that the exceptional
Fitting locus is proper, not that it is empty.

The projective source of that length-four fibre and the smaller
quadratic-remainder elimination are checked with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual.py \
  --through 7 --timeout 60 --projective-probe \
  --specialize s1=5 --specialize ell=7 --specialize u=2 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_projective_probe_s1_5_ell_7_u_2_exact.json
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_Q_residual_infinity.py \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_Q_residual_infinity_exact.json
```

The second command proves symbolically that the highest fibre-degree
parts of \(\mu_4,\mu_5\) share
\(6s_1us_5-s_6\), that any second infinity point is supported on the
inherited factor \(J_Q=6084\ell^2+4805u^2\), and that the tangent
determinant is a nonzero scalar times \(u^5D\), where \(D\) is the
residual border factor.  The first command confirms full projective
coprimality and infinity length two at the exact closed fibre.  It also
reduces \(\mu_5,\mu_6,\mu_7\) modulo the quadratic \(\mu_4\): their
\((s_6,s_5)\)-degrees become \((1,3),(1,3),(1,4)\), and the four
polynomials generate the unit ideal at that fibre.  These facts justify
the smaller univariate subresultant continuation; they do not yet exclude
its exceptional parameter locus globally.

The closed-point determinant-pencil oracle is reproduced by:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual_border_basis.py \
  --stage pencil --original-only --prime 43 --timeout 60 \
  --specialize s1=5 --specialize ell=7 --specialize u=2 \
  --write-moments-artifact \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_Q_corrected_moments_exact.json \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_pencil_s1_5_ell_7_u_2_mod43.json \
  --singular-output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_pencil_s1_5_ell_7_u_2_mod43.singular.log
```

It verifies the pivot open, \(\mu_3\), and the leading border before
constructing the length-twenty quotient and all coefficients of
\(\det(M_{\mu_6}+zM_{\mu_7})\).  This is a modular interpolation oracle,
not a characteristic-zero certificate.

One batched line reconstruction and the two transverse degree scouts are
reproduced by:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual_border_basis.py \
  --stage pencil --original-only --pivot-mode equation --prime 1009 \
  --scan-variable s1 --scan-values 0:500 \
  --specialize ell=7 --specialize u=2 \
  --moments-artifact \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_Q_corrected_moments_exact.json \
  --reconstruct-training-count 400 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_pencil_scan_s1_ell7_u2_500_mod1009.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual_border_basis.py \
  --stage pencil --original-only --pivot-mode equation --prime 1009 \
  --scan-variable ell --scan-values 0:1009 \
  --specialize s1=5 --specialize u=2 \
  --moments-artifact \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_Q_corrected_moments_exact.json \
  --reconstruct-training-count 800 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_pencil_scan_ell_s1_5_u2_full_mod1009.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_residual_border_basis.py \
  --stage pencil --original-only --pivot-mode equation --prime 1009 \
  --scan-variable u --scan-values 1:1009 \
  --specialize s1=5 --specialize ell=0 \
  --moments-artifact \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_Q_corrected_moments_exact.json \
  --reconstruct-training-count 800 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_residual_pencil_scan_u_s1_5_ell0_full_mod1009.json
```

The first fit validates all 21 rational pencil coefficients on 99
held-out points.  Representative numerator/denominator degrees are
\(209/91\) in \(s_1\), \(420/212\) in \(\ell\), and \(270/100\) in
\(v=u^2\) on the \(\ell=0\) chart.  These are exact finite-field line
certificates and interpolation estimates, not characteristic-zero
multivariate reconstruction.

The corrected sparse unspecialized ratio chart
\(\lambda=(s_1u-t_1)/u,\ v=u^2\) is profiled with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_flint.py \
  --stage evaluated --prime 43 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_unspecialized_flint_evaluated_mod43.json
```

The pivot and moment stages finish quickly and the displayed evaluation
finishes in about thirty seconds.  The next unspecialized basis stage was
stopped after four minutes under a \(3\)-GB cap, so it is not a
reproduction command for a completed artifact.

The bounded pre-pivot quadratic elimination is reproduced by:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_flint.py \
  --stage prepivot5 --prime 43 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_unspecialized_flint_prepivot5_mod43.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_flint.py \
  --stage prepivot6 --prime 43 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_unspecialized_flint_prepivot6_mod43.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_flint.py \
  --stage prepivot7 --prime 43 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_unspecialized_flint_prepivot7_mod43.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_t0_Q_flint.py \
  --stage raw_equations --prime 43 --timeout 60 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_t0_stratum_Q_unspecialized_flint_raw_equations_mod43.json
```

The first three commands pseudo-divide by the quadratic \(\mu_4\) before
substituting the dense \(s_3\)-pivot.  Their modular linear remainders
have respectively \(25,29,41\) raw monomials.  The last command forms
the three \(s_6\)-elimination equations in about five seconds; their
\((s_5,s_3)\)-degree pairs are \((5,12),(5,10),(6,12)\).
These are exact computations over the displayed finite-field residual
extension, not characteristic-zero component exclusions.  The guarded
`prepivot_cross6` stage confirms that expanding even the smallest
equation after the dense pivot exceeds 180 seconds; the intended next
consumer is therefore a resultant/subresultant implementation that keeps
the pivot linear.

The \(L=1\) trace/norm reconnaissance treats
\((s_3,s_5,t_4)\) as a rank-twelve finite fiber after
\(\mu_3,\mu_4,\mu_5\).  Export the rational-function-field generator
matrices, sample exact multiplication invariants at two good primes, and
replay every sampled joint-rank-drop point against the corrected later
moments with:

```bash
.venv/bin/python \
  scripts/explore_two_pair_sic_bidegree33_boundary_coefficients.py \
  --prime 47 --orders 2,3,4,5,6 --trace-norm --timeout 600
.venv/bin/python \
  scripts/explore_two_pair_sic_bidegree33_boundary_coefficients.py \
  --prime 47 --orders 2,3,4,5,6,7 \
  --trace-norm --trace-samples 1200 --timeout 600
.venv/bin/python \
  scripts/explore_two_pair_sic_bidegree33_boundary_coefficients.py \
  --prime 101 --orders 2,3,4,5,6,7 \
  --trace-norm --trace-samples 250 --timeout 600
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_boundary_trace_candidates.py \
  --prime 47 --timeout 600
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_boundary_trace_slice.py \
  --primes 47,101 --s1 1 --t0 1 --timeout 600
```

Every accepted specialization at both primes has leading ideal
\((s_3^2,s_5^2,s_5t_4^2,t_4^4)\).  The mod-\(47\) scan finds three
reduced length-one common \((\mu_6,\mu_7)\) fibers; corrected \(\mu_8\)
is nonzero at all three.  The commands write their exact finite-field
data under `artifacts/generated-results/`.  This is sampled modular
evidence only, not an exhaustion of the rank-drop locus or a
characteristic-zero nullcone certificate.  The final command is stronger:
on the complete slice \(L=s_1=t_0=1\), it finds length \(1128\) through
\(\mu_7\) modulo both primes and then uses Singular `modStd` with exactness
one to certify over \(\mathbb Q\) that adjoining corrected \(\mu_8\)
gives the unit ideal.  This is an exact slice exclusion, not a global
boundary certificate.  The same command keeps \(s_1\) free and checks the
larger \(L=t_0=1\) hyperslice with `msolve`, and separately keeps
\(t_0\) free on \(L=s_1=1\).  Corrected
\((\mu_3,\ldots,\mu_8)\) is the unit ideal modulo both primes and
directly over \(\mathbb Q\) on both hyperslices.  Each exact rational
`msolve` run uses deterministic sparse linear algebra (`-l 2`) and
outputs the one-element Gröbner basis \([1]\).  Singular's
separate verified modular reconstruction of the first larger hyperslice
hit the recorded \(600\)-second bound, but rational `msolve` supplies
both exact hyperslice certificates.  The corresponding unfixed
seven-variable \(L=1\) exact solve also hit its recorded \(600\)-second
bound, so the result is not a full \(L\)-open certificate.

Modular full-chart reconnaissance is available separately:

```bash
.venv/bin/python scripts/explore_two_pair_sic_bidegree33_full_anchor.py \
  --prime 43 \
  --orders 2,3,4,5,6,7,8,9,10,11,12,14 \
  --timeout 180 --backend msolve --charts s0 \
  --branch s0-A-open-sparse --msolve-linear-algebra 44
```

Recorded Singular and `msolve` runs reach the dense chart ideal quickly
but time out, including after the first pivot and on the common
third-moment pivot boundary. The fully substituted \(A\ne0\) branch
exports as eleven equations in ten variables, but the recorded `msolve`
run terminates inside the solver. Sparse principal-open encodings avoid
that expansion, but the full corrected \(A\)- and \(B\)-open systems still
exceed the recorded bounds. These are computational diagnostics, not
evidence for or against the anchor. Any later modular result remains
experimental until it has an exact characteristic-zero certificate.
The conceptual continuation introduces no new computed certificate:
[`TWO_PAIR_SIC_MOMENT_NULLCONE_PROGRAM.md`](extended-geometry/TWO_PAIR_SIC_MOMENT_NULLCONE_PROGRAM.md)
now records the bidegree-\((4,4)\) falsification of the formerly all-\(d\)
moment--nullcone conjecture, together with the surviving
Hilbert-series degree-selection layer, the quadratic-anchor target, and
the common-root synchronization induction. It also derives the all-\(d\)
nullcone dimension and the Krull-height lower bound of
\((d+1)^2-3\) global moments, and explains why the balanced higher-pair
analogue is already false for every \(n\geq3\) by the padded three-pair
counterexample.

The same witness has a separate factorial-functional translation and finite
prefix search:

```bash
make verify-factorial-moments
```

It extracts the torus diagonal of every checked power and verifies the
all-order binomial formula behind its zero factorial value.  The diagonal
sequence is not a power sequence: at order two its true translated middle
coefficient is `-8`, while squaring the first diagonal gives `-2`, with
factorial values zero and twelve respectively.  The target also reconstructs
cyclotomic linear forms in exact prime cyclotomic rings and verifies
witness-derived `2r`-term quartics whose first `2r-1` moments vanish for
`r=3,5,7,11`.  The displayed general identities and the homogeneous
linear-class minimality proof are in
[`FACTORIAL_MOMENT_WITNESSES.md`](extended-geometry/FACTORIAL_MOMENT_WITNESSES.md).
These are sharp finite-prefix examples, not counterexamples to the
Factorial Conjecture.

The first exact sparse frontier and fixed binary-form cutoffs are checked by

```bash
make verify-factorial-frontier
```

It exhausts all 3,276 three-monomial supports in two variables through total
degree six and finds no nonzero coefficient point with moments one through
three zero.  It also exhausts all 4,950 pairs of nontrivial monomial orbits
through degree six under the Dvorsky-aligned four-variable involution; odd
moments vanish automatically, but the exact second/fourth-moment gcd has no
nonzero root on any support.  Finally, projective rational Gröbner bases show
that for homogeneous binary forms of degrees one through four, the first
`d+1` moments force the zero form and the cutoff is sharp.  The formulas,
explicit sharp quadratic and cubic witnesses, and strict finite-search scope
are in
[`SPARSE_FACTORIAL_MOMENT_FRONTIER.md`](extended-geometry/SPARSE_FACTORIAL_MOMENT_FRONTIER.md).

The bounded four-variable descendant search can be replayed separately:

```bash
python3 scripts/search_dvorsky_gvc4_bounded.py
```

It exhausts the declared \(40\)-by-\(7{,}448\) normalized lattice slice,
checks pure contractions through order twelve, and screens the full space
of fixed linear multipliers on orders five through twelve.  It finds no
witness in that slice; this is a finite negative search, not a GVC(4) or
SIC(4) theorem.

One natural full-coefficient slice is now closed exactly:

```bash
.venv/bin/python scripts/verify_four_pair_dvorsky_slice_obstruction.py
```

For \(P=(t+a+b+d)(ad+bt)\) and the general ternary quadratic symbol
\(\Lambda=\partial_tR(\partial_a,\partial_b,\partial_d)\), the checker
constructs the first eight pure contractions, proves by exact Singular
reductions that their projective zero set consists of four rank-one square
directions, and verifies a strict weight gap for every direction.  The
written proof in
[`FOUR_PAIR_DVORSKY_SLICE_OBSTRUCTION.md`](extended-geometry/FOUR_PAIR_DVORSKY_SLICE_OBSTRUCTION.md)
upgrades this slice from a bounded lattice experiment to an all-order
no-counterexample theorem for arbitrary fixed multipliers.  It does not
prove GVC(4) or SIC(4).

The all-order nonvanishing proof is written in
[`IMAGE_VANISHING_COUNTEREXAMPLES.md`](extended-geometry/IMAGE_VANISHING_COUNTEREXAMPLES.md);
the generators check the finite artifacts and change-of-variable identities;
the dependency-free audit re-expands the 40-variable witness from scratch.
The dependency chain is `F1 -> LR1 -> IV1`: the foundational collision feeds
the essential cubic quotient, whose named inverse coordinate gives the direct
Image witness, identity slice, and the 40/42-variable Vanishing witnesses.
The parallel `LR1 -> GS1` branch is instead the nonexplicit route to
`not GMC(42)`.  The quantitative rank branch is
`LR1 -> LR2 -> LR3`, with `IV1 -> LR3` supplying the HN consequence
framework for the separate rank-37 realization.  These arrows record logical
or construction dependence.
The collision route retains witness sizes 20/40/42 and rank 37.
Historically, the independent Dvorsky--Long formulas lowered the certified
SIC and unrestricted GVC entries to 5/5.  The current repository witnesses
lower these to SIC pair-dimension 2 and unrestricted GVC dimension 3; binary
GVC makes the latter exact.  The ordinary-Laplacian 40 and homogeneous HN 42
entries are unchanged.  These are witness-ledger values, not literature-wide
minimality claims.
A local proof of the
fixed-dimensional DVEZ/Zhao implication, including Gaussian contraction, the
countable-union step, and formal inversion, completes the nonexplicit route to
`not GMC(42)`; `not GMC(158)` remains the exact conservative Long-route bound.
These high-dimensional bounds are retained as logical-transport regressions,
not active witness searches: Long's five-term three-real example already
settles all dimensions `n>=3`.
It also verifies the uniform weighted-seed Gaussian bridge:
first the standalone Gaussian--Lagrange identity for a nonlinear polynomial
map with nonzero constant terms, then the exact pencil branch, polynomial
determinant correction, and bounded Wick moments for canonical and split
seeds.  It also reverts the mixed-moment generating series to recover a
symbolic quartic and a concrete weighted quintic exactly, verifies the
optimal `N-3` normalized moment coordinates through degree eight and the
variable-scale `N-2` bound, and checks the determinantal reciprocal-series
equations, followed by a separate standard-library reconstruction.  The
all-order completed-ring and residue proof is
[`FORMAL_GAUSSIAN_LAGRANGE_LEMMA.md`](extended-geometry/FORMAL_GAUSSIAN_LAGRANGE_LEMMA.md);
the bounded exact script is explicitly a regression rather than a substitute
for that proof.  These checks are part of
`verify-regressions`, not `verify-minimal`.

The first nontrivial exceptional partition complex has a separate exact
moment-coordinate certificate:

```bash
.venv/bin/python scripts/verify_degree_six_gaussian_moment_geometry.py
.venv/bin/python scripts/verify_moment_prony_determinantal_geometry.py
```

The first command derives the irreducible sextic equation of the all-double
component, parametrizes the all-triple curve, and verifies that their four
all-six collision points have scheme-theoretic intersection length two.  It
also transports both degree-six vertical Ritt hypersurfaces into optimal
moment coordinates, verifies the second displayed sextic, and proves that the
`2 o 3` Ritt surface is exactly the all-double exceptional component.  The
second command constructs the equal-multiplicity loci and both degree-eight
Ritt orders from log-Prony and Krylov minors in the optimal moments, compares
their ideals scheme-theoretically, retains the degree-six dual-number
intersection, proves the degree-eight Ritt intersection reduced, and exhibits
the cubic collision thickness in naive mixed-weight Fitting minors.  It then
replaces that marked scheme by the saturated unmarked
Christoffel--Hankel/subresultant ideal, covers the minimal degree-five
`3+2` case and one-node collision strata, and verifies the length-two
degree-eight mixed/all-double
intersection.

The general primitive-merger theorem and its first failure have a separate
exact certificate:

```bash
.venv/bin/python scripts/verify_omitted_intersection_algebra.py
```

It constructs the allocation hypergraph and merger-cycle spaces for the
degree-twelve, degree-eighteen, and first degree-twenty-four faces.  It also
derives the primitive dual-number block, finds the first nonminimal failure
`k[t]/(t^3)` when its root meets a common double atom in degree eight, and
proves that coalescing two pure transfer blocks gives
`k[X,Y]/(X^3,XY,Y^2)`, which has the same length and Hilbert vector as two
dual numbers but a two-dimensional socle.

The underlying Hessian-incidence and Ritt-intersection calculation is replayed by

```bash
.venv/bin/python scripts/verify_hessian_ritt_degree_six.py
.venv/bin/python scripts/verify_degree_six_ritt_atlas.py
```

The second checker refines the Hessian-incidence result on the normalized
seed chart: it computes the `2^3` and `3^2` omitted-value intersections, the
four doubled type-`(6)` collision points, factored affine-sheet boundary cuts,
and clean rational witnesses for all open pieces.

The complete degree-six boundary atlas requires both SymPy and Singular:

```bash
make verify-ritt-boundary
```

It proves that the two Ritt surfaces have respectively two and three exact
affine-boundary curves, supplies a rational Hessian-clean witness on every
curve, and computes the common-curve deletions: one reduced sextic
zero-cluster orbit, two rational plus four conjugate extra-root points, and
four Hessian/type-`(6)` collisions disjoint from the affine boundary.

The first genuine braid of complete decompositions has a separate
scheme-theoretic certificate:

```bash
make verify-ritt-2-complex
```

It builds the Ritt Coxeter 2-complex with commuting-square and braid
relations, verifies the Dickson coefficient map at all six degree-thirty
vertices, and compares the two path ideals around the `S_3` hexagon.  Both
paths have the same smooth `A^2` reduction and normalization.  One path is
reduced; the other has nilpotence index four, with one excess tangent
direction and normalization-defect annihilator `(z^2)`, supported on the
monomial divisor.  It also identifies dual-number and length-five
curvilinear slices of the defect, computes the latter's `K`-adic length
filtration `2,4,5,5`, and verifies that the path tangent dimensions are
unchanged when computed directly in the ambient polynomial and Hessian
coefficient spaces.  The full ideal and doubled-annihilator comparison is
then repeated independently on the opposite `5 o 3 o 2` endpoint chart.
Restoring the omitted linear-coefficient residual leaves every endpoint,
path, and boundary ideal unchanged on both charts, proving exact
scheme-theoretic Hessian transfer for this braid component.
The checker then audits the four remaining vertex charts.  The three
composite-omission sectors `10`, `15`, and `6` have respectively
nilpotence/annihilator data
`(4,z^2)`, `(3,z^2)`, and `(4,z^4)`; the complementary prime-omission path
is reduced in every sector, and opposite endpoint charts agree.
Their annihilator slices have
`(length, embedding dimension, Hilbert vector)` equal to
`(5,1,(1,1,1,1,1))`, `(4,2,(1,2,1))`, and
`(8,2,(1,2,2,2,1))`.  All three have one-dimensional socle; the latter two
are codimension-two Artin Gorenstein complete intersections.  Exact
elimination identifies the three slice algebras as
`Q[u]/(u^5)`, `Q[u,v]/(u^2,v^2)`, and
`Q[u,v]/(u^4,v^2)`.  Their conormal ranks are `1,2,2`, with residue-field
Koszul Tor ranks `(1,1)`, `(1,2,1)`, and `(1,2,1)`.

The coefficient-decorated cellular and Postnikov generalization is replayed
by

```bash
.venv/bin/python scripts/verify_hessian_ritt_cellular_cotangent_prototype.py
.venv/bin/python scripts/verify_degree42_ritt_conormal_transitivity.py
.venv/bin/python scripts/verify_degree42_ritt_postnikov_overlap.py
.venv/bin/python scripts/verify_cellular_postnikov_transitivity.py
.venv/bin/python scripts/verify_hessian_ritt_cotangent_descent.py
.venv/bin/python scripts/verify_ritt_cellular_prototype_completion.py
```

The first command verifies the vertex, move, commuting-cell, braid-cell, and
relative-path totalizations.  The second and third use Singular to prove
degree-forty-two conormal non-splitting, overlap vanishing, and the separation
of non-flat base-change Tor.  The last command is fast exact rational linear
algebra: it validates arbitrary finite equivariant module towers, replays the
degree-thirty one-layer degeneration and the actual degree-forty-two
base-square action matrices, and writes
`artifacts/generated-results/cellular_postnikov_transitivity.json`.
The final command verifies the skeletal descent boundary: the filled braid
is complete in dimension two, while the oriented permutohedron three-cell
kills the topological `H2` line of the four-factor Coxeter two-skeleton
without changing `H0` or `H1`.  It also constructs the normalized
face-poset bars and canonical subdivision maps for the relative half-braid
and filled braid.  Their mapping cones are exactly acyclic before and after
tensoring with coefficient blocks of dimensions `2`, `4`, and `6`; the
same run computes the actual degree-forty-two tangent images of all six
factor charts in Hessian coefficient space.  The vertex ranks are all nine,
the adjacent-move intersection ranks are `(8,5,6,6,5,8)`, and the common
intersection is the Dickson tangent plane plus the Hessian projection of
`(W+1)^36-1`.  Intersecting the four vertex images along each half-braid
then verifies the conormal flag `(5,6,6,7)` for all three opposite-pair
sectors, with composite omissions `6`, `14`, and `21`.  The result is
written to
`artifacts/generated-results/hessian_ritt_cotangent_descent.json`.

The two nonlinear rotated degree-forty-two sectors are expensive,
specialized computations:

```bash
.venv/bin/python scripts/explore_degree42_ritt_rotated_conormal_flags.py --word 237
.venv/bin/python scripts/explore_degree42_ritt_rotated_conormal_flags.py --word 327
.venv/bin/python scripts/verify_degree42_ritt_cut14_postnikov_overlap.py
.venv/bin/python scripts/verify_degree42_ritt_cut21_postnikov_overlap.py
.venv/bin/python scripts/verify_degree42_ritt_cut14_tensor_split_q4.py
.venv/bin/python scripts/verify_degree42_ritt_inverse_limit_sections.py
.venv/bin/python scripts/verify_degree42_ritt_completed_splits.py
```

They reconstruct the cut-`14` and cut-`21` residual ideals, change to seven
normal plus two Dickson-base coordinates, and compute the exact quotient
modulo the fourth maximal-ideal power.  Both have conormal flag
`(5,6,6,7)`, thin path equal to boundary through order four, and common
spectator dimensions `(1,3,6)`.  Their sector dimensions are respectively
`(1,4,9)` and `(1,4,10)`, rather than the existing cut-`6` profile
`(1,5,13)`.  The commands write
`artifacts/generated-results/degree42_ritt_rotated_conormal_jet_237.json`
and
`artifacts/generated-results/degree42_ritt_rotated_conormal_jet_327.json`.
The first two reconstructions also write reusable compressed ordinary-ideal
caches
`artifacts/generated-results/degree42_ritt_rotated_source_ideals_237.json.gz`
and
`artifacts/generated-results/degree42_ritt_rotated_source_ideals_327.json.gz`.
To reconstruct only one cache from scratch, add
`--rebuild-source --build-source-only`; this replaces the selected cache
and does not run the local Singular audit.  The last two commands consume
those caches and perform exact Nakayama and Artin--Rees tests.  They prove
completed thin-path/boundary equality and vanishing of the completed
quadratic overlaps for cuts `14` and `21`, and write the corresponding
`degree42_ritt_cut14_postnikov_overlap.json` and
`degree42_ritt_cut21_postnikov_overlap.json` artifacts.
The two rotated-conormal JSON files remain finite-jet computations; the
separate completed certificates cover both new sectors.
The order-four command consumes the compressed order-four action-matrix cache
and proves that the tensor-presented cut-`14` conormal extension splits
over `B/(tau,zeta)^4`, with dimensions `9 -> 13 -> 4` and cocycle ranks
`32=32`.

The next command consumes the two compressed order-seven action-matrix
caches and constructs orders five and six as quotients of those single
order-seven presentations.  Thus its sections commute with truncation by
construction.  For both cuts `14` and `21`, the dimensions at orders
`5,6,7` are respectively

```text
12 -> 17 -> 5
15 -> 21 -> 6
18 -> 25 -> 7
```

and the cocycle/coboundary ranks are `55=55`, `84=84`, and `119=119`.
The section-difference restriction maps have two-dimensional cokernels, so
these finite splits alone do not imply an inverse-limit split.

The last command consumes the completed two-variable presentation caches
and verifies explicit polynomial sections.  With `u=1+tau`, their generator
images are

```text
cut 14: e4 + (-3*u^2 + 2*zeta)*e6
cut 21: e4 + (-4*u^3 + 8*u*zeta)*e7.
```

Every spectator relation maps to zero in the total presentation, and each
section followed by projection is the identity.  Hence both extensions
split over `Q[[tau,zeta]]`; the completed extension and inverse-limit torsor
obstruction classes vanish.  The earlier extra cut-`21` fourth-jet
quotient-ring dimension is a non-flat base-change/Tor discrepancy, while
the different correction polynomials retain the genuine labelled-sector
asymmetry.  Restriction coherence between different factor charts remains
open.

The last command writes the fully explicit degree-forty-two and
degree-thirty factor/move/labelled-cell diagrams, the totalized complexes,
all certified filtration cohomology rows, and the exact first failed split
reduction.  It verifies uniform cellular `H2=0` for the prototype while
locating the genuine obstruction at degree-forty-two filtration order three:
the sector--spectator extension is non-split and its completed cotangent
connecting morphism is nonzero.  The result is
`artifacts/generated-results/ritt_cellular_prototype_completion.json`.

The direct conductor-first node/cusp ansatz is replayed by

```bash
.venv/bin/python scripts/verify_conductor_first_one_chart_obstruction.py
```

It constructs finite quadratic marked-root algebras over
`Q+t(t-1)Q[t]` and `Q+t^2Q[t]`, verifies discriminant descent through the
node and cusp conductors, and solves the complete reconstruction
polynomiality equations against a prescribed conductor pole.  The systems
are incompatible in every pole order by divisibility; exact Gröbner
regressions cover the first four orders.  The independent unit-group
obstruction shows that the sole-conductor localization is not affine space
after any polynomial stabilization.  The theorem is scoped to the direct
separated one-chart ansatz; multi-chart and ambient-coupled conductors remain
open.  The result is
`artifacts/generated-results/conductor_first_one_chart_obstruction.json`.

The smallest symmetric ambient-coupled escape is replayed by

```bash
.venv/bin/python scripts/verify_conductor_three_boundary_cox_fill.py
```

It replaces the conductor localization by the three-prime fill
`x*y*z=c(t)`, with `c=t(t-1)` for the node and `c=t^2` for the cusp.
Singular verifies the descended hypersurface equations, smoothness of the
nodal normalization, and the exact three-axis singular locus of the
cuspidal normalization.  The remaining exact checks show that the marked
root and reconstruction equations descend polynomially, but the descended
dualizing form pulls back as `Omega_norm/(x*y*z)`.  Affine-space recognition
then fails: the smooth nodal fill has Hodge--Deligne polynomial
`(uv)^3+2(uv)^2-uv`, while the normal cuspidal fill is singular.  This is an
obstruction for the symmetric Cox-product ansatz, not for asymmetric affine
modifications or a distributed target conductor ledger.  The result is
`artifacts/generated-results/conductor_three_boundary_cox_fill.json`.

The general affine boundary-obstruction regressions are:

```bash
make verify-boundary-obstruction-theory
```

The checker compiles one saturated module and one genuine boundary-torsion
module, including regular-element and distinguished-class certificates.  It
then computes the finite normal jets `Q[x]/(x^n)` for `n=1,...,6`: every
transition is surjective, but the least boundary-annihilation exponent grows
as `1,2,3,4,5,6`.  This is the exact control showing why bounded jet
saturation does not imply formal saturation without a uniform exponent.
Independent rational matrices verify the node and cusp conductor pullbacks,
rank-three finite-free tensor descent, strict bounded lifting, and a
non-strict degree-loss example.  The command writes
`artifacts/generated-results/boundary_obstruction_theory.json`.

The conductor-first existence certificate is:

```bash
.venv/bin/python scripts/verify_conductor_first_foundational_cusp_keller.py
```

It begins with `Q[u^2,u^3] subset Q[u]`, derives the translated cusp
`4*S^3+27*V^2=0`, and reconstructs the cubic seed
`H(W)=W^2(1-W)`.  The cubic marked-root discriminant descends through the
conductor and has a triple root at its conductor point.  The checker then
solves the full weighted source/target ledger: all apparent divisions cancel
to `(F1,F2/2,F3/2)`, the reconstruction is `x=-C/E_W`, and the exact
Jacobian is `-1/2`.  Three rational source points have common target
`(-1/4,0,0)`.  This meets the conductor-first existence criterion but
recovers the known foundational weighted mechanism rather than a new stable
class.  The result is
`artifacts/generated-results/conductor_first_foundational_cusp_keller.json`.

The characteristic-labelled Hessian--Ritt linear complex and its
positive-characteristic Frobenius summand are replayed together by

```bash
make verify-ritt-deformation-complex
```

The first checker verifies the characteristic-zero tree differential and
cellular baseline.  The second exhausts normalized one-sided composition
tangents over small fields and verifies the exact Hessian-cutoff trichotomy:
the full \(r-1\)-dimensional kernel for \(H'=0\), the one-dimensional
\(kx\) kernel for \(H'=a\ne0\), and no invisible tangent when
\(\deg H'\ge1\).  The exhaustion is a regression for the written
all-degree degree-additivity proof; it is not itself the proof.

The chart-independent missing-linear-coefficient test is:

```bash
make verify-hessian-synchronization
```

It constructs the canonical lift \(\lambda_{a,b}\) using only Hessian
coefficients.  Exact ambient and canonical-factor-chart ideal membership
proves every multiple intersection through degree `18` is synchronized
scheme-theoretically.  In degree `24`, fourteen pairs reduce directly on
canonical factor charts.  The final outer-cut pair `{2,3}` is certified after
transporting the degree-six Dickson collision through a generic quartic and
changing to `4 normal | 5 base` coordinates; its exact Groebner basis has
size `63`.  Thus every degree-`24` multiple intersection is synchronized, and
each ordinary polynomial intersection is exactly one graph over its Hessian
intersection.  The same target verifies the augmentation-ideal lengths,
point-cotangent homology, and intrinsic Tor ranks of the degree-thirty
transverse sector models.  Finally, five exact degree-`30` pair reductions
with basis sizes `11,6,95,6,11` form the cut spanning tree
`2-6-3-15-5-10`.  Therefore the global all-six degree-thirty intersection is
scheme-theoretically synchronized.  A `4 normal | 7 base` common-refinement
calculation also closes the nested pair `{2,10}` with basis size `4`; five
incomparable two-cut subintersections remain uncertified.

The four larger non-tree pair certificates are intentionally split from the
fast spanning-tree regression:

```bash
make audit-degree30-hessian-synchronization-pairs
```

They certify `{5,6}`, `{6,10}`, `{6,15}`, and `{10,15}` with exact basis
sizes `502,189,12,96`.  Together with the default target, ten of the fifteen
degree-thirty pairs are therefore certified.

The same target runs the rank-two Poisson pre-audit and the independent
completion certificate.  The first verifies that the
single displayed output `R=x(2-3xq)` is exactly the foundational third output
after a polynomial source automorphism, and proves that the naive choices
`S=F_1/2`, `T=F_2` have no polynomial `D`-completion.  The second derives the
pole-cancelling shear `Z -> Z-9Q^2`, constructs exact polynomial `T,D,S`,
checks all six brackets and determinant one, and transports the complete
three-point fiber.  A dependency-free sparse-polynomial implementation then
rebuilds the formulas and separately checks all six brackets, the determinant,
term counts, and collision.  This proves a repository rank-two Poisson
theorem; it does not assert that these are the unavailable manuscript's
formulas.

The all-degree rigidity step behind the transported Hessian cases is
replayed by

```bash
make verify-common-right-factor-synchronization
```

It verifies the triangular top-jet reconstruction for every common-right
degree occurring in degrees `30` and `42`, checks that the two degree
censuses each have exactly three decorated incomparable pairs, and verifies
the characteristic-two dual-number counterexample when the total outer
degree is not invertible.  The theorem itself works over every ring in which
that outer degree is a unit.

The first degree-`42` primary transport certificate is:

```bash
make verify-degree42-hessian-normal-jets
```

On the `{2,7}` pair it constructs the `5 normal | 6 base` common-cubic
power chart and proves over `QQ` that the synchronization defect belongs to
the Hessian residual ideal plus the fifth power of the normal ideal.  Thus
the full six-parameter component synchronizes through normal order four.
The exact basis has size `88`.

The conceptual all-order upgrade on the dense power chart is replayed by

```bash
make verify-degree42-conormal-rees-synchronization
```

The normal Jacobian has maximal-minor ideal `(w0^2)`.  Hence away from
`w0=0` the residual conormal map is onto, complete Nakayama identifies the
completed residual and normal ideals, and the synchronization defect
vanishes at every Rees order.  On `w0=0` the conormal rank is exactly
three, so the existing fourth-order certificate remains the correct global
statement and the all-order primary frontier is confined to that divisor.

The divisor itself has a two-normal-variable Rees reduction:

```bash
make verify-degree42-divisor-rees-reduction
```

Three unit residual pivots eliminate `x3,x4,x5`.  The remaining binary
quadrics have resultant
`(81/256)*w1^4*((t+e1*e2)^2-4*e1^3)`.  Off this resultant their Hilbert
vector is `(1,2,1)`, so the completed normal ideal has cube zero and the
fourth-order defect certificate becomes exact.  The unresolved all-order
locus is reduced to `V(w0,w1)` together with
`V(w0,(t+e1*e2)^2-4*e1^3)`.

Dense opens of both residual branches are closed by

```bash
make verify-degree42-kuranishi-branches
```

On the discriminant branch, normalization exposes one common quadratic
tangent; its cubic obstruction is nonzero on `D(w1*w2*t)`, giving the
initial ideal `(ell^2,ell*s,s^3)`.  On `w1=0`, the two binary cubic
Kuranishi forms have resultant
`-15625/262144*w2^6*A*B`, with `A,B` displayed in the canonical note.
On `D(w2*A*B)` their complete-intersection Hilbert vector is
`(1,2,3,2,1)`.  In both cases the existing membership modulo the fifth
normal power is therefore exact.  The same checker computes the first
subresultant on the exceptional divisors: generically the cubics share
exactly one explicit linear factor on each of `A=0` and `B=0`.  Their
next obstruction is consequently a one-variable quartic restriction,
whereas, on this `w0=w1=0` branch, the further equation `w2=0` is the
sevenfold monomial collision.

The degenerate part of the discriminant branch is closed by

```bash
make verify-degree42-discriminant-quartics
```

When `t=0`, the common-tangent cubic vanishes but the terminal quartic
coefficient is `5*w2/64`.  This remains true at the cusp `e1=t=0`,
where the quadratic ideal becomes a single square.  Consequently the
whole discriminant branch synchronizes on `D(w1*w2)`.  The remaining
support is only `V(w0,w2)` together with `V(w0,w1,A*B)`.

Geometrically, for `W(z)=z^3+w2*z^2+w1*z+w0` and
`U(z)=z*W(z)^2`, the generic `V(w0,w2)` core is the odd polynomial
`z^3*(z^2+w1)^2`, while `V(w0,w1)` has the contact-five core
`z^5*(z+w2)^2`; only their deepest intersection is the sevenfold monomial
`z^7`.  On the `A` and `B` divisors the unique common tangent line is
handled by the common-line residual-intersection theorem: blow up the
normal plane, eliminate the transverse coordinate, and read the first
nonzero coefficient of the resulting one-variable residual series.

The generic `A/B` quartic restrictions are closed by

```bash
make verify-degree42-ab-residual-quartics
```

At the exact characteristic-zero point
`(e1,e2,t,w2)=(1,1,3/5,1)` on `A=0`, the residual scalar is
`-4203/1280`.  At the good-prime point `(1,1,21,1)` on `B=0` modulo
`103`, it is `47`; here `A=1` and the subresultant coefficient
`alphaB=9`.  The latter nonzero reduction excludes an identically zero
characteristic-zero restriction on the irreducible divisor `B`.  Hence
both generic resultant divisors synchronize.  Only the proper quartic-zero
subloci, together with the odd core `V(w0,w2)`, can remain.

The part of the proper quartic-zero analysis supported on the higher-gcd
locus is resolved by

```bash
make verify-degree42-higher-gcd-strata
```

The reduced higher-gcd locus on `D(w2)` is the union of four weighted
curves.  Quartic envelopes close every punctured curve: their Hilbert
vectors are `(1,2,3,2)`, except for one `(1,2,3,3,1)` curve.  The only
point of this locus not closed by the nilpotence cutoff is their common
contact-five vertex `e1=e2=t=0`; there the cubics vanish and the quartics
retain a common cubic factor.  This does not yet exclude additional zero
divisors of the scalar quartic residual away from the higher-gcd locus.

Those scalar residual divisors are factored by

```bash
make verify-degree42-ab-residual-factors
```

On the rational normalization of `A=0`, the residual is
`-75/512*e1^2*w2*(4*e1-e2^2)*P_A/
(e2*(6*e1-e2^2)^3)`, where `P_A` is affine-linear in `w2`.  On `e2=1`,
the normalization of `B=0` over `q^2=-3` is
`e1=(1-r+q*(r-1)*(2*r-1))/2`,
`t=-(1+q)*(r-1)^2*(2*r-1)/2`; there the residual is
`75/1024*w2*(q-1)*(r-1)^2*(2*r-1)*P_B`, with `P_B` also affine-linear
in `w2`.  The `e1=0`, `4*e1=e2^2`, `r=1`, and `r=1/2` factors are
exactly the already-classified higher-gcd branches `P4` and `P3`.
Thus the new degree-one-gcd support consists only of the residual graphs
`P_A=0` and `P_B=0`.  The checker also certifies the solved graph
identities for `w2` and verifies that the apparent coefficient-zero
values `5*e1=e2^2` and `r=3/5` do not add vertical components.

Dense opens of both residual graphs are closed by

```bash
make verify-degree42-ab-residual-quintics
```

The checker follows one terminal equation through its quartic transverse
correction, re-solves the three pivot equations through fourth order, and
computes the invariant fifth residual.  On the `A` graph its exact value is
`-250011279/8192000` at
`(e1,e2,t,w2)=(1,1,3/5,567/100)`.  On the geometric `B` normalization it
is `14 mod 103` at `(q,r)=(10,0)`.  The resulting homogeneous envelope has
Hilbert vector `(1,2,3,2,1)`, so the pre-existing defect membership modulo
the fifth normal power becomes exact on both dense opens.

The complete cutoff chain—from the conormal open through the higher-gcd
quartic strata and residual-graph quintics—can be replayed with

```bash
make verify-degree42-kuranishi-cutoff-chain
```

The two pieces are combined by the single support ideal

```text
k = (w0, w1*w2, A*B*w2).
```

The global non-jet target is `I:k^infinity = I`.  By the shared
[support-saturation principle](verified/SUPPORT_SATURATION_PRINCIPLE.md),
this is equivalent to excluding associated primes of the residual algebra
over `V(k)` and is sufficient after normal completion.  It is strictly
weaker than proving the residual algebra flat over the full Ritt base.

The first exact compression of this target is:

```bash
make verify-degree42-depth-reduction
```

Residuals 5, 11, and 17 are global unit-triangular pivots, not only pivots
on `w0=0`.  They eliminate `x3,x4,x5` exactly and present the same residual
algebra using only `x1,x2` over the six-dimensional base.  The checker also
has exploratory `--method height` and `--method colon` modes for

```text
f = w0 + w1*w2 + A*B*w2.
```

The height mode searches for a codimension-two perfect reduced ideal and a
one-step dimension drop after adjoining `f`; the colon mode tests `I:f=I`
directly.  Neither mode is part of the verified target until its
characteristic-zero computation completes.

The next normal order has an exact good-prime certificate, and one rational
point on the remaining divisor has an exact untruncated characteristic-zero
certificate:

```bash
.venv/bin/python scripts/verify_degree42_transported_27_normal_jets.py \
  --prime 32003 --normal-order 5 --timeout 240
.venv/bin/python scripts/verify_degree42_transported_27_normal_jets.py \
  --base-values 1,2,3,0,5,6 --normal-order 0 --timeout 240
```

They return basis sizes `179` and `8`, respectively.  To attempt the
remaining generic characteristic-zero calculation directly on the
conormal divisor, use

```bash
.venv/bin/python scripts/verify_degree42_transported_27_normal_jets.py \
  --w0-zero --normal-order 5
```

That function-field calculation currently exceeds 300 seconds; the timeout
is a performance boundary, not a failed reduction.

The generated certificates are stored as the conservative
[`79-variable artifact`](artifacts/generated-results/long_bcw_79_counterexample.json)
and optimized
[`33-variable artifact`](artifacts/generated-results/shared_bcw_33_counterexample.json),
together with the
[`24-variable rank-compressed artifact`](artifacts/generated-results/rank_compressed_bcw_24_counterexample.json)
and final
[`22-variable constant-kernel quotient`](artifacts/generated-results/constant_kernel_bcw_22_counterexample.json),
together with the new
[`21-variable essential quotient`](artifacts/generated-results/essential_bcw_21_counterexample.json)
and its
[`20/40-dimensional identity-slice witnesses`](artifacts/generated-results/image_vanishing_counterexamples_20_40.json)
and
[`homogeneous 21/42-dimensional witnesses`](artifacts/generated-results/image_vanishing_counterexamples_21_42.json).
They record the sparse cubic maps, every reduction-step choice, and the three
exact collision points, together with the expanded contraction and quartic
polynomials; regeneration is deterministic.

## Complete active suite

```bash
make verify
```

## Degreewise theorem audit

The former standalone manuscript has been retired. The retained theorem
statement, proof dependencies, and reproduction commands are in
[`DEGREEWISE_MULTIPLICITY_AUDIT.md`](DEGREEWISE_MULTIPLICITY_AUDIT.md).

To compile every standalone paper with the same discovery rule used by CI,
run:

```bash
make verify-papers
```

To retain an environment record and complete log under `artifacts/`, run:

```bash
make verify-logged
```

Generated outputs, bounded scans, and exploratory search programs are not
part of the public proof navigation.  Existing generated artifacts live under
`artifacts/generated-results/`; historical search tools are preserved under
`archive/tooling/`.

## Cyclic and dihedral absolute inverse-Galois audit

The first Programme 3 checker replays the Dickson power-sum recurrence,
primitive derivative, reduced-branch pullback, and odd/even discriminant
formulas through degree twelve.  It also checks the low-degree \(D_3,D_4,D_5\)
cards and the determinant-minus-one derivative-unit suspension:

```bash
.venv/bin/python scripts/verify_cyclic_dihedral_keller_audit.py
```

The all-degree proof and the distinction between geometric and arithmetic
monodromy are in
[`extended-geometry/ABSOLUTE_INVERSE_GALOIS_CYCLIC_DIHEDRAL_AUDIT.md`](extended-geometry/ABSOLUTE_INVERSE_GALOIS_CYCLIC_DIHEDRAL_AUDIT.md).
The bounded replay is a regression certificate, not an exhaustive proof over
all degrees.

## Absolute \(D_5\) affine-modification frontier

The degree-five precomputation checks the split derivative and branch
ledgers, the tangency of the two ramification colors, the singular product
and separated Cox fills, and the maximal-minor determinant identity for
affine-linear couplings with two and three new coordinates.  It also checks
the block determinant for arbitrary zero-section thickenings:

```bash
.venv/bin/python scripts/verify_absolute_dihedral_d5_modification_frontier.py
```

The all-degree affine-linear mask theorem, class-group calculations, and
the exact surviving nonlinear search conditions are in
[`extended-geometry/ABSOLUTE_DIHEDRAL_D5_MODIFICATION_FRONTIER.md`](extended-geometry/ABSOLUTE_DIHEDRAL_D5_MODIFICATION_FRONTIER.md).

## Nonlinear \(D_5\) obstruction classification

The nonlinear follow-up checks the complete branch-supported valuation
ledger, its primitive diagonal pole, the normalized target-cusp incidence
and tangent-rank drop, and the translated graph-section obstruction for two
and three auxiliary coordinates:

```bash
.venv/bin/python scripts/verify_d5_nonlinear_obstruction_classification.py
```

The eight-gate classification, proofs, and exact rank-five fibre
certificate requirements are in
[`extended-geometry/D5_NONLINEAR_MODIFICATION_OBSTRUCTION_CLASSIFICATION.md`](extended-geometry/D5_NONLINEAR_MODIFICATION_OBSTRUCTION_CLASSIFICATION.md).

## Canonical \(D_5\) two-mask blowdown

The first construction attempt verifies the determinant-\(\Delta\)
two-mask matrix and its adjugate inverse, proves the constant-linear
coefficient locks, checks the generic genus-two fibre used in the all-degree
automorphic rigidity theorem, and exhausts \(72\) coordinate assignments
from the unchanged and two ramification-incidence source charts.  It also
checks the first nonautomorphic contraction mismatch and solves the minimal
affine-normal tangential class:

```bash
.venv/bin/python scripts/verify_d5_two_mask_blowdown_obstructions.py
```

The chain-rule obstruction, all-degree genus-two proof, and precise
nonautomorphic continuation are in
[`extended-geometry/D5_TWO_MASK_BLOWDOWN_OBSTRUCTIONS.md`](extended-geometry/D5_TWO_MASK_BLOWDOWN_OBSTRUCTIONS.md).

## All-degree dihedral affine-completion obstructions

The uniform replay checks the Dickson branch identity, both components of
the even-degree branch, the determinant-\(\Delta_n\) two-mask blowdown, the
positive-genus automorphic-rigidity gate, the first nonautomorphic cusp
remainder, the affine-normal coefficient locks, and the nonlinear
normal-degree resonance, even factorization, and odd valuation-at-infinity
one-normal no-go gates for \(3\le n\le12\).  It also replays two
resonant-looking false positives and verifies that the bounded
one-normal search compiler has no open route:

```bash
.venv/bin/python scripts/verify_dihedral_all_degree_affine_completion_obstructions.py
```

The odd/even valuation ledgers and uniform proofs are in
[`extended-geometry/DIHEDRAL_ALL_DEGREE_AFFINE_COMPLETION_OBSTRUCTIONS.md`](extended-geometry/DIHEDRAL_ALL_DEGREE_AFFINE_COMPLETION_OBSTRUCTIONS.md).
The bounded replay is a regression certificate, not the proof of the
all-degree statements.

## \(S_4\) collision-frame Keller frontier

The collision-frame checker expands the decomposable absolute map
\(F\circ F\), verifies its determinant and degree-nine tower, and checks the
quartic factorization, six-edge action, discriminant, primitive conductor,
tiny normal form, and one-/two-normal determinant ledgers.  It also verifies
the rational determinant-one cotangent lift, the polynomial relative
logarithmic factorization, the natural resultant and Bezout obstruction
models, and finite-field point-count regressions for the displayed motivic
classes:

```bash
.venv/bin/python scripts/verify_s4_collision_frame_keller_frontier.py
```

The class-group calculations, the all-degree rank argument, and the
affine-space recognition requirements are written proofs in
[`extended-geometry/S4_COLLISION_FRAME_KELLER_FRONTIER.md`](extended-geometry/S4_COLLISION_FRAME_KELLER_FRONTIER.md).
The checker does not construct an ordinary polynomial Keller map for the
six-sheet collision cover.

## \(A_4\) Keller inverse-Galois frontier

The pure-target ledger, two-mask factorization, normalized-boundary
assembly, and root-incidence derivative-split checks are:

```bash
.venv/bin/python scripts/verify_a4_pure_target_ledger.py
.venv/bin/python scripts/verify_a4_two_mask_factorization.py
.venv/bin/python scripts/verify_a4_normalized_boundary_assembly.py
.venv/bin/python scripts/verify_a4_root_incidence_derivative_split.py
.venv/bin/python scripts/verify_a4_chart_unit_rank_four.py
.venv/bin/python scripts/verify_a4_two_mask_local_viability.py
.venv/bin/python scripts/verify_a4_affine_modification_obstruction.py
.venv/bin/python scripts/verify_a4_corrected_boundary_selector.py
.venv/bin/python scripts/verify_a4_boundary_coloring_surgery.py
Singular -q scripts/verify_a4_corrected_boundary_genus.sing
.venv/bin/python scripts/verify_a4_genus_zero_selector_search.py
Singular -q scripts/verify_a4_genus_zero_selector_search.sing
.venv/bin/python scripts/verify_a4_sharp_selector_plane.py
.venv/bin/python scripts/verify_a4_conic_principal_obstruction.py
.venv/bin/python scripts/verify_a4_degree_twenty_line_selectors.py
.venv/bin/python scripts/verify_a4_degree_twenty_line_selectors.py --conic-sieve
.venv/bin/python scripts/verify_a4_degree_twenty_line_selectors.py --cubic-sieve
.venv/bin/python scripts/verify_a4_degree_twenty_line_selectors.py --census-bound 6
.venv/bin/python scripts/verify_a4_degree_twenty_line_selectors.py --census-bound 3 --include-q5
Singular -q scripts/verify_a4_degree_twenty_near_selector.sing
```

The normalized-boundary command verifies the determinant-one ambient
completion and the exact obstruction to the resulting automorphic assembly.
It also verifies the first birational nonautomorphic solution of the corrected
log-Jacobian equation, namely homogeneous radial scaling, and the abstract
Jacobian identity used to prove that no choice of two polynomial mask outputs
can repair that radial base after the pure-target lift.  It also verifies
the determinant-boundary matrix factorization organizing the two inverse
masks and the direct-adjugate Jacobian obstruction.  The accompanying note
identifies the actual descended coupling space as the ideal-contraction
quotient \(\iota^{-1}(\iota(B_\pi)\Gamma)/(B_\pi)\).  The checker also rejects
the quadratic polar map at an explicit smooth boundary point and proves that
multiplying the cubic by any one hyperplane produces no linear Jacobian
syzygy and hence no Saito-free quartic.  It now also verifies the
etale-component sieve for coupling: a nonzero contraction class requires
the pulled boundary to be supported on the critical divisor \(WKL\).  For
the normalized nonradial triple the pulled boundary specializes at \(W=0\)
to \(-133z_2^3\), so its contraction module is zero.  For the old boundary,
the same checker proves that \(S^2(Q^2+3QR+9R^2)\) is the unique new
degree-four contraction class.  Elimination and the \(L=0\) dominance
test prove the full formula
\((\mathcal B,S^2(Q^2+3QR+9R^2),S^2P^3)\), which closes every base-fixed
use.  The displayed two-mask pair gives exact polynomial inverse quotients,
but the resulting composite Jacobian is \(W^3K^3L\rho_Vz_1/2\).  The
differential argument on \(\mathcal B=S=0\) then closes every incidence
whose new base boundary remains proportional to \(\mathcal B\), including
nonradial feedback.  The final mod-101 identity sieve is exhaustive through
base degree six and finds only \(\mathcal B\) and \(\mathcal B^2\), so a
genuinely different reduced exceptional boundary has degree at least seven.
For the nonreduced hit, the checker also verifies that \(\mathcal B\) has
ordinary node tangent cone \(3x(x+3y)\) and replays the resolved intersection
ledger whose residual ramification class is \(-E_2\) with intersection
\(-1\) against the pulled-back line.  Together with the cited
totally-invariant-curve theorem, the written argument excludes every
degree-at-most-two triple with
\(\mathcal B(p,q,\rho)=v\mathcal B^2\), including triples using the two mask
variables.  A nonreduced realization must therefore start in degree at
least three and cancel its higher terms.
Finally, the same command verifies the dominant node-chord rechart
\((\mu,\lambda,t)\mapsto
(\mu+\lambda f_1(t),\mu+\lambda f_2(t),\lambda f_3(t))\), its boundary
\(27\mu\lambda^2t^3(t-1)^3\), and its Jacobian
\(3\lambda t^2(t-1)^2\).  The irreducible etale witnesses
\(M,H,N_2,WN_2-1\) force its contraction module to vanish, so a successful
coupling must replace at least one node-chord factor by a selector supported
only on \(W,K,L\) after pullback.  The five coordinate-zero witnesses and
the five irreducible shifted-coordinate witnesses reject all sixty
injective placements of \((\mu,\lambda,t)\), so that selector must be
non-coordinate.
The same command then replaces \(\mu\) by the old boundary and verifies the
first nonzero different-boundary coupling.  For \((\lambda,t)=(P,Q)\), its
boundary is \(27\mathcal B P^2Q^3(Q-1)^3\); the class
\(27P^2Q^3(Q-1)^3S^2C_K\) gives two exact polynomial inverse masks and an
explicit polynomial composite, but its Jacobian remains nonconstant.  The
stronger exterior-form calculation proves that every admissible numerator
for a coordinate pair contains
\(h=\lambda^2t^3(t-1)^3\), whereas the log equation contains only
\(9\lambda t(t-1)\).  All twenty ordered coordinate pairs are therefore
closed for arbitrary polynomial masks, not just for the displayed pair.
The root-incidence command verifies a localized two-coordinate
representation of `1/P'(T)`, generic root-field recovery, and the residual
orientation pole obstructing target-only polynomial pullback.  It also
checks the selected rational root and proves that the ordinary `(U,V)`
pullback retains Jacobian `H^6/(2*Theta*K^6*L^3)`.  Neither command verifies
an ordinary polynomial Keller map.  The rank-four command expands the
correct reciprocal `H^3/(4*K^3*L)` in the quartic root basis, verifies its
common denominator `B^2*rho*sigma`, and checks the resulting localized
two-mask determinant-one suspension together with its three genuine
coefficient-ring pole divisors.  It also resolves their common cluster by
four point blowups, computes the branchwise Newton orders of the full
numerator, and solves the resulting local divisor-allocation intervals.
Finally it verifies the forced selector `T+a^3`, its exceptional Cartier
transforms, and the all-degree obstruction to realizing the exact divisor
with masks in the original polynomial root algebra.  It does not construct
the new affine modification needed to adjoin those exceptional quotients.
The local-viability command verifies that the three components have one
common nontransverse cusp/tangency cluster and that the simple second mask
has local order deficit two.  The affine-modification command follows the
forced quotients through the `E3` and `F` charts, verifies their nonnormal
and singular loci, and checks that the smooth negative-definite full-chain
resolution is nonaffine.  The corrected-selector command factors the first
exceptional quartic, derives the higher-order simple-branch correction,
verifies the exact rational two-mask divisor on every resolved branch, and
checks the unimodular relatively ample coarse deletion
`{F,E3,B,Hhat}`.  It then verifies the simple/triple normalized splitting
above retained `E1,E2`, proves that strict `B,Hhat` give only two
horizontal deleted primes, and certifies that the normalized open has
relative class rank at least two.  Finally it computes the normalized
seven-curve intersection matrix, contracts it to the `A3` chain, verifies
the cyclic order-four discriminant group and an explicit odd curvette, and
proves that the forced strict `B` class has content two in the full dual
lattice.  Thus no boundary basis can contain it; the `B,Hhat` pair also
has determinantal divisor two.  The same command then verifies the
`B`-free split `T*Lchi/Hhat | (T+a^3)/Lchi`, its two irreducible curvette
boundaries, the unique unimodular exceptional completion
`{Fs,R2,Ft,R1,Hhat,Dchi,D14}`, and an effective divisor on that support
whose seven exceptional intersections are all positive.  The final
Singular command computes the geometric genus of the irreducible
degree-16 `Hhat` norm as 13, verifies
`Norm(Lchi)=chi*q14`, and computes genus 20 for the degree-14 component.
Since a smooth completion of `A2` has only rational boundary components,
either positive-genus divisor rules out the surviving `B`-free complement
as an affine plane.  The bounded genus-zero search computes the complete
ordinary-total-degree-three valuation spaces in selected-root degrees one
and two.  It verifies four new exact root-linear selectors and proves that
every exact selector in this ansatz has horizontal norm degree at least
sixteen; root-quadratic terms start at norm degree twenty-four.  The bound
is sharp.  The accompanying Singular command verifies irreducible
degree-sixteen norms of genera twelve and fourteen for two sharp
representatives, genus twelve for one displayed `a^3` perturbation, and
genus ten for the near-selector `(b+6)*T-81*rho`.
Those samples do not classify the remaining degree-sixteen parameter plane,
so no rational selector or affine reconstruction is claimed.  The final
Python command constructs its birational degree-ten strict model, proves
generic absolute irreducibility, and verifies the complete fixed-infinity
tangent hierarchy.  Its two terminal members are absolutely irreducible of
genera ten and nine.  The same command factors the moving critical
determinant as a fixed conic squared times one irreducible degree-23 curve,
analyzes the exceptional conic rank drops, and verifies the first two
rational selectors found on that curve.  The parameter `[77:-16:-8]` is
absolutely irreducible of genus twelve with two additional nodes.  The
parameter `[103:-16:8]` has a rational line component, and its coefficient
norm has a rational conic factor, but the remaining component is absolutely
irreducible of genus ten.  Finally, the command checks all 864 primitive
parameters through height six for absolute irreducibility and computes a
144-member height-three genus census.  Those two censuses are bounded
experiments; the full parameter discriminant is not implicitized, and no
affine-space Keller map is asserted.  The conic-principal command purifies
the rational component to the prime ideal \((q_2,\ell)\).  It verifies
\(\operatorname{Norm}(\ell)=-3q_2R_{10}\), with \(R_{10}\) absolutely
irreducible of genus two, and the birational strict factorization
\(3UK^3G_5/H^3\), with \(G_5\) absolutely irreducible of genus two.  Its
conormal computation finds two quadratic closed points, hence four
geometric points, where the conic prime needs two local generators.  This
proves that no coefficient-one principal divisor isolates the reduced
conic.  It then identifies the cluster curvette with class two in the
local \`Z/4\` group, hence local index two, and resolves the other quadratic
pair as an ordinary two-branch conductor node.  The two exact root-chart
preimages have coefficient Jacobians \`-32/27\` and \`-4\`; the conic is
divisorial on the first branch and only a transverse codimension-two
incidence on the second.  Conductor matching and the principal ideal
theorem therefore exclude a support-only principal divisor at every
positive multiplicity.  The full moving-discriminant image, alternative
rational selectors, and a Keller map remain open.  The degree-twenty line
checker constructs the six strict pullbacks for `q0,...,q5`, proves the
complete fixed `K`, `M`, and chart-`rho` divisibility kernels, and gives an
exact four-minor resultant-gcd certificate for every affine line in the full
six-dimensional space.  Its only rational line is `U=0`; the kernel there
is generated by the known `[103:-16:0:8:0:0]` direction and the nonzero-`q5`
direction `T*(a^2-4*rho)`, which is the old rational conic multiplied by
`T`.  The affine-line problem is closed; nonlinear rational components
remain open.  The conic-sieve command exhausts all 3,875 projective
degree-two forms over \(\mathbb F _5\).  It finds only the reductions of
\(K,M,\rho_V\) and one exceptional point; exhaustive lifting gives 5, 25,
and 0 incidence points modulo 25, 125, and 625.  Thus the exception has no
characteristic-zero lift in the good-reduction chart, while the nonreduced
\(K/M/\rho_V\) neighborhoods and degree-dropping reductions remain open.
The cubic-sieve command factors all 3,906 projective selector members over
\(\mathbb F _5\).  The only irreducible cubic factors are the fixed
\(A_*,L,H\) factors and the known cubic component of
\(T(a^2-4\rho)\).  All 38 non-\(H\) incidence points have no moving-factor
tangent.  The artificial \(H\) plane is the sole nontransverse residue; on
the slice \(\widehat R_3+x\widehat R_4+y\widehat R_5\), its eight bounded
lift counts are 5, 5, 25, 25, 125, 125, 625, and 625.  Those counts are an
experiment, not a characteristic-zero existence or nonexistence theorem.
The next two commands are explicitly bounded factorization
experiments: all 175,680
primitive height-at-most-six parameters on `q5=0`, and all 58,095 primitive
height-at-most-three parameters in the full six-dimensional space, are
reducible exactly on the projectivized `K/M` kernels.  The final Singular
checker proves that the one-jet near-selector has an absolutely irreducible
degree-16 genus-10 norm.  Its explicit exact jet correction has an
absolutely irreducible degree-18 strict curve of genus 31.  These results
narrow the selector and coupling searches; they do not construct a Keller
map.

### Davenport alternating coefficient pencils

The first coefficient class beyond AS14--AS17, and the stronger
all-length-two obstruction, are replayed by:

```bash
.venv/bin/python scripts/verify_davenport_independent_marked_line.py
```

The checker forms the general quadratic--quadratic alternating Jung
coordinate with all translations and lower coefficients present, saturates
the four highest fiber coefficients at the outer quadratic leading
coefficient, and obtains the unit ideal over
`Q[a]/(a^2+a+2)`.  It then proves the unbounded length-two statement from
the three exact Newton-leading coefficients.  These are exact symbolic
obstructions, not bounded searches.  No unit-gate solution survives, so no
candidate is promoted to the nonlinear-`U` translated-incidence equations.

## Plane degree-frontier audit

The fixed-coordinate normalized sparse-support exclusions are replayed by:

```bash
.venv/bin/python plane-jc/cas/verify_sparse_support_exclusions.py
```

For `F=(x+P,y+Q)` with `P,Q` having no terms below degree two, the
arbitrary-degree proof classifies every exact support split `1+q` and
`q+1`, `q<=5`.  All charts are unit ideals except the quadratic, cubic, and
quartic
monomial shear chains; their Gröbner bases force compositions of two shears
with explicit polynomial inverses.  The same command checks explicit
Rabinowitsch unit identities for the balanced `2+2` class in arbitrary
degree: 256 determinant/divergence presence masks reduce to 15 linear
survivors and 20 canonical collision partitions, all exactly inconsistent.
The former `14,653,584`-support census through degree twelve remains as an
independent regression.  For `2+3`, the eleven Keller contributions give
2048 presence masks; 321 are compatible with the derivative/determinant
pattern, 98 survive the linear collision sieve, and the global exact
integer no-singleton formula is unsatisfiable.  Thus every `2+3` and `3+2`
chart has a singleton Rabinowitsch unit identity, and every normalized
support of cardinality at most five is invertible.  For support six, the
unbounded exponent formulas classify `1+5` as the quartic shear chain,
make every `2+4` and `4+2` chart a singleton unit ideal, and leave only the
common `{x^2,x*y,y^2}` exponent support in `3+3`; its saturated coefficient
ideal forces a directional quadratic shear with inverse `id-H`.  A separate
`5,290,000`-support census through degree six finds the same unique `3+3`
collision support.  Thus every normalized support of cardinality at most
six is invertible.  Minimizing support over
all tangent-to-identity affine normalizations shows that every noninvertible
plane Keller map has affine-normalized support complexity at least seven.
This is a support-cardinality theorem in fixed
normalized coordinates, not a new universal degree bound.  Support six is
the stopping point for sequential sparse-layer escalation; the next program
must connect affine support to Newton or boundary geometry.  The proof and
exact claim boundary are in
[`plane-jc/CERTIFIED_SPARSE_SUPPORT_EXCLUSIONS.md`](plane-jc/CERTIFIED_SPARSE_SUPPORT_EXCLUSIONS.md).
Intentional artifact regeneration uses `--refresh`.

The exact affine-support/Newton bridge audit is:

```bash
.venv/bin/python plane-jc/cas/verify_affine_support_newton_bridge.py
.venv/bin/python plane-jc/cas/classify_f2_75_125_layers.py
.venv/bin/python plane-jc/cas/reduce_f2_75_125_endpoint_system.py
.venv/bin/python plane-jc/cas/audit_f2_75_125_boundary_handoff.py
.venv/bin/python plane-jc/cas/test_f2_75_125_frontend.py
.venv/bin/python plane-jc/cas/generate_f2_modified_system.py --include-equations --output artifacts/generated-results/jc2_f2_modified_laurent_family.json
.venv/bin/python plane-jc/cas/verify_f2_modified_chart_bridge.py
.venv/bin/python plane-jc/cas/verify_f2_kummer_orbit_transfer.py
.venv/bin/python plane-jc/cas/verify_f2_terminal_residue_cover.py
.venv/bin/python plane-jc/cas/verify_f2_a6_simple_spectator_gluing.py
```

An optional exact characteristic-zero family sample requires Singular and
proves that the `r=5` P-only projected top-gap ideal is also `(1)`:

```bash
.venv/bin/python plane-jc/cas/verify_f2_modified_chart_bridge.py --extended-r5
```

The contact-only artifact was intentionally refreshed after its strategic
recommendation became historical, without changing the four-stratum census:

```bash
.venv/bin/python plane-jc/cas/audit_f2_75_125_boundary_handoff.py --refresh
```

Its current SHA-256 is
`77bffef9fed0ed9749f135a426de945dc27e226d43fce88bf6a45b79bb8a83e5`;
the refreshed checker SHA-256 recorded in `MATH_STATUS.json` is
`952d0955dee25eb96933d36f2510783cfe5610d04ac012bb0836749373ad6684`.
The refresh changes only the verdict context: the old recommended pivot is
now labelled historical and linked to `PF2KO1`, `PF2TR1`, and `PF2GC1`.

The first replay proves that coarse Newton vertices, geometric degree, and
nonproperness data cannot upper-bound affine-normalized support: a
Zariski-open family of triangular automorphisms has fixed coarse geometry
and support lower bound `d-2`.  It then verifies the Kummer chain-rule gate
and the live `(75,125)` F2 terminal characters `P={1,4}`,
`Q={0,1,3}` modulo five, which block constant-Jacobian descent.  The second
replay corrects the Laurent chart to `[t,z]=-z`.  Its upper window contains
all 35 zero layers (`39` through `5`), 665 band pairs, and 978 jet-reduced
linear parameters (973 after normalization).  The corrected post-jet
support-row upper bound is 5,344; exact reconstruction uses the jet equations
and all 165,980 compressed generators.  The same replay now carries the
corner bounds through the complete B0 tail: P bands `-75..15`, Q bands
`-125..25`, 2,418 jet-reduced parameters, 240 zero layers through `-200`,
13,741 band pairs, and 1,327,026 exact compressed generators.  This necessary
over-envelope does not by itself exclude F2.  The replay also proves that the
common-power top root is not an
arbitrary degree-18 polynomial: it has the exact two-parameter form
`H(t)=(1+u+...+u^4)^2*R(u^5)`, `u=1+t`, with `R` quadratic and
`R(1)=1/25`.  It also corrects the former upper-descent claim: the substitution
`p=C0^2*U`, `q=(-9/5)*C0^4*V` was an unproved divisibility restriction.  The
true exact source-band kernels at the first five descents have dimensions
`6,6,7,7,10`, because every P-band direction has the Q follower
`q=-3*C0^2*p`; the extra layer-35 direction is the commuting `C0^4` term.
The actual formal `lambda*C0^(-1)` resonance is at layer 10 and is not an
independent source-band kernel.  The intentionally refreshed layer artifact
also verifies the nonlinear resultant `1701*a^8`, forces root continuation
through descent 7, and leaves the first exact local residual
`27*y^2-9*y+1=0` at descent 8.  Its Q-band-one normalization excludes the
four fixed Kummer double-prime supports.  At this earliest spacing only the
nonzero double-root stratum of `R` remains; an exact normalized `P_3/Q_1` interpolation passes its first
local target jet.  At descent 40 the new 9- and 19-dimensional lower bands
show why `E5=0` is not a valid equation after the target: the correct object
is the layer-zero Fitting row.  The replay now reduces the complete target
operator to twelve local jets and two triangular residues.  It proves that
the local `P_3/Q_1` jet needs an off-grid lowest-`u` correction, derives the
edge equation `A'*D-B*C'=1/5`, and verifies an exact edge witness.  It also
proves that the all-`r` sparse witness completes only as an infinite formal
shear.  It classifies the polynomial order-two repair, proves that it never
terminates quadratically, and gives exact `r=3` unit eliminations for cubic
and quartic termination.  It now reconstructs the exact `v^5` Kummer packet:
the omitted fifth-binomial correction cancels the top conflict, its `8 x 10`
matrix has unit minor `3^5*5^16*e^13` over the rank-two branch, and a second
unit Bezout minor carries the edge recursion through `v^10`.  It also proves
the two all-`r` unit-minor formulas and traces their `18*r-1` pivots to exact
polynomial source combinations.  Their source-degree gaps are `12*r-8` on P
and `24*r-20` on Q, with minimum terminal-edge slack one; at `r=3` there are
`53` lifts with largest degrees `47<75` and `73<125`.  The two exact
confluent-CRT determinants at `0,1,w0` then
quotient the triangular `w=0` conditions and leave a rank-`24` global
Hermite coordinate module over the rank-two candidate algebra.  Its exact
fixed-endpoint cone has 22 target and 25 nonzero-weight layer-zero pairs.
The leading target coordinate is normalized, while the remaining ten
`w=1` rows have affine-linear determinant `75000`; degree-seven followers
divisible by `w^2` preserve the controlled `w=0` jets.  Eliminating those
ten variables leaves thirteen global Hermite coordinates.  It also
integrates layer zero into the length-15 algebra
`B[w]/(w^3*(w-1)^6*(w-w0)^6)`, whose quotient by constants is the exact
rank-14 Fitting residue.  Explicit target minors and a constructive
layer-zero span show why neither row alone contradicts the unrestricted old
B0 bands; their earlier triangular equations must be substituted first.
The new endpoint-reduction replay performs that substitution exactly.  Its
ten pivot circuits contain `1,489` terms and give a degree-eight upper bound
after forming brackets.  It derives `1,172` normalized zero-row coefficient
slots and retains all `1,061` active source coordinates, including the
P-minus-21/Q-minus-11 block which re-enters layer `3`.  Explicit tridiagonal
unit minors then eliminate `134` endpoint-disjoint new-Q coordinates on
layers `39..29`, leaving `219` upper Fitting slots and `927` active source
coordinates.  The surviving task is the coupled Schur/Fitting calculation
from descent `12` (layer `28`) through descent `37`, followed by the thirteen
`w=w0`/residue functionals, not another endpoint rank count.
The endpoint-reduction artifact has SHA-256
`9834ed2ba4e64a2b034a83cd0604140206f1a8192a561509c208d02e0a0ca189`;
its checker SHA-256 is
`762ff5e509abcf7701beea4a836b99853d33777f38cdf088eda356be23695858`.
<!-- status-consumer: PF2ER1 64378dad616fc3f2 -->
The artifact separately enumerates the still
open later first-defect spacings `9..90`.  The artifact
has SHA-256
`96e4fd2ff853fcba9d41a72973ea55a1acb2cd41eefad21e69efd6ae73df8b8b`.
The boundary-handoff replay then gives four exhaustive
contact partitions of 25.  It proves that their multiplicities do not
determine branch scales or finite-normalization rows; even the unsupported
strongest contact-to-row surrogate survives the degree-26 packet budget.
This stops the nonlinear triangular elimination, while the finite
character-resolved B0 system now includes every lower band.  Subsequent Kummer-orbit and
terminal-residue calculations, described below, bypass the contact surrogate
and reopen F2 at the global gluing stage; see
[`plane-jc/AFFINE_SUPPORT_NEWTON_BRIDGE.md`](plane-jc/AFFINE_SUPPORT_NEWTON_BRIDGE.md).
The frontend replay checks the forced chain and terminal normalization
independently.  The modified-system generator reconstructs the published
`r=2` systems and the 14- and 22-coefficient-function `r=3` windows under its
explicitly stated common-power ansatz.  It triangularly eliminates the power
rows, gives the remaining core in `B=A^r` coordinates, and presents the true
residue as an Artinian Fitting ideal: a `16 x 16` determinant for `d=2,h=2`
and the four-generator maximal-minor ideal for `d=3,h=2`.  It also checks the
certified `h=4` alternative, which has zero and two generic residue equations.
Every compact residue has an endpoint-binomial point.  More sharply, the
`d=2` residue is rational on a dense open, while the `d=3` residue has a
smooth fourfold meeting the coefficient torus and an exact cubic-invariant
subfamily.  The uniform congruence-support gate excludes all `d=2,3`
congruence sections under the certified `X^4` weight.  For the surviving
`d=3,h=2` section, the Laurent antiderivative constant is impossible and the
only monomial-bracket survivors lie on an exponent-seven, `lambda=0` ray;
its first member is stored as an exact eight-row formal residual point.
The calculation does not prove the ansatz, `d=2,3`, or the lower Laurent-`y`
support ledger, so it does not exclude the unrestricted family.
See
[`plane-jc/F2_MODIFIED_LAURENT_FAMILY.md`](plane-jc/F2_MODIFIED_LAURENT_FAMILY.md).
The generated JSON has SHA-256
`bca206498c153e41a2f31344015df2ce63890f8b12228b6d0ba2c0970eb87c85`;
its software assumptions are `.python-version`, `requirements.txt`, and exact
characteristic-zero SymPy arithmetic.
The modified-chart checker then derives `gamma=2`, the `d=3` monomial chart,
and every possible nonnegative-`xi` support position from the F2 corner
chain.  It preserves the translated binomial-jet relations instead of
treating the support box as independent: at `r=3` the exact P/Q image ranks
are `74/83` and `196/215`.  A primitive top-band relation excludes the formal
terminal point.  For the full projected top diagonal, the three P gap
equations define a length-27 Artinian algebra and the first Q gap has nonzero
resultant/multiplication determinant in it, so the combined ideal is `(1)`.
Thus every branch of the literal bracket-preserving polynomial projection is
excluded.  This also proves that naive deletion of the negative Laurent tail
cannot be the missing theorem.  The full seed has the explicit source lift
`x*(x*y^5-1)^2*R(x*y^5)`, and the corrected top tangent shows that its
first-five kernel dimensions are `6,6,7,7,10`; `lambda*C0^(-1)` is a
layer-10 formal resonance, not an independent source mode.  Deriving its
nonlinear `F`-tail cancellation remains open.  The
pinned chart artifact has SHA-256
`ac7dbc170cafbcf028079b9ccdb41afd78c333e3850abc0761f64cc056e7d7b8`.
<!-- status-consumer: PF2MCB1 6ff13314e0090f52 -->
The terminal checker first verifies the all-parameter degree-`2r` passport
`(2r-1,1)|(r,r)|(3,1^(2r-3))`, geometric monodromy `A_(2r)`, and discriminant
squareclass `(-1)^(r+1)*(2r-1)`.  It then specializes with the Kummer replay
to the terminal target row with transverse index one, residue degree six,
passport `(5,1)|(3,3)|(3,1,1,1)`, and geometric monodromy `A_6`.  Global
boundary gluing remains open.  The checker also proves that the natural
`A_6` action is four-transitive and primitive (so the residue cover has no
`2`-by-`3` factorization), that its target-fixed deck group is trivial, that
`e=1` gives zero transverse different, and that the residue formula is
parameter-free.  It emits residue-different packet `(4,2,2,2)` and verifies

```text
disc_s(125*s*(s+1)^5-r*(9*s^2+15*s+5)^3)
  = 5^17*r^4*(729*r-125)^2,
```

so the rational model has arithmetic `S_6` over `Q(r)`, geometric `A_6`, and
quadratic constant field `Q(sqrt(5))`.  After rescaling the third branch value
to one it is a Belyi map; its regular geometric `A_6` closure has signature
`(5,3,3)` and genus `25`.  The two target toric nodes have exactly three
preimages in the source-divisor interior, fixing three boundary-attachment
points with different contributions `(4,2,2)`.  The last contribution `2`
is at the source toric endpoint over the smooth third branch value.
The target valuation equality gives
geometric degree at least six, or at least twelve for two distinct packets
over the same target divisor.  Packets over distinct target divisors do not
add.  Since the certified target valuation is centered at infinity, the
affine-companion theorem supplies no `+1`; purity instead requires a separate
affine ramification row.  Global geometric monodromy has `A_6` as a
nonabelian simple composition factor.  Thus the remaining ledger has one
squarefree packet or two identical double-root packets attached to the same
versus distinct target components.  The final command exhausts the stated
two-transposition simple-spectator model: all six degree-seven genus-zero
classes survive and generate `S_7`, while paired-star witnesses show that the
coarse filters allow every larger remaining degree.  Adding the certified
endpoint/interior markings leaves three classes of signature `(5,3,1)` under
the strongest naive requirement that the connector anchor avoid both source
endpoints.  The same replay compares the two order-five structures exactly:
the terminal inertia normalizer is `AGL(1,5)`, its multiplier parity and the
arithmetic `S_6` sign both cut out `Q(sqrt(5))`, but `A_6` is perfect and the
full residue deck group is trivial.  Under the stronger assumption that one
Kummer orbit contributes five disjoint transpositions at one branch value,
Riemann--Hurwitz and a rational connected source boundary force degree eleven
and monodromy `S_11`; matching the five core anchors to the inertia support leaves one unoriented class, or four
Kummer-generator orientations.  Both gluing models are conditional because
sheet specialization is not a toroidal node-gluing theorem and the Kummer
spectators have not been assigned certified branch cycles.  These
calculations do not exclude `(75,125)`.

The expanded terminal checker is pinned in `MATH_STATUS.json` at SHA-256
`baa8fe7abdcf1652bc0a8636437b9505ebc6838bf2eebaf2492c03684fd63cbf`.
The software assumptions remain `.python-version` and `requirements.txt`;
the command above both recomputes the assertions and emits the final pass
marker.
<!-- status-consumer: PF2GC1 33dbc5ff48b5d064 -->

The small deterministic regression of the published candidate tables is:

```bash
python3 plane-jc/cas/frontier_125_150.py
```

The direct unibranch finite-flat attack is replayed by:

```bash
.venv/bin/python plane-jc/cas/verify_unibranch_spectator_models.py
```

It verifies a universal rank-`n+1` family with a clean singular unibranch
boundary collision of length `n` and a separate étale spectator.  The
quartic member realizes both exact `3+1` and `2+2` frontier fibers.  These
are countermodels to a purely local exclusion, not Keller maps: deleting
the principal ramification curve gives `A1 x G_m` with a nonconstant unit,
not the distinguished `A2` open.  See
[`plane-jc/UNIBRANCH_SPECTATOR_COUNTERMODELS.md`](plane-jc/UNIBRANCH_SPECTATOR_COUNTERMODELS.md).

The global quartic Cox-lattice continuation is included in:

```bash
.venv/bin/python plane-jc/cas/test_plane_boundary_exclusion.py
.venv/bin/python plane-jc/cas/test_degree_zero_endpoint_pairing.py
Singular -q plane-jc/cas/quartic_completed_deletion.sing
```

It separates the one-boundary row and the same-target/different-target
versions of the two-boundary row.  Their target-pullback lattices have
respectively index-two, rank-one, and index-two defects.  In every case the
single primitive ramified-boundary character saturates the exponent
lattice; it is also the canonical/different class.  This is a reduction,
not a quartic exclusion.  In the rank-one row, `g=a*s_E^2` defines the
finite target-side normalization input.  Base change `B/A` to that
hypersurface to obtain a rank-four finite-free order and normalize it in
`k(x,y)(s_E)`.  The distinguished source `A2` supplies the complementary
normal Gorenstein hypersurface `g(P,Q)=a*s_E^2`; its coordinate ring
contains the finite normalization and gives the graded Zariski--Main open
immersion restricting over `s_E!=0` to
`A2 x G_m -> X x G_m`.  The completed calculation is now part of the
replay: at the `3+1` cusp and each branch of a `2+2` connector the order is
`a*s^2=r^2*ell`, its normalization adjoins `z=r*ell/s`, and its conductor
and canonical module are `(r,s)`.  The source deletion `D(r)` is locally
compatible in every chart.  Its transitions are compatible too:
`r_i=u_ij*r_j` gives `ell_i=u_ij^-2*ell_j` and
`z_i=u_ij^-1*z_j`.  The revised target is therefore the two-generated
degree-zero global-section algebra and its cusp/connector endpoint
pairing, not nonprincipality alone.  The replay also checks the sharper
module form: the degree `-1` square map into degree `-2` has affine
companion cokernel `k[x,y]/(h)`.  A descended unit or exceptional curve is
now a possible witness, not the statement being assumed.  The second
Python command proves the proposed four-filter semigroup search is not yet
finite: for every `n>=1` the reduced divisor
`(y^2-x^3)*(x*y-y^(n+1)-1)` has cusp pole pair `(2,3)` and connector pole
pairs `(1,0),(n,1)`, while the connector matrix is unimodular and all
degree-zero, rank-one-piece, and odd-square-cokernel conditions survive.
The bounded loop through `n=20` is only a regression for the uniform
algebraic proof.  It also checks the triangular coordinates
`(x,y+x^m)`: the connector pole pairs become `(1,m),(n,m*n)` without
changing the affine plane, graded bridge, or packet data.  Hence no bound
on raw coordinate-generator poles can be intrinsic.  A genuine finite
compiler first needs the marked multivaluation semigroup or an
automorphism-minimal coordinate pair, a degree-four bound on its minimized
pole height, and the conductor equivalence relation pairing connector
endpoints.  On the displayed connector the inverse shear
`u=x-y^n=t^-1, v=y=t` gives the exact minimum height two.
Finally, the two-ended family
`X=t+t^-1, Y_c=t^2+c*t^-2+t` has the same pole row `(1,2)` at both
endpoints for every nonzero `c`, but one quadratic shear cancels both
leading poles only when their residue ratios agree (`c=1`).  The compiler
must therefore retain simultaneous initial-residue data, not only numerical
semigroup values.  The imported
`plane-jc/cas/endpoint_valuation_compiler.py` enumerates every monomial
triangular shear that strictly lowers the two-endpoint pole height, in both
orientations, and terminates by integer descent.  It reduces every displayed
connector to height two and distinguishes the residue-matched and
residue-mismatched rows.  It now also exhausts complete lowering polynomial
shears by recursively retaining forced cancellation terms in strictly
descending degree, even when a leading prefix is height-neutral or
height-increasing.  The witness
`u=t^-1+t^2`,
`v=u^3-t^6+t^5+u^2`
has no lowering monomial shear, but the compiler finds
`P(u)=u^3+u^2`, which lowers total height from eleven to nine.  This closes
the one-polynomial-shear gap.  Proposition 6.6 now closes the reduced
alternating-direction peak gate as well.  At each marked valuation, the pole
change made by the second of two opposite degree-at-least-two factors is at
least the change made by the first.  The factor-height increments along a
reduced Jung word are therefore nondecreasing, so a globally lowering word
starts with a lowering complete factor.  As a bounded regression, the
checker exhausts all 49
ordered nonempty-support pairs on Laurent exponents `{-1,0,1}` and every
alternating two-step monomial shear of degrees one or two with coefficients
`+/-1`.  It finds 16 paths which lower height after a nondecreasing first
step; every initial pair already has a lowering complete polynomial shear,
so no terminal peak counterexample occurs in this grid.
The expanded marked multi-pole experiment is:

```bash
.venv/bin/python plane-jc/cas/search_marked_multipole_peak.py \
  --max-degree 3 --max-length 2 --include-linear \
  --extended-seeds --max-seed-terms 2 --scan-all \
  --output artifacts/generated-results/marked_multipole_peak_search.json

.venv/bin/python plane-jc/cas/search_marked_multipole_peak.py \
  --max-degree 3 --max-length 3 \
  --extended-seeds --max-seed-terms 2 \
  --output artifacts/generated-results/marked_multipole_peak_reduced_degree27.json
```

The first command tracks exact valuations, initial coefficients, conductor
pairing, and every factor height in two- and three-pole rational charts.  It
finds 166 delayed lowering paths, including three-pole peaks of shape
`4 -> 5 -> 3`; all 166 initial states already admit a complete lowering
shear.  The second command tests reduced alternating degree-two/three words
through length three and polydegree 27 on the extended seed basis.  It checks
93,440 words from complete-shear-terminal states and finds no counterexample.
These exact bounded experiments are regressions for the uniform pole-change
proof, not its logical basis.  The proof, generated-orbit continuation, and
the complementary signed and complete-factor commands are recorded in
[`plane-jc/MARKED_MULTIPOLE_PEAK_EXPERIMENT.md`](plane-jc/MARKED_MULTIPOLE_PEAK_EXPERIMENT.md).
See
[`plane-jc/JC2_GLOBAL_COX_PACKET_ATTACK.md`](plane-jc/JC2_GLOBAL_COX_PACKET_ATTACK.md).
The normalization and conductor formulas are written algebraic proofs.  The
Python checker replays their cusp factorization, determinantal identities,
and monomial conductor quotient; the independent Singular command computes
the cusp normalization and conductor and verifies normality of the
determinantal overring.

The conductor-decorated endpoint-semigroup continuation is:

```bash
.venv/bin/python plane-jc/cas/experiment_quartic_endpoint_semigroups.py \
  --max-connectors 4 \
  --max-pole 8 \
  --max-contact 8 \
  --cutoff 12 \
  --output artifacts/generated-results/quartic_endpoint_semigroups.json
```

For each bounded connector row it computes the cusp semigroup
`<2,3>`, the displayed connector polar-bound monoid, the exact signed
two-endpoint valuation semigroup
`{(u,v) in Z^2 : u+v<=0}`, the residue completion, the conductor endpoint
pairing, and the odd-square contact vector.  It also replays the 24 cusp
braid pairs and three connector sheet matchings and checks the rank-one
graded bridge.  The output is a bounded feasibility report, not a quartic
exclusion.  The uniform carrier
`(y^2-x^3)*product_i(x*y-y^(n_i+1)-lambda_i)` proves that the listed inputs,
even with conductor pairing restored, do not bound connector count,
displayed pole parameters, or completed contact.  See
[`plane-jc/QUARTIC_ENDPOINT_SEMIGROUP_EXPERIMENT.md`](plane-jc/QUARTIC_ENDPOINT_SEMIGROUP_EXPERIMENT.md).

The exact quartic linear-pencil calibration is:

```bash
.venv/bin/python plane-jc/cas/experiment_quartic_keller_pencil.py \
  --output artifacts/generated-results/quartic_keller_pencil_calibration.json
```

It compares the finite-free packet model
`(y,x^4-x^3+x*y)` with its target shear
`(y,x^4-x^3+x*y+y^2)`. The maps have the same `3+1` cusp, the same
`2+2` connector, and the same Jacobian determinant. The script resolves
the symbolic linear pencil on the strata `beta!=0` and
`[alpha:beta]=[1:0]`. For the unsheared map the generic fiber is `C*` and
the zeta function is `1`; after the target shear it has genus one, two
punctures, and zeta `1/(1+t^2)`. Both maps are non-Keller packet
countermodels. The output proves that the packet does not determine the
pencil and that the linear pencil is not invariant under nonlinear target
equivalence; it does not exclude either quartic Keller packet. See
[`plane-jc/KELLER_PENCIL_AT_INFINITY_EXPERIMENT.md`](plane-jc/KELLER_PENCIL_AT_INFINITY_EXPERIMENT.md).

The exact chart-aware boundary localization/Smith-normal-form prefilter is:

```bash
.venv/bin/python plane-jc/cas/boundary_lattice_prefilter.py
```

The new pre-coefficient front ends are:

```bash
.venv/bin/python plane-jc/cas/test_intrinsic_a2_boundary.py
.venv/bin/python plane-jc/cas/test_log_boundary_compiler.py
.venv/bin/python plane-jc/cas/test_poisson_square_rigidity.py
```

The first reconstructs the canonical class of a complete `A2` boundary,
checks the adjunction/Noether identities, and audits target pole vectors,
ramification, and intrinsic dicriticals.  It proves that a nonproper Keller
resolution needs canonical free depth at least three.  The second compiles
certified branch scales to toroidal proximity and complete
boundary data.  It extracts the local `(2,1),(3,1),(4,1)` rays from the
published `(72,108)` case tree, distinguishes them from the longer adapted
map-base ideals `(t,x^4),(t,x^6),(t,x^8)`, compiles the isolated source chains
of lengths `4,6,8`, and verifies that the common order-four step collapses all
three cases to the same eight-blowup translation graph.  The `F_4` transition
and affine-plane fill then give a unimodular `10 x 10` source boundary
passing adjunction.  The factor-residue tree is encoded symbolically.  The
unselected order-three factor avoids both the common order-four center and
the filled divisor.  The complete common-graph pole vector has no dicritical
component.  A smooth point of `E3` is the unique one-blowup zero-pole
extension; exact two-step witnesses over `Yinf`, `E4`, `E7`, and `E8` delimit
that numerical minimality claim.  The first weighted Wronskian instead
forces the actual `E3 intersect E4` cluster with ten simple children.  At the
plane-return corner, the Poisson-square edge produces a quartic common
factor; all five root-partition fans compile with complete matrices,
differents, conductors, and ramification vectors.  The primary split-factor
formula and alternate-factor chart then select the quadruple-root package
and control its transverse terms.  Both terminal cases have the same
23-component boundary with one degree-twelve dicritical, so the
chain-to-boundary gap is closed.
The third classifies the entire
geometric reduced three-layer Poisson-square locus into the tangent closure
and the `C=0`, `A=0` components, with generic multiplicities `2,3,1`.  Its
exact principal-chart audit proves `I:d0^infinity=I`, excluding associated
primes on `d0=0`.  The normalized `d3,d2` colon filtration classifies the
complete associated-prime set: three minimal primes, three embedded
intersection surfaces, and two embedded core/intersection curves.
The four fast plane checks above run under:

```bash
make verify-plane-jc
```

The independent Singular scheme checks are:

```bash
make verify-plane-poisson-radical
make verify-plane-poisson-primary-charts
make verify-plane-poisson-separators
make verify-plane-poisson-primary-filtration
make verify-plane-poisson-filtered-modules
```

The exact 90 MB certificate archive and extracted replay source are pinned
under `plane-jc/external/zenodo-21479814/`.  Attachment hashes, environment
versions, the full replay command, and the independent hard-certificate
command are in
[`plane-jc/cas/README.md`](plane-jc/cas/README.md).  The independent checker
does not import the primary CAS or generation modules.

For the same archived exact replay with portable process-level CPU
parallelism and without changing the pinned external snapshot, run:

```bash
.venv/bin/python plane-jc/cas/verify_72_108_exact_fast.py --jobs 4
make verify-plane-72-108-exact-fast PYTHON=/absolute/path/to/venv/bin/python
```

GPU backends are intentionally deferred.  They may later be used for modular
or bounded discovery workloads only when their output has a separate portable
CPU verifier.

## Shared `JC_2`--`HC_4` isotropic boundary bridge

The combined programme is
[`JC2_HC4_SHARED_BOUNDARY_PROGRAM.md`](JC2_HC4_SHARED_BOUNDARY_PROGRAM.md).
Its exact first calculation identifies the cotangent determinant
\(\det\operatorname{Hess}(tP+mQ+H)=J(P,Q)^2\), the first isotropic Schur
remainder \(-\Phi_{mm}R(P)\), and the quartic packet's reduced conormal
residue \(2\ell\).  The completed continuation computes the \(3+1\) cusp
and both \(2+2\) connector initials, proves that the relevant positive
associated-graded conductor maps are isomorphisms, and finds
\(\operatorname{Obs}_{\rm pair}=0\) for all 72 monodromy-compatible
labellings.  Run:

```bash
.venv/bin/python scripts/verify_jc2_hc4_isotropic_boundary_bridge.py
.venv/bin/python scripts/verify_jc2_hc4_global_jet_transport.py
```

This does not verify isotropic-flag recognition for an arbitrary
four-variable potential.  It proves instead that the proposed local
paired initial-conormal cokernel cannot be nonzero without an additional
global transport between the two connector jet lines.  The second checker
shows that the marked affine-line normalization makes their projective jet
ratio intrinsic, with exact quartic value \([-1:1]\), but that the
dualizing-residue quotient is annihilated by the conductor.  The abstract
cusp-node family \(R_\lambda\subset k[T]\) realizes the varying ratio
\([-1:\lambda^2]\).  Its symmetric member is the exact two-generated plane
carrier \(k[T^2,T^3(T^2-1)]\), with equation
\(y^2=x^3(x-1)^2\), the anti-diagonal ratio \([-1:1]\), and the required
cusp/node conductor.  Its actual vertical braids are \(\sigma^3\) at the
cusp and \(\sigma^2\) at the node, so van Kampen imposes both the braid and
commutation relations.  These force equal meridian transpositions and
exclude a connected degree-four cover over this structured carrier.  Other
carriers still require their own global braid factorization.

## Direct Schur-descent audit for `HC_4`

The reusable six-to-five-variable construction is imported in
[`MENG_YANG_SCHUR_DESCENT_BRIDGE.md`](MENG_YANG_SCHUR_DESCENT_BRIDGE.md).
For \(\Phi=tA+B\), it separates the two exact hypotheses--constant bordered
Hessian determinant and identically singular reduced pencil
\(\operatorname{Hess}_w(B+sA)\)--from the automatic conclusion
\(\det\operatorname{Hess}(B+\lambda A^2/2+\mu A)=-\lambda c\).
The doubled-Keller row count proves the singular-pencil hypothesis directly,
so no homogeneous-face calculation is needed to reconstruct that descent.
The commands below concern the genuinely additional attempt to descend from
five variables to four.

Two algebraic strengthenings and the resulting research gates are in
[`SCHUR_DESCENT_CONTINUATIONS.md`](SCHUR_DESCENT_CONTINUATIONS.md).  The
exact scalar formula is
\[
\det\operatorname{Hess}\psi_{\lambda,\mu}
=\det M(\mu+\lambda A,w)-\lambda c,
\]
so a fixed descent needs only the specialized reduced determinant to be
constant.  The simultaneous \(r\)-pivot theorem uses the corank condition
\(\operatorname{rank}M\le m-r\); at \(n=3,r=2\) it gives a direct
six-to-four-variable bridge from any Keller collision admitting a jointly
affine two-dimensional source block.  The existing linear-direction audit
rules out such a block for the foundational gauge family, and the general
codimension-one block-affine theorem `SDX2` proves that every three-variable
Keller map with such a block is a polynomial automorphism, even after a
nonlinear source rechart.  Only mixed source--dual or coisotropic pivots
survive this route.

The parameter-uniform affine and low-degree graph audit of the v2
Meng--Yang family is:

```bash
.venv/bin/python scripts/verify_hc4_meng_yang_graph_obstructions.py
```

It verifies the three affine-normal chart coefficients
\(793152L^4,2160L^4,2160L^4\), the zero-dual-normal rank obstruction, and
the degree-at-most-three graph chain

\[
 \rho=-89/16,\qquad \sigma=-\tau,\qquad
 [t^5]D_R(0,t,t,0)=197LN^3/4.
\]

For degree four it retains every graph coefficient and computes the full
two-slope pencil \((0,t,ct,dt)\) only through degree eight.  The leading
square kills all ten \(x_1\)-free quartic jets containing \(y_2\), including
\([x_2^3y_2]R\).  The next coefficients reduce the remaining slice to

\[
 160\rho^2+1968\rho+6021=0,
 \qquad \operatorname{disc}=576\cdot34.
\]

At degree four the same truncated determinant then forces
\([y_1^3]R_3=[y_2^2]R_2=0\) and
\(8\rho^2+99\rho+279=0\).  The two quadratics have resultant
\(16959456\), so no branch survives over any characteristic-zero field.
These results are `HC4MYA1`, `HC4MYG3`, and `HC4MYG4`; degree five is the
first not fully classified single-graph degree.  See
[`HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md`](HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md).

The reverse scalar-pivot classifier and simultaneous matrix-pivot equation
builder are:

```bash
.venv/bin/python scripts/verify_hc4_reverse_schur_descent.py
```

The command generates
`artifacts/generated-results/hc4_reverse_schur_descent.json`. It verifies
the exact identities supporting `HC4RSD1`. An identically singular scalar
pencil whose generic kernel line is constant in the four reduced variables
has generic corank one, and the
bordered-unit equation forces a common constant kernel direction of `A`
and `B`.  Every reduced collision fiber is then a three-variable
constant-Hessian gradient fiber, hence a singleton by `HC3`.  This closes
the homogeneous scalar cone-pencil stratum and intersects it trivially with
all 318/306 live affine-degree-two/three projective-polar rows.  The same
checker constructs the exact corank-minor, integrability, collision, and
Schur equations for matrix pivots.  Nonhomogeneous scalar pencils with an
`x`-moving kernel line, nonsingular scalar pencils with determinant-term
cancellation, and moving matrix-pivot kernel planes remain open.  See
[`HC4_REVERSE_SCHUR_DESCENT.md`](HC4_REVERSE_SCHUR_DESCENT.md).

Continue from constant kernel directions to the first affine moving-kernel
stratum with:

```bash
.venv/bin/python scripts/verify_hc4_affine_moving_kernel_pencils.py
```

Together with the classification argument in the accompanying note, these
exact calculations establish `HC4RSD2`.  The adjugate Piola identity
classifies a unimodular affine kernel vector as either constant or, up to
affine coordinates,
`(z,1,0,0)`.  The checker integrates the latter kernel without a degree
bound, obtains
`f=y*C(z)+(x-y*z)*C'(z)+G(z,w)`, and imposes the complete pencil and
bordered-unit equations.  They force the normal form recorded in
[`HC4_AFFINE_MOVING_KERNEL_PENCILS.md`](HC4_AFFINE_MOVING_KERNEL_PENCILS.md),
whose every scalar Schur descendant has an explicit triangular polynomial
inverse.  The generated artifact is
`artifacts/generated-results/hc4_affine_moving_kernel_pencils.json`.
The affine branch left by `HC4RSD2` is parameter-moving, with primitive
kernel form `v0(x)+s*v1(x)`.

Close that parameter-moving affine branch with:

```bash
.venv/bin/python scripts/verify_hc4_parameter_moving_affine_kernel_pencils.py
```

Together with the proof in
[`HC4_PARAMETER_MOVING_AFFINE_KERNEL_PENCILS.md`](HC4_PARAMETER_MOVING_AFFINE_KERNEL_PENCILS.md),
this establishes `HC4RSD3`. The bordered adjugate bounds a primitive
kernel's parameter degree by one. Piola reduces its affine linear part to
rank-one compression pencils; proportional and common-image pencils collapse
to a fixed line, the constant-at-infinity corner is inconsistent, and the
sole moving common-covector integral has Hessian rank at most two. The
generated ledger is
`artifacts/generated-results/hc4_parameter_moving_affine_kernel_pencils.json`.

The first degree-unbounded nonlinear continuation is:

```bash
.venv/bin/python scripts/verify_hc4_univariate_shear_kernel_pencils.py
```

This supports `HC4RSD4`. For an arbitrary nonconstant `P(z,w)`, the fixed
kernel `v=(P(z,w),1,0,0)` integrates completely to
`f=x*a(z,w)+y*b(z,w)+G(z,w)` with `db=-P*da`. The bordered unit kills the
curvature of the common transverse composite, forcing it onto one linear
form. The resulting univariate normal form has an explicit triangular
polynomial inverse for every Schur descendant. See
[`HC4_UNIVARIATE_SHEAR_KERNEL_PENCILS.md`](HC4_UNIVARIATE_SHEAR_KERNEL_PENCILS.md).
The generated ledger is
`artifacts/generated-results/hc4_univariate_shear_kernel_pencils.json`.

Continue from one transverse-polynomial shear component to a unimodular
pair with:

```bash
.venv/bin/python scripts/verify_hc4_two_component_quasitranslation_kernels.py
```

Together with
[`HC4_TWO_COMPONENT_QUASITRANSLATION_KERNELS.md`](HC4_TWO_COMPONENT_QUASITRANSLATION_KERNELS.md),
this supports `HC4RSD5`. For a fixed primitive kernel
`v=(P,Q,0,0)`, Piola and unimodularity first force `P,Q` to be independent
of the two active variables. Exact mixed-partial integration then gives
`f=x*a(z,w)+y*b(z,w)+G(z,w)` and `P*da+Q*db=0`. The bordered unit makes
`P,Q` algebraically dependent and supplies a polynomial frame of constant
determinant over their closed common composite. Its two differential
coefficients cannot both be nonzero by polynomial degree, so `(P,Q)` lies
on an affine line. A constant active-coordinate change reduces the kernel
to the `HC4RSD4` shear form, whose descendants have triangular polynomial
inverses. The generated ledger is
`artifacts/generated-results/hc4_two_component_quasitranslation_kernels.json`.

Test whether a direct HC4 candidate is covered by an affine singular
scalar pivot with:

```bash
.venv/bin/python scripts/verify_hc4_affine_pivot_coverage_gate.py
```

Together with
[`HC4_AFFINE_PIVOT_COVERAGE_GATE.md`](HC4_AFFINE_PIVOT_COVERAGE_GATE.md),
this supports `HC4RSD6`. For a constant-Hessian potential with Hessian
`H`, a nonzero constant covector `ell` gives an affine singular lift if and
only if `ell^T*adj(H)*ell` is a nonzero constant. On an essential-rank-three
quintic top, the pivot must annihilate the constant top kernel. The next
metric face and the existing Schur face then force a constant relation
`a^T*adj(C)*grad(s3)=0`. Equivalently, all 3-by-3 minors of the degree-eight
coefficient matrix of the Schur vector vanish. The diagonal exact
calibration satisfies the Schur equation but has coefficient rank three,
showing that the coverage gate is a genuine additional restriction. More
generally on that diagonal top, Schur divisibility leaves coefficients
`alpha,beta,gamma`, and affine coverage is confined exactly to
`alpha*beta*gamma=0`. The generated ledger is
`artifacts/generated-results/hc4_affine_pivot_coverage_gate.json`.

Impose the marked collision on every affine-pivot coverage component with:

```bash
.venv/bin/python scripts/verify_hc4_affine_pivot_collision_fibers.py
```

Together with
[`HC4_AFFINE_PIVOT_COLLISION_FIBERS.md`](HC4_AFFINE_PIVOT_COLLISION_FIBERS.md),
this supports `HC4RSD7`. In coordinates adapted to a constant pivot
covector `ell`, the metric numerator `ell^T*adj(Hess(psi))*ell` is, up to a
nonzero constant square, the Hessian determinant of the restriction to
`ell.x=c`. If that numerator is a nonzero constant, every pivot fiber is a
three-variable constant-Hessian potential. `HC3` makes its tangential
gradient injective, so equal full gradients and equal pivot values force
the two points to coincide. Thus an affine zero-corner scalar parent may
represent a four-variable potential but cannot inherit its marked collision
at a common parent pivot value, whether its reduced Hessian is singular or
lies in the nonsingular exact-remainder branch. The generated ledger is
`artifacts/generated-results/hc4_affine_pivot_collision_fibers.json`.

Begin the nonlinear scalar-pivot branch with:

```bash
.venv/bin/python scripts/verify_hc4_quadratic_pivot_rank_obstruction.py
```

Together with
[`HC4_QUADRATIC_PIVOT_RANK_OBSTRUCTION.md`](HC4_QUADRATIC_PIVOT_RANK_OBSTRUCTION.md),
this supports `HC4RSD8`. If `A` is quadratic and
`det Hess(B+s*A)=0`, then a constant nonzero five-variable parent Hessian
determinant forces `rank Hess(A)<=2`. Rank four is excluded by the leading
pencil coefficient. In rank three, splitting off the null direction first
kills `D_z^2 B`; the cleared bordered identity then says that a polynomial
square equals a nonzero constant times `det(s*Q3+H)`, which has degree three
in `s`. The bordered unit makes the affine entries of `grad(A)` generate the
unit ideal, forcing a nonzero linear slice on `ker Hess(A)` and hence the
normal form `A=w+u^T*Qr*u/2`. The checker also gives an exact rank-two
fixed-kernel calibration, so the bound is sharp. The remaining quadratic
frontier is the rank-one and rank-two moving-kernel locus. In rank two, the
passive binary Hessian of `B` is singular. Its rank-zero stratum has
`det(M)=det(D)^2` and reduces to the fixed-support two-component kernel
theorem `HC4RSD5`; only passive rank one is genuinely new. The generated
ledger is
`artifacts/generated-results/hc4_quadratic_pivot_rank_obstruction.json`.

Close the rank-two quadratic branch with:

```bash
.venv/bin/python scripts/verify_hc4_quadratic_rank_two_pivots.py
```

Together with
[`HC4_QUADRATIC_RANK_TWO_PIVOTS.md`](HC4_QUADRATIC_RANK_TWO_PIVOTS.md),
this supports `HC4RSD9`. In the hyperbolic normal form `A=x*y+w`, the
leading pencil and parent faces make `B` affine in the other passive
variable. The next pencil coefficient leaves one active channel, and the
bordered unit makes that channel affine nonconstant. Exact integration gives
`B=x*z+rho*(y+h(x)*A)^2/2+beta(x)*y+gamma(x)*A+delta(x)`.
The checker verifies parent determinant `rho`, descendant determinant
`-kappa*rho`, and the triangular recovery of `x`, the displayed square
coordinate, `A`, `y`, `w`, and `z` from the descendant gradient. Thus every
rank-two quadratic-pivot descendant is a polynomial automorphism. The
generated ledger is
`artifacts/generated-results/hc4_quadratic_rank_two_pivots.json`.

Close the final passive three-by-three branch with:

```bash
.venv/bin/python scripts/verify_hc4_quadratic_rank_one_pivots.py
```

Together with
[`HC4_QUADRATIC_RANK_ONE_PIVOTS.md`](HC4_QUADRATIC_RANK_ONE_PIVOTS.md),
this supports `HC4RSD10`. Normalize the rank-one pivot to
`A=x^2/2+w`. The leading reduced-pencil and parent faces make the passive
three-variable Hessian `E` singular and impose
`a^T*adj(E)*a=0`. Passive ranks zero and two contradict the generic
corank-one bordered unit, so `rank(E)=1`. The rank-one polynomial-Hessian
normal form and the exact identity
`det Hess(Phi)=rho*det(a,d,ell)^2` turn the surviving factor into a unit
frame. Its Wronskian equation fixes the projective direction and gives
`B=x*z+rho*(y+h(x)*w)^2/2+alpha(x)*y+gamma(x)*w+delta(x)`.
The checker verifies the universal block faces, the frame identity, parent
determinant `rho`, descendant determinant `-kappa*rho`, and the triangular
recovery of all four variables. Hence all quadratic scalar pivots in the
identically singular reduced-pencil programme are collision-free. The
generated ledger is
`artifacts/generated-results/hc4_quadratic_rank_one_pivots.json`.

The same checker now continues into degree five.  On \(x_1=0\) it proves

\[
 D_R=\mathcal F(T,\partial T)-8LN^3S,
 \qquad T=R|_{x_1=0},\quad S=\partial_{x_1}R|_{x_1=0}.
\]

It follows that an \(x_1^2\)-divisible tail cannot repair any quartic graph
1-jet.  The leading quintic faces force
\(T_5=\kappa x_2^5\) and \(\partial_{y_2}T_4=0\).  For the complete v2 trace

\[
 T=\kappa x_2^5+d x_2^3y_1+\rho x_2^2y_2,
\]

exact resultants and first-transverse coefficients give a
characteristic-zero contradiction, even with every allowed off-plane term.
The checker also replays a rational graph that contains the marked collision
and has determinant \(17165601/25\) on all of \(x_1=0\), then detects its
nonzero coefficient \([x_1x_2^7]=22032/125\).  These are HC4MYGJ1,
HC4MYG5J, and HC4MYG5S; the general degree-five graph remains open.

The complementary rational top-cone and relative-linear checks are:

~~~bash
.venv/bin/python \
  scripts/verify_hc4_meng_yang_quintic_graph_normal_slice.py \
  --output \
  artifacts/generated-results/hc4_meng_yang_quintic_graph_normal_slice.json
.venv/bin/python \
  scripts/verify_hc4_meng_yang_quintic_q_kernel_slice.py \
  --output \
  artifacts/generated-results/hc4_meng_yang_quintic_q_kernel_slice.json
.venv/bin/python \
  scripts/verify_hc4_meng_yang_relative_linear_obstruction.py
~~~

The first forces the rational quintic top onto two constant-kernel charts;
the second substitutes the forced plane jet into the first
\(\partial_{y_2}\)-kernel slice with the complete degree-at-most-two lower
trace.  Its generic three immutable transverse equations generate the unit
ideal over \(\mathbb Q\); a square in the forced quartic normal jet kills the
\(y_2^2\) coefficient; and the exceptional nonzero \(x_2y_2\) branch ends in
two coprime transverse polynomials with resultant
\(986335129354383654912000\).  A separate resultant excludes the kernel
denominator chart.  Exact enumeration finds no admissible point modulo 101
or 103, hence no rational point to reconstruct.  As a diagnostic subfamily,
setting the \(y_1y_2\) coefficient to zero leaves two cubics with resultant
\(-108117004020524928=-2^7 3^{17}\cdot11\cdot13\cdot53\cdot863\).
The checker independently replays both transverse branches in the original
five-variable potential.  The third uses the external plane-Jacobian
degree-100 theorem to exclude collisions in every residual-linear correction
of degree at most 89.  The immediate graph target is the joint ideal obtained
from the forced plane normal jet, broader cubic and quartic traces, the
remaining top-cone charts, and first-transverse components above degree
three.  The generated JSON hashes for the projected normal slice and the
transverse kernel slice are,
respectively,
`cdf7fd3cb03dcaea616ce4e177ba87fddcffa423c5f29a0d5f6e4f5dc1e0fee5`
and
`dd28f7a44f6c813bdd422335133869ad7b4a513cc4426987930630c7cba859f9`.

The first bounded mixed canonical-pivot search is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/search_hc4_mixed_canonical_pivots.py \
  --output artifacts/generated-results/hc4_mixed_canonical_pivot_search.json
```

It searches 312 exact polynomial symplectic charts generated by quadratic
or cubic Hamiltonians in one mixed source--dual line, two commuting mixed
lines, or a cubic coisotropic-graph generator with one nonlinear constraint.
Pure source transformations are excluded by construction.  All 312 charts
have a scalar affine pivot and there are 258 jointly affine pairs.  Of 4320
specialized scalar-remainder trials, 240 are the exact inherited
`D=0` route and the other 4080 have modular nonconstancy witnesses.  Every
affine pair fails the simultaneous rank-at-most-two budget, and all 41796
complete descended determinants in the declared small repair box have
unequal values modulo `1000003`.  This is the finite-box result `HC4MCP1`,
not an exclusion of symbolic multi-parameter generators, mixed shear
compositions, coefficient-dependent repairs, or general coisotropic
embeddings.  See
[`HC4_MIXED_CANONICAL_PIVOT_SEARCH.md`](HC4_MIXED_CANONICAL_PIVOT_SEARCH.md).

The first genuinely compositional continuation is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/search_hc4_mixed_quadratic_words.py \
  --output artifacts/generated-results/hc4_mixed_quadratic_words.json
```

From 36 signed mixed quadratic Hamiltonian letters it forms all 1296
ordered words, keeps the 648 noncommuting words, removes 48 pure-source
cotangent maps, and deduplicates the remaining 600 exact linear symplectic
charts.  They have 1040 scalar affine pivots and 168 jointly affine pairs.
Exactly 864 scalar trials retain the inherited `D=0` mechanism.  Every pair
fails the rank-at-most-two budget and all 27216 complete determinants in the
same repair box are nonconstant by exact modular witnesses.  This is
`HC4MCP2`; words containing a cubic shear and general coefficient-dependent
repairs remain open.

The first word containing a cubic mixed shear is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/search_hc4_mixed_quadratic_words.py \
  --family quadratic-cubic \
  --output artifacts/generated-results/hc4_mixed_quadratic_cubic_words.json
```

It uses the 18 unit-time quadratic and 18 unit-time cubic mixed letters in
both orders.  Exact Poisson-bracket filtering and polynomial-map
deduplication leave 324 noncommuting nonlinear canonical words.  They have
576 scalar affine pivots and 108 jointly affine pairs, but no specialized
scalar remainder survives.  Every pair fails the corank budget and all
17496 complete repairs are nonconstant by exact modular witnesses.  The
post-gate audit also gives a nonconstant parent Hessian determinant in every
chart.  This is the normalized finite-box theorem `HC4MCP3`; signed flow
times, symbolic coefficient families, and coefficient-dependent repairs
remain open.

The fixed-order signed quadratic--cubic canonical box is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/search_hc4_mixed_quadratic_words.py \
  --family signed-quadratic-cubic \
  --output \
  artifacts/generated-results/hc4_canonical_signed_quadratic_cubic_words.json
```

It searches \(T_{H_2}\circ T_{H_1}\) with a signed quadratic \(H_1\) and
signed cubic \(H_2\).  Of 1296 raw words, 648 commute and are excluded.
The 648 noncommuting maps are distinct.  Exact support gives affine-pivot
dimension one for 432 words and dimension two for 216; the latter are
exactly the shared-dual words.  Their bracket-incidence census is 96 in
each one-sided type and 24 reciprocal words.  Every transformed reduced
Hessian pencil is generically rank four by an exact modular determinant
witness, every parent Hessian is nonconstant, and all 34992 complete
descended determinants are nonconstant in the declared repair box.  This
is the finite-box result `HC4MCP4`.  The box does not classify oblique
affine directions, symbolic coefficients, longer words, or
coefficient-dependent repairs.

The fixed-order symbolic coefficient closure is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/verify_hc4_symbolic_quadratic_cubic_words.py \
  --output \
  artifacts/generated-results/hc4_symbolic_quadratic_cubic_words.json
```

For each of the 54 noncommuting shared-dual support/sign patterns, it takes
\(H_1=aL_1^2\), \(H_2=bL_2^3\) over `Q[a,b]` and saturates by `a*b`.
Exact parent-Hessian determinant differences at integer probes give 14
localized monomial certificates and 40 unit standard bases in Singular.
Thus no pattern has a parent-constant specialization with `a*b != 0`.
This proves `HC4MCP5`, a coefficient-uniform parent obstruction for the
fixed degree order.  Zero coefficients, other Hamiltonian supports,
oblique directions, and longer words remain open.

The reverse-order symbolic coefficient family is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/verify_hc4_symbolic_quadratic_cubic_words.py \
  --order cubic-quadratic \
  --output \
  artifacts/generated-results/hc4_symbolic_cubic_quadratic_words.json
```

All 54 noncommuting shared-dual patterns now have an exact nonlinear
parent-preserving line.  The 48 one-sided cases force
`a = +/-1/2`; the six reciprocal cases force `a = +/-1/4`; in every case
`b` is arbitrary nonzero and the parent Hessian determinant is identically
`-16384`.  The checker verifies each complete six-variable identity in
Singular.  The retained coordinate affine-pivot dimensions are three in
24 patterns and two in 30 patterns.  All 102 constituent two-pivot pairs
have rank at least three and generic rank four for every `b != 0`, by exact
univariate minor and determinant certificates.  This is `HC4MCP6`.

The next unit quadratic--cubic commutator box is:

```bash
PYTHONHASHSEED=0 .venv/bin/python \
  scripts/search_hc4_mixed_commutator_words.py \
  --output \
  artifacts/generated-results/hc4_mixed_quadratic_cubic_commutators.json
```

It searches
`T_-H2 o T_-H1 o T_H2 o T_H1` for the 18 unit quadratic and 18 unit cubic
mixed-line letters.  Of 324 pairs, 162 commute.  The 162 noncommuting
commutator maps are distinct and split into 72 cases in each one-sided
Poisson-incidence type and 18 reciprocal cases.  An exact modular
chain-rule evaluation of the transformed Hessian gives unequal determinant
values at integer points for every word.  Thus every parent determinant is
nonconstant in characteristic zero and there are no survivors requiring
potential expansion.  This finite parent obstruction is `HC4MCP7`.

The smallest reciprocal mixed-line coisotropic graph box is:

```bash
.venv/bin/python \
  scripts/search_hc4_noncoordinate_coisotropic_scalar_gate.py \
  --output \
  artifacts/generated-results/hc4_noncoordinate_coisotropic_scalar_gate.json
```

For ordered `i != j`, it takes
`K=q_i+rho*p_j`, `L=q_j+rho*p_i`, and `H=tau*K*L^2`.
The commuting mixed forms make the time-one flow polynomial and send
`p_i=0` to a nonlinear reciprocal mixed graph.  In the boxes
`rho,tau in {-2,-1,1,2}`, `lambda in {-1,1}`, and
`mu in {-1,0,1}`, the 96 charts have 128 affine scalar pivots.  All 768
graph-specialized scalar Schur remainders have exact nonconstancy witnesses
modulo `1000003`.  The post-gate audit also proves that all 96 parent
Hessian determinants are nonconstant.  Thus no collision-transfer or full
descended-determinant calculation is reached.  This is the bounded result
`HC4MCP8`; arbitrary rational parameters and general coisotropic embeddings
remain open.  See
[`HC4_NONCOORDINATE_COISOTROPIC_GATE.md`](HC4_NONCOORDINATE_COISOTROPIC_GATE.md).

The direct one-variable calculation for the `PC(2)` graph is:

```bash
.venv/bin/python scripts/verify_hc4_direct_schur_descent.py
```

It first verifies that the nonsymmetric `PC(2)` Jacobian cannot literally be
a Hessian Schur complement.  It then classifies every coordinate chart and
omitted graph coordinate for an arbitrary polynomial auxiliary function of
`(X,Y,W,D)`, under transfer of a pair from the certified fiber.  Component
ideals, explicit fixed points, and a two-variable Jacobian-mate lemma leave
only charts `0010` and `0011`; in both, every slice is `W` up to scaling and
retained-coordinate gauge.  Their canonical generating families retain the
rational collision.  Irreducibility reduces every polynomial quadratic pivot
to three cases, and exact Hessian evaluations exclude constant nonzero
determinant in all six.  Finally, the checker proves that the Meng--Yang
five-variable potential has no further polynomial partial Legendre transform
along a constant linear direction.  Non-coordinate graph embeddings and
mixed coefficient-dependent critical equations are not tested.  Pure
univariate higher-degree repairs are excluded abstractly by `SDX1`.

The double-Schur audit for the parameterized quadratic-gauge families is:

```bash
.venv/bin/python scripts/verify_hc4_double_schur_gauge_obstruction.py
```

It tests two exact all-parameter routes from the six-variable Meng doubling
to four variables.  For elimination of two repaired dual variables, the
coefficient of the square of the retained dual variable is the bordered
Hessian invariant `K(L)`.  Three cubic coefficient equations and an
all-degree leading-layer calculation prove `K(L)` is nonzero for every
nonzero constant target linear form throughout the root-engineered gauge
family.  For a source-first descent, the common first coordinate forces the
only possible affine source direction to be `z`; higher-degree decorations
destroy that direction, while the cubic coefficient row has no constant
nonzero second pivot.  The calculation does not test nonlinear symplectic
changes, nonlinear retained dual coefficients, or nonconstant pivots with
exceptional divisibility.

The first nonlinear continuation, using triangular target coordinates, is:

```bash
.venv/bin/python scripts/verify_hc4_triangular_target_shears.py
```

After the canonical linear normalization of the full cubic gauge family,
it writes `L=W+H(U,V)` for each target permutation and a general polynomial
`H` of total degree at most three.  Sparse exact coefficients of the
bordered invariant `K(L)` first force the cubic part of `H` to vanish and
then give the quadratic contradiction.  In the third orientation the
coefficient `[x^8]K=9` is independent of every shear parameter.  Thus no
candidate reaches the full four-variable Hessian check.  The length-two
quadratic composition is treated next; degree at least four and general
nonlinear target or source--dual symplectic changes are not tested.

The length-two quadratic continuation is:

```bash
.venv/bin/python scripts/verify_hc4_two_quadratic_target_shears.py
```

For each of the six ordered words it writes
`A=X_i+Q(X_j,X_k)` and `L=X_j+R(A,X_k)`, with general positive-degree
quadratics `Q,R`.  Exact axis and transverse-line coefficients exclude the
`A^2` and `A*X_k` coupling strata over the reals.  If both couplings vanish,
the word is one quadratic triangular shear and is covered by the preceding
cubic single-shear checker.  This does not test degree-four single shears,
words with higher-degree factors, length at least three, general target
automorphisms, or nonlinear source--dual symplectic changes.

The preliminary bounded modular screen is:

```bash
.venv/bin/python scripts/search_hc4_two_quadratic_target_shears.py --prime 3
```

That command is an experiment, not a proof; the exact checker above is the
status-bearing calculation.

The degree-four single-shear closure is:

```bash
.venv/bin/python scripts/verify_hc4_quartic_target_shears.py
```

It includes every lower-degree shear coefficient and uses five extreme
spatial coefficients to force the homogeneous quartic layer to vanish in
the `P`- and `B`-retained orientations.  The cubic checker then applies.  In
the `C`-retained orientation, `[x^8]K=9` remains independent of every
coefficient through degree four.  This exact characteristic-zero result
does not test degree at least five, words with cubic-or-higher factors,
length at least three, general target automorphisms, or nonlinear
source--dual symplectic changes.

The all-order binary-form symbol underlying the two \(X\)-caustic normal-jet
recursions is documented in
[`HC4_PC2_GRAPH_POLARIZATION_AUDIT.md`](HC4_PC2_GRAPH_POLARIZATION_AUDIT.md).
Its finite exact regression in both charts is:

```bash
.venv/bin/python scripts/verify_hc4_logarithmic_normal_symbol.py
```

The checker verifies through tensor order ten that the order-\(n\) symbol is
\(L\) times a unimodular coefficient-extraction matrix.  The all-order
statement is proved by the triangular determinant formula in the note.  This
establishes a divisor-supported graded jet module, not yet a finite-rank
logarithmic connection or a Bernstein--Sato polynomial.

The nonlinear toric descent of the Meng--Yang potential is:

```bash
.venv/bin/python scripts/verify_hc5_nonlinear_toric_descent.py
```

It constructs an explicit determinant-one dual-coordinate change with
`t=A`, verifies the polynomial unit-pivot critical solution, and computes the
natural four-variable determinant `16*J(x1*x2)^2`.  A relative toric `SL(2)`
correction cancels the nonconstant factor and gives determinant `64`, but the
descended gradient then has an explicit polynomial inverse.  The checker
also proves the all-degree toric radial obstruction: constant nonzero
determinant forces both radial factors to be units, so this class cannot
retain the Meng--Yang collision.  Non-toric changes and non-coordinate
coisotropic embeddings remain open.

The bounded non-toric relative correction is:

```bash
.venv/bin/python scripts/verify_hc4_nontoric_sl2_correction_degree4.py
```

For a general matrix `C(x,y) in SL(2)` whose four entries have total degree
at most four, it sets `G=beta*C` for the natural complementary coefficient
row.  The four-variable identity is
`det Hess(psi)=16*det(DG)^2`.  Singular proves that the 45 determinant-one
equations and 218 nonconstant-Jacobian equations generate the unit ideal
over `QQ`; collision equality is not needed.  A second exact calculation
excludes arbitrary affine `SL(2)` perturbations of the known degree-ten
toric correction after the collision equations are imposed.  Raw degree at
least five, quadratic-or-higher perturbations of the toric correction,
mixed base--dual changes, and non-coordinate coisotropic embeddings remain
open.  This checker requires `Singular` on `PATH`.

The shortest exact replay of the canonical Meng descent chain
`HC5T1 -> HC4MQ1 -> HC4MCK` is:

```bash
.venv/bin/python scripts/verify_hc5_nonlinear_toric_descent.py
.venv/bin/python scripts/verify_hc4_meng_sparse_quartic_obstruction.py
.venv/bin/python scripts/verify_hc4_meng_full_cubic_kernel.py
```

The final command requires `Singular` on `PATH`.  The support-three and
support-four commands below are exact historical checkpoints and targeted
regressions, not logical prerequisites for `HC4MCK`.

The first collision-first non-toric Hamiltonian screen is:

```bash
.venv/bin/python scripts/verify_hc4_meng_sparse_quartic_obstruction.py
```

After the unit-pivot descent and a polynomial base gauge, the transported
points are antipodal for `psi_0=2*y*r+4*x*s`.  The checker exhausts every
homogeneous quartic correction supported on at most four monomials, imposing
the collision before determinant work.  It rejects 42,953 isolated collision
solutions and 515 one-parameter families exactly.  The degree-eight
principal part then leaves only 232 isolated quartics and two exact family
members; none admits a constant-determinant correction by one cubic monomial.
For two cubic monomials, degrees seven and one give a linear rank gate;
rank-zero cases pass through degree-six conic linearization and a
degree-four/two lift.  Only four bivariate families reach the terminal
calculation, and their full determinant ideals are units modulo `1000003`.
The complete homogeneous-quartic reduction is:

```bash
.venv/bin/python scripts/verify_hc4_meng_full_quartic_reduction.py
```

It verifies an explicit complex congruence taking
`2*y*r+4*x*s` to one half the sum of four squares, checks the gradient
chain rule on all 35 quartic monomials, and transports the antipodal
collision.  The transformed gradient is `z+grad(h4)`, with homogeneous
cubic nonlinear part and symmetric Jacobian.  The external
de Bondt--van den Essen dimension-four theorem therefore makes every
constant-Jacobian member a polynomial automorphism.  This closes the
complete homogeneous quartic chart, and the same congruence argument closes
every nondegenerate quadratic renormalization plus one homogeneous quartic.
Mixed nonlinear degrees are not covered.

The complete dense cubic--quartic continuation is:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_dense_cubic_quartic_reduction.py
```

For `psi=q2+h3+h4`, homogeneous determinant layers and
Gordan--Noether give a constant direction satisfying
`D_v(h4)=D_v^2(h3)=0` in every rank of `Hess(h4)`.  The rank-one quartic
case uses an additional exact ternary Hessian-pencil calculation.  A
nonisotropic direction descends to `HC(3)`; the isotropic bordered form
reduces either to a degree-at-most-three plane Keller map or to an `HC(2)`
block after a binary quadratic invariant forces `s=constant+x`.  This
excludes all 20 cubic and all 35 quartic coefficients simultaneously,
without a support bound.  The external inputs are Gordan--Noether,
`HC(3)`, `HC(2)`, and Moh's plane degree bound.

The independent dense mixed-quartic coefficient regression is:

```bash
.venv/bin/python scripts/verify_hc4_meng_dense_mixed_quartic.py
```

It treats all 25 homogeneous quartic monomials of base--dual bidegrees
`(1,3)`, `(2,2)`, and `(3,1)`, then separately adjoins every pure-base
quartic and every pure-dual quartic.  Collision and the linear determinant
layer have rank 14, leaving 11, 16, and 16 parameters.  Singular proves that
the remaining exact coefficient ideals, with 262, 273, and 273 generators,
are unit ideals over `QQ`.  The complete homogeneous-quartic theorem
subsumes this one-sided coefficient calculation, and the dense
cubic--quartic theorem subsumes the later sparse cubic-kernel chain.
Quartic--sextic, simultaneous cubic--quartic--sextic, and non-coordinate
embeddings remain open.  This command requires `Singular` on `PATH`.

The finite-field continuation through exactly three cubic monomials is:

```bash
.venv/bin/python scripts/verify_hc4_meng_three_cubic_rank_gate.py
```

It checks all `234*binomial(20,3)=266760` quartic/triple pairs.  The combined
degree-seven/degree-one rank census is `5480, 53364, 130508, 77408` in ranks
zero through three.  After delegating support-at-most-two boundary loci to
the preceding certificate, unit full-determinant gcds or ideals exclude all
920 genuine rank-two lines, 2,952 genuine rank-one planes, and 5,480
rank-zero three-parameter spaces over `F_1000003`.  This is an exact
finite-field computation.

The characteristic-zero promotion is:

```bash
.venv/bin/python scripts/verify_hc4_meng_three_cubic_characteristic_zero.py
```

It reconstructs all 234 quartics over `QQ`, reproduces the same odd-layer
rank census over `QQ`, and obtains unit gcds or Gröbner ideals for every
genuine line, plane, and three-parameter space.  The maximum evaluation
prefix lengths are respectively 5, 7, and 9 points.  Thus cubic support at
most three is excluded in characteristic zero.  This bounded-support
checkpoint is subsumed by the full cubic-kernel theorem below.

The finite-field continuation through exactly four cubic monomials is:

```bash
.venv/bin/python scripts/verify_hc4_meng_four_cubic.py
```

Its two targeted stages are:

```bash
.venv/bin/python scripts/verify_hc4_meng_four_cubic_rank_gate.py
.venv/bin/python scripts/verify_hc4_meng_four_cubic_rank_zero.py
```

The first checker exhausts all
`234*binomial(20,4)=1133730` quartic/quadruple pairs.  The odd-layer ranks
zero through four are `5430, 79396, 353740, 504818, 190346`.  After
delegating support-at-most-three boundaries to `HC4MC3`, unit gcds or ideals
exclude all 466 genuine rank-three lines, 6,082 genuine rank-two planes, and
7,956 genuine rank-one three-spaces.  The second checker isolates the 5,430
rank-zero four-parameter spaces and proves that every full determinant ideal
is a unit modulo `1000003`, using at most twelve evaluation points.

The characteristic-zero promotion is:

```bash
.venv/bin/python scripts/verify_hc4_meng_four_cubic_characteristic_zero.py
```

It reconstructs the same rank and boundary census over `QQ`.  Exact rational
gcds or Gröbner ideals exclude all 466 lines, 6,082 planes, 7,956
three-spaces, and 5,430 four-spaces.  The maximum evaluation-prefix lengths
for ranks three through zero are `5, 7, 8, 11`.  Thus cubic support at most
four is excluded in characteristic zero.  This bounded-support checkpoint is
also subsumed by the full cubic-kernel theorem below.

The full cubic-kernel characteristic-zero checker is:

```bash
.venv/bin/python scripts/verify_hc4_meng_full_cubic_kernel.py
```

It requires `Singular` on `PATH`.  For all 234 rational quartics, it
constructs the complete odd determinant kernel in the 20-dimensional cubic
space, symbolically extracts every spatial coefficient of the Hessian
determinant, clears denominators, and adds the coefficient ideals in
descending spatial degree.  The exact odd ranks are 8 for 229 quartics, 7
for one, and 4 for four exceptional `u^3*L` quartics, giving kernel
dimensions 12, 13, and 16.  Singular reaches the unit ideal first at degree
six for 62 quartics, degree five for 16, and degree four for 156.  Therefore
arbitrary homogeneous cubic corrections are excluded in characteristic
zero whenever the collision quartic has support at most four.

The parallel sparse sextic collision-carrier checker is:

```bash
.venv/bin/python scripts/verify_hc4_meng_sparse_sextic_obstruction.py
```

It exhausts all supports of at most four among the 84 homogeneous sextic
monomials.  Collision-first linear algebra leaves 1,725,838 isolated points,
7,566 lines, and one plane.  Five leading-Hessian evaluations reduce the
isolated points to 748 candidates, all rejected by the first two full
determinant evaluations.  Exact rational gcds leave only two principal
lines, at parameters `-81/8` and `9/2`; determinant degree twelve rejects
both.  The unique binary-cubic plane in `y*r` and `x*s` has unit principal
coefficient ideal over `QQ`.  Thus no sextic-only correction supported on
at most four monomials works.  The later joint theorem `HC4JQS4` treats
mixed quartic--sextic corrections of combined support at most four, and
the later support-free theorem `HC4E46` subsumes both calculations.

The characteristic-zero continuation over the 234 quartic principal parts
is:

```bash
.venv/bin/python scripts/verify_hc4_meng_mixed_quartic_sextic.py
```

Exactly four of the 234 quartic principal parts have zero immutable
determinant-degree-two signature.  In their homogeneous zero-gradient
sextic kernel, the checker finds 976 three-support lines, 205,494
four-support lines, and 519 planes.  Principal cancellation leaves 121,146
lines.  For each quartic, 52,686 lines have projectively complete modular
certificates and the remaining 68,460 have unit exact rational gcds.
Separate exact affine-matrix expansion proves that all `519*4=2076`
three-evaluation plane ideals are units over `QQ`.  Thus no zero-gradient
sextic supported on at most four monomials repairs any of the 234 quartic
principal parts.  The later `HC4JQS4` theorem treats quartics outside that
sextic-free principal screen when the combined support is at most four;
the larger mixed chart and simultaneous cubic corrections remain open.

The joint total-support-four continuation is:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_joint_quartic_sextic_total_support_four.py
```

It exhausts all 6,133,820 genuinely mixed quartic--sextic supports of
combined size at most four.  Collision-first linear algebra leaves
5,225,684 isolated points, 44,300 lines, and 34 planes.  Descending
determinant-layer evaluations reject every isolated point modulo the
rank-preserving prime `1000003`.  Every line is reconstructed over `QQ`
and has unit exact evaluation gcd; every plane has unit exact
three-evaluation ideal.  Together with the pure quartic theorem `HC4MQ1`
and pure sextic theorem `HC4MS6`, this proves `HC4JQS4`: no
quartic-plus-sextic correction of total support at most four works when
there is no cubic term.  The later `HC4E46` theorem removes the support
bound completely; simultaneous cubic and sextic corrections remain open.

The support-free common-kernel continuation is:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_dense_rank_three_sextic_reduction.py
```

It checks the determinant identities in theorem `HC4DCK`.  The theorem
excludes every support-free correction having a constant direction \(v\)
with `D_v h6=0` and `D_v^2 h4=0`.  Gordan--Noether and the
degree-fourteen determinant layer make this automatic for a rank-three
sextic Hessian.  A non-isotropic quadratic pivot reduces to `HC_3`.  In the
isotropic bordered chart, the next two coefficients reduce either to a
degree-at-most-five plane Keller map or to a two-variable constant-Hessian
block.  Moh's degree bound and `HC_2` exclude both.  Only sextic Hessian
rank at most two with a variable quartic null direction is left by
`HC4DCK` itself; theorem `HC4E46` below closes it.  Simultaneous cubic and
sextic corrections remain open.

The source/dual reorganization of that rank-two boundary is:

```bash
.venv/bin/python scripts/verify_hc4_source_dual_bigrading.py
```

It checks the complete \(2+2\) bidegree ledger, the corrected weighted
Hessian-face identity, and the cotangent determinant
`det Hess(t*F+m*G+H)=Jac(F,G)^2`.  The canonical note proves theorem
`HC4SDW`: the dual-linear stratum is exactly a degree-at-most-five
`JC(2)` locus, and every vanishing sequence of rank-two residual Hessian
faces synchronizes by a Schur recursion to one rational projective cone.
The rotating example `(x*t+y*m)^2` is synchronized but nonconstant, so
the synchronization theorem alone leaves a nonlinear moving-cone
algebraization problem.

The support-free closure of the full even quartic--sextic chart is:

```bash
.venv/bin/python scripts/verify_hc4_even_quartic_sextic_closure.py
```

The primitive cone-degree lemma shows that bihomogeneity leaves only
`c*(X^T*M*U)^2` as a nonconstant rank-two residual cone.  A rank-one `M`
has a constant dual direction and is `HC4DCK`; for invertible `M`, the
dual-degree-four part of the spatial `z^4` determinant coefficient is

```text
48*c^4*det(M)^2*(X^T*M*U)^4,
```

and cannot be cancelled by the quadratic block, a dual-linear/source-only
quartic block, or the source-only sextic block.  The rank-one sextic
boundary becomes constant by one-base scaling, rank three is `HC4DCK`,
and rank zero is `HC4HQ1`.  This proves `HC4E46`: no support-free
potential `q2+h4+h6` with nondegenerate `q2` has both constant nonzero
Hessian determinant and an antipodal collision.  Only simultaneous cubic
and sextic interaction, and non-coordinate coisotropic embeddings, remain.

The rank-three part of the simultaneous cubic--quartic--sextic chart is:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_triple_rank_three_reduction.py
```

For `psi=q2+h3+h4+h6`, the degree-fourteen and degree-thirteen determinant
layers on the one-dimensional sextic kernel force both
`D_t^2(h4)=0` and `D_t^2(h3)=0`.  The nonisotropic direction descends to
`HC(3)`.  In the isotropic chart, the checker reconstructs the complete
cubic bordered invariant.  The two rank-two binary-cubic orbits have
coefficient ideal `(qxm,qym,lm)^2`; the rank-one orbit has two explicit
radical branches, each with a constant missing direction.  The remaining
descent is either a plane Keller cotangent lift of degrees at most three
and five or an `HC(2)` block.  This proves `HC4T31` without a support
restriction.  Only simultaneous corrections with
`rank Hess(h6)<=2`, and non-coordinate coisotropic embeddings, remain.

The rank-two continuation is:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_triple_rank_two_reduction.py
```

The degree-twelve face makes the binary quartic Hessian on the constant
sextic kernel plane singular.  A constant cone immediately recovers the
common-direction reduction; if that binary Hessian is zero, the
degree-ten face makes the cubic binary Hessian singular, and total degree
three forbids a moving cone.  In the only moving quartic case, the
degree-eleven face forces the high-dual cubic to align as
`(X^T*M*U)*(alpha*t+beta*m)`.  The checker keeps this cubic and all
compatible lower blocks and proves that the later dual-degree-four face
is still the uncancellable
`48*c^4*det(M)^2*(X^T*M*U)^4`.  This proves `HC4T21`.

The rank-one continuation and complete coordinate-chart exhaustion are:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_triple_rank_one_reduction.py
```

The small-rank Hessian normal form makes the sextic a sixth power
`c*L^6`, with constant three-dimensional kernel `W`.  The degree-ten
face makes `Hess_W(h4)` singular.  Its rank-two branch aligns the cubic
at degree nine.  In rank one, the degree-eight face is exactly the binary
discriminant `a*d-b^2` of `Hess(h3)` on the constant quartic-kernel
plane; this is the `Sym^2` nullcone `L^2` from the SIC(2) binary-root
classification.  In rank zero, degree seven makes `Hess_W(h3)` singular.
One-base homogeneity makes every residual projective cone constant, so
each branch reaches the common-direction descent `HC4T31`.  This proves
`HC4T11`; ranks three, two, one, and zero then give `HC4TC1`, the complete
support-free obstruction for `q2+h3+h4+h6`.  Quintic and higher
homogeneous layers, as well as non-coordinate coisotropic embeddings,
remain outside this chart.

The first quintic extension of the common-direction descent is:

```bash
.venv/bin/python scripts/verify_hc4_quintic_common_direction.py
```

The checker proves the quartic bordered lemma `HC4BL4` by the binary-root
stratification used in SIC(2), including exact radical certificates for
the double-double and pure-fourth exceptional strata.  It then proves
`HC4CD5`: for `psi=q2+h3+h4+h5+h6`, a common direction satisfying
`D_v h6=0` and `D_v^2 h5=D_v^2 h4=D_v^2 h3=0` reduces by two Schur
steps to `HC(3)`, Moh's plane theorem, or `HC(2)`.  In sextic Hessian rank
three this closes the branch `D_v h5=0`.  The first remaining quintic
face is the exact divisibility

```text
det(Cbar) | grad(D_v h5)^T * adj(Cbar) * grad(D_v h5),
```

where `Cbar` is the nondegenerate ternary block of `Hess(h6)`.

The full Fermat-sextic diagonal Schur-norm stratum is checked by:

```bash
.venv/bin/python scripts/verify_hc4_quintic_diagonal_schur.py
cd formal/finite-etale-keller
lake env lean FiniteEtaleKeller/HC4QuinticDiagonal.lean
```

The first command retains an arbitrary base quintic and every lower
quartic, cubic, and quadratic coefficient.  Exact truncated determinant
expansion shows that the `lambda^13*t*x^4*y^4*z^4` face fixes the `t^3`
coefficient of `h3`, after which three `lambda^11*t^3` coefficients are
`1024*a^5`, `1024*b^5`, and `1024*c^5`.  Before that prolongation, an
exact 66-equation radical certificate proves that every quartic with
polynomial Fermat Schur norm is diagonal.  Thus `HC4QF1` closes the entire
Fermat-sextic quintic stratum.  Lean checks
the scalar Schur identity and the characteristic-zero fifth-power
conclusion; determinant coefficient extraction remains the Python
certificate.

The first full non-diagonal sextic pencil is:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_symmetric_sextic_schur.py
```

For `h6=(x^6+y^6+z^6)/30+mu*x^2*y^2*z^2`, the checker constructs the
complete 111-equation Schur-divisibility ideal for a generic 15-coefficient
quartic and six-coefficient quadratic quotient.  Saturation by `mu` has an
exact 261-element rational Gröbner basis.  The saturated ideal lies in the
21-variable coefficient origin, and the fourth powers of all 21 generators
reduce to zero.  Thus `HC4QS1` closes every `mu!=0` member; `HC4QF1` closes
`mu=0`.  The fourth-power reduced-ring endpoint is replayed in
`FiniteEtaleKeller/HC4QuinticDiagonal.lean`.

The two-parameter generic broadening is:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py
```

For
`h6=(x^6+y^6+z^6)/30+mu*x^2*y^2*z^2+nu*sum_(i!=j)x_i^4*x_j^2`,
six quotient pivots have determinant `4096*nu^12`.  After eliminating the
quadratic quotient and clearing `2*nu^2`, 114 intrinsic equations remain.
Over `Q(mu,nu)` their exact 117-element Gröbner basis contains the cubes of
all fifteen quartic coefficients.  This proves `HC4QSG2`, generic rigidity
on the two-parameter surface.  It does not exclude exceptional curves
inside `nu!=0`: exact uniform saturation timed out at 300 seconds, and
raw projective-chart runs timed out at 120 seconds.  Those timeouts are
diagnostics, not certificates.

The surface nevertheless contains an explicit exceptional Schur pair:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_radial_exceptional_schur.py
```

At `(mu,nu)=(1/5,1/10)`, the sextic is
`(x^2+y^2+z^2)^3/30`.  The checker proves that the nonzero quartic
`s4=(x^2+y^2+z^2)^2` has polynomial Schur quotient
`16*(x^2+y^2+z^2)` and that the sextic Hessian determinant is
`(x^2+y^2+z^2)^6/25`.  This verifies only the Schur face; it does not
by itself claim extension through the lower collision identities.  In
the invariant coordinates `R=x^2+y^2+z^2`,
`P2=x^2*y^2+x^2*z^2+y^2*z^2`, and `P3=x^2*y^2*z^2`, the surface is
`R^3/30+(nu-1/10)*R*P2+(mu-3*nu+1/10)*P3`.  The checker also proves that
the fixed radial quartic `R^2` has a polynomial Schur norm only when both
deformation coefficients vanish, so this pair occurs only at the stated
parameter point.

The lower-face prolongation is excluded by:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_radial_prolongation.py
```

The checker retains an arbitrary base quintic and every allowed quartic,
cubic, and quadratic coefficient.  Sparse exact determinant extraction
gives `2*(3*delta-32)/25` at `lambda^13*t*x^12`, forcing
`delta=32/3`, while `lambda^11*t^3*x^8` then equals `1024/25`.
Every arbitrary lower coefficient cancels.  Thus the radial Schur pair is
exceptional at the Schur face but cannot produce a constant-Hessian
collision.

The full even permutation-invariant quartic line is classified by:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_even_symmetric_schur.py
```

For `s4=a*R^2+b*P2`, the `a=1` chart has radical support at exactly two
points: the radial pair and the Fermat pair
`(mu,nu,s4)=(0,0,x^4+y^4+z^4)`.  The `a=0` chart is the unit ideal.
Consequently this complete symmetric-even slice contains no exceptional
parameter curve, and both isolated points are already excluded at lower
faces.

The component-directed exceptional-locus research transcript is:

```bash
.venv/bin/python scripts/verify_hc4_exceptional_schur_atlas.py
.venv/bin/python scripts/research_hc4_exceptional_schur_locus.py
.venv/bin/python scripts/research_hc4_exceptional_schur_locus.py \
  --exact-pure-chart --singular-timeout 900
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --extract-basis-denominators --basis-profile
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --cube-torsion-profile --cube-index 0 --cube-torsion-stage generic
```

The first command is the exact full reduced atlas `HC4QSE4`.  It builds
the 120-equation projective incidence system, expresses the Hessian
determinant in the seven symmetric invariants, and proves that its
nonsquarefree locus consists only of Fermat and radial.  The apparent
third discriminant point `(5/7,-1/14)` has two distinct coprime Hessian
factors.  Hence the reduced incidence is the Fermat projective plane of
diagonal quartics plus the single radial quartic `R^2`.  The lower-face
checkers `HC4QF1` and `HC4QSE2` exclude both before additional antipodal
collision equations can contribute.

The later commands are the earlier modular reconstruction,
special-fiber geometry, and denominator research recorded in
[`HC4_EXCEPTIONAL_SCHUR_LOCUS.md`](HC4_EXCEPTIONAL_SCHUR_LOCUS.md) and
[`hc4_exceptional_schur_locus_modular.json`](artifacts/generated-results/hc4_exceptional_schur_locus_modular.json).
The even-quartic charts reconstruct the Fermat and radial parameter
points modulo \(47,101,103\).  The transformation-aware
`--extract-denominators` calculation timed out at its 900-second Singular
bound and supplies no certificate.  The basis-only denominator is the
constant \(2\), so the desired exceptional divisor must occur in the lift
certificates.  The displayed first sign-character-block cube lift also
timed out at 900 seconds; it has 191 target cubic monomials and 441
multiplication columns.  The direct characteristic-zero `a=1`
even-quartic elimination likewise reached its 900-second bound before
returning a standard basis.  Those interrupted routes are historical
diagnostics and are not used by `HC4QSE4`.

The degree-three cube-torsion research modes are:

```bash
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --cube-torsion-profile --cube-index 2 \
  --cube-torsion-stage finite-field --cube-prime 19
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --cube-torsion-profile --cube-index 2 \
  --cube-torsion-stage specialize \
  --cube-mu-value=-5/3 --cube-nu-value=-1/6
.venv/bin/python \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --cube-torsion-profile --cube-index 2 \
  --cube-torsion-stage fiber \
  --cube-mu-value=-5/3 --cube-nu-value=-1/6
```

The construction uses the canonical multiplication presentation
`A^1710 -> A^680` and the fifteen cube targets, not a noncanonical
`A^114 -> A^15` compression.  Complete scans of the four
coefficient-monomial orbits modulo `11,13,17,19` reconstruct the radial
point and the additional point `(-5/3,-1/6)` on `nu!=0`.  Exact rational
specialization shows that precisely the three `x_i^2*x_j^2` coefficient
cubes survive at the new point, while every fourth power is zero; the
60-dimensional fiber is therefore supported at the coefficient origin.
This is a nilpotence-order jump, not a reduced exceptional Schur pair.
The even-block integral annihilator and function-field lift each reached a
900-second timeout; relation extraction and the zeroth Fitting ideal were
not reached.  See
[`HC4_FITTING_DENOMINATOR_EXTRACTION.md`](HC4_FITTING_DENOMINATOR_EXTRACTION.md)
and
[`hc4_fitting_denominator_extraction.json`](artifacts/generated-results/hc4_fitting_denominator_extraction.json).

The full-coefficient reduced-fiber scans are:

```bash
for prime in 7 11 13; do
  .venv/bin/python -u \
    scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
    --fourth-power-profile --fourth-prime "$prime"
done
```

The corresponding symbolic attempt is:

```bash
.venv/bin/python -u \
  scripts/verify_hc4_quintic_two_parameter_symmetric_schur.py \
  --fourth-power-profile --fourth-prime 7 \
  --fourth-stage annihilator --fourth-timeout 900
```

The zero sign-character block of degree four has 819 target monomials and
3474 equation-times-quadratic columns.  At all \(42+110+156\) parameter
points on `nu!=0`, membership of all fifteen coefficient fourth powers
certifies an empty reduced projective Schur fiber except at the radial
reduction.  The cube-torsion point `(-5/3,-1/6)` is certified
reduced-empty at every prime.  These scans cover \(\mathbb F_p\)-rational
parameters, not points over proper finite-field extensions and not the
characteristic-zero support ideal.  The symbolic \(\mathbb F_7[\mu,\nu]\)
annihilator command reached its 900-second bound before returning a
standard basis and supplies no factorization.  See
[`hc4_fourth_power_support.json`](artifacts/generated-results/hc4_fourth_power_support.json).

The direct collision-normalized finite-field experiment in degree bounds
five through eight is:

```bash
.venv/bin/python scripts/search_hc4_finite_field_potentials.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-bound 2 \
  --points 6 \
  --output artifacts/generated-results/hc4_finite_field_sparse_search.json
```

It exhausts one- and two-direction affine perturbations in the full linear
kernel of the normalized gradient-collision condition, totaling 45,181,194
coefficient choices.  No exact modular candidate survives.  This is a
bounded experiment, not evidence for unrestricted `HC_4`; its construction
and precise scope are in
[`HC4_FINITE_FIELD_SEARCH.md`](HC4_FINITE_FIELD_SEARCH.md).

The denser sampled-support coefficient solve is:

```bash
.venv/bin/python scripts/search_hc4_finite_field_dense_supports.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --strategies uniform homogeneous mixed \
  --trials 2 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_dense_support_search.json
```

It forms the complete determinant coefficient ideal for 96 deterministic
supports of 6, 8, 10, or 12 collision-kernel directions.  Adjoining
`a_i^p-a_i` makes each Singular calculation an exact existence test over the
selected prime field.  All 192 support-prime ideals are unit ideals, with no
timeouts.  The support selection is sampled, so this remains a bounded
experiment.

The structurally guided companion forces every sampled support to contain as
many directions as possible that can alter both base determinant defects on
the normalized collision axis:

```bash
.venv/bin/python scripts/search_hc4_finite_field_dense_supports.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --strategies axis \
  --trials 2 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_axis_support_search.json
```

All 64 ideals from these 32 additional supports are also unit ideals.

The principal-Hessian companion makes every degree-\(d\) correction a
three-variable cone, so its top homogeneous Hessian determinant vanishes
identically, then adds lower-degree monomials involving the omitted fourth
variable:

```bash
.venv/bin/python scripts/search_hc4_finite_field_dense_supports.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --strategies cone2 cone3 \
  --trials 4 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_cone_bridge_search.json
```

All 256 full coefficient ideals from these 128 additional supports are unit
ideals.  Thus the failure occurs below the principal homogeneous Hessian
gate, but the support search remains bounded.

The non-coordinate cone search uses
`u=x2+lambda*x3`, `v=x3` for `lambda=-1,1,2`.  In the adapted coordinates
the quadratic base is `x0*x1+u*v-lambda*v^2`, the top correction depends only
on `x0,x1,u`, and the lower bridges involve `v`:

```bash
.venv/bin/python scripts/search_hc4_oblique_cone_bridges.py \
  --degrees 5 6 7 8 \
  --primes 11 13 \
  --support-sizes 6 8 10 12 \
  --slopes -1 1 2 \
  --trials 3 \
  --timeout 30 \
  --output artifacts/generated-results/hc4_finite_field_oblique_cone_bridge_search.json
```

All 288 exact ideals from the 144 oblique families are unit ideals, with no
timeout.  The adapted change has determinant one, so the Hessian determinant
identity is equivalent to the tied-monomial identity in the original
coordinates.

## Proof-carrying arithmetic compiler

The active common-arithmetic-fibers paper has a separate correspondence
compiler for its displayed Berend--Bilu example. One compact JSON
specification generates the Lean-readable coefficients and map, TeX macros
used directly by the manuscript, runnable SymPy and Sage inputs, and an
expanded sparse JSON certificate. The verifier checks that every generated
view is current, runs the SymPy input, and builds the Lean theorems equating
the generated paper polynomial, map, inverse polynomial, output scalings, and
distinguished targets with the existing formal definitions:

```bash
make verify-common-arithmetic-fibers-correspondence
```

Intentional regeneration is:

```bash
make refresh-common-arithmetic-fibers-example
```

The canonical source is
`papers/common-arithmetic-fibers/data/explicit-quintic-spec.json`; generated
results are recorded under `artifacts/generated-results/`, with the
Lean-readable module at
`formal/finite-etale-keller/FiniteEtaleKeller/GeneratedPaperExample.lean`.

The pinned ramified-quintic local-field specification compiles to a portable
minimal-gauge JSON proof object and a Lean specialization of the formal
algebra-to-Keller theorem.  Two further JSON objects certify the same
ramified quintic in the power-shift family at `m=2` and the connected cubic
`T^3-T-1` in the cubic family at `n=7`.  Replay all three maps with two
independent arithmetic implementations and check all three generated Lean
specializations with:

```bash
python3 scripts/verify_arithmetic_keller_certificate.py
gp -q -f scripts/verify_arithmetic_keller_certificate.gp
python3 scripts/verify_arithmetic_keller_certificate.py \
  artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json
ARITHMETIC_CERTIFICATE=artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json \
  gp -q -f scripts/verify_arithmetic_keller_certificate.gp
python3 scripts/verify_arithmetic_keller_certificate.py \
  artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json
ARITHMETIC_CERTIFICATE=artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json \
  gp -q -f scripts/verify_arithmetic_keller_certificate.gp
cd formal/finite-etale-keller
lake env lean FiniteEtaleKeller/GeneratedArithmeticQuintic.lean
lake env lean FiniteEtaleKeller/GeneratedArithmeticQuinticStableM2.lean
lake env lean FiniteEtaleKeller/GeneratedArithmeticCubicStableN7.lean
```

From the repository root, the combined non-mutating command is:

```bash
make verify-arithmetic-compilation
```

Intentional regeneration is separate:

```bash
.venv/bin/python scripts/compile_arithmetic_keller_certificate.py
.venv/bin/python scripts/compile_arithmetic_keller_certificate.py \
  --stable-parameter 2 \
  --certificate artifacts/generated-results/arithmetic_keller_quintic_stable_m2.json \
  --lean-module FiniteEtaleKeller.GeneratedArithmeticQuinticStableM2 \
  --lean formal/finite-etale-keller/FiniteEtaleKeller/GeneratedArithmeticQuinticStableM2.lean
.venv/bin/python scripts/compile_arithmetic_keller_certificate.py \
  --spec arithmetic/specifications/connected_cubic_stable_n7.json \
  --certificate artifacts/generated-results/arithmetic_keller_cubic_stable_n7.json \
  --lean formal/finite-etale-keller/FiniteEtaleKeller/GeneratedArithmeticCubicStableN7.lean
```

The format and exact claim boundary are documented in
[`arithmetic/PROOF_CARRYING_COMPILATION.md`](arithmetic/PROOF_CARRYING_COMPILATION.md).

## Support-saturation compiler

The reusable module compiler computes
\((N:I^\infty)/N=H_I^0(F/N)\), module associated primes, exact regular
elements, distinguished-class annihilators and radicals, finite normal-jet
transitions, and finite-tower uniform-exponent tests.  The shared JSON schema
also records the completion ideal, parameter/base variables, normal
variables, and exact-versus-modular assurance.  Compile the checked example
with:

```bash
.venv/bin/python scripts/compile_support_saturation.py \
  schemas/support_saturation_example.json \
  --output /tmp/support_saturation_certificate.json
```

Its fast exact calibration is:

```bash
make verify-support-saturation-compiler
```

This target also runs the standard-library-only replay of the exact
characteristic-zero degree-42 \(c_6\) certificate.  It verifies
\[
c_6\notin J_6,\qquad w_0c_6,w_2c_6\in J_6
\]
on the specialized Ritt fiber, using explicit rational multipliers and a
finite-support Macaulay dual functional.  The narrower command is:

```bash
make verify-degree42-c6-macaulay
```

Intentional modular block-Wiedemann reconstruction over the pinned 31-bit
primes is:

```bash
make refresh-degree42-c6-macaulay
```

Intentional regeneration of the homogeneous cubic-symbol atlas, the exact
degree-42 finite-jet computation over `GF(32003)`, and the normalized
plane-JC cyclic `d3` boundary layer is:

```bash
make refresh-support-saturation-cases
```

Regenerate only the cubic search stratification imported from the proved
formal-gauge cokernel atlas with:

```bash
.venv/bin/python scripts/compile_support_saturation_cases.py --case cubic-frontier
```

This writes
`artifacts/generated-results/support_saturation_cubic_annihilator_frontier.json`.
It closes further smooth-symbol quartic saturation searches, queues the six
singular squarefree cases by annihilator type, and places the
generically-étale/Keller gate before saturation for the double-line,
triple-line, and zero symbols.  It is a routing certificate, not a new
singular saturation computation.

The older degree-42 and plane cases require Singular; the plane primary
decomposition is the longer run.  Their exact scopes, especially the
remaining full characteristic-zero saturation and order-seven gaps, the
unproved generic all-order degree-42 statement, and the still-undefined
plane Case-1 conductor/residue module, are recorded in
[`extended-geometry/SUPPORT_SATURATION_COMPILER.md`](extended-geometry/SUPPORT_SATURATION_COMPILER.md).

## HC4 projective polar atlas

Regenerate and verify the low-degree projective-degree/Segre-signature
atlas with:

```bash
.venv/bin/python scripts/verify_hc4_projective_polar_atlas.py
```

Independently compute the graph-compactification and full-polar
multidegrees for the quadratic and cubic constant-Hessian calibrations
with Macaulay2:

```bash
M2 --script scripts/verify_projective_polar_calibrations.m2
```

The Python command writes
`artifacts/generated-results/hc4_projective_polar_atlas.json`.  The
formula, the graph-versus-polar distinction, the cotangent and
Meng--Yang controls, Wang's exclusion of all sixteen quadratic-gradient
affine-degree-two/three rows, the `HC4CQ1` exclusion of all 139
cubic-gradient rows, the 319 and 307 quartic-gradient numerical rows, the
rank-one/two/three leading-quintic determinant faces, the exact
rank-three cubic Schur gap, its squarefree Hessian-discriminant
obstruction and exact witness, the resulting potential-degree lower bound
five, the generic essential-rank top-gradient/Rees support sieve, its exact
intersection with the three atlas codimension columns, and the exact
nonexistence scope are documented in
[`HC4_PROJECTIVE_POLAR_GEOMETRY.md`](HC4_PROJECTIVE_POLAR_GEOMETRY.md).

Construct the universal 56-coefficient quintic top part, verify its
gradient/Hessian/Euler/Koszul and midpoint-collision identities, build the
generic essential Hessian-rank strata, and intersect their exact support
codimensions with the atlas using:

```bash
.venv/bin/python scripts/analyze_hc4_quintic_infinity_rees.py
```

Independently certify that the generic smooth rank-one/two/three top ideals
are equal-degree complete intersections of linear type, with the stated
pure-top projective degrees, using:

```bash
M2 --script scripts/verify_hc4_quintic_infinity_rees_strata.m2
```

The Python command writes
`artifacts/generated-results/hc4_quintic_infinity_rees_strata.json`.
The pure-top Segre vectors are degeneration calibrations, not completed
constant-Hessian Segre classes.  The exact restrictions on the actual
atlas are the support-codimension/Segre-vanishing filters; lower-layer
normal-cone multiplicities remain open at this stage.  The next checker
below closes the smooth rank-three vertex packet.
The recorded replay environment is the repository Python lock together with
Macaulay2 1.22 and its `Cremona` and `ReesAlgebra` packages over
\(\mathbb Q\).

Close the smooth essential rank-three, codimension-four packet by verifying
the \(\epsilon\)-flat length-\(256\) local complete intersection, the
socle bound \(\dim(Bs_3)\ge2\), and the resulting affine-degree bound
\(\delta\ge6\):

```bash
.venv/bin/python scripts/verify_hc4_rank3_vertex_colength.py
M2 --script scripts/verify_hc4_rank3_vertex_colength.m2
```

The Python command writes
`artifacts/generated-results/hc4_rank3_vertex_colength.json` and intersects
the theorem with the atlas, excluding the signatures
`(1,4,16,64,2)` and `(1,4,16,64,3)`.  The Macaulay2 command independently
checks the complete-intersection Hilbert function and exact Fermat and
deformed local calibrations.  The universal conclusion comes from the
flatness/socle proof in
[`HC4_PROJECTIVE_POLAR_GEOMETRY.md`](HC4_PROJECTIVE_POLAR_GEOMETRY.md), not
from extrapolating those representatives.

Refine the two codimension-three packets with the rank-two
constant-kernel/Schur calculation and the rank-three ordinary-singularity
incidence:

```bash
.venv/bin/python scripts/verify_hc4_codim3_gradient_strata.py
M2 --script scripts/verify_hc4_codim3_gradient_strata.m2
```

The Python command writes
`artifacts/generated-results/hc4_codim3_gradient_strata.json`.  It proves
that a nonzero rank-two kernel restriction \(h_4|_K\) synchronizes a
constant direction, closes the squarefree binary-Hessian branch through
`HC4CD5`, and forces \(\sigma_3=16\) on the nonsquarefree remainder.  It
also checks that the rank-three Schur cubic vanishes at every isolated
singular point where the top Hessian has rank two.  The Macaulay2 replay
checks the radical synchronization powers, transverse lengths \(64\to16\),
and a nodal Hessian calibration.  These are packet restrictions, not
unconditional deletions of codimension-three atlas rows.

Apply `PGS3` to the essential-rank-two singular binary-quintic packet with:

```bash
.venv/bin/python scripts/verify_hc4_binary_root_partition_segre.py
M2 --script scripts/verify_hc4_binary_root_partition_segre.m2
```

On the open stratum where a redundant active gradient has
\(X_0\)-order one with unit coefficient at every repeated root, a root of
multiplicity \(e\) contributes exactly \(e-1\) to \(\sigma_2\).  Hence a
binary quintic with \(q\) distinct roots has \(\sigma_2=5-q\); the generic
double-root packet retains only 51 and 50 atlas rows for affine degrees two
and three.  The Macaulay2 replay checks multiplicities \(2,3,4\).
The higher-\(X_0\)-torsion failure locus remains open, so this is not an
unconditional row deletion.

## All-dimensional projective-gradient Segre machinery

Verify the canonical
\((g_0,\ldots,g_n)\leftrightarrow(\sigma_1,\ldots,\sigma_n)\) transform,
the actual affine-gradient and full-polar constructors, the leading
integrability/Euler reconstruction, and regenerate the typed family registry
with:

```bash
.venv/bin/python scripts/verify_projective_gradient_segre_machinery.py
.venv/bin/python scripts/verify_projective_gradient_normal_slices.py
.venv/bin/python scripts/verify_projective_gradient_singular_slices.py
```

Independently compute the exact plane-cotangent and
quadratic-stabilization multidegrees, and replay representative
smooth-essential normal slices, with:

```bash
M2 --script scripts/verify_projective_gradient_segre_families.m2
M2 --script scripts/verify_projective_gradient_normal_slices.m2
M2 --script scripts/verify_projective_gradient_singular_slices.m2
```

The Python commands write
`artifacts/generated-results/projective_gradient_segre_registry.json` and
the smooth and singular normal-slice ledgers
`artifacts/generated-results/projective_gradient_normal_slices.json` and
`artifacts/generated-results/projective_gradient_singular_slices.json`.
Complete multidegree/Segre vectors, top-degree-only transport controls, and
explicit families with uncomputed vectors are distinct record types.  The
normal-slice artifact records the dimension-free complete-intersection
Hilbert series, filtered missing-generator bound, exact unit-penultimate
law, and the HC4 specializations `HC4PPG7` and `HC4PPG8`.  The canonical
singular ledger records the kernel-vertex/singularity join, the exact
truncated DVR-module formula, and a repeated-root binary quintic whose
lower quartics realize active lengths \(8,3,2\).  This proves that the
singular support alone does not determine its Segre multiplicity.  The
canonical scope and the resulting restrictions on cotangent, Schur, HN,
coefficient-scheme, and boundary-normalization consumers are documented in
[`PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md`](PROJECTIVE_GRADIENT_SEGRE_MACHINERY.md).

## Coefficient-space and Kuranishi calculations

Certify the exact full-box tangent ranks for `F_4,F_5,F_6`, their visible
seed ranks and source gauges, and one nonzero quadratic obstruction in each
degree with:

```bash
.venv/bin/python scripts/verify_all_degree_coefficient_tangents.py
```

Certify the complete characteristic-zero quartic quadratic Kuranishi rank and
the explicit reduced-family tangent ranks with:

```bash
.venv/bin/python scripts/verify_quartic_full_box_kuranishi.py
.venv/bin/python scripts/verify_generic_coefficient_family_tangents.py
```

Regenerate the modular first-order source filtration with:

```bash
.venv/bin/python scripts/research_filtered_source_tangent_profile.py \
  --prime 32003 \
  --json-output artifacts/generated-results/filtered_source_tangent_profiles_mod32003.json
```

Regenerate the modular quartic slices and the Singular/Macaulay2 input files
with:

```bash
.venv/bin/python scripts/research_quartic_coefficient_kuranishi.py \
  --prime 32003 --jet-order 8 \
  --json-output artifacts/generated-results/quartic_coefficient_kuranishi_mod32003.json \
  --singular-output artifacts/generated-results/quartic_coefficient_kuranishi_mod32003.sing \
  --macaulay2-output artifacts/generated-results/quartic_coefficient_kuranishi_mod32003.m2

.venv/bin/python scripts/research_quartic_generic_component.py \
  --prime 32003 --greedy-jet-order 6 \
  --json-output artifacts/generated-results/quartic_generic_component_mod32003.json \
  --singular-output artifacts/generated-results/quartic_generic_component_mod32003.sing \
  --singular-order3-output artifacts/generated-results/quartic_generic_component_order3_mod32003.sing \
  --macaulay2-output artifacts/generated-results/quartic_generic_component_mod32003.m2
```

The optional cubic-layer compilation takes several minutes.  The CAS files
are research inputs.  Their full primary decompositions have not completed
within the available memory and are not certificate artifacts.
The theorem/computation boundary and the all-degree formal-versus-algebraic
statement are documented in
[`extended-geometry/JELONEK_COEFFICIENT_COMPONENTS.md`](extended-geometry/JELONEK_COEFFICIENT_COMPONENTS.md).

## Free-discriminant and Saito-matrix experiment

Verify the four first marked-root discriminants, the three full-target
Saito bases, the fixed-\(P\) quadratic-gauge Saito basis, all regular
marked-root lifts, and the type-\(A_3\) reflection control with:

```bash
.venv/bin/python scripts/verify_free_discriminant_saito.py
```

The weighted and cancellation branch surfaces, and the full
quadratic-gauge branch and ledger divisors, have a separate exact Singular
nonfreeness certificate:

```bash
Singular -q scripts/verify_free_discriminant_saito_nonfree.sing
```

The Singular command takes about one minute.  It verifies codimension two
and minimal-resolution length three for both Jacobian ideals; it is not a
bounded search.  The formulas, the corrected Saito--incidence proposition,
and the external-candidate gates are in
[`cancellation/FREE_DISCRIMINANT_SAITO_EXPERIMENT.md`](cancellation/FREE_DISCRIMINANT_SAITO_EXPERIMENT.md).

## Bidegree-\((3,3)\) Rodrigues survivor and sparse census

Verify the full-rank all-order pure-moment survivor, its beta factorization,
the Rodrigues identity, and the arbitrary-multiplier SIC cutoff with:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rodrigues_survivor.py
```

Reproduce its exact normalized null-quadratic local certificate and the
five-variable Singular slice with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_null_quadratic_s6.py \
  --orders 2,3,4,5,6,7,8,9,10,11 --skip-solver \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_null_quadratic_s6_local.json
```

The optional unsaturated full ten-variable modular solve is an experiment,
not part of the local theorem.  It can be run by omitting
`--skip-solver` and adding `--prime 43 --timeout 600`.

Reproduce the exact characteristic-zero anti-Weyl exclusion with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_anti_weyl.py \
  --prime 0 --through 14 --backend msolve \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_anti_weyl_normalized_msolve14_char0.json
```

Reproduce the exact isolated rank-two nine-moment component and its
nonlifting signs with:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_finite_prefix.py
```

This exact calculation builds the primitive anti-Weyl
moments, performs the rational Krawczyk inclusion in the radius-\(10^{-10}\)
box, proves coefficient rank two, certifies tangent rank eight on the smooth
rank-two chart, and bounds the primitive moments with signs
\(\mu_{10}>0,\mu_{12}<0,\mu_{14}>0\).  It writes
`artifacts/generated-results/two_pair_sic_bidegree33_rank_two_finite_prefix.json`.
It also derives the five-variable anti-Weyl square-invariant quotient and,
using characteristic-zero `msolve`, proves the unit ideal for the corrected
rank-two system on this chart.

Reproduce the two exact exclusions on the generic rank-two Hurwitz chart
with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_hurwitz.py \
  --characteristic-zero --backend msolve --minor 01 \
  --lambda-value 0 --orders 2,3,4,5,6,7,8 \
  --timeout 60 --memory-gb 3 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_lambda0_char0.json

.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_rank_two_hurwitz.py \
  --characteristic-zero --backend msolve --minor 01 \
  --orders 2,3,4,5,6,7,8 \
  --mu2-pivot-boundary-reduced secondary \
  --timeout 120 --memory-gb 3 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_rank_two_hurwitz_secondary_boundary_char0.json
```

The first command proves the unit ideal on the full \(\lambda=0\) fibre
after localizing the quadratic discriminant and the \(01\) channel minor.
The second imposes both successive \(\mu_2\)-pivot equations and proves
that branch unit.  The \(521\)-bit modular construction is an exact
integer recovery: the script checks the coefficient bound
\((3m)!\,52^m\) before invoking `msolve` over \(\mathbb Q\).

For a numerical-algebraic complexity estimate only, run:

```bash
julia --project=. \
  scripts/research_two_pair_sic_bidegree33_rank_two_homotopy.jl
```

It reports mixed volume \(74\,144\) for the unreduced square system
\(\mu_2,\ldots,\mu_8\).  This is not an exclusion, a solution count, or
a characteristic-zero certificate.

The complete six-entry coefficient-torus census is split into four
independent exact characteristic-zero shards:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 6 --through 12 --start 0 --limit 1897 \
  --output /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard0.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 6 --through 12 --start 1897 --limit 1897 \
  --output /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard1.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 6 --through 12 --start 3794 --limit 1897 \
  --output /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard2.json
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 6 --through 12 --start 5691 --limit 1897 \
  --output /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard3.json

.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 6 \
  --combine-input /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard0.json \
  --combine-input /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard1.json \
  --combine-input /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard2.json \
  --combine-input /tmp/two_pair_sic_bidegree33_sparse_six_support_screen_shard3.json \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_sparse_six_support_screen.json
```

For support size seven, replace the four starts by
`0,2800,5600,8400`, use `--limit 2800 --support-size 7`, name the
temporary files `two_pair_sic_bidegree33_sparse_support7_shard0.json`
through `shard3.json`, and combine them into
`artifacts/generated-results/two_pair_sic_bidegree33_sparse_support7_screen.json`.
The size-six census has exactly two normalized survivors, both on the
Rodrigues orbit; the size-seven census excludes all \(11{,}200\) mixed
coefficient tori.  Boundaries are covered by the separately verified
smaller-support results.

Certify that the two size-six nonunit systems each have exactly one
complex point, rather than merely one real box, with:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_survivor_rur.py
```

For support size eight, use four contiguous shards of length 3,195 with
starts `0,3195,6390,9585`, `--support-size 8 --through 12`, and combine
them into
`artifacts/generated-results/two_pair_sic_bidegree33_sparse_support8_screen.json`.
The sole timeout in that census is reproduced exactly with:

```bash
.venv/bin/python \
  scripts/research_two_pair_sic_bidegree33_sparse_six_counterexample.py \
  --support-size 8 --through 14 --start 8384 --limit 1 \
  --timeout 600 --threads 4 \
  --output \
  artifacts/generated-results/two_pair_sic_bidegree33_sparse_support8_parity_msolve14_char0.json
```

The recorded run takes about eight minutes.  Validate that unit
certificate, rerun full complex RURs for the fourteen nonunit systems,
and check their explicit one-sided normal forms with:

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_support8.py
```

Add `--rerun-parity` only when intentionally refreshing the pinned
long-running parity artifact.  The resulting theorem is a complete
coefficient-torus classification through support size eight: any actual
bidegree-\((3,3)\) SIC counterexample has at least nine nonzero standard
monomial coefficients.

The first complete nine-entry class, the sixteen \(3\times3\) coefficient
rectangles, is replayed by

```bash
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_rectangle9.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_two_row_fringe9.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_cross_two9.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_three_line9.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_three_line4329.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_sparse_full_line43119.py
.venv/bin/python \
  scripts/verify_two_pair_sic_bidegree33_rank_two_parity_channels.py
```

The checker uses six transpose/reversal orbits.  It forms exact finite
schemes over \(\mathbb Q\) through moment fourteen and proves that all 54
\(2\times2\)-minor localizations are unit ideals.  Thus every reduced
dense component has rank one; boundaries are covered by the size-eight
theorem.  This does not run or claim the full \(11{,}420\)-support
size-nine census.

The second checker closes 96 complete-two-row/column fringe supports
in 24 exact symmetry orbits.  It also verifies the full discrete census
of 11,420 mixed size-nine supports in 2,924 transpose/reversal orbits,
of which 30 are closed by the rectangle and fringe theorems.  The third
checker closes all 576 complete-row/complete-column cross-plus-two
supports in 156 exact symmetry orbits through \(\mu_{10}\); six of these
supports are outside the mixed census.  The fourth sparse checker closes
all 480 regular three-line supports in 120 exact symmetry orbits.  Its 114
unit ideals and six unique rational fixed-flag rank-two points close the
regular \(3+3+3\) class.  The fifth sparse checker closes all 1,148
\(4+3+2\) three-line supports in 287 exact symmetry orbits.  This finishes
every one of the 1,740 nine-entry supports with an empty row or column.
Overall, 2,310 mixed supports are closed in 591 orbits, with 9,110
supports in 2,333 orbits remaining.  The sixth sparse checker closes all
1,244 full-line \(4+3+1+1\) supports in 311 exact symmetry orbits through
\(\mu_{10}\).  Overall, 3,554 mixed supports are closed in 902 orbits,
with 7,866 supports in 2,022 orbits remaining.  The last checker uses the
reversal-centralizer orbit cover to classify the complete exact-rank-two
reversal-parity factor family.  Exact characteristic-zero msolve
calculations close both projective semistable charts through \(\mu_6\),
and a Singular minimal-prime decomposition finds exactly two components
on the invariant-zero boundary.  Fixed-flag factorization then certifies
their all-order recurrence, initial vanishing, nonzero degree-two mixed
values, and the mixed cutoff \(2m>e\).
