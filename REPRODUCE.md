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
map.
It also checks the normalized `(2,3)` factorization slice: the unimodular
boundary lattice, class `L^5-L^3`, direct counts `q^5-q^3` for four small
prime fields, and generic degree ten.  A separate two-chart certificate
checks the Euclidean quadratic norm, its affine-modification presentation
over `A^2 x SL_2`, the complementary-chart transition, and the integral
residue coefficients `1`, `2`, and the nonzero mod-two boundary used to prove
that integral cohomology is `Z` in degrees zero and three only.  The scripts
check the algebraic, arithmetic, and Gysin inputs; the written audit supplies
the localization sequences and homotopy argument.
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
`verify-foundations` adds the weighted construction and its clean-room checker.
It also runs the all-degree rational-fiber checker, whose symbolic odd/even
identities prove uniform admissibility and whose exact degrees `3,...,100`
remain as a regression:

```bash
.venv/bin/python scripts/verify_all_degree_rational_fibers.py
.venv/bin/python scripts/verify_real_fiber_spectrum.py
.venv/bin/python scripts/verify_adelic_fiber_engineering.py
.venv/bin/python scripts/verify_hasse_keller_fiber.py
.venv/bin/python scripts/verify_stratified_adelic_engineering.py
```

The Hasse-fiber command expands an explicit degree-eight weighted map, checks
its determinant `-38`, proves that its complete target fiber has no rational
point, and audits roots over `R` and every `Q_p` through the elementary
quadratic-residue covering and the two exceptional Hensel lifts.

The last command audits the constructive CRT/weak-approximation lift and an
explicit nonsurjective type-`(3,2)` quintic seed with trivial Hessian
stabilizer and complete fibers of all three quintic signatures, each with
cycle types `(5)` at `7` and `(2,2,1)` at `11`. The preceding adelic command
audits an explicit totally imaginary quartic complete fiber that is inert at
`7` and has unramified splitting type `(2,1,1)` at `11`.

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
The decisive exact cubic-field check is

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --exact-cubic --print-radical-basis --seventh-line
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --seventh-component-elimination
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
The characteristic-zero lift of that final unit identity remains open.

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

The separately authored Lean certificate is optional because it downloads a
pinned toolchain:

```bash
make verify-lean-foundational
```

GitHub Actions runs this target in the required `formal-lean` job using the
pinned upstream commit and Lean action.  The `papers` job compiles every
standalone paper discovered at `papers/*/main.tex`, and
`macaulay2-independent-check` runs the pinned
Macaulay2 comparison.  Together with the four Python matrix jobs, these are
the complete CI verification pipeline.  The final `verification-complete`
job is the single aggregation check intended for GitHub branch protection.

## Cancellation programme

```bash
make verify-master
```

This target includes the exact quadratic-gauge/cancellation intersection
regression.  To run its symbolic `N=4,5,6,7` discriminant and all-factorization
checks directly:

```bash
.venv/bin/python scripts/verify_quadratic_cancellation_intersection.py
```

The minimal-boundary gateway and classification program has a separate fast
cubic certificate:

```bash
make verify-minimal-boundary
```

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
is one.
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
for `r=5,6,7`.  At `r=7` it verifies a squarefree degree-42 branch polynomial,
excludes `c=0` and `z=infinity`, and reconstructs a unique finite `z`.  This
shows that the branch-at-infinity mechanism survives for the next fixed
column.  The first command is the quick limiting-system audit.  The second
also constructs the full bidegree-`(42,126)` endpoint eliminant, identifies
its complete top Newton edge, and proves eventual nonvanishing for `r=7`.
Neither command claims an effective `r=7` threshold or a uniform theorem in
`r`:

```bash
.venv/bin/python scripts/verify_contact_resultant_fixed_r_branch_schema.py
.venv/bin/python scripts/verify_contact_resultant_r7_asymptotic.py
```

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
the finite H2 factor formal over the corrected H1 graph.  The
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

Christopher D. Long's direct Gaussian-moment and `(xz)` identities, together
with the exact normalization of the foundational map used in his BCW
discussion, have a dedicated target:

```bash
make verify-external-consequences
```

The Gaussian and `(xz)` scripts use only the Python standard library.  Their
bounded exact regressions are distinguished from Long's written all-exponent
proofs.  A companion symbolic checker proves the `SU(2)=S^3` Haar density in
Hopf coordinates, completing the local integration proof.  The same target
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

The sharp remaining Gaussian dimension and the cross-conjecture minimum
ledger have their own fast exact target:

```bash
make verify-counterexample-scoreboard
```

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
target also recomputes three unit Groebner bases for the direct Long-style
collapse and writes the dimension/rank/index/degree scoreboard.  It does not
claim to settle GMC(2).

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

The all-order nonvanishing proof is written in
[`IMAGE_VANISHING_COUNTEREXAMPLES.md`](extended-geometry/IMAGE_VANISHING_COUNTEREXAMPLES.md);
the generators check the finite artifacts and change-of-variable identities;
the dependency-free audit re-expands the 40-variable witness from scratch.
A local proof of the
fixed-dimensional DVEZ/Zhao implication, including Gaussian contraction, the
countable-union step, and formal inversion, completes the nonexplicit route to
`not GMC(42)`; `not GMC(158)` remains the exact conservative Long-route bound.
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

## Canonical degreewise paper

The canonical statement and proof are in the standalone paper; the four-input
verification companion is
[`DEGREEWISE_MULTIPLICITY_AUDIT.md`](DEGREEWISE_MULTIPLICITY_AUDIT.md).  Build
the paper with:

```bash
cd papers/marked-root-multiplicity
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

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

## Plane degree-frontier audit

The small deterministic regression of the published candidate tables is:

```bash
python3 plane-jc/cas/frontier_125_150.py
```

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
