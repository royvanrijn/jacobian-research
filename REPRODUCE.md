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

This runs the construction, parameter arithmetic, boundary, monodromy, and
current-ansatz rigidity regressions.  It includes the endpoint-moment
reduction of the cancellation contact resultant: the general triangular
identity is checked exactly on a bounded grid, while the complete
`r=1,2,3,4` columns are proved uniformly in `m`.  It also checks the
irreducibility transfer proving every `1<=m<=6` column uniformly in `r` and
the eventual `r`-tail for each fixed `m`:

```bash
.venv/bin/python scripts/verify_contact_resultant_irreducible_ranges.py
```

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

The first not-yet-complete fixed-`r` column has an exact bounded-degree
reduction: the following Singular-backed checker constructs the
quintic--sextic endpoint equations and verifies that their residual eliminant
has degree 29 in `y` and degree 90 in `m`.

```bash
.venv/bin/python scripts/verify_contact_resultant_r6_reduction.py
```

The branch-at-infinity replay then proves eventual nonvanishing in that
column.  It checks the complete Newton edge after `y=1+c/m`, the squarefree
degree-29 edge polynomial, and the linear reconstruction of the limiting
`z`.  Lindemann--Weierstrass separates algebraic `z` from `exp(c)`.  The
argument does not yet provide an explicit threshold in `m`.

```bash
.venv/bin/python scripts/verify_contact_resultant_r6_asymptotic.py
```

For an effective certificate on the limiting edge itself, the following
atlas computes the first two `y` terms and first two `z` terms for all 29
branches (compressed to 15 complex-conjugation classes).  It uses 29 disjoint rational
Rouche disks and rational exponential enclosures to prove the strict modulus
gap `|z_0|!=|exp(c)|` branch by branch.  It still does not extract a positive
tail threshold in `t=1/m`.

```bash
.venv/bin/python scripts/explore_contact_resultant_r6_branch_atlas.py
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
.venv/bin/python scripts/verify_multicluster_ll_comparison.py
.venv/bin/python scripts/verify_rerooting_groupoid_boundary.py
.venv/bin/python scripts/verify_coarse_affine_mark_descent.py
.venv/bin/python scripts/verify_restricted_ll_degree.py
.venv/bin/python scripts/verify_caustic_maxwell_boundary.py
```

These checks support the generic affine-mark faithfulness theorem: the coarse
fiber is the exact rerooting orbit and every nontrivial rerooting moves the
unique unramified affine sheet into the reconstruction boundary.  The
selected root extends on the marked admissible-cover stack, and the
normalized-Stein, completed-chart, and conductor comparisons are complete at
arbitrary simultaneous collisions.  Coarse affine-mark descent is also
complete: the marked invariant ring is the universal monic-root incidence,
and the total-collision fiber `k[T]/(T^mu)` has one geometric point.  The
specialized restricted-LL
degree and caustic/Maxwell boundary-class calculations have no recorded
external review.

The companion affine-stratum audit verifies that the root-one component is
regular and that a nontrivial rerooting sends it to an extra-root boundary
component.  The multicluster audit checks distinct tangent lines, all pairwise
intersection numbers, the conductor exponent
`e_i(sum_j e_j-1)`, and regularity of the full marked-root incidence at
collisions.  The rerooting-groupoid audit separately checks the quotient degree `N-2`,
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
```

It derives the irreducible sextic equation of the all-double component,
parametrizes the all-triple curve, and verifies that their four all-six
collision points have scheme-theoretic intersection length two.  It now also
transports both degree-six vertical Ritt hypersurfaces into optimal moment
coordinates, verifies the second displayed sextic, and proves that the
`2 o 3` Ritt surface is exactly the all-double exceptional component.  The
underlying Hessian-incidence and Ritt-intersection calculation is replayed by

```bash
.venv/bin/python scripts/verify_hessian_ritt_degree_six.py
.venv/bin/python scripts/verify_degree_six_ritt_atlas.py
```

The second checker refines the Hessian-incidence result on the normalized
seed chart: it computes the `2^3` and `3^2` omitted-value intersections, the
four doubled type-`(6)` collision points, factored affine-sheet boundary cuts,
and clean rational witnesses for all open pieces.

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

The exact 90 MB certificate archive and extracted replay source are pinned
under `plane-jc/external/zenodo-21479814/`.  Attachment hashes, environment
versions, the full replay command, and the independent hard-certificate
command are in
[`plane-jc/cas/README.md`](plane-jc/cas/README.md).  The independent checker
does not import the primary CAS or generation modules.
