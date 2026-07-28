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
is the non-finite-generation proof.  The same note proves that the two
leading divisors in a tangent-normalized factorization slice are disjoint.
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
.venv/bin/python scripts/verify_padic_inverse_branches.py
.venv/bin/python scripts/verify_composite_degree_twelve.py
.venv/bin/python scripts/verify_degree_twelve_wreath_elimination.py
.venv/bin/python scripts/verify_all_degree_rational_fibers.py
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
.venv/bin/python scripts/verify_low_rank_multiplicity_boundaries.py
.venv/bin/python scripts/verify_real_fiber_spectrum.py
.venv/bin/python scripts/verify_adelic_fiber_engineering.py
.venv/bin/python scripts/verify_local_global_keller_fibers.py
.venv/bin/python scripts/verify_a5_grunwald_keller_fiber.py
.venv/bin/python scripts/verify_hasse_keller_fiber.py
.venv/bin/python scripts/verify_infinite_hasse_keller_fibers.py
.venv/bin/python scripts/verify_multiplicative_hasse_artifact.py
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

The separately authored Lean certificate is optional because it downloads a
pinned toolchain:

```bash
make verify-lean-foundational
```

GitHub Actions runs this target in the required `formal-lean` job using the
pinned upstream commit and Lean action.  The `papers` job compiles the
finalized and active manuscripts listed in `papers/README.md`, while parked
manuscripts remain available for direct local builds.  The
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
```

For the homogeneous tensor, all seven squarefree strata have saturated
cotangent presentation and a length-six `Ext_A^2(T,A)` support defect with
Hilbert function `3+3t` (three-dimensional top and zero `m^2` action);
double and triple lines instead have a one-dimensional support defect.
The zero homogeneous tensor passes both module tests but is nowhere
generically étale.  One explicit order-four kernel tensor makes the support
defect finite of length six in all ten strata, while cotangent saturation
still passes.  This is an exact leading-model computation.  It neither
proves lift-independence nor constructs a normal lift with a Keller open.
The second command works over `Q[t,x,y,z]`.  On each of the seven
squarefree lines `phi_h+t*psi_4`, it verifies uniform cotangent saturation,
no parameter torsion in relative `Ext_A^2(T,A)`, radical support equal to
the collision axis, multiplicity six, and equality of the relative
presentation with the scalar extension of its central specialization.
This proves constancy on those lines, not on the full 24-parameter
order-four space.

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

The companion one-profile Hopf classification is checked directly by

```bash
.venv/bin/python scripts/verify_hopf_lift_classification.py
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
.venv/bin/python scripts/audit_keller_near_invariant_backtrace.py
python3 scripts/audit_keller_observable_quotients.py
python3 scripts/audit_keller_provenance_compression.py
```

The first command gives characteristic-zero rank and Lie-image certificates
that the degree-at-most-five invariant space of the stored 21-variable map is
`Q[X_20]_{<=5}`. It also records the near-invariant
`Q=X_18*X_20-X_6*X_8` and its one-term pullback defect. The second uses two
exact torus gradings to exclude both sextic correction channels, while leaving
pure homogeneous sextic invariants open. The third reconstructs the frozen
17-step circuit and identifies `Q=c_4*s-v_3*v_5` as a determinantal
shared-factor gate residual whose stable-source restriction is `x^2*y*z`.
The fourth proves that any rational semiconjugate quotient carrying either
`X_0` or the restricted
quadratic observable has dimension at least 13; its longer rank plateau is
printed as experiment only, while a stacked rank-20 certificate excludes a
common constant translation direction behind that plateau. The fifth verifies that the normalized
three-variable canonical contraction fails at its first pure moment, checks
the full twenty-coordinate inverse-recurrence dependency closure after the
known identity slice, and reports the finite stored-circuit census.  Only
the degree-at-most-five invariant statement is a nonlinear quotient-class
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
the generalized repeated-pole beta recurrence with its order-one
factorial/double-factorial evaluation; it also formalizes the scalar chart
identity, the monomial contraction/coefficient-extraction equality, the
product-polynomial evaluation at natural numbers, and both normalized
displayed binomial-sum identities.  It does not formalize their linear
assembly for \(F\) or the derivation of the chart constant terms, so the
repository does not label SIC2C4 itself as formally verified.

The exact local geometry of this displayed \(F\) is checked by

```bash
.venv/bin/python scripts/verify_two_pair_counterexample_local_moduli.py
```

The minimum-degree separating invariant and the low-degree invariant-ring
calculation are checked independently by

```bash
.venv/bin/python scripts/verify_two_pair_counterexample_missing_invariant.py
```

This constructs the five Casimir projectors for
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

This computes the contraction-preserving stabilizer and orbit, uses the
all-order Hopf coefficient identity to obtain the thirteen-dimensional
all-moment tangent space, records the seven independent quadratic lifting
obstructions, proves that their quotient radical is a five-plane, constructs
a polynomial cubic lift for every direction on that plane, and verifies the
defect-preserving family
\[
F_{a,b}=\frac{aR+bZ}{2}
\left(2W(aR+bZ)^2-2abR^3-b^2R^2Z\right).
\]
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
The same replay checks all opposite odd-height monomial pairs and all
opposite even-height pairs of phase \(s\geq3\) in the displayed degree
range against the written all-degree formulas.  Finally it verifies the
degree-five phase-one elimination for
\(RF+aZT^4+bWT^4\): the normalized moments \(2,3,4\), their exact
lexicographic remainder, and the nonzero resultant
\(-418538718730248905250\).  Together these calculations exclude every
opposite monomial pair in degree five.
See
[`TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md`](extended-geometry/TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).
This proves `not MN_d` for every `d>=4`, `not SIC(2)`, and, with the known
one-pair theorem, the exact minimum failing pair dimension two.  The finite
replay is not being used as the all-order proof.

The coefficient-rank frontier inside bidegree \((4,4)\) is replayed by

```bash
python3 scripts/verify_two_pair_sic_bidegree44_rank_frontier.py
python3 scripts/verify_two_pair_sic_bidegree44_rank_two_invariants.py
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
necessary Hilbert test through degree \(100\), but their zero fiber remains
open. The exact rank of the extra truncated-moment point may be one or two,
and no all-order counterexample is claimed.
On the rank-one boundary, the same checker proves that moments one through
six have exact Jacobian rank six and computes their nonnegative Hilbert
numerator, of coefficient sum \(50\). This reduces the boundary gate to
the finitely many exceptional squarefree quartic cross-ratios and a
uniform finite cutoff on the low-root orbits; it does not claim that this
remaining fiber calculation is complete.

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

Over both \(\mathbb F_{47}\) and \(\mathbb F_{101}\), with the six base
variables treated as rational-function parameters, the checker obtains a
three-element basis and quotient length six for
\((\mu_4,\mu_5)\).  Both \(\mu_6,\mu_7\) reduce to six fiber coordinates.
The denominators reconstruct as products of
\(L=s_1t_0-t_1\) and
\(Q=s_1^2-s_2-(13/3)t_0^2\).  This identifies the principal open and two
degeneracy divisors for the next elimination; it remains modular
reconstruction evidence rather than a characteristic-zero certificate.

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
The collision route retains witness sizes 20/40/42 and rank 37.  The
independent Dvorsky--Long formulas lower the overall certified SIC and
unrestricted GVC entries to 5/5, but do not change the ordinary-Laplacian
40 or homogeneous HN 42 entries.  These are witness-ledger values, not
literature-wide minimality claims.
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
Singular -q scripts/verify_a4_corrected_boundary_genus.sing
```

The normalized-boundary command verifies the determinant-one ambient
completion and the exact obstruction to the resulting automorphic assembly.
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
as an affine plane.

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
.venv/bin/python plane-jc/cas/audit_f2_75_125_boundary_handoff.py
.venv/bin/python plane-jc/cas/test_f2_75_125_frontend.py
```

The first replay proves that coarse Newton vertices, geometric degree, and
nonproperness data cannot upper-bound affine-normalized support: a
Zariski-open family of triangular automorphisms has fixed coarse geometry
and support lower bound `d-2`.  It then verifies the Kummer chain-rule gate
and the live `(75,125)` F2 terminal characters `P={1,4}`,
`Q={0,1,3}` modulo five, which block constant-Jacobian descent.  The second
replay corrects the Laurent chart to `[t,z]=-z` and classifies the exact B0
degree/halfspace envelope of all 35 zero layers (`39` through `5`): 665 band
pairs, 978 jet-reduced linear parameters (973 after normalization), and
5,348 structurally active character-split Keller rows.  This does not
exclude F2.  The same replay proves that the common-power top root is not an
arbitrary degree-18 polynomial: it has the exact two-parameter form
`H(t)=(1+u+...+u^4)^2*R(u^5)`, `u=1+t`, with `R` quadratic and
`R(1)=1/25`.  The triangular recurrence forces layers `39,...,36` to
continue that common root.  Layer 35 is the first genuine branch: a
one-dimensional `lambda*C0^2` mode survives, so the upper five layers do not
give a unit ideal.  The boundary-handoff replay then gives four exhaustive
contact partitions of 25.  It proves that their multiplicities do not
determine branch scales or finite-normalization rows; even the unsupported
strongest contact-to-row surrogate survives the degree-26 packet budget.
This stops the degree-specific F2 route and replaces the failed coarse bridge
by a finite character-resolved B0 system and a precise B1 polygon/`gamma`
obligation; see
[`plane-jc/AFFINE_SUPPORT_NEWTON_BRIDGE.md`](plane-jc/AFFINE_SUPPORT_NEWTON_BRIDGE.md).
The final replay checks the forced chain and terminal normalization
independently.

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
now a possible witness, not the statement being assumed.  See
[`plane-jc/JC2_GLOBAL_COX_PACKET_ATTACK.md`](plane-jc/JC2_GLOBAL_COX_PACKET_ATTACK.md).
The normalization and conductor formulas are written algebraic proofs.  The
Python checker replays their cusp factorization, determinantal identities,
and monomial conductor quotient; the independent Singular command computes
the cusp normalization and conductor and verifies normality of the
determinantal overring.

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

## Direct Schur-descent audit for `HC_4`

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
higher-degree critical equations are not tested.

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
