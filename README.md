# Jacobian counterexample research workspace

This repository independently checks and studies the polynomial map announced
on 19--20 July 2026:

\[
F(x,y,z)=\left(
(1+xy)^3z+y^2(1+xy)(4+3xy),
y+3x(1+xy)^2z+3xy^2(4+3xy),
2x-3x^2y-x^3z
\right).
\]

It also implements a chart-first computational attack on the still-open plane
Jacobian conjecture. The 3D verification and the plane search have different
logical status:

- the displayed 3D map has a finite, exact counterexample certificate;
- the 2D work has found obstructions in explicitly bounded support families,
  but no plane counterexample and no general 2D impossibility theorem;
- modular solutions and solver timeouts are never treated as characteristic-
  zero evidence.

Start with [the verification, discovery, and analysis
checklist](notes/FIRST_CHECKLIST.md). The quickest supporting orientation
documents are [the facts ledger](notes/FACTS.md),
[the implementation ledger](notes/IMPLEMENTATION_STATUS.md), and
[the clustered obstruction report](results/OBSTRUCTION_CLUSTERS.md).

## Exact 3D certificate

Exact symbolic expansion gives

\[
\det DF=-2.
\]

Exact rational substitution gives three distinct points in one fiber:

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).
\]

Thus the map is everywhere etale over `C`, is not injective, and is not a
polynomial automorphism. Padding with identity coordinates gives examples in
every dimension at least three. The coordinate degrees are `(7,6,4)`.

Run the independent finite check:

```bash
python3 scripts/verify_counterexample_independent.py
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python scripts/verify_counterexample.py
```

The first verifier is dependency-free and implements its own exact sparse
polynomial arithmetic; the second uses SymPy.

## What this project has added

The repository goes substantially beyond checking the determinant and the
three-point collision.

### Construction and inverse geometry

- The map is rediscovered from the rational chart
  `(t,r,c)=(y+1/x,2/x,R(x,y,z))`; its two Jacobian factors cancel exactly.
- The inverse problem is reduced to the primitive cubic

  \[
  cT^3-2T^2+bT-2a=0,
  \]

  with rational reconstruction of `(x,y,z)` from a root. The generic
  function-field degree is exactly three and the extension is non-normal, with
  `S_3` monodromy/Galois closure.
- The exact image, fiber stratification, and nonproperness set are computed:

  \[
  F(\mathbb C^3)=\mathbb C^3\setminus\Gamma,
  \qquad S_F=V(Q),
  \]

  where

  \[
  Q=27a^2c^2-18abc+16a+b^3c-b^2,
  \quad
  \Gamma=V(3bc-4,12a-b^2).
  \]

  Fibers have respectively `3`, `1`, and `0` affine points on `Q!=0`, on
  `Q=0` away from `Gamma`, and on `Gamma`. The curve `Gamma` is exactly the
  singular locus of `V(Q)`.
- The boundary chart identifies one dicritical divisor, gives the normalization
  of `V(Q)`, and explains how length-two or length-three fibers escape to
  infinity without finite ramification.
- A general weighted seed construction is implemented. Within that family the
  quadratic seed, and hence the cubic-sheet mechanism, is minimal.

See [CONSTRUCTION.md](notes/CONSTRUCTION.md) and
[IMAGE_AND_NONPROPERNESS.md](notes/IMAGE_AND_NONPROPERNESS.md).

### Arithmetic, dynamics, and numerical benchmarks

- Exact finite-field fiber distributions are proved in every characteristic.
  For characteristic greater than three the only fiber sizes are `0,1,3`, and
  the image density tends to `2/3`. Characteristics two and three are handled
  separately, with complete factorization and collision refinements.
- The inverse-Jacobian polynomial frame consists of three commuting vector
  fields. Every nonzero constant linear combination is incomplete. A complete
  weighted Euler combination is identified, as are the affine target
  symmetries and degree-at-most-two tangency constraints along `Gamma`.
- The three-point fiber is converted into a degree-14 polynomial with exactly
  three nondegenerate global minima and an explicit Palais--Smale sequence at
  infinity.
- The first weighted compactification of its gradient flow is computed. It has
  a nonhyperbolic equilibrium line; a formal center calculation shows why the
  Palais--Smale curve is not itself a gradient orbit. This analysis is partial,
  not a basin-boundary theorem.
- An adversarial homotopy benchmark has two exact solution paths escaping to
  one projective point with multiplicity two while a third path remains finite.
  It records the expected `delta^-3` conditioning growth and observed path-loss
  threshold.

See [FINITE_FIELD_VALUE_DISTRIBUTION.md](notes/FINITE_FIELD_VALUE_DISTRIBUTION.md),
[COMMUTING_FLOWS.md](notes/COMMUTING_FLOWS.md),
[GRADIENT_INFINITY.md](notes/GRADIENT_INFINITY.md), and
[NONPROPER_FIBER_BENCHMARK.md](notes/NONPROPER_FIBER_BENCHMARK.md).

### Explicit normal-form counterexamples

The standard high-dimensional reductions have been executed rather than merely
cited:

- a 95-dimensional map `I+H`, with `H` cubic homogeneous, 148 nonzero cubic
  terms, determinant one, nilpotent `JH`, and a stored rational collision;
- a 510-dimensional Druzkowski map `X-(AX)^{*3}`, with `rank(A)=95`, exact
  pairing matrices, determinant one, and a transported rational collision.

The sparse JSON artifacts are generated by the construction scripts and
independently checked by separate verifier scripts. See
[CUBIC_HOMOGENEOUS_REDUCTION.md](notes/CUBIC_HOMOGENEOUS_REDUCTION.md) and
[CUBIC_LINEAR_REDUCTION.md](notes/CUBIC_LINEAR_REDUCTION.md).

### Audited consequences

Using published implications and contraposition, the verified 3D example also
implies failure of the Mathieu conjecture for `SU(3)` and existential failures
of the all-dimensional Gaussian Moments, quartic Vanishing, and Image
conjectures. Except for `SU(3)`, these implications do not identify the least
dimension or construct a compact explicit witness. The claim-by-claim audit is
in [DIRECT_CONSEQUENCES.md](notes/DIRECT_CONSEQUENCES.md).

## The plane search problem

The central computational formulation is to find rational charts

\[
\mathbb C^2\dashrightarrow\mathbb C^2\dashrightarrow\mathbb C^2
\]

whose Jacobian factors are reciprocal, whose composed coordinates have exact
pole cancellation, whose composition extends polynomially across the pole
divisor, and whose distinct boundary branches acquire the same finite image.

The search is organized as:

1. enumerate short Cremona/toric/blow-up charts and canonicalize their
   valuation signatures;
2. solve pole cancellation with exact sparse linear algebra;
3. impose the reciprocal-Jacobian equation and collision normalization;
4. use modular F4/F4SAT only as a solving and component-removal technique;
5. apply homotopy only to plausible surviving zero-dimensional systems;
6. Hensel-lift and rationally reconstruct modular survivors;
7. verify every final identity over characteristic zero;
8. reject birational, Galois, too-small field-degree, or incorrectly merged
   boundary constructions.

No candidate has reached steps 5--7. Homotopy and Hensel reconstruction are
installed and validated on planted/nonempty controls, but current plane
families are either exactly empty, birational controls, or unresolved due to
scaling.

See [METHODS.md](notes/METHODS.md),
[IMPLEMENTATION_STATUS.md](notes/IMPLEMENTATION_STATUS.md), and the executable
stage scripts listed below.

## What has been tried in 2D

The following table is deliberately precise about finite scope.

| Stage | Executed scope | Result | What it does not prove |
|---|---|---|---|
| 3D chart rediscovery | Sparse rational-chart ansatz | Recovers the announced map and two deformations exactly | Nothing about dimension two |
| Degree 5--20 pilot | 96 named sparse edge systems per pass | Unnormalized: 64 exact contradictions and 32 birational controls; collision-normalized: all 96 inconsistent | Not all polynomial maps of degrees 5--20 |
| One-divisor triangular chart | `deg p<=4`, `0<=i<=6`, `-8<=j<=8`, plus sparse/bound increases | Complete exact pole kernels; dimensions `39--63`; first raw forbidden pole `x^-1`; no tested bad prime | Not the full 119+119 nonlinear ideal or all one-divisor supports |
| Stage A supports | Five named 10--30-term collision-normalized families and degree controls | Exact characteristic-zero reduced bases `[1]`; saturated agreement at three large primes | Only those supports |
| Two-divisor toric charts | Three unimodular charts, 10- and 15-term cones | Six exact characteristic-zero bases `[1]` | Only those charts/cones |
| Translated two-divisor chart | Five 13--49-term families | Five exact bases `[1]`; proves symmetric Laurent boxes equal ordinary bounded-degree polynomial spaces | Does not exclude nonsymmetric valuation bands |
| Newton-directed Stage D | `(3,2)`, `(2,3)`, and both primitive `(2,7)` orientations | Four exact collision-normalized bases `[1]`; identity Keller control remains in every kernel | Prototypes, not the full admissible chain |
| Dense Stage D control | Total-degree edges 6/9, 52 parameters | First modular F4 run exceeded 30 seconds | Timeout is not emptiness or a survivor |
| Published `(9,27)` regression | Nine displayed Section 5 equations | Exact elimination reproduces equation (5.9); valuation/degree contradiction reproduced | Assumes the paper's preceding reductions |
| Reduced `(72,108)` target | Proposition 4.3 polygon systems | Broad sizes `186/302` and `72/92` variables/equations; smaller F4SAT run timed out | No candidate and no exclusion |
| Experimental `(72,108)` edge powers | Common square/cube boundary substitution | Sizes reduce to 163 and 49 variables; 49-variable modular run timed out | Exhaustiveness depends on an unfinished audit of leading-form implications |

Here “exact basis `[1]`” means a stored characteristic-zero msolve input and
output proving that the named algebraic system has no solution over an
algebraic closure of `Q`. Agreement at several primes is a cross-check, not the
proof.

The detailed finite-family ledger is
[IMPLEMENTATION_STATUS.md](notes/IMPLEMENTATION_STATUS.md). The current
`(72,108)` boundary is in
[DEGREE_72_108_STATUS.md](notes/DEGREE_72_108_STATUS.md).

## Important corrections learned during the search

Several tempting broad claims were falsified or refined by positive controls:

- A single exceptional divisor can polynomialize two algebraically independent
  reciprocal-Jacobian coordinates. For
  `Phi_p=(y+p(1/x),1/x)`, the planted inverse
  `(A,B)=(1/v,u-p(v))` pulls back to `(x,y)`. Any obstruction must use
  collision, nonbirationality, or stronger boundary data.
- Pole cancellation is generally not the final obstruction. The current thin
  Stage D kernels are nonempty, contain an identity Keller control, and can
  produce the required bracket exponent; inconsistency appears only after the
  normalized collision and full nonlinear bracket equations are imposed.
- Symmetric Laurent boxes in the translated chart do not exploit infinity.
  The box `[-n,n]^2` polynomializes exactly to ordinary polynomials of total
  degree at most `n`. Enlarging boxes merely recreates dense coefficient
  search.
- Finite characteristic creates fake Keller behavior. Modular results are used
  for ranks, F4 solving, and lifting candidates, never as standalone evidence.
- Collision normalization is coordinate-sensitive. The normalization
  `F(0)=F(1,0)=0`, `DF(0)=I` is valid in the original affine plane search, but
  cannot be copied after the rational transformations in Proposition 4.3
  without transporting those points through the charts.
- A solver timeout is only a complexity datum. The 52-parameter dense control
  and the 72- and 49-variable `(72,108)` systems remain mathematically
  unclassified.

These corrections are as important as the negative searches: they determine
which chart families and support enlargements are worth trying next.

## Installation

The exact Python code uses SymPy 1.14 and mpmath:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

The Gröbner backend is msolve 0.10.1 or a compatible recent version available
as `msolve` on `PATH`. On macOS with Homebrew:

```bash
brew install msolve
msolve -V
```

Julia and HomotopyContinuation.jl are optional and used only for numerical
continuation controls:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Project.toml` and `Manifest.toml` pin the tested Julia environment.

After setup, the main exact regression suite is:

```bash
make verify
```

## Reproduction guide

### Core exact geometry

```bash
python3 scripts/verify_counterexample_independent.py
.venv/bin/python scripts/verify_counterexample.py
.venv/bin/python scripts/rediscover_3d_chart.py
.venv/bin/python scripts/cubic_model.py
.venv/bin/python scripts/image_nonproperness.py
.venv/bin/python scripts/weighted_lift.py
```

### Arithmetic and dynamics

```bash
.venv/bin/python scripts/finite_field_distribution.py
.venv/bin/python scripts/finite_field_refinements.py
.venv/bin/python scripts/commuting_flows.py
.venv/bin/python scripts/low_degree_symmetries.py
.venv/bin/python scripts/analyze_gradient_infinity.py
.venv/bin/python scripts/nonproper_fiber_benchmark.py
julia --project=. scripts/nonproper_fiber_homotopy.jl
```

### Explicit normal forms

```bash
.venv/bin/python scripts/cubic_homogeneous_reduction.py
.venv/bin/python scripts/verify_cubic_homogeneous_counterexample.py
.venv/bin/python scripts/cubic_linear_reduction.py
.venv/bin/python scripts/verify_cubic_linear_counterexample.py
```

The construction scripts generate the large sparse JSON artifacts under
`results/`; the verifier scripts independently parse and check them.

### Plane-search validation ladder

```bash
.venv/bin/python scripts/validate_ladder.py
.venv/bin/python scripts/exhaust_pole_kernels.py
.venv/bin/python scripts/parameter_rank_scan.py
.venv/bin/python scripts/stage_a_support_families.py
.venv/bin/python scripts/stage_c_toric.py
.venv/bin/python scripts/stage_c_translated.py
.venv/bin/python scripts/verify_translated_box_theorem.py
.venv/bin/python scripts/stage_d_newton_bands.py --timeout 30
.venv/bin/python scripts/verify_newton_translation.py
.venv/bin/python scripts/newton_9_27_regression.py
.venv/bin/python scripts/verify_certificates.py
```

The degree 5--20 pilot is regenerated with:

```bash
.venv/bin/python scripts/scan_2d_charts.py \
  --output results/scan_2d_5_20.json
.venv/bin/python scripts/scan_2d_charts.py \
  --collision-normalized \
  --output results/scan_2d_5_20_collision.json
```

### Reduced `(72,108)` systems

Build the broad systems without launching F4SAT:

```bash
.venv/bin/python scripts/build_72_108_reduced_system.py
```

Build the explicitly experimental edge-power refinement:

```bash
.venv/bin/python scripts/build_72_108_reduced_system.py --edge-powers
```

The script supports bounded modular runs through `--run`, `--prime`, and
`--timeout`. A timeout must remain recorded as unclassified.

## Repository map

### Exact search library

`jcsearch/` contains reusable, auditable components:

- `laurent.py`: Laurent dictionaries, forbidden-pole equations, exact linear
  ranks, RREF witnesses, and modular rank checks;
- `charts.py`, `canonical.py`: elementary rational charts and invariant-based
  word deduplication;
- `triangular.py`: the one-divisor `Phi_p` family and polynomialization kernels;
- `toric.py`, `translated.py`: two-divisor charts and exact simultaneous
  cancellation;
- `newton.py`: Newton-edge/band generation, monomial-to-Laurent translation,
  polygon lattice masks, and the reduced `(72,108)` polygons;
- `msolve.py`: small F4/F4SAT/elimination wrapper with exact input export;
- `arithmetic.py`: CRT, Hensel/Newton lifting, and rational reconstruction;
- `filters.py`: conservative generic-degree, collision, line, Newton, and
  Puiseux screens;
- `weighted.py`: weighted-lift validation controls inherited from the 3D
  mechanism.

### Result artifacts

- `results/certificates/*.input` and matching `*.msolve` files are retained
  characteristic-zero proof artifacts. `verify_certificates.py` reruns both
  unit-ideal and elimination certificates.
- `results/systems/` contains generated msolve inputs for the reduced
  `(72,108)` systems; their final polynomial is the vertex saturator and they
  must be invoked with `-S`.
- [results/README.md](results/README.md) explains artifact retention and
  regeneration.
- `results/OBSTRUCTION_CLUSTERS.md` and `results/SCAN_SUMMARY.md` are compact
  human-readable summaries.
- Top-level `results/*.json` files are generated output and are intentionally
  ignored by Git. Regenerate them with the scripts named in this README. This
  includes the large 95D/510D sparse artifacts and detailed scan records.

### Additional scripts and controls

- `search_rational_pairs.py` and `solve_rational_ansatz.py` are the original
  small rational-composition enumerator and symbolic ansatz solver. Their
  survivors are deliberately planted triangular automorphisms.
- `adaptive_chart_generator.py` enumerates short Cremona/toric words and
  deduplicates them by necessary valuation/Jacobian/lattice invariants. The
  current depth-two run gives 20 sampled invariant classes; the signature is
  not a complete equivalence decision procedure.
- `pole_cancellation_basis.py`, `exhaust_pole_kernels.py`, and
  `parameter_rank_scan.py` build and stratify exact Laurent cancellation
  matrices.
- `minimal_unsat_core.py` greedily reduces the degree-three contradiction to a
  19-equation exact core. It is not claimed inclusion-minimal.
- `modular_groebner.py`, `export_homotopy.jl`, and `hensel_demo.py` validate the
  F4, numerical discovery, lifting, and reconstruction stages on controls.
- `puiseux_screen.py` is only a leading-term Puiseux control; full expansions
  and dicritical resolution are not implemented.
- `keller_tangent.py` explores the linearized degree-at-most-seven Keller
  locus and bounded coordinate-orbit directions. It is exploratory and is not
  used as an obstruction certificate.
- `cluster_obstructions.py` regenerates the compact human-readable obstruction
  report from the JSON records.

### Notes

- Audit order and research gates: `FIRST_CHECKLIST`.
- Verified map and geometry: `FACTS`, `CONSTRUCTION`,
  `IMAGE_AND_NONPROPERNESS`, `FINITE_FIELD_VALUE_DISTRIBUTION`.
- Dynamics and continuation: `COMMUTING_FLOWS`, `GRADIENT_INFINITY`,
  `NONPROPER_FIBER_BENCHMARK`.
- Normal forms and consequences: `CUBIC_HOMOGENEOUS_REDUCTION`,
  `CUBIC_LINEAR_REDUCTION`, `DIRECT_CONSEQUENCES`.
- Plane methods and finite searches: `METHODS`, `TRIANGULAR_CHART_STATUS`,
  `TRANSLATED_TWO_DIVISOR`,
  `NEWTON_BAND_STAGE_D`.
- Current large-degree boundary: `NEWTON_9_27_REGRESSION`,
  `DEGREE_72_108_STATUS`.
- Claim boundaries and provenance: `IMPLEMENTATION_STATUS`, `SOURCES`.

## Current research boundary

No plane counterexample candidate has been found. No exclusion of either
reduced `(72,108)` case has been proved. The next useful mathematical step is
to audit and encode the common-leading-form, valuation, Puiseux, and dicritical
constraints used in the published admissible-chain reductions. Those
constraints must shrink the current 49--186 generic coefficients into a thin
recursive system before more F4 or homotopy is likely to be informative.

The announcement and some structural commentary are same-day material. The
algebraic identities stored here are independently checkable; provenance,
historical priority, peer review, and broad claims not reduced to exact
certificates remain unsettled. See [SOURCES.md](notes/SOURCES.md) for the dated
source list.
