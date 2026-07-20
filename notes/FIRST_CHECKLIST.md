# What to verify, discover, and analyse first

This is the entry checklist for a fresh audit of the repository. It orders work
by logical dependency: establish the small counterexample certificate and the
reproducibility boundary before relying on the larger geometric claims,
derived normal forms, or plane-search results.

Labels used below:

- **Verify**: check a stated result independently or reproduce a stored proof.
- **Discover**: obtain information or a construction the repository does not
  yet contain.
- **Analyse**: explain structure already visible in the exact computations.

## Most promising next steps

1. **Independent minimal certificate.** Remove the single-algebra-engine
   failure mode from the determinant, collision, and degree checks. Implemented
   below with a dependency-free sparse-polynomial verifier.
2. **Provenance and formula audit.** Establish the primary source and compare
   every repository occurrence of the map mechanically, including signs,
   coordinate order, and normalization.
3. **Independent exceptional-stratum geometry.** Reprove the image,
   nonproperness, and fibre statements without silently restricting to the
   generic rational reconstruction chart.
4. **Normal-form artifact audit.** Add the 95D and 510D independent verifiers to
   the declared regression suite, then audit the cited stable-equivalence steps.
5. **Published 2D reduction audit.** Verify the reductions preceding the
   `(9,27)` terminal system and every leading-form condition behind the reduced
   `(72,108)` systems before increasing solver budgets.

This order prioritizes short checks that secure all downstream claims before
open-ended classification or larger computation.

## Current verification run: 20 July 2026

`make verify` completed successfully from a newly created ignored `.venv` on
macOS with Python 3.14.2, SymPy 1.14.0, mpmath 1.3.0, msolve 0.10.1, and Julia
1.12.6. The first invocation failed only because `.venv` did not yet exist;
after installing `requirements.txt`, the complete declared suite exited zero.

The run reproduced:

- Python compilation and all 19 README relative-link checks;
- `det DF = -2`, the three-point rational collision, coordinate degrees
  `(7,6,4)`, and the displayed Palais--Smale escape curve;
- the primitive cubic identities, rational reconstruction, discriminant
  identity, simple-root and `x=0` fibre checks, singular-locus calculation,
  boundary elimination/normalization, and explicit root-meridian checks;
- the validation ladder, including msolve F4/F4SAT controls and the planted
  chart and weighted-lift controls;
- translated-box checks through `n=4`, Newton translation through exponents
  eight, and the cached Stage D outcome classifications;
- the scoped `(9,27)` terminal elimination over `Q` and its three modular
  cross-checks; as documented by the script, this assumes the preceding
  published Newton reductions; and
- all 24 retained characteristic-zero msolve certificates: 23 reduced bases
  `[1]` and the exact `(9,27)` elimination relation.

After that run, `verify_counterexample_independent.py` was added and wired into
`make verify-core`. It restates the map using a standalone integer sparse-
polynomial implementation, computes the full Jacobian determinant, evaluates
the collision with standard-library rational arithmetic, and checks total
degrees without importing SymPy or `jcsearch`.

This run establishes reproducibility of the repository's declared regression
suite on one local environment. It is not a second-CAS verification, an audit
of the literature reductions, an independent ideal-membership proof, or a run
of scripts omitted from the `make verify` target (including the 95D/510D
artifact generators/verifiers, finite-field scripts, dynamics scripts, Julia
homotopy benchmarks, and broad `(72,108)` builders).

## 1. Establish a trustworthy baseline

- [ ] **Verify** the provenance of the displayed 3D map: locate the earliest
  primary source, original search code or prompt if available, timestamps,
  random seeds, and intermediate candidates. Keep same-day exposition distinct
  from a citable or peer-reviewed source.
- [ ] **Verify** repository integrity before treating results as archival.
  Record a commit, Python/Julia/msolve versions, platform, commands, runtimes,
  and hashes for retained certificates and generated large artifacts.
- [x] **Verify** file-retention claims. Top-level generated result JSON/TXT
  files and local noise are now ignored; certificates, solver inputs, Markdown
  summaries, source, and pinned dependency manifests remain trackable.
- [x] **Verify** a clean-environment run of `make verify`, then run the 3D core
  scripts without msolve or Julia to make the minimal dependency boundary
  explicit. The 20 July 2026 run above passed; `verify-core` uses only the
  pinned Python environment. Preserve complete logs in a future archival run.
- [ ] **Verify** that verifier scripts do not merely re-use the same formulas,
  serialization assumptions, or helper routines as their generators. Add
  independent parsers/checks where shared failure modes remain.

**Exit criterion:** another machine can reproduce the finite core certificate
from a pinned checkout and can identify exactly which later claims require
msolve, Julia, generated artifacts, or external literature.

## 2. Independently certify the 3D counterexample

- [x] **Verify** coefficient-by-coefficient that `det DF = -2` in two unrelated
  implementations: SymPy and the local dependency-free sparse-polynomial
  verifier.
- [ ] **Derive** a short hand-checkable determinant proof using the weighted
  chart identities rather than only full symbolic expansion.
- [x] **Verify** the three rational substitutions and their pairwise
  distinctness directly from the displayed formula. These two checks alone are
  the minimal counterexample certificate. Reproduced exactly by
  `verify_counterexample.py` in the current run.
- [ ] **Verify** that the formula used by every script and note is identical,
  including signs, output order, determinant normalization, and coordinate
  degrees `(7,6,4)`.
- [x] **Verify** the primitive cubic identities and rational reconstruction.
  These exact identities were reproduced by `cubic_model.py`.
- [ ] **Verify** generic degree three independently. Separate what follows from
  one explicit three-point fibre plus etaleness from what depends on
  elimination or irreducibility.
- [ ] **Analyse** the smallest robust certificate format: displayed map,
  determinant proof, collision, and optional reconstruction identities. It
  should be reviewable without trusting repository infrastructure.

**Exit criterion:** the claim that the usual complex Jacobian conjecture fails
in dimension three no longer depends on a single CAS, a same-day website, or a
large generated artifact.

## 3. Audit the exact geometry before downstream use

- [ ] **Verify** irreducibility of `Q`, its nonsquare discriminant class, and
  the resulting non-normal cubic extension with `S_3` Galois closure over the
  function field. Make every characteristic and base-field assumption explicit.
- [ ] **Verify** the fibre classification on all exceptional cases, especially
  `c=0`, degree drops of the cubic, repeated roots, and points where the rational
  reconstruction formula has a pole.
- [ ] **Verify** both inclusions in
  `F(C^3) = C^3 \ Gamma` and both inclusions in `S_F = V(Q)`; check that no
  affine or boundary stratum is lost by clearing a denominator.
- [ ] **Verify** that `Gamma` is exactly the singular locus of `V(Q)` and that
  the proposed boundary parameterization is its normalization.
- [ ] **Discover** a full resolution of the rational extension between suitable
  compactifications, listing every boundary divisor, dicritical component,
  valuation, multiplicity, and incidence relation.
- [ ] **Analyse** how the one-, two-, and three-point fibres lose sheets at
  infinity without finite ramification, and reconcile this with the resolved
  boundary model and monodromy.

**Exit criterion:** all image, nonproperness, fibre, and monodromy claims are
proved across denominator and degree-drop strata, not only on the generic
chart.

## 4. Audit derived constructions and claimed consequences

- [ ] **Verify** the 95-dimensional cubic-homogeneous artifact from its sparse
  data: determinant one, nilpotent nonlinear Jacobian, homogeneity, support
  count, and transported collision. Independently audit each stable-equivalence
  step against the cited reduction theorem.
- [ ] **Verify** the 510-dimensional Druzkowski artifact: matrix dimensions and
  rank, pairing identities, determinant one, and collision transport. Check
  that generator and verifier are logically independent.
- [ ] **Verify** each implication in `DIRECT_CONSEQUENCES.md` against the exact
  statement, quantifiers, field, and dimension hypotheses in the primary
  literature. Keep existential consequences separate from explicit witnesses.
- [ ] **Discover** whether the stable-equivalence and pairing choices can lower
  dimensions 95 and 510 or reduce support while retaining a transparent exact
  collision certificate.
- [ ] **Analyse** polynomial equivalence of the weighted seed family and search
  for lower degree, fewer monomials, or lower coefficient height under affine,
  triangular, and tame automorphisms.

## 5. Audit the plane-search machinery before enlarging it

- [ ] **Verify** the precise finite scope of every negative result from its
  generator, stored input, solver flags, and output. No timeout, modular unit
  ideal, or sampled support family should be summarized as a characteristic-zero
  exclusion.
- [ ] **Verify** all characteristic-zero msolve certificates with a pinned
  solver version and an independent ideal-membership check for representative
  `[1]` outputs. Confirm that saturators encode exactly the stated nonvanishing
  conditions and do not remove valid components.
- [ ] **Verify** the collision and derivative normalization wherever it is used.
  Prove its validity in the original affine coordinates and do not transfer it
  through rational charts without transporting the selected points.
- [ ] **Verify** the claim suggested by the modular sweeps that the triangular
  pole-matrix rank depends only on `deg p` in characteristic zero, or record an
  explicit counter-stratum.
- [ ] **Discover** a complete equivalence test, or controlled normal form, for
  chart words beyond the current necessary invariant signatures. Quantify what
  the depth-two catalogue can miss through over-identification or
  under-identification.
- [ ] **Verify** the full chain of published reductions leading to the `(9,27)`
  terminal system and the two `(72,108)` Proposition 4.3 polygons. The current
  executable `(9,27)` check begins only at the displayed terminal equations.
- [ ] **Verify** every common-leading-form and square/cube edge-power implication
  before treating the 49-variable `(72,108)` refinement as exhaustive.
- [ ] **Discover** the correct transported collision or noninjectivity condition
  in the Proposition 4.3 coordinates, including charts where a selected affine
  point moves to a boundary divisor.
- [ ] **Discover** full Newton--Puiseux and dicritical-component analysis for
  survivors; the current one-leading-term screen is only a filter.

**Exit criterion:** each plane-search statement has a machine-checkable scope,
a characteristic-zero evidence class, and a documented path from chart
coordinates back to a possible polynomial Keller map.

## 6. Analyse the next mathematical bottleneck

Do these after the audits above, in roughly this order:

- [ ] **Analyse** which ingredient of the 3D construction is impossible or
  constrained on a surface: boundary dual graph, reciprocal valuations,
  dicritical components, nonproper value-set geometry, or function-field
  monodromy.
- [ ] **Verify** the 2024 prime-degree claim before using it to exclude cubic
  generic degree in dimension two. If valid, reformulate searches around the
  first allowed composite non-Galois degree.
- [ ] **Analyse** the common obstruction behind the exact empty finite families
  after quotienting out normalization artifacts and birational controls.
- [ ] **Discover** the valuation and leading-form constraints that reduce the
  broad 49--186 parameter `(72,108)` systems to thin recursive systems. Encode
  these before spending larger F4 or homotopy budgets.
- [ ] **Analyse** the resolved infinity geometry of any survivor before
  numerical continuation; use homotopy and Hensel lifting only after dimension,
  saturation, and component structure are understood.

## Lower-priority follow-ups

- [ ] Resolve the nonhyperbolic equilibrium line in the weighted gradient chart
  and identify basin boundaries.
- [ ] Compute a presentation of the full discriminant-complement fundamental
  group and compare its action with the resolved boundary model.
- [ ] Classify the weighted seed family for generic degree and polynomial
  equivalence.
- [ ] Turn finite-field enumeration checks into theorem-level independent
  derivations for every characteristic, retaining them as arithmetic structure
  rather than evidence for the characteristic-zero certificate.

These are valuable, but none should delay the reproducibility audit, the
minimal 3D certificate, or the denominator-stratum checks in the exact geometry.
