# From one collision to marked-root Keller maps

This repository verifies a three-dimensional polynomial Keller map with a
three-point collision, explains it through a tangent-map normal form, and
develops weighted, cancellation, decorated-normalization, Hurwitz, and
rank-two symplectic consequences.  In generic degree `N>=4`, the coarse
decorated normalization already gives an `(N-3)`-dimensional family of stable
classes.  Adding one affine root sheet generically recovers the seed exactly.
The same parameters yield four-real-Gaussian witnesses with an optimal
`(N-3)`-moment algebraic fingerprint; in those moment coordinates, the full
affine nonsurjective partition-lattice geometry and every other functorially
defined seed locus are preserved in their native categories.  In degree six,
the two vertical Ritt loci become explicit sextic moment hypersurfaces, one
coinciding with the all-double nonsurjective component.  The same parameters
also yield inequivalent filtered Weyl endomorphisms.  The nonsurjective seed locus has
dimension `floor(N/2)-1`: quartics are generically exceptional, whereas
generic seeds of every degree `N>=5` are surjective.
The [exact degree-spectrum corollary](verified/GEOMETRIC_DEGREE_SPECTRUM.md)
shows that the geometric degrees of noninvertible Keller maps of
`A^3_C -> A^3_C` are precisely `3,4,5,...`.

## The foundational map

Put `u=1+xy` and define

\[
F(x,y,z)=\left(
u^3z+y^2u(4+3xy),
y+3xu^2z+3xy^2(4+3xy),
2x-3x^2y-x^3z
\right).
\]

On `x!=0`, set `t=y+1/x` and `r=2/x`.  For target coordinates `(a,b,c)`,

\[
a=t^2+rt/2-ct^3,\qquad b=r+4t-3ct^2.
\]

The two coordinate Jacobians give

\[
\det\frac{\partial(a,b,c)}{\partial(t,r,c)}
\det\frac{\partial(t,r,c)}{\partial(x,y,z)}
=\frac r2(-2x)=-2,
\]

so `det DF=-2` polynomially everywhere.  Exact substitution gives

\[
F(0,0,-1/4)=F(1,-3/2,13/2)=F(-1,3/2,13/2)=(-1/4,0,0).
\]

The minimal independent replay is:

```bash
python3 scripts/verify_counterexample_independent.py
```

## Canonical proof path

The common geometric normal form is

\[
(W,s)\longmapsto (s,Ws-H(W)).
\]

The [tangent-map core](verified/TANGENT_MAP_CORE.md) supplies the inverse
pencil, generic degree, critical divisor, discriminant normalization, full
Hessian Fitting divisor, incidence form, weighted suspension, and comparison
with cancellation suspension.  Later papers cite this theorem instead of
rederiving those calculations.

Its [all-degree rational-fiber corollary](verified/ALL_DEGREE_RATIONAL_FIBERS.md)
gives, for every `N>=3`, a Keller map with a complete fiber of exactly `N`
distinct rational points and an open real target neighborhood with `N` real
sheets.  The construction uses explicit integer roots and is uniform in `N`;
the former computation through degree 100 is retained only as a regression.
The [real-sheet spectrum theorem](verified/REAL_FIBER_SPECTRUM.md) sharpens
this to the exact chamber spectrum `N,N-2,...,N mod 2`, proves the minimum is
zero in even degree and one in odd degree, supplies rational targets for every
count, and exhibits the full parity chain by successive fold crossings.

The [standalone universal-monodromy theorem](verified/UNIVERSAL_SYMMETRIC_MONODROMY.md)
proves geometric and arithmetic `S_N` monodromy for every polynomial pencil
`H(W)-sW+t`, including monomial, Chebyshev, decomposable, symmetric, and
critical-value-collision cases; no generic-seed hypothesis is present.  The
[universal weighted-seed theorem](verified/WEIGHTED_SEED_THEOREM.md) applies
it to every admissible weighted map.
Its [effective finite-field Chebotarev corollary](extended-geometry/FINITE_FIELD_CHEBOTAREV.md)
realizes every prescribed cycle/factorization type over all sufficiently
large certified good fields.  It gives the split density `1/N!`, irreducible
density `1/N`, the random-permutation fixed-point law and its effective
`Poisson(1)` limit, plus a deterministic generator for modular Keller-fiber
witnesses and small rational lifts carrying the selected modular fingerprint.
Combining this with the real chambers by constructive weak approximation
gives the
[adelic complete-fiber theorem](verified/ADELIC_FIBER_ENGINEERING.md): every
degree and signature occurs as a complete fiber field of the explicit `F_N`,
and finitely many additional unramified splitting conditions may be imposed
simultaneously at sufficiently large good primes.
The [Hasse-principle fiber theorem](verified/HASSE_PRINCIPLE_KELLER_FIBER.md)
goes in the complementary arithmetic direction: one explicit degree-eight
complete regular fiber has points over `R` and every `Q_p` but no rational
point.  Its integral target is `(12138,-308652,1)`, and the full fiber is the
finite etale scheme cut out by an elementary intersective polynomial.
The [stratified extension](extended-geometry/STRATIFIED_ADELIC_FIBER_ENGINEERING.md)
moves the seed and target together on any rational chart of a selected locally
closed seed stratum. It isolates geometric nonemptiness as the only extra
compatibility condition and gives an exact nonsurjective quintic of omitted
type `2^1 3^1` and trivial Hessian symmetry whose three allowed real
signatures all occur with Frobenius types `(5)` at `7` and `(2,2,1)` at `11`.

The primary dependency chain is:

1. `F1` — foundational Keller collision.
2. `W1` — tangent-map core and weighted suspension.
3. `S1` — stable normalization functoriality.
4. `WB1` — weighted clean-locus intrinsic boundary reconstruction.
5. `C1` — universal cancellation construction.
6. `B1` — complete canonical boundary exhaustion.
7. `P1` — cancellation reconstruction residue and parameter faithfulness.
8. `M1` — finite degreewise stable multiplicity.
9. `D1` — marked Hessian-divisor moduli of dimension `N-3`.
10. `F2` — generic affine-mark faithfulness.
11. `H1` — internal Hurwitz--LL compactification theorem.
12. `R1`, `R2` — all-degree rank-two descent and parameter faithfulness.

The exact scopes, dependencies, checkers, and review states live only in [`MATH_STATUS.json`](MATH_STATUS.json).  [STATUS.md](STATUS.md) is generated
from it.  A checker is reproducibility evidence; it is not external review.

## Main papers

- [Foundational counterexample](papers/core-counterexample/main.pdf)
- [Discriminant pencils](papers/discriminant-pencils/main.pdf)
- [Decorated discriminant normalization](papers/decorated-discriminant-normalization/main.pdf)
- [Marked-root degreewise multiplicity](papers/marked-root-multiplicity/main.pdf)
- [Hurwitz--LL rerooting](papers/hurwitz-ll-rerooting/main.pdf)

The degree-five and degree-six calculations are worked regressions generated
from the all-degree results, not independent theorem authorities.  Quartic
models and external-map classifications are likewise examples or audits.

## Reproduction

Create the pinned Python environment in [REPRODUCE.md](REPRODUCE.md), then run:

```bash
make check verify-minimal verify-master verify-theorems verify-papers
```

`make check` validates links and the status graph.  `verify-minimal` needs only
the Python standard library.  `verify-master` covers the cancellation chain,
including the unconditional all-`(m,r)` thick-contact formula.  The conditional
contact-resultant refinements are separate from the `M1` proof path.
`verify-papers` discovers and compiles every `papers/*/main.tex`; CI uses the
same discovery rule and archives every resulting PDF.

The full command catalogue, heavier symbolic runs, optional Lean certificate,
and independent Macaulay2 comparison are documented in
[REPRODUCE.md](REPRODUCE.md).

## Research status

The H1--H3 chain has been adversarially audited.  The finite-cover valuative
lemma, quotient-stack map, and local collision algebra are proved, including
simultaneous clusters and relative stabilizers.  However, the global
comparison with the repository-specific admissible-cover contraction is not
yet constructed.  H1 and H2 are therefore partial, and H3 is plausible but
high-risk rather than foundational infrastructure.  See the standalone
[DVR marking audit](papers/hurwitz-ll-rerooting/dvr-marking-audit.tex).

The LR continuation is now split into three independent primary problems,
alongside the cancellation and Hessian--Ritt frontiers:

- `OP-CR`: cancellation contact resultants in the residual wedge `m>=7`,
  `r>=7`; the first six fixed-`r` columns, the six all-`r` columns `m<=6`,
  every fixed-`m` eventual tail, and all other parameter-irreducibility
  ranges are complete.  The `r=6` column is closed by 7424 rigorous Arb
  branch tubes on `0<=1/m<=1/41` and exact modular gcds for `m<=40`.
- `OP-LR-REES`: linear target-lift Rees strictness as a finite module/SAGBI
  problem, including a structural cutoff to finitely many torus weights.
- `OP-LR-II`: the minimal pair is now nonzero; classify the remaining
  generator matrices of `II_(F,p,-p)`, prove the finite-weight cutoff, and
  determine whether the minimal class generates the quadratic normal image.
- `OP-LR-NE`: valuative no escape through the finite proper closure of the
  marked normalized-cover `Isom` relation, using `F2`, `H2`, `H3`, and generic
  deck rigidity before reconstructing automorphisms.
- `OP-RITT`: `HR12` now proves the exact all-degree codimension formula and
  completes the degree-eight and degree-twelve diagrams, including the
  all-four Chebyshev intersection in degree twelve.  What remains is to extend
  the explicit diagrams beyond degree twelve and restrict the H3 boundary
  equations to their incidence components.

In particular, `b_m=34m+1` is an exact source-only profile in one
determinant-normalized target gauge.  It is not target-minimal and is not a
stable LR obstruction.  The better target torus produces a surviving
`v^(6m)S^(4m)` class in the invariant-ring-saturated equivariant target
quotient.  Modulo `gamma` its `v`-degree is `10m`, which already points to an
ordinary source-degree bound around `20m`; the exact torus-gauge slope `24` is
not needed.  The open algebraic step is to identify this as a universal
associated-graded obstruction after arbitrary lower gauges.  Formal
target-lifting now shows that all lower source jets are uniquely forced by the
target jets.  At order two, logarithmic coordinates reduce the new gauge
dependence to the affine-linear bracket class `Xi_2(Y_1)`, followed by one
explicit quadratic coordinate reconstruction.  The entire five-parameter
degree-25 homogeneous weight-zero kernel has now been audited: its order-two
quadratic residue never vanishes.  A filtered-coset BCH lemma now shows that
pairwise opposite-weight expansion is unnecessary at order two if the Rees
degeneration of the target-lift coset is strict.  The next implementation is
to construct the graded target-field, lifted, and normal modules and test
their generators by Gröbner/module membership, not to enlarge the symbolic
jet calculation.  Its quadratic obstruction is the explicit second
fundamental form
`II_F(Y_1,Y_2)=-(DF)^(-1)D^2F[ell_F(Y_1),ell_F(Y_2)]`; opposite weights reduce
to invariant-module maps `II_(F,p,-p)` into the weight-zero normal module.  The
smallest possible pair is already nonzero: the constant fields
`(partial_B,partial_C)` have weights `(1,-1)` and saturated normal symbol
`-146880u^5/7`.  Quadratic strict descent therefore fails, and this class is
the canonical quadratic LR invariant anticipated by `OP-LR-II`.  Valuative no
escape is a separate marked-cover problem, not an automorphism-coefficient
estimate.  None is a bottleneck for stable moduli, which already follow from
decorated normalization and the affine sheet.

Other questions—arithmetic Galois theory, controlled-boundary suspension
classification, wider quantization, coefficient-scheme gluing,
quadratic--cubic flexibility, the plane degree frontier, Gaussian-equivalence
and three-real-variable descent, and further Image/Vanishing consequences—are
retained as parked side programmes in
[STATUS.md](STATUS.md).

## Provenance

The earliest public item located by the repository audit is
[Levent Alpöge's post](https://x.com/__alpoge__/status/2079028340955197566).
Contemporaneous sources attribute the example to Alpöge and Fable; Dean
Cureton supplied a separate
[Lean certificate](https://github.com/deancureton/jacobian), and Alexis
Gallagher developed the weighted-lift viewpoint in a
[same-day explainer](https://jacobianfun.org/jacobian-explained).

The surviving record does not settle the complete discovery history, original
prompt, model conversation, or exact UTC timestamp.  This repository repeats
the contemporaneous attribution without making a priority claim.  See the
[provenance audit](archive/legacy-notes/PROVENANCE_AUDIT.md) for the detailed
record and the canonical sources for precise mathematical claims.

## Repository policy

Proofs belong in canonical sources, statuses belong in `MATH_STATUS.json`, and
open-problem IDs belong in `STATUS.md`.  Papers and notes may cite those IDs
but do not maintain independent continuation queues.  Superseded derivations
remain available under `archive/` and outside primary navigation.
