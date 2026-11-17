# Programme 8: complexity and compression theory

## Status

This is research infrastructure and an initial canonical census, not a new
reduction theorem.  It turns the repository's currently scattered size data
into a typed database and records seven already justified monotone relations.
The database has one compulsory slot for every requested measure, but the
initial fifteen rows do **not** yet enumerate every construction in the
repository.  Missing values are recorded as `uncomputed`, `not_defined`,
`not_applicable`, or `not_established`; none is inferred from a nearby map.

The machine-readable inputs and generated audit are:

- [`complexity_compression_database.json`](complexity_compression_database.json);
- [`verify_complexity_compression_database.py`](../scripts/verify_complexity_compression_database.py);
- [`programme8_complexity_compression_report.json`](../artifacts/generated-results/programme8_complexity_compression_report.json).

`MATH_STATUS.json` remains the sole mathematical status authority.  This
programme cites its entries and does not create a second theorem ledger.

## 1. Why a typed database is necessary

There is no useful total ordering on the current constructions.  The
42-variable homogeneous quartic Hessian witness is smaller than the
44-variable witness, while the latter has smaller Hessian rank.  The
21-variable sparse cubic is smaller than both the index-optimized
22-variable cubic and the rank-optimized 24-variable cubic, but has worse
nilpotency index and rank.  A scalar "complexity score" would conceal these
tradeoffs.

Each construction row therefore contains the following independent fields:

| requested datum | database interpretation |
|---|---|
| source dimension | dimension of the displayed source or potential-variable space |
| coordinate-degree vector | exact ordinary total degrees, in displayed coordinate order |
| support vector | exact expanded coordinate-monomial counts, in displayed coordinate order |
| coefficient height | maximum reduced rational height in the displayed expanded basis |
| geometric degree | function-field degree, never inferred merely from a collision |
| monodromy | named group together with the action on generic sheets |
| boundary-prime count | non-affine height-one primes in a stated canonical normalization ledger |
| puncture count | deleted points on the stated normalized one-dimensional root parameter |
| stable-moduli dimension | dimension of the construction family's image in coarse stable left--right moduli |
| rank and nilpotency indices | operator-labelled data such as `JH` or `Hess(P)` |
| formal-verification status | proof-assistant coverage, separated from exact CAS checking |

Three conventions prevent false comparisons.

1. A potential's monomial count is not substituted for the support vector of
   its gradient map.
2. A finite-field sampled nilpotency index is a diagnostic, not an exact
   polynomial-matrix index.
3. Boundary-prime and puncture counts name the normalization chart and its
   exclusions.  They are not guessed from informal unit-rank language.

## 2. Initial coverage

The first census covers:

- the foundational weighted cubic;
- the all-degree clean weighted and quadratic-gauge families;
- the naive, jointly height-balanced, and coefficient-only optimized sparse
  quartic representatives;
- the original and sparse-conjugate 21-variable cubic witnesses;
- the index-18 and rank-17 circuit witnesses;
- the 42- and 44-variable homogeneous quartic Hessian witnesses;
- the promoted universal atomic maps;
- the six-variable doubled potential and its five-variable Schur descent.

The generated report currently gives:

| measure | populated theorem/exact/formula rows | explicit gaps or inapplicable rows |
|---|---:|---:|
| source dimension | 15 | 0 |
| coordinate-degree vector | 9 | 6 |
| support vector | 8 | 7 |
| coefficient height | 8 | 7 |
| geometric degree | 7 | 8 |
| monodromy | 7 | 8 |
| boundary-prime count | 5 | 10 |
| puncture count | 6 | 9 |
| stable-moduli dimension | 6 | 9 |
| complete operator rank-and-index pair | 4 | 11 |
| proof-assistant status | 1 complete, 1 partial | 13 without proof-assistant coverage |

The large boundary and monodromy gaps are informative.  The BCW and Hessian
artifacts are excellent sparse algebra certificates but do not yet export
the canonical finite-normalization data needed to compare them with the
three-dimensional marked-root constructions.

## 3. A transformation is an invariant contract

A transformation record contains:

\[
(\text{source construction},\ \text{target construction},\
 \text{preserved invariants},\ \text{selected objectives},\
 \text{side effects}).
\]

The checker resolves the cited field paths in the source and target rows.  It
requires equality on every machine-checkable preserved field, weak
improvement in every selected objective, and strict improvement in at least
one objective.

It also distinguishes:

- **equivalences**, which transport all intrinsic invariants in their stated
  category;
- **eliminations**, which preserve only a proved structural packet;
- **construction replacements**, which retain properties such as constant
  Jacobian and noninjectivity but need not identify the maps.

This distinction is essential.  Circuit redesign is a legitimate way to
lower the best known rank or index, but it is not a monotone endomorphism of
one stable left--right class.

## 4. Certified monotone relations already present

| relation | selected strict improvement | preserved packet | cost or limitation |
|---|---|---|---|
| rational quadratic-gauge moduli scaling | coefficient height `2248704 -> 21875` | dimension, degree vector, support, geometric degree `4`, `S_4`, boundary count, punctures, stable class | no support improvement |
| coefficient-only quadratic-gauge scaling | coefficient height `21875 -> 4648` | same intrinsic and displayed-combinatorial packet | collision-coordinate height worsens `19856 -> 49640`; `4648` is optimal only on the declared 25281-point grid |
| five-step linear conjugation of the 21-variable cubic | total coordinate support `81 -> 74`; height `48 -> 36` | dimension `21`, cubic degree vector, `rank(JH)=18`, index `19`, conjugacy class | support vector is not componentwise monotone; nonzero `JH` entries increase |
| circuit redesign for index | `index(JH): 19 -> 18` | cubic homogeneous determinant-one noninjective class of problems; rank remains `18` | dimension `21 -> 22`; not an equivalence |
| circuit redesign for rank | `rank(JH): 18 -> 17` | cubic homogeneous determinant-one noninjective class of problems; index remains `18` | dimension `22 -> 24`; not an equivalence |
| cotangent circuit redesign | Hessian rank `38 -> 37` | homogeneous quartic HN noninjective-gradient packet | dimension `42 -> 44`; not an equivalence |
| one-pivot Schur elimination | dimension `6 -> 5` | nonzero constant Hessian and the selected equal-gradient collision | geometric degree, monodromy, boundary, support, and height are not proved preserved |

The generated frontiers expose two strict within-class dominations:

\[
\begin{aligned}
F_{\rm quartic,coeff}
 &\prec F_{\rm quartic,bal}\prec F_{\rm quartic,naive}
 &&\text{by coefficient height with all other required metrics equal},\\
F_{\rm BCW,21,sparse}
 &\prec F_{\rm BCW,21,original}
 &&\text{by total support and height}.
\end{aligned}
\]

The collision-coordinate height is an auxiliary objective and reverses the
first quartic comparison between `coeff` and `bal`.  This is exactly why
selected measures must be named rather than collapsed into one score.

The remaining four relations move between different Pareto points.  They are
useful search operators precisely because they expose which invariant packet
survives the move.

## 5. First search programme

### 5.1 Torus height descent

The quartic height reduction is not an isolated numerical trick.  On a
proved moduli-torus orbit, every expanded coefficient transforms by a
character.  Prime by prime, the logarithmic height objective is therefore a
maximum of finitely many affine functions of the scaling valuations.  This
suggests an exact finite convex minimax problem:

1. expand the coefficient characters;
2. solve the real piecewise-linear minimax problem;
3. enumerate the adjacent rational or integral valuation lattice points;
4. certify the winner by exact expansion;
5. transport and audit the collision points separately.

This search preserves support automatically while targeting coefficient and
point height.  The first all-family experiment should be the split
quadratic-gauge seeds through degrees four to eight.

### 5.2 Sparse triangular left--right descent

The 21-variable conjugation shows that strict greedy descent is inadequate:
a term-neutral fourth move exposes the cancellations in the fifth move.
The appropriate local search object is therefore a certified two-move
neighborhood, not a list of individually improving shears.

Candidate moves should be ranked lexicographically by:

\[
(\text{total support},\ \text{maximum component support},\
 \text{coefficient height},\ \text{nonzero Jacobian entries}),
\]

while exact polynomial automorphisms guarantee preservation of geometric
degree, monodromy, and canonical boundary data even when those common values
have not yet been computed.

### 5.3 Rank-index circuit rewriting

The current cubic witnesses do not dominate one another:

| witness | dimension | `rank(JH)` | index |
|---|---:|---:|---:|
| sparse 21-variable | 21 | 18 | 19 |
| index-directed 22-variable | 22 | 18 | 18 |
| rank-directed 24-variable | 24 | 17 | 18 |

The sharp next target is consequently not another single-objective record.
It is one exact witness satisfying

\[
\boxed{n\le22,\qquad \operatorname{rank}JH\le17,\qquad \nu(JH)\le18.}
\]

Search should rewrite the degree-lowering circuit before homogenization and
carry the full power-rank tuple as a Pareto objective.  Final basis changes
cannot shorten the existing generic Jordan chain.

### 5.4 Coordinate-pair and kernel quotients

Coordinate-pair restriction and constant-kernel quotient are direct
dimension-reducing operators.  They preserve constant Jacobian and a
collision under their exact hypotheses, but they do not automatically
preserve geometric degree, monodromy, or canonical boundary counts.
Programme 8 should therefore attach a post-reduction recertification stage
rather than silently copying those invariants from the parent.

The first targets are the 22- and 24-variable circuit artifacts, whose exact
kernel and identity-output relations already provide candidate coordinate
pairs.

### 5.5 Schur compression

Schur elimination is the model example of a theorem-level compression
operator: it lowers dimension and preserves a precisely stated
constant-Hessian collision packet.  The existing obstruction results rule
out the pure-source and constant-dual second pivots for the quadratic-gauge
families.  Any continued search should therefore begin with mixed
source--dual or coisotropic pivots, not repeat the closed ansatzes.

## 6. Database completion order

The next census additions should be made in this order.

1. Export exact coordinate degree, support, and height vectors for the
   promoted universal maps in ranks three through eight.
2. Add the cancellation and power-shifted gauge families, including their
   boundary-prime formulas, punctures, nilpotent boundary contacts, and
   stable lattice invariants.
3. Add the twelve-variable coordinate-pair map and each intermediate
   MacFarlane/BCW reduction as explicit transformation edges.
4. Build a canonical-boundary export for the cubic-homogeneous and Hessian
   witnesses before assigning them boundary or puncture values.
5. Add formalization scope at theorem granularity rather than using a single
   repository-wide Boolean.

An entry is useful even when most fields are unknown: it makes the missing
certificate visible and prevents an unjustified invariant-preservation
claim.

## 7. Reproduction

Validate the database, replay the sparse artifacts, check every declared
monotone relation, and refresh the generated report with:

```bash
python3 scripts/verify_complexity_compression_database.py --write-report
```

Inspect the two within-class Pareto frontiers with:

```bash
jq '.comparison_frontiers' \
  artifacts/generated-results/programme8_complexity_compression_report.json
```

The checker is dependency-free.  It does not run Singular, Macaulay2, Lean,
or any long symbolic search; it consumes the exact artifacts already
certified by their canonical reproduction commands.
