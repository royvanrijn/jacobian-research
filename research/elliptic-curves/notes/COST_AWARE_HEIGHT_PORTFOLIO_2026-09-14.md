# Cost-aware boxes on the retained three-model chart

The first-chart certificate admits a smaller **address-work proxy** when its
models receive different box heights. At the same certified elliptic-height
range as the best single model at height 125000, the best allocation uses
heights `(21129,37298,84441)` on the existing full minimal, 5-neighbour and
7-neighbour models. Its sum of squared heights is **0.573943** of the single
box's, compared with **0.912676** for two equal integer-height boxes.

This is a finite optimization of the existing sufficient bounds. It is not a
CPU measurement, a new bound, a prediction of target accessibility, or a point
search. No new model is constructed. The
[frozen input](../../artifacts/generated-results/elliptic-curves/height_portfolio_cost_v1/protocol.json),
[allocation](../../artifacts/generated-results/elliptic-curves/height_portfolio_cost_v1/allocation.json)
and [independent replay](../../artifacts/generated-results/elliptic-curves/height_portfolio_cost_v1/verified.json)
retain the exact data and scope.

The [cold three-arm benchmark](HEIGHT_MODEL_NEXT_DIRECTION_BENCHMARK_2026-09-13.md)
remains stopped: its six bounded arms timed out during discriminant-support
preparation and made no point calls. The separate target diagnostics show that
the best uniform bound selected a smallest actual coordinate in only one of
three control cases, with a tie there. None of those twelve target elliptic
heights lay in the corresponding guaranteed range. This postprocessor resolves
neither preparation cost nor that pointwise slack, and does not replace the
sealed benchmark's selector.

## What was already optimized

The [original three-model selector](POINTED_CHART_HEIGHT_BOUNDS_2026-09-13.md)
already minimized `|T|^2 D_T`, where

\[
D_T=\max_{s\in S}\min_{i\in T}D_i(s).
\]

At a fixed multiplicative elliptic-height limit `X`, ignoring integer rounding,
this is equivalent to minimizing `|T| sqrt(X D_T)`: the address-work proxy of
equal-height boxes. Thus the previous rule already charged the number of models
under equal per-address prices. It did not minimize the coverage constant alone.

The extensions considered here are model-specific costs, activation costs, and
unequal box heights. Cost calibration itself is not supplied by a height proof.

## A finite coverage optimization

Fix one elliptic equation, one anchor `Q`, and a complete certified joint-state
table `S`. Suppose for every point in state `s` the complete model `i` satisfies

\[
H(t_i(P))^4\le D_i(s)H_x(2P-Q),\qquad D_i(s)>0.
\]

Keep the same model at every place for a given point. States must come from a
valid joint partition, not independently minimized local contributions. Safe
extra states may remain even if their global attainability is unknown.

Given `X>=1`, let integer `H_i>=0` be a model's box height; `H_i=0` means off.
The sufficient coverage constraint is

\[
\boxed{\quad
\forall s\in S\quad\exists i:\ H_i>0\ \text{and}\ H_i^4\ge X D_i(s).
\quad}
\]

Every point with `H_x(2P-Q)<=X` is then in at least one active model's box.
This includes the projective parameters covered by the input certificate.
The search implementation must enumerate these boxes and exceptional parameters
completely before claiming executed coverage. The condition neither guarantees
a new independent point nor containment of all the original search boxes.

A simple supplied cost model is

\[
\widehat C(\mathbf H)=C_{\rm bank}+C_{\rm selector}
+\sum_{i:H_i>0}(a_i+c_iH_i^2),
\qquad a_i\ge0,\quad c_i>0.
\]

`C_bank` includes constructing and certifying every model used to make the
selection, including models later discarded. Shared cold support is charged
once at the appropriate interface. `a_i` is an additional activation cost, not
a second charge for that construction. Unknown preparation cost stays unknown.
For genuine runtime estimates, prices must come from separately frozen,
target-blind calibration and must account for sieve, map, duplicate and
verification work. The quadratic cost model may itself be inaccurate.

For a fixed subset with one common box height, the continuous expression is

\[
C_{\rm bank}+C_{\rm selector}+\sum_{i\in T}a_i
+\left(\sum_{i\in T}c_i\right)\sqrt{X D_T}.
\]

Equal prices and zero activation recover the original rule. Common already-paid
bank costs do not change the selected allocation, but must not disappear from a
complete-cost comparison with V3. Choosing which models to construct is a
different, upstream problem; this calculation starts with an existing bank.

**Finite reduction.** Define exact integer thresholds

\[
q_{i,s}=\left\lceil (X D_i(s))^{1/4}\right\rceil.
\]

There is an optimal allocation with each `H_i` either zero or one of that model's
thresholds. Indeed, shrink a positive box to the greatest threshold it currently
meets, or turn it off if it meets none. This preserves every state it covers and
does not increase a monotone cost. Consequently at most `(|S|+1)^m` allocations
suffice for `m` models. This argument also applies to other supplied monotone
cost functions, though the implementation here supports only the displayed
activation-plus-quadratic form.

An independent finite proof assigns every state to one model that covers it.
For each of the `m^|S|` assignments, set a model's height to the greatest assigned
threshold, or zero if none was assigned. Any feasible allocation admits such
an assignment at no greater cost. Minimizing over assignments therefore gives
the same optimum. The producer uses box combinations; the independent checker
uses state assignments and replays the original real/local bounds and maps.

Optimality is relative to the supplied table and cost function. It is not
optimality among all equivalent models, all valid sharper inequalities, or
actual rational points. It optimizes the cost of completed coverage. CPU to the
first quotient gain additionally depends on search order, early stopping and
where independent points occur; no such distribution follows from these bounds.

## Exact first-chart instance

Inputs are the original three models and four states, with the retained
first-chart table hash bound in the protocol. No target coordinates, ranks of
discovered points or winning chart addresses enter selection. The chosen range
is exactly

\[
X=125000^4/D_{\rm best\ single}.
\]

This retains the earlier reference exposure. Set every `c_i=1`, `a_i=0` and
the common term to zero **in address-proxy units**. These zeros exclude
preparation from this illustrative metric; they are not measurements of zero
preparation time. The integer thresholds are:

| Joint state `(inside 5, inside 7)` | Full minimal | 5-neighbour | 7-neighbour |
|---|---:|---:|---:|
| `(0,0)` | 21129 | 83399 | 99912 |
| `(0,1)` | 55902 | 220653 | 37763 |
| `(1,0)` | 47246 | 37298 | 223409 |
| `(1,1)` | 125000 | 98679 | 84441 |

There are 125 candidate allocations and 93 feasible ones. Independent enumeration
of 81 state assignments gives the same minimum, including the constraints of at
most one or two active models.

| Allocation | Full minimal height | 5-neighbour height | 7-neighbour height | `sum H_i^2` | Ratio to single |
|---|---:|---:|---:|---:|---:|
| Best single | 125000 | off | off | 15625000000 | 1 |
| Original equal pair | 84441 | off | 84441 | 14260564962 | 0.912676 |
| Best at most two | 47246 | off | 84441 | 9362466997 | 0.599198 |
| Best at most three | 21129 | 37298 | 84441 | 8967857926 | 0.573943 |

The other equal-height pair, the 5- and 7-neighbours, ties the original pair's
coverage cost. The new diagnostic's lexicographic height tie-break chooses that
pair; it does not alter the old benchmark's frozen tie-break. The small change
from the earlier continuous ratio `0.912664` is exact upward integer rounding.

In the three-model optimum the full model covers `(0,0)`, the 5-neighbour covers
`(1,0)`, and the 7-neighbour covers both states inside the 7-ball. Adding the
5-neighbour does not improve the old common-height `D_T`, but allows the full
model's box to shrink. This is why selecting only by the final common bound
would miss its usefulness for unequal allocations.

The third model saves only `394609071` address-proxy units over the unequal pair.
With equal address prices, an added activation charge greater than that erases
this advantage. Different address rates can also change the selection; exact
synthetic tests exercise both cases. No empirical model prices are inferred
from these arithmetic counts.

## Replay, cost and stopping boundary

```sh
python3 research/elliptic-curves/tests/test_height_portfolio_cost.py -v
sage -python research/elliptic-curves/cas/verify_height_portfolio_cost.py
```

To regenerate the allocation without overwriting the retained result:

```sh
python3 research/elliptic-curves/cas/height_portfolio_cost.py \
  --output /tmp/height-portfolio-allocation-replay.json
```

The four tests check exact root boundaries, sensitivity to both price and
activation charges, independent exhaustive integer boxes on small examples,
and rejection of incomplete inputs or an exceeded node cap. The retained
checker also binds the exact input and source bytes and replays the upstream
three-model mathematics. Its distinct-prime guard preserves the separate
same-prime UNKNOWN case from the benchmark.

The [meter](../../artifacts/generated-results/elliptic-curves/height_portfolio_cost_v1/meter.json)
records 0.025773 process-tree CPU seconds for allocation and 0.409999 for the
independent replay, including their startup. Together they used about 0.529
elapsed seconds. These are postprocessing/replay costs only, not historical
model preparation, development or point search.

This calculation stops at a certified allocation example. It establishes a
cheap way to optimize a declared coverage budget on an already prepared chart.
Measured cold preparation remains unresolved, and the three target diagnostics
warn against treating global sufficient bounds as pointwise predictors. A new
complete-cost benchmark would require a specifically frozen preparation and
pricing policy; the earlier arms have not been rerun, and no production search,
class-constructor or ancestry lane is opened.
