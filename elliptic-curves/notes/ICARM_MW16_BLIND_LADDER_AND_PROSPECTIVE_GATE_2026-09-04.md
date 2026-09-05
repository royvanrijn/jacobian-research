# MW16 blind ladder and prospective parent-pencil gate

Date: 2026-09-04; arithmetic-model audit updated 2026-09-05.

The later [sensitivity calibration](MW16_SENSITIVITY_RECOVERY_2026-09-05.md)
recovers 55/55 control directions with the new backend. The 54/55 historical
ladder and the original twenty-direction pointed-sieve calibration below
remain retained baselines.

Status: exact five-curve calibration complete; first prospective attempt
wholly timeout-censored; subsequent direct target-free search complete.
The exact model/section-size audit below separates coefficient reduction from
backend execution. The specialized pointed-quartic sieve also completes all
856 charts, using exact denominator and slope-lattice transforms.

## Exact presentation count

The nine atlas hit labels are not nine fibrations.  Exact rational base changes,
constant square Weierstrass scalings, target-parameter transport, and integral
specialized-MW16 comparisons give:

| target curve | atlas priorities | exact fibration classes |
|---:|---|---:|
| 398 | 16,875; 63,669 | 1 |
| 400 | 53,042; 62,992 | 1 |
| 401 | 57,487 | 1 |
| 542 | 30,486 | 1 |
| 548 | 31,627; 54,835; 63,647 | 1 |

Thus the correct statistical unit is five target curves.  The nine labels are
still useful as nine search coordinate systems: affine changes of the base do
not preserve a bounded projective-height box.  Results from those coordinate
systems are nested trials and are never counted as independent observations.

The five representative fibrations are mutually separated.  Nine of ten
pairs have frozen modular `j`-map landmark witnesses; the remaining
curve-398/curve-548 pair has different exact generic MW16 half-lattice coset
histograms, which is incompatible with a fibration isomorphism.

## Complement-blind five-curve ladder

The search fixture contains each target equation, its sixteen specialized
generic sections, and the exact generic MW16 Gram.  It contains no public
exceptional point, complement coordinate, target rank lower bound, or jump
label.  The first stage completely enumerates `M/2M` and searches every class
in the exact maximum-depth stratum.

| target | demonstrated atlas jump | deepest initial charts | blind initial gain | blind adaptive gain | best blind gain |
|---:|---:|---:|---:|---:|---:|
| 398 | +14 | 12 | +5 | +9 | +14 |
| 400 | +12 | 4 | +5 | +7 | +12 |
| 401 | +11 | 8 | +10 | not run | +10 |
| 542 | +10 | 10 | +10 | not needed | +10 |
| 548 | +8 | 8 | +8 | not needed | +8 |

Repeated coordinate presentations give the same specialized integral MW16
subgroup and the same initial exact recovered points.  They therefore do not
inflate the table.

At the five-curve level, the initial wave recovers 38 of the 55 demonstrated
quotient directions, completely recovers two curves, and comes within one
direction on three.  The previously completed curve-398 five-bit adaptive
wave and the new curve-400 five-bit adaptive wave raise this to 54 of 55.
The curve-400 run searches all `4(2^5-1)=124` adaptive charts at height
100,000; all 124 complete, and exact group-law plus finite-reduction
certificates give `M16 -> M21 -> M28`.

Only curve 401 remains one direction short of its demonstrated atlas lower
bound.  Its literal complete next wave has `8(2^10-1)=8,184` charts and has
not been run.  This purposive five-curve panel is detector calibration, not a
population success-rate estimate.  Matching an atlas lower bound also does
not prove that a target has no further rational directions.

## Prospective ladder

The enforced pipeline is:

```text
cheap Nagao/local ordering
    -> exact bounded half-lattice jump recovery
    -> complete same-minimal-curve residual 2-Selmer gate
    -> only then any expensive continuation
```

There is no Nagao-to-unrestricted-point-search edge.

The first prospective pass evaluates every primitive projective parameter of
height at most 300 in each of the nine coordinate presentations: 109,592
parameters per presentation.  Three disjoint prime blocks and frozen
within-bucket caps leave 104 finalists across the nine presentations.  Nagao
is used only for ordering; it contributes no rank evidence.  Exact
specialization transports all sixteen generic sections, rechecks the
candidate ledger, and finds 104 distinct `QQ`-isomorphism classes inside this
bounded selection.

The first half-lattice attempt covers all 104 finalists and all 856 exact
maximum-depth charts at rational height 10,000 with a two-second per-chart
budget.  Every specialized generic MW16 remains independent, but every one of
the 856 quartic processes times out.  Consequently:

- exact new quotient directions recovered: 0;
- fail-closed structural rejections: 0;
- completed height boxes: 0;
- timeout-censored chart attempts: 856;
- residual-Selmer calls authorized: 0.

This is a complete execution of the declared attempt, not a complete
height-10,000 search.  It supplies no negative rank evidence and does not
evaluate the Nagao ordering.  The raw specialized short models have maximum
rational coefficient size 5,871--6,302 bits, and the generated quartic
coefficients reach 10,140--11,302 bits.  Larger search budgets are therefore
not evidence for increasing the timeout on that representation. The completed
model/section audit below addresses the proposed arithmetic reduction; every
transformed candidate re-proves all sixteen section identities and
independence.

## Subsequent direct execution

The anonymous target-free candidate ledger reconstructs exactly the same set
of 104 raw equations as the earlier finalist ledger. Its worker bypasses
quartic minimalization/reduction and calls `hyperellratpoints` on the exact
integral quartic. The retained
`a1_mw16_target_free_parameter_search_h300_v1.json` records all 856 boxes
completed at height 100,000, with zero affine points or quotient gains and
no timeout. This supersedes the claim that no prospective box has ever
completed; it does not turn the historical timeouts into completed searches.
The result remains a bounded null experiment, with no upper rank bound.

## Exact arithmetic-model audit

`prepare_mw16_short_models.sage` consumes only the anonymous 104-fibre ledger
and its frozen deepest masks. It computes Sage/PARI's rational global minimal
model, explicitly stores `u,r,s,t` and all sixteen transported sections, and
checks every inverse map. It additionally measures the raw, rational short,
integral short, completed-square, and section-centered models. The selected
short model receives a fresh rank-16 finite-reduction certificate and an
explicit no-rational-2-torsion prime.

For every one of the 856 frozen masks, exact group addition replays the chart
base point in both models. A menu of at most 31 rational binary coordinate
changes tests denominator and real-size scalings, translations by slopes of
generic sections, reciprocal charts, and the projective change determined by
three generic-section slopes. Integer normalization removes exactly verified
square content using fixed small-prime valuations and a residual square-root
test, without factoring a discriminant. Content `2*k^2` must not cause the
whole `k^2` to be retained merely because the content itself is nonsquare.
Both model and chart
selection use only maximum/total coefficient bit sizes, with deterministic
ties; point counts, timings, Nagao values and exceptional points are absent
from selection. Candidate checkpoints precede all benchmarking.

All 104 global minimal models and all 1,664 section transports pass. Maximum
coefficient sizes, in bits, are:

| measured representation | minimum | median | maximum |
|---|---:|---:|---:|
| raw Weierstrass, rational coefficient height | 5,871 | 6,189 | 6,302 |
| global minimal Weierstrass, integral coefficients | 4,586 | 5,193 | 5,310 |
| raw pointed quartic, rational coefficient height | 10,140 | 11,005 | 11,302 |
| selected pointed quartic, integral coefficients | 6,262 | 10,287 | 10,693 |

The paired chart-by-chart reduction has median **6.14%**, with a best
reduction of **39.47%** (10,346 to 6,262 bits). Of 856 selected transforms,
756 retain the short coordinate, 23 use a nontrivial scale, and 77 use three
generic-section slopes. This does **not** meet the orders-of-magnitude size
target. The invariant argument below gives necessary integral-quartic
coefficient sizes of **1,522--1,768 bits**, depending on the fibre.

The frozen nine-chart benchmark completes all **36/36** calls at height
10,000 within the unchanged two-second cap. Observed total subprocess wall
times are **0.009--0.060 seconds**. The raw-direct and short-direct paths
report no affine points; the selected-direct and selected-reduced paths each
report four points. All eight reports have exact group-law expressions in
the supplied MW16 basis, stored in
`mw16_short_models_benchmark_verification_v1.json`; no quotient direction is
recovered. Raw-direct also completing quickly means model shortening is not
the explanation for executability. These timings describe this tiny panel,
not a performance guarantee for all 856 charts.

Pure `x` translation cannot shorten a fixed pointed quartic. On
`y^2=x^3+a2*x^2+a4*x+a6` its coefficients are

```text
t^4 - (6*xq+2*a2)*t^2 - 8*yq*t
    + a2^2 - 2*a2*xq - 3*xq^2 - 4*a4.
```

Substituting `x=x'+r` and transporting `Q` cancels every occurrence of `r`.
The audit checks this identity on every chart. Rational scaling and projective
changes can help; the same fixed fibre's `j` remains invariant even if its
family base parameter is re-expressed.

There is an exact obstruction to the requested orders-of-magnitude reduction
in **bit length**. For an integral binary quartic with coefficients of
absolute value at most `H<2^b`, its classical invariants satisfy
`|I|<=16*H^2` and `|J|<=137*H^3`, and

```text
j = 6912*I^3/(4*I^3-J^2).
```

Thus its reduced `j` numerator has at most `25+6b` bits and its denominator at
most `16+6b` bits: `6912*16^3<2^25` and
`4*16^3+137^2<2^16`. The certificate records the resulting necessary lower
bound for each fibre. It applies to every integral binary quartic representing
that curve, not merely the finite menu tested here. It does not assert that
the selected charts attain this bound or are optimal over `PGL2(Q)`.

## Specialized pointed-quartic execution (2026-09-05)

The new [pointed-quartic backend](POINTED_QUARTIC_SIEVE.md) makes no
`hyperellminimalmodel` or `hyperellred` call. It constructs the lattice
combination with GMP, removes the point's square/cube denominators by an
explicit congruence, balances an exact two-dimensional slope lattice, and
sieves rational slopes modulo small primes before large-integer square tests.
All transformations and returned points have exact rational replay.

The frozen 104-fibre ledger now completes all **856/856** maximum-depth boxes
at height **10,000**, with the unchanged **two-second search cap** and no
timeout or failed candidate. Integral quartic sizes are **1,537--1,789 bits**
(median **1,737**), close to the necessary invariant bounds in the preceding
audit. The modular workers take **53.08 seconds** in total and at most
**0.078 seconds** per chart. Including the exact group calculation and chart
transform, chart work totals **72.83 seconds**, with maximum **0.111 seconds**;
the separate generic-census and independence checks are additional work.
The sieve covers 171,208,560,000 integer numerator/denominator pairs and
performs only 28,134 large-integer square tests. It finds no finite point or
new quotient direction. This is a completed bounded null result in these
specific slope coordinates, not a rank bound.

A separate complement-blind **initial** control run completes 42/42 charts at
height 10,000 and proves gains **0, 0, 2, 10, 8** on curves
398, 400, 401, 542, 548 respectively. Independent group-law and finite-reduction
replay verifies these twenty directions. This does **not** reproduce the
historical 54/55 adaptive recovery: both the height bound and horizontal
coordinates differ from that calibration. The subsequent
[sensitivity calibration](MW16_SENSITIVITY_RECOVERY_2026-09-05.md) recovers
55/55 directions; speed alone did not establish that recovery.

The [compact summary](../../artifacts/generated-results/elliptic-curves/icarm_mw16_pointed_sieve_h10000_summary_v1.json)
pins the full deterministic gzip certificate, source hashes, software,
per-candidate coverage, and controls. The checker regenerates all 898
prospective/control chart base points and five-coefficient transformation
identities when invoked with `--replay-charts`. Small exhaustive tests check
the sieve, bad-prime cases, infinity, GMP group arithmetic, timeout semantics,
and checkpoint integrity. Incomplete chart checkpoints are replayed rather
than promoted to completed boxes.

## Replay

The current model/section audit and tiny benchmark are reproducible with:

```bash
sage -python elliptic-curves/cas/prepare_mw16_short_models.sage
sage -python elliptic-curves/cas/verify_mw16_short_models.sage
sage -python -m unittest discover -s elliptic-curves/tests -p test_mw16_model_size.py
```

The first command checkpoints each candidate under ignored
`artifacts/local/elliptic-curves/mw16-short-models-square-content/`, then writes
`mw16_short_models_h300_v1.json.gz` and its compact summary. It freezes the
median candidate by worst selected-chart size within each presentation and
the median chart by selected size on that candidate. Only those nine charts
are benchmarked, each on four paths: raw direct, short direct, selected
direct, and selected plus direct reduction. Each call has rational height
10,000 and two wall seconds, with no retry. Different coordinates define
different height boxes; their point counts are not interchangeable coverage
claims. `--prepare-only` omits search; `--benchmark-only` reuses a hash-matching
certificate; `--check` recomputes the deterministic full size audit without
search. The verifier replays the saved maps, sections, finite reductions and
quartic identities without rerunning the transform menu.

Historical calibration and prospective campaigns remain reproducible:

Build the complement-blind nine-presentation fixture and replay the exact
initial calibration:

```bash
sage -python elliptic-curves/cas/prepare_icarm_mw16_parent_ladder_inputs.sage --check
sage -python elliptic-curves/cas/run_icarm_mw16_parent_ladder_blind.sage
sage -python elliptic-curves/cas/audit_icarm_mw16_parent_presentations.sage --check
```

Replay the curve-400 adaptive calibration and the held-out five-curve
comparison:

```bash
sage -python elliptic-curves/cas/run_icarm_mw16_curve400_adaptive_calibration.sage
python3 elliptic-curves/cas/verify_icarm_mw16_blind_ladder_calibration.py --check
```

Replay the prospective local ordering and exact specializations:

```bash
python3 elliptic-curves/cas/sieve_icarm_mw16_parent_presentations_nagao.py --check
sage -python elliptic-curves/cas/specialize_icarm_mw16_nagao_finalists.sage --check
```

The current specialized runner supports `--candidate-start` and
`--maximum-candidates`, with additional atomic checkpoints after every chart.
Its replay command is:

```bash
sage -python elliptic-curves/cas/run_icarm_mw16_nagao_finalist_half_lattice.sage \
  --backend pointed-sieve --height-bound 10000 --timeout-seconds 2 \
  --output artifacts/local/elliptic-curves/mw16-pointed-sieve-h10000-all.json
sage -python elliptic-curves/cas/calibrate_icarm_mw16_pointed_sieve.sage
python3 elliptic-curves/cas/verify_icarm_mw16_pointed_sieve.py --check --replay-charts
```

The original eight contiguous 13-candidate timeout shards remain under
ignored `artifacts/local/`; their historical compact merged certificate is
unchanged. The specialized backend has separate outputs and source hashes.
See [the backend note](POINTED_QUARTIC_SIEVE.md) for packaging a fresh run and
replaying the retained gzip certificate without local checkpoints.

## Claim boundary

The five atlas jumps are rank lower-bound material, not exact target ranks.
The blind recovery certificates prove only the independent subgroup they
construct. The original prospective zero is entirely timeout-censored; the
subsequent direct run is a completed bounded null result. Model shortening
and benchmark timings imply no new independent direction. No candidate has
passed the residual-Selmer gate, no expensive continuation is authorized,
and rank 32 remains open.

<!-- status-consumer: EC-K3-ICARM-MW16-SENSITIVITY f88886c066d6cb45 -->
