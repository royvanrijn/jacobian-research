# Performance gate for the next V3 transfer

The preparation narrative below is retained historical evidence. The adapter
and productive-parent searches have since completed; see the dated assessment
immediately below for current operational guidance. Historical runs and frozen
sources remain unchanged.

## Current assessment — 2026-09-09

The [fresh-six cohort](FRESH6_RETAINED_SEED_COHORT_2026-09-09.md) provides a
new population branch after the bounded continuation misses. Six equations
are selected from retained scores after exact measured/queued-curve exclusions,
and their generic17 packets pass replay. Extra-point searching has not begun
on this cohort; prior high-rank suffixes remain separately checkpointed.

Parent choice also needs diversification. A new
[302 mask audit and prospective complement-bank experiment](PRODUCTIVE_PARENT_SPAN_REASSESSMENT_2026-09-09.md)
shows that many winning classes in two completed302 cascades lie outside the
span of earlier winning classes. Pairwise closure cannot reach them. The new
curve48 branch samples and certifies generic classes outside its old productive
span in2.986 seconds; it does not transfer302 points or labels. This is a
separate bounded parent-selection experiment, not a predictor or a reason to
discard the earlier productive-parent successes.

For new native-seed passes, use the
[lean preconditioned runner](../cas/run_lean_preconditioned_seed_v3.py). It retains the
100-invocation budget, exact admission and immediate rebuild, but bounds map
construction separately at5 seconds/1GiB before the10-second point-search
timer. Factor-free integral preconditioning precedes full PARI minimization;
the independent factor-free coordinate lane remains available. Exact map
identities are replayed without rerunning minimization. A map timeout remains
unresolved and does not stall the next lane.

The [exposure preflight](../cas/prepare_lean_preconditioned_exposure.py) freezes the
same landscape before any point search and compares its rational centres up
to sign against named same-equation histories. This measures new exposure,
not point existence. The motivating stall, four fixed winning-point controls,
and bounded validation runs are recorded in the
[high-rank roster](HIGH_RANK_SHORT_PASS_ROSTER_2026-09-09.md).
Older runners and their frozen receipts are preserved. The older cached
continuation wrapper does not yet accept this new protocol schema; use it
only for the earlier compatible short/cached runs.

Lean passes now have a separate
[cached continuation wrapper](../cas/run_lean_cached_seed_v3.py). It creates
a new directory, binds the parent's independent landscape verification,
retains original map-worker command/path receipts and adds100 actual point
invocations. Inherited maps are always replayed at their original paths;
missing or censored inherited maps cannot silently restart a worker. The
extended clouds and all point witnesses are replayed again; a new basis still
requires fresh landscape construction and independent rational-CVP replay.
Three focused inherited-map routing/rejection tests pass in0.155 seconds.

On curve48 the first continuation adds100 boxes,200 cumulative, with no gain
or timeout. Search takes88.352 seconds; independent replay takes2.414 seconds
using the bound previously verified landscape, versus117.099 seconds for the
initial full landscape replay. These phases do different verification work:
the earlier full CVP proof is inherited rather than recomputed. See the
[sealed continuation](../../artifacts/generated-results/elliptic-curves/curve48_lean_cached_v3_v1/result.json).

The lean variant initializes PARI directly through `cypari2`, avoiding the
full `sage.all` import in every isolated map worker. A
[fixed equivalence audit](../cas/audit_lean_map_workers.py) checks both policies
on curve71 centre indices0,1,2,3,46,47,48,49. All16 complete within5 seconds,
match their retained mapping dictionaries exactly, and pass independent
quartic-identity verification. Aggregate worker time is2.476 seconds versus
10.660 recorded for those same historical maps, about4.3x on this fixed
component comparison. The timings include process startup. This is not a
controlled end-to-end search comparison or a universal equivalence proof.
Evidence: `artifacts/local/elliptic-curves/lean-map-fixed-controls-v1/`.
The three [receipt regressions](../tests/test_lean_map_receipts.py) also pass:
a forced startup timeout replays without launching a worker, the factor-free
fallback matches a retained map, and an altered quartic is rejected even
after its container hashes are rebound. The test run takes1.412 seconds.

For cross-policy preflights, the
[projective-box key](../cas/pointed_box_equivalence.py) additionally identifies
signed coordinate permutations `t`, `-t`, `1/t`, `-1/t`, which preserve
`max(|a|,|b|)` for reduced rational `a/b`. Apply this only to completed boxes
at the same bound which include infinity. Three tests pass. The first16
preconditioned maps proposed for curve92's completed pairwise policy all
match completed boxes under this check, so no new point search is launched
for that roster. A separate audit verifier replays the292 parent charts and
checks the16 map comparisons. See the
[curve92 preflight](CURVE92_PRODUCTIVE_V3_SEARCH_2026-09-09.md).
This helper is available for future preflights; existing frozen runners keep
their declared projective-coordinate deduplication rules.

The future runners now use exact integer CVP, explicit generic dimension,
cached finite-group point admission, immutable chart receipts and one complete
cloud audit per epoch. Independent replay uses the original rational CVP.
Both coordinate maps are retained: the verified curve113 gain is inside the
factor-free height bound but outside the minimized map's bound at its winning
centre. Finite-column failure remains UNKNOWN.

A read-only [cost accountant](../cas/audit_future_search_costs.py) checks the
ordered receipt chains and aggregates recorded component times for seven
completed searches. Its frozen report is at
`artifacts/local/elliptic-curves/future-v3-cost-summary-v2/report.json`.
All seven searches have passed independent replay. The earlier v1 report
retains its pending-replay snapshot. This report is timing accounting, not an
arithmetic replay or a controlled comparison with the old implementation.

| Seed / parent bank | Search seconds | Point-search wall share | Landscape seconds |
| --- | ---: | ---: | ---: |
| 90 productive | 525.047 | 87.71% | 31.818 |
| 92 productive | 247.423 | 85.08% | 15.541 |
| 200 productive | 613.019 | 86.22% | 42.079 |
| 113 productive | 517.262 | 84.47% | 41.439 |
| 116 productive | 482.586 | 85.04% | 36.208 |
| 188 productive | 637.509 | 85.02% | 59.414 |
| 188 pairwise | 2387.337 | 84.79% | 239.795 |

The remaining time includes admission, maps, certificates, receipts and
supervisor overhead; it has not been separately attributed. Even removing all
non-backend work would improve these recorded totals by only about1.14–1.18x
with point-search time held fixed. This is an arithmetic ceiling for these
runs, not a forecast. The earlier3.093x measurement applies only to CVP.

The verified gains on200,113,116 occurred at charts81,12,27 respectively.
After immediate basis rebuilds, their respective734,750,652-chart epochs found
no additional certified gain. These three selected successes do not establish
a general success probability. They motivate the next scheduling experiment:
a frozen100-invocation first pass over prepared, deduplicated seeds, immediate
rebuild on each certified gain, and a separate queue for deeper continuation.
An unfinished pass must retain its exact cursor and be reported as budget
exhaustion, never as a completed policy or dependence proof.

The [short-pass runner](../cas/run_short_seed_v3.py) now freezes100 total
invocations across all adaptive epochs, retaining the unchanged parent
runner's centre/map order, gain checks and independent replay. It uses separate
1800-second/3-GiB search and replay limits and one worker. A cutoff preserves
the full landscape, ordered per-chart receipts and final seed; it does not
implement automatic promotion to the deeper queue. This validates a bounded
execution policy, not its population-level yield.

Its first prospective seed is inventory52, native R17 family08f72 at
`-164/2855`, from `next24_r17_results_v1.json`, entry `08f72-053`.
Preparation independently re-certifies M25 and the seven productive masks
19918,57793,64272,73128,80170,94624,94859, in11.247 seconds under120 seconds/1GiB.
It is the remaining M25 candidate in the earlier exact-root+1 scheduling
screen. This is a heuristic priority, not a parity theorem, an independent
point, or a prediction of rank26. The initial seed and parent selection use
only its existing certified history; neither a public higher point nor a
retrospective visibility oracle enters the run. Raw preparation and receipts
are under `artifacts/local/elliptic-curves/curve52-v3-preparation-v1/` and
`curve52-short-v3-discovery-v1/`. The
[completed first pass](CURVE52_SHORT_V3_SEARCH_2026-09-09.md) certifies25→26 on
chart73, then rebuilds and completes27 further invocations without another
gain. Search110.055 seconds and independent replay117.552 seconds both finish
successfully. A precise continuation cursor retains the partially visited
centre at the100-invocation cutoff. Full-landscape replay is now a material
cost relative to the shorter search; this is not an end-to-end speedup
benchmark against the longer policies.

The next two short passes use disjoint pairwise-derived parent classes on
[curve90](CURVE90_PRODUCTIVE_V3_SEARCH_2026-09-09.md) and
[curve92](CURVE92_PRODUCTIVE_V3_SEARCH_2026-09-09.md). Both independently replay
with no gain or timeout in100 invocations. Search/replay times are183.010 /
355.775 seconds on90 and95.416 /53.632 seconds on92. Curve90's28-parent
full-landscape replay costs more than its short search. Before extending a
queued prefix, reuse its sealed verified landscape with exact basis/source
bindings rather than paying to reconstruct and verify it on every batch.
The [cached continuation wrapper](../cas/run_cached_seed_v3.py) now implements
this separately from the queue recorder. Its first run extends curve92 from
100 inherited receipts to a cumulative cap of200, with at most100 new searches.
Exact file and source bindings preserve the previous full CVP verification;
all point receipts and the extended cloud are checked again. Gain-triggered
new landscapes still require independent rational-CVP replay. The parent
directory remains unchanged. The first continuation completes100 new charts
with no gain,200 cumulative, and passes replay. Search takes84.431 seconds
and replay3.599 seconds, versus95.416 /53.632 seconds for the earlier first
pass. The batches search different chart sets; this is operational evidence
of avoided landscape work, not a controlled end-to-end speedup benchmark.
The inherited full CVP proof remains an explicit dependency of the result.
The next cached batch completes the remaining92 invocations with no gain,
292 cumulative, in67.260 search seconds and4.068 replay seconds. Both batches
pass their source/file guards and fresh point-cloud replay; the original
queue is resolved by a separate completion receipt. This validates repeated
cache reuse on the unchanged basis while preserving the first full CVP proof.

The [queue recorder](../cas/queue_verified_short_pass.py) now binds verified
terminal records and preserves the next coordinate policy even when the cap
cuts a centre between its two maps. A gain at the cap is explicitly marked
as requiring a rebuild. It neither launches the queue nor labels the suffix
as searched. Completed90/92 prefixes and the curve52 gain remain intact.

Prioritize fewer unproductive invocations over another broad CVP rewrite.
Keep both maps, provenance, independent gain proofs and replay. Reuse known
productive parent classes without restricting away successful generic norm12
or half-integral shells. Test any new cutoff with separately masked recovery
controls before relying on it. No old checkout, old-search repair, higher box
height or additional campaign is required for this cost assessment.

## Historical preparation

The next candidate is the certified rank-27 MW16 curve
[`new-20260906-90`, at -1867/270](MW16_RANK27_VISIBILITY_AND_TRANSLATION_CLASSES_2026-09-06.md).
Preparation has independently reproduced its existing 27-point certificate.
No new point search or rank improvement is reported here. Historical runners,
frozen sources and searches are not repaired or restarted by this work.

## Measured costs

The completed `residual-strict-03` curve-302 cascade has 385 charts over 13
epochs. Its retained epoch wall times sum to 3937.978 seconds, of which backend
point searches account for 271.019 seconds (6.88%). The remaining 3666.958
seconds include landscape construction, maps, admission, checkpointing and
within-search verification; this audit does not attribute them individually.
Independent final replay is outside those epoch totals.

The individual chart files total 10,548,928 bytes. Cumulative cloud snapshots
total 299,061,547 bytes because each snapshot repeats preceding charts.
These are file sizes, not measured disk traffic or time attributable to I/O.
A future runner should retain each chart once and bind incremental receipts
to its ordered hash chain, with independent whole-epoch verification.

On curve90, the complete finite rank certificate took 0.112 seconds;
certified-state construction 0.115 seconds; the 27-point height matrix
0.026 seconds; and finite coset fingerprints 0.731 seconds. These single-run
measurements do not support prioritizing those preparation components.

## Implemented optimization, future runs only

[`IntegerExactParity`](../cas/visibility_lattice_fast.py) inherits the original
rational LDL decomposition and Babai implementation. For each LDL level choose
a common denominator b_i for its shift coefficients. Write its quadratic
contribution as

    d_i / b_i^2 * (b_i*w_i + a_i)^2.

One positive common scale makes every weight integral. Recursive used norms,
remaining radii and branch bounds can then use integers exclusively. Integer
square roots still give exact closed-ellipsoid bounds. Sorting by
`abs(b_i*w_i+a_i)` preserves the reference traversal and tie order.
At every leaf the scaled decomposition is checked against the original
quadratic form. Node limits retain their existing meaning.

The fixed retained benchmark takes the first anchor at M18/M26/M30 and the
first, middle and last distinct refined entries: eight cases in total (the
M18 anchor has only two entries). All eight reproduce the retained norm,
complete minimizing-vector list, node count and proof record exactly, as
well as a fresh run of the reference implementation. Aggregate solver time
improves by **3.093x** on this sample. This is not an end-to-end speed claim.

Four focused tests also cover all parities of 16 deterministic positive
definite matrices in dimensions 1–4, independent small-box enumeration,
ties, zero, malformed inputs and node-budget exhaustion.

## Next release gate

Build the future parent adapter with an explicit generic dimension. MW16 at
M27 has 32*2^11 = 65,536 extension cosets, twice the M27/MW17 count; silently
reusing the original hard-coded 17 would be incorrect. Bind its actual generic
basis and parent bank. Use the integer solver for search and the original
rational solver for independent verification. Verify masked recovery of the
retained 27th direction before opening the full-M27 search for a new direction.
No hidden known point may enter prospective centre selection.

## Adapter and masked-control result

The future adapter now reproduces the complete retained302 M18 selection,
including all114 centres, exact CVP records and score-array hashes. That
calibration completed in3.881 seconds. The exact MW16 norm8/10 bank contains
20,010 classes, computed in12.248 seconds; its labels agree with the historical
floating census. Preparation reconstructs the generic sixteen through the
native atlas and exact scale6 transport, then recertifies both M26 and M27.

The masked curve90 M26 landscape scores32,768 extension cosets and selects1511
centres in89.579 seconds. Its component timings are83.818 seconds for CVP plus
centre arithmetic,4.405 seconds for maps and0.271 seconds for scoring. All
selected CVPs subsequently replay with the original rational solver; the
combined reference-CVP and pointwise visibility audit takes272.609 seconds.

**The initial norm8/10 masked visibility gate fails:** neither sign of the
withheld point is inside height125000 on any selected chart. The smallest
coordinate height is1,398,686. This is not absence of the entire exceptional
direction or a prediction that a constructor must fail. No point search is
launched from this null calibration.

There is a parent-specific coverage omission: MW16 also has half-integral
generic shells9.5 and11.5. The historical winning generic mask585 is in9.5.
It was excluded before specialized scoring by the literal norm8/10 transfer.
The next preparation audits the complete generic ellipsoid through scaled
norm23 and retains all upper shells at least8, rather than adding that single
known mask. This is retrospective calibration informed by the null result,
not a blind prospective selector. The next masked roster still needs its own
freeze and verification before a discovery run.

That upper-shell enumeration completes in23.652 seconds and represents all
65,535 nonzero generic parities. The upper-shell counts are19,918 at norm8,
9,895 at9.5,92 at10 and12 at11.5:29,917 upper classes in total. Thus the literal
norm8/10 adapter omitted9,907 upper classes. These are exact minima in the
declared generic parity quotient, not a continuous covering-radius theorem.

`future_point_admission.py` prepares complete finite quotient tables once per
curve and caches point columns across charts. Its retained302 control finds
the same ordered19-point basis on the same13th chart, in2.516 seconds including
a fresh standalone rank certificate. A failed column test remains UNKNOWN.
The prepared constructor consumes sealed centres in order and never reads
withheld coordinates; it remains unlaunched pending the visibility gate.

Evidence directories under `artifacts/local/elliptic-curves/`:

- `curve90-v3-preparation-v1/`: native basis, certificates and exact bank;
- `future-v3-302-calibration-v1/`: complete retained-selection comparison;
- `future-v3-curve90-masked-v1/`: sealed1511-centre landscape and evaluation;
- `curve90-upper-shell-audit-v1/`: separately bounded upper-shell enumeration;
- `v3-curve90-cost-audit-v2/cached-admission-control.json`: retained admission check.

Each supervised operation has a sibling supervision directory. The old
campaigns remain unchanged. These are preparation/control results, not a
new rank or conductor record.

## Expanded control and productive-anchor continuation

The separately frozen upper-shell adapter scores61,440 cosets at60 anchors
and selects2,828 centres. Its supervisor completes in169.522 seconds. An exact
post-selection coordinate audit compares both the factor-free and historical
quartic-minimizing maps, without point searches or changing the selected
centres. Neither policy exposes either sign of the withheld point at125000.
Their smallest coordinate heights are1,398,686 and699,343 respectively.
The expanded CVPs have not been independently replayed; the coordinate
identities and pointwise audit do not require a CVP-optimality assertion.

For the fixed historical winning centre, the factor-free map gives height
162,862, versus79,466 for the historical map. The respective single-map
timings are0.00312 and0.00413 seconds. Across the full2,828-centre comparison,
map construction totals6.886 and9.523 seconds respectively. This argues for
retaining both coordinate constructions in a future bounded experiment;
it is not a general complexity bound on quartic minimization.

Adding shells does not fix the anchor ranking: the exact generic class585
has9,312 strictly deeper classes in its own9.5 shell, with no norm tie.
The policy retains only16 anchors and25 canonical centres per shell, so585
is absent from both lanes. This is a retrospective diagnosis of one point,
not proof that deeper anchors cannot find other directions.

The next bank is therefore defined by **all historical certified gains** on
this fibre, rather than another tuned depth cutoff. A fresh complete finite
admission replay of the43-chart generic attempt and301-chart adaptive attempt
reproduces lower bounds26 and27, producing eight gain-chart certificates.
It freezes exactly these eight generic masks:

    585, 2681, 10131, 15841, 15888, 41207, 52060, 59125

Every representative and generic norm comes from the exact upper-shell bank.
This is explicitly a productive subset, not a complete generic shell bank.
Preparation takes14.559 supervised seconds with zero point searches. At M27
these eight anchors would score16,384 extension cosets, versus122,880 at60
anchors. The reduced count is not an end-to-end speed measurement or a
prediction of a new direction. All existing27 seed points may be used in the
next discovery attempt; none is itself a new-rank success.

Additional evidence under `artifacts/local/elliptic-curves/`:

- `future-v3-curve90-upper-masked-v1/`: selection, coordinate comparison and
  `winning-anchor-coverage.json`;
- `curve90-upper-shell-audit-v1/historical-winning-centre-*.json`: the two
  fixed-centre map comparisons;
- `curve90-productive-anchor-bank-v1/`: source/input protocol, eight fresh
  gain proofs and `anchor-bank.json`.

The expanded implementation is in `visibility_future_upper.py`; it preserves
the first `visibility_future.py` and both null controls. No constructor or
discovery search has yet been launched from the productive subset.

## Reproduction and evidence

Raw evidence is retained under
`artifacts/local/elliptic-curves/v3-curve90-cost-audit-v2/`: `report.json`,
`supervisor.json`, `worker.log`, and `integer-cvp-benchmark-v2.json`.
The preparation supervisor completed in 2.032 seconds under a 300-second,
2-GiB cap. The failed first audit remains in the sibling `v1` directory;
it used the wrong loader call signature and made no search invocations.
The first benchmark is retained but superseded by v2, which corrects its
status text from nine cases to the actual eight-case roster.

From the repository root, using fresh output paths:

```sh
sage -python research/elliptic-curves/cas/audit_v3_continuation_costs.py --output /tmp/v3-cost-audit.json
python3 research/elliptic-curves/cas/benchmark_visibility_lattice_fast.py --output /tmp/v3-cvp-benchmark.json
python3 -m unittest discover -s research/elliptic-curves/tests -p test_visibility_lattice_fast.py
```

Reports bind read inputs and sources; reruns preserve earlier output files.
The cost audit consumes completed local logs; it does not independently
replay the full historical search. No mathematical-status promotion is needed
for these implementation and timing results.
