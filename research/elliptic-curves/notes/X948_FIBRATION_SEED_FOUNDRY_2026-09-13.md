# Eight new X948 fibrations: separate seed and amplification exposure

Eight new A1/MW16 fibrations passed exact admission. Each has a rational
equation, sixteen explicit sections generating the saturated generic group,
and a separating invariant proving inequivalence to the other admissions and
the frozen roster of eight previously searched X948 A1 fibrations. Generic
rank is exactly sixteen. Distinct abstract `O(NS)` orbits of the embedded `U`
are **not** claimed; the actual integral embeddings are retained.

Both positive detector controls passed using only equations and generic
sections. The detached first wave **stopped cleanly at its budget ceiling**:
9.923 charged elapsed hours and4.052 metered process-tree CPU hours. Stage1
is incomplete; Stage2 did not start. Mathematical admission is complete;
amplification and any useful advantage over the baseline remain unestablished.
Rank at least32 remains open.
The separately authorized short follow-up has also finished:128 further calls
on the strongest rank21 seed produced no new certified direction. It is one
targeted exposure, not the missing matched Stage2.
The [exported snapshot](../../artifacts/generated-results/elliptic-curves/x948_seed_foundry_v1/summary.json)
is dated data, not evidence that a process is currently running. Use the status
command below for that distinction.

## Budget-stop result

The [stopped-run receipt](../../artifacts/generated-results/elliptic-curves/x948_seed_foundry_v1/budget-stop.json)
retains554 terminal address jobs out of1,152 planned. The eight new fibrations
account for493 jobs:371 certified seeds,78 completed no-seed exposures, and44
unresolved starting-subgroup certificates. The baseline has61 jobs:49 seeds,
seven completed no-seed exposures and five unresolved starting certificates.
Thus505 jobs reached independently certified point-search outcomes;49 stopped
before point search because the finite independence gate could not certify all
sixteen starting points. Those failures are not proofs of dependence or rank
loss. No point-search or map timeout occurred among the505 certified outcomes.

The strongest seed certifies rank at least21 on priority1510 at `t=-16/15`,
from its first point-search call. Its five extra directions arrived in one
cloud; this is not a later amplification event. Packet hashes agree with the
retained independent-verification receipt.

Every fibration finished its first61 addresses; five new fibrations also
finished address62. The [matched61-address table](../../artifacts/generated-results/elliptic-curves/x948_seed_foundry_v1/matched-61-summary.csv)
compares only the common prefix. New-fibration seed counts range from39 to50,
versus49 for the baseline. Their metered CPU per seed ranges from29.66 to42.70
seconds, versus30.87 for the baseline. The best observed CPU enrichment is
1.041, below the predeclared1.5 threshold. This incomplete prefix shows good
seed production on new inputs, with no material advantage demonstrated over
the retained baseline. It does not evaluate the full-panel decision rule.

The remaining598 address slots are unrun. Neither Stage1's completion seal nor
the amplification roster exists. No continuation, parameter replacement or
automatic restart was launched at the stop.

## Separately authorized short rank21 amplification

After the first wave stopped, the user authorized a short amplification run.
The target is priority1510 at `t=-16/15`, selected as the strongest fresh seed.
This is a targeted follow-up, separate from the original matched Stage2. Its
selection cannot estimate an amplification rate or decide the fibration comparison.

The [thin launcher](../cas/run_x948_rank21_amplification.py) freezes the retained
worker and its4,064 source files, the exact parent and the complete21-point seed.
The new input replay independently verifies all21 points with two finite-group
implementations before any search. It passed. The starting packet SHA256 is
`1149d043b38174742d15cf23397fa15a2e5663bf050d55c517419c41d6183f0d`.

One worker uses one fixed logical CPU and one arithmetic thread, with3 GiB RSS.
It receives128 further dual-map point calls at height125000,10 seconds per point
worker and5 seconds per map. Search and its internal exact replay have900 seconds;
initial and final independent verification each have60 seconds, inside a1,200-second
aggregate ceiling. Preparation, failed attempts and verification are metered.
The new generic banks8..15 have exactly disjoint parity proposals from the old
bank0; the full rank21 basis receives fresh search landscapes. No old landscape
is reused across that basis change. Maps, full point returns and partial receipts
are checkpointed. A timeout retains the certified input bound and labels unreplayed
online gains as pending. It does not turn a bounded miss into an upper bound.

The detached run launched on2026-09-13 and completed all128 calls with no point
or map timeout. Full cloud reconciliation and separate independent output replay
passed. The certified subgroup stayed at rank21: no new direction was found.
Complete charged elapsed time was196.131 seconds; metered process-tree CPU was
173.998 seconds. The [sealed result](../../artifacts/local/elliptic-curves/x948-rank21-amplification-v1/seal.json)
binds the output packet, result and frozen plan. Its packet SHA256 is
`b57f86a7596384abb4a8d5293093d395949e42501cae8ad349edcfd8452f6af5`.

This completed bounded miss neither establishes exact rank21 nor estimates the
new fibrations' amplification rates. The matched Phase2 remains unrun. The
controller ended with `FINISHED_NO_AUTOMATIC_CONTINUATION`; the status command
checks the process start token as well as the retained result.

```sh
research/artifacts/local/elliptic-curves/x948-rank21-amplification-v1/status.sh
python3 research/elliptic-curves/cas/run_x948_rank21_amplification.py stop
```

The run folder contains `plan.json`, `input-verified.json`, `state.json`, worker
logs, the frozen runtime and `seal.json`, issued after the separate replay passed.
There is no automatic retry, additional parameter or budget increase.

## Admission and exact proof

The [complete norm-eight atlas](../../elkies-k3/ICARM_A1_MW16_ATLAS_2026-09-04.md)
contains63,917 translation classes, of which1,266 are genuine A1/MW16
candidates. The existing [residual-chord constructor](../cas/parent_foundry_geometry.sage)
supplies the equations and maps. No shell enumeration, target matching or
specialization scoring was repeated.

Previously searched priorities414,488,1070 are excluded before ordering. From
the complete retained A1 table, take the cheapest64 remaining trace words by
the stored group-addition bound, support, largest coefficient, coefficient
L1 norm and literal coordinates. A frozen24-entry order greedily maximizes
minimum source-lattice distance modulo sign, with construction cost breaking
ties. This trace-distance proxy supplies geometric variety without claiming
an orbit classification. Accept the first eight successful, inequivalent
constructions with compact coefficient size at most256 bits.

| New priority | Compact coefficient bits | Generic height determinant |
|---:|---:|---:|
|1191|144|474|
|1510|170|474|
|1792|143|474|
|1615|145|474|
|1956|142|474|
|4380|150|474|
|2021|142|474|
|2330|141|474|

The [roster](../../artifacts/generated-results/elliptic-curves/x948_seed_foundry_v1/roster.json)
links the complete equations, sections and admission packets. All eight first
attempts succeeded; ordering and construction took about20 seconds altogether.
The exact resource receipts remain in the exposure-state snapshot.

Write the saturated source lattice as `NS=U+MW17(-1)`, of determinant948.
For each norm-eight trace word `w`, the new fibre is `F=O_old+P_w`.
The compiler supplies the physical zero `O_new`, the nonidentity I2 component
`C=O_old`, and sixteen source-section classes that become new sections.
The admission verifier checks that

```
F, O_new, C, S_1, ..., S_16
```

is an integral19-by-19 change of basis with determinant±1. It also verifies
`F²=0`, `O_new²=C²=-2`, `F.O_new=1`, `F.C=O_new.C=0` and `F.S_i=1`.
Thus these classes generate the full source NS, not just a rational span.
Quotienting by the trivial lattice supplies the full MW group. The single
semistable I2 fibre allows at most a1/2 height correction: a nonzero section
has height `4+2(P.O)-correction >= 7/2`, so there is no geometric torsion.
The displayed basis is saturated and its height determinant is948/2=474.

Exact rational-function checks verify all sections and coordinate transports.
The raw coefficient degrees are8 and12; the finite discriminant is squarefree
of degree22 and coprime to A. Hence there are22 I1 fibres and one I2 fibre at
infinity. The source geometric Picard rank19 then gives the upper bound
`19-2-1=16`, matching sixteen rational sections.

An isomorphism of these fibrations must preserve the unique I2 fibre. In the
raw charts its base change is affine. The existing exact, depressed finite
discriminant weight ratios are affine invariants. Different vectors prove
inequivalence even over the algebraic closure. Equal vectors are conservatively
excluded, without claiming that equality proves equivalence. The nine old
presentations reproduce five keys; the three already searched new A1 families
add three more. Every admission has a new key.

The symbolic verifier shares the construction implementation. No independent
CAS reconstruction or formal proof assistant replay is claimed. The point
search independently checks specialized rank with two finite-group
implementations; on-curve tests alone never admit the starting subgroup.

## Frozen experiment

The [plan](../../artifacts/generated-results/elliptic-curves/x948_seed_foundry_v1/plan.json)
and source snapshot were frozen before the first new-fibration point call.
Allowlisted inputs contain generic equations, bases, lattice words and old
fibration invariants. Historical exceptional points and score tables do not
enter construction or search. This is input separation, not an operating-system
sandbox. Existing ranks and hits only select the separate commissioning controls.

The two controls are old MW16-01 at7/8 and previously searched X948 priority414
at1/8. Both factor-free and dual maps receive the same32-call allowance from
the generic sixteen. Both pass; factor-free requires1 and2 point calls, while
dual requires1 and1. Complete metered CPU totals are47.03 and37.81 seconds.
The frozen rule selects dual. These two controls calibrate execution and policy
choice; they do not establish a universal speed advantage or enter production
denominators.

Stage1 uses the same128 reduced, nonzero rational addresses in the height band
9..32 for every new fibration and old MW16-01 as a matched baseline. Addresses
are ordered by the declared SHA256 seed in the plan. The total roster is
1,024 new-fibration slots plus128 baseline slots. Jobs are interleaved by
address and fibration. Each receives32 actual point calls at height125000,
with10 seconds per call and5 seconds per coordinate-map worker. Searching
stops at the first gaining cloud; its full return is reconciled and certified
without searching another epoch. A failed generic rank certificate is UNKNOWN,
not a rank16 outcome, and does not cause replacement of the address.

Only after every Stage1 job finishes is its ledger sealed and Stage2's roster
written. The first sixteen certified seeds per fibration in address order,
or all seeds if fewer, receive128 further point calls each. Each starts with
the whole sealed Stage1 cloud. Selection never uses its rank or first-jump
size. Report another acquiring call, reaching23 through32, and gains occurring
after an earlier call already reached23 separately. Reaching23 in the first
cloud is recorded as such, not as a later deep amplification event.

One worker uses one fixed logical CPU and one arithmetic thread, with a3-GiB
process-tree RSS cap. The full first wave has a ten-hour charged elapsed
ceiling, including preparation, admission, calibration, unsuccessful attempts
and verification. Construction has a30-minute aggregate gate and180 seconds
per candidate; seed jobs have300 seconds and continuation jobs900 seconds,
with a separate60-second final replay cap inside the aggregate ceiling.
The controller does not start a job unless its declared bound fits. It drains
the current bounded job on a stop request. Unrun slots, failed preparations
and censored searches remain distinguishable. The ceiling may leave either
stage incomplete; it is not a promise to finish1,152 slots and all continuations.
CPU receipts measure reaped worker process trees and stage supervision.
Controller hashing and ledger work between stages are included in charged
elapsed time but are not separately allocated to per-fibre CPU. Both arms use
the same accounting scope; these ratios are not complete host-cost claims.

The practical seed gate is at least1.5 times the baseline's seeds per complete
CPU cost, with at least eight seeds in both arms and at most10% unresolved
slots. Construction/calibration costs are also reported, with common setup
allocated equally in the cold-cost comparison. For amplification, compare
the same first `min(n_new,n_baseline,16)` seeds; require at least eight fully
exposed in each arm and at least three extra rank23 hits with twice the baseline
count. Higher thresholds are descriptive secondary endpoints. These are
screening thresholds, not success-rate theorems. No winner automatically starts
another campaign: scaling requires a separately frozen fresh address interval.

The measured object is the fibration together with its compact chart and
detector exposure. Equal address boxes are not intrinsic under PGL2. A wide
range of observed rates would be a useful phenotype; similar rates would not
by themselves prove that the surface is the decisive variable.

## Follow the detached run

From the repository root:

```sh
research/artifacts/local/elliptic-curves/x948-seed-foundry-v1/status.sh
python3 research/elliptic-curves/cas/run_x948_seed_foundry.py stop
python3 research/elliptic-curves/cas/report_x948_seed_foundry.py
```

The status command checks the controller PID and its process-start token.
`state.json` records active jobs, complete process-tree CPU, elapsed time and
all terminal exposure records. Worker logs and immutable maps, clouds and
point certificates live under `runtime/research/jobs/` in the local run folder.
There is no automatic restart or continuing guardian.

To refresh the exported reporting snapshot without starting computation:

```sh
python3 research/elliptic-curves/cas/report_x948_seed_foundry.py \
  --export research/artifacts/generated-results/elliptic-curves/x948_seed_foundry_v1
```

To replay one exported admission without constructing another parent:

```sh
sage -python research/elliptic-curves/cas/x948_seed_foundry_geometry.sage verify \
  --folder research/artifacts/generated-results/elliptic-curves/x948_seed_foundry_v1 \
  --index 0
python3 -m unittest discover -s research/elliptic-curves/tests \
  -p test_x948_seed_foundry_reporting.py
```

The four endpoint regressions separate failed preparation, control rows,
same-call clouds and later deep gains. This branch reuses the successful point
finder and equation compiler. It does not reopen fixed-word continuation,
marked carriers, class-group construction, X1092 classification or deeper
search of the old stalled27/28 curves.
