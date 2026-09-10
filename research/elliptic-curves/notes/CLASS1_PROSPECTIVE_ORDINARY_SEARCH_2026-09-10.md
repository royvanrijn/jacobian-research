# Class1 ordinary prospective search

**The user has capped this commissioning run at64 ranked fibres and16 independent
controls. Scoring is unchanged; automatic window expansion is disabled.**
See the [conditional frame-breadth continuation](X1092_FRAME_BREADTH_COMMISSIONING_2026-09-10.md).
A strict-class selector is not a prerequisite. The certified class1 equation and generic basis remain
byte-frozen at `7c6ee40c46f5a1f3d1fc464b5685a0e4c77ed1f3347b67865d7c9f93b2acd862`.
The earlier [strict intake](X1092_CLASS1_FROZEN_ARITHMETIC_GATE_2026-09-10.md)
remains a separate construction lane. Its unavailable checker does not veto
ordinary parameter/point search. Class3 realization is conditionally authorized
after complete class1 commissioning without a certified ≥20 fibre.

## Fixed selection and cheap arithmetic

The active v2 window is **65,536** distinct signed Calkin--Wilf rational addresses,
generated from indices 65,536 through 131,071. The preceding window remains
the preserved v1 commissioning cohort described below. Address generation and the input
projection use no historical exceptional parameters, points, scores or ranks.
There is no record-based exclusion or targeting. The frozen runtime contains
code and the generic-only parent projection, not historical result packets.

At primes 5 through 997, coefficientwise common `p^4/p^6` scaling is removed
from the polynomial model before constructing the residue tables, including
the base-infinity chart. For every smooth reduction count the points and sum
the existing Mestre-style score
`round(10^12*(2-a_p)*log(p)/(p+1-a_p))`. Two table entries per prime are checked
by a separate scalar Legendre sum. Local singular reductions are recorded
and omitted from the sum; possible further parameter-specific minimization
is not attempted. This is a heuristic scheduling signal, not an exact rank
classifier or a fully minimized Nagao score at every parameter.

Each address also retains its exact discriminant-zero check, discriminant and
equation bit lengths, local singular-prime roster and smooth-prime count.
No conductor factorization is run; a discriminant is never labelled a conductor.
Only exact singularity prevents an elliptic fibre from passing specialization.

The first **8,192** indices in the fixed SHA256 order
`SHA256("class1-controls-v2/" + index)` form the score-independent control stream,
in that same hash order. This window has **4,153 positive and 4,039 negative**
controls. The other addresses are sorted by decreasing score, then equation
bit length, then original index. Dispatch interleaves **seven ranked fibres
and one control**. Control selection and order do not read scores. All addresses are
retained, including low scores and local-singularity flags. Consequently the
first window contains 8,192 controls and 57,344 ranked addresses.

The complete v1 commissioning scoring took approximately **3.55 CPU seconds**, including
its 0.30-second cold residue-table construction and recorded scoring overhead.
This measures triage only, not the much more expensive point search. The
residue table and every per-parameter score are retained for replay/audit.
The corrected v2 window took **3.70 CPU seconds** for the corresponding
triage work. Its sealed
[commissioning packet](../../artifacts/generated-results/elliptic-curves/class1_prospective_search_commissioning_v2.json)
binds the plan, source snapshot, all score rows and corrected queue.

## Rank-dependent exposure

| Certified lower bound | Cumulative point-call allowance |
|---|---:|
| 17–19 | 24 |
| 20–22 | 256 |
| 23–24 | 1,024 |
| 25–26 | 2,048 |
| 27–31 | 8,192 |

The controller runs two existing-style bounded worker slots. Each parameter
must first pass exact specialization and finite independence of its seventeen
generic sections. Failure remains UNKNOWN, not dependence or an upper bound.
The existing `parent_foundry_worker.py` then runs V3, varying generic parent
banks across batches and retaining all minimum representatives. Batches use
24, 64, 128 or 256 calls according to rank and stop at the remaining allowance.
The search height is 125,000; the existing per-map/per-point limits apply.
Each batch has a 7,200-second/3-GiB process-tree resource cap.

**Every batch reconciles the entire saved point cloud**, including below23,
and replays rank certificates with two finite implementations. Rank27 receives
the largest cumulative allowance and fresh banks through bounded continuations;
the controller does not grind it indefinitely after that allowance. A certified
rank32 packet ends that fibre's exposure. There is no exact-rank upper-bound,
worldwide novelty or conductor-record inference from a search endpoint.

## Detached operation and costs

Campaign directory:
`research/artifacts/local/elliptic-curves/class1-prospective-v2/`.
The controller PID is recorded in `controller.json` (detached launch PID
3024429). The launcher waits for v1 to drain before starting two worker slots.
It executes a source/hash-frozen runtime without model calls.
`STATUS.json` lists active indices and completed fibres; raw logs, all returned
clouds, bank proofs, batch certificates and terminal results remain underneath
the runtime's `ordinary-search/window-000/` directory.

Every parameter retains its score, control flag, direct scoring CPU, equal
share of cold table cost and scoring overhead, and search CPU. An isolated
subreaper job driver accounts for worker and descendant user/system CPU,
including bank construction, points, replay and full-cloud reconciliation.
Wall time and censorship remain separate. Open or abnormally interrupted
work has pending/UNKNOWN cost, not a fabricated zero. Preparation and global
controller overhead are reported separately; controller CPU includes scoring,
so those figures must not be double-counted.

The controller drains on a `STOP` file and refuses overlapping controllers.
It reuses completed receipts. A missing cost receipt after dispatch is an
interrupted UNKNOWN, not an automatic duplicate invocation. If a whole window
ends, the original controller would double the next disjoint address window up
to262,144 addresses. **That expansion policy is superseded by the64+16 cap.**
The frozen original runtime is retained; a separate capped controller reuses its
scored rows and certified worker without changing the scoring function or budgets.
There is no daily budget or AI decision loop.

Read-only reporting:

```sh
python3 research/elliptic-curves/cas/report_class1_prospective_search.py \
  --folder research/artifacts/local/elliptic-curves/class1-prospective-v2 \
  --output research/artifacts/local/elliptic-curves/class1-prospective-v2/REPORT.json
```

Compare certified threshold yield and CPU between ranked and control arms,
keeping unfinished/censored exposure visible. Adaptive escalation is identical
in both arms. Raw split counts, cheap +1 counts or a handful of successful
ranked fibres are not evidence that triage improves high-rank yield.

Validation before v2 launch: five policy regressions passed (thresholds, control
invariance, ranked ordering, prime roster, exact control fraction across both
signs), the source snapshot sealed, and
the residue tables passed their independent scalar checks. Both first
dispatched fibres passed the specialized generic-M17 gate, completed24 calls
each and replayed their full-cloud certificates. Both remain at lower bound17;
their total charged CPU was approximately47.5 and49.7 seconds. The controller
moved on to the next two fibres. These first two outcomes are not evidence
for or against triage effectiveness; no new record is asserted.

## Preserved commissioning correction

V1 launched with controls at address indices `7 mod 8`. The signed-address
enumerator assigns every odd index a negative rational, so that stride was
score-independent but sign-confounded. The audit caught this before the first
control point call; one control had already been dispatched and was allowed
to drain with the other bounded worker. A STOP marker prevented further v1
dispatch. Its scores, point packets, CPU costs and any paused continuations
remain intact under `class1-prospective-v1`, and its initial commissioning
snapshot remains historical. Do not combine that control cohort with v2 when
estimating triage efficacy.

V2 uses the next disjoint address window and the fixed hash rule above, with
no parameter or score from v1 selecting its candidates. This corrects both the
sign lock and the early-prefix height bias of address-ordered controls.
