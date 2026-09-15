# Basis-aware amplification with unfitted coordinate exploration

**The frozen performance gate fails.** Across24 new paired controls, rebuilding
recovers12 later-cloud directions in783.018782 CPU seconds, against13 in
706.898576 for the fixed bank. The later-direction rate ratio is0.833341,
with central97.5% paired interval **[0.685765, 0.928738]**. The candidate uses
10.77% more complete CPU and recovers one fewer direction. One interrupted
bank replay is retained and charged; it independently blocks promotion.

| Arm | Two-direction targets | All added directions | Later-cloud directions | Complete CPU seconds | Point calls |
|---|---:|---:|---:|---:|---:|
| Fixed bank | 14/24 | 34 | 13 | 706.898576 | 847 |
| Basis refresh | 13/24 | 33 | 12 | 783.018782 | 797 |

The [completed accounting](../../artifacts/generated-results/elliptic-curves/cancellation_basis_amplification_v1/summary.json)
includes every certified partial gain and both arms' first positive clouds,
including the one two-direction first cloud. Each arm has four zero-direction
cases. The candidate has seven one-direction cases to the baseline's six.
The all-direction rate ratio is0.876234. Charging all26.697880 shared preflight
CPU seconds to the candidate gives a later-direction ratio of0.805866.
These are arm-level counts on known controls, not67 distinct new discoveries.

The [mechanism audit](../../artifacts/generated-results/elliptic-curves/cancellation_basis_amplification_v1/mechanism-audit.json)
confirms **18 verified rank18 -> rank19 rebuilds**, exporting271 anchors that
use newly admitted generators. Only **7 of308 post-refresh calls** reach those
new anchors; all12 later gains arise on anchors already available in the
initial bank. The existing ordering spends nearly all subsequent calls on
that old bank. This rejects the specified policy and budget as an improvement;
the broader subgroup-amplification hypothesis remains open. Any next proposed
gate must explicitly budget early exposure of new anchors and handle the
observed Sage alarm at the bank boundary. No additional trial follows here.

All24 pairs have identical independently verified initial banks. Their
pre-first-gain point prefixes agree on23 pairs; the remaining zero-gain pair
stops at60 versus59 calls under the CPU cap. Across both arms there are1624
completed calls and20 explicit timeouts. There are1064 calls with exact
beyond-training coordinate witnesses, of which1053 complete. None of those
witnesses asserts a quartic square or a new elliptic direction.

This experiment tests the subgroup hypothesis left open by the
[failed fixed-bank cloud policy](CANCELLATION_CLOUD_SCHEDULER_2026-09-14.md).
The earlier 5.175586 CPU-second integration check established one compatible
rank18 bank, with no point search or performance claim. Its certificate and
the failed two-direction comparison remain unchanged.

## Frozen comparison

The [protocol](../../artifacts/generated-results/elliptic-curves/cancellation_basis_amplification_v1/protocol.json)
was sealed before any new control outcome. It selects **24 new whole-j CPU
groups**, four per retained R17 family, excluding all prior scheduler controls.
Each starts from an independently certified rank18 subgroup and has a separately
certified retained endpoint of rank at least20. Endpoints construct these
masked controls; the workers receive only the starting seed and parent bank.
These are retained-corpus controls, not externally pristine population data
or new rank discoveries. The rank18 scope keeps the extension calculation
near the preceding integration evidence; it does not repeat the failed
rank20 cohort.

Both arms reconstruct the same initial bank from the first16 generic anchors
of the retained compatible landscape and the complete rank18 seed. Every
bank is produced with integer CVP, compared against the independent rational
CVP solver, and checked by native Sage point arithmetic before use. Selection
uses the recorded rounded positive numerical height metric; the comparison
does not certify an exact canonical-height matrix or a full generic shell.

The two arms use an identical deterministic box order. Each anchor gets its
factor-free H125000 box followed by its retained prime-neighbour boxes, at
most two, at the same new-coordinate height. Then the next anchor opens.
There is **no residue-weight, radius, hazard or cloud-yield fit**. This ablation
measures the actual subgroup change without selecting another fitted policy.
It does not claim the box order is optimal or enumerate all possible neighbours.

The fixed-bank arm accumulates certified output while keeping its initial
search subgroup and bank. The basis-refresh arm reconciles the entire cloud,
independently certifies every admitted direction, and rebuilds immediately
when another direction is still required and CPU remains. The new landscape
uses the full enlarged basis. Old epoch maps, anchor indices and exposure
state are discarded. Only exact identities of completed boxes survive across
epochs. Directions beyond the stopping target are counted in full.

Finite-column admission uses the same places through500 in both arms.
An unresolved column in the known finite span remains UNKNOWN, not an exact
dependence. Every returned point is considered. Exact infinity witnesses
remain valid even when the finite backend times out; a timeout certifies no
completed finite prefix. Every rank increase receives both a portable proof
and independent Sage finite-group arithmetic before it can change the basis.

## Exploration outside the fitted range

The old empirical CDF stops at height125000. Completing its base box can
exhaust the fitted mass while leaving new old-coordinate addresses accessible
through a neighbour. The new policy has no model-support pruning. A distinct
neighbour remains eligible even when every old fitted marginal would be zero.
The only deduplication rule is a proved sufficient equivalence of **completed**
boxes: the same actual anchor and height, with coordinate matrices differing
by a projective signed coordinate permutation, including infinity.

A bounded exact probe records a primitive coordinate inside a proposed
neighbour box, beyond old H125000 and outside all completed boxes on that
anchor when it finds one. A failed probe is UNKNOWN and cannot remove the
job. The witness proves coordinate exposure only: it need not be a quartic
square or a rational elliptic point. The regression using
`[m:n] -> [3m:n]` exhibits this distinction after complete empirical-CDF
exhaustion. Model-support exhaustion is never an arithmetic exclusion.

## Accounting and decision rule

Each arm has **40 working CPU seconds**, including interpreter/imports,
initial construction, all rebuilds and rational reference replays, native
anchor checks, maps, full-cloud classification and independent rank proofs.
An atomic operation can overrun this launch allowance. Final exact
transcript/policy/rank replay is mandatory and also charged to the complete
isolated process CPU. Hard limits are100 CPU and150 wall seconds per arm.
Every failed or censored attempt remains in the denominator. Both arms start
from the same retained seed/parent interface; neither receives a free
initial or refreshed bank. Historical parent/seed construction is outside
that interface. Shared preflight is reported separately and charged entirely
to the candidate in a conservative sensitivity calculation.

The primary endpoint is **later-cloud directions per complete CPU**, where
later means after the first positive *full* cloud, not after its first point.
All cases, including zeros, remain in this rate. Promotion requires at least
10% improvement with a central97.5% paired family bootstrap interval entirely
above parity, using10000 draws and seed20260914. It also requires at least as
many two-direction target completions and uncapped total directions as the
baseline. Undefined reference rates or bootstrap samples make the gate
UNKNOWN; any preparation or infrastructure unknown blocks promotion.

A retained development cloud separately checks rank18 -> rank20 reconciliation
and an independently replayed replacement bank without point search. A runner
smoke check uses a previously tested group. The first development runner
recipe omitted the GP executable hash and failed before a point call; the
corrected recipe and the failed attempt are both retained. Neither development
attempt enters the 24-group performance comparison. Pre-seal development is
outside both arm meters; the failed wrapper's CPU was not measured. The
preflight-charged sensitivity is therefore not a total development-programme
cost estimate.

## Retained interruption and independent checks

All48 frozen arms ran once. Forty-seven completed ordinary final verification.
In the remaining candidate arm, a `cysignals.signals.AlarmInterrupt` escaped
the worker's `TimeoutError` handler while the rank19 rational-CVP reference
replay approached the working limit. Its rank19 point packet and59 earlier
point calls had already been sealed. No search used the interrupted bank.
The original failed supervisor, log and incomplete producer/reference files
remain unchanged; only unstarted arms continued under the original seal.

The [supplementary finisher](../cas/cancellation_basis_finish.py) independently
replays that literal point-search prefix on a temporary proof view, ending
before the unverified bank attempt. It neither completes the CVP nor retries
a point call. The original40.596087 CPU seconds and the separately metered
2.554188 seconds of final replay are both charged to this arm,43.150275 in
total. Its single certified direction is retained, its bank remains failed,
and the protocol cannot pass. This is explicit interruption accounting,
not a retroactively successful arm or a replacement source seal.

All67 arm-level added directions receive independent Sage rank certification.
The four focused exploration regressions pass, as do four retained cloud-policy
regressions. A further6.007090 CPU-second
[transition check](../../artifacts/generated-results/elliptic-curves/cancellation_basis_amplification_v1/transition-regression.json)
replays one positive arm and rejects two deliberately corrupted copies: one
omits a mandatory rebuild and the other drops the certified cloud from its
replacement basis. These supplementary tests add no point searches and are
outside the arm performance meters.

The mechanism diagnostic, interruption finisher and negative replay checks
are explicitly post-protocol supplements. Their source hashes are recorded
separately; the original109-source seal is unchanged. The bundle also retains
compact copies of the already verified endpoint packets, avoiding a dependency
on large historical campaign files for their rank proofs.

The [byte-checked replay bundle](../../artifacts/generated-results/elliptic-curves/cancellation_basis_amplification_v1/replay-bundle.tar.gz)
contains12301 files in18159391 compressed bytes. Its
[member manifest](../../artifacts/generated-results/elliptic-curves/cancellation_basis_amplification_v1/replay-manifest.json)
binds every retained file. The pre-execution protocol hash is
`ec5a5f84038166af87e075b03c7c5cb1defbce26bfe06d0393f1c1316807089a`.

## Reproduction and boundaries

The [entry point](../cas/cancellation_basis_amplification.py) separates roster,
preflight, development integration, sealing, execution, reporting and packaging.
The [worker](../cas/cancellation_basis_epoch.py) and
[independent verifier](../cas/verify_cancellation_basis_amplification.py) bind
every search to the correct epoch and replay every cloud and bank transition.
The [exploration module](../cas/cancellation_basis_exploration.py) contains no
oracle or fitted model access.

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
sage -python research/elliptic-curves/cas/cancellation_basis_amplification.py run
python3 research/elliptic-curves/cas/cancellation_basis_finish.py recover
python3 research/elliptic-curves/cas/cancellation_basis_finish.py report
python3 research/elliptic-curves/cas/audit_cancellation_basis_mechanism.py
sage -python research/elliptic-curves/cas/check_cancellation_basis_transitions.py
sage -python research/elliptic-curves/cas/cancellation_basis_amplification.py pack
sage -python -m pytest -q \
  research/elliptic-curves/tests/test_cancellation_basis_exploration.py
```

`run` resumes only unstarted arms of this exact sealed protocol. It does not
retry a failed arm or change an existing certificate's source hashes. A
completed run must not be restarted. Raw development and control records,
all calls and failed work are retained under
`artifacts/local/elliptic-curves/cancellation-basis-amplification-v1/`.
No fresh-fibre continuation or automatic budget increase is authorized by
this comparison. Exact rank, saturation, population-wide amplification and
rank32 remain outside its scope.
