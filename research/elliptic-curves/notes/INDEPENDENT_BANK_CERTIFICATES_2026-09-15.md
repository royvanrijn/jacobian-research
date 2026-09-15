# Independent certificates for compatible bank selection

**The bank-cost engineering gate passes.** All 34 paired rebuilds certify
byte-identical exported banks and identical complete selections. The new
certificate path costs **175.300703 complete child CPU seconds**, against
**251.028573** for full reference reconstruction: **30.167% less CPU**.
No point search runs, and sustained later-gain performance remains unproved.
Mathematical authority is `EC-CANCELLATION-BANK-CERTIFICATES-20260915` in
[MATH_STATUS.json](../../MATH_STATUS.json).

The [early-generator comparison](EARLY_GENERATOR_EXPOSURE_CONTROL_2026-09-14.md)
already verified 28 actual later gains through newly admitted generators, but
failed both charged performance blocks. Its 34 successful rebuilds spent
228.592993 component CPU seconds: 67.158259 in production, 144.486169 in the
independent reference selector, and 3.194898 checking native anchor identities,
with the remainder in subgroup/input preparation and serialization. This
motivates checking the supplied selection without repeating the optimizer.
It supplies no counterfactual recovery or timing result for an optimized search.

## Charged result

The [report](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2/summary.json)
and [supervisor receipts](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2/supervision.json)
record all 68 successful workers on 19 retained curves. There are no missing,
substituted or retried successor workers. The maximum complete worker costs
10.468125 CPU seconds; total child CPU is 426.329276, below the frozen 600 cap.

| Measurement | Full reference rebuild | Independent certificate |
|---|---:|---:|
| Certified identical banks | 34 | 34 |
| Complete child CPU seconds | 251.028573 | 175.300703 |
| M19 complete CPU seconds | 92.555923 | 70.347761 |
| M20 complete CPU seconds | 158.472650 | 104.952942 |
| Bank construction/checking component | 224.954360 | 149.447797 |
| Producer component | 65.472092 | 65.599001 |
| Reference selector / independent verifier component | 143.475402 | 71.011555 |
| Fresh native subgroup check component | 17.133865 | 16.584986 |

The paired baseline/candidate CPU ratio is **1.431988**, with central 97.5%
whole-curve, within-family bootstrap interval **[1.406525, 1.455576]**.
M19 saves 23.994% and M20 saves 33.772%. Thus all predeclared engineering
conditions pass. The nearly unchanged producer cost locates the saving in
independent verification; it does not demonstrate a faster discovery algorithm.

The [independent accounting check](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2/audit.json)
passes all 68 subgroup receipt bindings, 34 byte-identical bank pairs, 3264
independent minimum packets and 6528 direct integer norm/parity checks. Primary
candidate workers exhaust 3666070 verification nodes, versus 4020833 recorded
producer nodes, and perform 7062 native word checks and 3264 map recipes.
The accounting independently reproduces the paired timing and bootstrap result;
it does not rerun native rank proofs, maps or ellipsoids.

Supervisor CPU is 0.399504 seconds. Charging the candidate all **recorded**
development measurements and the failed predecessor adds 29.288524 seconds,
reducing the saving to **18.4996%**. This sensitivity is below the 25% primary
gate. Three early commissioning receipts measure component CPU and omit
interpreter startup; they must not be described as fully metered development.
The complete paired-worker comparison and this separate setup-cost sensitivity
are both retained.

## Exact minimum verification

The [minimum verifier](../cas/verify_parity_minimum.py) accepts an integer
symmetric positive definite Gram matrix `G`, a parity vector `p`, a proposed
minimum squared norm `R`, and the complete proposed set of minimizers.
All numeric certificate fields must be integers. A nonpositive or asymmetric
metric, missing tie, shorter vector, malformed packet or exhausted node cap
fails the certificate.

An independent rational symmetric Schur elimination produces `G = U^t D U`,
where `U` is upper unitriangular and every diagonal entry of `D` is positive.
The checker verifies this identity exactly. Thus

\[
 v^tGv=\sum_i D_i\left(v_i+\sum_{j>i}U_{ij}v_j\right)^2.
\]

It enumerates from the last coordinate to the first, retaining the fixed
closed radius `R`. If later coordinates contribute `S` and the current shift
is `a/b` in lowest terms, the admissible integer coordinates satisfy

\[
 |bv_i+a|\le m,
 \qquad m=\left\lfloor\sqrt{(R-S)b^2/D_i}\right\rfloor.
\]

The integer interval is intersected with `v_i = p_i (mod 2)`. Every omitted
branch has a sum of nonnegative exact quadratic terms exceeding `R`. Every
visited leaf is checked by direct multiplication with `G`. A leaf below `R`
rejects the proposed optimum, and the leaves at `R` must equal the entire
claimed set. Consequently a pass proves minimality and all ties in this
**rounded integer metric**. It does not certify numerical canonical heights.

This uses neither the producer's dynamic incumbent and traversal order nor
its denominator-scaled arithmetic. Recorded producer node counts remain
execution metadata; their exact values are not independently reconstructed.
The fixed-radius node count cannot exceed a correct dynamic traversal whose
incumbent is always at least the final optimum. The node cap supplies a
computational failure boundary, never an arithmetic exclusion.

## Whole-bank verification

The [selection checker](../cas/verify_cancellation_basis_bank.py) reconstructs
the numerical Gram recipe, exact rounded Gram, integral LLL transport,
fingerprint columns, complete extension census, heuristic shortlist and
canonical parent ordering. It checks every saved score array. Each selected
coset receives the independent minimum proof above; every minimizing word
is transported through LLL and evaluated using native Sage group arithmetic.

The prescribed factor-free map and reduction recipe is reproduced, including
its quartic profile. The final shortlist, point deduplication and complete
centre ordering must match. The stored map's historical `centre` cache object
is not mistaken for geometric input: its actual word is evaluated independently,
and all other map fields must match the prescribed recipe. Only exact-word
and exact-point repetitions within this verification are cached.

The [live producer adapter](../cas/verify_cancellation_basis_bank_v2.py)
serializes native Python tuples through the declared JSON representation
before checking equality with the saved selection. It changes no arithmetic,
minimum, anchor, map, bank membership or order. The certificate binds the
exact selection bytes read before verification. The outer rebuild separately
certifies the supplied subgroup before using this selection.

## Fixed cost comparison

The [successor protocol](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2/protocol.json)
has SHA256
`39b78c47a40bc3376b1eb2effe2647c093e3e18152636c07e2ad302b1892e4d8`.
It binds 121 source files and **all 34** successful early-exposure transitions,
17 at rank19 and 17 at rank20, without missing targets in construction.
This is a retained development population, not a fresh search-validation cohort.

Both arms freshly certify the complete subgroup with independent native finite
group arithmetic. The baseline runs the unchanged producer, rational reference
selector and native anchor check. The candidate runs the unchanged producer
and new whole-selection certificate. Both compare their complete selection
and exported bank with the retained output **after** rebuilding; neither uses
that output to build its bank. The final exported bank must be byte-identical.

There are 68 sequential isolated workers, alternating arm order by transition.
Each has 30 soft and 35 hard process CPU seconds and 90 wall seconds. The
campaign reserves both hard limits before starting a pair and has a 600
full-child-CPU ceiling. All imports, input/source checks, subgroup verification,
production, independent verification, comparison and writing are charged.
Supervisor CPU and preliminary development measurements are reported separately.
An incomplete worker stops the comparison without substitution or retry.

The engineering gate requires at least **25%** aggregate full-child-CPU saving,
a central 97.5% paired bootstrap interval for baseline/candidate CPU with lower
endpoint above one, and positive savings in both rank strata. Bootstrap draws
resample entire curve groups within each family, keeping both epochs and both
arms together: 10000 draws, seed 20260914. This gate is solely about bank cost.

## Retained commissioning and failure

The [initial commissioning](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v1/development.json)
checks retained M19/M20 selections. The [17 positive and corruption checks](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v1/corruption-checks.json)
include omitted ties, a false larger optimum, malformed metrics, a node cap,
omitted cosets, a wrong word, changed map and score payloads with updated hashes,
and changed final ordering. The
[bound-selection successor](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v1/binding-checks.json)
additionally rejects objects that differ from the saved selection. Original
source snapshots and receipts remain unchanged.

The [first timing protocol](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v1/protocol.json)
stops after its first baseline worker. It independently rebuilt a byte-identical
bank, then its harness compared Python minimum tuples with JSON arrays and
rejected that in-memory comparison. Its **4.834009 full child CPU seconds**,
source seal, log and failure are retained. It remains a failed protocol.

The successor uses JSON equality and explicit saved-bank byte equality.
Before its timing seal, a separately metered
[live-producer adapter check](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2/adapter-check.json)
freshly certifies and rebuilds both ranks and rejects mismatched selection
objects. It costs **9.259756 full child CPU seconds**. No bank policy, missing
point input, roster, cap or acceptance threshold changes. Earlier development
receipts report component CPU, not complete interpreter startup costs; the
recorded-cost sensitivity must retain that distinction.

## Retention and replay boundary

The [repaired archive](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2/replay-bundle-v2.tar.gz)
contains **11266 files**, compressed to **27443649 bytes**, with SHA256
`cf81d523dc0f48feb0a48d4dd7210b5b0ad39b7ac955e9edba95c2d7facb8f9d`.
Its [manifest](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2/replay-manifest-v2.json)
and [empty-root check](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2/portable-check-v2.json)
bind every file, source/input guard and independent accounting result.
The latter costs 1.109343 full child CPU seconds plus 2.075808 packaging-check
component seconds. It repeats integer norm/parity and timing accounting,
not the primary ellipsoids, native rank proofs, maps or point searches.

The first archive and its failed portability receipt remain at their original
paths. They omitted the `ecsearch` package imported indirectly by a legacy
rank helper. The [dependency supplement](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2/dependency-audit-v2.json)
adds and hashes the existing import dependencies after execution; these are
not retrospectively labelled part of the 121-source primary seal. Both archive
versions, the failed first timing protocol and all source snapshots remain
retained. The [completion record](../../artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2/completion.json)
identifies the passing archive and accounting files.

After extraction, the bank-cost accounting can be rerun without rewriting an
existing certificate:

```sh
cd research
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 sage -python - <<'PY'
import sys
sys.path.insert(0, 'elliptic-curves/cas')
from verify_cancellation_basis_bank_cost import verify
verify(write_result=False)
PY
```

This is a check of already retained evidence. The old `seal`, `run`, `report`
and commissioning commands preserve their used output paths; do not restart
them as an implicit campaign replay.

## Research boundary

This work performs **zero point searches**. It cannot change the failed
49-versus-49 early-exposure outcome, establish later-gain throughput, or certify
that a direction is absent from an entire subgroup coset. A future search worker
must consume the new certificate schema honestly and include its complete cost.
Its later-gain comparison needs separately frozen, disjoint whole-j validation
blocks. Target-independent prediction and an established external implementation
comparison remain open.

Complete-cloud admission and exact subgroup enlargement remain required.
Exploration must retain distinct neighbour boxes beyond the old fitted
height-125000 range. Zero model support is not an arithmetic exclusion; this
bank-cost change introduces no fitted-radius pruning.
