# Rank32: find the next independent point

The active objective is an unconditional certified rank32 subgroup, or a search
improvement that increases the chance of finding one. The next experiment measures
CPU until one independently certified direction beyond a supplied subgroup
`M_r`. A class-constructor transfer, parameter theorem or class-group upper bound
is not a prerequisite. The
[constructor pilot is closed after failed positive calibration](../rank-jump/FRESH_CONSTRUCTOR_TRANSFER_2026-09-12.md).

**Prepared, not launched.** The worker-time ceiling is awaiting the user's choice.
No old constructor allowance is transferred to this experiment. The executable
accepts a ceiling of at most eight worker-hours; that is a maximum, never a
minimum spend. A failed calibration, verification or validation gate stops it.

## Frozen controls and comparison

Reuse [the public28 recovery control](BLIND_FACTOR_FREE_CONTROL_AND_PROSPECTIVE_EXPOSURE_2026-09-07.md),
[Curve302's complete V3 ladder](ADAPTIVE_HALF_LATTICE_V3_2026-09-07.md) and
[the alternative seeded ladder](CURVE302_SEEDED_V3_RESULTS_2026-09-08.md).
The local packet preserves all fourteen pre-acquisition Curve302 states M17
through M30. The first comparison uses these four cases:

| Starting subgroup | Next-direction control | Centre bank |
|---|---|---|
| Curve302 M27 | 27→28 | First256 centres in the retained V3 order |
| Curve302 M29 | 29→30 | First256 centres in the retained V3 order |
| Curve302 M30 | 30→31 | First256 centres in the retained V3 order |
| Public inventory188 M27 | 27→28 | Complete original49-centre bank |

The coordinator projects only the equation, the starting basis, certificate
primes and centre coefficients. The worker reconstructs and certifies the
starting subgroup before searching. Higher-rank point suffixes, winning chart
indices and prior map receipts are not worker inputs. The original bank selection
and the choice of controls are retrospective. This is a fixed-bank representation
benchmark; it does not remeasure cold V3 landscape construction or prove a
prospective centre-selection improvement.

| Policy | Point search preparation and box |
|---|---|
| `v3_dual_125k` | Existing lean full minimization and factor-free reduction; height125000 |
| `factor_free_125k` | Factor-free reduction alone; height125000 |
| `v3_dual_500k` | Same two reductions; height500000 |
| `factor_free_shears_125k` | Factor-free map followed by the fixed coordinate changes below; height125000 |

The five changes are the identity, `u+1`, `u-1`, `u/(u+1)` and
`u/(1-u)`. Binary quartics, quadratic ordinate terms and the complete map to the
original equation are transported exactly. Infinity and poles are included.
The existing projective-box key removes certified signed-permutation duplicates;
arbitrary PGL2 changes are not assumed to preserve a finite height box.
Translating a centre, alternate fibrations and new centre-selection policies
remain possible later comparisons, not automatic additions to this run.

Each development arm has300 wall seconds, including starting-rank checks,
map construction and point search, followed by at most90 seconds for an exact
fresh-process replay. Each map is separately limited to5 seconds and1GiB;
each PARI point call to10 seconds. A failed map provides no point-box coverage.
One CPU is pinned, numerical-library thread counts are one, and the worker tree
has a3GiB RSS cap. Linux subreaping records descendant CPU separately from wall
time, including failed and interrupted attempts. Snapshot creation and input
certification are retained as shared preparation costs. Historical landscape
costs remain separately identified; absent complete cold costs remain UNKNOWN.

The score orders policies by the number of verified next-direction successes,
then summed CPU with twice the arm-plus-replay allowance charged for every miss.
Lane order rotates between cases. A gain stops that arm at the first certified
new direction; a point already in the initial subgroup is not a success.
Unknown finite-column admissions are not dependence proofs. Exact point maps,
square identities and finite-group independence are replayed in a fresh process.
This shares the maintained arithmetic implementation; it is not a claim of a
second independent CAS implementation. Certificate comparison uses canonical
JSON digests, retaining the tuple/list serialization regression exposed by the
old pilot.

At least two development successes are required. Freeze the provisional winner,
then compare it with the baseline on the alternative recovered-strict03 M30
subgroup, with450 seconds per arm plus90 for replay. Require a verified gain;
a nonbaseline winner that is slower than a successful baseline stops at this
validation gate. These are different subgroups on the same Curve302, not
independent samples of curves. The winner and both seals are frozen before
any production point call. A single timing comparison is conditional evidence,
not a universal speed theorem.

## Follow-up roster and exposure

The strongest locally generated inputs available with certified prepared banks
are rank27. Inventory188's rank28 includes reproduction of a known public point;
it is a control, not a newly generated rank28 candidate. The first follow-up
roster uses four existing rank27 candidates in stable inventory-ID order:

| Inventory ID | Family | Parameter |
|---|---|---|
| new-20260906-40 | 074d9 | 2818/1535 |
| new-20260906-48 | 11952 | 2828/2015 |
| new-20260906-71 | 103b2 | 3726/881 |
| new-20260906-90 | A1/MW16 | -1867/270 |

These are follow-ups on previously searched curves, not rank-blind fresh-fibre
selection. Subsequent public registrations do not by themselves establish prior
public provenance. This roster makes no worldwide novelty claim. Other rank27
rows, including41,72,186 and `r17-panel-103b2-low-02`, remain in the inventory;
the selected four reuse existing exact prepared banks without new lattice work.

The [retained exposure](HIGH_RANK_SHORT_PASS_ROSTER_2026-09-09.md) and
[parent-span reassessment](PRODUCTIVE_PARENT_SPAN_REASSESSMENT_2026-09-09.md)
contain substantial completed misses. In particular, earlier short passes and
the curve48 complementary-parent continuations are not unspent budgets. The
builder records completed same-basis coordinate boxes from every available
sealed native pass and skips exact matches. A larger height or a distinct
coordinate box is new bounded coverage, not a new curve or rank direction by
itself. This deduplication is scoped to the bound local passes, not all historical
exposure under every possible basis and isomorphism.

Only the frozen winning policy is applied. Remaining allowance is divided
equally between the four inputs with at most two hours per input, and a reserved
replay allowance. Each stops after one verified new direction. There is no
automatic continuation to a second point, new bank, parameter sweep, or larger
campaign. Any gain must subsequently enter the inventory through its usual
independent certificate and provenance review; a search receipt alone does not
publish a new rank record.

## Run and follow

The immutable local packet is
`research/artifacts/local/elliptic-curves/next-direction-benchmark-v1/`.
Its `plan.json` binds source snapshots, projected inputs, starting-rank
certificates, the complete control ladder and retained exposure. Workers have
an application artifact-read restriction; no operating-system isolation is
claimed. Source snapshots include the transitive Python packages and native
sources, avoiding the old launch omissions.

From the repository root, after choosing the ceiling:

```sh
# Example ceiling: eight worker-hours; launch detaches and returns immediately.
python3 research/elliptic-curves/cas/run_next_direction_benchmark.py launch \
  --folder research/artifacts/local/elliptic-curves/next-direction-benchmark-v1 --hours 8

research/elliptic-curves/next-direction-status.sh
watch -n 30 research/elliptic-curves/next-direction-status.sh
research/elliptic-curves/next-direction-status.sh stop
```

Status reports actual controller liveness, the current arm, progress and log
path. Every arm retains its maps, returned points, certificates and CPU/wall
receipts. The supervisor kills owned descendants if the controller exits.
Host sleep pauses local computation; the launcher does not change power settings.
An interrupted run is retained without automatic restart or reset of its budget.

Implementation: [input projection](../cas/prepare_next_direction_benchmark.py),
[bounded worker and exact replay](../cas/next_direction_benchmark.py),
[detached runner](../cas/run_next_direction_benchmark.py),
[narrow transport regressions](../tests/test_next_direction_benchmark.py).
The four transport, box-equivalence and finite-certificate regressions pass.
No new next-direction result or rank32 certificate is asserted by this protocol.
