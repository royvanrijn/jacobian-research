# Curve 385 sparse quotient rank-32 primary campaign (2026-09-04)

## Outcome

The precommitted primary campaign is complete.  Starting from the exactly
certified curve-385 `M29`, it searched the natural quotient-basis shells of
weight one and weight two in the structured cylinder built from the 43 frozen
generic-deep classes.  Both exact discovered-group classifications returned
the unchanged `M29` basis:

```text
M29 -- natural weight 1 --> M29 -- natural weight 2 --> M29.
```

This is a bounded negative point-search result, not a rank upper bound.  It
does not prove that `M29` is saturated in `E(Q)`, and rank at least 32 remains
open.

## Complete accounting

| stage | planned charts | fresh searches | exact prior-chart skips | charts with finite points | finite point occurrences | exact result |
|---|---:|---:|---:|---:|---:|---|
| natural weight 1 | 516 | 394 | 122 | 111 | 296 | no group growth |
| natural weight 2 | 2,838 | 2,722 | 116 | 554 | 1,376 | no group growth |
| total | 3,354 | 3,116 | 238 | 665 | 1,672 | rank lower bound remains 29 |

Every one of the 3,116 fresh quartic searches completed at reduced-coordinate
height bound `100000`, a 15-second per-chart wall limit, and a 1 GB GP stack.
There were zero timeouts, zero PARI failures, and no retries.  The 238 skipped
charts had exactly the same rational base-point key as a chart already present
in the frozen source ledger or an earlier primary stage.

The searches produced 1,445 distinct points carrying a primary-campaign
source tag, of which 1,360 were new to the source ledger.  After each complete
stage, the exact classifier supplied an integral relation for every discovered
nonbasis point.  It found neither a new independent direction nor a
finite-index enlargement inside the discovered group.

## Evidence

The compact audited manifest is
[`curve385_sparse_quotient_rank32_primary_v1.json`](../../artifacts/generated-results/elliptic-curves/curve385_sparse_quotient_rank32_primary_v1.json),
with SHA-256

```text
ff1b0a2e8dd29b9a34a3b81cf8db3bed5350d12ee0dcfd9dd3936f226d255d61
```

The complete 27,761,376-byte JSON ledger is retained as the deterministic
compressed artifact
[`curve385_sparse_quotient_rank32_primary_ledger_v1.json.gz`](../../artifacts/generated-results/elliptic-curves/curve385_sparse_quotient_rank32_primary_ledger_v1.json.gz).
Its compressed SHA-256 is
`08a2e416255910f733ef98283332e3a60a947350646329e4e2045cbc08d802c0`;
the decompressed ledger SHA-256 is
`17600ab552c8c4c5184d8ec02c6743c475424998e2672a8eeacc3ee75df5b77d`.

The manifest auditor checks the frozen protocol and source hashes, exact
old-plus-new chart-key union, all per-stage count identities, every search
status, the provenance of every returned rational point, and the complete
exact-relation accounting.  Replay it with:

```bash
python3 elliptic-curves/cas/promote_curve385_sparse_quotient_rank32_primary.py --check

python3 -m unittest -v \
  elliptic-curves/tests/test_curve385_sparse_quotient_rank32_primary.py
```

The original checkpointed search command was:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/run_curve385_sparse_quotient_rank32_search.sage \
  --height-bound 100000 \
  --timeout-seconds 15 \
  --stack-bytes 1000000000 \
  --max-stage 2 \
  --max-lattice-states 4 \
  --checkpoint-every 10
```

## Post-freeze restart-budget repair

That completed v1 command and its no-growth ledger remain immutable.  Its
single `--max-lattice-states 4` allowance is not sufficient for future
rank-32 work: two same-rank finite-index enlargements would consume the slack
needed to search all three unit rank gains from 29 to 32.

Future runs use the source-pinned v2 amendment and runner.  They have no
combined lattice-state cap.  Exact rank-changing and saturation-only group
changes instead use independent limits of three and four respectively; an
event with a rank gain is charged only as rank-changing.  Certified rank at
least 32 is checked before a budget stop.  The frozen regression path
`saturation, saturation, 29 -> 30 -> 31 -> 32` therefore succeeds with restart
accounting `(3,2)`.  The independent counters are operational limits only and
do not turn a stop into a rank or saturation theorem.

Replay the amendment and inspect the v2 plan without opening a search outcome:

```bash
python3 elliptic-curves/cas/build_curve385_sparse_restart_budget.py --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/run_curve385_sparse_quotient_rank32_search_v2.sage \
  --plan-only \
  --max-stage 2

python3 -m unittest -v \
  elliptic-curves/tests/test_curve385_sparse_restart_budget.py
```
<!-- status-consumer: EC-K3-R17-CURVE385-INDEPENDENT-RESTART-BUDGETS 39cfce110e3e494f -->

## Next gate

The next precommitted stage is
`alternate-a-weight-at-most-2`, containing 3,268 new charts after exact
physical-word deduplication.  The protocol intentionally requires explicit
stage-limit escalation before running it.  The primary result supplies
evidence that the natural low-weight heuristic does not transfer naively from
the three-bit recovery cylinder to this twelve-bit quotient; it does not by
itself decide whether the alternate-basis stage is worthwhile.
