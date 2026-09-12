# Constructor transfer: failed positive calibration, pilot closed

**Closed on 2026-09-12. No fresh fibre ran.** The controller reached
`STOPPED_COMMISSIONING_GATE`; both reference arms have terminated. The
[terminal receipt](../../artifacts/generated-results/elliptic-curves/fresh_constructor_transfer_v1/terminal-closure.json)
binds the retained outputs. The aggregate charged wall time is 543.228 seconds,
including failed launches. Wall time is not a process-tree CPU measurement.

| Stage | Retained result |
|---|---|
| Starting subgroup | Exact rank16 control certificate; cold field preparation passed |
| Class construction | 841 distinct atoms from 23,244,696 tested norm pairs; projected outside-S parity rank841, kernel dimension0 |
| Strict class and lift | Not reached; no new dependency was extracted |
| V3 reference | Rank22 search checkpoint after six gains; final verification failed at certificate equality, so the benchmark rank gain remains UNKNOWN |
| Fresh panel | Zero of eight inputs started; commissioning gate failed |

The [binary stop diagnostic](../../artifacts/generated-results/elliptic-curves/fresh_constructor_transfer_v1/reference-class-stop.json)
replays the parity statement from retained factorizations. It does not independently
recertify every factorization. The class arm cost321.670 charged wall seconds.
It tested a restricted empty-pool cold adapter: the long relation warm-up of the
successful reference constructor was not reproduced. This is a failed positive
calibration of that adapter, not a disproof of the earlier constructor or its
two soluble classes.

V3's final replay raised `AssertionError` at `transfer_v3.py:187`,
`assert fresh == old`. Its exact cause has not been certified. The saved
rank22 checkpoint is preserved, but no repaired or independently completed
benchmark result is claimed. The post-seal union routine fell back to the
initial16 when outputs were unverified; that number does not establish equality
with V3 or absence of complementary directions.

The initial import-package failure and subsequent missing native-source failure
remain charged. Source snapshots, all arithmetic chunks, failed receipts and
the unstarted eight-fibre roster are retained. The original preflight checker
still verifies the frozen selection and input identities only; it is not a
certificate of end-to-end transfer. The128 fresh point equations did not prove
sixteen independent points on each unstarted fibre.

**No further arithmetic-constructor experiment is scheduled.** Reopening needs
a new mathematical reason to expect success and a separately bounded proposal;
another calibration or larger bank alone is insufficient. Column7, fixed-word
continuation, the carrier bank and full class-group computations remain outside
the active search. The research priority is the
[rank32 objective](../README.md).

The complete prior design and launch narrative are preserved in the
[closure archive](../../archive/constructor-transfer-2026-09-12/FRESH_CONSTRUCTOR_TRANSFER.before.md.txt),
with [hashes and prior ledger entries](../../archive/constructor-transfer-2026-09-12/MANIFEST.json).
They are historical instructions, not an active queue. All existing scripts
and proof artifacts are retained for reproducibility.

To inspect the stopped run, from the repository root:

```sh
research/elliptic-curves/rank-jump/transfer-status.sh
```

The local checkpoint is `research/artifacts/local/elliptic-curves/fresh-constructor-transfer-v3/`.
Do not restart it to obtain a cleaner comparison ledger.
