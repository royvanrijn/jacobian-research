# Split-first det1092 search: implementation and handoff

Status: **experimental implementation; not a new theorem, rank result, or measured production speedup.** Reviewed against main `36d0b5a4315e67e6df817d6044a40c5febbe7be6`.

## Decision

Replace the det1092 **seed-discovery front end**, not the certified cascade and not the entire new-A1-parent programme. Initially reallocate two existing worker slots to this bounded pilot; do not add unbudgeted parallel workers. Neither active campaign has been stopped or modified by this change.

The current parent foundry already makes genuinely new A1 fibrations. Its checked-in report has six accepted parents, first panels of at most twelve slots, several certified +3 detections, no +5 detections, and explicit incomplete exposures. That is not enough evidence to prefer one family or declare that strategy exhausted. It remains a useful exploration lane.

The actionable theorem result is different: generic norm-ten traces now construct rational bisections by Euclidean arithmetic, without the 19-by-20 RR solve or polynomial factorization. Once constructed, incidence at a fibre is just an exact rational-square test. The closed fifteen-orbit panel at 302 is **not** a closure of the 40,917-orbit atlas, and its successful construction timing excluded the cost of obtaining the generic traces.

Sources:
- [Current foundry report](../../artifacts/generated-results/elliptic-curves/parent-foundry-v2/REPORT.md)
- [Euclidean formula, scope, and positive/dependent admission controls](DET1092_EUCLIDEAN_BISECTION_FORMULA_2026-09-09.md)
- [Existing exact specialization/cascade worker](../cas/parent_foundry_worker.py)
- [Existing finite admission gate](../cas/future_point_admission.py)

## Replacement pipeline

`generic norm-ten words -> cached exact conics -> cheap split sieve -> certified extra directions -> existing V3 cascade`

1. Freeze a generic-only projection of the reduced det1092 equation, its seventeen sections, their Gram matrix, and a target-independent hash-ordered window of rational norm-ten orbit words. Check each word's norm exactly. No saved conic, exceptional point, record parameter, or later search output is admitted as a construction input.
2. In bounded Sage subprocesses, construct the negative trace from the generic word. Independently repeat its group arithmetic in reversed order. The Sage-free Fraction kernel checks the trace equation, inverse and exact divisions, degree bounds, smoothness, and both polynomial map identities. A bad chart, timeout, or other construction failure remains unresolved, not excluded.
3. Scan a predeclared sequence of distinct signed Calkin-Wilf rational addresses against every successfully built conic. Exact integer square-root receipts distinguish nonsplit, ramified, and nonzero split cases. The two roots of one bisection are not counted as two directions.
4. On split fibres, retain the existing generic finite-independence gate. Admit candidates only when they increase the finite-column rank, then construct and replay the standalone rank certificate with both existing implementations. `UNKNOWN_FINITE_COLUMN_IN_SPAN` stays UNKNOWN: the current implementation does **not** substitute a halving-theorem conclusion for the finite certificate required by V3.
5. Only a certified gain starts the seeded cascade. Reuse `parent_foundry_worker.py` and its unchanged V3 search/full-cloud replay, with 64 point calls. Every fourth predeclared address also gets an unguided control, whether or not any conic splits. On those addresses the seeded/control comparison uses the same equation, parameter, bank index, height and point-call allowance. The seed construction/admission cost is additional and must be included in efficiency comparisons.

The implementation is deliberately scoped to det1092. It does not silently transfer the all-prime/halving theorems to the new MW16 A1 parents. It also does not claim a split always gives a new direction on a particular rational fibre.

## Files

- `../cas/euclidean_seed_sieve.py`: exact polynomial arithmetic, Euclidean construction and map replay, square receipts, deterministic addresses and orbit selection.
- `../cas/run_euclidean_seed_foundry.py`: immutable preparation, frozen-runtime dispatch, bounded parallel construction, exact admission, cascade/control bridge, resume and status reporting.
- `../cas/test_euclidean_seed_sieve.py`: arithmetic, independent symbolic oracle, integrity and controller integration tests.

Preparation never launches a process in the background. `run` is a foreground bounded campaign using the existing resource supervisor. It runs its controller and subprocesses from the frozen source snapshot, not a moving checkout. Completed receipts are reused on resume; resource-censored work is not silently retried or relabelled negative. Source, executable and completed-output changes fail closed. A `STOP` file prevents new jobs; already running bounded subprocesses drain.

## Search-agent handoff

First run the existing cold positive/dependent theorem replay and the new unit tests:

```sh
sage -python research/elliptic-curves/cas/verify_det1092_euclidean_seed_admission.sage
python -m unittest discover -s research/elliptic-curves/cas -p test_euclidean_seed_sieve.py -v
```

Run a small integration smoke campaign in a **new** folder. This exercises actual Sage imports, the cold trace builder, frozen-source dependencies and the unguided cascade path:

```sh
SMOKE=research/artifacts/local/elliptic-curves/euclidean-seed-smoke-v1
python research/elliptic-curves/cas/run_euclidean_seed_foundry.py prepare \
  --folder "$SMOKE" --workers 1 --orbits 2 --parameters 1 --control-every 1
python research/elliptic-curves/cas/run_euclidean_seed_foundry.py run --folder "$SMOKE"
```

Do not equate a process exit with mathematical success. Inspect `STATUS.json` and the retained worker logs. Require a completed control result with its replay certificate. Inspect construction failures individually: a chart or resource cap is not an atlas exclusion. If the two smoke traces are censored, test a larger independently predeclared window before interpreting the pilot.

Then reserve two slots from the existing allocation and prepare the prospective pilot:

```sh
PILOT=research/artifacts/local/elliptic-curves/euclidean-seed-foundry-v1
python research/elliptic-curves/cas/run_euclidean_seed_foundry.py prepare \
  --folder "$PILOT" --workers 2 --orbits 128 --parameters 48 --control-every 4
python research/elliptic-curves/cas/run_euclidean_seed_foundry.py run --folder "$PILOT"
python research/elliptic-curves/cas/run_euclidean_seed_foundry.py status --folder "$PILOT"
```

`--sage /absolute/path/to/sage` is supported during preparation. Each trace has a default 25-second cap, admission 900 seconds, cascade 7,200 seconds, and each supervised subprocess a 3-GiB RSS cap. `--trace-seconds` changes the declared trace cap before preparation. No oracle or model calls are made by the campaign.

Resume with the identical `run` command. For a new declared window, use a fresh folder and explicit `--orbit-offset` / `--parameter-offset`; the old manifest is immutable. Prepared plans cannot silently expand based on results. `touch "$PILOT/STOP"` stops further dispatch; remove that file only when deliberately resuming.

## What to measure before a larger pivot

Report attempted and successfully constructed orbits separately; construction failures retain their reason/log and supervisor receipt. Report square splits, certified new directions, UNKNOWN admission outcomes, final certified rank lower bounds, extra cascade gains beyond the seed rank, and wall time split among construction/admission/cascade/control. Charge the one-off trace construction cost rather than quoting only the warm square-test speed. Never aggregate the seeded and control evaluations of the same parameter as independent fibres.

The decision objective is certified high-rank output per total CPU/wall budget, not the number of split conics or easy +1 seeds. A frequent +1 with no subsequent amplification does not justify abandoning the A1 exploration lane. The finite admission gate can miss genuinely independent seeds; preserve those UNKNOWN cases for a separately certified halving/nonhalving bridge rather than declaring dependence. No general curve novelty or conductor claim is made by this runner; promote successful packets through the existing repository novelty/conductor intake.

## Validation performed for this change

23 tests pass under Python, including an independent SymPy reconstruction of the sixth-order cancellation, fifteen exact normal-form examples, randomized polynomial division, large rational square tests, denominator transport, ramification handling, immutable-output/source tamper rejection, and mocked-CAS controller tests for seeded/control dispatch, resume, timeouts and graceful stop. Python compilation and CLI parsing were also checked.

**Sage/GP are unavailable in the implementation environment.** The full Sage trace/admission/cascade integration has not been executed here, no prospective rank gain has been obtained, and no production speedup is asserted. The explicit smoke gate above is mandatory before reallocating the production search. Historical frozen checkers, mathematical status entries, and running-agent state were left untouched.
