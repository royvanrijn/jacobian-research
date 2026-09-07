# Shared arithmetic and search runtime

Engineering ledger for the technical-debt changes. Mathematical authority
remains [MATH_STATUS.json](../../MATH_STATUS.json) and its canonical proofs.
Existing certificates, their generation sources, and the frozen 100-field BNF
experiment are retained. The concurrent MW16→MW17/MW18 port was integrated last.

The production entry points are [`run_arithmetic_pipeline.py`](../cas/run_arithmetic_pipeline.py),
[`run_mw_search.py`](../cas/run_mw_search.py),
[`run_pointed_quartic_search.py`](../cas/run_pointed_quartic_search.py) and
[`run_surface_proof.py`](../cas/run_surface_proof.py). Their shared implementation
is [`research_runtime`](../cas/research_runtime/). Misses, timeouts, missing
hypotheses and ambiguous reduction columns remain UNKNOWN.

| Requested debt | Implementation and connected use |
|---|---|
| Repeated factorization and field setup | `ArithmeticContext` retains exact minimal-model transport and discriminant factors. `SageArithmetic` factors before minimalization, registers proved support, caches maximal orders and inverse `polredabs` maps, and transports the maximal basis without factoring a new defining polynomial. Relative descent, rank-28, exceptional probes and fixed-field arithmetic use these adapters. |
| Full BNF as the front door | `SubspaceDescent` accepts certified known squareclasses, intersects local matrices, builds norm covers, computes restricted Fisher CT and searches only its radical. `full_selmer` requires an explicit completeness/upper-bound purpose; its Simon backend consumes prepared BNF and local frames. The MW18 residual runner defaults to local features. |
| Late regulator obstruction | Section exports, direct brute-force/msolve and V4 campaigns query `function_field_gate_record` before constructing equations or compiling workers. The gate uses current theorem constraints and automatically replays available finite-reduction proof packets. Missing packets never imply a bound. |
| Product-specific regulator code | `SurfaceProofEngine` and `SurfaceProofRepository` accept exact surfaces, good reductions, fibre data and a proof verifier. The toric K3 adapter checks retained Frobenius output, the exact model, boundary factors, Hodge data, reciprocity and Weil roots. It has been checked on a singleton twist as well as the product certificates. |
| Repeated independence classification | `ReductionCache` retains finite quotient witnesses and per-point signatures. `IncrementalReductions` maintains a binary column basis; a new point evaluates new columns and requests extra primes only on ambiguity. Existing mod-2 helpers share this cache. |
| Reconstructed subgroup state | Immutable `MWState` holds the exact model/context, independent basis, labelled height Gram, Kummer classes, finite signatures, parity/CVP data and observations. Subspace searches, lazy searches and generic pointed jobs retain state transitions. Compatibility chart adapters accept historical generator lists while distinguishing them from a certified basis. |
| Nearest-first parity enumeration | `VoronoiIterator` lazily expands integer-child streams using exact LDL lower bounds. It emits cosets in increasing nearest-representative distance; its diversity window still samples an initially shallow pool. This is a different policy from selecting maximum-depth classes. Checkpoints resume after a node-budget stop, without a full parity table. The MW18 deep-centre experiment uses a separate finite census with exact maximum-stratum audits. The MW17 `--single-index` path uses the iterator; the frozen 43/301 census requires `--legacy-census-regression`. |
| MW16-only metric policy | `ChartPolicy` separates the centre Gram from the two-dimensional chart metric and its weight. The retained sweep covers sixteen R17/MW17 controls and five MW18 anchors, with two identities held out of policy ranking. See the measured result below. |
| Coupled normalization and enumeration | Raw search contexts need no minimal model, factorization or number field. Descent upgrades the context explicitly, including an integral labelled generator when needed. `RepresentationPipeline` times normalization, chart construction and enumeration separately. The concrete benchmark varies raw/minimal normalization, raw/metric coordinates and GMP/PARI enumeration. Fixed-field conic parameterization and quartic normalization are also separate choices. |
| Per-class local arithmetic | `LocalKummerBasis` prepares one local basis per place and evaluates squareclass matrices. Kernel/quotient operations are binary linear algebra. The fixed-field family and generic subspace backend use simultaneous intersections; individually inadmissible classes may combine to an admissible class. |
| Repeated 2-congruent field arithmetic | `TwoTorsionContext` keys the labelled algebra, independently of the curve. Orders, units/class data when known, completions, prime ideals, local squareclass bases and generator transports are shared. Isomorphism or `polredabs` equality alone never identifies labels. |
| Expensive verification | Local, norm-cover and CT discovery retains exact witnesses. Portable facts allow empty-cache replay. Pointed replay checks maps, square hits, exact points and finite independence without running a census. MW18 `--check` verifies all positive cover/section/quotient witnesses and the exhaustive negative census through retained polynomial group-law DAGs and chord identities; `--check --regenerate` explicitly repeats discovery. |
| Bespoke worker lifecycle | `supervisor.py` owns process groups and descendants, RSS/RLIMIT/cgroup options, PARI stack limits, TERM/KILL escalation, atomic checkpoints, retained inputs/logs and a parent-death watchdog. Active descent supervisors, pointed GMP workers, section solvers and MW17/MW18 workers use it. |
| Over-computing scheduling features | `fast_features` asks only for selected local polynomial factors, discriminant valuations and real signature, plus class-group 2-parts already known. It requests neither an order nor BNF. The frozen unconditional 100-field experiment is unchanged. |
| Finite-field outputs discarded | `FiniteFieldFacts` stores exact-model/prime/extension/label identities, traces and orders, quotient bases, signatures, surface fibre data, local maps and verified Frobenius packets. Family trace tables are shared with twist scoring. |
| Theorems not pruning searches | `PruningRegistry` binds typed scopes to current status and pinned proof witnesses. Six imported rules cover five rank-zero product twists and the fixed-field radical. Stale authority disables a rule. Search labels and heuristic height scores are never substituted for Kummer classes or rank theorems. |

## Reproducible controls

```sh
.venv/bin/python elliptic-curves/cas/run_mw_search.py \
  --request elliptic-curves/data/runtime_mw_search_control.json \
  --output artifacts/local/mw-runtime-example/search.json
.venv/bin/python elliptic-curves/cas/run_mw_search.py \
  --verify artifacts/local/mw-runtime-example/search.json \
  --output artifacts/local/mw-runtime-example/replayed.json

.venv/bin/python elliptic-curves/cas/run_arithmetic_pipeline.py \
  --request elliptic-curves/data/runtime_subspace_control.json \
  --output artifacts/local/subspace-example/discovery.json \
  --cache-dir artifacts/local/subspace-example/cache --wall-seconds 30
.venv/bin/python elliptic-curves/cas/run_arithmetic_pipeline.py \
  --verify artifacts/local/subspace-example/discovery.json \
  --output artifacts/local/subspace-example/replayed.json \
  --cache-dir artifacts/local/subspace-example/empty-cache --wall-seconds 30

.venv/bin/python elliptic-curves/cas/benchmark_search_representations.py \
  --request elliptic-curves/data/runtime_representation_control.json \
  --output artifacts/local/representation-example/benchmark.json

sage -python elkies-k3/scripts/certify_r17_extreme_anchored_mw18_covers.sage \
  --check --no-resume --replay-output artifacts/local/mw18-complete-replay.json
```

The tiny subspace control uses two independent point classes on
`y²=x³−7x+10`. All three nonzero covers and the restricted pairing replay
without conic solving, local-point search or BNF. The rational-model upgrade
control additionally checks a nonintegral raw generator and its integral
field presentation. The surface control is
[`runtime_surface_control.json`](../data/runtime_surface_control.json): its
two retained reductions give rank zero without relying on the theorem registry.

A lazy request fixes an initial basis, a positive definite **scoring** Gram,
`next_holes`, a node budget, chart height/time and optional extra reduction
primes. Supplying its `cvp_checkpoint` requests the next unseen holes under the
same state/metric binding. A basis change requires a new metric or an explicitly
retained old subspace; the implementation never invents heights for new points.

An arithmetic request defaults to `mode: features`. `mode: subspace` supplies
source points and their finite-independence places; an optional target has an
explicit labelled algebra and generator map. `search_masks` must belong to the
restricted radical. `mode: complete-selmer` explicitly requests a full upper
bound. Retained complete-Selmer integrity is labelled separately from an
independent upper-bound proof replay.

## Metric transfer result

The [MW18 deep-centre comparison](MW18_DEEP_CENTRE_CALIBRATION_2026-09-05.md)
separately tests centre selection with the exact generic height Gram, forty
charts per anchor presentation, and height 100,000. The older trial below used
initially shallow centres and an identity scoring Gram on MW18. Its null
MW18 result does not establish backend insensitivity or test the successful
maximum-depth policy.

All 600 new boxes complete: nearest-first certifies no additional directions,
deterministic deepest selection certifies 22, and diverse deep selection 21,
across five anchor presentations on four distinct curves. These improve on the
shallow baseline but fail the frozen 35/50, minimum-five-per-anchor gate, so the
balanced 27-candidate roster remains unrun.

The [retained summary](../../artifacts/generated-results/elliptic-curves/runtime_chart_policy_sweep_v1.json)
and [portable witness bundle](../../artifacts/generated-results/elliptic-curves/runtime_chart_policy_sweep_v1.zip)
contain 84 cells and 1,008 completed charts: 21 historical control fibres,
four policies, twelve lazy centres each, height 10,000 and 0.5 seconds per chart.
Nineteen fibres enter ranking; two preselected identities are held out. Policy
inputs omit public complements and missing-point outcomes. R17 uses its supplied
generic Gram; MW18 uses a declared Euclidean scoring Gram, not a height claim.

Weights **1/16, 1 and 16 each certify four additional directions on one R17
control**, with different point representatives. Raw coordinates certify none.
Other controls, including the MW18 anchors and held-out fibres, yield no new
certified directions at these budgets. The timing differences between the three
weights are too small to establish a unique best policy. Weight 16 remains
configurable; this result does not prove its MW18 optimality.

A preliminary diagnostic omitted extra-prime escalation. It was stopped, and
the corrected admission allowance was frozen before rerunning the unchanged
panel, centres and metric choices. Both diagnostic and corrected logs are kept
locally. This is historical control calibration, not the frozen scientific BNF
experiment and not a prospective population claim.

```sh
.venv/bin/python elliptic-curves/cas/retain_chart_policy_sweep.py \
  --verify-bundle artifacts/generated-results/elliptic-curves/runtime_chart_policy_sweep_v1.zip \
  --output artifacts/local/portable-policy-example
```

That command reconstructs retained facts into an empty cache and checks every
point and state transition without CVP generation or point enumeration. A fresh
protocol can be frozen with `calibrate_chart_policy.py --freeze PATH`; running
or regenerating a calibration is a separate command.

## Verification boundary

The integration regression passed 101 tests across arithmetic, local/CT,
surface proofs, constraints, CVP, worker lifecycle and the pointed adapters.
Additional controls replay the eight normalization/backend combinations, the
MW17 lazy entry, the large MW18 local-feature entry and all eight retained
refreshed MW18 covers. The calibration bundle has replayed from an empty cache. Five additional
group-law/census tamper tests pass, including replay with Sage elliptic
addition disabled.

A hash is an integrity check, not a mathematical proof. Retained CAS Frobenius
and full-Selmer computations keep their declared trust boundary. The MW18
checker additionally checks the complete negative census: 312,952 finite
nonresidue obstructions and 78,240 global validation entries, with eight
positive covers checked by its companion exact replayer. The retained
[proof bundle](../../artifacts/generated-results/elliptic-curves/runtime_mw18_census_witnesses_v1.zip)
contains 129 batches across 38 chart/prime pairs. The negative census replay
took 36.3 seconds from an empty arithmetic cache on the development host.
The [complete replay record](../../artifacts/generated-results/elliptic-curves/runtime_mw18_complete_replay_v1.json)
also retains the positive MW states and finite-field facts. Each group-law node is bound
to its exact equation-basis word; chord frames and nonresidues are checked
by polynomial identities. No elliptic additions or modular census run during
negative replay. Positive point relations still use exact elliptic addition.
The original census and its source remain available for explicit regeneration.
CVP is lazy, not polynomial-time in the worst case. Unsupported surface proof
hypotheses, uncached full class data and unresolved finite columns remain
UNKNOWN. Historical fixed-centre/native-CAS experiments remain reproducible
comparators rather than the production scheduling entry points.

The census witness backfill is independently reproducible and resumable:

```sh
sage -python elkies-k3/scripts/retain_r17_mw18_census_witnesses.sage \
  --all-directory artifacts/local/mw18-census-backfill --workers 2 \
  --wall-seconds 120 --rss-bytes 1073741824
sage -python elkies-k3/scripts/retain_r17_mw18_census_witnesses.sage \
  --pack-directory artifacts/local/mw18-census-backfill \
  --bundle artifacts/local/mw18-census-backfill/witnesses.zip
sage -python elkies-k3/scripts/retain_r17_mw18_census_witnesses.sage \
  --verify --bundle artifacts/local/mw18-census-backfill/witnesses.zip \
  --output artifacts/local/mw18-census-backfill/replayed.json
```

Discovery uses a verified transport to the sparse published trace basis.
Replay checks the original equation-basis labels directly. Packaging retains
its implementation sources for provenance, while original certificate input
hashes remain unchanged. A missing batch, altered word, invalid chord,
nonresidue failure or incomplete positive/negative partition fails closed.
