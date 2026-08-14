# Reproducing the elliptic-curve programme

Run commands from the repository root.  The exact lattice code uses only the
Python standard library.  The curve replay additionally requires PARI/GP;
version 2.15.4 produced the pinned manifest.

## Unit tests

```bash
python3 -m unittest discover -s elliptic-curves/tests -v
```

These include simple and singular Hensel lifts, generalized CRT, skew Gauss
reduction, the primitive-vector cancellation trap, complete asymmetric height
boxes, exact Fermigier quartic-to-Weierstrass identities, and finite-reduction
independence certificates.

## Benchmark arithmetic

```bash
python3 elliptic-curves/scripts/verify_benchmarks.py
```

This recomputes the exact conductor of the Fermigier E22 benchmark and the
literal integer cutoff.  It also pins the unresolved source normalization:
the printed shift `19754/39` gives a different model and conductor, while the
doubled shift and canonical adapter reproduce E22 exactly.  It does not
itself reproduce point independence; the exact lower-bound checker below does.
No command here supplies an unconditional rank upper bound.

Cross-check both family metadata files against the executable equations and
the calibration family's rank-two specialization with:

```bash
python3 elliptic-curves/scripts/verify_family_data.py
```

## CRT--lattice calibration

Replay the checked-in manifest:

```bash
python3 elliptic-curves/scripts/verify_crt_lattice_calibration.py
```

Generate an unpinned fresh copy in the ignored local cache:

```bash
python3 elliptic-curves/scripts/run_crt_lattice_calibration.py \
  --output artifacts/local/elliptic-curves/crt_lattice_calibration.json
```

The generator refuses to overwrite an existing file.  The pinned artifact was
created with:

```bash
python3 elliptic-curves/scripts/run_crt_lattice_calibration.py \
  --output artifacts/generated-results/elliptic-curves/crt_lattice_calibration_v1.json
```

Its SHA-256 is
`eb1543031e68026042c921ee2b93e765070b65340b8129b74f0629a9b3d5c8fa`.

## Fermigier high-family CRT seed

Replay the checked-in three-prime seed:

```bash
python3 elliptic-curves/scripts/verify_fermigier_crt_seed.py
```

Generate a fresh copy in the ignored cache:

```bash
python3 elliptic-curves/scripts/run_fermigier_crt_seed.py \
  --output artifacts/local/elliptic-curves/fermigier_crt_seed.json
```

The pinned command uses output
`artifacts/generated-results/elliptic-curves/fermigier_crt_seed_v1.json`.
Its SHA-256 is
`a4f2e27d63bbf2160cb8afaed1b171295bf941e99ac8db8f3d2bb85424edaf0c`.
The replay certifies the exhaustive CRT/height result and local reduction at
89, 131, and 137.  It intentionally does not compute a global conductor or a
rank certificate.

## Fermigier rank evaluator and certificates

Replay the reconstructed thirteenth point, the twelve independent generic
section differences, and all 22 independent published E22 points:

```bash
python3 elliptic-curves/scripts/verify_fermigier_rank_certificates.py
```

The replay uses exact finite-field arithmetic and independently checks every
stored group order, generator order, and discrete-log equality with PARI/GP.
The pinned artifact SHA-256 is
`94fc64d7f1744f6a20a0396d32914cd36330107db2538e03ee95cc3e32927051`.
Generate an unpinned copy with:

```bash
python3 elliptic-curves/scripts/run_fermigier_rank_certificates.py \
  --output artifacts/local/elliptic-curves/fermigier_rank_certificates.json
```

Evaluate another adapter parameter, optionally run PARI's bounded quartic
point search, and certify a modularly independent subset of the point cloud:

```bash
python3 elliptic-curves/scripts/evaluate_fermigier_specialization.py 19754/39 \
  --quartic-height 100000 --certify-searched-subset \
  --output artifacts/local/elliptic-curves/e22_evaluation.json
```

The `hyperellratpoints` height limit is a bounded search, not a completeness
claim.  The exact certificate proves only the rank of the selected subset.
For larger batches, an installed `ratpoints` executable can replace PARI and
apply a denominator cap:

```bash
python3 elliptic-curves/scripts/evaluate_fermigier_specialization.py 3251/16 \
  --search-engine ratpoints --quartic-height 2000000 \
  --denominator-bound 13000 --certify-searched-subset \
  --output artifacts/local/elliptic-curves/3251_16_evaluation.json
```

`ratpoints` is optional and is not vendored or installed by this repository.
An existing quiet abscissa-only output can be replayed without rerunning the
search by replacing the search options with
`--ratpoints-output artifacts/local/elliptic-curves/POINTS.out`.

## Rank-20 low-conductor near miss

Replay the pinned 58-abscissa search output, exact 20-point independence
certificate, global minimal model, and conductor:

```bash
python3 elliptic-curves/scripts/verify_fermigier_rank20_near_miss.py
```

Regenerate with an installed `ratpoints` 2.1.3:

```bash
python3 elliptic-curves/scripts/run_fermigier_rank20_near_miss.py \
  --output artifacts/local/elliptic-curves/fermigier_rank20_near_miss.json
```

The generator runs height `2000000` with denominator bound `13000`.  To replay
an already captured quiet output instead, add `--ratpoints-output PATH`.  The
pinned artifact SHA-256 is
`8416e835887236e9e4eafcb01384a710ce4f1be0628701a97f4a7d7a07fe63b1`.
It is a rank-at-least-20 near miss, not a target solution.

## Staged Fermigier score sweep

Build and run the dependency-light C++ ranking pass with:

```bash
g++ -O3 -march=native -fopenmp -std=c++20 \
  -o /tmp/fermigier-score-sweep \
  elliptic-curves/ecsearch/fermigier_score_sweep.cpp
OMP_NUM_THREADS=32 /tmp/fermigier-score-sweep 100000 500 20000 \
  > artifacts/local/elliptic-curves/fermigier_score_sweep.tsv
```

The three integer arguments are the maximum numerator, maximum denominator,
and retained output count.  The score is a search heuristic, not a rank or
conductor computation, and the command intentionally writes only to the
ignored local cache.

## Kihara rank-14 baseline

Replay the exact specialization and independence certificate with:

```bash
python3 elliptic-curves/scripts/verify_kihara_rank14.py
```

Generate an unpinned copy with:

```bash
python3 elliptic-curves/scripts/run_kihara_rank14.py \
  --output artifacts/local/elliptic-curves/kihara_rank14_t2.json
```

The pinned artifact SHA-256 is
`851ff6da6ccf4f4dca947048edd43846ff7da41161e83fde419747e715a0df46`.
This is a rank-at-least-14 family baseline, not a rank-30 candidate.

## Public rank-29 baseline

Replay all 29 published points and their exact finite-reduction independence
certificate with:

```bash
python3 elliptic-curves/scripts/verify_e29_independence.py
```

Generate an unpinned copy with:

```bash
python3 elliptic-curves/scripts/run_e29_independence.py \
  --output artifacts/local/elliptic-curves/elkies_klagsbrun_e29.json
```

The pinned artifact SHA-256 is
`a585a8bc081c67fc6314b7be8ea29721b465fcd8f147d170b534ecb52395891e`.
It proves the public lower bound 29 locally and exactly; it neither supplies a
thirtieth point nor replays the conditional upper bound.

## Combined gate

```bash
make verify-elliptic-curves PYTHON=python3
```

The repository's default `.venv/bin/python` is not present in every checkout;
the explicit override is sufficient for the dependency-free Python layer.
The curve and local-reduction replays additionally require `gp` on `PATH`.
