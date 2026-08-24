# Reproducing the elliptic-curve programme

Run commands from the repository root. This catalogue covers the active
certificates and current research gates. The pre-cleanup catalogue, including
commands for every bounded historical scan, is preserved as
[`REPRODUCE_2026-08-24.txt`](../archive/elliptic-curves/REPRODUCE_2026-08-24.txt).

## Environment

The dependency-free checks use the repository virtual environment:

```sh
.venv/bin/python --version
```

Some certificates additionally require PARI/GP, Sage/eclib, Singular, or
Magma. Those requirements are stated beside their commands. Raw output and
restart state belong under the ignored `artifacts/local/elliptic-curves/`
tree; do not overwrite pinned generated results during an exploratory run.

## Standard checks

Compile the active Python surface, validate links/status/layout, and run the
current elliptic-curve regression suite:

```sh
make check
make verify-elliptic-curves
```

Audit the evidence labels, JSON/gzip readability, generators, and coverage of
the compact current artifact directory:

```sh
.venv/bin/python elliptic-curves/scripts/audit_artifact_catalog.py
```

## Primary record certificates

### ICARM curve 302: rank at least 31

The fast checker verifies both pinned hashes and recomputes the complete exact
certificate:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/check_icarm_curve302_rank31_pinned.py
```

To generate an unpinned plain-JSON replay with optional primality checks:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_icarm_curve302_rank31.py \
  --output artifacts/local/elliptic-curves/icarm_curve302_rank31_v1.json \
  --verify-primality
```

This proves `rank E(Q) >= 31`, not an unconditional exact rank. See
[`ICARM_CURVE302_RANK31.md`](notes/ICARM_CURVE302_RANK31.md).

### ICARM curve 273: rank at least 30

```sh
.venv/bin/python elliptic-curves/cas/verify_icarm_curve273_rank30.py --check
```

The independent Sage replay is:

```sh
sage -python elliptic-curves/scripts/verify_icarm_curve273_rank30_sage.py
```

See [`ICARM_CURVE273_RANK30.md`](notes/ICARM_CURVE273_RANK30.md).

### Comparative height lattices: ranks 28--31

PARI/GP is required.  Compute the 100-digit canonical height matrices, LLL
transforms and reduced Grams for the public rank-28, rank-29, curve-273 and
curve-302 point lists:

```sh
PYTHONPATH=elliptic-curves/cas \
  python3 elliptic-curves/cas/compare_record_height_lattices.py --digits 100
```

The bounded additive short-vector search is run separately for each declared
height cutoff.  For example, the rank-29 control is:

```sh
PYTHONPATH=elliptic-curves/cas \
  python3 elliptic-curves/cas/search_record_rank17_core.py \
  rank29 --bound 60 --additive-pair-limit 1507 \
  --pool 300 --trials 2000 --seed 29017 \
  --output artifacts/local/elliptic-curves/rank29-r17-additive-ransac.txt
```

The matching commands for `rank28`, `curve273`, and `curve302` use respectively
`(--bound, --additive-pair-limit, --trials, --seed)=(60,2423,800,28017)`,
`(65,2500,800,27317)`, and `(70,3000,800,30217)`.  With those four local
search files present, exactly saturate the candidate spaces and generate their
100-digit core Grams and approximate 1,311-vector profiles:

```sh
PYTHONPATH=elliptic-curves/cas \
  python3 elliptic-curves/cas/analyze_record_rank17_candidates.py --digits 100
```

Calibrate the forced rank-17 fingerprint against ICARM curve 245, whose exact
Fermigier--Mestre rank-12 parent is already reconstructed.  This also fits
exact unimodular R17 shell bases and evaluates the out-of-sample integral-point
enrichment of each selected core.  It exactly transports the thirteen known
generic curve-245 points into the public basis, verifies their rank-12 span and
relation, measures its intersection with the forced rank-17 control, and
replays all 1,311 R17 minimal lines through every fitted basis:

```sh
PYTHONPATH=elliptic-curves/cas \
  python3 elliptic-curves/cas/calibrate_record_rank17_fingerprint.py \
  --digits 100 --restarts 64
```

The negative control fails for both the optimized R17 basis-entry score and
direct pairwise GL(17,Z)-plus-scale fitting.  The fitted full shells are also
highly dispersed.  Integrality enrichment survives as evidence for structured
cores, but the exact curve-245 replay shows that it need not recover the true
generic subgroup.  None of these calculations is a specialization
certificate.

Run the stronger full-shell coordinate descent on all four records and the
known negative control:

```sh
PYTHONPATH=elliptic-curves/cas \
  python3 elliptic-curves/cas/search_record_rank17_shell_embedding.py \
  rank28 rank29 curve273 curve302 curve245-negative-control \
  --restarts 16 --random-steps 80 --maximum-sweeps 25
```

This minimizes the variation of all 1,311 mapped R17 minimal-vector heights.
It also fails calibration: the curve-245 negative control scores better than
the known R17-positive rank-29 specialization, so the values for curves 273
and 302 are diagnostic only.

Extend the exact bounded Mestre fingerprint census to curve 302:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/analyze_icarm_construction_fingerprints.py \
  --include-curve302
```

These are numerical/bounded provenance calculations, not a K3 specialization
certificate.  See
[`RECORD_CURVES_28_29_273_302_HEIGHT_LATTICES.md`](notes/RECORD_CURVES_28_29_273_302_HEIGHT_LATTICES.md).

### ICARM curves 285 and 286: low-conductor candidates

```sh
.venv/bin/python elliptic-curves/cas/analyze_icarm_7fff_zip_sequence.py --check
```

The command exactly replays equation membership, finite-reduction
independence, torsion, invariants, and pairwise `j` comparisons for curves
281, 282, 285, and 286. For curves 285 and 286 it proves the 21-point rank
lower bound. It currently reads the conductor integers from the pinned public
source; it does not independently rerun global minimization or Tate's
algorithm. See [`ICARM_7FFF_ZIP_SEQUENCE.md`](notes/ICARM_7FFF_ZIP_SEQUENCE.md).

## Low-conductor exact baselines

### ICARM curve 245

```sh
.venv/bin/python elliptic-curves/cas/verify_icarm_curve245_rank20.py --check

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/explicit_formula_icarm_curve245_delta22.py --check
```

The first command proves `rank E(Q) >= 20` and independently reconstructs the
exact conductor. The second is a GRH-conditional fixed-fibre upper-bound
diagnostic and does not alter the unconditional status.

### Fermigier rank-20 near miss and `E22`

```sh
.venv/bin/python elliptic-curves/scripts/verify_fermigier_rank20_near_miss.py
.venv/bin/python elliptic-curves/scripts/verify_fermigier_rank_certificates.py
.venv/bin/python elliptic-curves/scripts/verify_benchmarks.py
```

These commands respectively replay the sub-cutoff rank-at-least-20 near miss,
the exact generic-rank/E22 independence certificates, and the family/model
normalization. The literal parameter factor-two discrepancy in the printed
Fermigier source remains open.

### Mestre frontiers

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_mestre_dsquare_four_u197.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/certify_mestre_dsquare_rank19_frontiers.py --check
```

The first is an exact rank-at-least-17 certificate. The second checks the two
rank-at-least-19 frontiers and their conditional explicit-formula diagnostics.

## Family and structural certificates

### New six-root family: exact rank 14 at `T=83/6`

This requires SageMath with eclib and PARI:

```sh
sage -python elliptic-curves/cas/newfamily/certify_rank_t83_6.py \
  --efforts 0 \
  --output artifacts/local/elliptic-curves/newfamily/rank_bounds_t83_6.json
```

The exact subgroup rank is 14 and PARI returns the unconditional interval
`[14,14]`. The pinned exact-rank and lower-bound artifacts are deliberately
separate. See [`NEWFAMILY_RANK14_T83_6.md`](notes/NEWFAMILY_RANK14_T83_6.md)
and [`cas/newfamily/README.md`](cas/newfamily/README.md).

### Kihara and Elkies--Klagsbrun baselines

```sh
.venv/bin/python elliptic-curves/scripts/verify_kihara_rank14.py
.venv/bin/python elliptic-curves/scripts/verify_e29_independence.py
```

These replay unconditional rank lower bounds 14 and 29. The public exact-rank
29 statement is conditional and is not used by the second command.

### Exact Nagao certificates

The retained exact/status entry points are:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_nagao_section7_picard_bound.py \
  --output artifacts/local/elliptic-curves/nagao_section7_picard_bound.json

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/certify_nagao_rank17_frontier.py \
  --output artifacts/local/elliptic-curves/nagao_rank17_frontier.json

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/certify_nagao_rank20_t5081.py \
  --output artifacts/local/elliptic-curves/nagao_rank20_t5081.json
```

The input discovery artifacts for these proofs are preserved and hash-pinned.
The much larger negative search history is archived.

### Fermigier exceptional transport and Mestre two-section geometry

The exact status checkers remain in `elliptic-curves/cas/`:

```sh
PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/analyze_fermigier_exceptional_transport.py \
  --output artifacts/local/elliptic-curves/fermigier_exceptional_transport.json

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_mestre_fermigier_two_section_generic_rank13.py

PYTHONPATH=elliptic-curves/cas \
  .venv/bin/python elliptic-curves/cas/verify_mestre_diameter235_eight_companion_component.py
```

The first command is an exact transport classification; the separate bounded
point searches retain their bounded label even though their enumeration is
exact.

## Current open computations

- The curve-273 residual 2-Selmer pipeline is under `cas/` with `bnf_free`,
  `residual_selmer`, and `curve273` in the filenames. Its intermediate local
  artifacts remain ignored until a complete certificate exists.
- The current low-conductor searches are the retained Fermigier rank-20,
  denominator-offset, mixed-small-prime, and six-root drivers.
- The H3/rootless-MW17 equation transport lives primarily in `elkies-k3/` and
  has its own reproduction catalogue.

No command in this section turns a partial Selmer calculation, timeout, score,
or bounded negative search into a rank theorem.

## Historical searches

The archived command snapshot and manifest are:

- [`REPRODUCE_2026-08-24.txt`](../archive/elliptic-curves/REPRODUCE_2026-08-24.txt);
- [`MANIFEST.tsv`](../archive/elliptic-curves/MANIFEST.tsv);
- [`archive/elliptic-curves/README.md`](../archive/elliptic-curves/README.md).

The manifest maps every old path to its archive path and records its SHA-256.
Use the historical Git revision named in the archive README when an old script
must be run exactly in its former directory layout.
