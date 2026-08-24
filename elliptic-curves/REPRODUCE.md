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
