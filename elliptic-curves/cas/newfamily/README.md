# Newfamily CAS and replay workflow

This directory contains the canonical Sage/C++ workflow for the six-root quartic family documented in:

- `../../notes/NEWFAMILY_QUARTIC_ROOTS_AND_CONSTANT_SECTIONS.md`
- `../../notes/NEWFAMILY_RANK14_T83_6.md`
- `../../notes/NEWFAMILY_REPLAY_CHECKLIST.md`
- `../../families/newfamily_symmetric_quartic_roots.json`

The primary root set used by the current high-rank specialization search is

```text
(-47,-43,-31,30,45,46)
```

## Repository policy

Exploratory calculations may use `/tmp` and ignored local artifacts, but reusable logic must move into this directory promptly. A result is not considered part of the research record merely because it exists in a terminal log.

Use these locations consistently:

```text
elliptic-curves/cas/newfamily/                  canonical replay/search code
elliptic-curves/notes/                          mathematical interpretation and status
artifacts/local/elliptic-curves/newfamily/     raw scans, logs, large tables, restart data (ignored)
artifacts/generated-results/elliptic-curves/   compact pinned certificates/results (versioned)
```

The current exact high-rank certificate is:

```text
artifacts/generated-results/elliptic-curves/newfamily_rank14_t83_6_v1.json
```

It records a baseline-first exact eclib verification at `T=83/6`:

```text
known hidden sections : rank 11
+ Q2                  : 11 -> 12
+ Q3                  : 12 -> 13
+ Q4                  : 13 -> 14
```

Hence the pinned claim is `rank >= 14`; no rank upper bound is claimed.

## Current canonical pipeline

### Hidden generic sections

```text
export_hidden_sections.py
verify_hidden_sections.py
```

The recovered hidden basis is the preferred generic rank-11 basis. The remaining cleanup task is to remove the last default dependency on `/tmp/newfamily_hidden_sections_complete.sobj` by committing the generated exact section data in a stable source representation.

Generate the deterministic source representation once with:

```bash
sage -python elliptic-curves/cas/newfamily/export_hidden_sections.py \
  --input /tmp/newfamily_hidden_sections_complete.sobj \
  --out elliptic-curves/cas/newfamily/hidden_sections_data.py

sage -python elliptic-curves/cas/newfamily/verify_hidden_sections.py
```

After successful replay, `hidden_sections_data.py` should be committed and the search tools should use it as the canonical source instead of the `.sobj` fallback.

### Rational Nagao search

```text
make_rational_nagao_tables.py
scan_rational_nagao_tables.cpp
```

The table builder creates projective local-symbol tables once. The C++ scanner then explores reduced rational `T=a/b` values without point counting in the hot loop.

### Fast specialization triage

```text
screen_height_rank_rational_candidates.py
measure_specialized_section_heights.py
```

The height-rank screen is numerical triage only. It is used to reject pathological specializations cheaply and to estimate whether ordinary point search heights are meaningful. It is not an independence proof by itself.

### Extra-point discovery

```text
search_unseeded_extra_points_v2.py
```

This runs unseeded eclib search on a global minimal model with `pp=0`, then tests discovered points against the known rank-11 height lattice using a two-precision Schur-complement residual. Stable Schur hits are only candidates; they are never promoted as rank gains without exact verification.

The `T=11` exact-rank-11 specialization is the numerical control: dependent points there give relative Schur residuals around `4e-76`, while the three exact new directions at `T=83/6` had residuals around `0.23--0.44`.

### Exact rank-gain verification

```text
batch_verify_v2_rank_gain_hits.py
verify_v2_rank_gain_hit.py
```

The batch verifier is baseline-first:

1. process all eleven known specialized sections;
2. record their exact processed subgroup rank;
3. process every distinct Schur-hit point;
4. record each exact rank increase.

A final processed subgroup rank `r` proves `rank(E(Q)) >= r`. It does not prove equality and does not imply full saturation.

### Pinning a verified result

```text
pin_batch_rank_result.py
```

Once a baseline-first batch result is complete, promote a compact record into the versioned generated-results directory, for example:

```bash
python3 elliptic-curves/cas/newfamily/pin_batch_rank_result.py \
  --input-json artifacts/local/elliptic-curves/newfamily/batch_exact_rank_gain_hits_recovered_v2.json \
  --parameter 83/6 \
  --output artifacts/generated-results/elliptic-curves/newfamily_rank14_t83_6_from_batch.json
```

The pinning step copies only exact verification fields. Raw search populations, payload directories, logs, and restart state remain local.

## Evidence discipline

Raw search logs and large intermediate tables belong under ignored `artifacts/local/elliptic-curves/newfamily/`. Once a computation gives a stable theorem-strength result, promote a compact JSON certificate into `artifacts/generated-results/elliptic-curves/` and add or update a note under `elliptic-curves/notes/`.

For new high-rank specializations, a promoted certificate should include at minimum:

- root set and rational parameter;
- exact specialized known-subgroup rank;
- exact rank gain over that subgroup;
- final processed subgroup rank/lower bound;
- exact coordinates of rank-increasing points;
- root number and minimal-discriminant size when available;
- the replay script used;
- a clear statement that no upper bound is claimed unless separately proved.

Do not infer a full Mordell--Weil rank from numerical height rank, Schur residuals, root number, or the rank of a displayed subgroup alone.
