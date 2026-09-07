# Six-root new-family workflow

This directory contains the active Sage/C++ workflow for the quartic family
with roots

```text
(-47,-43,-31,30,45,46).
```

The mathematical sources are
[`NEWFAMILY_QUARTIC_ROOTS_AND_CONSTANT_SECTIONS.md`](../../notes/NEWFAMILY_QUARTIC_ROOTS_AND_CONSTANT_SECTIONS.md),
[`NEWFAMILY_RANK14_T83_6.md`](../../notes/NEWFAMILY_RANK14_T83_6.md), and
[`newfamily_symmetric_quartic_roots.json`](../../families/newfamily_symmetric_quartic_roots.json).

## Proved checkpoint

At `T=83/6`, eleven specialized hidden sections and three additional points
give a rank-14 subgroup. PARI `ellrank` returns the unconditional interval
`[14,14]`, so this specialization has exact rank 14.

Two artifacts are intentionally kept distinct:

- `newfamily_rank14_t83_6_v1.json` proves only `rank >= 14` from exact subgroup
  processing;
- `newfamily_rank14_t83_6_pari_exact_rank_v1.json` adds the matching
  unconditional PARI upper bound and proves `rank = 14`.

They live in `artifacts/generated-results/elliptic-curves/`. Numerical Schur
residuals and search scores were discovery filters, not proof inputs.

## Active pipeline

The canonical files are:

```text
hidden_sections_data.py
hidden_sections.py
verify_hidden_sections.py
newfamily_rank11_common.py
newfamily_rank11_minimal_common.py

make_rational_nagao_tables.py
scan_rational_nagao_tables.cpp
screen_height_rank_rational_candidates.py
measure_specialized_section_heights.py

search_unseeded_extra_points_v3.py
batch_verify_v3_rank_gain_hits.py
pin_batch_rank_result.py
certify_rank_t83_6.py
batch_rank_bounds.py
recover_pari_missing_generators.py
```

The v3 search entry point still imports a small amount of v2 implementation
code. Those files are retained as internal compatibility modules, not as
separate recommended workflows. Superseded standalone drivers and tests are
in `archive/elliptic-curves/cas/newfamily/`.

## Replays

Verify all eleven generic hidden sections:

```sh
sage -python elliptic-curves/cas/newfamily/verify_hidden_sections.py
```

Replay the exact-rank-14 specialization without overwriting a pinned result:

```sh
sage -python elliptic-curves/cas/newfamily/certify_rank_t83_6.py \
  --efforts 0 \
  --output artifacts/local/elliptic-curves/newfamily/rank_bounds_t83_6.json
```

For discovery, use `search_unseeded_extra_points_v3.py`; for exact promotion,
use `batch_verify_v3_rank_gain_hits.py`, then `pin_batch_rank_result.py`.

## Storage and evidence

- Compact exact certificates: `artifacts/generated-results/elliptic-curves/`.
- Raw populations, logs, payloads, and restart state:
  `artifacts/local/elliptic-curves/newfamily/` (ignored).
- Historical scripts, tests, artifacts, and old commands:
  `archive/elliptic-curves/`.

A promoted rank claim must retain the normalized roots and parameter, exact
point coordinates, exact subgroup growth, software and command provenance,
and a separate upper bound if equality is claimed. Root number, score, Schur
residual, or a displayed point count alone is insufficient.
