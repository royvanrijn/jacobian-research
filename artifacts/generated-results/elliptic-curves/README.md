# Active elliptic-curve artifacts

This directory is the complete active artifact surface for
[`elliptic-curves/`](../../../elliptic-curves/README.md). Exploratory outputs and superseded results live in the
[`archive`](../../../archive/elliptic-curves/README.md); large or resumable
runs belong under the ignored `artifacts/local/elliptic-curves/` directory.

[`CATALOG.tsv`](CATALOG.tsv) is the index. Every row gives the exact SHA-256,
an evidence label, the governing `MATH_STATUS.json` identifier when one
exists, and a one-line scope statement. The evidence labels are deliberately
strict:

- `theorem-certificate`, `exact-rank-certificate`, and
  `exact-lower-bound-certificate` are proof-bearing;
- `conditional-bound` states its hypothesis;
- `exact-computation` proves only the calculation described;
- `bounded-experiment` never promotes a negative search to a theorem;
- `partial-reproduction` records precisely what remains dependent on a public
  source or missing local replay;
- `source-transcription`, `reproducibility-fixture`, and `search-plan` carry no
  mathematical conclusion by themselves.

Important distinctions made explicit by the catalogue:

- `icarm_curve302_rank31_v1.json.gz` proves rank at least 31, not exact rank 31.
- The ICARM 285/286 analysis exactly proves independence of 21 displayed
  points and now independently replays global minimality and every local
  conductor exponent.
- `conductor_first_near_miss_descent_targets_v1.json` pins exact known-subgroup
  Kummer inputs for four fixed fibres; it does not claim a complete Selmer
  group or rank upper bound.
- `newfamily_rank14_t83_6_v1.json` proves only rank at least 14; the separate
  `newfamily_rank14_t83_6_pari_exact_rank_v1.json` supplies the PARI interval
  `[14,14]` used for the exact-rank statement.
- Bounded Fermigier searches are indexed as experiments even when every
  calculation inside the declared box is exact.

Run the catalogue/archive integrity check with:

```sh
python3 elliptic-curves/scripts/audit_artifact_catalog.py
```

The pre-cleanup bytes and every provenance-only refresh are recorded under
[`archive/elliptic-curves/`](../../../archive/elliptic-curves/README.md).
