# Elliptic-curve computational archive

This directory preserves the exploratory and superseded computational history
removed from the active `elliptic-curves/` surface on 2026-08-24. Nothing was
discarded: 189 CAS/source programs, 146 matching regression tests, 177
generated artifacts, and two documentation snapshots are retained byte for
byte.

This archive is historical. It is not a source of current mathematical status;
use [`MATH_STATUS.json`](../../MATH_STATUS.json) and the active
[`elliptic-curves/README.md`](../../elliptic-curves/README.md).

## Index

[`MANIFEST.tsv`](MANIFEST.tsv) has one row per preserved item:

```text
original_path    archived_path    sha256    kind
```

The SHA-256 column identifies the bytes at the archived path. It makes the reorganization
auditable and distinguishes a moved file from a refreshed result.

The archive layout is:

- `cas/`: old Python, Sage, Singular, and C++ programs;
- `tests/`: regression tests associated with those programs;
- `artifacts/generated-results/`: bounded searches, scans, probes, and
  superseded certificates;
- `artifacts/snapshots/pre-cleanup-2026-08-24/`: the 52 active artifacts as
  they stood before provenance relocation, with a separate SHA-256 manifest;
- `ARTIFACT_MIGRATION_2026-08-24.tsv`: every current artifact byte change made
  solely to relocate or disambiguate provenance;
- `migrations/`: the deterministic one-time metadata migration used here;
- `REPRODUCE_2026-08-24.txt`: the full pre-cleanup command catalogue;
- `NEWFAMILY_README_2026-08-24.txt`: the pre-cleanup new-family workflow.

## Campaign map

Search by filename prefix rather than reading the directory linearly:

| Prefix or term | Historical campaign |
| --- | --- |
| `elkies_klagsbrun_rank30`, `curve273_rank31` | Searches around the former rank-29/30 frontiers |
| `fermigier_rank22`, `fermigier_global`, `multiple_root` | Fermigier CRT, residue-class, auxiliary-orbit, and global scans |
| `mestre_root_tuple`, `mestre_rank15`, `transverse` | Six-root censuses, specialization searches, and two-section branches |
| `nagao_rank13`, `nagao_rank21`, `nagao_u` | Nagao score, local-CRT, neighborhood, skew-height, and alternate-cover searches |
| `kihara`, `kumar_kuwata`, `elkies_rank18` | Literature/source reconstruction and older baselines |
| `anonymous`, `curve90` | Early anonymous-candidate and saturation investigations |
| `newfamily` | Superseded standalone new-family drivers; the active v3 pipeline remains under `elliptic-curves/cas/newfamily/` |

Examples:

```sh
rg -n "rank31|Selmer|alternate cover" archive/elliptic-curves

awk -F '\t' '$4 == "generated-artifact" && $1 ~ /nagao/' \
  archive/elliptic-curves/MANIFEST.tsv

sha256sum -c <(
  awk -F '\t' 'NR > 1 {print $3 "  " $2}' \
    archive/elliptic-curves/MANIFEST.tsv
)
```

## Exact historical replay

Many programs still contain their original relative paths. For an exact replay
in the former layout, use Git revision
`1ba81f31e98d8b0ef831197fae65fd26bbf4482a` and the original command from
`REPRODUCE_2026-08-24.txt`. For quick code archaeology, the archived files can
usually be imported with

```sh
PYTHONPATH=archive/elliptic-curves/cas \
  .venv/bin/python -m unittest discover -s archive/elliptic-curves/tests -v
```

That historical test command is best-effort: paths embedded in old scripts are
preserved as provenance, so exact execution should use the recorded Git
revision. No archived bounded miss is an upper bound or a present research
claim.
