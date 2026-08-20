from __future__ import annotations

from pathlib import Path
import argparse
import csv
import json
import shutil

import numpy as np

from rank_growth import (
    atomic_write_json,
    extension_metrics,
    jsonable_metrics,
    load_gram,
    matrix_numerical_rank,
)


parser = argparse.ArgumentParser(
    description=(
        "Score specialization candidates for the first rank-17 -> rank-18 "
        "ignition event using Schur-complement height."
    )
)
parser.add_argument(
    "--candidates",
    required=True,
    help=(
        "TSV with columns candidate_id,gram and optional parameter,points,"
        "new_point,source. Relative paths are resolved beside the TSV."
    ),
)
parser.add_argument("--out", required=True)
parser.add_argument("--base-rank", type=int, default=17)
parser.add_argument("--delta-min", type=float, default=1e-8)
parser.add_argument(
    "--max-hits",
    type=int,
    default=0,
    help="0 keeps every numerically independent candidate; otherwise keep best N.",
)
args = parser.parse_args()

candidate_file = Path(args.candidates).resolve()
out = Path(args.out).resolve()
out.mkdir(parents=True, exist_ok=True)
(out / "checkpoints").mkdir(exist_ok=True)


def read_candidates(path: Path):
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "candidate_id" not in reader.fieldnames or "gram" not in reader.fieldnames:
            raise SystemExit("candidate TSV requires candidate_id and gram columns")
        yield from reader


def resolve_optional(row: dict[str, str], key: str) -> Path | None:
    raw = (row.get(key) or "").strip()
    if not raw:
        return None
    p = Path(raw)
    if not p.is_absolute():
        p = candidate_file.parent / p
    return p.resolve()


records = []

for seq, row in enumerate(read_candidates(candidate_file)):
    candidate_id = (row.get("candidate_id") or f"candidate-{seq:06d}").strip()
    gram_path = resolve_optional(row, "gram")
    if gram_path is None or not gram_path.exists():
        print(f"SKIP|id={candidate_id}|reason=missing_gram|path={gram_path}", flush=True)
        continue

    try:
        gram = load_gram(gram_path)
    except Exception as exc:
        print(f"SKIP|id={candidate_id}|reason=bad_gram|error={exc}", flush=True)
        continue

    if gram.shape[0] <= args.base_rank:
        print(
            f"SKIP|id={candidate_id}|reason=no_extra_vector|shape={gram.shape}",
            flush=True,
        )
        continue

    numerical_rank = matrix_numerical_rank(gram)
    scored = []

    for j in range(args.base_rank, gram.shape[0]):
        m = extension_metrics(gram, args.base_rank, j)
        if np.isfinite(m.orthogonal_height) and m.orthogonal_height > args.delta_min:
            scored.append(m)

    if not scored:
        print(
            f"NO_IGNITION|id={candidate_id}|numerical_rank={numerical_rank}",
            flush=True,
        )
        continue

    # Nagao laboratory result: the invariant absolute transverse height is the
    # primary ignition quantity.  Keep the shallowest genuinely positive exit.
    best = min(scored, key=lambda m: m.orthogonal_height)

    record = {
        "candidate_id": candidate_id,
        "parameter": (row.get("parameter") or "").strip(),
        "source": (row.get("source") or "").strip(),
        "gram": str(gram_path),
        "points": str(resolve_optional(row, "points") or ""),
        "new_point": str(resolve_optional(row, "new_point") or ""),
        "numerical_rank": numerical_rank,
        **jsonable_metrics(best),
    }
    records.append((record, row, gram))

    print(
        f"IGNITION|id={candidate_id}"
        f"|rank={numerical_rank}"
        f"|j={best.candidate_index}"
        f"|delta={best.orthogonal_height:.12g}"
        f"|ratio={best.orthogonal_ratio:.12g}",
        flush=True,
    )

records.sort(key=lambda item: item[0]["orthogonal_height"])
if args.max_hits > 0:
    records = records[: args.max_hits]

hits_tsv = out / "hits.tsv"
fields = [
    "rank",
    "candidate_id",
    "parameter",
    "numerical_rank",
    "candidate_index",
    "candidate_height",
    "projection_height",
    "orthogonal_height",
    "orthogonal_ratio",
    "hit_dir",
    "source",
]

with hits_tsv.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()

    for rank, (record, source_row, gram) in enumerate(records, start=1):
        hit_dir = out / f"hit-{rank:06d}"
        hit_dir.mkdir(parents=True, exist_ok=True)

        np.savetxt(hit_dir / "height-gram.txt", gram, fmt="%.17g")
        np.savetxt(
            hit_dir / "pairing-vector.txt",
            np.asarray(record["pairing_vector"], dtype=float),
            fmt="%.17g",
        )
        np.savetxt(
            hit_dir / "projection-coefficients.txt",
            np.asarray(record["projection_coefficients"], dtype=float),
            fmt="%.17g",
        )

        parameter = record.get("parameter", "")
        (hit_dir / "parameter.txt").write_text(str(parameter) + "\n")

        for key, filename in (("points", "points.txt"), ("new_point", "new-point.txt")):
            src = resolve_optional(source_row, key)
            if src is not None and src.exists():
                shutil.copy2(src, hit_dir / filename)

        payload = {
            "version": 1,
            "kind": "rank18-ignition",
            **record,
            "hit_rank": rank,
            "hit_dir": str(hit_dir),
        }
        atomic_write_json(hit_dir / "ignition.json", payload)

        writer.writerow(
            {
                "rank": rank,
                "candidate_id": record["candidate_id"],
                "parameter": record["parameter"],
                "numerical_rank": record["numerical_rank"],
                "candidate_index": record["candidate_index"],
                "candidate_height": f'{record["candidate_height"]:.17g}',
                "projection_height": f'{record["projection_height"]:.17g}',
                "orthogonal_height": f'{record["orthogonal_height"]:.17g}',
                "orthogonal_ratio": f'{record["orthogonal_ratio"]:.17g}',
                "hit_dir": str(hit_dir),
                "source": record["source"],
            }
        )

atomic_write_json(
    out / "run.json",
    {
        "version": 1,
        "kind": "rank18-ignition-run",
        "candidate_file": str(candidate_file),
        "base_rank": args.base_rank,
        "delta_min": args.delta_min,
        "hits": len(records),
        "hits_tsv": str(hits_tsv),
    },
)

print()
print("DONE")
print("hits =", len(records))
print("saved =", hits_tsv)
