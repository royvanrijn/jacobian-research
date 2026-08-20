from __future__ import annotations

from pathlib import Path
import argparse
import csv
import json
import shutil

import numpy as np

from rank_growth import (
    atomic_write_json,
    cascade_metrics,
    extension_metrics,
    jsonable_metrics,
    load_gram,
    matrix_numerical_rank,
)


parser = argparse.ArgumentParser(
    description=(
        "Score rank-19/20/21 candidates after an ignition/cascade hit. "
        "Candidates are ranked by shallow transverse height and coupling to "
        "the newest orthogonal rank-growth increment."
    )
)
parser.add_argument(
    "--current-hit",
    required=True,
    help="hit directory from search_rank18_ignition.py or this script",
)
parser.add_argument(
    "--candidates",
    required=True,
    help=(
        "TSV with candidate_id,gram and optional parameter,points,new_point,source. "
        "Each Gram must begin with the current canonical basis, followed by extras."
    ),
)
parser.add_argument("--out", required=True)
parser.add_argument("--original-rank", type=int, default=17)
parser.add_argument("--delta-min", type=float, default=1e-8)
parser.add_argument("--min-alignment", type=float, default=0.85)
parser.add_argument("--basis-tol", type=float, default=1e-7)
parser.add_argument(
    "--max-hits",
    type=int,
    default=0,
    help="0 keeps all passing candidates; otherwise keep best N",
)
args = parser.parse_args()

hit_dir = Path(args.current_hit).resolve()
candidate_file = Path(args.candidates).resolve()
out = Path(args.out).resolve()
out.mkdir(parents=True, exist_ok=True)


def load_hit(path: Path):
    for name in ("cascade.json", "ignition.json"):
        p = path / name
        if p.exists():
            return json.loads(p.read_text()), p
    raise SystemExit(f"No cascade.json or ignition.json in {path}")


def canonical_current_gram(path: Path, meta: dict) -> np.ndarray:
    canonical = path / "current-height-gram.txt"
    if canonical.exists():
        return load_gram(canonical)

    # Backward-compatible construction for ignition hits: first original_rank
    # rows plus the selected ignition candidate from the full Gram.
    full_path = path / "height-gram.txt"
    if not full_path.exists():
        raw = meta.get("gram")
        if not raw:
            raise SystemExit("hit has neither current-height-gram.txt nor height-gram.txt")
        full_path = Path(raw)
    full = load_gram(full_path)
    base_rank = int(meta.get("base_rank", args.original_rank))
    j = int(meta["candidate_index"])
    indices = list(range(base_rank)) + [j]
    current = full[np.ix_(indices, indices)]
    np.savetxt(canonical, current, fmt="%.17g")
    return current


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


parent_meta, parent_meta_path = load_hit(hit_dir)
reference = canonical_current_gram(hit_dir, parent_meta)
base_rank = reference.shape[0]
target_rank = base_rank + 1

if base_rank <= args.original_rank:
    raise SystemExit(
        f"current hit rank {base_rank} must exceed original rank {args.original_rank}"
    )
if target_rank > 21:
    raise SystemExit(f"current basis already has rank {base_rank}; target {target_rank} > 21")

print(
    f"CASCADE_START|base_rank={base_rank}|target_rank={target_rank}"
    f"|min_alignment={args.min_alignment}",
    flush=True,
)

all_rows = []
hits = []

for seq, row in enumerate(read_candidates(candidate_file)):
    candidate_id = (row.get("candidate_id") or f"candidate-{seq:06d}").strip()
    gram_path = resolve_optional(row, "gram")
    if gram_path is None or not gram_path.exists():
        print(f"SKIP|id={candidate_id}|reason=missing_gram", flush=True)
        continue

    try:
        gram = load_gram(gram_path)
    except Exception as exc:
        print(f"SKIP|id={candidate_id}|reason=bad_gram|error={exc}", flush=True)
        continue

    if gram.shape[0] <= base_rank:
        print(f"SKIP|id={candidate_id}|reason=no_extra_vector", flush=True)
        continue

    if not np.allclose(
        gram[:base_rank, :base_rank],
        reference,
        rtol=args.basis_tol,
        atol=args.basis_tol,
    ):
        print(f"SKIP|id={candidate_id}|reason=current_basis_mismatch", flush=True)
        continue

    numerical_rank = matrix_numerical_rank(gram)

    for j in range(base_rank, gram.shape[0]):
        ext = extension_metrics(gram, base_rank, j)
        if not np.isfinite(ext.orthogonal_height) or ext.orthogonal_height <= args.delta_min:
            continue

        cascade = cascade_metrics(
            gram,
            args.original_rank,
            base_rank,
            j,
        )
        alignment = cascade.last_increment_corr
        passes = bool(np.isfinite(alignment) and alignment >= args.min_alignment)

        # Absolute transverse height is primary.  Alignment is a learned
        # cascade prior; dividing by corr^2 smoothly rewards the ~0.95 regime
        # without allowing alignment to turn a non-independent point into a hit.
        cascade_score = ext.orthogonal_height / max(alignment * alignment, 1e-12)

        record = {
            "candidate_id": candidate_id,
            "parameter": (row.get("parameter") or "").strip(),
            "source": (row.get("source") or "").strip(),
            "gram": str(gram_path),
            "points": str(resolve_optional(row, "points") or ""),
            "new_point": str(resolve_optional(row, "new_point") or ""),
            "numerical_rank": numerical_rank,
            "target_rank": target_rank,
            "passes_alignment": passes,
            "cascade_score": float(cascade_score),
            **{f"extension_{k}": v for k, v in jsonable_metrics(ext).items()},
            **{f"cascade_{k}": v for k, v in jsonable_metrics(cascade).items()},
        }
        all_rows.append((record, row, gram, j))
        if passes:
            hits.append((record, row, gram, j))

        print(
            f"CANDIDATE|id={candidate_id}|j={j}"
            f"|delta={ext.orthogonal_height:.12g}"
            f"|align={alignment:.6f}"
            f"|block={cascade.new_block_projection_share:.6f}"
            f"|pass={int(passes)}",
            flush=True,
        )

all_rows.sort(key=lambda item: item[0]["cascade_score"])
hits.sort(key=lambda item: item[0]["cascade_score"])
if args.max_hits > 0:
    hits = hits[: args.max_hits]

candidate_fields = [
    "rank",
    "candidate_id",
    "parameter",
    "target_rank",
    "numerical_rank",
    "candidate_index",
    "orthogonal_height",
    "orthogonal_ratio",
    "last_increment_corr",
    "new_block_projection_share",
    "cascade_score",
    "passes_alignment",
    "source",
]

with (out / "candidates.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=candidate_fields)
    writer.writeheader()
    for rank, (record, _, _, _) in enumerate(all_rows, start=1):
        writer.writerow(
            {
                "rank": rank,
                "candidate_id": record["candidate_id"],
                "parameter": record["parameter"],
                "target_rank": target_rank,
                "numerical_rank": record["numerical_rank"],
                "candidate_index": record["extension_candidate_index"],
                "orthogonal_height": f'{record["extension_orthogonal_height"]:.17g}',
                "orthogonal_ratio": f'{record["extension_orthogonal_ratio"]:.17g}',
                "last_increment_corr": f'{record["cascade_last_increment_corr"]:.17g}',
                "new_block_projection_share": f'{record["cascade_new_block_projection_share"]:.17g}',
                "cascade_score": f'{record["cascade_score"]:.17g}',
                "passes_alignment": int(record["passes_alignment"]),
                "source": record["source"],
            }
        )

with (out / "hits.tsv").open("w", newline="") as handle:
    fields = candidate_fields + ["hit_dir"]
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()

    for rank, (record, source_row, gram, j) in enumerate(hits, start=1):
        child = out / f"hit-{rank:06d}"
        child.mkdir(parents=True, exist_ok=True)

        # The input contract guarantees first base_rank rows are the current
        # canonical basis; append the chosen candidate as the next canonical row.
        indices = list(range(base_rank)) + [j]
        current = gram[np.ix_(indices, indices)]
        np.savetxt(child / "current-height-gram.txt", current, fmt="%.17g")
        np.savetxt(child / "height-gram.txt", gram, fmt="%.17g")
        np.savetxt(
            child / "pairing-vector.txt",
            np.asarray(record["extension_pairing_vector"], dtype=float),
            fmt="%.17g",
        )
        np.savetxt(
            child / "new-block-coefficients.txt",
            np.asarray(record["cascade_new_block_coefficients"], dtype=float),
            fmt="%.17g",
        )
        (child / "parameter.txt").write_text(str(record["parameter"]) + "\n")

        for key, filename in (("points", "points.txt"), ("new_point", "new-point.txt")):
            src = resolve_optional(source_row, key)
            if src is not None and src.exists():
                shutil.copy2(src, child / filename)

        payload = {
            "version": 1,
            "kind": "rank-growth-cascade",
            "parent_hit": str(hit_dir),
            "parent_meta": str(parent_meta_path),
            "base_rank": base_rank,
            "target_rank": target_rank,
            "hit_rank": rank,
            "hit_dir": str(child),
            **record,
        }
        atomic_write_json(child / "cascade.json", payload)

        writer.writerow(
            {
                "rank": rank,
                "candidate_id": record["candidate_id"],
                "parameter": record["parameter"],
                "target_rank": target_rank,
                "numerical_rank": record["numerical_rank"],
                "candidate_index": record["extension_candidate_index"],
                "orthogonal_height": f'{record["extension_orthogonal_height"]:.17g}',
                "orthogonal_ratio": f'{record["extension_orthogonal_ratio"]:.17g}',
                "last_increment_corr": f'{record["cascade_last_increment_corr"]:.17g}',
                "new_block_projection_share": f'{record["cascade_new_block_projection_share"]:.17g}',
                "cascade_score": f'{record["cascade_score"]:.17g}',
                "passes_alignment": 1,
                "source": record["source"],
                "hit_dir": str(child),
            }
        )

atomic_write_json(
    out / "run.json",
    {
        "version": 1,
        "kind": "rank-growth-cascade-run",
        "parent_hit": str(hit_dir),
        "candidate_file": str(candidate_file),
        "original_rank": args.original_rank,
        "base_rank": base_rank,
        "target_rank": target_rank,
        "delta_min": args.delta_min,
        "min_alignment": args.min_alignment,
        "candidates_scored": len(all_rows),
        "hits": len(hits),
    },
)

print()
print("DONE")
print("base rank =", base_rank)
print("target rank =", target_rank)
print("scored =", len(all_rows))
print("hits =", len(hits))
print("saved =", out / "hits.tsv")
