from __future__ import annotations

from pathlib import Path
import argparse
import csv


parser = argparse.ArgumentParser(
    description="Build the TSV input consumed by rank-growth ignition/cascade scorers."
)
parser.add_argument(
    "--gram-glob",
    required=True,
    help="Glob for candidate Gram matrices, e.g. 'results/specializations/*/height-gram.txt'",
)
parser.add_argument("--out", required=True)
parser.add_argument(
    "--root",
    default=None,
    help="Optional root used to make stored paths relative; defaults to output directory.",
)
parser.add_argument(
    "--parameter-file",
    default="parameter.txt",
    help="Sibling filename from which to read the specialization parameter when present.",
)
parser.add_argument("--points-file", default="points.txt")
parser.add_argument("--new-point-file", default="new-point.txt")
args = parser.parse_args()

out = Path(args.out).resolve()
out.parent.mkdir(parents=True, exist_ok=True)
root = Path(args.root).resolve() if args.root else out.parent

# pathlib does not support absolute patterns through Path.glob; split at the
# first wildcard by using the current working directory for ordinary repo paths.
import glob
paths = [Path(x).resolve() for x in sorted(glob.glob(args.gram_glob))]

if not paths:
    raise SystemExit(f"No Gram matrices matched: {args.gram_glob}")


def stored_path(path: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


fields = ["candidate_id", "gram", "parameter", "points", "new_point", "source"]

with out.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    writer.writeheader()

    for i, gram in enumerate(paths):
        directory = gram.parent
        parameter_path = directory / args.parameter_file
        points_path = directory / args.points_file
        new_point_path = directory / args.new_point_file
        parameter = parameter_path.read_text().strip() if parameter_path.exists() else ""

        writer.writerow(
            {
                "candidate_id": directory.name or f"candidate-{i:06d}",
                "gram": stored_path(gram),
                "parameter": parameter,
                "points": stored_path(points_path) if points_path.exists() else "",
                "new_point": stored_path(new_point_path) if new_point_path.exists() else "",
                "source": str(directory),
            }
        )

print("candidates =", len(paths))
print("saved =", out)
print("path root =", root)
