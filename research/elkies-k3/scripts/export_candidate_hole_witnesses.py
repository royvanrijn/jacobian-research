from pathlib import Path
import csv
from collections import defaultdict

BASE = Path(__file__).resolve().parents[1]
D = BASE / "data/holes"

candidate_file = D / "candidate_half_holes.tsv"
witness_file = D / "half_hole_witnesses.tsv"
out_file = D / "candidate_half_hole_witnesses.tsv"

candidates = {}

with candidate_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        candidates[int(row["mask"])] = {
            "priority": int(row["priority"]),
            "multiplicity": int(row["multiplicity"]),
            "hex": row["hex"],
        }

witnesses = defaultdict(list)

with witness_file.open() as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        m = int(row["mask"])
        if m in candidates:
            witnesses[m].append(row["representative"])

with out_file.open("w") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow([
        "priority",
        "mask",
        "hex",
        "multiplicity",
        "witness_index",
        "representative",
    ])

    for m, meta in sorted(
        candidates.items(),
        key=lambda kv: kv[1]["priority"]
    ):
        for i, rep in enumerate(witnesses[m], 1):
            w.writerow([
                meta["priority"],
                m,
                meta["hex"],
                meta["multiplicity"],
                i,
                rep,
            ])

print("candidate cosets =", len(candidates))
print("exported witnesses =", sum(len(v) for v in witnesses.values()))
print("saved", out_file)
