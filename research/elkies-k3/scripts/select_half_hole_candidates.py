from pathlib import Path
import csv
import random

BASE = Path(__file__).resolve().parents[1]
D = BASE / "data/holes"

src = D / "half_hole_multiplicity.tsv"
out = D / "candidate_half_holes.tsv"

rows=[]

with src.open() as f:
    reader=csv.DictReader(f,delimiter="\t")
    for r in reader:
        r["rank"]=int(r["rank"])
        r["mask"]=int(r["mask"])
        r["multiplicity"]=int(r["multiplicity"])
        rows.append(r)

selected=[]

# Keep all extreme classes.
selected += [
    r for r in rows
    if r["multiplicity"] >= 28
]

# Representative samples from the broader distribution.
rng=random.Random(1729)

for lo,hi,n in [
    (26,27,50),
    (24,25,50),
    (20,23,75),
    (12,19,25),
    (0,11,10),
]:
    pool=[
        r for r in rows
        if lo <= r["multiplicity"] <= hi
        and r not in selected
    ]
    rng.shuffle(pool)
    selected += pool[:n]

# Deduplicate, preserving priority order.
seen=set()
final=[]

for r in selected:
    if r["mask"] in seen:
        continue
    seen.add(r["mask"])
    final.append(r)

with out.open("w") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow([
        "priority",
        "mask",
        "hex",
        "multiplicity",
        "witness",
    ])

    for i,r in enumerate(final,1):
        w.writerow([
            i,
            r["mask"],
            r["hex"],
            r["multiplicity"],
            r["witness"],
        ])

print("selected =",len(final))
print("extreme >=28 =",sum(r["multiplicity"]>=28 for r in final))
print("saved",out)
