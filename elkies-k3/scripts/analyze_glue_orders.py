from pathlib import Path
from collections import Counter, defaultdict
from math import gcd, lcm
import re

BASE = Path(__file__).resolve().parents[1]
G = BASE / "data/glue"

mods = [2,2,2,2,2,20]

def order(label):
    o = 1
    for x,m in zip(label,mods):
        if x:
            o = lcm(o, m // gcd(x,m))
    return o

labels = {}

for line in (G/"minimal_vector_cosets.txt").read_text().splitlines():
    lhs,rhs = line.split(":")
    lab = tuple(map(int,lhs.split()))
    vectors = rhs.split()
    labels[lab] = len(vectors)

hist = Counter()
vectors_by_order = Counter()
occupancies = defaultdict(Counter)

for lab,n in labels.items():
    o = order(lab)
    hist[o] += 1
    vectors_by_order[o] += n
    occupancies[o][n] += 1

print("OCCUPIED COSETS BY GROUP ORDER")
for o in sorted(hist):
    print(
        f"order={o:2d}"
        f" occupied={hist[o]:3d}"
        f" vectors={vectors_by_order[o]:4d}"
    )

print()
print("OCCUPANCY HISTOGRAM BY ORDER")

for o in sorted(occupancies):
    print()
    print("order",o)
    for n,c in sorted(occupancies[o].items()):
        print(f"  {n:3d} vectors : {c:3d} cosets")

print()
print("SIGNED TOTAL =",sum(labels.values()))
print("PAIRS =",sum(labels.values())//2)
