from sage.all import *
from pathlib import Path
from collections import Counter
import argparse

# PARI's default stack is too small for rank-17 short-vector
# enumeration at these bounds.
pari.allocatemem(4 * 1024**3)

ap=argparse.ArgumentParser()
ap.add_argument("--max-norm",type=int,default=18)
args=ap.parse_args()

if args.max_norm % 2:
    raise SystemExit("--max-norm must be even")

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
D = BASE / "data/holes"
D.mkdir(parents=True,exist_ok=True)

def readmat(path):
    rows=[]
    for line in path.read_text().splitlines():
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        line=line.replace("["," ").replace("]"," ")
        try:
            row=[ZZ(x) for x in line.split()]
        except Exception:
            continue
        if row:
            rows.append(row)
    return matrix(ZZ,rows)

H=readmat(J/"short_vector_basis_gram.txt")
n=17

coeff=[]
for i in range(n):
    for j in range(i,n):
        coeff.append(
            H[i,i]//2 if i==j else H[i,j]
        )

Q=QuadraticForm(ZZ,n,coeff)

hole_file = D/"half_lattice_hole_masks.txt"

if not hole_file.exists():
    raise SystemExit(
        f"missing {hole_file}; run enumerate_half_lattice_holes.sage first"
    )

holes=[
    int(x)
    for x in hole_file.read_text().split()
]

hole_set=set(holes)

def vec_to_mask(v):
    m=0
    for i,x in enumerate(v):
        if int(x)&1:
            m |= 1<<i
    return m

# ------------------------------------------------------------
# Enumerate only relevant norm classes.
#
# Every qualifying hole has norm congruent to 2 mod 4 and has
# already been filtered to exclude norm 2 and norm 6.
#
# Therefore possible minimum norms are:
#
#     10, 14, 18, 22, ...
#
# We enumerate only those shells, one at a time.
# ------------------------------------------------------------

minimum = {}
unresolved = set(holes)

for N in range(10, args.max_norm + 1, 4):

    q = N // 2

    print(
        f"enumerating exact norm {N}",
        f"unresolved_before={len(unresolved)}",
        flush=True
    )

    # qfminim bound convention gives vectors with Q(v) <= q.
    # We enumerate through q, then select exact norm N.
    S = Q.short_vector_list_up_to_length(q + 1, True)

    shell = S[q]

    hits = 0

    for v in shell:
        m = vec_to_mask(v)

        if m in unresolved:
            minimum[m] = N
            unresolved.remove(m)
            hits += 1

    print(
        f"norm={N:2d}"
        f" pairs={len(shell):8d}"
        f" new_holes={hits:6d}"
        f" resolved={len(minimum):6d}"
        f" unresolved={len(unresolved):6d}",
        flush=True
    )

    if not unresolved:
        break

# ------------------------------------------------------------
# Highest minimum = deepest hole.
# Unresolved means min norm > search bound.
# ------------------------------------------------------------

hist=Counter(minimum.values())

print()
print("RESOLVED MINIMUM-NORM HISTOGRAM")
for N,c in sorted(hist.items()):
    print(
        f"minnorm={N:2d}"
        f" halfnorm={QQ(N)/4}"
        f" cosets={c}"
    )

unresolved=[
    m for m in holes
    if m not in minimum
]

print()
print(
    f"unresolved (> {args.max_norm}) =",
    len(unresolved)
)

ranked=sorted(
    holes,
    key=lambda m:(
        minimum.get(m,10**9),
        m
    ),
    reverse=True
)

with (D/"ranked_half_lattice_holes.tsv").open("w") as f:
    f.write(
        "rank\tmask\thex\tmin_norm\t"
        "half_hole_norm\tstatus\n"
    )

    for rank,m in enumerate(ranked,1):
        if m in minimum:
            N=minimum[m]
            f.write(
                f"{rank}\t{m}\t0x{m:05x}\t"
                f"{N}\t{QQ(N)/4}\texact-through-bound\n"
            )
        else:
            f.write(
                f"{rank}\t{m}\t0x{m:05x}\t"
                f">{args.max_norm}\t"
                f">{QQ(args.max_norm)/4}\tunresolved\n"
            )

print()
print("TOP 40 DEEPEST")

for rank,m in enumerate(ranked[:40],1):
    if m in minimum:
        N=minimum[m]
        print(
            rank,
            f"0x{m:05x}",
            "min",N,
            "half-hole",QQ(N)/4
        )
    else:
        print(
            rank,
            f"0x{m:05x}",
            f"min>{args.max_norm}"
        )

print()
print(
    "saved",
    D/"ranked_half_lattice_holes.tsv"
)
