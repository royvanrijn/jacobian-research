from sage.all import *
from pathlib import Path
from collections import Counter, defaultdict

# Large rank-17 norm-10 shell.
pari.allocatemem(4 * 1024**3)

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
D = BASE / "data/holes"
D.mkdir(parents=True, exist_ok=True)

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

holes=set(
    int(x)
    for x in (D/"half_lattice_hole_masks.txt")
        .read_text().split()
)

def vec_to_mask(v):
    m=0
    for i,x in enumerate(v):
        if int(x)&1:
            m |= 1<<i
    return m

print("holes =",len(holes),flush=True)
print("enumerating norm-10 vectors ...",flush=True)

S=Q.short_vector_list_up_to_length(6,True)
shell=S[5]

print("norm-10 +/- pairs =",len(shell),flush=True)

counts=Counter()
first={}
examples=defaultdict(list)

for v in shell:
    m=vec_to_mask(v)

    if m not in holes:
        continue

    counts[m]+=1

    if m not in first:
        first[m]=vector(ZZ,v)

    # Keep first few representatives for later K3/quartic work.
    if len(examples[m]) < 4:
        examples[m].append(vector(ZZ,v))

assert set(counts)==holes

hist=Counter(counts.values())

print()
print("MINIMAL-REPRESENTATIVE MULTIPLICITY")
for multiplicity,ncosets in sorted(hist.items()):
    print(
        f"{multiplicity:4d} norm10-pairs"
        f" : {ncosets:6d} cosets"
    )

print()
print("total cosets =",sum(hist.values()))
print(
    "total norm10 pairs in qualifying cosets =",
    sum(k*v for k,v in hist.items())
)

ordered=sorted(
    holes,
    key=lambda m:(-counts[m],m)
)

print()
print("TOP 50 BY MULTIPLICITY")

for rank,m in enumerate(ordered[:50],1):
    print(
        f"{rank:3d}"
        f" mask=0x{m:05x}"
        f" multiplicity={counts[m]:4d}"
        f" witness={first[m]}"
    )

with (D/"half_hole_multiplicity.tsv").open("w") as f:
    f.write(
        "rank\tmask\thex\tmultiplicity\twitness\n"
    )

    for rank,m in enumerate(ordered,1):
        f.write(
            f"{rank}\t{m}\t0x{m:05x}\t"
            f"{counts[m]}\t"
            + " ".join(map(str,first[m]))
            + "\n"
        )

with (D/"half_hole_witnesses.tsv").open("w") as f:
    f.write("mask\trepresentative\n")

    for m in ordered:
        for v in examples[m]:
            f.write(
                f"{m}\t"
                + " ".join(map(str,v))
                + "\n"
            )

print()
print("saved",D/"half_hole_multiplicity.tsv")
print("saved",D/"half_hole_witnesses.tsv")
