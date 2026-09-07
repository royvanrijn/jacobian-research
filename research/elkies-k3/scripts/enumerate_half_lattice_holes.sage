from sage.all import *
from pathlib import Path
from collections import Counter

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"
OUT = BASE / "data/holes"
OUT.mkdir(parents=True, exist_ok=True)

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

H = readmat(J/"short_vector_basis_gram.txt")
n = 17

assert H.det() == 948
assert H.is_positive_definite()
assert all(H[i,i] % 2 == 0 for i in range(n))

coeff=[]
for i in range(n):
    for j in range(i,n):
        coeff.append(
            H[i,i]//2 if i==j else H[i,j]
        )

Q = QuadraticForm(ZZ,n,coeff)

# ------------------------------------------------------------
# Helpers: M/2M represented as 17-bit masks.
# ------------------------------------------------------------

def mask_to_vec(mask):
    return vector(ZZ, [(mask >> i) & 1 for i in range(n)])

def vec_to_mask(v):
    m=0
    for i,x in enumerate(v):
        if int(x) & 1:
            m |= 1 << i
    return m

def norm(v):
    return ZZ((v * H * v.column())[0])

# ------------------------------------------------------------
# Norm mod 4 is constant on each coset modulo 2M.
# Candidate quadratic-section cosets require norm == 2 mod 4.
# ------------------------------------------------------------

residue2=[]

for mask in range(1 << n):
    v=mask_to_vec(mask)
    if norm(v) % 4 == 2:
        residue2.append(mask)

print("total M/2M cosets =",1<<n)
print("norm 2 mod 4 cosets =",len(residue2))

# ------------------------------------------------------------
# A qualifying coset has no representative of norm <10.
#
# For norm == 2 mod 4 the only possibilities below 10 are:
#   norm 2 and norm 6.
#
# Our lattice is rootless, but check both anyway.
# ------------------------------------------------------------

S = Q.short_vector_list_up_to_length(4, True)

bad={}
for q in (1,3):
    N=2*q
    for v in S[q]:
        m=vec_to_mask(v)
        bad[m]=min(bad.get(m,10**9),N)

print("norm-2 pairs =",len(S[1]))
print("norm-6 pairs =",len(S[3]))
print("bad residue cosets =",len(bad))

holes=[
    m for m in residue2
    if m not in bad
]

print()
print("QUALIFYING QUADRATIC-SECTION COSETS =",len(holes))

# ------------------------------------------------------------
# Basic representative statistics.
#
# The 0/1 representative is not generally shortest, but gives
# useful deterministic metadata.
# ------------------------------------------------------------

rep_norm_hist=Counter()

with (OUT/"half_lattice_holes.tsv").open("w") as f:
    f.write(
        "mask\thex\tbit_vector\trepresentative_norm\t"
        "half_representative_norm\n"
    )

    for m in holes:
        v=mask_to_vec(m)
        N=norm(v)

        rep_norm_hist[N]+=1

        f.write(
            f"{m}\t0x{m:05x}\t"
            + "".join(map(str,v))
            + f"\t{N}\t{QQ(N)/4}\n"
        )

print()
print("0/1 REPRESENTATIVE NORM HISTOGRAM")
for N,c in sorted(rep_norm_hist.items()):
    print(N,c)

(OUT/"half_lattice_hole_masks.txt").write_text(
    "\n".join(str(m) for m in holes)+"\n"
)

print()
print("saved",OUT/"half_lattice_hole_masks.txt")
print("saved",OUT/"half_lattice_holes.tsv")
