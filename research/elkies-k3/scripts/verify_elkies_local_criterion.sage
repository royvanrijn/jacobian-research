from sage.all import *
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
J = BASE / "checkpoints/24A1-JACKPOT-948-r0-n1311"

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

H = readmat(J/"rank17_gram.txt")

print("det =", H.det())

Hi = H.inverse()

# Dual lattice elements can be represented as z * H^{-1},
# with z integral. We seek examples whose norm has denominator p.
#
# Enumerate small z to exhibit witnesses to Elkies's criterion.

N = ZZ(474)

def chi(p,a):
    a %= p
    if a == 0:
        return 0
    return kronecker(a,p)

for p in [3,79]:

    expected = -chi(
        p,
        ZZ(2*N//p)
    )

    print()
    print("="*70)
    print("p =",p)
    print("required chi_p(c) =",expected)

    found = []

    vals=[-3,-2,-1,1,2,3]
    candidates=[]

    for i in range(17):
        for a in vals:
            z=vector(ZZ,17)
            z[i]=a
            candidates.append(z)

    for i in range(17):
        for j in range(i+1,17):
            for a in vals:
                for b in vals:
                    z=vector(ZZ,17)
                    z[i]=a
                    z[j]=b
                    candidates.append(z)

    for z in candidates:
        v = z * Hi
        q = (v * H * v.column())[0]

        qp = QQ(q) * p

        if qp.denominator() != 1:
            continue

        c = ZZ(qp)

        if gcd(c,p) != 1:
            continue

        if chi(p,c) == expected:
            found.append((q,c,z,v))

            print("witness:")
            print("  norm =",q)
            print("  c =",c)
            print("  chi =",chi(p,c))
            print("  z =",z)
            print("  dual vector =",v)
            break

    if not found:
        print("no sparse witness at bound",B)
