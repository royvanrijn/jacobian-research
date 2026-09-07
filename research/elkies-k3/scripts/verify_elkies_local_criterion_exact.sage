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
n = H.nrows()

print("rank =",n)
print("det =",H.det())

S,U,V = H.smith_form()

assert U*H*V == S

diag=[abs(ZZ(S[i,i])) for i in range(n)]

print("SNF =",diag)

# Last diagonal coordinate has modulus 948.
i=max(j for j,d in enumerate(diag) if d>1)
D=diag[i]

print("cyclic discriminant modulus =",D)
assert D == 948

Vi = V.inverse()

def canonical_order_p_class(p):
    """
    In SNF coordinates the last standard generator has order D.
    Multiplying it by D/p gives the unique subgroup of order p.

    The quotient is represented as row vectors modulo row(H).
    Since z -> z*V sends row(H) to row(S), transform back by V^-1.
    """
    y=vector(ZZ,n)
    y[i]=D//p

    z=y*Vi

    assert all(x.denominator()==1 for x in z)
    z=vector(ZZ,z)

    # Dual vector in original lattice coordinates.
    v=z*H.inverse()

    # Verify p*v is integral but v is not.
    assert all((p*x).denominator()==1 for x in v)
    assert any(x.denominator()!=1 for x in v)

    q=(v*H*v.column())[0]

    return z,v,q

# Elkies's criterion here concerns the quaternion discriminant D,
# not D*level.  Our rank-17 surface lies on X(6,79):
#
#     quaternion discriminant D = 6 = 2*3
#     Eichler level M = 79
#
# Therefore the only odd ramified prime to test is p=3.
Dquat = ZZ(6)
level = ZZ(79)

for p in [3]:

    print()
    print("="*72)
    print("p =",p)

    z,v,q=canonical_order_p_class(p)

    print("dual numerator z =",z)
    print("dual vector =",v)
    print("norm =",q)

    qp=QQ(q)*p

    print("p * norm =",qp)

    if qp.denominator()==1:
        c=ZZ(qp)

        print("c =",c)
        print("gcd(c,p) =",gcd(c,p))
        print("Legendre(c/p) =",kronecker(c,p))

        required = -kronecker(ZZ(2*Dquat//p),p)

        print("required =",required)
        print("criterion satisfied =",kronecker(c,p)==required)
    else:
        print("WARNING: norm does not have denominator p exactly")

    # Also inspect every nonzero multiple of the order-p generator.
    print()
    print("all nonzero classes in order-p subgroup:")

    seen=set()

    for a in range(1,p):
        va=a*v

        # reduce only conceptually modulo lattice; norm mod integers
        qa=(va*H*va.column())[0]

        cp=qa*p

        if cp.denominator()!=1:
            continue

        c=ZZ(cp)

        key=(c % p)

        if key in seen:
            continue
        seen.add(key)

        print(
            "a =",a,
            "norm =",qa,
            "c mod p =",c%p,
            "chi =",kronecker(c,p)
        )
