from sage.all import *
from pathlib import Path
import argparse, math, itertools

ap = argparse.ArgumentParser(description="Compare the recovered K3 transcendental lattice with its natural Clifford/Gross companion forms and test CM discriminants.")
ap.add_argument("--candidates", default="artifacts/local/elkies-k3/cm-t2-candidates.txt")
ap.add_argument("--T", type=int, default=0)
ap.add_argument("--targets", default="3,24", help="absolute CM order discriminants to test")
ap.add_argument("--bound", type=int, default=80)
ap.add_argument("--out", default="artifacts/local/elkies-k3/gross-companion-test.txt")
args = ap.parse_args()

def load_T(path, tid):
    rows=[]
    for ln in Path(path).read_text().splitlines():
        if ln.startswith(f"TGRAM|{tid}|"):
            rows.append([ZZ(x) for x in ln.split("|",2)[2].split()])
    if len(rows)!=3:
        raise SystemExit(f"could not load T={tid}")
    return matrix(ZZ,rows)

def primitive(v):
    return gcd([abs(ZZ(x)) for x in v]) == 1

def primitive_integral_form(A):
    dens = [QQ(x).denominator() for x in A.list()]
    scale = lcm(dens)
    B = (scale * A).change_ring(ZZ)

    coeff = [B[i,i] for i in range(3)]
    coeff += [2*B[0,1], 2*B[0,2], 2*B[1,2]]

    g = gcd([abs(ZZ(x)) for x in coeff if x])

    if g > 1:
        C = matrix(QQ, B) / g
        if all(QQ(x).denominator() == 1 for x in C.list()):
            B = C.change_ring(ZZ)
            scale = QQ(scale) / g

    return scale, B

def evalq(A,v):
    return QQ(v*A*v)

def reps(A,n,B):
    a00,a01,a02 = A[0,0], A[0,1], A[0,2]
    a11,a12 = A[1,1], A[1,2]
    a22 = A[2,2]

    out = []

    for x in range(-B,B+1):
        for y in range(-B,B+1):
            base = a00*x*x + 2*a01*x*y + a11*y*y
            lin = 2*(a02*x + a12*y)

            for z in range(-B,B+1):
                if x == 0 and y == 0 and z == 0:
                    continue
                if gcd([abs(x),abs(y),abs(z)]) != 1:
                    continue

                q = base + lin*z + a22*z*z
                if q == n:
                    out.append((x,y,z))
                    if len(out) >= 20:
                        return out

    return out

H=load_T(args.candidates,args.T)
assert H.det()==-948
Q=H.change_ring(QQ)/2  # q_T(x)=1/2 x^T H x

# For a ternary quadratic space q with Gram Q, the norm form on the
# trace-zero part of the EVEN CLIFFORD ALGEBRA is naturally the second
# exterior-power / adjugate form (up to basis orientation and integral order).
Adj=Q.adjugate()

# Sign chosen so that we can inspect positive values on the relevant companion.
forms=[
    ("T_q", Q),
    ("minus_T_q", -Q),
    ("clifford_adjugate", Adj),
    ("minus_clifford_adjugate", -Adj),
]

# Add inverse/dual presentations explicitly; useful for detecting normalization.
forms += [
    ("T_dual", Q.inverse()),
    ("minus_T_dual", -Q.inverse()),
]

targets=[ZZ(x) for x in args.targets.split(",") if x.strip()]
lines=[]
print(f"GROSS|stage=start|T={args.T}|detH={H.det()}|targets={targets}|bound={args.bound}",flush=True)

for name,A in forms:
    sc,Bint=primitive_integral_form(A)
    det=A.det()
    ev=A.eigenvalues(); sig=(sum(1 for x in ev if x>0), sum(1 for x in ev if x<0))
    print(f"GROSSFORM|name={name}|det={det}|signature={sig}|integral_scale={sc}",flush=True)
    lines.append(f"GROSSFORM|name={name}|det={det}|signature={sig}|integral_scale={sc}")
    lines += ["GROSSGRAM|%s|%s"%(name," ".join(map(str,row))) for row in A.rows()]

    # Test both raw rational normalization and primitive-integral normalization.
    for mode,M in [("raw",A),("primitive_integral",Bint)]:
        for n in targets:
            rr=reps(M,n,args.bound)
            print(f"GROSSREP|name={name}|mode={mode}|n={n}|count_shown={len(rr)}|first={rr[:5]}",flush=True)
            lines.append(f"GROSSREP|name={name}|mode={mode}|n={n}|count_shown={len(rr)}|first={rr[:20]}")

# Also test square classes in T directly. A special negative line has imaginary
# quadratic stabilizer with field determined rationally by the negative norm
# square-class, even though the ORDER conductor depends on the integral Gross lattice.
print("GROSS|square_class_scan",flush=True)
seen={}
B=min(args.bound,50)
for x in range(-B,B+1):
  for y in range(-B,B+1):
   for z in range(-B,B+1):
    v=vector(ZZ,[x,y,z])
    if not primitive(v): continue
    q=evalq(Q,v)
    if q>=0: continue
    num=ZZ(q.numerator()); den=ZZ(q.denominator())
    # squarefree representative of rational square class q.
    sf=ZZ(num*den)
    sign=-1 if sf<0 else 1
    sf=sign*prod(p for p,e in factor(abs(sf)) if e%2)
    if sf in (-3,-6):
        key=(int(sf),QQ(q))
        if key not in seen:
            seen[key]=tuple(map(int,v))

for (sf,q),v in sorted(seen.items(),key=lambda kv:(abs(kv[0][1]),kv[0][0]))[:50]:
    print(f"SQUARECLASS|field_rad={sf}|q={q}|v={v}",flush=True)
    lines.append(f"SQUARECLASS|field_rad={sf}|q={q}|v={v}")

Path(args.out).parent.mkdir(parents=True,exist_ok=True)
Path(args.out).write_text("\n".join(lines)+"\n")
print(f"GROSS|stage=done|out={args.out}",flush=True)
print("GROSS|interpretation=representation by the raw adjugate is rational-Clifford evidence; representation by the primitive-integral companion is stronger but still not a proof that this is the exact Gross lattice of the Eichler order.",flush=True)
