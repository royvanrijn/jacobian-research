from sage.all import *
from pathlib import Path
import argparse, ast, itertools, math

ap=argparse.ArgumentParser(description="Compare candidate ternary T discriminant forms with recovered rank-17 MW lattice M.")
ap.add_argument("--mw-gram", default="elkies-k3/data/lattice/rank17_gram.txt")
ap.add_argument("--candidates", default="artifacts/local/elkies-k3/cm-t2-candidates.txt")
ap.add_argument("--out", default="artifacts/local/elkies-k3/transcendental-discform-filter.txt")
a=ap.parse_args()

N=948

def read_gram(path):
    return matrix(ZZ, [[ZZ(x) for x in ln.split()] for ln in Path(path).read_text().splitlines()
                       if ln.strip() and not ln.startswith("#")])

def frac_mod_2(q):
    q=QQ(q)
    # canonical in [0,2)
    k=floor(q/2)
    r=q-2*k
    while r<0: r+=2
    while r>=2: r-=2
    return QQ(r)

def class_order(Hinv,z):
    x=Hinv*vector(ZZ,z)
    return lcm([QQ(xx).denominator() for xx in x])

def find_cyclic_generator(H):
    Hi=H.inverse()
    n=H.nrows()
    # First try standard basis; then tiny combinations.
    vecs=[]
    for i in range(n):
        z=[0]*n; z[i]=1; vecs.append(tuple(z))
    for B in (1,2):
        # enough for these cyclic SNF cases; stop at first full-order class.
        if n <= 3:
            vecs.extend(itertools.product(range(-B,B+1), repeat=n))
        else:
            # rank 17: combinations of pairs/triples of standard basis are enough in practice.
            for i in range(n):
                for j in range(i+1,n):
                    for ci in (-1,1):
                        for cj in (-1,1):
                            z=[0]*n; z[i]=ci; z[j]=cj; vecs.append(tuple(z))
    seen=set()
    for z in vecs:
        if z in seen or all(c==0 for c in z): continue
        seen.add(z)
        if class_order(Hi,z)==abs(H.det()):
            x=Hi*vector(ZZ,z)
            q=frac_mod_2((vector(QQ,z)*Hi*vector(QQ,z)))
            return vector(ZZ,z),x,q
    raise RuntimeError("could not find full-order discriminant generator")

def cyclic_qform_isometric(q1,q2,n=N):
    # For cyclic A=Z/n, generator may change by any unit u mod n.
    # Need q2 == u^2*q1 mod 2Z.
    hits=[]
    for u in range(1,n):
        if gcd(u,n)!=1: continue
        if frac_mod_2(QQ(u*u)*q1-q2)==0:
            hits.append(u)
    return hits

M=read_gram(a.mw_gram)
assert abs(M.det())==N
zg,xg,qM=find_cyclic_generator(M)
print(f"DFORM|M|rank={M.nrows()}|det={M.det()}|generator_z={tuple(zg)}|q={qM}")

# Parse TGRAM blocks from candidate output.
Ts={}
for ln in Path(a.candidates).read_text().splitlines():
    if ln.startswith("TGRAM|"):
        _,sid,row=ln.split("|",2)
        Ts.setdefault(int(sid),[]).append([ZZ(x) for x in row.split()])

results=[]
for tid in sorted(Ts):
    H=matrix(ZZ,Ts[tid])
    assert H.det()==-N
    z,x,qT=find_cyclic_generator(H)
    hits=cyclic_qform_isometric(qT,qM,N)
    ok=bool(hits)
    results.append((tid,ok,qT,tuple(z),hits[:12],H))
    print(f"DFORM|T={tid}|ok={ok}|q={qT}|generator_z={tuple(z)}|unit_hits={hits[:12]}")

surv=[r for r in results if r[1]]
lines=[
    f"DFORM|mw_q={qM}|candidate_count={len(results)}|survivors={len(surv)}"
]
for tid,ok,qT,z,hits,H in surv:
    lines.append(f"SURVIVE|T={tid}|q={qT}|z={z}|units={hits}")
    for row in H.rows():
        lines.append("SURVGRAM|T=%d|%s"%(tid," ".join(map(str,row))))
Path(a.out).write_text("\n".join(lines)+"\n")
print(f"DFORM|stage=done|survivors={len(surv)}|out={a.out}")
if not surv:
    print("DFORM|WARNING=no survivors: either q_T sign convention, lattice scaling, or candidate construction assumption is wrong")
elif len(surv)==1:
    print("DFORM|JACKPOT=single discriminant-form-compatible ternary candidate")
else:
    print("DFORM|next=classify survivors up to integral isometry, then restrict CM complements to survivors")
