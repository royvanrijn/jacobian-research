from sage.all import *
from pathlib import Path
import argparse, math
from fractions import Fraction

ap = argparse.ArgumentParser(description="Search candidate transcendental lattices for Elkies X(6,79), then enumerate CM rank-2 complements.")
ap.add_argument("--coeff-bound", type=int, default=8, help="search |a|,|b|,|d|,|e|,|f| <= B; c is solved from det=-948")
ap.add_argument("--c-bound", type=int, default=100)
ap.add_argument("--top", type=int, default=20, help="keep this many simplest ternary candidates")
ap.add_argument("--vbound", type=int, default=12, help="primitive negative-vector coordinate box")
ap.add_argument("--cm-disc-max", type=int, default=20000)
ap.add_argument("--out", default="artifacts/local/elkies-k3/cm-t2-candidates.txt")
args = ap.parse_args()

DET = -948
TARGET_RAM = (2,3)

def vpf(q,p):
    q=QQ(q); n=ZZ(q.numerator()); d=ZZ(q.denominator())
    v=0
    while n % p == 0: n//=p; v+=1
    while d % p == 0: d//=p; v-=1
    return v, QQ(n)/QQ(d)

def unit_mod(q,p,m=None):
    q=QQ(q)
    mod = p if m is None else m
    n=ZZ(q.numerator()) % mod
    d=ZZ(q.denominator()) % mod
    return (n * inverse_mod(d,mod)) % mod

def hilb(a,b,p):
    a=QQ(a); b=QQ(b)
    if p == -1:
        return -1 if a < 0 and b < 0 else 1
    va,ua=vpf(a,p); vb,ub=vpf(b,p)
    if p == 2:
        u=ZZ(unit_mod(ua,2,8)); v=ZZ(unit_mod(ub,2,8))
        ee = ((u-1)//2)*((v-1)//2) + va*((v*v-1)//8) + vb*((u*u-1)//8)
        return -1 if ee % 2 else 1
    lu=kronecker(unit_mod(ua,p),p)
    lv=kronecker(unit_mod(ub,p),p)
    s = -1 if (va*vb*((p-1)//2)) % 2 else 1
    if vb % 2: s *= lu
    if va % 2: s *= lv
    return int(s)

def diag_q(H):
    # Integral lattice quadratic form q(x)=1/2*x^T H*x.
    A=H.change_ring(QQ)/2
    d1=A[0,0]
    if d1 == 0: return None
    d2=A[1,1]-A[1,0]*A[0,1]/d1
    if d2 == 0: return None
    d3=A.det()/(d1*d2)
    return d1,d2,d3

def clifford_ramification(H):
    ds=diag_q(H)
    if ds is None: return None,None
    d1,d2,d3=ds
    # Even Clifford algebra of <d1,d2,d3>.
    aa=-d1*d2
    bb=-d1*d3
    ps={2}
    for z in (aa.numerator(),aa.denominator(),bb.numerator(),bb.denominator()):
        ps.update(prime_divisors(abs(ZZ(z))))
    ram=tuple(sorted(int(p) for p in ps if hilb(aa,bb,int(p)) == -1))
    if hilb(aa,bb,-1) == -1: ram=ram+(-1,)
    return (aa,bb),ram

def signature21(H):
    # Exact LDL signs, avoiding floating point.
    A=H.change_ring(QQ)
    # Try permutations until nonzero pivots.
    # Sage Permutations(3) uses values 1,2,3, while matrix indices are 0,1,2.
    for perm1 in Permutations(3):
        perm = [ZZ(k)-1 for k in perm1]
        P=matrix(QQ,3,3,lambda i,j: 1 if perm[i]==j else 0)
        M=P*A*P.transpose()
        d1=M[0,0]
        if d1==0: continue
        d2=M[1,1]-M[1,0]^2/d1
        if d2==0: continue
        d3=M.det()/(d1*d2)
        pos=sum(1 for x in (d1,d2,d3) if x>0)
        neg=sum(1 for x in (d1,d2,d3) if x<0)
        return (pos,neg)==(2,1)
    return False

def cyclic_disc(H):
    D=H.smith_form()[0]
    vals=[abs(ZZ(D[i,i])) for i in range(3)]
    return vals == [1,1,948]

def complexity(H):
    vals=[abs(ZZ(x)) for x in H.list()]
    return (max(vals),sum(vals),sum(x*x for x in vals))

def binary_reduce_key(G):
    # G=[[2a,b],[b,2c]], positive definite.
    aa=ZZ(G[0,0]//2); bb=ZZ(G[0,1]); cc=ZZ(G[1,1]//2)
    try:
        Q=BinaryQF([aa,bb,cc]).reduced_form()
        return tuple(map(int,Q))
    except Exception:
        # Safe fallback: enough for deduping most small cases.
        return (int(aa),int(bb),int(cc))

print(f"TSEARCH|stage=start|det={DET}|coeff_bound={args.coeff_bound}|target_ram={TARGET_RAM}",flush=True)

B=args.coeff_bound
found={}
tested=0
for a in range(-B,B+1):
    if a==0: continue
    for b in range(-B,B+1):
        if b==0: continue
        for d in range(-B,B+1):
            den=2*(4*a*b-d*d)
            if den==0: continue
            for e in range(-B,B+1):
                for f in range(-B,B+1):
                    # det [[2a,d,e],[d,2b,f],[e,f,2c]]
                    # = 2c(4ab-d^2)+2def-2af^2-2be^2.
                    rest=2*d*e*f-2*a*f*f-2*b*e*e
                    num=DET-rest
                    if num % den: continue
                    c=num//den
                    if abs(c)>args.c_bound: continue
                    H=matrix(ZZ,[[2*a,d,e],[d,2*b,f],[e,f,2*c]])
                    if H.det()!=DET or not signature21(H): continue
                    if not cyclic_disc(H): continue
                    tested+=1
                    quat,ram=clifford_ramification(H)
                    if ram!=TARGET_RAM: continue
                    # LLL is not applicable to indefinite forms; use sorted raw
                    # Gram as a cheap duplicate key under sign/permutation changes.
                    variants=[]
                    for perm1 in Permutations(3):
                        perm = [ZZ(k)-1 for k in perm1]
                        P=matrix(ZZ,3,3,lambda i,j: 1 if perm[i]==j else 0)
                        G=P*H*P.transpose()
                        for signs in cartesian_product_iterator([[1,-1]]*3):
                            S=diagonal_matrix(ZZ,signs)
                            X=S*G*S
                            variants.append(tuple(map(int,X.list())))
                    key=min(variants)
                    if key not in found or complexity(H)<complexity(found[key]):
                        found[key]=H

cands=sorted(found.values(),key=complexity)[:args.top]
print(f"TSEARCH|stage=ternary_done|tested_cyclic_sig21={tested}|matching_classes_cheap={len(found)}|kept={len(cands)}",flush=True)

lines=[]
allcm={}
for idx,H in enumerate(cands):
    quat,ram=clifford_ramification(H)
    lines.append(f"T|id={idx}|det={H.det()}|smith={list(H.smith_form()[0].diagonal())}|clifford={quat}|ram={ram}|complexity={complexity(H)}")
    lines.extend("TGRAM|%d|%s"%(idx," ".join(map(str,row))) for row in H.rows())

    # A CM specialization is a positive rank-2 sublattice orthogonal to a
    # primitive negative vector of the signature-(2,1) generic T.
    cm={}
    V=args.vbound
    for x in range(-V,V+1):
      for y in range(-V,V+1):
       for z in range(-V,V+1):
        if x==y==z==0 or gcd([x,y,z])!=1: continue
        v=vector(ZZ,[x,y,z])
        n=ZZ(v*H*v)
        if n>=0: continue
        row=matrix(ZZ,1,3,list(v*H))
        K=row.right_kernel_matrix()
        if K.nrows()!=2: continue
        G=K*H*K.transpose()
        if G.det()<=0 or G[0,0]<=0: continue
        disc=ZZ(G.det())
        if disc>args.cm_disc_max: continue
        key=binary_reduce_key(G)
        old=cm.get(key)
        score=(disc,abs(n),max(abs(q) for q in v))
        if old is None or score<old[0]:
            cm[key]=(score,v,n,G)

    cms=sorted(cm.items(),key=lambda kv:kv[1][0])
    print(f"TSEARCH|stage=cm|T={idx}|cm_classes={len(cms)}",flush=True)
    for j,(key,(score,v,n,G)) in enumerate(cms[:50]):
        rec=(idx,j,key,score,v,n,G)
        allcm[(idx,key)]=rec
        lines.append(f"CM|T={idx}|id={j}|binary={key}|disc={G.det()}|v={tuple(v)}|vnorm={n}|score={score}")
        lines.extend("CMGRAM|%d|%d|%s"%(idx,j," ".join(map(str,row))) for row in G.rows())

out=Path(args.out); out.parent.mkdir(parents=True,exist_ok=True)
out.write_text("\n".join(lines)+"\n")
print(f"TSEARCH|stage=done|out={out}|cm_total={len(allcm)}",flush=True)
print("TSEARCH|next=inspect smallest recurring CM binary forms across distinct ternary representatives; then match their discriminants to the CM orbit t=2 using Shimura/CM tables.",flush=True)
