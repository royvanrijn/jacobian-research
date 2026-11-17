#!/usr/bin/env python3
"""Mix quadratic outputs before completing squares; exact shared-cube exchanges."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,random,time
from itertools import combinations,combinations_with_replacement
from math import gcd,lcm
from functools import reduce
from power_search import DB,terms,decomposition,linear,sp,fmpq_mat,cubic12,homogeneous18,X
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def plus(a,b,c):
    out=dict(a)
    for j,v in b.items():out[j]=out.get(j,Q(0))+c*v
    return {j:v for j,v in out.items() if v}
def pivot(A,B,basis,p,q):
    value=A[p][q];ar={j:c/value for j,c in A[p].items()};br={j:c/value for j,c in B[p].items()}
    for k in range(len(A)):
        if k==p:continue
        factor=A[k].get(q,Q(0))
        if factor:
            A[k]=plus(A[k],ar,-factor)
            if br:B[k]=plus(B[k],br,-factor)
    A[p]=ar;B[p]=br;basis[p]=q
def signatures(B):
    groups={};scalars={}
    for p,b in enumerate(B):
        if not b:continue
        value=b[min(b)];key=tuple((j,c/value) for j,c in sorted(b.items()))
        groups.setdefault(key,[]).append(p);scalars[p]=value
    byrow={p:group for group in groups.values() for p in group}
    return byrow,scalars
def primitive(v):
    denominator=lcm(*(c.denominator for c in v));w=[int(c*denominator) for c in v]
    g=reduce(gcd,w,0)
    if not g:return None
    if next(c for c in w if c)<0:g=-g
    return tuple(c//g for c in w)
def square_factors(M,rng):
    A=[list(row) for row in M];answer=[]
    while any(c for row in A for c in row):
        choices=[(i,) for i in range(12) if A[i][i]]
        if not choices:choices=[(i,j) for i in range(12) for j in range(i+1,12) if A[i][j]]
        v=rng.choice(choices);a=sum(A[i][j] for i in v for j in v);assert a
        lin=[sum(A[i][j] for i in v)/a for j in range(12)]
        answer.append(primitive(lin))
        A=[[A[i][j]-a*lin[i]*lin[j] for j in range(12)] for i in range(12)]
        assert len(answer)<=12
    return answer
def main():
    out=HERE/"power-basis-pencils";out.mkdir(exist_ok=True);old=json.loads((DB/"polarization.json").read_text())
    oldforms=list(map(tuple,old["forms"]));pool=set()
    for _,_,factors in terms(True):
        for mode in range(5):pool.update(decomposition(factors,mode))
    _,base,_=cubic12(2)
    matrices=[]
    for h in base:
        M=[[Q(0) for _ in range(12)] for _ in range(12)]
        for e,c in sp.Poly(h,*X).terms():
            if sum(e)!=2:continue
            ids=[i for i,power in enumerate(e) for _ in range(power)];i,j=ids
            if i==j:M[i][i]+=Q(c)
            else:M[i][j]+=Q(c)/2;M[j][i]+=Q(c)/2
        if any(c for row in M for c in row):matrices.append(M)
    prng=random.Random(261118)
    pencils=[(i,j,sign) for i,j in combinations(range(len(matrices)),2) for sign in [-1,1]]
    prng.shuffle(pencils);selected=[(i,i,0) for i in range(len(matrices))]+pencils[:20]
    for i,j,sign in selected:
        M=[[matrices[i][a][b]+sign*matrices[j][a][b] for b in range(12)] for a in range(12)]
        for f in square_factors(M,prng):
            for sign_tau in [-1,1]:pool.add(primitive([Q(c) for c in f]+[Q(0)]*5+[Q(sign_tau)]))
    assert set(oldforms)<=pool
    forms=oldforms+sorted(pool-set(oldforms));m=len(forms)
    (out/"pool.json").write_text(json.dumps(dict(forms=forms,pencils=selected),indent=2)+"\n")
    print("POOL",m,flush=True)
    cubes=[]
    for f in forms:
        poly={}
        for i,j,k in combinations_with_replacement([a for a,c in enumerate(f) if c],3):
            e=[0]*18
            for a in (i,j,k):e[a]+=1
            mult=1 if i==k else (3 if i==j or j==k else 6)
            poly[tuple(e)]=mult*f[i]*f[j]*f[k]
        cubes.append(poly)
    _,K,_=cubic12(2);H,_=homogeneous18(K);H[-1]=X[-1]**3
    target=[dict(sp.Poly(f,*X).terms()) for f in H]
    monomials=sorted(set().union(*(set(p) for p in cubes),*(set(p) for p in target)))
    C=fmpq_mat([[str(p.get(e,0)) for p in cubes] for e in monomials])
    T=fmpq_mat([[str(p.get(e,0)) for p in target] for e in monomials])
    seed_path=HERE/"power-basis-anneal/candidate.json"
    seed=json.loads(seed_path.read_text());assert seed["size"]==89
    positions={f:i for i,f in enumerate(forms)}
    initial_indices=[positions[tuple(f)] for f in seed["forms"]]
    priority=initial_indices+[i for i in range(m) if i not in initial_indices]
    ordered=fmpq_mat([[str(C[i,j]) for j in priority] for i in range(C.nrows())])
    rr,rank=ordered.rref()
    base_basis=[priority[next(j for j in range(m) if rr[i,j])] for i in range(rank)]
    assert base_basis[:89]==initial_indices
    joint=fmpq_mat([[str(C[i,j]) for j in base_basis]+[str(C[i,j]) for j in range(m)]+[str(T[i,j]) for j in range(18)] for i in range(C.nrows())])
    rr,r=joint.rref();assert r==rank
    baseA=[{j:Q(str(rr[i,rank+j])) for j in range(m) if rr[i,rank+j]} for i in range(rank)]
    baseB=[{j:Q(str(rr[i,rank+m+j])) for j in range(18) if rr[i,rank+m+j]} for i in range(rank)]
    padding=next(i for i,f in enumerate(forms) if f==(0,)*17+(1,))
    print("RANK",rank,flush=True)
    rng=random.Random(261118);best=89;logs=[];seen=set();started=time.monotonic();steps=0;answer=None
    for restart in range(20):
        A=[dict(row) for row in baseA];B=[dict(row) for row in baseB];basis=list(base_basis);history=[]
        for step in range(1200):
            if time.monotonic()-started>60:break
            active=[p for p,b in enumerate(B) if b];zeros=[p for p,b in enumerate(B) if not b];count=len(active)
            if step%3==0 and zeros:
                # A zero-output basis change preserves the decomposition but changes circuit coordinates.
                candidates=[(p,q) for p in zeros for q in A[p] if q not in basis]
                if candidates:
                    p,q=rng.choice(candidates);pivot(A,B,basis,p,q);history.append([p,q])
            groups,scalars=signatures(B);basis_set=set(basis);current={basis[p] for p in active};options=[];mindelta=999
            for q in range(m):
                if q in basis_set:continue
                added=[k for k in zeros if A[k].get(q)]
                for p in active:
                    ap=A[p].get(q)
                    if not ap or basis[p]==padding:continue
                    removed=[k for k in groups[p] if k!=p and A[k].get(q,Q(0))*scalars[p]==ap*scalars[k]]
                    delta=len(added)-len(removed)
                    if count+delta>best+5:continue
                    support=(current-{basis[p]}-{basis[k] for k in removed})|{q}|{basis[k] for k in added}
                    sig=tuple(sorted(support))
                    if sig in seen:continue
                    mindelta=min(mindelta,delta)
                    options.append((delta,p,q,sig))
            if not options:break
            bucket=mindelta if rng.random()<0.72 else rng.choice(sorted({o[0] for o in options if o[0]<=mindelta+2}))
            chosen_delta,p,q,sig=rng.choice([o for o in options if o[0]==bucket]);before=count;pivot(A,B,basis,p,q);history.append([p,q]);steps+=1
            actual=sum(bool(b) for b in B);assert actual==before+chosen_delta
            seen.add(sig)
            improved=actual<best
            if improved:
                best=actual;print("BEST",best,"restart",restart,"step",step,flush=True)
            if steps%10==0:
                logs.append(dict(restart=restart,step=step,support=actual,best=best,seconds=round(time.monotonic()-started,3)))
                (out/"checkpoint.json").write_text(json.dumps(dict(status="running",source_sha256=sha(Path(__file__)),candidate_forms=m,rank=rank,best=best,steps=steps,log=logs),indent=2,sort_keys=True)+"\n")
            if improved and actual<89:
                rows={basis[k]:B[k] for k in range(rank) if B[k]};selected=sorted(rows)
                mat=fmpq_mat([[str(rows[i].get(o,Q(0))) for o in range(18)] for i in selected])
                selectedC=fmpq_mat([[str(C[a,i]) for i in selected] for a in range(C.nrows())])
                assert selectedC*mat==T
                assert [rows[i].get(17,Q(0)) for i in selected]==[Q(i==padding) for i in selected]
                answer=dict(source_sha256=sha(Path(__file__)),forms=[forms[i] for i in selected],
                    B=[[str(rows[i].get(o,Q(0))) for i in selected] for o in range(18)],
                    size=actual,padding_index=selected.index(padding),restart=restart,pivots=history,
                    scope="Exact shared-cube identity candidate. No smaller Keller lift claim until determinant and collision are constructed.")
                (out/"candidate.json").write_text(json.dumps(answer,indent=2,sort_keys=True)+"\n")
                baseA=[dict(row) for row in A];baseB=[dict(row) for row in B];base_basis=list(basis)
        if time.monotonic()-started>60:break
    receipt=dict(status="candidate" if answer else "bounded-no-improvement",source_sha256=sha(Path(__file__)),
       candidate_forms=m,rank=rank,best=best,steps=steps,seconds=round(time.monotonic()-started,3),log=logs,seed=261118,pool_sha256=sha(out/"pool.json"),seed_candidate_sha256=sha(seed_path),
       source_inputs={name:sha(DB/name) for name in ["verify.py","polarize.py","polarization.json"]},
       scope="Time-limited exact basis-exchange search allowing temporary support growth up to five above the current best in a pool enlarged by exact square completions of quadratic-output pencils; no lower-bound or optimality claim.")
    (out/"checkpoint.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n");print(json.dumps({k:v for k,v in receipt.items() if k!="log"}),flush=True)
if __name__=="__main__":main()

