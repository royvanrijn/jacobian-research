#!/usr/bin/env python3
"""All-pole strong-period audit for the new shear-envelope functional.

The X cutoff follows from relevance and operator shifts, not a pole search.
"""
from collections import defaultdict
from pathlib import Path
import json
import time
import hashlib
import verify_transport as a
from explore_degree_five_a2_subprincipal import poisson, pi_power, scale, d_z, delta
from verify_degree_five_laurent_quantum_obstruction import monomial_bracket_functional,split_coordinates

HERE=Path(__file__).resolve().parent

def product_period(left,right,lam):
    return sum(c*v*right.get((u[0]-m[0],u[1]-m[1],u[2]-m[2]),a.QQ.zero)
               for m,c in left.items() for u,v in lam.items())

def main():
    start=time.monotonic()
    certificate=json.loads((HERE/"combined_cocycle.json").read_text())
    lam={tuple(m):a.QQ(c) for m,c in certificate["functional"]}
    S,T=a.degree_five_family(a.QQ,-a.QQ.one/2,a.QQ.one)
    sm=a.laurent_monomials(28,3,15,4)
    tm=a.laurent_monomials(23,2,15,4)
    split=len(sm)
    assert min(m[0] for m in lam)==14 and max(m[0] for m in lam)==20
    assert (max(m[0] for m in S),max(m[0] for m in T))==(11,10)
    assert (max(m[2] for m in S),max(m[2] for m in T))==(5,4)
    assert (max(i+j+4*k for i,j,k in S),max(i+j+4*k for i,j,k in T))==(34,29)
    columns=[{m:c for m,c in poisson({m0:a.QQ.one},T).items() if m[0]>=-3} for m0 in sm]
    columns +=[{m:c for m,c in poisson(S,{m0:a.QQ.one}).items() if m[0]>=-3} for m0 in tm]
    rhs={m:c for m,c in scale(pi_power(S,T,3),-a.QQ.one/24).items() if m[0]>=-3}
    print("PROJECTED_COLUMNS",len(sm),len(tm),flush=True)
    base,kernel,rank=a.solve_affine(columns,rhs,a.QQ)
    assert (len(columns),rank,len(kernel))==(4910,3719,1191)
    print("PROJECTED_SOLVE",rank,len(kernel),"seconds",time.monotonic()-start,flush=True)
    bp=split_coordinates(base,split)
    kp=[split_coordinates(v,split) for v in kernel]
    # All pairs that can meet coefficient support, without a quadratic box scan.
    tindex={m:i for i,m in enumerate(tm)}
    B=defaultdict(list)
    for i,(x,q,z) in enumerate(sm):
        hits=set()
        for xx,qq,zz in lam:
            hits.update(((xx-x-1,qq-q,zz-z+1),(xx-x,qq-q+1,zz-z+1)))
        for m in hits:
            if m in tindex:
                v=monomial_bracket_functional((x,q,z),m,lam)
                if v: B[i].append((tindex[m],v))
    def bilinear(s,t):
        return sum(v*w*t.get(j,a.QQ.zero) for i,v in s.items() for j,w in B.get(i,()))
    fixedT=[d_z(delta(T,3-j),j) for j in range(4)]
    fixedS=[delta(d_z(S,3-j),j) for j in range(4)]
    L=[]
    for side,mons in ((0,sm),(1,tm)):
        for m in mons:
            p={m:a.QQ.one}
            val=a.QQ.zero
            for j in range(4):
                if side==0:
                    left=delta(d_z(p,3-j),j); right=fixedT[j]
                else:
                    left=d_z(delta(p,3-j),j); right=fixedS[j]
                val+=((-1)**j*a.comb(3,j))*product_period(left,right,lam)
            L.append(val/24)
    print("BILINEAR_ENTRIES",sum(map(len,B.values())),"linear done",time.monotonic()-start,flush=True)
    lin0=sum(L[i]*c for i,c in base.items())
    fixed=sum(c*lam.get(m,a.QQ.zero) for m,c in pi_power(S,T,5).items())/1920
    value=bilinear(*bp)+lin0+fixed
    print("PERIOD",str(value),flush=True)
    linear_bad=[]
    for i,(vec,(s,t)) in enumerate(zip(kernel,kp)):
        dv=bilinear(bp[0],t)+bilinear(s,bp[1])+sum(L[j]*c for j,c in vec.items())
        if dv: linear_bad.append((i,str(dv)))
    tinc=defaultdict(list)
    for i,(_,t) in enumerate(kp):
        for j,c in t.items(): tinc[j].append((i,c))
    Q=[]
    for s,_ in kp:
        image=defaultdict(lambda:a.QQ.zero)
        for i,c in s.items():
            for j,v in B.get(i,()): image[j]+=c*v
        row=defaultdict(lambda:a.QQ.zero)
        for j,c in image.items():
            if c:
                for k,v in tinc[j]: row[k]+=c*v
        Q.append({j:c for j,c in row.items() if c})
    quadratic_bad=[]
    for i,row in enumerate(Q):
        for j,c in row.items():
            if c+Q[j].get(i,a.QQ.zero): quadratic_bad.append((i,j,str(c+Q[j].get(i,a.QQ.zero))))
    current_bad=[]
    current_count=0
    for which,bound,zcap in ((0,22,1),(1,17,0)):
        for m in a.laurent_monomials(bound,zcap,0,4):
            col=poisson({m:a.QQ.one},T) if which==0 else poisson(S,{m:a.QQ.one})
            val=sum(c*lam.get(k,a.QQ.zero) for k,c in col.items())
            if val: current_bad.append((which,list(m),str(val)))
            current_count+=1
    # Certify every cubic potential basis vector: deg(delta^3 H)<=6,
    # deg(delta^5 H)<=8, and their induced seed shifts lie in the boxes.
    potentials=[]
    for i in range(4):
        for j in range(4-i):
            H={(i,j,0):a.QQ.one}
            B2=scale(delta(H,3),a.QQ.one/24)
            B4=scale(delta(H,5),a.QQ.one/1920)
            assert all(x+y<=6 for x,y,z in B2) and all(x+y<=8 for x,y,z in B4)
            from explore_degree_five_a2_subprincipal import multiply
            for seed,cap,zcap in ((S,28,3),(T,23,2)):
                assert all(x+y+4*z<=cap and z<=zcap for x,y,z in multiply(B2,d_z(seed,3)))
            for seed,cap,zcap in ((S,22,1),(T,17,0)):
                assert all(x+y+4*z<=cap and z<=zcap for x,y,z in multiply(B4,d_z(seed,5)))
                assert not d_z(seed,6)
            potentials.append([i,j])
    # The parametrization packet permits verification without rerunning a
    # matrix factorization or importing the discovery code.
    pack={"S":a.serial(S),"T":a.serial(T),"s_monomials":[list(m) for m in sm],"t_monomials":[list(m) for m in tm],
          "base":[[i,str(c)] for i,c in sorted(base.items())],
          "kernel":[[[i,str(c)] for i,c in sorted(v.items())] for v in kernel],
          "functional":certificate["functional"],"projected_rank":rank}
    (HERE/"combined_linear_packet.json").write_text(json.dumps(pack,separators=(",",":"))+"\n")
    out={"status":"PASS" if not(linear_bad or quadratic_bad or current_bad) and value==1 else "FAIL",
         "lower_columns":len(columns),"rank":rank,"dimension":len(kernel),"period":str(value),
         "linear_bad":linear_bad,"quadratic_bad":quadratic_bad,"current_bad":current_bad,
         "bilinear_entries":sum(map(len,B.values())),"quadratic_entries":sum(map(len,Q)),
         "current_polynomial_columns":current_count,"cubic_potential_basis_checks":potentials,
         "hashes":{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ("combined_cocycle.json","combined_linear_packet.json","combined_laurent.py","verify_transport.py")},
         "source_hashes":{str((a.OLD/n).relative_to(a.ROOT)):hashlib.sha256((a.OLD/n).read_bytes()).hexdigest() for n in ("explore_degree_five_quantum_residue.py","explore_degree_five_a2_subprincipal.py","verify_degree_five_laurent_quantum_obstruction.py")},
         "scope":"All finite-X-Laurent polynomials with lower caps (weight4-degree,Z-order)=(28,3),(23,2) and current caps (22,1),(17,0). Projection X>=-3 and input cutoff X>=-15 are proved by support shifts.",
         "elapsed_seconds":time.monotonic()-start}
    (HERE/"combined_laurent.json").write_text(json.dumps(out,indent=2)+"\n")
    print("RESULT",out["status"],"bad",len(linear_bad),len(quadratic_bad),len(current_bad),"seconds",out["elapsed_seconds"],flush=True)
    assert out["status"]=="PASS"

if __name__=="__main__":main()
