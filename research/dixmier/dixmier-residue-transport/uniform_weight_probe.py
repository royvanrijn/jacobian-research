#!/usr/bin/env python3
"""Boundary-weight reduction of the unbounded positional shear category.

This is a structural probe. A generic identity is not a theorem at exceptional
weights until every possible integer specialization is addressed.
"""
from pathlib import Path
from collections import defaultdict
import json
import time
import hashlib
import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices.sdm import sdm_irref,sdm_nullspace_from_rref
from sympy.polys.matrices import DomainMatrix
import verify_transport as a
from explore_degree_five_a2_subprincipal import poisson,pi_power,scale

HERE=Path(__file__).resolve().parent
A=sp.symbols("a")
K=QQ.frac_field(A)
N=10

def add(out,m,c):
    if c:
        out[m]=out.get(m,0)+c
        if not out[m]: del out[m]

def system(weight,top_bound=True,high=False,field=QQ):
    # s=X^weight*f(XQ,X^2Z), t=X^(weight+1)*g(XQ,X^2Z).
    sm=[(j,k) for k in range(3 if high else 4) for j in range(N+1)
        if not(k==3 and top_bound and weight+2*j>10)]
    tm=[(j,k) for k in range(2 if high else 3) for j in range(N+1)
        if not(k==2 and top_bound and weight+2*j>10)]
    cols=[]
    S,T=a.degree_five_family(QQ,-QQ.one/2,QQ.one)
    for side,mons in ((0,sm),(1,tm)):
        for j,k in mons:
            # The X exponent is symbolic; the Poisson bracket formula is
            # translated directly to (Q-degree,Z-degree) output coordinates.
            i=weight+side+j+2*k
            col={}
            for (ii,jj,kk),c0 in (T if side==0 else S).items():
                c=field.convert(c0)
                z=k+kk-1
                if z<0: continue
                first=3*(k*(ii-2*jj)-(i-2*j)*kk)
                second=2*(k*jj-j*kk)
                if side: first,second=-first,-second
                if j+jj<=N-1: add(col,(j+jj,z),c*first)
                if j+jj-1<=N-1: add(col,(j+jj-1,z),c*second)
            cols.append(col)
    rows={}
    keys=sorted(set().union(*(set(c) for c in cols)))
    keyidx={m:i for i,m in enumerate(keys)}
    for c,col in enumerate(cols):
        for m,v in col.items():
            if v: rows.setdefault(keyidx[m],{})[c]=field.convert(v)
    reduced,pivots,nz=sdm_irref(rows)
    kernel,_=sdm_nullspace_from_rref(reduced,field.one,len(cols),pivots,nz)
    return {"sm":sm,"tm":tm,"columns":cols,"rows":rows,"reduced":reduced,
            "pivots":pivots,"kernel":kernel,"rank":len(pivots),"field":field}

def pair_matrix(left,right,wa,wb,lam,field):
    # For different homogeneous sectors the scalar quadratic cross term is
    # Lambda({s_left,t_right}+{s_right,t_left}).
    out={}
    def pb(m,n,x,y):
        j,k=m; jj,kk=n
        if k+kk<1: return field.zero
        i=x+j+2*k; ii=y+1+jj+2*kk
        c1=3*(k*(ii-2*jj)-(i-2*j)*kk)
        c2=2*(k*jj-j*kk)
        z=k+kk-1
        return field.convert(c1)*field.convert(lam.get((j+jj,z),0))+field.convert(c2)*field.convert(lam.get((j+jj-1,z),0))
    sl,tl=len(left["sm"]),len(left["tm"])
    sr,tr=len(right["sm"]),len(right["tm"])
    for i,m in enumerate(left["sm"]):
        for j,n in enumerate(right["tm"]):
            c=pb(m,n,wa,wb)
            if c: out[i,sr+j]=c
    for i,n in enumerate(left["tm"]):
        for j,m in enumerate(right["sm"]):
            c=pb(m,n,wb,wa)
            if c: out[sl+i,j]=c
    return out

def restrict(matrix,left_kernel,right_kernel,field):
    bad=[]
    for i,u in enumerate(left_kernel):
        intermediate={}
        for (r,c),v in matrix.items():
            if r in u: add(intermediate,c,u[r]*v)
        for j,v in enumerate(right_kernel):
            z=sum(c*v.get(k,field.zero) for k,c in intermediate.items())
            if z: bad.append([i,j,str(z)])
    return bad

def main():
    start=time.monotonic()
    cert=json.loads((HERE/"combined_cocycle.json").read_text())
    lam={(j,k):QQ(v) for (i,j,k),v in cert["functional"]}
    assert all(i-j-2*k==12 for (i,j,k),v in cert["functional"])
    high=system(K.gens[0],False,True,K)
    low=system(K(8)-K.gens[0],False,False,K)
    print("GENERIC_RANKS",high["rank"],len(high["kernel"]),low["rank"],len(low["kernel"]),flush=True)
    M=pair_matrix(high,low,K.gens[0],K(8)-K.gens[0],lam,K)
    bad=restrict(M,high["kernel"],low["kernel"],K)
    # Final RREF denominators alone do NOT certify rank under specialization.
    # Freeze actual full-rank minors of the original polynomial matrices.
    R=QQ.poly_ring(A)
    minors=[]
    for sys0 in (high,low):
        pcs=sys0["pivots"]
        transpose={i:{row:values[c] for row,values in sys0["rows"].items() if c in values} for i,c in enumerate(pcs)}
        _,prs,_=sdm_irref(transpose)
        assert len(prs)==len(pcs)
        vals=[[R.from_sympy(K.to_sympy(sys0["rows"].get(r,{}).get(c,K.zero))) for c in pcs] for r in prs]
        det=DomainMatrix(vals,(len(pcs),len(pcs)),R).det()
        expr=R.to_sympy(det)
        minors.append({"rows":prs,"columns":pcs,"rank":len(pcs),"determinant":str(sp.factor(expr)),
                       "det_coefficients":[str(c) for c in sp.Poly(expr,A).all_coeffs()]})
        print("MINOR",len(pcs),sp.factor(expr),flush=True)
    lcm=sp.Poly(1,A)
    for m in minors: lcm*=sp.Poly(m["determinant"],A)
    roots=sp.polys.polytools.ground_roots(lcm)
    integer_roots=sorted(int(r) for r in roots if r.is_Integer)
    out={"generic_bad":bad,"ranks":[high["rank"],low["rank"]],"dimensions":[len(high["kernel"]),len(low["kernel"])],
         "rank_minors":minors,"rank_minor_product":str(sp.factor(lcm.as_expr())),"integer_exception_candidates":integer_roots}
    print("GENERIC_BAD",bad[:2],"EXCEPTIONS",integer_roots,flush=True)
    finite=[]
    for w in range(-2,5):
        x=system(w); y=system(8-w)
        m=pair_matrix(x,y,w,8-w,lam,QQ)
        fail=restrict(m,x["kernel"],y["kernel"],QQ)
        finite.append({"weight":w,"ranks":[x["rank"],y["rank"]],"dimensions":[len(x["kernel"]),len(y["kernel"])],"bad":fail})
    exceptional=[]
    for w in integer_roots:
        if w<=10: continue
        x=system(w); y=system(8-w)
        fail=restrict(pair_matrix(x,y,w,8-w,lam,QQ),x["kernel"],y["kernel"],QQ)
        exceptional.append({"weight":w,"bad":fail})
    out["finite_weights"]=finite;out["exceptions_above10"]=exceptional
    # The inhomogeneous weight-four sector: use an actual retained lower
    # lift as origin, and check its full projected linear variation.
    S,T=a.degree_five_family(QQ,-QQ.one/2,QQ.one)
    fam=a.third_order_family(S,T,QQ)
    bs,bt=fam.base
    assert {i-j-2*k for i,j,k in bs}=={4}
    assert {i-j-2*k for i,j,k in bt}=={5}
    lam3={tuple(m):QQ(v) for m,v in cert["functional"]}
    x=system(4)
    linear_bad=[]
    for v in x["kernel"]:
        ss={(4+j+2*k,j,k):c for i,c in v.items() if i<len(x["sm"]) for j,k in [x["sm"][i]]}
        tt={(5+j+2*k,j,k):c for i,c in v.items() if i>=len(x["sm"]) for j,k in [x["tm"][i-len(x["sm"])]]}
        p=poisson(bs,tt)
        for q,factor in ((poisson(ss,bt),1),(pi_power(ss,T,3),QQ.one/24),(pi_power(S,tt,3),QQ.one/24)):
            for m,c in q.items(): add(p,m,c*factor)
        val=sum(c*lam3.get(m,QQ.zero) for m,c in p.items())
        if val: linear_bad.append(str(val))
    current_bad=[];current_count=0
    for side,k,jmax in ((0,0,7),(0,1,3),(1,0,3)):
        for j in range(jmax+1):
            i=j+(12 if side==0 and k==1 else 10 if side==0 else 11)
            p=poisson({(i,j,k):QQ.one},T) if side==0 else poisson(S,{(i,j,k):QQ.one})
            v=sum(c*lam3.get(m,QQ.zero) for m,c in p.items());current_count+=1
            if v: current_bad.append([side,i,j,k,str(v)])
    out["weight4_linear_bad"]=linear_bad;out["current_bad"]=current_bad;out["current_relevance_checks"]=current_count
    # Retain the generic bases as rational expressions for a separate replay.
    out["generic_bases"]=[{"sm":z["sm"],"tm":z["tm"],"rank":z["rank"],
       "kernel":[[[i,str(K.to_sympy(c))] for i,c in sorted(v.items())] for v in z["kernel"]]} for z in (high,low)]
    # Clear denominators per basis vector. The independent checker can then
    # verify every generic identity using polynomial arithmetic alone.
    polynomial_bases=[]
    for z in (high,low):
        vectors=[]
        incidence={i:sum(i in v for v in z["kernel"]) for i in range(len(z["columns"]))}
        for v in z["kernel"]:
            D=sp.Poly(1,A)
            for c in v.values(): D=sp.lcm(D,sp.Poly(sp.denom(K.to_sympy(c)),A))
            vec={i:sp.Poly(sp.cancel(K.to_sympy(c)*D.as_expr()),A) for i,c in v.items()}
            free=min(i for i,c in v.items() if c==1 and incidence[i]==1)
            vectors.append({"denominator":[str(c) for c in reversed(D.all_coeffs())],"free_coordinate":free,
                            "vector":[[i,[str(c) for c in reversed(p.all_coeffs())]] for i,p in sorted(vec.items())]})
        polynomial_bases.append({"sm":z["sm"],"tm":z["tm"],"rank":z["rank"],"vectors":vectors})
    out["polynomial_bases"]=polynomial_bases
    all_denoms=sp.Poly(1,A)
    for z in polynomial_bases:
        for v in z["vectors"]:
            all_denoms*=sp.Poly(sum(sp.Rational(c)*A**i for i,c in enumerate(v["denominator"])),A)
    out["basis_denominator_integer_roots"]=sorted(int(r) for r in sp.polys.polytools.ground_roots(all_denoms) if r.is_Integer)
    assert set(out["basis_denominator_integer_roots"])<=set(integer_roots)
    out["seed_S"]=[[list(m),str(c)] for m,c in sorted(S.items())]
    out["seed_T"]=[[list(m),str(c)] for m,c in sorted(T.items())]
    out["functional"]=cert["functional"]
    out["status"]="PASS" if not bad and not linear_bad and not current_bad and all(not r["bad"] for r in finite+exceptional) else "FAIL"
    out["checker_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out["seed_source_sha256"]=hashlib.sha256((a.OLD/"explore_degree_five_quantum_residue.py").read_bytes()).hexdigest()
    out["elapsed_seconds"]=time.monotonic()-start
    (HERE/"uniform_weight_probe.json").write_text(json.dumps(out,indent=2)+"\n")
    print("FINITE_BAD",[(r["weight"],len(r["bad"])) for r in finite],"EXCEPTION_BAD",exceptional,flush=True)
    print("LINEAR_BAD",linear_bad,"CURRENT_BAD",current_bad,flush=True)
    assert not bad and not linear_bad and not current_bad
    assert all(not r["bad"] for r in finite+exceptional)

if __name__=="__main__":main()
