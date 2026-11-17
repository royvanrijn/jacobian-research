#!/usr/bin/env python3
"""Separate exact replay of the uniform boundary-weight certificate.

Fraction polynomial arithmetic and integer Bareiss determinants; no imports
from the uniform producer, discovery, legacy algebra, SymPy or FLINT.
"""
from fractions import Fraction as F
from pathlib import Path
from math import lcm
import hashlib
import json
import time
from check_combined_independent import bracket,pi,period,read

HERE=Path(__file__).resolve().parent
def trim(p):
    while p and not p[-1]: p.pop()
    return p
def plus(p,q):
    z=[F(0)]*max(len(p),len(q))
    for i,c in enumerate(p): z[i]+=c
    for i,c in enumerate(q): z[i]+=c
    return trim(z)
def times(p,q):
    if not p or not q:return []
    z=[F(0)]*(len(p)+len(q)-1)
    for i,c in enumerate(p):
        for j,d in enumerate(q):z[i+j]+=c*d
    return trim(z)
def scale(p,c):return trim([v*c for v in p])
def val(p,a):
    z=F(0)
    for c in reversed(p):z=z*a+c
    return z
def addpoly(out,key,p):
    z=plus(out.get(key,[]),p)
    if z:out[key]=z
    elif key in out:del out[key]

def matrix(S,T,wa,high=False,top=False):
    sm=[(j,k) for k in range(3 if high else 4) for j in range(11) if not(top and k==3 and wa[0]+2*j>10)]
    tm=[(j,k) for k in range(2 if high else 3) for j in range(11) if not(top and k==2 and wa[0]+2*j>10)]
    cols=[]
    for side,mons in ((0,sm),(1,tm)):
        for j,k in mons:
            col={}; ic=wa[0]+side+j+2*k; slope=wa[1]
            for (ii,jj,kk),c in (T if not side else S).items():
                z=k+kk-1
                if z<0:continue
                sign=-1 if side else 1
                first=scale([F(3*(k*(ii-2*jj)-(ic-2*j)*kk)),F(-3*slope*kk)],c*sign)
                second=[c*sign*2*(k*jj-j*kk)]
                if j+jj<=9:addpoly(col,(j+jj,z),first)
                if 0<=j+jj-1<=9:addpoly(col,(j+jj-1,z),second)
            cols.append(col)
    keys=sorted(set().union(*(set(c) for c in cols)))
    return {"sm":sm,"tm":tm,"cols":cols,"keys":keys}

def cross(left,right,wa,wb,lam):
    out={}
    def pb(m,n,x,y):
        j,k=m;jj,kk=n;z=k+kk-1
        if z<0:return []
        i=x[0]+j+2*k;ii=y[0]+1+jj+2*kk
        z1=scale([F(3*(k*(ii-2*jj)-(i-2*j)*kk)),F(3*(k*y[1]-x[1]*kk))],lam.get((j+jj,z),F(0)))
        return plus(z1,[2*(k*jj-j*kk)*lam.get((j+jj-1,z),F(0))])
    for i,m in enumerate(left["sm"]):
        for j,n in enumerate(right["tm"]):
            p=pb(m,n,wa,wb)
            if p:out[i,len(right["sm"])+j]=p
    for i,n in enumerate(left["tm"]):
        for j,m in enumerate(right["sm"]):
            p=pb(m,n,wb,wa)
            if p:out[len(left["sm"])+i,j]=p
    return out

def restrict(M,ks,ls):
    for u in ks:
        middle={}
        for (i,j),c in M.items():
            if i in u:addpoly(middle,j,times(u[i],c))
        for v in ls:
            z=[]
            for j,c in middle.items():z=plus(z,times(c,v.get(j,[])))
            assert not z,"nonzero restricted polynomial"

def bareiss(a):
    a=[r[:] for r in a];n=len(a);prev=1;sign=1
    for k in range(n-1):
        if not a[k][k]:
            p=next((i for i in range(k+1,n) if a[i][k]),None)
            if p is None:return 0
            a[k],a[p]=a[p],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            aik=a[i][k]
            for j in range(k+1,n):
                num=a[i][j]*pivot-aik*a[k][j]
                q,r=divmod(num,prev);assert not r
                a[i][j]=q
            a[i][k]=0
        prev=pivot
    return sign*a[-1][-1]

def root_factorization(p,roots):
    p=p[:]
    for r in roots:
        while len(p)>1 and not val(p,r):
            q=[F(0)]*(len(p)-1);q[-1]=p[-1]
            for i in range(len(q)-2,-1,-1):q[i]=p[i+1]+r*q[i+1]
            assert p[0]+r*q[0]==0
            p=trim(q)
    assert len(p)==1 and p[0],"uncertified polynomial roots"

def numeric_kernel(sys):
    cols=sys["cols"]; rows=[]
    for key in sys["keys"]:
        rows.append({i:c[key][0] for i,c in enumerate(cols) if key in c})
    piv={}
    for row in rows:
        row=dict(row)
        while row:
            k=min(row)
            if k not in piv:
                c=row[k];piv[k]={i:v/c for i,v in row.items()};break
            c=row[k]
            for i,v in piv[k].items():
                row[i]=row.get(i,F(0))-c*v
                if not row[i]:del row[i]
    free=[i for i in range(len(cols)) if i not in piv]
    kernel=[]
    for f in free:
        v={f:F(1)}
        for k,row in reversed(sorted(piv.items())):
            z=-sum(c*v.get(i,F(0)) for i,c in row.items() if i!=k)
            if z:v[k]=z
        kernel.append({i:[c] for i,c in v.items()})
    return kernel,len(piv)

def main():
    start=time.monotonic()
    d=json.loads((HERE/"uniform_weight_probe.json").read_text());assert d["status"]=="PASS"
    assert hashlib.sha256((HERE/"uniform_weight_probe.py").read_bytes()).hexdigest()==d["checker_sha256"]
    S,T=read(d["seed_S"]),read(d["seed_T"])
    original=json.loads((HERE/"inputs.json").read_text())
    assert S=={tuple(m[:3]):F(c) for m,c in original["S"]}
    assert T=={tuple(m[:3]):F(c) for m,c in original["T"]}
    assert bracket(S,T)=={(0,0,0):F(1)}
    assert {i-j-2*k for i,j,k in S}=={-2} and {i-j-2*k for i,j,k in T}=={-1}
    lam3=read(d["functional"]);assert all(i-j-2*k==12 for i,j,k in lam3)
    assert d["functional"]==json.loads((HERE/"combined_cocycle.json").read_text())["functional"]
    lam={(j,k):c for (i,j,k),c in lam3.items()}
    systems=[matrix(S,T,(0,1),True),matrix(S,T,(8,-1))]
    roots=d["integer_exception_candidates"];assert roots==list(range(-2,12))
    bases=[];interpolations=0
    for sy,b,minor in zip(systems,d["polynomial_bases"],d["rank_minors"]):
        assert sy["sm"]==list(map(tuple,b["sm"])) and sy["tm"]==list(map(tuple,b["tm"]))
        vecs=[{i:trim(list(map(F,p))) for i,p in v["vector"]} for v in b["vectors"]]
        for v,info in zip(vecs,b["vectors"]):
            D=trim(list(map(F,info["denominator"])))
            root_factorization(D,roots)
            assert v[info["free_coordinate"]]==D
            for other in vecs:
                if other is not v:assert info["free_coordinate"] not in other
            residual={}
            for i,p in v.items():
                for key,c in sy["cols"][i].items():addpoly(residual,key,times(p,c))
            assert not residual
        assert len(vecs)+minor["rank"]==len(sy["cols"])
        n=minor["rank"]; rows=minor["rows"];columns=minor["columns"]
        table=[[sy["cols"][c].get(sy["keys"][r],[]) for c in columns] for r in rows]
        D=lcm(*(c.denominator for row in table for p in row for c in p))
        detp=trim([F(c) for c in reversed(minor["det_coefficients"])])
        assert len(detp)<=n+1
        root_factorization(detp,roots)
        for x in range(20,21+n):
            ints=[[int(val(p,x)*D) for p in row] for row in table]
            assert bareiss(ints)==val(detp,x)*D**n
            interpolations+=1
        bases.append(vecs)
        print("PASS independent generic minor and basis",n,flush=True)
    restrict(cross(systems[0],systems[1],(0,1),(8,-1),lam),bases[0],bases[1])
    finite=[]
    for a in list(range(-2,5))+[r for r in roots if r>10]:
        x=matrix(S,T,(a,0),top=True);y=matrix(S,T,(8-a,0),top=True)
        kx,rx=numeric_kernel(x);ky,ry=numeric_kernel(y)
        restrict(cross(x,y,(a,0),(8-a,0),lam),kx,ky)
        finite.append({"weight":a,"ranks":[rx,ry],"dimensions":[len(kx),len(ky)]})
    # Inhomogeneous central sector; higher/lower weights cannot contribute
    # to its linear variation because the functional has pure weight12.
    bs={tuple(m[:3]):F(c) for m,c in original["s2"]};bt={tuple(m[:3]):F(c) for m,c in original["t2"]}
    assert {i-j-2*k for i,j,k in bs}=={4} and {i-j-2*k for i,j,k in bt}=={5}
    lower=bracket(bs,T)
    for m,c in bracket(S,bt).items():lower[m]=lower.get(m,F(0))+c
    for m,c in pi(S,T,3).items():lower[m]=lower.get(m,F(0))+c/24
    assert all(not c for c in lower.values())
    O=bracket(bs,bt)
    for p,f in ((pi(bs,T,3),F(1,24)),(pi(S,bt,3),F(1,24)),(pi(S,T,5),F(1,1920))):
        for m,c in p.items():O[m]=O.get(m,F(0))+c*f
    assert period(O,lam3)==1
    z=matrix(S,T,(4,0),top=True);ks,rank=numeric_kernel(z)
    for v in ks:
        ss={};tt={};split=len(z["sm"])
        for i,c in v.items():
            j,k=(z["sm"][i] if i<split else z["tm"][i-split])
            (ss if i<split else tt)[((4 if i<split else 5)+j+2*k,j,k)]=c[0]
        linear=period(bracket(bs,tt),lam3)+period(bracket(ss,bt),lam3)+period(pi(ss,T,3),lam3)/24+period(pi(S,tt,3),lam3)/24
        assert not linear
    count=0
    for side,k,top in ((0,0,7),(0,1,3),(1,0,3)):
        for j in range(top+1):
            i=j+(12 if side==0 and k else 10 if side==0 else 11)
            assert period(bracket({(i,j,k):F(1)},T) if side==0 else bracket(S,{(i,j,k):F(1)}),lam3)==0
            count+=1
    out={"status":"PASS","generic_kernel_polynomial_identities":True,"generic_cross_period_identity":True,
         "determinant_interpolation_points":interpolations,"generic_ranks":[50,70],"finite_weights":finite,
         "rank_minor_roots_complete":roots,"weight4_linear_variations":len(ks),"current_relevance_checks":count,"period":"1",
         "arithmetic":"Fraction polynomials, integer Bareiss determinants and separate rational Gaussian elimination",
         "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"certificate_sha256":hashlib.sha256((HERE/"uniform_weight_probe.json").read_bytes()).hexdigest(),
         "elapsed_seconds":time.monotonic()-start,"boundary":"Independent finite algebraic certificate; homogeneity, tail relevance and all-degree shear transport are written mathematical proofs, not formal verification."}
    (HERE/"uniform_independent.json").write_text(json.dumps(out,indent=2)+"\n")
    print("PASS independent uniform coefficient certificate; period1; seconds",out["elapsed_seconds"],flush=True)

if __name__=="__main__":main()
