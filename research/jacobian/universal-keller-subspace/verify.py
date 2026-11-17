#!/usr/bin/env python3
"""Bounded exact controls for the written subspace and dimension theorems."""
from pathlib import Path
from itertools import product
from collections import Counter
import importlib.util
import hashlib
import json
import platform
import time
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('retained_symmetry',ROOT/'jacobian/universal-keller-symmetry/verify.py')
retained=importlib.util.module_from_spec(spec)
spec.loader.exec_module(retained)
Field=retained.Field
checks=[]

def ck(label,condition):
    assert condition,label
    checks.append(label)

def norm(f,a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return a

def pa(f,a,b):
    return norm(f,[f.add(a[i] if i<len(a) else 0,b[i] if i<len(b) else 0)
                   for i in range(max(len(a),len(b)))])

def pm(f,a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]=f.add(c[i+j],f.mul(x,y))
    return norm(f,c)

def pp(f,a,n):
    out=[1]
    while n:
        if n&1:out=pm(f,out,a)
        a=pm(f,a,a); n>>=1
    return out

def pe(f,a,x):
    v=0
    for c in reversed(a):v=f.add(f.mul(v,x),c)
    return v

def compose(f,a,b):
    out=[0]
    for c in reversed(a):out=pa(f,pm(f,out,b),[c])
    return out

def subspace(f,W):
    out=[1]
    for w in sorted(W):out=pm(f,out,[f.mul(f.p-1,w),1])
    assert out[1]
    return [f.mul(c,f.power(out[1],f.q-2)) for c in out]

def deformation(f,a,m):
    assert not a[0]
    return {((j-1)//m,j):c for j,c in enumerate(a) if c}

def deform_compose(f,outer,inner):
    # Sparse bivariate coefficients, using additivity of the inner powers.
    out={}
    for (ui,xi),a in outer.items():
        for (uj,xj),b in inner.items():
            key=(ui+xi*uj,xi*xj)
            out[key]=f.add(out.get(key,0),f.mul(a,f.power(b,xi)))
    return {k:v for k,v in out.items() if v}

def quotient_dimension(gb,variables,bound):
    leading=[p.LM(order=gb.order).exponents for p in gb.polys]
    return sum(not any(all(e>=v for e,v in zip(exp,lm)) for lm in leading)
               for exp in product(range(bound),repeat=len(variables)))

def main():
    started=time.monotonic()
    u,x,z,t=s.symbols('u x z t')
    # Rational-function bases F_p(t) do not acquire F_(p^2) constants.
    rational_controls=[]
    for p,m in [(2,1),(3,2)]:
        K=s.GF(p).frac_field(t)
        Q=s.Poly(s.prod(x-i-j*t for i in range(p) for j in range(p)),x,domain=K)
        Q=Q.mul_ground(1/Q.diff().nth(0))
        h=s.Poly(x-x**p,x,domain=K)
        b=t-t**p
        outer=s.Poly(x-x**p/b**(p-1),x,domain=K)
        ck(f'F{p}(t) actual non-field subspace factorization',outer.compose(h)==Q)
        def deform(poly):
            return s.Poly(sum(poly.nth(j)*u**((j-1)//m)*x**j
                              for j in range(1,poly.degree()+1) if poly.nth(j)),x,domain=s.GF(p).frac_field(t,u))
        D=deform(Q); dh=deform(h); do=deform(outer)
        ck(f'F{p}(t) normalized subspace derivative',Q.diff()==s.Poly(1,x,domain=K))
        ck(f'F{p}(t) weighted global coefficient identity',do.compose(dh)==D)
        ck(f'F{p}(t) both deformed factors Jacobian one',dh.diff().as_expr()==1 and do.diff().as_expr()==1)
        ck(f'F{p}(t) complete degree and rank-one boundary',D.degree()==p*p and s.cancel(D.as_expr().subs(u,0)-x)==0)
        rational_controls.append({'base':f'F{p}(t)','m':m,'degree':p*p,
                                  'polynomial':str(D.as_expr()),'factor_degrees':[p,p]})

    # Three genuinely non-field subspaces inside explicit finite fields.
    finite_controls=[]
    cases=[(Field(2,(1,1,0,1)),2,1),
           (Field(3,(1,2,0,1)),3,2),
           (Field(2,(1,1,0,0,0,0,1)),4,3)]
    for f,B,m in cases:
        poly=sum(c*t**i for i,c in enumerate(f.modulus))
        ck(f'ambient F{f.q} irreducible modulus',s.Poly(poly,t,modulus=f.p).is_irreducible)
        scalars={c for c in range(f.q) if f.power(c,B)==c}
        ck(f'ambient F{f.q} scalar subfield F{B}',len(scalars)==B)
        b=next(v for v in range(f.q) if v not in scalars)
        W={f.add(c,f.mul(d,b)) for c in scalars for d in scalars}
        A=subspace(f,W); H=subspace(f,scalars)
        T={pe(f,H,w) for w in W}; C=subspace(f,T)
        ck(f'F{f.q} actual translation subspace size',len(W)==B*B)
        ck(f'F{f.q} W not a subfield',any(f.mul(a,b) not in W for a in W for b in W))
        support={j for j,c in enumerate(A) if c}
        ck(f'F{f.q} non-field normalized B-power support',support<={1,B,B*B} and A[1]==1)
        ck(f'F{f.q} exact full kernel',{v for v in range(f.q) if pe(f,A,v)==0}==W)
        ck(f'F{f.q} full subspace polynomial flag identity',compose(f,C,H)==A)
        D=deformation(f,A,m); DH=deformation(f,H,m); DC=deformation(f,C,m)
        ck(f'F{f.q} deformation respects actual coefficient composition',deform_compose(f,DC,DH)==D)
        ck(f'F{f.q} both deformed factor degrees B',len(H)-1==B and len(C)-1==B)
        scalar_roots={a for a in scalars if a and f.power(a,m)==1}
        ck(f'F{f.q} complement has full specified order',len(scalar_roots)==m)
        for a in sorted(scalar_roots):
            ck(f'F{f.q} scalar {a} equivariance of normalized map',all(
               pe(f,A,f.mul(a,v))==f.mul(a,pe(f,A,v)) for v in range(f.q)))
        remaining=set(W); scalar_orbits=[]
        while remaining:
            w0=min(remaining)
            orbit={f.mul(a,w0) for a in scalar_roots}
            remaining-=orbit; scalar_orbits.append(orbit)
        ck(f'F{f.q} normalization boundary orbit count',len(scalar_orbits)==1+(len(W)-1)//m)
        ck(f'F{f.q} normalization exact inertia cycle lengths',
           sorted(map(len,scalar_orbits))==[1]+[m]*((len(W)-1)//m))
        ck(f'F{f.q} normalization ramification degree sum',sum(map(len,scalar_orbits))==len(W))
        finite_controls.append({'ambient_field':f.q,'B':B,'m':m,'W':sorted(W),
                                'normalized_subspace_coefficients':A,
                                'inner_coefficients':H,'outer_coefficients':C,'degrees':[B,B]})

    # Exact normalization-to-source identity over the finite target boundary.
    tt,vv=s.symbols('tt vv')
    for B,m,p in [(3,2,3),(4,3,2),(4,1,2)]:
        A=lambda z:z-z**B
        fpoly=lambda uu,z:z-uu**((B-1)//m)*z**B
        ck(f'B{B} m{m} normalization equation divides by t to the actual source',
           s.Poly(s.expand(A(tt*x)-tt*fpoly(tt**m,x)),tt,x,modulus=p).is_zero)

    # Complete fibres of the sharp F4 degree16 example, not just root lists.
    ck('sharp three-coordinate example has global Jacobian one',
       s.Poly(s.det(s.Matrix([u,x-u*x**4,z-u*z**4]).jacobian([u,x,z]))-1,u,x,z,modulus=2).is_zero)
    f=Field(2,(1,1,1)); w,X,Y,L=s.symbols('w X Y L')
    def lift(a):return sum(v*w**i for i,v in enumerate(f.digits(a)))
    whole=[]; census=Counter(); classes=set(); source_points=0
    for uu,vv,ww in product(range(4),repeat=3):
        if uu==0:
            boundary=s.groebner([X-lift(vv),Y-lift(ww),w*w+w+1],X,Y,w,modulus=2,order='grevlex')
            boundary_rank=quotient_dimension(boundary,(X,Y,w),5)//2
            ck(f'F4 complete boundary quotient target {(uu,vv,ww)}',boundary_rank==1)
            source_points+=1
            whole.append({'target':[uu,vv,ww],'whole_rank':boundary_rank,'base_points':1,'framed':False})
            continue
        inv=f.power(uu,2)
        polys=[X**4+lift(inv)*X+lift(f.mul(inv,vv)),
               Y**4+lift(inv)*Y+lift(f.mul(inv,ww)),w*w+w+1]
        gb=s.groebner(polys,X,Y,w,modulus=2,order='grevlex')
        rank=quotient_dimension(gb,(X,Y,w),5)//2
        ck(f'F4 whole rank16 target {(uu,vv,ww)}',rank==16)
        counts={}
        for n in range(1,4):
            fixed=s.groebner(polys+[X**(4**n)-X,Y**(4**n)-Y],X,Y,w,modulus=2,order='grevlex')
            counts[n]=quotient_dimension(fixed,(X,Y,w),5)//2
        if uu==1:
            order=1 if vv==0 and ww==0 else 2
            if order==1:cl=('identity',)
            else:
                orbit={tuple(f.mul(a,v) for v in (vv,ww)) for a in range(1,4)}
                cl=('translation',)+min(orbit)
            expected={n:(16 if n%order==0 else 0) for n in range(1,4)}
        else:
            order=3; cl=('scalar',inv)
            expected={n:(16 if n%3==0 else 1) for n in range(1,4)}
        ck(f'F4 whole extension-point counts target {(uu,vv,ww)}',counts==expected)
        classes.add(cl); census[cl[0]]+=1; source_points+=counts[1]
        whole.append({'target':[uu,vv,ww],'whole_rank':rank,'framed':True,
                      'extension_point_counts':counts,'Frobenius_class':cl})
        # The full frame quotient: lambda^3=1/u and both root equations.
        frame=s.groebner(polys+[L**3+lift(inv)],L,X,Y,w,modulus=2,order='grevlex')
        ck(f'F4 whole specified frame rank48 target {(uu,vv,ww)}',
           quotient_dimension(frame,(L,X,Y,w),5)//2==48)
        ck(f'F4 quotient Frobenius multiplier target {(uu,vv,ww)}',
           frame.reduce(L**4-lift(inv)*L)[1]==0)
    ck('F4 exact eight framed conjugacy classes realized',len(classes)==8)
    ck('F4 complete off-boundary framed census',dict(census)=={'identity':1,'translation':15,'scalar':32})
    ck('F4 whole-map rational source accounting',source_points==4**3)
    ck('F4 sharp weighted source lower bound excludes dimension two',16+3-1>4**2)

    # Verify the weighted class formula independently in several natural actions.
    budget_controls=[]
    for p,m,dim in [(3,2,1),(3,2,2),(5,2,2),(5,4,1),(5,4,2)]:
        W=list(product(range(p),repeat=dim)); scalars=[a for a in range(1,p) if pow(a,m,p)==1]
        G=[(a,b) for a in scalars for b in W]
        def mul(g,h):
            a,b=g; c,d=h
            return (a*c%p,tuple((x+a*y)%p for x,y in zip(b,d)))
        def inverse(g):
            a,b=g; ia=pow(a,-1,p)
            return (ia,tuple(-ia*x%p for x in b))
        unseen=set(G); reps=[]
        while unseen:
            g=min(unseen); orbit={mul(mul(h,g),inverse(h)) for h in G}
            unseen-=orbit; reps.append(g)
        total=sum(sum(tuple((a*x+y)%p for x,y in zip(v,b))==v for v in W) for a,b in reps)
        q=p**dim
        ck(f'group F{p}^{dim} semidirect C{m} exact class count',len(reps)==(q-1)//m+m)
        ck(f'group F{p}^{dim} semidirect C{m} exact weighted budget',total==q+m-1)
        budget_controls.append({'p':p,'e':dim,'m':m,'classes':len(reps),'weighted_budget':total})

    for src in json.loads((HERE/'sources.json').read_text())['local_sources']:
        ck('retained source pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,
             'rational_function_subspace_controls':rational_controls,
             'finite_nonfield_subspace_controls':finite_controls,'whole_F4_fibres':whole,
             'framed_F4_census':dict(census),'weighted_budget_controls':budget_controls,
             'python':platform.python_version(),'sympy':s.__version__,
             'runtime_seconds':round(time.monotonic()-started,3),
             'assurance':'Written all-degree proofs with bounded exact non-field subspace identities, rational-function base controls, complete quotient algebras of all64 F4 target fibres and all48 frame fibres, and independently enumerated small-group conjugacy budgets. The finite-field helper is retained and reused. Not a formalization of Hilbert90, descent, generic monodromy or all-degree atomicity; no independent whole-proof replay, external review or priority certification.',
             'sha256':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']}}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'whole_F4_fibres':len(whole),
                      'whole_frame_fibres':48,'seconds':receipt['runtime_seconds']}))

if __name__=='__main__':main()
