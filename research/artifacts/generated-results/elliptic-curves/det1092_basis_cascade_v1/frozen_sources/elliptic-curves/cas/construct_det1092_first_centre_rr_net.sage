#!/usr/bin/env sage-python
"""Generic-input RR net one fibre above the historical norm10 bisection.

Limits: one20x23 kernel, one fixed member squarefree decomposition; no point
search, parameter sweep, exceptional points, or modification of pilot data.
"""
import hashlib,json,runpy
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,matrix,vector,lcm,gcd
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves';CAS=ROOT/'elliptic-curves/cas'
OUT=ART/'det1092_first_centre_rr_net_v1.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    pp=ART/'curve302_recovered_mw17_parent_v1.json';cp=ART/'det1092_historical_unlock_obstruction_v1.json';loader=CAS/'load_curve302_recovered_parent.sage'
    d=json.loads(cp.read_text());E,basis,_=runpy.run_path(str(loader))['load_curve302_recovered_parent'](pp)
    R=E.base_ring().ring();K=E.base_ring();t=R.gen()
    T=sum((n*P for n,P in zip(d['trace_word'],basis)),E(0))
    funcs=[K(t**j)*v for v,bound in [(K(1),10),(T[0],6),(T[1],4)] for j in range(bound+1)]
    den=R(1)
    for f in funcs:den=den.lcm(f.denominator())
    pol=[R(f*den) for f in funcs];deg=max(f.degree() for f in pol)
    M=matrix(QQ,deg+1,len(pol),lambda i,j:pol[j][i]);ker=M.right_kernel()
    assert M.dimensions()==(20,23) and M.rank()==20 and ker.dimension()==3
    A=[R(d['line_coefficients'][k]) for k in ['f0','f1','f2']]
    def flatten(row):return vector(QQ,[f[j] for f,bound in zip(row,[10,6,4]) for j in range(bound+1)])
    va=flatten(A);vta=flatten([t*f for f in A]);old=matrix(QQ,[va,vta]).row_space()
    raw=next(v for v in ker.basis() if v not in old)
    raw*=lcm(x.denominator() for x in raw);raw=vector(ZZ,raw);raw/=gcd(list(raw))
    if next(v for v in raw if v)<0:raw=-raw
    B=[R(list(raw[:11])),R(list(raw[11:18])),R(list(raw[18:23]))]
    assert M*va==0 and M*vta==0 and M*flatten(B)==0
    assert matrix(QQ,[va,vta,flatten(B)]).rank()==3
    assert A[0](0)+A[1](0)*T[0](0)+A[2](0)*T[1](0)==0
    assert B[0](0)+B[1](0)*T[0](0)+B[2](0)*T[1](0)==0
    restrict=matrix(QQ,[[f(0) for f in row] for row in [A,[t*f for f in A],B]])
    assert restrict.rank()==2
    # Deterministic new kernel vector, not selected by splitting at zero.
    X=PolynomialRing(K,'x');x=X.gen();f0,f1,f2=B;u=f0+f1*x;a1,a2,a3,a4,a6=E.a_invariants()
    elim=u*u-a1*x*u*f2-a3*u*f2-f2*f2*(x**3+a2*x*x+a4*x+a6)
    residual,rem=elim.quo_rem(x-T[0]);assert not rem and residual.degree()==2
    c,b,a=residual.list();disc=b*b-4*a*c
    # Squarefree polynomial model without integer or number-field factorization.
    num=R(disc.numerator());dd=R(disc.denominator());poly=num*dd
    fac=poly.squarefree_decomposition();q=R(fac.unit());h=R(1)
    for f,e in fac:q*=f**(int(e)%2);h*=f**(int(e)//2)
    assert disc==(K(h)/dd)**2*q and q.gcd(q.derivative()).degree()==0
    def rec(p):return list(map(str,p.list()))
    def rat(v):return {'numerator':rec(v.numerator()),'denominator':rec(v.denominator())}
    return {'classification':'verified application and new deduction','status':'PASS_EXACT_RR_NET_AND_SURJECTIVE_FIBRE_RESTRICTION',
        'trace_word':d['trace_word'],'matrix_shape':[20,23],'matrix_rank':20,'kernel_dimension':3,
        'A':[rec(f) for f in A],'B':[rec(f) for f in B],
        'net':'f_i(t;u,v)=B_i(t)+(u+v*t)*A_i(t); homogeneous span A,t*A,B',
        'restriction_at_zero_rank':2,'restriction_kernel':'span(t*A)',
        'slope_to_u':'u=-(B1(0)+m*B2(0))/(A1(0)+m*A2(0)); omitted denominator-zero value lies in the homogeneous chart',
        'generic_class':{'description':'D=C_min+F=2O+5F+phi(w)','self_intersection':2,'arithmetic_genus':2,'D_dot_C_min':0,'D_dot_O':1},
        'fixed_member':{'choice':'B, i.e. u=v=0, fixed by kernel row before examining splitting',
            'residual_coefficients':[rat(f) for f in residual.list()],'q_coefficients':rec(q),
            'discriminant_square_factor':rat(K(h)/dd),'squarefree_degree':int(q.degree()),
            'geometric_genus':int((q.degree()-1)//2)},
        'limits':{'RR_systems':1,'fixed_members':1,'parameter_sweeps':0,'point_searches':0,'exceptional_point_inputs':0},
        'boundary':'All fibrewise chord lines are represented; this is not a rational-solubility or independent-rank prediction.',
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [pp,cp,loader,Path(__file__)]}}
if __name__=='__main__':
    if OUT.exists():raise FileExistsError(OUT)
    d=build();OUT.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n')
    print(d['status'],'fixed squarefree degree',d['fixed_member']['squarefree_degree'],flush=True)
