#!/usr/bin/env python3
"""Exact algebra for COERCIVE_GRADIENT_MOMENTS.md, not an analytic certificate."""
from pathlib import Path
from functools import lru_cache
from math import factorial, comb
import hashlib, json, time
import sympy as S
from sympy.polys.rings import ring


def monomials(n, degree):
    if n == 1:
        yield (degree,)
    else:
        for k in range(degree+1):
            for rest in monomials(n-1,degree-k): yield (k,)+rest


def hesse_reduction(parameter, lower):
    K=S.QQ.frac_field(parameter) if isinstance(parameter,S.Symbol) else S.QQ
    R,x,y,z,u=ring('x,y,z,u',K); v=[x,y,z]
    a=K.convert(parameter); delta=1-a**3
    assert delta != 0
    F=[v[i]**2-v[(i+1)%3]*v[(i+2)%3]*a for i in range(3)]
    one=R.one; zero=R.zero
    def mono(m):
        out=one
        for w,k in zip(v,m):out*=w**k
        return out
    def jac(m):
        if max(m)<=1:return mono(m),[zero,zero,zero]
        supp=[i for i,k in enumerate(m) if k]
        i=next(i for i,k in enumerate(m) if k>=2)
        if len(supp)==1:
            h=[zero,zero,zero]; h[i]=v[i]**(m[i]-2)
            tail=v[i]**(m[i]-2)*v[(i+1)%3]*v[(i+2)%3]*a
            if m[i]<=3:return tail,h
            # Here tail is itself divisible by a mixed squared monomial.
            mm=list(m);mm[i]-=2;mm[(i+1)%3]+=1;mm[(i+2)%3]+=1
            jj,hh=jac(tuple(mm))
            return jj*a,[h[k]+hh[k]*a for k in range(3)]
        j=next(j for j in supp if j!=i); k=3-i-j
        mm=list(m);mm[i]-=2;mm[j]-=1; factor=mono(mm).quo_ground(delta)
        h=[zero,zero,zero]
        h[i]=factor*v[j];h[j]=factor*a*v[k];h[k]=factor*a*a*v[i]
        return zero,h
    checked=0
    for d in range(7):
        for m in monomials(3,d):
            normal,h=jac(m)
            assert mono(m)-normal==sum((F[i]*h[i] for i in range(3)),zero)
            assert all(max(e[:3])<=1 for e in normal)
            checked+=1
    if lower is None:return {'all_monomials_degree_at_most_6':checked,'parameter':'lambda; only denominator 1-lambda^3','division_witnesses':'PASS'}
    L=R.from_expr(lower); gi=[L.diff(w) for w in v]
    @lru_cache(None)
    def nf(m):
        normal,h=jac(m)
        if max(m)<=1:return normal,tuple(h)
        tail=sum(((u*(v[i]*h[i]-h[i].diff(v[i]))).quo_ground(3)-(gi[i]*h[i]).quo_ground(3) for i in range(3)),zero)
        out=normal; witness=[(u*hh).quo_ground(3) for hh in h]
        for exp,c in tail.items():
            mm=exp[:3]
            assert sum(mm)<sum(m)
            nn,ww=nf(mm); factor=u**exp[3]*c
            out+=factor*nn
            for i in range(3):witness[i]+=factor*ww[i]
        return out,tuple(witness)
    def reduce(poly):
        out=zero;witness=[zero,zero,zero]
        for exp,c in poly.items():
            nn,ww=nf(exp[:3]);factor=u**exp[3]*c;out+=factor*nn
            for i in range(3):witness[i]+=factor*ww[i]
        assert u*(poly-out)==sum((u*witness[i].diff(v[i])+(3*F[i]+gi[i]-u*v[i])*witness[i] for i in range(3)),zero)
        assert all(max(e[:3])<=1 for e in out)
        return out
    P=x**3+y**3+z**3-x*y*z*(3*a)+L
    basis=[mono(m) for m in [(0,0,0),(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1),(1,1,1)]]
    columns=[reduce(P*b) for b in basis]
    max_u=max(e[3] for col in columns for e in col)
    assert max_u<=5
    serialized=[[str(col) for col in columns],str(L),str(parameter)]
    return {'lambda':str(parameter),'lower_part':str(L),'dimension':8,'connection_columns':8,'max_inverse_t_degree':max_u,
            'exact_twisted_divergence_witnesses':'PASS','jacobian_monomials':checked,
            'matrix_sha256':hashlib.sha256(json.dumps(serialized).encode()).hexdigest(),
            'computed_normal_forms':nf.cache_info().currsize}


def even_gaussian(n):
    return 0 if n%2 else factorial(n)//(2**(n//2)*factorial(n//2))


def wick(poly,variables):
    return S.expand(sum(c*S.prod(even_gaussian(k) for k in m) for m,c in S.Poly(poly,*variables).terms()))


def main():
    started=time.monotonic();x,y,z,w=S.symbols('x y z w');lam=S.symbols('lambda')
    universal=hesse_reduction(lam,None)
    dense=2*x*x+3*x*y-5*x*z+7*y*y+11*y*z-13*z*z+17*x-19*y+23*z+29
    matrices=[hesse_reduction(k,lower) for k,lower in [(0,dense),(2,dense),(2,S.Integer(0))]]
    # Smooth Hesse gradient common zero: any zero coordinate forces all zero;
    # for xyz nonzero, multiplication of the three equations gives delta=0.
    assert S.expand((x*x)*(y*y)*(z*z)-(lam*y*z)*(lam*x*z)*(lam*x*y)-(1-lam**3)*x*x*y*y*z*z)==0
    assert all(S.diff(x**3+y**3+z**3-3*x*y*z,q).subs({x:1,y:1,z:1})==0 for q in [x,y,z])
    # Nodal examples, with exact affine projective patches.
    nodal=[]
    for f,points in [(z*(x*x+y*y)+x**3,[(0,0,1)]),(x*y*z,[(1,0,0),(0,1,0),(0,0,1)])]:
        grad=[S.diff(f,q) for q in [x,y,z]];H=S.hessian(f,[x,y,z])
        patches=[]
        for i in range(3):
            vv=[x,y,z];remaining=[q for q in vv if q!=vv[i]]
            equations=[g.subs(vv[i],1) for g in grad]
            # The patches overlap; expected points in that patch exactly exhaust the ideal.
            gb=S.groebner(equations,*remaining)
            expected=[p for p in points if p[i]!=0]
            if not expected:assert list(gb)==[1]
            else:
                assert len(expected)==1
                pp=expected[0];assert all(gb.reduce(q-S.Rational(pp[vv.index(q)],pp[i]))[1]==0 for q in remaining)
            patches.append(str(gb))
        for p in points:
            sub=dict(zip([x,y,z],p));assert all(g.subs(sub)==0 for g in grad)
            assert H.subs(sub).rank()==2
            assert (x*x+y*y+z*z).subs(sub)!=0
        nodal.append({'leading':str(f),'singular_points':points,'Hessian_rank':2,'nonzero_F2':'x^2+y^2+z^2','projective_patches':'PASS'})
    # F2(v)=0 is an actual boundary for the linear gradient lower bound.
    assert [S.diff(x*y*z,q).subs({x:0,y:0}) for q in [x,y,z]]==[0,0,0]
    # Exact full cubic prefix trap in Q(omega), omega^2+omega+1=0.
    cyclo=w*w+w+1;P=x**3+w*y**3+w*w*z**3
    moments=[S.rem(wick(S.expand(P**m),[x,y,z]),cyclo,w) for m in range(1,7)]
    assert moments==[0,0,0,0,0,96665400]
    # Pure-null negative control lies outside the smooth-leading hypothesis.
    null=(x+S.I*y)**3
    assert all(wick(S.expand(null**m),[x,y,z])==0 for m in range(1,9))
    vv={x:1,y:S.I,z:0};assert all(S.diff(null,q).subs(vv)==0 for q in [x,y,z])
    # Cumulant/root filters: exact regressions; the all-N proof is combinatorial.
    mu=[even_gaussian(3*m) for m in range(11)];kappa=[0]*11
    for m in range(1,11):kappa[m]=mu[m]-sum(comb(m-1,k-1)*kappa[k]*mu[m-k] for k in range(1,m))
    filters=[]
    for n in range(1,6):
        cp=S.cyclotomic_poly(2*n,w)
        assert all(S.rem(sum(w**(2*j*r) for j in range(n)),cp,w)==0 for r in range(1,n))
        assert S.rem(sum(w**(2*j*n) for j in range(n)),cp,w)==n
        assert all(kappa[m]==0 for m in range(1,2*n+1,2)) and kappa[2*n]>0
        filters.append({'variables':n,'first_nonzero_order':2*n,'first_nonzero_moment':n*kappa[2*n]})
    gates=[]
    for n,d in [(1,3),(2,3),(3,3),(4,3),(3,4)]:
        dim=(d-1)**n;B=n*(d-2)+2;gate=B*dim*(dim+1)//2+dim
        gates.append({'variables':n,'degree':d,'dimension':dim,'B':B,'sufficient_nonzero_prefix':gate})
    assert gates[2]['sufficient_nonzero_prefix']==188
    assert [v['sufficient_nonzero_prefix'] for v in gates]==[11,44,188,832,3051]
    p=Path(__file__).resolve().parent
    rec={'scope':'Exact Jacobian and twisted-connection algebra, nodal criteria and finite Wick regressions for the adjacent written proof. No analytic certification.',
         'universal_Hesse_division':universal,'connection_models':matrices,'nodal_examples':nodal,
         'ternary_smooth_prefix_trap':list(map(str,moments)),'all_N_root_filter_regressions':filters,'finite_gate_table':gates,
         'negative_controls':'lambda=1 singular Hesse; pure-null (x+i*y)^3; nodal xyz with F2=0',
         'analytic_boundary':'Tameness, local singularity resolution and full period comparison are written mathematics, not established by this checker. No formal verification or external review.',
         'sympy_version':S.__version__,'elapsed_seconds':time.monotonic()-started,
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (p/'tame_gaussian_verification.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps(rec,indent=2))


if __name__=='__main__':main()
