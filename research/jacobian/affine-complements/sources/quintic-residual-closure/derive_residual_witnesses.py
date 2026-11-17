"""Exact residual jet coefficients, keeping every lower trace parameter and f.
The two conjugates are treated together modulo Q(rho), not sampled.
"""
from pathlib import Path

import tempfile as _check_tempfile

_CHECK_OUTPUT_DIRECTORY = None

def _check_output(original):
    """Write new runs to a temporary directory; preserve deposited receipts."""
    global _CHECK_OUTPUT_DIRECTORY
    if _CHECK_OUTPUT_DIRECTORY is None:
        _CHECK_OUTPUT_DIRECTORY = Path(_check_tempfile.mkdtemp(prefix="math-check-"))
        print("Verification output:", _CHECK_OUTPUT_DIRECTORY)
    return _CHECK_OUTPUT_DIRECTORY / Path(original).name


from itertools import permutations
import json,time
import sympy as S
from sympy.polys.rings import ring

def derive_residual_witnesses():
    start=time.perf_counter()
    x,y,p,q=S.symbols('x y p q')
    rho,h,i,e,f,a,b,c,d,g,j,k,l,m,n,C=S.symbols('rho h i e f a b c d g j k l m n C')
    params=(rho,h,i,e,f,a,b,c,d,g,j,k,l,m,n,C)
    D=(16*rho+89)/2
    T=D*y**3*p+h*y**3+i*y*y*p+e*y*p*p+f*p**3+rho*y*y*q+a*p*q+b*y*q+c*q*q+d*p*p+g*y*p+j*y*y+k*p+l*q+m*y+n
    Ky=3*(2*p+13)*p+2*(21*y*y*p+24*y*q+2*S.diff(T,y))
    Kp=12*y*p+6*q+39*y+14*y**3+4*S.diff(T,p)
    Kq=6*p+39+24*y*y+4*S.diff(T,q)
    HH=30*y*y*p*p+60*y*p*q+18*q*q+78*y*y*p+156*y*q+12*y**4*p+36*y**3*q-12*y*T
    plane=S.Poly(S.expand(-8*HH+8*Ky*Kq+4*Kp*Kp-64*y*Kp*Kq-(32*p-256*y*y)*Kq*Kq),y,p,q)
    normal=S.expand(sum(co*y**ex[0]*p**ex[1]*q**ex[2] for ex,co in plane.terms() if sum(ex)<=4)/64-C/64)
    S4=S.expand(sum(co*y**ex[0]*p**ex[1]*q**ex[2] for ex,co in S.Poly(normal,y,p,q).terms() if sum(ex)==4))
    fs=[p*p+13*p+8*y*y*p+2*y*q,
        6*y*p*p+6*p*q+39*y*p+39*q+14*y**3*p+24*y*y*q+4*T,
        15*y*y*p*p+30*y*p*q+9*q*q+39*y*y*p+78*y*q+6*y**4*p+18*y**3*q-6*y*T+4*normal,
        20*y**3*p*p+60*y*y*p*q+36*y*q*q+13*y**3*p+39*y*y*q-(2*p+13)*T-6*y*normal]
    H0=S.zeros(4);H1=S.zeros(4);vs=(y,p,q)
    H0[0,0]=2*fs[2];H1[0,0]=6*fs[3]
    for u,vu in enumerate(vs,1):
        H0[0,u]=H0[u,0]=S.diff(fs[1],vu);H1[0,u]=H1[u,0]=2*S.diff(fs[2],vu)
        for v,vv in enumerate(vs,1):
            H0[u,v]=S.diff(fs[0],vu,vv);H1[u,v]=S.diff(fs[1],vu,vv)
    Domain=S.QQ.poly_ring(*params);Ring,Y=ring('Y',Domain)
    def extract(active,order,degree):
        def poly(expr):return Ring.from_expr(S.expand(expr).subs(y,S.Symbol('Y')))
        def mul(P,Q):
            out={}
            for (ii,),aa in P.items():
                for (jj,),bb in Q.items():
                    if ii+jj<=degree:out[(ii+jj,)]=out.get((ii+jj,),Domain.zero)+aa*bb
            return Ring.from_dict(out)
        def jets(H):
            return [[tuple(poly(S.diff(H[u,v],active,deg).subs({p:0,q:0})/S.factorial(deg)) for deg in range(order+1)) for v in range(4)] for u in range(4)]
        Hs=[jets(H0),jets(H1)];ans=[Ring.zero]*(order+1)
        for perm in permutations(range(4)):
            sign=(-1)**sum(perm[u]>perm[v] for u in range(4) for v in range(u+1,4))
            for chosen in range(4):
                acc=[Ring.one]+[Ring.zero]*order
                for u in range(4):
                    vals=Hs[int(u==chosen)][u][perm[u]]
                    acc=[sum((mul(acc[z],vals[deg-z]) for z in range(deg+1)),Ring.zero) for deg in range(order+1)]
                ans=[aa+sign*bb for aa,bb in zip(ans,acc)]
        return S.factor(Domain.to_sympy(ans[order].get((degree,),Domain.zero)))
    wc=extract(q,2,4);wp=extract(p,1,5)
    return dict(source_symbols=(x,y,p,q),params=params,trace=T,plane=plane,normal=normal,S4=S4,potential_coefficients=fs,
                wc=wc,wp=wp,elapsed_seconds=time.perf_counter()-start)

if __name__=='__main__':
    z=derive_residual_witnesses();rho,h,i,e,f,a,b,c,d,g,j,k,l,m,n,C=z['params']
    Q=160*rho*rho+1968*rho+6021;U=S.Rational(96,25)*(2898*rho+16821)
    wp=z['wp'];wpc0=S.expand(wp.subs(c,0));wpcf0=S.expand(wpc0.subs(f,0))
    co_c=S.cancel((wp-wpc0)/c);co_f=S.cancel((wpc0-wpcf0)/f);co_Q=S.cancel((wpcf0-U)/Q)
    for coef in (co_c,co_f,co_Q):assert S.denom(coef)==1 or S.denom(coef).is_Rational
    assert S.expand(wp-U-c*co_c-f*co_f-Q*co_Q)==0
    out=dict(parameters=[str(v) for v in z['params']],wc=str(z['wc']),wp=str(wp),unit_residue=str(U),
             coefficient_c=str(co_c),coefficient_f=str(co_f),coefficient_Q=str(co_Q),elapsed_seconds=z['elapsed_seconds'])
    _check_output(Path(__file__).with_name('residual_witness_certificate.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(wc=str(z['wc']),wp_terms=len(S.Poly(wp,*z['params']).terms()),elapsed_seconds=z['elapsed_seconds'],decomposition_exact=True)))
