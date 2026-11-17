"""Narrow exact diagnostic: three second-jet witnesses on the full reduced trace.
Field-valued leading-face normal form, f=[p3]T3=0; not a scheme radical claim.
"""
from pathlib import Path
from itertools import permutations
import json,time
import sympy as S
from sympy.polys.rings import ring

def derive_witnesses():
    start=time.perf_counter()
    x,y,p,q=S.symbols('x y p q')
    kap,A,D,h,i,e,rho,a,b,c,d,g,j,k,l,m,n,C=S.symbols('kap A D h i e rho a b c d g j k l m n C')
    params=(kap,A,D,h,i,e,rho,a,b,c,d,g,j,k,l,m,n,C)
    T=kap*y**5+A*y**4+D*y**3*p+h*y**3+i*y*y*p+e*y*p*p+rho*y*y*q+a*p*q+b*y*q+c*q*q+d*p*p+g*y*p+j*y*y+k*p+l*q+m*y+n
    Ky=3*(2*p+13)*p+2*(21*y*y*p+24*y*q+2*S.diff(T,y));Kp=12*y*p+6*q+39*y+14*y**3+4*S.diff(T,p);Kq=6*p+39+24*y*y+4*S.diff(T,q)
    HH=30*y*y*p*p+60*y*p*q+18*q*q+78*y*y*p+156*y*q+12*y**4*p+36*y**3*q-12*y*T
    F=S.Poly(S.expand(-8*HH+8*Ky*Kq+4*Kp*Kp-64*y*Kp*Kq-(32*p-256*y*y)*Kq*Kq),y,p,q)
    normal=sum(co*y**ex[0]*p**ex[1]*q**ex[2] for ex,co in F.terms() if sum(ex)<=4)/64-C/64
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
    def pol(expr):return Ring.from_expr(S.expand(expr).subs(y,S.Symbol('Y')))
    def mul(P,Q):
     out={}
     for (ii,),aa in P.items():
      for (jj,),bb in Q.items():
       if ii+jj<=7:out[(ii+jj,)]=out.get((ii+jj,),Domain.zero)+aa*bb
     return Ring.from_dict(out)

    def jets(H):
     return [[tuple(pol(S.diff(H[u,v],zz).subs({p:0,q:0})) if zz else pol(H[u,v].subs({p:0,q:0})) for zz in (None,p,q)) for v in range(4)] for u in range(4)]
    Hs=[jets(H0),jets(H1)]
    ans=[Ring.zero,Ring.zero,Ring.zero]
    for perm in permutations(range(4)):
     sign=(-1)**sum(perm[u]>perm[v] for u in range(4) for v in range(u+1,4))
     for chosen in range(4):
      acc=[Ring.one,Ring.zero,Ring.zero]
      for u in range(4):
       val=Hs[int(u==chosen)][u][perm[u]]
       acc=[mul(acc[0],val[0]),mul(acc[1],val[0])+mul(acc[0],val[1]),mul(acc[2],val[0])+mul(acc[0],val[2])]
      ans=[aa+sign*bb for aa,bb in zip(ans,acc)]
    coeffs=[S.factor(Domain.to_sympy(poly.get((deg,),Domain.zero))) for poly,deg in zip(ans,(7,5,4))]
    return dict(witnesses=coeffs,params=params,source_symbols=(x,y,p,q),trace=T,plane=F,normal=normal,potential_coefficients=fs,elapsed_seconds=time.perf_counter()-start)

if __name__=='__main__':
    result=derive_witnesses()
    names=('y7','y5p','y4q')
    out=dict(field_trace_f_p3_set_to_zero=True,parameters=[str(z) for z in result['params']],witnesses={name:str(coef) for name,coef in zip(names,result['witnesses'])},elapsed_seconds=result['elapsed_seconds'])
    Path(__file__).with_name('witnesses.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(outcome='EXACT_DERIVATION',elapsed_seconds=result['elapsed_seconds'])))
