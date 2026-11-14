"""Exact three-chart controls for the mixed-root pole-(3,5) layer.

The written coverage and generic-fibre arguments are separate. These
coefficient equations alone are neither a finite map nor a full proof.
"""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib
HERE=Path(__file__).resolve().parent
s,t,w,z=S.symbols('s t w z')
t0,x0,y0,k,lam=S.symbols('t0 x0 y0 k lam')
n1,delta,d0,d1,d2=S.symbols('n1 delta d0 d1 d2')
c1,c2,c3,c4=S.symbols('c1 c2 c3 c4')
D=d0*s*s+d1*s*t+d2*t*t

def lift(f,const=0):
    p=S.Poly(S.expand(f.subs({s:w,t:1})),w)
    return w*p.as_expr()/z+p.nth(2)*w+2*p.nth(3)*w*w+p.nth(3)*z+const

g=[w/z,w*w/z,(z+w*w)*w/z,(z+w*w)**2/z]
def case(name,L,M):
    tau=L*L*M;u=n1*L*M*M;v=delta*n1*n1*M**3+L*D
    T=lift(tau,t0);X=lift(u,x0);Y=lift(v,y0)
    W=lift(n1*n1*L*M*M)*lift(n1*M**3)
    C=sum(c*gg for c,gg in zip((c1,c2,c3,c4),g))
    P=T**3+3*T*X+3*Y
    G5=T**5+5*T**3*X+5*T*T*Y+5*T*X*X+10*X*Y
    G4=T**4+4*T*T*X+4*T*Y+2*X*X
    Q=G5+k*G4-S.Rational(5,3)*W+lam*T*T+C
    pp=S.Poly(S.expand(P*z**3),z);qq=S.Poly(S.expand(Q*z**5),z)
    coeff={}
    for j in range(5):
        coeff[j]=S.Poly(S.expand(sum((i-3)*pp.nth(i)*S.diff(qq.nth(j-i),w)
                 -(j-i-5)*S.diff(pp.nth(i),w)*qq.nth(j-i) for i in range(j+1))),w)
    for j in range(4):assert coeff[j].is_zero,(name,j)
    f5=coeff[4]
    if L==s:
        functional=f5.nth(7);label='coefficient_w7'
    elif L==t:
        functional=f5.nth(11);label='coefficient_w11'
    else:
        for j in range(7):
            val=S.factor(S.diff(f5.as_expr(),w,j).subs(w,1))
            print('at1 derivative',j,val)
            if val!=0:
                functional=val;label='derivative_'+str(j)+'_at_1';break
    print(name,label,S.factor(functional))
    return dict(name=name,leading_sections=list(map(str,(tau,u,v))),
        W=str(W),coefficient_equations=[(d[0],str(S.factor(v))) for d,v in f5.terms() if v!=0],
        functional=label,functional_value=str(S.factor(functional)))
start=monotonic()
rows=[case('double_at_marked',s,t),case('simple_at_marked',t,s),
      case('neither_at_marked',s-t,t)]
(HERE/'pole35_repeated_root_equations.json').write_text(json.dumps(dict(
    as_of='2026-11-14',outcome='NECESSARY_EQUATIONS_GENERATED',cases=rows,
    nonproportionality_open='n1 != 0',
    assurance=dict(coverage_formally_verified=False,layer_excluded=False,finite_map_constructed=False),
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    elapsed_seconds=round(monotonic()-start,6)),indent=2)+'\n')
print('seconds',monotonic()-start)
