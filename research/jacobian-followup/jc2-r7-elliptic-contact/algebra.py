"""Exact key relation and coefficient cubic, reconstructed from the two charts."""
from pathlib import Path
import sys,json
import sympy as S
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE.parent/'jc2-r7-complement'))
import geometry as G
A,M,s,W,Y=S.symbols('A M s W Y')

def family(name):
    row=next(v for v in json.loads((HERE.parent/'jc2-r7-quintic-normalizations/normalizations.json').read_text())['families'] if v['family']==name)
    parameter=G.F.d if name=='L' else G.K.z
    parse=lambda value:S.sympify(value,locals={'A':A,'M':M,'s':s,str(parameter):parameter})
    p=G.F.p.subs(G.F.line) if name=='L' else G.K.P
    q=G.F.q.subs(G.F.line) if name=='L' else G.K.Q
    gc,g,_=G.approximate_quintic(p,q,G.F.t)
    residual=S.Poly(S.expand(g*g),G.F.t);key=W*W
    for degree,target,restriction in [(18,G.U**3,p**3),(16,G.U*G.V,p*q),(15,G.U*W,p*g),(12,G.U**2,p*p),(10,G.V,q),(9,W,g),(6,G.U,p),(0,S.S.One,S.S.One)]:
        coefficient=S.cancel(residual.nth(degree)/S.Poly(restriction,G.F.t).LC())
        residual=S.Poly(S.expand(residual.as_expr()-coefficient*restriction),G.F.t);key-=coefficient*target
    key=S.expand(key)
    f=parse(row['groebner_basis'][0]);r8=parse(row['necessary_sufficient_equations'][1])
    pb=parse(row['normalization_p']);qb=parse(row['normalization_q'])
    delta=S.expand(S.resultant(S.diff(pb,s),S.diff(qb,s),s))
    coefficient=S.Poly(f,M).nth(1);rest=S.Poly(f,M).nth(0)
    cubic=S.expand(coefficient*coefficient+36*rest)
    return dict(name=name,parameter=parameter,p=p,q=q,G=gc,g=g,key=key,key_residual=residual.as_expr(),
        f=f,r8=r8,P=pb,Q=qb,delta=delta,cubic=cubic,M_coefficient=coefficient)
