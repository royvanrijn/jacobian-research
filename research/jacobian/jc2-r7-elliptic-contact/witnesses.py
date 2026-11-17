"""Exact contact34 witnesses on the coefficient cubic; no root embeddings."""
from pathlib import Path
import sys,json,hashlib
import sympy as S
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'jc2-r7-quintic-incidence'))
from finite import Algebra,rank,modular
sys.path.insert(0,str(HERE.parent/'jc2-r7-complement'))
import geometry as G
A,M,s,t,v=S.symbols('A M s t v')

def witness(family):
    prime=101
    row=next(v for v in json.loads((HERE.parent/'jc2-r7-quintic-normalizations/normalizations.json').read_text())['families'] if v['family']==family)
    par=S.Symbol(row['parameter']);loc=dict(A=A,M=M,s=s,**{str(par):par})
    parse=lambda value:S.sympify(value,locals=loc).subs({par:1,A:0})
    h=S.Poly(parse(row['groebner_basis'][0]),M).monic().as_expr()
    pb=parse(row['normalization_p']);qb=parse(row['normalization_q'])
    b,c,e,f,n=[S.Poly(qb,s).nth(j) for j in [4,3,2,1,0]]
    pc=(G.F.p.subs(G.F.line).subs(G.F.d,1) if family=='L' else G.K.P.subs(G.K.z,1)).subs(G.F.t,t)
    qc=(G.F.q.subs(G.F.line).subs(G.F.d,1) if family=='L' else G.K.Q.subs(G.K.z,1)).subs(G.F.t,t)
    other=S.Symbol('other');dp=S.cancel((pc-pc.subs(t,other))/(t-other));dq=S.cancel((qc-qc.subs(t,other))/(t-other))
    dc=S.Poly(S.resultant(dp,dq,other),t).monic().as_expr()
    _,rc,_=G.approximate_quintic(pc,qc,t);k={mon:parse(expr) for mon,expr in row['implicit_difference'].items()}
    finite=S.expand(rc-k['(1, 1)']*pc*qc-k['(2, 0)']*pc*pc-k['(0, 1)']*qc-k['(1, 0)']*pc-k['(0, 0)'])
    finite=S.rem(finite,h,M)
    R=v**4+b*v**3-e*v-f;ub=M-v**3;vb=v*v*(v**3-e)+n
    cases=[('two_reduced_coefficient_choices',(M,),(2,),[h],S.diff(h,M)),
        ('contact_exactly34',(M,),(2,),[h],k['(1, 1)']),
        ('B_collision_roots_simple',(M,v),(2,4),[h,R],S.diff(R,v)),
        ('B_immersive',(M,v),(2,4),[h,R],3*v*v),
        ('B_has_no_triple_value',(M,),(2,),[h],f-b*e),
        ('B_avoids_Sing_C',(M,t),(2,24),[h,dc],finite),
        ('C_avoids_Sing_B',(M,v,t),(2,4,6),[h,R,pc-ub],qc-vb)]
    outputs=[]
    for name,variables,degrees,relations,element in cases:
        algebra=Algebra(variables,degrees,relations,prime);matrix=algebra.matrix(element);r,_=rank(matrix,prime)
        outputs.append(dict(name=name,dimension=len(matrix),rank=r,prime=prime,matrix_sha256=hashlib.sha256(json.dumps(matrix,separators=(',',':')).encode()).hexdigest()))
    tofield=lambda expr:S.Poly(sum(modular(co,prime)*M**mon[0] for mon,co in S.Poly(expr,M).terms()),M,modulus=prime)
    hp=tofield(h);fp=S.Poly(finite,t);inverse=S.invert(tofield(fp.LC()),hp)
    normalized=sum((tofield(co)*inverse).rem(hp).as_expr()*t**mon[0] for mon,co in fp.terms())
    algebra=Algebra((M,t),(2,16),[h,normalized],prime);matrix=algebra.matrix(S.diff(normalized,t));r,_=rank(matrix,prime)
    outputs.append(dict(name='sixteen_transverse_intersections',dimension=32,rank=r,prime=prime,matrix_sha256=hashlib.sha256(json.dumps(matrix,separators=(',',':')).encode()).hexdigest()))
    return dict(family=family,first_parameter=1,A=0,M_equation=str(h),normalization_p=str(pb),normalization_q=str(qb),checks=outputs)

if __name__=='__main__':
    for family in ['L','K']:
        result=witness(family);print(result)
