"""Exact simultaneous incidence witnesses in monic finite-free algebras."""
from pathlib import Path
import sys,json,time,hashlib
import sympy as S
from finite import Algebra,rank,modular
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'jc2-r7-complement'))
import geometry as G
A,M,s,v,t=S.symbols('A M s v t')

def witness(family):
    start=time.monotonic();prime=101
    row=next(v for v in json.loads((HERE.parent/'jc2-r7-quintic-normalizations/normalizations.json').read_text())['families'] if v['family']==family)
    par=G.F.d if family=='L' else G.K.z
    loc=dict(A=A,M=M,s=s,**{str(par):par})
    parse=lambda x:S.sympify(x,locals=loc).subs(par,1)
    h=parse(row['witness_univariate']);mf=S.expand(M-parse(row['witness_linear']))
    def red(expr):return S.rem(S.expand(expr.subs(M,mf)),h,A)
    pc=(G.F.p.subs(G.F.line).subs(G.F.d,1) if family=='L' else G.K.P.subs(G.K.z,1));qc=(G.F.q.subs(G.F.line).subs(G.F.d,1) if family=='L' else G.K.Q.subs(G.K.z,1))
    pc=pc.subs(G.F.t,t);qc=qc.subs(G.F.t,t)
    other=S.Symbol('other')
    divided_p=S.cancel((pc-pc.subs(t,other))/(t-other));divided_q=S.cancel((qc-qc.subs(t,other))/(t-other))
    dc=S.Poly(S.resultant(divided_p,divided_q,other),t).monic().as_expr()
    assert S.degree(dc,t)==24 and S.gcd(dc,S.diff(dc,t))==1
    assert S.gcd(S.diff(pc,t),S.diff(qc,t))==1
    assert all(red(parse(e))==0 for e in row['necessary_sufficient_equations'])
    pb=parse(row['normalization_p']);qb=parse(row['normalization_q']);b,c,e,f,n=[S.Poly(qb,s).nth(j) for j in [4,3,2,1,0]]
    R=red(v**4+b*v**3+A*v**2+(2*A*b-e)*v-(A*A+f-c*A))
    ub=red(M-v*(v*v+A));vb=red((v*v+A)*(v**3+2*A*v+A*b-e)+c*(M-v*(v*v+A)-M)+n)
    k={mon:parse(x) for mon,x in row['implicit_difference'].items()};_,rc,_=G.approximate_quintic(pc,qc,t)
    finite=red(rc-k['(2, 0)']*pc**2-k['(0, 1)']*qc-k['(1, 0)']*pc-k['(0, 0)'])
    cases=[('seven_choices_reduced',(A,),(7,),[h],S.diff(h,A)),
           ('B_pair_roots_simple',(A,v),(7,4),[h,R],S.diff(R,v)),
           ('B_avoids_diagonal_pairs',(A,v),(7,4),[h,R],3*v*v+4*A),
           ('B_has_no_triple_value',(A,),(7,),[h],red(A*A+A*b*b-b*e-A*c+f)),
           ('B_avoids_Sing_C',(A,t),(7,24),[h,dc],finite),
           ('C_avoids_Sing_B',(A,v,t),(7,4,6),[h,R,pc-ub],qc-vb),
           ('intersection_degree12',(A,),(7,),[h],red(S.Poly(finite,t).LC()))]
    result=[]
    for name,variables,degrees,relations,element in cases:
        algebra=Algebra(variables,degrees,relations,prime);matrix=algebra.matrix(element);r,pivots=rank(matrix,prime)
        item=dict(name=name,dimension=len(matrix),rank=r,prime=prime,matrix_sha256=hashlib.sha256(json.dumps(matrix,separators=(',',':')).encode()).hexdigest());result.append(item)
    hp=S.Poly(sum(modular(co,prime)*A**m[0] for m,co in S.Poly(h,A).terms()),A,modulus=prime)
    field=lambda value:S.Poly(sum(modular(co,prime)*A**m[0] for m,co in S.Poly(value,A).terms()),A,modulus=prime)
    fp=S.Poly(finite,t);inverse=S.invert(field(fp.LC()),hp)
    normalized=sum((field(co)*inverse).rem(hp).as_expr()*t**m[0] for m,co in fp.terms())
    algebra=Algebra((A,t),(7,12),[h,normalized],prime);matrix=algebra.matrix(S.diff(normalized,t));r,pivots=rank(matrix,prime)
    item=dict(name='twelve_simple_intersections',dimension=len(matrix),rank=r,prime=prime,matrix_sha256=hashlib.sha256(json.dumps(matrix,separators=(',',':')).encode()).hexdigest());result.append(item)
    return dict(family=family,parameter=1,checks=result,elapsed_seconds=time.monotonic()-start,status='Exact modular witness; characteristic-zero interpretation uses explicitly monic finite-free algebras.')

def main():
    result=witness(sys.argv[1]);(HERE/('ranks-'+result['family']+'.json')).write_text(json.dumps(result,indent=2)+'\n');print(result,flush=True)

if __name__=='__main__':main()
