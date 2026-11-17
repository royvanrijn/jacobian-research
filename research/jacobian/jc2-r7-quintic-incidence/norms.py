"""Exact uniform finite-free norm of quintic nonimmersion."""
from pathlib import Path
import sys,json,time
import sympy as S
from sympy.polys.rings import ring
from sympy.polys.matrices import DomainMatrix
HERE=Path(__file__).resolve().parent
A,M,s=S.symbols('A M s')

def nonimmersion_norm(family):
    start=time.monotonic()
    row=next(r for r in json.loads((HERE.parent/'jc2-r7-quintic-normalizations/normalizations.json').read_text())['families'] if r['family']==family)
    par=S.Symbol(row['parameter']);loc={'A':A,'M':M,'s':s,str(par):par}
    parse=lambda x:S.sympify(x,locals=loc)
    q=parse(row['normalization_q']);b,c,e,f=[S.Poly(q,s).nth(j) for j in [4,3,2,1]]
    delta=S.expand(25*A**4+48*A**3*b**2-144*A**2*b*e+90*A**2*(f-c*A)+108*A*e**2+81*(f-c*A)**2)
    assert S.expand(S.resultant(3*s**2+A,S.diff(q,s),s)-delta)==0
    R,a,m=ring('A,M',S.QQ.poly_ring(par),order=lambda mon:(2*mon[0]+3*mon[1],mon[0]))
    basis=[R.from_expr(parse(v)) for v in row['groebner_basis']]
    standard=[tuple(v) for v in row['standard_monomials']];index={v:i for i,v in enumerate(standard)}
    element=R.from_expr(delta);matrix=[[S.S.Zero]*7 for _ in range(7)]
    for j,mon in enumerate(standard):
        rem=element.mul_monom(mon).div(basis)[1]
        for exponent,co in rem.items():matrix[index[exponent]][j]=co.as_expr()
    determinant=DomainMatrix.from_Matrix(S.Matrix(matrix)).det().as_expr()
    factored=S.factor(determinant)
    out=dict(family=family,parameter=str(par),delta=str(delta),matrix=[[str(v) for v in row] for row in matrix],norm=str(determinant),factorization=str(factored),elapsed_seconds=time.monotonic()-start,status='EXACT_DISCOVERY; determinant and interpretation require canonical replay')
    return out

def main():
    out=nonimmersion_norm(sys.argv[1]);(HERE/('nonimmersion-norm-'+out['family']+'.json')).write_text(json.dumps(out,indent=2)+'\n');print(out['factorization'],flush=True)

if __name__=='__main__':main()
