"""Bounded exact coefficient reconstruction for closest quintics."""
from pathlib import Path
import sys,json,time
import sympy as S
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'jc2-r7-complement'))
import geometry as G

def main():
    start=time.monotonic()
    s,A,M,b,c,e,f,n,al,be=S.symbols('s A M b c e f n alpha beta')
    family=sys.argv[1] if len(sys.argv)>1 else 'L'
    row=next(v for v in json.loads((G.HERE/'geometry.json').read_text())['closest_quintics'] if v['family']==family)
    parameter=G.F.d if family=='L' else G.K.z
    equation=S.sympify(row['equation'],locals={'U':G.U,'V':G.V,str(parameter):parameter})
    p=s**3+A*s+M;q=s**5+b*s**4+c*s**3+e*s**2+f*s+n
    value=S.Poly(S.expand(equation.subs({G.U:p,G.V:q})),s)
    recovered={};steps=[]
    for degree,var in [(14,b),(13,c),(12,e),(11,f),(10,n)]:
        coefficient=S.factor(value.nth(degree).subs(recovered))
        pol=S.Poly(coefficient,var)
        assert pol.degree()==1
        result=S.factor(-pol.nth(0)/pol.nth(1))
        recovered[var]=result;steps.append({'degree':degree,'variable':str(var),'value':str(result)})
        print('RECOVER',family,degree,var,result,flush=True)
    residual={str(k):str(S.factor(value.nth(k).subs(recovered))) for k in range(9,-1,-1)}
    out=dict(family=family,parameter=str(parameter),steps=steps,residual=residual,
             free=['A','M','alpha','beta'],status='EXACT_RECONSTRUCTION_ONLY; no exclusion asserted',elapsed_seconds=time.monotonic()-start)
    (HERE/('reconstruction-'+family+'.json')).write_text(json.dumps(out,indent=2)+'\n')
    for k,v in residual.items():print('RESIDUAL',k,v,flush=True)

if __name__=='__main__':main()
