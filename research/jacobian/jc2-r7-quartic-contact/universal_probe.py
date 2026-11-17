"""Weighted exact coefficient extraction for universal quartic contact minors.

The determinant coefficient has a fixed parameter weight. Ignore only
entry coefficients whose nonnegative weight already exceeds that total.
This is an exact extraction, not a parameter sample or series heuristic.
Discovery use: run with an external 60-second cap.
"""
from itertools import permutations,product
from pathlib import Path
import json,sys,time
import sympy as S
from sympy.polys.rings import ring

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'foundations/scripts'))
from verify_f2_degree_6_10_gap import construct_obstructions
from contact import R as REM,u as ru,v as rv,w as rw,z as rz,power_remainders

P,*gens=ring('a,b,c,d,B,D,E,u,v,w',S.QQ)
a,b,c,d,B,D,E,u,v,w=gens
t=S.Symbol('t')
weights=(2,3,4,5,2,6,8,1,2,3)
degrees=(10,12,16,20)
row_powers=(1,2,3,4,5,7)
_,rho,_,_,powers=construct_obstructions()
pp=t**6+a.as_expr()*t**4+b.as_expr()*t**3+c.as_expr()*t**2+d.as_expr()*t
qq=S.expand(sum(co*sum(powers[j][n].as_expr()*t**(2*j-n) for n in range(2*j+1))
                for j,co in [(5,1),(4,B.as_expr()),(2,D.as_expr()),(1,E.as_expr())]))
polynomials=[S.Poly(f,t) for f in (qq,pp*pp,pp*qq,qq*qq)]
columns=[[P.from_expr(f.nth(j)) for j in range(f.degree()+1)] for f in polynomials]
remainders=power_remainders(20)


def remainder_coefficient(n,r,k):
    result=P.zero
    for monomial,coefficient in remainders[n][r].items():
        if monomial[-1]!=k:continue
        exponents=(0,)*7+monomial[:3]
        result+=P.from_dict({exponents:coefficient})
    return result


def matrix_entry(column,r,k,total):
    f=columns[column]; deg=degrees[column]
    result=sum((f[n]*remainder_coefficient(n,r,k)
                for n in range(max(0,deg-total),deg+1)),P.zero)
    pr=P.from_expr(S.Poly(pp,t).nth(r))
    if pr:
        result-=pr*sum((f[n]*remainder_coefficient(n,6,k)
                       for n in range(max(0,deg-total),deg+1)),P.zero)
    return result


def coefficient(rows,power,substitutions=None):
    total=sum(degrees)-sum(row_powers[r] for r in rows)-4*power
    entries=[]
    for row in rows:
        r=row_powers[row]; row_entries=[]
        for col,deg in enumerate(degrees):
            values={k:matrix_entry(col,r,k,total) for k in range(6)
                    if 0<=deg-r-4*k<=total}
            if substitutions:
                values={k:f.compose(substitutions) for k,f in values.items()}
            row_entries.append({k:f for k,f in values.items() if f})
        entries.append(row_entries)
    answer=P.zero
    for perm in permutations(range(4)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        factors=[entries[i][j] for i,j in enumerate(perm)]
        for chosen in product(*(list(f.items()) for f in factors)):
            if sum(k for k,f in chosen)!=power:continue
            term=P.ground_new(sign)
            for k,f in chosen:term*=f
            answer+=term
    preserves_weight=not substitutions or all(
        all(sum(weight*exponent for weight,exponent in zip(weights,m))==weights[gens.index(variable)]
            for m in value)
        for variable,value in substitutions.items())
    if preserves_weight:
        assert all(sum(weight*exponent for weight,exponent in zip(weights,m))==total for m in answer)
    return answer


if __name__=='__main__':
    start=time.monotonic()
    targets=[((0,2,3,4),10,None),((0,1,2,3),10,{u:P.zero}),((0,1,2,3),9,{u:P.zero})]
    results=[]
    for rows,power,subs in targets:
        out=coefficient(rows,power,subs)
        factored=S.factor(out.as_expr())
        print('rows',rows,'Z',power,'terms',len(out),'elapsed',round(time.monotonic()-start,3),flush=True)
        print(factored,flush=True)
        results.append(dict(rows=list(rows),z_power=power,u_zero=bool(subs),expression=str(factored)))
        (HERE/'universal_discovery.json').write_text(json.dumps(results,indent=2)+'\n')
