"""Small exact one-variable branch replay after the quartic leading minors."""
from itertools import combinations
from pathlib import Path
import json,time
import sympy as S
from quartic_probe import R,u,v,w,z,t,remainder,determinant

HERE=Path(__file__).resolve().parent
start=time.monotonic(); output=[]
cases=json.loads((HERE.parent/'jc2-normal-r7/classification.json').read_text())['cases']
for i,row in enumerate(cases,1):
    p,q=S.sympify(row['p']),S.sympify(row['q'])
    raw=[remainder(f) for f in (q,p*p,p*q,q*q)]
    matrix=[[f[j]-f[6]*R.ground_new(S.Poly(p,t).nth(j)) for f in raw] for j in [1,2,3,4,5,7]]
    fixed_v=[S.Rational(-31,90),S.Rational(-157,90)][i-1]
    fixed_w=[S.Rational(-4,9),S.Rational(-8,45)][i-1]
    for label,variable,value in [('fixed_v',v,fixed_v),('fixed_w',w,fixed_w)]:
        parameter=S.Symbol('w' if label=='fixed_v' else 'v')
        submat=[[f.compose(u,R.zero).compose(variable,R.ground_new(value)) for f in row] for row in matrix]
        polys=[]; info=[]
        for rows in combinations(range(6),4):
            det=determinant([[submat[j][k] for k in range(4)] for j in rows])
            if not det:continue
            degree=max(m[-1] for m in det)
            for power in range(degree,-1,-1):
                co=R.from_dict({m[:-1]+(0,):c for m,c in det.items() if m[-1]==power})
                if co:
                    poly=S.Poly(co.as_expr(),parameter).clear_denoms()[1].primitive()[1]
                    polys.append(poly);info.append(dict(rows=list(rows),z_power=power,coefficient=str(co.as_expr())))
        common=polys[0]
        for f in polys[1:]:common=S.gcd(common,f)
        print('case',i,label,'gcd',S.factor(common.as_expr()),'polynomials',len(polys),'elapsed',round(time.monotonic()-start,2),flush=True)
        # Retain a pair whose gcd already equals the whole gcd if possible.
        pair=None
        for j in range(1,len(polys)):
            if S.gcd(polys[0],polys[j]).monic()==common.monic():
                pair=[info[0],info[j]];break
        output.append(dict(case=i,branch=label,constant=str(value),parameter=str(parameter),gcd=str(common.as_expr()),pair=pair))
(HERE/'quartic_branch_discovery.json').write_text(json.dumps(output,indent=2)+'\n')
