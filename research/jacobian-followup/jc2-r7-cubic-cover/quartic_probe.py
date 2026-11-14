"""Capped discovery probe: doubled quartic-fibre contact minors.

Necessary algebra only; no rank conclusion is promoted without replay.
Run under an external 60-second timeout. Output is a small coefficient
ideal for later exact certificates, never a Keller-source assertion.
"""
from pathlib import Path
from itertools import combinations,permutations
import json,time
import sympy as S
from sympy.polys.rings import ring

HERE=Path(__file__).resolve().parent
t=S.Symbol('t')
R,u,v,w,z=ring('u,v,w,z',S.QQ)
lower=[z*z,-2*w*z,w*w-2*v*z,2*v*w-2*u*z,
       v*v+2*u*w-2*z,2*u*v+2*w,u*u+2*v,2*u]
remainders=[]
for n in range(21):
    if n<8: remainders.append([R.one if i==n else R.zero for i in range(8)])
    else:
        remainders.append([-sum((lower[i]*remainders[n-8+i][j] for i in range(8)),R.zero) for j in range(8)])

def remainder(f):
    p=S.Poly(f,t)
    return [sum((R.ground_new(p.nth(n))*remainders[n][j] for n in range(p.degree()+1)),R.zero) for j in range(8)]

def determinant(M):
    result=R.zero
    for perm in permutations(range(4)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(4) for j in range(i+1,4))
        term=R.ground_new(sign)
        for i,j in enumerate(perm): term*=M[i][j]
        result+=term
    return result

if __name__=='__main__':
    start=time.monotonic()
    cases=json.loads((HERE.parent/'jc2-normal-r7/classification.json').read_text())['cases']
    output=[]
    for i,row in enumerate(cases,1):
        p,q=S.sympify(row['p']),S.sympify(row['q'])
        raw=[remainder(f) for f in (q,p*p,p*q,q*q)]
        matrix=[[f[j]-f[6]*R.ground_new(S.Poly(p,t).nth(j)) for f in raw] for j in [1,2,3,4,5,7]]
        coefficients=[]
        for rows in combinations(range(6),4):
            det=determinant([[matrix[j][k] for k in range(4)] for j in rows])
            if not det: continue
            degree=max(m[-1] for m in det)
            for power in range(degree,max(-1,degree-2),-1):
                co=R.from_dict({m[:-1]+(0,):c for m,c in det.items() if m[-1]==power})
                if co:
                    expression=S.factor(co.as_expr())
                    coefficients.append(dict(rows=list(rows),z_power=power,expression=str(expression)))
            print('case',i,'rows',rows,'zdegree',degree,'terms',len(det),'elapsed',round(time.monotonic()-start,2),flush=True)
        output.append(dict(case=i,coefficients=coefficients))
        (HERE/'quartic_discovery.json').write_text(json.dumps(output,indent=2)+'\n')
    print('Saved necessary coefficient data; no exclusion inferred.',flush=True)
