"""Prepare one independently generated complete passive-block certificate.

The whole 35-coefficient class is invariant under GL2(z,w), so nondegenerate
A1 may be put into I2 without restricting a previously normalized stratum.
"""
from pathlib import Path
import json
import time
from sage.all import PolynomialRing,QQ

HERE=Path(__file__).resolve().parent
data=json.loads((HERE/'zquad-input.json').read_text())
names=data['variables'];R=PolynomialRing(QQ,names,order='degrevlex');rv=dict(zip(names,R.gens()))
indices={tuple(e):i for i,e in enumerate(data['monomials'])}
fixed={f'c{indices[1,0,2,0]}':QQ(1)/2,
       f'c{indices[1,0,1,1]}':QQ(0),
       f'c{indices[1,0,0,2]}':QQ(1)/2}
N=PolynomialRing(QQ,[n for n in names if n not in fixed and n not in ('T','U')],order='degrevlex')
sub={rv[k]:R(v) for k,v in fixed.items()}
equations=[]
for raw in data['equations'][:-2]:
    p=N(str(R(raw).subs(sub)))
    if p!=0 and p not in equations:equations.append(p)
c=N(str(R(data['constant_determinant']).subs(sub)))
payload=dict(variables=list(N.variable_names()),equations=list(map(str,equations)),target=str(c),normalization={k:str(v) for k,v in fixed.items()},scope='Complete GL2-invariant unary-top Z class, A1=I2; no leaf-specific substitutions')
(HERE/'zquad-normal-input.json').write_text(json.dumps(payload,indent=2)+'\n')
start=time.monotonic()
I=N.ideal(equations);G=I.groebner_basis(algorithm='libsingular:slimgb')
powers=[i for i in range(1,5) if (c**i).reduce(G)==0]
result=dict(basis_size=len(G),target_powers_checked=list(range(1,5)),zero_reductions=powers,elapsed_seconds=round(time.monotonic()-start,3))
(HERE/'zquad-normal-result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
if powers:
    power=min(powers)
    text='ring R=0,('+','.join(N.variable_names())+'),dp;\n'
    text+='ideal I='+','.join(map(str,equations))+';\n'
    text+='poly target=('+str(c)+')^'+str(power)+';\n'
    text+='matrix L=lift(I,ideal(target));\n'
    text+='print("LIFT_BEGIN");\n'
    for i in range(len(equations)):
        text+=f'print("ROW {i}");print(L[{i+1},1]);\n'
    text+='print("LIFT_END");quit;\n'
    (HERE/'zquad-normal-lift.sing').write_text(text)
