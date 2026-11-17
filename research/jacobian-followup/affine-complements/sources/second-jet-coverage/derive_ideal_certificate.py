"""Exact tracked polynomial division; produces V in a four-generator ideal.
No completion, locus coverage or nilpotence is inferred from a search prefix.
"""
from pathlib import Path
from sympy.polys.rings import ring
import sympy as S
import json,time,signal
HERE=Path(__file__).resolve().parent
start=time.perf_counter()
data=json.loads((HERE/'non_q_kernel_basis.json').read_text())
R,a,V,r=ring('a,V,rho',S.QQ,order='grevlex')
inputs=[R.from_expr(S.sympify(z,locals={'a':S.Symbol('a'),'V':S.Symbol('V'),'rho':S.Symbol('rho')})) for z in data['input']]
n=len(inputs);polys=[];reps=[];pairs=set()
def add(poly,rep):
 if not poly:return
 quot,rem=poly.div(polys)
 for j,quo in enumerate(quot):
  if quo:rep=[v-quo*w for v,w in zip(rep,reps[j])]
 if not rem:return
 lc=rem.LC;rem=rem/lc;rep=[v/lc for v in rep]
 assert sum((v*w for v,w in zip(rep,inputs)),R.zero)==rem
 idx=len(polys);polys.append(rem);reps.append(rep)
 for j in range(idx):
  lcm=R.monomial_lcm(rem.LM,polys[j].LM)
  if R.monomial_mul(rem.LM,polys[j].LM)!=lcm:pairs.add((j,idx))
def certificate():
 quot,rem=V.div(polys)
 if rem:return None
 rep=[R.zero]*n
 for j,quo in enumerate(quot):
  if quo:rep=[v+quo*w for v,w in zip(rep,reps[j])]
 assert sum((v*w for v,w in zip(rep,inputs)),R.zero)==V
 return rep
for i,poly in enumerate(inputs):add(poly,[R.one if i==j else R.zero for j in range(n)])
def alarm(*args):raise TimeoutError('45-second tracked-division cap')
signal.signal(signal.SIGALRM,alarm);signal.alarm(45)
steps=0;out=certificate()
while out is None and pairs:
 ij=min(pairs,key=lambda pair:R.order(R.monomial_lcm(polys[pair[0]].LM,polys[pair[1]].LM)));pairs.remove(ij);i,j=ij
 lcm=R.monomial_lcm(polys[i].LM,polys[j].LM)
 mi=R.from_dict({R.monomial_div(lcm,polys[i].LM):S.QQ.one});mj=R.from_dict({R.monomial_div(lcm,polys[j].LM):S.QQ.one})
 add(mi*polys[i]-mj*polys[j],[mi*u-mj*v for u,v in zip(reps[i],reps[j])]);steps+=1
 out=certificate()
signal.alarm(0)
assert out is not None
result=dict(target='V',variables=['a','V','rho'],input=data['input'],multipliers=[str(z.as_expr()) for z in out],identity_verified=True,intermediate_basis_size=len(polys),pair_steps=steps,elapsed_seconds=time.perf_counter()-start)
(HERE/'V_membership_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(identity_verified=True,elapsed_seconds=result['elapsed_seconds'],pair_steps=steps,basis_size=len(polys),certificate_bytes=(HERE/'V_membership_certificate.json').stat().st_size)))
