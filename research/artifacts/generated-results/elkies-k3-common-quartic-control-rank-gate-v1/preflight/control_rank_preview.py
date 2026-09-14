import json,time,resource
from pathlib import Path
import sympy as S
resource.setrlimit(resource.RLIMIT_CPU,(30,35))
start=time.process_time(); t=S.symbols('t')
pkt=json.load(open('research/artifacts/generated-results/elkies-k3-common-quartic-singularities-v1/input.json'))
rows=[]
for c in pkt['controls']:
 D=sum(v*t**i for i,v in enumerate(c['D'])); s=sum(v*t**i for i,v in enumerate(c['s']))
 A=S.Poly(D*(2*s+1)-1,t); B=S.Poly(D*s*s,t); delta=4*A**3+27*B**2
 for p in [5,7,11,13,17,19,23,29,31,37,41,43]:
  d=S.Poly(delta,t,modulus=p)
  if d.degree()!=24 or S.gcd(d,d.diff()).degree()!=0: continue
  aa=[int(A.nth(i))%p for i in range(9)]; bb=[int(B.nth(i))%p for i in range(13)]
  count=[]
  for u in list(range(p))+[None]:
   a=aa[-1] if u is None else sum(v*pow(u,i,p) for i,v in enumerate(aa))%p
   b=bb[-1] if u is None else sum(v*pow(u,i,p) for i,v in enumerate(bb))%p
   n=1+sum(sum((y*y-x*x*x-a*x-b)%p==0 for y in range(p)) for x in range(p))
   count.append(n)
  trace=sum(count)-1-p*p
  bound=(22*p+trace)//(2*p)
  rows.append({'control':c['name'],'p':p,'fibre_counts':count,'surface_count':sum(count),'trace':trace,'rational_picard_upper':bound,'any_Q_fibration_rank_upper':bound-2})
  break
assert len(rows)==2
result={'scope':'Two fixed control parents only; first squarefree degree24 discriminant prime from5,7,11,13,17,19,23,29,31,37,41,43. CPU30 seconds, no rank or point search. This preview is retained separately from the subsequent frozen independent certificate.','rows':rows,'cpu_seconds':time.process_time()-start}
path=Path('research/artifacts/generated-results/elkies-k3-common-quartic-control-rank-gate-v1')
path.mkdir(parents=True,exist_ok=True)
with (path/'preview.json').open('x') as f: json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result))
