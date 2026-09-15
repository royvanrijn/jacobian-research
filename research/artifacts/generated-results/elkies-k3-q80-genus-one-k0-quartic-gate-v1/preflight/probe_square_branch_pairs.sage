from sage.all import GF, PolynomialRing
import json, resource, time
from pathlib import Path
resource.setrlimit(resource.RLIMIT_CPU,(30,35));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
R=Path(__file__).resolve().parents[4]
points=json.loads((R/'artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1/norm4-sections.json').read_text())['records']
F=GF(131);P=PolynomialRing(F,'t');t=P.gen()
x=[P(r['x']) for r in points];y=[P(r['y']) for r in points]
start=time.process_time();found=[];shared2=0;nodal=[]
q0=t*t+62*t+88
for i in range(len(points)):
 if y[i]%q0==0:nodal.append({'index':i,'q0_order':int(y[i].valuation(q0)),'x_mod_q0':[int(c) for c in (x[i]%q0).list()]})
 for j in range(i):
  g=y[i].gcd(y[j]);inf=min(6-y[i].degree(),6-y[j].degree())
  if g.degree()+inf<2:continue
  shared2+=1;diff=x[i]-x[j]
  while True:
   c=g.gcd(diff)
   if c.degree()<=0:break
   g//=c
  if x[i][4]==x[j][4]:inf=0
  if g.degree()+inf<2:continue
  fac=[(f,e) for f,e in g.factor() if f.degree()<=2]
  quadratics=[]
  for f,e in fac:
   if f.degree()==2:quadratics.append([int(c) for c in f.list()])
   elif e>=2:quadratics.append([int(c) for c in (f*f).list()])
  linear=[f for f,e in fac if f.degree()==1]
  for a in range(len(linear)):
   for b in range(a):quadratics.append([int(c) for c in (linear[a]*linear[b]).list()])
  if inf:
   quadratics += [[int(c) for c in f.list()]+[0] for f in linear]
   if inf>=2:quadratics.append([1,0,0])
  if quadratics:found.append({'pair':[j,i],'H':quadratics})
print(json.dumps({'pairs_examined':len(points)*(len(points)-1)//2,'shared_degree_at_least2':shared2,
                  'distinct_root_contact_pairs':found,'sections_with_q0_ordinate_factor':nodal,
                  'cpu_seconds':time.process_time()-start},indent=2))
