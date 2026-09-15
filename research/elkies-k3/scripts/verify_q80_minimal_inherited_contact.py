"""Exact bounded contact check; no higher-height or all-trace exclusion."""
import json,time,resource,hashlib,sys
from pathlib import Path
from flint import fmpq,fmpq_poly
resource.setrlimit(resource.RLIMIT_CPU,(60,60))
resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
root=Path(__file__).resolve().parents[2]; source=root/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'; roster=root/'artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1/norm4-sections.json'
s=json.loads(source.read_text()); rows=json.loads(roster.read_text())['records'];G=s['sections']['height_gram'];p=lambda cs:fmpq_poly([fmpq(str(c)) for c in cs]); A=p(s['weierstrass_model']['A_coefficients_low_to_high']);B=p(s['weierstrass_model']['B_coefficients_low_to_high'])
def coord(r,k):return tuple(p(r[k][d+'_coefficients_low_to_high']) for d in ['numerator','denominator'])
basis=[(coord(r,'X'),coord(r,'Y')) for r in s['sections']['records']]; TX,TD=basis[3][0]
def add(P,Q,a):
 if P is None:return Q
 if Q is None:return P
 x,y=P;u,v=Q
 if x==u:
  if y==-v:return None
  m=(3*x*x+a)/(2*y)
 else:m=(v-y)/(u-x)
 z=m*m-x-u;return z,m*(x-z)-y
sites=[]
for t in range(30):
 if all(d(t) for co in basis for n,d in co) and 4*A(t)**3+27*B(t)**2:sites.append(t)
 if len(sites)==7:break
assert len(sites)==7
values=[[(nx(t)/dx(t),ny(t)/dy(t)) for (nx,dx),(ny,dy) in basis] for t in sites]
def interp(vs):
 z=fmpq_poly([])
 for i,t in enumerate(sites):
  term=fmpq_poly([vs[i]])
  for j,u in enumerate(sites):
   if i!=j:term*=fmpq_poly([-u,1])/fmpq(t-u)
  z+=term
 return z
def witness(f):
 for prime in range(5,998):
  if any(prime % d == 0 for d in range(2,int(prime**.5)+1)):continue
  cs=[]
  for c in f.coeffs():
   if int(c.denominator)%prime==0:break
   cs.append(int(c.numerator)*pow(int(c.denominator),-1,prime)%prime)
  else:
   if not cs[-1]:continue
   if all(sum(c*pow(v,j,prime) for j,c in enumerate(cs))%prime for v in range(prime)):return prime
 raise AssertionError('No rational-root obstruction prime found')
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def phash(f):return hashlib.sha256(json.dumps([str(v) for v in f.coeffs()]).encode()).hexdigest()
out=[];start=time.monotonic()
for row in rows:
 w=row['word'];h=22-4*sum(G[3][j]*w[j] for j in range(17)); sign=1 if h==14 else -1 if h==30 else 0
 if not sign:continue
 w=[sign*v for v in w]; assert sum(w[i]*G[i][j]*w[j] for i in range(17) for j in range(17))==4; vv=[]
 for t,bb in zip(sites,values):
  R=None
  for c,P in zip(w,bb):
   if c<0:P=(P[0],-P[1])
   for _ in range(abs(c)):R=add(R,P,A(t))
  assert R is not None;vv.append(R)
 x=interp([v[0] for v in vv]);y=interp([v[1] for v in vv]);assert x.degree()<=4 and y.degree()<=6 and y*y==x*x*x+A*x+B
 n=(3*x*x+A)**2-8*x*y*y;d=4*y*y;contact=TX*d-TD*n
 while True:
  g=contact.gcd(TD*d)
  if g.degree()==0:break
  contact=contact//g
 unit,fac=contact.factor(); product=fmpq_poly([unit])
 for f,e in fac:product*=f**e
 assert product==contact and sorted((f.degree(),e) for f,e in fac)==[(5,1),(13,1)]
 assert TD.gcd(y).degree()==0, 'Possible contact over the trace pole'
 assert (TX.degree()-TD.degree(),x.degree(),y.degree(),n.degree(),d.degree())==(4,4,6,16,12)
 assert TX.leading_coefficient()/TD.leading_coefficient()!=n.leading_coefficient()/d.leading_coefficient(), 'Possible contact at infinity'
 out.append({'index':row['index'],'word':w,'height':14,'factor_checks':[{'degree':f.degree(),'sha256':phash(f),'no_rational_root_prime':witness(f)} for f,e in fac],'common_pole_contact':False,'infinity_contact':False,'x_sha256':phash(x),'y_sha256':phash(y)})
assert len(rows)==1313 and len(out)==230
result={'status':'PASS','source_sha256':digest(source),'roster_sha256':digest(roster),'checker_sha256':digest(Path(__file__)),'trace_index':3,'inherited_height':4,'distance_height':14,'sites':sites,'rows':out,'rational_contacts':0,'correlated_target_complete':False}
packet=root/'artifacts/generated-results/elkies-k3-q80-minimal-inherited-contact-v1'
if '--write' in sys.argv:
 packet.mkdir(exist_ok=True);(packet/'result.json').write_text(json.dumps(result,indent=2)+'\n')
else:assert result==json.loads((packet/'result.json').read_text())
print('PASS: all 230 height-four sections at distance fourteen; no rational contact, including poles and infinity')
