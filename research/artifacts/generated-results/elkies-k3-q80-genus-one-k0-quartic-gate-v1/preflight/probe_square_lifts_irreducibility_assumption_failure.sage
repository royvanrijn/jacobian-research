from sage.all import GF, Integers, PolynomialRing, QQ, matrix, vector
from pathlib import Path
import json, resource, time
resource.setrlimit(resource.RLIMIT_CPU,(30,35));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
R=Path(__file__).resolve().parents[4];out=Path(__file__).resolve().parent
source=json.loads((R/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json').read_text())['weierstrass_model']
points=json.loads((R/'artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1/norm4-sections.json').read_text())['records']
pairs=json.loads((out/'square-branch-pairs-preview.json').read_text())['distinct_root_contact_pairs']
p=131;F=GF(p);P=PolynomialRing(F,'t');t=P.gen();Z=Integers(p*p);P2=PolynomialRing(Z,'t')
A=P([F(QQ(v)) for v in source['A_coefficients_low_to_high']]);B=P([F(QQ(v)) for v in source['B_coefficients_low_to_high']])
A2=P2([Z(QQ(v).numerator())/Z(QQ(v).denominator()) for v in source['A_coefficients_low_to_high']])
B2=P2([Z(QQ(v).numerator())/Z(QQ(v).denominator()) for v in source['B_coefficients_low_to_high']])
lift=lambda f:P2([int(v) for v in f.list()])
start=time.process_time();rows=[]
for row in pairs:
 for hh in row['H']:
  H=P(hh);D=H*H
  assert H.degree()==2 and H.is_irreducible()
  J=matrix(F,26,24);error=[]
  for i,n in enumerate(row['pair']):
   X,Y=P(points[n]['x']),P(points[n]['y']);r,rem=Y.quo_rem(H);assert not rem
   xp=[(3*X*X+A)*t**k for k in range(5)];rp=[-2*D*r*t**k for k in range(5)];dp=[-r*r*t**k for k in range(4)]
   for j,f in enumerate(xp+rp):
    for k in range(13):J[13*i+k,10*i+j]=f[k]
   for j,f in enumerate(dp):
    for k in range(13):J[13*i+k,20+j]=f[k]
   err=lift(X)**3+A2*lift(X)+B2-lift(D)*lift(r)**2
   assert all(int(err[k])%p==0 for k in range(13))
   error.extend([int(F(-int(err[k])//p)) for k in range(13)])
  rank=int(J.rank());aug=int(J.augment(vector(F,error).column()).rank())
  rows.append({'pair':row['pair'],'H':hh,'rank':rank,'augmented_rank':aug,'error_target':error})
print(json.dumps({'rows':rows,'cpu_seconds':time.process_time()-start},indent=2))
