import json,sympy as S,time
from pathlib import Path
r=Path('research');s=json.load(open(r/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'));t=S.symbols('t')
def pol(cs):return S.Poly(sum(S.Rational(c)*t**j for j,c in enumerate(cs)),t)
def coord(i,key):
 a=s['sections']['records'][i][key];return pol(a['numerator_coefficients_low_to_high']),pol(a['denominator_coefficients_low_to_high'])
A=pol(s['weierstrass_model']['A_coefficients_low_to_high']);x,dx=coord(16,'X');y,dy=coord(16,'Y');assert dx.degree()==dy.degree()==0
x=x.mul_ground(1/dx.LC());y=y.mul_ground(1/dy.LC());den=4*y*y;num=(3*x*x+A)**2-8*x*y*y
for i in [1,2]:
 start=time.monotonic();X,D=coord(i,'X');N=X*den-D*num
 print('trace',i,'degree',N.degree(),'seconds',time.monotonic()-start,flush=True)
 assert S.gcd(D,y).degree()==0
 for p in list(S.primerange(5,200)):
  try:
   cs=[int(c.p%p)*pow(int(c.q),-1,p)%p for c in N.all_coeffs()]
  except ValueError:continue
  if cs[0]==0:continue
  F=S.Poly.from_list(cs,t,modulus=p);z=S.Poly(t,t,modulus=p)
  xp=S.Poly(1,t,modulus=p)
  def power(a,n):
   o=S.Poly(1,t,modulus=p)
   while n:
    if n&1:o=(o*a)%F
    a=(a*a)%F;n//=2
   return o
  if S.gcd(F,power(z,p*p)-z).degree()==0:
   print('no degree1/2 factors at',p,flush=True);break
 else:print('NO WITNESS')
 infinity_X=X.LC()/D.LC();infinity_double=num.LC()/den.LC()
 print('infinity distinct',infinity_X!=infinity_double,flush=True)
 Path(f'/tmp/contact-trace-{i}.json').write_text(json.dumps({'trace':i,'N':[str(c) for c in reversed(N.all_coeffs())],'witness_prime':int(p),'infinity_distinct':bool(infinity_X!=infinity_double)},indent=2)+'\n')
