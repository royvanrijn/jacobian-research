import json,sympy as S,time,resource
from pathlib import Path
resource.setrlimit(resource.RLIMIT_CPU,(30,30));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
r=Path('research');src=r/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json';s=json.load(open(src));t=S.symbols('t');out=Path('/tmp/q80-signed-contact-bank');out.mkdir(exist_ok=True)
primes=[181,191,193,197,199,211,223,227,229]
def pol(cs):return S.Poly(sum(S.Rational(c)*t**j for j,c in enumerate(cs)),t)
def coord(i,key):
 v=s['sections']['records'][i][key];return pol(v['numerator_coefficients_low_to_high']),pol(v['denominator_coefficients_low_to_high'])
def power(a,n,f):
 o=S.Poly(1,t,modulus=f.get_modulus())
 while n:
  if n&1:o=(o*a)%f
  a=(a*a)%f;n//=2
 return o
A=pol(s['weierstrass_model']['A_coefficients_low_to_high']);H=s['sections']['height_gram']
for i in [1,2]:
 X,D=coord(i,'X')
 for j in range(17):
  hs={str(sign):H[i][i]+4*H[j][j]-4*sign*H[i][j] for sign in [1,-1]}
  if max(hs.values())<20:continue
  x,dx=coord(j,'X');y,dy=coord(j,'Y')
  # x(2R)=((3x^2+A)^2/(4y^2))-2x.
  NN=(3*x*x+A*dx*dx)**2*dy*dy-8*x*dx**3*y*y
  DD=4*y*y*dx**4
  common=S.gcd(NN,DD);NN=NN.exquo(common);DD=DD.exquo(common)
  N=X*DD-D*NN
  # Remove coordinate-pole factors. Common poles are audited separately.
  exceptional=D*DD
  while True:
   g=S.gcd(N,exceptional)
   if g.degree()==0:break
   N=N.exquo(g)
  common_poles=S.gcd(D,DD).monic()
  # At a common pole T=2R=O, the half is either O or nonzero2-torsion.
  # Existing basis pole/torsion gate excludes nonzero2-torsion; O is the ruled-section degeneracy.
  result={'trace':i,'basis_section':j,'signed_heights':hs,'polynomial':[str(c) for c in reversed(N.all_coeffs())],'degree':N.degree(),'common_poles':[str(c) for c in reversed(common_poles.all_coeffs())],'status':'UNRESOLVED'}
  for p in primes:
   try:cs=[int(c.p%p)*pow(int(c.q),-1,p)%p for c in N.all_coeffs()]
   except ValueError:continue
   if not cs[0]:continue
   f=S.Poly.from_list(cs,t,modulus=p);z=S.Poly(t,t,modulus=p)
   dg=S.gcd(f,power(z,p*p,f)-z).degree()
   if dg==0:result.update(status='NO_DEGREE_LE_TWO_NONPOLE',prime=p);break
  (out/f'{i}-{j}.json').write_text(json.dumps(result,indent=2)+'\n')
  print(i,j,hs,N.degree(),result['status'],result.get('prime'),flush=True)
