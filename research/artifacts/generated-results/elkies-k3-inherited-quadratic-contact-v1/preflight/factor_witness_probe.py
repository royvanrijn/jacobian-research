import json,sympy as S,resource
from pathlib import Path
resource.setrlimit(resource.RLIMIT_CPU,(30,30));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
t=S.symbols('t');out=Path('/tmp/q80-signed-contact-factors');out.mkdir(exist_ok=True)
def power(a,n,f):
 o=S.Poly(1,t,modulus=f.get_modulus())
 while n:
  if n&1:o=(o*a)%f
  a=(a*a)%f;n//=2
 return o
for path in sorted(Path('/tmp/q80-signed-contact-bank').glob('*.json')):
 r=json.loads(path.read_text())
 if r['status']!='UNRESOLVED':continue
 N=S.Poly(sum(S.Rational(c)*t**i for i,c in enumerate(r['polynomial'])),t);factors=[]
 for f,e in S.factor_list(N)[1]:
  f=f.monic();rec={'degree':f.degree(),'multiplicity':e,'coefficients':[str(c) for c in reversed(f.all_coeffs())]}
  if f.degree()<=2:rec['candidate']=True
  else:
   for p in S.primerange(5,998):
    try:cs=[int(c.p%p)*pow(int(c.q),-1,p)%p for c in f.all_coeffs()]
    except ValueError:continue
    ff=S.Poly.from_list(cs,t,modulus=p);z=S.Poly(t,t,modulus=p)
    if S.gcd(ff,power(z,p*p,ff)-z).degree()==0:rec['prime']=int(p);break
  factors.append(rec)
 row={'trace':r['trace'],'basis_section':r['basis_section'],'factors':factors}
 (out/path.name).write_text(json.dumps(row,indent=2)+'\n')
 print(path.name,[(f['degree'],f.get('prime'),f.get('candidate')) for f in factors],flush=True)
