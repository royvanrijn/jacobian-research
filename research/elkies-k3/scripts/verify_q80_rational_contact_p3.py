import json,sympy as S,resource,time
from pathlib import Path
resource.setrlimit(resource.RLIMIT_CPU,(25,25));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
root=Path(__file__).resolve().parents[2]
source=root/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
s=json.loads(source.read_text());t=S.symbols('t');H=s['sections']['height_gram']
def pol(cs):return S.Poly(sum(S.Rational(c)*t**i for i,c in enumerate(cs)),t,domain=S.QQ)
def coord(i,k):return [pol(s['sections']['records'][i][k][z+'_coefficients_low_to_high']) for z in ['numerator','denominator']]
A=pol(s['weierstrass_model']['A_coefficients_low_to_high']);B=pol(s['weierstrass_model']['B_coefficients_low_to_high']);i=3;X,D=coord(i,'X');Y,DY=coord(i,'Y');rows=[]
for j in range(17):
 x,dx=coord(j,'X');y,dy=coord(j,'Y');NN=(3*x*x+A*dx*dx)**2*dy*dy-8*x*dx**3*y*y;DD=4*y*y*dx**4;g=S.gcd(NN,DD);NN=NN.exquo(g);DD=DD.exquo(g);cp=S.gcd(D,DD).monic();N=X*DD-D*NN
 while True:
  g=S.gcd(N,D*DD)
  if g.degree()==0:break
  N=N.exquo(g)
 fs=S.factor_list(N)[1];witnesses=[]
 for f,e in fs:
  assert f.degree()>1
  witness=None
  for p in S.primerange(5,998):
   try:cs=[int(c.p)*pow(int(c.q),-1,p)%p for c in f.all_coeffs()]
   except ValueError:continue
   if not cs[0]:continue
   if all(sum(c*pow(v,len(cs)-1-k,p) for k,c in enumerate(cs))%p for v in range(p)):
    witness=int(p);break
  assert witness is not None
  witnesses.append({'degree':f.degree(),'multiplicity':e,'no_rational_root_prime':witness})
 assert S.prod(f**e for f,e in fs).monic()==N.monic()
 roots=[-f.nth(0)/f.nth(1) for f,e in fs if f.degree()==1]
 inf_equal=(NN.degree()-DD.degree()==4 and X.degree()-D.degree()==4 and NN.LC()/DD.LC()==X.LC()/D.LC())
 assert X.degree()-D.degree()==4 and NN.degree()-DD.degree()==4 and not inf_equal
 assert cp.degree()==0 or (j==i and cp==D.monic() and S.sqf_part(cp).degree()==1)
 row={'factor_witnesses':witnesses,'infinity_x_equality':inf_equal,'basis':j,'degree':N.degree(),'heights':[int(H[i][i]+4*H[j][j]-4*sign*H[i][j]) for sign in [1,-1]],'factors':[[f.degree(),e] for f,e in fs],'rational_roots':[str(v) for v in roots],'common_pole_factors':[[f.degree(),e] for f,e in S.factor_list(cp)[1]]};rows.append(row)
import hashlib,sys
result={'status':'PASS','trace_index':3,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'cases':rows,'nonpole_rational_contacts':0,'infinity_contacts':0,'only_common_pole':'P3=O at its unique rational pole','positive_correlated_target_complete':False}
packet=root/'artifacts/generated-results/elkies-k3-q80-rational-contact-p3-v1'
if '--write' in sys.argv:
 packet.mkdir(exist_ok=True);(packet/'result.json').write_text(json.dumps(result,indent=2)+'\n')
else:assert result==json.loads((packet/'result.json').read_text())
print('PASS: first height-six trace, signed basis contact test; no eligible rational contacts')
