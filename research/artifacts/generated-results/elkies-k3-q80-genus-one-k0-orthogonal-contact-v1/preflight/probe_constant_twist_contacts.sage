from sage.all import GF, PolynomialRing, QQ
from pathlib import Path
import json,resource,time
resource.setrlimit(resource.RLIMIT_CPU,(20,25));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
ROOT=Path(__file__).resolve().parents[4];source=json.loads((ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json').read_text())['weierstrass_model']
F=GF(131);P=PolynomialRing(F,'t');t=P.gen();A=P([F(QQ(x)) for x in source['A_coefficients_low_to_high']]);B=P([F(QQ(x)) for x in source['B_coefficients_low_to_high']]);q=t*t+62*t+88;H=q*q;S=73+6*t
started=time.process_time();x0=(S-(S**3+A*S+B)*(3*S*S+A).inverse_mod(H))%H;rows=[]
for c in F:
 x=x0+c*H;f,r=(x**3+A*x+B).quo_rem(H);assert not r
 lead=f.leading_coefficient();monic=f/lead
 if not monic.is_square():continue
 y=monic.sqrt();rows.append({'c':int(c),'x':[int(z) for z in x.list()],'r':[int(z) for z in y.list()],'scalar':int(lead),'scalar_square':bool(lead.is_square())})
print(json.dumps({'x0':[int(z) for z in x0.list()],'tried':131,'rows':rows,'cpu_seconds':time.process_time()-started},indent=2))
