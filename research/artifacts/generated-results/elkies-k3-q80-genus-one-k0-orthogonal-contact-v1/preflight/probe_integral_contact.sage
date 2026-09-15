from sage.all import Qp, PolynomialRing, QQ, ZZ
from pathlib import Path
import json,resource,time
resource.setrlimit(resource.RLIMIT_CPU,(30,35));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
ROOT=Path(__file__).resolve().parents[4]
source=json.loads((ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json').read_text())
point=json.loads((ROOT/'artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1/norm4-sections.json').read_text())['records'][732]
started=time.process_time();base=Qp(131,prec=20);R=PolynomialRing(base,'z');z=R.gen();K=base.extension(z*z+62*z+88,names='theta');theta=K.gen();P=PolynomialRing(K,'t');t=P.gen()
poly=lambda values:P([K(QQ(v)) for v in values])
model=source['weierstrass_model'];A=poly(model['A_coefficients_low_to_high']);B=poly(model['B_coefficients_low_to_high']);delta=4*A**3+27*B**2
alpha=theta
for i in range(6):alpha-=delta(alpha)/delta.derivative()(alpha)
e=-3*B(alpha)/(2*A(alpha));s=(3*e).sqrt();torus=K(1);rows=[]
for i,w in enumerate(point['word']):
 if not w:continue
 row=source['sections']['records'][i]
 X=poly(row['X']['numerator_coefficients_low_to_high'])(alpha)/poly(row['X']['denominator_coefficients_low_to_high'])(alpha)
 Y=poly(row['Y']['numerator_coefficients_low_to_high'])(alpha)/poly(row['Y']['denominator_coefficients_low_to_high'])(alpha)
 num=Y+s*(X-e);den=Y-s*(X-e);value=num/den
 torus*=value**w
 rows.append({'basis_index':i,'word':w,'x_valuation':str(X.valuation()),'y_valuation':str(Y.valuation()),'torus_value':str(value)})
u=s*(torus+1)/(torus-1);xt=u*u-2*e;yt=u*(u*u-3*e)
print(json.dumps({'status':'BOUNDED_LOCAL_PREVIEW','precision':20,'alpha':str(alpha),'nodal_root':str(e),'sqrt_3e':str(s),
 'basis_rows':rows,'T_torus':str(torus),'T_u':str(u),'T_y':str(yt),'T_y_valuation':str(yt.valuation()),
 'T_x_minus_simple_valuation':str((xt+2*e).valuation()),'cpu_seconds':time.process_time()-started},indent=2))
