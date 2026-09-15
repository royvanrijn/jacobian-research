from sage.all import *
from pathlib import Path
import json, argparse
P=Path(__file__).resolve().parents[2]/'artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1'
x=json.loads((P/'rank1-frames.json').read_text());hit=next(h for h in x['hits'] if h['P_dot_O']==0)
F=matrix(ZZ,hit['frame_gram']);R=F[:16,:16];v=F[:16,16];h=F[16,16]-(v.transpose()*R.inverse()*v)[0,0]
assert R==block_diagonal_matrix([CartanMatrix(['A',n]) for n in [2,3,3,8]])
assert h==QQ(85)/36 and F[16,16]==4 and F.det()==1020 and F.is_positive_definite()
assert F.elementary_divisors()==[1]*16+[1020]
D,U,V=F.smith_form();g=V.inverse().row(16)*F.inverse();g=vector(QQ,[z-floor(z) for z in g]);q=g*F*g
assert lcm(z.denominator() for z in g)==1020
units=[u for u in range(1020) if gcd(u,1020)==1 and (q+u*u*QQ(403)/1020)/2 in ZZ];assert units
# Norm decomposition n^2*h + (r+n*R^-1*v)^2 proves all roots are in R since h>2.
z=F.__pari__().qfminim(2,10000,2);assert z[0]==2*(3+6+6+36)
result={'status':'PASS_RANK16_ROOT_FRAME_MW1','root_type':'A2+2A3+A8','frame_gram':[[int(z) for z in r] for r in F.rows()], 'height':str(h),'root_lattice_determinant':int(R.det()),'root_vectors':int(z[0]),'determinant':1020,'generator_frame_norm':4,'P_dot_O':0,'component_indices':[0,0,1,1], 'discriminant_generator':[str(z) for z in g],'q_generator':str(q),'isometry_multiplier_to_T':units[0], 'scope':'Exact abstract NS frame. Indefinite uniqueness and the already admitted full rational NS imply arithmetic MW1 fibration existence. No physical nef basis or equation is supplied.'}
a=argparse.ArgumentParser();a.add_argument('--check',action='store_true');args=a.parse_args()
if args.check: assert json.loads((P/'rank1-certificate.json').read_text())==result
else: (P/'rank1-certificate.json').write_text(json.dumps(result,indent=2)+'\n')
print(result['status'])
