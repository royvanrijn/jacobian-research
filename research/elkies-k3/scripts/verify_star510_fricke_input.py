#!/usr/bin/env python3
"""Verify a retained modular point in the correct cached differential basis."""
import argparse,hashlib,json,resource
from itertools import product
from pathlib import Path
import sympy as S
from flint import fmpz_poly
ROOT=Path(__file__).resolve().parents[2]
PACKET=ROOT/'artifacts/generated-results/elkies-k3-star510-fricke-input-v1'
def run():
 resource.setrlimit(resource.RLIMIT_CPU,(20,20));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
 d=json.loads((PACKET/'input.json').read_text());assert d['level']==510
 index=510*S.Rational(3,2)*S.Rational(4,3)*S.Rational(6,5)*S.Rational(18,17);bound=int(8*index/12)
 assert index==1296 and bound==d['sturm_bound']==864
 fs=[fmpz_poly(f) for f in d['forms']];assert len(fs)==3 and all(len(f)==865 for f in d['forms'])
 mons=d['monomials'];cs=d['coefficients'];ps=[fs[0]**a*fs[1]**b*fs[2]**c for a,b,c in mons]
 assert S.Matrix([[int(f[j]) for f in ps] for j in range(80)]).rank()==14
 assert all(sum(k*int(f[j]) for k,f in zip(cs,ps))==0 for j in range(bound+1))
 x,y,z=S.symbols('x y z');F=sum(k*x**a*y**b*z**c for k,(a,b,c) in zip(cs,mons));P=dict(zip([x,y,z],d['point']))
 assert F.subs(P)==0
 grad=[S.diff(F,v) for v in [x,y,z]];gradient=[int(g.subs(P)) for g in grad];assert gradient==[72,-376,552]
 assert S.expand(F.subs(y,3*z)-2*z*(x-8*z)*(x*x-8*x*z+36*z*z))==0
 # Smoothness over one finite field is enough for the literal plane quartic.
 witness=None
 for p in [3,5,7,11,13,19,23,29,31]:
  affine=S.groebner([g.subs(z,1) for g in grad],x,y,modulus=p)
  infinity=S.groebner([g.subs({y:1,z:0}) for g in grad],x,modulus=p)
  endpoint=any(int(g.subs({x:1,y:0,z:0}))%p for g in grad)
  if list(affine)==[1] and list(infinity)==[1] and endpoint:witness=p;break
 assert witness is not None
 # Distinguish isotropic Inose T from the existing determinant1020 frame.
 frame_path=ROOT/'artifacts/generated-results/elkies-k3-det1020-root-gate-v1/certificate.json'
 frame=json.loads(frame_path.read_text());assert S.Rational(frame['q_frame'])==S.Rational(10877,1020)
 assert 10877%5==2 and {u*u%5 for u in range(1,5)}=={1,4}
 return {'status':'PASS','input_sha256':hashlib.sha256((PACKET/'input.json').read_bytes()).hexdigest(),'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'upstream_commit':d['upstream_commit'],'curve':str(F),'point':d['point'],'gradient':gradient,'smoothness_prime':witness,'geometric_smoothness_via_unit_jacobian_ideals':True,'weight8_sturm_bound':bound,'q_expansion_identity':True,'modularity_and_nonCM_status':'inherited from authors, not independently rederived','fricke_to_star_degree':8,'rational_fricke_lift':'UNKNOWN','isotropic_Inose_MW17_parent':'UNKNOWN','existing_det1020_frame_source_sha256':hashlib.sha256(frame_path.read_bytes()).hexdigest(),'existing_det1020_frame_wrong_discriminant_form_modulus':5,'positive_correlated_target_complete':False}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--write',action='store_true');args=p.parse_args();out=run()
 if args.write:(PACKET/'result.json').write_text(json.dumps(out,indent=2)+'\n')
 else:assert out==json.loads((PACKET/'result.json').read_text())
 print(json.dumps(out,sort_keys=True))
