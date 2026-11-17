#!/usr/bin/env python3
"""Exact application/normal-form controls for a written small-component theorem.
Does not independently reprove the universal linear algebra or imported Moh theorem.
"""
from pathlib import Path
from itertools import product
import hashlib,json,sys,time
import sympy as S
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 started=time.monotonic();ledger=json.loads((ROOT/'MATH_STATUS.json').read_text());claims={e['id']:e for e in ledger['entries']};kp=claims['KPAIR1'];base=ROOT/'jacobian/sparse-eight';source=base/'verify_coordinate_gates.py';old=json.loads((base/'coordinate-gates-validation.json').read_text());assert kp['state']=='proved' and kp['artifact_hash']=='sha256:'+sha(source) and old['status']=='pass' and old['source_sha256']==sha(source) and old['plane_degree_import']==100
 assert sha(source)=='252f9f217444f8d300f9dd559c714043dd90bf56172e63500b28b08cdd419a7e';assert sha(base/'coordinate-gates-validation.json')=='9d625c89cbc41fbc3b9e0001b418bb41e77d442f4d0b19eef5678b842303854a'
 primary=HERE/'moh-primary-account.html';assert sha(primary)=='29591bbc082ec140b2dd85b9a34c891af021f336a9ac4cb7d84fc81690c01ed3';assert b'degrees of the two polynomials' in primary.read_bytes() and b'one hundred' in primary.read_bytes()
 x,y,z,t,u,v=S.symbols('x y z t u v');xyz=(x,y,z);controls=[]
 for rank in range(3):
  q=S.Add(*(a*a for a in xyz[:rank]));R=z+q;E=(u,v,t-q.subs({x:u,y:v},simultaneous=True));assert S.expand(R.subs(dict(zip(xyz,E)),simultaneous=True)-t)==0;assert S.Matrix(E).jacobian((t,u,v)).det()==1;assert S.diff(R,z)==1;assert S.hessian(R,xyz).rank()==rank;controls.append(dict(hessian_rank=rank,normal_form=str(R),inverse_straightening=list(map(str,E)),directional_derivative='1',jacobian_determinant='1'))
 ell=S.symbols('ell1:4');A=S.diag(2,2,2);R=sum(a*a for a in xyz)+sum(a*b for a,b in zip(ell,xyz));critical=[-a/2 for a in ell];assert all(S.expand(S.diff(R,a).subs(dict(zip(xyz,critical)),simultaneous=True))==0 for a in xyz)
 # All exponent triples through each degree cap: substitution exponent weight
 # a+b+s*c is an upper bound for every plane monomial after z=t-H(u,v).
 degree_bounds=[]
 for s,D in [(1,100),(2,50)]:
  count=0;maximum=0
  for a in range(D+1):
   for b in range(D+1-a):
    for c in range(D+1-a-b):
     weight=a+b+s*c;assert weight<=s*D<=100;maximum=max(maximum,weight);count+=1
  assert maximum==100;degree_bounds.append(dict(component_degree=s,maximum_map_degree=D,monomials_checked=count,maximum_plane_degree=maximum))
 profiles=sorted({(a,b,c) for a,b,c in product(range(1,8),repeat=3) if a>=b>=c and (a,b,c)<(7,6,4)});sys.path.insert(0,str(ROOT/'foundations'));from jcsearch.low_degree_census import degree_profiles_below
 assert tuple(profiles)==degree_profiles_below((7,6,4));closed=[p for p in profiles if p[2]<=2];remaining=[p for p in profiles if p[2]>=3];assert len(profiles)==74 and len(closed)==47 and len(remaining)==27 and set(closed).isdisjoint(remaining) and sorted(closed+remaining)==profiles;assert sum(p[2]==1 for p in closed)==27 and sum(p[2]==2 for p in closed)==20;prior=[(1,1,1),(2,1,1),(2,2,1),(2,2,2)];assert set(prior)<=set(closed) and len(set(closed)-set(prior))==43;assert (7,6,4) not in profiles
 cubic=[p for p in remaining if p[0]<=3];quartic=[p for p in remaining if p[0]<=4];assert cubic==[(3,3,3)] and quartic==[(3,3,3),(4,3,3),(4,4,3),(4,4,4)]
 result=dict(status='pass',source_sha256=sha(Path(__file__)),proof_sha256=sha(HERE/'PROOF.md'),prior_checker_sha256=sha(source),prior_validation_sha256=sha(base/'coordinate-gates-validation.json'),prior_gate_claim='KPAIR1',universal_linear_algebra_proof='Written annihilator proof, not formalized or inferred from normal-form controls.',symbolic_controls=controls,rank3_unique_critical_point=list(map(str,critical)),all_hessian_rank_strata_covered_in_written_proof=True,degree_bounds=degree_bounds,plane_import=dict(author='T.T.Moh',bound=100,url=old['primary_plane_source'],primary_snapshot_sha256=sha(primary),independently_reproved=False),profile_compiler_sha256=sha(ROOT/'foundations/jcsearch/low_degree_census.py'),profile_enumeration_independent=True,total_profiles=74,excluded_profiles=closed,excluded_count=47,previous_dense_quadratic_profiles=prior,additional_dense_profiles=43,remaining_profiles=remaining,remaining_count=27,remaining_smallest_entry_at_least=3,conditional_cubic_import_adds=cubic,conditional_cubic_remaining=26,conditional_quartic_import_adds=quartic,conditional_quartic_remaining=23,new_external_theorem_imported=False,global_degree_or_sparsity_minimum_proved=False,positive_construction_bounds_changed=False,normalized_support_interval=[10,48],scope='All characteristic-zero3D Keller collisions with a degree1 output combination and D<=100, or degree2 output combination and D<=50, excluded by written linear algebra and plane reduction importing Moh100. All47 invariant profiles below764 with smallest entry<=2 are excluded without support bound. Checker verifies exact application, degree bookkeeping and controls; it is not an independent proof of the general theorem or Moh.',seconds=round(time.monotonic()-started,3));(HERE/'validation.json').write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print('PASS47 dense profiles excluded;27 remain; degree caps100/50; controls only, Moh imported; seconds',result['seconds'])
if __name__=='__main__':main()
