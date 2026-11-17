#!/usr/bin/env python3
"""Exact controls for the written intrinsic-period and unary-descent proofs.

Does NOT machine-verify the topology, units theorem, or curve normalization.
"""
from pathlib import Path
from fractions import Fraction as Q
from collections import deque
from math import comb
import hashlib,json
from verify_degree12_four_weights import groups

ROOT=Path(__file__).resolve().parent

def group_derived(name):
 field,gens,values=groups(name);add,mul,scale=field.add,field.mul,field.scale
 identity=(field.one,field.zero,field.zero,field.one)
 def mm(a,b):
  return (add(mul(a[0],b[0]),mul(a[1],b[2])),add(mul(a[0],b[1]),mul(a[1],b[3])),
          add(mul(a[2],b[0]),mul(a[3],b[2])),add(mul(a[2],b[1]),mul(a[3],b[3])))
 def inv(a):return (a[3],scale(a[1],-1),scale(a[2],-1),a[0])
 def comm(a,b):return mm(mm(mm(a,b),inv(a)),inv(b))
 generating=set(comm(a,b) for a in gens for b in gens)-{identity}
 while True:
  subgroup={identity};queue=deque([identity])
  while queue:
   a=queue.popleft()
   for g in generating:
    b=mm(a,g)
    if b not in subgroup:subgroup.add(b);queue.append(b)
  extra={mm(mm(g,a),inv(g)) for a in generating for g in gens}-subgroup
  if not extra:break
  generating|=extra
 assert subgroup<=values
 # Normality on subgroup generators suffices; commutators of the original
 # group generators in the subgroup make the quotient abelian.
 assert all(mm(mm(g,a),inv(g)) in subgroup for a in generating for g in gens)
 assert all(comm(a,b) in subgroup for a in gens for b in gens)
 expected=24 if name=='2O' else 120
 assert len(subgroup)==expected
 minus_identity=tuple(scale(x,-1) for x in identity)
 assert minus_identity in subgroup
 return {'order':len(values),'derived_order':len(subgroup),
         'abelianization_order':len(values)//len(subgroup),
         'center_minus_one_in_derived':True,'normal_and_abelian_quotient_checked':True}

def poly_mul(a,b):
 out={}
 for i,c in a.items():
  for j,d in b.items():out[i+j]=out.get(i+j,Q(0))+c*d
 return {i:c for i,c in out.items() if c}

def compose(a,b):
 out={};power={0:Q(1)}
 for i in range(max(a,default=0)+1):
  for j,c in power.items():out[j]=out.get(j,Q(0))+a.get(i,0)*c
  power=poly_mul(power,b)
 return {i:c for i,c in out.items() if c}

def main():
 # For t=tr(g), SU2 has h(t^(2k))=Catalan(k), odd moments zero.
 catalan=lambda k:Q(comb(2*k,k),k+1)
 haar=lambda p:sum(c*catalan(i//2) for i,c in p.items() if not i%2)
 t={1:Q(1)};square_trace={2:Q(1),0:Q(-2)};jacobian={2:Q(2)}
 assert haar(t)==0 and haar(square_trace)==-1
 weighted=[]
 for r in range(25):
  f={r:Q(1)}
  left=haar(poly_mul(jacobian,compose(f,square_trace)))
  right=2*haar(f)
  assert left==right
  weighted.append(str(left))
 # This is a ramified map, not an etale counterexample: its Jacobian is
 # 2 tr(g)^2 and vanishes on tr(g)=0. Unweighted transport really fails.
 assert compose(square_trace,t)==square_trace and haar(jacobian)==2
 # Affine centering in the unary theorem cannot be omitted.
 f={3:Q(1),1:Q(3),0:Q(5)}
 r={2:Q(1),1:Q(-10),0:Q(26)}  # (T-5)^2+1.
 composed=compose(r,f)
 assert all(i%2==0 for i in composed)
 assert any(i%2 for i in f) and f[0]==5
 centered={i:c for i,c in f.items() if i}
 assert all(i%2 for i in centered)
 records={name:group_derived(name) for name in ['2O','2I']}
 assert records['2O']['abelianization_order']==2
 assert records['2I']['abelianization_order']==1
 receipt={
  'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  'imported_arithmetic_sha256':{n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in [
   'verify_degree12_four_weights.py','verify_haar_shear_witnesses.py','verify_octahedral_minimal_type.py']},
  'finite_groups':records,
  'square_map_negative_control':{'h_trace':0,'h_trace_after_square':-1,'jacobian':'2*tr(g)^2','mean_jacobian':2},
  'weighted_square_map_controls_orders0_to24':weighted,
  'unary_affine_center_control':{'f':'t^3+3t+5','R':'(T-5)^2+1','invariant_composition':True,'center':5},
  'completed':True,'universal_period_proof_machine_verified':False,
  'universal_curve_normalization_proof_machine_verified':False,
  'exceptional_MZ_decided':False,'independent_replay':False,'external_review':False,'formal_verification':False}
 (ROOT/'intrinsic_haar_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))

if __name__=='__main__':main()
