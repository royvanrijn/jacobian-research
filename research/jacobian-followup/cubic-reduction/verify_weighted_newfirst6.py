#!/usr/bin/env python3
"""Independent audit of the marked weighted degree-six exclusion via theta^7."""
from pathlib import Path
from fractions import Fraction as Q
import hashlib,json,time
from verify import add,mul,variable,derivative,determinant3
from verify_power_single_compress import decode
from weighted_verify_utils import substitution,expression,equal_generators,rank,coefficient_buckets
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 started=time.monotonic();source=HERE/'weighted-newfirst6/input.json';d=json.loads(source.read_text());free=d['variables'];m=len(free);assert m==13 and d['degree']==6 and d['source_sha256']==sha(HERE/'weighted_newfirst.py')
 atoms=d['atoms'];assert len(atoms)==17;expected=set()
 for row,target in enumerate([1,-1,-2]):
  for a in range(7):
   for b in range(7):
    for c in range(7):
     if 2<=a+b+c<=6 and a-b-2*c==target:expected.add((row,(a,b,c)))
 assert {(a['row'],tuple(a['exponent'])) for a in atoms}==expected
 names=list(d['substitution']);n=len(names);assert n==15;index={name:i for i,name in enumerate(names)};N=3+n;F=[variable(i,N) for i in range(3)]
 for atom in atoms:
  e=list(atom['exponent'])+[0]*n
  if atom['coefficient'] in index:e[3+index[atom['coefficient']]]=1;c=Q(1)
  else:c=Q(atom['coefficient']);assert atom['row']==0 and c in [Q(-3),Q(-1)]
  F[atom['row']]=add(F[atom['row']],{tuple(e):c})
 det=add(determinant3([[derivative(f,j) for j in range(3)] for f in F],N),{(0,)*N:Q(1)},Q(-1));buckets=coefficient_buckets(det,3)
 assert all(e[0]==e[1]+2*e[2] for e in buckets)
 def point(f):
  out={}
  for e,c in f.items():out=add(out,{e[3:]:c*(Q(-1)**e[0])*(Q(4)**e[2])})
  return out
 original=list(buckets.values())+[point(F[0]),add(point(F[1]),{(0,)*n:Q(1)},Q(-1))]
 linear=[p for p in original if p and max(sum(e) for e in p)<=1];assert all(sum(e)==1 for p in linear for e in p)
 values=[expression(d['substitution'][name],free) for name in names]
 for i,name in enumerate(names):
  if name in free:assert values[i]==variable(free.index(name),m)
 assert all(not substitution(p,values,m) for p in linear)
 assert rank([[p.get(tuple(int(k==j) for k in range(n)),Q(0)) for j in range(n)] for p in linear],n)==n-m==2
 native=[substitution(p,values,m) for p in original];native=[p for p in native if p];prepared=[decode(p,m) for p in d['equations']];equal_generators(native,prepared)
 out=HERE/'weighted-newfirst6/parameter-field';certpath=out/'certificate.json';cert=json.loads(certpath.read_text());assert cert['status']=='exact-primary-parameter-certificate' and cert['source_sha256']==sha(HERE/'weighted_newfirst6_lift.py') and cert['input_sha256']==sha(source)
 assert cert['singular_source_sha256']==sha(out/'lift.sing');P=cert['variables'];assert len(P)==m and P[0]=='theta' and set(P[1:])==set(free)-{'f2'}
 rename=[variable(P.index('theta' if name=='f2' else name),m) for name in free];derived=[substitution(p,rename,m) for p in prepared];generators=[decode(p,m) for p in cert['generators']];equal_generators(derived,generators)
 target=decode(cert['target'],m);e=[0]*m;e[0]=7;assert target=={tuple(e):Q(1)}
 multipliers=[decode(p,m) for p in cert['multipliers']];assert len(multipliers)==len(generators);result={}
 for a,b in zip(generators,multipliers):result=add(result,mul(a,b))
 assert result==target
 # The exceptional theta0 slice embeds in the independently classified degree7 tails.
 tailpath=HERE/'weighted-tail7/validation.json';tail=json.loads(tailpath.read_text());assert tail['status']=='pass' and tail['source_sha256']==sha(HERE/'verify_weighted_tail7.py')
 wide=json.loads((HERE/'weighted-tail7/input.json').read_text());wide_atoms={(a['row'],tuple(a['exponent'])) for a in wide['tail_atoms']}
 for atom in atoms:
  if atom['row']>0:assert (atom['row'],tuple(atom['exponent'])) in wide_atoms
 theta_index=free.index('f2');zero_values=[variable(i,m) if i!=theta_index else {} for i in range(m)]
 first=[a for a in atoms if a['row']==0 and a['coefficient'] in index];assert len(first)==2
 for atom in first:assert not substitution(values[index[atom['coefficient']]],zero_values,m)
 raw=json.loads((HERE/'raw13-specialization/input.json').read_text());forbidden=next(a for a in raw['atoms'] if a['row']==2 and a['exponent']==[3,3,1]);assert Q(forbidden['coefficient'])==8
 assert (2,(3,3,1)) not in expected and tail['unique_nonlinear_total']==13
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/'verify.py'),utility_sha256=sha(HERE/'weighted_verify_utils.py'),input_sha256=sha(source),parameter_certificate_sha256=sha(certpath),tail_classification_sha256=sha(tailpath),
  nonlinear_atoms=17,coefficient_parameters=15,eliminated_linear_parameters=2,certificate_target='theta^7',unit_identities=0,cleared_parameter_identities=1,
  multiplier_monomials=sum(len(p) for p in multipliers),max_multiplier_degree=max(sum(e) for p in multipliers for e in p),
  universal_three_variable_determinant_reconstructed=True,linear_parameterization_complete=True,cleared_identity_multiplied=True,exceptional_parameter_roots_exactly_zero=True,zero_slice_excluded_by_tail_classification=True,
  seconds=round(time.monotonic()-started,3),scope='No degree<=6 weight-compatible map with source weights(1,-1,-2), identity linear part, first-output a^2b/a^3c coefficients-3,-1, and marked collision between(0,1,c0) and(-1,1,4) exists over a characteristic-zero field. All17 possible nonlinear atoms are allowed. The certificate forces theta0, whose tails would have to contain a nonzero degree7 term. This is not a global degree-six or raw12-term exclusion.')
 (HERE/'weighted-newfirst6/validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
