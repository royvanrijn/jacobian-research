#!/usr/bin/env python3
"""Independent raw13 coefficient-pattern classification and certificate replay."""
from pathlib import Path
from fractions import Fraction as Q
from collections import Counter
import hashlib,json,time
from verify import add,mul,variable,derivative,determinant3,evaluate
from verify_power_single_compress import decode,scale
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def power(p,k,n):
 out={(0,)*n:Q(1)}
 for i in range(k):out=mul(out,p)
 return out
def canonical(p):return tuple(sorted(p.items()))
def substitute(p,values,n):
 result={};cache={}
 for e,c in p.items():
  term={(0,)*n:c}
  for i,k in enumerate(e):
   if not k:continue
   if (i,k) not in cache:cache[i,k]=power(values[i],k,n)
   term=mul(term,cache[i,k],)
  result=add(result,term)
 return result

def main():
 started=time.monotonic();out=HERE/'raw13-rigidity';rawpath=HERE/'raw13-specialization/input.json';normalpath=HERE/'raw13-normalized/input.json';raw=json.loads(rawpath.read_text());normal=json.loads(normalpath.read_text());atoms=raw['atoms']
 assert raw['source_sha256']==sha(HERE/'raw13_specialization.py') and raw['base_source_sha256']==sha(HERE.parent/'de-bondt-2610.24889/verify.py')
 assert normal['source_sha256']==sha(HERE/'raw13_normalized.py') and normal['input_sha256']==sha(rawpath)
 # Reconstruct the established representative directly, independently of SymPy.
 x,y,z=[variable(i,3) for i in range(3)];one={(0,)*3:Q(1)};v=add(one,mul(x,y),Q(2));v2=mul(v,v);v3=mul(v2,v)
 base=[add(add(x,mul(mul(x,x),y),Q(-3)),mul(mul(mul(x,x),x),z),Q(-1)),
       add(scale(mul(mul(v2,x),z),Q(3)),mul(add(add(scale(v2,Q(9)),v,Q(-6)),one,Q(-2)),y)),
       add(mul(v3,z),mul(add(scale(v2,Q(6)),v,Q(2)),mul(y,y)))]
 assert [sum(sum(e)>=2 for e in f) for f in base]==[2,5,6]
 expected={(i,e):c for i,f in enumerate(base) for e,c in f.items() if sum(e)>=2};assert len(expected)==13
 assert len(atoms)==13 and {(a['row'],tuple(a['exponent'])):Q(a['coefficient']) for a in atoms}==expected
 assert determinant3([[derivative(f,j) for j in range(3)] for f in base],3)==one
 points=[[Q(0),Q(1),Q(-8)],[Q(-1),Q(1),Q(4)]];images=[[evaluate(f,p) for f in base] for p in points];assert images[0]==images[1]==[Q(0),Q(1),Q(0)]
 # Universal coefficient pattern:3 source variables plus14 parameter variables.
 N=17;F=[variable(i,N) for i in range(3)]
 for j,atom in enumerate(atoms):
  e=list(atom['exponent'])+[0]*14;e[3+j]=1;F[atom['row']]=add(F[atom['row']],{tuple(e):Q(1)})
 det=add(determinant3([[derivative(f,j) for j in range(3)] for f in F],N),{(0,)*N:Q(1)},Q(-1));buckets={}
 for e,c in det.items():buckets.setdefault(e[:3],{})[e[3:]]=c
 def at_second(f):
  p={}
  for e,c in f.items():
   a=list(e[3:]);a[13]+=e[2];key=tuple(a);p[key]=p.get(key,Q(0))+c*((-1)**e[0])
  return {e:c for e,c in p.items() if c}
 actual=list(buckets.values())+[at_second(F[0]),add(at_second(F[1]),{(0,)*14:Q(1)},Q(-1))]
 rawgens=[decode(f,14) for f in raw['equations']];assert Counter(map(canonical,actual))==Counter(map(canonical,rawgens))
 # Independently derive the normalized system by substitution in every raw generator.
 one13={(0,)*13:Q(1)};values=[scale(one13,Q(-1)),scale(one13,Q(-3))]+[variable(i,13) for i in range(11)]+[add(one13,variable(11,13),Q(3))]
 derived=[substitute(f,values,13) for f in rawgens[:-2]];derived=[f for f in derived if f]
 normalized_point=[scale(one13,Q(-1)),variable(11,13),add(one13,variable(11,13),Q(3))]+values
 assert not substitute(F[0],normalized_point,13)
 derived.append(add(substitute(F[1],normalized_point,13),variable(11,13),Q(-1)))
 derived.append(add(one13,mul(variable(11,13),variable(12,13)),Q(-1)))
 ngens=[decode(f,13) for f in normal['equations']];assert Counter(map(canonical,derived))==Counter(map(canonical,ngens))
 baseline=[Q(a['coefficient']) for a in atoms[2:]]+[Q(1),Q(1)];assert all(evaluate(f,baseline)==0 for f in ngens)
 manifest=json.loads((out/'certification-direct.json').read_text());assert manifest['status']=='exact-primary-identities' and manifest['source_sha256']==sha(HERE/'parse_raw13_certificates.py')
 assert manifest['raw_input_sha256']==sha(rawpath) and manifest['normalized_input_sha256']==sha(normalpath)
 identities=0;maxdegree=0;monomials=0
 for tag,record in manifest['certificates'].items():
  path=out/record['path'];assert record['sha256']==sha(path);certificate=json.loads(path.read_text())
  assert certificate['source_sha256']==sha(HERE/'parse_raw13_certificates.py') and certificate['exporter_sha256']==sha(HERE/'certify_raw13_rigidity.py')
  assert certificate['singular_source_sha256']==sha(out/(tag+'.sing')) and certificate['singular_log_sha256']==sha(out/(tag+'.log'))
  if tag.startswith('raw-unit-'):
   i=int(tag[-1]);assert i in [0,1];variables=14;generators=rawgens+[variable(i,14)];targets=[{(0,)*14:Q(1)}]
  else:
   assert tag=='normalized-squares';variables=13;generators=ngens;targets=[mul(add(variable(i,13),one13,-b),add(variable(i,13),one13,-b)) for i,b in enumerate(baseline)]
  assert len(certificate['certificates'])==len(targets)
  for cert,target in zip(certificate['certificates'],targets):
   assert decode(cert['target'],variables)==target;multipliers=[decode(f,variables) for f in cert['multipliers']];assert len(multipliers)==len(generators)
   result={}
   for f,h in zip(generators,multipliers):result=add(result,mul(f,h))
   assert result==target,(tag,identities);identities+=1;monomials+=sum(len(h) for h in multipliers);maxdegree=max(maxdegree,max((sum(e) for h in multipliers for e in h),default=0))
 assert identities==15
 # Verify the one-parameter Laurent family and its complete collision algebraically.
 family=[variable(i,4) for i in range(3)]
 for atom in atoms:
  row=atom['row'];sourcee=atom['exponent'];k=sourcee[2]-int(row==2);e=tuple(sourcee+[k]);coefficient=Q(atom['coefficient'])*(Q(-1) if k%2 else Q(1));family[row]=add(family[row],{e:coefficient})
 assert determinant3([[derivative(f,j) for j in range(3)] for f in family],4)=={(0,)*4:Q(1)}
 def laurent_image(f,first,third):
  result={}
  for e,c in f.items():
   if first==0 and e[0]:continue
   coefficient=c*(Q(first)**e[0])*(Q(third)**e[2]);k=e[3]-e[2];result[k]=result.get(k,Q(0))+coefficient
  return {k:c for k,c in result.items() if c}
 for f,expected in zip(family,[{}, {0:Q(1)},{}]):assert laurent_image(f,0,8)==laurent_image(f,-1,-4)==expected
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/'verify.py'),helper_sha256=sha(HERE/'verify_power_single_compress.py'),raw_input_sha256=sha(rawpath),normalized_input_sha256=sha(normalpath),certification_sha256=sha(out/'certification-direct.json'),
  baseline_nonlinear_counts=[2,5,6],baseline_total=13,baseline_dimension=3,baseline_degree=7,raw_generators=len(rawgens),normalized_generators=len(ngens),identities=identities,max_multiplier_degree=maxdegree,multiplier_monomials=monomials,
  universal_determinant_reconstructed=True,normalized_generators_reconstructed=True,baseline_collision_checked=True,all_polynomial_identities_multiplied=True,laurent_family_determinant_checked=True,laurent_family_collision_checked=True,seconds=round(time.monotonic()-started,3),
  scope='In the exact13-atom raw three-variable support pattern, every determinant-one map with a collision between(0,1,c) and(-1,1,d) is the one-parameter diagonal conjugate of the established representative. All13 coefficients are nonzero. Thus no coefficient deletion in this pattern and collision family gives12 or fewer terms. This is not a global raw support minimum or a normalized-axis sparsity bound; the representative is pre-existing, not a new construction.')
 (out/'validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
