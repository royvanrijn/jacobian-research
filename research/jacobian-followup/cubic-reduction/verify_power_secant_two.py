#!/usr/bin/env python3
"""Independent Fraction replay of the unrestricted two-new-cube obstruction."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations,combinations_with_replacement
from collections import defaultdict
import hashlib,json,time
from verify import add,mul,derivative
from verify_power_single_compress import decode,echelon,reduce,lin,scale
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 started=time.monotonic();out=HERE/'power-secant-two';source=HERE/'power-basis-anneal/candidate.json';d=json.loads(source.read_text());n=d['size'];assert n==89
 spanpath=HERE/'power-span-split/validation.json';span=json.loads(spanpath.read_text());assert span['status']=='pass' and span['source_sha256']==sha(HERE/'verify_power_span_optimality.py') and span['input_sha256']==sha(source)
 B=[[Q(c) for c in row] for row in d['B']];forms=d['forms'];cubes=[mul(mul(lin(f),lin(f)),lin(f)) for f in forms];basis={}
 for j,cube in enumerate(cubes):
  p=dict(cube);coordinate={j:Q(1)}
  for e,(b,bc) in sorted(basis.items()):
   c=p.get(e,Q(0))
   if c:p=add(p,b,-c);coordinate=add(coordinate,bc,-c)
  assert p;e=min(p);c=p[e];basis[e]=(scale(p,1/c),scale(coordinate,1/c))
 cache={}
 def project_unit(e):
  if e in cache:return cache[e]
  p={e:Q(1)};coords={}
  for pivot,(b,bc) in sorted(basis.items()):
   c=p.get(pivot,Q(0))
   if c:p=add(p,b,-c);coords=add(coords,bc,c)
  cache[e]=(p,coords);return p,coords
 groups=defaultdict(list)
 for j in range(n):
  v=tuple(row[j] for row in B);c=next(a for a in v if a);groups[tuple(a/c for a in v)].append((j,c))
 triples=[tuple(t) for ids in groups.values() for t in combinations(ids,3)];assert len(triples)==214
 supports=[json.loads((HERE/'power-single-compress'/('case-%02d.json'%i)).read_text())['support'] for i in range(15)]
 adj=[0]*18
 for f in forms:
  for i,j in combinations([i for i,c in enumerate(f) if c],2):adj[i]|=1<<j;adj[j]|=1<<i
 nonedges=[(i,j) for i,j in combinations(range(18),2) if not adj[i]&(1<<j)]
 for i,j in nonedges:assert all(not derivative(derivative(cube,i),j) for cube in cubes)
 masks=[]
 for S in supports:
  assert all(adj[i]&(1<<j) for i,j in combinations(S,2))
  assert all(any(not adj[v]&(1<<u) for u in S) for v in range(18) if v not in S)
  masks.append(sum(1<<i for i in S))
 clique=bytearray(1<<18);clique[0]=1;covered=0
 for mask in range(1,1<<18):
  bit=mask&-mask;i=bit.bit_length()-1;rest=mask^bit
  if clique[rest] and rest&~adj[i]==0:clique[mask]=1;covered+=1;assert any(mask&~c==0 for c in masks)
 primary=json.loads((out/'checkpoint.json').read_text());manifest=json.loads((out/'certification.json').read_text())
 assert primary['status']=='primary-all-unit' and primary['source_sha256']==sha(HERE/'power_secant_two.py')
 assert manifest['status']=='exact-primary-identities' and manifest['source_sha256']==sha(HERE/'certify_power_secant_two.py') and manifest['primary_checkpoint_sha256']==sha(out/'checkpoint.json')
 expected_pairs=list(combinations_with_replacement(range(15),2));assert [tuple(r['pair']) for r in primary['cases']]==expected_pairs
 cert_by_pair={tuple(r['pair']):r for r in manifest['cases']};identities=0;maxdegree=0;allcases=[];zero_cases=0;terms_count=0
 for entry in primary['cases']:
  i,j=entry['pair'];S,T=supports[i],supports[j];m=len(S);variables=m+len(T);residual={};actual=[{} for k in range(89)]
  for support,offset in [(S,0),(T,m)]:
   for ids in combinations_with_replacement(range(len(support)),3):
    e=[0]*18;u=[0]*variables
    for a in ids:e[support[a]]+=1;u[offset+a]+=1
    mult=1 if ids[0]==ids[2] else (3 if ids[0]==ids[1] or ids[1]==ids[2] else 6);r,c=project_unit(tuple(e));u=tuple(u)
    for exponent,value in r.items():residual.setdefault(exponent,{})[u]=Q(mult)*value
    for k,value in c.items():actual[k][u]=Q(mult)*value
  J=echelon(residual.values())
  if entry['status']=='coordinate-zero-exclusion':
   reduced=[reduce(f,J) for f in actual]
   assert not any(all(reduced[a] for a,c in triple) for triple in triples);zero_cases+=1;allcases.append(dict(pair=[i,j],status='coordinate-zero-exclusion'));continue
  assert entry['status']=='unit';tag='pair-%02d-%02d'%(i,j);casepath=out/(tag+'.json');case=json.loads(casepath.read_text());assert entry['input_sha256']==sha(casepath)
  assert case['pair']==[i,j] and case['supports']==[S,T] and len(case['variables'])==variables
  constraints=[decode(f,variables) for f in case['membership']];coords=[decode(f,variables) for f in case['coordinates']];assert len(coords)==89
  assert all(sum(e)==3 for f in constraints+coords for e in f)
  I=echelon(constraints);assert len(I)==len(J)
  assert all(not reduce(f,I) for f in residual.values()) and all(not reduce(f,J) for f in constraints)
  assert all(not reduce(add(a,b,Q(-1)),I) for a,b in zip(coords,actual))
  selected=[triple for triple in triples if all(coords[a] for a,c in triple)]
  encoded=[tuple((a,Q(c)) for a,c in triple) for triple in case['triples']]
  assert set(encoded)==set(selected) and len(encoded)==len(selected)==entry['systems']
  certpath=out/(tag+'-cert.json');certificate=json.loads(certpath.read_text());assert cert_by_pair[(i,j)]['certificate_sha256']==sha(certpath)
  assert certificate['input_sha256']==sha(casepath) and certificate['source_sha256']==sha(HERE/'certify_power_secant_two.py') and certificate['singular_sha256']==sha(out/(tag+'-cert.sing'))
  assert len(certificate['certificates'])==len(selected);one={(0,)*variables:Q(1)}
  for triple,cert in zip(encoded,certificate['certificates']):
   assert cert['indices']==[a for a,c in triple]
   generators=constraints+[add(coords[a],one,-c) for a,c in triple];multipliers=[decode(f,variables) for f in cert['multipliers']];assert len(generators)==len(multipliers)
   result={}
   for g,h in zip(generators,multipliers):result=add(result,mul(g,h))
   assert result==one,(i,j,triple);identities+=1;terms_count+=sum(len(f) for f in multipliers);maxdegree=max(maxdegree,max((sum(e) for f in multipliers for e in f),default=0))
  allcases.append(dict(pair=[i,j],identities=len(encoded),membership_rank=len(I),certificate_sha256=sha(certpath)))
  if len(allcases)%20==0:print('VERIFIED',len(allcases),'SUPPORT PAIRS',identities,'IDENTITIES',flush=True)
 assert len(allcases)==120 and zero_cases==7 and identities==primary['systems']==manifest['identities']==607
 assert len(cert_by_pair)==113
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/'verify.py'),helper_sha256=sha(HERE/'verify_power_single_compress.py'),
  input_sha256=sha(source),fixed_span_validation_sha256=sha(spanpath),primary_checkpoint_sha256=sha(out/'checkpoint.json'),certification_sha256=sha(out/'certification.json'),
  source_cliques=15,source_support_pairs=120,coordinate_zero_pairs=7,nonedge_pairs=len(nonedges),covered_nonempty_cliques=covered,eligible_proportional_triples=214,
  unit_identities=identities,max_multiplier_degree=maxdegree,multiplier_monomials=terms_count,cases=allcases,
  all_nonedge_mixed_derivatives_zero=True,all_support_pairs_covered=True,all_membership_spaces_reconstructed=True,all_coordinate_maps_verified=True,all_unit_identities_multiplied=True,
  seconds=round(time.monotonic()-started,3),scope='The saved89-form identity cannot be shortened by retaining any subset, adding at most two arbitrary new forms and re-solving all retained output coefficients, over characteristic-zero extensions. Fixed-span optimality handles independent new output columns; the mixed second derivative clique lemma and607 identities handle proportional columns. This is not a global minimum89 for the tensor or for Keller maps; at least three new forms remain possible.')
 (out/'validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps({k:v for k,v in receipt.items() if k!='cases'}),flush=True)
if __name__=='__main__':main()
