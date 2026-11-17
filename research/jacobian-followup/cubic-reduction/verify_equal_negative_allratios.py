#!/usr/bin/env python3
"""Independent exact finite checks for the all-ratio equal-negative degree7 gate."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
from itertools import product
from collections import defaultdict
import hashlib,json
from verify import variable,add,scale,mul,compose,derivative,determinant3,decode
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def total(*polys):
 out={}
 for p in polys:out=add(out,p)
 return out

def power(p,k,n):
 out={(0,)*n:Q(1)}
 for j in range(k):out=mul(out,p)
 return out

def jac(f,g):return add(mul(derivative(f,1),derivative(g,2)),mul(derivative(f,2),derivative(g,1)),Q(-1))
def main():
 folder=HERE/'equal-negative-smallweights-v2';spth=folder/'structural-models.json';struct=json.loads(spth.read_text());assert struct['status']=='complete-primary-structural-identities' and struct['source_sha256']==sha(HERE/'equal_negative_structural.py')
 pairs=[(p,q) for d in range(3,7) for p in range(1,d) for q in [d-p] if gcd(p,q)==1];assert len(pairs)==10
 assert [(r['weights'][0],-r['weights'][1]) for r in struct['pairs']]==pairs
 identities=0;atoms=0;control_degrees=[]
 for rec,(pw,qw) in zip(struct['pairs'],pairs):
  layers=2 if pw+qw==3 else 1;degs={'Q':pw,'T':2*pw,'R':pw+1,'S':pw+1,'U':2*pw+1,'V':2*pw+1};names=[name+str(j) for name in ['Q','T','R','S','U','V'] for j in range(degs[name]+1) if layers==2 or name not in ['T','U','V']];assert rec['parameter_names']==names and rec['layers']==layers and rec['weights']==[pw,-qw,-qw]
  n=3+len(names);X=[variable(j,n) for j in range(3)];params={name:variable(3+j,n) for j,name in enumerate(names)};F=list(X);allowed=[]
  # Enumerate every weighted monomial independently of the producer's layer formula.
  for i,wgt in enumerate([pw,-qw,-qw]):
   for e in product(range(8),repeat=3):
    if not 2<=sum(e)<=7 or pw*e[0]-qw*(e[1]+e[2])!=wgt:continue
    if i==0:prefix='Q' if e[0]==qw+1 else 'T';assert e[0] in [qw+1,2*qw+1]
    else:prefix=('R' if i==1 else 'S') if e[0]==qw else ('U' if i==1 else 'V');assert e[0] in [qw,2*qw]
    name=prefix+str(e[2]);assert name in params;mon={e+(0,)*len(names):Q(1)};F[i]=add(F[i],mul(mon,params[name]));allowed.append((i,e,name))
  assert len(allowed)==len(names)==rec['nonlinear_atoms'];atoms+=len(allowed)
  assert F==[decode(p,n) for p in rec['outputs']]
  B={}
  for prefix,degree in degs.items():B[prefix]=total(*(mul(params[prefix+str(j)],mul(power(X[1],degree-j,n),power(X[2],j,n))) for j in range(degree+1))) if layers==2 or prefix not in ['T','U','V'] else {}
  assert B=={k:decode(p,n) for k,p in rec['blocks'].items()}
  Qb,T,R,S,U,V=[B[k] for k in ['Q','T','R','S','U','V']];J=determinant3([[derivative(p,j) for j in range(3)] for p in F],n);assert J.get((0,)*n)==1
  expected_names=['first','one-layer-top'] if layers==1 else ['first','one-layer-boundary-top','full-top','quartic-first-full-top','no-high-tail-top','one-high-tail-quartic-first','one-high-tail-mixed'];assert [d['name'] for d in rec['identities']]==expected_names
  expected=[total(derivative(R,1),derivative(S,2),scale(Qb,qw+1))]
  if layers==1:expected.append(scale(mul(Qb,jac(R,S)),Q(pw+qw+1,pw+1)))
  else:expected.extend([scale(mul(Qb,jac(R,S)),Q(pw+qw+1,pw+1)),scale(mul(T,jac(U,V)),Q(2*pw+2*qw+1,2*pw+1)),scale(mul(Qb,jac(U,V)),Q(2*pw+qw+1,2*pw+1)),scale(mul(T,jac(R,S)),Q(pw+2*qw+1,pw+1)),scale(add(scale(mul(Qb,jac(U,S)),pw+qw+1),scale(mul(U,jac(Qb,S)),-qw)),Q(1,pw+1)),scale(add(scale(mul(T,jac(U,S)),pw+2*qw+1),scale(mul(U,jac(T,S)),-qw)),Q(1,pw+1))])
  for d,wanted in zip(rec['identities'],expected):
   zero_indices={3+names.index(name) for prefix in d['zero_blocks'] for name in names if name.startswith(prefix)};actual={}
   for e,c in J.items():
    if e[0]!=d['a_power'] or any(e[j] for j in zero_indices):continue
    key=(0,)+e[1:];actual[key]=actual.get(key,Q(0))+c
   actual={e:c for e,c in actual.items() if c};assert actual==decode(d['coefficient'],n)==wanted;identities+=1
  assert rec['degree']==1+layers*(pw+qw)<=7
  # Valid non-affine triangular controls ensure the gate is not asserting all maps are linear.
  xx=[variable(j,3) for j in range(3)];shear=mul(power(xx[0],qw,3),power(xx[2],pw+1,3))
  if layers==2:shear=add(shear,mul(power(xx[0],2*qw,3),power(xx[2],2*pw+1,3)))
  H=[xx[0],add(xx[1],shear),xx[2]];inv=[xx[0],add(xx[1],shear,Q(-1)),xx[2]];assert [compose(p,inv,3) for p in H]==xx and [compose(p,H,3) for p in inv]==xx
  assert determinant3([[derivative(p,j) for j in range(3)] for p in H],3)=={(0,)*3:Q(1)};control_degrees.append(max(sum(e) for p in H for e in p))
 assert atoms==147 and identities==30
 # Complete UFD valuation possibilities over the algebraic closure, including coincident lines.
 valuation_specs=[('211-quartic',(2,3,5),(1,6,4),[(0,2,3),(2,1,2),(2,3,5)]),('211-mixed',(4,3,5),(1,7,5),[(1,2,3),(3,1,2),(4,3,5)]),('122-quartic',(1,2,3),(2,5,4),[(1,2,3)]),('122-mixed',(2,2,3),(2,7,6),[(2,2,3)])]
 valuation_records=[]
 for name,bounds,weights,expected_states in valuation_specs:
  states=[e for e in product(*(range(k+1) for k in bounds)) if any(e) and weights[0]*e[0]+weights[1]*e[1]==weights[2]*e[2]];assert states==expected_states
  partitions=[counts for counts in product(range(max(bounds)+1),repeat=len(states)) if all(sum(k*state[j] for k,state in zip(counts,states))==bounds[j] for j in range(3))]
  assert len(partitions)==(2 if name.startswith('211') else 1);valuation_records.append({'name':name,'states':states,'degree_partitions':partitions})
 checkpoint=json.loads((folder/'checkpoint.json').read_text());liftcheck=json.loads((folder/'lift-checkpoint.json').read_text());assert checkpoint['status']==liftcheck['status']=='complete' and checkpoint['source_sha256']==sha(HERE/'equal_negative_smallweights_v2.py') and liftcheck['source_sha256']==sha(HERE/'certify_equal_negative_smallweights_v2.py') and liftcheck['producer_checkpoint_sha256']==sha(folder/'checkpoint.json')
 cases={r['name']:r for r in checkpoint['cases']};lifts={r['name']:r for r in liftcheck['cases']};five={'bb':('b','b'),'bc':('b','c'),'b-bc':('b','b+c'),'cc':('c','c'),'cb':('c','b')};expected_cases={(p,q,family,label) for p,q in [(2,1),(1,2)] for family in ['quartic-first','mixed'] for label in (five if p==2 else ['b','c'])};assert len(expected_cases)==len(cases)==len(lifts)==14
 records=[];hashes={};terms_total=0;generators_total=0;maxdeg=0
 for pw,qw,family,label in sorted(expected_cases):
  name='p%dq%d-%s-%s'%(pw,qw,family,label);meta=cases[name];mp=folder/meta['model'];m=json.loads(mp.read_text());assert m['name']==name and m['weights']==[pw,-qw,-qw] and m['family']==family and m['chart']==label and m['source_sha256']==sha(HERE/'equal_negative_smallweights_v2.py') and meta['model_sha256']==sha(mp)
  names=(['q','s','t','r','rho'] if family=='quartic-first' else (['q0','q1','q2','s','t','h','r','rho'] if pw==2 else ['q0','q1','s','t','h','r','rho']));assert m['parameter_names']==names;n=3+len(names);xx=[variable(j,n) for j in range(3)];a,b,c=xx;pv={v:variable(3+i,n) for i,v in enumerate(names)};LL={'b':b,'c':c,'b+c':add(b,c)};L,M=(LL[v] for v in (five[label] if pw==2 else (label,label)));s,t=pv['s'],pv['t']
  if family=='quartic-first':
   T={}
   if pw==2:Qb=mul(pv['q'],power(L,2,n));S=mul(s,mul(L,power(M,2,n)));U=mul(t,mul(power(L,2,n),power(M,3,n)))
   else:Qb=mul(pv['q'],L);S=mul(s,power(L,2,n));U=mul(t,power(L,3,n))
   saturated=mul(pv['rho'],mul(pv['q'],mul(s,t)))
  else:
   if pw==2:Qb=total(mul(pv['q0'],power(b,2,n)),mul(pv['q1'],mul(b,c)),mul(pv['q2'],power(c,2,n)));S=mul(s,mul(power(L,2,n),M));T=mul(t,mul(L,power(M,3,n)));U=mul(pv['h'],mul(power(L,3,n),power(M,2,n)))
   else:Qb=add(mul(pv['q0'],b),mul(pv['q1'],c));S=mul(s,power(L,2,n));T=mul(t,power(L,2,n));U=mul(pv['h'],power(L,3,n))
   saturated=mul(pv['rho'],mul(s,mul(t,pv['h'])))
  rb=scale(add(derivative(S,2),scale(Qb,qw+1)),-1);R={}
  for e,coef in rb.items():f=list(e);f[1]+=1;R[tuple(f)]=coef/Q(f[1])
  R=add(R,mul(pv['r'],power(c,pw+1,n)));F=[total(a,mul(power(a,qw+1,n),Qb),mul(power(a,2*qw+1,n),T)),total(b,mul(power(a,qw,n),R),mul(power(a,2*qw,n),U)),add(c,mul(power(a,qw,n),S))]
  assert F==[decode(p,n) for p in m['outputs']]
  for block,p in [('Q',Qb),('S',S),('T',T),('U',U),('R',R)]:assert p=={(0,)+e:coef for e,coef in decode(m[block],n-1).items()}
  J=add(determinant3([[derivative(p,j) for j in range(3)] for p in F],n),{(0,)*n:Q(1)},Q(-1));top=(4 if family=='quartic-first' else 5)*qw;assert all(e[0]!=top for e in J)
  groups=defaultdict(dict)
  for e,coef in J.items():groups[e[:3]][e[3:]]=coef
  keys=sorted(groups);gens=[groups[k] for k in keys];keys=[list(k) for k in keys];sat={e[3:]:coef for e,coef in add(saturated,{(0,)*n:Q(1)},Q(-1)).items()};gens.append(sat);keys.append('saturation');assert keys==m['generator_monomials'] and gens==[decode(row,len(names)) for row in m['generators']]
  cp=folder/(name+'-unit.json');cert=json.loads(cp.read_text());script=folder/(name+'.sing');assert cert['status']=='exact-primary-unit-identity' and cert['source_sha256']==sha(HERE/'certify_equal_negative_smallweights_v2.py') and cert['model_sha256']==sha(mp) and cert['singular_source_sha256']==sha(script) and lifts[name]['certificate_sha256']==sha(cp)
  multipliers=[decode(row,len(names)) for row in cert['multipliers']];assert len(multipliers)==len(gens);unit=total(*(mul(x,y) for x,y in zip(multipliers,gens)));assert unit=={(0,)*len(names):Q(1)}
  terms=sum(map(len,multipliers));degree=max(sum(e) for row in multipliers for e in row);assert terms==lifts[name]['multiplier_terms'] and degree==lifts[name]['max_multiplier_degree'];terms_total+=terms;generators_total+=len(gens);maxdeg=max(maxdeg,degree);hashes[name]={'model_sha256':sha(mp),'unit_sha256':sha(cp),'singular_sha256':sha(script)};records.append({'name':name,'generators':len(gens),'multiplier_terms':terms,'max_multiplier_degree':degree});print('UNIT PASS',name,len(gens),terms,flush=True)
 assert (generators_total,terms_total,maxdeg)==(145,376,13)
 # Explicit prior premises retain their written geometric/plane trust boundaries.
 premises={}
 for claim,rel in [('GLDC5','small-component-gate'),('KCUB23','equal-negative-completion')]:
  base=HERE.parent/rel;vp=base/'validation.json';v=json.loads(vp.read_text());assert v['status']=='pass' and v['source_sha256']==sha(base/'verify.py') and v['proof_sha256']==sha(base/'PROOF.md');premises[claim]={'validation_sha256':sha(vp),'proof_sha256':sha(base/'PROOF.md'),'checker_sha256':sha(base/'verify.py')}
 receipt={'status':'pass','claim':'KCUB32','checker_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'structural_producer_sha256':sha(HERE/'equal_negative_structural.py'),'structural_models_sha256':sha(spth),'normal_producer_sha256':sha(HERE/'equal_negative_smallweights_v2.py'),'unit_producer_sha256':sha(HERE/'certify_equal_negative_smallweights_v2.py'),'checkpoint_sha256':sha(folder/'checkpoint.json'),'lift_checkpoint_sha256':sha(folder/'lift-checkpoint.json'),'primitive_nonlinear_pairs_checked':pairs,'complete_weighted_nonlinear_atoms':atoms,'generic_top_and_first_identities':identities,'triangular_inverse_controls':len(control_degrees),'control_degrees':control_degrees,'valuation_cases':valuation_records,'normal_charts':len(records),'coefficient_generators':generators_total,'exact_unit_multiplier_terms':terms_total,'maximum_multiplier_degree':maxdeg,'cases':records,'case_hashes':hashes,'imported_premises':premises,'all_primitive_ratios_except_one_one_closed_through7':True,'one_one_strict_branch_imported_from_KCUB23':True,'one_one_full_degree7_not_closed':True,'first_pass_three_off_locus_models_not_used':True,'global_profile_counts_unchanged':[62,12],'positive_bounds_unchanged':True,'goal_complete':False,'trust_boundary':'Finite exact reconstruction and unit identities plus written weight-lattice, homogeneous-power/UFD, Borel-cover and plane-fibre arguments. GLDC5/Moh100 and KCUB23 strict-only are explicit imports; no formal or external-review claim.'}
 proof=HERE/'EQUAL_NEGATIVE_ALL_RATIOS.md'
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (folder/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS all primitive equal-negative ratios except1/1 throughdegree7;14 units,145 generators,376 multiplier terms; all-ratio strict sector closes only by importing KCUB23; global goal open',flush=True)
if __name__=='__main__':main()
