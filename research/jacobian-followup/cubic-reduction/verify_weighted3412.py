#!/usr/bin/env python3
"""Independent source, complete-stratum, unit and matching upper-control replay."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import json,hashlib
from verify import decode,variable,add,scale,mul,compose,derivative,determinant3,evaluate
HERE=Path(__file__).resolve().parent;OUT=HERE/'weighted-312-412';ROOT=HERE.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def total(*ps):
 out={}
 for p in ps:out=add(out,p)
 return out

def con(c,n):return {(0,)*n:Q(c)} if c else {}
def power(p,k,n):
 out=con(1,n)
 for _ in range(k):out=mul(out,p)
 return out

def cof(p,i,k):return {e[:i]+e[i+1:]:c for e,c in p.items() if e[i]==k}
def bracket(f,g):return add(mul(derivative(f,0),derivative(g,1)),mul(derivative(f,1),derivative(g,0)),Q(-1))
def dictionary(f,g,h,p):return total(mul(f,bracket(g,h)),scale(mul(g,bracket(f,h)),Q(1,p)),scale(mul(h,bracket(g,f)),Q(2,p)))
def main():
 cp=OUT/'checkpoint.json';check=json.loads(cp.read_text());assert check['status']=='complete' and check['source_sha256']==sha(HERE/'weighted3412.py');models={};model_hashes={};atoms=0;ngens=0
 for rec in check['models']:
  p=rec['positive_weight'];ip=OUT/rec['input'];d=json.loads(ip.read_text());assert rec['input_sha256']==sha(ip) and d['source_sha256']==sha(HERE/'weighted3412.py');basis=[]
  for i in range(3):
   for degree in range(2,8):
    for a in range(degree+1):
     for b in range(degree-a+1):
      e=(a,b,degree-a-b)
      if p*e[0]-e[1]-2*e[2]==[p,-1,-2][i]:basis.append((i,e))
  basis.sort();N=len(basis);assert N=={3:15,4:14}[p] and d['basis']==[{'output':i,'exponent':list(e)} for i,e in basis];assert d['parameter_names']==['c'+str(j) for j in range(N)];n=3+N;x=[variable(j,n) for j in range(3)];F=list(x)
  for j,(i,e) in enumerate(basis):F[i]=add(F[i],mul({e+(0,)*N:Q(1)},variable(3+j,n)))
  assert F==[decode(z,n) for z in d['source_outputs']];ni=2+N;u,v=variable(0,ni),variable(1,ni);sub=[con(1,ni),u,v]+[variable(2+j,ni) for j in range(N)];forms=[compose(z,sub,ni) for z in F];assert forms==[decode(z,ni) for z in d['cover_outputs']];D=dictionary(*forms,p);assert D==decode(d['quotient_jacobian'],ni);lift={}
  for e,c in D.items():assert (e[0]+2*e[1])%p==0;lift[((e[0]+2*e[1])//p,e[0],e[1])+e[2:]]=c
  J=determinant3([[derivative(z,j) for j in range(3)] for z in F],n);assert lift==J;eqs=[{e[2:]:c for e,c in add(D,con(1,ni),Q(-1)).items() if e[:2]==m} for m in sorted({e[:2] for e in add(D,con(1,ni),Q(-1))},reverse=True)];assert eqs==[decode(z,N) for z in d['native_generators']];parent=HERE/'weight-resonance7/mixed-cases'/('p'+str(p)+'q1r2.json');assert d['parent_timeout_model_sha256']==sha(parent) and json.loads(parent.read_text())['basis']==d['basis'];models[p]=(d,F,forms,D,J,eqs);model_hashes[ip.name]=sha(ip);atoms+=N;ngens+=len(eqs);print('SOURCE PASS',p,N,len(eqs),flush=True)
 assert (atoms,ngens)==(29,52)
 # Two exact b=0 boundary factorizations, with L=0 imposed only for p3.
 for p in [3,4]:
  d,F,forms,D,J,eqs=models[p];N=len(d['parameter_names']);n=3+N;x=[variable(j,n) for j in range(3)];nr=1+N;t=variable(0,nr);P=[variable(1+j,nr) for j in range(N)];one=con(1,nr)
  if p==4:f0=total(one,mul(P[0],t),mul(P[3],power(t,2,nr)));g0=total(one,mul(P[4],t),mul(P[7],power(t,2,nr)));h0=total(one,mul(P[9],t),mul(P[13],power(t,2,nr)));aa,bb=1,2;T=mul(x[0],power(x[2],2,n));actual=J
  else:
   f0=add(one,mul(P[2],t));g0=add(one,mul(P[7],t));h0=add(one,mul(P[13],t));aa,bb=2,3;T=mul(power(x[0],2,n),power(x[2],3,n));S=[variable(j,n) for j in range(n)];S[3+4]={};actual=compose(J,S,n)
  expected=mul(g0,add(mul(h0,add(f0,scale(mul(t,derivative(f0,0)),aa))),scale(mul(t,mul(f0,derivative(h0,0))),bb)));assert compose(expected,[T]+[variable(3+j,n) for j in range(N)],n)=={e:c for e,c in actual.items() if e[1]==0}
 sp=OUT/'structural-identities.json';struct=json.loads(sp.read_text());assert struct['source_sha256']==sha(HERE/'weighted3412.py') and check['structural_identities_sha256']==sha(sp)
 # P4 reduced source model, complete specialization, all seven active supports.
 names=['A','B','C','D','E','G','L','M'];d=struct['p4_boundary_reduced'];assert d['parameters']==names;n=11;a,b,c=[variable(j,n) for j in range(3)];P={z:variable(3+j,n) for j,z in enumerate(names)};U=mul(a,power(b,4,n));V=mul(a,mul(power(b,2,n),c));one=con(1,n);F=[mul(a,total(one,mul(P['A'],V),mul(P['B'],U))),mul(b,total(one,mul(P['C'],V),mul(P['D'],U))),total(mul(c,total(one,mul(P['E'],V),mul(P['G'],U))),mul(power(b,2,n),add(P['L'],mul(P['M'],U))))];assert F==[decode(z,n) for z in d['source_outputs']];J=determinant3([[derivative(z,j) for j in range(3)] for z in F],n);JJ={}
 for e,co in J.items():assert e[0]>=e[2] and e[1]==4*(e[0]-e[2])+2*e[2];JJ[(e[0]-e[2],e[2])+e[3:]]=co
 assert JJ==decode(d['invariant_jacobian'],10);_,full,_,_,_,_=models[4];nf=17;S=[variable(j,nf) for j in range(nf)]
 for j in [0,3,4,7,9,13]:S[3+j]={}
 RF=[compose(z,S,nf) for z in full];assert RF==[compose(z,[variable(j,nf) for j in range(3)]+[variable(3+j,nf) for j in [1,2,5,6,10,11,8,12]],nf) for z in F]
 diag=[{e[2:]:co for e,co in JJ.items() if e[:2]==(0,j)} for j in [1,2,3]];pp={z:variable(i,8) for i,z in enumerate(names)};expected=[total(scale(pp['A'],2),scale(pp['C'],3),scale(pp['E'],2)),total(scale(mul(pp['A'],pp['C']),4),scale(mul(pp['A'],pp['E']),3),scale(mul(pp['C'],pp['E']),4)),scale(mul(pp['A'],mul(pp['C'],pp['E'])),5)];assert diag==expected==[decode(z,8) for z in d['diagonal_coefficients']]
 for mask in range(1,8):
  bits=[int(bool(mask&(1<<j))) for j in range(3)];S=[variable(j,8) for j in range(8)];target=con(1+bits[0]+2*bits[1]+bits[2],8)
  for z,on in zip(['A','C','E'],bits):
   if on:target=mul(target,pp[z])
   else:S[names.index(z)]={}
  assert compose(diag[sum(bits)-1],S,8)==target
 ni=10;S=[variable(j,ni) for j in range(ni)]
 for z in ['A','C','E']:S[2+names.index(z)]={}
 pure=compose(JJ,S,ni);assert pure=={(e[0],0)+e[1:]:co for e,co in decode(d['pure_U_factorization'],9).items()};uu=variable(0,ni);BB,DD,GG=[variable(2+names.index(z),ni) for z in ['B','D','G']];rhs=mul(add(con(1,ni),mul(GG,uu)),total(con(1,ni),mul(add(scale(BB,2),scale(DD,5)),uu),scale(mul(mul(BB,DD),power(uu,2,ni)),6)));assert pure==rhs
 # P3 L0 full reduction and product derivative at U=0.
 names=['mu','P','Q','M','N','B','S','H','J','K','X'];d=struct['p3_Lzero_boundary_reduced'];assert d['parameters']==names;n=14;a,b,c=[variable(j,n) for j in range(3)];P={z:variable(3+j,n) for j,z in enumerate(names)};U=mul(a,power(b,3,n));V=mul(a,mul(b,c));one=con(1,n);f=total(one,mul(P['mu'],V),mul(P['P'],U),mul(P['Q'],power(V,2,n)));g=total(one,mul(P['M'],V),mul(P['N'],U),mul(P['B'],power(V,2,n)));h=total(one,mul(P['H'],V),mul(P['J'],U),mul(P['X'],power(V,2,n)));F=[mul(a,f),mul(b,g),add(mul(c,h),mul(power(b,2,n),add(P['S'],mul(P['K'],U))))];assert F==[decode(z,n) for z in d['source_outputs']];J=determinant3([[derivative(z,j) for j in range(3)] for z in F],n);JJ={}
 for e,co in J.items():assert e[0]>=e[2] and e[1]==3*(e[0]-e[2])+e[2];JJ[(e[0]-e[2],e[2])+e[3:]]=co
 assert JJ==decode(d['invariant_jacobian'],13);_,full,_,_,_,_=models[3];nf=18;S=[variable(j,nf) for j in range(nf)]
 for j in [4,2,7,13]:S[3+j]={}
 assert [compose(z,S,nf) for z in full]==[compose(z,[variable(j,nf) for j in range(3)]+[variable(3+j,nf) for j in [0,1,3,5,6,8,9,10,11,12,14]],nf) for z in F];ni=12;t=variable(0,ni);pp={z:variable(1+j,ni) for j,z in enumerate(names)};prod=con(1,ni)
 for A,B in [('mu','Q'),('M','B'),('H','X')]:prod=mul(prod,total(con(1,ni),mul(pp[A],t),mul(pp[B],power(t,2,ni))))
 assert cof(JJ,0,0)==decode(d['U_zero_derivative_product'],ni)==derivative(mul(t,prod),0);S=[variable(j,13) for j in range(13)]
 for z in ['mu','Q','M','B','H','X']:S[2+names.index(z)]={}
 pure=compose(JJ,S,13);assert pure=={(e[0],0)+e[1:]:co for e,co in decode(d['pure_U_factorization'],12).items()};uu=variable(0,13);PP,NN,HH=[variable(2+names.index(z),13) for z in ['P','N','J']];rhs=mul(add(con(1,13),mul(HH,uu)),total(con(1,13),mul(add(scale(PP,2),scale(NN,4)),uu),scale(mul(mul(PP,NN),power(uu,2,13)),5)));assert pure==rhs
 # The b=0 cubic for the hard312 stratum and its complete four-branch cover.
 _,_,_,D,_,original=models[3];zero=cof(D,0,0);bc=[cof(zero,0,3*j) for j in [1,2,3]];assert all(e[0]%3==0 for e in zero);n=15;C=[variable(j,n) for j in range(n)];mu,T,L,H,A,W=C[0],C[2],C[4],C[10],C[7],C[13];one=con(1,n);expected=[total(scale(T,3),A,scale(W,4),scale(mul(L,H),-2),scale(mul(mu,L),-1)),total(scale(mul(A,W),4),scale(mul(T,A),3),scale(mul(T,W),6),scale(mul(T,mul(L,H)),-3)),scale(mul(T,mul(A,W)),6)];assert bc==expected;S=[variable(j,n) for j in range(n)];S[0]=one;S[4]=one;S[10]=scale(total(scale(T,3),A,scale(W,4),scale(one,-1)),Q(1,2));assert compose(bc[0],S,n)=={};red=compose(bc[1],S,n);wanted=add(mul(A,total(scale(W,8),scale(T,3))),scale(mul(T,add(scale(T,3),one,Q(-1))),3),Q(-1));assert scale(red,2)==wanted
 # Residual diagonal conjugacy (a,b,c)->(a/t,t*b,t*c) keeps L=1 and scales(mu,T) by t,(P,Q) by t^2.
 basis=json.loads((OUT/'p3.json').read_text())['basis'];chars=[(1 if z['output']==0 else -1)-z['exponent'][0]+z['exponent'][1]+z['exponent'][2] for z in basis];assert [chars[j] for j in [4,0,2,1,3]]==[0,1,1,2,2]
 # Reconstruct every constant-pivot model and independently multiply all seven unit identities.
 identity=con(1,15);model_pins={};unit_pins={};units=0;terms=0;unitgens=0;pivots=0;successful_pivots=0;maxdegree=0;timeout_records=[]
 for cpname,producer in [('strata-checkpoint.json','weighted312_strata.py'),('boundary-checkpoint.json','weighted312_boundary.py')]:
  checkpoint=json.loads((OUT/cpname).read_text());assert checkpoint['status']=='complete' and checkpoint['source_sha256']==sha(HERE/producer) and checkpoint['parent_model_sha256']==sha(OUT/'p3.json')
  for rec in checkpoint['cases']:
   name=rec['case'];ip=OUT/(name+'.json');d=json.loads(ip.read_text());assert rec['input_sha256']==sha(ip) and d['source_sha256']==sha(HERE/producer);fixed={int(j):(decode(v,15) if isinstance(v,list) else con(v,15)) for j,v in d['fixed_parameters'].items()};expected_fixed={'mu1':{4:1,0:1},'mu0T1':{4:1,0:0,2:1},'muT0P1':{4:1,0:0,2:0,1:1},'muTP0Q1':{4:1,0:0,2:0,1:0,3:1},'mu1-T0W0':{4:1,0:1,2:0,13:0},'mu1-T0A0':{4:1,0:1,2:0,7:0},'mu1-W0A3Tminus1':{4:1,0:1,13:0,7:'linear'},'mu1-A0Tthird':{4:1,0:1,7:0,2:Q(1,3)}}[name];assert set(fixed)==set(expected_fixed)
   for j,value in expected_fixed.items():assert fixed[j]==(add(scale(C[2],3),identity,Q(-1)) if value=='linear' else con(value,15))
   S=[fixed.get(j,variable(j,15)) for j in range(15)];active=[j for j in range(15) if j not in fixed]
   for step in d['pivot_steps']:
    j=step['parameter'];assert j in active;eq=compose(original[step['native_equation_index']],S,15);pivot=Q(step['constant_pivot']);assert pivot and derivative(eq,j)==con(pivot,15);expr=add(variable(j,15),scale(eq,-1/pivot));assert expr==decode(step['expression'],15);R=[variable(k,15) for k in range(15)];R[j]=expr;S=[compose(z,R,15) for z in S];active.remove(j);pivots+=1
   nf=len(active);assert d['free_parameter_names']==['c'+str(j) for j in active];R=[variable(active.index(j),nf) if j in active else {} for j in range(15)];SS=[compose(z,R,nf) for z in S];assert SS==[decode(z,nf) for z in d['substitutions']];allgens=[compose(z,SS,nf) for z in original];indices=[j for j,z in enumerate(allgens) if z];gens=[allgens[j] for j in indices];assert indices==d['native_generator_indices'] and gens==[decode(z,nf) for z in d['generators']];script=OUT/(name+'.sing');assert sha(script)==rec['script_sha256'];model_pins[name]=sha(ip)
   if rec['status']=='changed-reduced-timeout10':assert name=='mu1';timeout_records.append({'case':name,'status':rec['status'],'checkpoint_sha256':sha(OUT/cpname)});continue
   assert rec['status']=='exact-primary-unit';certp=OUT/(name+'-unit.json');cert=json.loads(certp.read_text());assert sha(certp)==rec['certificate_sha256'] and cert['source_sha256']==sha(HERE/producer) and cert['input_sha256']==sha(ip) and cert['script_sha256']==sha(script);mm=[decode(z,nf) for z in cert['multipliers']];assert len(mm)==len(gens) and total(*(mul(a,b) for a,b in zip(mm,gens)))==con(1,nf);tt=sum(map(len,mm));assert tt==rec['multiplier_terms'];units+=1;terms+=tt;unitgens+=len(gens);successful_pivots+=len(d['pivot_steps']);maxdegree=max(maxdegree,max((sum(e) for z in mm for e in z),default=0));unit_pins[name]=sha(certp);print('UNIT PASS',name,len(d['pivot_steps']),tt,flush=True)
 assert (units,terms,unitgens,pivots,successful_pivots,maxdegree)==(7,1096,123,48,42,7)
 # Generic shear inverse controls for the two directly classified boundary families.
 for p in [3,4]:
  n=5;a,b,c,L,M=[variable(j,n) for j in range(n)];shear=add(mul(L,power(b,2,n)),mul(M,mul(a,power(b,p+2,n))));F=[a,b,add(c,shear)];inv=[a,b,add(c,shear,Q(-1))];extra=[L,M];assert [compose(z,inv+extra,n) for z in F]==[a,b,c] and [compose(z,F+extra,n) for z in inv]==[a,b,c] and determinant3([[derivative(z,j) for j in range(3)] for z in F],n)==con(1,n)
 # Reconstruct the existing matching upper control directly; no new witness claim.
 a,b,c=[variable(j,3) for j in range(3)];one=con(1,3);v=add(one,scale(mul(a,b),2));v2=power(v,2,3);v3=power(v,3,3);upper=[total(a,scale(mul(power(a,2,3),b),-3),scale(mul(power(a,3,3),c),-1)),add(scale(mul(v2,mul(a,c)),3),mul(total(scale(v2,9),scale(v,-6),scale(one,-2)),b)),add(mul(v3,c),mul(add(scale(v2,6),scale(v,2)),power(b,2,3)))];assert [{e:co for e,co in z.items() if sum(e)<=1} for z in upper]==[a,b,c];assert [sum(sum(e)>=2 for e in z) for z in upper]==[2,5,6] and max(sum(e) for z in upper for e in z)==7;assert all(e[0]-e[1]-2*e[2]==w for z,w in zip(upper,[1,-1,-2]) for e in z);assert determinant3([[derivative(z,j) for j in range(3)] for z in upper],3)==one;points=[[Q(0),Q(1),Q(-8)],[Q(-1),Q(1),Q(4)]];assert [[evaluate(z,p) for z in upper] for p in points]==[[Q(0),Q(1),Q(0)]]*2;oldcontrol=HERE/'raw13-specialization/input.json';old=json.loads(oldcontrol.read_text());assert {(z['row'],tuple(z['exponent'])):Q(z['coefficient']) for z in old['atoms']}=={(i,e):co for i,z in enumerate(upper) for e,co in z.items() if sum(e)>=2}
 parent=HERE/'weighted-124-125/validation.json';prev=json.loads(parent.read_text());assert prev['status']=='pass' and len(prev['remaining_diagonal_weight_improvement_types'])==12;remaining=[z for z in prev['remaining_diagonal_weight_improvement_types'] if z['mixed_positive_negative_magnitudes'] not in [[3,1,2],[4,1,2]]];assert len(remaining)==10 and sum(z['mixed_positive_negative_magnitudes']==[1,1,1] for z in remaining)==1;ledger=json.loads((ROOT/'MATH_STATUS.json').read_text());entries={e['id']:e for e in ledger['entries']};imports={}
 for k in ['KCUB36','KCUB37','KCUB15','GLDC10']:
  e=entries[k];assert e['state']=='proved';imports[k]={'scope_fingerprint':hashlib.sha256((e['state']+'\0'+e['title']+'\0'+e['scope']).encode()).hexdigest()[:16],'proof_sha256':sha(ROOT/e['canonical_source'])}
 oldcp=HERE/'weight-resonance7/mixed-lift-checkpoint.json';old=json.loads(oldcp.read_text());assert len([z for z in old['cases'] if z['status']=='timeout'])==4;proof=HERE/'DISTINCT_NONZERO_WEIGHT_MINIMUM.md';receipt={'status':'pass','claim':'KCUB38','checker_sha256':sha(Path(__file__)),'model_producer_sha256':sha(HERE/'weighted3412.py'),'strata_producer_sha256':sha(HERE/'weighted312_strata.py'),'boundary_producer_sha256':sha(HERE/'weighted312_boundary.py'),'arithmetic_sha256':sha(HERE/'verify.py'),'checkpoint_sha256':sha(cp),'structural_identities_sha256':sha(sp),'model_hashes':model_hashes,'stratum_model_hashes':model_pins,'unit_hashes':unit_pins,'source_atoms':atoms,'native_coefficient_generators':ngens,'boundary_factorizations':2,'p4_diagonal_identities':3,'p4_active_support_cases':7,'p3_Lzero_product_derivative':True,'complete_normalization_strata':4,'mu1_complete_boundary_branches':4,'all_constant_pivots_verified':pivots,'successful_unit_constant_pivots':successful_pivots,'exact_unit_identities':units,'unit_generators':unitgens,'unit_multiplier_monomials':terms,'max_multiplier_degree':maxdegree,'changed_reduced_timeouts_preserved':timeout_records,'old_native_timeout_checkpoint_sha256':sha(oldcp),'generic_direct_shear_inverse_controls':2,'matching_existing_upper_control':{'source_sha256':sha(oldcontrol),'degree':7,'nonlinear_counts':[2,5,6],'weight':[1,-1,-2],'full_determinant_one':True,'points':[[str(z) for z in p] for p in points],'common_image':['0','1','0'],'new_witness':False},'parent_frontier_receipt_sha256':sha(parent),'remaining_diagonal_weight_improvement_types':remaining,'imports':imports,'sharp_class_scope':'Characteristic zero,displayed pairwise-distinct nonzero diagonal Gm weights: minimum collision dimension3 in every degree by the written rank-one consequence in dimensions<=2 plus the verified3D control. In dimension3,minimum ordinary collision degree7; minimum raw nonlinear occurrences13 with degree<=7. General source-frame changes,repeated/zero weights and larger-degree support are outside the term result.','minimum_dimension_for_this_additional_form':3,'dimension_lower_bound_is_written_rank_one_consequence':True,'goal_complete':False,'global_positive_bounds_and_unrestricted_profile_counts_unchanged':True,'trust_boundary':'Independent rational full-source/pivot/unit/witness reconstruction plus written complete normalization/field covers,boundary units,degree and imported GLDC10/KCUB classification arguments. Not formal or externally reviewed.'}
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (OUT/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS312/412 closure:29 atoms,52 native equations,48 constant pivots,7 exact units/1096 terms; existing13-term collision rechecked; distinct nonzero diagonal-weight minimum13 throughdegree7; full frontier12->10; goal open',flush=True)
if __name__=='__main__':main()
