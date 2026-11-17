#!/usr/bin/env python3
"""Independent rational reconstruction and universal indicator certificates."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from collections import Counter,defaultdict
import hashlib,json
from verify import decode,variable,add,scale,mul,compose,identity,madd,mmul,derivative,determinant3,evaluate
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def psum(*ps):
 out={}
 for p in ps:out=add(out,p)
 return out

def main():
 path=HERE/'cubic12-integer-witness.json';raw=json.loads(path.read_text());F=[decode(p,12) for p in raw['coefficients']];mp=HERE/'joint-unipotent5/model.json';m=json.loads(mp.read_text());assert m['status']=='complete-symbolic-family' and m['source_sha256']==sha(HERE/'joint_unipotent5.py') and m['input_sha256']==sha(path)
 P=[variable(j,5) for j in range(5)];u,v,w,z,t=P;one={(0,)*5:Q(1)};const=lambda q:scale(one,Q(q));I=identity(12,5);S=[[dict(p) for p in row] for row in I];Si=[[dict(p) for p in row] for row in I]
 for i,j,p in [(1,6,u),(1,9,v),(2,8,w),(2,7,add(z,mul(w,t))),(8,7,t)]:S[i][j]=p
 for i,j,p in [(1,6,u),(1,9,v),(2,8,w),(2,7,z),(8,7,t)]:Si[i][j]=scale(p,-1)
 assert S==[[decode(p,5) for p in row] for row in m['target_matrix']] and Si==[[decode(p,5) for p in row] for row in m['source_matrix']]
 assert mmul(S,Si)==I and mmul(Si,S)==I;A=madd(S,I,Q(-1));assert mmul(mmul(A,A),A)==[[{} for j in range(12)] for i in range(12)]
 # Also verify an elementary determinant-one factorization, without a determinant oracle.
 E=I
 for i,j,p in [(1,6,u),(1,9,v),(2,8,w),(2,7,z),(8,7,t)]:
  T=[[dict(q) for q in row] for row in I];T[i][j]=p;E=mmul(E,T)
 assert E==S
 nv=17;X=[variable(j,nv) for j in range(12)];embedp=lambda p:{(0,)*12+e:c for e,c in p.items()};embedf=lambda p:{e+(0,)*5:c for e,c in p.items()};sub=[]
 for row in Si:sub.append(psum(*(mul(embedp(c),x) for c,x in zip(row,X))))
 sub.extend(variable(12+j,nv) for j in range(5));R=[compose(embedf(p),sub,nv) for p in F];G=[psum(*(mul(embedp(c),p) for c,p in zip(row,R))) for row in S]
 assert R==[decode(p,nv) for p in m['right_source_coefficients']] and G==[decode(p,nv) for p in m['coefficients']]
 assert all({e:c for e,c in p.items() if sum(e[:12])<=1}==X[i] for i,p in enumerate(G))
 assert all(sum(e[:12])<=3 for p in G for e in p)
 cube=[0]*17;cube[1]=3;assert G[4][tuple(cube)]==1
 factors=[decode(f['polynomial'],5) for f in m['factors']];assert len(factors)==19 and [f['id'] for f in m['factors']]==list(range(19))
 for j,p in [(0,w),(1,z),(2,v),(8,u),(12,t)]:assert factors[j]==p
 def factored(rec):
  out=const(rec['factor_scalar']);assert out
  for a in rec['factors']:
   assert 0<=a['factor']<19 and a['power']>=1
   for k in range(a['power']):out=mul(out,factors[a['factor']])
  return out
 def coefficient_groups(polys):
  out={}
  for i,p in enumerate(polys):
   groups=defaultdict(dict)
   for e,c in p.items():
    if sum(e[:12])>=2:groups[e[:12]][e[12:]]=c
   for e,q in groups.items():out[(i,e)]=q
  return out
 slots=coefficient_groups(G);assert len(slots)==m['generic_nonlinear_total']==90
 actual_slots={(r['output']-1,tuple(r['source_exponent'])):r for r in m['nonlinear_slots']};assert set(actual_slots)==set(slots)
 inv=Counter()
 for key,rec in actual_slots.items():
  assert slots[key]==decode(rec['parameter_polynomial'],5)==factored(rec);inv[tuple((a['factor'],a['power']) for a in rec['factors'])]+=1
 assert inv==Counter({tuple((a['factor'],a['power']) for a in r['factors']):r['count'] for r in m['factor_inventory']})
 rightslots=coefficient_groups(R);columns=defaultdict(dict)
 for (i,e),q in rightslots.items():columns[e][i]=q
 source_inv=Counter();source_records={tuple(r['source_exponent']):r for r in m['source_columns']};assert set(columns)==set(source_records) and len(columns)==68
 for e,rec in source_records.items():
  gcd=decode(rec['gcd_polynomial'],5);assert gcd==factored(rec);quotients={int(j)-1:decode(q,5) for j,q in rec['quotients'].items()};assert set(quotients)==set(columns[e])
  for i,q in quotients.items():assert mul(gcd,q)==columns[e][i]
  unit={};assert rec['quotient_unit'] is not None
  assert len(rec['quotient_unit'])==len({a['output'] for a in rec['quotient_unit']})
  for a in rec['quotient_unit']:assert Q(a['coefficient']);unit=add(unit,quotients[a['output']-1],Q(a['coefficient']))
  assert unit==one;source_inv[tuple((a['factor'],a['power']) for a in rec['factors'])]+=1
 assert source_inv==Counter({tuple((a['factor'],a['power']) for a in r['factors']):r['count'] for r in m['source_factor_inventory']})
 def collapse(inv):
  out=Counter()
  for ids,count in inv.items():out[tuple(j for j,k in ids)]+=count
  return out
 ci=collapse(inv);cs=collapse(source_inv)
 assert [ci[k] for k in [(),(0,),(1,),(2,),(8,)]]==[24,5,5,10,10]
 assert [cs[k] for k in [(),(0,),(1,),(2,),(8,)]]==[24,4,4,9,9]
 def cost(inventory,flags):return sum(count for ids,count in inventory.items() if all(flags[j] for j in ids))
 def flagmap(mapping,bits):return [a if isinstance(a,bool) else all(bits[j] for j in a) for a in mapping]
 def verify_mapping(sub,bases,mapping):
  for f,a in zip(factors,mapping):
   q=compose(f,sub,5);wanted=(one if a else {}) if isinstance(a,bool) else one
   if not isinstance(a,bool):
    for j in a:wanted=mul(wanted,bases[j])
   if not wanted:assert not q
   else:assert q and scale(wanted,q[min(q)]/wanted[min(wanted)])==q
 # Complete coefficient identities on u=v=0; all2^9 indicator assignments are checked algebraically.
 ba=[w,z,add(z,const(-4)),psum(mul(t,w),scale(w,-3),z,const(8)),t,add(w,const(-2)),add(scale(t,2),z),add(w,const(4)),add(t,const(-3))]
 ma=[(0,),(1,),False,True,True,True,(2,),True,False,True,True,(3,),(4,),(5,),(6,),True,(4,),(7,),(8,)]
 verify_mapping([{}, {},w,z,t],ba,ma)
 assert add(ba[5],ba[7],Q(-1))==const(-6) and add(ba[4],ba[8],Q(-1))==const(3)
 zero_wz=[u,v,{}, {},t]
 assert compose(ba[3],zero_wz,5)==const(8) and compose(ba[6],zero_wz,5)==scale(t,2)
 assert compose(ba[5],[u,v,{},z,t],5)==const(-2) and compose(ba[7],[u,v,{},z,t],5)==const(4) and compose(ba[2],[u,v,w,{},t],5)==const(-4)
 accepted=0;minimum=999;next_minimum=999
 for bits in product([False,True],repeat=9):
  flags=flagmap(ma,bits);N=cost(ci,flags);C=cost(cs,flags)
  assert N==31+6*bits[0]+5*bits[1]+bits[2]+bits[3]+3*bits[4]+bits[5]+bits[6]+bits[7]+bits[8]
  assert C==29+4*bits[0]+4*bits[1]+bits[5]+2*bits[4]+bits[6]
  if not(bits[5] or bits[7]) or not(bits[4] or bits[8]):continue
  if not bits[0] and not(bits[5] and bits[7]):continue
  if not bits[1] and not bits[2]:continue
  if not bits[0] and not bits[1] and (not bits[3] or bits[6]!=bits[4]):continue
  accepted+=1;minimum=min(minimum,N)
  if bits[0] or bits[1] or bits[4]:assert N>=39 and C>=33;next_minimum=min(next_minimum,N)
  else:assert N==36 and C==30
 assert minimum==36 and next_minimum==39
 # The other two axes are entire polynomial families, not sampled parameter points.
 bb=[u,add(u,const(Q(1,3))),add(u,const(3)),add(u,const(Q(1,6))),t,add(t,const(-3))]
 mb=[False,False,False,(1,),True,(2,),True,(1,),(0,),True,(3,),True,(4,),True,(4,),(1,),(4,1),True,(5,)]
 verify_mapping([u,{}, {},{},t],bb,mb)
 bc=[v,add(v,const(-1)),add(v,const(12)),add(v,const(Q(-1,2))),t,add(t,const(-3))]
 mc=[False,False,(0,),True,(1,),True,(2,),(1,),False,(3,),True,True,(4,),True,(4,),True,(4,),True,(5,)]
 verify_mapping([{},v,{}, {},t],bc,mc)
 for rest in product([False,True],repeat=5):
  bits=(True,)+rest;assert cost(ci,flagmap(mb,bits))==42+4*bits[1]+bits[2]+bits[3]+3*bits[4]+bits[4]*bits[1]+3*bits[5]
  assert cost(ci,flagmap(mc,bits))==44+3*bits[1]+bits[2]+bits[3]+5*bits[4]+3*bits[5]
 # Global parameter cases: both u,v nonzero, or one plus w/z, are excluded by disjoint mandatory groups.
 for bu,bv,bw,bz in product([False,True],repeat=4):
  lb=24+10*bu+10*bv+5*bw+5*bz
  if bu and bv:assert lb>=44
  elif (bu or bv) and (bw or bz):assert lb>=39
 # For v!=0,u=0, these three extra source columns are fixed and disjoint from the mandatory24+9.
 for j in [3,10,15]:assert cs[(j,)]>=1 and set(compose(factors[j],[{},v,w,z,t],5))=={(0,)*5}
 # At v=0, three source columns always form a pairwise-nonparallel projective triangle.
 def exponent(**kwargs):
  e=[0]*12
  for j,k in kwargs.items():e[int(j[1:])-1]=k
  return tuple(e)
 triangle=[exponent(x1=1,x2=1),exponent(x1=2,x3=1),exponent(x1=1,x10=1)]
 wanted=[{3:const(3),11:one},{3:one},{11:one}]
 for e,expect in zip(triangle,wanted):
  actual={i:compose(q,[u,{},w,z,t],5) for i,q in columns[e].items()};actual={i:q for i,q in actual.items() if q};assert actual==expect
 points=[[Q(v) for v in p] for p in raw['points']];image=list(map(Q,raw['common_image']));assert points[0][0]!=points[1][0]
 for point in points:
  pp=[embedp(psum(*(scale(c,a) for c,a in zip(row,point)))) for row in S];pp.extend(variable(12+j,17) for j in range(5))
  for i,p in enumerate(G):assert compose(p,pp,17)==({(0,)*17:image[i]} if image[i] else {})
 # Two matching controls, each with a full all-variable determinant check.
 def specialize(p,values):
  out={}
  for e,c in p.items():
   for a,k in zip(values,e[12:]):c*=a**k
   out[e[:12]]=out.get(e[:12],Q(0))+c
  return {e:c for e,c in out.items() if c}
 def numeric_mm(A,B):return [[sum(a*B[k][j] for k,a in enumerate(row)) for j in range(len(B[0]))] for row in A]
 def determinant_check(H):
  J=[[derivative(p,j) for j in range(12)] for p in H];D=[row[3:] for row in J[3:]];Id=identity(9,12);N=madd(D,Id,Q(-1));assert all(not N[i][j] for i in range(9) for j in range(i,9));inverse=Id;power=Id
  for k in range(1,9):power=mmul(power,N);inverse=madd(inverse,power,Q((-1)**k))
  assert mmul(power,N)==[[{} for j in range(9)] for i in range(9)] and mmul(D,inverse)==Id and mmul(inverse,D)==Id
  schur=madd([row[:3] for row in J[:3]],mmul(mmul([row[3:] for row in J[:3]],inverse),[row[:3] for row in J[3:]]),Q(-1));assert determinant3(schur,12)=={(0,)*12:Q(1)}
 values=[Q(0)]*4+[Q(3)];spoint=[[evaluate(p,values) for p in row] for row in S];ipoint=[[evaluate(p,values) for p in row] for row in Si];target=[[Q(i==j) for j in range(12)] for i in range(12)]
 for i,j,c in [(1,6,-3),(1,9,-12),(2,8,-4),(2,7,-20),(8,7,3),(10,9,3)]:target[i][j]=Q(c)
 matching=[]
 assert len(m['examples'])==2 and {a['label'] for a in m['examples']}=={'identity-linear-t3','free-output-t3'}
 for example in m['examples']:
  assert list(map(Q,example['parameters']))==values and [[Q(v) for v in row] for row in example['source_matrix']]==ipoint
  if example['label']=='identity-linear-t3':H=[specialize(p,values) for p in G];Tout=spoint;expected_count=39
  else:
   assert example['label']=='free-output-t3';rp=[specialize(p,values) for p in R];H=[psum(*(scale(p,c) for p,c in zip(rp,row))) for row in target];Tout=target;expected_count=34
  assert H==[decode(p,12) for p in example['coefficients']] and [[Q(v) for v in row] for row in example['target_matrix']]==Tout
  linear=numeric_mm(Tout,ipoint);assert [[derivative(p,j).get((0,)*12,Q(0)) for j in range(12)] for p in H]==linear
  nonlinear=[sum(sum(e)>=2 for e in p) for p in H];linearc=[sum(sum(e)==1 for e in p) for p in H]
  assert nonlinear==example['nonlinear_counts'] and sum(nonlinear)==example['nonlinear_total']==expected_count
  assert linearc==example['linear_counts'] and sum(linearc)==example['linear_total'] and sum(linearc)+sum(nonlinear)==example['total']==51
  assert len(set(e for p in H for e in p if sum(e)>=2))==example['nonlinear_source_columns']==33
  pp=[[sum(a*b for a,b in zip(row,point)) for row in spoint] for point in points];assert pp==[list(map(Q,p)) for p in example['points']] and pp[0]!=pp[1]
  assert all([evaluate(p,point) for p in H]==list(map(Q,example['common_image'])) for point in pp)
  assert all(c.denominator==1 for p in H for c in p.values());determinant_check(H);matching.append({'label':example['label'],'nonlinear':expected_count,'linear':sum(linearc),'full_symbolic_determinant_one':True,'integral_collision':True})
 # Previously established controls/premises are pinned, not new searches.
 assert [specialize(p,[Q(0)]*5) for p in G]==F
 prevpath=HERE/'output-frame31/validation.json';prev=json.loads(prevpath.read_text());assert prev['status']=='pass' and prev['checker_sha256']==sha(HERE/'verify_output31.py') and prev['proof_sha256']==sha(HERE/'OUTPUT_FRAME_TRADEOFF.md') and prev['input_sha256']==sha(path) and prev['arithmetic_sha256']==sha(HERE/'verify.py')
 prevw=HERE/'output-frame31/witness.json';assert prev['witness_sha256']==sha(prevw);pw=json.loads(prevw.read_text());assert [specialize(p,list(map(Q,[-3,-12,-4,-8,3]))) for p in G]==[decode(p,12) for p in pw['identity_normalized_conjugate']['coefficients']]
 oldpath=HERE/'integer-validation12.json';old=json.loads(oldpath.read_text());assert old['status']=='pass' and old['witness_sha256']==sha(path) and old['source_sha256']==sha(HERE/'verify_integer12.py') and old['arithmetic_sha256']==sha(HERE/'verify.py')
 rigidpath=HERE/'target-derivation12/validation.json';rigid=json.loads(rigidpath.read_text());assert rigid['status']=='pass' and rigid['input_sha256']==sha(path) and rigid['checker_sha256']==sha(HERE/'verify_target_rigidity12.py') and rigid['proof_sha256']==sha(HERE/'TARGET_CUBIC_RIGIDITY.md') and rigid['arithmetic_sha256']==sha(HERE/'verify.py')
 receipt={'status':'pass','claim':'KCUB31','checker_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'model_sha256':sha(mp),'producer_sha256':sha(HERE/'joint_unipotent5.py'),'input_sha256':sha(path),'base_validation_sha256':sha(oldpath),'output31_validation_sha256':sha(prevpath),'target_rigidity_validation_sha256':sha(rigidpath),'all_parameter_matrices_inverse_and_elementary':True,'all_generic_coefficients_and_collisions':True,'generic_nonlinear_slots':90,'coefficient_factor_products':90,'factor_polynomials':19,'generic_nonlinear_source_columns':68,'source_column_factorizations':68,'source_column_quotient_units':68,'branch_A_boolean_assignments_checked':512,'branch_A_necessary_compatible_assignments':accepted,'branch_B_boolean_assignments_checked':32,'branch_C_boolean_assignments_checked':32,'normalized_minimum':36,'normalized_minimum_unique_identity':True,'normalized_nonidentity_minimum':39,'normalized_37_38_excluded':True,'nonlinear_source_column_minimum':30,'nonidentity_source_column_minimum':33,'free_output_minimum_importing_KCUB30':31,'nonidentity_source_free_output_minimum':34,'matching_controls':matching,'known_identity36_and_clearing74_controls':True,'degree_at_most5_target_corollary_imports_KCUB26':True,'goal_complete':False,'boundary':'Only the specified five-parameter source family, with prescribed inverse normalization or arbitrary affine target mixing. No full source-frame/global optimum or improved headline witness.'}
 proof=HERE/'JOINT_UNIPOTENT_LIMITS.md'
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (HERE/'joint-unipotent5/validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS complete five-parameter family: normalized36 only at identity, next39; source columns30/33; arbitrary output minima31/34; two full determinants and integral collisions',flush=True)
if __name__=='__main__':main()
