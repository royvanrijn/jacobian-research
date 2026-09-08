#!/usr/bin/env python3
"""Post-barrier evaluator. Full parent is forbidden until all17 output hashes pass."""
import hashlib,json,time,resource
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[4];PKG=ROOT/'research/artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1'
def resolve(p):
 if p.exists():return p
 try:rel=p.relative_to(PKG/'runs')
 except ValueError:return p
 return ROOT/'research/artifacts/local/elliptic-curves/det1092_blind_mw16_v1/runs'/rel
def sha(p):return hashlib.sha256(resolve(p).read_bytes()).hexdigest()
def put(p,d):
 with p.open('x') as f:json.dump(d,f,indent=2,sort_keys=True);f.write('\n')
 p.chmod(0o444)
barrier=json.loads((PKG/'primary_complete.json').read_text());assert len(barrier['runs'])==17
for arm in barrier['runs']:
 for p,h in arm['files'].items():assert sha(PKG/p)==h
protocol=json.loads((PKG/'protocol.json').read_text());commit=protocol['starting_commit']
parent_rel='research/artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json'
pp=PKG/'frozen_sources'/parent_rel
assert sha(pp)==protocol['source_bindings'][parent_rel]
# First permitted read of the unredacted parent in this stage.
d=json.loads(pp.read_text());R=PolynomialRing(QQ,'t');K=R.fraction_field()
def dec(x):return K(R(x['numerator']))/R(x['denominator'])
def ht(P):return ZZ(0) if P.is_zero() else ZZ(max(P[0].numerator().degree(),P[0].denominator().degree()+4))
def pairing(P,Q):return (ht(P)+ht(Q)-ht(P-Q))/2
E=EllipticCurve(K,list(map(dec,d['a_invariants'])))
B=[E(list(map(dec,p))) for p in d['basis_weierstrass_coordinates']]
G=matrix(QQ,17,lambda i,j:pairing(B[i],B[j]));assert G==matrix(QQ,d['generic_height_gram']) and G.det()==1092
assert G.is_positive_definite()
roster=json.loads((PKG/'roster.json').read_text())['arms'];rows=[];start=time.monotonic()
for arm in roster:
 aid=arm['arm_id'];i=arm['omitted_basis_index_one_based']-1
 result=json.loads(resolve(PKG/'runs'/aid/'result.json').read_text());fixture=json.loads((PKG/arm['fixture']).read_text())
 keep=[j for j in range(17) if j!=i]
 assert fixture['sections']==[d['basis_weierstrass_coordinates'][j] for j in keep]
 assert matrix(QQ,fixture['gram'])==G.matrix_from_rows_and_columns(keep,keep)
 row={'arm_id':aid,'omitted_basis_index_one_based':i+1,'terminal_status':result['terminal_status'],'candidates':[],
      'retained_determinant':str(matrix(QQ,fixture['gram']).det()),'missing_orthogonal_height':str(QQ(1092)/matrix(QQ,fixture['gram']).det()),
      'resource':json.loads(resolve(PKG/'runs'/aid/'accounting.json').read_text()),'primary_result_sha256':sha(PKG/'runs'/aid/'result.json')}
 for ci,c in enumerate(result.get('candidates',[])):
  Q=E(list(map(dec,c['coordinates'])));v=vector(QQ,[pairing(P,Q) for P in B]);w=G.solve_right(v)
  residual=ht(Q)-w*G*w
  assert residual==0,('UNEXPECTED_GENERIC_OVERGROUP',aid,ci,residual)
  assert all(x.denominator()==1 for x in w),('SATURATED_BASIS_INCONSISTENCY',aid,ci)
  assert sum((n*P for n,P in zip(w,B)),E(0))==Q
  gain=bool(w[i]);assert c['rank']==(17 if gain else 16)
  relation=('EQUALS_OMITTED_SECTION' if Q==B[i] else 'NEGATIVE_OMITTED_SECTION' if Q==-B[i] else 'DIFFERS_BY_RETAINED_SUBGROUP' if w[i]==1 else 'NEGATIVE_QUOTIENT_PLUS_RETAINED_SUBGROUP' if w[i]==-1 else 'NONPRIMITIVE_QUOTIENT_MULTIPLE' if gain else 'RETAINED_SUBGROUP')
  row['candidates'].append({'candidate_id':ci,'full_basis_coordinates':list(map(str,w)),'omitted_coefficient':str(w[i]),'quotient_classification':'Generic quotient recovery' if gain else 'No generic gain','representative_relation':relation,'exact_group_identity':True,'orthogonal_residual_height':str(residual)})
 gains=[x for x in row['candidates'] if x['quotient_classification']=='Generic quotient recovery']
 row['success']=bool(gains);row['recovered_generic_rank']=17 if gains else 16
 row['quotient_classification']='Generic quotient recovery' if gains else 'No generic gain'
 row['first_success_candidate_id']=gains[0]['candidate_id'] if gains else None
 if gains:
  c=result['candidates'][gains[0]['candidate_id']];ci,n,li=c['producer']
  selection=json.loads(resolve(PKG/'runs'/aid/'selection.json').read_text())
  row['observed_geometry']={'first_success_centre_index_zero_based':ci,'centre_height':selection['centres'][ci]['height'],'RR_degree_n':n,'kernel_member_index':li,'schur_complement':c['schur_complement'],'quotient_coefficient':gains[0]['omitted_coefficient']}
 else:row['observed_geometry']={'completed_RR_stages':sum(s['status'].startswith('COMPLETE') for s in result.get('stages',[])),'returned_sections':len(result.get('candidates',[]))}
 rows.append(row);print(aid,row['quotient_classification'],flush=True)
count=sum(r['success'] for r in rows);failures=[r['arm_id'] for r in rows if r['terminal_status'] not in ['GENERIC_QUOTIENT_RECOVERY','BOUNDED_NO_RECOVERY']]
conclusion='BROAD_16_TO_17_RECOVERY' if count==17 else 'CORE_DEPENDENT_16_TO_17_RECOVERY' if count else 'UNRESOLVED_DUE_TO_EXPERIMENT_FAILURE' if failures else 'BOUNDED_NO_RECOVERY'
put(PKG/'evaluation.json',{'starting_commit':commit,'primary_barrier_sha256':sha(PKG/'primary_complete.json'),'parent_sha256':sha(pp),'evaluator_sha256':sha(Path(__file__)),'successes':count,'conclusion':conclusion,'failed_or_censored_arms':failures,'rows':rows,'resource':{'wall_seconds':time.monotonic()-start,'user_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime,'max_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},'boundary':'Known generic17 positive-control reconstruction only. Alternative representatives and all misses retained. No exceptional302 point input; no rank change.'})
print(conclusion,count,'of17',flush=True)
