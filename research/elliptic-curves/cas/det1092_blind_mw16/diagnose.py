#!/usr/bin/env python3
"""Secondary diagnostic only after primary freeze, evaluation AND independent replay.

No point search, exceptional point input or feedback to primary arms. The finite
atlas is17 singleton and272 signed-pair centres,153 parity classes per basis.
"""
import hashlib,json,time,statistics,resource
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[4];PKG=ROOT/'research/artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1'
def resolve(p):
 if p.exists():return p
 try:rel=p.relative_to(PKG/'runs')
 except ValueError:return p
 return ROOT/'research/artifacts/local/elliptic-curves/det1092_blind_mw16_v1/runs'/rel
def sha(p):return hashlib.sha256(resolve(p).read_bytes()).hexdigest()
replay=json.loads((PKG/'replay.json').read_text());assert replay['status']=='PASS_INDEPENDENT_GENERIC_REPLAY'
ev=json.loads((PKG/'evaluation.json').read_text());assert sha(PKG/'evaluation.json')==replay['evaluation_sha256']
protocol=json.loads((PKG/'protocol.json').read_text());commit=protocol['starting_commit']
parent_path=PKG/'frozen_sources/research/artifacts/generated-results/elliptic-curves/curve302_recovered_mw17_parent_v1.json'
parent=json.loads(parent_path.read_text());R=PolynomialRing(QQ,'t');K=R.fraction_field()
def dec(v):return K(R(v['numerator']))/R(v['denominator'])
E=EllipticCurve(QQ,[dec(a)(0) for a in parent['a_invariants']]);B=[E([dec(c)(0) for c in p]) for p in parent['basis_weierstrass_coordinates']]
G=matrix(QQ,parent['generic_height_gram']);a1,a2,a3,a4,a6=E.a_invariants();b2=a1*a1+4*a2
A=a4+a1*a3/2-b2*b2/48
identity=matrix.identity(ZZ,17)
def bits(q):return max(abs(q.numerator()).nbits(),q.denominator().nbits())
def atlas(points,U):
 rows=[];parities=set()
 words=[]
 for j in range(17):words.append((identity.row(j),points[j]))
 for j in range(17):
  for k in range(j+1,17):
   for s in [-1,1]:words.append((identity.row(j)+s*identity.row(k),points[j]+s*points[k]))
 for w,P in words:
  v=w*U;parity=tuple(int(x%2) for x in v);parities.add(parity)
  X=P[0]+b2/12;Y=P[1]+(a1*P[0]+a3)/2
  quartic=[-3*X*X-4*A,-8*Y,-6*X,QQ(0),QQ(1)]
  # Raw slope quartic; no coordinate reduction, searched point or visibility
  # result is inferred from coefficient bit size.
  rows.append({'basis_word':list(map(int,w)),'full_basis_word':list(map(int,v)),'generic_centre_height':str(v*G*v),'raw_quartic_max_rational_bits':int(max(map(bits,quartic)))})
 costs=[r['raw_quartic_max_rational_bits'] for r in rows]
 return {'centres':rows,'parity_masks':sorted(sum(v<<i for i,v in enumerate(p)) for p in parities),'cost_min':min(costs),'cost_median':statistics.median(costs),'cost_max':max(costs)}
start=time.monotonic();original=atlas(B,identity);rows=[]
for row in ev['rows']:
 aid=row['arm_id'];i=row['omitted_basis_index_one_based']-1;diag={'arm_id':aid,'generic_recovery':row['success']}
 if not row['success']:
  diag['status']='NO_RECOVERED_GENERIC_QUOTIENT_TO_SPECIALIZE';rows.append(diag);continue
 result=json.loads(resolve(PKG/'runs'/aid/'result.json').read_text());j=row['first_success_candidate_id'];candidate=result['candidates'][j]
 Q=E([dec(c)(0) for c in candidate['coordinates']]);w=vector(ZZ,row['candidates'][j]['full_basis_coordinates'])
 assert sum((n*P for n,P in zip(w,B)),E(0))==Q
 embedding=matrix(ZZ,parent['basis_embedding_in_public_D']);publicword=embedding*w
 # The original public D independence and injective specialization remain
 # bound to their canonical parent proof; no public point coordinates used.
 diag.update(status='PASS_EXACT_SPECIALIZATION_QUOTIENT',specialized_point=list(map(str,Q.xy())),full_generic_basis_word=list(map(str,w)),public_D_word=list(map(str,publicword)),expected_missing_generic_quotient=str(w[i]),subgroup_index_in_full_generic_specialization=int(abs(w[i])))
 if abs(w[i])==1:
  U=matrix(ZZ,identity);U.set_row(i,w);assert abs(U.det())==1
  alternative=list(B);alternative[i]=Q;data=atlas(alternative,U)
  old=set(original['parity_masks']);new=set(data['parity_masks'])
  diag.update(same_integral_MW17_specialization=True,intrinsic_half_lattice_unchanged=True,finite_atlas=data,parity_overlap_with_original=len(old&new),new_parity_classes=len(new-old),lost_parity_classes=len(old-new))
 else:diag.update(same_integral_MW17_specialization=False,intrinsic_half_lattice_unchanged=None)
 rows.append(diag)
# Post-run observational stretch only: no trained selector and no additional
# recovery arm. Determinant and execution cost are available generically.
ordered=sorted([r for r in ev['rows'] if r['success']],key=lambda r:r['resource']['wall_seconds'])
largest=max(ev['rows'],key=lambda r:QQ(r['retained_determinant']))
maxdet_diagnostic={'rule':'Choose the largest exact retained Gram determinant, breaking ties by arm_id; evaluated post-run only, without a new reconstruction.','selected_arm':largest['arm_id'],'retained_determinant':largest['retained_determinant'],'success_in_frozen_primary':largest['success'],'wall_seconds':largest['resource']['wall_seconds'],'cost_rank_among_successes':next((i+1 for i,r in enumerate(ordered) if r['arm_id']==largest['arm_id']),None),'interpretation':'For this known parent the missing orthogonal height is1092/det(core). The score uses only core data; its apparent predictive utility is post-hoc and unvalidated on other parents.'}
out={'starting_commit':commit,'status':'PASS_SECONDARY_DIAGNOSTIC','checker_sha256':sha(Path(__file__)),'primary_replay_sha256':sha(PKG/'replay.json'),'evaluation_sha256':sha(PKG/'evaluation.json'),'diagnostic_policy':'Fixed289 raw chord centres,17 singleton and272 signed pairs,153 parity classes; equation at t=0, no point search, no quartic reduction, no exceptional points. This post-run footprint compares representations, not empirical discovery visibility.','original_atlas':original,'arms':rows,'stretch_observation':{'fastest_recovered_core':ordered[0]['arm_id'] if ordered else None,'max_determinant_posthoc_diagnostic':maxdet_diagnostic,'generic_only_selector_validated':False,'reason':'Execution cost and core determinants can be compared after the fixed experiment; fitting a chooser on this known parent would be training data. No further optimization or additional arms were run.'},'resource':{'wall_seconds':time.monotonic()-start,'user_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime,'max_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},'boundary':'Same integral subgroup implies identical intrinsic half-lattice. Finite atlas coverage and raw coefficient sizes can change with basis. Whether changed charts discover exceptional directions is UNMEASURED.'}
with (PKG/'diagnostic.json').open('x') as f:json.dump(out,f,indent=2,sort_keys=True);f.write('\n')
(PKG/'diagnostic.json').chmod(0o444)
print(out['status'],flush=True)
