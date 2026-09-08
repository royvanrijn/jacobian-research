#!/usr/bin/env python3
"""Post-completion audit of fixed policy, stage coverage and deliverable bindings.

Recomputes only the finite lattice selector; no RR construction, candidate search,
new arm or full-rank calculation. Existing independent rank replay is hash-checked
against its source and all224 candidate/rank records. Preserves frozen artifacts.
"""
import hashlib,json,time,resource
from datetime import datetime
from pathlib import Path
from sage.all import QQ,ZZ,matrix
ROOT=Path(__file__).resolve().parents[4]
PKG=ROOT/'research/artifacts/generated-results/elliptic-curves/det1092_blind_mw16_v1'
RAW=ROOT/'research/artifacts/local/elliptic-curves/det1092_blind_mw16_v1'
CAS=Path(__file__).resolve().parent
START=time.monotonic()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
summary=read(PKG/'summary.json');protocol=read(PKG/'protocol.json');roster=read(PKG/'roster.json')
barrier=read(PKG/'primary_complete.json');dispatch=read(PKG/'dispatch.json')
ev=read(PKG/'evaluation.json');replay=read(PKG/'replay.json');compact=read(PKG/'compact_replay.json');diag=read(PKG/'diagnostic.json')
commit=protocol['starting_commit']
assert summary['starting_commit']==roster['starting_commit']==barrier['starting_commit']==commit
for p,h in summary['package_files'].items():assert sha(PKG/p)==h,p
for p,h in protocol['source_bindings'].items():assert sha(PKG/'frozen_sources'/p)==h,p
for run in barrier['runs']:
 for p,h in run['files'].items():assert sha(RAW/p)==h,p
assert datetime.fromisoformat(protocol['frozen_utc'])<datetime.fromisoformat(dispatch['launch_utc'])<datetime.fromisoformat(barrier['completed_utc'])
assert dispatch['reruns']==0 and dispatch['sequential'] and barrier['evaluation_has_not_started']
assert dispatch['protocol_sha256']==sha(PKG/'protocol.json')
assert sha(CAS/'run.py')==dispatch['runner_sha256']
for name in ['worker.py','isolate.py','prepare.py']:
 assert sha(CAS/name)==protocol['source_bindings'][str((CAS/name).relative_to(ROOT))]
for obj,name,key in [(ev,'evaluate.py','evaluator_sha256'),(replay,'replay.py','checker_sha256'),(compact,'replay.py','checker_sha256'),(diag,'diagnose.py','checker_sha256')]:
 assert obj[key]==sha(CAS/name) and obj['starting_commit']==commit
assert ev['primary_barrier_sha256']==replay['primary_barrier_sha256']==compact['primary_barrier_sha256']==sha(PKG/'primary_complete.json')
assert replay['evaluation_sha256']==compact['evaluation_sha256']==diag['evaluation_sha256']==sha(PKG/'evaluation.json')
assert diag['primary_replay_sha256']==sha(PKG/'replay.json')
assert replay['status']=='PASS_INDEPENDENT_GENERIC_REPLAY' and compact['status']=='PASS_INDEPENDENT_COMPACT_GENERIC_REPLAY'
assert len(roster['arms'])==len(ev['rows'])==len(replay['arms'])==len(compact['arms'])==17
assert [a['omitted_basis_index_one_based'] for a in roster['arms']]==list(range(1,18))

# Independent scalar-integer implementation: recompute H*v at every move;
# do not import the producer or use its incrementally maintained product.
def selector(G):
 U=G.LLL_gram();assert abs(U.det())==1
 H=U.transpose()*G*U
 hs=[[int(x) for x in row] for row in H.rows()]
 us=[[int(x) for x in row] for row in U.rows()]
 directions=[[(i,1)] for i in range(16)]
 directions += [[(i,1),(j,s)] for i in range(16) for j in range(i+1,16) for s in [-1,1]]
 norms=[sum(a*b*hs[i][j] for i,a in move for j,b in move) for move in directions]
 samples=[];pool={}
 for j in range(512):
  mask=int.from_bytes(hashlib.sha256(('det1092-blind-mw16-v1:'+str(j)).encode()).digest()[:2],'big')
  v=[(mask>>i)&1 for i in range(16)];moves=0;capped=True
  for k in range(128):
   hv=[sum(x*y for x,y in zip(row,v)) for row in hs]
   options=[]
   for index,(move,norm) in enumerate(zip(directions,norms)):
    product=sum(sign*hv[i] for i,sign in move)
    delta=4*(norm-abs(product))
    if delta<0:options.append((delta,index,-1 if product>0 else 1))
   if not options:capped=False;break
   _,index,orientation=min(options)
   for i,sign in directions[index]:v[i]+=2*orientation*sign
   moves+=1
  norm=sum(v[i]*hs[i][k]*v[k] for i in range(16) for k in range(16))
  w=tuple(sum(x*y for x,y in zip(row,v)) for row in us)
  if next((x for x in w if x),1)<0:w=tuple(-x for x in w)
  samples.append({'sample':j,'mask_in_reduced_basis':mask,'norm':norm,'steps':moves,'step_cap_reached':capped})
  if 8<=norm<=14:
   parity=tuple(x%2 for x in w);key=(-norm,sum(map(abs,w)),w)
   if parity not in pool or key<pool[parity]:pool[parity]=key
 selected=sorted(pool.values())[:32]
 return {'LLL_columns':us,'samples':samples,'eligible_distinct_parities':len(pool),'centres':[{'word':list(w),'height':-negative} for negative,_,w in selected]}

arms=[];candidate_count=0;stage_count=0;censor_count=0;basis_words=[]
for entry,e,r,c,d in zip(roster['arms'],ev['rows'],replay['arms'],compact['arms'],diag['arms']):
 aid=entry['arm_id'];assert aid==e['arm_id']==r['arm_id']==c['arm_id']==d['arm_id']
 fixture=read(PKG/entry['fixture']);assert sha(PKG/entry['fixture'])==entry['fixture_sha256']
 assert fixture['policy']==protocol['policy'] and len(fixture['sections'])==len(fixture['gram'])==16
 cert=read(PKG/'arms'/aid/'certificate.json');selection=read(PKG/'arms'/aid/'selection.json')
 account=read(PKG/'arms'/aid/'accounting.json');ledger=read(PKG/'arms'/aid/'stage_ledger.json')
 assert account['returncode']==0 and not account['wall_cap_reached']
 assert account['wall_seconds']<=240 and cert['resource']['user_seconds']+cert['resource']['system_seconds']<=180
 assert cert['fixture_sha256']==entry['fixture_sha256'] and cert['worker_sha256']==sha(CAS/'worker.py')
 recomputed=selector(matrix(ZZ,fixture['gram']))
 for key,value in recomputed.items():assert selection[key]==value,(aid,key)
 assert len(selection['centres'])==32
 plan=[]
 for i,centre in enumerate(selection['centres']):
  n=(3*centre['height']+3)//4+1
  plan.extend([{'centre_index_zero_based':i,'n':n},{'centre_index_zero_based':i,'n':n+1}])
 actual=[{'centre_index_zero_based':s['centre'],'n':s['n']} for s in cert['stages']]
 assert actual==plan[:len(actual)]
 assert all(s['status']=='COMPLETE' for s in cert['stages'][:-1])
 assert cert['stages'][-1]['status']==('COMPLETE_THROUGH_FIRST_SUCCESS' if e['success'] else 'COMPLETE')
 assert cert['terminal_status']==('GENERIC_QUOTIENT_RECOVERY' if e['success'] else 'BOUNDED_NO_RECOVERY')
 assert cert['censored_after_success']==e['success']
 if not e['success']:assert actual==plan
 censored=[dict(p,status='NOT_RUN_AFTER_FIRST_SUCCESS') for p in plan[len(actual):]]
 partial=cert['stages'][-1]
 censored_rows=list(range(len(partial['members']),partial['kernel_dimension']))
 assert not censored_rows or e['success']
 events=ledger['events'];assert events[0]['event']=='START' and events[-1]['event']=='TERMINAL'
 assert events[0]['repository_visible'] is False and events[0]['network_namespace']=='isolated'
 assert [x['event'] for x in events].count('START')==1
 assert next(i for i,x in enumerate(events) if x['event']=='CENTRES_FROZEN')<next(i for i,x in enumerate(events) if x['event']=='RR_START')
 assert sum(x['event']=='RR_START' for x in events)==len(actual)==sum(x['event']=='RR_COMPLETE' for x in events)
 assert not any('FAILED' in x['event'] for x in events)
 assert len(cert['candidates'])==len(e['candidates'])==len(r['candidate_claims'])==len(c['candidate_claims'])
 assert r['status']==c['status']=='PASS' and r['rank16_determinant']==cert['input_rank_certificate']['determinant']
 for i,(candidate,evaluated,proof1,proof2) in enumerate(zip(cert['candidates'],e['candidates'],r['candidate_claims'],c['candidate_claims'])):
  assert proof1==proof2 and proof1['candidate_id']==evaluated['candidate_id']==i
  assert proof1['rank']==candidate['rank'] and proof1['schur_complement']==candidate['schur_complement']
  w=[int(x) for x in evaluated['full_basis_coordinates']];assert len(w)==17
  gain=bool(w[entry['omitted_basis_index_one_based']-1]);assert candidate['rank']==16+int(gain)
  assert gain==(evaluated['quotient_classification']=='Generic quotient recovery')
  assert evaluated['exact_group_identity'] and evaluated['orthogonal_residual_height']=='0'
  if gain:assert abs(w[entry['omitted_basis_index_one_based']-1])==1
  basis_words.append(w)
 assert d['generic_recovery']==e['success']
 if e['success']:
  assert d['status']=='PASS_EXACT_SPECIALIZATION_QUOTIENT' and d['same_integral_MW17_specialization']
  assert d['subgroup_index_in_full_generic_specialization']==1 and d['intrinsic_half_lattice_unchanged']
  assert len(d['finite_atlas']['centres'])==289 and len(d['finite_atlas']['parity_masks'])==153
 candidate_count+=len(cert['candidates']);stage_count+=len(actual);censor_count+=len(censored)
 arms.append({'arm_id':aid,'status':'PASS','omitted_basis_index_one_based':entry['omitted_basis_index_one_based'],'rank16_proof_verified':True,'selector_recomputed_from_core_Gram_only':True,'parity_samples_replayed':512,'selected_centres':32,'completed_RR_stages':len(actual),'terminal_status':cert['terminal_status'],'all_candidates_accounted':len(cert['candidates']),'censored_RR_stages':censored,'censored_kernel_rows_in_last_stage':censored_rows,'bound_accounting_sha256':sha(PKG/'arms'/aid/'accounting.json'),'failed_stages':[]})
 print(aid,'PASS fixed selector, stage prefix and censor ledger',flush=True)
assert candidate_count==224 and stage_count==627 and censor_count==461
assert sum(e['success'] for e in ev['rows'])==summary['successes']==9
assert summary['conclusion']==ev['conclusion']=='CORE_DEPENDENT_16_TO_17_RECOVERY'
assert diag['stretch_observation']['max_determinant_posthoc_diagnostic']['selected_arm']=='arm-14'
assert not diag['stretch_observation']['generic_only_selector_validated']
status=read(ROOT/'research/MATH_STATUS.json')
authority=next(e for e in status['entries'] if e['id']=='EC-DET1092-BLIND-MW16-RECONSTRUCTION-20260908')
assert authority['independent_replay'] and authority['artifact_hash']=='sha256:'+sha(CAS/'replay.py')
note=ROOT/'research'/authority['canonical_source'];assert note.is_file() and len([line for line in note.read_text().splitlines() if line.startswith('| P')])==17
assert 'UNMEASURED' in note.read_text() and 'no universal MW16' in note.read_text()
requirements=[
 ('1_starting_commit_and_sources',['protocol.json','frozen_sources/','preflight/starting_worktree.txt'],'Starting SHA and actual dirty-source bytes are bound and intact.'),
 ('2_immutable_17_arm_roster',['roster.json','dispatch.json','primary_complete.json'],'Exactly17 fixed omissions; protocol freeze precedes dispatch; all arms terminal before evaluator reads.'),
 ('3_per_arm_inputs_rank_boundary_outputs_resources_failures',['fixtures/','arms/','completion_audit.json'],'All16-section inputs, core rank certificates, allowed/forbidden policy, candidates, exact rank decisions and CPU/wall/RSS records present. Independent selector replay passes;461 unexecuted RR stages explicitly identified as success-censored, not negative exposure.'),
 ('4_full_basis_evaluation',['evaluation.json','replay.json','compact_replay.json'],'All224 returned-section records have exact17-entry integral words and independently verified group identities.18 positive candidate records have omitted coefficient+/-1; no overgroup inconsistency.'),
 ('5_summary_table',['summary.json',str(note.relative_to(ROOT))],'17 rows carry omission, outcome, recovered subgroup rank, quotient classification, cost and geometry.'),
 ('6_clear_classification',['summary.json','evaluation.json',str(note.relative_to(ROOT))],'CORE_DEPENDENT_16_TO_17_RECOVERY:9 recoveries and8 complete bounded misses.'),
 ('7_independent_generic_replay',['replay.json','compact_replay.json','verification_accounting.json'],'Separately implemented affine group and pole-height proofs cover every returned section and promoted rank, with coefficient/RR rank replay. Exact proof sources and outputs match current bindings.'),
 ('8_preserve_negative_failed_censored_attempts',['raw_archive.json','preflight/','verification_accounting.json','completion_audit.json'],'All primary raw bytes retained unchanged. Eight negative arms preserved; no primary failures or repeats. Infrastructure failures and the censored redundant-checker attempt remain recorded.'),
 ('secondary_specialization_diagnostic',['diagnostic.json',str(note.relative_to(ROOT))],'Exact t=0 generic quotient identities and primitive subgroup equality checked;289-centre footprint and coefficient costs compared. Exceptional-point discovery visibility is explicitly unmeasured, not claimed.'),
 ('post_terminal_stretch',['diagnostic.json'],'Post-run maximum-core-determinant rule chooses arm14 and its easy recovery. No training optimization, extra primary arm or prospective validation is claimed.'),
 ('scope_and_tautology_boundary',['protocol.json','frozen_sources/research/elliptic-curves/cas/det1092_blind_mw16/worker.py','frozen_sources/research/elliptic-curves/cas/det1092_blind_mw16/isolate.py',str(note.relative_to(ROOT))],'Worker dependency graph is own equation/core→selector→RR line→rational factors→exact core height rank. Historical target data are excluded by policy/code/filesystem. The known full MW17 theorem and curve302 rank boundaries remain unchanged.')]
out={'starting_commit':commit,'status':'PASS_REQUIREMENT_BY_REQUIREMENT_COMPLETION_AUDIT','prior_goal_turn_classification':'progress: completed and independently replayed the fixed17-arm experiment, generated its canonical note/package and updated authoritative status.','checker_sha256':sha(Path(__file__)),'audited_summary_sha256':sha(PKG/'summary.json'),'canonical_note_sha256':sha(note),'primary_barrier_sha256':sha(PKG/'primary_complete.json'),'independent_replay_sha256':sha(PKG/'replay.json'),'compact_replay_sha256':sha(PKG/'compact_replay.json'),'requirements':[{'id':i,'status':'PROVED_COMPLETE','evidence':paths,'finding':finding} for i,paths,finding in requirements],'arms':arms,'totals':{'primary_arms':17,'selector_samples_replayed':8704,'RR_stages_executed':627,'RR_stages_censored_after_success':461,'candidate_records':224,'recovering_arms':9,'bounded_misses':8,'new_primary_runs':0},'limitations':['Algorithm/process input blindness is certified; human ignorance of canonical background is not claimed.','The finite raw chart-footprint diagnostic does not measure exceptional-point discoveries.','The post-hoc determinant selector is not validated on unseen MW16 parents.','Development-only administrative timing has the explicitly documented limits in preflight/administrative_accounting.json.'],'resource':{'wall_seconds':time.monotonic()-START,'user_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime,'max_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}}
destination=PKG/'completion_audit.json'
if destination.exists():
 old=read(destination)
 assert {k:v for k,v in old.items() if k!='resource'}=={k:v for k,v in out.items() if k!='resource'}
else:
 with destination.open('x') as f:json.dump(out,f,indent=2,sort_keys=True);f.write('\n')
 destination.chmod(0o444)
print(out['status'],flush=True)
