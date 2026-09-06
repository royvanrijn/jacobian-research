#!/usr/bin/env python3
"""Post-terminal exact equation novelty; point ranks remain separately certified."""
import argparse,sys,time
from pathlib import Path
import certify_compact_r17_candidates as cert
import strata60_mw16_pari_batch as batch
from research_runtime.store import checkpoint
from research_runtime.supervisor import run,Limits
ROOT=batch.ROOT;ART=batch.ART;D=batch.extension.D/'novelty-controller'
UPSTREAM=batch.extension.D/'accounting-controller'
OUT=ART/'strata60_mw16_novelty_v1.json'
PREVIOUS=ART/'corrected60_mw16_results_v1.json'
REPORT=ART/'retained_mw16_score_strata_experiment_v1.json'
AUDIT=ART/'retained_mw16_score_strata_accounting_replay_v1.json'


def jkey(model):
    v=cert.weierstrass_invariants(tuple(map(cert.F,model)))
    if not v['discriminant']:raise ArithmeticError('singular equation')
    return v['c4']**3/v['discriminant']


def index(rows,key):
    result={}
    for row in rows:result.setdefault(jkey(row[key]),[]).append(row)
    return result


def expected():
    if cert.read(UPSTREAM/'ledger.json')['status']!='PASS' or cert.read(batch.extension.D/'point-controller/ledger.json')['status']!='COMPLETE_FIXED_COMPARISON_AND_ACCOUNTING':
        raise ArithmeticError('terminal point/proof/accounting gates required before novelty')
    p=batch.protocol();report=cert.read(REPORT);audit=cert.read(AUDIT)
    if audit['status']!='PASS':raise ArithmeticError('independent accounting required')
    for source in (report,audit):
        if any(cert.hashed(ROOT/n)!=h for n,h in source['sources'].items()):raise ArithmeticError('terminal source binding changed')
    previous=cert.read(PREVIOUS)
    old=previous['previous_equations']+[{'address':PREVIOUS.name+':'+r['family']+':'+r['parameter'],'curve':r['curve']} for r in previous['curves']]
    catalogue=previous['catalogue']['equations']
    if len(old)!=1345 or len(catalogue)!=593:raise ArithmeticError('fixed1345 prior and593 catalogue equations required')
    ci=index(catalogue,'ainvs');oi=index(old,'curve');rows=[];pairs=[]
    if [r['id'] for r in report['rows']]!=[r['id'] for r in p['rows']] or len(p['rows'])!=60:raise ArithmeticError('all60 allocations required')
    sources={str(path.relative_to(ROOT)):cert.hashed(path) for path in (Path(__file__).resolve(),Path(cert.__file__),PREVIOUS,REPORT,AUDIT,batch.D/'protocol.json',batch.extension.OUT)}
    for i,(chosen,outcome) in enumerate(zip(p['rows'],report['rows'])):
        model=chosen['model'];j=jkey(model);point_binding=None
        if outcome['rank_lower_bound'] is not None:
            folder=batch.D/chosen['id'];v=cert.read(folder/'verification.json')
            if v['status']!='PASS' or v['rank_lower_bound']!=outcome['rank_lower_bound']:raise ArithmeticError('verified point rank differs')
            cloud=cert.read(ROOT/v['cloud_path']);odd=cert.read(ROOT/v['odd_path'])
            if cert.hashed(ROOT/v['cloud_path'])!=v['cloud_sha256'] or cert.hashed(ROOT/v['odd_path'])!=v['odd_sha256'] or not cert.isomorphic(model,cloud['curve']) or cloud['curve']!=odd['curve']:
                raise ArithmeticError('selected equation/point certificate differs')
            lower=max(cloud['rank_lower_bound'],*(a['finite_column_rank'] for a in odd['audits']))
            if lower!=outcome['rank_lower_bound']:raise ArithmeticError('point lower bound differs')
            point_binding={k:v[k] for k in ('cloud_path','cloud_sha256','odd_path','odd_sha256')}
            for path in (folder/'verification.json',ROOT/v['cloud_path'],ROOT/v['odd_path']):sources[str(path.relative_to(ROOT))]=cert.hashed(path)
        row={k:chosen[k] for k in ('id','family','parameter','arm','band','sign','block')}
        row.update(curve=model,j_invariant=str(j),rank_lower_bound=outcome['rank_lower_bound'],
                   point_proof_binding=point_binding,
                   catalogue_matches=[q['id'] for q in ci.get(j,[]) if cert.isomorphic(model,q['ainvs'])],
                   previous_matches=[q['address'] for q in oi.get(j,[]) if cert.isomorphic(model,q['curve'])])
        for other in rows:
            if other['j_invariant']==str(j) and cert.isomorphic(model,other['curve']):pairs.append([other['id'],row['id']])
        rows.append(row)
    high=[r['id'] for r in rows if r['rank_lower_bound'] is not None and r['rank_lower_bound']>=22 and not r['catalogue_matches'] and not r['previous_matches']]
    return {'schema':'elliptic-curves.strata60-novelty.v1','status':'PASS','sources':sources,
        'rows':rows,'within_cohort_isomorphic_pairs':pairs,'unmatched_high_rank_ids':high,
        'catalogue':previous['catalogue'],'previous_equations':old,
        'claim_boundary':'Exact rational-isomorphism comparisons for all60 frozen selected equations against593 pinned catalogue equations and1345 previously measured equations, after terminal proofs and independent computation accounting. Certified lower bounds are bound to separate point proofs and exactly isomorphic search equations; this check does not rerun point independence. Unresolved ranks remain null. A catalogue absence is not universal novelty; no exact rank, minimal model, conductor, new search or production-policy change is asserted.'}


def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve terminal novelty protocol')
    paths=[Path(__file__).resolve(),Path(cert.__file__),PREVIOUS,batch.extension.OUT,UPSTREAM/'protocol.json']
    checkpoint(D/'protocol.json',{'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
        'wait_seconds':172800,'seconds_per_stage':300,'rss_bytes':2147483648,
        'scope':'After the fixed comparison and independent accounting finish, exactly compare all60 prospective equations with593 pinned catalogue and1345 prior equations. Report certified>=22 unmatched rows separately; leave unknown ranks unknown. Build and read-only replay, no new point or parameter search.'})


def launch():
    p=cert.read(D/'protocol.json')
    if (D/'ledger.json').exists():raise FileExistsError('preserve novelty launch')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('novelty source changed')
    ledger={'status':'WAITING_FOR_TERMINAL_ACCOUNTING','rows':[]};checkpoint(D/'ledger.json',ledger);deadline=time.monotonic()+p['wait_seconds']
    try:
        while True:
            state=cert.read(UPSTREAM/'ledger.json')['status']
            if state=='PASS':break
            if state=='NOT_RUN_MATCHING_INCOMPLETE':
                ledger['status']='NOT_RUN_MATCHING_INCOMPLETE';checkpoint(D/'ledger.json',ledger);return
            if state not in ('WAITING_FOR_TERMINAL_COMPARISON','RUNNING') or time.monotonic()>deadline:raise ArithmeticError('upstream failed or fixed wait expired')
            time.sleep(5)
        if cert.read(batch.extension.D/'point-controller/ledger.json')['status']!='COMPLETE_FIXED_COMPARISON_AND_ACCOUNTING':
            ledger['status']='NOT_RUN_PRESEARCH_GATE_FAILURE';checkpoint(D/'ledger.json',ledger);return
        if any(cert.hashed(ROOT/n)!=h for n,h in p['sources'].items()):raise ArithmeticError('novelty source changed during wait')
        ledger['status']='RUNNING';checkpoint(D/'ledger.json',ledger)
        for stage in ('build','check'):
            s=run([sys.executable,str(Path(__file__).resolve()),stage],limits=Limits(p['seconds_per_stage'],p['rss_bytes']),cwd=ROOT,log_path=D/(stage+'.log'),checkpoint_path=D/(stage+'.supervisor.json'))
            ok=s['outcome']=='completed' and s['returncode']==0
            ledger['rows'].append({'name':stage,'status':'PASS' if ok else 'FAILED_OR_CENSORED','supervision':s});checkpoint(D/'ledger.json',ledger)
            if not ok:raise ArithmeticError('novelty check failed; no retry')
        ledger['status']='PASS';checkpoint(D/'ledger.json',ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'ledger.json',ledger);raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','launch','build','check']);a=p.parse_args()
    if a.stage in ('prepare','launch'):globals()[a.stage]()
    else:
        d=expected()
        if a.stage=='check':
            if cert.read(OUT)!=d:raise ArithmeticError('exact novelty comparison replay differs')
        else:
            if OUT.exists():raise FileExistsError('preserve novelty certificate')
            checkpoint(OUT,d)
        print('STRATA60 NOVELTY',d['status'],'UNMATCHED>=22',d['unmatched_high_rank_ids'],flush=True)
