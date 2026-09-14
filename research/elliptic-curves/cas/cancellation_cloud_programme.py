#!/usr/bin/env sage -python
"""Independent two-direction control stage after the first-recovery programme."""
import argparse
from collections import defaultdict
import gzip
import json
from pathlib import Path
import time

import cancellation_scheduler_prepare as preparation
import cancellation_scheduler_round4 as supervisor
from cancellation_scheduler_fresh import source_closure
from cancellation_cloud_training import OUT,RAW,CAS
from finite_cancellation_corpus import ROOT,OUT as CORPUS,canonical,digest,write

PRIOR=ROOT/'artifacts/generated-results/elliptic-curves/cancellation_scheduler_v4'


def read(path):return json.loads(path.read_text())
def sha(path):return digest(path.read_bytes())
def need(value,message):
    if not value:raise ArithmeticError(message)


def prepare():
    need(read(PRIOR/'summary.json')['high_rank_transfer_gate'] is True,'first-direction policy not validated')
    # The stopping objective and release gates are fixed before new outcomes.
    design={'status':'DESIGN_BEFORE_TWO_DIRECTION_OUTCOMES','arms':['factor_free','adaptive_cloud'],
        'cases':24,'target_directions':2,'maximum_centres':48,'maximum_calls':432,
        'heights':[8000,32000,125000],'search_cpu_seconds':40,'point_wall_seconds':5,
        'arm_wall_seconds':120,'hard_process_cpu_seconds':100,
        'selection':'24 new whole-j CPU controls, four from each retained family. Exclude all preceding CPU groups. Require a certified input of rank19 or20 and retained endpoint at least two ranks larger. Within each j choose highest starting rank then lexical bank; within each family prefer rank20 then SHA256(cancellation-cloud-v1/ + j). Known endpoint eligibility is only masked-control construction, never prospective selection.',
        'fit':'Refit the identical q-only residue mixture excluding all new and previous CPU groups. Reconcile all545 old V2 factor-free development calls at the same fixed finite places through500, independently certify their full clouds and use all selected independent witness radii. Gamma prior strength4; smoothed bucket rate=(directions+1)/(calls+10). No new control outcome enters a fit.',
        'policy':'V4 exact box geometry with expected capped Poisson direction yield per predicted complete cost. Completed boxes remove overlap and update per-anchor Gamma intensity by their new independently admissible column count. Timeouts change neither intensity nor covered exposure. All points in every returned cloud are considered before stopping at two added directions.',
        'subgroup':'Every point call is defined on the original certified subgroup and fixed compatible bank. Accumulated output columns do not become a new search basis within this arm. A later enlarged-basis epoch requires a newly verified landscape and a separate cost allowance.',
        'primary_gate':'Adaptive must achieve at least as many two-direction completions AND at least as many uncapped total independent directions as V3. Its uncapped directions per complete CPU must be at least1.1 times V3 with central97.5-percent paired family/stratum bootstrap interval wholly above1;10000 draws seed20260914. Extra cloud directions count; they are not discarded to improve the comparison. Zero reference yield makes the gate UNKNOWN.',
        'cost':'Complete isolated process CPU includes setup, every model, scoring, full-cloud admission, failed work, exact transcript/policy replay and independent Sage rank certification. One-time fit and shared input checks are reported separately and conservatively charged to candidate in campaign sensitivity. Retained compatible banks are common inputs; this does not measure cold landscape creation.',
        'continuation_gate':'Only a passing new control stage permits a separately frozen generic-only fresh-fibre two-direction test with complete cold costs. No automatic timing retries, expanded case allowance or restart of an old campaign.',
        'boundary':'This is a finite multi-direction exposure test, not a theorem about Poisson rational-point supply, a new basis landscape, population speedup or rank32.'}
    preparation.new_write(OUT/'design.json',design)
    prior=read(PRIOR/'roster.json');excluded=set(prior['excluded_training_j_groups'])
    eligible=json.loads(gzip.decompress((preparation.V2/'eligible_banks.json.gz').read_bytes()))
    by_j={}
    for row in sorted(eligible,key=lambda r:(-r['initial_rank'],r['source_landscape'])):
        if row['j_group'] not in excluded and row['initial_rank'] in (19,20) and row['endpoint_rank']>=row['initial_rank']+2:
            by_j.setdefault(row['j_group'],row)
    chosen=[]
    for family in preparation.FAMILIES:
        candidates=sorted((r for r in by_j.values() if r['family']==family),key=lambda r:(-r['initial_rank'],digest(('cancellation-cloud-v1/'+r['j_group']).encode())))
        need(len(candidates)>=4,'insufficient control family; no replacement');chosen.extend(dict(r) for r in candidates[:4])
    corpus={r['id']:r for r in json.loads(gzip.decompress((CORPUS/'corpus.json.gz').read_bytes()))};inputs=[];oracle=[]
    for row in chosen:
        need(sha(ROOT/row['source_seed'])==row['seed_sha256'] and sha(ROOT/row['source_landscape'])==row['landscape_sha256'],'retained bank changed')
        seed=read(ROOT/row['source_seed']);landscape=read(ROOT/row['source_landscape']);need(seed['points']==landscape['basis'],'bank basis differs')
        ident=digest(canonical(['cloud-v1',row['j_group'],row['source_landscape']]))[:20];row['id']=ident
        inputs.append({'id':ident,'j_group':row['j_group'],'family':row['family'],'stratum':'M'+str(row['initial_rank']),
            'seed':{k:seed[k] for k in ('curve','points','proof')},'centres':[c['representative'] for c in landscape['centres'][:48]]})
        endpoint=corpus[row['corpus_id']]
        oracle.append({'id':ident,'curve':endpoint['curve'],'points':endpoint['generic_points']+endpoint['targets'],
            'source':endpoint['source'],'proof_sha256':endpoint['proof_sha256']})
    preparation.new_write(OUT/'inputs.json',inputs);preparation.new_write(OUT/'oracle.json',oracle)
    preparation.new_write(OUT/'roster.json',{'chosen':chosen,'prior_cpu_j_groups':sorted(excluded),
        'excluded_training_j_groups':sorted(excluded|{r['j_group'] for r in inputs}),
        'previous_roster_sha256':sha(PRIOR/'roster.json'),'inputs_sha256':sha(OUT/'inputs.json'),'oracle_sha256':sha(OUT/'oracle.json')})
    print(json.dumps({'status':'FROZEN_TWO_DIRECTION_ROSTER','cases':len(inputs),'ranks':{r:sum(len(c['seed']['points'])==r for c in inputs) for r in (19,20)}}),flush=True)


def train():preparation.OUT=OUT;preparation.train()
def preflight():
    preparation.OUT=OUT;preparation.preflight()
    need(all(r['endpoint']['rank']>=r['seed']['rank']+2 for r in read(OUT/'preflight.json')['rows']),'two-direction input not certified')
def training_replay():
    import verify_cancellation_scheduler_training as verifier
    verifier.OUT=OUT;verifier.main()
def cloud_train():
    import cancellation_cloud_training as module
    module.train()
def cloud_replay():
    import cancellation_cloud_training as module
    module.verify()


def smoke():
    """Software control on a known two-point DEVELOPMENT cloud, not a holdout."""
    import cancellation_cloud_cpu as module
    training=read(OUT/'cloud-fit.json');row=next(r for r in training['records'] if r['directions']>=2)
    case=next(c for c in read(preparation.V2/'inputs.json') if c['id']==row['case'])
    case={**case,'id':'DEVELOPMENT-'+case['id'],'centres':[read(ROOT/row['source'])['centre']]}
    plan={**read(OUT/'design.json'),'fit':training['fit'],'search_cpu_seconds':10,
        'gp_sha256':read(PRIOR/'protocol.json')['gp_sha256'],
        'source_sha256':source_closure([Path(__file__),CAS/'cancellation_cloud_cpu.py',CAS/'verify_cancellation_cloud.py'])}
    dest=RAW/'development';preparation.new_write(dest/'recipe.json',{'case':case,'plan':plan,'development_source':row['source'],
        'boundary':'Known winning development anchor validates multiple-point software semantics only. Never included in the new control comparison.'})
    results=[]
    for arm in plan['arms']:
        target=dest/arm;module.run(case['id'],arm,development_context=(plan,case,target))
        result=read(target/'result.json');need(result['success'] and result['new_directions']>=2,'multi-point development control failed')
        results.append({'arm':arm,'directions':result['new_directions'],'result_sha256':sha(target/'result.json')})
    preparation.new_write(OUT/'development-smoke.json',{'status':'PASS_DEVELOPMENT_ONLY','rows':results,'recipe_sha256':sha(dest/'recipe.json')})


def seal():
    plan=read(OUT/'design.json');base=read(OUT/'fit.json');cloud=read(OUT/'cloud-fit.json')
    need(read(OUT/'training-replay.json')['fit_sha256']==sha(OUT/'fit.json'),'base fit replay missing')
    need(read(OUT/'cloud-training-replay.json')['cloud_fit_sha256']==sha(OUT/'cloud-fit.json'),'cloud fit replay missing')
    need(read(OUT/'preflight.json')['status']=='PASS','preflight missing')
    need(read(OUT/'development-smoke.json')['status']=='PASS_DEVELOPMENT_ONLY','software smoke missing')
    old=read(PRIOR/'protocol.json')
    sources=source_closure([ROOT/p for p in old['source_sha256']]+[Path(__file__),CAS/'cancellation_cloud_cpu.py',CAS/'verify_cancellation_cloud.py',CAS.parent/'tests/test_cancellation_cloud_policy.py'])
    for name,h in old['source_sha256'].items():need(sources[name]==h,'ancestor implementation changed')
    test={c['j_group'] for c in read(OUT/'inputs.json')}
    need(not test&set(base['training_j_groups']) and not test&set(cloud['development_j_groups']),'fit leakage')
    need(cloud['base_fit_sha256']==sha(OUT/'fit.json'),'cloud fit refers to another base')
    plan.update(status='FROZEN_TWO_DIRECTION_COMPARISON',source_sha256=sources,gp_sha256=old['gp_sha256'],fit=cloud['fit'],
        input_sha256={n:sha(OUT/n) for n in ('design.json','inputs.json','oracle.json','roster.json','fit.json','preflight.json','training-replay.json','cloud-fit.json','cloud-training-replay.json','development-smoke.json')})
    preparation.new_write(OUT/'protocol.json',plan);print('SEALED_TWO_DIRECTION_COMPARISON',sha(OUT/'protocol.json'),flush=True)


def worker(case,arm):
    import cancellation_cloud_cpu as module
    module.OUT=OUT;module.RAW=RAW;module.run(case,arm)
def run():
    supervisor.OUT=OUT;supervisor.RAW=RAW;supervisor.ENTRY=Path(__file__).resolve();supervisor.run()


def compare(rows):
    import numpy as np
    a,b='adaptive_cloud','factor_free';table={(r['case'],r['arm']):r for r in rows};groups=defaultdict(list)
    for row in rows:
        if row['arm']==b:groups[row['family'],row['stratum']].append(row['case'])
    groups=[sorted(ids) for _,ids in sorted(groups.items())]
    totals={arm:{'cases':sum(r['arm']==arm for r in rows),'target_completions':sum(r['success'] for r in rows if r['arm']==arm),
        'directions':sum(r['new_directions'] for r in rows if r['arm']==arm),'cpu_seconds':sum(r['cpu_seconds'] for r in rows if r['arm']==arm),
        'calls':sum(r['calls'] for r in rows if r['arm']==arm)} for arm in (a,b)}
    A,B=totals[a],totals[b];ratio=(A['directions']/A['cpu_seconds'])/(B['directions']/B['cpu_seconds']) if B['directions'] else None
    rng=np.random.default_rng(20260914);ratios=[];undefined=0
    for _ in range(10000):
        ids=[x for group in groups for x in rng.choice(group,len(group),replace=True)]
        ca=sum(table[x,a]['cpu_seconds'] for x in ids);cb=sum(table[x,b]['cpu_seconds'] for x in ids)
        ga=sum(table[x,a]['new_directions'] for x in ids);gb=sum(table[x,b]['new_directions'] for x in ids)
        if gb:ratios.append((ga/ca)/(gb/cb))
        else:undefined+=1
    interval=list(map(float,np.quantile(ratios,[.0125,.9875]))) if ratios and not undefined else None
    gate=bool(A['target_completions']>=B['target_completions'] and A['directions']>=B['directions'] and ratio>=1.1 and interval[0]>1) if ratio is not None and interval else None
    return totals,{'uncapped_direction_rate_ratio':ratio,'paired_stratified_ci975':interval,'undefined_samples':undefined,'gate_passed':gate,
        'candidate_only_target':sum(table[x,a]['success'] and not table[x,b]['success'] for g in groups for x in g),
        'reference_only_target':sum(table[x,b]['success'] and not table[x,a]['success'] for g in groups for x in g)}


def report():
    plan=read(OUT/'protocol.json');supervision=read(OUT/'supervision.json');rows=[];unknowns=0;statuses=defaultdict(int)
    need(supervision['status']=='COMPLETE' and len(supervision['records'])==48,'incomplete control stage')
    for receipt in supervision['records']:
        dest=RAW/'arms'/receipt['case']/receipt['arm'];result=read(dest/'result.json');verification=read(dest/'independent-verification.json')
        need(receipt['status']=='COMPLETE' and receipt['result_sha256']==sha(dest/'result.json'),'process receipt differs')
        need(result['independent_verification_sha256']==sha(dest/'independent-verification.json') and verification['status']=='PASS','independent proof differs')
        need(result['rank_lower_bound']==verification['rank']['rank']==result['initial_rank']+result['new_directions'],'rank count differs')
        unknowns+=verification['preparation_unknowns']
        for call in result['calls']:statuses[call['status']]+=1
        rows.append({'case':receipt['case'],'family':receipt['family'],'stratum':receipt['stratum'],'arm':receipt['arm'],
            'success':result['success'],'new_directions':result['new_directions'],'initial_rank':result['initial_rank'],
            'rank_lower_bound':result['rank_lower_bound'],'calls':len(result['calls']),'cpu_seconds':receipt['charged_cpu_seconds'],
            'status':result['status']})
    totals,comparison=compare(rows)
    offline={n:read(OUT/n)['cpu_seconds'] for n in ('fit.json','preflight.json','training-replay.json','cloud-fit.json','cloud-training-replay.json')}
    a,b=totals['adaptive_cloud'],totals['factor_free']
    conservative=(a['directions']/(a['cpu_seconds']+sum(offline.values())))/(b['directions']/b['cpu_seconds']) if b['directions'] else None
    write(OUT/'summary.json',{'status':'COMPLETE_TWO_DIRECTION_COMPARISON','rows':rows,'totals':totals,'comparison':comparison,
        'point_call_statuses':dict(statuses),'preparation_unknowns':unknowns,'offline_cpu_seconds':offline,
        'conservative_campaign_rate_ratio':conservative,'fresh_fibre_gate':comparison['gate_passed'] is True and not unknowns,
        'protocol_sha256':sha(OUT/'protocol.json'),'rank32':'UNKNOWN',
        'boundary':'24 new whole-j CPU controls with two known withheld dimensions; no external-pristine corpus claim. Counts include every certified cloud gain and failed-target partial gain. All searches use the fixed original subgroup and bank. Exact ranks, new record, fresh transfer and improved later amplification are not asserted.'})
    print(json.dumps({'totals':totals,'comparison':comparison,'offline_cpu_seconds':offline,'conservative_campaign_rate_ratio':conservative,'point_call_statuses':dict(statuses)},indent=2),flush=True)


def pack():
    import pack_cancellation_scheduler as module
    module.OUT=OUT;module.RAW=RAW;module.pack()
def replay():
    import verify_cancellation_cloud as module
    module.OUT=OUT;module.RAW=RAW;module.main()


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['prepare','train','preflight','training_replay','cloud_train','cloud_replay','smoke','seal','worker','run','report','pack','replay'])
    p.add_argument('--case');p.add_argument('--arm');a=p.parse_args();worker(a.case,a.arm) if a.command=='worker' else globals()[a.command]()
