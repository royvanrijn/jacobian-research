#!/usr/bin/env sage -python
"""Disjoint two-block later-gain control with certified banks in both arms."""
import argparse
from collections import Counter
import gzip
import json
from pathlib import Path
import shutil

import cancellation_scheduler_prepare as preparation
import cancellation_scheduler_round4 as supervisor
from cancellation_basis_amplification import parent_subset
from cancellation_basis_certificate_receipts import CHECKERS
from cancellation_scheduler_fresh import source_closure
from finite_cancellation_corpus import ROOT, LOCAL, OUT as CORPUS, canonical, digest, write

CAS = Path(__file__).resolve().parent
OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_certified_v1'
RAW = LOCAL/'cancellation-basis-certified-v1'
PRIOR = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_early_v1'
INTEGRATION = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_certificate_integration_v1'
COST = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v2'


def read(path): return json.loads(path.read_text())
def sha(path): return digest(path.read_bytes())
def need(condition, message):
    if not condition: raise ArithmeticError(message)
def new_write(path, data): preparation.new_write(path, data)


def prepare():
    from check_cancellation_basis_certificate import guard as integration_guard
    integration_guard()
    integration = read(INTEGRATION/'summary.json')
    need(integration['status'] == 'PASS_COMPLETE_CERTIFICATE_ENGINE_INTEGRATION', 'complete engine integration missing')
    need(read(COST/'completion.json')['engineering_gate'] is True, 'identical-bank cost gate missing')
    prior = read(PRIOR/'design.json')
    keys = ('arms','cases','initial_rank','target_directions','validation_blocks','cases_per_block',
        'bank_policy','maximum_extension_dimension','maximum_centres','prime_bound','bank_wall_seconds',
        'maximum_calls','height','point_wall_seconds','working_cpu_seconds','hard_process_cpu_seconds',
        'arm_wall_seconds','early_order','exploration','subgroup','primary_gate','boundary')
    plan = {k:prior[k] for k in keys}
    plan.update(status='DESIGN_BEFORE_CERTIFIED_BANK_SEARCH_OUTCOMES',
        selection='Exclude every whole-j group in the completed early-exposure roster exclusion set. Keep retained M18 seeds with certified endpoints at least21; lexical first compatible source landscape per j. Start each family allocation at min(4,2*floor(eligible_count/2)); redistribute two slots at a time to the family with most unused eligible groups, ties by fixed family order, until24. All allocations are even. Within each family order by SHA256(cancellation-basis-certified-v1/ + j); first half assigned block0, second half block1. Both twelve-pair blocks fixed before outcomes. No replacement, interim tuning or reuse of prior controls.',
        hypothesis='After equalizing the bank constructor, complete-cloud subgroup enlargement and latest-generator-first compatible anchors yield a repeatable increase in later directions per complete CPU sufficient to repay every rebuild.',
        ablation='Both arms use the same independently certified bank constructor, complete-cloud admission, native rank checks, exploration and final replay. Both reconstruct the same initial M18 bank from16 parent anchors and export at most48 centres. fixed_bank retains its initial search subgroup and bank while retaining all admitted points. basis_refresh reconstructs from every admitted point after each positive cloud below target3, then visits latest-generator anchors first. No missing point, radius fit, residue fit or Schur diagnostic enters selection. The common initial-bank engineering saving cannot be credited as subgroup amplification.',
        cost='The unchanged40 working CPU allowance includes interpreter/imports, initial bank, every rebuild, all exact minimum and selection certificates, native arithmetic, maps, point calls and complete-cloud rank certificates. Atomic operations may overrun the launch allowance. Mandatory final replay is also charged in complete isolated child-process CPU. Both arms have100 CPU/150 wall hard limits. Seed/parent packets are the shared input interface. Shared eligibility preflight, transcript-only development and postmortem processes are reported separately, with conservative candidate-only charge sensitivities. No retrospective timing substitution or retry.',
        mechanism_gate='After both blocks finish, audit all actual later candidate gains using new generators. Independently bind native subgroup/anchor identities and square maps. Prepare old-bank and coefficient-deletion maps without targets. Evaluate both P and -P, then both signs of P + sum k_i G_i for all generators admitted above M18 before the call, each k_i in {-1,0,1}; at most18 representatives. Promotion requires at least one old-dictionary-exclusive, strict coefficient-deletion witness in each block surviving this larger finite dictionary. Numerical rank_growth diagnostics are retrospective only. No full-coset exclusion.',
        postmortem={'mechanism_cpu_limit':300,'representative_cpu_limit':90,'wall_seconds':420,
            'point_search_calls':0,'selection':'Every qualifying actual later gain in both fixed blocks, including target overshoots; never a selected winning subset.'},
        continuation='Exactly48 sequential arms, two fixed12-pair blocks, maximum4800 hard child CPU and7200 nominal worker wall seconds. All started failures and prefixes retained; no failed arm retried. Then one bounded mechanism audit and one bounded finite representative audit. No automatic fresh-fibre search or larger cap follows.',
        prior_protocol_sha256=sha(PRIOR/'protocol.json'), bank_cost_protocol_sha256=sha(COST/'protocol.json'),
        integration_protocol_sha256=sha(INTEGRATION/'protocol.json'),
        bank_certificate_sources={'elliptic-curves/cas/'+n:sha(CAS/n) for n in CHECKERS})
    new_write(OUT/'design.json', plan)
    excluded = set(read(PRIOR/'roster.json')['excluded_training_j_groups'])
    eligible = json.loads(gzip.decompress((preparation.V2/'eligible_banks.json.gz').read_bytes())); banks = {}
    for row in sorted(eligible,key=lambda r:r['source_landscape']):
        if row['j_group'] not in excluded and row['initial_rank']==18 and row['endpoint_rank']>=21:
            banks.setdefault(row['j_group'],row)
    families = preparation.FAMILIES
    available = {f:sum(r['family']==f for r in banks.values()) for f in families}
    allocation = {f:min(4,2*(available[f]//2)) for f in families}
    need(all(allocation.values()), 'all six families must remain represented')
    while sum(allocation.values())<24:
        choices = [f for f in families if available[f]-allocation[f]>=2]
        need(choices,'insufficient disjoint controls; no replacement')
        f = max(choices,key=lambda f:(available[f]-allocation[f],-families.index(f)))
        allocation[f] += 2
    need(sum(allocation.values())==24,'invalid family allocation')
    chosen = []
    for family in families:
        rows = sorted((r for r in banks.values() if r['family']==family),
            key=lambda r:digest(('cancellation-basis-certified-v1/'+r['j_group']).encode()))
        n = allocation[family]
        chosen.extend({**r,'validation_block':i//(n//2)} for i,r in enumerate(rows[:n]))
    chosen.sort(key=lambda r:(r['validation_block'],families.index(r['family'])))
    corpus = {c['id']:c for c in json.loads(gzip.decompress((CORPUS/'corpus.json.gz').read_bytes()))}
    inputs=[]; oracle=[]
    for row in chosen:
        need(sha(ROOT/row['source_seed'])==row['seed_sha256'] and sha(ROOT/row['source_landscape'])==row['landscape_sha256'],'retained input changed')
        seed=read(ROOT/row['source_seed']); landscape=read(ROOT/row['source_landscape'])
        need(seed['points']==landscape['basis'] and landscape['generic_rank']==17,'incompatible seed')
        ident=digest(canonical(['basis-certified-v1',row['j_group'],row['source_landscape']]))[:20]; row['id']=ident
        inputs.append({'id':ident,'j_group':row['j_group'],'family':row['family'],'stratum':'M18',
            'validation_block':row['validation_block'],'seed':{k:seed[k] for k in ('curve','points','proof')},
            'parent_bank':parent_subset(landscape)})
        endpoint=corpus[row['corpus_id']]
        oracle.append({'id':ident,'curve':endpoint['curve'],'points':endpoint['generic_points']+endpoint['targets'],
            'source':endpoint['source'],'proof_sha256':endpoint['proof_sha256']})
    new_write(OUT/'inputs.json',inputs); new_write(OUT/'oracle.json',oracle)
    new_write(OUT/'roster.json',{'chosen':chosen,'eligible_whole_j_counts':available,'family_allocation':allocation,
        'prior_cpu_j_groups':sorted(excluded),'excluded_training_j_groups':sorted(excluded|{c['j_group'] for c in inputs}),
        'previous_roster_sha256':sha(PRIOR/'roster.json'),'eligible_banks_sha256':sha(preparation.V2/'eligible_banks.json.gz'),
        'inputs_sha256':sha(OUT/'inputs.json'),'oracle_sha256':sha(OUT/'oracle.json')})
    new_write(OUT/'integration-prerequisite.json',{'summary':integration,'summary_sha256':sha(INTEGRATION/'summary.json'),
        'negative_checks':read(INTEGRATION/'negative-checks.json'), 'cost_completion':read(COST/'completion.json'),
        'boundary':'Already completed known-transcript development and identical-bank engineering prerequisites; zero new point searches. Neither supplies a validation outcome.'})
    print(json.dumps({'status':'FROZEN_CERTIFIED_BANK_ROSTER','cases':len(inputs),'allocation':allocation,
        'blocks':Counter(c['validation_block'] for c in inputs)}),flush=True)


def preflight():
    import resource
    resource.setrlimit(resource.RLIMIT_CPU,(120,125))
    preparation.OUT=OUT; preparation.preflight()
    need(all(r['seed']['rank']==18 and r['endpoint']['rank']>=21 for r in read(OUT/'preflight.json')['rows']),'three-direction eligibility missing')


def seal():
    from sage.all import pari
    from sage.env import SAGE_VERSION
    from pari_pointed_backend import sources
    import sys
    plan=read(OUT/'design.json'); inputs=read(OUT/'inputs.json'); pre=read(OUT/'preflight.json')
    need(pre['status']=='PASS' and len(pre['rows'])==24,'input preflight missing')
    need(len({c['j_group'] for c in inputs})==24 and Counter(c['validation_block'] for c in inputs)=={0:12,1:12},'invalid disjoint blocks')
    need(not {c['j_group'] for c in inputs}&set(read(OUT/'roster.json')['prior_cpu_j_groups']),'prior control reuse')
    new_write(OUT/'environment.json',{'sage':SAGE_VERSION,'pari':str(pari('version()')),'python':sys.version,
        'threads':{'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}})
    names=[Path(__file__),CAS/'cancellation_basis_certificate_epoch.py',CAS/'verify_cancellation_basis_certificate.py',
        CAS/'audit_cancellation_basis_certified.py',CAS/'audit_cancellation_basis_certified_representatives.py',
        CAS/'retain_cancellation_basis_certified.py']+[ROOT/n for n in sources()]
    names+=list((ROOT/'elliptic-curves/ecsearch').glob('*.py'))
    plan.update(status='FROZEN_TWO_BLOCK_CERTIFIED_BANK_CONTROL',source_sha256=source_closure(names),
        gp_sha256=sha(Path('/usr/bin/gp')),
        input_sha256={n:sha(OUT/n) for n in ('design.json','inputs.json','oracle.json','roster.json','preflight.json','integration-prerequisite.json','endpoint-packets.json','environment.json')})
    new_write(OUT/'protocol.json',plan)
    for name in plan['source_sha256']:
        p=RAW/'sources'/name; p.parent.mkdir(parents=True,exist_ok=True); shutil.copyfile(ROOT/name,p)
    print(json.dumps({'status':plan['status'],'protocol_sha256':sha(OUT/'protocol.json'),'sources':len(plan['source_sha256'])}),flush=True)


def guard():
    plan=read(OUT/'protocol.json')
    for n,h in plan['source_sha256'].items():need(sha(ROOT/n)==h,'sealed source changed: '+n)
    for n,h in plan['input_sha256'].items():need(sha(OUT/n)==h,'sealed input changed: '+n)
    need(sha(Path('/usr/bin/gp'))==plan['gp_sha256'],'GP executable changed')
    return plan,read(OUT/'inputs.json')


def worker(case_id,arm):
    from cancellation_basis_certificate_epoch import run
    plan,inputs=guard(); case=next(c for c in inputs if c['id']==case_id)
    run(case,arm,plan,RAW/'arms'/case_id/arm,sha(OUT/'protocol.json'))


def run():
    guard(); supervisor.OUT=OUT; supervisor.RAW=RAW; supervisor.ENTRY=Path(__file__).resolve(); supervisor.run()


def report():
    from cancellation_cloud_programme import compare
    from retain_cancellation_basis_certified import bindings
    audit=bindings(); plan,inputs=guard(); cases={c['id']:c for c in inputs}; rows=[]
    for receipt in read(OUT/'supervision.json')['records']:
        r=read(RAW/'arms'/receipt['case']/receipt['arm']/'result.json')
        rows.append({**{k:r[k] for k in ('case','family','stratum','arm','success','new_directions','later_cloud_directions','first_cloud_gain','unknowns','status')},
            'validation_block':cases[r['case']]['validation_block'],'calls':len(r['calls']),'epochs':len(r['epochs']),
            'cpu_seconds':receipt['charged_cpu_seconds'],'components':r['components']})
    def comparison(selected,later):
        renamed=[{**r,'arm':'factor_free' if r['arm']=='fixed_bank' else 'adaptive_cloud',
            'new_directions':r['later_cloud_directions'] if later else r['new_directions']} for r in selected]
        totals,outcome=compare(renamed)
        return {'totals':{('fixed_bank' if k=='factor_free' else 'basis_refresh'):v for k,v in totals.items()},'comparison':outcome}
    all_gain=comparison(rows,False); later=comparison(rows,True)
    blocks={str(b):{'all':comparison([r for r in rows if r['validation_block']==b],False),
                   'later':comparison([r for r in rows if r['validation_block']==b],True)} for b in (0,1)}
    unknowns=sum(r['unknowns'] for r in rows); a,b=all_gain['totals']['basis_refresh'],all_gain['totals']['fixed_bank']
    performance=bool(later['comparison']['gate_passed'] and a['directions']>=b['directions'] and a['target_completions']>=b['target_completions'] and not unknowns)
    repeat=all(x['later']['comparison']['uncapped_direction_rate_ratio'] is not None and x['later']['comparison']['uncapped_direction_rate_ratio']>1 for x in blocks.values())
    mechanism=read(OUT/'mechanism.json'); representatives=read(OUT/'representative-audit.json')
    causal=bool(mechanism['status']=='PASS' and representatives['status']=='PASS_FINITE_REPRESENTATIVE_SENSITIVITY' and
        all(representatives['blocks'][str(i)]['both']>0 for i in (0,1)))
    shared=read(OUT/'preflight.json')['cpu_seconds']; diagnostic=sum(read(OUT/n)['full_child_cpu_seconds'] for n in ('mechanism-process.json','representative-process.json'))
    development=read(OUT/'integration-prerequisite.json')['summary']
    dev_cpu=development['full_child_cpu_seconds']+development['negative_check_component_cpu_seconds']
    la,lb=later['totals']['basis_refresh'],later['totals']['fixed_bank']
    def ratio(extra):return (la['directions']/(la['cpu_seconds']+extra))/(lb['directions']/lb['cpu_seconds']) if lb['directions'] else None
    summary={'status':'COMPLETE_FIXED_TWO_BLOCK_CERTIFIED_BANK_CONTROL','rows':rows,'all':all_gain,'later':later,'blocks':blocks,
        'aggregate_performance_gate':performance,'repeatability_rate_gate':repeat,'actual_mechanism_gate':causal,
        'promotion_gate':performance and repeat and causal,'preparation_unknowns':unknowns,
        'shared_preflight_component_cpu_seconds':shared,'postmortem_full_child_cpu_seconds':diagnostic,
        'transcript_development_recorded_cpu_seconds':dev_cpu,
        'conservative_later_rate_ratio':ratio(shared),'conservative_ratio_including_postmortem':ratio(shared+diagnostic),
        'conservative_ratio_including_recorded_integration':ratio(shared+diagnostic+dev_cpu),
        'exposure':{k:audit[k] for k in ('totals','paired_initial_exposure','point_call_statuses')},
        'protocol_sha256':sha(OUT/'protocol.json'),'fresh_fibres_run':0,'rank32':'UNKNOWN','boundary':plan['boundary']}
    write(OUT/'summary.json',summary)
    print(json.dumps({k:summary[k] for k in ('all','later','aggregate_performance_gate','repeatability_rate_gate','actual_mechanism_gate','promotion_gate','preparation_unknowns')},indent=2),flush=True)


def pack():
    import pack_cancellation_scheduler as p
    p.OUT=OUT; p.RAW=RAW; p.pack()


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('command',choices=['prepare','preflight','seal','worker','run','report','pack']); p.add_argument('--case'); p.add_argument('--arm'); a=p.parse_args()
    worker(a.case,a.arm) if a.command=='worker' else globals()[a.command]()
