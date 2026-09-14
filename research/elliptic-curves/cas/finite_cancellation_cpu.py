#!/usr/bin/env sage -python
"""Prepare or execute a sealed next-direction model-policy CPU comparison.

The worker reads only its known subgroup, retained V3 centre bank and fitted
coefficients. No final point packet, target height, or discovery chart is read.
"""
import argparse
from copy import deepcopy
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import resource
import signal
import time

from finite_cancellation_corpus import ROOT, OUT, canonical, digest, write, short

CPU=OUT/'cpu'
FAMILIES=['074d9','07ca9','08234','08f72','103b2','11952']

def prepare():
    from analyze_finite_cancellation import FINITE, ridge
    import numpy as np
    start=time.process_time();corpus=json.loads(gzip.decompress((OUT/'corpus.json.gz').read_bytes()))
    chosen=[];bindings={}
    for family in FAMILIES:
        candidates=sorted((r for r in corpus if r['family']==family and 'broad-cases/' in r['source']),key=lambda r:r['j_group'])
        for c in candidates:
            root=ROOT/c['source'].split('/batch-')[0]
            folder=root/'batch-000/search-00';sp=folder/'seed.json';lp=folder/'epoch-00/landscape/selection.json'
            if not sp.exists() or not lp.exists():continue
            seed=json.loads(sp.read_text());selection=json.loads(lp.read_text())
            if seed['rank_lower_bound']>=c['rank_lower_bound']:continue
            if selection['basis']!=seed['points']:continue
            cm,cp=short(seed['curve'],seed['points'])
            assert cm==c['curve'] and cp[:len(c['generic_points'])]==c['generic_points']
            centres=[r['representative'] for r in selection['centres'][:48]]
            if len(centres)<16:continue
            for p in [sp,lp]:bindings[str(p.relative_to(ROOT))]=digest(p.read_bytes())
            chosen.append({'id':c['id'],'j_group':c['j_group'],'family':family,'seed':seed,'centres':centres,
                'source_seed':str(sp.relative_to(ROOT)),'source_landscape':str(lp.relative_to(ROOT)),
                'known_larger_subgroup_rank':c['rank_lower_bound']})
            if sum(x['family']==family for x in chosen)==2:break
    assert len(chosen)==12 and len({x['j_group'] for x in chosen})==12
    excluded={x['j_group'] for x in chosen};train=[]
    for path in sorted((OUT/'cases').glob('*.json.gz')):
        packet=json.loads(gzip.decompress(path.read_bytes()))
        if packet['status']!='PASS_EXACT_ACCESSIBILITY' or packet['j_group'] in excluded or packet['family']=='Curve302-development':continue
        models=packet['prepared']['models'];base=models[0]['features'];obs=packet['observations']
        for ti in {r['target_index'] for r in obs}:
            hs=[min(r['logH'] for r in obs if r['target_index']==ti and r['model_index']==mi) for mi in range(len(models))]
            for mi,m in enumerate(models):train.append({'x':{k:m['features'][k]-base[k] for k in FINITE},'y':hs[mi]-hs[0]})
    _,fit=ridge(train,[train[0]],FINITE)
    write(CPU/'inputs.json',chosen)
    plan={'status':'FROZEN_BEFORE_CPU_RUN','cases':12,'selection':'Two lowest j hashes per six R17 families among retained broad next-direction banks with a separately certified larger endpoint and at least16 centres. No target coordinate or target accessibility used.',
        'arms':['factor_free','finite_selector'],'height':125000,'cpu_seconds_per_arm':40,'point_wall_seconds':5,
        'maximum_centres':48,'stop':'First independent finite column beyond the input subgroup, followed by standalone exact rank certificate; otherwise cap or bank exhaustion.',
        'order':'Alternate arm order by case ordinal; each arm in its own process.',
        'cost':'Charge Python plus child-process CPU from worker entry, including initial finite tables, exact seed verification, all candidate maps/features, discarded models, point calls and output rank proof. Cached input landscapes are common to both arms. Offline training cost reported separately.',
        'blinding':'Worker reads inputs.json and this fit only. All12 j groups and Curve302 excluded from training. Legacy final endpoints establish withholding eligibility only. This is a newly scheduled retrospective control, not previously unknown curves.',
        'fit':fit,'training_model_rows':len(train),'training_cpu_seconds':time.process_time()-start,
        'excluded_training_j_groups':sorted(excluded),'input_bindings':bindings,'inputs_sha256':digest((CPU/'inputs.json').read_bytes()),
        'feature_protocol_sha256':digest((OUT/'protocol.json').read_bytes()),
        'source_sha256':{p.name:digest(p.read_bytes()) for p in [Path(__file__),Path(__file__).with_name('finite_cancellation_features.py'),Path(__file__).with_name('pari_pointed_backend.py'),Path(__file__).with_name('lean_factor_free_pari_mapping.sage')]},
        'gp_sha256':digest(Path('/usr/bin/gp').read_bytes()),
        'success_gate':'At least as many next-direction recoveries as factor-free at the identical per-case CPU cap, and at least10 percent less aggregate charged CPU to those terminal outcomes. Report paired disagreements; this small fixed control does not establish a population success probability.'}
    write(CPU/'protocol.json',plan)
    print(json.dumps({'status':plan['status'],'cases':[(x['id'],x['family'],len(x['centres']),x['seed']['rank_lower_bound']) for x in chosen],'training_cpu_seconds':plan['training_cpu_seconds']}))

def cpu_now():
    r=resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time()+r.ru_utime+r.ru_stime

def run(case_id,arm):
    start=cpu_now();wall=time.monotonic()
    from finite_cancellation_features import prepare as chart_prepare, alarm
    from half_lattice_pointed_sieve import linear_combination
    from importlib.machinery import SourceFileLoader
    from v3_warm_engine import certified_state
    from future_point_admission import FinitePointAdmission
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import execute
    from memory_rank_certificate import checked_rank
    plan=json.loads((CPU/'protocol.json').read_text());assert digest((CPU/'inputs.json').read_bytes())==plan['inputs_sha256']
    for name,h in plan['source_sha256'].items():assert digest(Path(__file__).with_name(name).read_bytes())==h
    case=next(r for r in json.loads((CPU/'inputs.json').read_text()) if r['id']==case_id)
    seed=case['seed'];curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']]
    state=certified_state(curve,basis,seed['proof'])
    admission=FinitePointAdmission(curve,basis,prime_bound=500)
    mapper=SourceFileLoader('cpu_factor_free_mapper',str(Path(__file__).with_name('lean_factor_free_pari_mapping.sage'))).load_module()
    mapper.pari.allocatemem(256000000,silent=True)
    folder=CPU/'arms'/case_id/arm;folder.mkdir(parents=True,exist_ok=True)
    if (folder/'result.json').exists():raise FileExistsError('preserve prior completed timing')
    signal.signal(signal.SIGALRM,alarm);records=[];success=False;proof=None;status='BANK_EXHAUSTED';map_cpu=0.;point_cpu=0.;fit=plan['fit']
    for ci,word in enumerate(case['centres']):
        if cpu_now()-start>=plan['cpu_seconds_per_arm']:status='CPU_CAP';break
        before=cpu_now();signal.alarm(10)
        try:
            centre={'representative':word}
            if arm=='factor_free':
                mapping=mapper.mapping(curve,basis,centre);selection={'selected':'factor_free','models_built':1}
            else:
                anchor=linear_combination(curve,basis,word)
                prepared=chart_prepare(list(map(str,curve)),[list(map(str,anchor))],0,mapper)
                base=prepared['models'][0]['features'];scores=[]
                for model in prepared['models']:
                    z=[(model['features'][k]-base[k]-mu)/sd for k,mu,sd in zip(fit['features'],fit['means'],fit['scales'])]
                    scores.append(fit['coefficients'][0]+sum(a*b for a,b in zip(fit['coefficients'][1:],z)))
                mi=min(range(len(scores)),key=lambda i:scores[i]);mapping=prepared['models'][mi]['mapping']
                selection={'selected':prepared['models'][mi]['name'],'models_built':len(scores),'scores':scores,
                    'prepared_sha256':digest(canonical(prepared)),'features':[m['features'] for m in prepared['models']]}
            search=PointedQuarticSearch(state=state,centre={'coefficients':word},coordinate_policy=mapping['coordinate_policy'])
        except TimeoutError:
            records.append({'index':ci,'status':'MAP_TIMEOUT','cpu_seconds':cpu_now()-before});continue
        finally:signal.alarm(0)
        used=cpu_now()-before;map_cpu+=used
        remaining=plan['cpu_seconds_per_arm']-(cpu_now()-start)
        if remaining<=0:status='CPU_CAP_AFTER_PREPARATION';break
        beforepoint=cpu_now()
        transcript,points=execute(search,mapping,plan['height'],min(plan['point_wall_seconds'],remaining),plan['gp_sha256'])
        point_cpu+=cpu_now()-beforepoint
        row={'index':ci,'centre':word,'selection':selection,'mapping':mapping,'search':transcript,'preparation_cpu_seconds':used}
        for point in points:
            if admission.consider(point)['status']=='INDEPENDENT_FINITE_COLUMN':
                success=True;row['new_point']=list(map(str,point));break
        write(folder/f'chart-{ci:03d}.json',row);records.append({'index':ci,'status':transcript['status'],
            'selection':selection['selected'],'new_direction':success,'cpu_after':cpu_now()-start})
        if success:
            proof=checked_rank(curve,admission.points,admission.primes,seed['proof']['no_rational_2_torsion_prime'])
            write(folder/'rank-input.json',{'curve':seed['curve'],'points':[list(map(str,p)) for p in admission.points],'proof':proof})
            status='CERTIFIED_NEXT_DIRECTION';break
    result={'status':status,'case':case_id,'family':case['family'],'arm':arm,'success':success,
        'initial_rank':len(basis),'rank_lower_bound':len(admission.points),'charts':records,
        'cpu_seconds':cpu_now()-start,'wall_seconds':time.monotonic()-wall,
        'map_and_feature_cpu_seconds':map_cpu,'point_backend_cpu_seconds':point_cpu,
        'protocol_sha256':digest((CPU/'protocol.json').read_bytes()),'point_search_calls':sum('selection' in r for r in records),
        'boundaries':'Known-curve control. Exact independence beyond input subgroup; no novelty or rank upper bound. Cached V3 centre bank reused, all per-worker costs charged.'}
    write(folder/'result.json',result);print(json.dumps({k:v for k,v in result.items() if k!='charts'}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['prepare','run']);p.add_argument('--case');p.add_argument('--arm',choices=['factor_free','finite_selector']);a=p.parse_args()
    if a.command=='prepare':prepare()
    else:run(a.case,a.arm)
