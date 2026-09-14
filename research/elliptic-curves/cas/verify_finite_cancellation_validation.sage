#!/usr/bin/env sage -python
"""Replay the clean ablation, point transcripts, independent ranks and costs."""
import argparse
from collections import Counter,defaultdict
from fractions import Fraction as F
import gzip
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import time

from finite_cancellation_corpus import OUT as V1,LOCAL,canonical,digest,pointkey,short,write
from finite_cancellation_validation_audit import OUT
from finite_cancellation_validation_features import FEATURES,model_features,prepare,select
from finite_cancellation_validation_replay import verify_prepared

RAW=LOCAL/'finite-cancellation-validation-v2'


def training():
    from analyze_finite_cancellation import ridge
    start=time.process_time();fit=json.loads((OUT/'fit.json').read_text());roster=json.loads((OUT/'roster.json').read_text())
    assert digest((OUT/'training-features.json.gz').read_bytes())==fit['training_features_sha256']
    stored=json.loads(gzip.decompress((OUT/'training-features.json.gz').read_bytes()));excluded=set(roster['excluded_training_j_groups'])
    rows=[];leaves=0;models=0
    for ri,record in enumerate(stored):
        path=V1/'cases'/record['source'];raw=path.read_bytes();assert digest(raw)==fit['source_cases_sha256'][record['source']]
        packet=json.loads(gzip.decompress(raw));assert packet['j_group']==record['j_group'] and packet['j_group'] not in excluded
        assert packet['family']!='Curve302-development'
        prep=packet['prepared'];fs=[];ms=[]
        for mi,m in enumerate(prep['models']):
            n,d,q=([int(v) for v in m[k]] for k in ('N','D','q'))
            f,local=model_features(n,d,q,True)
            assert f==record['features'][mi]
            fs.append(f);ms.append({**m,'features':f,'local':local})
        leaves+=verify_prepared({**prep,'models':ms});models+=len(ms)
        for ti in sorted({r['target_index'] for r in packet['observations']}):
            hs=[min(o['logH'] for o in packet['observations'] if o['target_index']==ti and o['model_index']==mi) for mi in range(len(fs))]
            for mi,f in enumerate(fs):rows.append({'group':packet['j_group'],'x':{k:f[k]-fs[0][k] for k in FEATURES['real_q_g']},'y':hs[mi]-hs[0]})
        if (ri+1)%500==0:print(json.dumps({'training_anchors_replayed':ri+1,'leaves':leaves,'cpu_seconds':time.process_time()-start}),flush=True)
    assert len(rows)==fit['training_model_rows']
    assert sorted({r['group'] for r in rows})==fit['training_j_groups']
    for arm,names in FEATURES.items():
        _,actual=ridge(rows,[rows[0]],names)
        assert actual==fit['fits'][arm]
    result={'status':'PASS_TRAINING_FEATURES_PARTITIONS_AND_EXCLUSION','anchors':len(stored),'models':models,
        'leaves':leaves,'training_rows':len(rows),'excluded_groups':len(excluded),'cpu_seconds':time.process_time()-start,
        'fit_sha256':digest((OUT/'fit.json').read_bytes()),'source_sha256':digest(Path(__file__).read_bytes())}
    write(OUT/'training-verified.json',result);print(json.dumps(result),flush=True)


def cpu():
    from half_lattice_pointed_sieve import linear_combination
    from v3_warm_engine import certified_state
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import replay
    from memory_rank_certificate import checked_rank
    rankcheck=SourceFileLoader('validation_independent_rank',str(Path(__file__).with_name('verify_finite_cancellation_cpu.sage'))).load_module().sage_rank
    start=time.process_time();plan=json.loads((OUT/'protocol.json').read_text());ph=digest((OUT/'protocol.json').read_bytes())
    supervision=json.loads((OUT/'supervision.json').read_text());inputs=json.loads((OUT/'inputs.json').read_text())
    assert supervision['status']=='COMPLETE' and len(supervision['records'])==3*len(inputs)
    for name,h in plan['source_sha256'].items():assert digest(Path(__file__).with_name(name).read_bytes())==h
    for name,h in plan['input_sha256'].items():assert digest((OUT/name).read_bytes())==h
    roster=json.loads((OUT/'roster.json').read_text())
    assert digest((OUT/'oracle.json.gz').read_bytes())==roster['oracle_sha256']
    oracle={o['id']:o for o in json.loads(gzip.decompress((OUT/'oracle.json.gz').read_bytes()))}
    mapper=SourceFileLoader('validation_replay_mapper',str(Path(__file__).with_name('lean_factor_free_pari_mapping.sage'))).load_module()
    mapper.pari.allocatemem(256000000,silent=True)
    rows=[];leaves=0;proofs=0;calls=0;shared_choices=Counter()
    for ci,case in enumerate(inputs):
        seed=case['seed'];curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']]
        state=certified_state(curve,basis,seed['proof']);cm,cp=short(seed['curve'],seed['points'])
        assert cm==oracle[case['id']]['curve']
        withheld={pointkey(p) for p in oracle[case['id']]['points']}-{pointkey(p) for p in cp}
        assert withheld
        q_inputs={}
        for arm in plan['arms']:
            dest=RAW/'arms'/case['id']/arm;supervisor=json.loads((dest/'supervisor.json').read_text())
            assert supervisor in supervision['records'] and supervisor['protocol_sha256']==ph
            if supervisor['status']!='COMPLETE':
                rows.append({'case':case['id'],'family':case['family'],'stratum':case['stratum'],'arm':arm,
                    'status':supervisor['status'],'success':False,'verification':'UNKNOWN_WORKER_FAILURE',
                    'charged_cpu_seconds':supervisor['charged_cpu_seconds'],'point_calls':None,
                    'initial_rank':len(basis),'literal_withheld_points_recovered':None})
                continue
            result=json.loads((dest/'result.json').read_text())
            assert digest((dest/'result.json').read_bytes())==supervisor['result_sha256'] and result['protocol_sha256']==ph
            hits=set();new=[];pointcalls=0;statuses=Counter();selected=Counter()
            for item in result['charts']:
                i=item['index'];path=dest/f'chart-{i:03d}.json';assert digest(path.read_bytes())==item['row_sha256']
                row=json.loads(path.read_text());word=case['centres'][i]
                assert row['centre']==word and row['status']==item['status']
                statuses[row['status']]+=1
                if not row.get('mapping'):
                    assert not row['point_call'] and row['status'] in ('MAP_TIMEOUT','MAP_FAILURE');continue
                mapping=row['mapping'];selected[row['selection']['selected']]+=1
                if arm=='factor_free':
                    assert mapping==mapper.mapping(curve,basis,{'representative':word})
                else:
                    pp=dest/f'prepared-{i:03d}.json.gz';assert digest(pp.read_bytes())==row['selection']['prepared_file_sha256']
                    raw=gzip.decompress(pp.read_bytes());assert digest(raw)==row['selection']['prepared_sha256']
                    saved=json.loads(raw);anchor=linear_combination(curve,basis,word)
                    actual=prepare(seed['curve'],list(map(str,anchor)),mapper,arm=='real_q_g')
                    assert actual==saved
                    leaves+=verify_prepared(saved)
                    mi,scores=select(saved,plan['fits'][arm])
                    assert scores==row['selection']['scores'] and mapping==saved['models'][mi]['mapping']
                    common=[{k:m['features'][k] for k in FEATURES['real_q']} for m in saved['models']]
                    if i in q_inputs:
                        assert common==q_inputs[i];shared_choices['identical_common_feature_checks']+=1
                    else:q_inputs[i]=common
                if not row['point_call']:
                    assert row['status']=='CPU_CAP_AFTER_PREPARATION';continue
                search=PointedQuarticSearch(state=state,centre={'coefficients':word},coordinate_policy=mapping['coordinate_policy'])
                pts=replay(search,mapping,row['search']);_,normalized=short(seed['curve'],pts)
                hits|={pointkey(p) for p in normalized}&withheld
                if row.get('new_point'):
                    assert tuple(map(F,row['new_point'])) in pts
                    new.append(row['new_point'])
                pointcalls+=1
            assert pointcalls==result['point_search_calls']
            rank_result=None
            if result['success']:
                packet=json.loads((dest/'rank-input.json').read_text())
                assert packet['points'][:len(basis)]==seed['points'] and packet['points'][len(basis):]==new and len(new)==1
                proof=packet['proof'];exact=checked_rank(curve,[tuple(map(F,p)) for p in packet['points']],
                    [s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'])
                assert json.loads(json.dumps(exact))==proof
                rank_result=rankcheck(packet);proofs+=1
            else:assert not new
            r={'case':case['id'],'family':case['family'],'stratum':case['stratum'],'arm':arm,'initial_rank':len(basis),
                'status':result['status'],'success':result['success'],'verification':'PASS_EXACT_REPLAY',
                'rank_lower_bound':result['rank_lower_bound'],'rank_replay':rank_result,
                'literal_withheld_points_recovered':len(hits),'point_calls':pointcalls,
                'charged_cpu_seconds':supervisor['charged_cpu_seconds'],'wall_seconds':supervisor['wall_seconds'],
                'result_sha256':supervisor['result_sha256'],'models_selected':dict(selected),'chart_statuses':dict(statuses),
                'cpu_components':{k:result[k] for k in ['initial_verification_cpu_seconds','map_and_feature_cpu_seconds','point_backend_cpu_seconds','admission_cpu_seconds','rank_proof_cpu_seconds']}}
            rows.append(r);calls+=pointcalls
        if (ci+1)%5==0:print(json.dumps({'cases_replayed':ci+1,'calls':calls,'rank_proofs':proofs,'cpu_seconds':time.process_time()-start}),flush=True)
    result={'status':'PASS_COMPLETED_ARITHMETIC_WITH_FAILURES_EXPLICIT' if any(r['verification']!='PASS_EXACT_REPLAY' for r in rows) else 'PASS_ALL_ARMS_MAPS_TRANSCRIPTS_AND_INDEPENDENT_RANKS',
        'rows':rows,'cases':len(inputs),'point_calls':calls,'independent_gain_proofs':proofs,'residue_leaves':leaves,
        **shared_choices,'protocol_sha256':ph,'source_sha256':digest(Path(__file__).read_bytes()),'cpu_seconds':time.process_time()-start,
        'boundary':'Point transcript completeness trusts the pinned GP enumerator. Exact maps and leaves independently verified; every positive subgroup checked using portable signatures and complete Sage finite quotient arithmetic. Unsuccessful exposure is not a rank upper bound.'}
    write(OUT/'verified.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('command',choices=['training','cpu']);a=p.parse_args()
    globals()[a.command]()
