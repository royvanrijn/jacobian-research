#!/usr/bin/env sage -python
"""Fail-closed policy/map/transcript replay plus independent Sage rank proof."""
from fractions import Fraction as F
import gzip
from importlib.machinery import SourceFileLoader
import json
from math import gcd
import time

from finite_cancellation_corpus import ROOT, digest, write
from cancellation_scheduler_prepare import OUT, RAW, CAS


def verify_arm(case, arm, dest, plan, result):
    from cancellation_scheduler import choose, prepare, AnchorExposure, CostModel
    from finite_cancellation_validation_replay import verify_tree
    from half_lattice_pointed_sieve import binary_transform, linear_combination
    from v3_warm_engine import certified_state
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import replay
    from finite_cancellation_features import xmap
    start=time.process_time()
    seed=case['seed'];curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']]
    state=certified_state(curve,basis,seed['proof'])
    mapper=SourceFileLoader('scheduler_replay_mapper',str(CAS/'lean_factor_free_pari_mapping.sage')).load_module()
    mapper.pari.allocatemem(256000000,silent=True)
    events=json.loads((dest/'events.json').read_text())
    if digest((dest/'events.json').read_bytes())!=result['events_sha256']:
        raise ArithmeticError('event log changed')
    anchors={};maps={};cost=CostModel(plan['fit']);next_index=0;new=[];calls=0;leaves=0;unknown=0
    for event in events:
        decision=event['decision'];ci=decision['index'];word=case['centres'][ci]
        if arm!='factor_free':
            expected=choose(anchors,cost,plan['fit'],next_index,len(case['centres']))
            if expected!=decision:raise ArithmeticError('scheduler decision changed')
        if decision['kind']=='prepare':
            if ci!=next_index:raise ArithmeticError('anchor order changed')
            next_index+=1
            if event['status']=='PREPARATION_UNKNOWN':
                unknown+=1;continue
            path=dest/f'prepared-{ci:03d}.json.gz'
            if digest(path.read_bytes())!=event['sha256']:raise ArithmeticError('prepared artifact changed')
            prepared=json.loads(gzip.decompress(path.read_bytes()))
            if arm=='factor_free':
                expected={'models':[{'name':'factor_free','mapping':mapper.mapping(curve,basis,{'representative':word})}]}
            else:
                anchor=linear_combination(curve,basis,word)
                expected=prepare(seed['curve'],list(map(str,anchor)),mapper,
                    plan['fit']['local_weight'] if arm=='adaptive_local' else 0.)
                n,d,q=xmap(curve[3],curve[4],anchor,prepared['models'][0]['mapping'])
                for partition in prepared['partitions']:
                    if partition['tree'] is not None:
                        leaves+=verify_tree(partition['tree'],n,d,q)
                    if abs(sum(c['probability'] for c in partition['categories'])-1)>1e-12:
                        raise ArithmeticError('inconsistent categorical probabilities')
                for model in prepared['models'][1:]:
                    T=tuple(map(F,model['mapping']['second_matrix']))
                    values=binary_transform(n,T)+binary_transform(d,T)
                    content=gcd(*(int(z) for z in values))
                    if content!=model['prime']**model['content_exponent']:
                        raise ArithmeticError('independent content replay failed')
                anchors[ci]=AnchorExposure(prepared,plan['fit'],ci)
            if expected!=prepared:raise ArithmeticError('target-blind model/feature reconstruction differs')
            maps[ci]=prepared
        elif decision['kind']=='search':
            if arm=='factor_free' and (ci!=next_index-1 or decision!={'kind':'search','index':ci,'model':0,'height':125000}):
                raise ArithmeticError('baseline schedule changed')
            path=dest/event['call'];row=json.loads(path.read_text())
            expected_record=result['calls'][calls]
            if expected_record['file']!=path.name or digest(path.read_bytes())!=expected_record['sha256']:
                raise ArithmeticError('call record changed')
            mapping=maps[ci]['models'][decision['model']]['mapping']
            if row['mapping']!=mapping or row['decision']!=decision or row['centre']!=word:
                raise ArithmeticError('map/centre decision binding failed')
            if row['search']['height_bound']!=decision['height'] or row['call_cpu_seconds']!=event['cpu_seconds']:
                raise ArithmeticError('box/cost binding failed')
            search=PointedQuarticSearch(state=state,centre={'coefficients':word},coordinate_policy=mapping['coordinate_policy'])
            points=replay(search,mapping,row['search'])
            if row['new_point']:
                if tuple(map(F,row['new_point'])) not in points:raise ArithmeticError('new point missing from transcript')
                new.append(row['new_point'])
            complete=row['search']['status']=='bounded_search_complete'
            if arm!='factor_free':
                anchors[ci].observe(decision['job'],complete)
                cost.observe(decision['height'],event['cpu_seconds'],complete)
            calls+=1
        else:raise ArithmeticError('unrecognized event')
    if calls!=len(result['calls']) or len(new)!=int(result['candidate']):
        raise ArithmeticError('candidate/call accounting differs')
    if new:
        packet=json.loads((dest/'rank-input.json').read_text())
        if packet['curve']!=seed['curve'] or packet['points']!=seed['points']+new:
            raise ArithmeticError('certificate subgroup does not match discovered extension')
    else:packet=seed
    # Different implementation: complete Sage E(Fp)/2E(Fp) cosets and point
    # arithmetic, not the producer's portable finite-column tables.
    sage_rank=SourceFileLoader('scheduler_independent_rank',str(CAS/'verify_finite_cancellation_cpu.sage')).load_module().sage_rank
    rank=sage_rank(packet)
    if rank['rank']!=len(basis)+len(new):raise ArithmeticError('independent rank mismatch')
    return {'status':'PASS','rank':rank,'calls':calls,'residue_leaves':leaves,
            'preparation_unknowns':unknown,'cpu_seconds':time.process_time()-start,
            'boundary':'Every selected action, polynomial map, square witness and subgroup extension replayed. Sage finite groups independently certify rank lower bounds; timeout or finite-column ambiguity is not an exclusion.'}


def main():
    plan=json.loads((OUT/'protocol.json').read_text())
    inputs={c['id']:c for c in json.loads((OUT/'inputs.json').read_text())}
    supervision=json.loads((OUT/'supervision.json').read_text());rows=[];start=time.process_time()
    if supervision['status']!='COMPLETE':raise ArithmeticError('benchmark incomplete')
    for receipt in supervision['records']:
        dest=RAW/'arms'/receipt['case']/receipt['arm']
        if receipt['status']!='COMPLETE':raise ArithmeticError('incomplete process accounting')
        result=json.loads((dest/'result.json').read_text())
        if digest((dest/'result.json').read_bytes())!=receipt['result_sha256']:raise ArithmeticError('result changed')
        verification=verify_arm(inputs[receipt['case']],receipt['arm'],dest,plan,result)
        rows.append({'case':receipt['case'],'arm':receipt['arm'],**verification})
    write(OUT/'replay.json',{'status':'PASS','rows':rows,'cpu_seconds':time.process_time()-start,
        'protocol_sha256':digest((OUT/'protocol.json').read_bytes())})
    print(json.dumps({'status':'PASS','arms':len(rows),'cpu':time.process_time()-start}),flush=True)


if __name__=='__main__':main()
