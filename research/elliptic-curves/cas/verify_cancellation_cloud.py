"""Independent arithmetic replay of full-cloud admission and yield scheduling."""
from fractions import Fraction as F
import gzip
from importlib.machinery import SourceFileLoader
import json
from math import gcd
import time

from cancellation_cloud_training import OUT,RAW,CAS
from finite_cancellation_corpus import ROOT,digest,write


def verify_arm(case,arm,dest,plan,result):
    from cancellation_scheduler import prepare,CostModel
    from cancellation_cloud_policy import choose,CloudExposure
    from half_lattice_pointed_sieve import binary_transform,linear_combination
    from finite_cancellation_features import xmap
    from finite_cancellation_validation_replay import verify_tree
    from v3_warm_engine import certified_state
    from future_point_admission import FinitePointAdmission
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import replay
    start=time.process_time();seed=case['seed'];curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']]
    state=certified_state(curve,basis,seed['proof']);admission=FinitePointAdmission(curve,basis,prime_bound=500)
    mapper=SourceFileLoader('cloud_replay_mapper',str(CAS/'lean_factor_free_pari_mapping.sage')).load_module()
    mapper.pari.allocatemem(256000000,silent=True)
    events=json.loads((dest/'events.json').read_text())
    if digest((dest/'events.json').read_bytes())!=result['events_sha256']:raise ArithmeticError('event bytes changed')
    anchors={};maps={};cost=CostModel(plan['fit']);next_index=0;calls=0;leaves=0;unknown=0;new=[]
    for event in events:
        remaining=plan['target_directions']-len(new)
        if remaining<=0:raise ArithmeticError('search continues beyond sealed target')
        decision=event['decision'];ci=decision['index'];word=case['centres'][ci]
        if arm!='factor_free' and decision!=choose(anchors,cost,plan['fit'],next_index,len(case['centres']),remaining):raise ArithmeticError('policy choice differs')
        if decision['kind']=='prepare':
            if ci!=next_index:raise ArithmeticError('anchor order differs')
            next_index+=1
            if event['status']=='PREPARATION_UNKNOWN':unknown+=1;continue
            path=dest/f'prepared-{ci:03d}.json.gz';prepared=json.loads(gzip.decompress(path.read_bytes()))
            if digest(path.read_bytes())!=event['sha256']:raise ArithmeticError('preparation bytes changed')
            if arm=='factor_free':expected={'models':[{'name':'factor_free','mapping':mapper.mapping(curve,basis,{'representative':word})}]}
            else:
                anchor=linear_combination(curve,basis,word);expected=prepare(seed['curve'],list(map(str,anchor)),mapper,plan['fit']['local_weight'])
                n,d,q=xmap(curve[3],curve[4],anchor,prepared['models'][0]['mapping'])
                for part in prepared['partitions']:
                    if part['tree'] is not None:leaves+=verify_tree(part['tree'],n,d,q)
                for model in prepared['models'][1:]:
                    T=tuple(map(F,model['mapping']['second_matrix']))
                    content=gcd(*(int(z) for z in binary_transform(n,T)+binary_transform(d,T)))
                    if content!=model['prime']**model['content_exponent']:raise ArithmeticError('content replay differs')
                anchors[ci]=CloudExposure(prepared,plan['fit'],ci)
            if prepared!=expected:raise ArithmeticError('prepared feature/map differs')
            maps[ci]=prepared
        elif decision['kind']=='search':
            if arm=='factor_free' and (ci!=next_index-1 or decision!={'kind':'search','index':ci,'model':0,'height':125000}):raise ArithmeticError('baseline schedule differs')
            record=result['calls'][calls];path=dest/event['call'];row=json.loads(path.read_text())
            if record['file']!=path.name or digest(path.read_bytes())!=record['sha256']:raise ArithmeticError('point receipt differs')
            mapping=maps[ci]['models'][decision['model']]['mapping']
            if row['mapping']!=mapping or row['decision']!=decision or row['centre']!=word:raise ArithmeticError('search binding differs')
            if row['search']['height_bound']!=decision['height'] or row['call_cpu_seconds']!=event['cpu_seconds']:raise ArithmeticError('box/cost binding differs')
            search=PointedQuarticSearch(state=state,centre={'coefficients':word},coordinate_policy=mapping['coordinate_policy'])
            returned=replay(search,mapping,row['search']);found=[]
            for point in returned:
                if admission.consider(point)['status']=='INDEPENDENT_FINITE_COLUMN':found.append(list(map(str,point)))
            if found!=row['new_points'] or len(found)!=event['new_directions'] or len(found)!=record['new_directions']:raise ArithmeticError('cloud admission differs')
            new.extend(found)
            if row['rank_after']!=len(basis)+len(new):raise ArithmeticError('cloud rank trace differs')
            complete=row['search']['status']=='bounded_search_complete'
            if not complete and found:raise ArithmeticError('unexpected partial-cloud inference')
            if arm!='factor_free':
                anchors[ci].observe_cloud(decision['job'],complete,len(found));cost.observe(decision['height'],event['cpu_seconds'],complete)
            calls+=1
        else:raise ArithmeticError('unrecognized event')
    if calls!=len(result['calls']) or len(new)!=result['candidate_directions']:raise ArithmeticError('candidate accounting differs')
    packet=json.loads((dest/'rank-input.json').read_text())
    if packet['curve']!=seed['curve'] or packet['points']!=seed['points']+new:raise ArithmeticError('certificate subgroup differs')
    checker=SourceFileLoader('cloud_independent_sage',str(CAS/'verify_finite_cancellation_cpu.sage')).load_module().sage_rank
    independent=checker(packet)
    if independent['rank']!=len(basis)+len(new):raise ArithmeticError('independent cloud rank differs')
    return {'status':'PASS','rank':independent,'new_directions':len(new),'calls':calls,
        'residue_leaves':leaves,'preparation_unknowns':unknown,'cpu_seconds':time.process_time()-start,
        'boundary':'Full clouds, all incremental finite-column decisions, maps and residue policy replayed. Complete Sage finite-group cosets independently certify the enlarged subgroup. The search subgroup and bank stay fixed throughout the arm.'}


def main():
    plan=json.loads((OUT/'protocol.json').read_text());cases={c['id']:c for c in json.loads((OUT/'inputs.json').read_text())}
    supervision=json.loads((OUT/'supervision.json').read_text());rows=[];start=time.process_time()
    if supervision['status']!='COMPLETE':raise ArithmeticError('stage incomplete')
    for receipt in supervision['records']:
        dest=RAW/'arms'/receipt['case']/receipt['arm'];result=json.loads((dest/'result.json').read_text())
        if receipt['result_sha256']!=digest((dest/'result.json').read_bytes()):raise ArithmeticError('result changed')
        rows.append({'case':receipt['case'],'arm':receipt['arm'],**verify_arm(cases[receipt['case']],receipt['arm'],dest,plan,result)})
    write(OUT/'replay.json',{'status':'PASS','rows':rows,'cpu_seconds':time.process_time()-start,'protocol_sha256':digest((OUT/'protocol.json').read_bytes())})


if __name__=='__main__':main()
