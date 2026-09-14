#!/usr/bin/env sage -python
"""Costed multi-direction exposure on a fixed certified subgroup and bank."""
import argparse
from fractions import Fraction as F
import gzip
from importlib.machinery import SourceFileLoader
import json
from math import ceil
import resource
import signal

from cancellation_cloud_training import OUT, RAW, CAS
from cancellation_scheduler_cpu import cpu,atomic
from finite_cancellation_corpus import ROOT, canonical, digest, write


def run(case_id,arm,*,development_context=None):
    start=cpu()
    from cancellation_scheduler import prepare,CostModel
    from cancellation_cloud_policy import CloudExposure,choose
    from finite_cancellation_features import alarm
    from half_lattice_pointed_sieve import linear_combination
    from v3_warm_engine import certified_state
    from future_point_admission import FinitePointAdmission
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import execute
    from memory_rank_certificate import checked_rank
    if development_context is None:
        plan=json.loads((OUT/'protocol.json').read_text());ph=digest((OUT/'protocol.json').read_bytes())
        if digest((OUT/'inputs.json').read_bytes())!=plan['input_sha256']['inputs.json']:raise ArithmeticError('input changed')
        case=next(c for c in json.loads((OUT/'inputs.json').read_text()) if c['id']==case_id)
        dest=RAW/'arms'/case_id/arm
    else:
        plan,case,dest=development_context;ph='DEVELOPMENT_ONLY:'+digest(canonical(plan))
    if arm not in plan['arms']:raise ValueError('unsealed arm')
    for name,h in plan['source_sha256'].items():
        if digest((ROOT/name).read_bytes())!=h:raise ArithmeticError('source changed: '+name)
    dest.mkdir(parents=True,exist_ok=True)
    if (dest/'result.json').exists() or (dest/'events.json').exists():raise FileExistsError('preserve prior arm')
    resource.setrlimit(resource.RLIMIT_CPU,(plan['hard_process_cpu_seconds']-5,plan['hard_process_cpu_seconds']))
    signal.signal(signal.SIGALRM,alarm)
    seed=case['seed'];curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']]
    state=certified_state(curve,basis,seed['proof']);admission=FinitePointAdmission(curve,basis,prime_bound=500)
    mapper=SourceFileLoader('cloud_factor_free',str(CAS/'lean_factor_free_pari_mapping.sage')).load_module()
    mapper.pari.allocatemem(256000000,silent=True)
    cost=CostModel(plan['fit']);anchors={};maps={};events=[];calls=[];next_index=0
    status='SURROGATE_EXHAUSTED';target=plan['target_directions'];initial=len(basis)
    components={'initial':cpu()-start,'preparation':0.,'search':0.,'admission':0.}
    while len(calls)<plan['maximum_calls']:
        remaining_time=plan['search_cpu_seconds']-(cpu()-start)
        if remaining_time<=0:status='CPU_CAP';break
        remaining_gain=target-(len(admission.points)-initial)
        if remaining_gain<=0:status='PENDING_INDEPENDENT_TARGET_CERTIFICATION';break
        if arm=='factor_free':
            decision={'kind':'prepare','index':next_index} if next_index<len(case['centres']) else {'kind':'exhausted'}
        else:decision=choose(anchors,cost,plan['fit'],next_index,len(case['centres']),remaining_gain)
        if decision['kind']=='exhausted':break
        if decision['kind']=='prepare':
            ci=decision['index'];word=case['centres'][ci];before=cpu();signal.alarm(max(1,min(10,ceil(remaining_time))))
            try:
                if arm=='factor_free':prepared={'models':[{'name':'factor_free','mapping':mapper.mapping(curve,basis,{'representative':word})}]}
                else:
                    anchor=linear_combination(curve,basis,word)
                    prepared=prepare(seed['curve'],list(map(str,anchor)),mapper,plan['fit']['local_weight'])
                    anchors[ci]=CloudExposure(prepared,plan['fit'],ci)
                maps[ci]=prepared;path=dest/f'prepared-{ci:03d}.json.gz'
                with gzip.GzipFile(str(path),'wb',mtime=0) as f:f.write(canonical(prepared))
                event={'decision':decision,'status':'PREPARED','sha256':digest(path.read_bytes())}
            except (TimeoutError,ArithmeticError,ValueError,AssertionError) as e:
                event={'decision':decision,'status':'PREPARATION_UNKNOWN','error':type(e).__name__+': '+str(e)}
            finally:signal.alarm(0)
            spent=cpu()-before;components['preparation']+=spent;event['cpu_seconds']=spent
            events.append(event);next_index+=1;atomic(dest/'events.json',events)
            if arm!='factor_free' or event['status']!='PREPARED':continue
            decision={'kind':'search','index':ci,'model':0,'height':125000}
        ci,mi,height=decision['index'],decision['model'],decision['height']
        remaining_time=plan['search_cpu_seconds']-(cpu()-start)
        if remaining_time<=0:status='CPU_CAP_AFTER_PREPARATION';break
        mapping=maps[ci]['models'][mi]['mapping'];word=case['centres'][ci];before=cpu()
        search=PointedQuarticSearch(state=state,centre={'coefficients':word},coordinate_policy=mapping['coordinate_policy'])
        transcript,points=execute(search,mapping,height,min(plan['point_wall_seconds'],remaining_time),plan['gp_sha256'])
        elapsed=cpu()-before;components['search']+=elapsed;complete=transcript['status']=='bounded_search_complete'
        before=cpu();new=[]
        for point in points:
            if admission.consider(point)['status']=='INDEPENDENT_FINITE_COLUMN':new.append(list(map(str,point)))
        components['admission']+=cpu()-before
        if not complete and new:raise ArithmeticError('unexpected partial cloud; no completed-prefix inference')
        row={'decision':decision,'centre':word,'mapping':mapping,'search':transcript,
            'new_points':new,'rank_after':len(admission.points),'call_cpu_seconds':elapsed}
        path=dest/f'call-{len(calls):03d}.json';write(path,row)
        calls.append({'file':path.name,'sha256':digest(path.read_bytes()),'decision':decision,
            'status':transcript['status'],'cpu_seconds':elapsed,'new_directions':len(new)})
        event={'decision':decision,'status':transcript['status'],'call':path.name,'cpu_seconds':elapsed,'new_directions':len(new)}
        events.append(event)
        if arm!='factor_free':
            anchors[ci].observe_cloud(decision['job'],complete,len(new));cost.observe(height,elapsed,complete)
        atomic(dest/'events.json',events)
        if len(admission.points)>=initial+target:status='PENDING_INDEPENDENT_TARGET_CERTIFICATION';break
    if len(calls)>=plan['maximum_calls'] and len(admission.points)<initial+target:status='CALL_CAP'
    search_end=cpu()
    proof=checked_rank(curve,admission.points,admission.primes,seed['proof']['no_rational_2_torsion_prime'])
    write(dest/'rank-input.json',{'curve':seed['curve'],'points':[list(map(str,p)) for p in admission.points],'proof':proof})
    result={'status':status,'case':case_id,'family':case['family'],'stratum':case['stratum'],'arm':arm,
        'initial_rank':initial,'rank_lower_bound':initial,'candidate_directions':len(admission.points)-initial,
        'success':False,'new_directions':0,'target_directions':target,'calls':calls,
        'development_only':development_context is not None,'protocol_sha256':ph,
        'events_sha256':digest((dest/'events.json').read_bytes()),'components':components,
        'search_phase_cpu_seconds':search_end-start,
        'boundary':'All searches use the original certified subgroup and its fixed bank. Admission accumulates every cloud direction for the final certificate; it does not invent an enlarged search landscape.'}
    from verify_cancellation_cloud import verify_arm
    verification=verify_arm(case,arm,dest,plan,result);write(dest/'independent-verification.json',verification)
    if verification['status']!='PASS':raise ArithmeticError('independent replay failed')
    actual=verification['rank']['rank']-initial;success=actual>=target
    result.update(success=success,new_directions=actual,rank_lower_bound=initial+actual,
        status='CERTIFIED_MULTI_DIRECTION_TARGET' if success else status,
        independent_verification_sha256=digest((dest/'independent-verification.json').read_bytes()),
        proof_and_independent_replay_cpu_seconds=cpu()-search_end,total_internal_cpu_seconds=cpu()-start)
    write(dest/'result.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('calls','components','boundary')}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',required=True);p.add_argument('--arm',required=True);a=p.parse_args();run(a.case,a.arm)
