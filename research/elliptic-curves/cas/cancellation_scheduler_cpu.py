#!/usr/bin/env sage -python
"""One sealed adaptive-policy arm, ending with independent rank certification."""
import argparse
from fractions import Fraction as F
import gzip
from importlib.machinery import SourceFileLoader
import json
from math import ceil
import os
import resource
import signal
import time

from cancellation_scheduler_prepare import OUT, RAW, CAS
from finite_cancellation_corpus import ROOT, canonical, digest, write


def cpu():
    r=resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time()+r.ru_utime+r.ru_stime


def atomic(path, value):
    temp=path.with_suffix(path.suffix+'.tmp')
    write(temp,value)
    os.replace(temp,path)


def run(case_id, arm, *, development_context=None):
    start=cpu()
    from cancellation_scheduler import prepare, choose, AnchorExposure, CostModel
    from finite_cancellation_features import alarm
    from half_lattice_pointed_sieve import linear_combination
    from v3_warm_engine import certified_state
    from future_point_admission import FinitePointAdmission
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import execute
    from memory_rank_certificate import checked_rank
    if development_context is None:
        plan=json.loads((OUT/'protocol.json').read_text())
        ph=digest((OUT/'protocol.json').read_bytes())
        if digest((OUT/'inputs.json').read_bytes())!=plan['input_sha256']['inputs.json']:
            raise ArithmeticError('worker input changed')
        case=next(c for c in json.loads((OUT/'inputs.json').read_text()) if c['id']==case_id)
        dest=RAW/'arms'/case_id/arm
    else:
        # Explicit API for isolated development controls, unavailable on the
        # sealed production CLI. These destinations never enter supervision.
        plan,case,dest=development_context
        ph='DEVELOPMENT_ONLY:'+digest(canonical(plan))
    if arm not in plan['arms']:raise ValueError('unsealed arm')
    for name,h in plan['source_sha256'].items():
        if digest((ROOT/name).read_bytes())!=h:raise ArithmeticError(f'source changed: {name}')
    dest.mkdir(parents=True,exist_ok=True)
    if (dest/'result.json').exists() or (dest/'events.json').exists():
        raise FileExistsError('preserve earlier or interrupted arm')
    resource.setrlimit(resource.RLIMIT_CPU,(plan['hard_process_cpu_seconds']-5,plan['hard_process_cpu_seconds']))
    signal.signal(signal.SIGALRM,alarm)
    seed=case['seed'];curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']]
    state=certified_state(curve,basis,seed['proof'])
    admission=FinitePointAdmission(curve,basis,prime_bound=500)
    mapper=SourceFileLoader('scheduler_factor_free',str(CAS/'lean_factor_free_pari_mapping.sage')).load_module()
    mapper.pari.allocatemem(256000000,silent=True)
    anchors={}; maps={}; events=[]; calls=[]; next_index=0; status='SURROGATE_EXHAUSTED'; success=False
    cost=CostModel(plan['fit']); components={'initial':cpu()-start, 'preparation':0., 'search':0., 'admission':0.}
    while len(calls)<plan['maximum_calls']:
        remaining=plan['search_cpu_seconds']-(cpu()-start)
        if remaining<=0:status='CPU_CAP';break
        if arm=='factor_free':
            decision={'kind':'prepare','index':next_index} if next_index<len(case['centres']) else {'kind':'exhausted'}
        else:
            decision=choose(anchors,cost,plan['fit'],next_index,len(case['centres']))
        if decision['kind']=='exhausted':break
        if decision['kind']=='prepare':
            ci=decision['index'];word=case['centres'][ci];before=cpu()
            signal.alarm(max(1,min(10,ceil(remaining))))
            try:
                if arm=='factor_free':
                    mapping=mapper.mapping(curve,basis,{'representative':word})
                    prepared={'models':[{'name':'factor_free','mapping':mapping}]}
                else:
                    anchor=linear_combination(curve,basis,word)
                    prepared=prepare(seed['curve'],list(map(str,anchor)),mapper,
                                     plan['fit']['local_weight'] if arm=='adaptive_local' else 0.)
                    anchors[ci]=AnchorExposure(prepared,plan['fit'],ci)
                maps[ci]=prepared
                path=dest/f'prepared-{ci:03d}.json.gz'
                with gzip.GzipFile(str(path),'wb',mtime=0) as f:f.write(canonical(prepared))
                event={'decision':decision,'status':'PREPARED','sha256':digest(path.read_bytes())}
            except (TimeoutError,ArithmeticError,ValueError,AssertionError) as e:
                event={'decision':decision,'status':'PREPARATION_UNKNOWN','error':type(e).__name__+': '+str(e)}
            finally:signal.alarm(0)
            components['preparation']+=cpu()-before
            event['cpu_seconds']=cpu()-before;events.append(event);next_index+=1
            atomic(dest/'events.json',events)
            if arm!='factor_free' or event['status']!='PREPARED':continue
            decision={'kind':'search','index':ci,'model':0,'height':125000}
        ci,mi,height=decision['index'],decision['model'],decision['height']
        mapping=maps[ci]['models'][mi]['mapping'];word=case['centres'][ci]
        remaining=plan['search_cpu_seconds']-(cpu()-start)
        if remaining<=0:status='CPU_CAP_AFTER_PREPARATION';break
        before=cpu()
        search=PointedQuarticSearch(state=state,centre={'coefficients':word},coordinate_policy=mapping['coordinate_policy'])
        transcript,points=execute(search,mapping,height,min(plan['point_wall_seconds'],remaining),plan['gp_sha256'])
        elapsed=cpu()-before;components['search']+=elapsed
        complete=transcript['status']=='bounded_search_complete'
        before=cpu(); new_point=None
        for point in points:
            if admission.consider(point)['status']=='INDEPENDENT_FINITE_COLUMN':
                new_point=list(map(str,point));success=True;break
        components['admission']+=cpu()-before
        row={'decision':decision,'centre':word,'mapping':mapping,'search':transcript,
             'new_point':new_point,'call_cpu_seconds':elapsed}
        path=dest/f'call-{len(calls):03d}.json';write(path,row)
        calls.append({'file':path.name,'sha256':digest(path.read_bytes()),'decision':decision,
                      'status':transcript['status'],'cpu_seconds':elapsed})
        events.append({'decision':decision,'status':transcript['status'],'call':path.name,'cpu_seconds':elapsed})
        if arm!='factor_free':
            anchors[ci].observe(decision['job'],complete)
            cost.observe(height,elapsed,complete)
        atomic(dest/'events.json',events)
        if success:status='PENDING_INDEPENDENT_CERTIFICATION';break
    if not success and len(calls)>=plan['maximum_calls']:
        status='CALL_CAP'
    search_end=cpu()
    if success:
        proof=checked_rank(curve,admission.points,admission.primes,seed['proof']['no_rational_2_torsion_prime'])
        write(dest/'rank-input.json',{'curve':seed['curve'],'points':[list(map(str,p)) for p in admission.points],'proof':proof})
    result={'status':status,'case':case_id,'arm':arm,'success':False,'candidate':success,
        'development_only':development_context is not None,
        'family':case['family'],'stratum':case['stratum'],'initial_rank':len(basis),
        'rank_lower_bound':len(basis),'calls':calls,'events_sha256':digest((dest/'events.json').read_bytes()),
        'protocol_sha256':ph, 'components':components, 'search_phase_cpu_seconds':search_end-start}
    # Every arm pays for transcript and policy replay, including failed arms.
    # Independent Sage group arithmetic is inside the outer timing boundary.
    from verify_cancellation_scheduler import verify_arm
    verification=verify_arm(case,arm,dest,plan,result)
    write(dest/'independent-verification.json',verification)
    if verification['status']!='PASS':raise ArithmeticError('independent verification failed')
    result.update(success=success, status='CERTIFIED_NEXT_DIRECTION' if success else status,
        rank_lower_bound=len(admission.points), independent_verification_sha256=digest((dest/'independent-verification.json').read_bytes()),
        proof_and_independent_replay_cpu_seconds=cpu()-search_end, total_internal_cpu_seconds=cpu()-start)
    write(dest/'result.json',result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('calls','components')}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',required=True);parser.add_argument('--arm',required=True)
    args=parser.parse_args();run(args.case,args.arm)
