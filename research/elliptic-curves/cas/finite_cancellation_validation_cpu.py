#!/usr/bin/env sage -python
"""One isolated validation arm; worker inputs contain no endpoint or target."""
import argparse
from fractions import Fraction as F
import gzip
from importlib.machinery import SourceFileLoader
import json
from math import ceil
from pathlib import Path
import resource
import signal
import time

from finite_cancellation_corpus import ROOT, LOCAL, canonical, digest, write
from finite_cancellation_validation_audit import OUT

RAW=LOCAL/'finite-cancellation-validation-v2'


def cpu_now():
    r=resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time()+r.ru_utime+r.ru_stime


def run(case_id,arm):
    start=cpu_now();wall=time.monotonic()
    from finite_cancellation_features import alarm
    from finite_cancellation_validation_features import prepare,select
    from half_lattice_pointed_sieve import linear_combination
    from v3_warm_engine import certified_state
    from future_point_admission import FinitePointAdmission
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import execute
    from memory_rank_certificate import checked_rank
    plan=json.loads((OUT/'protocol.json').read_text());ph=digest((OUT/'protocol.json').read_bytes())
    assert arm in plan['arms'] and digest((OUT/'inputs.json').read_bytes())==plan['input_sha256']['inputs.json']
    for name,h in plan['source_sha256'].items():assert digest(Path(__file__).with_name(name).read_bytes())==h
    case=next(c for c in json.loads((OUT/'inputs.json').read_text()) if c['id']==case_id)
    dest=RAW/'arms'/case_id/arm;dest.mkdir(parents=True,exist_ok=True)
    if (dest/'result.json').exists():raise FileExistsError('Preserve the prior arm.')
    resource.setrlimit(resource.RLIMIT_CPU,(plan['hard_process_cpu_seconds']-5,plan['hard_process_cpu_seconds']))
    seed=case['seed'];curve=tuple(map(F,seed['curve']));basis=[tuple(map(F,p)) for p in seed['points']]
    state=certified_state(curve,basis,seed['proof'])
    admission=FinitePointAdmission(curve,basis,prime_bound=500)
    mapper=SourceFileLoader('validation_factor_free_mapper',str(Path(__file__).with_name('lean_factor_free_pari_mapping.sage'))).load_module()
    mapper.pari.allocatemem(256000000,silent=True)
    signal.signal(signal.SIGALRM,alarm)
    records=[];success=False;status='BANK_EXHAUSTED';prep_cpu=0.;point_cpu=0.;admission_cpu=0.;proof_cpu=0.
    initial_cpu=cpu_now()-start
    for ci,word in enumerate(case['centres']):
        remaining=plan['search_cpu_seconds_per_arm']-(cpu_now()-start)
        if remaining<=0:status='CPU_CAP';break
        before=cpu_now();signal.alarm(max(1,min(10,ceil(remaining))))
        row={'index':ci,'centre':word,'point_call':False}
        try:
            if arm=='factor_free':
                mapping=mapper.mapping(curve,basis,{'representative':word})
                selection={'selected':'factor_free','models_built':1}
            else:
                anchor=linear_combination(curve,basis,word)
                prepared=prepare(seed['curve'],list(map(str,anchor)),mapper,arm=='real_q_g')
                mi,scores=select(prepared,plan['fits'][arm]);mapping=prepared['models'][mi]['mapping']
                raw=canonical(prepared)
                path=dest/f'prepared-{ci:03d}.json.gz'
                with gzip.GzipFile(str(path),'wb',mtime=0) as f:f.write(raw)
                selection={'selected':prepared['models'][mi]['name'],'models_built':len(prepared['models']),
                    'scores':scores,'prepared_sha256':digest(raw),'prepared_file_sha256':digest(path.read_bytes())}
            search=PointedQuarticSearch(state=state,centre={'coefficients':word},coordinate_policy=mapping['coordinate_policy'])
            row.update(mapping=mapping,selection=selection)
        except TimeoutError:
            row['status']='MAP_TIMEOUT'
        except (ArithmeticError,ValueError,AssertionError) as e:
            row.update(status='MAP_FAILURE',error=type(e).__name__+': '+str(e))
        finally:signal.alarm(0)
        spent=cpu_now()-before;prep_cpu+=spent;row['preparation_cpu_seconds']=spent
        if row.get('status') not in ('MAP_TIMEOUT','MAP_FAILURE'):
            remaining=plan['search_cpu_seconds_per_arm']-(cpu_now()-start)
            if remaining<=0:
                row['status']='CPU_CAP_AFTER_PREPARATION';status=row['status']
            else:
                before=cpu_now()
                transcript,points=execute(search,mapping,plan['height'],min(plan['point_wall_seconds'],remaining),plan['gp_sha256'])
                point_cpu+=cpu_now()-before
                row.update(search=transcript,point_call=True,status=transcript['status'])
                before=cpu_now()
                for point in points:
                    if admission.consider(point)['status']=='INDEPENDENT_FINITE_COLUMN':
                        row['new_point']=list(map(str,point));success=True;break
                admission_cpu+=cpu_now()-before
        path=dest/f'chart-{ci:03d}.json'
        write(path,row)
        records.append({'index':ci,'status':row['status'],'point_call':row['point_call'],
            'selection':row.get('selection',{}).get('selected'),'new_direction_candidate':success,
            'cpu_after':cpu_now()-start,'row_sha256':digest(path.read_bytes())})
        if success:
            before=cpu_now()
            proof=checked_rank(curve,admission.points,admission.primes,seed['proof']['no_rational_2_torsion_prime'])
            write(dest/'rank-input.json',{'curve':seed['curve'],'points':[list(map(str,p)) for p in admission.points],'proof':proof})
            proof_cpu+=cpu_now()-before;status='CERTIFIED_NEXT_DIRECTION';break
        if row['status']=='CPU_CAP_AFTER_PREPARATION':break
    result={'status':status,'case':case_id,'family':case['family'],'stratum':case['stratum'],'arm':arm,
        'success':success,'initial_rank':len(basis),'rank_lower_bound':len(admission.points),'charts':records,
        'point_search_calls':sum(r['point_call'] for r in records),'cpu_seconds':cpu_now()-start,
        'wall_seconds':time.monotonic()-wall,'initial_verification_cpu_seconds':initial_cpu,
        'map_and_feature_cpu_seconds':prep_cpu,'point_backend_cpu_seconds':point_cpu,
        'admission_cpu_seconds':admission_cpu,'rank_proof_cpu_seconds':proof_cpu,'protocol_sha256':ph,
        'boundary':'Retained known-curve control with targets absent from worker input. Exact subgroup lower bounds only. All candidate construction and discarded features charged; common centre-bank creation and offline fitting reported separately.'}
    write(dest/'result.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='charts'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--case',required=True);p.add_argument('--arm',choices=['factor_free','real_q','real_q_g'],required=True)
    a=p.parse_args();run(a.case,a.arm)
