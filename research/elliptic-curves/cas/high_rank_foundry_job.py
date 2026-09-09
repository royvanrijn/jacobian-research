"""One recoverable, bounded job. Only its sealed result can advance the ledger."""
import fcntl
from fractions import Fraction as F
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback

from research_runtime.supervisor import Limits, run as supervise
from v3_warm_support import read, sha, require, atomic as _atomic, process_info
from high_rank_foundry_intake import conductor_gate, novelty

CAS=Path(__file__).resolve().parent
ROOT=CAS.parents[1]


def save(path,data,immutable=False):
    _atomic(path,json.loads(json.dumps(data)),immutable=immutable)


class PhaseFailure(RuntimeError):
    def __init__(self,phase,record):
        self.phase,self.record=phase,record
        super().__init__(phase+': '+record['outcome'])


def step(job,phase,index=0,repeat=0,seconds=None,allow_resource=False):
    config=read(job/'request.json')['config']
    if phase in ('search','replay'):
        suffix='-supervision' if phase=='search' else '-replay-supervision'
        directory=job/(f'search-{index:02d}'+suffix)
    else:
        directory=job/'phases'/f'{phase}-{index:02d}-{repeat}'
    receipt=directory/'supervisor.json'
    if receipt.exists():
        old=read(receipt)
        if old['outcome']=='completed' and old['returncode']==0:
            return old
        if allow_resource and old['outcome'] in ('strict_wall_timeout','strict_rss_limit'):
            return old
        # Crash recovery keeps all prior supervisor logs. Arithmetic checkpoints
        # are replayed by the same frozen worker, never replaced by a new search.
        retained=directory.with_name(directory.name+'-attempt-'+str(time.time_ns()))
        directory.rename(retained)
    directory.mkdir(parents=True,exist_ok=True)
    save(job/'progress.json',{'phase':phase,'index':index,'updated_at':time.time()})
    seconds=seconds or (config['phase_seconds'] if phase in ('search','replay','seed-search','seed-replay') else 600)
    cmd=[config['sage'],'-python',str(CAS/'high_rank_foundry_arithmetic.py'),phase,'--job',str(job),'--index',str(index)]
    record=supervise(cmd,limits=Limits(wall_seconds=seconds,rss_bytes=config['rss_bytes']),
                     log_path=directory/'worker.log',checkpoint_path=receipt,cwd=ROOT,
                     env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    if record['outcome']=='completed' and record['returncode']==0:
        return record
    if allow_resource and record['outcome'] in ('strict_wall_timeout','strict_rss_limit'):
        return record
    raise PhaseFailure(phase,record)


def segment_events(terminal, reconciled, inherited, used):
    events=[]; offsets={}; subtotal=0
    for s in terminal['stages']:
        offsets[s['epoch']]=subtotal;subtotal+=s['charts']
        if s['after']>s['before'] and subtotal>inherited:
            events.append({'phase':'amplify','call':used+subtotal-inherited,'before':s['before'],'after':s['after']})
    for g in reconciled['gains']:
        path=Path(g['chart']);epoch=int(path.parts[0].split('-')[1]);i=int(path.stem.split('-')[1])
        origin=offsets[epoch]+i+1
        if origin<=inherited:
            raise ArithmeticError('a previously reconciled prefix hid an additional direction')
        events.append({'phase':'cloud','call':used+origin-inherited,'before':g['after']-1,'after':g['after']})
    return events


def conductor_job(job,request):
    packet=read(ROOT/request['candidate']['packet'])
    require(sha(ROOT/request['candidate']['packet'])==request['candidate']['packet_sha256'],'conductor input changed')
    model=list(map(F,packet['curve']));d=math.lcm(*(a.denominator for a in model))
    source={'id':request['candidate']['id'],'original_curve':packet['curve'],
            'curve':[str(a*d**w) for a,w in zip(model,(1,2,3,4,6))],
            'integral_model_scale':str(d),'rank_lower_bound':packet['rank_lower_bound'],
            'known_primes':[],'factor_hints':[]}
    save(job/'conductor-input.json',source,immutable=True)
    build=step(job,'conductor-build',seconds=request['config']['conductor_seconds'],allow_resource=True)
    if not (job/'conductor.json').exists():
        raise PhaseFailure('conductor-build-no-checkpoint',build)
    step(job,'conductor-replay',seconds=120)
    result=read(job/'conductor.json');threshold=request['conductor_gate']['threshold']
    comparison=('BEATS_PINNED_REPORTED_THRESHOLD' if result['status']=='EXACT' and threshold and int(result['exact_conductor'])<int(threshold)
                else 'ABOVE_PINNED_THRESHOLD' if threshold and int(result['conductor_divisor'])>=int(threshold) else 'UNKNOWN')
    return {'status':'PASS_CONDUCTOR_REPLAY','conductor':result,'comparison':comparison,
            'gate':request['conductor_gate'],'calls':0,'rank_lower_bound':packet['rank_lower_bound'],
            'claim_boundary':'Exact conductor only if factorization and local proofs close; benchmark is a dated catalogue, not a world-record assertion.'}


def search_job(job,request):
    row=request['candidate'];config=request['config'];case=row['id']
    folder=job/'cases'/case
    used=0;events=[];segments=[];timeouts=0;map_timeouts=0
    head=None;packet_path=ROOT/row['packet'] if row.get('packet') else None
    if request['kind']=='fresh' or (request.get('revival') and row['rank']==17):
        step(job,'generic')
        gate=read(job/'generic-seed-gate.json')
        if gate['status']=='UNRESOLVED_GENERIC_SEED':
            require(request['kind']=='fresh' and not row.get('packet'),
                    'a previously certified seed failed its input gate')
            return {'status':'UNRESOLVED_GENERIC_SEED','id':case,'family':row['family'],
                    'parameter':row['parameter'],'rank_lower_bound':None,'calls':0,
                    'point_invocation_intents':0,'completed_point_invocations':0,
                    'seed_gate':gate,'seed_gate_sha256':sha(job/'generic-seed-gate.json'),
                    'claim_boundary':gate['claim_boundary']}
        require(gate['status']=='PASS_GENERIC_SEED_GATE','unexpected generic seed gate result')
        for phase in ('seed-prepare','seed-search','seed-replay','seed-cloud'):
            step(job,phase)
        step(job,'seed-cloud',repeat=1)
        packet_path=folder/'seed-reconciled.json';packet=read(packet_path)
        seed=read(folder/'seed-search/terminal.json');used=seed['charts']
        timeouts+=seed['point_timeouts']
        map_timeouts+=sum(a['skip_reason'] in ('strict_wall_timeout','strict_rss_limit') for a in seed['map_attempts'])
        if seed['rank_lower_bound']>17:
            events.append({'phase':'seed','call':used,'before':17,'after':18})
        for g in packet['reconciliation_gains']:
            events.append({'phase':'seed-cloud','call':g['seed_call'],'before':g['after']-1,'after':g['after']})
        seed_calls=used
    else:
        require(sha(packet_path)==row['packet_sha256'],'continuation packet changed')
        packet=read(packet_path);seed_calls=0
    unresolved=False
    if 17<packet['rank_lower_bound']<32:
        remaining=request['allowance']
        for index in range(16):
            if remaining<=0: break
            dest=job/f'search-{index:02d}'
            parent=None
            if index==0 and request.get('use_cached'):
                parent=row['head']
            spec={'packet':str(packet_path.relative_to(ROOT)),'allowance':remaining,'parent':parent}
            save(job/f'segment-{index:02d}.json',spec,immutable=True)
            if parent is None: step(job,'bank',index)
            for phase in ('search','replay','cloud'): step(job,phase,index)
            step(job,'cloud',index,repeat=1)
            terminal=read(dest/'terminal.json');reconciled=read(job/f'reconciled-{index:02d}.json')
            inherited=read(dest/'protocol.json').get('inherited_charts',0)
            calls=terminal['charts']-inherited
            require(0<=calls<=remaining,'point call allowance exceeded')
            events.extend(segment_events(terminal,reconciled,inherited,used))
            new_stages=terminal['stages']
            timeouts+=sum(s['censored'] for s in new_stages)
            map_timeouts+=sum(s['censored_maps'] for s in new_stages)
            if inherited:
                old=read(ROOT/parent/'terminal.json')
                timeouts-=sum(s['censored'] for s in old['stages'])
                map_timeouts-=sum(s['censored_maps'] for s in old['stages'])
            used+=calls;remaining-=calls
            segments.append({'folder':str(dest.relative_to(ROOT)),'calls':calls,'inherited_calls':inherited,
                             'rank':reconciled['rank_lower_bound'],'stop':terminal['stop_reason']})
            packet_path=job/f'reconciled-{index:02d}.json';packet=reconciled
            # Audited odd-prime evidence cannot be silently dismissed if mod2
            # admission has not yet resolved it. Preserve the certified packet.
            strongest=max([terminal['rank_lower_bound']]+[r for s in new_stages for r in s['odd_ranks'].values()])
            if packet['rank_lower_bound']<strongest:
                unresolved=True;break
            if packet['rank_lower_bound']>=32: break
            if terminal['stop_reason']=='ADDITIONAL_FINITE_RANK_REQUIRES_RECONCILIATION':
                if packet['rank_lower_bound']<=terminal['rank_lower_bound']:
                    unresolved=True;break
                continue
            if terminal['stop_reason']=='CHART_BUDGET_EXHAUSTED' and packet['points']==terminal['points']:
                head=str(dest.relative_to(ROOT))
            break
    require(used<=request['allowance']+98,'overall job point allowance exceeded')
    from high_rank_foundry_certificate import compact_packet
    packet=compact_packet(packet)
    save(job/'packet.json',packet,immutable=True)
    step(job,'verify-packet')
    catalogue=read(ROOT/request['catalogue'])
    result={'status':'PASS_CERTIFIED_SEARCH','id':case,'family':row['family'],'parameter':row['parameter'],
            'rank_lower_bound':packet['rank_lower_bound'],'calls':used,'seed_calls':seed_calls,
            'point_invocation_intents':len(list((job/'point-invocations').glob('*/attempt-*.json'))),
            'completed_point_invocations':len(list((job/'point-invocations').glob('*/result.json'))),
            'gain_timeline':events,'segments':segments,'point_timeouts':timeouts,'map_timeouts':map_timeouts,
            'unresolved_cloud':unresolved,'bank_index':request['bank_index'],'head':head,
            'packet_path':str((job/'packet.json').relative_to(ROOT)),'packet_sha256':sha(job/'packet.json'),
            'novelty':novelty(packet['curve'],catalogue),
            'conductor_gate':conductor_gate(packet['curve'],packet['rank_lower_bound'],catalogue),
            'claim_boundary':'Certified subgroup lower bound, never an exact rank. Every returned cloud reconciled. Bounded misses and unresolved conductors remain UNKNOWN.'}
    return result


def run(job):
    job.mkdir(parents=True,exist_ok=True)
    with (job/'job.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        save(job/'lease.json',{'pid':os.getpid(),'token':process_info(os.getpid())['start_token'],'started_at':time.time()})
        if (job/'result.json').exists(): return
        started=time.time();request=read(job/'request.json')
        try:
            result=conductor_job(job,request) if request['kind']=='conductor' else search_job(job,request)
        except Exception as e:
            result={'status':'UNRESOLVED_JOB_FAILURE','error':traceback.format_exc(),'calls':None,
                    'phase':getattr(e,'phase',None),'outcome':getattr(e,'record',{}).get('outcome'),
                    'claim_boundary':'No rank upper bound or no-point conclusion. Raw checkpoints preserved.'}
        result.update(request_sha256=sha(job/'request.json'),wall_seconds=time.time()-started,
                      evidence={str(p.relative_to(job)):sha(p) for p in sorted(job.rglob('*.json'))
                                if p not in (job/'progress.json',job/'result.json',job/'lease.json')})
        save(job/'result.json',result,immutable=True)
        print(json.dumps({k:result.get(k) for k in ('status','rank_lower_bound','calls','wall_seconds')}),flush=True)


if __name__=='__main__':
    run(Path(sys.argv[1]).resolve())
