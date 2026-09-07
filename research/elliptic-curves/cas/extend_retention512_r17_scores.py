#!/usr/bin/env python3
"""Checkpointed extended scoring of4608 newly retained fixed-box candidates."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
import benchmark_r17_extended_prime_traces as engine
import extend_retained_r17_prime_scores as scoring
import retain512_compact_r17 as retention
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=retention.ROOT;D=retention.LOCAL/'r17-retention512-extended-v1';OLD=(retention.LOCAL/'r17-retained-extended-primes-v1/result.json',retention.LOCAL/'r17-discarded-shards-extended-v1/result.json')
def sources():
    return {**engine.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(scoring.__file__).resolve(),Path(retention.__file__).resolve())}}
def choose(rows):
    order=sorted(rows,key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],r['numerator']))
    return {'extended_top_four':[r['retained_index'] for r in order[:4]]}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve wider trace protocol')
    retained=cert.read(retention.D/'result.json');replay=cert.read(retention.D/'replay.supervisor.json')
    if retained['status']!='COMPLETE_FIXED_RETENTION' or replay['outcome']!='completed' or replay['returncode']!=0:raise ArithmeticError('wider retention replay incomplete')
    excluded={(r['family'],r['parameter']) for path in OLD for r in cert.read(path)['rows']}
    if len(excluded)!=1536:raise ArithmeticError('previous scored population differs')
    families={f['family']:f for f in cert.read(spec.ATLAS)['families']};rows=[]
    for family in sorted(families):
        pool=[r for shard in retained['shards'] if shard['family']==family for r in shard['rows']]
        if len(pool)!=1024 or len({r['parameter'] for r in pool})!=1024 or len([r for r in pool if (family,r['parameter']) in excluded])!=256:raise ArithmeticError('signed retained1024/old256 inclusion differs')
        remaining=[r for r in pool if (family,r['parameter']) not in excluded];remaining.sort(key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],r['numerator']))
        for index,r in enumerate(remaining):
            model,points=spec.specialize(families[family],r['parameter'])
            if any(q.denominator!=1 for q in model):raise ArithmeticError('integral specialization required')
            rows.append({**r,'family':family,'retained_index':index,'model':list(map(str,model))})
    if len(rows)!=4608:raise ArithmeticError('fixed4608 extension roster differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.r17-retention512-extended.v1','sources':sources(),'rows':rows,'retention_sha256':cert.hashed(retention.D/'result.json'),'retention_replay_sha256':cert.hashed(retention.D/'replay.supervisor.json'),'old_scored_hashes':{str(p.relative_to(ROOT)):cert.hashed(p) for p in OLD},'gp_sha256':cert.hashed(engine.GP),'prime_roster':engine.PRIMES,'selection_extension':[4099,32749],'disjoint_validation':[32771,65521],'per_curve_seconds':5,'outer_seconds':1500,'rss_bytes_per_curve':536870912,'maximum_workers':2,'checkpoint_block_size':16,'score':'Original quantized562-prime score plus the same quantized extension4099..32749. Validation32771..65521 does not enter selection. Ties use combined good-prime count, denominator, signed numerator.','gate':'The exact512-per-sign retention campaign preserves every old128 prefix and independently checks all6144 returned short-prime scores. Exclude all1536 already extended addresses, leaving4608 fixed addresses. Broader retention addresses an experimentally demonstrated cutoff loss; no rank theorem is inferred.','selection':'Exactly four per family by the frozen extended selection score after all4608 rows complete and replay. No public points, catalogue labels, known-record parameters or validation scores enter ordering. No point search or automatic next stage here.','point_scope':'At most24 new point attempts after a separate protocol, with all43/49 exact maximum generic parity classes and uniformly125000 height under the prior completed control gate.','checkpoints':'Every raw invocation is immutable and checkpointed. Aggregate rows commit in16-address blocks. Explicit resume validates all committed rows and reuses successful raw records; failed raw invocations are never silently retried or overwritten.','boundaries':'This is still a short-prime-truncated fixed H4096 population. No global extended-score optimum, exact rank, upper rank, point absence or universal novelty is asserted.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(engine.GP) or p['retention_sha256']!=cert.hashed(retention.D/'result.json') or p['retention_replay_sha256']!=cert.hashed(retention.D/'replay.supervisor.json'):raise ArithmeticError('frozen trace inputs differ')
    for name,h in p['old_scored_hashes'].items():
        if cert.hashed(ROOT/name)!=h:raise ArithmeticError('prior trace population changed')
    return p
def raw_path(r):return D/r['family']/f"candidate-{r['retained_index']:04}"/'raw.json'
def checked_row(r,p,create):
    path=raw_path(r);command=engine.program(r['model'])
    if not path.exists():
        if not create:raise ArithmeticError('missing immutable raw trace')
        cap=capture([str(engine.GP),'-q','-s','256000000'],input_text=command,limits=Limits(p['per_curve_seconds'],p['rss_bytes_per_curve']),log_path=path.parent/'gp.log',separate_stderr=True,check=False);checkpoint(path,{'program':command,'stdout':cap.stdout,'stderr':cap.stderr,'supervision':cap.supervision})
    raw=cert.read(path)
    if raw['program']!=command or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('failed/censored raw trace retained; no automatic retry')
    traces,ms=engine.parse(raw['stdout'],r['model']);scores=scoring.sums(traces)
    return {**r,**scores,'combined_selection_units':r['score_units']+scores['extension_selection_units'],'combined_good':r['good_primes']+scores['extension_good'],'trace_cpu_ms':ms,'raw_sha256':cert.hashed(path)}
def run(resume=False):
    p=protocol();out=D/'result.json'
    if out.exists():
        if not resume:raise FileExistsError('use explicit resume for an existing trace checkpoint')
        data=cert.read(out)
        if data['protocol_hash']!=digest(p) or data['status']!='RUNNING':raise ArithmeticError('only a valid running checkpoint can resume')
        for original,r in zip(p['rows'],data['rows']):
            if checked_row(original,p,False)!=r:raise ArithmeticError('committed trace row changed')
    else:
        if resume:raise FileNotFoundError('no checkpoint to resume')
        data={'status':'RUNNING','protocol_hash':digest(p),'rows':[]};checkpoint(out,data)
    with ThreadPoolExecutor(max_workers=2) as pool:
        for start in range(len(data['rows']),len(p['rows']),p['checkpoint_block_size']):
            block=p['rows'][start:start+p['checkpoint_block_size']];rows=list(pool.map(lambda r:checked_row(r,p,True),block));data['rows'].extend(rows);checkpoint(out,data)
            if len(data['rows'])%128==0:print('RETENTION512 EXTENDED',len(data['rows']),'of4608',flush=True)
    data['selection']={f:choose([r for r in data['rows'] if r['family']==f]) for f in sorted({r['family'] for r in p['rows']})};data['status']='COMPLETE_FROZEN_TRACE_EXTENSION';checkpoint(out,data)
def replay():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='COMPLETE_FROZEN_TRACE_EXTENSION' or d['protocol_hash']!=digest(p) or len(d['rows'])!=4608:raise ArithmeticError('complete4608 traces required')
    for original,r in zip(p['rows'],d['rows']):
        if checked_row(original,p,False)!=r:raise ArithmeticError('trace row replay differs')
    selection={f:choose([r for r in d['rows'] if r['family']==f]) for f in sorted({r['family'] for r in p['rows']})}
    if selection!=d['selection']:raise ArithmeticError('frozen top-four selection differs')
    print('REPLAYED4608 TRACE ROSTERS AND FIXED24 SELECTION',selection,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','resume','replay']);a=p.parse_args();run(True) if a.stage=='resume' else globals()[a.stage]()
