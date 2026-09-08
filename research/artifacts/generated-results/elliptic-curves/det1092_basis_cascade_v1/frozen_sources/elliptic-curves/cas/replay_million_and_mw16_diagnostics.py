#!/usr/bin/env python3
"""Read-only arithmetic replays, using temporary destinations for frozen emitters."""
import argparse,tempfile
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import digest
import benchmark_native29_million_chart as native
import benchmark_pointed_sieve_v2 as sieve
import audit_native29_cvp_neighbours as neighbours
import benchmark_mw16_extended_prime_traces as traces

def emitted(module):
    original=module.OUT;expected=cert.read(original)
    with tempfile.TemporaryDirectory(prefix='elliptic-replay-') as name:
        module.OUT=Path(name)/'replay.json'
        try:
            module.replay()
            if cert.read(module.OUT)!=expected:raise ArithmeticError('frozen emitted certificate differs')
        finally:module.OUT=original

def neighbour_replay():
    expected=cert.read(neighbours.ART/'native29_cvp_neighbours_replay_v1.json')
    with tempfile.TemporaryDirectory(prefix='elliptic-neighbours-') as name:
        path=Path(name)/'replay.json';neighbours.replay(path)
        if cert.read(path)!=expected:raise ArithmeticError('frozen neighbour certificate differs')

def trace_replay():
    p=cert.read(traces.D/'protocol.json');d=cert.read(traces.D/'result.json');old=cert.read(traces.D/'replay.json')
    if p['sources']!=traces.sources() or d['status']!='PASS' or d['protocol_hash']!=digest(p) or len(d['rows'])!=5:raise ArithmeticError('five trace benchmark bindings differ')
    atlas={r['fibration_id']:r for r in cert.read(traces.spec.ATLAS)['families']};count=0
    for source,row in zip(p['rows'],d['rows']):
        pop=traces.PARENT/source['family']/'population.json'
        if cert.hashed(pop)!=source['population_sha256'] or cert.read(pop)['retained_candidates'][4]['parameter']!=source['parameter']:raise ArithmeticError('fifth saved address differs')
        model,_=traces.spec.specialize(atlas[source['family']],source['parameter'])
        if list(map(str,model))!=source['model']:raise ArithmeticError('MW16 benchmark model differs')
        rawpath=traces.D/source['family']/'raw.json';raw=cert.read(rawpath)
        if raw['program']!=traces.program(source['model']) or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('trace execution failed or differs')
        values,ms=traces.parse(raw['stdout'],source['model']);checks=[]
        for i in p['direct_validation_indices']:
            q,t=values[i];expected=traces.direct(source['model'],q)
            if t!=expected:raise ArithmeticError('independent character sum differs')
            checks.append([q,expected]);count+=1
        expected={**source,'status':'PASS','trace_count':len(values),'cpu_ms':ms,'wall_seconds':raw['supervision']['wall_seconds'],'direct_checks':checks,'raw_sha256':cert.hashed(rawpath)}
        if row!=expected:raise ArithmeticError('trace benchmark row differs')
    if count!=40 or old['status']!='PASS' or old['exact_character_sums']!=40 or old['protocol_sha256']!=cert.hashed(traces.D/'protocol.json') or old['result_sha256']!=cert.hashed(traces.D/'result.json'):raise ArithmeticError('recorded40-check replay differs')
    print('REPLAYED FIVE MW16 TRACE ROSTERS AND40 EXACT CHARACTER SUMS',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['native','sieve','neighbours','traces']);a=p.parse_args();{'native':lambda:emitted(native),'sieve':lambda:emitted(sieve),'neighbours':neighbour_replay,'traces':trace_replay}[a.stage]()
