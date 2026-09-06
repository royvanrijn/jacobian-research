#!/usr/bin/env python3
"""One full-size repeated-fixture benchmark; no new candidate or point search."""
import argparse,json
from pathlib import Path
import certify_compact_r17_candidates as cert
import benchmark_retained_extended_cache as extended
import benchmark_retained_short_cache as short
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=extended.ROOT;ART=extended.ART;D=ROOT/'artifacts/local/elliptic-curves/million-retained-cache-benchmark-v1';OUT=ART/'million_retained_cache_benchmark_v1.json';INPUT=D/'repeated-candidates.txt'

def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),extended.OUT,short.OUT,extended.D/'protocol.json',short.D/'protocol.json',extended.BINARY,short.BINARY]}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve full-size performance fixture')
    a=cert.read(extended.OUT);b=cert.read(short.OUT);ep=extended.protocol();sp=short.protocol()
    if a['status']!='PASS' or b['status']!='PASS' or [(r['numerator'],r['denominator']) for r in ep['rows']]!=[(r['numerator'],r['denominator']) for r in sp['rows']]:raise ArithmeticError('same967 exact two-band fixtures required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.million-retained-cache-benchmark.v1','sources':sources(),'rows':1048576,'reference_rows':967,'pattern':'Repeat the967 immutable old candidate coordinates cyclically, preserving their existing exact short and extended expected scores.','bands':{'short':{'binary':str(short.BINARY.relative_to(ROOT)),'cache':str(short.TABLE.relative_to(ROOT)),'cache_sha256':b['cache_sha256'],'expected':sp['rows'],'primes':562},'extended':{'binary':str(extended.BINARY.relative_to(ROOT)),'cache':str(extended.TABLE.relative_to(ROOT)),'cache_sha256':a['cache_sha256'],'expected':ep['rows'],'primes':2948}},'seconds_per_band':300,'cost_gate_seconds_per_band':120,'rss_bytes':2147483648,'maximum_workers':1,'scope':'One full-size1048576-row repeated performance fixture and one cached invocation per score band, with exact comparison of every output row to its existing scalar-score certificate. These are duplicates for software cost testing, not new rational candidates, a parameter population, a rank experiment or a throughput guarantee on unseen inputs. No parameter scan or point search follows automatically.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen million-fixture sources changed')
    return p

def expected_input(p):
    yield f"R17-CANDIDATES-V1 {p['rows']}\n"
    rows=p['bands']['short']['expected']
    for i in range(p['rows']):
        r=rows[i%p['reference_rows']];yield f"{r['numerator']} {r['denominator']}\n"

def run(check=False):
    p=protocol()
    if check:
        with INPUT.open() as f:
            for expected in expected_input(p):
                if f.readline()!=expected:raise ArithmeticError('repeated input differs')
            if f.read(1):raise ArithmeticError('trailing repeated input')
    else:
        if INPUT.exists() or OUT.exists():raise FileExistsError('preserve one full-size benchmark')
        with INPUT.open('x') as f:f.writelines(expected_input(p))
    rows=[]
    for label,band in p['bands'].items():
        path=D/(label+'.json');cache=ROOT/band['cache'];cmd=[str(ROOT/band['binary']),str(cache),str(INPUT)]
        if cert.hashed(cache)!=band['cache_sha256']:raise ArithmeticError('encoded score cache changed')
        if not check:
            if path.exists():raise FileExistsError('preserve full-size cached invocation')
            c=capture(cmd,limits=Limits(p['seconds_per_band'],p['rss_bytes']),log_path=D/(label+'.log'),separate_stderr=True,check=False)
            checkpoint(path,{'command':cmd,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
        raw=cert.read(path);lines=raw['stdout'].splitlines()
        if raw['command']!=cmd or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0 or len(lines)!=p['rows']+1 or lines[-1]!=f"S {p['rows']} {band['primes']}":raise ArithmeticError('complete million-row invocation failed or censored')
        for i,line in enumerate(lines[:-1]):
            r=band['expected'][i%p['reference_rows']]
            if line!=f"R {i} {r['expected_units']} {r['expected_good']}":raise ArithmeticError('full-size cached output differs from exact reference')
        rows.append({'band':label,'rows':p['rows'],'primes':band['primes'],'wall_seconds':raw['supervision']['wall_seconds'],'raw_sha256':cert.hashed(path),'cost_gate_passed':raw['supervision']['wall_seconds']<=p['cost_gate_seconds_per_band']});print('FULL-SIZE CACHE BENCHMARK',label,rows[-1]['wall_seconds'],'SECONDS',flush=True)
    result={'schema':'elliptic-curves.million-retained-cache-benchmark-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'input_sha256':cert.hashed(INPUT),'rows':rows,'cost_gate_passed':all(r['cost_gate_passed'] for r in rows),'claim_boundary':p['scope']}
    if check:
        if cert.read(OUT)!=result:raise ArithmeticError('full-size benchmark report differs')
    else:checkpoint(OUT,result)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','run','check']);v=a.parse_args();prepare() if v.stage=='prepare' else run(v.stage=='check')
