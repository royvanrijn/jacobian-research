#!/usr/bin/env python3
"""Original short-score cache agreement on all967 saved11952 outer rows."""
import argparse,json,struct,sys
from array import array
from math import log
from pathlib import Path
import certify_compact_r17_candidates as cert
import build_extended_projective_trace_cache_11952 as parent
import extend_outer131072_r17 as extension
import benchmark_periodic_nagao_scanner as legacy
from mod2_reduction_independence import _primes_up_to
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=parent.ROOT;CAS=parent.CAS;ART=parent.ART;D=ROOT/'artifacts/local/elliptic-curves/retained-short-cache-benchmark-v1'
CPP=CAS/'newfamily/score_retained_projective_cache_v2.cpp';BINARY=D/'scorer';TABLE=D/'selection-cache.bin';INPUT=D/'candidates.txt';OUT=ART/'retained_short_cache_benchmark_v1.json'

def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),CPP,ROOT/'elliptic-curves/tests/test_retained_projective_cache_scorer_v2.py',parent.OUT,parent.D/'protocol.json',parent.D/'check.supervisor.json',extension.D/'result.json',ART/'retained_extended_cache_benchmark_v1.json',CAS/'extend_retained_r17_prime_scores.py',Path(legacy.__file__).resolve()]}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve short-cache benchmark')
    completed=cert.read(ART/'retained_extended_cache_benchmark_v1.json')
    if completed['status']!='PASS' or completed['candidates_checked']!=967:raise ArithmeticError('exact extended cached-score gate required')
    old=cert.read(extension.D/'result.json');rows=[r for r in old['rows'] if r['family']=='11952']
    primes=[p for p in _primes_up_to(4093) if p>=5];tables=[]
    f=next(r for r in cert.read(legacy.spec.ATLAS)['families'] if r['family']=='11952')
    model={k:f[k]+['0']*(n-len(f[k])) for k,n in [('A_coefficients_low_to_high',9),('B_coefficients_low_to_high',13)]}
    for prime in primes:
        path=legacy.old.PARENT/'11952/trace-tables'/f'{prime}.json';tables.append({'prime':prime,'table_sha256':cert.hashed(path)})
    if len(rows)!=967 or len(primes)!=562:raise ArithmeticError('fixed967/562 score comparison required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.retained-short-cache-benchmark.v1','sources':sources(),'rows':[{'retained_index':r['retained_index'],'numerator':r['numerator'],'denominator':r['denominator'],'expected_units':r['score_units'],'expected_good':r['good_primes']} for r in rows],'prime_rows':tables,'model_hash':digest(model),'compiler':['g++','-O3','-std=c++17'],'compile_seconds':60,'encode_seconds':60,'replay_seconds':60,'score_seconds':60,'rss_bytes':2147483648,'quantization':'Preserve the original short expression round((2-a_p)/(p+1-a_p)*log(p)*10**12), including Python operation order. Singular rows have zero units and good=0.','scope':'Encode the existing562-prime11952 short cache and compare every byte and all967 saved short sums with a compiled retained-list scorer. Version two permits primes5through32749 with unchanged integer lookup logic. This enables exact returned-short-score replay on a larger retained set without changing the original scanner or its retention policy. No parameter enumeration, ranking change, new candidate, point search or automatic wider campaign.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen cached-score benchmark changed')
    return p

def frame(prime,p):
    row=next(r for r in p['prime_rows'] if r['prime']==prime);path=legacy.old.PARENT/'11952/trace-tables'/f'{prime}.json'
    if cert.hashed(path)!=row['table_sha256']:raise ArithmeticError('certified projective table changed')
    data=cert.read(path)
    if data['input']!={'family':'11952','model_hash':p['model_hash'],'prime':prime} or len(data['traces'])!=prime+1 or len(data['good'])!=prime+1:raise ArithmeticError('projective cache shape differs')
    units=array('q',(round((2-t)/(prime+1-t)*log(prime)*10**12) if good else 0 for t,good in zip(data['traces'],data['good'])))
    if units.itemsize!=8:raise ArithmeticError('64-bit integer array required')
    if sys.byteorder!='little':units.byteswap()
    return struct.pack('<II',prime,prime+1)+units.tobytes()+bytes(data['good'])

def encode(check=False):
    p=protocol();header=b'R17XS001'+struct.pack('<I',len(p['prime_rows']));tail=b'ENDXSC01'
    if not check and TABLE.exists():raise FileExistsError('preserve encoded selection cache')
    with TABLE.open('rb' if check else 'xb') as out:
        def consume(data):
            if check:
                if out.read(len(data))!=data:raise ArithmeticError('encoded cache bytes differ')
            else:out.write(data)
        consume(header)
        for row in p['prime_rows']:consume(frame(row['prime'],p))
        consume(tail)
        if check and out.read(1):raise ArithmeticError('trailing cache bytes')
    text='R17-CANDIDATES-V1 '+str(len(p['rows']))+'\n'+''.join(f"{r['numerator']} {r['denominator']}\n" for r in p['rows'])
    if check:
        if INPUT.read_text()!=text:raise ArithmeticError('candidate bytes differ')
    else:
        if INPUT.exists():raise FileExistsError('preserve candidate fixture')
        INPUT.write_text(text)
    print('EXACT ENCODED SHORT CACHE',TABLE.stat().st_size,'BYTES',flush=True)

def compile():
    p=protocol()
    if BINARY.exists():raise FileExistsError('preserve scorer binary')
    cmd=[*p['compiler'],str(CPP),'-o',str(BINARY)];s=capture(cmd,limits=Limits(p['compile_seconds'],p['rss_bytes']),log_path=D/'compile.log',separate_stderr=True,check=False)
    if s.supervision['outcome']!='completed' or s.supervision['returncode']!=0:raise ArithmeticError('scorer compilation failed')
    checkpoint(D/'build.json',{'command':cmd,'binary_sha256':cert.hashed(BINARY),'supervision':s.supervision,'stderr':s.stderr})

def score(check=False):
    p=protocol();build=cert.read(D/'build.json')
    if build['binary_sha256']!=cert.hashed(BINARY):raise ArithmeticError('cached scorer binary changed')
    path=D/'raw.json';cmd=[str(BINARY),str(TABLE),str(INPUT)]
    if not check:
        if path.exists() or OUT.exists():raise FileExistsError('preserve one cached score benchmark')
        s=capture(cmd,limits=Limits(p['score_seconds'],p['rss_bytes']),log_path=D/'score.log',separate_stderr=True,check=False)
        checkpoint(path,{'command':cmd,'stdout':s.stdout,'stderr':s.stderr,'supervision':s.supervision})
    raw=cert.read(path);lines=raw['stdout'].splitlines()
    if raw['command']!=cmd or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0 or len(lines)!=len(p['rows'])+1 or lines[-1]!=f"S {len(p['rows'])} {len(p['prime_rows'])}":raise ArithmeticError('cached score invocation/frame differs')
    for i,(row,line) in enumerate(zip(p['rows'],lines[:-1])):
        if line!=f"R {i} {row['expected_units']} {row['expected_good']}":raise ArithmeticError('cached sum differs from scalar GP proof')
    result={'schema':'elliptic-curves.retained-short-cache-benchmark-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'binary_sha256':cert.hashed(BINARY),'cache_sha256':cert.hashed(TABLE),'cache_bytes':TABLE.stat().st_size,'candidate_sha256':cert.hashed(INPUT),'raw_sha256':cert.hashed(path),'candidates_checked':len(p['rows']),'selection_primes':len(p['prime_rows']),'lookup_count':len(p['rows'])*len(p['prime_rows']),'wall_seconds':raw['supervision']['wall_seconds'],'claim_boundary':p['scope']}
    if check:
        if cert.read(OUT)!=result:raise ArithmeticError('cached-score benchmark report differs')
    else:checkpoint(OUT,result)
    print('CACHED SHORT SUMS MATCH ALL967 ORIGINAL SCORES',raw['supervision']['wall_seconds'],'SECONDS',flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','encode','encoding-check','compile','score','score-check']);v=a.parse_args()
    if v.stage=='encoding-check':encode(True)
    elif v.stage=='score-check':score(True)
    else:globals()[v.stage]()
