#!/usr/bin/env python3
"""Binary-cache replay and exact score agreement on all967 saved11952 outer rows."""
import argparse,json,struct,sys
from array import array
from math import log
from pathlib import Path
import certify_compact_r17_candidates as cert
import build_extended_projective_trace_cache_11952 as parent
import extend_outer131072_r17 as extension
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=parent.ROOT;CAS=parent.CAS;ART=parent.ART;D=ROOT/'artifacts/local/elliptic-curves/retained-extended-cache-benchmark-v1'
CPP=CAS/'newfamily/score_retained_projective_cache.cpp';BINARY=D/'scorer';TABLE=D/'selection-cache.bin';INPUT=D/'candidates.txt';OUT=ART/'retained_extended_cache_benchmark_v1.json'

def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),CPP,ROOT/'elliptic-curves/tests/test_retained_projective_cache_scorer.py',parent.OUT,parent.D/'protocol.json',parent.D/'check.supervisor.json',extension.D/'result.json',CAS/'extend_retained_r17_prime_scores.py']}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve cached-score benchmark')
    cache=cert.read(parent.OUT);s=cert.read(parent.D/'check.supervisor.json');old=cert.read(extension.D/'result.json')
    if cache['status']!='PASS' or s['outcome']!='completed' or s['returncode']!=0 or old['status']!='COMPLETE_FROZEN_TRACE_EXTENSION':raise ArithmeticError('complete cache replay and old exact score cohort required')
    rows=[r for r in old['rows'] if r['family']=='11952']
    if len(rows)!=967 or cache['selection_primes']!=2948:raise ArithmeticError('fixed comparison population differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.retained-extended-cache-benchmark.v1','sources':sources(),'rows':[{'retained_index':r['retained_index'],'numerator':r['numerator'],'denominator':r['denominator'],'expected_units':r['extension_selection_units'],'expected_good':r['extension_good']} for r in rows],'prime_rows':cache['rows'],'compiler':['g++','-O3','-std=c++17'],'compile_seconds':60,'encode_seconds':600,'replay_seconds':600,'score_seconds':60,'rss_bytes':2147483648,'quantization':'Use exactly the existing extension expression round((2-a_p)*log(p)/(p+1-a_p)*10**12), preserving Python operation order; singular rows have zero units and good=0. No validation prime enters this file.','scope':'Encode the complete11952 selection-prime cache as a fixed little-endian integer table, independently replay its bytes against every source table and invoke the compiled retained-list scorer once on all967 previously scored outer candidates. Require every quantized sum and good-prime count to match its existing scalar-GP certificate. This benchmarks cached score lookup only; it performs no parameter enumeration, ranking change, point search, future survivor selection or automatic wider campaign.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen cached-score benchmark changed')
    return p

def frame(prime,p):
    row=next(r for r in p['prime_rows'] if r['prime']==prime);path=parent.D/str(prime)/'table.json'
    if cert.hashed(path)!=row['table_sha256']:raise ArithmeticError('certified projective table changed')
    data=cert.read(path)
    if data['input']['prime']!=prime or len(data['traces'])!=prime+1 or len(data['good'])!=prime+1:raise ArithmeticError('projective cache shape differs')
    units=array('q',(round((2-t)*log(prime)/(prime+1-t)*10**12) if good else 0 for t,good in zip(data['traces'],data['good'])))
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
    print('EXACT ENCODED EXTENDED CACHE',TABLE.stat().st_size,'BYTES',flush=True)

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
    result={'schema':'elliptic-curves.retained-extended-cache-benchmark-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'binary_sha256':cert.hashed(BINARY),'cache_sha256':cert.hashed(TABLE),'cache_bytes':TABLE.stat().st_size,'candidate_sha256':cert.hashed(INPUT),'raw_sha256':cert.hashed(path),'candidates_checked':len(p['rows']),'selection_primes':len(p['prime_rows']),'lookup_count':len(p['rows'])*len(p['prime_rows']),'wall_seconds':raw['supervision']['wall_seconds'],'claim_boundary':p['scope']}
    if check:
        if cert.read(OUT)!=result:raise ArithmeticError('cached-score benchmark report differs')
    else:checkpoint(OUT,result)
    print('CACHED SUMS MATCH ALL967 SCALAR GP SCORES',raw['supervision']['wall_seconds'],'SECONDS',flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','encode','encoding-check','compile','score','score-check']);v=a.parse_args()
    if v.stage=='encoding-check':encode(True)
    elif v.stage=='score-check':score(True)
    else:globals()[v.stage]()
