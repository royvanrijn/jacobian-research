#!/usr/bin/env python3
"""Bounded trace-engine calibration on six prospective retained R17 addresses."""
import argparse
from pathlib import Path
from math import isqrt
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as spec
from mod2_reduction_independence import _is_prime
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=Path(__file__).resolve().parents[2];LOCAL=ROOT/'artifacts/local/elliptic-curves';D=LOCAL/'r17-extended-prime-benchmark-v1';PARENT=LOCAL/'compact-six-r17-h4096-v1';GP=Path('/usr/bin/gp')
PRIMES=[p for p in range(4099,65522) if _is_prime(p)]

def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(spec.__file__),spec.ATLAS,ROOT/'elliptic-curves/cas/research_runtime/supervisor.py')}

def program(model):
    a,b=model[3:];return f'A={a};B={b};gettime();forprime(p=4099,65521,if((4*A^3+27*B^2)%p==0,print("B|",p),print("T|",p,"|",p+1-ellcard(ellinit([Mod(A,p),Mod(B,p)])))));print("MS|",gettime());print("DONE");quit\n'

def parse(stdout,model):
    lines=stdout.splitlines()
    if len(lines)!=len(PRIMES)+2 or lines[-1]!='DONE' or not lines[-2].startswith('MS|'):raise ArithmeticError('trace output framing differs')
    rows=[];a,b=map(int,model[3:])
    for p,line in zip(PRIMES,lines):
        v=line.split('|');bad=(4*a*a*a+27*b*b)%p==0
        if len(v)!=(2 if bad else 3) or v[0]!=('B' if bad else 'T') or int(v[1])!=p:raise ArithmeticError('prime/trace roster differs')
        t=None if bad else int(v[2])
        if t is not None and t*t>4*p:raise ArithmeticError('Hasse bound failed')
        rows.append([p,t])
    return rows,int(lines[-2][3:])

def direct(model,p):
    a,b=(int(q)%p for q in model[3:])
    if (4*a*a*a+27*b*b)%p==0:return None
    squares=bytearray(p)
    for x in range(1,p):squares[x*x%p]=1
    return -sum(0 if (z:=(x*x*x+a*x+b)%p)==0 else (1 if squares[z] else -1) for x in range(p))

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve trace benchmark')
    atlas=cert.read(spec.ATLAS);rows=[]
    for f in atlas['families']:
        path=PARENT/f['family']/'population.json';t=cert.read(path)['retained_candidates'][4]['parameter'];model,pts=spec.specialize(f,t)
        if any(q.denominator!=1 for q in model):raise ArithmeticError('integral family model required')
        rows.append({'family':f['family'],'parameter':t,'model':list(map(str,model)),'population_sha256':cert.hashed(path),'retained_index':4})
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.r17-extended-prime-benchmark.v1','sources':sources(),'gp_sha256':cert.hashed(GP),'gp_version':'2.15.4','rows':rows,'primes':PRIMES,'seconds_per_curve':20,'rss_bytes':536870912,'outer_seconds':150,'gate':'The current562-prime selector has incomplete retrospective recovery. Before another selector campaign, measure cost and exact sample agreement of prime extension on six predetermined retained addresses, index4 in each family. No record labels, exceptional points or catalogue equations enter this benchmark.','direct_validation_indices':[0,1,2,len(PRIMES)//4,len(PRIMES)//2,3*len(PRIMES)//4,len(PRIMES)-2,len(PRIMES)-1],'claim_boundary':'Trace-engine benchmark only; no selected candidate, point search or rank prediction. Singular displayed reductions are omitted, not mistaken for good trace data.'})

def run():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(GP):raise ArithmeticError('benchmark binding changed')
    if (D/'result.json').exists():raise FileExistsError('preserve benchmark result')
    output={'status':'RUNNING','protocol_hash':digest(p),'rows':[]};checkpoint(D/'result.json',output)
    for r in p['rows']:
        model=r['model'];cmd=program(model);folder=D/r['family'];cap=capture([str(GP),'-q','-s','256000000'],input_text=cmd,limits=Limits(p['seconds_per_curve'],p['rss_bytes']),log_path=folder/'gp.log',separate_stderr=True,check=False)
        raw={'program':cmd,'stdout':cap.stdout,'stderr':cap.stderr,'supervision':cap.supervision};checkpoint(folder/'raw.json',raw)
        if cap.supervision['outcome']!='completed' or cap.supervision['returncode']!=0 or cap.stderr:raise ArithmeticError('trace invocation failed or censored')
        traces,ms=parse(cap.stdout,model);checks=[]
        for i in p['direct_validation_indices']:
            q,t=traces[i];expected=direct(model,q)
            if t!=expected:raise ArithmeticError('direct character sum disagrees')
            checks.append([q,expected])
        row={**r,'status':'PASS','trace_count':len(traces),'cpu_ms':ms,'wall_seconds':cap.supervision['wall_seconds'],'direct_checks':checks,'raw_sha256':cert.hashed(folder/'raw.json')};output['rows'].append(row);checkpoint(D/'result.json',output);print('EXTENDED TRACE BENCHMARK',r['family'],len(traces),ms,'ms',flush=True)
    output['status']='PASS';checkpoint(D/'result.json',output)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run']);a=p.parse_args();globals()[a.stage]()
