#!/usr/bin/env python3
"""Exact primitive-hit and bounded cost comparison for the optional v2 GMP sieve."""
import argparse,json,subprocess,random
from pathlib import Path
from math import gcd,isqrt
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture,Limits
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/pointed-sieve-v2-benchmark-v1';CONTROL=ROOT/'artifacts/local/elliptic-curves/native11952-pari49-control-v1/candidate-00/result.json';OUT=ART/'pointed_sieve_v2_benchmark_v1.json'
def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),CAS/'pointed_quartic_sieve.cpp',CAS/'pointed_quartic_sieve_v2.cpp',CONTROL)}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve sieve benchmark')
    rng=random.Random(20260906);cases=[[0,0,0,0,0],[1,0,0,0,1],[-1,0,0,0,-1],[1,0,-2,0,1],[0,0,0,0,1],[9,0,0,0,-1]]+[[rng.randrange(-12,13) for _ in range(5)] for _ in range(42)]
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.pointed-sieve-v2-benchmark.protocol.v1','sources':sources(),'test_cases':cases,'test_height':17,'test_intervals':[[1,17],[5,13]],'test_seconds':3,'benchmark_indices':[0,12,48],'benchmark_height':100000,'benchmark_seconds':20,'rss_bytes':1073741824,'gate':'At million height the unchanged GMP control chart completed only117838 denominators in60seconds, while PARI completed the full box and recovered29. The frozen scalar GMP loop performs a remainder for each filter/64-numerator block and scans locally impossible primitive denominator rows. Test exact primitive-square filtering and empty-row rejection plus division-free residue stepping in a separate v2 source.','proof':'A primitive pair cannot have n=d=0 modulo any prime. Every retained rational point must pass the square condition at every prime. If a denominator row has no such numerator residue, it contains no primitive rational point. For r in[0,p), replacing (r+64)modp by r+(64modp) and one conditional subtraction is exact.','boundaries':'No production migration or new mathematical rank claim. Exhaustive Python enumeration on48 fixed quartics and two denominator intervals, then full-hit agreement against the original frozen worker on three real control charts. Both code versions, compiler output, raw transcripts and timeouts remain. Full-box execution still trusts the pinned worker.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('benchmark sources changed')
    return p
def build():
    p=protocol()
    for version,name in [('v1','pointed_quartic_sieve.cpp'),('v2','pointed_quartic_sieve_v2.cpp')]:
        binary=D/version
        if binary.exists():raise FileExistsError('preserve benchmark executable')
        args=['g++','-O3','-std=c++17',str(CAS/name),'-lgmpxx','-lgmp','-o',str(binary)];r=capture(args,limits=Limits(60,p['rss_bytes']),log_path=D/(version+'-build.log'),check=True)
        v=subprocess.run([str(binary),'--version'],capture_output=True,text=True,check=True)
        checkpoint(D/(version+'-build.json'),{'command':args,'supervision':r.supervision,'binary_sha256':cert.hashed(binary),'version':v.stdout})
def parse(stdout,height,first,last,coefficients):
    lines=stdout.splitlines();primes=[];hits=[];done=None
    for line in lines:
        fields=line.split()
        if fields[0]=='PRIMES':primes=list(map(int,fields[1:]))
        elif fields[0]=='POINT':hits.append(tuple(map(int,fields[1:])))
        elif fields[0]=='DONE':done=fields[1:]
        else:raise ArithmeticError('unknown raw sieve line')
    if done is None or int(done[0])!=last or len(hits)!=int(done[4]) or len(set(hits))!=len(hits):raise ArithmeticError('complete finite interval required')
    for n,d,r in hits:
        if not -height<=n<=height or not first<=d<=last or gcd(abs(n),d)!=1 or r<0 or r*r!=sum(c*n**i*d**(4-i) for i,c in enumerate(coefficients)):raise ArithmeticError('incorrect primitive square hit')
    return {'hits':[list(h) for h in hits],'completed_denominator':int(done[0]),'word_survivors':int(done[1]),'all_prime_survivors':int(done[2]),'exact_tests':int(done[3]),'seconds':float(done[5]),'primes':primes}
def call(version,coefficients,height,first,last,seconds,folder):
    binary=D/version;build=cert.read(D/(version+'-build.json'))
    if cert.hashed(binary)!=build['binary_sha256']:raise ArithmeticError('compiled binary differs')
    path=folder/(version+'.json')
    if path.exists():raise FileExistsError('preserve raw sieve attempt')
    command=f'{height} {first} {last} {seconds}\n'+'\n'.join(map(str,coefficients))+'\n';cap=capture([str(binary)],input_text=command,limits=Limits(seconds+5,1073741824),log_path=folder/(version+'.log'),separate_stderr=True,check=False)
    raw={'input':command,'stdout':cap.stdout,'stderr':cap.stderr,'supervision':cap.supervision,'binary_sha256':build['binary_sha256']};checkpoint(path,raw)
    if cap.stderr or cap.supervision['outcome']!='completed' or cap.supervision['returncode']!=0:raise ArithmeticError('sieve call failed/censored')
    return parse(cap.stdout,height,first,last,coefficients)
def exhaustive(c,h,first,last):
    points=[]
    for d in range(first,last+1):
        for n in range(-h,h+1):
            if gcd(abs(n),d)!=1:continue
            v=sum(f*n**i*d**(4-i) for i,f in enumerate(c))
            if v>=0 and isqrt(v)**2==v:points.append([n,d,isqrt(v)])
    return points
def run():
    p=protocol();out=D/'result.json'
    if out.exists():raise FileExistsError('preserve benchmark result')
    data={'status':'RUNNING','sources':sources(),'tests':[],'real_charts':[]};checkpoint(out,data)
    for i,c in enumerate(p['test_cases']):
        for first,last in p['test_intervals']:
            expected=exhaustive(c,p['test_height'],first,last);r=call('v2',c,p['test_height'],first,last,p['test_seconds'],D/'tests'/f'{i:03}-{first}-{last}')
            if r['hits']!=expected:raise ArithmeticError('exhaustive primitive enumeration differs')
            data['tests'].append({'case':i,'interval':[first,last],'hits':len(expected)});checkpoint(out,data)
    control=cert.read(CONTROL)
    for index in p['benchmark_indices']:
        c=list(map(int,control['charts'][index]['search']['coefficients']));results={version:call(version,c,p['benchmark_height'],1,p['benchmark_height'],p['benchmark_seconds'],D/'charts'/str(index)) for version in ('v1','v2')}
        if results['v1']['hits']!=results['v2']['hits']:raise ArithmeticError('old/new full-box hits differ')
        data['real_charts'].append({'index':index,'results':results,'speedup':results['v1']['seconds']/results['v2']['seconds']});checkpoint(out,data);print('SIEVE V2',index,results['v1']['seconds'],results['v2']['seconds'],'speedup',data['real_charts'][-1]['speedup'],flush=True)
    data['status']='COMPLETE_FIXED_BENCHMARK';checkpoint(out,data)
def replay():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='COMPLETE_FIXED_BENCHMARK' or d['sources']!=sources() or len(d['tests'])!=96:raise ArithmeticError('complete96 exhaustive checks required')
    for r in d['tests']:
        i=r['case'];first,last=r['interval'];c=p['test_cases'][i];raw=cert.read(D/'tests'/f'{i:03}-{first}-{last}'/'v2.json');parsed=parse(raw['stdout'],p['test_height'],first,last,c)
        if parsed['hits']!=exhaustive(c,p['test_height'],first,last) or r['hits']!=len(parsed['hits']):raise ArithmeticError('exhaustive replay differs')
    control=cert.read(CONTROL)
    if [r['index'] for r in d['real_charts']]!=p['benchmark_indices']:raise ArithmeticError('real chart roster differs')
    for row in d['real_charts']:
        c=list(map(int,control['charts'][row['index']]['search']['coefficients']))
        for version in ('v1','v2'):
            raw=cert.read(D/'charts'/str(row['index'])/(version+'.json'));actual=parse(raw['stdout'],p['benchmark_height'],1,p['benchmark_height'],c)
            if actual!=row['results'][version]:raise ArithmeticError('raw comparison replay differs')
        if row['results']['v1']['hits']!=row['results']['v2']['hits']:raise ArithmeticError('complete primitive hit sets differ')
    if OUT.exists():raise FileExistsError('preserve benchmark certificate')
    checkpoint(OUT,{'schema':'elliptic-curves.pointed-sieve-v2-benchmark.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'result_sha256':cert.hashed(D/'result.json'),'tests_checked':96,'real_charts':d['real_charts'],'claim_boundary':p['boundaries'],'argument':p['proof']});print('REPLAYED96 EXHAUSTIVE CASES AND THREE REAL FULL BOXES',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','build','run','replay']);a=p.parse_args();globals()[a.stage]()
