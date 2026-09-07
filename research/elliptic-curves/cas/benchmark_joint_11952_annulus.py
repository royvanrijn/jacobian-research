#!/usr/bin/env python3
"""Complete longer-score enumeration on one frozen new outer slice, not a sweep."""
import argparse,json
from math import gcd
from pathlib import Path
import certify_compact_r17_candidates as cert
import scan_11952_new_annulus as prior
import benchmark_retained_short_cache as short
import benchmark_retained_extended_cache as extended
import benchmark_11952_annulus_cache_v3 as reader
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture,Limits
ROOT=prior.ROOT;CAS=prior.CAS;ART=prior.ART;D=ROOT/'artifacts/local/elliptic-curves/joint11952-annulus-benchmark-v1'
CPP=CAS/'newfamily/scan_joint_cache_annulus.cpp';TEST=ROOT/'elliptic-curves/tests/test_joint_cache_annulus.py';BINARY=D/'scanner';OUT=ART/'joint11952_annulus_benchmark_v1.json'
def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),CPP,TEST,short.TABLE,short.OUT,extended.TABLE,extended.OUT,reader.BINARY,reader.OUT,prior.D/'protocol.json',prior.D/'result.json',prior.D/'replay.json',CAS/'verify_periodic_nagao_scanner.py',CAS/'research_runtime/supervisor.py',CAS/'research_runtime/store.py']}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve complete-score benchmark')
    old=prior.protocol();row=old['rows'][0]
    for module in (short,extended):
        proof=cert.read(module.OUT)
        if proof['status']!='PASS' or cert.hashed(module.TABLE)!=proof['cache_sha256']:raise ArithmeticError('pinned exact score cache required')
    if cert.read(reader.OUT)['status']!='PASS' or cert.hashed(reader.BINARY)!=cert.read(reader.OUT)['binary_sha256'] or cert.read(prior.D/'replay.json')['status']!='PASS':raise ArithmeticError('new domain/annulus replay gates required')
    cases=[]
    for sign in (-1,1):
        cases.extend([{'id':str(sign)+'-small','sign':sign,'N':31,'M':29,'inner':7,'shard':0,'shards':1},{'id':str(sign)+'-actual-frame','sign':sign,'N':31,'M':2*row['shards'],'inner':row['shards'],'shard':row['shard'],'shards':row['shards']}])
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.joint11952-annulus-benchmark.v1','sources':sources(),'cases':cases,'real_case':{**row,'N':row['outer'],'M':row['outer']},'real_keep':512,'prime_count':3510,'compile_seconds':60,'test_seconds':30,'seconds_per_call':120,'maximum_workers':1,'real_cost_gate_seconds':45,'rss_bytes':1610612736,'gate':'The new11952 annular short scanner needs only about2.6 seconds per roughly9.6-million-address slice. Retained candidates are then scored through32749, but any short-stage truncation can omit the true longer-score leaders. Existing exact short and extended projective caches make complete combined-score periodic accumulation a concrete alternative. Benchmark its mathematical output and cost before another parameter campaign.','scope':'One fixed already short-scanned outer slice, the first negative quarter, all of whose parameters lie beyond131072 and through524288. Apply all3510 quantized selection primes before any heap admission; keep512 by combined units, good count, denominator and absolute numerator. No public target, catalogue, old point outcome or validation value is read. Exhaustive synthetic tests and four complete real-cache signed frames precede the single full-slice cost gate. Check all returned scores with a separately compiled retained-list reader on both caches; compare their parameter overlap with the previous4096 short survivors. No new original-fibre point search, automatic broader sweep, rank estimate, density or complete65521-score claim.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen joint benchmark inputs differ')
    return p

def execute(name,cmd,seconds,check=False):
    path=D/(name+'.json')
    if not check:
        if path.exists():raise FileExistsError('preserve exact benchmark call')
        c=capture(cmd,limits=Limits(seconds,1610612736),log_path=D/(name+'.log'),separate_stderr=True,check=False)
        checkpoint(path,{'command':cmd,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    if raw['command']!=cmd or raw['supervision']['command']!=cmd or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('joint benchmark call failed/censored')
    return raw

def scan_case(c,keep,name,check=False):
    values=[c['sign'],c['N'],c['M'],keep,c['shard'],c['shards'],c['inner'],3510];cmd=[str(BINARY),str(short.TABLE),str(extended.TABLE),*map(str,values)]
    raw=execute(name,cmd,120,check);lines=raw['stdout'].splitlines()
    population=prior.annulus.counts.population(c['N'],c['M'],c['shard'],c['shards'])-prior.annulus.counts.population(min(c['N'],c['inner']),min(c['M'],c['inner']),c['shard'],c['shards'])
    count=min(keep,population)
    if raw['stderr'] or len(lines)!=count+3 or lines[:2]!=['JOINT_NAGAO_ANNULUS_V1','P '+' '.join(map(str,values))] or lines[-1]!=f'S {population} {count}':raise ArithmeticError('exact joint population frame differs')
    rows=[]
    for line in lines[2:-1]:
        fields=line.split()
        if len(fields)!=5 or fields[0]!='C':raise ArithmeticError('joint candidate frame differs')
        n,d,s,g=map(int,fields[1:]);r={'numerator':n,'denominator':d,'combined_selection_units':s,'combined_good':g,'parameter':str(cert.F(n,d))}
        if n*c['sign']<=0 or gcd(n,d)!=1 or not 1<=abs(n)<=c['N'] or not 1<=d<=c['M'] or not c['inner']<max(abs(n),d) or (d-1)%c['shards']!=c['shard'] or not 0<=g<=3510:raise ArithmeticError('joint candidate outside fixed primitive frame')
        rows.append(r)
    if len({r['parameter'] for r in rows})!=len(rows) or rows!=sorted(rows,key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],abs(r['numerator']))):raise ArithmeticError('joint uniqueness/order differs')
    return rows,raw,population

def references(rows,name,check=False):
    path=D/(name+'-candidates.txt');text='R17-CANDIDATES-V1 '+str(len(rows))+'\n'+''.join(f"{r['numerator']} {r['denominator']}\n" for r in rows)
    if check:
        if path.read_text()!=text:raise ArithmeticError('independent reader roster changed')
    else:
        if path.exists():raise FileExistsError('preserve independent reader input')
        path.write_text(text)
    values=[]
    for label,module,count in [('short',short,562),('extended',extended,2948)]:
        raw=execute(name+'-'+label,[str(reader.BINARY),str(module.TABLE),str(path),'524288'],120,check);lines=raw['stdout'].splitlines();pairs=[]
        if raw['stderr'] or len(lines)!=len(rows)+1 or lines[-1]!=f'S {len(rows)} {count}':raise ArithmeticError('independent cache frame differs')
        for i,line in enumerate(lines[:-1]):
            v=line.split()
            if len(v)!=4 or v[:2]!=['R',str(i)]:raise ArithmeticError('independent cache order differs')
            pairs.append(tuple(map(int,v[2:])))
        values.append(pairs)
    return [{**r,'combined_selection_units':a[0]+b[0],'combined_good':a[1]+b[1]} for r,a,b in zip(rows,*values)]

def run(check=False):
    p=protocol();execute('compile',['g++','-O3','-std=c++17',str(CPP),'-o',str(BINARY)],60,check);test=execute('tests',['/usr/bin/python3',str(TEST)],30,check)
    if 'Ran 3 tests' not in test['stderr'] or not test['stderr'].rstrip().endswith('OK'):raise ArithmeticError('all three exhaustive test groups required')
    tests=[]
    for c in p['cases']:
        bare=[{'numerator':c['sign']*n,'denominator':d,'parameter':str(cert.F(c['sign']*n,d))} for d in range(c['shard']+1,c['M']+1,c['shards']) for n in range(1,c['N']+1) if gcd(n,d)==1 and max(n,d)>c['inner']]
        expected=references(bare,c['id']+'-reference',check);expected.sort(key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],abs(r['numerator'])))
        for k in (1000,7):
            got,raw,count=scan_case(c,k,c['id']+'-k'+str(k),check)
            if got!=expected[:k] or count!=len(expected):raise ArithmeticError('entire exact real-cache frame/top7 differs')
        tests.append({'id':c['id'],'population':len(expected),'complete_order_and_top7':'PASS'})
    rows,raw,count=scan_case(p['real_case'],p['real_keep'],'real-slice',check)
    if references(rows,'real-reference',check)!=rows:raise ArithmeticError('some full-slice retained combined score differs')
    old=next(s for s in cert.read(prior.D/'result.json')['shards'] if s['id']==p['real_case']['id']);prior_parameters={r['parameter'] for r in old['rows']};missed=[r['parameter'] for r in rows if r['parameter'] not in prior_parameters]
    elapsed=raw['supervision']['wall_seconds'];result={'schema':'elliptic-curves.joint11952-annulus-benchmark-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'binary_sha256':cert.hashed(BINARY),'tests':tests,'primitive_population':count,'selection_primes':3510,'retained_rows':rows,'wall_seconds':elapsed,'cost_gate_passed':elapsed<=p['real_cost_gate_seconds'],'prior_short_retained':len(prior_parameters),'joint_top512_outside_prior_short4096':missed,'claim_boundary':p['scope']}
    if check:
        if cert.read(OUT)!=result:raise ArithmeticError('joint benchmark replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve joint benchmark certificate')
        checkpoint(OUT,result)
    print('JOINT11952 FULL-SLICE BENCHMARK',count,'ADDRESSES',elapsed,'SECONDS; COST',result['cost_gate_passed'],'MISSED BY SHORT4096',len(missed),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','check']);a=p.parse_args();prepare() if a.stage=='prepare' else run(a.stage=='check')
