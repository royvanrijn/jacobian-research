#!/usr/bin/env python3
"""Strict command/binary/population binding and exhaustive small signed-frame checks."""
import argparse
from pathlib import Path
import benchmark_periodic_nagao_scanner as b
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture,Limits
D=b.ROOT/'artifacts/local/elliptic-curves/periodic-nagao-strict-v1';OUT=b.ART/'periodic_nagao_scanner_strict_v1.json'
def population(N,M,shard,shards):
    spf=list(range(M+1))
    for q in range(2,M+1):
        if spf[q]==q:
            for k in range(q*q,M+1,q):
                if spf[k]==k:spf[k]=q
    total=0
    for d in range(1+shard,M+1,shards):
        n=d;factors=[]
        while n>1:
            q=spf[n];factors.append(q)
            while n%q==0:n//=q
        divisors=[(1,1)]
        for q in factors:divisors += [(v*q,-s) for v,s in list(divisors)]
        total+=sum(s*(N//v) for v,s in divisors)
    return total

def main(check):
    p=b.protocol();build=cert.read(b.D/'build.json');binary=b.D/'scanner-v2'
    if build['binary_sha256']!=cert.hashed(binary) or build['supervision']['outcome']!='completed' or build['supervision']['returncode']!=0 or build['command']!=['g++','-O3','-std=c++17',str(b.CPP),'-o',str(binary)]:raise ArithmeticError('prototype build binding differs')
    bindings={}
    for case in p['cases']:
        count=population(case['numerator_bound'],case['denominator_bound'],case['shard'],case['shards'])
        for v,path in [('v1',b.old.BINARY),('v2',binary)]:
            rawpath=b.D/case['name']/(v+'.json');raw=cert.read(rawpath);expected=[str(path),str(b.ROOT/p['table']['path']),*[str(case[k]) for k in ('numerator_bound','denominator_bound','keep','shard','shards')]];rows,summary=b.old.parse(raw['stdout'],-1)
            if raw['command']!=expected or raw['supervision']['command']!=expected or summary[3:5]!=[count,count]:raise ArithmeticError('benchmark invocation or exact count differs')
            bindings[str(rawpath.relative_to(b.ROOT))]=cert.hashed(rawpath)
    parent=cert.read(b.old.DIRECTORY/'protocol.json');results=[]
    for family in parent['families']:
        for sign in (-1,1):
            table=parent['trace_tables'][family][str(sign)];input_path=b.ROOT/table['path']
            if cert.hashed(input_path)!=table['sha256']:raise ArithmeticError('signed table differs')
            values=[]
            for version,path in [('v1',b.old.BINARY),('v2',binary)]:
                folder=D/family/str(sign);out=folder/(version+'.json');command=[str(path),str(input_path),'31','29','1000','0','1']
                if not check:
                    if out.exists():raise FileExistsError('preserve small exhaustive scanner call')
                    c=capture(command,limits=Limits(15,1073741824),log_path=folder/(version+'.log'),separate_stderr=True,check=False);checkpoint(out,{'command':command,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
                raw=cert.read(out)
                if raw['command']!=command or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('small exact call failed/censored')
                rows,summary=b.old.parse(raw['stdout'],sign);count=population(31,29,0,1)
                if summary!=[31,29,1000,count,count,count] or {(r['numerator'],r['denominator']) for r in rows}!={(sign*n,d) for n in range(1,32) for d in range(1,30) if b.gcd(n,d)==1}:raise ArithmeticError('small exhaustive primitive population differs')
                values.append((rows,raw['stdout']));bindings[str(out.relative_to(b.ROOT))]=cert.hashed(out)
            if values[0]!=values[1]:raise ArithmeticError('entire small ordered population differs')
            hashes=b.exact_scores(values[1][0],{'family':family});bindings.update(hashes);results.append({'family':family,'sign':sign,'curves_checked':count});print('EXACT SMALL SIGNED FRAME',family,sign,count,flush=True)
    result={'schema':'elliptic-curves.periodic-nagao-strict.v1','status':'PASS','sources':{str(q.relative_to(b.ROOT)):cert.hashed(q) for q in (Path(__file__).resolve(),Path(b.__file__).resolve(),b.D/'protocol.json',b.D/'build.json',binary,b.ART/'periodic_nagao_scanner_benchmark_v1.json')},'bindings':bindings,'small_populations':results,'claim_boundary':'Exact old/new invocation, binary and primitive-count bindings for two real benchmark populations; all primitive scores and complete ordered outputs on31-by29 boxes for all twelve signed compactR17 frames. No broad scan or point search. Performance evidence remains the two fixed benchmark timing pairs.'}
    if check:
        if cert.read(OUT)!=result:raise ArithmeticError('strict benchmark replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve strict benchmark proof')
        checkpoint(OUT,result)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();main(a.check)
