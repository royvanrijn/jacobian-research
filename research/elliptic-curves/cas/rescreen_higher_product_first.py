#!/usr/bin/env python3
"""Product-first retention on the unchanged twelve higher-parameter slices."""
import argparse
from math import log,gcd
from pathlib import Path
import certify_compact_r17_candidates as cert
import scan_higher_r17_stratified as parent
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=parent.ROOT;CAS=parent.CAS;ART=parent.ART;D=parent.LOCAL/'higher32768-product-first-v1';BINARY=parent.BINARY
def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(parent.__file__).resolve(),Path(parent.engine.__file__).resolve(),Path(parent.strict.__file__).resolve(),BINARY,ROOT/'elliptic-curves/cas/compare_bounded_prime_selectors.py')}
def freeze():
    if (D/'protocol.json').exists():raise FileExistsError('preserve product-first protocol')
    gates=[ART/'product22_comparison_portable_replay_v1.json',ART/'higher_r17_masked_relations_v1.json'];a,b=map(cert.read,gates)
    if a['status']!='PASS' or b['status']!='COMPLETE_BOUNDED_RELATION_AUDIT' or any(r['recovery_status']!='RECOVERED_KNOWN_DIRECTION' for r in b['rows']):raise ArithmeticError('completed paired and masked diagnostic gates required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.higher-product-first.v1','sources':sources(),'gate_hashes':{str(p.relative_to(ROOT)):cert.hashed(p) for p in gates},'parent_protocol_sha256':cert.hashed(parent.D/'protocol.json'),'rows':parent.roster(),'height':32768,'keep_per_slice':512,'small_box':[31,29],'small_keep':1000,'seconds_per_call':60,'rss_bytes':1073741824,'maximum_workers':1,'outer_seconds':600,'formula':'At each good displayed reduction replace the existing table contribution by round(log((p+1-a_p)/p)*10^12). Bad displayed reductions contribute zero. Keep the identical binary, local traces, primitive population, denominator slices, sign convention and ties. The duplicate held prime in this legacy table format is output only, never a selection tie breaker.','gate':'The paired product-versus-S1 comparison was restricted to6144 addresses already truncated by short-prime S1. It cannot test product-first retention of discarded addresses. Six exact masked controls recover their withheld generic directions on the larger fibres, without proving exceptional sensitivity. Test the earlier retention-stage gap on exactly the same122368792 primitive addresses, using only cached traces and no larger height or new slice.','stages':'Freeze before reweighting. Before any large slice call, check all576 primitive scores and full ordered output on the31-by29 box for every signed family. Then retain512 product-first rows from each unchanged signed slice. Independently recompute every returned score from canonical residue traces. Record overlap with the old6144 pool only after selection.','future_scope':'This protocol contains no additional finite-field traces, point searches or automatic extension. A separate protocol may use completed population and replay evidence to decide a bounded next step.','boundaries':'Reproducible numerical selection, not a rank bound or exact transcendental ordering. Full large-box top retention trusts the unchanged pinned scanner; the small populations are exhaustive checks. No calibrated superiority, point absence, density, new curve or novelty claim.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['parent_protocol_sha256']!=cert.hashed(parent.D/'protocol.json') or p['rows']!=parent.roster():raise ArithmeticError('frozen product-first inputs changed')
    return p
def table(row,check):
    source=ROOT/row['table']['path'];out=D/row['id']/'product-table.txt'
    if cert.hashed(source)!=row['table']['sha256']:raise ArithmeticError('original signed table changed')
    if not check and out.exists():raise FileExistsError('preserve product table')
    lines=[];prime=None;remaining=0
    for line in source.read_text().splitlines():
        fields=line.split()
        if remaining:
            good,trace,old=map(int,fields)
            if good not in (0,1):raise ArithmeticError('invalid local good flag')
            units=round(log((prime+1-trace)/prime)*10**12) if good else 0;lines.append(f'{good} {trace} {units}');remaining-=1
        else:
            lines.append(line)
            if fields[0]=='P':prime=int(fields[1]);remaining=prime+1
    if remaining:raise ArithmeticError('truncated original table')
    text='\n'.join(lines)+'\n'
    if check:
        if out.read_text()!=text:raise ArithmeticError('product table differs')
    else:out.parent.mkdir(parents=True,exist_ok=True);out.write_text(text)
    return out
def exact(rows,family):
    scores=[0]*len(rows);goods=[0]*len(rows);bindings={}
    f=next(f for f in cert.read(parent.engine.spec.ATLAS)['families'] if f['family']==family);model={k:f[k]+['0']*(n-len(f[k])) for k,n in [('A_coefficients_low_to_high',9),('B_coefficients_low_to_high',13)]}
    for q in (q for q in range(5,4094) if parent.engine._is_prime(q)):
        path=parent.engine.old.PARENT/family/'trace-tables'/f'{q}.json';t=cert.read(path)
        if t['input']!={'family':family,'model_hash':digest(model),'prime':q}:raise ArithmeticError('canonical local model differs')
        bindings[str(path.relative_to(ROOT))]=cert.hashed(path)
        for i,r in enumerate(rows):
            x=r['numerator']*pow(r['denominator'],-1,q)%q if r['denominator']%q else q
            if t['good'][x]:scores[i]+=round(log((q+1-t['traces'][x])/q)*10**12);goods[i]+=1
    if any((r['score_units'],r['good_primes'])!=(s,g) for r,s,g in zip(rows,scores,goods)):raise ArithmeticError('canonical product score differs')
    return bindings
def call(row,p,small,check):
    N,M,keep,shard,shards=(31,29,1000,0,1) if small else (32768,32768,512,row['shard'],64);folder=D/row['id'];path=folder/('small.json' if small else 'large.json');tp=folder/'product-table.txt';command=[str(BINARY),str(tp),*map(str,(N,M,keep,shard,shards))]
    if not check:
        if path.exists():raise FileExistsError('preserve scanner result')
        c=capture(command,limits=Limits(p['seconds_per_call'],p['rss_bytes']),log_path=path.with_suffix('.log'),separate_stderr=True,check=False);checkpoint(path,{'command':command,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    if raw['command']!=command or raw['supervision']['command']!=command or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('product scanner invocation failed or differs')
    rows,summary=parent.engine.old.parse(raw['stdout'],row['sign']);count=parent.strict.population(N,M,shard,shards)
    if summary!=[N,M,keep,count,count,min(keep,count)]:raise ArithmeticError('exact primitive population differs')
    if len({r['parameter'] for r in rows})!=len(rows) or any(gcd(abs(r['numerator']),r['denominator'])!=1 or not 1<=abs(r['numerator'])<=N or not 1<=r['denominator']<=M or r['numerator']*row['sign']<=0 or (r['denominator']-1)%shards!=shard for r in rows):raise ArithmeticError('returned address differs')
    if small and {(r['numerator'],r['denominator']) for r in rows}!={(row['sign']*n,d) for n in range(1,32) for d in range(1,30) if gcd(n,d)==1}:raise ArithmeticError('small exhaustive roster differs')
    if rows!=sorted(rows,key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],abs(r['numerator']))):raise ArithmeticError('complete returned order differs')
    hashes=exact(rows,row['family']);return {'id':row['id'],'family':row['family'],'sign':row['sign'],'rows':rows,'summary':summary,'canonical_table_hashes':hashes,'table_sha256':cert.hashed(tp),'raw_sha256':cert.hashed(path),'wall_seconds':raw['supervision']['wall_seconds']}
def run(check=False):
    p=protocol();out=D/'result.json'
    if not check and out.exists():raise FileExistsError('preserve product-first result')
    data={'status':'RUNNING','protocol_hash':digest(p),'small_checks':[],'shards':[]}
    for row in p['rows']:
        table(row,check);data['small_checks'].append(call(row,p,True,check))
        if not check:checkpoint(out,data)
        print('PRODUCT-FIRST SMALL PASS',row['id'],flush=True)
    for row in p['rows']:
        data['shards'].append(call(row,p,False,check))
        if not check:checkpoint(out,data)
        print('PRODUCT-FIRST LARGE PASS',row['id'],flush=True)
    data['rows']=parent.merge(data['shards']);old=cert.read(parent.D/'result.json');oldset={(r['family'],r['parameter']) for r in old['rows']};data['previous_pool_overlap']=sum((r['family'],r['parameter']) in oldset for r in data['rows']);data['status']='COMPLETE_FIXED_PRODUCT_FIRST_POPULATION'
    if check:
        if cert.read(out)!=data:raise ArithmeticError('full product-first replay differs')
    else:checkpoint(out,data)
    print('PRODUCT-FIRST6144 PASS; OLD POOL OVERLAP',data['previous_pool_overlap'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['freeze','run','replay']);a=p.parse_args();freeze() if a.stage=='freeze' else run(a.stage=='replay')
