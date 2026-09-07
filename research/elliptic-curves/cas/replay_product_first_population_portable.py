#!/usr/bin/env python3
"""Read-only product-first replay with recorded argv bound to its original root."""
from pathlib import Path
from math import gcd
import rescreen_higher_product_first as original
from rescreen_higher_product_first import cert,parent,ROOT,D,BINARY,exact,digest
def checked(row,p,small):
    N,M,keep,shard,shards=(31,29,1000,0,1) if small else (32768,32768,512,row['shard'],64);folder=D/row['id'];path=folder/('small.json' if small else 'large.json');tp=folder/'product-table.txt';command=[str(BINARY),str(tp),*map(str,(N,M,keep,shard,shards))]
    raw=cert.read(path);recorded_root=Path(raw['command'][0]).parents[4];command=[str(recorded_root/BINARY.relative_to(ROOT)),str(recorded_root/tp.relative_to(ROOT)),*map(str,(N,M,keep,shard,shards))]
    if raw['command']!=command or raw['supervision']['command']!=command or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('product scanner invocation failed or differs')
    rows,summary=parent.engine.old.parse(raw['stdout'],row['sign']);count=parent.strict.population(N,M,shard,shards)
    if summary!=[N,M,keep,count,count,min(keep,count)]:raise ArithmeticError('exact primitive population differs')
    if len({r['parameter'] for r in rows})!=len(rows) or any(gcd(abs(r['numerator']),r['denominator'])!=1 or not 1<=abs(r['numerator'])<=N or not 1<=r['denominator']<=M or r['numerator']*row['sign']<=0 or (r['denominator']-1)%shards!=shard for r in rows):raise ArithmeticError('returned address differs')
    if small and {(r['numerator'],r['denominator']) for r in rows}!={(row['sign']*n,d) for n in range(1,32) for d in range(1,30) if gcd(n,d)==1}:raise ArithmeticError('small exhaustive roster differs')
    if rows!=sorted(rows,key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],abs(r['numerator']))):raise ArithmeticError('complete returned order differs')
    hashes=exact(rows,row['family']);return {'id':row['id'],'family':row['family'],'sign':row['sign'],'rows':rows,'summary':summary,'canonical_table_hashes':hashes,'table_sha256':cert.hashed(tp),'raw_sha256':cert.hashed(path),'wall_seconds':raw['supervision']['wall_seconds']}
def main():
    p=original.protocol();data={'status':'RUNNING','protocol_hash':digest(p),'small_checks':[],'shards':[]}
    for row in p['rows']:
        original.table(row,True);data['small_checks'].append(checked(row,p,True))
    for row in p['rows']:data['shards'].append(checked(row,p,False))
    data['rows']=parent.merge(data['shards']);old=cert.read(parent.D/'result.json');oldset={(r['family'],r['parameter']) for r in old['rows']};data['previous_pool_overlap']=sum((r['family'],r['parameter']) in oldset for r in data['rows']);data['status']='COMPLETE_FIXED_PRODUCT_FIRST_POPULATION'
    if cert.read(D/'result.json')!=data:raise ArithmeticError('relocated product-first replay differs')
    print('RELOCATED PRODUCT-FIRST SCORE AND POPULATION REPLAY PASS',flush=True)
if __name__=='__main__':main()
