#!/usr/bin/env python3
"""Read-only replay with recorded argv rooted at the original pinned checkout."""
from pathlib import Path
import scan_higher_r17_stratified as original
import certify_compact_r17_candidates as cert
from research_runtime.store import digest

def checked(row,p):
    path=original.D/row['id']/'raw.json';raw=cert.read(path)
    recorded_root=Path(raw['command'][0]).parents[4]
    cmd=[str(recorded_root/original.BINARY.relative_to(original.ROOT)),str(recorded_root/row['table']['path']),'32768','32768','512',str(row['shard']),'64']
    if raw['command']!=cmd or raw['supervision']['command']!=cmd or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('immutable scanner failed/censored or differs')
    rows,summary=original.engine.old.parse(raw['stdout'],row['sign']);count=row['primitive_population']
    if summary!=[32768,32768,512,count,count,512] or raw['stdout'].splitlines().count(f"R {row['shard']} 64")!=1:raise ArithmeticError('population framing differs')
    if len(rows)!=512 or len({r['parameter'] for r in rows})!=512 or any(original.engine.gcd(abs(r['numerator']),r['denominator'])!=1 or not 1<=abs(r['numerator'])<=32768 or not 1<=r['denominator']<=32768 or r['numerator']*row['sign']<=0 or (r['denominator']-1)%64!=row['shard'] for r in rows):raise ArithmeticError('primitive signed slice roster differs')
    if rows!=sorted(rows,key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],abs(r['numerator']))):raise ArithmeticError('ordered slice differs')
    hashes=original.engine.exact_scores(rows,{'family':row['family']})
    return {**row,'status':'PASS','rows':rows,'summary':summary,'raw_sha256':cert.hashed(path),'canonical_table_hashes':hashes,'wall_seconds':raw['supervision']['wall_seconds']}

def main():
    p=original.protocol();d=cert.read(original.D/'result.json');shards=[checked(r,p) for r in p['rows']]
    if d['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or d['protocol_hash']!=digest(p) or d['shards']!=shards or d['rows']!=original.merge(shards):raise ArithmeticError('exact stratified population replay differs')
    print('REPLAYED RELOCATED12 SLICES,6144 SCORES AND',p['total_primitive_population'],'PRIMITIVE ADDRESSES',flush=True)
if __name__=='__main__':main()
