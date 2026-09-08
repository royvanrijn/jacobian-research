#!/usr/bin/env python3
"""Frozen equal-area skew rectangles chosen only by exact family coefficient bounds."""
import argparse
from math import gcd
from pathlib import Path
import certify_compact_r17_candidates as cert
import scan_higher_r17_stratified as parent
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=parent.ROOT;CAS=parent.CAS;ART=parent.ART;D=parent.LOCAL/'skew-r17-boxes-v1';BINARY=parent.BINARY;BOUNDS=ART/'r17_parameter_box_skew_v1.json'

def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(parent.__file__).resolve(),Path(parent.engine.__file__).resolve(),Path(parent.strict.__file__).resolve(),BINARY,BOUNDS,CAS/'audit_r17_parameter_box_skew.py')}
def roster():
    bounds=cert.read(BOUNDS)
    if bounds['status']!='PASS_EXACT_FINITE_BOUND_OPTIMIZATION':raise ArithmeticError('exact coefficient-bound audit required')
    chosen={r['family']:r for r in bounds['rows'] if r['selected_k']!=0 and cert.F(r['weighted_bound_improvement'])>=256};rows=[]
    for r in parent.roster():
        if r['family'] not in chosen:continue
        q=chosen[r['family']];N=q['selected_numerator_bound'];M=q['selected_denominator_bound'];rows.append({**r,'k':q['selected_k'],'numerator_bound':N,'denominator_bound':M,'primitive_population':parent.strict.population(N,M,r['shard'],64),'small_numerator_bound':int(32*cert.F(2)**q['selected_k']),'small_denominator_bound':int(32/cert.F(2)**q['selected_k'])})
    if len(rows)!=4 or {r['family'] for r in rows}!={'08234','08f72'}:raise ArithmeticError('fixed two improving families required')
    return rows

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve skew population protocol')
    gate=ART/'specialized_parity_portable_replay_v1.json'
    if cert.read(gate)['status']!='PASS':raise ArithmeticError('completed specialized exposure test required')
    rows=roster();checkpoint(D/'protocol.json',{'schema':'elliptic-curves.skew-r17-boxes.v1','sources':sources(),'original_root':str(ROOT),'gate_sha256':cert.hashed(gate),'rows':rows,'keep_per_slice':512,'small_keep':2048,'prime_bound':4093,'seconds_per_call':60,'rss_bytes':1073741824,'maximum_workers':1,'outer_seconds':600,'total_primitive_population':sum(r['primitive_population'] for r in rows),'gate':'The exact finite equal-area rectangle audit finds coefficient-bound improvements only for08234 and08f72 among the nine declared axis scalings. They meet the fixed improvement factor256. For08234 the weighted triangle bound improves by more than2^54;08f72 improves by more than2^9. This is a conservative coefficient bound, not a rank or visibility predictor. Test this unmeasured population choice without changing S1 scores, signed denominator residue classes or retained counts.','population':'Use the two signs and the same pre-existing SHA-selected denominator residues modulo64 for each of the two improving families. The08234 rectangle is |n|<=262144,d<=4096;08f72 is |n|<=65536,d<=16384. Each has the same N*D as32768 squared; exact primitive counts are computed separately. No public parameter or point enters shape selection. Do not rerun the four unchanged-family square boxes.','validation':'Before any large call, enumerate all primitive pairs in each signed small rectangle obtained from base height32 with the same aspect ratio, without slicing. Retain every pair and verify complete score/order from canonical traces. Then retain512 rows per large signed slice and recompute every returned562-prime score. No score modification or newly computed trace.','future_scope':'After population replay, a separately frozen protocol may extend the2048 retained addresses using cached overlaps and at most four point candidates per family. This protocol itself runs no point search and no automatic refill.','boundaries':'Exact finite coefficient inequalities, returned scores and primitive counts. Whole large-box ordering trusts the pinned scanner. Equal area is not equal primitive population or arithmetic height, and bounds are not global minimality, conductor, rank prediction or point solubility. No new rank or novelty claim.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['rows']!=roster():raise ArithmeticError('frozen skew inputs differ')
    return p

def call(row,p,small,check):
    N,M,keep,shard,shards=(row['small_numerator_bound'],row['small_denominator_bound'],2048,0,1) if small else (row['numerator_bound'],row['denominator_bound'],512,row['shard'],64)
    folder=D/row['id'];path=folder/('small.json' if small else 'large.json');origin=Path(p['original_root']);command=[str(origin/BINARY.relative_to(ROOT)),str(origin/row['table']['path']),*map(str,(N,M,keep,shard,shards))]
    if cert.hashed(ROOT/row['table']['path'])!=row['table']['sha256']:raise ArithmeticError('signed score table changed')
    if not check:
        if ROOT!=origin or path.exists():raise FileExistsError('preserve original scanner invocation')
        c=capture(command,limits=Limits(p['seconds_per_call'],p['rss_bytes']),log_path=path.with_suffix('.log'),separate_stderr=True,check=False);checkpoint(path,{'command':command,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    if raw['command']!=command or raw['supervision']['command']!=command or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('skew scanner invocation differs or failed')
    rows,summary=parent.engine.old.parse(raw['stdout'],row['sign']);count=parent.strict.population(N,M,shard,shards)
    if summary!=[N,M,keep,count,count,min(keep,count)] or raw['stdout'].splitlines().count(f'R {shard} {shards}')!=1:raise ArithmeticError('exact primitive population framing differs')
    if len({r['parameter'] for r in rows})!=len(rows) or any(gcd(abs(r['numerator']),r['denominator'])!=1 or not 1<=abs(r['numerator'])<=N or not 1<=r['denominator']<=M or r['numerator']*row['sign']<=0 or (r['denominator']-1)%shards!=shard for r in rows):raise ArithmeticError('returned parameter differs')
    if small and {(r['numerator'],r['denominator']) for r in rows}!={(row['sign']*n,d) for n in range(1,N+1) for d in range(1,M+1) if gcd(n,d)==1}:raise ArithmeticError('small exhaustive roster differs')
    if rows!=sorted(rows,key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],abs(r['numerator']))):raise ArithmeticError('ordered output differs')
    hashes=parent.engine.exact_scores(rows,{'family':row['family']})
    return {'id':row['id'],'family':row['family'],'sign':row['sign'],'rows':rows,'summary':summary,'canonical_table_hashes':hashes,'raw_sha256':cert.hashed(path),'wall_seconds':raw['supervision']['wall_seconds']}

def run(check=False):
    p=protocol();out=D/'result.json'
    if not check and out.exists():raise FileExistsError('preserve skew population')
    data={'status':'RUNNING','protocol_hash':digest(p),'small_checks':[],'shards':[]}
    for row in p['rows']:
        data['small_checks'].append(call(row,p,True,check))
        if not check:checkpoint(out,data)
        print('SKEW SMALL PASS',row['id'],flush=True)
    for row in p['rows']:
        data['shards'].append(call(row,p,False,check))
        if not check:checkpoint(out,data)
        print('SKEW LARGE PASS',row['id'],flush=True)
    data['rows']=parent.merge(data['shards']);old={(r['family'],r['parameter']) for r in cert.read(parent.D/'result.json')['rows']};data['old_S1_retained_overlap']=sum((r['family'],r['parameter']) in old for r in data['rows']);data['status']='COMPLETE_FROZEN_SKEW_POPULATION'
    if check:
        if cert.read(out)!=data:raise ArithmeticError('exact skew population replay differs')
    else:checkpoint(out,data)
    print('SKEW2048 REPLAY PASS; OLD S1 OVERLAP',data['old_S1_retained_overlap'],'POPULATION',p['total_primitive_population'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args();prepare() if a.stage=='prepare' else run(a.stage=='replay')
