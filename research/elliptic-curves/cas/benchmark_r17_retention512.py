#!/usr/bin/env python3
"""Cost and exact-output gate for retaining512 per sign on the existing H4096 box."""
import argparse
from pathlib import Path
from math import log,gcd
import certify_compact_r17_candidates as cert
import select_compact_r17_wide as original
import compact_atlas_specialization as spec
from mod2_reduction_independence import _is_prime
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=Path(__file__).resolve().parents[2];LOCAL=ROOT/'artifacts/local/elliptic-curves';ART=ROOT/'artifacts/generated-results/elliptic-curves';D=LOCAL/'r17-retention512-benchmark-v1';PARENT=original.DIRECTORY;FAMILY='103b2';SIGN=-1;KEEP=512

def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),Path(original.__file__).resolve(),original.BINARY,ROOT/'elliptic-curves/cas/newfamily/scan_rational_nagao_tables.cpp',spec.ATLAS,ROOT/'elliptic-curves/cas/research_runtime/supervisor.py')}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve retention benchmark')
    gate=ART/'height_and_discarded_portable_replay_v1.json';g=cert.read(gate);parent=cert.read(PARENT/'protocol.json');table=parent['trace_tables'][FAMILY][str(SIGN)];old=PARENT/FAMILY/f'scan-{SIGN}.json'
    if g['status']!='PASS' or parent['scanner_binary_sha256']!=cert.hashed(original.BINARY) or cert.hashed(ROOT/table['path'])!=table['sha256']:raise ArithmeticError('prior discovery or scanner gate differs')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.r17-retention512-benchmark.v1','sources':sources(),'family':FAMILY,'sign':SIGN,'keep':KEEP,'height':4096,'table':table,'old_shard_sha256':cert.hashed(old),'gate_sha256':cert.hashed(gate),'seconds':120,'rss_bytes':1073741824,'maximum_workers':1,'gate':'Eight new curves, including a globally minimal26, came from saved addresses excluded by the global128 short-prime merge. Before another full retention campaign, benchmark512 per sign on the predetermined first atlas family103b2 and negative sign. Same parameter box, table bytes and scanner. No record labels or public points.','validation':'Require all512 scores recomputed from562 canonical trace tables, identical old128 prefix, exact primitive parameter count, and deterministic absolute-numerator tie order within the negative-sign scanner.','boundaries':'Cost and finite-output replay only, not independent full-population ranking or rank prediction. A wider campaign needs its own fixed protocol after this gate passes.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['keep']!=512 or p['family']!=FAMILY or p['sign']!=SIGN or cert.hashed(ROOT/p['table']['path'])!=p['table']['sha256'] or cert.hashed(PARENT/FAMILY/f'scan-{SIGN}.json')!=p['old_shard_sha256']:raise ArithmeticError('frozen retention benchmark differs')
    return p
def verify_rows(rows,summary,family,sign):
    if len(rows)!=512 or len({r['parameter'] for r in rows})!=512 or summary!=[4096,4096,512,10200039,10200039,512]:raise ArithmeticError('fixed signed population/roster differs')
    if rows!=sorted(rows,key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],abs(r['numerator']))):raise ArithmeticError('scanner tie order differs')
    if any(gcd(abs(r['numerator']),r['denominator'])!=1 or not 1<=abs(r['numerator'])<=4096 or not 1<=r['denominator']<=4096 or r['numerator']*sign<=0 for r in rows):raise ArithmeticError('primitive signed address differs')
    old=cert.read(PARENT/family/f'scan-{sign}.json')
    if rows[:128]!=old['rows']:raise ArithmeticError('old128 prefix differs')
    totals=[0]*512;goods=[0]*512;table_hashes={};primes=[p for p in range(5,4094) if _is_prime(p)]
    f=next(f for f in cert.read(spec.ATLAS)['families'] if f['family']==family);model={key:f[key]+['0']*(n-len(f[key])) for key,n in [('A_coefficients_low_to_high',9),('B_coefficients_low_to_high',13)]}
    for p in primes:
        path=original.PARENT/family/'trace-tables'/f'{p}.json';table=cert.read(path)
        if table['input']!={'family':family,'model_hash':digest(model),'prime':p}:raise ArithmeticError('canonical trace table differs')
        table_hashes[str(path.relative_to(ROOT))]=cert.hashed(path)
        for i,r in enumerate(rows):
            n,d=r['numerator'],r['denominator'];t=n*pow(d,-1,p)%p if d%p else p;ap,good=table['traces'][t],table['good'][t]
            if good:totals[i]+=round((2-ap)/(p+1-ap)*log(p)*10**12);goods[i]+=1
    if any(r['score_units']!=s or r['good_primes']!=g for r,s,g in zip(rows,totals,goods)):raise ArithmeticError('full562-prime score recomputation differs')
    return table_hashes
def run():
    p=protocol();out=D/'result.json'
    if out.exists():raise FileExistsError('preserve retention output')
    cmd=[str(original.BINARY),str(ROOT/p['table']['path']),'4096','4096','512','0','1'];cap=capture(cmd,limits=Limits(p['seconds'],p['rss_bytes']),log_path=D/'scanner.log',separate_stderr=True,check=False);checkpoint(D/'raw.json',{'command':cmd,'stdout':cap.stdout,'stderr':cap.stderr,'supervision':cap.supervision})
    if cap.supervision['outcome']!='completed' or cap.supervision['returncode']!=0 or cap.stderr:raise ArithmeticError('benchmark incomplete/censored')
    rows,summary=original.parse(cap.stdout,SIGN);hashes=verify_rows(rows,summary,FAMILY,SIGN);checkpoint(out,{'schema':'elliptic-curves.r17-retention512-benchmark-result.v1','status':'PASS','protocol_hash':digest(p),'raw_sha256':cert.hashed(D/'raw.json'),'rows':rows,'summary':summary,'canonical_table_hashes':hashes,'all562_scores_recomputed':True,'old128_prefix_identical':True,'scanner_seconds':cap.supervision['wall_seconds'],'claim_boundary':'Fixed512 returned scores and old128 prefix verified; full10200039-address ranking still trusts the unchanged pinned scanner. No rank claim.'});print('RETENTION512 PASS',cap.supervision['wall_seconds'],'seconds',flush=True)
def replay():
    p=protocol();d=cert.read(D/'result.json');raw=cert.read(D/'raw.json')
    if d['status']!='PASS' or d['protocol_hash']!=digest(p) or d['raw_sha256']!=cert.hashed(D/'raw.json') or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('benchmark binding differs')
    rows,summary=original.parse(raw['stdout'],SIGN)
    if rows!=d['rows'] or summary!=d['summary'] or verify_rows(rows,summary,FAMILY,SIGN)!=d['canonical_table_hashes']:raise ArithmeticError('retention replay differs')
    print('REPLAYED RETENTION512 AND562-PRIME SCORES',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args();globals()[a.stage]()
