#!/usr/bin/env python3
"""Fixed wider retention on the existing six-family H4096 parameter box."""
import argparse,sys
from concurrent.futures import ThreadPoolExecutor,as_completed
from pathlib import Path
import benchmark_r17_retention512 as bench
import certify_compact_r17_candidates as cert
import select_compact_r17_wide as original
import compact_atlas_specialization as spec
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=bench.ROOT;LOCAL=bench.LOCAL;D=LOCAL/'compact-six-r17-retention512-v1'
def sources():
    return {**bench.sources(),str(Path(__file__).resolve().relative_to(ROOT)):cert.hashed(Path(__file__).resolve())}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve retention campaign')
    b=cert.read(bench.D/'result.json');check=cert.read(bench.D/'replay.supervisor.json');parent=cert.read(original.DIRECTORY/'protocol.json')
    if b['status']!='PASS' or b['scanner_seconds']>120 or check['outcome']!='completed' or check['returncode']!=0:raise ArithmeticError('retention benchmark gate incomplete')
    rows=[{'family':f['family'],'sign':s,'id':f['family']+('-negative' if s<0 else '-positive'),'table':parent['trace_tables'][f['family']][str(s)],'reuse_benchmark':f['family']=='103b2' and s==-1} for f in cert.read(spec.ATLAS)['families'] for s in (-1,1)]
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.compact-six-r17-retention512.v1','sources':sources(),'rows':rows,'height':4096,'keep_per_sign':512,'retained_total':6144,'benchmark_sha256':cert.hashed(bench.D/'result.json'),'benchmark_raw_sha256':cert.hashed(bench.D/'raw.json'),'scanner_seconds':120,'outer_seconds':600,'rss_bytes_per_worker':1073741824,'maximum_workers':2,'gate':'The predetermined512-per-sign benchmark passes with identical old128 prefix and exact562-prime scores. Eight discoveries from the saved discarded pool demonstrate an operational short-prime retention loss. Increase only retained width on the same122400468-address parameter box, preserving all existing search/certificate files.','policy':'Twelve fixed family/sign shards,512 per shard. Reuse the passed benchmark for103b2 negative and run the other eleven with the unchanged scanner and table bytes. Verify every returned562-prime score and all twelve old128 prefixes. Keep all6144 returned addresses, with no global short-prime merge cutoff. Scanner ties use absolute numerator; later extended-score selection uses signed numerator. No catalogue, public point or validation score is read.','prospective_scope':'After complete retention replay, freeze an extended-prime cohort excluding all1536 already scored addresses. At most4608 new extended scores, with the same disjoint validation band; at most24 equally exposed point candidates under a separate protocol. No point search or automatic next stage here.','boundaries':'Full population ranking trusts the unchanged pinned scanner. Returned scores and addresses are independently replayed; short-prime retention remains finite and does not prove global extended-score optimality or rank.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['benchmark_sha256']!=cert.hashed(bench.D/'result.json') or p['benchmark_raw_sha256']!=cert.hashed(bench.D/'raw.json'):raise ArithmeticError('frozen retention inputs differ')
    for r in p['rows']:
        if cert.hashed(ROOT/r['table']['path'])!=r['table']['sha256']:raise ArithmeticError('retained trace bytes changed')
    return p
def one(r,p):
    folder=D/r['id'];out=folder/'result.json'
    if out.exists():raise FileExistsError('preserve retention shard')
    if r['reuse_benchmark']:
        raw=cert.read(bench.D/'raw.json');checkpoint(folder/'raw.json',raw)
    else:
        cmd=[str(original.BINARY),str(ROOT/r['table']['path']),'4096','4096','512','0','1'];cap=capture(cmd,limits=Limits(p['scanner_seconds'],p['rss_bytes_per_worker']),log_path=folder/'scanner.log',separate_stderr=True,check=False);raw={'command':cmd,'stdout':cap.stdout,'stderr':cap.stderr,'supervision':cap.supervision};checkpoint(folder/'raw.json',raw)
    if raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0 or raw['stderr']:raise ArithmeticError('retention shard incomplete/censored')
    rows,summary=original.parse(raw['stdout'],r['sign']);hashes=bench.verify_rows(rows,summary,r['family'],r['sign']);result={'status':'PASS','id':r['id'],'family':r['family'],'sign':r['sign'],'rows':rows,'summary':summary,'canonical_table_hashes':hashes,'raw_sha256':cert.hashed(folder/'raw.json'),'protocol_hash':digest(p),'reused_benchmark':r['reuse_benchmark']};checkpoint(out,result);print('RETAINED512',r['id'],raw['supervision']['wall_seconds'],'seconds',flush=True);return result

def run():
    p=protocol();out=D/'result.json'
    if out.exists():raise FileExistsError('preserve retention aggregate')
    result={'status':'RUNNING','protocol_hash':digest(p),'shards':[]};checkpoint(out,result);got={}
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(one,r,p):i for i,r in enumerate(p['rows'])}
        for f in as_completed(pending):
            got[pending[f]]=f.result();result['shards']=[got[i] for i in sorted(got)];checkpoint(out,result)
    if sum(len(r['rows']) for r in result['shards'])!=6144:raise ArithmeticError('retention total differs')
    result['status']='COMPLETE_FIXED_RETENTION';checkpoint(out,result)
def replay():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='COMPLETE_FIXED_RETENTION' or d['protocol_hash']!=digest(p) or len(d['shards'])!=12:raise ArithmeticError('incomplete retention')
    for r,s in zip(p['rows'],d['shards']):
        raw=cert.read(D/r['id']/'raw.json');rows,summary=original.parse(raw['stdout'],r['sign'])
        expected=[str(original.BINARY),str(ROOT/r['table']['path']),'4096','4096','512','0','1']
        # Recorded commands are absolute original paths; compare their repository-relative operands.
        if Path(raw['command'][0]).name!=Path(expected[0]).name or Path(raw['command'][1]).parts[-3:]!=Path(expected[1]).parts[-3:] or raw['command'][2:]!=expected[2:] or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('retention execution differs')
        if r['id']!=s['id'] or rows!=s['rows'] or summary!=s['summary'] or s['raw_sha256']!=cert.hashed(D/r['id']/'raw.json') or bench.verify_rows(rows,summary,r['family'],r['sign'])!=s['canonical_table_hashes']:raise ArithmeticError('retention replay differs')
    print('REPLAYED6144 RETAINED SCORES AND ALL OLD128 PREFIXES',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','run','replay']);a=p.parse_args();globals()[a.stage]()
