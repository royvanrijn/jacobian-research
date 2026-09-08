#!/usr/bin/env python3
"""Fixed-output/cost gate for periodic score accumulation and one higher-height shard."""
import argparse
from pathlib import Path
from math import gcd,log
import certify_compact_r17_candidates as cert
import select_compact_r17_wide as old
import benchmark_r17_retention512 as prior
import compact_atlas_specialization as spec
from mod2_reduction_independence import _is_prime
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=old.ROOT;CAS=old.CAS;ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/periodic-nagao-scanner-benchmark-v1';OUT=ART/'periodic_nagao_scanner_benchmark_v1.json';CPP=CAS/'newfamily/scan_rational_nagao_tables_v2.cpp'
def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in (Path(__file__).resolve(),CPP,CAS/'newfamily/scan_rational_nagao_tables.cpp',old.BINARY,Path(old.__file__).resolve(),spec.ATLAS)}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve periodic benchmark protocol')
    gate=ART/'public_compact_parameter_heights_v1.json';h=cert.read(gate);parent=cert.read(old.DIRECTORY/'protocol.json');table=parent['trace_tables']['103b2']['-1']
    if h['status']!='PASS' or h['reported_rank_coverage']['28']['inside_box_counts']['16384']!=0:raise ArithmeticError('exact higher-height population diagnostic required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.periodic-nagao-scanner-benchmark.v1','sources':sources(),'table':table,'height_audit_sha256':cert.hashed(gate),'family':'103b2','sign':-1,'cases':[{'name':'old4096','numerator_bound':4096,'denominator_bound':4096,'keep':512,'shard':0,'shards':1},{'name':'higher32768-one-shard','numerator_bound':32768,'denominator_bound':32768,'keep':512,'shard':0,'shards':64}],'compiler':'g++ -O3 -std=c++17','build_seconds':60,'seconds_per_call':120,'rss_bytes':1073741824,'maximum_workers':1,'gate':'The exact retrospective coordinate audit places all eight mapped public>=28 successes outside16384; their sample has unknown historical exposure and does not predict success density. Before a new higher-height population, test a score-engine cost improvement on the first atlas family and negative sign, retaining identical512 output at4096 and comparing both engines on one fixed denominator shard at32768. No public parameter, point or score enters this benchmark selection.','proof':'For fixed d not divisible by p, the local symbol for n/d is periodic in n mod p. Build one length-p cycle using d inverse and add identical copies to exact signed64-bit score/int good-count arrays. If p divides d, all n use the same infinity symbol. Apply the unchanged gcd filter and deterministic heap order only after all prime contributions. Held scores remain computed by the original survivor routine.','boundaries':'Two finite fixed parameter populations and two unchanged score policies; no new finalist promoted, no full32768 scan, no point search or automatic campaign. Exact output agreement and returned-score replay do not replace trust in the pinned enumeration workers.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or cert.hashed(ROOT/p['table']['path'])!=p['table']['sha256']:raise ArithmeticError('benchmark bindings changed')
    return p
def build():
    p=protocol();binary=D/'scanner-v2'
    if binary.exists():raise FileExistsError('preserve prototype binary')
    cmd=['g++','-O3','-std=c++17',str(CPP),'-o',str(binary)];c=capture(cmd,limits=Limits(60,p['rss_bytes']),log_path=D/'build.log',check=True);checkpoint(D/'build.json',{'command':cmd,'supervision':c.supervision,'binary_sha256':cert.hashed(binary)})
def exact_scores(rows,p):
    scores=[0]*len(rows);goods=[0]*len(rows);hashes={};family=p['family'];f=next(r for r in cert.read(spec.ATLAS)['families'] if r['family']==family);model={k:f[k]+['0']*(n-len(f[k])) for k,n in [('A_coefficients_low_to_high',9),('B_coefficients_low_to_high',13)]}
    for q in (q for q in range(5,4094) if _is_prime(q)):
        path=old.PARENT/family/'trace-tables'/f'{q}.json';t=cert.read(path)
        if t['input']!={'family':family,'model_hash':digest(model),'prime':q}:raise ArithmeticError('canonical score table differs')
        hashes[str(path.relative_to(ROOT))]=cert.hashed(path)
        for i,r in enumerate(rows):
            x=r['numerator']*pow(r['denominator'],-1,q)%q if r['denominator']%q else q
            if t['good'][x]:ap=t['traces'][x];scores[i]+=round((2-ap)/(q+1-ap)*log(q)*10**12);goods[i]+=1
    if any((r['score_units'],r['good_primes'])!=(s,g) for r,s,g in zip(rows,scores,goods)):raise ArithmeticError('exact returned score differs')
    return hashes
def checked(case,p):
    values={}
    for version in ('v1','v2'):
        raw=cert.read(D/case['name']/(version+'.json'))
        if raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('censored/failed benchmark')
        rows,summary=old.parse(raw['stdout'],p['sign']);values[version]=(rows,summary,raw)
    rows,summary,raw=values['v2']
    if values['v1'][:2]!=values['v2'][:2] or values['v1'][2]['stdout']!=raw['stdout']:raise ArithmeticError('complete ordered old/new output differs')
    if summary[:3]!=[case['numerator_bound'],case['denominator_bound'],512] or summary[3]!=summary[4] or summary[5]!=512:raise ArithmeticError('complete population framing differs')
    if len(rows)!=512 or len({r['parameter'] for r in rows})!=512 or any(gcd(abs(r['numerator']),r['denominator'])!=1 or not 1<=abs(r['numerator'])<=case['numerator_bound'] or not 1<=r['denominator']<=case['denominator_bound'] or (r['denominator']-1)%case['shards']!=case['shard'] or r['numerator']>=0 for r in rows):raise ArithmeticError('fixed primitive roster differs')
    if rows!=sorted(rows,key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],abs(r['numerator']))):raise ArithmeticError('ordered scores differ')
    if case['name']=='old4096' and rows!=cert.read(prior.D/'result.json')['rows']:raise ArithmeticError('full old512 prefix differs')
    hashes=exact_scores(rows,p)
    return {'case':case,'rows':rows,'summary':summary,'canonical_table_hashes':hashes,'wall_seconds':{v:values[v][2]['supervision']['wall_seconds'] for v in values},'raw_hashes':{v:cert.hashed(D/case['name']/(v+'.json')) for v in values},'speedup':values['v1'][2]['supervision']['wall_seconds']/raw['supervision']['wall_seconds']}
def run():
    p=protocol();out=D/'result.json'
    if out.exists():raise FileExistsError('preserve benchmark outcomes')
    data={'status':'RUNNING','protocol_hash':digest(p),'cases':[]};checkpoint(out,data)
    for case in p['cases']:
        for version,binary in [('v1',old.BINARY),('v2',D/'scanner-v2')]:
            cmd=[str(binary),str(ROOT/p['table']['path']),*[str(case[k]) for k in ('numerator_bound','denominator_bound','keep','shard','shards')]];c=capture(cmd,limits=Limits(120,p['rss_bytes']),log_path=D/case['name']/(version+'.log'),separate_stderr=True,check=False);checkpoint(D/case['name']/(version+'.json'),{'command':cmd,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
        row=checked(case,p);data['cases'].append(row);checkpoint(out,data);print('PERIODIC NAGAO',case['name'],row['wall_seconds'],row['speedup'],flush=True)
    data['status']='PASS';checkpoint(out,data)
def replay():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='PASS' or d['protocol_hash']!=digest(p) or d['cases']!=[checked(c,p) for c in p['cases']]:raise ArithmeticError('full benchmark replay differs')
    result={'schema':'elliptic-curves.periodic-nagao-scanner-benchmark-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'result_sha256':cert.hashed(D/'result.json'),'cases':[{'name':r['case']['name'],'summary':r['summary'],'wall_seconds':r['wall_seconds'],'speedup':r['speedup']} for r in d['cases']],'claim_boundary':p['boundaries']}
    if OUT.exists():
        if cert.read(OUT)!=result:raise ArithmeticError('benchmark certificate differs')
    else:checkpoint(OUT,result)
    print('REPLAYED BOTH PERIODIC SCANNER BOXES AND1024 EXACT562-PRIME SCORES',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','build','run','replay']);a=p.parse_args();globals()[a.stage]()
