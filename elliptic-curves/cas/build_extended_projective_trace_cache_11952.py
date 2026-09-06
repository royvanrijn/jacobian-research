#!/usr/bin/env python3
"""One complete11952 projective cache at selection primes4099through32749."""
import argparse,json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import certify_compact_r17_candidates as cert
import benchmark_extended_projective_trace_cache_v2 as benchmark
import benchmark_r17_extended_prime_traces as scalar
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=benchmark.ROOT;ART=benchmark.ART;CAS=ROOT/'elliptic-curves/cas'
D=ROOT/'artifacts/local/elliptic-curves/extended-projective-trace-cache-11952-v1'
OUT=ART/'extended_projective_trace_cache_11952_v1.json'

def sources():
    return {**benchmark.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),benchmark.OUT,benchmark.D/'protocol.json',benchmark.D/'check.supervisor.json']}}

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve complete cache protocol')
    p=benchmark.protocol();b=cert.read(benchmark.OUT);check=cert.read(benchmark.D/'check.supervisor.json')
    if b['status']!='PASS' or not b['cost_gate_passed'] or b['character_sum_checks']!=40 or check['outcome']!='completed' or check['returncode']!=0:raise ArithmeticError('eight-prime exact and cost gates required')
    reused={str(r['prime']):{'raw_sha256':r['raw_sha256'],'table_sha256':r['table_sha256']} for r in b['rows']}
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.extended-projective-trace-cache-11952.v1','sources':sources(),'family':'11952','model':p['model'],'model_hash':p['model_hash'],'primes':p['full_selection_prime_roster'],'projective_rows':p['full_selection_projective_rows'],'gp_sha256':p['gp_sha256'],'gp_seconds_per_prime':20,'rss_bytes':536870912,'maximum_workers':2,'checkpoint_block_size':16,'build_wall_seconds':3600,'replay_wall_seconds':1200,'reused_benchmark_tables':reused,'gate':'The exact known29 retention contrast is435th short versus2nd extended among the same967 saved candidates. Eight deterministic projective-table benchmarks through32749 pass all40 independent character sums, all residue/discriminant/Hasse frames and the frozen1800-second single-worker projection gate. This complete one-family cache removes repeated finite-field work when a later experiment retains more short-stage candidates. The projection is not a runtime guarantee; enforce the actual finite limits.','scope':'Exactly the fixed11952 raw homogeneous projective fibres at every prime4099through32749, with allp+1 residues. Reuse the eight immutable benchmark tables; every remaining prime gets one bounded PARI call. Check the complete residue/discriminant/Hasse frame and five independent character sums per prime, preserve every raw call and table, and replay the whole cache without executing new GP calls. Stop after the current sixteen-prime block on any failure or censoring; no retry, broader prime range or second family. No parameter scan, candidate selection, point search, inferred rank, good-reduction restoration or new curve is part of this cache build.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(scalar.GP):raise ArithmeticError('frozen complete cache inputs changed')
    return p

def evaluate(c,t,p):
    z=0
    for a in reversed(c):z=(z*t+a)%p
    return z

def program(model,p):
    a=[int(v)%p for v in model['A_coefficients_low_to_high']];b=[int(v)%p for v in model['B_coefficients_low_to_high']]
    return f'A={a};B={b};p={p};ev(v,t)={{my(z=0);forstep(j=#v,1,-1,z=(z*t+v[j])%p);z}};\ngettime();for(i=0,p,a=if(i==p,A[#A],ev(A,i));b=if(i==p,B[#B],ev(B,i));if((4*a^3+27*b^2)%p==0,print("B|",i),print("T|",i,"|",p+1-ellcard(ellinit([Mod(a,p),Mod(b,p)])))));print("MS|",gettime());print("DONE");quit\n'

def table(prime,p,create):
    folder=D/str(prime);path=folder/'raw.json';code=program(p['model'],prime)
    if create:
        if path.exists():raise FileExistsError('preserve prime-table benchmark call')
        c=capture([str(scalar.GP),'-q','-s','256000000'],input_text=code,limits=Limits(p['gp_seconds_per_prime'],p['rss_bytes']),log_path=folder/'gp.log',separate_stderr=True,check=False)
        checkpoint(path,{'program':code,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    if raw['program']!=code or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('table call failed or censored')
    lines=raw['stdout'].splitlines();a=[int(v)%prime for v in p['model']['A_coefficients_low_to_high']];b=[int(v)%prime for v in p['model']['B_coefficients_low_to_high']];values=[]
    if len(lines)!=prime+3 or lines[-1]!='DONE' or not lines[-2].startswith('MS|'):raise ArithmeticError('complete projective table framing differs')
    for i,line in enumerate(lines[:-2]):
        aa=a[-1] if i==prime else evaluate(a,i,prime);bb=b[-1] if i==prime else evaluate(b,i,prime);bad=(4*aa**3+27*bb**2)%prime==0;v=line.split('|')
        if len(v)!=(2 if bad else 3) or v[0]!=('B' if bad else 'T') or int(v[1])!=i:raise ArithmeticError('residue/discriminant framing differs')
        value=None if bad else int(v[2])
        if value is not None and value*value>4*prime:raise ArithmeticError('projective Hasse bound differs')
        values.append(value)
    checks=[]
    for i in [0,prime//3,2*prime//3,prime-1,prime]:
        aa=a[-1] if i==prime else evaluate(a,i,prime);bb=b[-1] if i==prime else evaluate(b,i,prime);direct=scalar.direct([0,0,0,aa,bb],prime)
        if direct!=values[i]:raise ArithmeticError('independent table character sum differs')
        checks.append([i,direct])
    saved={'input':{'family':p['family'],'model_hash':p['model_hash'],'prime':prime},'traces':[v or 0 for v in values],'good':[v is not None for v in values]};tp=folder/'table.json'
    if create:
        if tp.exists():raise FileExistsError('preserve benchmark table')
        checkpoint(tp,saved)
    elif cert.read(tp)!=saved:raise ArithmeticError('saved projective table differs')
    return {'prime':prime,'projective_rows':prime+1,'cpu_ms':int(lines[-2][3:]),'wall_seconds':raw['supervision']['wall_seconds'],'raw_sha256':cert.hashed(path),'table_sha256':cert.hashed(tp),'direct_checks':checks}

def one(prime,p):
    if str(prime) in p['reused_benchmark_tables']:
        folder=D/str(prime);folder.mkdir(parents=True,exist_ok=True)
        for name,key in [('raw.json','raw_sha256'),('table.json','table_sha256')]:
            original=benchmark.D/str(prime)/name;path=folder/name
            if path.exists():raise FileExistsError('preserve cached benchmark copy')
            if cert.hashed(original)!=p['reused_benchmark_tables'][str(prime)][key]:raise ArithmeticError('benchmark table changed')
            path.write_bytes(original.read_bytes())
        return table(prime,p,False)
    return table(prime,p,True)

def result(p,rows):
    if [r['prime'] for r in rows]!=p['primes'] or sum(r['projective_rows'] for r in rows)!=p['projective_rows']:raise ArithmeticError('complete selection-prime product differs')
    return {'schema':'elliptic-curves.extended-projective-trace-cache-11952-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows,'projective_rows':p['projective_rows'],'selection_primes':len(rows),'character_sum_checks':5*len(rows),'reused_benchmark_tables':len(p['reused_benchmark_tables']),'claim_boundary':p['scope']}

def run():
    p=protocol();path=D/'ledger.json'
    if path.exists() or OUT.exists():raise FileExistsError('preserve one full cache attempt')
    ledger={'status':'RUNNING','protocol_hash':digest(p),'rows':[]};checkpoint(path,ledger)
    try:
        with ThreadPoolExecutor(max_workers=p['maximum_workers']) as pool:
            for start in range(0,len(p['primes']),p['checkpoint_block_size']):
                ledger['rows']+=list(pool.map(lambda q:one(q,p),p['primes'][start:start+p['checkpoint_block_size']]))
                checkpoint(path,ledger)
                if len(ledger['rows'])%256==0:print('COMPLETE CACHE PROGRESS',len(ledger['rows']),'of',len(p['primes']),flush=True)
        certificate=result(p,ledger['rows']);checkpoint(OUT,certificate);ledger['status']='PASS';checkpoint(path,ledger)
    except Exception as error:
        ledger['status']='FAILED_OR_CENSORED';ledger['reason']=str(error);checkpoint(path,ledger);raise

def check():
    p=protocol();ledger=cert.read(D/'ledger.json')
    if ledger['status']!='PASS' or ledger['protocol_hash']!=digest(p):raise ArithmeticError('complete cache terminal ledger required')
    rows=[table(q,p,False) for q in p['primes']];r=result(p,rows)
    if ledger['rows']!=rows or cert.read(OUT)!=r:raise ArithmeticError('complete cache replay differs')
    print('REPLAYED COMPLETE11952 CACHE',len(rows),'PRIMES',p['projective_rows'],'PROJECTIVE RESIDUES',flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','run','check']);v=a.parse_args();globals()[v.stage]()
