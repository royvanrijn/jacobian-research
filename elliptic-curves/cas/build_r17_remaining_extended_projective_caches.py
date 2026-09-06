#!/usr/bin/env python3
"""Five complete target-free R17 projective caches for outer-region selection."""
import argparse,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import benchmark_r17_remaining_extended_projective_tables as benchmark
import benchmark_extended_projective_trace_cache_v2 as engine
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=benchmark.ROOT;CAS=benchmark.CAS;ART=benchmark.ART
D=ROOT/'artifacts/local/elliptic-curves/r17_remaining-extended-projective-caches-v1'
OUT=ART/'r17_remaining_extended_projective_caches_v1.json'
def sources():
    paths=[Path(__file__).resolve(),benchmark.OUT,benchmark.D/'protocol.json',benchmark.D/'controller/ledger.json']
    return {**benchmark.sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}}
def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve five-cache protocol')
    ledger=cert.read(benchmark.D/'controller/ledger.json')
    if ledger['status']!='PASS' or [r['name'] for r in ledger['rows']]!=['run','check'] or any(r['status']!='PASS' for r in ledger['rows']):raise ArithmeticError('complete fifteen-table benchmark and read-only replay required')
    p=benchmark.protocol();gate=cert.read(benchmark.OUT)
    if gate['status']!='PASS' or not gate['cost_gate_passed'] or gate['direct_character_sums']!=75:raise ArithmeticError('balanced exact and cost gate required')
    reused={r['family']+':'+str(r['prime']):{k:r[k] for k in ('raw_sha256','table_sha256')} for r in gate['rows']}
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.r17_remaining-extended-projective-caches.v1','sources':sources(),'families':p['families'],'primes':p['full_prime_roster'],'gp_sha256':p['gp_sha256'],'gp_seconds_per_case':20,'rss_bytes':536870912,'maximum_workers':5,'checkpoint_block_size':80,'build_wall_seconds':7200,'replay_wall_seconds':3600,'reused':reused,'gate':'All five remaining R17 families pass the same fixed three-prime projective-table cost gate and75 independent character sums. All score-range prime-scaling exclusions through131071 independently pass. The user requests broader earlier selection in new territory; the completed broader MW16 outcome supports a separately bounded extension to the five R17 fibrations still lacking complete extended projective tables. Height-independent finite-field tables allow later untouched parameter populations to receive all selection primes before short-score truncation. Equal treatment of all five missing families;11952 already has its complete cache and is excluded. No known target enters generation.','scope':'Exactly the five existing R17 atlas models other than11952, every prime4099through32749 and allp+1 raw projective residues. Reuse exactly15 benchmark tables and make at most14725 new bounded PARI calls, five workers and80-case checkpoints. Verify all frames/discriminants/Hasse bounds and five independent character sums per table. Replay all tables without new GP calls. No retry or larger prime range. The15 original benchmark raw/table pairs are reused once with their exact hashes. This cache build includes no parameter scan, point search, new rank inference, restoration of bad raw-model reduction or catalogue input. A subsequent fresh-parameter trial requires its own frozen protocol after this replay passes.'})
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(engine.scalar.GP):raise ArithmeticError('frozen cache sources changed')
    return p
def table(f, prime, p, create):
    folder = D/f['family']/str(prime)
    rawpath = folder/'raw.json'
    code = engine.program(f['model'], prime)
    if create:
        if rawpath.exists():
            raise FileExistsError('preserve attempted table')
        result = capture([str(engine.scalar.GP),'-q','-s','256000000'],
                         input_text=code, limits=Limits(p['gp_seconds_per_case'],p['rss_bytes']),
                         log_path=folder/'gp.log', separate_stderr=True, check=False)
        checkpoint(rawpath, {'program':code,'stdout':result.stdout,'stderr':result.stderr,
                             'supervision':result.supervision})
    raw = cert.read(rawpath)
    if raw['program'] != code or raw['stderr'] or raw['supervision']['outcome'] != 'completed' or raw['supervision']['returncode'] != 0:
        raise ArithmeticError('table failed or censored')
    lines = raw['stdout'].splitlines()
    if len(lines) != prime+3 or lines[-1] != 'DONE' or not lines[-2].startswith('MS|'):
        raise ArithmeticError('incomplete projective table')
    a,b = ([int(c)%prime for c in f['model'][k]] for k in
           ('A_coefficients_low_to_high','B_coefficients_low_to_high'))
    values = []
    for i,line in enumerate(lines[:-2]):
        aa = a[-1] if i == prime else engine.evaluate(a,i,prime)
        bb = b[-1] if i == prime else engine.evaluate(b,i,prime)
        bad = (4*aa**3+27*bb**2)%prime == 0
        fields = line.split('|')
        if len(fields) != (2 if bad else 3) or fields[0] != ('B' if bad else 'T') or int(fields[1]) != i:
            raise ArithmeticError('residue or discriminant mismatch')
        value = None if bad else int(fields[2])
        if value is not None and value*value > 4*prime:
            raise ArithmeticError('Hasse bound violated')
        values.append(value)
    direct = []
    for i in [0,prime//3,2*prime//3,prime-1,prime]:
        aa = a[-1] if i == prime else engine.evaluate(a,i,prime)
        bb = b[-1] if i == prime else engine.evaluate(b,i,prime)
        value = engine.scalar.direct([0,0,0,aa,bb],prime)
        if value != values[i]:
            raise ArithmeticError('independent character sum differs')
        direct.append([i,value])
    saved = {'input':{'family':f['family'],'model_hash':f['model_hash'],'prime':prime},
             'traces':[v or 0 for v in values], 'good':[v is not None for v in values]}
    path = folder/'table.json'
    if create:
        if path.exists(): raise FileExistsError('preserve prime table')
        checkpoint(path,saved)
    elif cert.read(path) != saved:
        raise ArithmeticError('table replay differs')
    return {'family':f['family'],'prime':prime,'projective_rows':prime+1,
            'wall_seconds':raw['supervision']['wall_seconds'],'direct_checks':direct,
            'raw_sha256':cert.hashed(rawpath),'table_sha256':cert.hashed(path)}

def one(job,p):
    f,prime=job;key=f['family']+':'+str(prime)
    if key in p['reused']:
        folder=D/f['family']/str(prime);folder.mkdir(parents=True,exist_ok=True)
        for name,keyhash in [('raw.json','raw_sha256'),('table.json','table_sha256')]:
            source=benchmark.D/f['family']/str(prime)/name;target=folder/name
            if target.exists():raise FileExistsError('preserve benchmark copy')
            if cert.hashed(source)!=p['reused'][key][keyhash]:raise ArithmeticError('benchmark table changed')
            target.write_bytes(source.read_bytes())
        return table(f,prime,p,False)
    return table(f,prime,p,True)
def jobs(p):return [(f,q) for q in p['primes'] for f in p['families']]
def result(p,rows):
    if [(r['family'],r['prime']) for r in rows]!=[(f['family'],q) for f,q in jobs(p)]:raise ArithmeticError('full five-family product differs')
    return {'schema':'elliptic-curves.r17_remaining-extended-projective-caches-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows,'families':len(p['families']),'primes_per_family':len(p['primes']),'projective_rows':sum(r['projective_rows'] for r in rows),'direct_character_sums':5*len(rows),'reused_tables':len(p['reused']),'claim_boundary':p['scope']}
def run():
    p=protocol();path=D/'ledger.json'
    if path.exists() or OUT.exists():raise FileExistsError('preserve cache attempt')
    ledger={'status':'RUNNING','protocol_hash':digest(p),'rows':[]};checkpoint(path,ledger);start=time.monotonic();roster=jobs(p)
    try:
        with ThreadPoolExecutor(max_workers=p['maximum_workers']) as pool:
            for i in range(0,len(roster),p['checkpoint_block_size']):
                if time.monotonic()-start>p['build_wall_seconds']:raise TimeoutError('fixed build deadline')
                ledger['rows']+=list(pool.map(lambda j:one(j,p),roster[i:i+p['checkpoint_block_size']]))
                checkpoint(path,ledger)
                if len(ledger['rows'])%400==0:print('FIVE R17 CACHE TABLES',len(ledger['rows']),'OF',len(roster),flush=True)
        checkpoint(OUT,result(p,ledger['rows']));ledger['status']='PASS';checkpoint(path,ledger)
    except Exception as exc:
        ledger.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(path,ledger);raise

def check():
    p=protocol();ledger=cert.read(D/'ledger.json');start=time.monotonic()
    if ledger['status']!='PASS' or ledger['protocol_hash']!=digest(p):raise ArithmeticError('terminal full cache required')
    rows=[]
    for f,q in jobs(p):
        if time.monotonic()-start>p['replay_wall_seconds']:raise TimeoutError('fixed replay deadline')
        rows.append(table(f,q,p,False))
    if rows!=ledger['rows'] or result(p,rows)!=cert.read(OUT):raise ArithmeticError('full five-cache replay differs')
    checkpoint(D/'replay.json',{'status':'PASS','sources':sources(),'certificate_sha256':cert.hashed(OUT),'tables':len(rows),'direct_character_sums':len(rows)*5})
    print('FIVE COMPLETE R17 CACHES REPLAY PASS',len(rows),flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','run','check']);a=parser.parse_args();globals()[a.stage]()
