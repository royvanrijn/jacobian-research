#!/usr/bin/env python3
"""Untouched higher MW16 slices scored through32749 before any truncation."""
import argparse
from pathlib import Path
from hashlib import sha256
from concurrent.futures import ThreadPoolExecutor,as_completed
import benchmark_joint_11952_annulus as engine
import scan_mw16_outer_bands as old
import encode_mw16_joint_caches as binary
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture,Limits
ROOT=old.ROOT;CAS=old.CAS;ART=old.ART
D=ROOT/'artifacts/local/elliptic-curves/joint-mw16-higher-annuli-v1'
SALT='target-free-mw16-complete3510-higher-annuli-v1';KEEP=4096

def sources():
    paths=[Path(__file__).resolve(),Path(old.__file__),old.D/'protocol.json',
           Path(binary.__file__),binary.OUT,binary.D/'score-replay.json',
           ROOT/'artifacts/local/elliptic-curves/mw16-joint-cache-controller-v1/ledger.json',
           engine.BINARY,engine.CPP,engine.TEST,engine.OUT,engine.D/'ledger.json',
           engine.reader.BINARY,engine.reader.OUT,CAS/'verify_periodic_nagao_scanner.py',
           CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def roster():
    rows=[]
    for r in cert.read(old.D/'protocol.json')['rows']:
        if r['band'] not in (2,3):continue
        modulus=r['shards'];choices=[v for v in range(r['shard']%2,modulus,2) if v!=r['shard']]
        h=sha256(f"{SALT}|{r['band']}|{r['family']}|{r['sign']}".encode()).hexdigest();shard=choices[int(h,16)%len(choices)]
        count=old.counts.population(r['outer'],r['outer'],shard,modulus)-old.counts.population(r['inner'],r['inner'],shard,modulus)
        rows.append({'id':r['id'],'family':r['family'],'band':r['band'],'sign':r['sign'],
                     'N':r['outer'],'M':r['outer'],'inner':r['inner'],'shard':shard,'shards':modulus,
                     'excluded_previous_shard':r['shard'],'sha256_choice':h,'primitive_population':count})
    return rows

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve higher full-score scanner')
    proof=cert.read(binary.OUT);gate=cert.read(ROOT/'artifacts/local/elliptic-curves/mw16-joint-cache-controller-v1/ledger.json')
    compiled=cert.read(engine.OUT)
    if proof['status']!='PASS' or proof['candidates_checked']!=40960 or gate['status']!='PASS':raise ArithmeticError('all five cache proofs and40960 scalar agreements required')
    if compiled['status']!='PASS' or cert.hashed(engine.BINARY)!=compiled['binary_sha256'] or cert.read(engine.D/'ledger.json')['status']!='PASS':raise ArithmeticError('exact full-score binary gate required')
    rows=roster()
    if len(rows)!=20 or sum(r['shard']%2 for r in rows)!=10 or any(r['shard']==r['excluded_previous_shard'] for r in rows):raise ArithmeticError('balanced fresh20 required')
    bindings={}
    for f in {r['family'] for r in rows}:
        for label in ('short','extended'):
            path=binary.cache(f,label);meta=next(r for r in proof['encoding']['rows'] if r['family']==f and r['label']==label)
            if cert.hashed(path)!=meta['sha256']:raise ArithmeticError('encoded cache changed')
            bindings[str(path.relative_to(ROOT))]=meta['sha256']
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.joint-mw16-higher-annuli.v1','sources':sources(),
        'rows':rows,'cache_bindings':bindings,'keep_per_slice':KEEP,'prime_count':3510,
        'maximum_workers':2,'seconds_per_call':120,'first_slice_per_band_cost_gate_seconds':45,
        'run_wall_seconds':1200,'replay_wall_seconds':600,'rss_bytes':1610612736,
        'total_primitive_population':sum(r['primitive_population'] for r in rows),
        'gate':'Fresh outer MW16 discoveries demonstrate loss beyond a signed512 short prefix, while full-score11952 enumeration measures substantial loss even beyond4096. All five complete projective caches and40960 existing scalar agreements now pass. Revisit the higher two parameter bands on different denominator slices with all3510 primes scored before truncation. Earlier shallow-retention low yields were finite observations, not exclusions of these bands.',
        'population':'Exactly twenty fresh signed denominator slices: five families times two signs times the two annuli16384<H<=65536 and65536<H<=262144. SHA256 chooses a same-parity denominator residue modulo256 or4096, excluding the corresponding previous outer-band residue. Every address is outside compact4096, outside the fresh4096..16384 trial and disjoint from all previously scanned outer-band addresses of that family.',
        'selection':'Every primitive parameter receives all3510 quantized selection primes before heap admission. Retain4096 per signed slice by units,good count,denominator,absolute numerator. There is no short-score prefilter. Public record equations,parameters,points,ranks,j-invariants and jump labels remain outside selection and execution.',
        'validation':'Complete signed actual-modulus frames and top7 ordering for all20 cases against the independent retained-list reader on both family caches. Check every retained score component and exact primitive population. One first slice per band has a45-second cost gate and is reused once. No retry, adaptive expansion or replacement slices.',
        'future_scope':'Only after complete scan replay, a separate protocol may select1024 distinct prospective equations per band/family from the8192-address signed pool,10240 total, score with fresh scalar traces through65521 and independently check the cached first extension. Then select six per band/family,60 total, with disjoint65537..131071 validation. A later generic16-only point protocol may use43 charts at125000/10seconds, all maps before points. No catalogue or prior-record exclusion enters MW16 selection; post-terminal exact novelty comparisons remain separate.',
        'boundaries':'Finite full-score selection inside20 declared slices, not complete annulus coverage, a rank predictor, a rank upper bound or universal novelty. No original point search under this scanner.'})
    print('FROZEN20 FRESH HIGHER MW16 SLICES',sum(r['primitive_population'] for r in rows),'FULL3510-PRIME ADDRESSES',flush=True)

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['rows']!=roster():raise ArithmeticError('frozen higher MW16 scanner changed')
    if any(cert.hashed(ROOT/n)!=h for n,h in p['cache_bindings'].items()):raise ArithmeticError('frozen family binary cache changed')
    return p

def execute(name,cmd,check=False):
    path=D/(name+'.json')
    if not check:
        if path.exists():raise FileExistsError('preserve one fresh joint call')
        c=capture(cmd,limits=Limits(120,1610612736),log_path=D/(name+'.log'),separate_stderr=True,check=False)
        checkpoint(path,{'command':cmd,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path)
    if raw['command']!=cmd or raw['supervision']['command']!=cmd or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('fresh joint call failed/censored')
    return raw

def scan_case(c,keep,name,check=False):
    values=[c['sign'],c['N'],c['M'],keep,c['shard'],c['shards'],c['inner'],3510]
    raw=execute(name,[str(engine.BINARY),str(binary.cache(c['family'],'short')),str(binary.cache(c['family'],'extended')),*map(str,values)],check);lines=raw['stdout'].splitlines()
    count=old.counts.population(c['N'],c['M'],c['shard'],c['shards'])-old.counts.population(min(c['N'],c['inner']),min(c['M'],c['inner']),c['shard'],c['shards']);kept=min(count,keep)
    if len(lines)!=kept+3 or lines[:2]!=['JOINT_NAGAO_ANNULUS_V1','P '+' '.join(map(str,values))] or lines[-1]!=f'S {count} {kept}':raise ArithmeticError('fresh joint population framing differs')
    rows=[]
    for line in lines[2:-1]:
        v=line.split()
        if len(v)!=5 or v[0]!='C':raise ArithmeticError('candidate row framing differs')
        n,d,u,g=map(int,v[1:])
        if n*c['sign']<=0 or engine.gcd(n,d)!=1 or not 1<=abs(n)<=c['N'] or not 1<=d<=c['M'] or not c['inner']<max(abs(n),d) or (d-1)%c['shards']!=c['shard'] or not 0<=g<=3510:raise ArithmeticError('nonprimitive/out-of-frame candidate')
        rows.append({'numerator':n,'denominator':d,'parameter':str(cert.F(n,d)),'combined_selection_units':u,'combined_good':g})
    if len({r['parameter'] for r in rows})!=len(rows) or rows!=sorted(rows,key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],abs(r['numerator']))):raise ArithmeticError('exact slice ordering/uniqueness differs')
    return rows,raw,count

def components(rows,family,name,check=False):
    path=D/(name+'-candidates.txt');text='R17-CANDIDATES-V1 '+str(len(rows))+'\n'+''.join(f"{r['numerator']} {r['denominator']}\n" for r in rows)
    if check:
        if path.read_text()!=text:raise ArithmeticError('component reader roster changed')
    else:
        if path.exists():raise FileExistsError('preserve component reader input')
        path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text)
    sums=[]
    for label,count in [('short',562),('extended',2948)]:
        raw=execute(name+'-'+label,[str(engine.reader.BINARY),str(binary.cache(family,label)),str(path),'262144'],check);lines=raw['stdout'].splitlines();values=[]
        if len(lines)!=len(rows)+1 or lines[-1]!=f'S {len(rows)} {count}':raise ArithmeticError('independent component frame differs')
        for i,line in enumerate(lines[:-1]):
            v=line.split()
            if len(v)!=4 or v[:2]!=['R',str(i)]:raise ArithmeticError('independent component index differs')
            values.append(tuple(map(int,v[2:])))
        sums.append(values)
    result=[]
    for row,a,b in zip(rows,*sums):
        if 'combined_selection_units' in row and (row['combined_selection_units'],row['combined_good'])!=(a[0]+b[0],a[1]+b[1]):raise ArithmeticError('complete combined score differs from independent components')
        result.append({**row,'score_units':a[0],'good_primes':a[1],'extension_selection_units':b[0],'extension_good':b[1],'combined_selection_units':a[0]+b[0],'combined_good':a[1]+b[1]})
    return result

def benchmark(check=False):
    p=protocol();records=[]
    for row in p['rows']:
        c={**row,'N':31,'M':2*row['shards'],'inner':row['shards']};name='benchmark/'+row['id']
        bare=[{'numerator':c['sign']*n,'denominator':d,'parameter':str(cert.F(c['sign']*n,d))} for d in range(c['shard']+1,c['M']+1,c['shards']) for n in range(1,32) if engine.gcd(n,d)==1 and max(n,d)>c['inner']]
        expected=components(bare,row['family'],name+'-reference',check);expected.sort(key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],abs(r['numerator'])))
        keys=['numerator','denominator','parameter','combined_selection_units','combined_good'];target=[{k:r[k] for k in keys} for r in expected]
        for keep in (1000,7):
            got,raw,count=scan_case(c,keep,name+'-k'+str(keep),check)
            if got!=target[:keep] or count!=len(target):raise ArithmeticError('complete new signed frame/top7 differs')
        records.append({'id':row['id'],'primitive_population':len(target),'status':'PASS'})
    result={'status':'PASS','protocol_sha256':cert.hashed(D/'protocol.json'),'rows':records}
    if check:
        if cert.read(D/'benchmark.json')!=result:raise ArithmeticError('fresh joint frame replay differs')
    else:
        if (D/'benchmark.json').exists():raise FileExistsError('preserve fresh joint regression')
        checkpoint(D/'benchmark.json',result)
    print('ALL20 FRESH MW16 FULL-SCORE FRAMES PASS',flush=True)

def checked(row,check=False):
    rows,raw,count=scan_case(row,KEEP,row['id']+'/scan',check);verified=components(rows,row['family'],row['id']+'/components',check)
    if count!=row['primitive_population'] or len(rows)!=KEEP:raise ArithmeticError('full retained slice/count differs')
    return {**row,'rows':verified,'raw_sha256':cert.hashed(D/row['id']/'scan.json'),'wall_seconds':raw['supervision']['wall_seconds']}

def run():
    p=protocol();benchmark(True)
    if (D/'result.json').exists():raise FileExistsError('preserve higher joint population')
    data={'status':'RUNNING','protocol_sha256':cert.hashed(D/'protocol.json'),'shards':[]};checkpoint(D/'result.json',data);done={}
    try:
        for band in (2,3):
            i=next(i for i,r in enumerate(p['rows']) if r['band']==band)
            done[i]=checked(p['rows'][i]);data['shards']=[done[j] for j in sorted(done)];checkpoint(D/'result.json',data)
            if done[i]['wall_seconds']>p['first_slice_per_band_cost_gate_seconds']:raise ArithmeticError('first higher-band slice cost gate failed')
        with ThreadPoolExecutor(max_workers=2) as pool:
            pending={pool.submit(checked,row):i for i,row in enumerate(p['rows']) if i not in done}
            for future in as_completed(pending):
                i=pending[future];done[i]=future.result();data['shards']=[done[j] for j in sorted(done)];checkpoint(D/'result.json',data)
                print('FRESH HIGHER MW16',done[i]['id'],done[i]['wall_seconds'],flush=True)
        data['status']='COMPLETE_FIXED_STRATIFIED_POPULATION';checkpoint(D/'result.json',data)
    except Exception as exc:
        data.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'result.json',data);raise

def replay():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or d['protocol_sha256']!=cert.hashed(D/'protocol.json') or len(d['shards'])!=20:raise ArithmeticError('all20 higher slices required')
    for row,actual in zip(p['rows'],d['shards']):
        if checked(row,True)!=actual:raise ArithmeticError('full higher slice replay differs')
    out=D/'replay.json';result={'status':'PASS','protocol_sha256':cert.hashed(D/'protocol.json'),'result_sha256':cert.hashed(D/'result.json'),'retained_scores':81920,'selection_primes':3510,'primitive_population':p['total_primitive_population'],'disjoint_previous_outer_slices':True}
    if out.exists():
        if cert.read(out)!=result:raise ArithmeticError('saved higher scanner replay differs')
    else:checkpoint(out,result)
    print('ALL81920 HIGHER MW16 FULL-SCORE SURVIVORS REPLAY',flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','benchmark','benchmark-check','run','replay']);a=parser.parse_args();benchmark(True) if a.stage=='benchmark-check' else globals()[a.stage]()
