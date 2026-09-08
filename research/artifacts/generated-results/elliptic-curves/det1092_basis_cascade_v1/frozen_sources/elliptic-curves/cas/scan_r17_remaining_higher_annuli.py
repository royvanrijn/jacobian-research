#!/usr/bin/env python3
"""Untouched higher R17 slices scored through32749 before any truncation."""
import argparse
from pathlib import Path
from hashlib import sha256
from concurrent.futures import ThreadPoolExecutor,as_completed
import benchmark_joint_11952_annulus as engine
import verify_periodic_nagao_scanner as counts
import scan_outer131072_r17_v2 as prior
import scan_skew_r17_boxes as skew
import encode_r17_remaining_joint_caches as binary
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture,Limits
ROOT=prior.ROOT;CAS=prior.CAS;ART=prior.ART
D=ROOT/'artifacts/local/elliptic-curves/r17-remaining-higher-annuli-v1'
SALT='user-fill-five-r17-gaps-new-higher-16-slices-v1';KEEP=4096

def sources():
    paths=[Path(__file__).resolve(),Path(prior.__file__),prior.D/'protocol.json',
           Path(skew.__file__),skew.D/'protocol.json',Path(counts.__file__),
           Path(binary.__file__),binary.OUT,binary.D/'score-replay.json',binary.D/'controller/ledger.json',
           engine.BINARY,engine.CPP,engine.TEST,engine.OUT,engine.D/'ledger.json',
           engine.reader.BINARY,engine.reader.OUT,ART/'r17_extended_score_prime_minimality_v1.json',
           ART/'r17_extended_score_prime_minimality_sage_replay_v1.json',CAS/'research_runtime/store.py',
           CAS/'research_runtime/supervisor.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def roster():
    rows=[]
    old_rows=[r for r in cert.read(prior.D/'protocol.json')['rows'] if r['family']!='11952']
    skew_rows=cert.read(skew.D/'protocol.json')['rows']
    if len(old_rows)!=10 or len({r['family'] for r in old_rows})!=5:
        raise ArithmeticError('all five remaining signed families required')
    for band,inner,outer,modulus in [(2,32768,131072,1024),(3,131072,524288,16384)]:
        for r in old_rows:
            exclusions=[{'modulus':1024,'shard':r['shard'],'origin':'outer131072'}]
            exclusions += [{'modulus':q['shards'],'shard':q['shard'],'origin':'skew'}
                           for q in skew_rows if (q['family'],q['sign'])==(r['family'],r['sign'])]
            choices=[]
            for shard in range(r['shard']%2,modulus,2):
                if any(modulus%e['modulus'] or shard%e['modulus']==e['shard'] for e in exclusions):continue
                h=sha256(f"{SALT}|{band}|{r['family']}|{r['sign']}|{shard}".encode()).hexdigest()
                choices.append((h,shard))
            if len(choices)<16:raise ArithmeticError('insufficient untouched residues')
            for slot,(h,shard) in enumerate(sorted(choices)[:16]):
                count=counts.population(outer,outer,shard,modulus)-counts.population(inner,inner,shard,modulus)
                parent_id=f"b{band}-{r['family']}-"+('neg' if r['sign']<0 else 'pos')
                rows.append({'id':f"r17-{parent_id}-{slot:02}",'parent_slice_id':parent_id,'slot':slot,
                    'family':r['family'],'band':band,'sign':r['sign'],'N':outer,'M':outer,'inner':inner,
                    'shard':shard,'shards':modulus,'excluded_previous_residue_classes':exclusions,
                    'sha256_choice':h,'primitive_population':count})
    return rows

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve higher full-score scanner')
    proof=cert.read(binary.OUT);gate=cert.read(binary.D/'controller/ledger.json')
    compiled=cert.read(engine.OUT)
    if any(cert.read(ART/name)['status']!='PASS' for name in ['r17_extended_score_prime_minimality_v1.json','r17_extended_score_prime_minimality_sage_replay_v1.json']):raise ArithmeticError('complete independently replayed R17 prime range required')
    if proof['status']!='PASS' or proof['candidates_checked']!=4831 or gate['status']!='PASS':raise ArithmeticError('all five cache proofs and4831 scalar agreements required')
    if compiled['status']!='PASS' or cert.hashed(engine.BINARY)!=compiled['binary_sha256'] or cert.read(engine.D/'ledger.json')['status']!='PASS':raise ArithmeticError('exact full-score binary gate required')
    rows=roster()
    if len(rows)!=320 or sum(r['shard']%2 for r in rows)!=160 or any(any(r['shard']%e['modulus']==e['shard'] for e in r['excluded_previous_residue_classes']) for r in rows):raise ArithmeticError('balanced fresh320 required')
    bindings={}
    for f in {r['family'] for r in rows}:
        for label in ('short','extended'):
            path=binary.cache(f,label);meta=next(r for r in proof['encoding']['rows'] if r['family']==f and r['label']==label)
            if cert.hashed(path)!=meta['sha256']:raise ArithmeticError('encoded cache changed')
            bindings[str(path.relative_to(ROOT))]=meta['sha256']
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.broad-r17_remaining-higher-annuli.v1','sources':sources(),
        'rows':rows,'cache_bindings':bindings,'keep_per_slice':KEEP,'prime_count':3510,
        'maximum_workers':4,'seconds_per_call':120,'first_slice_per_band_cost_gate_seconds':45,
        'run_wall_seconds':7200,'replay_wall_seconds':1800,'rss_bytes':1610612736,
        'total_primitive_population':sum(r['primitive_population'] for r in rows),
        'gate':'The user explicitly requests filling the five remaining R17 extended-score coverage gaps and broadening early search in new high-height territory. All five complete extended projective caches and4831 independent saved scalar-score agreements must pass first. The established joint scanner and exact R17 score-prime minimality theorem support longer scoring before any truncation. Give each family sixteen fresh denominator slices per sign and band, then the same10240 scalar and sixty finalist limits used in the recent MW16 higher-band trials. No generic-rank-only family choice or rank-density inference.',
        'population':'Exactly320 fresh signed denominator slices on103b2,074d9,07ca9,08234,08f72 in32768<H<=131072 and131072<H<=524288. Sixteen SHA256-selected same-parity residues per family/sign/band, modulo1024 and16384. Every slice excludes the previous outer131072 residue and any prior skew-rectangle residue by its exact modulus. Earlier full squares through32768 lie below the inner cut.11952 is excluded because its complete cache and separate higher campaigns already exist. All new slices are disjoint within each family; mathematical equation novelty is checked separately.',
        'selection':'Every primitive parameter receives all3510 quantized selection primes before heap admission. Retain4096 per signed slice by units,good count,denominator,absolute numerator. There is no short-score prefilter. Public record equations,parameters,points,ranks,j-invariants and jump labels remain outside selection and execution.',
        'validation':'Complete signed actual-modulus frames and top7 ordering for all320 cases against the independent retained-list reader on both family caches. Check every retained score component and exact primitive population. One first slice per band has a45-second cost gate and is reused once. No retry, adaptive expansion or replacement slices.',
        'future_scope':'Only after complete scan replay, a separate protocol may select1024 distinct prospective equations per band/family from the131072-address combined signed pool,10240 total, score with fresh scalar traces through65521 and independently check the cached first extension. Then select six per band/family,60 total, with disjoint65537..131071 validation. A later generic17-only point protocol may use at most49 charts per fibre at125000/10seconds and2940 boxes overall, all maps before points. No catalogue or prior-record labels enter this prospective selection; post-terminal exact novelty comparisons remain separate.',
        'boundaries':'Finite full-score selection inside320 declared slices, not complete annulus coverage, a rank predictor, a rank upper bound or universal novelty. No original point search under this scanner.'})
    print('FROZEN320 BROAD HIGHER R17 SLICES',sum(r['primitive_population'] for r in rows),'FULL3510-PRIME ADDRESSES',flush=True)

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['rows']!=roster():raise ArithmeticError('frozen higher R17 scanner changed')
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
    count=counts.population(c['N'],c['M'],c['shard'],c['shards'])-counts.population(min(c['N'],c['inner']),min(c['M'],c['inner']),c['shard'],c['shards']);kept=min(count,keep)
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
        raw=execute(name+'-'+label,[str(engine.reader.BINARY),str(binary.cache(family,label)),str(path),'524288'],check);lines=raw['stdout'].splitlines();values=[]
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
    print('ALL320 BROAD R17 FULL-SCORE FRAMES PASS',flush=True)

def load_shard(manifest):
    path=ROOT/manifest['retained_path']
    if path!=D/manifest['id']/'retained.json' or cert.hashed(path)!=manifest['retained_sha256']:
        raise ArithmeticError('immutable broad retained shard changed')
    result=cert.read(path)
    if result['id']!=manifest['id'] or len(result['rows'])!=KEEP:
        raise ArithmeticError('broad retained shard roster differs')
    return result

def checked(row,check=False):
    rows,raw,count=scan_case(row,KEEP,row['id']+'/scan',check)
    verified=components(rows,row['family'],row['id']+'/components',check)
    if count!=row['primitive_population'] or len(rows)!=KEEP:raise ArithmeticError('full broad slice/count differs')
    saved={**row,'rows':verified,'raw_sha256':cert.hashed(D/row['id']/'scan.json'),'wall_seconds':raw['supervision']['wall_seconds']}
    path=D/row['id']/'retained.json'
    if check:
        if cert.read(path)!=saved:raise ArithmeticError('broad retained point/score vectors differ')
    else:
        if path.exists():raise FileExistsError('preserve broad retained shard')
        checkpoint(path,saved)
    return {**row,'retained_path':str(path.relative_to(ROOT)),'retained_sha256':cert.hashed(path),
            'retained_count':len(verified),'raw_sha256':saved['raw_sha256'],'wall_seconds':saved['wall_seconds']}

def run():
    p=protocol();benchmark(True)
    if (D/'result.json').exists():raise FileExistsError('preserve higher joint population')
    data={'status':'RUNNING','protocol_sha256':cert.hashed(D/'protocol.json'),'shards':[]};checkpoint(D/'result.json',data);done={}
    try:
        for band in (2,3):
            i=next(i for i,r in enumerate(p['rows']) if r['band']==band)
            done[i]=checked(p['rows'][i]);data['shards']=[done[j] for j in sorted(done)];checkpoint(D/'result.json',data)
            if done[i]['wall_seconds']>p['first_slice_per_band_cost_gate_seconds']:raise ArithmeticError('first higher-band slice cost gate failed')
        with ThreadPoolExecutor(max_workers=4) as pool:
            pending={pool.submit(checked,row):i for i,row in enumerate(p['rows']) if i not in done}
            for future in as_completed(pending):
                i=pending[future];done[i]=future.result();data['shards']=[done[j] for j in sorted(done)];checkpoint(D/'result.json',data)
                print('BROAD HIGHER R17',done[i]['id'],done[i]['wall_seconds'],flush=True)
        data['status']='COMPLETE_FIXED_STRATIFIED_POPULATION';checkpoint(D/'result.json',data)
    except Exception as exc:
        data.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'result.json',data);raise

def replay():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or d['protocol_sha256']!=cert.hashed(D/'protocol.json') or len(d['shards'])!=320:raise ArithmeticError('all320 higher slices required')
    for row,actual in zip(p['rows'],d['shards']):
        if checked(row,True)!=actual:raise ArithmeticError('full higher slice replay differs')
    out=D/'replay.json';result={'status':'PASS','protocol_sha256':cert.hashed(D/'protocol.json'),'result_sha256':cert.hashed(D/'result.json'),'retained_scores':1310720,'selection_primes':3510,'primitive_population':p['total_primitive_population'],'disjoint_previous_outer_slices':True}
    if out.exists():
        if cert.read(out)!=result:raise ArithmeticError('saved higher scanner replay differs')
    else:checkpoint(out,result)
    print('ALL1310720 BROAD HIGHER R17 FULL-SCORE SURVIVORS REPLAY',flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','benchmark','benchmark-check','run','replay']);a=parser.parse_args();benchmark(True) if a.stage=='benchmark-check' else globals()[a.stage]()
