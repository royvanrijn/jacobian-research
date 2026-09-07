#!/usr/bin/env python3
"""Fresh outer slices scored through32749 before any survivor truncation."""
import argparse
from pathlib import Path
from hashlib import sha256
from concurrent.futures import ThreadPoolExecutor,as_completed
import benchmark_joint_11952_annulus as engine
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=engine.ROOT;CAS=engine.CAS;ART=engine.ART;D=ROOT/'artifacts/local/elliptic-curves/joint11952-fresh-annulus-v1'
SALT='fresh11952-complete3510-prime-annular-ranking-v1';KEEP=4096

def sources():return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),Path(engine.__file__),engine.OUT,engine.BINARY,engine.CPP,engine.TEST,engine.short.TABLE,engine.short.OUT,engine.extended.TABLE,engine.extended.OUT,engine.reader.BINARY,engine.reader.OUT,engine.prior.D/'protocol.json',engine.D/'ledger.json',CAS/'verify_periodic_nagao_scanner.py',CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py']}
def roster():
    rows=[]
    for old in cert.read(engine.prior.D/'protocol.json')['rows']:
        quarter=old['quarter'];parity=old['shard']%2;choices=[quarter*4096+2*k+parity for k in range(2048) if quarter*4096+2*k+parity!=old['shard']]
        h=sha256(f"{SALT}|{old['sign']}|{quarter}".encode()).hexdigest();shard=choices[int(h,16)%len(choices)]
        count=engine.prior.annulus.counts.population(524288,524288,shard,16384)-engine.prior.annulus.counts.population(131072,131072,shard,16384)
        rows.append({'id':old['id'],'family':'11952','sign':old['sign'],'quarter':quarter,'N':524288,'M':524288,'inner':131072,'shard':shard,'shards':16384,'excluded_previous_shard':old['shard'],'sha256_choice':h,'primitive_population':count})
    return rows

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve fresh joint parameter protocol')
    proof=cert.read(engine.OUT)
    if proof['status']!='PASS' or not proof['cost_gate_passed'] or cert.read(engine.D/'ledger.json')['status']!='PASS' or cert.hashed(engine.BINARY)!=proof['binary_sha256']:raise ArithmeticError('complete joint-score output/cost/replay gate required')
    rows=roster()
    if len(rows)!=8 or sum(r['shard']%2 for r in rows)!=4 or any(r['shard']==r['excluded_previous_shard'] for r in rows):raise ArithmeticError('fresh balanced signed slices required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.joint11952-fresh-annulus.v1','sources':sources(),'rows':rows,'keep_per_slice':KEEP,'prime_count':3510,'maximum_workers':2,'seconds_per_call':120,'first_slice_cost_gate_seconds':45,'run_wall_seconds':600,'rss_bytes':1610612736,'total_primitive_population':sum(r['primitive_population'] for r in rows),'gate':'The full longer-score benchmark enumerates12690811 new-annulus addresses in20.755 seconds and passes independent whole-frame/top7 and returned-score replays. Of its top512 combined-score parameters,366 were omitted by the prior4096 short-score retention. This exact loss and the passed45-second cost gate justify computing the stronger score before truncation on a fresh prospective population. No rank improvement is inferred from score retention.','population':'Eight new signed denominator slices modulo16384 in131072<max(abs(n),d)<=524288, one in each residue quarter per sign and balanced parity. SHA256 chooses among2047 same-parity alternatives excluding the earlier slice. Every address is disjoint from the previous complete131072 square and from the earlier eight-slice outer trial.','selection':'Every primitive address receives all3510 quantized scores at5..4093 and4099..32749 before heap admission. Retain4096 per slice by combined units, good count, denominator and absolute numerator. No short-score heap precedes this. Quantization remains exactly the two pinned cache expressions. Scores through65521 and points are not computed here.','validation':'Reuse the proven joint binary. Before real slices, check entire signed large-modulus frames and top7 order for all eight new residues against the independently compiled retained-list scorer on both caches. Require exact population counts and recheck all32768 retained combined scores, retaining their separate short and extension components. First real slice has a45-second cost gate and is reused once. No retries or replacement slices.','future_scope':'After all scan/score replays, a separately frozen protocol may choose4096 equation-distinct candidates under combined selection score, apply scalar traces through65521 and wholly disjoint validation on64 finalists. Previous and catalogue equations may only be exclusions. Generic17-only point exposure requires a later fixed protocol with all maps before points, at most49 charts per fibre at125000/10seconds. No automatic point search or adaptive wave follows this scanner.','boundaries':'Complete combined scoring within eight declared slices, not complete annulus coverage, an unrestricted global optimum, complete65521 scoring, rank density, point absence or universal novelty.'})
    print('FROZEN8 FRESH JOINT11952 SLICES',sum(r['primitive_population'] for r in rows),'FULL3510-PRIME ADDRESSES',flush=True)
def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['rows']!=roster():raise ArithmeticError('frozen fresh joint sources/population differ')
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
    raw=execute(name,[str(engine.BINARY),str(engine.short.TABLE),str(engine.extended.TABLE),*map(str,values)],check);lines=raw['stdout'].splitlines()
    count=engine.prior.annulus.counts.population(c['N'],c['M'],c['shard'],c['shards'])-engine.prior.annulus.counts.population(min(c['N'],c['inner']),min(c['M'],c['inner']),c['shard'],c['shards']);kept=min(count,keep)
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

def components(rows,name,check=False):
    path=D/(name+'-candidates.txt');text='R17-CANDIDATES-V1 '+str(len(rows))+'\n'+''.join(f"{r['numerator']} {r['denominator']}\n" for r in rows)
    if check:
        if path.read_text()!=text:raise ArithmeticError('component reader roster changed')
    else:
        if path.exists():raise FileExistsError('preserve component reader input')
        path.parent.mkdir(parents=True,exist_ok=True);path.write_text(text)
    sums=[]
    for label,module,count in [('short',engine.short,562),('extended',engine.extended,2948)]:
        raw=execute(name+'-'+label,[str(engine.reader.BINARY),str(module.TABLE),str(path),'524288'],check);lines=raw['stdout'].splitlines();values=[]
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
        c={**row,'N':31,'M':32768,'inner':16384};name='benchmark/'+row['id']
        bare=[{'numerator':c['sign']*n,'denominator':d,'parameter':str(cert.F(c['sign']*n,d))} for d in range(c['shard']+1,c['M']+1,c['shards']) for n in range(1,32) if engine.gcd(n,d)==1 and max(n,d)>16384]
        expected=components(bare,name+'-reference',check);expected.sort(key=lambda r:(-r['combined_selection_units'],-r['combined_good'],r['denominator'],abs(r['numerator'])))
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
    print('ALL8 FRESH FULL-SCORE FRAMES PASS',flush=True)

def checked(row,check=False):
    rows,raw,count=scan_case(row,KEEP,row['id']+'/scan',check);verified=components(rows,row['id']+'/components',check)
    if count!=row['primitive_population'] or len(rows)!=KEEP:raise ArithmeticError('full retained slice/count differs')
    return {**row,'rows':verified,'raw_sha256':cert.hashed(D/row['id']/'scan.json'),'wall_seconds':raw['supervision']['wall_seconds']}

def run():
    p=protocol();benchmark(True)
    if (D/'result.json').exists():raise FileExistsError('preserve fresh joint population')
    data={'status':'RUNNING','protocol_sha256':cert.hashed(D/'protocol.json'),'shards':[]};checkpoint(D/'result.json',data);done={}
    try:
        done[0]=checked(p['rows'][0]);data['shards']=[done[0]];checkpoint(D/'result.json',data)
        if done[0]['wall_seconds']>p['first_slice_cost_gate_seconds']:raise ArithmeticError('fresh full-score first-slice runtime gate failed')
        with ThreadPoolExecutor(max_workers=2) as pool:
            pending={pool.submit(checked,row):i for i,row in enumerate(p['rows']) if i}
            for future in as_completed(pending):
                i=pending[future];done[i]=future.result();data['shards']=[done[j] for j in sorted(done)];checkpoint(D/'result.json',data);print('FRESH JOINT11952',done[i]['id'],done[i]['wall_seconds'],flush=True)
        data['status']='COMPLETE_FIXED_STRATIFIED_POPULATION';checkpoint(D/'result.json',data)
    except Exception as exc:
        data.update(status='FAILED_OR_CENSORED',reason=str(exc));checkpoint(D/'result.json',data);raise

def replay():
    p=protocol();d=cert.read(D/'result.json')
    if d['status']!='COMPLETE_FIXED_STRATIFIED_POPULATION' or d['protocol_sha256']!=cert.hashed(D/'protocol.json') or len(d['shards'])!=8:raise ArithmeticError('all eight completed full-score slices required')
    for row,actual in zip(p['rows'],d['shards']):
        if checked(row,True)!=actual:raise ArithmeticError('fresh full-score slice replay differs')
    out=D/'replay.json';result={'status':'PASS','protocol_sha256':cert.hashed(D/'protocol.json'),'result_sha256':cert.hashed(D/'result.json'),'retained_scores':32768,'selection_primes':3510,'primitive_population':p['total_primitive_population'],'disjoint_previous_outer_slices':True}
    if out.exists():
        if cert.read(out)!=result:raise ArithmeticError('saved joint replay differs')
    else:checkpoint(out,result)
    print('ALL32768 FRESH FULL-SCORE SURVIVORS REPLAY',flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','benchmark','benchmark-check','run','replay']);a=p.parse_args();benchmark(True) if a.stage=='benchmark-check' else globals()[a.stage]()
