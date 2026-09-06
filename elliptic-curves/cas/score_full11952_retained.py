#!/usr/bin/env python3
"""Verify a million retained short scores and select64 by unchanged extended S1."""
import argparse,heapq,json,struct
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import certify_compact_r17_candidates as cert
import scan_full11952_h131072 as scan
import extend_outer131072_r17 as extension
import benchmark_r17_extended_prime_traces as scalar
import extend_retained_r17_prime_scores as scoring
from fixed_cubic_geometry import poly,power,scale,add
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
ROOT=scan.ROOT;CAS=scan.CAS;ART=scan.ART;D=ROOT/'artifacts/local/elliptic-curves/full11952-h131072-retained-v1'
INPUT=D/'candidates.txt';META=D/'candidate-metadata.bin';EXT=D/'extended-scores.bin';OUT=ART/'full11952_h131072_selection_v1.json';OLD=ART/'outer48_r17_results_v1.json'
RECORD=struct.Struct('<iiqI');EXTRA=struct.Struct('<qI')

def sources():
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in [Path(__file__).resolve(),Path(scan.__file__).resolve(),scan.D/'protocol.json',scan.D/'ledger.json',scan.D/'scan.supervisor.json',OLD,Path(scalar.__file__).resolve(),Path(scoring.__file__).resolve(),Path(extension.__file__).resolve(),Path(cert.__file__).resolve(),extension.spec.ATLAS,CAS/'fixed_cubic_geometry.py']}

def previous():
    old=cert.read(OLD);rows=old['previous_equations']+[{'address':OLD.name+':'+r['id'],'curve':r['curve']} for r in old['curves']]
    if len(rows)!=789:raise ArithmeticError('fixed789 previous equations required')
    return rows

def prepare():
    if (D/'protocol.json').exists():raise FileExistsError('preserve full retained-score protocol')
    p=scan.protocol();ledger=cert.read(scan.D/'ledger.json');s=cert.read(scan.D/'scan.supervisor.json')
    if ledger['status']!='COMPLETE_FIXED_SCAN_PENDING_SCORE_REPLAY' or s['outcome']!='completed' or s['returncode']!=0 or len(ledger['rows'])!=2048 or any(r['status']!='COMPLETE_PENDING_CACHED_SCORE_REPLAY' for r in ledger['rows']):raise ArithmeticError('all2048 terminal scanner invocations required')
    previous();f=next(r for r in cert.read(extension.spec.ATLAS)['families'] if r['family']=='11952');a=poly(f['A_coefficients_low_to_high']);b=poly(f['B_coefficients_low_to_high']);n=scale(power(a,3),6912);d=add(scale(power(a,3),4),scale(power(b,2),27))
    size=max(len(n),len(d));n=list(n)+[0]*(size-len(n));d=list(d)+[0]*(size-len(d))
    if not any(d) or size>25 or all(n[i]*d[j]==n[j]*d[i] for i in range(size) for j in range(size)):raise ArithmeticError('nonconstant degree-at-most24 j-map required')
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.full11952-h131072-retained.v1','sources':sources(),'recorded_root':str(ROOT),'family':'11952','rows':1048576,'cached_score_bands':p['cached_score_bands'],'short_primes':562,'extension_primes':2948,'bank_seconds':600,'cached_seconds_per_band':300,'selection_seconds':600,'selected_curves':64,'retained_sorted_prefix':32768,'previous_equations':789,'j_degree_upper_bound':24,'validation_seconds_per_curve':20,'validation_workers':2,'validation_wall_seconds':900,'direct_validation_indices':[0,1,len(scalar.PRIMES)//2,len(scalar.PRIMES)-1],'gp_sha256':cert.hashed(scalar.GP),'rss_bytes':2147483648,'ordering':'Combined exact quantized S1 through32749 descending, combined good-prime count descending, denominator ascending, signed numerator ascending. No validation prime or public rank enters ties.','deduplication':'Exactly64 rational-isomorphism-distinct nonsingular equations after excluding all789 previously measured equations and earlier selected models, using equations only. At most32768 leading score rows need be retained: a nonconstant j-map of degree at most24 has at most24 preimages of each prior or selected j, plus at most24 singular finite parameters. The conservative prefix exceeds24*(789+63+1)+64. No result-dependent prefix refill occurs.','gate':'The full11952 short population and all2048 raw calls are terminal. Cached score lookup matches the scalar reference on967 old inputs and every row of the full-size repeated fixture. Keep every one of the1048576 short survivors until extended S1 is available, then apply exact equation exclusions. This tests the larger retention population motivated by the known-control ordering contrast.','scope':'Verify every raw slice frame and primitive address, write a hash-bound million-row bank, and require all562-prime short sums to match the original scanner. Apply the unchanged cached selection extension, freeze64 new equations, then check fresh scalar traces on all64, including exact agreement of cached selection sums and a disjoint validation band. No global short-prefix truncation, public-target selection, point search, inferred rank or automatic adaptive wave. A separate point protocol must freeze all64 maps before any points.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources() or p['gp_sha256']!=cert.hashed(scalar.GP):raise ArithmeticError('frozen retained-score inputs changed')
    return p

def bank(check=False):
    p=protocol();parent=scan.protocol();ledger=cert.read(scan.D/'ledger.json');manifest=D/'bank.json';header=b'R17BANK1'+struct.pack('<II',p['rows'],RECORD.size)
    if not check and any(q.exists() for q in (INPUT,META,manifest)):raise FileExistsError('preserve retained parameter bank')
    with META.open('rb' if check else 'xb') as meta,INPUT.open('r' if check else 'x') as text:
        def binary(v):
            if check:
                if meta.read(len(v))!=v:raise ArithmeticError('candidate metadata differs')
            else:meta.write(v)
        def line(v):
            if check:
                if text.readline()!=v:raise ArithmeticError('candidate text differs')
            else:text.write(v)
        binary(header);line(f"R17-CANDIDATES-V1 {p['rows']}\n");count=0;bindings=[]
        for row,kept in zip(parent['rows'],ledger['rows']):
            actual,rows=scan.checked(row,parent,False)
            if actual!=kept:raise ArithmeticError('terminal slice ledger differs')
            bindings.append({'id':row['id'],'raw_sha256':kept['raw_sha256']})
            for r in rows:
                binary(RECORD.pack(r['numerator'],r['denominator'],r['score_units'],r['good_primes']));line(f"{r['numerator']} {r['denominator']}\n");count+=1
        if count!=p['rows']:raise ArithmeticError('fixed million-row bank incomplete')
        if check and (meta.read(1) or text.read(1)):raise ArithmeticError('trailing bank bytes')
    result={'status':'PASS','protocol_hash':digest(p),'rows':count,'metadata_sha256':cert.hashed(META),'input_sha256':cert.hashed(INPUT),'slice_bindings':bindings}
    if check:
        if cert.read(manifest)!=result:raise ArithmeticError('retained bank report differs')
    else:checkpoint(manifest,result)
    print('FULL11952 EXACT RETAINED BANK',count,'ROWS',flush=True)

def iter_bank(p):
    with META.open('rb') as f:
        if f.read(16)!=b'R17BANK1'+struct.pack('<II',p['rows'],RECORD.size):raise ArithmeticError('bank header differs')
        for i in range(p['rows']):
            value=f.read(RECORD.size)
            if len(value)!=RECORD.size:raise ArithmeticError('short metadata bank')
            yield i,RECORD.unpack(value)
        if f.read(1):raise ArithmeticError('trailing metadata bank')

def cached(label,check=False):
    p=protocol();bank=cert.read(D/'bank.json');band=p['cached_score_bands'][label];folder=D/label;path=folder/'raw.json';output=folder/'result.json';cache=ROOT/band['cache'];binary=ROOT/band['binary']
    if bank['status']!='PASS' or bank['protocol_hash']!=digest(p) or cert.hashed(INPUT)!=bank['input_sha256'] or cert.hashed(META)!=bank['metadata_sha256'] or cert.hashed(cache)!=band['cache_sha256'] or cert.hashed(binary)!=band['binary_sha256']:raise ArithmeticError('cached-score input binding differs')
    cmd=[str(binary),str(cache),str(INPUT)]
    if not check:
        if path.exists() or output.exists():raise FileExistsError('preserve one retained cached-score call')
        c=capture(cmd,limits=Limits(p['cached_seconds_per_band'],p['rss_bytes']),log_path=folder/'score.log',separate_stderr=True,check=False)
        checkpoint(path,{'command':cmd,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
    raw=cert.read(path);expected=[str(Path(p['recorded_root'])/band['binary']),str(Path(p['recorded_root'])/band['cache']),str(Path(p['recorded_root'])/INPUT.relative_to(ROOT))];primes=p['short_primes'] if label=='short' else p['extension_primes'];lines=raw['stdout'].splitlines()
    if raw['command']!=expected or raw['supervision']['command']!=expected or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0 or len(lines)!=p['rows']+1 or lines[-1]!=f"S {p['rows']} {primes}":raise ArithmeticError('cached retained invocation failed/censored or differs')
    f=None
    try:
        if label=='extended':
            if not check and EXT.exists():raise FileExistsError('preserve extended-score vector')
            f=EXT.open('rb' if check else 'xb');header=b'R17EXT01'+struct.pack('<II',p['rows'],primes)
            if check:
                if f.read(16)!=header:raise ArithmeticError('extended-score header differs')
            else:f.write(header)
        for (i,(_,_,units,good)),line in zip(iter_bank(p),lines[:-1]):
            v=line.split()
            if len(v)!=4 or v[0]!='R' or int(v[1])!=i:raise ArithmeticError('retained score index differs')
            score,count=int(v[2]),int(v[3])
            if not 0<=count<=primes or abs(score)>primes*10**13 or (not count and score):raise ArithmeticError('retained score range differs')
            if label=='short':
                if (score,count)!=(units,good):raise ArithmeticError('cached short score differs from original scanner')
            else:
                data=EXTRA.pack(score,count)
                if check:
                    if f.read(EXTRA.size)!=data:raise ArithmeticError('extended-score bytes differ')
                else:f.write(data)
        if check and f and f.read(1):raise ArithmeticError('trailing extended-score bytes')
    finally:
        if f:f.close()
    result={'status':'PASS','protocol_hash':digest(p),'band':label,'rows':p['rows'],'primes':primes,'bank_sha256':cert.hashed(D/'bank.json'),'raw_sha256':cert.hashed(path),'wall_seconds':raw['supervision']['wall_seconds'],'extended_vector_sha256':cert.hashed(EXT) if label=='extended' else None}
    if check:
        if cert.read(output)!=result:raise ArithmeticError('retained cached-score report differs')
    else:checkpoint(output,result)
    print('FULL11952 CACHED',label,p['rows'],'ROWS',raw['supervision']['wall_seconds'],'SECONDS',flush=True)

def ranked_prefix(p):
    bank=cert.read(D/'bank.json');a=cert.read(D/'short/result.json');b=cert.read(D/'extended/result.json')
    if any(r['status']!='PASS' or r['protocol_hash']!=digest(p) for r in (bank,a,b)) or cert.hashed(META)!=bank['metadata_sha256'] or cert.hashed(EXT)!=b['extended_vector_sha256']:raise ArithmeticError('complete cached-score proofs required')
    def rows():
        with EXT.open('rb') as f:
            if f.read(16)!=b'R17EXT01'+struct.pack('<II',p['rows'],p['extension_primes']):raise ArithmeticError('extended vector header differs')
            for i,(n,d,s,g) in iter_bank(p):
                raw=f.read(EXTRA.size)
                if len(raw)!=EXTRA.size:raise ArithmeticError('short extended score vector')
                e,h=EXTRA.unpack(raw);yield (s+e,g+h,-d,-n,i,s,g,e,h)
            if f.read(1):raise ArithmeticError('trailing extended score vector')
    return heapq.nlargest(p['retained_sorted_prefix'],rows())

def selection():
    p=protocol();prefix=ranked_prefix(p);f=next(r for r in cert.read(extension.spec.ATLAS)['families'] if r['family']==p['family']);seen={};rows=[];skips=[]
    def jvalue(model):
        inv=cert.weierstrass_invariants(model);return inv['c4']**3/inv['discriminant']
    for row in previous():
        model=tuple(map(cert.F,row['curve']));seen.setdefault(jvalue(model),[]).append((model,row['address']))
    for u,g,negative_d,negative_n,i,s,sg,e,eg in prefix:
        n,d=-negative_n,-negative_d;t=str(cert.F(n,d))
        try:model=tuple(map(cert.F,extension.model_at(f,t)))
        except ArithmeticError:
            value=cert.F(t);aa=extension.spec.polynomial(f['A_coefficients_low_to_high'],value)*d**8;bb=extension.spec.polynomial(f['B_coefficients_low_to_high'],value)*d**12
            if aa.denominator!=1 or bb.denominator!=1 or 4*aa**3+27*bb**2!=0:raise
            skips.append({'index':i,'parameter':t,'reason':'SINGULAR_DISPLAYED_FIBRE'});continue
        j=jvalue(model);matches=[name for q,name in seen.get(j,[]) if cert.isomorphic(model,q)]
        if matches:
            skips.append({'index':i,'parameter':t,'reason':'PREVIOUS_OR_SELECTED_Q_ISOMORPHISM','matches':matches});continue
        row={'id':f"11952-{i:07}",'family':'11952','retained_index':i,'parameter':t,'numerator':n,'denominator':d,'score_units':s,'good_primes':sg,'extension_selection_units':e,'extension_good':eg,'combined_selection_units':u,'combined_good':g,'model':list(map(str,model))};rows.append(row);seen.setdefault(j,[]).append((model,row['id']))
        if len(rows)==p['selected_curves']:break
    if len(rows)!=64:raise ArithmeticError('fixed64 eligible roster not obtained; no refill')
    return {'schema':'elliptic-curves.full11952-h131072-selection.v1','status':'PASS_FROZEN64_SELECTION','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'bank_sha256':cert.hashed(D/'bank.json'),'short_replay_sha256':cert.hashed(D/'short/result.json'),'extended_replay_sha256':cert.hashed(D/'extended/result.json'),'sorted_prefix_sha256':digest(prefix),'sorted_prefix_rows':len(prefix),'selected':rows,'skipped_before_completion':skips,'previous_equations':789,'claim_boundary':p['scope']}

def select(check=False):
    result=selection()
    if check:
        if cert.read(OUT)!=result:raise ArithmeticError('fixed64 selection replay differs')
    else:
        if OUT.exists():raise FileExistsError('preserve fixed64 selection')
        checkpoint(OUT,result)
    print('FROZEN64 FULL11952 SELECTION',len(result['skipped_before_completion']),'EXACT EQUATION SKIPS',flush=True)

def validate_selected(check=False):
    p=protocol();selected=cert.read(OUT);path=D/'selected-validation/result.json'
    if selected['status']!='PASS_FROZEN64_SELECTION' or selected['protocol_sha256']!=cert.hashed(D/'protocol.json'):raise ArithmeticError('fixed64 roster required before validation')
    if not check and path.exists():raise FileExistsError('preserve fixed64 scalar validation')
    def one(row):
        rawpath=D/'selected-validation'/row['id']/'raw.json';code=scalar.program(row['model']);cmd=[str(scalar.GP),'-q','-s','256000000']
        if not check:
            if rawpath.exists():raise FileExistsError('preserve scalar finalist call')
            c=capture(cmd,input_text=code,limits=Limits(p['validation_seconds_per_curve'],p['rss_bytes']),log_path=rawpath.with_suffix('.log'),separate_stderr=True,check=False)
            checkpoint(rawpath,{'program':code,'stdout':c.stdout,'stderr':c.stderr,'supervision':c.supervision})
        raw=cert.read(rawpath)
        if raw['program']!=code or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0:raise ArithmeticError('selected scalar trace failed or censored')
        values,ms=scalar.parse(raw['stdout'],row['model']);scores=scoring.sums(values);checks=[]
        if scores['extension_selection_units']!=row['extension_selection_units'] or scores['extension_good']!=row['extension_good']:raise ArithmeticError('new finalist scalar selection score differs from cache')
        for i in p['direct_validation_indices']:
            q,value=values[i]
            if scalar.direct(row['model'],q)!=value:raise ArithmeticError('selected direct character sum differs')
            checks.append([q,value])
        return {'id':row['id'],'status':'PASS','scores':scores,'direct_checks':checks,'raw_sha256':cert.hashed(rawpath),'cpu_ms':ms}
    rows=[]
    with ThreadPoolExecutor(max_workers=p['validation_workers']) as pool:
        for start in range(0,64,8):
            rows+=list(pool.map(one,selected['selected'][start:start+8]))
            if not check:checkpoint(path,{'status':'RUNNING','selection_sha256':cert.hashed(OUT),'rows':rows})
    result={'status':'PASS','selection_sha256':cert.hashed(OUT),'protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows,'direct_character_sum_checks':256,'claim_boundary':'All64 frozen finalists have fresh scalar traces agreeing exactly with the cached selection sums, plus a disjoint validation band and256 independent character sums. No validation value changes selection and no point search occurs.'}
    if check:
        if cert.read(path)!=result:raise ArithmeticError('selected64 scalar validation replay differs')
    else:checkpoint(path,result)
    print('ALL64 FRESH SCALAR TRACE AND DISJOINT VALIDATION CHECKS PASS',flush=True)
if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('stage',choices=['prepare','bank','bank-check','short','short-check','extended','extended-check','select','selection-check','validate-selected','validation-check']);v=a.parse_args()
    if v.stage=='prepare':prepare()
    elif v.stage.startswith('bank'):bank(v.stage.endswith('-check'))
    elif v.stage.startswith('short'):cached('short',v.stage.endswith('-check'))
    elif v.stage.startswith('extended'):cached('extended',v.stage.endswith('-check'))
    elif v.stage in ('select','selection-check'):select(v.stage=='selection-check')
    else:validate_selected(v.stage=='validation-check')
