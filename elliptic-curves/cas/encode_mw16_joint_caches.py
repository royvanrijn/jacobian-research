#!/usr/bin/env python3
"""Exact binary caches and all40960 target-free saved scalar-score comparisons."""
import argparse, struct, sys, time
from array import array
from math import log
from pathlib import Path
import build_mw16_extended_projective_caches as parent
import scan_mw16_outer_bands as old
import benchmark_11952_annulus_cache_v3 as reader
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture, Limits

ROOT=parent.ROOT; CAS=parent.CAS; ART=parent.ART
D=ROOT/'artifacts/local/elliptic-curves/mw16-joint-binary-caches-v1'
OUT=ART/'mw16_joint_binary_caches_v1.json'
SAVED=ROOT/'artifacts/local/elliptic-curves/mw16-fresh-retention-scores-v1'

def sources():
    paths=[Path(__file__).resolve(),Path(parent.__file__),parent.OUT,parent.D/'replay.json',
           parent.D/'protocol.json',old.spec.ATLAS,SAVED/'protocol.json',SAVED/'result.json',
           SAVED/'controller/ledger.json',reader.BINARY,reader.OUT,
           CAS/'newfamily/score_retained_projective_cache_v3.cpp',
           CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py']
    return {str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}

def prepare():
    if (D/'protocol.json').exists(): raise FileExistsError('preserve cache encoding protocol')
    p=parent.protocol();proof=cert.read(parent.OUT);replay=cert.read(parent.D/'replay.json')
    if proof['status']!='PASS' or replay['status']!='PASS' or replay['certificate_sha256']!=cert.hashed(parent.OUT):
        raise ArithmeticError('all five full projective cache proofs required')
    if cert.read(SAVED/'controller/ledger.json')['status']!='PASS': raise ArithmeticError('complete saved target-free scalar trial required')
    saved=cert.read(SAVED/'protocol.json');scores=cert.read(SAVED/'result.json')
    if scores['status']!='COMPLETE_FROZEN40960' or len(saved['rows'])!=40960 or len(scores['rows'])!=40960:
        raise ArithmeticError('fixed40960 prior scalar cases required')
    extended={(r['family'],r['prime']):r for r in proof['rows']};families=[]
    for f in p['families']:
        tables={}
        for label,primes in [('short',old.PRIMES),('extended',p['primes'])]:
            rows=[]
            for prime in primes:
                path=(old.TABLES/f['family']/'trace-tables'/f'{prime}.json') if label=='short' else (parent.D/f['family']/str(prime)/'table.json')
                h=cert.hashed(path)
                if label=='extended' and h!=extended[(f['family'],prime)]['table_sha256']: raise ArithmeticError('proved table changed')
                rows.append({'prime':prime,'path':str(path.relative_to(ROOT)),'sha256':h})
            tables[label]=rows
        fixtures=[]
        for row,score in zip(saved['rows'],scores['rows']):
            if row['id']!=score['id']: raise ArithmeticError('saved scalar ordering differs')
            if row['family']==f['family']:
                fixtures.append({'id':row['id'],'numerator':row['numerator'],'denominator':row['denominator'],
                    'short':[row['score_units'],row['good_primes']],
                    'extended':[score['scores']['extension_selection_units'],score['scores']['extension_good']]})
        if len(fixtures)!=8192: raise ArithmeticError('balanced saved score fixtures required')
        families.append({'family':f['family'],'model_hash':f['model_hash'],'tables':tables,'fixtures':fixtures})
    checkpoint(D/'protocol.json',{'schema':'elliptic-curves.mw16-joint-binary-caches.v1','sources':sources(),
        'families':families,'encode_seconds':2400,'replay_seconds':2400,'score_seconds_per_call':60,
        'rss_bytes':1610612736,'maximum_height':262144,
        'quantization':'Short: round((2-ap)/(p+1-ap)*log(p)*10**12). Extended: round((2-ap)*log(p)/(p+1-ap)*10**12). Preserve Python operation order and zero units for singular raw-model residues.',
        'scope':'Encode all562 short and2948 extended prime tables for each of the five MW16 families in the existing little-endian R17XS001 format. Replay every byte against the frozen exact tables. Require every short/extended sum and good count on all40960 previously scored target-free candidates to agree with their saved canonical/scalar values using the independent compiled retained-list reader. No public-record input, new parameter scan, point exposure, selector or rank claim. The explicit reader bound262144 prepares later higher-parameter inputs.'})

def protocol():
    p=cert.read(D/'protocol.json')
    if p['sources']!=sources(): raise ArithmeticError('fixed encoding sources changed')
    return p

def cache(family,label): return D/family/(label+'.bin')

def frame(f,row,label):
    path=ROOT/row['path'];prime=row['prime']
    if cert.hashed(path)!=row['sha256']: raise ArithmeticError('frozen prime table changed')
    data=cert.read(path)
    if data['input']!={'family':f['family'],'model_hash':f['model_hash'],'prime':prime} or len(data['traces'])!=prime+1 or len(data['good'])!=prime+1:
        raise ArithmeticError('canonical projective table shape/model differs')
    if label=='short': values=[round((2-t)/(prime+1-t)*log(prime)*10**12) if good else 0 for t,good in zip(data['traces'],data['good'])]
    else: values=[round((2-t)*log(prime)/(prime+1-t)*10**12) if good else 0 for t,good in zip(data['traces'],data['good'])]
    units=array('q',values)
    if units.itemsize!=8: raise ArithmeticError('64-bit integers required')
    if sys.byteorder!='little': units.byteswap()
    return struct.pack('<II',prime,prime+1)+units.tobytes()+bytes(data['good'])

def encode(check=False):
    p=protocol();start=time.monotonic();rows=[]
    if not check and (D/'encoding.json').exists(): raise FileExistsError('preserve encoding attempt')
    for f in p['families']:
        for label,table_rows in f['tables'].items():
            path=cache(f['family'],label);path.parent.mkdir(parents=True,exist_ok=True)
            with path.open('rb' if check else 'xb') as out:
                def consume(value):
                    if check:
                        if out.read(len(value))!=value: raise ArithmeticError('binary encoding differs')
                    else: out.write(value)
                consume(b'R17XS001'+struct.pack('<I',len(table_rows)))
                for row in table_rows:
                    if time.monotonic()-start>p['replay_seconds' if check else 'encode_seconds']: raise TimeoutError('finite encoding deadline')
                    consume(frame(f,row,label))
                consume(b'ENDXSC01')
                if check and out.read(1): raise ArithmeticError('trailing cache data')
            rows.append({'family':f['family'],'label':label,'primes':len(table_rows),'sha256':cert.hashed(path),'bytes':path.stat().st_size})
            if not check: checkpoint(D/'encoding-progress.json',{'status':'RUNNING','rows':rows})
            print('MW16 BINARY CACHE',f['family'],label,'BYTE REPLAY' if check else 'ENCODED',flush=True)
    result={'status':'PASS','protocol_sha256':cert.hashed(D/'protocol.json'),'rows':rows}
    if check:
        if cert.read(D/'encoding.json')!=result: raise ArithmeticError('cache metadata differs')
        checkpoint(D/'encoding-replay.json',result)
    else: checkpoint(D/'encoding.json',result)

def score(check=False):
    p=protocol();encoded=cert.read(D/'encoding.json');replay=cert.read(D/'encoding-replay.json')
    if encoded!=replay or encoded['status']!='PASS': raise ArithmeticError('all binary bytes must replay first')
    if not check and OUT.exists(): raise FileExistsError('preserve complete score comparison')
    records=[]
    for f in p['families']:
        folder=D/f['family'];path=folder/'candidates.txt'
        value='R17-CANDIDATES-V1 '+str(len(f['fixtures']))+'\n'+''.join(f"{r['numerator']} {r['denominator']}\n" for r in f['fixtures'])
        if check:
            if path.read_text()!=value: raise ArithmeticError('fixed saved score roster changed')
        else:
            with path.open('x') as stream: stream.write(value)
        for label in ('short','extended'):
            table_path=cache(f['family'],label)
            meta=next(r for r in encoded['rows'] if r['family']==f['family'] and r['label']==label)
            if cert.hashed(table_path)!=meta['sha256']: raise ArithmeticError('encoded cache changed')
            cmd=[str(reader.BINARY),str(table_path),str(path),str(p['maximum_height'])];rawpath=folder/(label+'-raw.json')
            if not check:
                if rawpath.exists(): raise FileExistsError('preserve score call')
                result=capture(cmd,limits=Limits(p['score_seconds_per_call'],p['rss_bytes']),log_path=folder/(label+'.log'),separate_stderr=True,check=False)
                checkpoint(rawpath,{'command':cmd,'stdout':result.stdout,'stderr':result.stderr,'supervision':result.supervision})
            raw=cert.read(rawpath);lines=raw['stdout'].splitlines()
            expected=[f"R {i} {r[label][0]} {r[label][1]}" for i,r in enumerate(f['fixtures'])]+[f"S 8192 {meta['primes']}"]
            if raw['command']!=cmd or raw['stderr'] or raw['supervision']['outcome']!='completed' or raw['supervision']['returncode']!=0 or lines!=expected:
                raise ArithmeticError('compiled cached sums disagree with saved exact scalar/canonical proofs')
            records.append({'family':f['family'],'label':label,'candidates_checked':8192,'raw_sha256':cert.hashed(rawpath),'wall_seconds':raw['supervision']['wall_seconds']})
    result={'schema':'elliptic-curves.mw16-joint-binary-caches-result.v1','status':'PASS','sources':sources(),'protocol_sha256':cert.hashed(D/'protocol.json'),'encoding':encoded,'comparisons':records,'candidates_checked':40960,'component_sums_checked':81920,'claim_boundary':p['scope']}
    if check:
        if cert.read(OUT)!=result: raise ArithmeticError('complete score comparison replay differs')
        checkpoint(D/'score-replay.json',{'status':'PASS','certificate_sha256':cert.hashed(OUT),'sources':sources()})
    else: checkpoint(OUT,result)
    print('ALL40960 SAVED MW16 SHORT/EXTENDED SCORES MATCH',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['prepare','encode','encoding-check','score','score-check']);a=parser.parse_args()
    if a.stage=='encoding-check': encode(True)
    elif a.stage=='score-check': score(True)
    else: globals()[a.stage]()
