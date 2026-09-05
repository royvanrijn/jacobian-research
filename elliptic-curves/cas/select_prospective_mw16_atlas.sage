#!/usr/bin/env sage-python
"""Target-free finite generic census and H1024 selector for five compact MW16 families."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import lcm, log
from pathlib import Path
import sys
import numpy as np
from sage.all import EllipticCurve, GF, prime_range

ROOT=Path(__file__).resolve().parents[2]; CAS=ROOT/'elliptic-curves/cas'
sys.path.insert(0,str(CAS))
import compact_mw16_specialization as spec
from research_runtime.store import checkpoint, digest
from research_runtime.supervisor import capture, Limits
geometry=SourceFileLoader('prospective_geometry',str(CAS/'prospective_half_lattice.sage')).load_module()
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/prospective-mw16-h1024-v1'
BINARY=ROOT/'artifacts/local/elliptic-curves/compact-r17-wide-v1/scanner'
SCANNER=CAS/'newfamily/scan_rational_nagao_tables.cpp'

def read(p): return json.loads(p.read_text())
def hashed(p): return sha256(p.read_bytes()).hexdigest()
def sources():
    paths=[Path(__file__).resolve(),spec.ATLAS,Path(spec.__file__),CAS/'certify_compact_r17_candidates.py',
           CAS/'prospective_half_lattice.sage',SCANNER,CAS/'research_runtime/store.py',CAS/'research_runtime/supervisor.py']
    return {str(p.relative_to(ROOT)):hashed(p) for p in paths}
def bindings(directory):
    p=read(directory/'protocol.json')
    if p['sources']!=sources() or p['scanner_binary_sha256']!=hashed(BINARY): raise ArithmeticError('frozen source changed')
    return p
def family_record(family): return next(f for f in read(spec.ATLAS)['families'] if f['fibration_id']==family)

def prepare(directory):
    if (directory/'protocol.json').exists(): raise FileExistsError('selection protocol already frozen')
    if hashed(spec.ATLAS)!='88a8c61c2722394bb812da601e24ea2c19e54361eef4b0fe9c7924e78c828801': raise ArithmeticError('atlas handoff changed')
    if hashed(Path(spec.__file__))!='f6aeab97b980bc6e9a5db71675d4f47fbf3ac15cac5f35f2ea310a77e2a8a7c0': raise ArithmeticError('helper handoff changed')
    checkpoint(directory/'protocol.json',{
        'schema':'elliptic-curves.prospective-mw16-selection.v1','sources':sources(),'scanner_binary_sha256':hashed(BINARY),
        'families':[f['fibration_id'] for f in read(spec.ATLAS)['families']],
        'height':1024,'prime_bound':4093,'retained_per_family':128,'finalists_per_family':4,
        'score':'All562 primes5..4093, sum round(1e12*(2-a_p)*log(p)/(p+1-a_p)) at good residues.',
        'order':'score descending, good prime count descending, denominator ascending, signed numerator ascending',
        'generic_census':{'classes_per_family':65536,'dimension':16,'retained_nonzero_classes':43,
                          'order':'computed representative norm descending, mask ascending','wall_seconds_per_family':300,
                          'rss_bytes':1073741824,'maximum_parallel_workers':2,'checkpoint_every':2048},
        'selection_worker_wall_seconds':300,'selection_worker_rss_bytes':1073741824,'maximum_selection_workers':2,
        'scanner_wall_seconds':120,'population':'Signed primitive nonzero n/d with abs(n),d<=1024 and d>=1; zero and infinity excluded.',
        'endpoint':'Candidate incidence by score and generic geometry only; no point or rank conclusion.',
        'point_followup_scope':'Separate frozen protocol, fixed20 finalists,43 charts each at height100000 and4 seconds; workers300seconds/1.5GiB, at most4 concurrent.',
        'target_free_boundary':'No catalogue equations, known parameters, target j values, public points, ranks or jump labels are read by selection or prospective point workers. Novelty checks occur only after the fixed point batch.',
        'failure_semantics':'Timeouts and incomplete generic independence are censored attempts, never rank upper bounds.',
        'software':{'sage':str(__import__('sage.env',fromlist=['SAGE_VERSION']).SAGE_VERSION),'numpy':np.__version__}})
    print('FROZEN FIVE-FAMILY TARGET-FREE SELECTOR',flush=True)

def census(directory,family):
    protocol=bindings(directory);f=family_record(family);path=directory/family/'generic-census.json'
    gram=[[F(x) for x in row] for row in f['generic_height_gram']];scale=lcm(*(x.denominator for row in gram for x in row))
    integers=[[int(x*scale) for x in row] for row in gram];oracle=geometry.CosetOracle(integers)
    data=read(path) if path.exists() else {'protocol_hash':digest(protocol),'family':family,'gram':f['generic_height_gram'],
        'integer_gram_scale':scale,'records':[],'status':'RUNNING','claim':'All parities enumerated; floating CVP optimality is not an exact covering-radius theorem.'}
    if data['protocol_hash']!=digest(protocol) or data['gram']!=f['generic_height_gram']: raise ArithmeticError('census binding mismatch')
    for mask in range(len(data['records']),65536):
        norm,rep,error=oracle.solve(tuple((mask>>j)&1 for j in range(16)))
        if len(rep)!=16 or any((rep[j]-(mask>>j))%2 for j in range(16)): raise ArithmeticError('generic parity mismatch')
        data['records'].append({'mask':mask,'norm':str(F(norm,scale)),'representative':list(rep),'cvp_error':error})
        if (mask+1)%2048==0: checkpoint(path,data);print('MW16 CENSUS',family,mask+1,flush=True)
    selected=sorted(data['records'][1:],key=lambda r:(-F(r['norm']),r['mask']))[:43]
    data.update(status='COMPLETE_DECLARED_CENSUS',selected=selected,norm_histogram=dict(Counter(r['norm'] for r in data['records'])))
    checkpoint(path,data);print('MW16 GENERIC CLASSES',family,data['norm_histogram'],flush=True)

def trace_table(model,p):
    if not 3<p<=4093: raise ValueError('prime outside int64-safe protocol')
    t=np.arange(p,dtype=np.int64)
    def evaluate(key):
        cs=[int(F(c).numerator*pow(F(c).denominator,-1,p)%p) for c in model[key]]
        v=np.zeros(p,dtype=np.int64)
        for c in reversed(cs): v=(v*t+c)%p
        return np.append(v,cs[-1])
    aa,bb=evaluate('A_coefficients_low_to_high'),evaluate('B_coefficients_low_to_high')
    chi=np.full(p,-1,dtype=np.int64);chi[(t*t)%p]=1;chi[0]=0;cube=(t*t%p)*t%p
    traces=np.array([-int(chi[(cube+a*t+b)%p].sum()) for a,b in zip(aa,bb)],dtype=np.int64)
    good=(4*((aa*aa%p)*aa%p)+27*(bb*bb%p))%p!=0
    return traces,good

def parse(text,sign):
    rows=[];summary=None
    for line in text.splitlines():
        parts=line.split()
        if not parts: continue
        if parts[0]=='C':
            _,n,d,a,b,g,h=parts;n=sign*int(n);d=int(d)
            rows.append({'numerator':n,'denominator':d,'parameter':str(F(n,d)),
                         'score_units':int(F(a)*10**12),'good_primes':int(g)})
        if parts[0]=='S': summary=list(map(int,parts[1:]))
    if summary is None: raise ArithmeticError('scanner did not finish')
    return rows,summary

def select(directory,family):
    protocol=bindings(directory);folder=directory/family;output=folder/'population.json'
    if output.exists(): raise FileExistsError('population already retained')
    f=family_record(family);model={k:f[k]+['0']*(n-len(f[k])) for k,n in [('A_coefficients_low_to_high',9),('B_coefficients_low_to_high',13)]}
    records=[];checks=[]
    for p in map(int,prime_range(5,4094)):
        path=folder/'trace-tables'/f'{p}.json';key={'family':family,'model_hash':digest(model),'prime':p}
        if path.exists():
            table=read(path)
            if table['input']!=key: raise ArithmeticError('trace table binding mismatch')
        else:
            traces,good=trace_table(model,p);table={'input':key,'traces':traces.tolist(),'good':good.tolist()};checkpoint(path,table)
        if p in (5,7,11,13):
            for t in range(p+1):
                a,b=[sum((F(c).numerator*pow(F(c).denominator,-1,p)%p)*pow(t,i,p) for i,c in enumerate(model[k]))%p if t<p else int(F(model[k][-1]))%p for k in model]
                good=(4*a**3+27*b*b)%p!=0
                if good!=table['good'][t]: raise ArithmeticError('singularity flag failed')
                if good and p+1-int(EllipticCurve(GF(p),[a,b]).cardinality())!=table['traces'][t]: raise ArithmeticError('independent cardinality check failed')
                checks.append([p,t,good])
        records.append(table)
        if p in (997,1999,2999,4093): print('MW16 TABLE',family,p,flush=True)
    shards=[]
    for sign in (-1,1):
        table_path=folder/f'tables-{sign}.txt'
        if not table_path.exists():
            with table_path.open('x') as out:
                out.write('RATIONAL_NAGAO_LOCAL_TABLE_V1\nF PROSPECTIVE_MW16 8 12\n')
                for label,band in [('D',records),('H',records[:1])]:
                    out.write(f'B {label} {len(band)}\n')
                    for r in band:
                        p=r['input']['prime'];out.write(f'P {p}\n')
                        for i in range(p+1):
                            j=p if i==p else sign*i%p;a=r['traces'][j];good=r['good'][j]
                            units=round((2-a)/(p+1-a)*log(p)*10**12) if good else 0
                            out.write(f'{int(good)} {a} {units}\n')
                out.write('END\n')
        path=folder/f'scan-{sign}.json'
        if path.exists(): shard=read(path)
        else:
            result=capture([str(BINARY),str(table_path),'1024','1024','128','0','1'],limits=Limits(120,536870912),log_path=folder/f'scan-{sign}.log')
            rows,summary=parse(result.stdout,sign)
            shard={'rows':rows,'summary':summary,'supervision':result.supervision,'table_sha256':hashed(table_path),'protocol_hash':digest(protocol)};checkpoint(path,shard)
        if shard['protocol_hash']!=digest(protocol) or shard['table_sha256']!=hashed(table_path): raise ArithmeticError('scanner binding changed')
        shards.append(shard)
    rows=[r for s in shards for r in s['rows']];rows.sort(key=lambda r:(-r['score_units'],-r['good_primes'],r['denominator'],r['numerator']))
    # Signed numerator tie ordering differs within the negative C++ shard.
    # Strict score separation makes that irrelevant to the fixed final four.
    if not all(rows[3]['score_units']>s['rows'][-1]['score_units'] for s in shards): raise ArithmeticError('finalist boundary needs a tie-complete rescan')
    checkpoint(output,{'protocol_hash':digest(protocol),'family':family,'candidate_count':sum(s['summary'][3] for s in shards),
                       'retained_candidates':rows[:128],'finalists':rows[:4],'independent_small_prime_checks':checks,
                       'unused_H_band':'duplicate p5; unused and not validation','target_free':True})
    print('MW16 SELECTED',family,[(r['parameter'],r['score_units']/10**12) for r in rows[:4]],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','census','select']);p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--family');a=p.parse_args()
    if a.stage=='prepare': prepare(a.directory)
    elif a.stage=='census': census(a.directory,a.family)
    else: select(a.directory,a.family)
