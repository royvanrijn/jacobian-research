#!/usr/bin/env sage-python
"""Frozen balanced H1024 pilot on six newly compactified existing R17 families."""
import argparse
from importlib.machinery import SourceFileLoader
from pathlib import Path
from fractions import Fraction as F
from math import log
import sys
from sage.all import EllipticCurve,QQ,prime_range

ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path[:0]=[str(CAS),str(ROOT/'elliptic-curves')]
import certify_compact_r17_candidates as cert
import compact_atlas_specialization as specialization
from research_runtime.store import checkpoint,digest
from research_runtime.supervisor import capture,Limits
from research_runtime.finite_reduction import ReductionCache
from research_runtime.memory_store import MemoryFactStore
from research_runtime.search_state import raw_state
from research_runtime.mw_state import MWState
from pointed_quartic_search import PointedQuarticSearch,sources as pointed_sources
old=SourceFileLoader('atlas_search_common',str(CAS/'compact_r17_prospective.sage')).load_module()
wide=SourceFileLoader('atlas_scanner_parse',str(CAS/'compact_r17_wide.sage')).load_module()
engine=SourceFileLoader('atlas_height_geometry',str(CAS/'half_lattice_fake_descent_replay.sage')).load_module()
DIRECTORY=ROOT/'artifacts/local/elliptic-curves/compact-six-r17-h1024-v2'
HISTORY=ROOT/'artifacts/generated-results/elliptic-curves/r17_refresh_jump_ladder_blind_v2.json'
BLIND_INPUT=ROOT/'elliptic-curves/data/r17_refresh_jump_ladder_blind_inputs_v1.json'
BINARY=ROOT/'artifacts/local/elliptic-curves/compact-r17-wide-v1/scanner'


def sources():
    paths=(Path(__file__).resolve(),specialization.ATLAS,Path(specialization.__file__).resolve(),
           CAS/'compact_r17_prospective.sage',CAS/'compact_r17_wide.sage',CAS/'half_lattice_fake_descent_replay.sage',wide.SCANNER,
           CAS/'research_runtime/memory_store.py',Path(cert.__file__).resolve(),cert.DATABASE)
    return {**pointed_sources(),**{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths}}


def prepare(directory):
    if (directory/'protocol.json').exists():raise FileExistsError('pilot already frozen')
    atlas=cert.read(specialization.ATLAS);history=cert.read(HISTORY)['results'];blind=cert.read(BLIND_INPUT)['cases'];masks={}
    for f in atlas['families']:
        entry=next(r for r in history if r['native_chart']=='norm12-orbit-'+f['family'])
        basis=next(r for r in blind if r['curve_id']==entry['curve_id'])
        if basis['generic_height_gram']!=f['generic_height_gram']:raise ArithmeticError('generic mask basis mismatch')
        records=entry['initial']['cover_records'];selected=sorted(r['mask'] for r in records)
        if len(set(selected))!=43:raise ArithmeticError('generic deep mask count changed')
        oracle=engine.CosetOracle(f['generic_height_gram'])
        for mask in selected:
            norm,rep,error=oracle.solve(mask)
            exact=sum(rep[i]*f['generic_height_gram'][i][j]*rep[j] for i in range(17) for j in range(17))
            if norm!=12 or exact!=12 or any((rep[i]-(mask>>i))%2 for i in range(17)):raise ArithmeticError('generic deep representative failed')
        masks[f['family']]=selected
    protocol={'schema':'elliptic-curves.compact-six-r17-h1024.v1','sources':sources(),
        'families':[f['family'] for f in atlas['families']],'height':1024,'prime_bound':4093,'retained_per_family':128,'finalists_per_family':4,
        'score':'all 562 primes 5..4093: rounded integer 1e12*(2-a_p)*log(p)/(p+1-a_p) at good residues; no prefix truncation',
        'order':'score descending, good-prime count descending, denominator, signed numerator',
        'generic_masks':masks,'generic_mask_source_hashes':{str(p.relative_to(ROOT)):cert.hashed(p) for p in (HISTORY,BLIND_INPUT)},
        'generic_mask_boundary':'Only the generic Gram and unordered 43 generic norm-12 masks are reused; no exceptional point, rank label, or specialization-dependent ordering is a selection input.',
        'scanner_binary_sha256':cert.hashed(BINARY),'table_worker_wall_seconds':300,'scanner_wall_seconds':120,
        'charts_per_fibre':43,'chart_height':100000,'seconds_per_chart':4,'admission_prime_bound':251,
        'coordinate_policy':{'kind':'metric','weight':'16'},'point_worker_wall_seconds':300,'point_worker_rss_bytes':1610612736,
        'maximum_table_workers':2,'maximum_point_workers':4,'maximum_point_workers_total':24,
        'target_rank_lower_bound':28,'record_target_rank_lower_bound':32,
        'scope':'Balanced pilot on newly compactified models. Post-freeze catalogue matching and curve deduplication precede point work. Known equations are skipped without refilling the fixed roster. Every miss stays bounded.'}
    checkpoint(directory/'protocol.json',protocol);print('FROZEN SIX-FAMILY PILOT',flush=True)


def bindings(directory):
    p=cert.read(directory/'protocol.json')
    if p['sources']!=sources():raise ArithmeticError('frozen search sources changed')
    return p


def select(directory,family):
    protocol=bindings(directory);folder=directory/family;output=folder/'population.json'
    if output.exists():print('RETAINED POPULATION',family,flush=True);return
    f=next(r for r in cert.read(specialization.ATLAS)['families'] if r['family']==family)
    model={k:f[k]+['0']*(n-len(f[k])) for k,n in [('A_coefficients_low_to_high',9),('B_coefficients_low_to_high',13)]}
    records=[]
    for p in prime_range(5,4094):
        p=int(p);path=folder/'trace-tables'/f'{p}.json';key={'family':family,'model_hash':digest(model),'prime':p}
        if path.exists():
            table=cert.read(path)
            if table['input']!=key:raise ArithmeticError('trace cache binding mismatch')
        else:
            traces,good=old.trace_table(model,p);table={'input':key,'traces':traces.tolist(),'good':good.tolist()};checkpoint(path,table)
        records.append(table)
        if p in (997,1999,2999,4093):print('ATLAS TABLE',family,p,flush=True)
    shards=[]
    for sign in (-1,1):
        table_path=folder/f'tables-{sign}.txt'
        if not table_path.exists():
            with table_path.open('x') as out:
                out.write('RATIONAL_NAGAO_LOCAL_TABLE_V1\nF COMPACT_ATLAS_'+family+' 8 12\n')
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
        if path.exists():shard=cert.read(path)
        else:
            if cert.hashed(BINARY)!=protocol['scanner_binary_sha256']:raise ArithmeticError('scanner binary changed')
            result=capture([str(BINARY),'%s'%table_path,'1024','1024','128','0','1'],limits=Limits(120,536870912),log_path=folder/f'scan-{sign}.log')
            rows,summary=wide.parse(result.stdout,sign)
            for r in rows:r['score_units']=r.pop('prefix_units');r.pop('extension_units');r.pop('extension_good')
            shard={'rows':rows,'summary':summary,'supervision':result.supervision,'table_sha256':cert.hashed(table_path),'protocol_hash':digest(protocol)};checkpoint(path,shard)
        if shard['protocol_hash']!=digest(protocol) or shard['table_sha256']!=cert.hashed(table_path):raise ArithmeticError('scanner binding changed')
        shards.append(shard)
    rows=[r for s in shards for r in s['rows']];rows.sort(key=lambda r:(-r['score_units'],-r['prefix_good'],r['denominator'],r['numerator']))
    if len(set(r['parameter'] for r in rows))!=len(rows):raise ArithmeticError('duplicate population row')
    checkpoint(output,{'family':family,'protocol_hash':digest(protocol),'candidate_count':sum(s['summary'][3] for s in shards),
        'retained_candidates':rows[:128],'finalists':rows[:4],
        'unused_H_band':'duplicate p=5; unused, not validation','public_equations_or_exceptional_points_used_for_selection':False})
    print('SELECTED',family,[(r['parameter'],r['score_units']/10**12) for r in rows[:4]],flush=True)


def run(directory,family,index):
    protocol=bindings(directory);population=cert.read(directory/family/'population.json')
    if population['protocol_hash']!=digest(protocol):raise ArithmeticError('population binding changed')
    t=population['finalists'][index]['parameter'];f=next(r for r in cert.read(specialization.ATLAS)['families'] if r['family']==family)
    folder=directory/family/f'candidate-{index:02}';output=folder/'result.json'
    if output.exists() and cert.read(output)['status']!='RUNNING':print('RETAINED',family,t,flush=True);return
    original,generic=specialization.specialize(f,t)
    matches=[r['id'] for r in cert.read(cert.DATABASE)['curves'] if cert.isomorphic(original,r['ainvs'])]
    if matches:
        checkpoint(output,{'status':'SKIPPED_CATALOGUED','family':family,'parameter':t,'matches':matches,'protocol_hash':digest(protocol)})
        print('KNOWN',family,t,matches,flush=True);return
    curve=EllipticCurve(QQ,list(map(lambda q:QQ(str(q)),original)));pts=[curve([QQ(str(x)),QQ(str(y))]) for x,y in generic]
    minimal=curve.local_data(2).minimal_model();iso=curve.isomorphism_to(minimal)
    a1,a2,a3,a4,a6=minimal.a_invariants();c2,c4,c6=a2+a1*a1/4,a4+a1*a3/2,a6+a3*a3/4
    model=tuple(F(str(v)) for v in (0,0,0,c4-c2*c2/3,c6-c2*c4/3+2*c2**3/27))
    points=[]
    for P in pts:
        x,y=iso(P).xy();points.append((F(str(x+c2/3)),F(str(y+(a1*x+a3)/2))))
    specialization.family_check(f,t,model,points)
    cache=ReductionCache(MemoryFactStore());initial=raw_state(model,points,cache=cache,prime_bound=1000)
    if initial.rank!=17:raise ArithmeticError('generic independence not certified')
    if output.exists():
        result=cert.read(output)
        if result['initial_state']!=initial.record() or result['protocol_hash']!=digest(protocol):raise ArithmeticError('resume mismatch')
        cache.store.import_snapshot(result['arithmetic_facts']);state=MWState.from_record(result['final_state'],cache=cache)
    else:
        gram=engine.canonical_height_gram(model,points);rounded=[[int((x*1000000).to_integral_value()) for x in row] for row in gram]
        oracle=engine.CosetOracle(rounded);centres=[]
        for mask in protocol['generic_masks'][family]:
            norm,rep,error=oracle.solve(mask);centres.append({'mask':mask,'representative':list(rep),'metric_norm':norm})
        centres.sort(key=lambda r:(-r['metric_norm'],r['mask']));state=initial
        result={'schema':'elliptic-curves.compact-atlas-point-search.v1','family':family,'parameter':t,'protocol_hash':digest(protocol),
            'curve':list(map(str,model)),'generic_points':[list(map(str,p)) for p in points],
            'initial_state':initial.record(),'final_state':initial.record(),'arithmetic_facts':cache.store.snapshot(),
            'metric_gram':[[str(x) for x in row] for row in gram],'centres':centres,'charts':[],'status':'RUNNING','rank_lower_bound':17}
        checkpoint(output,result)
    primes=tuple(map(int,prime_range(3,252)))
    for i in range(len(result['charts']),43):
        centre=result['centres'][i];rep=centre['representative']+[0]*(state.rank-17)
        search=PointedQuarticSearch(state=state,centre={'coefficients':rep},coordinate_policy=protocol['coordinate_policy'])
        outcome=search.search(100000,4,checkpoint_dir=folder/'charts'/state.key)
        for point in outcome.curve_points:state=state.adjoin(point,cache=cache,extra_primes=primes)
        result['charts'].append({'centre':centre,'search':outcome.record,'rank_lower_bound':state.rank,'admission_prime_bound':251})
        result.update(final_state=state.record(),arithmetic_facts=cache.store.snapshot(),rank_lower_bound=state.rank);checkpoint(output,result)
        print('ATLAS POINTS',family,t,'chart',i+1,'rank',state.rank,flush=True)
        if state.rank>=28:result['status']='TARGET_REACHED_PENDING_INDEPENDENT_REPLAY';checkpoint(output,result);return
    result['status']='COMPLETE_DECLARED_PILOT';checkpoint(output,result)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','select','run']);p.add_argument('--directory',type=Path,default=DIRECTORY);p.add_argument('--family',default='103b2');p.add_argument('--index',type=int,default=0);a=p.parse_args()
    if a.stage=='prepare':prepare(a.directory)
    elif a.stage=='select':select(a.directory,a.family)
    else:run(a.directory,a.family,a.index)
