#!/usr/bin/env sage-python
"""Backfill/replay exact retained obstructions for the frozen MW18 census.

Only the recorded obstruction/validation primes are visited. Discovery retains
shared group-law DAGs and chord numerator identities. Replay uses polynomial
cross multiplication and the stated finite nonresidues, never a new sieve.
Each prime has an atomic resumable checkpoint and a bounded supervisor.
"""
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import runpy
import sys
import time
from functools import lru_cache

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
from research_runtime.function_field_witness import GroupLawDiscovery,chord_witness,replay_group_law,replay_chord
from research_runtime.store import atomic_write,checkpoint
from research_runtime.supervisor import Limits,run
from research_runtime.witnesses import retained_source

SOURCE=ROOT/'elkies-k3/scripts/certify_r17_extreme_anchored_mw18_covers.sage'
CERTIFICATE=ROOT/'artifacts/generated-results/elkies-k3-r17-extreme-anchored-mw18-covers-v1.json'
SCHEMA='elliptic-curves.mw18-census-prime-witness.v3'
BATCH_SIZE=2000


def save(path,row):
    atomic_write(path,gzip.compress((json.dumps(row,sort_keys=True,separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))

def read(path):return json.loads(gzip.decompress(path.read_bytes()))


def prime_indices(chart):
    sieve=chart['sieve'];result={}
    arrays=[sieve['global_validation_primes_in_priority_order'],*sieve['obstruction_primes_in_priority_order'].values()]
    if any(len(a)!=39120 for a in arrays):raise ArithmeticError('wrong retained census dimension')
    for i in range(39120):
        for p in {a[i] for a in arrays}-{0}:result.setdefault(p,[]).append(i)
    return result


@lru_cache(maxsize=2)
def common_setup(chart_key):
    h=runpy.run_path(str(SOURCE));certificate=json.loads(CERTIFICATE.read_text())
    chart=next(c for c in certificate['charts'] if c['chart']=='norm12-orbit-'+chart_key)
    for name,expected in chart['inputs'].items():
        pinned=retained_source(ROOT,name,expected)
        if Path(name).suffix not in ('.py','.sage') and pinned.resolve()!=(ROOT/name).resolve():
            raise ArithmeticError('active mathematical input differs from the pinned census')
    direct=json.loads(h['CHARTS'][chart_key]['direct'].read_text())
    priority=h['parse_priority_rows'](h['PRIORITY']);words=h['transported_words'](direct,priority)[-1]
    helper=runpy.run_path(str(h['CHORD_HELPER']))
    return h,chart,direct,words,helper


def setup(chart_key,prime):
    h,chart,direct,words,helper=common_setup(chart_key)
    finite,R,K,A,B,delta,curve,_=h['build_modular_context'](direct,[],prime)
    basis=[curve(h['parse_rational_function'](r['X'],R,K),h['parse_rational_function'](r['Y'],R,K))
           for r in direct['sections']['records']]
    parameters={f['curve_id']:h['reduce_rational'](f['native_parameter'],finite) for f in chart['fibres']}
    return h,chart,words,R,K,A,B,delta,curve,basis,helper,parameters


@lru_cache(maxsize=2)
def short_transport(chart_key):
    from sage.all import matrix,ZZ
    h,chart,direct,words,helper=common_setup(chart_key)
    priority=h['parse_priority_rows'](h['PRIORITY'])
    coefficients=[tuple(map(int,r['published_basis_w'].split())) for r in priority]
    U=matrix(ZZ,coefficients);V=matrix(ZZ,words);rows=U.transpose().pivots()
    if len(rows)!=17:raise ArithmeticError('short trace coordinates do not span')
    T=U.matrix_from_rows(rows).solve_right(V.matrix_from_rows(rows))
    if any(c not in ZZ for c in T.list()) or U*T!=V:raise ArithmeticError('short basis transport failed')
    return coefficients,[tuple(map(int,r)) for r in T.rows()]


def witness_identity(chart_key,prime,chart,batch_index=None):
    return {'schema':SCHEMA,'certificate_sha256':sha256(CERTIFICATE.read_bytes()).hexdigest(),
            'chart':chart_key,'prime':prime,'inputs':chart['inputs'],'batch_index':batch_index,'batch_size':BATCH_SIZE}


def discovery(chart_key,prime,path,max_roots=None,batch_index=None):
    h,chart,words,R,K,A,B,delta,curve,basis,helper,parameters=setup(chart_key,prime)
    indices=selected_indices(chart,prime,batch_index);identity=witness_identity(chart_key,prime,chart,batch_index)
    old=read(path) if path.exists() else None
    if old and old['identity']!=identity:raise ArithmeticError('checkpoint inputs differ')
    graph=GroupLawDiscovery(curve,basis,retained=old['nodes'] if old else None)
    roots=old['roots'] if old else []
    short_words,transport=short_transport(chart_key);short_basis=[graph.trace(w) for w in transport]
    if [r['index'] for r in roots]!=indices[:len(roots)]:raise ArithmeticError('checkpoint roots do not form the required prefix')
    started=time.monotonic();last=started;new=0
    for index in indices[len(roots):]:
        if max_roots is not None and new>=max_roots:break
        node=graph.trace_in_basis(short_words[index],short_basis);frame=chord_witness(graph.points[node],A,R,K,helper)
        # Discovery checks the same explicit identities that portable replay uses.
        from research_runtime.function_field_witness import read_point,point_record
        q=replay_chord(R,A,delta,read_point(R,point_record(graph.points[node])),frame)
        ids=[int(cid) for cid,values in chart['sieve']['obstruction_primes_in_priority_order'].items() if values[index]==prime]
        for cid in ids:
            value=q(parameters[cid])
            if not value or value.is_square():raise ArithmeticError('the retained prime is not a nonresidue witness')
        roots.append({'index':index,'node':node,'chord':frame,'obstructed_curves':ids});new+=1
        if len(roots)%1000==0 or time.monotonic()-last>=40:
            save(path,{'identity':identity,'nodes':graph.nodes,'roots':roots,'status':'PARTIAL'})
            print(f'CENSUS_WITNESS|chart={chart_key}|p={prime}|roots={len(roots)}/{len(indices)}|nodes={len(graph.nodes)}|seconds={time.monotonic()-started:.2f}',flush=True)
            last=time.monotonic()
    status='COMPLETE' if len(roots)==len(indices) else 'PARTIAL'
    save(path,{'identity':identity,'nodes':graph.nodes,'roots':roots,'status':status})
    report={'status':status,'root_count':len(roots),'required_root_count':len(indices),'node_count':len(graph.nodes),
            'wall_seconds':time.monotonic()-started,'witness_sha256':sha256(path.read_bytes()).hexdigest()}
    checkpoint(path.with_suffix('.summary.json'),report)
    print(f'CENSUS_WITNESS|chart={chart_key}|p={prime}|status={status}|roots={len(roots)}|nodes={len(graph.nodes)}|seconds={report["wall_seconds"]:.2f}',flush=True)
    return report


def replay(chart_key,prime,path,allow_partial=False,batch_index=None):
    h,chart,words,R,K,A,B,delta,curve,basis,helper,parameters=setup(chart_key,prime)
    row=path if isinstance(path,dict) else read(path);indices=selected_indices(chart,prime,batch_index)
    if row['identity']!=witness_identity(chart_key,prime,chart,batch_index):raise ArithmeticError('witness inputs differ')
    if row['status']!='COMPLETE' and not allow_partial:raise ArithmeticError('incomplete prime witness')
    expected=indices[:len(row['roots'])] if allow_partial else indices
    if [r['index'] for r in row['roots']]!=expected:raise ArithmeticError('missing or duplicate retained census roots')
    started=time.monotonic();points,labelled_words=replay_group_law(R,A,B,basis,row['nodes'])
    count=0
    for root in row['roots']:
        index,node=root['index'],root['node']
        if type(node) is not int or not 0<=node<len(points) or labelled_words[node]!=words[index]:
            raise ArithmeticError('trace root differs from the transported priority word')
        q=replay_chord(R,A,delta,points[node],root['chord'])
        ids=[int(cid) for cid,values in chart['sieve']['obstruction_primes_in_priority_order'].items() if values[index]==prime]
        if ids!=root['obstructed_curves']:raise ArithmeticError('obstruction scope differs')
        for cid in ids:
            value=q(parameters[cid])
            if not value or value.is_square():raise ArithmeticError('retained finite nonresidue failed')
            count+=1
    result={'chart':chart_key,'prime':prime,'batch_index':batch_index,'verified_roots':len(expected),'verified_nonresidues':count,
            'group_nodes':len(points),'wall_seconds':time.monotonic()-started,'status':'PASS' if not allow_partial else 'PARTIAL_CONTROL_PASS'}
    print('CENSUS_REPLAY|'+json.dumps(result,sort_keys=True),flush=True);return result


def selected_indices(chart,prime,batch_index):
    indices=prime_indices(chart)[prime]
    if batch_index is None:return indices
    if type(batch_index) is not int or not 0<=batch_index<(len(indices)+BATCH_SIZE-1)//BATCH_SIZE:
        raise ValueError('invalid census batch index')
    return indices[batch_index*BATCH_SIZE:(batch_index+1)*BATCH_SIZE]


def census_jobs(certificate):
    return [(len(indices),chart['chart'].split('-')[-1],prime,batch)
            for chart in certificate['charts'] for prime,indices in prime_indices(chart).items()
            for batch in range((len(indices)+BATCH_SIZE-1)//BATCH_SIZE)]


def batch_path(directory,key,prime,batch):return directory/f'{key}-{prime}-{batch:03d}.json.gz'


def backfill_all(directory,workers,wall,rss):
    from concurrent.futures import ThreadPoolExecutor,as_completed
    certificate=json.loads(CERTIFICATE.read_text());jobs=census_jobs(certificate)
    jobs.sort(reverse=True)
    directory.mkdir(parents=True,exist_ok=True);records=[]
    def dispatch(job):
        count,key,prime,batch=job;path=batch_path(directory,key,prime,batch);summary=path.with_suffix('.summary.json')
        if summary.exists() and path.exists():
            old=json.loads(summary.read_text())
            if old['status']=='COMPLETE' and old['witness_sha256']==sha256(path.read_bytes()).hexdigest():
                return {'chart':key,'prime':prime,'batch_index':batch,'status':'COMPLETE','summary':old,'witness':str(path)}
        log=path.with_suffix('.discovery.log')
        receipt=run([sys.executable,str(Path(__file__).resolve()),'--worker','--chart',key,'--prime',str(prime),'--batch-index',str(batch),'--output',str(path)],
                    limits=Limits(wall,rss),log_path=log,checkpoint_path=log.with_suffix('.supervisor.json'))
        summary_value=json.loads(summary.read_text()) if receipt['outcome']=='completed' and summary.exists() else None
        return {'chart':key,'prime':prime,'batch_index':batch,'status':summary_value['status'] if summary_value else 'INCOMPLETE',
                'summary':summary_value,'supervisor':receipt,'witness':str(path)}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(dispatch,job) for job in jobs]
        for future in as_completed(futures):
            row=future.result();records.append(row)
            checkpoint(directory/'backfill.json',{'records':sorted(records,key=lambda r:(r['chart'],r['prime'],r['batch_index'])),
                'status':'COMPLETE' if len(records)==len(jobs) and all(r['status']=='COMPLETE' for r in records) else 'PARTIAL'})
            print(f"CENSUS_BACKFILL|jobs={len(records)}/{len(jobs)}|chart={row['chart']}|p={row['prime']}|status={row['status']}",flush=True)
    if any(r['status']!='COMPLETE' for r in records):raise SystemExit('incomplete prime witnesses remain checkpointed')


BUNDLE=ROOT/'artifacts/generated-results/elliptic-curves/runtime_mw18_census_witnesses_v1.zip'
BUNDLE_SCHEMA='elliptic-curves.mw18-census-witness-bundle.v1'


def check_inventory(certificate,selected):
    counts={'global_validations':0,'nonresidues':0,'positive_covers':0}
    if len({c['chart'] for c in certificate['charts']})!=len(certificate['charts']):raise ArithmeticError('duplicate chart')
    for key in selected:
        h,chart,direct,words,helper=common_setup(key);sieve=chart['sieve']
        if sieve['exact_global_validation_indices_zero_based'] or any(v==0 for v in sieve['global_validation_primes_in_priority_order']):
            raise ArithmeticError('an exact-only validation needs a separate retained witness')
        counts['global_validations']+=len(sieve['global_validation_primes_in_priority_order'])
        for fibre in chart['fibres']:
            cid=str(fibre['curve_id']);obstructions=sieve['obstruction_primes_in_priority_order'][cid]
            if sieve['exact_nonsquare_survivor_indices_zero_based'][cid] or sieve['exact_ramified_indices_zero_based'][cid]:
                raise ArithmeticError('an exact-only survivor needs a separate retained witness')
            positives=[cover['priority_index_zero_based'] for cover in fibre['covers']]
            if len(set(positives))!=len(positives) or sorted(positives)!=[i for i,p in enumerate(obstructions) if not p]:
                raise ArithmeticError('positive and nonresidue witnesses do not partition the census')
            for cover in fibre['covers']:
                if tuple(cover['equation_basis_trace_word'])!=words[cover['priority_index_zero_based']]:
                    raise ArithmeticError('positive cover differs from its inventory trace word')
            counts['nonresidues']+=sum(bool(p) for p in obstructions);counts['positive_covers']+=len(positives)
    return counts


def pack_bundle(directory,output):
    import zipfile
    certificate=json.loads(CERTIFICATE.read_text());jobs=sorted(census_jobs(certificate),key=lambda j:j[1:]);entries=[]
    sources={}
    for source in (Path(__file__).resolve(),SOURCE,ROOT/'elliptic-curves/cas/research_runtime/function_field_witness.py'):
        data=source.read_bytes();sources[str(source.relative_to(ROOT))]={'sha256':sha256(data).hexdigest(),'data':data}
    output.parent.mkdir(parents=True,exist_ok=True);temporary=output.with_suffix('.tmp')
    try:
        with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_STORED) as bundle:
            def write(name,data):
                info=zipfile.ZipInfo(name,date_time=(2026,9,5,0,0,0));bundle.writestr(info,data)
            for count,key,prime,batch in jobs:
                path=batch_path(directory,key,prime,batch);data=path.read_bytes();summary=json.loads(path.with_suffix('.summary.json').read_text())
                if summary['status']!='COMPLETE' or summary['witness_sha256']!=sha256(data).hexdigest():raise ArithmeticError('incomplete or changed discovery witness')
                name='witnesses/'+path.name;write(name,data)
                entries.append({'chart':key,'prime':prime,'batch_index':batch,'member':name,'sha256':sha256(data).hexdigest()})
            for name,row in sources.items():write('sources/'+row['sha256']+'/'+Path(name).name,row['data'])
            manifest={'schema':BUNDLE_SCHEMA,'certificate_sha256':sha256(CERTIFICATE.read_bytes()).hexdigest(),
                      'batch_size':BATCH_SIZE,'entries':entries,'implementation_sources_at_packaging':{k:v['sha256'] for k,v in sources.items()}}
            write('manifest.json',(json.dumps(manifest,sort_keys=True,separators=(',',':'))+'\n').encode())
        temporary.replace(output)
    finally:
        if temporary.exists():temporary.unlink()
    checkpoint(output.with_suffix('.json'),{'schema':BUNDLE_SCHEMA,'status':'RETAINED_AWAITING_EXACT_REPLAY',
        'bundle':str(output.relative_to(ROOT)) if output.is_relative_to(ROOT) else str(output),
        'bundle_sha256':sha256(output.read_bytes()).hexdigest(),'batch_count':len(entries)})


def replay_bundle(path=BUNDLE,selected=('07ca9','08234')):
    import gc
    import zipfile
    certificate=json.loads(CERTIFICATE.read_text());counts=check_inventory(certificate,selected);started=time.monotonic();reports=[]
    expected_jobs=sorted((key,prime,batch) for count,key,prime,batch in census_jobs(certificate))
    with zipfile.ZipFile(path) as bundle:
        manifest=json.loads(bundle.read('manifest.json'))
        if manifest['schema']!=BUNDLE_SCHEMA or manifest['certificate_sha256']!=sha256(CERTIFICATE.read_bytes()).hexdigest() or manifest['batch_size']!=BATCH_SIZE:
            raise ArithmeticError('census bundle identity differs')
        actual=[(r['chart'],r['prime'],r['batch_index']) for r in manifest['entries']]
        if actual!=expected_jobs:raise ArithmeticError('missing or duplicate prime/root batches')
        names=['manifest.json']+[r['member'] for r in manifest['entries']]
        for name,expected in manifest['implementation_sources_at_packaging'].items():
            member='sources/'+expected+'/'+Path(name).name;names.append(member)
            if sha256(bundle.read(member)).hexdigest()!=expected:raise ArithmeticError('generation source integrity failed')
        if sorted(bundle.namelist())!=sorted(names):raise ArithmeticError('unexpected or duplicate bundle members')
        for entry in manifest['entries']:
            key,prime,batch=entry['chart'],entry['prime'],entry['batch_index']
            if entry['member']!='witnesses/'+batch_path(Path('.'),key,prime,batch).name:raise ArithmeticError('unexpected witness member')
            if key not in selected:continue
            data=bundle.read(entry['member'])
            if sha256(data).hexdigest()!=entry['sha256']:raise ArithmeticError('retained witness integrity failed')
            row=json.loads(gzip.decompress(data))
            reports.append(replay(key,prime,row,batch_index=batch));del row,data;gc.collect()
    if sum(r['verified_nonresidues'] for r in reports)!=counts['nonresidues']:raise ArithmeticError('nonresidue coverage differs')
    return {'schema':BUNDLE_SCHEMA,'status':'PASS_EXACT_CENSUS_NONRESIDUE_AND_VALIDATION_WITNESSES',
            'certificate_sha256':sha256(CERTIFICATE.read_bytes()).hexdigest(),'bundle_sha256':sha256(path.read_bytes()).hexdigest(),
            'counts':counts,'batches':reports,'wall_seconds':time.monotonic()-started,
            'claim_boundary':'Every recorded finite validation and negative obstruction is checked by polynomial group-law and chord identities; positive covers require the companion exact cover replayer.'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--chart',choices=['07ca9','08234']);p.add_argument('--prime',type=int)
    p.add_argument('--bundle',type=Path);p.add_argument('--pack-directory',type=Path)
    p.add_argument('--all-directory',type=Path);p.add_argument('--workers',type=int,default=2)
    p.add_argument('--output',type=Path);p.add_argument('--verify',action='store_true')
    p.add_argument('--batch-index',type=int);p.add_argument('--max-roots',type=int);p.add_argument('--allow-partial-control',action='store_true')
    p.add_argument('--wall-seconds',type=int,default=1200);p.add_argument('--rss-bytes',type=int,default=1073741824)
    p.add_argument('--worker',action='store_true',help=argparse.SUPPRESS)
    a=p.parse_args()
    if a.pack_directory:
        pack_bundle(a.pack_directory,a.bundle or BUNDLE);return
    if a.bundle and a.verify:
        if a.worker:
            result=replay_bundle(a.bundle)
            if a.output:checkpoint(a.output,result)
            print(result['status'],flush=True)
        else:
            log=(a.output or a.bundle).with_suffix('.replay.log')
            result=run([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:],'--worker'],limits=Limits(a.wall_seconds,a.rss_bytes),log_path=log,checkpoint_path=log.with_suffix('.supervisor.json'))
            print(log.read_text(),end='')
            if result['outcome']!='completed':raise SystemExit(result['outcome'])
        return
    if a.all_directory:
        if not 1<=a.workers<=2:p.error('one or two bounded workers are supported')
        backfill_all(a.all_directory,a.workers,a.wall_seconds,a.rss_bytes);return
    if not a.chart or not a.prime or not a.output:p.error('--chart, --prime and --output are required')
    if a.worker:
        if a.verify:replay(a.chart,a.prime,a.output,a.allow_partial_control,a.batch_index)
        else:discovery(a.chart,a.prime,a.output,a.max_roots,a.batch_index)
        return
    log=a.output.with_suffix('.replay.log' if a.verify else '.discovery.log')
    result=run([sys.executable,str(Path(__file__).resolve()),*sys.argv[1:],'--worker'],
               limits=Limits(a.wall_seconds,a.rss_bytes),log_path=log,checkpoint_path=log.with_suffix('.supervisor.json'))
    print(log.read_text(),end='')
    if result['outcome']!='completed':raise SystemExit(result['outcome'])

if __name__=='__main__':main()
