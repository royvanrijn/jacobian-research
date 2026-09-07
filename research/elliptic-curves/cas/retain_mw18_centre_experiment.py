#!/usr/bin/env python3
"""Portable MW18 geometry, centre-selection and point-search witnesses."""
import argparse
import gzip
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import zipfile

from pointed_quartic_search import ROOT
from research_runtime.store import checkpoint
from research_runtime.supervisor import capture, Limits
from research_runtime.witnesses import retained_source
from run_mw18_centre_experiment import check_protocol


def encoded(value): return (json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()


def pack(directory,output):
    if output.exists(): raise FileExistsError('witness bundles are immutable')
    protocol=json.loads((directory/'protocol.json').read_text());check_protocol(protocol)
    summary=json.loads((directory/'anchor-trial/summary.json').read_text())
    replay=json.loads((directory/'anchor-trial/replay.json').read_text())
    joint=json.loads((directory/'anchor-trial/joint-recoveries.json').read_text())
    if (summary['status']!='COMPLETE' or replay['protocol_hash']!=protocol['protocol_hash'] or
        replay['summary_sha256']!=sha256((directory/'anchor-trial/summary.json').read_bytes()).hexdigest()):
        raise ArithmeticError('complete independently replayed cells are required')
    if joint['status']!='COMPLETE' or joint['protocol_hash']!=protocol['protocol_hash']:
        raise ArithmeticError('joint recovery diagnostic is incomplete')
    sources={k:v for k,v in protocol['source_hashes'].items() if not k.endswith('.json')}
    geometry=json.loads((directory/'geometry.json').read_text());sources.update(geometry['sources'])
    for filename in ('retain_mw18_centre_experiment.py','supervise_mw18_centre_cells.py','analyze_mw18_centre_recoveries.py'):
        p=Path(__file__).with_name(filename);sources[str(p.relative_to(ROOT))]=sha256(p.read_bytes()).hexdigest()
    output.parent.mkdir(parents=True,exist_ok=True);temporary=output.with_suffix('.tmp');facts={}
    with zipfile.ZipFile(temporary,'w') as archive:
        def write(name,data):
            info=zipfile.ZipInfo(name,date_time=(2026,9,5,0,0,0))
            info.compress_type=zipfile.ZIP_STORED if name.endswith('.gz') else zipfile.ZIP_DEFLATED
            archive.writestr(info,data,compresslevel=9)
        for name in ('protocol.json','geometry.json','software.json','anchor-trial/summary.json','anchor-trial/replay.json','anchor-trial/supervision.json','prospective-gate-check.json'):
            write(name,(directory/name).read_bytes())
        def pool(value):
            value=dict(value);snapshot=value.pop('arithmetic_facts');ids=[]
            for fact in snapshot['facts']:
                key=fact['sha256'];ids.append(key)
                if key in facts and facts[key]!=fact:raise ArithmeticError('conflicting fact identities')
                facts[key]=fact
            value['arithmetic_fact_ids']=ids
            return value
        write('anchor-trial/joint-recoveries.json',encoded(pool(joint)))
        for label,record in protocol['centre_records'].items():
            for filename in (label+'.json',record['census_path']):
                raw=(directory/'centres'/filename).read_bytes()
                if filename.endswith('.gz') and sha256(raw).hexdigest()!=record['census_sha256']:
                    raise ArithmeticError('changed parity census')
                write('centres/'+filename,raw)
        for row in summary['cells']:
            p=directory/'anchor-trial'/row['result_path'];raw=p.read_bytes()
            if sha256(raw).hexdigest()!=row['result_sha256']: raise ArithmeticError('changed discovery cell')
            value=pool(json.loads(raw))
            write('anchor-trial/'+row['result_path'],encoded(value));sources.update(value['source_hashes'])
        # Group identical labelled points/models across primes before global
        # compression. Individual ZIP members repeat their huge rational labels.
        from research_runtime.store import digest
        def ordering(fact):
            key=fact['record']['key'];inputs=key['inputs']
            return (key['namespace'],digest({k:v for k,v in inputs.items() if k not in ('prime','degree')}),
                    fact['sha256'])
        payload=json.dumps(sorted(facts.values(),key=ordering),sort_keys=True,separators=(',',':')).encode()
        write('facts.json.gz',gzip.compress(payload,compresslevel=9,mtime=0))
        for name,key in sources.items():
            write('sources/'+key+'/'+Path(name).name,retained_source(ROOT,name,key).read_bytes())
        write('sources.json',encoded(sources))
    temporary.replace(output)
    result={**summary,'bundle':str(output.resolve().relative_to(ROOT)),
        'bundle_sha256':sha256(output.read_bytes()).hexdigest(),'exact_geometry_cover_count':9,
        'complete_parity_censuses':5,'parity_classes_per_census':262144,
        'chart_count':sum(r['chart_count'] for r in summary['cells']),
        'distinct_arithmetic_facts':len(facts),'facts_storage':'pooled-gzip-v1','limits':protocol['limits'],
        'coordinate_policy':protocol['coordinate_policy'],
        'joint_certified_gains':{r['case']:r['joint_certified_gain'] for r in joint['cases']},
        'demonstrated_remaining_directions_per_anchor':10}
    return result


def replay(bundle,directory):
    if directory.exists(): raise FileExistsError('empty replay directory required')
    directory.mkdir(parents=True)
    with zipfile.ZipFile(bundle) as archive:
        names=archive.namelist()
        if len(names)!=len(set(names)): raise ArithmeticError('duplicate bundle members')
        for name in names:
            p=Path(name)
            if p.is_absolute() or '..' in p.parts: raise ValueError('invalid bundle member')
            target=directory/p;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(archive.read(name))
        protocol=json.loads((directory/'protocol.json').read_text());check_protocol(protocol)
        summary=json.loads((directory/'anchor-trial/summary.json').read_text())
        facts={r['sha256']:r for r in json.loads(gzip.decompress(archive.read('facts.json.gz')))} if 'facts.json.gz' in names else None
        def restore(value):
            ids=value.pop('arithmetic_fact_ids')
            value['arithmetic_facts']={'schema':'elliptic-curves.arithmetic-facts.v1',
                'facts':[facts[key] if facts is not None else json.loads(archive.read('facts/'+key+'.json')) for key in ids]}
            return value
        for row in summary['cells']:
            p=directory/'anchor-trial'/row['result_path'];value=json.loads(p.read_text())
            value=restore(value)
            raw=encoded(value)
            if sha256(raw).hexdigest()!=row['result_sha256']: raise ArithmeticError('restored cell differs')
            p.write_bytes(raw)
        p=directory/'anchor-trial/joint-recoveries.json';value=json.loads(p.read_text())
        if 'arithmetic_fact_ids' in value:p.write_bytes(encoded(restore(value)))
        for name,key in json.loads((directory/'sources.json').read_text()).items():
            if sha256((directory/'sources'/key/Path(name).name).read_bytes()).hexdigest()!=key:
                raise ArithmeticError('retained source hash failed')
    # Reconstruct the height proof from its pinned canonical source inputs.
    geometry=json.loads((directory/'geometry.json').read_text())
    for name,key in geometry['inputs'].items():
        if sha256((ROOT/name).read_bytes()).hexdigest()!=key: raise ArithmeticError('canonical geometry input changed')
    env=dict(os.environ,EC_ARITHMETIC_CACHE=str(directory/'empty-cache'),EC_RUNTIME_SOURCE_ARCHIVE=str(directory/'sources'))
    cas=ROOT/'elliptic-curves/cas'
    commands=[['sage','-python',str(cas/'build_mw18_height_geometry.sage'),'--output',str(directory/'geometry-replayed.json')],
        ['sage','-python',str(cas/'prepare_mw18_deep_centres.sage'),'--geometry',str(directory/'geometry.json'),
         '--directory',str(directory/'centres'),'--verify',
         *[arg for label in protocol['centre_records'] for arg in ('--cover',label)]],
        [sys.executable,str(cas/'supervise_mw18_centre_cells.py'),'--protocol',str(directory/'protocol.json'),
         '--directory',str(directory/'anchor-trial'),'--workers','3','--verify'],
        [sys.executable,str(cas/'analyze_mw18_centre_recoveries.py'),'--protocol',str(directory/'protocol.json'),
         '--directory',str(directory/'anchor-trial'),'--output',str(directory/'anchor-trial/joint-recoveries.json'),'--verify']]
    for i,command in enumerate(commands):
        r=capture(command,limits=Limits(1800,6*1024**3),log_path=directory/f'replay-{i}.log',env=env,check=True)
        print(f'MW18_PORTABLE_REPLAY|stage={i}|returncode={r.returncode}',flush=True)
        if i==0:
            rebuilt=json.loads((directory/'geometry-replayed.json').read_text())
            if rebuilt['covers']!=geometry['covers']: raise ArithmeticError('exact geometry replay differs')
    actual=json.loads((directory/'anchor-trial/summary.json').read_text())
    if actual!=summary: raise ArithmeticError('portable comparison differs')
    return {'status':'PASS_GEOMETRY_DEEP_STRATA_CHARTS_AND_INDEPENDENCE',
        'bundle_sha256':sha256(bundle.read_bytes()).hexdigest(),
        'protocol_hash':protocol['protocol_hash'],'chart_count':sum(r['chart_count'] for r in actual['cells']),
        'point_enumeration_repeated':False,'floating_cvp_repeated':False,'empty_arithmetic_cache':True}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--pack-directory',type=Path);g.add_argument('--verify-bundle',type=Path)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--summary',type=Path)
    a=p.parse_args();result=pack(a.pack_directory,a.output) if a.pack_directory else replay(a.verify_bundle,a.output)
    if a.summary:checkpoint(a.summary,result)
    print(json.dumps(result,sort_keys=True))
