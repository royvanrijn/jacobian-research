#!/usr/bin/env python3
"""Pack portable chart witnesses, or replay them without CVP/point enumeration."""
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
import zipfile

from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits,captured_run
from research_runtime.witnesses import retained_source
from pointed_quartic_search import ROOT


def encoded(value):return (json.dumps(value,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
def write(archive,name,data):
    info=zipfile.ZipInfo(name,date_time=(2026,9,5,0,0,0));info.compress_type=zipfile.ZIP_DEFLATED
    archive.writestr(info,data,compresslevel=9)


def pack(directory,output):
    if output.exists():raise FileExistsError('retain immutable witness bundles')
    result=json.loads((directory/'result.json').read_text())
    verified=json.loads((directory/'verified.json').read_text())
    if result['status']!='COMPLETE' or verified['status']!='COMPLETE' or result['cells']!=verified['cells']:
        raise ArithmeticError('the full sweep must replay before publication')
    sources={};seen=set();output.parent.mkdir(parents=True,exist_ok=True)
    temporary=output.with_suffix('.tmp')
    with zipfile.ZipFile(temporary,'w') as archive:
        write(archive,'protocol.json',(directory/'protocol.json').read_bytes())
        write(archive,'result.json',encoded(result))
        for row in result['cells']:
            path=directory/row['result'];raw=path.read_bytes()
            if sha256(raw).hexdigest()!=row['result_sha256']:raise ArithmeticError('changed discovery cell')
            value=json.loads(raw);snapshot=value.pop('arithmetic_facts')
            ids=[]
            for fact in snapshot['facts']:
                key=fact['sha256'];ids.append(key)
                if key not in seen:
                    write(archive,'facts/'+key+'.json',encoded(fact));seen.add(key)
            value['arithmetic_fact_ids']=ids
            write(archive,row['result'],encoded(value))
            sources.update(value['source_hashes'])
        for name in ('run_mw_search.py','calibrate_chart_policy.py','research_runtime/cvp.py','research_runtime/chart_policy.py'):
            path=ROOT/'elliptic-curves/cas'/name;sources[str(path.relative_to(ROOT))]=sha256(path.read_bytes()).hexdigest()
        for name,key in sources.items():
            write(archive,'sources/'+key+'/'+Path(name).name,retained_source(ROOT,name,key).read_bytes())
        write(archive,'sources.json',encoded(sources))
    temporary.replace(output)
    return summarize(output)


def summarize(output):
    with zipfile.ZipFile(output) as archive:
        result=json.loads(archive.read('result.json'));frozen=json.loads(archive.read('protocol.json'))
        facts=sum(name.startswith('facts/') for name in archive.namelist())
    absolute=output.resolve()
    return {'schema':'elliptic-curves.portable-chart-policy-sweep.v1','bundle':str(absolute.relative_to(ROOT)) if absolute.is_relative_to(ROOT) else str(absolute),
            'bundle_sha256':sha256(output.read_bytes()).hexdigest(),'protocol_hash':result['protocol_hash'],
            'chart_count':len(result['cells'])*frozen['protocol']['limits']['next_holes'],'cell_count':len(result['cells']),
            'ranking':result['ranking'],'distinct_arithmetic_facts':facts,
            'calibration_fibre_count':len(frozen['protocol']['panel']),'held_out_fibre_count':len(frozen['protocol']['held_out_controls']),
            'rank_groups':{rank:len(ids) for rank,ids in frozen['panel_roles'].items()},
            'claim_boundary':'Historical control calibration with retained exact positive point witnesses. Timing differences and equal recovery counts do not establish a unique optimal metric.'}


def replay(bundle,directory):
    if directory.exists():raise FileExistsError('choose a fresh replay directory')
    directory.mkdir(parents=True)
    with zipfile.ZipFile(bundle) as archive:
        protocol=json.loads(archive.read('protocol.json'));result=json.loads(archive.read('result.json'))
        sources=json.loads(archive.read('sources.json'))
        for name,key in sources.items():
            if len(key)!=64 or any(c not in '0123456789abcdef' for c in key):raise ValueError('invalid source hash')
            data=archive.read('sources/'+key+'/'+Path(name).name)
            if sha256(data).hexdigest()!=key:raise ArithmeticError('corrupt retained source')
            target=directory/'sources'/key/Path(name).name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data)
        checkpoint(directory/'protocol.json',protocol);checkpoint(directory/'result.json',result)
        for row in result['cells']:
            name=Path(row['result'])
            if name.is_absolute() or name.parts[0]!='cells' or '..' in name.parts:raise ValueError('invalid cell member')
            value=json.loads(archive.read(str(name)))
            ids=value.pop('arithmetic_fact_ids')
            value['arithmetic_facts']={'schema':'elliptic-curves.arithmetic-facts.v1',
                'facts':[json.loads(archive.read('facts/'+key+'.json')) for key in ids]}
            raw=encoded(value)
            if sha256(raw).hexdigest()!=row['result_sha256']:raise ArithmeticError('restored cell differs from discovery witness')
            path=directory/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw)
    env=dict(os.environ,EC_ARITHMETIC_CACHE=str(directory/'empty-arithmetic-cache'),
             EC_RUNTIME_SOURCE_ARCHIVE=str(directory/'sources'))
    captured_run([sys.executable,str(ROOT/'elliptic-curves/cas/calibrate_chart_policy.py'),
        '--protocol',str(directory/'protocol.json'),'--verify','--output',str(directory/'verified.json')],
        limits=Limits(1200,2_147_483_648),env=env,text=True,check=True)
    verified=json.loads((directory/'verified.json').read_text())
    if verified['cells']!=result['cells']:raise ArithmeticError('portable sweep replay differs')
    return verified['status']


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--pack-directory',type=Path);g.add_argument('--verify-bundle',type=Path);g.add_argument('--summarize-bundle',type=Path)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--summary',type=Path)
    a=p.parse_args()
    if a.summarize_bundle:
        checkpoint(a.output,summarize(a.summarize_bundle));return
    if a.pack_directory:
        result=pack(a.pack_directory,a.output)
        if a.summary:checkpoint(a.summary,result)
        print(json.dumps(result,sort_keys=True))
    else:print('PORTABLE_POLICY_REPLAY|'+replay(a.verify_bundle,a.output))

if __name__=='__main__':main()
