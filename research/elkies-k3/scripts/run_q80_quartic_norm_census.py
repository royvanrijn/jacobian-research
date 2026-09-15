#!/usr/bin/env python3
"""Freeze and execute the exact degree-four residue gate in ten resumable ranges."""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import resource
import subprocess
import tempfile
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-k0-quartic-gate-v1'
SOURCE='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
PRODUCER='elkies-k3/scripts/certify_q80_quartic_norm_census.cpp'
REPLAY='elkies-k3/scripts/replay_q80_quartic_norm_census.cpp'
INPUT_HELPER='elkies-k3/scripts/verify_q80_genus_one_k1_reciprocity.py'
RANGES=[[i,min(i+429,4290)] for i in range(0,4290,429)]
def digest(p):return sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def require(b,message):
    if not b:raise ValueError(message)
def write(p,obj):
    with p.open('x') as f:json.dump(obj,f,indent=2,sort_keys=True);f.write('\n')
def input_text():
    spec=importlib.util.spec_from_file_location('q80_k1_input',ROOT/INPUT_HELPER)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    # The inherited helper has its own frozen identity; only rational polynomial arithmetic is reused.
    source=read(ROOT/SOURCE);model=source['weierstrass_model'];lines=[]
    def emit(a):lines.append(str(len(a))+' '+' '.join(map(str,a)))
    emit(m.coeff(model['A_coefficients_low_to_high']));emit(m.coeff(model['B_coefficients_low_to_high']));lines.append('17')
    require(len(source['sections']['records'])==17,'literal basis size')
    for r in source['sections']['records']:
        x=r['X'];n=m.coeff(x['numerator_coefficients_low_to_high']);d=m.coeff(x['denominator_coefficients_low_to_high']);g=m.gcd(n,d)
        n,rn=m.divrem(n,g);d,rd=m.divrem(d,g);require(not rn and not rd,'exact reduced coordinates');emit(n);emit(d)
    return '\n'.join(lines)+'\n'
def freeze(out):
    paths=[SOURCE,PRODUCER,REPLAY,INPUT_HELPER,str(Path(__file__).relative_to(ROOT))]
    packet={'schema':1,'scope':'All genuine degree-four base orbits at p=131 with fully split cubic and all three inherited-root norm characters zero; finite necessary gate only.',
            'prime':131,'orbits':73620690,'tower':'i^2=-1, j^2=1+i','monomial_modulus':[2,0,-2,0,1],
            'pair_ranges':RANGES,'orbits_per_pair':17161,'cpu_seconds_per_range':90,'memory_bytes_per_range':512*1024**2,
            'driver_cpu_seconds':30,'driver_memory_bytes':2*1024**3,'process_wall_timeout':95,
            'bindings':{p:digest(ROOT/p) for p in paths},
            'preserved_preflight':{str(p.relative_to(ROOT)):digest(p) for p in sorted((out/'preflight').glob('*')) if p.is_file()}}
    text=input_text();packet['coefficient_input_sha256']=sha256(text.encode()).hexdigest()
    with (out/'quartic-input.txt').open('x') as f:f.write(text)
    write(out/'quartic-input.json',packet)
def validate_packet(out):
    p=read(out/'quartic-input.json')
    require(p['pair_ranges']==RANGES and p['orbits']==73620690 and p['orbits_per_pair']==17161,'complete partition')
    require(p['prime']==131 and p['monomial_modulus']==[2,0,-2,0,1],'field scope')
    for path,h in {**p['bindings'],**p['preserved_preflight']}.items():require(digest(ROOT/path)==h,'frozen binding '+path)
    text=input_text();require(text==(out/'quartic-input.txt').read_text(),'literal coefficient input')
    require(sha256(text.encode()).hexdigest()==p['coefficient_input_sha256'],'coefficient hash')
    return p,text
def run(out,kind):
    packet,text=validate_packet(out);source=PRODUCER if kind=='producer' else REPLAY
    directory=out/('quartic-'+kind);directory.mkdir(exist_ok=True);started=time.monotonic()
    with tempfile.TemporaryDirectory(prefix='q80-quartic-'+kind+'-') as tmp:
        exe=Path(tmp)/'census';subprocess.run(['g++','-O3','-std=c++17',str(ROOT/source),'-o',str(exe)],check=True,timeout=30)
        for begin,end in RANGES:
            path=directory/f'{begin:04d}-{end:04d}.json'
            if path.exists():
                previous=read(path);require(previous['source_sha256']==digest(ROOT/source) and previous['input_sha256']==digest(out/'quartic-input.json'),'checkpoint binding')
                continue
            step=time.monotonic()
            try:
                completed=subprocess.run([str(exe),str(begin),str(end)],input=text,text=True,capture_output=True,check=True,timeout=95)
                result=json.loads(completed.stdout)
                require(result['status']=='COMPLETE_RANGE' and result['pair_begin']==begin and result['pair_end']==end and result['orbits']==(end-begin)*17161,'complete range')
                require(sum(result['first_failed_character'])==result['smooth_split'],'exhaustive split histogram')
            except Exception as exc:
                write(directory/f'failure-{begin:04d}-{time.time_ns()}.json',{'status':'FAILED','range':[begin,end],'error':repr(exc)})
                raise
            write(path,{'source_sha256':digest(ROOT/source),'input_sha256':digest(out/'quartic-input.json'),'result':result,'wall_seconds':time.monotonic()-step})
            print(json.dumps({'kind':kind,'completed_pair_range':[begin,end],'orbits':result['orbits'],'candidates':len(result['candidate_rows']),'cpu_seconds':result['cpu_seconds']}),flush=True)
    rows=[read(directory/f'{a:04d}-{b:04d}.json')['result'] for a,b in RANGES]
    aggregate={'status':'COMPLETE','input_sha256':digest(out/'quartic-input.json'),'source_sha256':digest(ROOT/source),
               'orbits':sum(r['orbits'] for r in rows),'smooth_split':sum(r['smooth_split'] for r in rows),
               'first_failed_character':[sum(r['first_failed_character'][i] for r in rows) for i in range(18)],
               'candidate_rows':[c for r in rows for c in r['candidate_rows']],
               'cpu_seconds':sum(r['cpu_seconds'] for r in rows),'session_wall_seconds':time.monotonic()-started,
               'checkpoints':{str(path.relative_to(out)):digest(path) for path in sorted(directory.glob('[0-9]*.json'))}}
    require(aggregate['orbits']==73620690,'full exact orbit count')
    target=out/('quartic-'+kind+'.json')
    if not target.exists():write(target,aggregate)
    else:
        previous=read(target)
        for key in ['input_sha256','source_sha256','orbits','smooth_split','first_failed_character','candidate_rows','checkpoints']:require(previous[key]==aggregate[key],'preserved aggregate '+key)
    print(json.dumps(aggregate),flush=True)
def compare(out):
    validate_packet(out);a=read(out/'quartic-producer.json');b=read(out/'quartic-replay.json')
    for lo,hi in RANGES:
        name=f'{lo:04d}-{hi:04d}.json';ar=read(out/'quartic-producer'/name);br=read(out/'quartic-replay'/name)
        for obj,source in [(ar,PRODUCER),(br,REPLAY)]:
            require(obj['input_sha256']==digest(out/'quartic-input.json') and obj['source_sha256']==digest(ROOT/source),'range identities')
        for key in ['status','pair_begin','pair_end','orbits','smooth_split','first_failed_character','split_character_checksum','candidate_rows']:require(ar['result'][key]==br['result'][key],'independent exact range comparison '+key)
    for obj,source in [(a,PRODUCER),(b,REPLAY)]:
        require(obj['source_sha256']==digest(ROOT/source) and obj['input_sha256']==digest(out/'quartic-input.json'),'aggregate identities')
        for path,h in obj['checkpoints'].items():require(digest(out/path)==h,'checkpoint content binding')
    for key in ['status','orbits','smooth_split','first_failed_character','candidate_rows']:require(a[key]==b[key],'full independent comparison '+key)
    require(a['orbits']==73620690 and not a['candidate_rows'],'empty full quartic gate')
    receipt={'status':'PASS','input_sha256':digest(out/'quartic-input.json'),'producer_sha256':digest(out/'quartic-producer.json'),'replay_sha256':digest(out/'quartic-replay.json'),
             'orbits':a['orbits'],'smooth_split':a['smooth_split'],'norm_candidates':0,'producer_cpu_seconds':a['cpu_seconds'],'replay_cpu_seconds':b['cpu_seconds'],
             'full_characters_tested_only_until_first_failure':True,'nodal_denominator_boundary_excluded':False,'positive_mw17_target_complete':False}
    target=out/'quartic-comparison.json'
    if not target.exists():write(target,receipt)
    else:require(read(target)==receipt,'preserved comparison receipt')
    print(json.dumps(receipt),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['freeze','producer','replay','compare']);p.add_argument('--directory',type=Path,default=DEFAULT);args=p.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(30,95));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    if args.mode=='freeze':freeze(args.directory)
    elif args.mode=='compare':compare(args.directory)
    else:run(args.directory,args.mode)
