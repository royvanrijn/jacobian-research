#!/usr/bin/env python3
"""One fixed unconditional online descent; immutable transcript and fail-closed parser."""
import argparse
from pathlib import Path
import re
import subprocess
import sys
import retrospective as r

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'GENUS_TWO_RANK_AND_SIMPLICITY_PROTOCOL.json'
SOURCE=r.OUT/'rank_jump_native_genus_two_lift_gate_v1.json'
INPUT=r.OUT/'rank_jump_genus_two_rank_probe_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_genus_two_rank_probe_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-genus-two-rank-probe-v1'
ONLINE=r.ROOT/'elliptic-curves/cas/run_fixed_field_radical_covers.py'


def program():
    d=r.read(SOURCE);coef=','.join(d['sextic_coefficients'])
    return f'''SetColumns(0); SetSeed(1); SetNthreads(1);
print "RJG2_VERSION",GetVersion();
Q:=Rationals(); R<t>:=PolynomialRing(Q);
f:=R![{coef}]; C:=HyperellipticCurve(f); J:=Jacobian(C);
assert Genus(C) eq 2;
assert Evaluate(f,{d['retained_genus_two_point']['t']}) eq ({d['retained_genus_two_point']['y']})^2;
T:=TwoTorsionSubgroup(J); assert #T eq 4;
print "RJG2_TORSION",#T;
SetVerbose("Selmer",1); SetVerbose("ClassGroup",1);
print "RJG2_DESCENT_BEGIN";
S:=TwoSelmerGroup(J);
print "RJG2_SELMER_ORDER",#S;
print "RJG2_DONE";
'''


def prepare():
    WORK.mkdir(parents=True,exist_ok=True)
    p=program();(WORK/'probe.m').write_text(p)
    r.write_new(INPUT,{'schema':'rank-jump.genus-two-rank-probe-inputs.v1','program':p,
        'bindings':{str(x.relative_to(r.ROOT)):r.digest(x.read_bytes()) for x in (SOURCE,PROTOCOL,ONLINE)}})


def worker():
    inp=r.read(INPUT)
    for path,sha in inp['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    assert (WORK/'probe.m').read_text()==inp['program']==program()
    sys.path.insert(0,str(ONLINE.parent))
    import run_fixed_field_radical_covers as existing
    existing.run_online(WORK/'probe.m')


def analyze(raw):
    import xml.etree.ElementTree as ET
    if not raw:return {'status':'UNKNOWN','reason':'No completed server response'}
    root=ET.fromstring(raw);lines='\n'.join(''.join(x.itertext()) for x in root.findall('.//results/line'))
    warnings=[''.join(x.itertext()) for x in root.findall('.//warning')]
    errors=[s for s in ('Runtime error','User error','System Error','Syntax error','memory limit','Internal error') if s in lines]
    m=re.search(r'RJG2_SELMER_ORDER\s+(\d+)',lines)
    if warnings or errors or not m or 'RJG2_DONE' not in lines:
        return {'status':'UNKNOWN','warnings':warnings,'errors':errors,'transcript':lines}
    n=int(m[1]);assert n>=4 and n&(n-1)==0
    upper=n.bit_length()-1-2
    return {'status':'COMPLETE','Selmer_order':n,'rational_2_torsion_dimension':2,'rank_upper_bound':upper,
        'ordinary_Chabauty_gate':'OPEN' if upper<=1 else 'UNKNOWN','transcript':lines}


def collect():
    inp=r.read(INPUT);assert inp['program']==program()
    raw=(WORK/'probe.xml').read_text() if (WORK/'probe.xml').exists() else None
    return {'schema':'rank-jump.genus-two-rank-probe.v1','analysis':analyze(raw),'raw_xml':raw,
        'execution':r.read(WORK/'execution.json'),'worker_log':(WORK/'worker.log').read_text(),
        'transport_error':(WORK/'probe.error').read_text() if (WORK/'probe.error').exists() else None,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (INPUT,PROTOCOL,SOURCE,Path(__file__),ONLINE,HERE/'retrospective.py')}}


def run():
    log=WORK/'worker.log'
    if not log.exists():
        with log.open('x') as out:
            try:
                p=subprocess.run(['/home/royvanrijn/.local/bin/sage','-python',str(Path(__file__).resolve()),'worker'],stdout=out,stderr=out,timeout=120)
                status={'status':'COMPLETE' if p.returncode==0 else 'FAILED','returncode':p.returncode}
            except subprocess.TimeoutExpired:status={'status':'TIMEOUT'}
        r.write_new(WORK/'execution.json',status)
    d=collect();r.write_new(OUTPUT,d);print(d['analysis'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','worker','run','check']);a=p.parse_args()
    if a.mode=='prepare':prepare()
    elif a.mode=='worker':worker()
    elif a.mode=='run':run()
    else:
        d=r.read(OUTPUT)
        for path,sha in d['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
        assert d['analysis']==analyze(d['raw_xml']);assert r.read(INPUT)['program']==program()
        print('PASS offline transcript replay:',d['analysis']['status'])
