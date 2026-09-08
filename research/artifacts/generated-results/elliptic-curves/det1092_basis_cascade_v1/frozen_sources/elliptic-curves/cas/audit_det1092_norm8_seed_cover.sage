#!/usr/bin/env sage-python
"""Checkpoint one pre-existing generic norm8 construction; no point search.

Stage zero only. Preserve failure and inputs; 25-second cap. The original
producer and any existing artifacts are never overwritten.
"""
import hashlib,json,runpy,signal,traceback
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
CAS=ROOT/'elliptic-curves/cas'
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_norm8_seed_cover_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def main():
    paths=[CAS/'construct_det1092_norm8_elliptic_pencil.sage',
           CAS/'load_curve302_recovered_parent.sage',CAS/'visibility_lattice_v2.py',
           ART/'curve302_recovered_mw17_parent_v1.json',
           ART/'curve302_parent_degree2_multisection_orbits_v1.tsv',Path(__file__)]
    protocol={'classification':'bounded audit of pre-existing generic construction',
       'rule':'Execute build() of the unchanged generic norm8 pencil script once. Its fixed choice is minimum(l1,linfinity,orbit) from the old norm8 census, its auxiliary address is z=0, and its finite certificate pool is unchanged. Do not select a replacement on failure. No exceptional point is an input to this stage.',
       'limits':{'wall_seconds':25,'pencils':1,'parity_node_limit':200000,
                 'fixed_auxiliary_addresses':1,'point_searches':0,
                 'exceptional_point_inputs':0,'V3_inputs':0,'pilot_changes':0},
       'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    retain(DIR/'protocol.json',protocol)
    def expired(signum,frame):raise TimeoutError('FROZEN_25_SECOND_LIMIT')
    signal.signal(signal.SIGALRM,expired);signal.alarm(25)
    try:
        result=runpy.run_path(str(paths[0]))['build']()
        retain(DIR/'generic.json',result)
        print(result['status'],'orbit',result['selection']['orbit'],flush=True)
    except BaseException as exc:
        retain(DIR/'failure.json',{'status':'UNRESOLVED_FIXED_CONSTRUCTION',
          'error':str(exc),'traceback':traceback.format_exc(),
          'inputs':protocol['inputs']})
        raise
    finally:signal.alarm(0)

if __name__=='__main__':
    DIR.mkdir(parents=True,exist_ok=True);main()
