#!/usr/bin/env python3
"""Check finite premises of the written Q80 k0 pole proof; no analytic proof automation."""
import argparse
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import resource
import time

ROOT=Path(__file__).resolve().parents[2]
DEFAULT=ROOT/'artifacts/generated-results/elkies-k3-q80-genus-one-k0-pole-boundary-v1'
NOTE='elkies-k3/Q80_GENUS_ONE_K0_POLE_BOUNDARY_2026-09-14.md'
PARENT='artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
K0='artifacts/generated-results/elkies-k3-q80-genus-one-k0-quartic-gate-v1'
NODE='artifacts/generated-results/elkies-k3-q80-genus-one-k2-orientation-v1'
HELPER='elkies-k3/scripts/verify_q80_genus_one_k2_orientation.py'
HEIGHT='elkies-k3/COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md'
def read(p):return json.loads(p.read_text())
def digest(p):return sha256(p.read_bytes()).hexdigest()
def require(b,message):
    if not b:raise ValueError(message)
def bindings():
    paths={NOTE,PARENT,HELPER,HEIGHT,str(Path(__file__).relative_to(ROOT)),
          *[K0+'/'+n for n in ['input.json','result.json','independent-replay.json','quartic-input.json','quartic-producer.json','quartic-replay.json','quartic-comparison.json']],
          *[NODE+'/'+n for n in ['input.json','result.json','independent-replay.json','nodal-hypotheses.json']]}
    return {p:digest(ROOT/p) for p in sorted(paths)}
def freeze(out):
    out.mkdir(parents=True,exist_ok=True)
    with (out/'input.json').open('x') as f:
        json.dump({'schema':1,'prime':131,'cpu_seconds':20,'memory_bytes':2*1024**3,
                   'scope':'Source hashes, exact nodal hypotheses and the finite height/pole-degree table only. Analytic, Gauss-reduction and divisor-descent arguments require the written proof.',
                   'bindings':bindings()},f,indent=2,sort_keys=True);f.write('\n')
def check_prior_packet(base,input_name='input.json',result_name='result.json'):
    packet=read(base/input_name);result=read(base/result_name)
    require(result['status']=='PASS' and result['input_sha256']==digest(base/input_name),'prior result/input binding')
    for path,h in {**packet['bindings'],**packet.get('preserved_preflight',{})}.items():require(digest(ROOT/path)==h,'prior immutable source '+path)
    for path,h in result.get('records',{}).items():require(digest(base/path)==h,'prior exact record '+path)
    return result
def verify(out):
    started=time.process_time();packet=read(out/'input.json')
    require(packet['schema']==1 and packet['prime']==131 and packet['cpu_seconds']==20 and packet['memory_bytes']==2*1024**3,'finite check scope')
    require(packet['bindings']==bindings(),'written source and premise bytes')
    gate=check_prior_packet(ROOT/K0);prior=read(ROOT/K0/'independent-replay.json')
    require(prior['status']=='PASS' and prior['result_sha256']==digest(ROOT/K0/'result.json') and prior['input_sha256']==digest(ROOT/K0/'input.json'),'completed k0 independent replay binding')
    require(gate['necessary_D_mod131']==[15,39,90,124,1] and gate['at_least_one_negative_abscissa_required'],'precise inherited k0 premise')
    require(gate['remaining_nodal_denominator_boundary']=='UNKNOWN' and gate['genus1_k0_allocations_remaining']==5,'historical component assurance unchanged')
    # Check the completed quartic receipts and their checkpoint hashes, without repeating that census.
    qp=read(ROOT/K0/'quartic-input.json')
    for path,h in {**qp['bindings'],**qp['preserved_preflight']}.items():require(digest(ROOT/path)==h,'quartic frozen source '+path)
    totals=[]
    for kind in ['producer','replay']:
        q=read(ROOT/K0/('quartic-'+kind+'.json'));require(q['input_sha256']==digest(ROOT/K0/'quartic-input.json'),'quartic input binding')
        for path,h in q['checkpoints'].items():require(digest(ROOT/K0/path)==h,'quartic checkpoint content')
        totals.append((q['orbits'],q['smooth_split'],q['candidate_rows']))
    require(totals[0]==totals[1]==(73620690,12267623,[]),'complete retained quartic conclusion')
    spec=importlib.util.spec_from_file_location('q80_pole_nodal_premises',ROOT/HELPER);helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    nodal=helper.verify(ROOT/NODE)
    require(nodal['status']=='PASS' and nodal['written_analytic_proof_required'],'independently rechecked finite nodal premises')
    q0=[88,62,1];require(all((t*t+62*t+88)%131 for t in range(131)),'q0 has degree two and no rational residue root')
    table=[]
    for j in range(5):
        cross=4-2*j;minus=16-2*cross;plus=16+2*cross
        require((minus,plus)==(8+4*j,24-4*j),'polarized height identities')
        require((minus-8)%4==0 and (plus-8)%4==0,'integral descended pole degrees')
        sminus=(minus-8)//4;splus=(plus-8)//4
        table.append({'j':j,'cross_height':cross,'height_minus':minus,'height_plus':plus,
                      'pole_degree_minus':sminus,'pole_degree_plus':splus,'smaller_pole_degree':min(sminus,splus)})
    excluded=[r['j'] for r in table if r['smaller_pole_degree']<2]
    require(excluded==[0,1,3,4] and table[2]['pole_degree_minus']==table[2]['pole_degree_plus']==2,'exact four-case versus orthogonal boundary')
    return {'status':'PASS','input_sha256':digest(out/'input.json'),'checker_sha256':digest(Path(__file__)),
            'written_proof_sha256':digest(ROOT/NOTE),'scope':'Finite premises only; the exclusions and pole condition use the written proof.',
            'nodal_premise_replay':nodal,'height_pole_table':table,'excluded_k0_j_by_written_proof':excluded,
            'remaining_k0_j':[2],'remaining_pole_quadratic_reduction':q0,'remaining_boundary':'UNKNOWN',
            'new_analytic_proof_independently_replayed':False,'formal_verification':False,'external_review':False,
            'complete_quartic_or_contact_censuses_rerun':False,'positive_mw17_target_complete':False,
            'cpu_seconds':time.process_time()-started}
if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mode',choices=['freeze','check']);p.add_argument('--directory',type=Path,default=DEFAULT);p.add_argument('--receipt',type=Path);args=p.parse_args()
    resource.setrlimit(resource.RLIMIT_CPU,(20,25));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    if args.mode=='freeze':freeze(args.directory)
    else:
        record=verify(args.directory)
        if args.receipt:
            with args.receipt.open('x') as f:json.dump(record,f,indent=2,sort_keys=True);f.write('\n')
        print(json.dumps(record),flush=True)
