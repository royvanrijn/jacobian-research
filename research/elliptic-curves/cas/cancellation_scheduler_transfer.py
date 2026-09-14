#!/usr/bin/env sage -python
"""Gated M27/M29/M30 transfer test; retained controls, never new records."""
import argparse
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import time

import cancellation_scheduler_round4 as round4
from cancellation_scheduler_prepare import new_write
from finite_cancellation_corpus import ROOT, LOCAL, canonical, digest, write

PRIOR=round4.OUT
OUT=ROOT/'artifacts/generated-results/elliptic-curves/cancellation_scheduler_transfer_v1'
RAW=LOCAL/'cancellation-scheduler-transfer-v1'
CAS=Path(__file__).resolve().parent


def prepare():
    from memory_rank_certificate import checked_rank
    from finite_cancellation_validation_audit import j_group
    prior=json.loads((PRIOR/'protocol.json').read_text())
    design={'status':'DESIGN_BEFORE_HIGH_RANK_CONTROL_SEARCH','cases':3,
        'arms':prior['arms'],'candidate':prior['candidate'],'heights':[8000,32000,125000],
        'maximum_centres':256,'maximum_calls':2304,'search_cpu_seconds':150,
        'point_wall_seconds':5,'arm_wall_seconds':240,'hard_process_cpu_seconds':220,
        'selection':'All three existing height-model benchmark inputs: Curve302 M27, M29 and M30, with their complete retained 256-centre order. No winning centre or target location selects inputs or scheduling. These are one known curve, excluded from all corpus fitting; no independent-population bootstrap.',
        'gate':'Run only after V4 independent withheld gate passes. Fresh-fibre gate then requires all three independently certified control recoveries under both arms, at least 10 percent lower aggregate complete CPU for the candidate, and no candidate case exceeding 1.5 times baseline CPU. No gain or curve is a new record.',
        'fit':prior['fit'],'gp_sha256':prior['gp_sha256'],
        'cost':'Complete isolated process CPU through independent Sage certification, including all maps and discarded work. The fixed retained centre bank is a common input; its historical construction is not a cold search measurement.'}
    new_write(OUT/'design.json',design)
    old=ROOT/'artifacts/generated-results/elliptic-curves/height_model_next_direction_v1'
    source_plan=json.loads((old/'protocol.json').read_text());inputs=[];bindings={};preflight=[];start=time.process_time()
    checker=SourceFileLoader('transfer_independent_rank',str(CAS/'verify_finite_cancellation_cpu.sage')).load_module().sage_rank
    for case in source_plan['cases']:
        path=LOCAL/'height-model-next-direction-v1'/case['input']
        if digest(path.read_bytes())!=case['sha256']:raise ArithmeticError('retained high-rank input changed')
        envelope=json.loads(path.read_text());source=envelope['payload']
        seed={'curve':source['curve'],'points':source['points']}
        seed['proof']=checked_rank(tuple(map(F,seed['curve'])),[tuple(map(F,p)) for p in seed['points']],source['primes'],source['torsion_prime'])
        result=checker(seed)
        inputs.append({'id':source['id'],'j_group':j_group(source['curve']),
            'family':'Curve302-known-control','stratum':'high_rank_transfer','seed':seed,
            'centres':[c['representative'] for c in source['centres']]})
        if len(inputs[-1]['centres'])!=256:raise ArithmeticError('retained centre count differs')
        bindings[str(path.relative_to(ROOT))]=digest(path.read_bytes());preflight.append(result)
    new_write(OUT/'inputs.json',inputs)
    new_write(OUT/'preflight.json',{'status':'PASS','rows':preflight,'source_inputs_sha256':bindings,
        'cpu_seconds':time.process_time()-start,'inputs_sha256':digest((OUT/'inputs.json').read_bytes())})
    print(json.dumps({'status':'PREPARED_NOT_SEARCHED','ranks':[len(c['seed']['points']) for c in inputs]}),flush=True)


def seal():
    summary=json.loads((PRIOR/'summary.json').read_text())
    if summary['high_rank_transfer_gate'] is not True:raise ArithmeticError('withheld validation did not pass; no transfer search')
    prior=json.loads((PRIOR/'protocol.json').read_text());plan=json.loads((OUT/'design.json').read_text())
    if json.loads((OUT/'preflight.json').read_text())['status']!='PASS':raise ArithmeticError('preflight incomplete')
    sources=dict(prior['source_sha256']);sources[str(Path(__file__).relative_to(ROOT))]=digest(Path(__file__).read_bytes())
    plan.update(status='FROZEN_HIGH_RANK_TRANSFER',source_sha256=sources,
        prerequisite_summary_sha256=digest((PRIOR/'summary.json').read_bytes()),
        input_sha256={n:digest((OUT/n).read_bytes()) for n in ('design.json','inputs.json','preflight.json')})
    new_write(OUT/'protocol.json',plan);print('SEALED_HIGH_RANK_TRANSFER',flush=True)


def configure():
    round4.OUT=OUT;round4.RAW=RAW;round4.ENTRY=Path(__file__).resolve()


def run():
    configure();round4.run()


def worker(case,arm):
    configure();round4.worker(case,arm)


def report():
    plan=json.loads((OUT/'protocol.json').read_text());supervision=json.loads((OUT/'supervision.json').read_text())
    if supervision['status']!='COMPLETE' or len(supervision['records'])!=6:raise ArithmeticError('incomplete transfer')
    rows=[];table={};candidate=plan['candidate']
    for receipt in supervision['records']:
        path=RAW/'arms'/receipt['case']/receipt['arm'];result=json.loads((path/'result.json').read_text())
        verification=json.loads((path/'independent-verification.json').read_text())
        if receipt['status']!='COMPLETE' or receipt['result_sha256']!=digest((path/'result.json').read_bytes()):raise ArithmeticError('bad receipt')
        if verification['status']!='PASS' or result['independent_verification_sha256']!=digest((path/'independent-verification.json').read_bytes()):raise ArithmeticError('bad independent proof')
        row={'case':receipt['case'],'arm':receipt['arm'],'success':result['success'],
             'rank_lower_bound':result['rank_lower_bound'],'calls':len(result['calls']),
             'cpu_seconds':receipt['charged_cpu_seconds'],'status':result['status'],
             'preparation_unknowns':verification['preparation_unknowns']}
        rows.append(row);table[row['case'],row['arm']]=row
    totals={arm:sum(r['cpu_seconds'] for r in rows if r['arm']==arm) for arm in plan['arms']}
    cases=sorted({r['case'] for r in rows})
    ratios={c:table[c,candidate]['cpu_seconds']/table[c,'factor_free']['cpu_seconds'] for c in cases}
    gate=all(r['success'] and not r['preparation_unknowns'] for r in rows) and totals[candidate]<=.9*totals['factor_free'] and max(ratios.values())<=1.5
    write(OUT/'summary.json',{'status':'COMPLETE_KNOWN_HIGH_RANK_TRANSFER','rows':rows,'cpu_totals':totals,
        'cpu_ratios':ratios,'fresh_fibre_gate':gate,'rank32':'UNKNOWN',
        'protocol_sha256':digest((OUT/'protocol.json').read_bytes()),
        'boundary':'Three subgroups on one known curve. Exact lower-bound recoveries only; no new rank, no population inference, no automatic fresh-fibre input construction or larger allowance.'})
    print(json.dumps({'cpu_totals':totals,'cpu_ratios':ratios,'fresh_fibre_gate':gate,'rows':rows},indent=2),flush=True)


def pack():
    configure();round4.pack()


def replay():
    configure();round4.replay()


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['prepare','seal','run','worker','report','pack','replay'])
    parser.add_argument('--case');parser.add_argument('--arm');args=parser.parse_args()
    worker(args.case,args.arm) if args.command=='worker' else globals()[args.command]()
