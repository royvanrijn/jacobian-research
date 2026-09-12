#!/usr/bin/env python3
"""Project retained pre-acquisition inputs; no new point search or enumeration."""
import argparse
import json
from pathlib import Path
import shutil
import time

from next_direction_benchmark import POLICIES, read, sha
from research_runtime.store import checkpoint, digest
from pointed_box_equivalence import box_key
from fractions import Fraction as F
from memory_rank_certificate import checked_rank

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
LOCAL = ROOT/'artifacts/local/elliptic-curves'


def prepare(output):
    started = time.monotonic()
    started_cpu = time.process_time()
    output.mkdir(exist_ok=False)
    bindings, cases, ladder = {}, [], []
    def bound(path):
        data = read(path)
        bindings[str(path.relative_to(ROOT))] = sha(path)
        return data
    def projection(identifier, role, curve, points, primes, torsion, centres, source, excluded=()):
        packet = {'id':identifier, 'curve':curve, 'points':points, 'primes':primes,
            'torsion_prime':torsion, 'centres':[{'representative':c['representative']} for c in centres],
            'completed_boxes':sorted(excluded)}
        if any(len(c['representative']) != len(points) for c in packet['centres']):
            raise ArithmeticError('centre uses a different starting basis')
        path = output/'inputs'/(identifier+'.json')
        checkpoint(path, {'sha256':digest(packet), 'payload':packet})
        cases.append({'id':identifier, 'role':role, 'initial_rank':len(points),
            'input':'inputs/'+path.name, 'sha256':sha(path), 'centres':len(centres), 'source':source})

    base = LOCAL/'adaptive-visibility-cascade-v3/replay-M17'
    stages = bound(base/'stages.json')
    # Preserve the complete ladder's pre-acquisition banks. Execute only the declared late rungs.
    for stage in stages:
        ep = base/f"epoch-{stage['epoch']:02d}"
        selection = bound(ep/'selection.json')
        cloud = bound(ep/'cloud-000.json')
        state = cloud['final_state']['state']
        if state['reductions']['points'] != selection['basis']:
            raise ArithmeticError('historical state is not the pre-acquisition basis')
        r = stage['before']
        packet = {'curve':cloud['curve'], 'points':selection['basis'],
            'primes':state['reductions']['primes'], 'torsion_prime':state['no_two_torsion_prime'],
            'centres':[{'representative':c['representative']} for c in selection['centres']]}
        path = output/'ladder'/f'curve302-M{r}.json'
        checkpoint(path,packet)
        ladder.append({'initial_rank':r,'path':str(path.relative_to(output)), 'sha256':sha(path),
            'centres':len(packet['centres'])})
        if r in (27,29,30):
            projection(f'curve302-M{r}','development',packet['curve'],packet['points'],packet['primes'],
                packet['torsion_prime'],packet['centres'][:256],str(ep.relative_to(ROOT)))

    path = LOCAL/'blind-factor-free-28-control-v1'
    seed, bank = bound(path/'seed.json'), bound(path/'maps.json')
    proof = seed['rank_certificate']
    projection('public188-M27','development',seed['curve'],seed['points'],[s['prime'] for s in proof['signatures']],
        proof['no_rational_2_torsion_prime'],bank['centres'],str(path.relative_to(ROOT)))

    # An existing alternative subgroup on the same control curve; not a new independent curve sample.
    path = LOCAL/'curve302-seeded-v3-amplifier-v1/recovered-strict-03/replay-M17/epoch-12'
    selection, cloud = bound(path/'selection.json'), bound(path/'cloud-000.json')
    state = cloud['final_state']['state']
    if selection['basis'] != state['reductions']['points'] or len(selection['basis']) != 30:
        raise ArithmeticError('holdout pre-acquisition basis differs')
    projection('curve302-alternative-M30','validation',cloud['curve'],selection['basis'],
        state['reductions']['primes'],state['no_two_torsion_prime'],selection['centres'][:256],str(path.relative_to(ROOT)))

    # Existing highest locally generated rank27 inputs with prepared, certified banks.
    # Stable inventory IDs determine order; no known public rank28 point is a production input.
    for number in (40,48,71,90):
        suffix = 'productive-v3' if number == 90 else 'exact-maximum'
        path = LOCAL/f'curve{number}-{suffix}-discovery-v1'
        seed, selection = bound(path/'seed.json'), bound(path/'epoch-00/landscape/selection.json')
        proof = seed['proof']
        if seed['points'] != selection['basis'] or len(seed['points']) != 27:
            raise ArithmeticError('production seed is not the declared rank27 subgroup')
        excluded = set()
        # Read every locally retained, independently sealed pass with exactly this basis.
        for prior in sorted(LOCAL.glob(f'curve{number}-*discovery*')):
            if not (prior/'verified.json').is_file() or not (prior/'terminal.json').is_file():
                continue
            verified = bound(prior/'verified.json')
            if verified.get('terminal_sha256') != sha(prior/'terminal.json'):
                raise ArithmeticError('prior pass terminal seal differs')
            bound(prior/'terminal.json')
            for selection_path in prior.glob('epoch-*/landscape/selection.json'):
                prior_selection = bound(selection_path)
                if prior_selection['basis'] != seed['points']:
                    continue
                for chart_path in selection_path.parent.parent.glob('chart-*.json'):
                    chart = bound(chart_path)
                    if chart['search']['status'] != 'bounded_search_complete':
                        continue
                    m = chart['mapping']
                    excluded.add((str(m['centre']['representative']),*map(str,box_key(m['matrix'])),str(chart['search']['height_bound'])))
        projection(f'new-20260906-{number}','production',seed['curve'],seed['points'],
            [s['prime'] for s in proof['signatures']],proof['no_rational_2_torsion_prime'],
            selection['centres'],str(path.relative_to(ROOT)),excluded)

    # Exact starting-rank preflight, not just on-curve checks. No search is invoked.
    rank_receipts = {}
    for case in cases:
        seed = read(output/case['input'])['payload']
        proof = checked_rank(tuple(map(F,seed['curve'])),
            [tuple(map(F,p)) for p in seed['points']],seed['primes'],seed['torsion_prime'])
        target = output/'starting-ranks'/(case['id']+'.json')
        checkpoint(target,proof)
        rank_receipts[str(target.relative_to(output))] = sha(target)

    snapshot = output/'runtime/elliptic-curves/cas'
    sources = {}
    for p in sorted(CAS.rglob('*')):
        if not p.is_file() or p.suffix not in ('.py','.sage','.cpp','.h') or '__pycache__' in p.parts:
            continue
        target = snapshot/p.relative_to(CAS)
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(p,target)
        sources[str(target.relative_to(output))] = sha(target)
    # Every packet and source is frozen before execution; no winning-point suffix is projected.
    plan = {'schema':'next-direction-benchmark.v1','status':'FROZEN_NOT_STARTED',
        'cases':sorted(cases,key=lambda c:(['development','validation','production'].index(c['role']),c['id'])),
        'complete_curve302_ladder':ladder,'inputs':bindings,'sources':sources,
        'starting_rank_certificates':rank_receipts,
        'policies':POLICIES,'benchmark_arm_seconds':300,'validation_arm_seconds':450,
        'replay_seconds':90,'map_python':shutil.which('sage'),'seconds_per_map':5,
        'map_rss_bytes':1024**3,'gp_sha256':sha('/usr/bin/gp'),'rss_bytes':3*1024**3,
        'sage_launcher_sha256':sha(shutil.which('sage')),
        'selection':'Unchanged pre-acquisition centre order, first256 for control rungs; complete49 public188 control bank. No winning chart indices or withheld point coordinates enter worker inputs.',
        'metric':'Verified next-direction successes first, then summed process-tree CPU with twice the arm allowance charged for every miss. Charge map preparation, point calls, failures and replay.',
        'gate':'Require at least two development successes and a separately replayed validation gain. Compare winner and baseline on the alternative subgroup before any production point call. If a nonbaseline winner loses that validation comparison, stop.',
        'production':'Only the frozen winning policy, once per frozen rank27 input, stop at its first independently certified gain. Deduplicate completed same-basis boxes. No automatic bank enlargement or successive rank ladder.',
        'boundary':'Retrospective fixed-bank representation comparison. Historical landscape costs are retained separately and unknown complete cold cost remains UNKNOWN. No fresh parameter construction, rank upper bound, speed theorem or success rate.',
        'ceiling':'Execution requires an explicit --hours value, at most8; the previous constructor allowance is not reused.',
        'preparation_wall_seconds':time.monotonic()-started,
        'preparation_cpu_seconds':time.process_time()-started_cpu}
    checkpoint(output/'plan.json',plan)
    print('FROZEN',len(cases),'cases;',len(ladder),'ladder states; no point search',flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args = p.parse_args()
    output = args.output.resolve()
    wall,cpu = time.monotonic(),time.process_time()
    status = 'FAILED_INPUT_PROJECTION'
    try:
        prepare(output)
        status = 'PASS_INPUT_PROJECTION_AND_STARTING_RANKS'
    finally:
        if output.exists():
            checkpoint(output/'preparation-meter.json',{'status':status,
                'wall_seconds':time.monotonic()-wall,'cpu_seconds':time.process_time()-cpu,'point_calls':0})
