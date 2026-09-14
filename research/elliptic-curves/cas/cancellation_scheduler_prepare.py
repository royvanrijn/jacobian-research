#!/usr/bin/env sage -python
"""Freeze new controls and train the residue/exposure scheduler on retained data."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import gzip
import json
from pathlib import Path
import statistics
import sys
import time

from finite_cancellation_corpus import ROOT, OUT as CORPUS, LOCAL, canonical, digest, write
from finite_cancellation_validation_audit import OUT as V2, FAMILIES

OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_scheduler_v3'
RAW = LOCAL/'cancellation-scheduler-v3'
CAS = Path(__file__).resolve().parent
ARMS = ['factor_free', 'adaptive_uniform', 'adaptive_local']


def new_write(path, value):
    if path.exists():
        raise FileExistsError(f'Preserve existing evidence: {path}')
    write(path, value)


def roster():
    design = {
        'status': 'DESIGN_BEFORE_NEW_CONTROL_OUTCOMES',
        'arms': ARMS, 'heights': [8000, 32000, 125000], 'maximum_centres': 48,
        'search_cpu_seconds': 40, 'maximum_calls': 144, 'point_wall_seconds': 5,
        'arm_wall_seconds': 120, 'hard_process_cpu_seconds': 100,
        'selection': 'Reuse V2 eligible bank audit; exclude all 53 prior CPU curves and Curve302. All remaining deep j groups, at most six per family, then six shallow groups per family by SHA256(cancellation-scheduler-v3/ + j). Largest starting rank then lexical bank for deep; lexical bank for shallow. No new eligibility census or missing-artifact reconstruction.',
        'train': 'V1 corpus trains one clipped mixture of uniform and differential root-ball probabilities; weight one per j group. V2 factor-free completed calls train radial detection thresholds and anchor hazards; all failures retained as censored at H125000. No new control j group may enter either fit. Fixed 32-point midpoint real quadrature and H levels. The full corpus already informed development.',
        'scheduler': 'Greedy estimated new union-of-boxes detection mass per complete CPU. Fit empirical radial thresholds from independent development witnesses. Joint categorical states at the same prime; product surrogate only across distinct primes. Update cost from completed calls and remove overlap after complete misses. Timeout removes no exposure. New anchors compete with prepared jobs. Baseline retains its order and H125000.',
        'cancellation_ablation': 'adaptive_uniform builds the identical prime-neighbour bank and uses uniform P1 residue mass. adaptive_local replaces it with the training-calibrated q-only differential mass. No full gcd residue tree and no height-regression score.',
        'primary_gate': 'adaptive_local must recover at least as many directions as each comparator, have at least 10 percent higher recoveries per CPU, and a paired family/stratum bootstrap lower 97.5 percent bound above one against BOTH comparators. 10000 draws with seed 20260914. Each of the two comparisons uses the central 97.5 percent interval (Bonferroni two-comparison family). No post-outcome tuning/repeats. Zero reference gains make the gate UNKNOWN.',
        'cost': 'Outer process CPU includes interpreter, all candidate models, policy scoring, failed calls, exact admission, transcript replay and independent Sage finite-quotient rank certification. Shared retained landscapes are the common input interface, not free cold construction. Offline fitting and corpus replay are retained separately and charged once in campaign-inclusive totals.',
        'fresh_gate': 'Fresh-fibre work is conditional on BOTH primary comparisons passing with no infrastructure failure. Before any fresh execution freeze generic-only inputs, previously unsearched addresses, source hashes, common cold construction and a separate finite paired CPU budget. Failure does not trigger expansion.',
        'scope': 'This tests a finite adaptive policy. It cannot prove universal non-predictiveness, rank32, a rank upper bound or saturation. Exact new points prove only independently certified lower bounds.'}
    new_write(OUT/'design.json', design)
    old = json.loads((V2/'inputs.json').read_text())+json.loads((CORPUS/'cpu/inputs.json').read_text())
    corpus = json.loads(gzip.decompress((CORPUS/'corpus.json.gz').read_bytes()))
    excluded = {c['j_group'] for c in old} | {c['j_group'] for c in corpus if c['family']=='Curve302-development'}
    candidates = json.loads(gzip.decompress((V2/'eligible_banks.json.gz').read_bytes()))
    chosen, used = [], set(excluded)
    for stratum in ('deep', 'shallow'):
        key = (lambda r: (-r['initial_rank'], r['source_landscape'])) if stratum=='deep' else (lambda r:r['source_landscape'])
        banks = {}
        for row in sorted((r for r in candidates if r['stratum']==stratum and r['j_group'] not in used), key=key):
            banks.setdefault(row['j_group'], row)
        for family in FAMILIES:
            eligible = sorted((r for r in banks.values() if r['family']==family),
                              key=lambda r:digest(('cancellation-scheduler-v3/'+r['j_group']).encode()))
            for row in eligible[:6]:
                chosen.append(dict(row)); used.add(row['j_group'])
    inputs, oracle = [], []
    for row in chosen:
        for field, hashfield in [('source_seed','seed_sha256'), ('source_landscape','landscape_sha256')]:
            if digest((ROOT/row[field]).read_bytes()) != row[hashfield]:
                raise ArithmeticError('retained bank changed')
        seed = json.loads((ROOT/row['source_seed']).read_text())
        landscape = json.loads((ROOT/row['source_landscape']).read_text())
        if seed['points'] != landscape['basis']:
            raise ArithmeticError('basis binding changed')
        ident = digest(canonical([row['j_group'], row['source_landscape']]))[:20]
        row['id'] = ident
        inputs.append({'id':ident, 'j_group':row['j_group'], 'family':row['family'], 'stratum':row['stratum'],
                       'seed':{k:seed[k] for k in ('curve','points','proof')},
                       'centres':[x['representative'] for x in landscape['centres'][:48]]})
        endpoint = next(c for c in corpus if c['id']==row['corpus_id'])
        oracle.append({'id':ident, 'curve':endpoint['curve'], 'points':endpoint['generic_points']+endpoint['targets'],
                       'source':endpoint['source'], 'proof_sha256':endpoint['proof_sha256']})
    new_write(OUT/'inputs.json', inputs)
    new_write(OUT/'oracle.json', oracle)
    new_write(OUT/'roster.json', {'chosen':chosen, 'excluded_training_j_groups':sorted(used),
        'prior_cpu_j_groups':sorted(excluded), 'design_sha256':digest((OUT/'design.json').read_bytes()),
        'eligible_banks_sha256':digest((V2/'eligible_banks.json.gz').read_bytes()),
        'inputs_sha256':digest((OUT/'inputs.json').read_bytes()), 'oracle_sha256':digest((OUT/'oracle.json').read_bytes())})
    print(json.dumps({'status':'FROZEN_NEW_ROSTER', 'cases':len(inputs),
                      'deep':sum(c['stratum']=='deep' for c in inputs)}), flush=True)


def coordinate_height(point, anchor, matrix):
    from search_observability import primitive
    x, y = map(F, point); a, b = map(F, anchor)
    A, B, C, D = map(F, matrix)
    heights = []
    for sign in (1, -1):
        slope = (sign*y+b)/(x-a)
        u, v = primitive(D*slope-B, -C*slope+A)
        heights.append(max(abs(u), abs(v)))
    return min(heights)


def train():
    from cancellation_scheduler import root_distribution
    from half_lattice_pointed_sieve import linear_combination
    start = time.process_time()
    roster_data = json.loads((OUT/'roster.json').read_text())
    excluded = set(roster_data['excluded_training_j_groups'])
    groups = defaultdict(list); withheld = defaultdict(list); bindings = {}
    for path in sorted((CORPUS/'cases').glob('*.json.gz')):
        packet = json.loads(gzip.decompress(path.read_bytes()))
        if packet['status'] != 'PASS_EXACT_ACCESSIBILITY':
            continue
        group = packet['j_group']
        if group in roster_data['prior_cpu_j_groups']:
            continue
        prep = packet['prepared']; models = prep['models']; roots = defaultdict(set)
        for model in models[1:]:
            roots[model['prime']].add(model['residue'])
        distributions = {p:root_distribution(list(map(int, models[0]['q'])), p, sorted(rs))[0] for p,rs in roots.items()}
        target_rows = [o for o in packet['observations'] if o['model_index']==0]
        observations = []
        for model in models[1:]:
            p, r = model['prime'], model['residue']
            for obs in target_rows:
                m, n = map(int, obs['coordinate'])
                actual = int(n%p != 0 and (m-r*n)%p == 0)
                observations.append([1/(p+1), distributions[p][r], actual])
        (withheld if group in excluded else groups)[group].extend(observations)
        bindings[path.name] = digest(path.read_bytes())
        if len(bindings)%500 == 0:
            print(json.dumps({'corpus_anchors':len(bindings), 'cpu':time.process_time()-start}), flush=True)
        if time.process_time()-start > 600:
            raise TimeoutError('offline 600 CPU second cap')
    def fit(gs):
        numerator = denominator = 0.
        for rows in gs.values():
            for u, d, y in rows:
                numerator += (d-u)*(y-u)/len(rows)
                denominator += (d-u)**2/len(rows)
        # Keep at least 2% uniform exploration; no censored cell is excluded.
        return min(.98, max(0., numerator/denominator)) if denominator else 0.
    def brier(gs, beta):
        return statistics.mean(sum((u+beta*(d-u)-y)**2 for u,d,y in rs)/len(rs) for rs in gs.values() if rs)
    beta = fit(groups)
    folds=[]
    for fold in range(5):
        testing={g:r for g,r in groups.items() if int(g[:8],16)%5==fold}
        training={g:r for g,r in groups.items() if int(g[:8],16)%5!=fold}
        weight=fit(training)
        folds.append({'fold':fold, 'training_j_groups':len(training), 'testing_j_groups':len(testing),
                      'weight':weight, 'uniform_brier':brier(testing,0), 'local_brier':brier(testing,weight)})
    # Completed V2 calls are development only. Censoring is retained instead
    # of manufacturing actual target heights for unsuccessful charts.
    radial, backend, rows = [], [], []
    calls=[0,0,0]; successes=[0,0,0]
    inputs=json.loads((V2/'inputs.json').read_text())
    for case in inputs:
        if case['j_group'] in {c['j_group'] for c in json.loads((OUT/'inputs.json').read_text())}:
            raise ArithmeticError('development/validation leakage')
        folder=LOCAL/'finite-cancellation-validation-v2/arms'/case['id']/'factor_free'
        result=json.loads((folder/'result.json').read_text())
        for record in result['charts']:
            path=folder/f'chart-{record["index"]:03d}.json'; row=json.loads(path.read_text())
            if row['search']['status']!='bounded_search_complete':
                continue
            ci=row['index']; bucket=0 if ci<4 else 1 if ci<16 else 2
            calls[bucket]+=1; backend.append(row['search']['search_cpu_ms']/1000)
            h=None
            if row.get('new_point'):
                anchor=linear_combination(tuple(map(F,case['seed']['curve'])),
                    [tuple(map(F,p)) for p in case['seed']['points']],row['centre'])
                h=coordinate_height(row['new_point'],anchor,row['mapping']['matrix'])
                if h>125000:
                    raise ArithmeticError('development independent witness not in searched box')
                radial.append(h);successes[bucket]+=1
            rows.append({'case':case['id'],'index':ci,'independent_height':h,
                         'censor_height':125000,'source':str(path.relative_to(ROOT)), 'sha256':digest(path.read_bytes())})
    fitted={'local_weight':beta, 'anchor_priors':[(s+1)/(n+10) for s,n in zip(successes,calls)],
        'radial_heights':sorted(radial), 'backend_cpu_at_125000':statistics.median(backend),
        'per_call_overhead_cpu':.015, 'independent_certificate_cpu_prior':1.0,
        'calls_by_anchor_bucket':calls, 'successes_by_anchor_bucket':successes,
        'training_j_groups':sorted(groups), 'withheld_j_groups':sorted(withheld),
        'crossfit':folds, 'withheld_brier':{'uniform':brier(withheld,0), 'local':brier(withheld,beta)},
        'corpus_sha256':digest((CORPUS/'corpus.json.gz').read_bytes()), 'corpus_cases_sha256':bindings,
        'development_rows':rows, 'cpu_seconds':time.process_time()-start,
        'boundary':'Calibration uses target membership, not target input at execution. Radius observations are conservative first-certified witnesses on 39 positive completed development charts; other completed charts are detection-censored, not absence of points. Brier scores cannot promote the policy.'}
    new_write(OUT/'fit.json', fitted)
    print(json.dumps({k:v for k,v in fitted.items() if k in ('local_weight','anchor_priors','radial_heights','withheld_brier','crossfit','cpu_seconds')}), flush=True)


def preflight():
    from importlib.machinery import SourceFileLoader
    from finite_cancellation_corpus import short, pointkey
    from v3_warm_engine import certified_state
    checker=SourceFileLoader('scheduler_preflight_sage',str(CAS/'verify_finite_cancellation_cpu.sage')).load_module().sage_rank
    start=time.process_time(); rows=[]
    oracle={o['id']:o for o in json.loads((OUT/'oracle.json').read_text())}
    for case in json.loads((OUT/'inputs.json').read_text()):
        seed=case['seed']; o=oracle[case['id']]
        certified_state(tuple(map(F,seed['curve'])),[tuple(map(F,p)) for p in seed['points']],seed['proof'])
        seed_check=checker(seed)
        path=o['source']; split=path.index('.json')+5
        source=ROOT/path[:split]; packet=json.loads(source.read_text())
        for k in path[split+1:].split('/'):
            if k:packet=packet[int(k)] if isinstance(packet,list) else packet[k]
        proof=packet.get('rank_certificate',packet.get('proof'))
        if digest(canonical(proof))!=o['proof_sha256']:
            raise ArithmeticError('endpoint proof changed')
        model=packet.get('curve',packet.get('ainvs',packet.get('model')))
        curve,points=short(model,packet['points']); basecurve,basepoints=short(seed['curve'],seed['points'])
        if curve!=basecurve or not {pointkey(p) for p in basepoints} < {pointkey(p) for p in points}:
            raise ArithmeticError('starting subgroup is not a strict certified subset')
        endpoint_check=checker({'curve':model,'points':packet['points'],'proof':proof})
        rows.append({'id':case['id'],'seed':seed_check,'endpoint':endpoint_check,
                     'source_sha256':digest(source.read_bytes())})
    new_write(OUT/'preflight.json',{'status':'PASS', 'rows':rows,'cpu_seconds':time.process_time()-start,
        'inputs_sha256':digest((OUT/'inputs.json').read_bytes()),'oracle_sha256':digest((OUT/'oracle.json').read_bytes())})
    print(json.dumps({'status':'PASS', 'cases':len(rows),'cpu':time.process_time()-start}),flush=True)


def seal():
    from pari_pointed_backend import sources
    plan=json.loads((OUT/'design.json').read_text())
    fit=json.loads((OUT/'fit.json').read_text()); pre=json.loads((OUT/'preflight.json').read_text())
    inputs=json.loads((OUT/'inputs.json').read_text())
    if pre['status']!='PASS' or len(pre['rows'])!=len(inputs):raise ArithmeticError('preflight incomplete')
    training_replay=json.loads((OUT/'training-replay.json').read_text())
    if training_replay['status']!='PASS' or training_replay['fit_sha256']!=digest((OUT/'fit.json').read_bytes()):
        raise ArithmeticError('training replay missing or stale')
    if set(fit['training_j_groups'])&{c['j_group'] for c in inputs}:raise ArithmeticError('training leakage')
    names=['cancellation_scheduler.py','cancellation_scheduler_prepare.py','cancellation_scheduler_cpu.py',
        'run_cancellation_scheduler.py','verify_cancellation_scheduler.py','report_cancellation_scheduler.py',
        'verify_cancellation_scheduler_training.py',
        'finite_cancellation_features.py','finite_cancellation_validation_features.py',
        'lean_factor_free_pari_mapping.sage','future_point_admission.py','memory_rank_certificate.py',
        'v3_warm_engine.py','verify_finite_cancellation_cpu.sage','pointed_box_equivalence.py']
    plan.update(status='FROZEN_BEFORE_NEW_CPU', cases=len(inputs),
        fit={k:fit[k] for k in ('local_weight','anchor_priors','radial_heights','backend_cpu_at_125000',
            'per_call_overhead_cpu','independent_certificate_cpu_prior')},
        input_sha256={n:digest((OUT/n).read_bytes()) for n in ('design.json','roster.json','inputs.json','oracle.json','fit.json','preflight.json','training-replay.json')},
        source_sha256={**sources(), **{str((CAS/n).relative_to(ROOT)):digest((CAS/n).read_bytes()) for n in names}},
        gp_sha256=digest(Path('/usr/bin/gp').read_bytes()))
    new_write(OUT/'protocol.json',plan)
    print('FROZEN',digest((OUT/'protocol.json').read_bytes()),flush=True)


if __name__=='__main__':
    if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['roster','train','preflight','seal'])
    globals()[parser.parse_args().command]()
