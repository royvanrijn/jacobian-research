#!/usr/bin/env python3
"""Bounded post-execution accessibility audit, with no rational-point search.

Old experiment inputs and results remain immutable. All chart choices are
sealed before withheld points enter the per-case diagnostic. This is not a
performance experiment and does not train or promote a scheduler.
"""
import argparse
from collections import Counter
from copy import copy
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

from finite_cancellation_corpus import ROOT, LOCAL, canonical, digest, write, short
from cancellation_basis_amplification import OUT as PRIOR, guard as prior_guard
from cancellation_basis_finish import effective
from cancellation_scheduler_cpu import cpu
from search_observability import primitive

CAS = Path(__file__).resolve().parent
OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_accessibility_v1'
RAW = LOCAL/'cancellation-basis-accessibility-v1'


def read(path): return json.loads(path.read_text())
def sha(path): return digest(path.read_bytes())
def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def new_write(path, value):
    need(not path.exists(), 'preserve existing audit artifact: '+str(path))
    write(path, value)


def freeze():
    plan, inputs = prior_guard()
    cases = {r['id']: r for r in inputs}
    rows = []
    for receipt in read(PRIOR/'supervision.json')['records']:
        if receipt['arm'] != 'basis_refresh':
            continue
        dest, result, _ = effective(receipt)
        for epoch in result['epochs']:
            if epoch['epoch'] == 0:
                continue
            need(epoch['epoch'] == 1 and epoch['rank'] == 19, 'unexpected extension; do not expand scope')
            paths = {key: dest/value for key, value in {
                'initial_bank': 'epoch-00/bank.json', 'new_bank': 'epoch-01/bank.json',
                'bank_verification': 'epoch-01/verification.json', 'result': 'result.json'}.items()}
            need(sha(paths['new_bank']) == epoch['bank_sha256'], 'old bank hash differs')
            need(sha(paths['bank_verification']) == epoch['verification_sha256'], 'old verification differs')
            calls = [r for r in result['calls'] if r['epoch'] == 1]
            rows.append({'id': receipt['case'], 'family': cases[receipt['case']]['family'],
                'inputs': {k: {'path': str(p.relative_to(ROOT)), 'sha256': sha(p)} for k, p in paths.items()},
                'executed_calls': [{'path': str((dest/r['file']).relative_to(ROOT)), 'sha256': r['sha256']} for r in calls]})
    need(len(rows) == 18, 'all eighteen successful rebuilds required')
    names = [Path(__file__), CAS/'verify_cancellation_basis_accessibility.py',
             ROOT/'elkies-k3/scripts/rank_growth.py']
    sources = {**plan['source_sha256'], **{str(p.relative_to(ROOT)): sha(p) for p in names}}
    new_write(OUT/'protocol.json', {
        'status': 'FROZEN_POST_EXECUTION_DIAGNOSTIC', 'cases': rows,
        'source_sha256': sources, 'prior_protocol_sha256': sha(PRIOR/'protocol.json'),
        'endpoint_packets_sha256': sha(PRIOR/'endpoint-packets.json'),
        'height': 125000, 'translations_in_new_generator': [-1, 0, 1], 'signs': [1, -1],
        'prime_bound': 500, 'hard_cpu_seconds_per_case': 90, 'wall_seconds_per_case': 120,
        'point_search_calls': 0,
        'selection': 'All eighteen successfully verified M18 to M19 rebuilds of the sealed ablation, without replacement. All exported initial anchors; all exported refreshed anchors with nonzero coefficient in the new generator; and each such anchor after setting that coefficient to zero. Factor-free plus the same at-most-two retained prime-neighbour maps. All maps are prepared and byte-sealed before opening endpoint points in each case. No CVP rebuild or point search.',
        'oracle': 'After maps are sealed, normalize every retained certified endpoint point onto the search model. Classify each separately above the complete M19 seed at the same finite places. UNKNOWN is not dependence. For every independently certified extra point P, evaluate both signs of P+kG for k=-1,0,1, identically on old and new charts. These six representatives do not exhaust a rational coset.',
        'explanation': 'Report exact minima within the finite chart/representative dictionary, square and point-map witnesses, the new generator and anchor word, and the coefficient-deletion counterfactual. Distinguish exported potential visibility from actually executed coverage and actual gained directions. Existing rank_growth.py supplies numerical Schur/cascade diagnostics only; no numerical score proves rank or visibility.',
        'boundary': 'Retrospective finite accessibility diagnostic on already executed controls. No discovery, performance advantage, optimal representation, arithmetic exclusion, fitting, or prospective selection validation. A positive witness motivates a separately frozen exposure test; it cannot pass that test. Interrupted or missing calculations remain UNKNOWN. Original experiment and all costs stay unchanged.'})
    print(json.dumps({'status': 'FROZEN', 'cases': len(rows), 'protocol_sha256': sha(OUT/'protocol.json')}), flush=True)


def guard():
    plan = read(OUT/'protocol.json')
    for name, value in plan['source_sha256'].items():
        need(sha(ROOT/name) == value, 'audit source changed: '+name)
    need(sha(PRIOR/'protocol.json') == plan['prior_protocol_sha256'], 'prior protocol changed')
    need(sha(PRIOR/'endpoint-packets.json') == plan['endpoint_packets_sha256'], 'endpoint packet changed')
    return plan


def loaded_inputs(row):
    result = {}
    for name, item in row['inputs'].items():
        path = ROOT/item['path']; need(sha(path) == item['sha256'], 'retained input changed')
        result[name] = read(path)
    for item in row['executed_calls']:
        need(sha(ROOT/item['path']) == item['sha256'], 'executed call changed')
    return result


def coordinate(A, anchor, target, matrix):
    """Exact fast height scan; the independent checker uses point_visibility."""
    a, b = map(F, anchor); x, y = map(F, target)
    if x == a and y == b:
        return None
    t = -(3*a*a+A)/(2*b) if x == a else (y+b)/(x-a)
    u, v, w, z = map(F, matrix)
    return primitive(z*t-v, -w*t+u)


def prepare_maps(row, plan, folder):
    from half_lattice_pointed_sieve import linear_combination
    from cancellation_basis_exploration import prepare
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import validate_map
    data = loaded_inputs(row)
    initial, current = data['initial_bank'], data['new_bank']
    seed = current['seed']; model = tuple(map(F, seed['curve']))
    basis = [tuple(map(F, p)) for p in seed['points']]
    need(seed['points'][:18] == initial['seed']['points'] and len(basis) == 19, 'basis prefix differs')
    mapper = SourceFileLoader('basis_accessibility_mapper', str(CAS/'lean_factor_free_pari_mapping.sage')).load_module()
    mapper.pari.allocatemem(256000000, silent=True)
    entries = []
    for i, centre in enumerate(initial['centres']):
        entries.append({'id': f'old-{i:03d}', 'group': 'old', 'bank_index': i,
                        'word': centre['representative']+[0], 'point': centre['point']})
    for i, centre in enumerate(current['centres']):
        word = centre['representative']
        if not word[18]:
            continue
        entries.append({'id': f'new-{i:03d}', 'group': 'new', 'bank_index': i,
                        'word': word, 'point': centre['point']})
        old_word = word[:18]+[0]; point = linear_combination(model, basis, old_word)
        entries.append({'id': f'drop-{i:03d}', 'group': 'drop', 'bank_index': i,
                        'word': old_word, 'point': None if point is None else list(map(str, point))})
    cache = {}
    for entry in entries:
        need(cpu() < plan['hard_cpu_seconds_per_case']-15, 'map preparation allowance exhausted')
        q = entry['point']
        if q is None:
            entry['status'] = 'POINTED_CENTRE_AT_INFINITY'; entry['models'] = []; continue
        key = tuple(q)
        if key not in cache:
            maps = prepare(model, [tuple(map(F, q))], [1], mapper)['models']
            for item in maps:
                search = PointedQuarticSearch(curve=model, subgroup=[], centre={'point': q},
                    coordinate_policy=item['mapping']['coordinate_policy'])
                validate_map(search, item['mapping']); item['chart'] = search.chart_record()
            cache[key] = maps
        entry['models'] = cache[key]; entry['status'] = 'EXACT_CHARTS_PREPARED'
    packet = {'case': row['id'], 'seed': seed, 'entries': entries,
              'protocol_sha256': sha(OUT/'protocol.json'), 'oracle_opened': False,
              'cpu_seconds': cpu(), 'point_search_calls': 0}
    new_write(folder/'maps.json', packet)
    return packet


def scan(maps, representatives, A):
    minima = {}
    for entry in maps['entries']:
        hits = []
        for mi, item in enumerate(entry['models']):
            for ri, rep in enumerate(representatives):
                pair = coordinate(A, entry['point'], rep['point'], item['mapping']['matrix'])
                if pair is not None:
                    hits.append((max(map(abs, pair)), mi, ri, pair))
        if hits:
            h, mi, ri, pair = min(hits)
            minima[entry['id']] = {'height': str(h), 'model': mi, 'representative': ri,
                                    'coordinate': list(map(str, pair))}
    return minima


def best(minima, group):
    rows = [(int(v['height']), key, v) for key, v in minima.items() if key.startswith(group+'-')]
    if not rows:
        return None
    _, key, value = min(rows)
    return {'entry': key, **value}


def evaluate(row, plan, folder, maps):
    from future_point_admission import FinitePointAdmission
    from memory_rank_certificate import checked_rank
    from finite_cancellation_features import add
    import numpy as np
    seed = maps['seed']; curve = tuple(map(F, seed['curve'])); basis = [tuple(map(F, p)) for p in seed['points']]
    endpoint = next(r for r in read(PRIOR/'endpoint-packets.json')['rows'] if r['case'] == row['id'])
    model, points = short(endpoint['packet']['curve'], endpoint['packet']['points'])
    need(model == seed['curve'], 'endpoint model differs')
    base = FinitePointAdmission(curve, basis, prime_bound=plan['prime_bound'])
    geo = SourceFileLoader('accessibility_height_geometry', str(CAS/'prospective_half_lattice_v3.sage')).load_module()
    growth = SourceFileLoader('accessibility_rank_growth', str(ROOT/'elkies-k3/scripts/rank_growth.py')).load_module()
    targets = []; statuses = []
    G = basis[18]
    for index, point in enumerate(points):
        need(cpu() < plan['hard_cpu_seconds_per_case']-10, 'oracle allowance exhausted')
        admission = copy(base)
        admission.pivots = dict(base.pivots); admission.known_keys = set(base.known_keys)
        admission.points = list(base.points); admission.columns = dict(base.columns)
        status = admission.consider(tuple(map(F, point)))
        statuses.append({'endpoint_index': index, **status})
        if status['status'] != 'INDEPENDENT_FINITE_COLUMN':
            continue
        proof = checked_rank(curve, admission.points, admission.primes, seed['proof']['no_rational_2_torsion_prime'])
        target_packet = {'curve': seed['curve'], 'points': [list(map(str, p)) for p in admission.points], 'proof': proof}
        representatives = []
        for k in plan['translations_in_new_generator']:
            P = tuple(map(F, point))
            if k:
                P = add(P, (G[0], k*G[1]), curve[3])
            need(P is not None, 'independent target translated to infinity')
            for sign in plan['signs']:
                representatives.append({'k': k, 'sign': sign, 'point': [str(P[0]), str(sign*P[1])]})
        minima = scan(maps, representatives, curve[3])
        gram, asymmetry = geo.canonical_height_gram(curve, basis+[tuple(map(F, point))])
        g = np.asarray(gram, dtype=float)
        geometry = {'gram': [[str(v) for v in r] for r in gram], 'maximum_asymmetry': str(asymmetry),
            'before': growth.jsonable_metrics(growth.extension_metrics(g, 18, 19)),
            'after': growth.jsonable_metrics(growth.extension_metrics(g, 19, 19)),
            'cascade': growth.jsonable_metrics(growth.cascade_metrics(g, 18, 19, 19)),
            'boundary': 'PARI numerical canonical heights at 384-bit precision; not an exact rank or chart-visibility proof.'}
        targets.append({'endpoint_index': index, 'point': point, 'rank_packet': target_packet,
            'representatives': representatives, 'minima': minima, 'best_old': best(minima, 'old'),
            'best_new': best(minima, 'new'), 'geometry': geometry})
    result = {'case': row['id'], 'family': row['family'], 'status': 'EXACT_SCAN_AWAITING_INDEPENDENT_REPLAY',
        'maps_sha256': sha(folder/'maps.json'), 'protocol_sha256': sha(OUT/'protocol.json'),
        'endpoint_statuses': statuses, 'targets': targets, 'cpu_seconds': cpu(), 'point_search_calls': 0}
    new_write(folder/'scan.json', result)
    return result


def worker(case_id):
    plan = guard(); row = next(r for r in plan['cases'] if r['id'] == case_id)
    resource.setrlimit(resource.RLIMIT_CPU, (plan['hard_cpu_seconds_per_case'], plan['hard_cpu_seconds_per_case']))
    folder = RAW/case_id; new_write(folder/'start.json', {'protocol_sha256': sha(OUT/'protocol.json')})
    maps = prepare_maps(row, plan, folder)
    evaluate(row, plan, folder, maps)
    from verify_cancellation_basis_accessibility import verify
    new_write(folder/'verification.json', verify(row, plan, folder))
    print(json.dumps({'case': case_id, 'status': 'PASS', 'cpu_seconds': cpu()}), flush=True)


def run():
    plan = guard(); records = []
    for row in plan['cases']:
        folder = RAW/row['id']; need(not folder.exists(), 'no audit retry or refill')
        folder.mkdir(parents=True)
        t = time.monotonic(); before = resource.getrusage(resource.RUSAGE_CHILDREN)
        with (folder/'worker.log').open('w') as log:
            try:
                child = subprocess.run(['sage', '-python', str(Path(__file__)), 'worker', '--case', row['id']],
                    stdout=log, stderr=subprocess.STDOUT, timeout=plan['wall_seconds_per_case'])
                code = child.returncode; status = 'COMPLETE' if code == 0 else 'WORKER_FAILURE'
            except subprocess.TimeoutExpired:
                code = None; status = 'WALL_LIMIT'
        after = resource.getrusage(resource.RUSAGE_CHILDREN)
        item = {'case': row['id'], 'status': status, 'returncode': code,
            'charged_cpu_seconds': after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
            'wall_seconds': time.monotonic()-t,
            'files': {p.name: sha(p) for p in folder.iterdir() if p.is_file()}}
        new_write(folder/'receipt.json', item); records.append(item)
        write(OUT/'supervision.json', {'status': 'RUNNING', 'records': records})
        print(json.dumps(item), flush=True)
        if status != 'COMPLETE':
            write(OUT/'supervision.json', {'status': 'STOPPED_RETAINED_FAILURE', 'records': records})
            return
    write(OUT/'supervision.json', {'status': 'COMPLETE', 'records': records})


def report():
    plan = guard(); supervision = read(OUT/'supervision.json')
    need(supervision['status'] == 'COMPLETE' and len(supervision['records']) == len(plan['cases']), 'incomplete audit')
    counts = Counter(); rows = []; examples = []
    for receipt in supervision['records']:
        folder = RAW/receipt['case']
        for name, value in receipt['files'].items():
            need(sha(folder/name) == value, 'audit receipt changed')
        result = read(folder/'scan.json'); verified = read(folder/'verification.json'); maps = read(folder/'maps.json')
        need(verified['status'] == 'PASS', 'unverified audit')
        counts.update(verified['counts'])
        for target in result['targets']:
            old, new = target['best_old'], target['best_new']
            h0, h1 = int(old['height']), int(new['height'])
            strict = h1 <= plan['height'] < h0
            counts['targets'] += 1; counts['strict_new_bank_visibility_targets'] += strict
            counts['old_visible_targets'] += h0 <= plan['height']; counts['new_visible_targets'] += h1 <= plan['height']
            drop = target['minima'].get(new['entry'].replace('new-', 'drop-'))
            direct = strict and drop is not None and int(drop['height']) > plan['height']
            counts['strict_coefficient_deletion_targets'] += direct
            item = {'case': result['case'], 'family': result['family'], 'endpoint_index': target['endpoint_index'],
                'old_minimum_height': str(h0), 'new_minimum_height': str(h1),
                'best_new_anchor_deleted_coefficient_minimum_height': None if drop is None else drop['height'],
                'strict_new_bank_visibility': strict, 'strict_coefficient_deletion_visibility': direct,
                'orthogonal_height_before': target['geometry']['before']['orthogonal_height'],
                'orthogonal_height_after': target['geometry']['after']['orthogonal_height']}
            rows.append(item)
            if h1 <= plan['height']:
                entry = next(e for e in maps['entries'] if e['id'] == new['entry'])
                examples.append({**item, 'new_generator': maps['seed']['points'][18],
                    'anchor_word': entry['word'], 'anchor': entry['point'], 'target': target['point'],
                    'representative': target['representatives'][new['representative']],
                    'best_new': new, 'verification_path': str((folder/'verification.json').relative_to(ROOT))})
    summary = {'status': 'PASS_FINITE_POST_EXECUTION_ACCESSIBILITY_AUDIT', 'cases': len(supervision['records']),
        'counts': dict(counts), 'rows': rows, 'visible_examples': examples,
        'charged_cpu_seconds': sum(r['charged_cpu_seconds'] for r in supervision['records']),
        'point_search_calls': 0, 'protocol_sha256': sha(OUT/'protocol.json'),
        'boundary': plan['boundary']}
    new_write(OUT/'summary.json', summary)
    print(json.dumps({k: v for k, v in summary.items() if k not in ('rows', 'visible_examples')}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('command', choices=['freeze', 'worker', 'run', 'report'])
    parser.add_argument('--case'); args = parser.parse_args()
    if args.command == 'worker': worker(args.case)
    else: globals()[args.command]()
