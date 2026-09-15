#!/usr/bin/env sage -python
"""Post-protocol finite translation audit; cached maps and no point search."""
import argparse
import itertools
from pathlib import Path
import resource
import time

from cancellation_basis_early import CAS, OUT, RAW, ROOT, guard, need, new_write, read, sha
from cancellation_scheduler_fresh import source_closure
from finite_cancellation_corpus import write


def freeze():
    main, _ = guard(); mechanism = read(OUT/'mechanism.json')
    need(mechanism['status'] == 'PASS' and len(mechanism['rows']) == 28, 'actual-gain audit missing')
    paths = {OUT/'mechanism.json', OUT/'mechanism-process.json', OUT/'summary.json'}
    for row in mechanism['rows']:
        paths.add(RAW/'mechanism-maps'/f'{row["case"]}.json')
        paths.add(RAW/'arms'/row['case']/'basis_refresh'/f'epoch-{row["epoch"]:02d}'/'bank.json')
    plan = {'status': 'FROZEN_POST_PROTOCOL_REPRESENTATIVE_SUPPLEMENT',
        'selection': 'All 28 actual later candidate gains in the completed primary comparison; no omitted or replacement rows.',
        'representatives': 'Both signs of P + sum k_i G_i for every generator admitted above the initial M18 before the call, with each k_i in {-1,0,1}. At most two added generators and 18 representatives per target. Use the same dictionary on the original bank, each single-coefficient-deletion bank, and the executed new chart.',
        'maps': 'Use only exact maps retained by the completed mechanism audit. Prepare no new models or banks; run no point search.',
        'purpose': 'Test whether the original P,-P visibility differences persist under a larger fixed finite dictionary for the same quotient direction. This post-protocol sensitivity cannot change the predeclared performance or mechanism gates.',
        'cpu_limit': 90, 'primary_protocol_sha256': sha(OUT/'protocol.json'),
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)},
        'source_sha256': source_closure([Path(__file__)]),
        'boundary': 'A finite representative dictionary, not a minimum over the full subgroup coset. A surviving or disappearing witness does not prove arithmetic exclusion or a performance advantage.'}
    new_write(OUT/'representative-protocol.json', plan)
    old = OUT/'dependency-audit.json'
    preserved = OUT/'dependency-audit-before-representatives.json'
    need(not preserved.exists(), 'preserve original dependency record')
    preserved.write_bytes(old.read_bytes())
    dependency = read(old)
    for name, value in plan['source_sha256'].items():
        if name in main['source_sha256']:
            need(value == main['source_sha256'][name], 'frozen primary source differs')
        else:
            dependency['extra_source_sha256'][name] = value
    dependency.update(previous_dependency_audit_sha256=sha(preserved),
        representative_protocol_sha256=sha(OUT/'representative-protocol.json'),
        boundary='Post-protocol retention plus finite representative sensitivity. The preceding dependency record is preserved verbatim; no primary source seal or decision gate changes.')
    write(old, dependency)
    print({'status': plan['status'], 'protocol_sha256': sha(OUT/'representative-protocol.json')}, flush=True)


def run():
    plan = read(OUT/'representative-protocol.json')
    resource.setrlimit(resource.RLIMIT_CPU, (plan['cpu_limit'], plan['cpu_limit']+5))
    from sage.all import EllipticCurve, QQ
    from search_observability import point_visibility
    from pointed_quartic_search import PointedQuarticSearch
    from fractions import Fraction as F
    need(not (OUT/'representative-audit.json').exists(), 'preserve completed supplement')
    guard()
    for name, value in {**plan['input_sha256'], **plan['source_sha256']}.items():
        need(sha(ROOT/name) == value, 'supplement input changed: '+name)
    start = time.process_time(); rows = []; evaluations = 0; blocks = {str(i): {'targets': 0, 'old_visible': 0, 'old_exclusive': 0, 'strict_deletion': 0, 'both': 0} for i in (0, 1)}
    mechanism = read(OUT/'mechanism.json')
    for source in mechanism['rows']:
        folder = RAW/'arms'/source['case']/'basis_refresh'/f'epoch-{source["epoch"]:02d}'
        bank = read(folder/'bank.json'); seed = bank['seed']; initial = 18
        need(len(seed['points']) == source['pre_call_rank'], 'pre-call rank differs')
        need(sha(folder/'bank.json') == source['bank_sha256'], 'bank changed')
        E = EllipticCurve(QQ, seed['curve']); P = E(source['point'])
        generators = [E(p) for p in seed['points'][initial:]]
        need(1 <= len(generators) <= 2, 'dictionary exceeds declared dimension')
        representatives = []
        for coefficients in itertools.product((-1, 0, 1), repeat=len(generators)):
            R = P + sum((k*G for k, G in zip(coefficients, generators)), E(0))
            for sign in (1, -1):
                Q = sign*R
                need(not Q.is_zero(), 'independent representative became zero')
                representatives.append({'coefficients': list(coefficients), 'sign': sign, 'point': list(map(str, Q.xy()))})
        old = read(RAW/'mechanism-maps'/f'{source["case"]}.json')['old_models']
        need(old, 'original map dictionary missing')

        def minimum(anchors):
            nonlocal evaluations
            best = None
            for ai, anchor in enumerate(anchors):
                for mi, model in enumerate(anchor['models']):
                    for ri, rep in enumerate(representatives):
                        witness = point_visibility(model['chart'], rep['point']); evaluations += 1
                        need(witness['status'] == 'OBSERVABLE_WITHOUT_TRANSCRIPT', 'known endpoint unexpectedly independent')
                        height = max(abs(int(v)) for v in witness['coordinate'])
                        key = (height, ai, mi, ri)
                        if best is None or key < best[0]:
                            best = (key, witness, model['chart'])
            if best is None:
                return None
            (height, ai, mi, ri), witness, chart = best
            search = PointedQuarticSearch(**chart['input'])
            need(search.chart_record() == chart, 'winning chart reconstruction differs')
            pair = list(map(int, witness['coordinate'])); root = int(witness['square_root_absolute'])
            need(tuple(map(F, representatives[ri]['point'])) in {search.map_hit(*pair, root), search.map_hit(*pair, -root)}, 'winning square map differs')
            return {'height': str(height), 'anchor_index': ai, 'model': mi, 'representative': ri, 'witness': witness}

        original = minimum(old)
        deleted = [minimum([drop]) for drop in source['counterfactuals']]
        actual = minimum([{'models': [{'chart': source['actual_chart']}]}])
        need(int(actual['height']) <= 125000, 'actual search witness lost')
        exclusive = int(original['height']) > 125000
        strict = any(d is not None and int(d['height']) > 125000 for d in deleted)
        block = blocks[str(source['validation_block'])]
        block['targets'] += 1; block['old_visible'] += not exclusive; block['old_exclusive'] += exclusive
        block['strict_deletion'] += strict; block['both'] += exclusive and strict
        rows.append({'case': source['case'], 'validation_block': source['validation_block'],
            'call': source['call'], 'call_sha256': source['call_sha256'], 'point': source['point'],
            'generators': seed['points'][initial:], 'representatives': representatives,
            'original_minimum': original, 'deleted_minima': deleted, 'actual_minimum': actual,
            'old_exclusive': exclusive, 'strict_deletion': strict})
        write(RAW/'representative-prefix.json', {'status': 'RUNNING', 'rows': rows, 'evaluations': evaluations, 'cpu_seconds': time.process_time()-start})
        print({'targets': len(rows), 'evaluations': evaluations, 'old_exclusive': exclusive, 'strict_deletion': strict}, flush=True)
    result = {'status': 'PASS_FINITE_REPRESENTATIVE_SENSITIVITY', 'rows': rows, 'blocks': blocks,
        'coordinate_square_evaluations': evaluations, 'cpu_seconds': time.process_time()-start,
        'protocol_sha256': sha(OUT/'representative-protocol.json'), 'point_search_calls': 0,
        'boundary': plan['boundary']}
    new_write(OUT/'representative-audit.json', result)
    print({k: v for k, v in result.items() if k != 'rows'}, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['freeze', 'run'])
    globals()[parser.parse_args().command]()
