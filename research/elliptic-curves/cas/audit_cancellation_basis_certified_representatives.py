#!/usr/bin/env sage -python
"""Predeclared finite translation audit; cached maps and no point search."""
import argparse
import itertools
from pathlib import Path
import resource
import time

from cancellation_basis_certified import CAS, OUT, RAW, ROOT, guard, need, new_write, read, sha
from cancellation_scheduler_fresh import source_closure
from finite_cancellation_corpus import write


def run():
    plan, _ = guard()
    cap = plan['postmortem']['representative_cpu_limit']
    resource.setrlimit(resource.RLIMIT_CPU, (cap, cap+5))
    from sage.all import EllipticCurve, QQ
    from search_observability import point_visibility
    from pointed_quartic_search import PointedQuarticSearch
    from fractions import Fraction as F
    need(not (OUT/'representative-audit.json').exists(), 'preserve completed supplement')
    need(read(OUT/'mechanism.json')['status'] == 'PASS', 'complete mechanism audit missing')
    start = time.process_time(); rows = []; evaluations = 0; blocks = {str(i): {'targets': 0, 'old_visible': 0, 'old_exclusive': 0, 'strict_deletion': 0, 'both': 0} for i in (0, 1)}
    mechanism = read(OUT/'mechanism.json')
    need(mechanism['protocol_sha256'] == sha(OUT/'protocol.json'), 'mechanism protocol differs')
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
        need(sha(RAW/'mechanism-maps'/f'{source["case"]}.json') == mechanism['map_sha256'][source['case']], 'cached maps changed')
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
        'protocol_sha256': sha(OUT/'protocol.json'), 'mechanism_sha256': sha(OUT/'mechanism.json'),
        'source_sha256': sha(Path(__file__)), 'point_search_calls': 0,
        'boundary': 'Predeclared finite representatives for every actual later gain using new generators; at most18 representatives per point. No full-coset minimum, arithmetic exclusion, or performance implication.'}
    new_write(OUT/'representative-audit.json', result)
    print({k: v for k, v in result.items() if k != 'rows'}, flush=True)


if __name__ == '__main__':
    run()
