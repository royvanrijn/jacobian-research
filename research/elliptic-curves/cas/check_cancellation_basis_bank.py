#!/usr/bin/env sage -python
"""Bounded positive and corruption checks before the bank-cost protocol."""
import copy
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import resource
import shutil
import time
import traceback

from cancellation_basis_early import guard, OUT as EARLY, RAW as RETAINED
from cancellation_basis_epoch import read, sha, need
from cancellation_scheduler_fresh import source_closure
from cancellation_scheduler_prepare import new_write
from finite_cancellation_corpus import ROOT, LOCAL, write
from verify_cancellation_basis_bank import verify_selection
from verify_parity_minimum import ParityMinimumVerifier

OUT = ROOT/'artifacts/generated-results/elliptic-curves/cancellation_basis_bank_certificate_v1'
RAW = LOCAL/'cancellation-basis-bank-certificate-v1'


def run():
    import numpy as np
    resource.setrlimit(resource.RLIMIT_CPU, (40, 45))
    started = time.process_time(); plan, cases = guard()
    case = next(c for c in cases if c['id'] == 'fcf6ce8e6a518eb715d2')
    epoch = RETAINED/'arms'/case['id']/'basis_refresh/epoch-01'
    packet = read(epoch/'seed.json'); original = read(epoch/'producer/selection.json')
    names = source_closure([Path(__file__)])
    recipe = {'status': 'FROZEN_CORRUPTION_CHECKS', 'primary_protocol_sha256': sha(EARLY/'protocol.json'),
        'case': case['id'], 'epoch': 1, 'source_sha256': names, 'cpu_soft_limit': 40, 'cpu_hard_limit': 45,
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in sorted((epoch/'producer').iterdir()) if p.is_file()},
        'seed_sha256': sha(epoch/'seed.json'), 'point_search_calls': 0,
        'bank_corruptions': ['omitted_minimum_tie', 'omitted_selected_coset', 'wrong_anchor_word',
                             'changed_map_recipe_with_updated_hash', 'changed_full_census_with_updated_hash',
                             'changed_complete_centre_order'],
        'boundary': 'Checks of supplied certificates only; no optimizer, point search, or arithmetic exclusion.'}
    new_write(RAW/'corruption-checks/recipe.json', recipe)
    for name in names:
        target = RAW/'corruption-checks/sources'/name
        target.parent.mkdir(parents=True, exist_ok=True); shutil.copyfile(ROOT/name, target)
    rows = []

    def rejects(name, action, expected):
        tick = time.process_time()
        try:
            action()
        except (ArithmeticError, ValueError, RuntimeError, TypeError) as error:
            need(expected in str(error), 'unexpected rejection for '+name+': '+str(error))
            rows.append({'name': name, 'status': 'REJECTED_AS_REQUIRED', 'error': str(error),
                         'cpu_seconds': time.process_time()-tick})
        else:
            raise ArithmeticError('corruption accepted: '+name)

    try:
        # Tiny independently exhaustive controls include an off-diagonal
        # negative entry, an exact zero, and the twofold sign tie.
        for gram, parity in [([[2,1],[1,2]], [1,1]), ([[2,-3],[-3,8]], [1,0]), ([[1]], [0])]:
            check = ParityMinimumVerifier(gram)
            vectors = [v for v in product(range(-6,7), repeat=len(parity)) if [x % 2 for x in v] == parity]
            radius = min(check.norm(v) for v in vectors)
            minima = [v for v in vectors if check.norm(v) == radius]
            result = check.verify(parity, radius, minima)
            rows.append({'name': 'tiny_complete_parity_'+str(len(rows)), 'status': result['status'],
                         'gram': gram, 'parity': parity, 'certificate': result})
        one = ParityMinimumVerifier([[1]])
        rejects('strictly_shorter_vector', lambda: one.verify([1], 9, [[-3],[3]]), 'strictly shorter')
        rejects('missing_sign_tie', lambda: one.verify([1], 1, [[1]]), 'complete minimum set')
        rejects('duplicate_minimum', lambda: one.verify([1], 1, [[1],[1]]), 'duplicate minimum')
        rejects('wrong_parity', lambda: one.verify([0], 1, [[-1],[1]]), 'outside the parity sphere')
        rejects('noninteger_radius', lambda: one.verify([1], 1.0, [[-1],[1]]), 'integer')
        rejects('node_cap', lambda: one.verify([1], 1, [[-1],[1]], 1), 'node limit')
        rejects('asymmetric_metric', lambda: ParityMinimumVerifier([[1,1],[0,1]]), 'asymmetric')
        rejects('nonpositive_metric', lambda: ParityMinimumVerifier([[1,2],[2,1]]), 'positive definite')

        for name in recipe['bank_corruptions']:
            folder = RAW/'corruption-checks'/name
            shutil.copytree(epoch/'producer', folder)
            selection = copy.deepcopy(original); expected = None
            if name == 'omitted_minimum_tie':
                selection['anchors'][0]['refined'][0]['cvp']['minima'].pop()
                expected = 'complete minimum set'
            elif name == 'omitted_selected_coset':
                selection['anchors'][0]['refined'].pop(); expected = 'shortlist omitted'
            elif name == 'wrong_anchor_word':
                selection['anchors'][0]['refined'][0]['representative'][-1] += 2
                expected = 'refined word, profile or multiplicity'
            elif name == 'changed_map_recipe_with_updated_hash':
                p = folder/'anchor-00-maps.json'; maps = read(p)
                maps[0]['matrix'][0] = str(F(maps[0]['matrix'][0])+1)
                write(p, maps); selection['anchors'][0]['maps_sha256'] = sha(p)
                expected = 'prescribed map recipe'
            elif name == 'changed_full_census_with_updated_hash':
                p = folder/'anchor-00-full.npz'
                with np.load(p, allow_pickle=False) as saved:
                    arrays = {key: saved[key] for key in saved.files}
                arrays['residues'][0,0] ^= 1; np.savez_compressed(p, **arrays)
                selection['anchors'][0]['full_scores_sha256'] = sha(p); expected = 'full score array'
            elif name == 'changed_complete_centre_order':
                selection['centres'][0], selection['centres'][1] = selection['centres'][1], selection['centres'][0]
                expected = 'complete exported centre order'
            write(folder/'selection.json', selection)
            rejects(name, lambda: verify_selection(packet, case, plan, folder, selection), expected)
            rows[-1]['selection_sha256'] = sha(folder/'selection.json')
        status = 'PASS_POSITIVE_AND_CORRUPTION_CHECKS'
    except BaseException as error:
        new_write(OUT/'corruption-checks.json', {'status':'RETAINED_CHECK_FAILURE','rows':rows,
            'error':type(error).__name__+': '+str(error),'traceback':traceback.format_exc(),
            'cpu_seconds':time.process_time()-started,'recipe_sha256':sha(RAW/'corruption-checks/recipe.json')})
        raise
    result = {'status':status,'rows':rows,'cpu_seconds':time.process_time()-started,
              'recipe_sha256':sha(RAW/'corruption-checks/recipe.json'),'point_search_calls':0}
    new_write(OUT/'corruption-checks.json', result)
    print({'status':status,'checks':len(rows),'cpu_seconds':result['cpu_seconds']}, flush=True)


if __name__ == '__main__':
    run()
