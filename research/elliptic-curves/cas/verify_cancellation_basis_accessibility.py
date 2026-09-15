#!/usr/bin/env sage -python
"""Independent exact point/coordinate replay; no map reduction or point search."""
from collections import Counter
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from math import isqrt

from cancellation_basis_accessibility import CAS, ROOT, read, sha, need, loaded_inputs
from pointed_box_equivalence import box_key
from search_observability import point_visibility, prepare_chart


def verify(row, plan, folder):
    from sage.all import EllipticCurve, QQ
    import numpy as np
    from pointed_quartic_search import PointedQuarticSearch
    from pari_pointed_backend import validate_map
    data = loaded_inputs(row); maps = read(folder/'maps.json'); result = read(folder/'scan.json')
    need(result['maps_sha256'] == sha(folder/'maps.json'), 'unsealed maps')
    need(maps['oracle_opened'] is False, 'map seal is not target blind')
    need(maps['seed'] == data['new_bank']['seed'], 'wrong enlarged subgroup')
    E = EllipticCurve(QQ, maps['seed']['curve']); basis = [E(p) for p in maps['seed']['points']]
    need(len(basis) == 19, 'unexpected subgroup dimension')
    rank = SourceFileLoader('accessibility_native_rank', str(CAS/'verify_finite_cancellation_cpu.sage')).load_module()
    growth = SourceFileLoader('accessibility_native_growth', str(ROOT/'elkies-k3/scripts/rank_growth.py')).load_module()
    counts = Counter(); proof_rows = []; entries = {e['id']: e for e in maps['entries']}
    expected = {}
    for i, centre in enumerate(data['initial_bank']['centres']):
        expected[f'old-{i:03d}'] = centre['representative']+[0]
    for i, centre in enumerate(data['new_bank']['centres']):
        if centre['representative'][18]:
            expected[f'new-{i:03d}'] = centre['representative']
            expected[f'drop-{i:03d}'] = centre['representative'][:18]+[0]
    need(set(entries) == set(expected), 'chart dictionary omits an anchor or counterfactual')
    for key, word in expected.items():
        entry = entries[key]; need(entry['word'] == word, 'anchor word differs')
        q = sum((int(c)*p for c, p in zip(word, basis)), E(0))
        if q.is_zero():
            need(entry['point'] is None and not entry['models'], 'infinite centre used as affine')
            counts['counterfactual_infinite_centres'] += 1; continue
        need(list(map(str, q.xy())) == entry['point'], 'native anchor identity failed')
        counts[entry['group']+'_anchors'] += 1
        need(1 <= len(entry['models']) <= 3 and entry['models'][0]['name'] == 'factor_free', 'chart vocabulary differs')
        for item in entry['models']:
            chart = item['chart']; prepare_chart(chart)
            search = PointedQuarticSearch(curve=E.a_invariants(), subgroup=[],
                centre={'point': entry['point']}, coordinate_policy=item['mapping']['coordinate_policy'])
            validate_map(search, item['mapping'])
            need(chart == search.chart_record(), 'chart record differs')
            counts['exact_maps'] += 1
    calls = [read(ROOT/r['path']) for r in row['executed_calls']]
    executed = {}
    for call in calls:
        anchor = call['search']['base_point']; q = (F(anchor['x']), F(anchor['y']))
        key = q, box_key(call['mapping']['matrix'])
        executed.setdefault(key, []).append(call)
    for target in result['targets']:
        native = rank.sage_rank(target['rank_packet'])
        need(native['rank'] == 20 and target['rank_packet']['points'][:19] == maps['seed']['points'], 'target independence differs')
        need(target['rank_packet']['points'][19] == target['point'], 'target proof mismatch')
        P = E(target['point']); G = basis[18]
        expected_reps = [(k, sign) for k in plan['translations_in_new_generator'] for sign in plan['signs']]
        need([(r['k'], r['sign']) for r in target['representatives']] == expected_reps, 'oracle dictionary differs')
        for rep in target['representatives']:
            q = rep['sign']*(P+rep['k']*G)
            need(not q.is_zero() and list(map(str, q.xy())) == rep['point'], 'target translation identity failed')
        minima = {}
        for key, entry in entries.items():
            hits = []
            for mi, item in enumerate(entry['models']):
                for ri, rep in enumerate(target['representatives']):
                    witness = point_visibility(item['chart'], rep['point'])
                    need(witness['status'] == 'OBSERVABLE_WITHOUT_TRANSCRIPT', 'independent target is a known anchor')
                    pair = tuple(map(int, witness['coordinate']))
                    hits.append((max(map(abs, pair)), mi, ri, pair, witness))
                    counts['exact_coordinate_square_checks'] += 1
            if not hits:
                continue
            h, mi, ri, pair, witness = min(hits, key=lambda v: v[:4])
            value = {'height': str(h), 'model': mi, 'representative': ri, 'coordinate': list(map(str, pair))}
            minima[key] = value
            if h <= plan['height']:
                item = entry['models'][mi]; rep = target['representatives'][ri]
                search = PointedQuarticSearch(curve=E.a_invariants(), subgroup=[], centre={'point': entry['point']},
                    coordinate_policy=item['mapping']['coordinate_policy'])
                root = int(witness['square_root_absolute'])
                expected_point = tuple(map(F, rep['point']))
                need(expected_point in {search.map_hit(*pair, root), search.map_hit(*pair, -root)}, 'square does not map to target')
                key0 = tuple(map(F, entry['point'])), box_key(item['mapping']['matrix'])
                matching = executed.get(key0, [])
                completed = [c for c in matching if c['search']['status'] == 'bounded_search_complete' and c['search']['height_bound'] >= h]
                recorded = any(any((F(p['x']), F(p['y'])) == expected_point for p in c['search']['finite_curve_points']) for c in completed)
                need(not completed or recorded, 'visible target missing from executed completed box')
                counts[entry['group']+'_visible_anchor_target_pairs'] += 1
                if entry['group'] == 'new':
                    counts['new_visible_pairs_executed_complete'] += bool(completed)
                    counts['new_visible_pairs_unexecuted'] += not bool(matching)
                    counts['new_visible_pairs_executed_incomplete_only'] += bool(matching) and not bool(completed)
                proof_rows.append({'endpoint_index': target['endpoint_index'], 'entry': key,
                    'model': mi, 'representative': ri, 'height': h, 'witness': witness,
                    'matching_executed_calls': len(matching), 'completed_executed_calls': len(completed),
                    'recorded_in_completed_call': recorded})
        need(minima == target['minima'], 'independent full finite minima differ')
        for name in ('old', 'new'):
            candidates = [(int(v['height']), k, v) for k, v in minima.items() if k.startswith(name+'-')]
            _, key, value = min(candidates)
            need(target['best_'+name] == {'entry': key, **value}, 'best bank witness differs')
        geometry = target['geometry']; gram = np.asarray(geometry['gram'], dtype=float)
        need(gram.shape == (20, 20) and np.allclose(gram, gram.T), 'numerical Gram malformed')
        need(geometry['before'] == growth.jsonable_metrics(growth.extension_metrics(gram, 18, 19)), 'old Schur calculation differs')
        need(geometry['after'] == growth.jsonable_metrics(growth.extension_metrics(gram, 19, 19)), 'new Schur calculation differs')
        need(geometry['cascade'] == growth.jsonable_metrics(growth.cascade_metrics(gram, 18, 19, 19)), 'cascade calculation differs')
        counts['independent_target_rank_packets'] += 1
    counts['endpoint_classification_unknowns'] += sum(r['status'] == 'UNKNOWN_FINITE_COLUMN_IN_SPAN' for r in result['endpoint_statuses'])
    return {'status': 'PASS', 'case': row['id'], 'maps_sha256': sha(folder/'maps.json'),
        'scan_sha256': sha(folder/'scan.json'), 'counts': dict(counts), 'visible_witnesses': proof_rows,
        'point_search_calls': 0,
        'boundary': 'Native Sage anchor/translation identities and complete finite-group independence; existing exact observability and point map replay for the finite dictionary. Numerical Gram and Schur metrics are diagnostics only. Original completed-search coverage retains its pinned PARI trust boundary.'}
