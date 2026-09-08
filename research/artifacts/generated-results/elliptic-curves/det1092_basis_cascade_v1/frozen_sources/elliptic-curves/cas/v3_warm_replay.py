"""Independent replay for the warm transfer, including schedules beyond index 999.

The numerical landscape checks are adapted from check_visibility_cascade_v3.
The original checker and all original calibration evidence remain unchanged.
No dynamic source rewriting or assertion suppression is used.
"""
from __future__ import annotations

from fractions import Fraction as F
from pathlib import Path
import csv
import json

from v3_warm_support import (atomic, bindings, check_chart, indexed_paths, read,
                             require, sha, point_tuple, curve_tuple, within)


def verify_landscape(engine, model, basis, tested, wd, policy):
    """Recompute the complete cheap landscape and EVERY retained exact CVP."""
    import numpy as np
    from sage.all import ZZ, matrix
    from visibility_lattice_v2 import ExactParity
    from visibility_selection_v3 import cheap_shortlist, final_shortlist, chart_profile

    selection = read(wd/'selection.json')
    require(selection['basis'] == [list(map(str, p)) for p in basis], 'landscape seed mismatch')
    g, u = matrix(ZZ, selection['rounded_gram']), matrix(ZZ, selection['LLL'])
    require(abs(u.det()) == 1 and g.is_positive_definite(), 'invalid lattice decision form')
    exact = ExactParity((u*g*u.transpose()).rows())
    count = 1 << (len(basis)-17)
    require(count == selection['extensions_per_anchor'], 'extension count mismatch')
    columns = engine.fingerprints(model, basis)
    require(columns == selection['fingerprint_columns'], 'fingerprints changed')
    ext_columns = columns[17:]
    rebased = ext_columns[:]
    if len(rebased) > 1:
        rebased[0] ^= rebased[1]
        rebased.reverse()

    def fp_set(cols):
        values = {0}
        for col in cols:
            values |= {value ^ col for value in values}
        return values

    require(fp_set(ext_columns) == fp_set(rebased) and len(fp_set(ext_columns)) == count,
            'extension coverage failed')
    buckets = {8: [], 10: []}
    with engine.ORBITS.open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            shell = int(row['minimum_norm'])
            if shell in buckets:
                word = list(map(int, row['parent_MW17_w'].split())) + [0]*(len(basis)-17)
                buckets[shell].append({'orbit': int(row['orbit_mask']), 'shell': shell,
                                      'representative': word, 'fingerprint': engine.fp(word, columns)})
    chosen, anchors = [], []
    gg = np.asarray(g.rows(), dtype=np.int64)
    for rows in buckets.values():
        require(bool(rows), 'empty generic norm shell')
        words = np.asarray([r['representative'] for r in rows], dtype=np.int64)
        require(len(basis)**2*int(abs(words).max())**2*int(abs(gg).max()) < 2**62,
                'canonical norm overflow')
        norms = np.einsum('ij,jk,ik->i', words, gg, words, optimize=True)
        for row, norm in zip(rows, norms):
            row['metric_norm'] = int(norm)
        rows.sort(key=lambda r: (-r['metric_norm'], r['fingerprint']))
        anchors += rows[:policy['anchors_per_shell']]
        canonical_count = 0
        for row in rows:
            key, word = engine.point_key(model, basis, row['representative'])
            if tuple(map(str, key)) in tested:
                continue
            chosen.append({**row, 'representative': word, 'point': list(map(str, key)), 'lane': 'canonical'})
            canonical_count += 1
            if canonical_count == policy['canonical_per_shell']:
                break
    require(anchors == [r['anchor'] for r in selection['anchors']], 'anchor selection differs')
    mapper = engine.load('factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    exact_count, cosets = 0, 0
    for ai, record in enumerate(selection['anchors']):
        path = wd/f'anchor-{ai:02d}-full.npz'
        require(sha(path) == record['full_scores_sha256'], 'full-score checkpoint hash mismatch')
        with np.load(path, allow_pickle=False) as data:
            residues = data['residues']
            require(residues.shape == (count, len(basis)), 'residue shape mismatch')
            require(bool(np.all((residues == 0) | (residues == 1))), 'nonbinary residues')
            require(bool(np.all(residues[:, :17] == np.asarray(record['anchor']['representative'][:17]) % 2)),
                    'base parity changed')
            require(len(set(map(tuple, residues[:, 17:]))) == count, 'extension masks omitted or duplicated')
            rp = (residues @ (np.asarray(u.inverse().rows(), dtype=np.int64) % 2)) % 2
            require(np.array_equal(rp, data['reduced_residues']), 'wrong reduced parity')
            words, norms = exact.babai(rp)
            require(np.array_equal(words, data['babai_words']) and np.array_equal(norms, data['babai_norms']),
                    'cheap landscape changed')
        keys = [engine.fp(row, columns) for row in residues]
        order = cheap_shortlist(norms, keys)
        require(set(order) == {r['extension'] for r in record['refined']} and
                len(order) == len(record['refined']), 'exact-CVP shortlist differs')
        for row in record['refined']:
            i = row['extension']
            result = exact.solve(rp[i], words[i], policy['exact_cvp_node_limit'])
            require(json.dumps(result, sort_keys=True) == json.dumps(row['cvp'], sort_keys=True),
                    'exact CVP certificate changed')
            require(row['metric_norm'] == result['norm'] and row['fingerprint'] == keys[i], 'CVP identity changed')
            word = matrix(ZZ, 1, len(basis), row['representative'])
            require(int((word*g*word.transpose())[0, 0]) == result['norm'], 'representative norm differs')
            require([int(x) % 2 for x in word.row(0)] == list(residues[i]), 'representative parity differs')
            key, _ = engine.point_key(model, basis, row['representative'])
            require(list(map(str, key)) == row['point'], 'centre point differs')
            minima = []
            for minimum in result['minima']:
                w = list(map(int, (matrix(ZZ, 1, len(basis), minimum)*u).row(0)))
                minima.append(engine.point_key(model, basis, w)[0])
            require(key == min(minima), 'semantic minimum-centre tie changed')
            profile = chart_profile(mapper.mapping(model, basis, row))
            require(all(row[k] == v for k, v in profile.items()), 'quartic ranking features changed')
            require(row['multiplicity'] == len(result['minima'])//2, 'minimum multiplicity differs')
            exact_count += 1
        ranked = sorted(record['refined'], key=lambda r: (-r['metric_norm'], r['fingerprint']))
        chosen += final_shortlist([r for r in ranked if tuple(r['point']) not in tested])
        cosets += count
        print('REPLAY_LANDSCAPE', len(basis), ai+1, '/', len(anchors), flush=True)
    require(cosets == selection['full_cosets_scored'], 'full coset total differs')
    chosen.sort(key=lambda r: (-r['metric_norm'], r['fingerprint'], r['lane'], tuple(map(F, r['point']))))
    unique, seen = [], set(tested)
    for row in chosen:
        key = tuple(row['point'])
        if key not in seen:
            unique.append(row)
            seen.add(key)
    require(unique == selection['centres'], 'final frozen schedule differs')
    # Independently recompute the numerical decision form as well, rather than
    # trusting a saved Gram matrix or advertising it as an exact height theorem.
    geo = engine.load('prospective_half_lattice_v3.sage')
    from sage.all import pari
    hg, asym = geo.canonical_height_gram(model, basis)
    actual_g = matrix(ZZ, geo.rounded_gram(hg, 1000000))
    actual_u = matrix(ZZ, pari(actual_g).qflllgram()).transpose()
    require(actual_g == g and actual_u == u and str(asym) == selection['height_asymmetry'],
            'height metric/LLL replay differs')
    return {'cosets': cosets, 'exact_cvps': exact_count, 'centres': len(unique)}


def replay_case(ctx):
    """Independently check the retained terminal; never re-execute a point search."""
    from v3_warm_engine import certified_state
    from v3_warm_support import terminal_structure
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    import audit_recorded_point_mod2_rank_v3 as mod2
    import audit_retained_cloud_modl as modl

    engine, model, state = ctx.engine, ctx.model, ctx.state
    folder, policy, root = ctx.folder, ctx.policy, ctx.root
    terminal = terminal_structure(folder)
    require(terminal is not None, 'no sealed search terminal to replay')
    cache_dir = folder/'warm-verification'
    cache_dir.mkdir(exist_ok=True)
    reports, tested = [], set()
    mapper = engine.load('factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    for stage in terminal['stages']:
        wd = folder/f"replay-M17/epoch-{stage['epoch']:02d}"
        vd = cache_dir/f"epoch-{stage['epoch']:02d}"
        vd.mkdir(exist_ok=True)
        basis = point_tuple(state.basis)
        selection = read(wd/'selection.json')
        require(point_tuple(selection['basis']) == basis, 'replay basis chain differs')
        landscape_inputs = [wd/'selection.json', folder/'protocol.json', folder/'orbits.tsv']
        landscape_inputs += list(wd.glob('anchor-*-full.npz'))
        landscape_bindings = {str(p.relative_to(root)): sha(p) for p in landscape_inputs}
        cached = vd/'landscape.json'
        if cached.exists():
            record = read(cached)
            require(record['sources'] == ctx.sources and record['bindings'] == landscape_bindings,
                    'independent landscape cache has changed dependencies')
            bindings(root, record['bindings'])
            metrics = record['metrics']
        else:
            metrics = verify_landscape(engine, model, basis, tested, wd, policy)
            atomic(cached, {'sources': ctx.sources, 'bindings': landscape_bindings, 'metrics': metrics}, immutable=True)
        # Numeric file indices are checked against payload indices and the frozen
        # schedule. In particular the 999->1000 boundary is not an ordering rule.
        charts = indexed_paths(wd, expected=stage['charts'])
        if charts:
            # The transcript binds the FULL MWState key, not just its points.
            # Restore and verify the historical state rather than regenerating
            # a mathematically equivalent state with different observations.
            from v3_warm_engine import restore_state
            state = restore_state(read(wd/'cloud-000.json')['final_state'], model, basis)
        prefix = []
        expected_points = list(basis)
        seen = {(x, abs(y)) for x, y in basis}
        for j, path in enumerate(charts):
            chart = read(path)
            check_chart(chart, selection, j)
            require(chart['search']['height_bound'] == policy['height'] and
                    chart['search']['timeout_seconds'] == policy['seconds_per_chart'] and
                    chart['search']['gp_binary_sha256'] == policy['gp_sha256'], 'chart used a different frozen budget/backend')
            prefix.append(chart)
            for raw in chart['search']['finite_curve_points']:
                point = F(raw['x']), F(raw['y'])
                key = point[0], abs(point[1])
                if key not in seen:
                    seen.add(key)
                    expected_points.append(point)
            cloud_path, audit_path = wd/f'cloud-{j:03d}.json', wd/f'mod2-{j:03d}.json'
            snapshot, audit = read(cloud_path), read(audit_path)
            require(audit.get('status') == 'COMPLETE_DECLARED_FINITE_AUDIT', 'partial mod-2 audit')
            require(curve_tuple(snapshot['curve']) == model and point_tuple(snapshot['final_state']['state']['reductions']['points']) == basis,
                    'point snapshot seed/model changed')
            require(snapshot['charts'] == prefix, f'cumulative snapshot changed at chart {j}')
            require(audit['input_sha256'] == sha(cloud_path) and within(root, audit['input_path']) == cloud_path,
                    'audit not bound to the executed snapshot')
            require(curve_tuple(audit['curve']) == model and point_tuple(audit['points']) == tuple(expected_points),
                    'certified cloud contains unreturned/missing points')
            chart_bindings = {str(p.relative_to(root)): sha(p) for p in (path, cloud_path, audit_path)}
            record_path = vd/f'chart-{j:06d}.json'
            if record_path.exists():
                old = read(record_path)
                require(old['sources'] == ctx.sources and old['bindings'] == chart_bindings and
                        old['selection_sha256'] == sha(wd/'selection.json'), 'chart replay cache changed')
            else:
                require(mapper.mapping(model, basis, chart['centre']) == chart['mapping'], 'saved quartic map changed')
                search = PointedQuarticSearch(state=state, centre={'coefficients': chart['centre']['representative']},
                                             coordinate_policy=chart['mapping']['coordinate_policy'])
                backend.replay(search, chart['mapping'], chart['search'])
                mod2.check(audit_path)
                atomic(record_path, {'sources': ctx.sources, 'bindings': chart_bindings,
                                     'selection_sha256': sha(wd/'selection.json')}, immutable=True)
            if j < len(charts)-1:
                require(audit['rank_lower_bound'] == len(basis), 'stale chart executed after a rank gain')
            tested.add(tuple(chart['centre']['point']))
            if j % 25 == 0 or j+1 == len(charts):
                print('REPLAY_CHART', ctx.case, stage['epoch'], j+1, '/', len(charts), flush=True)
        last = read(wd/stage['audit'])
        enlarged = point_tuple(last['independent_points'])
        require(enlarged[:len(basis)] == basis and len(enlarged) == stage['after'], 'certified prefix/rank changed')
        state = certified_state(model, enlarged, last['rank_certificate'])
        # Odd-prime audits are checked ALSO for no-gain terminals.
        odd_data = read(wd/'modl.json')
        require(odd_data['status'] == 'COMPLETE_BOUNDED_QUOTIENT_AUDIT' and len(odd_data['audits']) == 2 and
                all(a['status'] == 'COMPLETE_BOUNDED_QUOTIENT_AUDIT' for a in odd_data['audits']), 'partial odd-prime audit')
        require(odd_data['input_sha256'] == sha(wd/stage['audit']) and odd_data['points'] == last['points'] and
                curve_tuple(odd_data['curve']) == model, 'odd-prime audit is not bound to the final returned point cloud')
        modl.check(wd/'modl.json')
        odd = {str(a['modulus']): a['finite_column_rank'] for a in odd_data['audits']}
        require(odd == {'3': state.rank, '5': state.rank}, 'mod-3/5 replay disagrees')
        all_inputs = {str(p.relative_to(root)): sha(p) for p in wd.iterdir() if p.suffix in ('.json', '.npz')}
        reports.append({'epoch': stage['epoch'], 'before': len(basis), 'after': state.rank,
                        'charts_replayed': len(charts), 'independent_modl_ranks': odd,
                        'landscape': metrics, 'checkpoint_hashes': all_inputs})
        atomic(cache_dir/'progress.json', {'sources': ctx.sources, 'stages': reports})
    require(state.rank == terminal['final_rank_lower_bound'], 'final rank disagrees')
    result = {'schema': 'warm-transfer-independent-replay.v2', 'case': ctx.case,
              'sources': ctx.sources, 'protocol_sha256': sha(folder/'protocol.json'),
              'terminal_sha256': sha(folder/'replay-M17/terminal.json'),
              'rank_lower_bound': state.rank, 'stages': reports, 'charts': terminal['charts'],
              'scope': 'Exact point/rank and finite-schedule replay. No upper bound, new-curve or record claim.'}
    atomic(folder/'warm-replay.json', result, immutable=True)
    final_paths = [folder/'warm-replay.json', folder/'protocol.json', folder/'seed-input.json',
                   folder/'seed-proof.json', folder/'replay-M17/terminal.json']
    verified = {'schema': 'v3-warm-start-transfer-result.v2', 'status': 'PASS_INDEPENDENT_WARM_REPLAY',
                'case': ctx.case, 'initial_rank': policy['initial_rank'], 'rank_lower_bound': state.rank,
                'gain': state.rank-policy['initial_rank'], 'stop_reason': terminal['stop_reason'],
                'charts': terminal['charts'], 'sources': ctx.sources,
                'bindings': {str(p.relative_to(root)): sha(p) for p in final_paths},
                'claim_boundary': 'Warm-start extension of an existing curve; generic 17->17 control stays a no-gain result.'}
    atomic(folder/'warm-verified.json', verified, immutable=True)
    if verified['gain'] > 0:
        atomic(folder/'warm-discovery.json', {'curve': list(map(str, model)),
               'points': [list(map(str, p)) for p in state.basis], 'rank_lower_bound': state.rank,
               'verified_sha256': sha(folder/'warm-verified.json'),
               'status': 'VERIFIED_INCREASE_FOR_EXISTING_CURVE_NOT_A_CATALOGUE_RECORD_CLAIM'}, immutable=True)
    print('INDEPENDENT_WARM_REPLAY_COMPLETE', ctx.case, state.rank, flush=True)
    return verified
