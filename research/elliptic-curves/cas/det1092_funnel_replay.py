"""Prospective M18 replay adapted from the retained generic-M17 replay.

The prior replay and numerical V3 source files remain untouched. This adapter
checks initial rank18, retains exact map/cloud/finite-group verification and
publishes queue labels only after independent replay.
"""
from __future__ import annotations

from fractions import Fraction as F
from v3_warm_support import (atomic, bindings, check_chart, curve_tuple, indexed_paths,
                            point_tuple, read, require, sha, terminal_structure, within)
import det1092_funnel as funnel


def replay(ctx):
    from v3_warm_engine import certified_state, restore_state
    from v3_warm_replay import verify_landscape
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    import audit_recorded_point_mod2_rank_v3 as mod2
    import audit_retained_cloud_modl as modl

    engine, model, state = ctx.engine, ctx.model, ctx.state
    folder, policy, root = ctx.folder, ctx.policy, ctx.root
    terminal = terminal_structure(folder)
    require(terminal is not None, 'no sealed terminal to replay')
    require(terminal['initial_rank'] == policy['initial_rank'] == state.rank == 18, 'replay requires certified M18')
    require(terminal.get('execution_sources') == ctx.sources, 'search implementation binding differs')
    cache_dir = folder/'trial-verification'; cache_dir.mkdir(exist_ok=True)
    reports, tested = [], set()
    mapper = engine.load('factor_free_pari_mapping.sage'); mapper.pari.allocatemem(256000000, silent=True)
    for stage in terminal['stages']:
        wd = folder/f"replay-M17/epoch-{stage['epoch']:02d}"
        vd = cache_dir/f"epoch-{stage['epoch']:02d}"; vd.mkdir(exist_ok=True)
        basis = point_tuple(state.basis); selection = read(wd/'selection.json')
        require(point_tuple(selection['basis']) == basis, 'replay subgroup chain differs')
        landscape_paths = [wd/'selection.json', folder/'protocol.json', folder/'orbits.tsv']
        landscape_paths += list(wd.glob('anchor-*-full.npz'))
        landscape_bindings = {str(p.relative_to(root)):sha(p) for p in landscape_paths}
        cached = vd/'landscape.json'; tested_record = sorted(map(list,tested))
        if cached.exists():
            record = read(cached)
            require(record['sources'] == ctx.sources and record['bindings'] == landscape_bindings and
                    record['previously_tested'] == tested_record, 'landscape replay cache differs')
            bindings(root, record['bindings']); metrics = record['metrics']
        else:
            metrics = verify_landscape(engine, model, basis, tested, wd, policy)
            atomic(cached, {'sources':ctx.sources,'bindings':landscape_bindings,'metrics':metrics,
                           'previously_tested':tested_record}, immutable=True)
        charts = indexed_paths(wd, expected=stage['charts'])
        require(bool(charts), 'completed epoch without charts')
        state = restore_state(read(wd/'cloud-000.json')['final_state'], model, basis)
        require(read(wd/'epoch-seed-state.json') == state.record(), 'epoch seed state binding differs')
        prefix, expected_points = [], list(basis)
        seen = {(x,abs(y)) for x,y in basis}
        for j, path in enumerate(charts):
            chart = read(path); check_chart(chart, selection, j)
            require(chart['search']['height_bound'] == policy['height'] and
                    chart['search']['timeout_seconds'] == policy['seconds_per_chart'] and
                    chart['search']['gp_binary_sha256'] == policy['gp_sha256'], 'frozen chart budget/backend changed')
            prefix.append(chart)
            for raw in chart['search']['finite_curve_points']:
                point = F(raw['x']), F(raw['y']); key = point[0], abs(point[1])
                if key not in seen: seen.add(key); expected_points.append(point)
            cloud_path, audit_path = wd/f'cloud-{j:03d}.json', wd/f'mod2-{j:03d}.json'
            snapshot, audit = read(cloud_path), read(audit_path)
            require(audit.get('status') == 'COMPLETE_DECLARED_FINITE_AUDIT', 'partial mod-2 audit')
            require(curve_tuple(snapshot['curve']) == model and
                    snapshot['final_state'] == state.record(), 'snapshot seed or model changed')
            require(snapshot['charts'] == prefix, f'cumulative snapshot differs at chart {j}')
            require(audit['input_sha256'] == sha(cloud_path) and
                    within(root,audit['input_path']) == cloud_path, 'audit snapshot binding differs')
            require(curve_tuple(audit['curve']) == model and
                    point_tuple(audit['points']) == tuple(expected_points), 'audit contains missing/unreturned points')
            chart_bindings = {str(p.relative_to(root)):sha(p) for p in (path,cloud_path,audit_path)}
            record_path = vd/f'chart-{j:06d}.json'
            if record_path.exists():
                old = read(record_path)
                require(old['sources'] == ctx.sources and old['bindings'] == chart_bindings and
                        old['selection_sha256'] == sha(wd/'selection.json'), 'chart replay cache differs')
            else:
                require(mapper.mapping(model,basis,chart['centre']) == chart['mapping'], 'quartic map differs')
                search = PointedQuarticSearch(state=state, centre={'coefficients':chart['centre']['representative']},
                                              coordinate_policy=chart['mapping']['coordinate_policy'])
                backend.replay(search,chart['mapping'],chart['search']); mod2.check(audit_path)
                atomic(record_path, {'sources':ctx.sources,'bindings':chart_bindings,
                       'selection_sha256':sha(wd/'selection.json')}, immutable=True)
            if j < len(charts)-1:
                require(audit['rank_lower_bound'] == len(basis), 'stale chart executed after gain')
            tested.add(tuple(chart['centre']['point']))
            if j % 25 == 0 or j+1 == len(charts):
                print('REPLAY_CHART',ctx.case,stage['epoch'],j+1,'/',len(charts),flush=True)
        last = read(wd/stage['audit']); enlarged = point_tuple(last['independent_points'])
        require(enlarged[:len(basis)] == basis and len(enlarged) == stage['after'], 'certified prefix/rank differs')
        state = certified_state(model,enlarged,last['rank_certificate'])
        odd_data = read(wd/'modl.json')
        require(odd_data['status'] == 'COMPLETE_BOUNDED_QUOTIENT_AUDIT' and len(odd_data['audits']) == 2 and
                all(a['status'] == 'COMPLETE_BOUNDED_QUOTIENT_AUDIT' for a in odd_data['audits']), 'partial odd-prime audit')
        require(odd_data['input_sha256'] == sha(wd/stage['audit']) and
                odd_data['points'] == last['points'] and curve_tuple(odd_data['curve']) == model,
                'odd-prime audit not bound to actual returned cloud')
        modl.check(wd/'modl.json')
        odd = {str(a['modulus']):a['finite_column_rank'] for a in odd_data['audits']}
        require(odd == {'3':state.rank,'5':state.rank}, 'mod-3/5 replay disagrees')
        checkpoint_hashes = {str(p.relative_to(root)):sha(p) for p in wd.iterdir() if p.suffix in ('.json','.npz')}
        reports.append({'epoch':stage['epoch'],'before':len(basis),'after':state.rank,
                        'charts_replayed':len(charts),'independent_modl_ranks':odd,'landscape':metrics,
                        'checkpoint_hashes':checkpoint_hashes})
        atomic(cache_dir/'progress.json', {'sources':ctx.sources,'stages':reports})
    require(state.rank == terminal['final_rank_lower_bound'], 'terminal rank differs')
    result = {'schema':'det1092-funnel-independent-replay.v1','case':ctx.case,'sources':ctx.sources,
              'protocol_sha256':sha(folder/'protocol.json'),'terminal_sha256':sha(folder/'replay-M17/terminal.json'),
              'rank_lower_bound':state.rank,'stages':reports,'charts':terminal['charts'],'claim_boundary':policy['scope']}
    atomic(folder/'trial-replay.json',result,immutable=True)
    final_paths = [folder/'trial-replay.json',folder/'protocol.json',folder/'seed-input.json',
                   folder/'seed-proof.json',folder/'replay-M17/terminal.json']
    verified = {'schema':'det1092-funnel-result.v1','status':'PASS_INDEPENDENT_DET1092_REPLAY',
        'case':ctx.case,'parameter':read(folder/'seed-input.json')['parameter'],'initial_rank':18,
        'rank_lower_bound':state.rank,'gain':state.rank-18,'stop_reason':terminal['stop_reason'],
        'charts':terminal['charts'],'sources':ctx.sources,
        'bindings':{str(p.relative_to(root)):sha(p) for p in final_paths},'claim_boundary':policy['scope']}
    atomic(folder/'trial-verified.json',verified,immutable=True)
    if verified['gain'] > 0:
        atomic(folder/'trial-discovery.json', {'curve':list(map(str,model)),
            'points':[list(map(str,p)) for p in state.basis],'rank_lower_bound':state.rank,
            'parameter':verified['parameter'],'verified_sha256':sha(folder/'trial-verified.json'),
            'status':'VERIFIED_GAIN_NOT_A_NOVELTY_OR_RECORD_CLAIM'}, immutable=True)
    atomic(folder/'queue.json', {'certified_rank_lower_bound':state.rank,
        'queue':funnel.queue(state.rank),'verified_sha256':sha(folder/'trial-verified.json'),
        'extra_allocation_seconds':0,'scope':'Queue classification; no implicit extra run.'}, immutable=True)
    print('INDEPENDENT_DET1092_REPLAY_COMPLETE',ctx.case,state.rank,flush=True)
    return verified
