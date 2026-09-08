"""Explicit search/replay services for the warm runner; no runtime source patches.

Only this module (inside a Sage subprocess) imports the numerical stack. The
frozen V3 selector and chart backend are reused unchanged. A completed old
terminal is replayed, NEVER searched again. No result is called a new curve or
record without a separate catalogue/novelty check.
"""
from __future__ import annotations

from dataclasses import dataclass
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
from pathlib import Path
from fractions import Fraction as F
import contextlib
import csv
import os
import sys
import time
import uuid

from v3_warm_support import (WARM, atomic, assert_basis, bindings, check_chart,
    curve_tuple, indexed_paths, point_tuple, read, require, sha, terminal_structure)

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
D = ROOT/'artifacts/local/elliptic-curves/v3-transfer-11952-v3'


def own_sources():
    names = ('v3_warm_engine.py', 'v3_warm_replay.py', 'v3_warm_support.py', 'run_v3_warm_start_overnight.py')
    return {str((CAS/n).relative_to(ROOT)): sha(CAS/n) for n in names}


def load(name, path):
    loader = SourceFileLoader(name, str(path))
    spec = spec_from_loader(name, loader)
    module = module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


def certified_state(model, points, proof):
    """Actually verify finite images and torsion, not only JSON rank/cardinality."""
    from research_runtime.arithmetic import ArithmeticContext, CurveModel
    from research_runtime.mw_state import MWState
    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    model, points = curve_tuple(model), point_tuple(points)
    primes = tuple(int(r['prime']) for r in proof['signatures'])
    require(bool(primes) and len(set(primes)) == len(primes), 'invalid seed certificate primes')
    require(proof['rank_lower_bound'] == len(points), 'certificate does not describe supplied basis')
    cache = Cache(MemoryFactStore())
    state = MWState.empty(ArithmeticContext.for_search(CurveModel(model)), cache=cache, primes=primes,
                          no_two_torsion_prime=int(proof['no_rational_2_torsion_prime']))
    for i, point in enumerate(points):
        state = state.adjoin(point, cache=cache, extra_primes=())
        require(state.rank == i+1, f'finite independence failed at seed column {i}')
    assert_basis(state, model, points, len(points))
    return state


def restore_state(record, model, points):
    from research_runtime.mw_state import MWState
    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    state = MWState.from_record(record, cache=Cache(MemoryFactStore()))
    assert_basis(state, model, points, len(points))
    return state


@dataclass
class Context:
    case: str
    folder: Path
    root: Path
    policy: dict
    engine: object
    model: tuple
    state: object
    sources: dict


def setup(case):
    from sage.all import pari
    import sage.version
    require(case in WARM, 'case not in the authorized warm roster')
    control = read(D/'control-native11952/verified.json')
    require(control['status'] == 'PASS_INDEPENDENT_TRANSFER_REPLAY' and
            control['initial_rank'] == control['rank_lower_bound'] == 17 and control['gain'] == 0,
            'expected completed no-gain control, not a re-labelled success')
    bindings(ROOT, control['bindings'])
    roster = read(D/'roster.json')
    for filename, field in (('parent-bank.json','parent_bank_sha256'), ('gate.json','gate_sha256'),
                            ('preparation.json','preparation_sha256')):
        require(sha(D/filename) == roster[field], f'changed {filename}')
    entry = next(r for r in roster['cases'] if r['id'] == case)
    folder = D/case
    require(sha(folder/'protocol.json') == entry['protocol_sha256'], 'warm protocol changed')
    policy = read(folder/'protocol.json')
    bindings(ROOT, policy['inputs'])
    bindings(ROOT, policy['sources'])
    bindings(ROOT, policy['driver_sources'])
    require(policy['initial_rank'] == entry['initial_rank'] == 27, 'warm intake must start at certified 27')
    software = policy['software']
    require(software['sage'] == sage.version.version and software['pari'] == str(pari.version()),
            'Sage/PARI versions differ from frozen preparation')
    require(sha(Path('/usr/bin/gp')) == policy['gp_sha256'], 'PARI executable differs from frozen search')
    seed, proof = read(folder/'seed-input.json'), read(folder/'seed-proof.json')
    model, points = curve_tuple(seed['curve']), point_tuple(seed['points'])
    state = certified_state(model, points, proof)
    assert_basis(state, model, points, policy['initial_rank'])
    # Load the numerical engine directly, not the v1->v2->v3 transfer wrapper
    # chain. Its source bytes remain exactly those in the frozen protocol.
    engine = load('warm_frozen_numeric', CAS/'adaptive_visibility_cascade_v3.sage')
    engine.D, engine.ORBITS = folder, folder/'orbits.tsv'
    engine.v1.D, engine.v1.ORBITS = folder, engine.ORBITS
    require(engine.sources() == policy['sources'], 'numerical engine source set differs')
    sources = {**policy['sources'], **policy['driver_sources'], **own_sources()}
    engine.v1.guard()
    try:
        (D/'control-native11952/seed-input.json').open()
    except PermissionError:
        pass
    else:
        raise ValueError('worker artifact guard failed to isolate other-case input')
    return Context(case, folder, ROOT, policy, engine, model, state, sources)


def preflight(ctx):
    """Seed -> exact quartic -> backend map validation, without point search."""
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    with ctx.engine.ORBITS.open() as stream:
        row = next(csv.DictReader(stream, delimiter='\t'))
    word = list(map(int, row['parent_MW17_w'].split())) + [0]*(ctx.state.rank-17)
    centre = {'representative': word}
    mapper = ctx.engine.load('factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    mapping = mapper.mapping(ctx.model, point_tuple(ctx.state.basis), centre)
    search = PointedQuarticSearch(state=ctx.state, centre={'coefficients': word},
                                 coordinate_policy=mapping['coordinate_policy'])
    backend.validate_map(search, mapping)
    print('PREFLIGHT_PASS', ctx.case, 'seed_rank', ctx.state.rank, 'exact_chart_map', flush=True)


def _audit(snapshot, path):
    import audit_recorded_point_mod2_rank_v3 as mod2
    if path.exists() and read(path).get('status') != 'COMPLETE_DECLARED_FINITE_AUDIT':
        saved = path.parent/'incomplete-audits'/str(uuid.uuid4())
        saved.mkdir(parents=True)
        path.replace(saved/path.name)
    with (path.parent/(path.stem+'.check.log')).open('a') as log, contextlib.redirect_stdout(log):
        if not path.exists():
            mod2.build(snapshot, path, 1000, sha(snapshot))
        mod2.check(path)
    result = read(path)
    require(result['input_sha256'] == sha(snapshot), 'audit snapshot binding differs')
    return result


def _odd(audit, path):
    import audit_retained_cloud_modl as modl
    if path.exists() and read(path).get('status') != 'COMPLETE_BOUNDED_QUOTIENT_AUDIT':
        archived = path.parent/'incomplete-audits'/str(uuid.uuid4())
        archived.mkdir(parents=True)
        path.replace(archived/path.name)
    # The checker is always called, even on a no-gain/censored epoch.
    with (path.parent/'odd-check.log').open('a') as log, contextlib.redirect_stdout(log):
        if not path.exists():
            modl.build(audit, path)
        modl.check(path)
    data = read(path)
    require(data['input_sha256'] == sha(audit) and data['points'] == read(audit)['points'],
            'odd-prime audit is bound to a different point cloud')
    require(data['status'] == 'COMPLETE_BOUNDED_QUOTIENT_AUDIT' and len(data['audits']) == 2 and
            all(a['status'] == 'COMPLETE_BOUNDED_QUOTIENT_AUDIT' for a in data['audits']), 'partial odd audit')
    return {str(a['modulus']): a['finite_column_rank'] for a in data['audits']}


def run_search(ctx):
    from pointed_quartic_search import PointedQuarticSearch
    from v3_warm_replay import verify_landscape
    import pari_pointed_backend as backend
    engine, model, state = ctx.engine, ctx.model, ctx.state
    folder, policy = ctx.folder, ctx.policy
    out = folder/'replay-M17'
    out.mkdir(exist_ok=True)
    if (out/'terminal.json').exists():
        terminal_structure(folder)
        print('SEALED_SEARCH_REUSED', ctx.case, 'no point search rerun', flush=True)
        return
    mapper = engine.load('factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    total, stages, tested = 0, [], set()
    reason = 'EPOCH_BUDGET_EXHAUSTED'
    for epoch in range(policy['max_epochs']):
        if total >= policy['max_charts']:
            reason = 'CHART_BUDGET_EXHAUSTED'
            break
        wd = out/f'epoch-{epoch:02d}'
        wd.mkdir(exist_ok=True)
        basis, before = point_tuple(state.basis), state.rank
        started = time.monotonic()
        existed = (wd/'selection.json').exists()
        existing_charts = indexed_paths(wd)
        if existed:
            selection = read(wd/'selection.json')
            require(point_tuple(selection['basis']) == basis, 'resume basis differs from sealed schedule')
            # A resumed selection is replayed, not trusted or reranked.
            verify_landscape(engine, model, basis, tested, wd, policy)
            centres = selection['centres']
        else:
            require(not existing_charts and not (wd/'stage.json').exists(), 'charts/stage without a selection')
            unfinished = list(wd.glob('anchor-*'))
            if unfinished:
                archived = wd/'incomplete-landscapes'/str(uuid.uuid4())
                archived.mkdir(parents=True)
                for path in unfinished:
                    path.replace(archived/path.name)
                print('PRESERVED_PARTIAL_LANDSCAPE', archived, flush=True)
            centres = engine.landscape(model, basis, tested, wd, policy)
            selection = read(wd/'selection.json')
        require(bool(centres), 'no new finalists; no search coverage is claimed for this engineering stop')
        require(len(existing_charts) <= len(centres), 'extra chart checkpoints outside frozen schedule')
        if (wd/'epoch-seed-state.json').exists():
            state = restore_state(read(wd/'epoch-seed-state.json'), model, basis)
        elif existing_charts and (wd/'cloud-000.json').exists():
            state = restore_state(read(wd/'cloud-000.json')['final_state'], model, basis)
        else:
            # Commit the exact state BEFORE any chart; a crash between chart
            # publication and cloud publication must not lose its binding.
            atomic(wd/'epoch-seed-state.json', state.record(), immutable=True)
        stage_seed = state
        charts, censored, audit_path = [], 0, None
        for j, centre in enumerate(centres):
            if total >= policy['max_charts']:
                break
            path = existing_charts[j] if j < len(existing_charts) else wd/f'chart-{j:03d}.json'
            if path.exists():
                chart = read(path)
                check_chart(chart, selection, j)
                search = PointedQuarticSearch(state=stage_seed, centre={'coefficients': centre['representative']},
                                             coordinate_policy=chart['mapping']['coordinate_policy'])
                backend.replay(search, chart['mapping'], chart['search'])
            else:
                mapping = mapper.mapping(model, basis, centre)
                search = PointedQuarticSearch(state=stage_seed, centre={'coefficients': centre['representative']},
                                             coordinate_policy=mapping['coordinate_policy'])
                transcript, points = backend.execute(search, mapping, policy['height'],
                                                     policy['seconds_per_chart'], policy['gp_sha256'])
                require(backend.replay(search, mapping, transcript) == points, 'returned point replay changed')
                chart = {'index': j, 'centre': centre, 'mapping': mapping, 'search': transcript}
                atomic(path, chart, immutable=True)
            charts.append(chart)
            tested.add(tuple(centre['point']))
            total += 1
            censored += chart['search']['status'] != 'bounded_search_complete'
            snapshot, audit_path = wd/f'cloud-{j:03d}.json', wd/f'mod2-{j:03d}.json'
            if snapshot.exists():
                saved = read(snapshot)
                require(saved['charts'] == charts and curve_tuple(saved['curve']) == model and
                        point_tuple(saved['final_state']['state']['reductions']['points']) == basis,
                        'existing cumulative point snapshot differs')
            else:
                atomic(snapshot, {'status': 'INCREMENTAL_RETAINED_CLOUD', 'family': '11952-warm-transfer',
                       'parameter': ctx.case, 'curve': list(map(str, model)), 'charts': charts,
                       'final_state': stage_seed.record(), 'rank_lower_bound': before}, immutable=True)
            cloud = _audit(snapshot, audit_path)
            print(ctx.case, 'epoch', epoch, 'chart', j+1, '/', len(centres),
                  'certified', cloud['rank_lower_bound'], flush=True)
            if cloud['rank_lower_bound'] > before:
                enlarged = point_tuple(cloud['independent_points'])
                require(enlarged[:before] == basis, 'independent generic-first prefix lost')
                state = certified_state(model, enlarged, cloud['rank_certificate'])
                require(len(existing_charts) <= j+1, 'retained stale charts after rank gain')
                break
        require(audit_path is not None, 'no chart audit in epoch')
        odd = _odd(audit_path, wd/'modl.json')
        require(odd == {'3': state.rank, '5': state.rank}, 'odd-prime rank disagreement; not a no-gain result')
        if state.rank >= policy['target_rank']:
            reason = 'TARGET_LOWER_BOUND_REACHED'
        elif total >= policy['max_charts']:
            reason = 'CHART_BUDGET_EXHAUSTED'
        elif censored:
            reason = 'CENSORED_SEARCH'
        elif state.rank > before:
            reason = 'REBUILD_AFTER_CERTIFIED_GAIN'
        else:
            reason = 'COMPLETE_FINITE_NO_GAIN' if len(charts) == len(centres) else 'INCOMPLETE_EPOCH'
        stage = {'epoch': epoch, 'before': before, 'after': state.rank, 'charts': len(charts),
                 'stale_charts_cancelled': len(centres)-len(charts) if state.rank > before else 0,
                 'audit': audit_path.name, 'audit_sha256': sha(audit_path),
                 'full_cosets_scored': selection['full_cosets_scored'], 'stop_reason': reason,
                 'censored_charts': censored, 'wall_seconds': time.monotonic()-started}
        if (wd/'stage.json').exists():
            old = read(wd/'stage.json')
            require({k:v for k,v in old.items() if k != 'wall_seconds'} ==
                    {k:v for k,v in stage.items() if k != 'wall_seconds'}, 'sealed stage differs on resume')
            stage = old
        else:
            atomic(wd/'stage.json', stage, immutable=True)
        stages.append(stage)
        atomic(out/'stages.json', stages)
        if reason != 'REBUILD_AFTER_CERTIFIED_GAIN':
            break
    else:
        reason = 'EPOCH_BUDGET_EXHAUSTED'
    atomic(out/'terminal.json', {'status': 'TERMINAL_BOUNDED_TRANSFER', 'stop_reason': reason,
           'initial_rank': policy['initial_rank'], 'final_rank_lower_bound': state.rank,
           'stages': stages, 'charts': total, 'protocol_sha256': sha(folder/'protocol.json'),
           'read_paths': sorted(engine.v1.READS), 'execution_sources': ctx.sources,
           'scope': 'Search terminal only. Independent replay required; no upper bound or record claim.'}, immutable=True)
    terminal_structure(folder)
    print('WARM_SEARCH_TERMINAL', ctx.case, state.rank, reason, flush=True)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['preflight', 'search', 'replay'])
    parser.add_argument('--case', required=True, choices=WARM)
    parser.add_argument('--session', type=Path, required=True)
    args = parser.parse_args()
    require(os.environ.get('V3_WARM_SUPERVISED') == '1', 'use the supervised warm runner')
    session = read(args.session)
    require(session['sources'] == own_sources(), 'controller/worker implementation snapshot mismatch')
    bindings(ROOT, session['sources'])
    ctx = setup(args.case)
    if args.action == 'preflight':
        preflight(ctx)
    elif args.action == 'search':
        run_search(ctx)
    else:
        from v3_warm_replay import replay_case
        replay_case(ctx)


if __name__ == '__main__':
    main()
