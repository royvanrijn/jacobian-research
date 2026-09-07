#!/usr/bin/env sage-python
"""First bounded cross-parent V3 transfer: one control, three warm rank-27 fibres.

Existing V1/V2/V3 sources and evidence are never edited. The V3 landscape and
selector are imported, not copied or retuned. This driver supplies a separately
certified det-948 orbit bank and exact generic-first seed packages. Commands:
  prepare --v3-replay PATH
  next
  status
Run from research/, using `sage -python`. `next` executes at most ONE case and
its independent replay through the shared process-group supervisor.
"""
from __future__ import annotations
import argparse
import contextlib
import csv
import fcntl
import io
import json
import os
import sys
import time
from fractions import Fraction as F
from importlib.machinery import SourceFileLoader
from pathlib import Path

from v3_transfer_contract import (BASE_COMMIT, CASE_SPECS, LIMITS, require, sha,
    read_json, resolve_json, write_new, within, check_bindings, on_curve,
    extract_sage_literals, transport_to_short, check_orbit_rows,
    validate_gate, classify_stop)

ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT/'elliptic-curves/cas'
ART = ROOT/'artifacts/generated-results/elliptic-curves'
D = ROOT/'artifacts/local/elliptic-curves/v3-transfer-11952-v1'
V3_D = ROOT/'artifacts/local/elliptic-curves/adaptive-visibility-cascade-v3'
ATLAS = ART/'compact_six_r17_atlas_v1.json'
CONTROL = ROOT/'elliptic-curves/data/half_lattice_rank29_control_inputs_v1.json'


def module(name, path):
    return SourceFileLoader(name, str(path)).load_module()


def driver_sources():
    files = [Path(__file__), CAS/'v3_transfer_contract.py', CAS/'v3_transfer_orbits.py',
             CAS/'exact_parity_ellipsoid.py', CAS/'check_visibility_cascade_v3.sage',
             CAS/'check_visibility_metric_v3.sage', CAS/'research_runtime/supervisor.py',
             CAS/'export_compact_r17_atlas.sage']
    return {str(p.relative_to(ROOT)): sha(p) for p in files}


def engine_for(folder):
    engine = module('transfer_frozen_v3', CAS/'adaptive_visibility_cascade_v3.sage')
    engine.D = folder
    engine.ORBITS = folder/'orbits.tsv'
    engine.v1.D = folder
    engine.v1.ORBITS = engine.ORBITS
    # The legacy replay path is called replay-M17 even for a warm input.
    # The actual seed rank is explicit and never inferred from that folder name.
    def seed(_row):
        s = read_json(folder/'seed-input.json')
        return tuple(map(F, s['curve'])), tuple(tuple(map(F, p)) for p in s['points'])
    engine.v1.seed = seed
    return engine


def prepare(v3_replay):
    """No point searches. All four inputs and policy bytes freeze together."""
    require(not D.exists(), 'Preserve existing preparation, including failed preparation. No automatic retry/expansion.')
    gate = validate_gate(ROOT, V3_D, v3_replay, V3_D/'metric-replay-M17.json')
    from sage.all import ZZ, QQ, matrix, PolynomialRing, pari
    import sage.version as sage_version
    from v3_transfer_orbits import shell_bank
    import certify_compact_r17_candidates as cert
    from research_runtime.search_state import raw_state
    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    from research_runtime.store import atomic_write

    D.mkdir(parents=True)
    write_new(D/'gate.json', gate)
    write_new(D/'preparation.json', {'base_commit': BASE_COMMIT, 'limits': LIMITS,
        'sources': driver_sources(), 'case_specs': list(CASE_SPECS),
        'scope': 'Predeclared det948 transfer. Only a completed positive-control gain releases warm jobs. No new parameter search or automatic submission.'})
    # Check each compact generic section and its source transport before use.
    atlas_path = resolve_json(ATLAS)
    atlas = read_json(atlas_path)
    family = next(f for f in atlas['families'] if f['family'] == '11952')
    atlas_helper = module('transfer_atlas', CAS/'export_compact_r17_atlas.sage')
    sections, gram = atlas_helper.transport(family)
    require(sections == family['sections'] and gram == family['generic_height_gram'], 'Generic section/Gram replay failed')
    G = matrix(ZZ, gram)
    require(G.nrows() == 17 and G.det() == 948 and G.is_positive_definite(), 'Wrong parent lattice')
    require(all(G[i,i] % 2 == 0 for i in range(17)), 'Non-even generic lattice')
    U = matrix(ZZ, pari(G).qflllgram()).transpose()
    require(abs(U.det()) == 1, 'Nonunimodular lattice adapter')
    H = U*G*U.transpose()
    bank = shell_bank([list(map(int,r)) for r in H.rows()],
                      [list(map(int,r)) for r in U.rows()], node_limit=LIMITS['orbit_node_limit'])
    require(bank['norm_counts'].get(2, 0) == 0, 'Generic lattice is not rootless')
    check_orbit_rows(gram, bank['rows'])
    bank.update(generic_gram=gram, LLL=[list(map(int,r)) for r in U.rows()],
                reduced_gram=[list(map(int,r)) for r in H.rows()], family='11952', determinant=948)
    write_new(D/'parent-bank.json', bank)
    buf = io.StringIO()
    w = csv.writer(buf, delimiter='\t', lineterminator='\n')
    w.writerow(['orbit_mask', 'minimum_norm', 'parent_MW17_w'])
    for row in bank['rows']:
        w.writerow([row['mask'], row['norm'], ' '.join(map(str,row['word']))])
    table = buf.getvalue()
    boundary_path = resolve_json(CONTROL)
    boundary = read_json(boundary_path)
    require(boundary['boundary']['output_contains_exceptional_point_coordinates'] is False, 'Control input is not redacted')
    control = next(r for r in boundary['cases'] if r['label'] == 'curve12-2024-rank29')
    require(control['generic_height_gram'] == gram, 'Native control uses a different ordered generic basis')
    ring = PolynomialRing(QQ, 't')
    policy0 = read_json(V3_D/'protocol.json')
    original_engine = engine_for(D/'unused-engine-binding')
    require(policy0['sources'] == original_engine.sources(), 'V3 sources changed since the completed calibration')
    roster = []
    for spec in CASE_SPECS:
        provenance = {str(atlas_path.relative_to(ROOT)): sha(atlas_path)}
        if spec['kind'] == 'positive-control':
            source = list(map(F, control['short_model']))
            model = (F(0), F(0), F(0), source[3]/6**4, source[4]/6**6)
            generic = tuple((F(x)/36, F(y)/216) for x,y in control['generic_points'])
            supplied, transport = generic, None
            provenance[str(boundary_path.relative_to(ROOT))] = sha(boundary_path)
        else:
            t = QQ(spec['parameter']); q = t.denominator()
            A = ring(family['A_coefficients_low_to_high'])(t)*q**8
            B = ring(family['B_coefficients_low_to_high'])(t)*q**12
            model = (F(0), F(0), F(0), F(str(A)), F(str(B)))
            def ev(rec):
                return ring(rec['numerator_coefficients_low_to_high'])(t)/ring(rec['denominator_coefficients_low_to_high'])(t)
            generic = tuple((F(str(ev(s['X'])*q**4)), F(str(ev(s['Y'])*q**6))) for s in sections)
            exported = ART/spec['export']
            old_model, old_points = extract_sage_literals(exported.read_text())
            require(len(old_points) == 27, 'Unexpected retained export size')
            extras, transport = transport_to_short(old_model, model, old_points)
            supplied = generic+extras
            provenance[str(exported.relative_to(ROOT))] = sha(exported)
        require(len(generic) == 17 and all(on_curve(model,p) for p in supplied), 'Bad specialization or point input')
        cache = Cache(MemoryFactStore())
        state = raw_state(model, supplied, cache=cache, prime_bound=1000)
        basis = tuple(tuple(map(F,p)) for p in state.basis)
        require(basis[:17] == generic and len(basis) == spec['initial_rank'], 'Generic prefix/rank intake failed; do not adjust the target silently')
        proof = cert.checked_rank(model, basis)
        folder = D/spec['id']; folder.mkdir()
        atomic_write(folder/'orbits.tsv', table.encode())
        seed = {'curve': list(map(str,model)), 'points': [list(map(str,p)) for p in basis],
                'generic_rank': 17, 'initial_rank': len(basis), 'kind': spec['kind'],
                'family': '11952', 'generic_height_gram': gram}
        write_new(folder/'seed-input.json', seed)
        write_new(folder/'seed-proof.json', proof)
        # These preparation records are outside solver inputs.
        write_new(D/(spec['id']+'-provenance.json'), {'spec': spec, 'inputs': provenance,
            'exact_model_transport': transport, 'supplied_points': [list(map(str,p)) for p in supplied]})
        engine = engine_for(folder)
        policy = {**policy0, **{k:LIMITS[k] for k in ('max_charts','max_epochs','target_rank')}}
        policy.update(schema='v3-transfer-11952.v1', sources=engine.sources(), driver_sources=driver_sources(),
            software={'sage':sage_version.version,'pari':str(pari.version()),'python':sys.version},
            calibration_only=spec['kind']=='positive-control',
            inputs={str(p.relative_to(ROOT)):sha(p) for p in (folder/'seed-input.json',folder/'orbits.tsv',folder/'seed-proof.json')},
            family='11952', generic_determinant=948, initial_rank=len(basis),
            oracle_boundary='Known-positive control or explicitly warm-started inventory. Selection cannot read exceptional/public targets, post-hoc diagnostics, catalogue ranks, or other jobs. Not an untouched prospective population.',
            invariance='V3 selector unchanged. Full finite extension coverage does NOT establish LLL/Babai, anchor-bank or full-policy basis invariance.',
            scope='Only this four-case transfer is authorized. A bounded miss is not a rank upper bound; no publication or full rollout is automatic.')
        write_new(folder/'protocol.json', policy)
        roster.append({'id':spec['id'], 'kind':spec['kind'], 'initial_rank':len(basis),
                       'protocol_sha256':sha(folder/'protocol.json')})
    write_new(D/'roster.json', {'status':'FROZEN_TRANSFER_READY', 'cases':roster,
        'parent_bank_sha256':sha(D/'parent-bank.json'), 'gate_sha256':sha(D/'gate.json'),
        'preparation_sha256':sha(D/'preparation.json')})
    print('PREPARED: one redacted native11952 control and three independent warm rank27 inputs. No search has run.', flush=True)


def bind_job(case):
    require(case in {s['id'] for s in CASE_SPECS}, 'Unknown case')
    roster = read_json(D/'roster.json')
    require(sha(D/'parent-bank.json') == roster['parent_bank_sha256'], 'Parent bank changed')
    require(sha(D/'gate.json') == roster['gate_sha256'], 'Calibration gate changed')
    require(sha(D/'preparation.json') == roster['preparation_sha256'], 'Preparation changed')
    entry = next(r for r in roster['cases'] if r['id'] == case)
    folder = D/case
    require(sha(folder/'protocol.json') == entry['protocol_sha256'], 'Job protocol changed')
    engine = engine_for(folder)
    policy = read_json(folder/'protocol.json')
    require(policy['sources'] == engine.sources(), 'Frozen V3 engine changed')
    require(policy['driver_sources'] == driver_sources(), 'Transfer implementation changed')
    check_bindings(ROOT, policy['inputs'])
    if entry['kind'] != 'positive-control':
        control = read_json(D/'control-native11952/verified.json')
        require(control.get('status') == 'PASS_INDEPENDENT_TRANSFER_REPLAY' and control.get('initial_rank') == 17, 'Control has not been independently replayed')
        replay = read_json(D/'control-native11952/full-replay.json')
        require(control['rank_lower_bound'] == replay['rank_lower_bound'], 'Control report/replay rank mismatch')
        require(control['rank_lower_bound'] > 17, 'Warm stage blocked: no independently verified positive-control gain')
        require(control['stop_reason'] not in ('CENSORED_SEARCH','INCOMPLETE_EPOCH'), 'Control censored; review before transfer')
        check_bindings(ROOT, control['bindings'])
    return engine, folder, policy


def run_case(case):
    engine, folder, p = bind_job(case)
    model, basis = engine.v1.seed(None)
    # Install the same artifact-read guard AFTER gate/own-seed validation.
    engine.v1.guard()
    import audit_recorded_point_mod2_rank_v3 as mod2
    import audit_retained_cloud_modl as modl
    import pari_pointed_backend as backend
    from research_runtime.search_state import raw_state
    from research_runtime.memory_store import MemoryFactStore
    from research_runtime.quotient_only_reduction import QuotientOnlyReductionCache as Cache
    from pointed_quartic_search import PointedQuarticSearch
    from research_runtime.store import checkpoint
    cache = Cache(MemoryFactStore())
    state = raw_state(model,basis,cache=cache,prime_bound=1000)
    require(state.rank == p['initial_rank'] and tuple(state.basis) == basis, 'Seed rank replay failed')
    out = folder/'replay-M17'; out.mkdir(exist_ok=True)
    require(not (out/'terminal.json').exists(), 'Case already terminal')
    mapper = engine.load('factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000,silent=True)
    tested, stages, total, reason = set(), [], 0, 'EPOCH_BUDGET_EXHAUSTED'
    for epoch in range(p['max_epochs']):
        wd = out/f'epoch-{epoch:02d}'; wd.mkdir(exist_ok=True)
        basis = tuple(tuple(map(F,q)) for q in state.basis); before = len(basis)
        # Resume only an already sealed epoch; an unfinished job needs explicit review.
        if (wd/'stage.json').exists():
            stage = read_json(wd/'stage.json'); mod2.check(wd/stage['audit']); modl.check(wd/'modl.json')
            require(stage['before'] == before, 'Resume prefix mismatch')
            cloud = read_json(wd/stage['audit'])
            require(tuple(tuple(map(F,q)) for q in cloud['independent_points'][:before]) == basis, 'Resume basis changed')
            for path in sorted(wd.glob('chart-*.json')):
                tested.add(tuple(read_json(path)['centre']['point'])); total += 1
            state = raw_state(model,tuple(tuple(map(F,q)) for q in cloud['independent_points']),cache=cache,prime_bound=1000)
            stages.append(stage); reason = stage['stop_reason']
            if reason != 'REBUILD_AFTER_CERTIFIED_GAIN': break
            continue
        require(not any(wd.iterdir()), 'Unsealed partial epoch retained. Do not overwrite or auto-expand its budget.')
        started = time.monotonic()
        centres = engine.landscape(model,basis,tested,wd,p)
        require(bool(centres), 'No new finalists: retain as an engineering stop; do not create a partially replayable terminal')
        charts, censored, audit = [], 0, None
        for j,c in enumerate(centres):
            if total >= p['max_charts']: break
            mapping = mapper.mapping(model,basis,c)
            search = PointedQuarticSearch(state=state,centre={'coefficients':c['representative']},coordinate_policy=mapping['coordinate_policy'])
            transcript, points = backend.execute(search,mapping,p['height'],p['seconds_per_chart'],p['gp_sha256'])
            require(backend.replay(search,mapping,transcript) == points, 'Exact map replay changed')
            chart = {'index':j,'centre':c,'mapping':mapping,'search':transcript}
            write_new(wd/f'chart-{j:03d}.json',chart)
            charts.append(chart); total += 1; tested.add(tuple(c['point']))
            censored += int(transcript['status'] != 'bounded_search_complete')
            snapshot, audit = wd/f'cloud-{j:03d}.json', wd/f'mod2-{j:03d}.json'
            write_new(snapshot, {'status':'INCREMENTAL_RETAINED_CLOUD','family':'11952-v3-transfer',
                'parameter':case,'curve':list(map(str,model)),'charts':charts,
                'final_state':state.record(),'rank_lower_bound':before})
            with (wd/f'audit-{j:03d}.log').open('x') as log, contextlib.redirect_stdout(log):
                mod2.build(snapshot,audit,1000,sha(snapshot)); mod2.check(audit)
            cloud = read_json(audit)
            print(case,'epoch',epoch,'chart',j+1,'/',len(centres),'certified',cloud['rank_lower_bound'],flush=True)
            if cloud['rank_lower_bound'] > before:
                enlarged = tuple(tuple(map(F,q)) for q in cloud['independent_points'])
                require(enlarged[:before] == basis, 'Certified prefix lost')
                state = raw_state(model,enlarged,cache=cache,prime_bound=1000)
                require(state.rank == len(enlarged), 'Admission/certificate disagreement')
                break
        if audit is None:
            reason = 'CHART_BUDGET_EXHAUSTED' if total >= p['max_charts'] else 'NO_NEW_FINALISTS'
            break
        # Always include terminal no-gain/censored clouds in the odd-prime checks.
        with (wd/'odd-audit.log').open('x') as log, contextlib.redirect_stdout(log):
            modl.build(audit,wd/'modl.json'); modl.check(wd/'modl.json')
        odd = read_json(wd/'modl.json')
        require(all(a['finite_column_rank'] == state.rank for a in odd['audits']), 'Odd-prime cloud disagrees: preserve for rank-admission review, not a no-gain result')
        reason = classify_stop(rank=state.rank,target=p['target_rank'],charts=total,max_charts=p['max_charts'],
            executed=len(charts),scheduled=len(centres),censored=censored,gain=state.rank>before)
        stage = {'epoch':epoch,'before':before,'after':state.rank,'charts':len(charts),
            'stale_charts_cancelled':len(centres)-len(charts) if state.rank>before else 0,
            'audit':audit.name,'audit_sha256':sha(audit),'full_cosets_scored':read_json(wd/'selection.json')['full_cosets_scored'],
            'stop_reason':reason,'censored_charts':censored,'wall_seconds':time.monotonic()-started}
        write_new(wd/'stage.json',stage); stages.append(stage)
        checkpoint(out/'stages.json',stages)
        if reason != 'REBUILD_AFTER_CERTIFIED_GAIN': break
    write_new(out/'terminal.json', {'status':'TERMINAL_BOUNDED_TRANSFER','stop_reason':reason,
        'initial_rank':p['initial_rank'],'final_rank_lower_bound':state.rank,'stages':stages,'charts':total,
        'protocol_sha256':sha(folder/'protocol.json'),'read_paths':sorted(engine.v1.READS),
        'scope':'New lower bounds require independent replay; no exact-rank, novelty, saturation-at-odd-primes, or general sensitivity claim.'})


def replay_case(case):
    engine, folder, p = bind_job(case)
    terminal = read_json(folder/'replay-M17/terminal.json')
    require(bool(terminal['stages']), 'No completed search epoch to replay')
    # Reuse the separately implemented upstream checker, with only its input
    # adapter changed. It recomputes all full scores, exact CVP, shortlists,
    # maps, finite-rank proofs, and protocol/checkpoint hashes.
    checker = module('transfer_upstream_checker',CAS/'check_visibility_cascade_v3.sage')
    # Loading the legacy checker reloads its V1 module. Rebind AFTER that load.
    engine, folder, p = bind_job(case)
    checker.v2 = engine
    result = checker.replay(start=17)
    require(len(result['stages']) == len(terminal['stages']), 'Partial replay')
    require(result['rank_lower_bound'] == terminal['final_rank_lower_bound'], 'Terminal rank mismatch')
    geo = engine.load('prospective_half_lattice_v3.sage')
    from sage.all import ZZ, matrix, pari
    model, _ = engine.v1.seed(None)
    metrics = []
    for s in terminal['stages']:
        path = folder/f"replay-M17/epoch-{s['epoch']:02d}/selection.json"
        selection = read_json(path)
        wd = path.parent
        prefix = []
        for j, cp in enumerate(sorted(wd.glob('chart-*.json'))):
            chart = read_json(cp); prefix.append(chart)
            snapshot = wd/f'cloud-{j:03d}.json'
            audit = read_json(wd/f'mod2-{j:03d}.json')
            saved = read_json(snapshot)
            require(audit['input_sha256'] == sha(snapshot) and within(ROOT,audit['input_path']) == snapshot, 'Audit/transcript binding failed')
            require(saved['charts'] == prefix and saved['curve'] == list(map(str,model)), 'Snapshot includes unexecuted or wrong-curve charts')
            seed_points = saved['final_state']['state']['reductions']['points']
            require(seed_points == selection['basis'], 'Cloud seed differs from epoch basis')
            expected = [list(map(str,q)) for q in seed_points]
            seen = {(F(x),abs(F(y))) for x,y in expected}
            for old in prefix:
                for q in old['search']['finite_curve_points']:
                    point = [str(F(q['x'])),str(F(q['y']))]
                    key = F(point[0]),abs(F(point[1]))
                    if key not in seen: expected.append(point);seen.add(key)
            require(expected == audit['points'], 'Certified cloud is not exactly input plus returned witnesses')
        basis = tuple(tuple(map(F,q)) for q in selection['basis'])
        hg, asym = geo.canonical_height_gram(model,basis)
        G = matrix(ZZ,geo.rounded_gram(hg,1000000))
        U = matrix(ZZ,pari(G).qflllgram()).transpose()
        require([list(map(int,r)) for r in G.rows()] == selection['rounded_gram'], 'Metric replay failed')
        require([list(map(int,r)) for r in U.rows()] == selection['LLL'] and str(asym) == selection['height_asymmetry'], 'LLL/asymmetry changed')
        metrics.append({'epoch':s['epoch'],'selection_sha256':sha(path)})
    write_new(folder/'full-replay.json',result)
    write_new(folder/'metric-replay.json',metrics)
    paths = [folder/'protocol.json',folder/'full-replay.json',folder/'metric-replay.json',folder/'replay-M17/terminal.json']
    write_new(folder/'verified.json', {'status':'PASS_INDEPENDENT_TRANSFER_REPLAY', 'case':case,
        'initial_rank':p['initial_rank'],'rank_lower_bound':result['rank_lower_bound'],
        'gain':result['rank_lower_bound']-p['initial_rank'],'stop_reason':terminal['stop_reason'],
        'charts':terminal['charts'],'bindings':{str(x.relative_to(ROOT)):sha(x) for x in paths},
        'scope':'Finite curve-specific recovery/extension. Check current catalogue after any gain; do not claim a new curve or record automatically.'})


def supervise(action, case=None, v3_replay=None):
    from research_runtime.supervisor import Limits, run
    if action == 'prepare-worker':
        logdir = D.parent/'v3-transfer-11952-v1-preparation'
        seconds = LIMITS['prepare_wall_seconds']
    else:
        logdir = D/case
        seconds = LIMITS['case_wall_seconds'] if action == 'run-worker' else LIMITS['replay_wall_seconds']
    logdir.mkdir(parents=True,exist_ok=True)
    stem = action.replace('-worker','')
    require(not (logdir/(stem+'.supervisor.json')).exists(), 'Preserve previous supervised attempt. Review before any resume.')
    command = [sys.executable,str(Path(__file__).resolve()),action]
    if case: command += ['--case',case]
    if v3_replay: command += ['--v3-replay',str(v3_replay.resolve())]
    env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',V3_TRANSFER_SUPERVISED='1')
    report = run(command, limits=Limits(seconds,LIMITS['rss_bytes']),cwd=ROOT,env=env,
        log_path=logdir/(stem+'.log'),checkpoint_path=logdir/(stem+'.supervisor.json'))
    require(report['outcome'] == 'completed', f"{action} stopped: {report['outcome']}; retain as a method/resource outcome, not a negative rank result")


def next_case():
    require((D/'roster.json').exists(), 'Run prepare first')
    with (D/'controller.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for spec in CASE_SPECS:
            folder = D/spec['id']
            if (folder/'verified.json').exists():
                check_bindings(ROOT,read_json(folder/'verified.json')['bindings']); continue
            bind_job(spec['id'])
            if not (folder/'replay-M17/terminal.json').exists(): supervise('run-worker',spec['id'])
            supervise('replay-worker',spec['id'])
            print(json.dumps(read_json(folder/'verified.json'),indent=2))
            return
        print('All four fixed cases finished. No automatic parameter/family expansion.')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=['prepare','next','status','prepare-worker','run-worker','replay-worker'])
    ap.add_argument('--v3-replay',type=Path);ap.add_argument('--case',choices=[s['id'] for s in CASE_SPECS])
    a=ap.parse_args()
    if a.action.endswith('-worker'):
        require(os.environ.get('V3_TRANSFER_SUPERVISED') == '1', 'Use prepare/next: worker commands must run under the shared resource supervisor')
    if a.action in ('prepare','prepare-worker'):
        require(a.v3_replay is not None,'Supply the FINAL V3 M17 full-replay JSON, not its progress file')
        if a.action=='prepare': supervise('prepare-worker',v3_replay=a.v3_replay)
        else: prepare(a.v3_replay)
    elif a.action=='next': next_case()
    elif a.action=='run-worker': run_case(a.case)
    elif a.action=='replay-worker': replay_case(a.case)
    else:
        for s in CASE_SPECS:
            p=D/s['id']/'verified.json'
            print(s['id'],json.dumps(read_json(p)) if p.exists() else 'NOT_YET_VERIFIED')

if __name__=='__main__': main()
