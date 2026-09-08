"""Explicit Sage services for the same-parent pilot; no wrapper monkeypatching.

The numerical V3 landscape/selector and backend are imported unchanged. The
checkpoint loop follows the verified warm engine, with determinant-1092 labels
and generic-17 intake. Every case is a separate supervised process.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction as F
import hashlib
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time
import uuid

import det1092_v3_contract as c
from v3_warm_support import (atomic, assert_basis, bindings, check_chart, curve_tuple,
    indexed_paths, point_tuple, read, require, sha, terminal_structure)


def copy_immutable(source, target):
    target = Path(target); target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        require(sha(source) == sha(target), 'immutable source copy differs')
        return
    fd, name = tempfile.mkstemp(prefix='.copy-', dir=target.parent)
    try:
        with os.fdopen(fd, 'wb') as out, Path(source).open('rb') as inp:
            shutil.copyfileobj(inp, out); out.flush(); os.fsync(out.fileno())
        try:
            os.link(name, target)
        except FileExistsError:
            require(sha(source) == sha(target), 'concurrent immutable source copy differs')
    finally:
        os.unlink(name)


def prepare():
    from sage.all import ZZ, matrix, pari, QQ, PolynomialRing
    import sage.version
    from v3_transfer_contract import validate_gate
    from v3_warm_engine import certified_state, load
    import certify_compact_r17_candidates as cert

    selection, original_path = c.selection_input()
    pilot, reserve = c.choose(selection)
    for path, expected in ((c.REDUCED_PARENT,c.PARENT_BLOB), (c.ORBITS,c.ORBITS_BLOB),
                           (c.LATTICE,c.LATTICE_BLOB)):
        require(c.git_blob(path) == expected, 'retained generic input changed: '+str(path))
    original = read(c.CALIBRATION/'protocol.json')
    gate = validate_gate(c.ROOT, c.CALIBRATION, c.CALIBRATION/'replay-M17.json',
                         c.CALIBRATION/'metric-replay-M17.json')
    parent = read(c.REDUCED_PARENT)
    generic = read(c.CALIBRATION/'generic.json')
    gram = generic['generic_height_gram']
    G = matrix(ZZ, gram)
    require(G.nrows() == 17 and G.ncols() == 17 and G.det() == 1092 and
            G.is_positive_definite(), 'wrong ordered generic lattice')
    lattice = read(c.LATTICE)
    require(lattice['status'] == 'PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT' and
            lattice['orbits_tsv_sha256'] == sha(c.ORBITS), 'orbit/lattice binding differs')
    # Exact polynomial identities, independent of a chosen specialization.
    R = PolynomialRing(QQ, 't'); K = R.fraction_field()
    def rational(row):
        return K(R(row['numerator']))/K(R(row['denominator']))
    aa = [rational(row) for row in parent['a_invariants']]
    require(aa[:3] == [0,0,0], 'nonshort generic model')
    require(len(parent['basis_weierstrass_coordinates']) == 17, 'generic marking length differs')
    for x, y in parent['basis_weierstrass_coordinates']:
        x, y = rational(x), rational(y)
        require(y*y == x**3+aa[3]*x+aa[4], 'generic section identity failed')
    numerical = load('det1092_prepare_numeric', c.CAS/'adaptive_visibility_cascade_v3.sage')
    require(numerical.sources() == original['sources'], 'V3 numerical sources changed since calibration')
    impl = c.own_sources()
    source_paths = [original_path, c.REDUCED_PARENT, c.ORBITS, c.LATTICE,
                    c.CALIBRATION/'protocol.json', c.CALIBRATION/'generic.json']
    inputs = {str(p.relative_to(c.ROOT)):sha(p) for p in source_paths}
    inputs.update(gate['bindings'])
    c.D.mkdir(parents=True, exist_ok=True)
    # Plan is committed before ANY seed-rank or point-search outcome. No refill.
    atomic(c.D/'plan.json', {'pilot':pilot,'reserve':reserve, 'inputs':inputs,
           'implementation_sources':impl, 'scope':c.CLAIM}, immutable=True)
    projected = {'status':'PASS', 'selected':[{k:r[k] for k in c.FIELDS} for r in selection['selected']]}
    atomic(c.D/'frozen-selection.json', projected, immutable=True)
    atomic(c.D/'gate.json', gate, immutable=True)
    atomic(c.D/'parent-sections.json', {k:parent[k] for k in
           ('a_invariants','basis_weierstrass_coordinates')}, immutable=True)
    jobs = {}
    for row in pilot:
        model, points = c.specialize(parent, row)
        proof = cert.checked_rank(model, points)
        state = certified_state(model, points, proof)
        assert_basis(state, model, points, 17)
        folder = c.D/row['id']; folder.mkdir(exist_ok=True)
        copy_immutable(c.ORBITS, folder/'orbits.tsv')
        seed = {'family':'det1092', 'parameter':row['parameter'], 'curve':list(map(str,model)),
                'points':[list(map(str,p)) for p in points], 'initial_rank':17, 'generic_rank':17,
                'generic_height_gram':gram, 'parent_blob':c.PARENT_BLOB}
        atomic(folder/'seed-input.json', seed, immutable=True)
        atomic(folder/'seed-proof.json', proof, immutable=True)
        policy = c.trial_policy(original)
        policy.update(sources=numerical.sources(), implementation_sources=impl,
            software={'sage':sage.version.version, 'pari':str(pari.version()), 'python':sys.version},
            inputs={str(p.relative_to(c.ROOT)):sha(p) for p in
                    (folder/'seed-input.json',folder/'seed-proof.json',folder/'orbits.tsv')},
            frozen_302_protocol_sha256=sha(c.CALIBRATION/'protocol.json'),
            selection_role=row['pilot_role'], resource_limits=c.RESOURCE)
        atomic(folder/'protocol.json', policy, immutable=True)
        jobs[row['id']] = sha(folder/'protocol.json')
        print('PREPARED_GENERIC17', row['id'], row['pilot_role'], flush=True)
    for name in ('plan.json','frozen-selection.json','gate.json','parent-sections.json'):
        inputs[str((c.D/name).relative_to(c.ROOT))] = sha(c.D/name)
    atomic(c.D/'roster.json', {'status':'READY_EIGHT_FIBRES', 'pilot':pilot,'reserve':reserve,
           'inputs':inputs,'implementation_sources':impl,'jobs':jobs,'claim_boundary':c.CLAIM}, immutable=True)
    c.validate_roster()
    print('PREPARED_EIGHT_NO_SEARCH', ','.join(r['id'] for r in pilot), flush=True)


def setup(case):
    from sage.all import pari
    import sage.version
    from v3_warm_engine import Context, certified_state, load
    roster = c.validate_roster()
    row = next((r for r in roster['pilot'] if r['id'] == case), None)
    require(row is not None, 'case is not in the frozen eight-fibre pilot')
    folder = c.D/case; p = read(folder/'protocol.json')
    require(p['implementation_sources'] == c.own_sources(), 'implementation differs from frozen trial')
    for key in ('inputs','sources','implementation_sources'):
        bindings(c.ROOT, p[key])
    require(p['software']['sage'] == sage.version.version and
            p['software']['pari'] == str(pari.version()), 'Sage/PARI changed since preparation')
    require(sha(Path('/usr/bin/gp')) == p['gp_sha256'], 'GP binary changed')
    seed, proof = read(folder/'seed-input.json'), read(folder/'seed-proof.json')
    model, points = c.specialize(read(c.D/'parent-sections.json'), row)
    require(curve_tuple(seed['curve']) == model and point_tuple(seed['points']) == points,
            'frozen seed no longer equals ordered generic specialization')
    state = certified_state(model, points, proof); assert_basis(state, model, points, 17)
    engine = load('det1092_frozen_v3_numeric', c.CAS/'adaptive_visibility_cascade_v3.sage')
    engine.D, engine.ORBITS = folder, folder/'orbits.tsv'
    engine.v1.D, engine.v1.ORBITS = folder, engine.ORBITS
    require(engine.sources() == p['sources'], 'numerical engine source-set differs')
    sources = {**p['sources'], **p['implementation_sources']}
    engine.v1.guard()
    # This generic-only worker may not read other jobs, original score inputs,
    # or the calibration points. Binding happens BEFORE the guarded search.
    for forbidden in (c.D/'frozen-selection.json', c.CALIBRATION/'fixed-M30.json'):
        try:
            forbidden.open()
        except PermissionError:
            pass
        else:
            raise ValueError('artifact boundary failed: '+str(forbidden))
    return Context(case, folder, c.ROOT, p, engine, model, state, sources)



def preflight(ctx):
    """Generic seed -> one nonzero shell chart -> exact map check; no search."""
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    word = c.preflight_word(ctx.engine.ORBITS, ctx.state.rank)
    mapper = ctx.engine.load('factor_free_pari_mapping.sage')
    mapper.pari.allocatemem(256000000, silent=True)
    mapping = mapper.mapping(ctx.model, point_tuple(ctx.state.basis), {'representative':word})
    search = PointedQuarticSearch(state=ctx.state, centre={'coefficients':word},
                                 coordinate_policy=mapping['coordinate_policy'])
    backend.validate_map(search, mapping)
    print('PREFLIGHT_PASS', ctx.case, 'seed_rank', ctx.state.rank, 'exact_nonzero_chart', flush=True)



def run_search(ctx):
    """Checkpointed V3 loop; numerical choices are those of the frozen engine."""
    from v3_warm_engine import certified_state, restore_state, _audit, _odd
    from v3_warm_replay import verify_landscape
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    engine, model, state = ctx.engine, ctx.model, ctx.state
    folder, policy = ctx.folder, ctx.policy
    out = folder/'replay-M17'; out.mkdir(exist_ok=True)
    if (out/'terminal.json').exists():
        terminal_structure(folder)
        print('SEALED_SEARCH_REUSED', ctx.case, flush=True); return
    mapper = engine.load('factor_free_pari_mapping.sage'); mapper.pari.allocatemem(256000000, silent=True)
    total, stages, tested = 0, [], set()
    reason = 'EPOCH_BUDGET_EXHAUSTED'
    for epoch in range(policy['max_epochs']):
        if total >= policy['max_charts']:
            reason = 'CHART_BUDGET_EXHAUSTED'; break
        wd = out/f'epoch-{epoch:02d}'; wd.mkdir(exist_ok=True)
        basis, before = point_tuple(state.basis), state.rank
        started = time.monotonic(); existing = indexed_paths(wd)
        if (wd/'selection.json').exists():
            selection = read(wd/'selection.json')
            require(point_tuple(selection['basis']) == basis, 'resume seed differs from schedule')
            verify_landscape(engine, model, basis, tested, wd, policy)
            centres = selection['centres']
        else:
            require(not existing and not (wd/'stage.json').exists(), 'charts/stage without selection')
            unfinished = list(wd.glob('anchor-*'))
            if unfinished:
                archive = wd/'incomplete-landscapes'/str(uuid.uuid4()); archive.mkdir(parents=True)
                for path in unfinished: path.replace(archive/path.name)
            centres = engine.landscape(model, basis, tested, wd, policy)
            selection = read(wd/'selection.json')
        require(bool(centres), 'no finalists: engineering stop, not an arithmetic negative')
        require(len(existing) <= len(centres), 'extra retained charts outside schedule')
        if (wd/'epoch-seed-state.json').exists():
            state = restore_state(read(wd/'epoch-seed-state.json'), model, basis)
        elif existing and (wd/'cloud-000.json').exists():
            state = restore_state(read(wd/'cloud-000.json')['final_state'], model, basis)
        else:
            atomic(wd/'epoch-seed-state.json', state.record(), immutable=True)
        stage_seed = state; charts, censored, audit_path = [], 0, None
        for j, centre in enumerate(centres):
            if total >= policy['max_charts']: break
            path = existing[j] if j < len(existing) else wd/f'chart-{j:03d}.json'
            if path.exists():
                chart = read(path); check_chart(chart, selection, j)
                search = PointedQuarticSearch(state=stage_seed, centre={'coefficients':centre['representative']},
                                              coordinate_policy=chart['mapping']['coordinate_policy'])
                backend.replay(search, chart['mapping'], chart['search'])
            else:
                mapping = mapper.mapping(model, basis, centre)
                search = PointedQuarticSearch(state=stage_seed, centre={'coefficients':centre['representative']},
                                              coordinate_policy=mapping['coordinate_policy'])
                transcript, points = backend.execute(search, mapping, policy['height'],
                                                      policy['seconds_per_chart'], policy['gp_sha256'])
                require(backend.replay(search, mapping, transcript) == points, 'exact witness replay changed')
                chart = {'index':j,'centre':centre,'mapping':mapping,'search':transcript}
                atomic(path, chart, immutable=True)
            charts.append(chart); tested.add(tuple(centre['point'])); total += 1
            censored += chart['search']['status'] != 'bounded_search_complete'
            snapshot, audit_path = wd/f'cloud-{j:03d}.json', wd/f'mod2-{j:03d}.json'
            if snapshot.exists():
                saved = read(snapshot)
                require(saved['charts'] == charts and curve_tuple(saved['curve']) == model and
                        point_tuple(saved['final_state']['state']['reductions']['points']) == basis,
                        'retained cumulative snapshot differs')
            else:
                atomic(snapshot, {'status':'INCREMENTAL_RETAINED_CLOUD','family':'det1092-v3-pilot',
                    'parameter':ctx.case,'curve':list(map(str,model)),'charts':charts,
                    'final_state':stage_seed.record(),'rank_lower_bound':before}, immutable=True)
            cloud = _audit(snapshot, audit_path)
            print(ctx.case, 'epoch', epoch, 'chart', j+1, '/', len(centres), 'certified', cloud['rank_lower_bound'], flush=True)
            if cloud['rank_lower_bound'] > before:
                enlarged = point_tuple(cloud['independent_points'])
                require(enlarged[:before] == basis, 'generic-first prefix lost')
                state = certified_state(model, enlarged, cloud['rank_certificate'])
                require(len(existing) <= j+1, 'stale charts retained after gain'); break
        require(audit_path is not None, 'epoch has no executed chart')
        odd = _odd(audit_path, wd/'modl.json')
        require(odd == {'3':state.rank,'5':state.rank}, 'independence audits disagree')
        if state.rank >= policy['target_rank']: reason = 'TARGET_LOWER_BOUND_REACHED'
        elif total >= policy['max_charts']: reason = 'CHART_BUDGET_EXHAUSTED'
        elif censored: reason = 'CENSORED_SEARCH'
        elif state.rank > before: reason = 'REBUILD_AFTER_CERTIFIED_GAIN'
        else: reason = 'COMPLETE_FINITE_NO_GAIN' if len(charts) == len(centres) else 'INCOMPLETE_EPOCH'
        stage = {'epoch':epoch,'before':before,'after':state.rank,'charts':len(charts),
            'stale_charts_cancelled':len(centres)-len(charts) if state.rank > before else 0,
            'audit':audit_path.name,'audit_sha256':sha(audit_path),'full_cosets_scored':selection['full_cosets_scored'],
            'stop_reason':reason,'censored_charts':censored,'wall_seconds':time.monotonic()-started}
        if (wd/'stage.json').exists():
            old = read(wd/'stage.json')
            require({k:v for k,v in old.items() if k != 'wall_seconds'} ==
                    {k:v for k,v in stage.items() if k != 'wall_seconds'}, 'sealed stage differs on resume')
            stage = old
        else: atomic(wd/'stage.json', stage, immutable=True)
        stages.append(stage); atomic(out/'stages.json', stages)
        if reason != 'REBUILD_AFTER_CERTIFIED_GAIN': break
    else: reason = 'EPOCH_BUDGET_EXHAUSTED'
    atomic(out/'terminal.json', {'status':'TERMINAL_BOUNDED_TRANSFER','stop_reason':reason,
        'initial_rank':policy['initial_rank'],'final_rank_lower_bound':state.rank,'stages':stages,'charts':total,
        'protocol_sha256':sha(folder/'protocol.json'),'read_paths':sorted(engine.v1.READS),
        'execution_sources':ctx.sources,'scope':'Search terminal pending independent replay. '+c.CLAIM}, immutable=True)
    terminal_structure(folder)
    print('DET1092_SEARCH_TERMINAL', ctx.case, state.rank, reason, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare','preflight','search','replay'))
    parser.add_argument('--case'); parser.add_argument('--session', type=Path, required=True)
    a = parser.parse_args()
    require(os.environ.get('DET1092_V3_SUPERVISED') == '1', 'use the supervised controller')
    session = read(a.session)
    require(session['sources'] == c.own_sources(), 'worker/controller source snapshot differs')
    bindings(c.ROOT, session['sources'])
    if a.action == 'prepare': prepare(); return
    ctx = setup(a.case)
    if a.action == 'preflight':
        preflight(ctx)
    elif a.action == 'search': run_search(ctx)
    else:
        from det1092_v3_replay import replay_case
        replay_case(ctx)


if __name__ == '__main__': main()
