#!/usr/bin/env python3
"""Bounded retained-data cost audit; constructs no new point-search boxes.

Run with Sage Python. The withheld curve-90 point is used only for its existing
27-point proof and timing, never to select a prospective chart.
"""
import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import sys
import time

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
LOCAL = ROOT / 'artifacts/local/elliptic-curves'
ART = ROOT / 'artifacts/generated-results/elliptic-curves'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(output):
    if output.exists():
        raise FileExistsError('preserve the previous cost audit')
    from research_runtime.store import checkpoint
    from v3_warm_engine import load, certified_state
    from memory_rank_certificate import checked_rank
    inputs = {}

    def read(path):
        raw = path.read_bytes()
        inputs[str(path.relative_to(ROOT))] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    report = {'status': 'RUNNING', 'point_searches': 0, 'timings': {},
              'inputs': inputs, 'epochs': [],
              'claim_boundary': 'Retained timings and bounded preparation benchmarks only; '
              'no new points, speedup claim, or exact-rank assertion.'}

    def timed(name, fn):
        start = time.monotonic()
        result = fn()
        report['timings'][name] = time.monotonic() - start
        checkpoint(output, report)
        print(name, report['timings'][name], flush=True)
        return result

    folder = LOCAL / 'curve302-seed-universality-panel-v1/residual-strict-03'
    verified = read(folder / 'seeded-verified.json')
    if verified['status'] != 'PASS_INDEPENDENT_SEEDED_V3_REPLAY':
        raise ArithmeticError('completed replay required')
    terminal = folder / 'replay-M17/terminal.json'
    if digest(terminal) != verified['terminal_sha256']:
        raise ArithmeticError('sealed terminal changed')
    read(terminal)
    for old in verified['stages']:
        wd = folder / 'replay-M17' / f"epoch-{old['epoch']:02d}"
        stage = read(wd / 'stage.json')
        count = stage['charts']
        if count != old['charts_replayed']:
            raise ArithmeticError('replayed chart count differs')
        paths = [wd / f'chart-{j:03d}.json' for j in range(count)]
        charts = [read(p) for p in paths]
        if [c['index'] for c in charts] != list(range(count)):
            raise ArithmeticError('chart ordering differs')
        backend = sum(c['search']['wall_seconds'] for c in charts)
        report['epochs'].append({
            'rank_before': stage['before'], 'charts': count,
            'epoch_wall_seconds': stage['wall_seconds'],
            'backend_wall_seconds': backend,
            'non_backend_wall_seconds': stage['wall_seconds'] - backend,
            'chart_file_bytes': sum(p.stat().st_size for p in paths),
            'cumulative_cloud_file_bytes': sum((wd / f'cloud-{j:03d}.json').stat().st_size
                                              for j in range(count)),
            'chart_occurrences_in_cumulative_snapshots': count * (count + 1) // 2})
    checkpoint(output, report)

    data = read(ART / 'mw16_rank27_visibility_input_v1.json')
    proof = data['point_proof']
    basis = tuple(tuple(map(Fraction, p)) for p in proof['discovery_points'])
    model = tuple(map(Fraction, proof['discovery_curve']))
    if (len(basis) != 27 or proof['parameter'] != '-1867/270'
            or data['generic_points'] != proof['discovery_points'][:16]
            or data['initial_points'] != proof['discovery_points'][:26]):
        raise ArithmeticError('curve-90 basis identity failed')
    old = proof['rank_certificate']
    fresh = timed('curve90_complete_rank27_certificate', lambda: checked_rank(
        model, basis, [s['prime'] for s in old['signatures']],
        old['no_rational_2_torsion_prime']))
    if json.loads(json.dumps(fresh)) != old:
        raise ArithmeticError('independent rank certificate differs')
    state = timed('curve90_certified_state', lambda: certified_state(model, basis, fresh))
    report['curve90_certified_rank'] = state.rank
    engine = load('cost_audit_v3', CAS / 'adaptive_visibility_cascade_v3.sage')
    geo = load('cost_audit_geometry', CAS / 'prospective_half_lattice_v3.sage')
    gram, asym = timed('curve90_height_gram27', lambda: geo.canonical_height_gram(model, basis))
    report['curve90_rounded_gram'] = geo.rounded_gram(gram, 1000000)
    report['height_asymmetry'] = str(asym)
    columns = timed('curve90_fingerprints27', lambda: engine.fingerprints(model, basis))
    report['curve90_fingerprint_columns'] = columns
    report['extension_scaling'] = [
        {'rank': r, 'generic_rank': g, 'cosets_for_32_anchors': 32 * 2**(r-g)}
        for g in (16, 17) for r in (26, 27, 30, 31)]
    report['sources'] = {str(p.relative_to(ROOT)): digest(p)
        for p in sorted({Path(m.__file__).resolve() for m in list(sys.modules.values())
                         if getattr(m, '__file__', None)
                         and Path(m.__file__).resolve().is_relative_to(CAS)})
        if p.is_file()}
    if any(digest(ROOT / p) != h for p, h in inputs.items()):
        raise ArithmeticError('input changed during audit')
    report['status'] = 'PASS_RETAINED_COST_AND_CURVE90_PREPARATION_AUDIT'
    checkpoint(output, report)
    print(report['status'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    audit(parser.parse_args().output)
