#!/usr/bin/env python3
"""Exact preparation and bounded phase adapters for the new 60-fibre panel."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import json

from v3_warm_support import atomic as _atomic, read, require, sha
from select_r17_60_panel import ROOT, ART

D = ROOT / 'artifacts/local/elliptic-curves/r17-60-panel-v1'


def atomic(path, value, **kwargs):
    return _atomic(path, json.loads(json.dumps(value)), **kwargs)


def native_check(packet, row):
    import compact_atlas_specialization as atlas
    from memory_rank_certificate import checked_rank
    family = next(f for f in read(atlas.ATLAS)['families'] if f['family'] == row['family'])
    model, generic = atlas.specialize(family, row['parameter'])
    points = tuple(tuple(map(F, p)) for p in packet['points'])
    require(model == tuple(map(F, packet['curve'])) and points[:17] == generic,
            'native generic prefix differs')
    old = json.loads(json.dumps(packet['proof']))
    fresh = checked_rank(model, points, [s['prime'] for s in old['signatures']], old['no_rational_2_torsion_prime'])
    require(json.loads(json.dumps(fresh)) == old and len(points) == packet['rank_lower_bound'], 'finite independence differs')


def generic():
    import compact_atlas_specialization as atlas
    import certify_compact_r17_candidates as cert
    from memory_rank_certificate import checked_rank
    from future_point_admission import FinitePointAdmission
    from v3_warm_engine import certified_state
    root = D / 'generic'
    rows = read(D / 'roster.json')['rows']
    families = {f['family']: f for f in read(atlas.ATLAS)['families']}
    protocol = {'point_searches': 0, 'inputs': {str(p.relative_to(ROOT)): sha(p) for p in
                (D / 'roster.json', D / 'protocol.json', atlas.ATLAS)}}
    atomic(root / 'protocol.json', protocol, immutable=True)
    records = []
    for row in rows:
        path = root / row['id'] / 'seed-M17.json'
        if not path.exists():
            model, points = atlas.specialize(families[row['family']], row['parameter'])
            require(model == tuple(map(F, row['model'])), 'selected equation differs')
            admission = FinitePointAdmission(model, points, prime_bound=1000)
            torsion = cert.find_two_torsion_certificate_prime(model, prime_bound=200)
            proof = checked_rank(model, points, admission.primes, torsion)
            packet = {k: row[k] for k in ('id', 'family', 'parameter')}
            packet.update(curve=list(map(str, model)), points=[list(map(str, p)) for p in points],
                          generic_rank=17, rank_lower_bound=17, proof=proof)
            certified_state(model, points, proof)
            atomic(path, packet, immutable=True)
        native_check(read(path), row)
        records.append({'id': row['id'], 'path': str(path.relative_to(ROOT)), 'sha256': sha(path)})
        print('R17_60_GENERIC_VERIFIED', row['id'], flush=True)
    atomic(root / 'prepared.json', {'status': 'PASS_SIXTY_GENERIC17_PACKETS', 'records': records,
                                   'protocol_sha256': sha(root / 'protocol.json')}, immutable=True)
    atomic(root / 'verified.json', {'status': 'PASS_INDEPENDENT_GENERIC17_REPLAY',
                                   'prepared_sha256': sha(root / 'prepared.json')}, immutable=True)


def reconcile_seed(case):
    from memory_rank_certificate import checked_rank
    from future_point_admission import FinitePointAdmission
    from v3_warm_engine import certified_state
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend
    folder = D / 'cases' / case
    seedrun = folder / 'seed-search'
    terminal, verified = read(seedrun / 'terminal.json'), read(seedrun / 'verified.json')
    require(verified['status'] == 'PASS_INDEPENDENT_FIRST_SEED_REPLAY' and
            verified['terminal_sha256'] == sha(seedrun / 'terminal.json'), 'seed replay required')
    initial = read(seedrun / 'seed.json')
    model = tuple(map(F, initial['curve']))
    state = certified_state(model, initial['points'], initial['proof'])
    admission = FinitePointAdmission(model, tuple(tuple(map(F, p)) for p in terminal['points']))
    bindings = {n: sha(seedrun / n) for n in ('terminal.json', 'verified.json', 'protocol.json')}
    gains = []
    for i in range(terminal['charts']):
        path = seedrun / f'chart-{i:04d}.json'; chart = read(path); mapping = chart['mapping']
        bindings[path.name] = sha(path)
        search = PointedQuarticSearch(state=state, centre={'coefficients': mapping['centre']['representative']},
                                      coordinate_policy=mapping['coordinate_policy'])
        for point in backend.replay(search, mapping, chart['search']):
            before = len(admission.points)
            admission.consider(point)
            if len(admission.points) > before:
                gains.append({'seed_call': i+1, 'after': len(admission.points), 'point': list(map(str, point))})
    proof = checked_rank(model, admission.points, admission.primes, initial['proof']['no_rational_2_torsion_prime'])
    packet = {'status': 'PASS_REPLAYED_SEED_CLOUD', 'curve': initial['curve'],
              'points': [list(map(str, p)) for p in admission.points], 'proof': proof,
              'rank_lower_bound': len(admission.points), 'generic_rank': 17,
              'reconciliation_gains': gains, 'bindings': bindings, 'point_searches': 0}
    row = next(r for r in read(D / 'roster.json')['rows'] if r['id'] == case)
    native_check(packet, row)
    atomic(folder / 'seed-reconciled.json', packet, immutable=True)
    print('R17_60_SEED_CLOUD', case, terminal['rank_lower_bound'], len(admission.points), flush=True)


def complement(case):
    """Sixteen exact maximum parents outside the single seed-winning mask span."""
    from sage.all import ZZ, matrix, pari
    import numpy as np
    from visibility_lattice_fast import IntegerExactParity
    from visibility_lattice_v2 import ExactParity
    folder = D / 'cases' / case; out = folder / 'complement-preparation'
    packet = read(folder / 'seed-reconciled.json')
    row = next(r for r in read(D / 'roster.json')['rows'] if r['id'] == case)
    native_check(packet, row)
    require(packet['rank_lower_bound'] >= 18, 'no first seed')
    terminal = read(folder / 'seed-search/terminal.json')
    require(terminal['rank_lower_bound'] == 18, 'standalone first M18 missing')
    winning = read(folder / 'seed-search' / f"chart-{terminal['charts']-1:04d}.json")
    word = winning['mapping']['centre']['representative']
    original = sum((int(v) % 2) << j for j, v in enumerate(word))
    catalogue = ART / 'r17_exact_maximum_parity_classes_v1.json'
    family = next(f for f in read(catalogue)['families'] if f['family'] == row['family'])
    g = 2 * matrix(ZZ, family['gram'])
    u = matrix(ZZ, pari(g).qflllgram()).transpose(); inv = u.inverse()
    require(abs(u.det()) == 1, 'generic LLL not unimodular')
    fast = IntegerExactParity((u*g*u.transpose()).rows())
    ref = ExactParity((u*g*u.transpose()).rows())
    candidates, checks = [], []
    for mask in sorted(r['mask'] for r in family['classes'] if r['mask'] not in (0, original)):
        w = matrix(ZZ, 1, 17, [(mask >> j) & 1 for j in range(17)])
        residue = tuple(int(x) % 2 for x in (w*inv).row(0))
        starts, _ = fast.babai(np.asarray([residue], dtype=np.int64)); start = tuple(map(int, starts[0]))
        proof = fast.solve(residue, start, 2000000)
        require(proof == ref.solve(residue, start, 2000000), 'generic reference CVP differs')
        require(proof['norm'] == 2*family['exact_maximum_parity_minimum'], 'not a maximum class')
        words = []
        for v in proof['minima']:
            z = tuple(map(int, (matrix(ZZ, 1, 17, v)*u).row(0)))
            if next(x for x in z if x) < 0:
                z = tuple(-x for x in z)
            require(sum((x % 2) << j for j,x in enumerate(z)) == mask, 'parity transport differs')
            words.append(z)
        candidates.append({'mask': mask, 'word': list(min(words)), 'norm': proof['norm']})
        checks.append({'mask': mask, 'reduced_seed': start, 'proof': proof})
    pivots = {}
    def insert(mask):
        for k in sorted(pivots, reverse=True):
            if (mask >> k) & 1:
                mask ^= pivots[k]
        if mask:
            pivots[mask.bit_length()-1] = mask
        return bool(mask)
    insert(original); chosen = []
    for c in candidates:
        if insert(c['mask']):
            chosen.append(c)
        if len(chosen) == 16:
            break
    for c in candidates:
        if len(chosen) == 16:
            break
        if c not in chosen:
            chosen.append(c)
    require(len(chosen) == 16, 'sixteen complement parents required')
    protocol = {'family': row['family'], 'parameter': row['parameter'], 'initial_rank': packet['rank_lower_bound'],
                'generic_rank': 17, 'point_searches': 0,
                'inputs': {str(p.relative_to(ROOT)): sha(p) for p in
                           (D/'protocol.json', folder/'seed-reconciled.json', catalogue, Path(__file__))}}
    atomic(out / 'protocol.json', protocol, immutable=True)
    atomic(out / f"seed-M{packet['rank_lower_bound']}.json", packet, immutable=True)
    atomic(out / 'generic-cvp-proofs.json', {'gram': [list(map(int,r)) for r in g.rows()],
          'LLL': [list(map(int,r)) for r in u.rows()], 'checks': checks}, immutable=True)
    atomic(out / 'derivation.json', {'original_masks': [original], 'new_masks': [c['mask'] for c in chosen],
          'extended_span_rank': len(pivots), 'all_candidates': candidates,
          'rule': 'Ascending exact maximum masks outside the winning-mask span; extend that span then fill to16. Seed cloud is reconciled before exposure.'}, immutable=True)
    bank = {'status': 'COMPLETE_FROZEN_COMPLEMENT_PARENT_SUBSET', 'dimension': 17, 'generic_gram_scale': 2,
            'rows': chosen, 'shells': sorted({c['norm'] for c in chosen}),
            'protocol_sha256': sha(out/'protocol.json'), 'derivation_sha256': sha(out/'derivation.json'),
            'generic_cvp_sha256': sha(out/'generic-cvp-proofs.json')}
    atomic(out / 'anchor-bank.json', bank, immutable=True)
    atomic(out / 'prepared.json', {'status': 'PASS_COMPLEMENT_PARENT_PREPARATION', 'point_searches': 0,
           'files': {p.name: sha(p) for p in sorted(out.glob('*.json')) if p.name != 'prepared.json'}}, immutable=True)


def segment(case, index, remaining, phase):
    """Reconciliation may restart a landscape, but never reset the 100-call cap."""
    import run_complement_seed_v3 as runner
    folder = D / 'cases' / case
    prep = folder / ('complement-preparation' if index == 0 else f'preparation-{index:02d}')
    out = folder / f'complement-{index:02d}'
    runner.PREP = runner.BANK = prep
    if phase == 'search':
        require(0 < remaining <= 100, 'invalid remaining point allowance')
        if not out.exists():
            policy = runner.freeze(out)
            policy['max_charts'] = remaining
            policy['panel'] = {'total_complementary_allowance': 100, 'segment_index': index,
                               'used_before_segment': 100-remaining}
            policy['inputs'][str((D/'protocol.json').relative_to(ROOT))] = sha(D/'protocol.json')
            policy['sources'][str(Path(__file__).relative_to(ROOT))] = sha(Path(__file__))
            atomic(out/'protocol.json', policy)
        require(read(out/'protocol.json')['max_charts'] == remaining, 'segment budget changed')
        runner.search(out)
    else:
        runner.replay(out)


def main():
    import run_r17_60_panel as control
    control.check_protocol()
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=('generic', 'seed-prepare', 'seed-search', 'seed-replay',
                                   'seed-cloud', 'complement-prepare', 'search', 'replay', 'reconcile', 'rebuild'))
    p.add_argument('--case'); p.add_argument('--index', type=int, default=0)
    p.add_argument('--remaining', type=int, default=100)
    a = p.parse_args(); case = a.case
    if case is not None:
        require(case in {r['id'] for r in read(D/'roster.json')['rows']}, 'unknown case')
    folder = D / 'cases' / str(case)
    if a.phase == 'generic':
        generic()
    elif a.phase.startswith('seed-') and a.phase != 'seed-cloud':
        import run_fresh6_seed_confirmation_v2 as seed
        if a.phase == 'seed-prepare':
            seed.prepare(D/'generic'/case/'seed-M17.json', folder/'seed-search')
        else:
            seed.execute(folder/'seed-search', a.phase == 'seed-replay')
    elif a.phase == 'seed-cloud':
        reconcile_seed(case)
    elif a.phase == 'complement-prepare':
        complement(case)
    elif a.phase in ('search', 'replay'):
        segment(case, a.index, a.remaining, a.phase)
    elif a.phase == 'reconcile':
        import reconcile_verified_v3_cloud as rc
        rc.run(folder/f'complement-{a.index:02d}', folder/f'reconciled-{a.index:02d}.json')
    else:
        import prepare_reconciled_v3_continuation as rebuild
        prep = folder/('complement-preparation' if a.index == 0 else f'preparation-{a.index:02d}')
        rebuild.run(prep, folder/f'reconciled-{a.index:02d}.json', folder/f'preparation-{a.index+1:02d}')
    control.check_protocol()


if __name__ == '__main__':
    main()
