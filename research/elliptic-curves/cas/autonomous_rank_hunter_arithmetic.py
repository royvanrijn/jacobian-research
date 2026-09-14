#!/usr/bin/env python3
"""Exact arithmetic adapters used by the autonomous R17 rank hunter.

All expensive point construction remains in the repository's existing frozen
seed/V3 engines.  This file only creates fresh native packets, reconciles full
returned clouds, builds new exact parent banks, and adapts bounded V3 segments.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from v3_warm_support import atomic as _atomic, read, require, sha

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
CATALOGUE = ART / 'r17_exact_maximum_parity_classes_v1.json'


def atomic(path, value, **kwargs):
    return _atomic(path, json.loads(json.dumps(value)), **kwargs)


def candidate(case):
    return read(case / 'candidate.json')


def _native_packet(row):
    import compact_atlas_specialization as atlas
    import certify_compact_r17_candidates as cert
    from memory_rank_certificate import checked_rank
    from future_point_admission import FinitePointAdmission
    from v3_warm_engine import certified_state

    family = next(f for f in read(atlas.ATLAS)['families'] if f['family'] == row['family'])
    model, points = atlas.specialize(family, row['parameter'])
    require(model == tuple(map(F, row['model'])), 'selected equation differs from native specialization')
    admission = FinitePointAdmission(model, points, prime_bound=1000)
    torsion = cert.find_two_torsion_certificate_prime(model, prime_bound=200)
    proof = checked_rank(model, points, admission.primes, torsion)
    certified_state(model, points, proof)
    return {
        'id': row['id'], 'family': row['family'], 'parameter': row['parameter'],
        'curve': list(map(str, model)), 'points': [list(map(str, p)) for p in points],
        'generic_rank': 17, 'rank_lower_bound': 17, 'proof': proof,
    }


def generic(case, replay=False):
    row = candidate(case)
    dataset = case / 'generic'
    seed_dir = dataset / 'candidate'
    seed = seed_dir / 'seed-M17.json'
    if not replay:
        require(not dataset.exists(), 'preserve existing generic dataset')
        seed_dir.mkdir(parents=True)
        packet = _native_packet(row)
        protocol = {
            'schema': 'autonomous-r17-generic17.v1', 'point_searches': 0,
            'candidate_sha256': sha(case / 'candidate.json'),
            'sources': {str(Path(__file__).relative_to(ROOT)): sha(Path(__file__))},
        }
        atomic(dataset / 'protocol.json', protocol, immutable=True)
        atomic(seed, packet, immutable=True)
        prepared = {'status': 'PASS_SINGLE_GENERIC17_PACKET', 'records': [{
            'id': row['id'], 'path': str(seed.relative_to(ROOT)), 'sha256': sha(seed)}],
            'protocol_sha256': sha(dataset / 'protocol.json')}
        atomic(dataset / 'prepared.json', prepared, immutable=True)
        return

    packet = _native_packet(row)
    require(read(seed) == json.loads(json.dumps(packet)), 'independent native17 replay differs')
    prepared = read(dataset / 'prepared.json')
    require(prepared['records'][0]['sha256'] == sha(seed), 'generic packet seal differs')
    atomic(dataset / 'verified.json', {
        'status': 'PASS_INDEPENDENT_GENERIC17_REPLAY',
        'prepared_sha256': sha(dataset / 'prepared.json'),
        'seed_sha256': sha(seed), 'point_searches': 0,
    }, immutable=True)


def reconcile_seed(case):
    from memory_rank_certificate import checked_rank
    from future_point_admission import FinitePointAdmission
    from v3_warm_engine import certified_state
    from pointed_quartic_search import PointedQuarticSearch
    import pari_pointed_backend as backend

    run = case / 'seed-search'
    terminal, verified = read(run / 'terminal.json'), read(run / 'verified.json')
    require(verified['status'] == 'PASS_INDEPENDENT_FIRST_SEED_REPLAY' and
            verified['terminal_sha256'] == sha(run / 'terminal.json'), 'seed replay required')
    initial = read(run / 'seed.json')
    model = tuple(map(F, initial['curve']))
    state = certified_state(model, initial['points'], initial['proof'])
    admission = FinitePointAdmission(model, tuple(tuple(map(F, p)) for p in terminal['points']))
    bindings, gains = {}, []
    for name in ('terminal.json', 'verified.json', 'protocol.json'):
        bindings[name] = sha(run / name)
    for i in range(terminal['charts']):
        path = run / f'chart-{i:04d}.json'
        chart = read(path); mapping = chart['mapping']; bindings[path.name] = sha(path)
        search = PointedQuarticSearch(state=state,
            centre={'coefficients': mapping['centre']['representative']},
            coordinate_policy=mapping['coordinate_policy'])
        for point in backend.replay(search, mapping, chart['search']):
            before = len(admission.points); admission.consider(point)
            if len(admission.points) > before:
                gains.append({'seed_call': i + 1, 'after': len(admission.points),
                              'point': list(map(str, point))})
    proof = checked_rank(model, admission.points, admission.primes,
                         initial['proof']['no_rational_2_torsion_prime'])
    result = {
        'status': 'PASS_REPLAYED_SEED_CLOUD', 'curve': initial['curve'],
        'points': [list(map(str, p)) for p in admission.points], 'proof': proof,
        'rank_lower_bound': len(admission.points), 'generic_rank': 17,
        'reconciliation_gains': gains, 'bindings': bindings, 'point_searches': 0,
    }
    out = case / 'seed-reconciled.json'
    if out.exists():
        require(read(out) == json.loads(json.dumps(result)), 'independent seed-cloud reconciliation differs')
    else:
        atomic(out, result, immutable=True)


def _mask_from_winning_seed(case):
    terminal = read(case / 'seed-search/terminal.json')
    if terminal['rank_lower_bound'] < 18:
        return None
    chart = read(case / 'seed-search' / f"chart-{terminal['charts'] - 1:04d}.json")
    word = chart['mapping']['centre']['representative']
    return sum((int(v) % 2) << j for j, v in enumerate(word[:17]))


def _insert(mask, pivots):
    q = int(mask)
    for k in sorted(pivots, reverse=True):
        if (q >> k) & 1:
            q ^= pivots[k]
    if q:
        pivots[q.bit_length() - 1] = q
    return bool(q)


def _exact_rows(family, masks, allow_failures=False):
    from sage.all import ZZ, matrix, pari
    import numpy as np
    from visibility_lattice_fast import IntegerExactParity
    from visibility_lattice_v2 import ExactParity

    g = 2 * matrix(ZZ, family['gram'])
    u = matrix(ZZ, pari(g).qflllgram()).transpose(); inv = u.inverse()
    require(abs(u.det()) == 1 and inv.denominator() == 1, 'generic LLL not unimodular')
    reduced = u * g * u.transpose()
    fast, ref = IntegerExactParity(reduced.rows()), ExactParity(reduced.rows())
    rows, checks = [], []
    for mask in masks:
        w = matrix(ZZ, 1, 17, [(int(mask) >> j) & 1 for j in range(17)])
        residue = tuple(int(x) % 2 for x in (w * inv).row(0))
        starts, _ = fast.babai(np.asarray([residue], dtype=np.int64)); start = tuple(map(int, starts[0]))
        try:
            proof = fast.solve(residue, start, 2000000)
            require(proof == ref.solve(residue, start, 2000000), 'generic exact CVP solvers differ')
        except BaseException:
            if allow_failures:
                continue
            raise
        words = []
        for v in proof['minima']:
            z = tuple(map(int, (matrix(ZZ, 1, 17, v) * u).row(0)))
            if next(x for x in z if x) < 0:
                z = tuple(-x for x in z)
            require(sum((x % 2) << j for j, x in enumerate(z)) == int(mask), 'parity transport differs')
            if z not in words:
                words.append(z)
        rows.append({'mask': int(mask), 'word': list(min(words)), 'alternatives': [list(x) for x in sorted(words)],
                     'norm': proof['norm']})
        checks.append({'mask': int(mask), 'reduced_seed': start, 'proof': proof})
    return g, u, rows, checks


def prepare_bank(case, packet_path, output, generation):
    """Build a fresh <=16-parent bank. First exhaust exact maximum classes, then sampled shells."""
    import random

    row = candidate(case); packet = read(packet_path)
    require(packet['rank_lower_bound'] == len(packet['points']) and packet['rank_lower_bound'] >= 18,
            'certified seeded subgroup required')
    family = next(x for x in read(CATALOGUE)['families'] if x['family'] == row['family'])
    maximum = sorted(int(x['mask']) for x in family['classes'])
    prior = []
    for path in sorted((case / 'banks').glob('bank-*/preparation/anchor-bank.json')):
        prior.extend(int(x['mask']) for x in read(path)['rows'])
    used = set(prior)
    winning = _mask_from_winning_seed(case)

    if generation <= 2:
        pool = [m for m in maximum if m not in used and (generation or m != winning)]
    else:
        # Deterministic broad parity sample. Scores/ranks never enter this construction.
        rng = random.Random(int(hashlib.sha256(f'{row["family"]}/{generation}'.encode()).hexdigest(), 16))
        sample = set()
        while len(sample) < 256:
            q = rng.randrange(1, 1 << 17)
            if q not in used:
                sample.add(q)
        pool = sorted(sample)
    require(pool, 'no fresh parent classes in this bank generation')
    g, u, exact, checks = _exact_rows(family, pool, allow_failures=generation >= 3)
    require(len(exact) >= min(16, len(pool)), 'too few exact sampled parent classes')
    if generation >= 3:
        exact.sort(key=lambda r: (-r['norm'], r['mask']))
    else:
        exact.sort(key=lambda r: r['mask'])

    pivots = {}
    if winning is not None:
        _insert(winning, pivots)
    chosen = []
    for r in exact:
        if _insert(r['mask'], pivots):
            chosen.append(r)
        if len(chosen) == 16:
            break
    for r in exact:
        if len(chosen) == 16:
            break
        if r not in chosen:
            chosen.append(r)
    require(chosen and len({r['mask'] for r in chosen}) == len(chosen) <= 16, 'invalid parent subset')

    output.mkdir(parents=True, exist_ok=False)
    protocol = {
        'schema': 'autonomous-r17-parent-bank.v1', 'family': row['family'],
        'parameter': row['parameter'], 'generation': int(generation),
        'initial_rank': packet['rank_lower_bound'], 'generic_rank': 17, 'point_searches': 0,
        'selection': 'Exact maximum classes are partitioned before any sampled-shell bank. Within a bank, independent masks are preferred then filled to <=16. Generation>=3 uses a deterministic 256-mask broad sample ranked by exact generic CVP norm. No specialized point outcome selects a mask.',
        'inputs': {str(p.relative_to(ROOT)): sha(p) for p in
                   (case / 'candidate.json', packet_path, CATALOGUE, Path(__file__))},
    }
    atomic(output / 'protocol.json', protocol, immutable=True)
    seed = {k: packet[k] for k in ('curve', 'points', 'proof', 'rank_lower_bound')}; seed['generic_rank'] = 17
    atomic(output / f"seed-M{packet['rank_lower_bound']}.json", seed, immutable=True)
    atomic(output / 'derivation.json', {
        'generation': int(generation), 'winning_seed_mask': winning, 'prior_masks': sorted(used),
        'candidate_count': len(exact), 'new_masks': [r['mask'] for r in chosen],
        'mode': 'exact-maximum-partition' if generation <= 2 else 'deterministic-sampled-shell',
    }, immutable=True)
    atomic(output / 'generic-cvp-proofs.json', {
        'gram': [list(map(int, r)) for r in g.rows()], 'LLL': [list(map(int, r)) for r in u.rows()],
        'checks': [c for c in checks if c['mask'] in {r['mask'] for r in chosen}],
    }, immutable=True)
    bank = {
        'status': 'COMPLETE_FROZEN_COMPLEMENT_PARENT_SUBSET', 'dimension': 17,
        'generic_gram_scale': 2,
        'rows': [{k: r[k] for k in ('mask', 'word', 'norm')} for r in chosen],
        'shells': sorted({r['norm'] for r in chosen}),
        'protocol_sha256': sha(output / 'protocol.json'),
        'derivation_sha256': sha(output / 'derivation.json'),
        'generic_cvp_sha256': sha(output / 'generic-cvp-proofs.json'),
    }
    atomic(output / 'anchor-bank.json', bank, immutable=True)
    files = {p.name: sha(p) for p in sorted(output.glob('*.json'))}
    atomic(output / 'prepared.json', {'status': 'PASS_AUTONOMOUS_PARENT_PREPARATION',
                                      'point_searches': 0, 'files': files}, immutable=True)


def reseed(case, preparation, packet_path, output):
    """Retain one exact parent bank but replace its certified subgroup seed."""
    from memory_rank_certificate import checked_rank

    sealed = read(preparation / 'prepared.json')
    require(all(sha(preparation / name) == digest for name, digest in sealed['files'].items()),
            'old preparation changed')
    prior_path = next(preparation.glob('seed-M*.json')); prior = read(prior_path); packet = read(packet_path)
    require(packet['curve'] == prior['curve'] and packet['points'][:len(prior['points'])] == prior['points'],
            'reseed packet does not extend prepared subgroup')
    points = tuple(tuple(map(F, p)) for p in packet['points']); proof = packet['proof']
    fresh = checked_rank(tuple(map(F, packet['curve'])), points,
                         [r['prime'] for r in proof['signatures']], proof['no_rational_2_torsion_prime'])
    require(json.loads(json.dumps(fresh)) == proof and len(points) == packet['rank_lower_bound'],
            'reseed certificate differs')
    output.mkdir(parents=True, exist_ok=False)
    for name in sealed['files']:
        if name in ('protocol.json', 'anchor-bank.json') or name.startswith('seed-M'):
            continue
        (output / name).write_bytes((preparation / name).read_bytes())
    protocol = read(preparation / 'protocol.json')
    protocol['initial_rank'] = len(points)
    protocol['reseed_inputs'] = {str(p.relative_to(ROOT)): sha(p) for p in
                                 (packet_path, preparation / 'prepared.json', Path(__file__))}
    protocol['reseed_rule'] = 'Exact certified subgroup extension; retain identical frozen generic parent bank.'
    atomic(output / 'protocol.json', protocol, immutable=True)
    bank = read(preparation / 'anchor-bank.json'); bank['protocol_sha256'] = sha(output / 'protocol.json')
    atomic(output / 'anchor-bank.json', bank, immutable=True)
    seed = {k: packet[k] for k in ('curve', 'points', 'proof', 'rank_lower_bound')}; seed['generic_rank'] = 17
    atomic(output / f"seed-M{len(points)}.json", seed, immutable=True)
    files = {p.name: sha(p) for p in sorted(output.glob('*.json'))}
    atomic(output / 'prepared.json', {'status': 'PASS_AUTONOMOUS_RESEED_PREPARATION',
                                      'point_searches': 0, 'files': files}, immutable=True)


def v3(preparation, folder, max_calls, replay=False):
    import run_complement_seed_v3 as runner
    runner.PREP = runner.BANK = preparation
    if not replay and not folder.exists():
        policy = runner.freeze(folder)
        require(0 < max_calls <= 100, 'segment allowance must be 1..100')
        policy['max_charts'] = int(max_calls)
        policy['autonomous_controller'] = {'segment_allowance': int(max_calls),
            'adapter_sha256': sha(Path(__file__))}
        atomic(folder / 'protocol.json', policy)
    require(read(folder / 'protocol.json')['max_charts'] == int(max_calls), 'segment allowance changed')
    (runner.replay if replay else runner.search)(folder)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('phase', choices=('generic-prepare', 'generic-replay', 'seed-cloud', 'bank', 'reseed', 'v3-search', 'v3-replay'))
    p.add_argument('--case', type=Path, required=True)
    p.add_argument('--packet', type=Path); p.add_argument('--output', type=Path); p.add_argument('--preparation', type=Path)
    p.add_argument('--generation', type=int, default=0); p.add_argument('--max-calls', type=int, default=100)
    a = p.parse_args(); case = a.case.resolve()
    if a.phase == 'generic-prepare': generic(case, False)
    elif a.phase == 'generic-replay': generic(case, True)
    elif a.phase == 'seed-cloud': reconcile_seed(case)
    elif a.phase == 'bank': prepare_bank(case, a.packet.resolve(), a.output.resolve(), a.generation)
    elif a.phase == 'reseed': reseed(case, a.preparation.resolve(), a.packet.resolve(), a.output.resolve())
    elif a.phase == 'v3-search': v3(a.packet.resolve(), a.output.resolve(), a.max_calls, False)
    else: v3(a.packet.resolve(), a.output.resolve(), a.max_calls, True)


if __name__ == '__main__':
    main()
