#!/usr/bin/env python3
"""Exact retrospective arithmetic comparison of the two successful curve-302 V3 seeds.

This audit answers whether recovered-strict-02 and recovered-strict-03 are merely
two visibility representatives of one exceptional direction, or two genuinely
independent arithmetic directions over the generic M17 subgroup.

It consumes only completed immutable amplifier evidence plus the generic/visibility
artifacts.  It performs no point search and makes no prospective rank claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from dataclasses import asdict
from pathlib import Path
from fractions import Fraction as F

CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
LOCAL = ROOT / 'artifacts/local/elliptic-curves'
AMP = LOCAL / 'curve302-seeded-v3-amplifier-v1'
SUMMARY = AMP / 'summary.json'
VISIBILITY = ART / 'curve302_residual_visibility_geometry_v1.json'
OUT = ART / 'curve302_seed_pair_arithmetic_v1.json'
SEEDS = ('recovered-strict-02', 'recovered-strict-03')
PRIME_BOUND = 1000


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def point_tuple(rows):
    return tuple(tuple(F(x) for x in row) for row in rows)


def curve_tuple(row):
    values = tuple(F(x) for x in row)
    if len(values) == 2:
        values = (F(0), F(0), F(0), *values)
    require(len(values) == 5, 'expected Weierstrass model')
    return values


def atomic_immutable(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(value, sort_keys=True, indent=2) + '\n'
    if path.exists():
        require(read(path) == value, 'preserve immutable seed-pair arithmetic audit')
        return
    temp = path.with_name('.' + path.name + '.tmp-' + str(os.getpid()))
    temp.write_text(encoded)
    os.replace(temp, path)


def crop_rank(signatures, columns):
    from mod2_reduction_independence import gf2_rank
    rows = ([row[i] for i in columns] for signature in signatures for row in signature.rows)
    return gf2_rank(rows, len(columns))


def visibility_fingerprint(entry):
    """Retain small descriptive orbit/parity metadata; this is not a proof layer."""
    terms = ('orbit', 'parity', 'shell', 'extension', 'sector', 'strict', 'local', 'class')
    out = {}

    def visit(value, path=''):
        if isinstance(value, dict):
            for key, child in value.items():
                visit(child, f'{path}.{key}' if path else key)
        elif isinstance(value, (str, int, float, bool)) or value is None:
            key = path.rsplit('.', 1)[-1].lower()
            if any(term in key for term in terms):
                out[path] = value
    visit(entry)
    return dict(sorted(out.items()))


def build():
    from mod2_reduction_independence import (
        combined_mod2_rank,
        find_mod2_reduction_certificate,
        find_two_torsion_certificate_prime,
    )

    require(SUMMARY.exists(), 'completed seeded amplifier summary is missing')
    summary = read(SUMMARY)
    require(summary.get('status') == 'COMPLETE_TWO_SEED_AMPLIFIER', 'amplifier is not complete')
    result_by_seed = {row['seed_direction']: row for row in summary['results']}
    require(set(result_by_seed) == set(SEEDS), 'amplifier roster differs from the two frozen seeds')
    for seed in SEEDS:
        row = result_by_seed[seed]
        require(row['initial_rank'] == 18 and row['rank_lower_bound'] == 31,
                f'{seed} did not independently replay 18->31')

    evidence = {}
    seed_rows = {}
    seed_proofs = {}
    for seed in SEEDS:
        folder = AMP / seed
        seed_path = folder / 'seed-input.json'
        proof_path = folder / 'seed-proof.json'
        verified_path = folder / 'seeded-verified.json'
        for path in (seed_path, proof_path, verified_path):
            require(path.exists(), f'missing amplifier evidence: {path}')
            evidence[str(path.relative_to(ROOT))] = sha(path)
        seed_rows[seed] = read(seed_path)
        seed_proofs[seed] = read(proof_path)
        verified = read(verified_path)
        require(verified['status'] == 'PASS_INDEPENDENT_SEEDED_V3_REPLAY' and
                verified['seed_direction'] == seed and verified['rank_lower_bound'] == 31,
                f'{seed} independent replay endpoint differs')

    evidence[str(SUMMARY.relative_to(ROOT))] = sha(SUMMARY)
    if VISIBILITY.exists():
        evidence[str(VISIBILITY.relative_to(ROOT))] = sha(VISIBILITY)

    first, second = (seed_rows[s] for s in SEEDS)
    model = curve_tuple(first['curve'])
    require(curve_tuple(second['curve']) == model, 'seed curves differ')
    p1, p2 = point_tuple(first['points']), point_tuple(second['points'])
    require(len(p1) == len(p2) == 18, 'seed inputs are not M18 bases')
    require(p1[:17] == p2[:17], 'generic M17 prefixes differ')
    generic = p1[:17]
    seed_points = (p1[17], p2[17])
    require(seed_points[0] != seed_points[1], 'two labels resolve to the same rational point')

    # Recheck each sealed M18 certificate before combining the seeds.
    for seed, points in ((SEEDS[0], p1), (SEEDS[1], p2)):
        proof = seed_proofs[seed]
        require(proof['rank_lower_bound'] == 18, f'{seed} seed proof rank differs')
        sigs = find_mod2_reduction_certificate(model, points, prime_bound=PRIME_BOUND)
        require(combined_mod2_rank(sigs, 18) == 18,
                f'{seed} no longer has a deterministic rank-18 mod-2 certificate')

    joint = (*generic, *seed_points)
    signatures = find_mod2_reduction_certificate(model, joint, prime_bound=PRIME_BOUND)
    joint_rank = combined_mod2_rank(signatures, 19)
    require(joint_rank == 19,
            'bounded exact reductions do not certify the two seeds jointly independent; no collapse claim allowed')
    torsion_prime = find_two_torsion_certificate_prime(model, prime_bound=200)

    idx_h = tuple(range(17))
    idx_h1 = (*idx_h, 17)
    idx_h2 = (*idx_h, 18)
    idx_h12 = tuple(range(19))
    profile = {
        'generic_M17': crop_rank(signatures, idx_h),
        'M17_plus_recovered_strict_02': crop_rank(signatures, idx_h1),
        'M17_plus_recovered_strict_03': crop_rank(signatures, idx_h2),
        'M17_plus_both': crop_rank(signatures, idx_h12),
    }
    require(profile == {
        'generic_M17': 17,
        'M17_plus_recovered_strict_02': 18,
        'M17_plus_recovered_strict_03': 18,
        'M17_plus_both': 19,
    }, 'unexpected finite mod-2 rank profile')

    fingerprints = {}
    if VISIBILITY.exists():
        vis = read(VISIBILITY)
        entries = {row['id']: row for row in vis.get('directions', [])}
        require(all(seed in entries for seed in SEEDS), 'visibility metadata lacks a seeded direction')
        fingerprints = {seed: visibility_fingerprint(entries[seed]) for seed in SEEDS}

    result = {
        'schema': 'curve302-seed-pair-arithmetic.v1',
        'status': 'PASS_EXACT_RETROSPECTIVE_SEED_PAIR_AUDIT',
        'seeds': list(SEEDS),
        'generic_rank': 17,
        'joint_rank_lower_bound': 19,
        'finite_mod2_rank_profile': profile,
        'no_rational_2_torsion_prime': torsion_prime,
        'joint_rank_certificate': {
            'rank_lower_bound': 19,
            'signatures': [asdict(s) for s in signatures],
            'argument': ('The 19 displayed columns are independent in a product of exact finite '
                         'E(F_p)/2E(F_p) quotients and E(Q)[2]=0; infinite descent proves '
                         'their Z-linear independence.'),
        },
        'conclusions': {
            'distinct_over_generic_rational_span': True,
            'distinct_mod2_cosets_over_generic_subgroup': True,
            'single_seed_mechanism_not_forced_by_group_relation': True,
            'meaning': ('The two successful amplifier seeds define two independent exceptional '
                        'directions modulo the generic M17 rational span, and their detected '
                        '2-Kummer classes are independent modulo the generic subgroup. Any common '
                        'geometric construction must therefore produce at least two distinct classes.'),
        },
        'visibility_fingerprints_descriptive_only': fingerprints,
        'bindings': dict(sorted(evidence.items())),
        'prime_bound': PRIME_BOUND,
        'claim_boundary': ('Retrospective arithmetic comparison of two known curve-302 seeds. '
                           'This is not a prospective selector, exact-rank upper bound, or new-curve claim.'),
    }
    return result


def report(result):
    p = result['finite_mod2_rank_profile']
    print('CURVE302_SEED_PAIR|status=PASS|joint_rank=19|seeds=' + ','.join(SEEDS))
    print('CURVE302_SEED_PAIR_MOD2|H={}|H+02={}|H+03={}|H+02+03={}'.format(
        p['generic_M17'], p['M17_plus_recovered_strict_02'],
        p['M17_plus_recovered_strict_03'], p['M17_plus_both']))
    print('CURVE302_SEED_PAIR_CONCLUSION|distinct_rational_directions=True|distinct_mod2_classes=True')
    for seed, fp in result.get('visibility_fingerprints_descriptive_only', {}).items():
        print('CURVE302_SEED_PAIR_VISIBILITY|seed={}|{}'.format(seed, json.dumps(fp, sort_keys=True)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    require(not (args.write and args.check), 'choose --write or --check')
    result = build()
    if args.check:
        require(OUT.exists(), 'seed-pair audit artifact is missing')
        require(read(OUT) == result, 'seed-pair audit replay differs')
    elif args.write:
        atomic_immutable(OUT, result)
    report(result)


if __name__ == '__main__':
    main()
