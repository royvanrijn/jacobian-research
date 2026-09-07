#!/usr/bin/env sage-python
"""Exact genus gate for pairs of the 6417 previously frozen line covers.

No new point enumeration. Reconstruct every branch polynomial from the saved
residual line equations, certify its irreducibility with a small prime, and
deduplicate within each pencil. Distinct irreducible quadratics have disjoint
geometric branch sets, so their V4 fibre products have genus one.
"""
import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import signal
from sage.all import PolynomialRing, QQ, prime_range

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'artifacts/generated-results/elkies-k3-curve302-cubic-pencil-overlap-v1.json'
OUTPUT = ROOT / 'artifacts/generated-results/elkies-k3-curve302-shared-branch-gate-v1.json'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def square(q):
    return q >= 0 and isqrt(q.numerator)**2 == q.numerator and isqrt(q.denominator)**2 == q.denominator


def positive_control():
    ring = PolynomialRing(QQ, 'k')
    field = ring.fraction_field()
    k = field.gen()
    r = (k*k-6*k+12)/(k*k-12)
    w = k*r+3-k
    u = (4*r*r-1)/(r*r-1)
    v2 = w/(r*r-1)
    v1 = r*v2
    assert w*w == 12*r*r-3
    assert v1*v1 == u*(u-1) and v2*v2 == u*(u-4)
    assert max(u.numerator().degree(), u.denominator().degree()) == 4
    evaluate = lambda f: f.numerator()(1)/f.denominator()(1)
    assert tuple(evaluate(f) for f in [u, v1, v2]) == (QQ(-25)/24, QQ(35)/24, QQ(-55)/24)
    return {'branch_sets': [[0, 1], [0, 4]], 'union_size': 3,
            'normalization_genus': 0, 'map_degree': 4,
            'r': str(r), 'w': str(w), 'u': str(u), 'v1': str(v1), 'v2': str(v2),
            'rational_point_at_k1': ['-25/24', '35/24', '-55/24'],
            'identities': ['v1^2=u(u-1)', 'v2^2=u(u-4)'],
            'boundary': 'Base-curve positive control only; no elliptic-surface rank claim.'}


def compute():
    source = json.loads(SOURCE.read_text())
    stamps = {str(p.relative_to(ROOT)): digest(p) for p in [Path(__file__), SOURCE]}
    for path, expected in source['inputs'].items():
        assert digest(ROOT/path) == expected, path
        stamps[path] = expected
    assert source['anchor_pencil_count'] == len(source['pencils']) == 31
    primes = list(map(int, prime_range(3, 2000)))
    rows = []
    histogram = Counter()
    cover_count = 0
    for idx, pencil in enumerate(source['pencils']):
        path = ROOT/pencil['checkpoint']
        assert digest(path) == pencil['checkpoint_sha256'], str(path)
        stamps[str(path.relative_to(ROOT))] = pencil['checkpoint_sha256']
        checkpoint = json.loads(path.read_text())
        assert checkpoint['anchors_one_based'] == pencil['anchors_one_based']
        groups = checkpoint['groups']
        assert len(groups) == len(pencil['classes']) == pencil['cover_class_count'] == 207
        distinct = set()
        witnesses = []
        for group, old in zip(groups, pencil['classes']):
            strings = group['branch_polynomial_monic_low_to_high']
            assert sha256(json.dumps(strings).encode()).hexdigest() == old['branch_sha256']
            c, b, a = map(Fraction, strings)
            assert a == 1 and c != 0
            assert (c, b) not in distinct
            distinct.add((c, b))
            assert len(group['covers']) == 1
            assert group['extra_public_indices_one_based'] == old['extra_public_indices_one_based']
            assert [v['basepoint_index'] for v in group['covers']] == old['basepoint_indices']
            for cover in group['covers']:
                c0, b0, a0 = map(Fraction, cover['residual_F0_coefficients_low_to_high'])
                c1, b1, a1 = map(Fraction, cover['residual_F1_coefficients_low_to_high'])
                dc = b0*b0-4*a0*c0
                db = 2*b0*b1-4*(a0*c1+a1*c0)
                da = b1*b1-4*a1*c1
                assert da == Fraction(cover['discriminant_scalar']) != 0
                assert (dc/da, db/da) == (c, b)
                assert c0+b0+a0 == 0 and c1+b1+a1 != 0
                assert dc != 0 and square(dc)
                cover_count += 1
            delta = b*b-4*c
            assert delta != 0
            # If delta were a rational square, its reduced numerator times
            # denominator would be an integer square. A nonsquare residue
            # therefore certifies irreducibility without any factorization.
            product = delta.numerator*delta.denominator
            witness = next(((p, product % p) for p in primes
                            if pow(product % p, (p-1)//2, p) == p-1), None)
            assert witness is not None, 'No certificate within declared prime bound'
            witnesses.append(list(witness))
            histogram[witness[0]] += 1
        pairs = len(groups)*(len(groups)-1)//2
        rows.append({'pencil_zero_based': idx, 'classes': len(groups),
                     'distinct_irreducible_quadratics': len(distinct),
                     'nonsquare_witness_prime_and_residue_in_class_order': witnesses,
                     'disjoint_branch_pairs': pairs, 'genus_zero_pairs': 0,
                     'genus_one_pairs': pairs})
        print('SHARED_BRANCH_GATE', idx+1, len(groups), pairs, flush=True)
    assert cover_count == source['attempted_line_covers'] == 6417
    pair_count = sum(r['disjoint_branch_pairs'] for r in rows)
    assert pair_count == 660951
    return {'schema': 'curve302.shared-branch-gate.v1',
            'status': 'ALL_FROZEN_WITHIN_PENCIL_DISTINCT_PAIRS_HAVE_GENUS_ONE',
            'input_sha256': stamps, 'positive_control': positive_control(),
            'pencil_count': len(rows), 'cover_count': cover_count,
            'pair_count': pair_count, 'genus_zero_pairs': 0,
            'genus_one_pairs': pair_count,
            'nonsquare_prime_histogram': dict(sorted(histogram.items())),
            'records': rows,
            'proof': 'Each saved cover has two finite branch points and no infinity branch. Distinct irreducible monic quadratics over Q have disjoint roots. The geometrically connected degree-four V4 compositum has four branch points of inertia two: 2g-2=-8+4*2=0. No nonconstant map from P1 to this genus-one normalization exists in characteristic zero.',
            'boundary': 'Only pairs within each of the 31 frozen pencils; no comparison across different base maps, no exclusion of other covers or pencils, no exclusion of all Q(t) parents. This certificate does not construct an alternative parent.',
            'replay': 'sage -python elkies-k3/scripts/certify_curve302_shared_branch_gate.sage --check'}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    signal.alarm(180)
    text = json.dumps(compute(), indent=2, sort_keys=True)+'\n'
    if args.check:
        assert OUTPUT.read_text() == text, 'Certificate changed'
    else:
        assert not OUTPUT.exists(), 'Preserve existing certificate; use --check'
        OUTPUT.write_text(text)
    print('PASS_SHARED_BRANCH_GATE', flush=True)
