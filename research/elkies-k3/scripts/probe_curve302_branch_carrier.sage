#!/usr/bin/env sage-python
"""Frozen two-chart, height512 search on the genus-three branch carrier.

Select charts by coefficients, after exact square scaling; never by hits.
PARI searches are replayed by complete integer-pair enumeration with modular
sieves and integer square tests. Also compute two good-reduction Jacobian
orders and independently replay their Frobenius polynomials by point counts.
One worker, 300 seconds, no unbounded factorization or point search.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import permutations
import json
from math import gcd, isqrt
from pathlib import Path
import signal
from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, gcd as sgcd, lcm, pari, prime_range

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT/'artifacts/generated-results/elkies-k3-curve302-rational-branch-carrier-v1.json'
PROTOCOL = ROOT/'artifacts/generated-results/elkies-k3-curve302-branch-carrier-protocol-v1.json'
OUTPUT = ROOT/'artifacts/generated-results/elkies-k3-curve302-branch-carrier-probe-v1.json'
R = PolynomialRing(QQ, 'z')
z = R.gen()
SMALL_PRIMES = list(prime_range(10000))
HEIGHT = 512


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def square_scale(f):
    assert all(c.denominator() == 1 for c in f)
    content = sgcd([ZZ(c) for c in f])
    remainder, scale = content, ZZ(1)
    for p in SMALL_PRIMES:
        v = remainder.valuation(p)
        if v:
            remainder //= p**v
            scale *= p**(v//2)
    if remainder.is_square():
        scale *= remainder.sqrt()
    reduced = f/scale**2
    assert reduced*scale**2 == f and all(c.denominator() == 1 for c in reduced)
    return reduced, scale


def chart_pool(f, roots):
    f, first_scale = square_scale(f)
    best = {}
    for i, j, k in permutations(range(8), 3):
        ri, rj, rk = [roots[t] for t in [i, j, k]]
        entries = [ri*(rk-rj), rj*(ri-rk), rk-rj, ri-rk]
        denominator = lcm([q.denominator() for q in entries])
        entries = [ZZ(q*denominator) for q in entries]
        content = sgcd(entries)
        entries = [q//content for q in entries]
        if next(v for v in entries if v) < 0:
            entries = [-q for q in entries]
        a, b, c, d = entries
        assert a*d-b*c != 0
        transformed = sum(f[t]*(a*z+b)**t*(c*z+d)**(8-t) for t in range(9))
        reduced, scale = square_scale(transformed)
        assert reduced.degree() == 7
        score = max(abs(ZZ(q)).nbits() for q in reduced)
        row = {'root_triple_infinity_zero_one': [i, j, k],
               'matrix_a_b_c_d': list(map(str, entries)),
               'square_scale': str(first_scale*scale),
               'polynomial_low_to_high': list(map(str, reduced.list())),
               'coefficient_max_bits': int(score)}
        key = tuple(sorted((i, j, k)))
        if key not in best or (score, i, j, k) < (best[key]['coefficient_max_bits'], *best[key]['root_triple_infinity_zero_one']):
            best[key] = row
    assert len(best) == 56
    selected = sorted(best.values(), key=lambda e: (e['coefficient_max_bits'], e['root_triple_infinity_zero_one']))[:2]
    return {'ordered_charts_compared': 336, 'unordered_triples': 56,
            'selection': 'Best coefficient-height chart per unordered root triple; first two distinct triples. No target points enter selection.',
            'original_max_bits': int(max(abs(ZZ(q)).nbits() for q in (f*first_scale**2))),
            'square_scaled_original_max_bits': int(max(abs(ZZ(q)).nbits() for q in f)),
            'initial_square_scale': str(first_scale), 'selected': selected}


def modular_sieve(f):
    rows = []
    for p0 in prime_range(101, 2000):
        p = int(p0)
        fp = f.change_ring(GF(p))
        if fp.degree() != 7 or fp.gcd(fp.derivative()).degree() != 0:
            continue
        squares = {v*v % p for v in range(p)}
        allowed = [int(fp(v)) in squares for v in range(p)]
        inv = [0]+[pow(v, -1, p) for v in range(1, p)]
        rows.append((p, allowed, inv))
        if len(rows) == 8:
            return rows
    raise ArithmeticError('Insufficient good sieve primes in declared bound')


def integer_replay(f, height):
    coefficients = list(map(int, f.list()))
    sieves = modular_sieve(f)
    counts = Counter()
    xs = []
    for denominator in range(1, height+1):
        for numerator in range(-height, height+1):
            if gcd(numerator, denominator) != 1:
                continue
            counts['reduced_pairs'] += 1
            for p, allowed, inv in sieves:
                dm = denominator % p
                if dm and not allowed[(numerator % p)*inv[dm] % p]:
                    counts[f'excluded_mod_{p}'] += 1
                    break
            else:
                counts['exact_square_tests'] += 1
                # Degree-eight homogenization has a square denominator d^8.
                value = sum(a*numerator**i*denominator**(8-i) for i, a in enumerate(coefficients))
                if value >= 0 and isqrt(value)**2 == value:
                    xx = QQ(numerator)/denominator
                    yy = QQ(isqrt(value))/denominator**4
                    assert yy*yy == f(xx)
                    xs.append((xx, yy))
    return xs, dict(sorted(counts.items())), [p for p, _, _ in sieves]


def run_chart(row, height, original, D4=None):
    f = R(row['polynomial_low_to_high'])
    a, b, c, d = map(QQ, row['matrix_a_b_c_d'])
    scale = QQ(row['square_scale'])
    assert sum(original[t]*(a*z+b)**t*(c*z+d)**(8-t) for t in range(9)) == scale**2*f
    found = pari(f).hyperellratpoints([height, height])
    pari_xs = {QQ(str(point[0])) for point in found}
    xs, counts, primes = integer_replay(f, height)
    assert pari_xs == {xx for xx, _ in xs}
    results = []
    for xx, yy in xs:
        assert c*xx+d != 0
        m = (a*xx+b)/(c*xx+d)
        Y = scale*yy/(c*xx+d)**4
        assert Y*Y == original(m)
        point = {'z': str(xx), 'v_nonnegative': str(yy), 'm': str(m), 'Y': str(Y),
                 'degenerate': bool(Y == 0)}
        if D4 is not None:
            point['D4_is_square'] = bool(D4(m).is_square())
            if point['D4_is_square']:
                point['H'] = str(D4(m).sqrt())
        results.append(point)
    # z=infinity is the first designated branch root, hence degenerate.
    assert original(a/c) == 0
    return {'root_triple_infinity_zero_one': row['root_triple_infinity_zero_one'],
            'height_inclusive': height, 'replay_counts': counts, 'sieve_primes': primes,
            'pari_and_integer_replay_agree': True, 'finite_points_up_to_hyperelliptic_sign': results,
            'point_at_infinity': {'m': str(a/c), 'degenerate': True}}


def control():
    f = R(40320)
    for i in range(8):
        f *= z-i
    assert f(8) == 40320**2
    pool = chart_pool(f, list(map(QQ, range(8))))
    results = [run_chart(row, 64, f) for row in pool['selected']]
    assert all(any(p['m'] == '8' and not p['degenerate'] for p in row['finite_points_up_to_hyperelliptic_sign']) for row in results)
    return {'curve': 'y^2=40320*product(x-i,i=0..7)', 'known_point': ['8', '40320'],
            'selection_and_search_shared_with_target': True, 'height_inclusive': 64,
            'results': results, 'boundary': 'Constructed search/transport control, not a calibration of point heights on302.'}


def frobenius(f, a, b):
    tr = PolynomialRing(ZZ, 'T')
    t = tr.gen()
    rows = []
    for p in [47, 53]:
        fp = f.change_ring(GF(p))
        assert fp.degree() == 8 and fp.gcd(fp.derivative()).degree() == 0
        cp = tr(pari(fp).hyperellcharpoly())
        counts = []
        for n in [1, 2, 3]:
            field = GF(p**n, name='w')
            poly = f.change_ring(field)
            # Even degree: there are 1+chi(leading coefficient) points at infinity.
            chi = lambda v: 0 if v == 0 else (1 if v.is_square() else -1)
            count = p**n + sum(chi(poly(v)) for v in field) + 1 + chi(poly.leading_coefficient())
            counts.append(int(count))
        sums = [ZZ(p**n+1-counts[n-1]) for n in [1, 2, 3]]
        c1 = -sums[0]
        c2 = -(sums[1]+c1*sums[0])/2
        c3 = -(sums[2]+c1*sums[1]+c2*sums[0])/3
        independent = tr([p**3, p**2*c1, p*c2, c3, c2, c1, 1])
        assert independent == cp
        elliptic = EllipticCurve(GF(p), [a, b])
        ec = t*t-elliptic.trace_of_frobenius()*t+p
        assert cp.gcd(ec) == 1
        rows.append({'prime': p, 'extension_point_counts_n1_n2_n3': counts,
                     'frobenius_low_to_high': list(map(str, cp.list())),
                     'jacobian_order': str(cp(1)), 'curve302_frobenius_low_to_high': list(map(str, ec.list())),
                     'coprime_to_curve302_frobenius': True})
        print('DIRECT_FROBENIUS_REPLAY', p, counts, flush=True)
    assert sgcd([ZZ(r['jacobian_order']) for r in rows]) == 64
    return {'records': rows, 'jacobian_order_gcd': 64,
            'rational_two_torsion_order_from_eight_rational_branch_points': 64,
            'rational_torsion': '(Z/2)^6',
            'no_curve302_Q_isogeny_factor': True,
            'boundary': 'No Mordell-Weil rank upper bound. The genus-three Jacobian rank remains UNKNOWN.'}


def main(check):
    source = json.loads(SOURCE.read_text())
    for path, h in source['input_sha256'].items():
        assert digest(ROOT/path) == h
    stamps = {str(p.relative_to(ROOT)): digest(p) for p in [Path(__file__), SOURCE]}
    f = R(source['R8_low_to_high'])
    D4 = R(source['D4_low_to_high'])
    roots = list(map(QQ, source['R8_rational_roots']))
    calibration = control()
    print('PASS_SHARED_CHART_CONTROL', flush=True)
    pool = chart_pool(f, roots)
    protocol = {'schema': 'curve302.branch-carrier-protocol.v1', 'input_sha256': stamps,
                'limits': {'target_charts': 2, 'height_inclusive': HEIGHT, 'workers': 1, 'seconds': 300,
                           'square_content_trial_primes_below': 10000, 'sieve_primes': 8},
                'control': calibration, 'chart_pool': pool}
    rendered_protocol = json.dumps(protocol, indent=2, sort_keys=True)+'\n'
    if check:
        assert PROTOCOL.read_text() == rendered_protocol
    else:
        assert not PROTOCOL.exists() and not OUTPUT.exists(), 'Use --check; preserve frozen artifacts'
        PROTOCOL.write_text(rendered_protocol)
    results = []
    for row in pool['selected']:
        result = run_chart(row, HEIGHT, f, D4)
        results.append(result)
        print('TARGET_CHART', row['root_triple_infinity_zero_one'], 'points',len(result['finite_points_up_to_hyperelliptic_sign']), flush=True)
    reduced, _ = square_scale(f)
    arithmetic = frobenius(reduced, *map(QQ, source['short302_coefficients_a_b']))
    data = {'schema': 'curve302.branch-carrier-probe.v1',
            'status': 'COMPLETE_TWO_CHART_PROBE_AND_JACOBIAN_TORSION_CERTIFICATE',
            'input_sha256': {**stamps, str(PROTOCOL.relative_to(ROOT)): digest(PROTOCOL)},
            'results': results, 'jacobian_arithmetic': arithmetic,
            'nondegenerate_points': sum(not p['degenerate'] for r in results for p in r['finite_points_up_to_hyperelliptic_sign']),
            'nondegenerate_lifts_to302': sum(not p['degenerate'] and p['D4_is_square'] for r in results for p in r['finite_points_up_to_hyperelliptic_sign']),
            'boundary': 'Only the complete inclusive height512 boxes in two frozen rational charts, with chart infinities accounted for. No global rational-point exclusion, Jacobian rank upper bound, new parent or full generic MW basis.'}
    rendered = json.dumps(data, indent=2, sort_keys=True)+'\n'
    if check:
        assert OUTPUT.read_text() == rendered
    else:
        OUTPUT.write_text(rendered)
    print('PASS_BRANCH_CARRIER_PROBE', data['nondegenerate_points'], data['nondegenerate_lifts_to302'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    signal.alarm(300)
    main(args.check)
