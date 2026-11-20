#!/usr/bin/env python3
"""Independent complete finite actual-pair replay for the relative-gap cone.

For m>=1 the actual normalized factorial moment is
T_m(t)=sum_k binom(m,k)*(a*m+1)_(d*k)*t^k.
Each new terminal polynomial is constructed from its OWN shape. The preceding
actual polynomial is reused only after it has been constructed at the preceding
terminal. Coefficients use sliding rising products and binomial recurrences,
without factorial tables. Every attempted modular xgcd identity is multiplied
out. A completed finite PASS needs a separate analytical tail for an all-order
theorem. This program makes no all-binomial assertion.
"""
if not __debug__:
    raise SystemExit('Unoptimized Python is required for exact proof checks.')

from array import array
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
from math import comb, factorial, isqrt
from pathlib import Path
import argparse
import json
import struct
import sys
import time
import flint
from flint import nmod_poly

PRIMES = (1000033, 1000037, 1000039, 1000003)
MAX_FACTORIAL = 903084
_inverse_cache = {}
_rising_cache = {}


def is_prime(p):
    return p >= 2 and all(p % q for q in range(2, isqrt(p) + 1))


def domain():
    rows = [(a, d, max(4000, 128 * (a + d)) - 1)
            for d in range(9, 64) for a in range(3, d // 3 + 1)]
    assert len(rows) == len(set(rows)) == 532
    assert sum(last - 1 for a, d, last in rows) == 3694040
    assert max(last for a, d, last in rows) == 10751
    assert max((a + d) * last for a, d, last in rows) == MAX_FACTORIAL
    assert all(a >= 3 and d >= 3 * a and 2 <= last < PRIMES[0]
               for a, d, last in rows)
    for p in PRIMES:
        assert is_prime(p) and p > MAX_FACTORIAL
    return rows


def make_tables(d, limit, p):
    assert 1 <= d <= limit < p and is_prime(p)
    inv = array('I', [0]) * (limit + 1)
    inv[1] = 1
    for j in range(2, limit + 1):
        inv[j] = -(p // j) * inv[p % j] % p
    return inv, make_rising(d, limit, p, inv)


def make_rising(d, limit, p, inv):
    rise = array('I', [0]) * (limit + 1)
    q = 1
    for j in range(1, d + 1):
        q = q * j % p
    rise[1] = q
    for t in range(1, limit - d + 1):
        q = q * (t + d) % p * inv[t] % p
        rise[t + 1] = q
    return rise


def cached_tables(d, p):
    # The inverse table is shared across all supports handled by a worker.
    # Keep only the most recent gap's rising table for each prime.
    if p not in _inverse_cache:
        inv, rise = make_tables(d, MAX_FACTORIAL, p)
        _inverse_cache[p] = inv
        _rising_cache[p] = (d, rise)
    if _rising_cache[p][0] != d:
        _rising_cache[p] = (d, make_rising(d, MAX_FACTORIAL, p,
                                        _inverse_cache[p]))
    return _inverse_cache[p], _rising_cache[p][1]


def coefficients(a, d, m, p, tables):
    assert a >= 1 and d >= 1 and m >= 1 and (a + d) * m < p
    inv, rise = tables
    out = [1]
    value = 1
    first = a * m + 1
    for k in range(m):
        value = value * (m - k) % p * inv[k + 1] % p
        value = value * rise[first + d * k] % p
        assert value
        out.append(value)
    assert len(out) == m + 1 and out[0] == 1 and out[-1]
    return out


def controls():
    p = PRIMES[0]
    count = 0
    for a, d in [(1, 1), (2, 6), (3, 9), (3, 63), (4, 12),
                 (7, 21), (10, 30), (20, 60), (21, 63)]:
        tab = make_tables(d, 6 * (a + d), p)
        for m in range(1, 7):
            actual = coefficients(a, d, m, p, tab)
            literal = [comb(m, k) * factorial(a * m + d * k)
                       // factorial(a * m) % p for k in range(m + 1)]
            assert actual == literal
            count += 1
    # An actual wrong-shape negative control at a=3,d=9,N=2.
    a, d, N = 3, 9, 2
    actual = coefficients(a, d, N, p, make_tables(d, 24, p))
    wrong = [comb(N, k) * factorial(a * (N - 1) + d * k)
             // factorial(a * (N - 1)) % p for k in range(N + 1)]
    assert actual != wrong and actual[1] != wrong[1]
    f = nmod_poly([1, 1], p)
    g = f * f
    q, s, t = f.xgcd(g)
    assert q.degree() == 1 and s * f + t * g == q
    assert (s + 1) * f + t * g != q
    return {'literal_integer_coefficient_vectors': count,
            'wrong_actual_shape_rejected': True,
            'common_factor_detected': True,
            'corrupted_Bezout_identity_rejected': True}


def worker(args):
    a, d, last = args
    started = time.monotonic()
    previous = {}
    counts = {}
    fallbacks = []
    h = sha256()
    bezout_checks = 0
    for N in range(2, last + 1):
        attempts = []
        for p in PRIMES:
            tab = cached_tables(d, p)
            if p not in previous or previous[p][0] != N - 1:
                prev = nmod_poly(coefficients(a, d, N - 1, p, tab), p)
            else:
                prev = previous[p][1]
            current = nmod_poly(coefficients(a, d, N, p, tab), p)
            previous[p] = (N, current)
            assert prev.degree() == N - 1 and current.degree() == N
            assert prev[0] == current[0] == 1 and prev[N - 1] and current[N]
            q, s, t = current.xgcd(prev)
            assert s * current + t * prev == q
            bezout_checks += 1
            attempts.append([p, q.degree()])
            if q.degree() == 0:
                assert q == nmod_poly([1], p)
                counts[str(p)] = counts.get(str(p), 0) + 1
                h.update(struct.pack('<IIII', a, d, N, p))
                if len(attempts) > 1:
                    fallbacks.append({'terminal': N, 'attempts': attempts})
                break
        else:
            raise AssertionError(('Unresolved ACTUAL finite pair', a, d, N,
                                  attempts))
    assert sum(counts.values()) == last - 1
    return {'lower_exponent': a, 'gap': d, 'first_terminal': 2,
            'last_terminal': last, 'tested_cases': last - 1,
            'maximum_factorial_argument': (a + d) * last,
            'verified_full_Bezout_identities': bezout_checks,
            'successful_primes': counts, 'fallbacks': fallbacks,
            'maximum_upper_degree': last, 'maximum_previous_degree': last - 1,
            'case_stream_sha256': h.hexdigest(),
            'seconds': time.monotonic() - started}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--output', type=Path,
                    default=Path(__file__).resolve().with_name(
                        'independent_relative_gap_three_finite_box.json'))
    ap.add_argument('--benchmark-last', type=int)
    args = ap.parse_args()
    assert args.jobs >= 1
    assert args.benchmark_last is None or args.benchmark_last >= 2
    started = time.monotonic()
    coefficient_controls = controls()
    expected = domain()
    full = args.benchmark_last is None
    selected = expected if full else [(a, d, min(last, args.benchmark_last))
                                     for a, d, last in expected]
    rows = []
    args.output.parent.mkdir(parents=True, exist_ok=True)
    running = args.output.with_suffix('.running.json')

    def progress():
        data = {'status': 'RUNNING', 'complete_finite_box': False,
                'requested_complete_box': full,
                'expected_support_count': 532,
                'expected_terminal_pairs': 3694040,
                'completed_support_count': len(rows),
                'completed_terminal_pairs': sum(r['tested_cases'] for r in rows),
                'seconds': time.monotonic() - started,
                'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
                'rows': sorted(rows, key=lambda r: (r['gap'], r['lower_exponent']))}
        temp = running.with_suffix('.tmp')
        temp.write_text(json.dumps(data, indent=2) + '\n')
        temp.replace(running)

    progress()
    with ProcessPoolExecutor(max_workers=args.jobs) as pool:
        futures = {pool.submit(worker, item): item for item in selected}
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            progress()
            print(json.dumps({k: v for k, v in row.items() if k != 'fallbacks'}),
                  flush=True)
    rows.sort(key=lambda r: (r['gap'], r['lower_exponent']))
    assert [(r['lower_exponent'], r['gap'], r['last_terminal'])
            for r in rows] == selected
    total = sum(r['tested_cases'] for r in rows)
    if full:
        assert total == 3694040 and len(rows) == 532
        assert max(r['last_terminal'] for r in rows) == 10751
        assert max(r['maximum_factorial_argument'] for r in rows) == MAX_FACTORIAL
    out = {'status': 'PASS', 'complete_finite_box': full,
           'tested_supports': len(rows), 'tested_terminal_pairs': total,
           'verified_full_Bezout_identities': sum(
               r['verified_full_Bezout_identities'] for r in rows),
           'maximum_factorial_argument_for_full_domain': MAX_FACTORIAL,
           'actual_shapes': ['a*N+1', 'a*(N-1)+1'],
           'coefficient_construction':
               'Separate binomial/sliding-rising recurrences; preceding actual '
               'polynomial reused at next terminal; no factorial tables',
           'primary_prime': PRIMES[0], 'alternative_primes': list(PRIMES[1:]),
           'controls': coefficient_controls, 'rows': rows,
           'source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
           'versions': {'python': sys.version.split()[0],
                        'python_flint': flint.__version__,
                        'FLINT': flint.__FLINT_VERSION__},
           'seconds': time.monotonic() - started,
           'scope': 'Only the completed ACTUAL finite domain. Infinite '
                    'relative-gap conclusions require the separately proved '
                    'uniform analytical tail, sparse overlap and lower-two '
                    'antecedent. No assertion for all ordinary binomials.'}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'complete_finite_box': full,
                      'tested_supports': len(rows),
                      'tested_terminal_pairs': total,
                      'verified_full_Bezout_identities':
                          out['verified_full_Bezout_identities'],
                      'seconds': out['seconds']}), flush=True)


if __name__ == '__main__':
    main()
