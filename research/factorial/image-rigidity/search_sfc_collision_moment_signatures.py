#!/usr/bin/env python3
"""Find integer row blocks close to the collision in normalized log moments."""

import argparse
from functools import cache
from itertools import product
from math import comb, factorial, lgamma, prod

from flint import fmpz_poly, nmod_mpoly_ctx, nmod_poly
import mpmath as mp
import numpy as np


SOURCE_PAIR = ((5, 2, 3), (2, 5, 3))
CORE = ((5, 4, 4), (4, 5, 4), (5, 5, 4))
POWERS = tuple((n-j-k, j, k) for n in (2, 3, 4) for k in range(n+1)
               for j in range((n-k)//2+1))


def dot(row, power):
    return sum(a*b for a, b in zip(row, power))


def normalize(values):
    log = dict(zip(POWERS, values))
    slope = (log[2, 0, 0] - log[0, 0, 2])/2
    return [log[i, j, k] - log[i+j+k, 0, 0] + k*slope for i, j, k in POWERS]


def signature(rows):
    return normalize([sum(lgamma(dot(row, power)+1) for row in rows) for power in POWERS])


def wick(p, r, q):
    return sum(comb(p, h)*comb(r, h)*q**h for h in range(min(p, r)+1))


def physical_logs(q, rows):
    return [sum(mp.loggamma(dot(row, power)+1) for row in rows)
            + mp.log(wick(*(dot(row, power) for row in SOURCE_PAIR), q))
            for power in POWERS]


def reflected(pair):
    return (pair, (pair[1], pair[0], pair[2])) if pair != (0, 0, 0) else ()


def tensor(values):
    entries = dict(zip(POWERS, values))
    return lambda i, j, k: entries[max(i, j), min(i, j), k]


def symmetric_moment(moment, even, odd, third):
    return sum((-1)**(odd-b)*comb(even, a)*comb(odd, b)
               * moment(a+b, even+odd-a-b, third)
               for a in range(even+1) for b in range(odd+1))


def field_moment(moment, even, odd, field):
    return sum(comb(even, k)*field**k*symmetric_moment(moment, even-k, odd, k)
               for k in range(even+1))


def equations(values, field):
    moment = tensor(values)
    variance = field_moment(moment, 0, 2, field)
    s_squared = field_moment(moment, 2, 0, field)/variance
    cubic = field_moment(moment, 3, 0, field) - 3*s_squared*field_moment(moment, 1, 2, field)
    fourth = (field_moment(moment, 4, 0, field) - 6*s_squared*field_moment(moment, 2, 2, field)
              + s_squared**2*field_moment(moment, 0, 4, field))
    return cubic, fourth


def cubic_field(values, target):
    moment = tensor(values)
    second = [comb(2, k)*symmetric_moment(moment, 2-k, 0, k) for k in range(3)]
    mixed = [symmetric_moment(moment, 1-k, 2, k) for k in range(2)]
    variance = symmetric_moment(moment, 0, 2, 0)
    cubic = [comb(3, k)*symmetric_moment(moment, 3-k, 0, k)
             - 3*sum(second[k-j]*mixed[j] for j in range(2) if 0 <= k-j < 3)/variance
             for k in range(4)][::-1]
    roots = np.roots(np.array(cubic, dtype=float))
    seed = min((float(root.real) for root in roots if abs(root.imag)<1e-7), key=lambda r: abs(r-target))
    return mp.findroot(lambda field: mp.polyval(cubic, field), seed)


def exact_chart(rows):
    @cache
    def moment(i, j, k):
        return prod(factorial(dot(row, (i, j, k))) for row in rows)

    def form(i, j):
        return fmpz_poly([comb(i, k)*symmetric_moment(moment, i-k, j, k) for k in range(i+1)])

    variance = form(0, 2)[0]
    fourth = form(0, 4)[0]
    second = form(2, 0)
    cubic = variance*form(3, 0) - 3*second*form(1, 2)
    quartic = variance**2*form(4, 0) - 6*variance*second*form(2, 2) + fourth*second**2
    return cubic.gcd(quartic).degree(), cubic.degree() < 3 and quartic.degree() < 4


def modular_windows(rows, starts, prime):
    terminal = max(starts)+2
    factorials = [1]
    for i in range(1, terminal*max(map(max, rows))+1):
        factorials.append(factorials[-1]*i % prime)

    @cache
    def moment(i, j, k):
        return prod(factorials[dot(row, (i, j, k))] for row in rows) % prime

    context = nmod_mpoly_ctx.get(('y', 't'), prime)
    affine, infinity, pure_difference = {}, {}, {}
    for n in range(min(starts), terminal+1):
        coeffs = {(j//2, k): comb(n, k)*comb(n-k, j)*symmetric_moment(moment, n-j-k, j, k) % prime
                  for k in range(n+1) for j in range(0, n-k+1, 2)}
        affine[n] = context.from_dict(coeffs)
        infinity[n] = nmod_poly([coeffs[j, n-2*j] for j in range(n//2+1)], prime)
        pure_difference[n] = coeffs[n//2, 0] if n % 2 == 0 else 0

    survivors = []
    for start in starts:
        orders = range(start, start+3)
        anchor, left, right = (2, 1, 3) if start == 1 else orders
        left_eliminant = affine[anchor].resultant(affine[left], 'y')
        right_eliminant = affine[anchor].resultant(affine[right], 'y')
        if (affine[anchor].degrees()[0] < 1
                or not any(pure_difference[n] for n in orders)
                or infinity[start].gcd(infinity[start+1]).gcd(infinity[start+2]).degree() != 0
                or not left_eliminant.gcd(right_eliminant).is_constant()
                or left_eliminant.is_zero() and right_eliminant.is_zero()):
            survivors.append(start)
    return survivors


def closest_blocks(bound, target, keep):
    pairs = [(0, 0, 0)] + [(a, b, c) for a in range(1, bound+1) for b in range(a)
                            for c in range(bound+1)]
    shared = [(d, d, e) for d, e in product(range(bound+1), repeat=2)]
    pair_signatures = np.array([signature(reflected(pair)) for pair in pairs])
    shared_signatures = np.array([signature([row]) for row in shared])
    first, second = np.triu_indices(len(pairs))
    sums = pair_signatures[first] + pair_signatures[second] - target
    norm = np.einsum('ij,ij->i', sums, sums)
    best = []
    for index, shared_signature in enumerate(shared_signatures):
        distances = norm + 2*sums @ shared_signature + shared_signature @ shared_signature
        indices = np.argpartition(distances, keep-1)[:keep]
        best.extend((float(distances[i]), first[i], second[i], index) for i in indices)
        best = sorted(best)[:keep]
    print('Blocks scored:', len(sums)*len(shared), flush=True)
    return [(distance, reflected(pairs[a])+reflected(pairs[b])+(shared[c],))
            for distance, a, b, c in best]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bound', type=int, default=10)
    parser.add_argument('--keep', type=int, default=256)
    parser.add_argument('--max-start', type=int, default=12)
    args = parser.parse_args()
    mp.mp.dps = 55

    def physical(q, field):
        return equations([mp.exp(x) for x in normalize(physical_logs(q, SOURCE_PAIR+CORE))], field)

    start_q = mp.mpf('0.845')
    start_t = -2531*mp.exp((physical_logs(start_q, SOURCE_PAIR+CORE)[POWERS.index((0, 0, 2))]
                           -physical_logs(start_q, SOURCE_PAIR+CORE)[0])/2)
    q, field = mp.findroot(physical, (start_q, start_t))
    print('Collision q, normalized t:', mp.nstr(q, 22), mp.nstr(field, 22), flush=True)
    target = np.array(normalize(physical_logs(q, SOURCE_PAIR)), dtype=float)
    closest = closest_blocks(args.bound, target, args.keep)
    reports = []
    remaining = []
    starts = list(range(1, args.max_start+1))
    starts.remove(2)
    first_prime_survivors = 0
    for distance, block in closest:
        rows = tuple(row for row in block+CORE if any(row))
        logs = normalize([sum(mp.loggamma(dot(row, power)+1) for row in rows) for power in POWERS])
        values = [mp.exp(x) for x in logs]
        t = cubic_field(values, field)
        quartic = equations(values, t)[1]
        gcd_degree, infinity = exact_chart(rows)
        reports.append((distance, t, quartic, gcd_degree, infinity, block))
        if starts:
            survivors = modular_windows(rows, starts, 65521)
            first_prime_survivors += len(survivors)
            if survivors:
                survivors = modular_windows(rows, survivors, 65519)
                remaining.extend((rows, start) for start in survivors)
    for label, report in [('closest log moments', reports[0]),
                          ('smallest absolute quartic', min(reports, key=lambda r: abs(r[2]))),
                          ('closest log moments with negative quartic', next((r for r in reports if r[2]<0), None)),
                          ('closest log moments with positive quartic', next((r for r in reports if r[2]>0), None))]:
        if report is not None:
            distance, t, quartic, gcd_degree, infinity, block = report
            print(label, block, 'norm', np.sqrt(max(0, distance)), 't', mp.nstr(t, 18),
                  'Q4', mp.nstr(quartic, 18), 'gcd', gcd_degree, 'infinity', infinity)
    print('Exact chart checks:', len(reports), 'affine nontrivial:', sum(r[3]>0 for r in reports),
          'infinity survivors:', sum(r[4] for r in reports))
    print('Other starts:', starts, 'checked:', len(reports)*len(starts),
          'first-prime survivors:', first_prime_survivors, 'final survivors:', remaining)


if __name__ == '__main__':
    main()
