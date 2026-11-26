"""Finite conditional k=2 test of Pouly's first-excluded-prefix bound.

Run: python3 research/arithmetic/verify_semigroup_k2_stopped_expansion.py
The fixed target and eight events are a diagnostic, not an asymptotic theorem.
"""

from itertools import combinations, product
from math import comb, prod, sqrt

target = 100
cutoff = 5
exposed_selected = {2, 3, 5}
lefts = (1, 2, 3, 5, 7, 9, 11, 12)


def factors(n):
    ans = set()
    p = 2
    while p * p <= n:
        if n % p == 0:
            ans.add(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        ans.add(n)
    return ans


supports = []
for a in lefts:
    support = factors(a) | factors(target - a)
    assert all(p in exposed_selected for p in support if p <= cutoff)
    supports.append(frozenset(p for p in support if p > cutoff))

prime_set = sorted(set().union(*supports))
prob = {p: 1 / sqrt(p) for p in prime_set}
n = len(supports)
indices = tuple(range(n))


def active_components(I):
    unseen = set(I)
    count = 0
    while unseen:
        seed = unseen.pop()
        component = {seed}
        frontier = [seed]
        while frontier:
            i = frontier.pop()
            neighbors = {j for j in unseen if supports[i] & supports[j]}
            unseen -= neighbors
            component |= neighbors
            frontier.extend(neighbors)
        count += len(component) > 1
    return count


def q(I):
    union = set().union(*(supports[i] for i in I)) if I else set()
    return prod(prob[p] for p in union)


subsets = [frozenset(I) for k in range(n + 1) for I in combinations(indices, k)]
retained = {I for I in subsets if active_components(I) <= 1}
assert all(J in retained for I in retained for k in range(len(I) + 1)
           for J in map(frozenset, combinations(I, k)))
boundary = [I for I in subsets if I not in retained and I - {max(I)} in retained]
bound = sum(2 ** len(I) * q(I) for I in boundary)

retained_law = [0.0] * (n + 1)
for I in retained:
    for j in range(len(I) + 1):
        retained_law[j] += q(I) * comb(len(I), j) * (-1) ** (len(I) - j)

actual_law = [0.0] * (n + 1)
for bits in product((0, 1), repeat=len(prime_set)):
    selected = {p for p, b in zip(prime_set, bits) if b}
    weight = prod(prob[p] if p in selected else 1 - prob[p] for p in prime_set)
    z = sum(support <= selected for support in supports)
    actual_law[z] += weight

error = sum(abs(a - b) for a, b in zip(actual_law, retained_law))
pair_overlap_mass = sum(q((i, j)) for i, j in combinations(indices, 2)
                        if supports[i] & supports[j])
mu = sum(q((i,)) for i in indices)
one_copy_second = mu + 2 * sum(q((i, j)) for i, j in combinations(indices, 2))
two_copy_second = mu * mu
assert error <= bound + 1e-10
assert abs(sum(actual_law) - 1) < 1e-12

print('target, cutoff, exposed selected', target, cutoff, sorted(exposed_selected))
print('events', [(a, target-a, sorted(s)) for a, s in zip(lefts, supports)])
print('unexposed primes', prime_set)
print('retained subsets, first-excluded prefixes', len(retained), len(boundary))
print('actual L1 error, Pouly boundary bound', round(error, 8), round(bound, 8))
print('single-copy overlap pair mass', round(pair_overlap_mass, 8))
print('conditional mean, one-copy E[Z^2], two-copy E[Z1 Z2]',
      round(mu, 8), round(one_copy_second, 8), round(two_copy_second, 8))
