#!/usr/bin/env python3
"""Exact support for effective discriminant -23 primitive norm sumsets.

Builds actual Brandt correspondences at the joint 2/3 levels and ramified
23 axis level, compares complete coordinate shells, and certifies every
rational/numerical constant in the all-large proof. Infinite classical
Hecke imports and written correspondence/tail proofs remain necessary.
"""
if not __debug__:
    raise RuntimeError("This certificate requires Python assertions enabled.")
from fractions import Fraction as F
from itertools import combinations, product
from math import isqrt
from pathlib import Path
from types import SimpleNamespace
import argparse, hashlib, importlib.util, json

SUPPORT_NAME = "verify_primitive_norm_minus23_split_targets.py"
SUPPORT_SHA256 = "fb60a06b0b7bc1aea011dbdba43dc5cc2932da4cc20f89cbff9dac205b78ac16"

def load_support():
    path=Path(__file__).with_name(SUPPORT_NAME)
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != SUPPORT_SHA256:
        raise RuntimeError("Split-target support source does not match the pinned edition.")
    spec=importlib.util.spec_from_file_location("minus23_split_support",path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def joint_model(z):
    params = [(1, 0), (3, 1), (2, 0)]
    bases = [[tuple((int(i == j) for i in range(4))) for j in range(4)]]
    bases += [z.ideal_basis(n, s) for n, s in params[1:]]
    units = [[(1, 0, 0, 0), (-1, 0, 0, 0), (0, 0, 1, 0), (0, 0, -1, 0)]]
    units += [[tuple((F(v, n) for v in a)) for a in z.left_order_norm_elements(n, s, 1)] for n, s in params[1:]]

    def inI(a, i):
        if any((F(x).denominator != 1 for x in a)):
            return False
        return i == 0 or z.in_ideal(tuple((int(x) for x in a)), *params[i])

    def roots(p, h=7):
        modulus = p
        r = 0
        for _ in range(1, h):
            r = next((r + k * modulus for k in range(p) if ((r + k * modulus) ** 2 - (r + k * modulus) + 6) % (modulus * p) == 0))
            modulus *= p
        assert (r * r - r + 6) % modulus == 0
        return (r, 1 - r)

    def mul(A, B):
        return z.multiply(A, B)

    def inv(A):
        a, b = A[0]
        c, d = A[1]
        det = a * d - b * c
        return [[F(d, det), F(-b, det)], [F(-c, det), F(a, det)]]

    def embedding(a, p):
        r, s = roots(p)
        x, y, u, v = a
        return [[x + y * r, u + v * r], [-u - v * s, x + y * s]]

    def g(i, p):
        n, s = params[i]
        return [[1, 0], [s, p]] if n == p else [[1, 0], [0, 1]]

    def local(a, i, j, p):
        A = mul(mul(inv(g(i, p)), embedding(a, p)), g(j, p))
        assert all((F(x).denominator % p for row in A for x in row)), (a, i, j, p, A)
        return [[int(F(x).numerator * pow(F(x).denominator, -1, p) % p) for x in row] for row in A]

    def lines(p):
        return [(1, s) for s in range(p)] + [(0, 1)]

    def normalize(v, p):
        a, b = [x % p for x in v]
        assert a or b
        return (1, b * pow(a, -1, p) % p) if a else (0, 1)

    def act(A, v, p):
        return normalize([sum((x * y for x, y in zip(row, v))) for row in A], p)

    def image(A, p):
        return next((normalize(v, p) for v in zip(*A) if any(v)))

    def kernel(A, p):
        return next((v for v in lines(p) if all((sum((x * y for x, y in zip(row, v))) % p == 0 for row in A))))

    def inverse_mod(A, p):
        a, b = A[0]
        c, d = A[1]
        det = (a * d - b * c) % p
        return [[d * pow(det, -1, p) % p, -b * pow(det, -1, p) % p], [-c * pow(det, -1, p) % p, a * pow(det, -1, p) % p]]

    def neighbors(i, p):
        found = {}
        for j, (n, _) in enumerate(params):
            bound = p * params[i][0] * n
            rows = z.norm_vectors(bound)
            for u in range(bound + 1):
                for a, b in product(rows[u], rows[bound - u]):
                    alpha = tuple((F(x, n) for x in a + b))
                    if not all((inI(z.qmul(alpha, v), i) for v in bases[j])):
                        continue
                    Ap = local(alpha, i, j, p)
                    ell = image(Ap, p)
                    assert all((sum((x * y for x, y in zip(row, kernel(Ap, p)))) % p == 0 for row in Ap))
                    entry = (j, kernel(Ap, p), local(alpha, i, j, 5 - p), alpha)
                    if ell in found:
                        assert found[ell][0] == j
                    else:
                        found[ell] = entry
        assert set(found) == set(lines(p)), (i, p, found)
        return found
    unitmats = {(i, p): [local(u, i, i, p) for u in units[i]] for i in range(3) for p in (2, 3)}

    def canonical(i, l2, l3):
        return min(((i, act(A, l2, 2), act(B, l3, 3)) for A, B in zip(unitmats[i, 2], unitmats[i, 3])))
    states = sorted(set((canonical(i, l2, l3) for i in range(3) for l2 in lines(2) for l3 in lines(3))))
    assert len(states) == 22
    lookup = {state: k for k, state in enumerate(states)}
    neigh = {(i, p): neighbors(i, p) for i in range(3) for p in (2, 3)}
    flips = {}
    forwards = {}
    for p in (2, 3):
        flip = []
        forward = []
        for i, l2, l3 in states:
            ell = l2 if p == 2 else l3
            j, back, othermat, _ = neigh[i, p][ell]
            other = act(inverse_mod(othermat, 5 - p), l3 if p == 2 else l2, 5 - p)
            target = canonical(j, back, other) if p == 2 else canonical(j, other, back)
            k = lookup[target]
            flip.append([int(v == k) for v in range(22)])
            row = [0] * 22
            for newell in lines(p):
                nextstate = canonical(j, newell, other) if p == 2 else canonical(j, other, newell)
                row[lookup[nextstate]] += 1
            row[k] -= 1
            assert min(row) >= 0 and sum(row) == p
            forward.append(row)
        assert mul(flip, flip) == z.identity(22)
        assert mul(mul(flip, forward), flip) == [list(row) for row in zip(*forward)]
        flips[p] = flip
        forwards[p] = forward
    assert mul(flips[2], flips[3]) == mul(flips[3], flips[2])
    assert mul(forwards[2], forwards[3]) == mul(forwards[3], forwards[2])
    assert mul(forwards[2], flips[3]) == mul(flips[3], forwards[2])
    axes = sorted(set((lookup[canonical(0, u, v)] for u in [(1, 0), (0, 1)] for v in [(1, 0), (0, 1)])))
    assert len(axes) == 2

    def coordinate_count(N, rows):
        total = 0
        for n in range(N + 1):
            for x, y in rows[n]:
                for u, v in rows[N - n]:
                    if all((q % 2 == 0 for q in (x, y, u, v))) or all((q % 3 == 0 for q in (x, y, u, v))):
                        continue
                    if not (x % 2 == y % 2 == 0 or u % 2 == v % 2 == 0):
                        continue
                    if not (x % 3 == y % 3 == 0 or u % 3 == v % 3 == 0):
                        continue
                    total += 1
        return total
    rows = z.norm_vectors(2 ** 6 * 3 ** 4)
    powers = {p: [z.identity(22)] for p in (2, 3)}
    for p in (2, 3):
        for _ in range(6):
            powers[p].append(mul(powers[p][-1], forwards[p]))
    controls = []
    for a in range(1, 7):
        for b in range(1, 5):
            V = mul(mul(powers[2][a - 1], flips[2]), mul(powers[3][b - 1], flips[3]))
            expected = 8 * sum((V[i][j] for i in axes for j in axes))
            actual = coordinate_count(2 ** a * 3 ** b, rows)
            assert expected == actual, (a, b, expected, actual)
            controls.append({'a': a, 'b': b, 'N': 2 ** a * 3 ** b, 'predicted': expected, 'actual': actual})
    result = {'scope': 'Exploratory exact joint Iwahori2x3 return check, not the complete all-prime sieve.', 'states': states, 'axis_states': axes, 'uniform_weight': '1/22', 'flip2': flips[2], 'flip3': flips[3], 'forward2': forwards[2], 'forward3': forwards[3], 'coordinate_controls': controls}
    return SimpleNamespace(**locals())

def mixed_model(j):
    z = j.z

    def choices(p, div):
        return j.lines(p) if div else list(combinations(j.lines(p), 2))

    def action(A, v, p, div):
        return j.act(A, v, p) if div else tuple(sorted((j.act(A, w, p) for w in v)))

    def canon(i, v2, v3, mode):
        return min(((i, action(A, v2, 2, mode[0]), action(B, v3, 3, mode[1])) for A, B in zip(j.unitmats[i, 2], j.unitmats[i, 3])))

    def state_data(mode):
        raw = [(i, v2, v3) for i in range(3) for v2 in choices(2, mode[0]) for v3 in choices(3, mode[1])]
        states = sorted(set((canon(*v, mode) for v in raw)))
        lookup = {v: k for k, v in enumerate(states)}
        c = len(choices(2, mode[0])) * len(choices(3, mode[1]))
        coarse = [F(3, 11), F(2, 11), F(6, 11)]
        weights = [coarse[i] * sum((canon(*v, mode) == s for v in raw if v[0] == i)) / c for s in states for i in [s[0]]]
        assert sum(weights) == 1
        axes = ((1, 0), (0, 1))
        p2 = (1, 0) if mode[0] else tuple(sorted(axes))
        p3 = (1, 0) if mode[1] else tuple(sorted(axes))
        principal = lookup[canon(0, p2, p3, mode)]
        return (states, lookup, weights, principal)

    def joint_coordinates(N, rows):
        total = 0
        for n in range(N + 1):
            for x, y in rows[n]:
                for u, v in rows[N - n]:
                    if all((q % 2 == 0 for q in (x, y, u, v))) or all((q % 3 == 0 for q in (x, y, u, v))):
                        continue
                    if (x % 2 == y % 2 == 0 or u % 2 == v % 2 == 0) and (x % 3 == y % 3 == 0 or u % 3 == v % 3 == 0):
                        total += 1
        return total
    results = {}
    for p in (2, 3):
        mode = (p == 2, p == 3)
        states, lookup, weights, principal = state_data(mode)
        size = len(states)
        flip = []
        B = []
        for i, v2, v3 in states:
            ell = v2 if p == 2 else v3
            vertex, back, othermat, _ = j.neigh[i, p][ell]
            other = action(j.inverse_mod(othermat, 5 - p), v3 if p == 2 else v2, 5 - p, False)
            target = canon(vertex, back, other, mode) if p == 2 else canon(vertex, other, back, mode)
            k = lookup[target]
            flip.append([int(t == k) for t in range(size)])
            row = [0] * size
            for v in j.lines(p):
                state = canon(vertex, v, other, mode) if p == 2 else canon(vertex, other, v, mode)
                row[lookup[state]] += 1
            row[k] -= 1
            assert min(row) >= 0 and sum(row) == p
            B.append(row)
        assert z.multiply(flip, flip) == z.identity(size)
        adj = z.multiply(z.multiply(flip, B), flip)
        assert all((weights[i] * adj[i][k] == weights[k] * B[k][i] for i in range(size) for k in range(size)))
        expected_mass = F(1, 33) if p == 2 else F(1, 22)
        assert weights[principal] == expected_mass
        maximum = 10 if p == 2 else 7
        rows = z.norm_vectors(p ** maximum)
        power = z.identity(size)
        controls = []
        for t in range(1, maximum + 1):
            V = z.multiply(power, flip)
            predicted = 8 * V[principal][principal]
            actual = joint_coordinates(p ** t, rows)
            assert predicted == actual, (p, t, predicted, actual)
            controls.append({'t': t, 'N': p ** t, 'predicted': predicted, 'actual': actual})
            power = z.multiply(power, B)
        results[str(p)] = {'states': states, 'weights': [str(w) for w in weights], 'principal': principal, 'principal_mass': str(expected_mass), 'flip': flip, 'forward': B, 'controls': controls}

    def embed5(a):
        x, y, u, v = a
        W = [[0, -6], [1, 1]]
        C = [[1, 1], [0, -1]]
        A = [[x + W[0][0] * y, W[0][1] * y], [W[1][0] * y, x + W[1][1] * y]]
        E = [[u + W[0][0] * v, W[0][1] * v], [W[1][0] * v, u + W[1][1] * v]]
        J = [[2 * C[r][s] for s in range(2)] for r in range(2)]
        BJ = z.multiply(E, J)
        return [[int(F(A[r][s] + BJ[r][s]).numerator * pow(F(A[r][s] + BJ[r][s]).denominator, -1, 5) % 5) for s in range(2)] for r in range(2)]
    neigh5 = {}
    for i in range(3):
        found = {}
        for k, (n, _) in enumerate(j.params):
            bound = 5 * j.params[i][0] * n
            rows = z.norm_vectors(bound)
            for u in range(bound + 1):
                for a, b in product(rows[u], rows[bound - u]):
                    alpha = tuple((F(x, n) for x in a + b))
                    if not all((j.inI(z.qmul(alpha, v), i) for v in j.bases[k])):
                        continue
                    ell = j.image(embed5(alpha), 5)
                    if ell in found:
                        assert found[ell][0] == k
                    else:
                        found[ell] = (k, alpha)
        assert set(found) == set(j.lines(5))
        neigh5[i] = found
    mode = (False, False)
    states, lookup, weights, principal = state_data(mode)
    size = len(states)
    T5 = []
    for i, v2, v3 in states:
        row = [0] * size
        for k, alpha in neigh5[i].values():
            v2new = action(j.inverse_mod(j.local(alpha, i, k, 2), 2), v2, 2, False)
            v3new = action(j.inverse_mod(j.local(alpha, i, k, 3), 3), v3, 3, False)
            row[lookup[canon(k, v2new, v3new, mode)]] += 1
        assert sum(row) == 6
        T5.append(row)
    assert all((weights[i] * T5[i][k] == weights[k] * T5[k][i] for i in range(size) for k in range(size)))
    assert weights[principal] == F(1, 66)
    rows = z.norm_vectors(5 ** 4)
    ops = [z.identity(size), T5]
    controls = []
    for t in range(1, 5):
        if t >= 2:
            product_matrix = z.multiply(T5, ops[-1])
            previous = ops[-2]
            ops.append([[product_matrix[i][k] - 5 * previous[i][k] for k in range(size)] for i in range(size)])
        assert all((x >= 0 for row in ops[t] for x in row))
        predicted = 4 * ops[t][principal][principal]
        actual = joint_coordinates(5 ** t, rows)
        assert predicted == actual, (t, predicted, actual)
        controls.append({'t': t, 'N': 5 ** t, 'predicted': predicted, 'actual': actual})
    results['both_unit'] = {'states': states, 'weights': [str(w) for w in weights], 'principal': principal, 'principal_mass': '1/66', 'T5': T5, 'controls': controls}
    results['scope'] = 'Exact actual mixed and both-unit Cartan-normalizer return checks, finite only.'
    return SimpleNamespace(**locals())

def ramified_model(j):
    z = j.z
    p = 23

    def reduce_fraction(x):
        x = F(x)
        assert x.denominator % 23
        return x.numerator * pow(x.denominator, -1, 23) % 23

    def residue(q):
        x, y, u, v = q
        return (reduce_fraction(x + 12 * y), reduce_fraction(u + 12 * v))

    def cmul(a, b):
        return ((a[0] * b[0] - a[1] * b[1]) % 23, (a[0] * b[1] + a[1] * b[0]) % 23)

    def cinv(a):
        n = (a[0] * a[0] + a[1] * a[1]) % 23
        assert n
        return (a[0] * pow(n, -1, 23) % 23, -a[1] * pow(n, -1, 23) % 23)

    def quotient(a):
        u = j.normalize(a, 23)
        v = j.normalize((-a[1], a[0]), 23)
        return min(u, v)
    localpoints = sorted(set((quotient(v) for v in j.lines(23))))
    assert len(localpoints) == 12

    def action(u, v):
        return quotient(cmul(u, v))

    def canonical(i, v):
        return min(((i, action(residue(u), v)) for u in j.units[i]))
    raw = [(i, v) for i in range(3) for v in localpoints]
    states = sorted(set((canonical(i, v) for i, v in raw)))
    lookup = {v: k for k, v in enumerate(states)}
    coarse = [F(3, 11), F(2, 11), F(6, 11)]
    weights = [coarse[i] * sum((canonical(*v) == s for v in raw if v[0] == i)) / 12 for s in states for i in [s[0]]]
    assert sum(weights) == 1 and len(states) == 28
    principal = lookup[canonical(0, (1, 0))]
    assert weights[principal] == F(1, 44)
    character = [1 if pow((v[0] ** 2 + v[1] ** 2) % 23, 11, 23) == 1 else -1 for i, v in states]
    assert sum((w * c for w, c in zip(weights, character))) == 0
    assert character[principal] == 1
    operators = {}
    for ell in (2, 3):
        rows = []
        for i, v in states:
            row = [0] * len(states)
            for k, _, _, alpha in j.neigh[i, ell].values():
                point = canonical(k, action(cinv(residue(alpha)), v))
                row[lookup[point]] += 1
            assert sum(row) == ell + 1
            rows.append(row)
        assert all((weights[i] * rows[i][k] == weights[k] * rows[k][i] for i in range(28) for k in range(28)))
        assert all((sum((x * y for x, y in zip(row, character))) == (ell + 1) * character[i] for i, row in enumerate(rows)))
        operators[ell] = rows
    assert z.multiply(operators[2], operators[3]) == z.multiply(operators[3], operators[2])

    def hecke_powers(ell, maximum):
        ops = [z.identity(28), operators[ell]]
        for t in range(2, maximum + 1):
            product = z.multiply(operators[ell], ops[-1])
            ops.append([[product[i][k] - ell * ops[-2][i][k] for k in range(28)] for i in range(28)])
        return ops
    powers2 = hecke_powers(2, 6)
    powers3 = hecke_powers(3, 4)

    def primitive(ops, t):
        if t < 2:
            return ops[t]
        return [[ops[t][i][k] - ops[t - 2][i][k] for k in range(28)] for i in range(28)]
    vectors = z.norm_vectors(2 ** 6 * 3 ** 4)

    def coordinates(N):
        count = 0
        for n in range(N + 1):
            for x, y in vectors[n]:
                for u, v in vectors[N - n]:
                    if all((t % 2 == 0 for t in (x, y, u, v))) or all((t % 3 == 0 for t in (x, y, u, v))):
                        continue
                    if n % 23 == 0 or (N - n) % 23 == 0:
                        count += 1
        return count
    controls = []
    for a in range(7):
        for b in range(5):
            A = z.multiply(primitive(powers2, a), primitive(powers3, b))
            predicted = 4 * A[principal][principal]
            actual = coordinates(2 ** a * 3 ** b)
            assert predicted == actual, (a, b, predicted, actual)
            controls.append({'a': a, 'b': b, 'N': 2 ** a * 3 ** b, 'predicted': predicted, 'actual': actual})
    result = {'scope': 'Finite actual level-23 axis-compact Brandt returns, both norm characters, and simultaneous D2/D3 primitive selection. Not a proof of the unrestricted spectral theorem.', 'states': states, 'weights': [str(w) for w in weights], 'principal': principal, 'principal_mass': '1/44', 'norm_character': character, 'T2': operators[2], 'T3': operators[3], 'controls': controls}
    return SimpleNamespace(**locals())

def ramified_nonsquare_model(r, m):
    z = r.z
    T5 = []
    for i, v in r.states:
        row = [0] * len(r.states)
        for k, alpha in m.neigh5[i].values():
            target = r.canonical(k, r.action(r.cinv(r.residue(alpha)), v))
            row[r.lookup[target]] += 1
        assert sum(row) == 6
        T5.append(row)
    assert all((r.weights[i] * T5[i][k] == r.weights[k] * T5[k][i] for i in range(28) for k in range(28)))
    assert all((sum((a * b for a, b in zip(row, r.character))) == -6 * r.character[i] for i, row in enumerate(T5)))
    ops = [z.identity(28), T5]
    vectors = z.norm_vectors(625)
    controls = []
    for t in range(1, 5):
        if t >= 2:
            product = z.multiply(T5, ops[-1])
            ops.append([[product[i][k] - 5 * ops[-2][i][k] for k in range(28)] for i in range(28)])
        N = 5 ** t
        actual = 0
        for n in range(N + 1):
            if n % 23 and (N - n) % 23:
                continue
            actual += len(vectors[n]) * len(vectors[N - n])
        predicted = 4 * ops[t][r.principal][r.principal]
        assert predicted == actual
        if t % 2:
            assert actual == 0
        controls.append({'t': t, 'N': N, 'chi23': (-1) ** t, 'predicted': predicted, 'actual': actual})
    result = {'scope': 'Additional actual ramified-axis return check covering both norm-character signs; finite only.', 'T5': T5, 'character_eigenvalue': -6, 'controls': controls}
    return result

def primes_to(R):
    flag = bytearray(b'\x01') * (R + 1)
    flag[:2] = b'\x00\x00'
    for p in range(2, isqrt(R) + 1):
        if flag[p]:
            flag[p * p:R + 1:p] = b'\x00' * ((R - p * p) // p + 1)
    return [p for p in range(2, R + 1) if flag[p]]

def divisor_bound():
    k = 7
    maximum = 4 * k
    C = 1424
    value = F(1)
    packet = []
    ratio = F(maximum + 2, maximum + 1) ** k
    assert ratio < 2
    for p in primes_to(127):
        best, a = max(((F((a + 1) ** k, p ** a), a) for a in range(maximum + 1)))
        value *= best
        packet.append({'p': p, 'max_exponent': a, 'max_ratio': str(best)})
    assert value <= C ** k
    assert value > (C - 1) ** k
    return {'C': C, 'exponent': '1/7', 'prime_threshold': 128, 'last_exponent': maximum, 'decreasing_ratio_at_endpoint': str(ratio), 'Euler_product': str(value), 'maxima': packet, 'scope': 'For p>=128 use a+1<=2^a<=p^(a/7); for smaller p the displayed decreasing adjacent ratio proves the complete maximum. Thus tau(n)<=C*n^(1/7) for every positive n.'}

def mass_bound():
    R = 65536
    den = 2 ** 80
    rounded = 0
    packet = []
    for p in primes_to(R):
        if p < 5 or p == 23:
            continue
        split = pow(p % 23, 11, 23) == 1
        term = F(5 * p - 3, p * p * (p + 1)) if split else F(2, p * (p - 1))
        up = (term.numerator * den + term.denominator - 1) // term.denominator
        assert F(up, den) >= term
        rounded += up
        packet.append({'p': p, 'split': split, 'ceil80': up})
    tail = F(5, 2 * R) + F(5, R * R)
    upper = F(rounded, den) + tail
    assert upper < F(6, 25)
    assert F(5, 12) - F(6, 25) - F(1, 6) == F(1, 100)
    return {'R': R, 'denominator': den, 'finite_ceil_numerator': rounded, 'infinite_tail': str(tail), 'upper': str(upper), 'strict_bound': '6/25', 'conservative_local_survival': '1/100', 'terms': packet, 'scope': 'Exact rational masses rounded upward, plus a complete decreasing-series bound over every remaining odd integer.'}

def ramified_gate():
    controls = []
    for N in range(1, 23):
        points = [(x, y) for x in range(23) for y in range(23) if (x * x + y * y) % 23 == N]
        bad = [(x, y) for x, y in points if x == 0 or y == 0]
        chi = 1 if pow(N, 11, 23) == 1 else -1
        assert len(points) == 24
        assert F(len(bad), len(points)) == F(1 + chi, 12)
        controls.append({'target': N, 'chi23': chi, 'full': len(points), 'bad': len(bad)})
    assert (23 ** 2 - 1) // (2 * (23 - 1)) == 12
    assert pow(2, 11, 23) == pow(3, 11, 23) == 1
    return {'compact_index': 12, 'projective_residue_order': 24, 'axis_subgroup_order': 2, 'norm_characters': ['1', 'chi23'], 'principal_mass': 'w/12', 'character_coefficients': ['w/12', 'w/12'], 'strong_unit_loss_max': '1/6', 'controls': controls, 'scope': 'The residue model Z23[j]+Z23[j]*sqrt(-23) is established in the proof. This enumerates its norm fibers and biased two-character axis gate; scalar Hensel squares supply all principal units.'}

def onset(C):
    assert 240 + 48 == 288
    assert F(5, 12) - F(1, 6) - F(6, 25) == F(1, 100)
    assert F(1, 100) * F(12, 11) == F(3, 275)
    assert F(288, 52800) == F(3, 550)
    threshold = 52800 * C
    assert 2 ** (245 * 3) >= threshold ** 28
    assert 2 ** (249 * 3) >= 23 ** 3 * threshold ** 28
    assert 23 ** 5 < 10 ** 7
    assert 22000 * C * C < 45000000000
    return {'signed_lower': '3N/275 - 288C*N^(25/28)', 'half_main_condition': 'N^(3/28)>=52800C', 'half_main_signed_lower': '3N/550', 'unit_cutoff': '2^245', 'all_admissible_cutoff': '2^249', 'threshold_integer': threshold, 'unit_cutoff_cube_margin': 2 ** 735 - threshold ** 28, 'all_cutoff_cube_margin': 2 ** 747 - 23 ** 3 * threshold ** 28, 'unit_value_pair_lower': '3N^(5/7)/(2200C^2)', 'all_admissible_value_pair_lower': 'N^(5/7)/15000000000', 'scope': 'Exact integer cutoff comparisons and conversion constant; the written proof supplies uniform spectral error and complete tail estimates.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("primitive_norm_minus23_effective_checks.json"))
    args=parser.parse_args()
    z=load_support()
    joint=joint_model(z)
    mixed=mixed_model(joint)
    ramified=ramified_model(joint)
    divisor=divisor_bound()
    result={
        "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "support_checker":SUPPORT_NAME,"support_sha256":SUPPORT_SHA256,
        "joint_divisible_2_3":joint.result,
        "mixed_and_unit_2_3":mixed.results,
        "ramified_23_actual_brandt":ramified.result,
        "ramified_23_nonsquare_actual_brandt":ramified_nonsquare_model(ramified,mixed),
        "divisor_bound":divisor,"mass_bound":mass_bound(),
        "ramified_residue_gate":ramified_gate(),"onset":onset(divisor["C"]),
        "scope":"Exact finite Brandt/Iwahori/Cartan-normalizer returns, ramified norm characters, and complete rational/integer bounds supporting the written effective theorem. These checks do not replace its infinite correspondence, spectral and lattice arguments, or establish external/formal verification."
    }
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"passed":True,"output":str(args.output),"all_admissible_cutoff":"2^249"}))

if __name__=="__main__":main()
