#!/usr/bin/env python3
"""Pure rational/finite arithmetic replay of the unpointed norm cochain."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from pathlib import Path
import retrospective as r
import verify_explicit_governing_octic as old

INPUT = r.OUT/'rank_jump_unpointed_governing_norm_v1.json'
CONTROL = r.OUT/'rank_jump_scalar_cup_control_v1.json'
OUTPUT = r.OUT/'rank_jump_unpointed_governing_norm_verification_v1.json'


class Algebra:
    def __init__(self, modulus):
        self.f = tuple(map(Q, modulus)); self.n = len(modulus)-1
    def elt(self, a):
        a = list(map(Q, a))
        while len(a) > self.n:
            c = a.pop(); offset = len(a)-self.n
            for i in range(self.n): a[offset+i] -= c*self.f[i]
        return tuple(a+[Q(0)]*(self.n-len(a)))
    def add(self, a, b): return tuple(x+y for x, y in zip(a, b))
    def neg(self, a): return tuple(-x for x in a)
    def mul(self, a, b):
        out = [Q(0)]*(2*self.n-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b): out[i+j] += x*y
        return self.elt(out)
    def power(self, a, n):
        out = self.elt([1])
        while n:
            if n & 1: out = self.mul(out, a)
            a = self.mul(a, a); n //= 2
        return out
    def evaluate(self, coeff, a):
        out = self.elt([])
        for c in reversed(coeff): out = self.add(self.mul(out, a), self.elt([c]))
        return out
    def matrix(self, a):
        cols = [self.mul(a, self.elt([0]*i+[1])) for i in range(self.n)]
        return [list(row) for row in zip(*cols)]
    def norm(self, a):
        matrix = self.matrix(a)
        # Small determinant by permutations; used only for the cubic.
        assert self.n == 3
        return sum((-1)**sum(p[i] > p[j] for i in range(3) for j in range(i+1, 3))
                   *matrix[0][p[0]]*matrix[1][p[1]]*matrix[2][p[2]] for p in permutations(range(3)))
    def inverse(self, a):
        rows = [row+[Q(i == 0)] for i, row in enumerate(self.matrix(a))]
        for i in range(self.n):
            j = next(j for j in range(i, self.n) if rows[j][i])
            rows[i], rows[j] = rows[j], rows[i]
            q = rows[i][i]; rows[i] = [x/q for x in rows[i]]
            for j in range(self.n):
                if j != i:
                    q = rows[j][i]; rows[j] = [x-q*y for x, y in zip(rows[j], rows[i])]
        out = tuple(row[-1] for row in rows)
        assert self.mul(a, out) == self.elt([1])
        return out


def symbolic():
    zero = (0,)*6
    def add(a, b):
        out = dict(a)
        for m, c in b.items(): out[m] = out.get(m, 0)+c
        return {m: c for m, c in out.items() if c}
    def scale(a, c): return {m: v*c for m, v in a.items() if v*c}
    def mul(a, b):
        out = {}
        for m, c in a.items():
            for n, d in b.items():
                k = tuple(x+y for x, y in zip(m, n)); out[k] = out.get(k, 0)+c*d
        return {m: c for m, c in out.items() if c}
    def prod(items):
        out = {zero: 1}
        for item in items: out = mul(out, item)
        return out
    variables = [{tuple(int(i == j) for i in range(6)): 1} for j in range(6)]
    xs, bs = variables[:3], variables[3:]
    values = [prod(add(xs[i], scale(bs[i], -1 if v >> i & 1 else 1)) for i in range(3)) for v in (0, 3, 5, 6)]
    elementary = []
    for k in range(1, 5):
        total = {}
        for terms in combinations(values, k): total = add(total, prod(terms))
        elementary.append(total)
    P, B = prod(xs), prod(bs)
    delta = [add(mul(x, x), scale(mul(b, b), -1)) for x, b in zip(xs, bs)]
    D = prod(delta); middle = {}
    for k in range(3):
        term = prod([delta[i] for i in range(3) if i != k]+[add(mul(xs[k], xs[k]), mul(bs[k], bs[k]))])
        middle = add(middle, scale(term, 2))
    expected = [scale(add(P, B), 4), middle, scale(mul(D, add(P, scale(B, -1))), 4), mul(D, D)]
    assert elementary == expected
    return {'variables': 6, 'elementary_symmetric_identities': 4, 'method': 'integer sparse polynomial expansion'}


def group_certificate():
    V = (0, 3, 5, 6); G = list(permutations(range(3))); roots = list(product(V, range(2)))
    def act(g, a): return sum(((a >> i) & 1) << g[i] for i in range(3))
    def dot(a, b): return (a & b).bit_count() % 2
    def compose(g, h): return tuple(g[h[i]] for i in range(3))
    def multiply(x, y):
        e, f, g, k = x; ee, ff, h, l = y
        return e ^ act(g, ee), f ^ act(g, ff), compose(g, h), k ^ l ^ dot(e, act(g, ff))
    def action(x):
        e, f, g, k = x
        return tuple(roots.index((act(g, v) ^ f, s ^ k ^ dot(e, act(g, v)))) for v, s in roots)
    full = list(product(V, V, G, range(2)))
    actions = {x: action(x) for x in full}; assert len(set(actions.values())) == 192
    for x, y in product(full, repeat=2):
        assert actions[multiply(x, y)] == tuple(actions[x][actions[y][i]] for i in range(8))
    cycle = (1, 2, 0); inverse_cycle = (2, 0, 1); I = (0, 1, 2)
    small = [x for x in full if x[2] in (I, cycle, inverse_cycle) and x[1] == act(inverse_cycle, x[0])]
    assert len(small) == 24
    assert all(multiply(x, y) in small for x, y in product(small, repeat=2))
    unit = (0, 0, I, 0); order_counts = Counter(); prime_types = Counter()
    for x in small:
        y = unit; order = 0
        while True:
            y = multiply(y, x); order += 1
            if y == unit: break
            assert order <= 24
        order_counts[order] += 1
        if x[2] == I: continue
        u = next(v for v in V if act(x[2], v) ^ v == x[0])
        psi = dot(u, x[1]) ^ x[3]
        perm = actions[x]; seen = set(); lengths = []
        for i in range(8):
            if i in seen: continue
            j = i; size = 0
            while j not in seen: seen.add(j); size += 1; j = perm[j]
            lengths.append(size)
        typ = tuple(sorted(lengths)); assert typ == ((1, 1, 3, 3) if psi == 0 else (2, 6))
        prime_types[(psi, typ)] += 1
    assert dict(order_counts) == {1: 1, 2: 1, 3: 8, 6: 8, 4: 6}
    assert len({actions[x][0] for x in small}) == 8
    return {'universal_cup_action_checks': 192**2, 'universal_group_order': 192,
            'small_group_order': 24, 'small_order_counts': {str(k): v for k, v in sorted(order_counts.items())},
            'small_inert_cycle_counts': {str(k): v for k, v in sorted(prime_types.items())},
            'small_group_structure': 'Q8 semidirect C3 with nontrivial order-three action; SL(2,3)',
            'small_galois_relation': 'beta_i=alpha_(i+1); f_i=e_(i+1)'}


def modular(row, entry):
    p = entry['prime']; trim, sub, mul, power, gcd = old.finite_polynomials(p)
    f = trim(row['cubic_ascending']); h = trim(list(map(int, row['integral_octic_ascending'])))
    assert len(gcd(sub(power([0, 1], p, f), [0, 1]), f)) == 1
    xpk = [0, 1]; counts = {}; degrees = []
    for k in range(1, 9):
        xpk = power(xpk, p, h); total = len(gcd(sub(xpk, [0, 1]), h))-1
        old_count = sum(d*n for d, n in counts.items() if k % d == 0)
        assert (total-old_count) % k == 0
        counts[k] = (total-old_count)//k; degrees += [k]*counts[k]
    assert degrees == entry['factor_degrees']
    # Since the cubic extension has odd degree, a rational nonsquare suffices.
    q = p**3; odd = q-1; s = 0
    while odd % 2 == 0: odd //= 2; s += 1
    nonsquare = next(z for z in range(2, p) if pow(z, (p-1)//2, p) == p-1)
    b = trim(row['beta_ascending']); t = power(b, odd, f); x = power(b, (odd+1)//2, f)
    c = power([nonsquare], odd, f); m = s
    while t != [1]:
        u = t; i = 0
        while u != [1]: u = mul(u, u, f); i += 1
        assert i < m
        v = power(c, 2**(m-i-1), f); x = mul(x, v, f); c = mul(v, v, f); t = mul(t, c, f); m = i
    assert mul(x, x, f) == b
    norm_exponent = p*p+p+1; norm_beta_root = int(row['beta_norm_square_root']) % p
    if power(x, norm_exponent, f) != [norm_beta_root]: x = trim([-v for v in x])
    assert power(x, norm_exponent, f) == [norm_beta_root]
    xv = trim([r.mod(v, p) for v in row['X_ascending']])
    value = power(sub(xv, [-v for v in x]), norm_exponent, f)
    assert value == [entry['norm_mod_prime']]
    assert int(pow(value[0], (p-1)//2, p) == p-1) == entry['psi']
    return {'prime': p, 'psi': entry['psi']}


def compute():
    row = r.read(INPUT); control = r.read(CONTROL)
    assert row['status'] == control['status'] == 'PASS'
    for path, sha in row['bindings'].items(): assert r.digest((r.ROOT/path).read_bytes()) == sha
    assert control['polynomial_ascending'] == row['cubic_ascending']
    assert [list(map(Q, v['beta'])) for v in control['norm_witnesses']] == [list(map(Q, row[k])) for k in ('alpha_ascending', 'beta_ascending')]
    K = Algebra(row['cubic_ascending'])
    a, b, X, Y = [K.elt(row[k]) for k in ('alpha_ascending', 'beta_ascending', 'X_ascending', 'Y_ascending')]
    assert K.add(K.mul(X, X), K.neg(K.mul(a, K.mul(Y, Y)))) == b
    assert K.norm(a) == K.norm(b) == 625
    tau = K.elt([-2, -12, 1]); theta = K.elt([0, 1])
    assert K.evaluate(row['cubic_ascending'], tau) == K.elt([0])
    assert K.evaluate(tau, K.evaluate(tau, tau)) == theta
    assert K.evaluate(a, tau) == b
    root = K.elt([Q(-29, 25), Q(-23, 25), Q(2, 25)])
    assert K.evaluate(b, tau) == K.mul(K.mul(a, b), K.mul(root, root))
    assert K.mul(K.mul(a, b), root) == K.elt([-25])
    split_roots = [z for z in range(37) if sum(c*z**i for i, c in enumerate(row['cubic_ascending'])) % 37 == 0]
    signatures = [sum(int(pow(sum(int(c)*z**i for i, c in enumerate(v)) % 37, 18, 37) == 36) << j
                      for j, z in enumerate(split_roots)) for v in (a, b)]
    assert len(split_roots) == 3 and r.rank(signatures) == 2
    delta = K.add(K.mul(X, X), K.neg(b)); D = K.norm(delta); P = K.norm(X)
    quotient = K.mul(K.add(K.mul(X, X), b), K.inverse(delta))
    matrix = K.matrix(quotient); trace = sum(matrix[i][i] for i in range(3))
    h = [D*D, 0, -4*D*(P-25), 0, 2*D*trace, 0, -4*(P+25), 0, 1]
    assert list(map(Q, row['rational_octic_ascending'])) == h
    integral = list(map(int, row['integral_octic_ascending']))
    scale = Q(row['integral_root_scale'])
    assert [v*scale**(8-i) for i, v in enumerate(h)] == integral
    assert str(old.discriminant(integral)) == row['integral_octic_discriminant']
    reduced = [16, 0, -112, 0, 92, 0, -18, 0, 1]
    A = Algebra(reduced)
    change = A.elt([0, Q(135, 2), 0, Q(-435, 4), 0, Q(195, 8), 0, Q(-45, 32)])
    assert A.evaluate(h, change) == A.elt([0])
    assert old.discriminant(reduced) == 2**36*163**4
    # The residue-basis powers have rank eight: the map generates the whole algebra.
    powers = [A.power(change, i) for i in range(8)]
    matrix = [list(v) for v in zip(*powers)]
    for i in range(8):
        j = next(j for j in range(i, 8) if matrix[j][i]); matrix[i], matrix[j] = matrix[j], matrix[i]
        c = matrix[i][i]; matrix[i] = [v/c for v in matrix[i]]
        for j in range(i+1, 8):
            c = matrix[j][i]; matrix[j] = [x-c*y for x, y in zip(matrix[j], matrix[i])]
    replays = [modular(row, e) for e in row['inert_prime_table']]
    return {'schema': 'rank-jump.unpointed-governing-norm-verification.v1', 'status': 'PASS',
            'bindings': {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in (INPUT, CONTROL, Path(__file__), Path(old.__file__), Path(r.__file__))},
            'norm_identity': True, 'symbolic': symbolic(), 'group': group_certificate(),
            'cubic_automorphism': ['-2', '-12', '1'], 'conjugate_square_multiplier': list(map(str, root)),
            'split_prime_independence': {'prime': 37, 'roots': split_roots, 'signatures': signatures, 'rank': 2},
            'reduced_octic_ascending': reduced, 'old_root_in_reduced_algebra': list(map(str, change)),
            'reduction_map_generates_degree': 8, 'reduced_octic_discriminant': str(2**36*163**4),
            'finite_ramification_support_contained_in': [2, 163], 'prime_replays': replays,
            'boundary': 'Exact norm, cochain, Galois and ramification construction for class-selected Sha controls; no rational-solubility conclusion.'}


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['build', 'check']); args = p.parse_args()
    result = compute()
    if args.mode == 'build': r.write_new(OUTPUT, result)
    else: assert r.read(OUTPUT) == result
    print('PASS: unpointed norm, universal cup action, order24 small group, reduced octic, and 29 Frobenius replays')
