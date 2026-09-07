#!/usr/bin/env python3
"""Explicit pair cochain; bounded retrospective arithmetic, no point search."""
import argparse
from collections import Counter
from itertools import permutations, product
from pathlib import Path
import subprocess
import sys
import retrospective as r

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / 'EXPLICIT_GOVERNING_OCTIC_PROTOCOL.json'
INPUT = r.OUT / 'rank_jump_bad_prime_inputs_v1.json'
OUTPUT = r.OUT / 'rank_jump_explicit_governing_octic_v1.json'
CHECKPOINT = r.ROOT / 'artifacts/local/rank-jump-explicit-governing-octic-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes())
            for p in (Path(__file__), PROTOCOL, INPUT, r.INPUT, HERE / 'retrospective.py')}


def finite_certificate():
    V = (0, 3, 5, 6)
    G = list(permutations(range(3)))
    def act(g, v): return sum(((v >> i) & 1) << g[i] for i in range(3))
    def dot(a, b): return (a & b).bit_count() % 2
    def compose(g, h): return tuple(g[h[i]] for i in range(3))
    base = list(product(V, V, G))
    def multiply(x, y):
        e, f, g = x; ee, ff, h = y
        return e ^ act(g, ee), f ^ act(g, ff), compose(g, h)
    def beta(x, y): return dot(x[1], act(x[2], y[0] ^ y[1]))
    def correction(x): return dot(x[0], x[1]) ^ int(x[1] != 0)
    roots = list(product(V, range(2)))
    def root_action(x, k):
        e, f, g = x
        return tuple(roots.index((act(g, v) ^ e ^ f, s ^ k ^ dot(f, act(g, v)))) for v, s in roots)
    def cycles(perm):
        seen = set(); sizes = []
        for i in range(len(perm)):
            if i in seen: continue
            n = 0; j = i
            while j not in seen: seen.add(j); n += 1; j = perm[j]
            sizes.append(n)
        return tuple(sorted(sizes))
    for x, y in product(base, repeat=2):
        xy = multiply(x, y)
        assert beta(x, y) ^ correction(xy) ^ correction(x) ^ correction(y) == dot(x[0], act(x[2], y[1]))
        for k, l in product(range(2), repeat=2):
            a = root_action(x, k); b = root_action(y, l)
            assert tuple(a[b[i]] for i in range(8)) == root_action(xy, k ^ l ^ beta(x, y))
    actions = {root_action(x, k) for x in base for k in range(2)}
    assert len(actions) == 192
    assert all(len({p[i] for p in actions}) == 8 for i in range(8))
    rows = []
    for x in base:
        e, f, g = x
        if any(act(g, v) == v for v in V[1:]): continue
        a = next(v for v in V if act(g, v) ^ v == e)
        for k in range(2):
            psi = dot(a, f) ^ k ^ correction(x)
            typ = cycles(root_action(x, k))
            assert typ == ((1, 1, 3, 3) if psi == 0 else (2, 6))
            rows.append((psi, typ))
    return {'base_order': 96, 'extension_order': 192, 'cup_identity_checks': 96**2,
            'root_action_composition_checks': 4*96**2, 'faithful_transitive_degree': 8,
            'inert_cycle_type_counts': {str(k): v for k, v in sorted(Counter(rows).items())},
            'cocycle': 'beta(x,y)=f_x dot g_x(e_y+f_y)',
            'correction': 'r(x)=e_x dot f_x + [f_x!=0]',
            'gamma': 'kappa+r', 'psi': '((g-1)^-1 e) dot f + gamma'}


def symbolic_certificate():
    from sage.all import PolynomialRing, QQ
    R = PolynomialRing(QQ, names=('a0', 'a1', 'a2', 'c'))
    a = R.gens()[:3]; c = R.gen(3)
    squares = [z*z+c for z in a]
    def add(x, y): return {m: x.get(m, R(0))+y.get(m, R(0)) for m in set(x)|set(y)}
    def scale(x, u): return {m: z*u for m, z in x.items()}
    def mul(x, y):
        out = {}
        for m, v in x.items():
            for n, w in y.items():
                z = v*w
                for i in range(3):
                    if (m & n) >> i & 1: z *= squares[i]
                out[m ^ n] = out.get(m ^ n, R(0))+z
        return out
    F = {0: R(1)}
    for i in range(3): F = mul(F, {0: a[i], 1 << i: R(1)})
    p = {0: a[0]*a[1]*a[2]}; q = {7: R(1)}
    F2 = mul(F, F); F3 = mul(F2, F); F4 = mul(F3, F)
    residual = add(F4, scale(mul(add(p, q), F3), -4))
    residual = add(residual, scale(F2, c*c*(4*sum(z*z for z in a)+6*c)))
    residual = add(residual, scale(mul(add(q, scale(p, -1)), F), -4*c**3))
    residual = add(residual, {0: c**6})
    assert all(z == 0 for z in residual.values())
    return {'free_algebra_rank': 8, 'zero_residual_coefficients': 8,
            'relations': 'b_i^2=a_i^2+c', 'radical': 'F=product(a_i+b_i)',
            'identity': 'F^4-4(p+q)F^3+c^2(4*sum(a_i^2)+6c)F^2-4c^3(q-p)F+c^6=0'}


def compute_case(index):
    from sage.all import QQ, ZZ, PolynomialRing, GF, lcm
    from sage.version import version
    protocol = r.read(PROTOCOL)
    spec = protocol['cases'][index]
    old = r.read(INPUT)['cases'][spec['case_index']]
    source = next(x for x in r.read(r.INPUT)['rows'] if x['id'] == old['id'])
    short, points = r.short(source['model'], source['generic_points']+source['points'])
    assert short == old['short_model']
    d = QQ(old['elliptic_scaling_d'])
    cubic_coeff = list(map(QQ, old['integral_cubic_ascending']))
    B, A, zero, one = cubic_coeff
    assert zero == 0 and one == 1
    assert A == QQ(short[3])*d**4 and B == QQ(short[4])*d**6
    selected = [old['selected_input_indices'][i] for i in spec['basis_indices']]
    scaled = [(QQ(points[i][0])*d**2, QQ(points[i][1])*d**3) for i in selected]
    for x, y in scaled: assert y*y == x**3+A*x+B
    xP, yP = scaled[0]; xQ, yQ = scaled[1]; c = xQ-xP
    assert c and yP and yQ
    blocks = [(p, r.roots_at(str(A), str(B), p)) for p in source['prime_list']]
    assert all(roots is not None for p, roots in blocks)
    model = ['0', '0', '0', str(A), str(B)]
    sigs = [r.point_signature(model, list(map(str, P)), blocks) for P in scaled]
    assert r.rank(sigs) == 2
    generic_sigs = [r.point_signature(model, [str(QQ(points[i][0])*d**2), str(QQ(points[i][1])*d**3)], blocks)
                    for i in old['selected_input_indices'][:old['generic_dimension']]]
    assert r.rank(generic_sigs) == old['generic_dimension']
    relative = r.rank(generic_sigs+sigs)-len(generic_sigs)
    assert relative == (0 if index == 2 else 2)
    R = PolynomialRing(QQ, 'T'); T = R.gen()
    h = T**8-4*(yP+yQ)*T**6+6*c*c*(xP+xQ)*T**4-4*c**3*(yQ-yP)*T**2+c**6
    denominator = lcm(z.denominator() for z in h.list())
    integral = R([h[i]*denominator**(8-i) for i in range(9)])
    assert all(z.denominator() == 1 for z in integral.list())
    discriminant = ZZ(integral.discriminant()); assert discriminant
    disc_cubic = -4*A**3-27*B**2
    assert not disc_cubic.is_square()
    table = []; excluded = []
    for p in r.primes(protocol['prime_bound']):
        values = (denominator, c.numerator(), c.denominator(), yP.numerator(), yP.denominator(),
                  yQ.numerator(), yQ.denominator(), xP.denominator(), xQ.denominator(), disc_cubic)
        if any(z % p == 0 for z in values) or discriminant % p == 0:
            excluded.append({'prime': p, 'reason': 'explicit denominator, radical, or polynomial-discriminant support'}); continue
        Rp = PolynomialRing(GF(p), 'z'); fp = Rp(cubic_coeff)
        if not fp.is_irreducible(): continue
        degrees = sorted(int(g.degree()) for g, m in Rp(integral).factor() for _ in range(m))
        assert degrees in ([1, 1, 3, 3], [2, 6])
        psi = int(degrees == [2, 6])
        # Independent finite-field radical route: choose A,B in F_(p^3),
        # each with prescribed norm yP,yQ. Then F=Norm(A+B) is in F_p.
        K = GF(p**3, name='theta', modulus=fp); theta = K.gen()
        alpha = (K(xP)-theta).sqrt(); beta = (K(xQ)-theta).sqrt()
        norm_exponent = p*p+p+1
        if alpha**norm_exponent != K(yP): alpha = -alpha
        if beta**norm_exponent != K(yQ): beta = -beta
        assert alpha**norm_exponent == K(yP) and beta**norm_exponent == K(yQ)
        F = (alpha+beta)**norm_exponent
        assert F**p == F and F
        scalar = next(j for j in range(p) if K(j) == F)
        independent_psi = int(pow(scalar, (p-1)//2, p) == p-1)
        assert independent_psi == psi
        table.append({'prime': p, 'octic_factor_degrees': degrees, 'psi': psi,
                      'radical_norm_mod_p': scalar, 'independent_radical_psi': independent_psi})
    assert table
    return {'id': old['id'], 'role': spec['role'], 'case_index': spec['case_index'],
            'basis_indices': spec['basis_indices'], 'input_point_indices': selected,
            'scaled_points': [[str(v) for v in P] for P in scaled],
            'elliptic_scaling_d': str(d), 'cubic_ascending': [str(z) for z in cubic_coeff],
            'cubic_discriminant': str(disc_cubic), 'Galois_group': 'S3',
            'pair_Kummer_signatures': sigs, 'pair_Kummer_dimension': 2,
            'independence_prime_blocks': [[p, list(roots)] for p, roots in blocks],
            'generic_Kummer_dimension': len(generic_sigs), 'pair_dimension_mod_generic': relative,
            'rational_octic_ascending': [str(z) for z in h.list()],
            'integral_octic_root_scale': str(denominator),
            'integral_octic_ascending': [str(z) for z in integral.list()],
            'integral_octic_discriminant': str(discriminant),
            'ramification_statement': 'Every finite ramified prime of the octic splitting field divides this nonzero integral polynomial discriminant; this is a support superset, not an exact field discriminant.',
            'inert_prime_table': table, 'excluded_primes': excluded,
            'psi_counts': dict(sorted(Counter(str(z['psi']) for z in table).items())),
            'software': {'sage': version}, 'bindings': bindings()}


def build(check):
    protocol = r.read(PROTOCOL); records = []
    CHECKPOINT.mkdir(parents=True, exist_ok=True)
    for i in range(len(protocol['cases'])):
        dest = CHECKPOINT/f'case-{i}.json'
        if not dest.exists():
            proc = subprocess.run([sys.executable, str(Path(__file__)), 'worker', '--case', str(i), '--destination', str(dest)],
                                  capture_output=True, text=True, timeout=protocol['worker_timeout_seconds'])
            with (CHECKPOINT/f'case-{i}.log').open('x') as f: f.write(proc.stdout+proc.stderr)
            if proc.returncode or not dest.exists(): raise RuntimeError(f'worker {i}: inspect checkpoint log')
        row = r.read(dest); assert row['bindings'] == bindings()
        if check: assert compute_case(i) == row
        records.append(row)
        print('checkpoint', row['id'], row['psi_counts'], flush=True)
    result = {'schema': 'rank-jump.explicit-governing-octic.v1', 'status': 'PASS',
              'bindings': bindings(), 'finite_certificate': finite_certificate(),
              'symbolic_certificate': symbolic_certificate(), 'rows': records,
              'boundary': 'Explicit retrospective cochain for pairs of supplied rational points; no new twist or original parameter, full Selmer basis, CT radical, rational-solubility theorem, or original-family rank predictor.'}
    if check: assert r.read(OUTPUT) == result
    else: r.write_new(OUTPUT, result)
    print('PASS: exact radical identity, cup correction, faithful octic action, and production Frobenius replay')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('mode', choices=['build', 'check', 'worker'])
    parser.add_argument('--case', type=int); parser.add_argument('--destination', type=Path)
    args = parser.parse_args()
    if args.mode == 'worker': r.write_new(args.destination, compute_case(args.case))
    else: build(args.mode == 'check')
