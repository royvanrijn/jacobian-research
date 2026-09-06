#!/usr/bin/env python3
"""One point-free cubic norm computation for a retained strict pair."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE/'UNPOINTED_GOVERNING_NORM_PROTOCOL.json'
OUTPUT = r.OUT/'rank_jump_unpointed_governing_norm_v1.json'
WORK = r.ROOT/'artifacts/local/rank-jump-unpointed-governing-norm-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in (Path(__file__), PROTOCOL, Path(r.__file__))}


def worker():
    from sage.all import QQ, PolynomialRing, pari, NumberField, lcm, ZZ, GF
    from sage.version import version
    spec = r.read(PROTOCOL)
    pari.allocatemem(64000000, spec['limits']['pari_stack_bytes'], silent=True)
    pari.setrand(spec['seed'])
    R = PolynomialRing(QQ, 'y'); f = R(spec['cubic_ascending'])
    nf = pari.nfinit(f)
    alpha = pari.Mod(pari(R(spec['alpha_ascending'])), pari(f))
    beta = pari.Mod(pari(R(spec['beta_ascending'])), pari(f))
    relative = pari('x')**2-alpha
    setup = pari.rnfisnorminit(nf, relative, 1)
    z, remainder = pari.rnfisnorm(setup, beta)
    assert remainder == 1
    X, Y = [z.lift().polcoef(i) for i in range(2)]
    assert Y and X*X-alpha*Y*Y == beta
    enc = lambda v: [str(pari.lift(v).polcoef(i)) for i in range(3)]
    K = NumberField(f, 'theta')
    a, b, x, y = [K(list(map(QQ, row))) for row in (spec['alpha_ascending'], spec['beta_ascending'], enc(X), enc(Y))]
    assert x*x-a*y*y == b
    norm_a, norm_b = a.norm(), b.norm()
    assert norm_a.is_square() and norm_b.is_square()
    p = x.norm(); b0 = norm_b.sqrt(); delta = x*x-b; D = delta.norm()
    assert D and D == norm_a*y.norm()**2
    S = PolynomialRing(QQ, 'T'); T = S.gen()
    middle = 2*D*((x*x+b)/delta).trace()
    h = T**8-4*(p+b0)*T**6+middle*T**4-4*D*(p-b0)*T**2+D**2
    d = lcm(c.denominator() for c in h.list())
    integral = S([h[i]*d**(8-i) for i in range(9)])
    assert all(v.denominator() == 1 for v in integral.list())
    disc = ZZ(integral.discriminant()); assert disc
    # Bounded worker includes this small exact Galois identification.
    gal = pari.polgalois(integral)
    table = []; excluded = []
    denominator_support = lcm(v.denominator() for obj in (a, b, x, y) for v in obj.list())
    for ell in r.primes(spec['limits']['prime_bound']):
        if any(v % ell == 0 for v in (disc, d, denominator_support, norm_a, norm_b, D, f.discriminant())):
            excluded.append(ell); continue
        Fp = GF(ell); Rp = PolynomialRing(Fp, 's'); fp = Rp(f.list())
        if not fp.is_irreducible(): continue
        F = GF(ell**3, name='u', modulus=fp); u = F.gen()
        eval_at = lambda v: sum(F(c)*u**i for i, c in enumerate(v.list()))
        bv, xv = eval_at(b), eval_at(x)
        radical = bv.sqrt(); exponent = ell*ell+ell+1
        if radical**exponent != F(b0): radical = -radical
        assert radical**exponent == F(b0)
        norm = (xv+radical)**exponent
        assert norm and norm**ell == norm
        scalar = next(i for i in range(ell) if F(i) == norm)
        psi = int(pow(scalar, (ell-1)//2, ell) == ell-1)
        degrees = sorted(int(g.degree()) for g, m in Rp(integral).factor() for _ in range(m))
        assert degrees == ([1, 1, 3, 3] if psi == 0 else [2, 6])
        table.append({'prime': ell, 'factor_degrees': degrees, 'norm_mod_prime': scalar, 'psi': psi})
    assert table
    return {'schema': 'rank-jump.unpointed-governing-norm.v1', 'status': 'PASS', 'bindings': bindings(),
            'cubic_ascending': spec['cubic_ascending'], 'alpha_ascending': spec['alpha_ascending'],
            'beta_ascending': spec['beta_ascending'], 'X_ascending': enc(X), 'Y_ascending': enc(Y),
            'norm_alpha': str(norm_a), 'norm_beta': str(norm_b), 'norm_X': str(p),
            'norm_X_squared_minus_beta': str(D), 'beta_norm_square_root': str(b0),
            'rational_octic_ascending': [str(v) for v in h.list()], 'integral_root_scale': str(d),
            'integral_octic_ascending': [str(v) for v in integral.list()],
            'integral_octic_discriminant': str(disc), 'pari_polgalois': str(gal),
            'galois_order': int(gal[0]), 'inert_prime_table': table, 'excluded_primes': excluded,
            'software': {'sage': version, 'pari': str(pari.version())},
            'boundary': 'No elliptic points read or searched. The octic encodes a cochain, not rational solubility; the selected classes are already known strict Sha controls.'}


def capture():
    WORK.mkdir(parents=True, exist_ok=True); path = WORK/'witness.json'
    if not path.exists():
        with (WORK/'worker.log').open('x') as log:
            try:
                process = subprocess.run([sys.executable, str(Path(__file__)), 'worker'], stdout=log, stderr=log,
                                         timeout=r.read(PROTOCOL)['limits']['worker_seconds'])
                failure = None if process.returncode == 0 else 'worker failed'
            except subprocess.TimeoutExpired: failure = 'bounded worker timed out'
        if failure: raise RuntimeError(f'{failure}; inspect {WORK}/worker.log')
    row = r.read(path); assert row['bindings'] == bindings()
    r.write_new(OUTPUT, row)
    print('PASS norm cochain:', row['X_ascending'], row['Y_ascending'], row['pari_polgalois'], len(row['inert_prime_table']))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('mode', choices=['worker', 'capture']); args = parser.parse_args()
    if args.mode == 'worker': r.write_new(WORK/'witness.json', worker())
    else: capture()
