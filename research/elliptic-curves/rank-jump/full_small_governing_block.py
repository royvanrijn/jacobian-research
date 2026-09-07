#!/usr/bin/env python3
"""Two bounded norm equations completing the retained small Selmer basis."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE/'FULL_SMALL_GOVERNING_BLOCK_PROTOCOL.json'
PRIOR = r.OUT/'rank_jump_unpointed_governing_norm_v1.json'
OUTPUT = r.OUT/'rank_jump_full_small_governing_block_inputs_v1.json'
WORK = r.ROOT/'artifacts/local/rank-jump-full-small-governing-block-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in
            (Path(__file__), PROTOCOL, PRIOR, Path(r.__file__), HERE/'NORM_LIFTS_CAN_BE_ENTIRELY_SHA.md')}


def worker(index):
    from sage.all import QQ, ZZ, PolynomialRing, pari, NumberField
    from sage.version import version
    spec = r.read(PROTOCOL); indices = spec['norm_pairs'][index]
    pari.allocatemem(64000000, spec['limits']['pari_stack_bytes'], silent=True); pari.setrand(spec['seed'])
    R = PolynomialRing(QQ, 'y'); f = R(spec['cubic_ascending']); nf = pari.nfinit(f)
    left, right = [spec['classes_ascending'][i] for i in indices]
    alpha = pari.Mod(pari(R(left)), pari(f)); beta = pari.Mod(pari(R(right)), pari(f))
    setup = pari.rnfisnorminit(nf, pari('x')**2-alpha, 1)
    z, remainder = pari.rnfisnorm(setup, beta); assert remainder == 1
    X, Y = [z.lift().polcoef(i) for i in range(2)]
    assert Y and X*X-alpha*Y*Y == beta
    enc = lambda v: [str(pari.lift(v).polcoef(i)) for i in range(3)]
    K = NumberField(f, 'theta'); a, b, x, y = [K(list(map(QQ, v))) for v in (left, right, enc(X), enc(Y))]
    assert x*x-a*y*y == b and a.norm().is_square() and b.norm().is_square()
    p = x.norm(); b0 = b.norm().sqrt(); delta = x*x-b; D = delta.norm(); assert D
    S = PolynomialRing(QQ, 'T'); T = S.gen()
    h = T**8-4*(p+b0)*T**6+2*D*((x*x+b)/delta).trace()*T**4-4*D*(p-b0)*T**2+D**2
    reduced, change = pari.polredabs(h, 1)
    g = S(reduced); mapping = S(pari.lift(change)); assert h(mapping) % g == 0
    disc = ZZ(g.discriminant()); assert disc
    group = pari.polgalois(g)
    return {'schema': 'rank-jump.full-small-governing-pair.v1', 'bindings': bindings(),
            'class_indices': indices, 'alpha_ascending': left, 'beta_ascending': right,
            'X_ascending': enc(X), 'Y_ascending': enc(Y), 'rational_octic_ascending': [str(v) for v in h.list()],
            'reduced_octic_ascending': [str(v) for v in g.list()],
            'old_root_in_reduced_algebra': [str(v) for v in mapping.list()],
            'reduced_octic_discriminant': str(disc), 'discriminant_factorization': [[str(p), int(e)] for p, e in disc.factor()],
            'pari_polgalois': str(group), 'galois_order': int(group[0]),
            'software': {'sage': version, 'pari': str(pari.version())}}


def capture():
    WORK.mkdir(parents=True, exist_ok=True); rows = []
    for i in range(2):
        path = WORK/f'pair-{i}.json'
        if not path.exists():
            with (WORK/f'pair-{i}.log').open('x') as log:
                try:
                    p = subprocess.run([sys.executable, str(Path(__file__)), 'worker', '--index', str(i)],
                                       stdout=log, stderr=log, timeout=r.read(PROTOCOL)['limits']['seconds_per_worker'])
                    failure = None if p.returncode == 0 else 'worker failed'
                except subprocess.TimeoutExpired: failure = 'worker timed out'
            if failure: raise RuntimeError(f'{failure}: inspect {WORK}/pair-{i}.log')
        row = r.read(path); assert row['bindings'] == bindings(); rows.append(row)
        print('checkpoint', row['class_indices'], row['X_ascending'], row['Y_ascending'], row['pari_polgalois'], flush=True)
    r.write_new(OUTPUT, {'schema': 'rank-jump.full-small-governing-block-inputs.v1', 'status': 'PASS',
                        'bindings': bindings(), 'rows': rows,
                        'boundary': 'Two norm cochains for the fixed full Selmer basis; no new elliptic points or descents.'})


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('mode', choices=['worker', 'capture']); p.add_argument('--index', type=int)
    args = p.parse_args()
    if args.mode == 'worker': r.write_new(WORK/f'pair-{args.index}.json', worker(args.index))
    else: capture()
