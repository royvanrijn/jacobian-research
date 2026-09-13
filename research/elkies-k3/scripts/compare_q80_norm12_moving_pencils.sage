#!/usr/bin/env sage-python
"""Compare full moving-zero-contact norm12 pencils with the norm8 bank.

All inputs are generic. The retained product-character targets are not read:
only the parity inventory's 49 trace words are projected when freezing.
Each completed prime stage is immutable and separately reproducible.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import resource
import runpy
import time

import numpy as np
from sage.all import GF, PolynomialRing, matrix

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elkies-k3-q80-norm12-moving-pencils-v1'
BANK = ROOT / 'artifacts/generated-results/elkies-k3-q80-complete-genus-one-pairs-v1'
PARITY = ROOT / 'artifacts/generated-results/elkies-k3-r17-norm12-11952-product-tate-parity-v1.json'
WORKER = ROOT / 'elkies-k3/scripts/compare_q80_complete_genus_one_pencils.sage'
W = runpy.run_path(str(WORKER))
read, write_new = W['read'], W['write_new']


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def freeze(path):
    old = read(BANK / 'input.json.gz')
    parity = read(PARITY)
    records = parity['invariant_trace_parity']['deep_norm12_classes']
    words = [r['section_basis_w'] for r in records]
    gram = np.asarray(old['gram'], dtype=np.int64)
    ww = np.asarray(words, dtype=np.int64)
    assert len(words) == 49 and np.all(np.sum((ww @ gram) * ww, axis=1) == 12)
    assert len({tuple(c % 2 for c in w) for w in words}) == 49
    assert not ({tuple(c % 2 for c in w) for w in words} &
                {tuple(c % 2 for c in w) for w in old['words']})
    packet = {k: old[k] for k in ['A', 'B', 'basis', 'gram']}
    packet.update({
        'schema': 'q80-norm12-moving-pencils-input-v1',
        'words': words, 'orbit_masks': [r['orbit_mask'] for r in records],
        'norm8_class_count': len(old['words']),
        'scope': 'All49 retained minimum-norm12 pencils, full moving intersection with O; compare all1176 internal pairs and3131933 pairs with the complete63917 norm8 pencils, over the fixed original parameter.',
        'formula': 'For monic deg4 h and deg7 M0, c=coeff(M0,t^7), g=v*t-u, M=g*M0-v*c*h^2. Q=(M^4-6*M^2*Nx*g^2-8*M*Ny*g^3-3*Nx^2*g^4-4*A*h^4*g^4)/h^6; Q has bidegree at most(4,4).',
        'primes': [521, 523, 541],
        'limits': {'cpu_seconds_per_prime': 60, 'address_space_gib': 4, 'checkpoint_size': 512},
        'quarantine': 'Generic model, generic basis and complete parity words only. No product targets, exceptional points or specialization ranks.',
        'acceptance': old['acceptance'],
        'source_hashes': {str(p.relative_to(ROOT)): digest(p) for p in
                          [PARITY, WORKER, BANK / 'input.json.gz', BANK / 'frame-manifest.json',
                           BANK / 'independent-replay-final.json']},
        'goal_complete': False,
    })
    write_new(path / 'input.json.gz', packet)
    print(json.dumps({'status': 'FROZEN', 'norm12_pencils': 49,
                      'new_pairs': 3133109, 'primes': packet['primes']}), flush=True)


def compile_norm12(ctx, word, helper, inverse):
    R, K, A, B, E, basis, tracer = ctx
    t = R.gen()
    P = tracer(word)
    assert W['degree_height'](P) == 12, 'trace height changed'
    X, Y = K(P[0]), K(P[1])
    for shift in (None, 0, 1, 2, 3, 4):
        if shift is None:
            Xc, Yc, Ac = X, Y, A
        else:
            Xc = helper['invert_rational'](X(t + shift), 4, R, K)
            Yc = helper['invert_rational'](Y(t + shift), 6, R, K)
            Ac = inverse['reciprocal_polynomial'](A(t + shift), 8, R)
        frame = helper['trace_chord_frame'](Xc, Yc, R)
        h, nx, ny, m = (frame[k] for k in ['h', 'Nx', 'Ny', 'M0'])
        if h.degree() == 4:
            break
    else:
        raise ArithmeticError('no complete finite-pole chart')
    assert h.is_monic() and m.degree() == 7, 'norm12 regular frame changed'
    U = PolynomialRing(R, 'u')
    u = U.gen()
    g = U(t) - u
    moving = g * m - m[7] * h**2
    numerator = (moving**4 - 6*moving**2*nx*g**2 - 8*moving*ny*g**3
                 - 3*nx**2*g**4 - 4*Ac*h**4*g**4)
    qs = []
    for i in range(5):
        q, rem = R(numerator[i]).quo_rem(h**6)
        assert not rem and q.degree() <= 4, 'moving branch formula failed'
        qs.append(q)
    # The diagonal gives the two points at the moving contact with O and tau.
    diagonal = sum(q*t**i for i, q in enumerate(qs))
    assert diagonal == m[7]**4*h**2, 'moving rational-point identity failed'
    if shift is not None:
        qs = [R(sum(q[i]*(t-shift)**(4-i) for i in range(5))) for q in qs]
    branch = [[int(q[i]) for q in qs] for i in range(5)]
    return {'chart_shift': shift, 'h': list(map(int, h)), 'Nx': list(map(int, nx)),
            'Ny': list(map(int, ny)), 'M0': list(map(int, m)), 'branch_matrix': branch,
            'full_rank': bool(matrix(GF(R.characteristic()), branch).det())}


def load_old_frames(prime):
    result = read(BANK / ('prime-%d.json.gz' % prime))
    frames = {}
    for filename in result['checkpoint_files']:
        for r in read(BANK / ('prime-%d' % prime) / filename)['records']:
            frames[r['index']] = r
    assert sorted(frames) == result['compiled_indices']
    return frames


def cross_matches(old_records, deep_records, p):
    """Intersect packed projective images, retaining every matching member."""
    deep = W['packed_images']([r['branch_matrix'] for r in deep_records],
                              [r['index'] for r in deep_records], p)
    order = np.argsort(deep // np.uint64(65536))
    deep = deep[order]
    keys = deep // np.uint64(65536)
    old = W['packed_images']([r['branch_matrix'] for r in old_records],
                             [r['index'] for r in old_records], p)
    queries = old // np.uint64(65536)
    left, right = np.searchsorted(keys, queries, 'left'), np.searchsorted(keys, queries, 'right')
    matches = set()
    for k in np.flatnonzero(left < right):
        i = int(old[k] % np.uint64(65536))
        for value in deep[left[k]:right[k]]:
            matches.add((i, int(value % np.uint64(65536))))
    return matches


def run_prime(path, p):
    packet = read(path / 'input.json.gz')
    assert p in packet['primes']
    limit = packet['limits']['cpu_seconds_per_prime']
    resource.setrlimit(resource.RLIMIT_CPU, (limit, limit+5))
    resource.setrlimit(resource.RLIMIT_AS, (4*1024**3, 4*1024**3))
    start = time.process_time()
    earlier = [q for q in packet['primes'][:packet['primes'].index(p)]]
    for q in earlier:
        assert (path / ('prime-%d.json.gz' % q)).exists(), 'missing earlier stage'
    prior = read(path / ('prime-%d.json.gz' % earlier[-1])) if earlier else None
    old_indices = sorted({a for a, b in prior['surviving_cross_pairs']}) if prior else list(range(packet['norm8_class_count']))
    helper, inverse = runpy.run_path(str(W['HELPER'])), runpy.run_path(str(W['INVERSION']))
    ctx = W['context'](packet, p)
    checkpoint = path / ('frames-%d.json.gz' % p)
    if checkpoint.exists():
        deep = read(checkpoint)['records']
    else:
        deep = []
        for i, word in enumerate(packet['words']):
            row = compile_norm12(ctx, word, helper, inverse)
            row['index'] = i
            deep.append(row)
        write_new(checkpoint, {'input_sha256': digest(path / 'input.json.gz'), 'prime': p, 'records': deep})
    assert [r['index'] for r in deep] == list(range(49))
    print(json.dumps({'prime': p, 'norm12_compiled': 49, 'norm8_selected': len(old_indices),
                      'cpu_seconds': time.process_time()-start}), flush=True)
    old_frames = load_old_frames(p)
    missing = [i for i in old_indices if i not in old_frames]
    supplemental = []
    if missing:
        old_packet = read(BANK / 'input.json.gz')
        # The full old word set also pins the cached scalar-multiple bound.
        old_ctx = W['context'](old_packet, p)
        for i in missing:
            row = W['compile_frame'](old_ctx, old_packet['words'][i], helper, inverse)
            row['index'] = i
            old_frames[i] = row
            supplemental.append(row)
    found = set()
    size = packet['limits']['checkpoint_size']
    for offset in range(0, len(old_indices), size):
        chosen = old_indices[offset:offset+size]
        found.update(cross_matches([old_frames[i] for i in chosen], deep, p))
    if prior:
        found &= set(map(tuple, prior['surviving_cross_pairs']))
    packed = W['packed_images']([r['branch_matrix'] for r in deep], list(range(49)), p)
    buckets = W['collision_buckets'](packed)
    internal = {(a, b) for bucket in buckets for a in bucket['members'] for b in bucket['members'] if a < b}
    if prior:
        internal &= set(map(tuple, prior['surviving_internal_pairs']))
    result = {
        'schema': 'q80-norm12-moving-pencils-prime-v1', 'prime': p,
        'input_sha256': digest(path / 'input.json.gz'),
        'frames_sha256': digest(checkpoint), 'compiled_norm8_indices': old_indices,
        'supplemental_norm8_frames': supplemental,
        'surviving_cross_pairs': sorted(found), 'surviving_internal_pairs': sorted(internal),
        'norm12_internal_buckets': buckets,
        'norm12_full_rank_indices': [r['index'] for r in deep if r['full_rank']],
        'projective_image_count': (len(old_indices)+49)*(p+1),
        'status': 'ALL_NEW_PAIRS_EXCLUDED' if not found and not internal else 'FINITE_CANDIDATES_REQUIRE_LATER_PRIME_OR_EXACT_RESOLUTION',
        'cpu_seconds': time.process_time()-start, 'goal_complete': False,
    }
    write_new(path / ('prime-%d.json.gz' % p), result)
    print(json.dumps({k: result[k] for k in ['prime', 'status', 'cpu_seconds', 'projective_image_count']} |
                     {'cross_pairs': len(found), 'internal_pairs': len(internal),
                      'full_rank': len(result['norm12_full_rank_indices'])}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['freeze', 'prime'])
    parser.add_argument('--output', type=Path, default=DEFAULT)
    parser.add_argument('--prime', type=int)
    args = parser.parse_args()
    if args.mode == 'freeze':
        freeze(args.output)
    else:
        run_prime(args.output, args.prime)
