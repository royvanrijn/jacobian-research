#!/usr/bin/env python3
"""Independent arithmetic replay of the full norm12 moving-pencil extension.

The previous norm8 theorem certifies its retained trace frames. Their hashes
are checked here and their projective images are evaluated again. The 98 new
norm12 frames and the one supplemental norm8 frame are independently checked.
"""
import argparse
from collections import Counter
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT / 'artifacts/generated-results/elkies-k3-q80-norm12-moving-pencils-v1'
BANK = ROOT / 'artifacts/generated-results/elkies-k3-q80-complete-genus-one-pairs-v1'
PARITY = ROOT / 'artifacts/generated-results/elkies-k3-r17-norm12-11952-product-tate-parity-v1.json'
SPEC = importlib.util.spec_from_file_location('norm8_replay', ROOT / 'elkies-k3/scripts/verify_q80_complete_genus_one_pairs.py')
V = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(V)
read, require, ev, F = V.read, V.require, V.ev, V.F


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def require_completed_stages(path):
    for p in (521, 523):
        for prefix in ('prime', 'frames'):
            require((path / ('%s-%d.json.gz' % (prefix, p))).is_file(), 'missing completed stage %s-%d' % (prefix, p))


def add(P, Q, a, p):
    if P is None:
        return Q
    if Q is None:
        return P
    x, y = P
    u, v = Q
    if x == u and (y+v) % p == 0:
        return None
    m = ((3*x*x+a)*pow(2*y, -1, p) if x == u else (v-y)*pow((u-x) % p, -1, p)) % p
    z = (m*m-x-u) % p
    return z, (m*(x-z)-y) % p


def trace_context(packet, p):
    # This independently checks the basis equations and generic height Gram.
    ctx = V.context(packet, p)
    basis = [(V.primitive_rf(r['X'], p), V.primitive_rf(r['Y'], p)) for r in packet['basis']]
    sites = [t for t in range(p) if (4*ev(ctx['A'], t, p)**3+27*ev(ctx['B'], t, p)**2) % p][:27]
    require(len(sites) == 27, '27 smooth trace-identification fibres required')
    maximum = max(abs(c) for w in packet['words'] for c in w)
    sums = []
    for t in sites:
        a = ev(ctx['A'], t, p)
        multiples = []
        for x, y in basis:
            xd, yd = ev(x[1], t, p), ev(y[1], t, p)
            if xd:
                require(yd, 'inconsistent affine basis specialization')
                point = ev(x[0], t, p)*pow(xd, -1, p) % p, ev(y[0], t, p)*pow(yd, -1, p) % p
            else:
                require(not yd, 'inconsistent zero-section specialization')
                point = None
            row = {0: None}
            for k in range(1, maximum+1):
                row[k] = add(row[k-1], point, a, p)
                row[-k] = None if row[k] is None else (row[k][0], -row[k][1] % p)
            multiples.append(row)
        word_sums = []
        for word in packet['words']:
            point = None
            for i, n in enumerate(word):
                if n:
                    point = add(point, multiples[i][n], a, p)
            word_sums.append(point)
        sums.append(word_sums)
    return ctx | {'identification_sites': sites, 'identified_word_values': sums}


def validate_norm12(records, ctx, identify=True):
    p = ctx['p']
    for r in records:
        # Every emitted frame used this chart. No unchecked chart fallback.
        require(r['chart_shift'] is None, 'unimplemented norm12 replay chart')
        require(len(r['h']) == 5 and r['h'][-1] == 1 and len(r['Nx']) <= 13 and
                len(r['Ny']) <= 19 and len(r['M0']) == 8 and r['M0'][-1], 'norm12 degree bounds')
        require(len(r['branch_matrix']) == 5 and all(len(row) == 5 for row in r['branch_matrix']), 'branch matrix shape')
        values = [c for k in ['h', 'Nx', 'Ny', 'M0'] for c in r[k]] + [c for row in r['branch_matrix'] for c in row]
        require(all(isinstance(c, int) and 0 <= c < p for c in values), 'canonical finite coefficients')
        require(len(F.ff_gcd(r['h'], r['Nx'], p)) == 1, 'cancelled norm12 trace poles')
        h2 = F.ff_mul(r['h'], r['h'], p)
        congruence = F.ff_sub(F.ff_mul(r['M0'], r['Nx'], p), [-x % p for x in r['Ny']], p)
        require(not F.ff_rem(congruence, h2, p), 'norm12 chord congruence')
        # Degree <=36 in t: 37 distinct values certify the exact identity.
        for t in range(37):
            h, nx, ny = (ev(r[k], t, p) for k in ['h', 'Nx', 'Ny'])
            a, b = ev(ctx['A'], t, p), ev(ctx['B'], t, p)
            require((ny*ny-nx**3-a*nx*h**4-b*h**6) % p == 0, 'degree36 trace Weierstrass identity')
        # Bidegree <=(28,4): the complete 29 by 5 interpolation grid is exact.
        for t in range(29):
            h, nx, ny, m0 = (ev(r[k], t, p) for k in ['h', 'Nx', 'Ny', 'M0'])
            a, c = ev(ctx['A'], t, p), r['M0'][7]
            q_in_u = [ev([row[j] for row in r['branch_matrix']], t, p) for j in range(5)]
            for u in range(5):
                g = (t-u) % p
                m = (g*m0-c*h*h) % p
                rhs = m**4-6*m*m*nx*g*g-8*m*ny*g**3-3*nx*nx*g**4-4*a*h**4*g**4
                require((ev(q_in_u, u, p)*h**6-rhs) % p == 0, 'degree28 moving branch identity')
            # Degree <=8 after restricting the bidegree(4,4) curve to t=u.
            if t < 9:
                require((ev(q_in_u, t, p)-c**4*h*h) % p == 0, 'moving rational-point identity')
        if identify:
            matched = 0
            for t, values in zip(ctx['identification_sites'], ctx['identified_word_values']):
                h = ev(r['h'], t, p)
                if not h:
                    continue
                inv = pow(h, -1, p)
                point = ev(r['Nx'], t, p)*inv*inv % p, ev(r['Ny'], t, p)*inv*inv*inv % p
                require(point == values[r['index']], 'height-bounded norm12 trace identification')
                matched += 1
            require(matched >= 23, '23 affine trace comparisons required')
    ranks = V.full_rank([r['branch_matrix'] for r in records], p)
    require([bool(x) for x in ranks] == [r['full_rank'] for r in records], 'norm12 embedding witnesses')
    return ranks


def cross_pairs(old_records, deep_records, p):
    deep = V.images(deep_records, p)
    membership = {}
    for value in deep:
        key, i = divmod(int(value), 65536)
        membership.setdefault(key, set()).add(i)
    keys = np.asarray(list(membership), dtype=np.uint64)
    old = V.images(old_records, p)
    relevant = old[np.isin(old // np.uint64(65536), keys)]
    found = set()
    for value in relevant:
        key, i = divmod(int(value), 65536)
        found.update((i, j) for j in membership[key])
    return found


def bank_frames(p, manifest):
    result = read(BANK / ('prime-%d.json.gz' % p))
    out = {}
    for name in result['checkpoint_files']:
        relative = 'prime-%d/%s' % (p, name)
        require(digest(BANK / relative) == manifest['files'][relative], 'inherited norm8 frame hash')
        for r in read(BANK / relative)['records']:
            require(r['index'] not in out, 'duplicate inherited norm8 frame')
            out[r['index']] = r
    require(sorted(out) == result['compiled_indices'], 'inherited norm8 checkpoint coverage')
    return out


def verify(path, progress=False):
    require_completed_stages(path)
    packet = read(path / 'input.json.gz')
    require(packet['primes'] == [521, 523, 541], 'frozen prime panel')
    for name, expected in packet['source_hashes'].items():
        require(digest(ROOT / name) == expected, 'frozen generic source hash: ' + name)
    old = read(BANK / 'input.json.gz')
    direct = read(V.DIRECT)
    require(digest(V.DIRECT) == old['source_hashes'][str(V.DIRECT.relative_to(ROOT))], 'inherited direct model hash')
    require(packet['A'] == direct['weierstrass_model']['A_coefficients_low_to_high'] and
            packet['B'] == direct['weierstrass_model']['B_coefficients_low_to_high'], 'generic model projection')
    require(packet['basis'] == [{'X': r['X'], 'Y': r['Y']} for r in direct['sections']['records']], 'generic basis projection')
    require(all(packet[k] == old[k] for k in ['A', 'B', 'basis', 'gram']), 'inherited generic context')
    parity = read(PARITY)['invariant_trace_parity']
    records = parity['deep_norm12_classes']
    require(packet['words'] == [r['section_basis_w'] for r in records] and
            packet['orbit_masks'] == [r['orbit_mask'] for r in records], 'complete deep parity projection')
    require(parity['isotropic_minimum_norm_partition'] == {'4': 1313, '8': 63917, '12': 49}, 'inherited complete isotropic inventory')
    words = packet['words']
    require(len(words) == 49 and len({tuple(c % 2 for c in w) for w in words}) == 49, 'complete distinct deep class coverage')
    ww, gg = np.asarray(words, dtype=np.int64), np.asarray(packet['gram'], dtype=np.int64)
    require(np.all(np.sum((ww @ gg)*ww, axis=1) == 12), 'generic norm12 trace heights')
    require(not ({tuple(c % 2 for c in w) for w in words} & {tuple(c % 2 for c in w) for w in old['words']}), 'deep and norm8 parity overlap')
    require(packet['norm8_class_count'] == len(old['words']) == 63917, 'complete norm8 comparison size')
    inherited = read(BANK / 'independent-replay-final.json')
    require(inherited['status'] == 'PASS_COMPLETE_Q80_GENUS_ONE_PENCIL_INJECTIVITY' and inherited['pair_count'] == 2042659486, 'inherited norm8 theorem receipt')
    manifest = read(BANK / 'frame-manifest.json')
    prior_cross, prior_internal, embeddings = None, None, set()
    total_images, frames, supplemental_count, stages = 0, 0, 0, []
    for p in (521, 523):
        result = read(path / ('prime-%d.json.gz' % p))
        checkpoint = read(path / ('frames-%d.json.gz' % p))
        frozen = digest(path / 'input.json.gz')
        require(result['input_sha256'] == checkpoint['input_sha256'] == frozen and result['prime'] == checkpoint['prime'] == p, 'new prime input binding')
        require(digest(path / ('frames-%d.json.gz' % p)) == result['frames_sha256'], 'new norm12 frame binding')
        deep = checkpoint['records']
        require([r['index'] for r in deep] == list(range(49)), 'missing or repeated norm12 frame')
        ctx = trace_context(packet, p)
        ranks = validate_norm12(deep, ctx)
        full = [r['index'] for r, yes in zip(deep, ranks) if yes]
        require(full == result['norm12_full_rank_indices'], 'per-prime embeddings')
        embeddings.update(full)
        frames += len(deep)
        old_indices = sorted({i for i, j in prior_cross}) if prior_cross is not None else list(range(63917))
        require(result['compiled_norm8_indices'] == old_indices, 'adaptive cross-pair coverage')
        old_frames = bank_frames(p, manifest)
        missing = [i for i in old_indices if i not in old_frames]
        extra = result['supplemental_norm8_frames']
        require([r['index'] for r in extra] == missing, 'missing supplemental norm8 frames')
        if extra:
            V.validate_records(extra, V.context(old, p))
            supplemental_count += len(extra)
            old_frames.update({r['index']: r for r in extra})
        found = set()
        for offset in range(0, len(old_indices), 512):
            chosen = old_indices[offset:offset+512]
            found.update(cross_pairs([old_frames[i] for i in chosen], deep, p))
        if prior_cross is not None:
            found &= prior_cross
        internal_buckets = V.buckets(V.images(deep, p))
        require(Counter(internal_buckets) == Counter(tuple(b['members']) for b in result['norm12_internal_buckets']), 'internal collision bucket replay')
        internal = {(a, b) for members in internal_buckets for a in members for b in members if a < b}
        if prior_internal is not None:
            internal &= prior_internal
        require(sorted(map(list, found)) == result['surviving_cross_pairs'] and sorted(map(list, internal)) == result['surviving_internal_pairs'], 'independent projective collision replay')
        count = (len(old_indices)+49)*(p+1)
        require(count == result['projective_image_count'], 'complete projective parameter coverage')
        total_images += count
        stages.append({'prime': p, 'norm12_pencils': 49, 'norm8_pencils': len(old_indices),
                       'projective_parameters': count, 'cross_pairs_remaining': len(found),
                       'internal_pairs_remaining': len(internal), 'full_rank_norm12': len(full)})
        if progress:
            print(json.dumps(stages[-1]), flush=True)
        prior_cross, prior_internal = found, internal
    require(not prior_cross and not prior_internal, 'unresolved shared branch candidates')
    require(embeddings == set(range(49)), 'incomplete norm12 within-pencil injectivity')
    return {'status': 'PASS_COMPLETE_Q80_SMOOTH_GENUS_ONE_BISECTION_INJECTIVITY',
            'norm12_pencils': 49, 'norm8_pencils': 63917, 'total_pencils': 63966,
            'new_pairs_excluded': 3133109, 'total_pairs_excluded_with_inherited_theorem': 2045792595,
            'new_trace_frames_verified': frames, 'supplemental_norm8_frames_verified': supplemental_count,
            'projective_parameters_replayed': total_images, 'prime_stages': stages,
            'scope': 'The full norm12 pencils have disjoint rational projective branch images from each other and the full norm8 bank; all individual pencils are injective. The written nef/translation classification exhausts smooth genus-one bisections modulo inherited section translation. Singular higher-arithmetic-genus images remain outside this theorem.',
            'inherited_dependencies': ['direct11952 rational MW17 model and saturated height lattice',
                                       'complete M/2M minimum-norm inventory',
                                       'prior complete norm8 branch-injectivity theorem and hashed trace frames'],
            'goal_complete': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=DEFAULT)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--progress', action='store_true')
    args = parser.parse_args()
    start = time.monotonic()
    result = verify(args.input_dir, args.progress)
    result['elapsed_seconds'] = round(time.monotonic()-start, 3)
    if args.output:
        with args.output.open('x') as out:
            json.dump(result, out, indent=2, sort_keys=True)
            out.write('\n')
    print(json.dumps(result, sort_keys=True))
