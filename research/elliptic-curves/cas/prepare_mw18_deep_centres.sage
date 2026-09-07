#!/usr/bin/env sage -python
"""Bounded, replayable parity census and three equal-exposure MW18 policies."""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys
from time import monotonic

import numpy as np
from fpylll import GSO, IntegerMatrix, Enumeration

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'elliptic-curves/cas'))
from research_runtime.cvp import Hole, VoronoiIterator
from research_runtime.deep_centres import exact_coset_minimum, diverse_deep, parity
from research_runtime.store import checkpoint, digest, atomic_write


def census(geometry, output, count=40, verify=False, proposal_only=False):
    start = monotonic()
    gram = geometry['reduced_gram']; n = len(gram)
    if n != 18 or not 1 <= count <= 100: raise ValueError('rank 18 and at most 100 charts required')
    path = output.with_suffix('.census.json.gz')
    if path.exists():
        retained = json.loads(gzip.decompress(path.read_bytes()))
        if retained['geometry_hash'] != digest(geometry): raise ArithmeticError('census geometry changed')
        reps = np.array(retained['representatives'], dtype=np.int64)
        norms = np.array(retained['norm_upper_bounds'], dtype=np.int64)
    elif verify:
        raise FileNotFoundError('replay needs the retained parity witnesses')
    else:
        gso = GSO.Mat(IntegerMatrix.from_matrix(gram), gram=True, float_type='dd', update=True)
        mu = np.array([[gso.get_mu(i,j) if i > j else float(i == j) for j in range(n)] for i in range(n)])
        reps = np.zeros((1 << n, n), dtype=np.int64); norms = np.zeros(1 << n, dtype=np.int64)
        g = np.array(gram, dtype=np.int64)
        # Each residue / 2 itself supplies a uniform initial search radius.
        radius = sum(abs(v) for row in gram for v in row)/4 + 1
        for mask in range(1 << n):
            residue = np.array([(mask >> i) & 1 for i in range(n)], dtype=np.int64)
            target = tuple(-residue.dot(mu)/2)
            distance, coordinates = Enumeration(gso).enumerate(0, n, radius, 0, target=target)[0]
            rounded = np.rint(coordinates)
            if np.max(np.abs(rounded-coordinates)) > 1e-7: raise ArithmeticError('nonintegral CVP proposal')
            z = residue + 2*rounded.astype(np.int64)
            value = int(z.dot(g).dot(z))
            if abs(4*distance-value) > 1e-6: raise ArithmeticError('proposal norm discrepancy')
            reps[mask] = z; norms[mask] = value
        retained = {'geometry_hash': digest(geometry), 'representatives': reps.tolist(),
                    'norm_upper_bounds': norms.tolist(), 'boundary': 'Upper bounds only until exact maximum-stratum audit.'}
        atomic_write(path, gzip.compress(json.dumps(retained, separators=(',', ':')).encode(), mtime=0))
    # Independent exact integral upper-bound replay; object arithmetic forbids overflow.
    if reps.shape != (1 << n,n) or norms.shape != (1 << n,): raise ArithmeticError('incomplete census')
    if not np.array_equal((reps % 2).dot(1 << np.arange(n)), np.arange(1 << n)):
        raise ArithmeticError('parity census is not a complete partition')
    exact_norms = np.sum(reps.astype(object).dot(np.array(gram, dtype=object))*reps, axis=1)
    if not np.array_equal(exact_norms, norms): raise ArithmeticError('retained norm failed')
    maximum = int(max(norms)); maximum_masks = np.flatnonzero(norms == maximum).tolist()
    if proposal_only:
        print(f"MW18_PROPOSAL|{geometry['cover_label']}|max={maximum}|count={len(maximum_masks)}",flush=True)
        return
    if len(maximum_masks) < count: raise ArithmeticError('too few deepest classes for frozen exposure')
    audits = {}; exact_distances = {}
    def distance(mask):
        if mask not in exact_distances:
            value, rep, nodes = exact_coset_minimum(gram, mask, reps[mask].tolist())
            audits[str(mask)] = {'minimum_norm': value, 'representative': list(rep), 'nodes': nodes}
            exact_distances[mask] = Fraction(value, 4)
        return exact_distances[mask]
    for mask in maximum_masks:
        if distance(mask) != Fraction(maximum,4):
            raise ArithmeticError('floating census overestimated a maximum; repair proposals before selecting')
    deep = [Hole(mask, tuple(map(int,reps[mask])), Fraction(maximum,4)) for mask in maximum_masks]
    diverse = diverse_deep(deep, distance, count)
    if verify:
        old = json.loads(output.read_text())
        # No nearest enumeration in witness replay; validate its exact frontier
        # and every representative below. Ordering discovery remains retained.
        iterator = VoronoiIterator.resume(old['nearest_checkpoint'], binding={'geometry': digest(geometry)})
        near = [Hole.from_record(r) for r in old['policies']['nearest_first']]
    else:
        iterator = VoronoiIterator(gram, binding={'geometry': digest(geometry)})
        near = iterator.next_holes(count, diversity_window=2, node_budget=1000000)
    # The iterator chooses the baseline masks. All arms use the same retained
    # CVP representative table for chart construction, avoiding a tie-solver
    # change between the shallow and deep arms.
    near = [Hole(h.mask, tuple(map(int,reps[h.mask])), h.squared_distance) for h in near]
    policies = {'nearest_first': near, 'deepest': deep[:count], 'diverse_deep': diverse}
    transform = geometry['reduced_basis_rows_in_original_basis']
    original = {}
    for name, holes in policies.items():
        if len(holes) != count: raise ArithmeticError('unequal exposure')
        for h in holes:
            if h.squared_distance != distance(h.mask) or h.squared_distance != Fraction(int(norms[h.mask]),4):
                raise ArithmeticError('selected centre is not shortest')
        original[name] = [dict(coefficients=[sum(h.doubled_coordinates[i]*transform[i][j] for i in range(n))
                                             for j in range(n)]) for h in holes]
    result = {'schema': 'elliptic-curves.mw18-deep-centres.v1', 'geometry_hash': digest(geometry),
        'cover_label': geometry['cover_label'], 'census_path': path.name,
        'census_sha256': sha256(path.read_bytes()).hexdigest(), 'coset_count': 1 << n,
        'upper_bound_histogram': {str(k): v for k,v in sorted(Counter(map(int,norms)).items())},
        'maximum_minimum_norm': maximum, 'half_lattice_depth': str(Fraction(maximum,4)),
        'complete_deepest_stratum': maximum_masks, 'deepest_count': len(deep),
        'exact_cvp_audits': audits, 'nearest_checkpoint': iterator.checkpoint(),
        'policies': {k: [h.record() for h in hs] for k,hs in policies.items()}, 'original_basis_centres': original,
        'selection': {'count': count, 'nearest_first': 'current iterator masks, diversity_window=2, exact reduced generic Gram',
            'representatives': 'common retained CVP table for all three arms; selected minima independently audited exactly',
            'deepest': f'first {count} maximum-depth classes by reduced-basis mask',
            'diverse_deep': f'{count} maximum-depth classes, greedy farthest-first using exact torus distances, then depth and mask'},
        'claim_boundary': 'All 262144 upper bounds and every proposed maximum are checked exactly. This certifies the complete deepest stratum, not all other proposed minima; every distance used in selection is separately exact.'}
    if verify:
        if result != old: raise ArithmeticError('centre-policy replay differs')
    else: checkpoint(output, result)
    print(f"MW18_CENTRES|{geometry['cover_label']}|max={maximum}|deepest={len(deep)}|audits={len(audits)}|seconds={monotonic()-start:.1f}", flush=True)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--geometry',type=Path,required=True); p.add_argument('--directory',type=Path,required=True)
    p.add_argument('--cover',action='append',required=True); p.add_argument('--verify',action='store_true')
    p.add_argument('--count',type=int,default=40); p.add_argument('--proposal-only',action='store_true')
    a=p.parse_args(); d=json.loads(a.geometry.read_text())
    for label in a.cover:
        g=next(r for r in d['covers'] if r['cover_label']==label)
        output=a.directory/(label+'.json')
        census(g, output, count=a.count, verify=a.verify, proposal_only=a.proposal_only)
