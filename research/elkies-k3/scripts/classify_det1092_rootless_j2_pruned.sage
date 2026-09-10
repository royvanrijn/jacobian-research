#!/usr/bin/env sage-python
"""Exact orthogonal-norm pruning and resumable per-anchor X1092 census.

For prescribed pairings b against B, the projection has squared norm
b (B G B^t)^(-1) b^t. The residual root components are mutually orthogonal
and orthogonal to B. A component with Dynkin labels l costs l C^(-1) l^t.
Their SUM, not independent full norm budgets, bounds all admissible labels.
Every pruning decision is rational and exhaustive; no heuristic cutoff.
"""
import argparse
from collections import Counter
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path
import runpy
import time

from sage.all import QQ, ZZ, matrix, vector, pari

LEGACY = Path(__file__).with_name('classify_det1092_rootless_j2_niemeier.sage')
old = runpy.run_path(str(LEGACY))
ROOT = old['ROOT']
DEFAULT_WORK = ROOT / 'artifacts/local/elkies-k3/det1092-pruned-v1'


def key(G):
    return tuple(tuple(map(int, r)) for r in G.rows())


@lru_cache(maxsize=24)
def root_data(g):
    G = matrix(ZZ, g)
    roots = old['signed_roots'](G)
    return roots, [r * G for r in roots]


@lru_cache(maxsize=4096)
def residual_data(g, indices):
    G = matrix(ZZ, g)
    roots = [root_data(g)[0][i] for i in indices]
    simple = [old['simple_roots'](G, c) for c in old['root_components'](G, roots)]
    inverses = [(c * G * c.transpose()).inverse() for c in simple]
    return simple, inverses


@lru_cache(maxsize=4096)
def labels(g, bound, strict):
    return tuple(old['labels_up_to'](matrix(QQ, g), bound, strict))


def norm_choices(options, budget, stats):
    """Enumerate all products within the shared budget, in original order."""
    minima = [min((cost for _, cost in opts), default=None) for opts in options]
    if any(x is None for x in minima):
        return
    suffix = [QQ(0)] * (len(options) + 1)
    for i in range(len(options) - 1, -1, -1):
        suffix[i] = suffix[i + 1] + minima[i]
    chosen = []

    def visit(i, spent):
        if spent + suffix[i] > budget:
            stats['pruned_label_branches'] += 1
            return
        if i == len(options):
            yield tuple(chosen), spent
            return
        for item in options[i]:
            if spent + item[1] + suffix[i + 1] > budget:
                stats['pruned_label_branches'] += 1
                continue
            chosen.append(item)
            yield from visit(i + 1, spent + item[1])
            chosen.pop()
    yield from visit(0, QQ(0))


def chamber_vectors(G, B, pairings, norm, strict):
    g = key(G)
    roots, covectors = root_data(g)
    bt = B.transpose()
    indices = tuple(i for i, r in enumerate(covectors) if r * bt == 0)
    simple, inverses = residual_data(g, indices)
    rb = matrix(ZZ, [r for c in simple for r in c.rows()]) if simple else matrix(ZZ, 0, G.nrows())
    dimension = G.nrows() - B.nrows() - rb.nrows()
    base_cost = pairings * (B * G * bt).inverse() * pairings
    budget = QQ(norm) - base_cost
    stats = Counter()
    info = {'residual_root_rank': rb.nrows(), 'fixed_dimension': dimension,
            'base_projection_norm': str(base_cost), 'available_root_norm': str(budget),
            'strict': bool(strict)}
    # In a strict chamber l_i >= 1. ADE inverse Cartan entries are positive,
    # so the all-one labels minimize the component norm exactly.
    minimum_root_cost = sum((sum(C.list()) for C in inverses), QQ(0)) if strict else QQ(0)
    info['minimum_root_norm'] = str(minimum_root_cost)
    if minimum_root_cost > budget:
        return [], dict(info, exact_root_norm_obstruction=True, label_choices=0)
    options = [labels(tuple(tuple(r) for r in C.rows()), budget, strict) for C in inverses]
    info['component_label_counts'] = list(map(len, options))
    info['cartesian_label_products_before_sum_pruning'] = math.prod(map(len, options))
    constraints = (G * bt).augment(G * rb.transpose())
    for coordinate in range(G.nrows()):
        if constraints.ncols() == G.nrows():
            break
        unit = matrix(ZZ, G.nrows(), 1, [int(i == coordinate) for i in range(G.nrows())])
        if constraints.augment(unit).rank() > constraints.rank():
            constraints = constraints.augment(unit)
    inverse = constraints.inverse()
    coordinate_norm = inverse * G * inverse.transpose()
    size = B.nrows() + rb.nrows()
    Q = coordinate_norm[size:, size:]
    qi = Q.inverse() if dimension else matrix(QQ, 0, 0)
    # Verify the orthogonal norm formula against the Schur complement once.
    schur = coordinate_norm[:size, :size]
    cross = coordinate_norm[:size, size:]
    if dimension:
        schur = schur - cross * qi * cross.transpose()
    from sage.all import block_diagonal_matrix
    assert schur == block_diagonal_matrix([(B * G * bt).inverse()] + inverses)
    # Reuse LDL for every translated ellipsoid in this chamber.
    lower_diagonal = old['ldl'](Q) if dimension else None
    candidates = []
    for choice, cost in norm_choices(options, budget, stats):
        stats['label_choices'] += 1
        labs = [entry for component, _ in choice for entry in component]
        prefix = vector(QQ, list(pairings) + labs)
        remaining = budget - cost
        if not dimension and remaining != 0:
            stats['wrong_exact_norm'] += 1
            continue
        centre = -prefix * cross * qi if dimension else vector(QQ, [])
        fixeds = ellipsoid(Q, centre, remaining, lower_diagonal)
        for fixed in fixeds:
            stats['ellipsoid_shell_vectors'] += 1
            candidate = vector(QQ, list(prefix) + list(fixed)) * inverse
            if any(v.denominator() != 1 for v in candidate):
                stats['nonintegral'] += 1
                continue
            candidate = vector(ZZ, candidate)
            assert candidate * G * candidate == norm
            assert candidate * G * bt == pairings
            if not old['primitive'](B.stack(matrix(ZZ, [candidate]))):
                stats['nonprimitive'] += 1
                continue
            candidates.append({'vector': candidate, 'labels': labs, 'fixed_dimension': dimension})
    return candidates, dict(info, **dict(stats))


def ellipsoid(Q, centre, norm, decomposition):
    if Q.nrows() == 0:
        return [tuple()] if norm == 0 else []
    # Use the same exact recursion with the factorization supplied from cache.
    fn = old['ellipsoid_shell']
    namespace = dict(fn.__globals__, ldl=lambda unused: decomposition)
    from types import FunctionType
    return FunctionType(fn.__code__, namespace)(Q, centre, norm)


def bindings():
    paths = [Path(__file__), LEGACY] + [old[n] for n in ('CATALOG', 'ANCHORS', 'AUXILIARY', 'PARENT')]
    return {str(p.relative_to(ROOT)): old['digest'](p) for p in paths}


def write(path, obj):
    raw = old['canonical'](obj)
    temporary = path.with_suffix('.tmp')
    temporary.write_bytes(raw)
    temporary.replace(path)


def run_anchor(number, work):
    work.mkdir(parents=True, exist_ok=True)
    path = work / ('anchor-%02d.json' % number)
    bound = bindings()
    if path.exists():
        saved = json.loads(path.read_text())
        assert saved['bindings'] == bound, 'stale checkpoint; use a fresh work directory'
        if saved['status'] == 'PASS_COMPLETE_ANCHOR':
            print('CACHED', number, flush=True)
            return saved
    else:
        saved = None
    anchors = json.loads(old['ANCHORS'].read_text())['anchors']
    a = anchors[number - 1]
    catalog = json.loads(old['CATALOG'].read_text())['rooted_niemeier_lattices']
    G = matrix(ZZ, next(x for x in catalog if x['label'] == a['niemeier'])['gram'])
    B = matrix(ZZ, a['D5_basis_in_ambient'])
    auxiliary = matrix(ZZ, json.loads(old['AUXILIARY'].read_text())['auxiliary']['gram'])
    assert B * G * B.transpose() == old['D5']
    start = time.monotonic()
    sixths, sixth_info = chamber_vectors(G, B, old['SIXTH_PAIRINGS'], old['SIXTH_NORM'], False)
    print('SIXTHS', number, len(sixths), 'seconds', round(time.monotonic() - start, 3), flush=True)
    encoded_sixths = [{'vector': list(map(int, x['vector'])), 'labels': x['labels']} for x in sixths]
    state = saved or {'schema': 'det1092.pruned-anchor.v1', 'bindings': bound,
        'anchor_number': number, 'anchor': a, 'status': 'RUNNING',
        'sixths': encoded_sixths, 'sixth_accounting': sixth_info,
        'seventh_chambers': [], 'embeddings': []}
    assert state['sixths'] == encoded_sixths
    write(path, state)
    for index in range(len(state['seventh_chambers']), len(sixths)):
        x = sixths[index]
        six = B.stack(matrix(ZZ, [x['vector']]))
        sevenths, info = chamber_vectors(G, six, old['SEVENTH_PAIRINGS'], old['SEVENTH_NORM'], True)
        for y in sevenths:
            basis = six.stack(matrix(ZZ, [y['vector']]))
            assert old['primitive'](basis)
            assert basis * G * basis.transpose() == auxiliary
            complement = (basis * G).right_kernel_matrix()
            H = complement * G * complement.transpose()
            assert H.nrows() == 17 and H.det() == 1092
            assert int(pari(H).qfminim(2)[0]) == 0
            state['embeddings'].append({'sixth_index': index, 'seventh_labels': y['labels'],
                'auxiliary_basis_in_ambient': old['rows'](basis),
                'complement_basis_in_ambient': old['rows'](complement), 'gram': old['rows'](H)})
        state['seventh_chambers'].append(dict(info, rootless_embeddings=len(sevenths)))
        if (index + 1) % 25 == 0 or sevenths:
            write(path, state)
            print('PROGRESS', number, index + 1, '/', len(sixths), 'embeddings', len(state['embeddings']), flush=True)
    state['status'] = 'PASS_COMPLETE_ANCHOR'
    write(path, state)
    print('PASS_ANCHOR', number, a['niemeier'], 'sixths', len(sixths), 'rootless', len(state['embeddings']), 'seconds', round(time.monotonic() - start, 3), flush=True)
    return state


def regression():
    anchors = json.loads(old['ANCHORS'].read_text())['anchors']
    catalog = json.loads(old['CATALOG'].read_text())['rooted_niemeier_lattices']
    def compare(G, B, b, norm, strict):
        start = time.monotonic()
        slow, _ = old['chamber_vectors'](G, B, b, norm, strict)
        slow_seconds = time.monotonic() - start
        start = time.monotonic()
        fast, _ = chamber_vectors(G, B, b, norm, strict)
        assert sorted(tuple(x['vector']) for x in slow) == sorted(tuple(x['vector']) for x in fast)
        print('MATCH', len(fast), 'slow', round(slow_seconds, 3), 'fast', round(time.monotonic() - start, 3), flush=True)
        return fast
    for number in (1, 2):
        a = anchors[number - 1]
        G = matrix(ZZ, next(x for x in catalog if x['label'] == a['niemeier'])['gram'])
        B = matrix(ZZ, a['D5_basis_in_ambient'])
        xs = compare(G, B, old['SIXTH_PAIRINGS'], old['SIXTH_NORM'], False)
        for x in xs:
            compare(G, B.stack(matrix(ZZ, [x['vector']])), old['SEVENTH_PAIRINGS'], old['SEVENTH_NORM'], True)
    # Positive exact shells with a root-free residual space, and strict/non-strict A1.
    for G in (matrix(ZZ, [[2,0,0],[0,4,0],[0,0,6]]), matrix(ZZ, [[2,0,0],[0,2,0],[0,0,4]])):
        for strict in (False, True):
            compare(G, matrix(ZZ, [[1,0,0]]), vector(ZZ, [0]), ZZ(6), strict)
    print('PASS_REGRESSION', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--anchors', type=int, nargs='+')
    p.add_argument('--work', type=Path, default=DEFAULT_WORK)
    p.add_argument('--regression', action='store_true')
    args = p.parse_args()
    if args.regression:
        regression()
    if args.anchors:
        assert all(1 <= n <= 16 for n in args.anchors)
        for n in args.anchors:
            run_anchor(n, args.work)
