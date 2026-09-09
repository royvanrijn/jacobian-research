"""Bounded exact generic shell bank for future dimension-aware V3 runs."""
from math import isqrt
from visibility_lattice_fast import IntegerExactParity


def enumerate_bank(reduced_gram, transform, *, bound, shells, node_limit):
    oracle = IntegerExactParity(reduced_gram)
    n = oracle.n
    if len(transform) != n or any(len(row) != n for row in transform):
        raise ValueError('wrong transformation dimension')
    if bound < 0 or node_limit < 1:
        raise ValueError('invalid finite enumeration budget')
    words = [0]*n
    best, counts = {}, {}
    nodes = 0
    limit = bound * oracle.norm_scale
    masks = [sum((int(x) % 2) << j for j, x in enumerate(row)) for row in transform]

    def visit(i, used):
        nonlocal nodes
        nodes += 1
        if nodes > node_limit:
            raise RuntimeError('generic shell node budget exhausted; no complete bank')
        if i < 0:
            q, rem = divmod(used, oracle.norm_scale)
            if rem:
                raise ArithmeticError('nonintegral generic norm')
            counts[q] = counts.get(q, 0) + 1
            mask = 0
            for j, x in enumerate(words):
                if x % 2:
                    mask ^= masks[j]
            if mask == 0 or (mask in best and best[mask][0] < q):
                return
            w = tuple(sum(words[j]*transform[j][k] for j in range(n)) for k in range(n))
            if next(x for x in w if x) < 0:
                w = tuple(-x for x in w)
            if mask not in best or (q, w) < best[mask]:
                best[mask] = q, w
            return
        b = oracle.shift_denominators[i]
        a = sum(c * words[j] for j, c in enumerate(oracle.shift_coefficients[i], i + 1))
        weight = oracle.integer_weights[i]
        rad = isqrt((limit - used) // weight)
        lo, hi = -((rad + a) // b), (rad - a) // b
        for x in range(lo, hi + 1):
            words[i] = x
            visit(i - 1, used + weight*(b*x + a)**2)

    visit(n - 1, 0)
    return {'status': 'COMPLETE_EXACT_GENERIC_SHELL_BANK', 'dimension': n,
            'bound': bound, 'shells': list(shells), 'nodes': nodes,
            'norm_counts': counts,
            'rows': [{'mask': m, 'norm': q, 'word': list(w)}
                     for m, (q, w) in sorted(best.items()) if q in shells],
            'represented_parity_minima': {str(m): q for m, (q, w) in sorted(best.items())}}
