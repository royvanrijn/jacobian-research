#!/usr/bin/env python3
"""Exact local braid/Gram controls; stdlib only, not a knot-ideal solver."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
N = 4
PAIRS = [(i, j) for i in range(N) for j in range(i + 1, N)]
ZERO = (0,) * len(PAIRS)


def c(value):
    return {ZERO: value} if value else {}


def add(*ps):
    result = {}
    for p in ps:
        for e, a in p.items():
            result[e] = result.get(e, 0) + a
            if not result[e]:
                del result[e]
    return result


def neg(p):
    return {e: -a for e, a in p.items()}


def mul(p, q):
    result = {}
    for e, a in p.items():
        for f, b in q.items():
            g = tuple(x + y for x, y in zip(e, f))
            result[g] = result.get(g, 0) + a * b
            if not result[g]:
                del result[g]
    return result


def unit(i):
    return [c(int(i == j)) for j in range(N)]


def dot(a, b, gram):
    return add(*(mul(mul(x, gram[i][j]), y)
                 for i, x in enumerate(a) for j, y in enumerate(b)))


def reflect(normal, vector, gram):
    factor = dot(normal, vector, gram)
    return [add(x, neg(mul(factor, y)))
            for x, y in zip(vector, normal)]


def gram_matrix():
    gram = [[c(2 if i == j else 0) for j in range(N)] for i in range(N)]
    for index, (i, j) in enumerate(PAIRS):
        e = list(ZERO)
        e[index] = 1
        gram[i][j] = gram[j][i] = {tuple(e): 1}
    return gram


def root_move(roots, k, sign, gram):
    answer = [list(v) for v in roots]
    a, b = roots[k:k + 2]
    if sign == 1:
        answer[k:k + 2] = [reflect(a, b, gram), a]
    else:
        assert sign == -1
        answer[k:k + 2] = [b, reflect(b, a, gram)]
    return answer


def cord_move(g, k, sign):
    """Independent symmetric specialization of the source's phi rules."""
    h = [[dict(x) for x in row] for row in g]
    other = [i for i in range(N) if i not in (k, k + 1)]
    h[k][k + 1] = h[k + 1][k] = neg(g[k][k + 1])
    for i in other:
        if sign == 1:
            h[k][i] = h[i][k] = add(g[k + 1][i],
                                    neg(mul(g[k + 1][k], g[k][i])))
            h[k + 1][i] = h[i][k + 1] = g[k][i]
        else:
            h[k][i] = h[i][k] = g[k + 1][i]
            h[k + 1][i] = h[i][k + 1] = add(g[k][i],
                                            neg(mul(g[k][k + 1], g[k + 1][i])))
    return h


def word_roots(word, gram):
    roots = [unit(i) for i in range(N)]
    for letter in word:
        roots = root_move(roots, abs(letter) - 1,
                          1 if letter > 0 else -1, gram)
    return roots


def run():
    source = HERE / 'sources/cornwell-kch-1310.7526v3.pdf'
    expected = '991940273daed005b8af270233cceaa58c11d69168efb7006f3bfa9bd598bcbd'
    assert hashlib.sha256(source.read_bytes()).hexdigest() == expected
    gram = gram_matrix()
    base = [unit(i) for i in range(N)]
    moves = 0
    for k in range(N - 1):
        for sign in (-1, 1):
            roots = root_move(base, k, sign, gram)
            measured = [[dot(a, b, gram) for b in roots] for a in roots]
            assert measured == cord_move(gram, k, sign)
            assert all(dot(v, v, gram) == c(2) for v in roots)
            assert root_move(roots, k, -sign, gram) == base
            # Reflection conjugacy is checked on every ambient basis vector.
            a, b = base[k:k + 2]
            for vector in base:
                assert reflect(reflect(a, b, gram), vector, gram) == reflect(
                    a, reflect(b, reflect(a, vector, gram), gram), gram)
            moves += 1
    for sign in (-1, 1):
        for k in (1, 2):
            assert word_roots([sign*k, sign*(k+1), sign*k], gram) == word_roots(
                [sign*(k+1), sign*k, sign*(k+1)], gram)
        assert word_roots([sign, sign*3], gram) == word_roots([sign*3, sign], gram)
    singular = [[c(2) for _ in range(N)] for _ in range(N)]
    for i, root in enumerate(base):
        displacement = [add(x, neg(y))
                        for x, y in zip(reflect(root, root, singular), root)]
        assert displacement == [c(-2 if i == j else 0) for j in range(N)]
    assert singular[0] == singular[1] == singular[2] == singular[3]
    return {
        'status': 'PASS',
        'scope': 'local symmetric braid/Gram identities and singular-form control',
        'generic_roots': N,
        'generic_gram_variables': len(PAIRS),
        'positive_negative_moves': moves,
        'inverse_cancellations': moves,
        'adjacent_braid_relations': 4,
        'disjoint_commutations': 2,
        'singular_gram_rank': 1,
        'independent_displacement_directions': N,
        'source_sha256': expected,
        'knot_specific_obstruction_replayed': False,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


if __name__ == '__main__':
    receipt = run()
    (HERE / 'symmetric_augmentation_transport_verification.json').write_text(
        json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))
