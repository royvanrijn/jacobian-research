#!/usr/bin/env python3
"""Exact compact controls for arbitrary-cycle branched-cover rank bounds.

Standard library only. The generic builder accepts a finite group
presentation, right permutation actions, and one specified meridian word.
It is not a knot-presentation or meridian-provenance recognizer.
"""
from itertools import product
from fractions import Fraction
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def inverse(p):
    return tuple(p.index(i) for i in range(len(p)))


def cycles(p):
    seen, answer = set(), []
    for start in range(len(p)):
        if start in seen:
            continue
        cycle, v = [], start
        while v not in seen:
            seen.add(v)
            cycle.append(v)
            v = p[v]
        assert v == start
        answer.append(cycle)
    return answer


def covering_matrix(images, relators, meridian):
    """Lift all relators and fill each cycle of one meridian word.

    Words are traversed left to right; images are right actions on sheets.
    Edge ng+v is the positive generator g edge from v to images[g][v].
    Delete a spanning tree and retain all other directed edge coordinates.
    """
    assert images and meridian
    n = len(images[0])
    assert n >= 2
    for p in images:
        assert sorted(p) == list(range(n)), 'invalid sheet permutation'
    inverses = [inverse(p) for p in images]
    edge_count = n * len(images)

    def walk(word, start):
        row, v = [0] * edge_count, start
        for letter in word:
            assert letter and abs(letter) <= len(images)
            g = abs(letter) - 1
            if letter < 0:
                v = inverses[g][v]
            row[g*n+v] += 1 if letter > 0 else -1
            if letter > 0:
                v = images[g][v]
        return v, row

    parent = list(range(n))
    def root(v):
        while parent[v] != v:
            v = parent[v]
        return v
    tree = []
    for g, p in enumerate(images):
        for v in range(n):
            if root(v) != root(p[v]):
                parent[root(v)] = root(p[v])
                tree.append(g*n+v)
    assert len(tree) == n-1, 'nontransitive sheet action'
    retained = [e for e in range(edge_count) if e not in tree]
    rows = []
    for word in relators:
        for start in range(n):
            end, row = walk(word, start)
            assert end == start, 'permutation relator fails'
            rows.append([row[e] for e in retained])
    meridian_image = tuple(walk(meridian, v)[0] for v in range(n))
    peripheral_cycles = cycles(meridian_image)
    for cycle in peripheral_cycles:
        end, row = walk(meridian * len(cycle), cycle[0])
        assert end == cycle[0]
        rows.append([row[e] for e in retained])
    return {'matrix': rows, 'tree_edges': tree, 'retained_edges': retained,
            'meridian_image': list(meridian_image),
            'meridian_cycles': peripheral_cycles,
            'degree': n, 'delta': n-len(peripheral_cycles)}


def determinant(matrix):
    """Exact rational elimination, independent of integer Smith routines."""
    a = [[Fraction(x) for x in row] for row in matrix]
    n, answer = len(a), Fraction(1)
    assert all(len(row) == n for row in a)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            answer = -answer
        value = a[j][j]
        answer *= value
        for i in range(j+1, n):
            coefficient = a[i][j] / value
            for k in range(j+1, n):
                a[i][k] -= coefficient * a[j][k]
    assert answer.denominator == 1
    return int(answer)


def binary_characters(matrix):
    return [list(v) for v in product(range(2), repeat=len(matrix[0]))
            if all(sum(x*y for x,y in zip(v,row)) % 2 == 0
                   for row in matrix)]


def matrix_action(a):
    points = [tuple((i >> j) & 1 for j in range(3)) for i in range(1,8)]
    return tuple(points.index(tuple(sum(a[i][j]*v[j] for j in range(3)) % 2
                                    for i in range(3))) for v in points)


def multiply(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def power(p, exponent):
    answer = tuple(range(len(p)))
    if exponent < 0:
        return power(inverse(p), -exponent)
    for _ in range(exponent):
        answer = multiply(answer, p)
    return answer


def expect_failure(callback, message):
    try:
        callback()
    except AssertionError as error:
        assert message in str(error)
        return True
    raise AssertionError('corrupted input accepted')


def main():
    x = matrix_action([[1,0,0],[1,1,0],[0,1,1]])
    y = matrix_action([[0,1,0],[0,0,1],[1,1,0]])
    images = [inverse(x), inverse(y)]
    identity = tuple(range(7))
    m = multiply(inverse(x), power(y, 2))
    assert power(x,4) == power(y,7) == power(m,2) == identity
    z = multiply(m, multiply(y,m))
    assert multiply(y,z) != multiply(z,y)
    controls = []
    for t in (0,1,2,3):
        q, a = 28*t+7, 7*t+2
        relators, meridian = [[1]*4+[-2]*q], [-1]+[2]*a
        cover = covering_matrix(images, relators, meridian)
        assert cover['degree'] == 7 and cover['delta'] == 2
        chars = binary_characters(cover['matrix'])
        assert len(chars) == 2
        assert cover['meridian_image'] == list(m)
        # Independent nonabelian quotient of the double-cover presentation.
        yz = multiply(power(y,a),power(z,a))
        zy = multiply(power(z,a),power(y,a))
        assert power(yz,2) == power(y,q)
        assert power(zy,2) == power(z,q)
        controls.append({'t':t, 'q':q, 'a':a,
                         'H1_F2_dimension':1,
                         'binary_characters':chars,
                         'matrix_mod2':[[x % 2 for x in row]
                                        for row in cover['matrix']]})
        if t == 0:
            first = cover
    assert controls[0]['matrix_mod2'] == controls[2]['matrix_mod2']
    assert controls[1]['matrix_mod2'] == controls[3]['matrix_mod2']
    rows = first['matrix']
    assert first['tree_edges'] == [0,1,2,4,8,10]
    minor_rows = [0,1,3,7,8,9,10,11]
    minor = determinant([rows[i] for i in minor_rows])
    assert minor == 2
    assert binary_characters(rows) == [[0]*8,[1,1,0,1,1,0,1,1]]
    # A determinant-two full-rank minor bounds the abelian group order by2;
    # the explicit nonzero binary character forces order2, hence H1=C2.
    # Deliberately omit fixed-sheet fillings: this must change homology.
    moved_rows = rows[:7] + [rows[7+i] for i,c in
                            enumerate(first['meridian_cycles']) if len(c)>1]
    assert len(binary_characters(moved_rows)) > 2
    controls_rejected = {
        'invalid_permutation': expect_failure(
            lambda: covering_matrix([tuple([0]*7),images[1]],[[1]*4+[-2]*7],[-1,2,2]),
            'invalid sheet permutation'),
        'invalid_relator': expect_failure(
            lambda: covering_matrix(images,[[1]*3+[-2]*7],[-1,2,2]),
            'permutation relator fails'),
        'nontransitive_action': expect_failure(
            lambda: covering_matrix([identity],[],[1]),
            'nontransitive sheet action'),
        'omitting_fixed_fillings_changes_homology': True}
    inputs = ['verify_general_branched_covers.py',
              'GENERAL_BRANCHED_COVER_RANK_BOUNDS.tex']
    receipt = {
        'status':'PASS',
        'scope':'Exact compact torus controls, integral H1=C2 at T(4,7), binary parity controls, and nonabelian double-cover quotients. Uniform graph and parity proofs remain written; no novelty or full MRC claim.',
        'T4_7':first,
        'integral_H1_certificate':{'minor_rows':minor_rows,'minor_determinant':minor,
                                 'H1_free_rank':0,'H1_torsion_factors':[2]},
        'known_bound_comparison':{'simple_cover_rank_invariant':3,
                                 'same_action_index_only':3,
                                 'general_cover_H1_bound':4,
                                 'bridge_upper_bound':4},
        'parity_controls':controls,
        'negative_controls':controls_rejected,
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest()
                        for p in inputs}}
    output = HERE/'general_branched_cover_verification.json'
    output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','T4_7_H1':'Z/2',
                      'bound':4,'parity_controls':4,
                      'negative_controls':len(controls_rejected)},indent=2))


if __name__ == '__main__':
    main()
