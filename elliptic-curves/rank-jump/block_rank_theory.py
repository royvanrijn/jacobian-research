"""Exact rank-addition and partial-pairing gates; arithmetic provenance is required.

Matrices here do not by themselves prove point membership, Selmer membership,
or an upper envelope. Callers must supply the certificates stated in J7-J10.
"""
from math import isqrt

from mechanism_theory import qrank, quotient_dimension


def mod_rref(rows, p, width=None):
    if type(p) is not int or p < 2 or any(p % d == 0 for d in range(2, isqrt(p)+1)):
        raise ValueError('prime modulus required')
    rows = [list(row) for row in rows]
    width = len(rows[0]) if rows and width is None else (width or 0)
    if any(len(row) != width or any(type(x) is not int for x in row) for row in rows):
        raise ValueError('rectangular integer matrix required')
    rows = [[x % p for x in row] for row in rows]
    pivots, k = [], 0
    for j in range(width):
        i = next((i for i in range(k, len(rows)) if rows[i][j]), None)
        if i is None:
            continue
        rows[k], rows[i] = rows[i], rows[k]
        inv = pow(rows[k][j], -1, p)
        rows[k] = [x*inv % p for x in rows[k]]
        for i in range(len(rows)):
            if i != k:
                c = rows[i][j]
                rows[i] = [(x-c*y) % p for x, y in zip(rows[i], rows[k])]
        pivots.append(j)
        k += 1
    return rows[:k], pivots


def mod_rank(rows, p=2, width=None):
    return len(mod_rref(rows, p, width)[1])


def nullspace(rows, p, width):
    reduced, pivots = mod_rref(rows, p, width)
    result = []
    for j in range(width):
        if j in pivots:
            continue
        v = [int(i == j) for i in range(width)]
        for row, k in zip(reduced, pivots):
            v[k] = -row[j] % p
        result.append(v)
    return result


def signature_sandwich(generic_rank, torsion_dimension, generic, points, relations, p=2):
    if any(type(x) is not int or x < 0 for x in (generic_rank, torsion_dimension)):
        raise ValueError('nonnegative certified rank and torsion dimension required')
    if not generic and not points:
        width = 0
    else:
        width = len((generic+points)[0])
    a = mod_rank(generic, p, width)
    b = mod_rank(generic+points, p, width)
    if a > generic_rank+torsion_dimension:
        raise ValueError('generic signature rank contradicts supplied rank/torsion bound')
    lower = max(0, b-generic_rank-torsion_dimension)
    upper = len(points)-qrank(relations, len(points))
    if lower > upper:
        raise ValueError('independence and relation certificates contradict each other')
    return {'generic_signature_rank': a, 'combined_signature_rank': b,
            'observed_signature_increment': b-a,
            'unresolved_generic_signature_defect_budget': generic_rank+torsion_dimension-a,
            'quotient_lower_bound': lower, 'quotient_upper_bound': upper,
            'exact_quotient_rank': lower if lower == upper else None}


def block_overlap(generic, first, second):
    a = quotient_dimension(generic, first)
    b = quotient_dimension(generic, second)
    union = quotient_dimension(generic, first+second)
    return {'first_rank': a, 'second_rank': b, 'union_rank': union,
            'intersection_rank': a+b-union, 'second_increment_after_first': union-a}


def alternating(matrix):
    n = len(matrix)
    if any(len(row) != n or any(type(x) is not int or x not in (0, 1) for x in row) for row in matrix):
        raise ValueError('square binary matrix required')
    if any(matrix[i][i] or matrix[i][j] != matrix[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Cassels-Tate input must be alternating')


def radical_partner_bound(old, cross, partner_count=None):
    """J9 uses old-old and old-new pairings, never new-new pairings."""
    alternating(old)
    d = len(old)
    if len(cross) != d:
        raise ValueError('one cross row per old basis vector required')
    e = len(cross[0]) if cross else (partner_count or 0)
    if partner_count is not None and partner_count != e:
        raise ValueError('partner width mismatch')
    if any(len(row) != e or any(type(x) is not int or x not in (0, 1) for x in row) for row in cross):
        raise ValueError('binary cross-pairing matrix required')
    radical = nullspace(old, 2, d)
    obstruction = [[sum(v[i]*cross[i][j] for i in range(d)) % 2 for j in range(e)] for v in radical]
    a, c = mod_rank(obstruction, 2, e), mod_rank(old, 2, d)
    return {'old_dimension': d, 'old_pairing_rank': c, 'old_radical_basis': radical,
            'radical_partner_matrix': obstruction, 'radical_partner_rank': a,
            'certified_pairing_rank_lower_bound': c+2*a,
            'old_subspace_soluble_dimension_upper_bound': d-c-a,
            'new_independent_selmer_classes_mod_old_at_least': a}


def rank_exclusion(selmer_upper, torsion_dimension, pairing_lower, target, known_rank_lower=0):
    if any(type(x) is not int or x < 0 for x in (torsion_dimension, pairing_lower, target, known_rank_lower)):
        raise ValueError('nonnegative integer bounds required')
    if pairing_lower % 2:
        raise ValueError('use an even alternating-rank lower bound')
    if selmer_upper is None:
        return {'status': 'UNKNOWN_NO_SELMER_UPPER_ENVELOPE', 'rank_upper_bound': None}
    if type(selmer_upper) is not int or selmer_upper < 0:
        raise ValueError('a certified nonnegative integer Selmer upper bound is required')
    upper = selmer_upper-torsion_dimension-pairing_lower
    if upper < known_rank_lower:
        raise ValueError('upper inputs contradict a certified rank lower bound')
    needed = max(0, selmer_upper-torsion_dimension-target+1)
    needed += needed % 2
    return {'status': 'TARGET_EXCLUDED' if upper < target else 'NOT_EXCLUDED',
            'rank_upper_bound': upper,
            'even_pairing_rank_sufficient_for_exclusion': needed,
            'additional_pairing_rank_needed': max(0, needed-pairing_lower)}
