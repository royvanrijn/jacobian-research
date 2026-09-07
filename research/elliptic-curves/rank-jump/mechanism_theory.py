"""Exact small arithmetic gates for RANK_JUMP_MECHANISM_THEOREMS.md.

These functions certify algebra conditional on supplied point relations; they
never certify that a caller's coordinates are points or an independent basis.
No parameter enumeration, factorization, descent, or numerical heights.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, isqrt, prod


def qrank(rows, width=None):
    rows = [list(map(F, row)) for row in rows]
    width = len(rows[0]) if rows and width is None else (width or 0)
    if any(len(row) != width for row in rows):
        raise ValueError("ragged matrix")
    pivots = {}
    for row in rows:
        for j, pivot in sorted(pivots.items()):
            a = row[j]
            row = [x - a*y for x, y in zip(row, pivot)]
        if any(row):
            j = next(j for j, x in enumerate(row) if x)
            a = row[j]
            pivots[j] = [x/a for x in row]
    return len(pivots)


def quotient_dimension(generic_rows, point_rows):
    """Ranks of rational coordinates in a separately certified ambient basis."""
    rows = list(generic_rows) + list(point_rows)
    width = len(rows[0]) if rows else 0
    return qrank(rows, width) - qrank(generic_rows, width)


def signed_relation_bound(n, edges):
    """Relations a*P_i+b*P_j in the generic rational span, a,b in {+1,-1}."""
    if not isinstance(n, int) or n < 0:
        raise ValueError("invalid vertex count")
    graph = [[] for _ in range(n)]
    rows = []
    for i, j, a, b in edges:
        if not (0 <= i < n and 0 <= j < n and i != j and a in (-1, 1) and b in (-1, 1)):
            raise ValueError("invalid signed edge")
        row = [0]*n
        row[i], row[j] = a, b
        rows.append(row)
        graph[i].append((j, -a*b))
        graph[j].append((i, -a*b))
    signs = {}
    components = []
    for root in range(n):
        if root in signs:
            continue
        signs[root] = 1
        todo, vertices, balanced = [root], [], True
        while todo:
            i = todo.pop()
            vertices.append(i)
            for j, sign in graph[i]:
                expected = sign*signs[i]
                if j not in signs:
                    signs[j] = expected
                    todo.append(j)
                elif signs[j] != expected:
                    balanced = False
        components.append({"vertices": sorted(vertices), "balanced": balanced})
    bound = sum(c["balanced"] for c in components)
    rank = qrank(rows, n)
    if bound != n-rank:
        raise ArithmeticError("graph and rational matrix disagree")
    return {"relation_rank": rank, "quotient_upper_bound": bound, "components": components}


def determinant(rows):
    a = [list(map(F, row)) for row in rows]
    n = len(a)
    if any(len(row) != n for row in a):
        raise ValueError("nonsquare matrix")
    out = F(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            return F(0)
        if pivot != j:
            a[j], a[pivot] = a[pivot], a[j]
            out = -out
        p = a[j][j]
        out *= p
        for i in range(j+1, n):
            c = a[i][j]/p
            for k in range(j+1, n):
                a[i][k] -= c*a[j][k]
    return out


def resultant(f, g):
    """Sylvester determinant, coefficients from constant to leading."""
    if len(f) < 2 or len(g) < 2 or not f[-1] or not g[-1]:
        raise ValueError("positive exact polynomial degrees required")
    m, n = len(f)-1, len(g)-1
    rows = [[0]*i + list(reversed(f)) + [0]*(n-1-i) for i in range(n)]
    rows += [[0]*i + list(reversed(g)) + [0]*(m-1-i) for i in range(m)]
    value = determinant(rows)
    if value.denominator != 1:
        raise ValueError("integral coefficients required")
    return int(value)


def quadratic_support(forms):
    if len(forms) < 2:
        raise ValueError("at least two quadratic forms required")
    for q in forms:
        if len(q) != 3 or any(type(x) is not int for x in q) or q[2] == 0:
            raise ValueError("exact integral quadratics required; retain all constants")
        if q[1]*q[1] == 4*q[0]*q[2]:
            raise ValueError("branch divisor is not squarefree")
    resultants = [resultant(forms[i], forms[j]) for i, j in combinations(range(len(forms)), 2)]
    if not all(resultants):
        raise ValueError("branch divisors collide over Qbar")
    return abs(prod(resultants))


def strip_support(value, support):
    """Return (supported factor, coprime factor), without prime factorization."""
    if value <= 0 or support <= 0:
        raise ValueError("positive integers required")
    remaining, supported = value, 1
    while (d := gcd(remaining, support)) != 1:
        supported *= d
        remaining //= d
    return supported, remaining


def lift_at_parameter(forms, T, Z):
    """Exact homogeneous lift at one primitive projective address.

    BRANCH_REQUIRES_SEPARATE_CHART means this nonbranch theorem is inapplicable.
    NOT_PRODUCT_POINT and NONTRIVIAL_NATIVE_CLASS exclude only this address.
    """
    support = quadratic_support(forms)
    if type(T) is not int or type(Z) is not int or gcd(T, Z) != 1 or Z < 0:
        raise ValueError("primitive projective integer parameter with Z >= 0 required")
    values = [q[0]*Z*Z + q[1]*T*Z + q[2]*T*T for q in forms]
    if 0 in values:
        return {"status": "BRANCH_REQUIRES_SEPARATE_CHART"}
    product = prod(values)
    if product < 0 or isqrt(product)**2 != product:
        return {"status": "NOT_PRODUCT_POINT"}
    representatives, outside_roots = [], []
    for value in values:
        supported, outside = strip_support(abs(value), support)
        if isqrt(outside)**2 != outside:
            raise ArithmeticError("resultant support theorem violated")
        representatives.append(supported if value > 0 else -supported)
        outside_roots.append(isqrt(outside))
    split = all(v > 0 and isqrt(v)**2 == v for v in representatives)
    return {"status": "NATIVE_LIFT" if split else "NONTRIVIAL_NATIVE_CLASS",
            "signed_supported_representatives": representatives,
            "outside_square_roots": outside_roots,
            "homogeneous_roots": [isqrt(v) for v in values] if split else None}


def minimum_square_tests(m, clusters):
    """Minimum vertex covers of collision cliques; valuation support only.

    This small structural enumeration is capped at 16 labels. It is not a
    parameter search. A separate real-sign condition remains necessary.
    """
    if not 2 <= m <= 16:
        raise ValueError("declared structural limit is 2..16 labels")
    if any(c < 0 or c >= 1 << m for c in clusters):
        raise ValueError("invalid collision mask")
    edges = {sum(1 << i for i in pair) for c in clusters
             for pair in combinations([i for i in range(m) if c >> i & 1], 2)}
    for size in range(m):
        covers = [sum(1 << i for i in subset) for subset in combinations(range(m), size)
                  if all(any(edge >> i & 1 for i in subset) for edge in edges)]
        if covers:
            return {"minimum_tests": size, "index_masks": covers, "pair_masks": sorted(edges)}
    raise ArithmeticError("a graph on m vertices has a vertex cover of size <= m-1")
