#!/usr/bin/env python3
"""Exact helpers for the Curve302 short-vector/core experiment.

No point search, no floating ranking and no Sage dependency.  The only external
math dependency is SymPy for exact integer/rational lattice normal forms.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from functools import reduce
import math
from pathlib import Path
from typing import Iterable, Iterator, Sequence

from sympy import Matrix, QQ, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp


def require(test, message):
    if not test:
        raise ValueError(message)


def fraction(value) -> F:
    if isinstance(value, F):
        return value
    if isinstance(value, int):
        return F(value)
    return F(str(value))


def primitive(values: Sequence[int]) -> tuple[int, ...]:
    values = tuple(int(v) for v in values)
    g = 0
    for value in values:
        g = math.gcd(g, abs(value))
    if not g:
        return values
    sign = 1 if next(v for v in values if v) > 0 else -1
    return tuple(sign * v // g for v in values)


def qnorm(form: Sequence[Sequence[F]], values: Sequence[int]) -> F:
    values = tuple(int(v) for v in values)
    return sum(F(values[i]) * form[i][j] * values[j]
               for i in range(len(values)) for j in range(len(values)))


def exact_ldl(form: Sequence[Sequence[F]]):
    """Return unit-lower L and positive diagonal D with Q=L D L^T."""
    q = [[fraction(v) for v in row] for row in form]
    n = len(q)
    require(n and all(len(row) == n for row in q), "form is not square")
    require(all(q[i][j] == q[j][i] for i in range(n) for j in range(n)), "form is asymmetric")
    L = [[F(int(i == j)) for j in range(n)] for i in range(n)]
    D = []
    for i in range(n):
        d = q[i][i] - sum(L[i][k] * L[i][k] * D[k] for k in range(i))
        require(d > 0, "form is not positive definite")
        D.append(d)
        for j in range(i + 1, n):
            L[j][i] = (q[j][i] - sum(L[j][k] * L[i][k] * D[k] for k in range(i))) / d
    return L, D


def floor_div(a: int, b: int) -> int:
    require(b > 0, "positive denominator required")
    return a // b


def ceil_div(a: int, b: int) -> int:
    require(b > 0, "positive denominator required")
    return -((-a) // b)


@dataclass
class EnumerationStats:
    nodes: int = 0
    leaves: int = 0
    directions: int = 0


def enumerate_primitive_directions(form: Sequence[Sequence[F]], bound: F, *,
                                   max_directions: int = 2_000_000,
                                   max_nodes: int = 100_000_000) -> tuple[list[tuple[F, tuple[int, ...]]], EnumerationStats]:
    """Complete exact Fincke-Pohst enumeration modulo sign through ``bound``.

    The recursion uses the exact rational LDL decomposition, so interval bounds
    are obtained by integer square roots and cannot lose a boundary vector to
    floating roundoff.  Every emitted vector is primitive and sign-canonical.
    ``max_*`` are fail-closed resource guards, not truncation limits.
    """
    bound = fraction(bound)
    require(bound >= 0, "negative enumeration bound")
    L, D = exact_ldl(form)
    n = len(D)
    x = [0] * n
    answer: list[tuple[F, tuple[int, ...]]] = []
    stats = EnumerationStats()

    def recurse(i: int, remaining: F):
        stats.nodes += 1
        require(stats.nodes <= max_nodes, "enumeration node cap reached; result UNKNOWN")
        if i < 0:
            stats.leaves += 1
            if not any(x):
                return
            g = 0
            first = 0
            for value in x:
                g = math.gcd(g, abs(value))
                if first == 0 and value:
                    first = value
            if g != 1 or first < 0:
                return
            vec = tuple(x)
            norm = qnorm(form, vec)
            require(norm <= bound, "enumerator emitted vector beyond bound")
            answer.append((norm, vec))
            stats.directions += 1
            require(stats.directions <= max_directions, "enumeration direction cap reached; result UNKNOWN")
            return

        c = sum(L[j][i] * x[j] for j in range(i + 1, n))
        radius2 = remaining / D[i]
        if radius2 < 0:
            return
        # |k + a/b|^2 <= p/q  => |k*b+a| <= floor(sqrt(p*b^2/q)).
        a, b = c.numerator, c.denominator
        p, q = radius2.numerator, radius2.denominator
        t2 = (p * b * b) // q
        T = math.isqrt(t2)
        lo = ceil_div(-T - a, b)
        hi = floor_div(T - a, b)
        for value in range(lo, hi + 1):
            x[i] = value
            term = D[i] * (F(value) + c) ** 2
            if term <= remaining:
                recurse(i - 1, remaining - term)
        x[i] = 0

    recurse(n - 1, bound)
    answer.sort(key=lambda row: (row[0], row[1]))
    require(len({v for _, v in answer}) == len(answer), "duplicate primitive direction")
    return answer, stats


def _matrix_rows(rows: Sequence[Sequence[int]], n: int) -> Matrix:
    if not rows:
        return Matrix.zeros(0, n)
    M = Matrix([[int(v) for v in row] for row in rows])
    require(M.cols == n, "wrong lattice ambient dimension")
    return M


def rational_rank(rows: Sequence[Sequence[int]], n: int) -> int:
    return int(_matrix_rows(rows, n).rank())


def saturation_index(rows: Sequence[Sequence[int]], n: int) -> int:
    """Index of the row lattice in its saturation inside its rational span."""
    A = _matrix_rows(rows, n)
    if not A.rows:
        return 1
    rank = int(A.rank())
    D = smith_normal_form(A, domain=ZZ)
    diagonal = [abs(int(D[i, i])) for i in range(min(D.rows, D.cols)) if D[i, i] != 0]
    require(len(diagonal) == rank, "Smith rank mismatch")
    index = math.prod(diagonal)
    require(index >= 1, "invalid saturation index")
    return index


def _primitive_integer_row(values) -> list[int]:
    den = 1
    for value in values:
        den = math.lcm(den, int(value.q)) if hasattr(value, "q") else math.lcm(den, fraction(value).denominator)
    ints = [int(value * den) for value in values]
    g = 0
    for value in ints:
        g = math.gcd(g, abs(value))
    if g:
        ints = [value // g for value in ints]
    if next((v for v in ints if v), 1) < 0:
        ints = [-v for v in ints]
    return ints


def _integer_kernel_basis(equations: Matrix, n: int) -> Matrix:
    """Saturated integer right kernel via Smith transformations."""
    if equations.rows == 0:
        return Matrix.eye(n)
    require(equations.cols == n, "kernel equation width mismatch")
    domain = DomainMatrix.from_Matrix(equations)
    D, left, right = smith_normal_decomp(domain)
    check = left.to_Matrix() * equations * right.to_Matrix()
    require(check == D.to_Matrix(), "Smith transformation identity failed")
    rank = sum(D.to_Matrix()[i, i] != 0 for i in range(min(D.shape)))
    T = right.to_Matrix()
    # S*N*T=D, so the columns of T indexed by zero Smith columns form the
    # saturated integer right kernel of N. Return them as row generators.
    return Matrix([list(T[:, j]) for j in range(rank, n)])


def saturation_basis(rows: Sequence[Sequence[int]], n: int) -> tuple[tuple[int, ...], ...]:
    """Canonical HNF row basis of span_Q(rows) intersect Z^n."""
    A = _matrix_rows(rows, n)
    if not A.rows or A.rank() == 0:
        return ()
    null = A.nullspace()
    if null:
        annihilator = Matrix([_primitive_integer_row(list(v)) for v in null])
        kernel = _integer_kernel_basis(annihilator, n)
    else:
        kernel = Matrix.eye(n)
    require(kernel.cols == n, "integer kernel width mismatch")
    H = hermite_normal_form(kernel.T).T
    basis = tuple(tuple(int(v) for v in H.row(i)) for i in range(H.rows)
                  if any(H.row(i)))
    require(rational_rank(basis, n) == A.rank(), "saturation basis rank mismatch")
    require(saturation_index(basis, n) == 1, "computed saturation basis is not saturated")
    return basis


def lattice_contains(outer: Sequence[Sequence[int]], inner: Sequence[Sequence[int]], n: int) -> bool:
    """Containment for saturated lattices, equivalent to rational-span containment."""
    if not inner:
        return True
    r = rational_rank(outer, n)
    return rational_rank((*outer, *inner), n) == r


def intersection_saturated(lattices: Sequence[Sequence[Sequence[int]]], n: int) -> tuple[tuple[int, ...], ...]:
    """Intersection of saturated sublattices of Z^n, returned saturated/canonical."""
    require(lattices, "no lattices to intersect")
    annihilator_rows = []
    for rows in lattices:
        A = _matrix_rows(rows, n)
        require(saturation_index(rows, n) == 1, "intersection input is not saturated")
        for v in A.nullspace():
            annihilator_rows.append(_primitive_integer_row(list(v)))
    if not annihilator_rows:
        return tuple(tuple(int(i == j) for j in range(n)) for i in range(n))
    annihilator = Matrix(annihilator_rows)
    kernel = _integer_kernel_basis(annihilator, n)
    return saturation_basis(kernel.tolist(), n)


def shell_counts(rows: Sequence[tuple[F, tuple[int, ...]]]) -> list[tuple[F, int, int]]:
    """Return (norm, first 1-based rank, last 1-based rank) for exact norm shells."""
    answer = []
    i = 0
    while i < len(rows):
        j = i + 1
        while j < len(rows) and rows[j][0] == rows[i][0]:
            j += 1
        answer.append((rows[i][0], i + 1, j))
        i = j
    return answer


def rank_interval(rows: Sequence[tuple[F, tuple[int, ...]]], vector: Sequence[int]) -> tuple[int, int, F]:
    target = primitive(vector)
    lookup = {v: i for i, (_, v) in enumerate(rows)}
    require(target in lookup, "observed vector absent from complete enumeration")
    i = lookup[target]
    norm = rows[i][0]
    lo = i
    while lo and rows[lo - 1][0] == norm:
        lo -= 1
    hi = i + 1
    while hi < len(rows) and rows[hi][0] == norm:
        hi += 1
    return lo + 1, hi, norm


def annihilator_rows(rows: Sequence[Sequence[int]], n: int) -> tuple[tuple[F, ...], ...]:
    A = _matrix_rows(rows, n)
    return tuple(tuple(F(int(x.p), int(x.q)) if hasattr(x, "p") else fraction(x) for x in v)
                 for v in A.nullspace())


def in_rational_span(vector: Sequence[int], rows: Sequence[Sequence[int]], n: int) -> bool:
    if not rows:
        return not any(vector)
    return rational_rank(rows, n) == rational_rank((*rows, tuple(vector)), n)


def one_step_core_hit(prefix: Sequence[Sequence[int]], candidate: Sequence[int],
                      core: Sequence[Sequence[int]], n: int) -> bool:
    """Whether Sat(prefix + candidate) contains ``core``.

    Prefixes and cores are saturated in the actual experiment, so testing the
    rational span is exact and avoids one Smith/HNF computation per candidate.
    """
    require(saturation_index(prefix, n) == 1, "prefix must be saturated")
    require(saturation_index(core, n) == 1, "core must be saturated")
    if in_rational_span(candidate, prefix, n):
        return False
    return lattice_contains(saturation_basis((*prefix, tuple(candidate)), n), core, n)
