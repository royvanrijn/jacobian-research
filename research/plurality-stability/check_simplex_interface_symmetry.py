#!/usr/bin/env python3
"""Exact finite checks supporting SIMPLEX_INTERFACE_SYMMETRY.md.

No all-k induction theorem, analytic kernel inequality or global optimum
is certified by this script. The optional quadrature is exploration only.
"""

import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import itertools
import json
import math
from pathlib import Path
from statistics import NormalDist


def parity(p):
    return (-1) ** sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))


def cycle_type(p):
    seen, lengths = set(), []
    for i in range(len(p)):
        if i in seen:
            continue
        j, length = i, 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = p[j]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def edge_trace(p, complement_sign=False):
    """Induced action on an oriented reference edge and its complement."""
    n, result = len(p), 0
    for i, j in itertools.combinations(range(n), 2):
        if {p[i], p[j]} != {i, j}:
            continue
        sign = 1 if p[i] < p[j] else -1
        if complement_sign:
            remaining = [a for a in range(n) if a not in (i, j)]
            sign *= parity([remaining.index(p[a]) for a in remaining])
        result += sign
    return result


def character_checks():
    # Actual edge traces are compared with independently specified irreducibles.
    s4 = {
        (1, 1, 1, 1): {"31": 3, "22": 2, "211": 3, "1111": 1},
        (2, 1, 1): {"31": 1, "22": 0, "211": -1, "1111": -1},
        (2, 2): {"31": -1, "22": 2, "211": -1, "1111": 1},
        (3, 1): {"31": 0, "22": -1, "211": 0, "1111": 1},
        (4,): {"31": -1, "22": 0, "211": 1, "1111": -1},
    }
    checked = 0
    for p in itertools.permutations(range(4)):
        chi = s4[cycle_type(p)]
        assert edge_trace(p) == chi["31"] + chi["211"]
        assert edge_trace(p, True) == chi["22"] + chi["211"] + chi["1111"]
        checked += 1
    for p in itertools.permutations(range(3)):
        standard = sum(p[i] == i for i in range(3)) - 1
        assert edge_trace(p) == standard + parity(p)
        checked += 1
    # Negative control: forgetting the orientation is already wrong at (12).
    p = (1, 0, 2, 3)
    unsigned = sum({p[i], p[j]} == {i, j} for i, j in itertools.combinations(range(4), 2))
    assert unsigned != edge_trace(p)
    return {"permutations_checked": checked, "orientation_negative_control": "passed"}


def higher_cell_reflection_checks():
    """Compare the geometric reflection sum and determinant as formal sums.

    A term is exp(tau * a^T Q b), represented by its integer matrix Q.
    No transcendental evaluation or spectral conclusion is involved.
    """
    cases = []
    for k in range(3, 9):
        p = k - 1

        def coordinate(i):
            return 0 if i < 2 else i - 1

        def geometric_exponent(permutation):
            matrix = [0] * (p * p)
            for i, j in enumerate(permutation):
                matrix[p * coordinate(i) + coordinate(j)] += 1
            return tuple(matrix)

        def add_term(terms, matrix, coefficient):
            terms[matrix] = terms.get(matrix, 0) + coefficient
            if terms[matrix] == 0:
                del terms[matrix]

        folded, omitted, wrong_sign, standard_folded = {}, {}, {}, {}
        for rest in itertools.permutations(range(2, k)):
            permutation = [0, 1, *rest]
            sign = parity(rest)
            exponent = geometric_exponent(permutation)
            for terms in (folded, omitted, wrong_sign):
                add_term(terms, exponent, 2 * sign)
            add_term(standard_folded, exponent, 2)
            for l in range(2, k):
                reflected = permutation.copy()
                reflected[1], reflected[l] = reflected[l], reflected[1]
                exponent = geometric_exponent(reflected)
                add_term(folded, exponent, -2 * sign)
                if l != k - 1:
                    add_term(omitted, exponent, -2 * sign)
                add_term(wrong_sign, exponent, 2 * sign)
                add_term(standard_folded, exponent, 1)

        determinant, permanent = {}, {}
        for permutation in itertools.permutations(range(p)):
            matrix = [0] * (p * p)
            matrix[0] = 1  # The additional exp(tau*x_2*y_2) factor.
            for i, j in enumerate(permutation):
                matrix[p * i + j] += 1
            add_term(determinant, tuple(matrix), 2 * parity(permutation))
            add_term(permanent, tuple(matrix), 1 + (permutation[0] == 0))

        assert folded == determinant
        assert standard_folded == permanent
        assert standard_folded != determinant
        assert len(determinant) == math.factorial(p)
        assert omitted != determinant
        assert wrong_sign != determinant
        # The prefolding kernel's value at the origin, divided by Gamma.
        origin_coefficient = 2 * (3 - k)
        assert origin_coefficient == 2 - 2 * (k - 2)
        if k >= 4:
            assert origin_coefficient < 0
        cases.append({"k": k, "determinant_terms": len(determinant),
                      "prefolding_origin_coefficient": origin_coefficient})

    return {
        "cases": cases,
        "omitted_adjacent_and_wrong_sign_controls": "passed",
        "standard_permanent_identity": "passed",
        "scope": "Exact bilinear-exponent monomial comparisons for the sign determinant S44/S46 and standard permanent S49/S50 at 3<=k<=8. The all-k identities, determinant positivity, operator norm and dimension lifting require the written proofs; no all-k spectral inequality is certified.",
    }


def explicit_gap_check():
    t = Q(15, 16)
    assert t - t**3 / 6 > Q(4, 5)
    assert 6 * 20**2 < 49**2
    assert 21 * 11**2 > 50**2
    upper = Q(64, 9) * (Q(13, 24) - Q(1240, 2695) + Q(451, 2750))
    assert 1 < upper < Q(7, 4)
    assert not upper < Q(3, 2)  # Detect an unsupported stronger bound.
    return {
        "H_half_rational_upper": str(upper),
        "H_half_upper_less_than_7_over_4": True,
        "scope": "Rational arithmetic in S21; uses analytic Gaussian integral and classical 3<pi<22/7 from the proof.",
        "gap": "ell*(2/sqrt(3)-1), sign sector at rho=1/2 only",
    }


def hook_cone_checks():
    """Exact bounded controls for S52-S61, separate from their universal proof."""
    def dot(v, w):
        return sum(a * b for a, b in zip(v, w))

    def swap(v, i, j):
        w = list(v)
        w[i], w[j] = w[j], w[i]
        return w

    def wedge(v, w):
        return [v[i] * w[j] - v[j] * w[i]
                for i, j in itertools.combinations(range(len(v)), 2)]

    def symmetric(v, w):
        return [v[i] * w[j] + v[j] * w[i]
                for i in range(len(v)) for j in range(i, len(v))]

    def adjacent(v, l):
        # I - 2*w_l*w_l^T on H12; all selected vectors are integral.
        assert v[0] == v[1]
        s = v[0] - v[l]
        w = list(v)
        w[0] -= s
        w[1] -= s
        w[l] += 2 * s
        return w

    cases, two_row_cases = [], []
    wedge_pairs = 0
    for k in range(3, 8):
        rays = [[k-m] * m + [-m] * (k-m) for m in range(2, k)]
        root = [1, -1] + [0] * (k-2)
        for ray in rays:
            for l in range(2, k):
                first = wedge(swap(root, 1, l), swap(ray, 1, l))
                second = wedge(swap(root, 0, l), swap(ray, 0, l))
                assert [a+b for a, b in zip(first, second)] == wedge(root, adjacent(ray, l))
                wedge_pairs += 1
        for m, ray in enumerate(rays, 2):
            for j, other in enumerate(rays, 2):
                assert dot(ray, other) == k * min(j, m) * (k-max(j, m)) > 0

        x = [sum(ray[i] for ray in rays) for i in range(k)]
        y = [sum((m+1) * ray[i] for m, ray in enumerate(rays)) for i in range(k)]
        assert x[0] == x[1] and y[0] == y[1]
        assert all(x[i] > x[i+1] and y[i] > y[i+1] for i in range(1, k-1))
        terms = []
        for rest in itertools.permutations(range(2, k)):
            perm = [0, 1, *rest]
            z = [y[i] for i in perm]
            terms.append((dot(x, z), perm, None))
            for l in range(2, k):
                terms.append((dot(x, swap(z, 1, l)), perm, l))
        minimum = min(e for e, _, _ in terms)
        powers = {e: 1 << (e-minimum) for e, _, _ in terms}
        matrix = [[0] * (k-2) for _ in range(k-2)]
        for exponent, perm, l in terms:
            for m, ray in enumerate(rays):
                v = [ray[i] for i in perm]
                v = [2*a for a in v] if l is None else adjacent(v, l)
                for j in range(k-2):
                    difference = v[j+1] - v[j+2]
                    assert difference % k == 0
                    matrix[j][m] += (difference // k) * powers[exponent]
        assert all(entry > 0 for row in matrix for entry in row)
        cases.append({"k": k, "matrix_size": k-2, "positive_entries": (k-2)**2,
                      "tau": "log(2)", "x_gap_coordinates": [1] * (k-2),
                      "y_gap_coordinates": list(range(1, k-1))})
        if k >= 4:
            q = k-2
            rays0 = [[0, 0] + [q-m]*m + [-m]*(q-m) for m in range(1, q)]
            for m, ray in enumerate(rays0, 1):
                for j, other in enumerate(rays0, 1):
                    assert dot(ray, other) == q*min(j, m)*(q-max(j, m)) > 0
                for l in range(2, k):
                    first = symmetric(swap(root, 1, l), swap(ray, 1, l))
                    second = symmetric(swap(root, 0, l), swap(ray, 0, l))
                    assert [a+b for a, b in zip(first, second)] == symmetric(root, ray)
            matrix0 = [[0]*(q-1) for _ in range(q-1)]
            for exponent, perm, l in terms:
                for m, ray in enumerate(rays0):
                    v = [ray[i] for i in perm]
                    factor = 2 if l is None else 1
                    for j in range(q-1):
                        difference = v[j+2] - v[j+3]
                        assert difference % q == 0
                        matrix0[j][m] += factor*(difference//q)*powers[exponent]
            assert all(entry > 0 for row in matrix0 for entry in row)
            # Exact exponent identity forcing the ordered residual R_j-R_(j+1).
            for rest in itertools.permutations(range(2, k)):
                z = [y[i] for i in [0, 1, *rest]]
                for j in range(2, k-1):
                    first = dot(x, swap(z, 1, j))
                    second = dot(x, swap(swap(z, j, j+1), 1, j+1))
                    assert first-second == (x[j]-x[j+1])*(z[1]-z[j+1]) > 0
            if k == 4:
                # Scalar sign and (2,2) kernels have a strict positive difference.
                sign_sum = sum((2 if l is None else -2)*parity(perm)*powers[e]
                               for e, perm, l in terms)
                assert 0 < sign_sum < matrix0[0][0]
            two_row_cases.append({"k": k, "matrix_size": q-1,
                                  "positive_entries": (q-1)**2,
                                  "ordered_residual_exponent_identity": "passed"})

    # Domain control: the three-cell hook kernel is negative for x=-r2,y=r2.
    assert 2 - 2 * 2**9 < 0
    root, ray = [1, -1, 0], [1, 1, -2]
    assert wedge(root, adjacent(ray, 2)) != wedge(root, [-a for a in ray])
    return {
        "wedge_transport_pairs": wedge_pairs,
        "strict_integer_matrix_controls": cases,
        "two_row_symmetric_transport_and_matrices": two_row_cases,
        "four_cell_sign_domination_control": "passed",
        "outside_chamber_and_wrong_reflection_coefficient_controls": "passed",
        "scope": "Exact finite exterior/symmetric transport, Gram and ordered-residual exponent identities, selected strictly positive folded matrices for 3<=k<=7, and one exact four-cell sign domination control. Exponentials are exact powers of two after factoring a common positive term. These controls do not certify all x,y,k,rho; the universal pairing, supersolution and spectral arguments S54-S61 remain written mathematics.",
    }


def elementary_half_point_and_reflections():
    """The finite polynomial/rational arithmetic used in S32 and S36."""
    x = Q(93, 100)
    assert x - x**3/6 + x**5/120 - x**7/5040 > Q(4, 5)
    x = Q(101, 500)
    assert x - x**3/6 > Q(1, 5)
    x = Q(411, 1000)
    assert x - x**3/6 + x**5/120 < Q(2, 5)
    x = Q(253, 1000)
    assert x - x**3/6 > Q(1, 4)
    assert 16 * (Q(1, 5) - Q(1, 5)**3/3) - Q(4, 239) > Q(25, 8)
    assert 3 * Q(577, 1000)**2 < 1
    assert 15 * Q(258, 1000)**2 < 1
    assert 6 * Q(409, 1000)**2 > 1
    assert 21 * Q(219, 1000)**2 > 1
    R = Q(4, 3)*(1+2*Q(93, 100)/Q(25, 8))
    R += 2*Q(409, 1000)*(1+2*Q(101, 500)/Q(25, 8))
    R += Q(219, 1000)*(1-2*Q(411, 1000)/Q(22, 7))
    L = Q(4, 3)*Q(577, 1000) + Q(258, 1000)*(1-2*Q(253, 1000)/Q(25, 8))
    upper = Q(16, 9)*R - Q(16, 3)*L*L
    assert R == Q(2650234883, 825000000)
    assert L == Q(4619803, 4687500)
    assert upper == Q(32058163928093, 60424804687500)
    assert 0 < upper < Q(9, 16)
    # With u=w-1>=0, compare the exact two difference polynomials.
    plus = [Q(math.comb(4, j)) + (2 if j <= 1 else 0) - (3 if j == 0 else 0) for j in range(5)]
    minus = [Q(math.comb(4, j)) - (4 if j <= 1 else 0) + (3 if j == 0 else 0) for j in range(5)]
    assert plus == [0, 6, 6, 4, 1]
    assert minus == [0, 0, 6, 4, 1]
    assert all(x >= 0 for x in plus + minus)
    # Replacing 3/2 by 2 in K- >= c*D is false for every w>1.
    assert 4*(1-Q(2)) < 0
    return {
        "R_half_upper": str(R), "L_half_lower": str(L),
        "G_half_upper": str(upper), "G_half_less_than_9_over_16": True,
        "reflection_difference_coefficients_in_u": {"plus": [str(x) for x in plus], "minus": [str(x) for x in minus]},
        "scope": "Exact finite algebra in S32/S36. The all-correlation theorem additionally uses the written OU parity spectra, ground-state identities, interaction monotonicity and extra-dimensional Hermite decomposition; this script does not verify those analytic arguments.",
    }


def sign_vandermonde_controls():
    """Bounded exact normalization and minor controls for S63-S66."""
    def determinant(rows):
        matrix = [[Q(value) for value in row] for row in rows]
        value = Q(1)
        for i in range(len(matrix)):
            pivot = next((j for j in range(i, len(matrix)) if matrix[j][i]), None)
            if pivot is None:
                return Q(0)
            if pivot != i:
                matrix[i], matrix[pivot] = matrix[pivot], matrix[i]
                value = -value
            p = matrix[i][i]
            value *= p
            for j in range(i+1, len(matrix)):
                factor = matrix[j][i] / p
                for l in range(i+1, len(matrix)):
                    matrix[j][l] -= factor * matrix[i][l]
        return value

    def gaussian_moment(order, mean, variance):
        return sum(Q(math.factorial(order),
                     2**r * math.factorial(r) * math.factorial(order-2*r))
                   * variance**r * mean**(order-2*r)
                   for r in range(order//2+1))

    rho, variance = Q(1, 2), Q(3, 4)
    cases = []
    for q in range(1, 7):
        degree = q*(q-1)//2
        b = list(range(q, 0, -1))
        moments = [[gaussian_moment(q-1-j, rho*value, variance)
                    for j in range(q)] for value in b]
        vandermonde = math.prod(b[i]-b[j] for i in range(q) for j in range(i+1, q))
        expected = rho**degree * vandermonde
        assert determinant(moments) == expected > 0
        assert expected != rho**(degree-1) * vandermonde
        theta = rho**(degree-1)
        assert (1-theta > 0) == (q >= 3)

        p = q+1
        a, c = list(range(p, 0, -1)), [j*j for j in range(p, 0, -1)]
        matrix = [[2**(x*y) for y in c] for x in a]
        full, tail = determinant(matrix), determinant([row[1:] for row in matrix[1:]])
        assert 0 < full < matrix[0][0]*tail
        cases.append({"q_untied": q, "k": q+2, "vandermonde_degree": degree,
                      "rho": str(rho), "variance": str(variance),
                      "gaussian_vandermonde_expectation": str(expected),
                      "sign_schur_factor": str(theta),
                      "strict_explicit_gap_factor": str(1-theta) if q >= 3 else None,
                      "exponential_minor_comparison": "passed"})

    # All entries positive is insufficient for the totally positive comparison.
    bad = [[1, 1, 100], [100, 10, 1], [1, 1, 10]]
    assert determinant(bad) == 8100
    assert determinant(bad) > bad[0][0] * determinant([row[1:] for row in bad[1:]])
    return {
        "cases": cases,
        "positive_entries_without_total_positivity_control": "passed",
        "wrong_rho_degree_and_k3_k4_scope_controls": "passed",
        "scope": "Exact selected Gaussian Vandermonde moments and exponential-matrix minor comparisons for 1<=q<=6, and rejection of insufficient hypotheses or the wrong rho exponent. This does not certify the imported Hermitian transition theorem, the universal Schur-complement argument, the Gaussian facet lift or the all-parameter spectral estimate.",
    }


def universal_block_controls():
    """Exact bounded controls for the group-algebra and Young-basis reduction."""
    group_cases = []
    for k in range(3, 8):
        counts = {}
        def left(g):
            return tuple(1-a if a < 2 else a for a in g)
        def right(g):
            return (g[1], g[0], *g[2:])
        for rest in itertools.permutations(range(1, k)):
            g = (0, *rest)
            for term, coefficient in ((g, 1), (left(right(g)), 1),
                                      (left(g), -1), (right(g), -1)):
                counts[term] = counts.get(term, 0) + coefficient
        for g in itertools.permutations(range(k)):
            expected = int(g[0] == 0)-int(g[0] == 1)-int(g[1] == 0)+int(g[1] == 1)
            assert counts.get(g, 0) == expected
        assert counts[(1, 0, *range(2, k))] == -2
        group_cases.append({"k": k, "permutations": math.factorial(k)})

    @lru_cache(None)
    def tableaux(shape):
        if not shape:
            return ((),)
        result = []
        for row, length in enumerate(shape):
            if row+1 < len(shape) and shape[row+1] == length:
                continue
            smaller = list(shape)
            smaller[row] -= 1
            smaller = tuple(x for x in smaller if x)
            result.extend(t + ((row, length-1),) for t in tableaux(smaller))
        return tuple(result)

    def partitions(n, cap=None):
        if n == 0:
            yield ()
        else:
            for first in range(min(n, n if cap is None else cap), 0, -1):
                for rest in partitions(n-first, first):
                    yield (first, *rest)

    def identity(d):
        return [[Q(i == j) for j in range(d)] for i in range(d)]

    def product(a, b):
        return [[sum(value*b[r][j] for r, value in enumerate(row) if value)
                 for j in range(len(b))] for row in a]

    young_cases, triple_sign_cases, triple_frame_cases = [], [], []
    shapes = [shape for k in range(3, 6) for shape in partitions(k)]
    shapes += [(3, 2, 1), (2, 2, 2), (3, 1, 1, 1), (2, 2, 1, 1), (3, 3), (4, 2)]
    for shape in shapes:
        k, tabs = sum(shape), tableaux(shape)
        d, index = len(tabs), {t: i for i, t in enumerate(tabs)}
        lengths = []
        for t in tabs:
            word = sorted(range(k), key=lambda i: t[i])
            lengths.append(sum(word[i] > word[j] for i in range(k) for j in range(i+1, k)))
        generators = []
        for a in range(k-1):
            matrix = [[Q(0) for _ in range(d)] for _ in range(d)]
            for col, t in enumerate(tabs):
                axial = (t[a+1][1]-t[a+1][0])-(t[a][1]-t[a][0])
                assert axial
                matrix[col][col] = Q(1, axial)
                swapped = list(t)
                swapped[a], swapped[a+1] = swapped[a+1], swapped[a]
                other = index.get(tuple(swapped))
                if other is not None:
                    matrix[other][col] = Q(1) if lengths[other] > lengths[col] else 1-Q(1, axial**2)
                else:
                    assert abs(axial) == 1
            assert all(value+(i == j) >= 0 for i, row in enumerate(matrix) for j, value in enumerate(row))
            assert product(matrix, matrix) == identity(d)
            generators.append(matrix)
        for i in range(k-2):
            a, b = generators[i:i+2]
            assert product(a, product(b, a)) == product(b, product(a, b))
        for i in range(k-1):
            for j in range(i+2, k-1):
                assert product(generators[i], generators[j]) == product(generators[j], generators[i])
        minus = [i for i in range(d) if generators[0][i][i] == -1]
        assert all(generators[0][i][j] == 0 for i in range(d) for j in range(d) if i != j)
        first, second = generators[:2]
        fs, sf = product(first, second), product(second, first)
        third = product(fs, first)
        antisym = [[(Q(i == j)-first[i][j]-second[i][j]-third[i][j]+fs[i][j]+sf[i][j])/6
                    for j in range(d)] for i in range(d)]
        assert product(antisym, antisym) == antisym
        for i in range(d):
            for j in range(d):
                corner = 2*second[i][j] if i in minus and j in minus else Q(0)
                assert corner == Q(i == j and i in minus)-3*antisym[i][j]
        rank = sum(antisym[i][i] for i in range(d))
        assert rank.denominator == 1
        assert (rank == 0) == (len(shape) <= 2)
        triple_sign_cases.append({"partition": list(shape), "minus_dimension": len(minus),
                                  "three_label_sign_rank": int(rank)})
        if minus:
            def swap_rep(i, j):
                current = identity(d)
                word = [*range(i, j), *range(j-2, i-1, -1)]
                for a in word:
                    current = product(current, generators[a])
                return [[current[a][b] for b in minus] for a in minus]

            size = len(minus)
            projections = []
            for l in range(2, k):
                transport = swap_rep(1, l)
                projections.append([[(Q(i == j)-2*transport[i][j])/3 for j in range(size)]
                                    for i in range(size)])
            frame = [[sum(a[i][j] for a in projections) for j in range(size)] for i in range(size)]
            remaining_center = [[Q(0) for _ in minus] for _ in minus]
            for i, j in itertools.combinations(range(2, k), 2):
                matrix = swap_rep(i, j)
                for a in range(size):
                    for b in range(size):
                        remaining_center[a][b] += matrix[a][b]
            content = lambda beta: sum(length*(length-1)//2-i*length for i, length in enumerate(beta))
            for i in range(size):
                for j in range(size):
                    assert 3*frame[i][j] == (k-3-content(shape))*int(i == j)+remaining_center[i][j]
            roots = set()
            for i, j in itertools.combinations(range(len(shape)), 2):
                beta = list(shape)
                beta[i] -= 1
                beta[j] -= 1
                if all(a >= b for a, b in zip(beta, beta[1:])):
                    roots.add(Q(k-3-content(shape)+content(beta), 3))
            assert roots
            polynomial = identity(size)
            for value in sorted(roots):
                polynomial = product(polynomial, [[frame[i][j]-value*int(i == j) for j in range(size)]
                                                  for i in range(size)])
            assert polynomial == [[Q(0) for _ in minus] for _ in minus]
            if len(shape) >= 3:
                assert min(roots) >= Q(1, 3)
            if len(shape) == 3:
                assert roots <= {Q(shape[2], 3), Q(shape[1]+1, 3), Q(shape[0]+2, 3)}
            if k >= 4:
                four_sign = [[Q(i == j and all(tabs[minus[i]][a] == (a, 0) for a in range(4)))
                              for j in range(size)] for i in range(size)]
                a3, a4 = projections[:2]
                overlap = product(product(a3, a4), a3)
                assert overlap == [[a3[i][j]/9+Q(8, 9)*four_sign[i][j] for j in range(size)]
                                   for i in range(size)]
            triple_frame_cases.append({"partition": list(shape), "frame_eigenvalue_controls": [str(v) for v in sorted(roots)],
                                       "projection_frame_identity": "passed"})
        z = [[Q(0) for _ in range(d)] for _ in range(d)]
        for rest in itertools.permutations(range(1, k)):
            g = [0, *rest]
            inverse = [g.index(i) for i in range(k)]
            exponent = sum((k-i)*(k-inverse[i])**2 for i in range(1, k))
            weight = 2**exponent
            current, word, work = identity(d), [], list(g)
            for target in range(k):
                pos = work.index(target)
                while pos > target:
                    work[pos-1], work[pos] = work[pos], work[pos-1]
                    word.append(pos-1)
                    pos -= 1
            for a in word:
                current = product(generators[a], current)
            for i in range(d):
                for j in range(d):
                    z[i][j] += weight*current[i][j]
        assert all(value >= 0 for row in z for value in row)
        assert all(z[i][j] >= 0 for i in minus for j in minus)
        young_cases.append({"partition": list(shape), "k": k, "dimension": d,
                            "minus_dimension": len(minus), "selected_TP_matrix": "passed"})
    assert 1*1-2*3 < 0  # A positive-entry, non-TN matrix has a negative sign immanant.
    return {
        "free_group_algebra_identity": group_cases,
        "rational_Young_Coxeter_and_selected_TP_controls": young_cases,
        "triple_sign_corner_controls": triple_sign_cases,
        "triple_projection_frame_controls": triple_frame_cases,
        "orientation_and_non_TN_negative_controls": "passed",
        "scope": "Exact S68 coefficient identity for 3<=k<=7; exact rational seminormal Coxeter relations and nonnegative I+s_i matrices for every partition through k=5 and six selected k=6 partitions; selected nonsymmetric TP exponential-matrix sums and their (12)-minus corners. The same finite models verify S85's triple-sign corner identity, projection idempotence and height-two vanishing, plus S88-S89's S4 overlap, transposition-content frame identity and frame spectral polynomials. The all-k immanant/transport formulas, positivity, branching and scalar spectral reduction are written arguments; no analytic spectral gap is certified.",
    }


def triple_sign_reflection_controls():
    """Bounded exact half-facet controls; all-k spectral conclusions are written."""
    factor_count = 0
    for alpha in (0, 1, 3):
        for beta in (1, 2):
            for u in (0, 2):
                for t in (1, 3):
                    a, b, c = Q(2)**(alpha*t), Q(2)**(beta*t), Q(2)**(-beta*u)
                    actual = a*b+c-1-a*c
                    assert actual == (a-1)*(b-c)+(b-1) > 0
                    factor_count += 1

    def swap(value, i, j):
        value = list(value)
        value[i], value[j] = value[j], value[i]
        return value

    def center(value):
        mean = Q(sum(value), len(value))
        return [Q(v)-mean for v in value]

    def dot(x, y):
        return sum(a*b for a, b in zip(x, y))

    half_cases = []
    for k in range(4, 9):
        # The other coordinates deliberately need not be in decreasing order.
        x = center([10, 10, 5, 1]+[6 if j % 2 else -3 for j in range(k-4)])
        y = center([8, 8, 4, -2]+[-4 if j % 2 else 6 for j in range(k-4)])
        sy = swap(y, 2, 3)
        assert x[0] == x[1] > max(x[2:])
        assert y[0] == y[1] > max(y[2:])
        assert x[2] > x[3] and y[2] > y[3]
        terms = [(2, k*dot(x, y)), (-2, k*dot(x, sy))]
        for l in range(2, k):
            terms.extend(((1, k*dot(x, swap(y, 1, l))),
                          (-1, k*dot(x, swap(sy, 1, l)))))
        minimum = min(exponent for _, exponent in terms)
        assert all(exponent.denominator == 1 for _, exponent in terms)
        actual = sum(coefficient*(1 << int(exponent-minimum)) for coefficient, exponent in terms)
        assert actual > 0
        # Exchanging the y-half reverses the odd kernel sign exactly.
        assert sum(-coefficient*(1 << int(exponent-minimum)) for coefficient, exponent in terms) < 0
        # The row-rate ordering holds even when y3<y4 after the substitution.
        for z in (y, sy):
            sz = swap(z, 2, 3)
            exponent_gap = dot(x, swap(z, 1, 2))-dot(x, swap(sz, 1, 3))
            assert exponent_gap == (x[2]-x[3])*(z[1]-z[3]) > 0
        half_cases.append({"k": k, "tau": "k*log(2)", "strict_odd_kernel": "passed",
                           "other_coordinates_fully_ordered": all(x[j] > x[j+1] for j in range(2, k-1))})

    # A positive scalar exponential kernel with negative bilinear exponent is
    # not a positive-semidefinite operator: tau=4*log(2), x=s*w_l, s=1,2.
    determinant = Q(2)**(-3-12)-Q(2)**(-6-6)
    assert determinant == -Q(7, 32768) < 0
    assert Q(1)/Q(1, 2) > 1  # The invariant translation mode is excluded.
    return {"exceptional_pair_factorizations": factor_count,
            "half_facet_and_rate_order_controls": half_cases,
            "wrong_half_facet_rejection": "passed",
            "adjacent_kernel_negative_two_point_determinant": str(determinant),
            "invariant_translation_exclusion": "passed",
            "scope": "Exact exceptional-pair factorization controls and selected odd-kernel/rate exponent comparisons for4<=k<=8 without ordering all remaining coordinates, plus a negative adjacent-kernel determinant. Universal Young-model triple-sign identities are in exact_checks.universal_blocks. The all-k strict supersolution, nontrivial scalar type transfer and all two-row block gaps S85-S87 are written analytic and representation proofs, not certified by these finite controls."}


def antisymmetry_height_controls():
    """Finite formal-group and witness controls for S88-S89."""
    def multiply(a, b):
        result = {}
        for p, first in a.items():
            for q, second in b.items():
                key = tuple(p[q[i]] for i in range(len(p)))
                result[key] = result.get(key, Q(0))+first*second
        return {p: v for p, v in result.items() if v}

    def combine(*terms):
        result = {}
        for coefficient, value in terms:
            for p, v in value.items():
                result[p] = result.get(p, Q(0))+coefficient*v
        return {p: v for p, v in result.items() if v}

    def transposition(n, i, j):
        value = list(range(n))
        value[i], value[j] = value[j], value[i]
        return {tuple(value): Q(1)}

    def antisym(n, size):
        return {(*p, *range(size, n)): Q(parity(p), math.factorial(size))
                for p in itertools.permutations(range(size))}

    group_cases = []
    for m in range(2, 6):
        n = m+1
        e, larger = antisym(n, m), antisym(n, n)
        t = transposition(n, m-1, m)
        compressed = multiply(multiply(e, t), e)
        assert compressed == combine((Q(1, m), e), (-Q(m+1, m), larger))
        if m <= 4:
            other = multiply(multiply(t, e), t)
            overlap = multiply(multiply(e, other), e)
            assert overlap == combine((Q(1, m*m), e), (1-Q(1, m*m), larger))
        ident = {tuple(range(n)): Q(1)}
        p = combine((Q(1, 2), ident), (-Q(1, 2), transposition(n, 0, 1)))
        transport = combine((Q(2), multiply(multiply(p, transposition(n, 1, m)), p)))
        actual = multiply(e, transport)
        terms = [(Q(1), e), (-Q(m+1), larger)]
        terms.extend((-Q(1), multiply(e, transposition(n, r, m))) for r in range(2, m))
        assert actual == combine(*terms)
        group_cases.append({"antisymmetric_labels": m, "ambient_group_order": math.factorial(n),
                            "compressed_corner_and_transport": "passed",
                            "overlap_directly_expanded": m <= 4})

    def swap(v, i, j):
        v = list(v)
        v[i], v[j] = v[j], v[i]
        return v

    def scalar_terms(x, y):
        dot = lambda v: sum(a*b for a, b in zip(x, v))
        result = [(2, dot(y)), (-2, dot(swap(y, 1, 2)))]
        for l in range(3, len(x)):
            result.extend(((1, dot(swap(y, 1, l))),
                           (-1, dot(swap(swap(y, 2, l), 1, l)))))
        return result

    x = list(map(Q, (7, 7, -6, 0, -8)))
    y = [Q(v, 5) for v in (27, 27, 7, -23, -38)]
    assert sum(x) == sum(y) == 0
    assert x[0] == x[1] > max(x[2:]) and y[0] == y[1] > max(y[2:])
    assert x[3] > x[4] and y[3] > y[4]
    terms = scalar_terms(x, y)+[(-c, e) for c, e in scalar_terms(x, swap(y, 3, 4))]
    coefficients = {}
    for c, exponent in terms:
        exponent *= 5
        assert exponent.denominator == 1
        exponent = int(exponent)
        coefficients[exponent] = coefficients.get(exponent, 0)+c
    coefficients = {e: c for e, c in coefficients.items() if c}
    assert max(coefficients) == 680 and coefficients[680] == -1
    assert max(e for e, c in coefficients.items() if c > 0) == 650
    assert sum(c for c in coefficients.values() if c > 0) == 8
    assert -Q(1)+8*Q(2)**(-30) < 0
    smallest = min(coefficients)
    assert sum(c*2**(e-smallest) for e, c in coefficients.items()) < 0
    return {"nested_group_identities": group_cases,
            "negative_three_row_scalar_reflection": {
                "k": 5, "x": [str(v) for v in x], "y": [str(v) for v in y],
                "tau": "5*log(2)", "integer_exponential_coefficients":
                [[e, c] for e, c in sorted(coefficients.items(), reverse=True)],
                "leading_exponent": 680, "positive_exponent_upper": 650,
                "positive_coefficient_sum": 8, "sign": "negative"},
            "scope": "Exact nested antisymmetrizer corner/transport identities in S3 through S6, direct overlap expansions through S5, and the exact five-cell negative reflection witness. The all-size coset hierarchy and all-k three-row closed-image spectral reduction S88-S89 are written proofs, using the retained content-spectrum input. No new all-correlation gap follows from these controls."}


def standard_hook_gaussian_controls():
    """Finite exact algebra behind the Gaussian integration identities S71."""
    cases = []
    for k in range(3, 9):
        vectors = {l: [Q(1, 2) if i < 2 else Q(-int(i == l)) for i in range(k)]
                   for l in range(2, k)}
        dot = lambda x, y: sum(a*b for a, b in zip(x, y))
        outer = lambda w: [[a*b for b in w] for a in w]
        apply = lambda matrix, x: [sum(a*b for a, b in zip(row, x)) for row in matrix]
        matrices = {}
        for l, w in vectors.items():
            assert dot(w, w) == Q(3, 2)
            ww = outer(w)
            c = [[Q(i == j)-2*ww[i][j] for j in range(k)] for i in range(k)]
            d = [[Q(i == j)-ww[i][j] for j in range(k)] for i in range(k)]
            matrices[l] = c
            for i in range(k):
                for j in range(k):
                    assert sum(c[i][r]*d[r][j] for r in range(k)) == (i == j)
                    assert d[i][j]-c[i][j] == ww[i][j]
            assert apply(c, w) == [-2*a for a in w]
            wrong_reflection = [[Q(i == j)-Q(4, 3)*ww[i][j] for j in range(k)]
                                for i in range(k)]
            assert apply(wrong_reflection, apply(d, w)) != w
        for l, w in vectors.items():
            for r, other in vectors.items():
                if l == r:
                    continue
                assert dot(w, other) == Q(1, 2)
                assert apply(matrices[r], w) == [a-b for a, b in zip(w, other)]
                assert apply(matrices[l], other) == [b-a for a, b in zip(w, other)]
                # A basis of vectors with y1=y2=yl; equality holds before centering.
                for index in [0, *[i for i in range(2, k) if i != l]]:
                    y = [Q(int(i == index)) for i in range(k)]
                    if index == 0:
                        y[1] = y[l] = Q(1)
                    left = y[:]
                    left[1], left[r] = left[r], left[1]
                    right = y[:]
                    right[l], right[r] = right[r], right[l]
                    right[1], right[l] = right[l], right[1]
                    assert left == right
        cases.append({"k": k, "adjacent_labels": k-2, "ordered_boundary_pairs": (k-2)*(k-3)})
    return {"cases": cases, "wrong_orthogonal_reflection_control": "passed",
            "scope": "Exact C_l D_l=I, D_l-C_l=w_l w_l^T, same-face annihilation and all paired triple-junction permutation identities for 3<=k<=8. The unrestricted Gaussian integration-by-parts and differentiation identities S71 are written analytic proofs; no favorable residual sign or standard gap is certified."}


def higher_exterior_hook_checks():
    """Bounded exact controls for S72 and the complete third-hook gap S73."""
    def dot(x, y):
        return sum(a*b for a, b in zip(x, y))

    def swap(v, i, j):
        result = list(v)
        result[i], result[j] = result[j], result[i]
        return result

    def wedge(columns):
        rank, size = len(columns), len(columns[0])
        return [sum(parity(p)*math.prod(columns[j][indices[p[j]]] for j in range(rank))
                    for p in itertools.permutations(range(rank)))
                for indices in itertools.combinations(range(size), rank)]

    def transport(value, w):
        return [a-2*b*dot(w, value) for a, b in zip(value, w)]

    def helmert(value, j):
        return sum(value[:j])-j*value[j]

    cases, folded_cases, transport_count, boundary_count = [], [], 0, 0
    for k in range(4, 9):
        root = [Q(1), Q(-1)]+[Q(0)]*(k-2)
        rays = [[Q(k-m)]*m+[Q(-m)]*(k-m) for m in range(2, k)]
        p0 = [Q(1, 2)-Q(1, k)]*2+[Q(-1, k)]*(k-2)
        x = [sum(ray[i] for ray in rays) for i in range(k)]
        y = [sum((m+1)*ray[i] for m, ray in enumerate(rays)) for i in range(k)]
        vectors = {l: [Q(1, 2)]*2+[Q(-int(i == l)) for i in range(2, k)]
                   for l in range(2, k)}
        rates = {l: Q(k-l) for l in vectors}
        v = [sum(rates[l]*vectors[l][i] for l in vectors) for i in range(k)]
        for l, w in vectors.items():
            assert transport(p0, w) == [a-b for a, b in zip(p0, w)]
            for first, second in itertools.combinations(rays, 2):
                actual = [a+b for a, b in zip(
                    wedge([swap(z, 1, l) for z in (root, first, second)]),
                    wedge([swap(z, 0, l) for z in (root, first, second)]))]
                assert actual == wedge([root, transport(first, w), transport(second, w)])
                transport_count += 1
            for r, other in vectors.items():
                if l == r:
                    continue
                difference = [a-b for a, b in zip(w, other)]
                for index in range(k):
                    basis = [Q(i == index) for i in range(k)]
                    first = transport(basis, other)
                    second = transport(swap(basis, l, r), w)
                    assert [a-b for a, b in zip(first, second)] == [a*(w[index]+other[index]) for a in difference]
                boundary_count += 1
        xx, vv = [helmert(x, j) for j in range(2, k)], [helmert(v, j) for j in range(2, k)]
        for j in range(2, k):
            assert helmert(p0, j) == 1
            assert vv[j-2] == (j+1)*rates[j]+sum(rates[l] for l in range(j+1, k))
        assert all(0 < a < b for a, b in zip(xx, xx[1:]))
        assert all(a > b > 0 for a, b in zip(vv, vv[1:]))
        assert all(vv[i]*xx[j]-vv[j]*xx[i] > 0 for i, j in itertools.combinations(range(k-2), 2))
        cases.append({"k": k, "multiplicity_dimension": math.comb(k-2, 2),
                      "positive_linear_supersolution_coordinates": "passed",
                      "ordered_residual_coordinates": "passed"})

        if k <= 6:
            basis = [[Q(1)]*j+[Q(-j)]+[Q(0)]*(k-j-1) for j in range(2, k)]
            pairs = list(itertools.combinations(range(k-2), 2))
            terms = []
            for rest in itertools.permutations(range(2, k)):
                perm = [0, 1, *rest]
                z = [y[i] for i in perm]
                terms.append((dot(x, z), perm, None))
                terms.extend((dot(x, swap(z, 1, l)), perm, l) for l in vectors)
            minimum = min(e for e, _, _ in terms)
            powers = {e: 1 << (int(e-minimum)) for e, _, _ in terms}
            matrix = [[Q(0) for _ in pairs] for _ in pairs]
            for exponent, perm, l in terms:
                linear = []
                for column in basis:
                    value = [column[i] for i in perm]
                    if l is not None:
                        value = transport(value, vectors[l])
                    linear.append([helmert(value, j)/Q(j*(j+1)) for j in range(2, k)])
                factor = (2 if l is None else 1)*powers[exponent]
                for a, (i, j) in enumerate(pairs):
                    for b, (r, s) in enumerate(pairs):
                        matrix[a][b] += factor*(linear[r][i]*linear[s][j]-linear[s][i]*linear[r][j])
            assert all(entry > 0 for row in matrix for entry in row)
            folded_cases.append({"k": k, "matrix_size": len(pairs),
                                 "strict_positive_entries": len(pairs)**2, "tau": "log(2)"})

    # Mere positivity of abstract row rates does not prove the ordered residual.
    assert 105*203-108*200 == -285 < 0
    return {"third_exterior_transport_pairs": transport_count,
            "rank_one_junction_pairs": boundary_count,
            "supersolution_coordinate_controls": cases,
            "selected_positive_folded_matrices": folded_cases,
            "unordered_positive_rates_control": "passed",
            "scope": "Exact third-exterior transport, rank-one junction-pair and selected ordered coordinate identities for 4<=k<=8; selected strict positive folded exterior matrices for 4<=k<=6. A control rejects positivity of abstract row rates without their required ordering. The all-degree Gaussian creation identity S72 and complete third-hook spectral gap S73 are written proofs, not certified by these finite controls."}


def representation_hermite_floor_checks():
    """Finite character controls for the written all-k floor and Schur bound."""
    @lru_cache(None)
    def partitions(total, cap=None):
        if total == 0:
            return ((),)
        cap = total if cap is None else min(cap, total)
        return tuple((first,)+rest for first in range(cap, 0, -1)
                     for rest in partitions(total-first, first))

    def cells(shape):
        return {(i, j) for i, length in enumerate(shape) for j in range(length)}

    @lru_cache(None)
    def character(shape, cycle):
        # Murnaghan--Nakayama, enumerating every contained remainder shape.
        if not cycle:
            return int(not shape)
        length, rest = cycle[0], cycle[1:]
        outer, answer = cells(shape), 0
        for remainder in partitions(sum(shape)-length):
            inner = cells(remainder)
            if not inner <= outer:
                continue
            strip = outer-inner
            if any({(i, j), (i+1, j), (i, j+1), (i+1, j+1)} <= strip
                   for i, j in strip):
                continue
            reached = {next(iter(strip))}
            pending = list(reached)
            while pending:
                i, j = pending.pop()
                for neighbor in ((i-1, j), (i+1, j), (i, j-1), (i, j+1)):
                    if neighbor in strip and neighbor not in reached:
                        reached.add(neighbor)
                        pending.append(neighbor)
            if reached == strip:
                answer += (-1)**(len({i for i, j in strip})-1)*character(remainder, rest)
        return answer

    def centralizer(cycle):
        counts = {v: cycle.count(v) for v in set(cycle)}
        return math.prod(v**count*math.factorial(count) for v, count in counts.items())

    def symmetric_powers(cycle, cap):
        # det(I-t*g|E)^(-1) = (1-t)*product_l (1-t^l)^(-1).
        values = [1]+[0]*cap
        for length in cycle:
            for degree in range(length, cap+1):
                values[degree] += values[degree-length]
        return [values[i]-(values[i-1] if i else 0) for i in range(cap+1)]

    cases, multiplicity_tests, character_norms = [], 0, 0
    for k in range(3, 10):
        shapes = partitions(k)
        cap = (k-1)*(k-2)//2
        powers = {cycle: symmetric_powers(cycle, cap) for cycle in shapes}
        for shape in shapes:
            norm = sum(Q(character(shape, cycle)**2, centralizer(cycle)) for cycle in shapes)
            assert norm == 1
            character_norms += 1
            if len(shape) == 1:
                continue
            degree = sum(i*length for i, length in enumerate(shape))-(len(shape)-1)
            if degree == 0:
                assert shape == (k-1, 1)
            if degree == 1:
                assert shape in ((k-2, 1, 1), (k-2, 2))
            multiplicities = [sum(Q(character(shape, cycle)*(cycle.count(1)-1)*powers[cycle][m],
                                      centralizer(cycle)) for cycle in shapes)
                              for m in range(degree+1)]
            assert all(value == 0 for value in multiplicities[:degree])
            assert multiplicities[degree] > 0
            assert all(value.denominator == 1 for value in multiplicities)
            multiplicity_tests += degree+1
            cases.append({"k": k, "partition": list(shape), "first_target_degree": degree,
                          "target_multiplicity": int(multiplicities[degree])})

    # The standard target contains degree zero: ordinary balance is separate.
    assert next(c for c in cases if c["k"] == 5 and c["partition"] == [4, 1])["first_target_degree"] == 0
    assert next(c for c in cases if c["k"] == 6 and c["partition"] == [3, 2, 1])["first_target_degree"] == 2
    assert next(c for c in cases if c["k"] == 8 and c["partition"] == [2, 2, 2, 2])["first_target_degree"] == 9

    # Exact rational squares of the midpoint-b and optimized-b errors.
    midpoint_cases = []
    for k, rho, degree in ((8, Q(1, 2), 9), (12, Q(1, 2), 15), (20, Q(3, 4), 50)):
        b = (1+rho)/2
        direct = 4*(k-1)**2*(1-rho*rho)*rho**(2*degree-2)/(b**(2*degree)*(1-b*b))
        simplified = Q(64*(k-1)**2, 1)/((1+rho)*(3+rho))*(2*rho/(1+rho))**(2*degree-2)
        assert direct == simplified
        assert simplified < 25*(k-1)**2*(2*rho/(1+rho))**(2*degree-2)
        midpoint_cases.append({"k": k, "rho": str(rho), "degree": degree,
                               "relative_error_square": str(simplified)})
    k, degree, rho = 8, 9, Q(1, 2)
    optimized = 4*(k-1)**2*(1-rho*rho)*Q((degree+1)**(degree+1), degree**degree)*rho**(2*degree-2)
    coarse = Q(9*(k-1)**2*(degree+1), 4**(degree-1))
    assert optimized < Q(1, 16)
    assert optimized < coarse == Q(4410, 65536) < Q(1, 9)
    # A shallow-sector input must not be declared positive by the same bound.
    shallow = 4*(6-1)**2*(1-rho*rho)*Q(3**3, 2**2)*rho**2
    assert shallow > 1
    # Squared balanced-remainder coefficient for b^2=2/3 in S83.
    assert Q(16, 1)/(Q(2, 3)**2*Q(1, 3)) == 108
    remainder_cases = []
    for rho, b in ((Q(1, 4), Q(1, 2)), (Q(1, 2), Q(3, 4))):
        for degree in range(2, 13):
            assert rho**(degree-1) <= rho/b*b**(degree-1)
        # The constant and linear terms must be removed by balance/subtraction.
        assert rho**(-1) > rho/b*b**(-1)
        assert Q(1) > rho/b
        remainder_cases.append({"rho": str(rho), "b": str(b), "degrees": [2, 12]})
    endpoint_cases = []
    for k in range(3, 21):
        # rho_endpoint^2, and the error divided by phi0 after sqrt(3) cancels.
        endpoint_square = Q(1, 108*k*k*(k-1)**2)
        error = Q(1, k*(k-1))
        ell_lower = Q(2, k-1)
        assert endpoint_square < Q(1, 4) < Q(2, 3)
        assert error == ell_lower/(2*k)
        endpoint_cases.append({"k": k, "endpoint_square": str(endpoint_square),
                               "error_divided_phi0": str(error)})
    return {"cell_range": [3, 9], "character_norm_controls": character_norms,
            "graded_multiplicity_tests": multiplicity_tests, "partition_cases": cases,
            "midpoint_bound_controls": midpoint_cases,
            "two_column_first_instance_optimized_error_square": str(optimized),
            "two_column_first_instance_coarse_error_square": str(coarse),
            "shallow_sector_rejection": "passed",
            "balanced_remainder_power_controls": remainder_cases,
            "explicit_initial_endpoint_controls": endpoint_cases,
            "unbalanced_and_unsubtracted_remainder_rejections": "passed",
            "scope": "Exact graded-character controls through k=9 for the first degree in E tensor Sym^m(E) and the D=0/1 classification, rational Schur-bound algebra, the first two-column instance, and bounded Hermite-power/explicit initial-endpoint controls. The m=0/1 rejection controls require balance and subtraction of the first moment. These controls do not assert interface attainment of the polynomial target degree. The all-k Hermite vanishing, Schur bound, infinite-family positive gap and convex-covariance first-moment/explicit initial full-gap deductions in S80-S83 are written proofs."}


def five_cell_mixed_certificate():
    """Exhaustive formal positivity certificate for the five-cell (2,2,1) block."""
    k, shape = 5, (2, 2, 1)

    @lru_cache(None)
    def tableaux(partition):
        if not partition:
            return ((),)
        result = []
        for row, length in enumerate(partition):
            if row+1 < len(partition) and partition[row+1] == length:
                continue
            smaller = list(partition)
            smaller[row] -= 1
            smaller = tuple(v for v in smaller if v)
            result.extend(t+((row, length-1),) for t in tableaux(smaller))
        return tuple(result)

    def identity(n):
        return [[Q(i == j) for j in range(n)] for i in range(n)]

    def product(a, b):
        return [[sum(v*b[r][j] for r, v in enumerate(row) if v)
                 for j in range(len(b[0]))] for row in a]

    def apply(matrix, vector):
        return [sum(a*b for a, b in zip(row, vector)) for row in matrix]

    tabs = tableaux(shape)
    lookup, dimension = {t: i for i, t in enumerate(tabs)}, len(tabs)
    lengths = []
    for t in tabs:
        word = sorted(range(k), key=lambda i: t[i])
        lengths.append(sum(word[i] > word[j] for i in range(k) for j in range(i+1, k)))
    generators = []
    for index in range(k-1):
        matrix = [[Q(0) for _ in tabs] for _ in tabs]
        for col, t in enumerate(tabs):
            axial = (t[index+1][1]-t[index+1][0])-(t[index][1]-t[index][0])
            matrix[col][col] = Q(1, axial)
            other = list(t)
            other[index], other[index+1] = other[index+1], other[index]
            row = lookup.get(tuple(other))
            if row is not None:
                matrix[row][col] = Q(1) if lengths[row] > lengths[col] else 1-Q(1, axial**2)
        assert product(matrix, matrix) == identity(dimension)
        generators.append(matrix)
    for i in range(k-2):
        a, b = generators[i:i+2]
        assert product(a, product(b, a)) == product(b, product(a, b))
    for i in range(k-1):
        for j in range(i+2, k-1):
            assert product(generators[i], generators[j]) == product(generators[j], generators[i])

    @lru_cache(None)
    def representation(g):
        result, work = identity(dimension), list(g)
        for target in range(k):
            pos = work.index(target)
            while pos > target:
                work[pos-1], work[pos] = work[pos], work[pos-1]
                result = product(generators[pos-1], result)
                pos -= 1
        return result

    minus = [i for i in range(dimension) if generators[0][i][i] == -1]
    assert minus == [0, 1, 3]
    restrict = lambda a: [[a[i][j] for j in minus] for i in minus]
    average = [[(Q(i == j)+restrict(generators[3])[i][j])/2 for j in range(3)] for i in range(3)]
    a3 = apply(average, [Q(1)]*3)
    linear, transports = [[Q(0)]*5 for _ in range(3)], {}
    for l in range(2, 5):
        g = list(range(5))
        g[2], g[l] = g[l], g[2]
        vector = apply(restrict(representation(tuple(g))), a3)
        for row, value in enumerate(vector):
            linear[row][l] = value
        g = list(range(5))
        g[1], g[l] = g[l], g[1]
        transports[l] = [[2*v for v in row] for row in restrict(representation(tuple(g))) ]
    assert linear == [[Q(0), Q(0), Q(0), Q(5, 4), Q(-5, 4)],
                      [Q(0), Q(0), Q(5, 4), Q(-5, 12), Q(-5, 6)],
                      [Q(0), Q(0), Q(5, 8), Q(-5, 8), Q(0)]]
    metric = [Q(1), Q(9, 8), Q(3, 2)]
    for matrix in [*transports.values(), restrict(generators[2]), restrict(generators[3])]:
        assert all(metric[i]*matrix[i][j] == metric[j]*matrix[j][i] for i in range(3) for j in range(3))
    for index in (2, 3):
        swapped = [row[:] for row in linear]
        for row in swapped:
            row[index], row[index+1] = row[index+1], row[index]
        assert product(restrict(generators[index]), linear) == swapped
    rays = [[Q(5-m)]*m+[Q(-m)]*(5-m) for m in range(2, 5)]
    images = [apply(linear, ray) for ray in rays]
    assert all(v >= 0 for image in images for v in image)
    assert all(sum(image[row] for image in images) > 0 for row in range(3))

    nodes = [(0, *rest) for rest in itertools.permutations(range(1, 5))]
    node_index = {g: i for i, g in enumerate(nodes)}
    coefficients = [[[Q(0) for _ in nodes] for _ in rays] for _ in range(3)]
    for rest in itertools.permutations(range(2, 5)):
        pi = (0, 1, *rest)
        inverse = [pi.index(i) for i in range(5)]
        for l in range(2, 5):
            g = tuple(l if z == 1 else 1 if z == l else z for z in pi)
            for m, ray in enumerate(rays):
                z = [ray[i] for i in inverse]
                field = apply(linear, z)
                moved = apply(transports[l], field)
                for row in range(3):
                    coefficients[row][m][node_index[g]] += field[row]-moved[row]+2*linear[row][l]*(z[1]-z[l])

    # Formal polynomials in twelve Jacobi parameters and four positive diagonal entries.
    variables, degree_zero = 16, (0,)*16
    def variable(index):
        exponent = list(degree_zero)
        exponent[index] = 1
        return {tuple(exponent): Q(1)}

    def plus(a, b):
        result = a.copy()
        for exponent, value in b.items():
            result[exponent] = result.get(exponent, Q(0))+value
        return {e: c for e, c in result.items() if c}

    def multiply(a, b):
        result = {}
        for first, c in a.items():
            for second, d in b.items():
                exponent = tuple(u+v for u, v in zip(first, second))
                result[exponent] = result.get(exponent, Q(0))+c*d
        return {e: c for e, c in result.items() if c}

    def polynomial_identity():
        return [[{degree_zero: Q(1)} if i == j else {} for j in range(4)] for i in range(4)]

    lower, upper = polynomial_identity(), polynomial_identity()
    word = (1, 2, 1, 3, 2, 1)
    for index, i in enumerate(word):
        for row in range(4):
            lower[row][i-1] = plus(lower[row][i-1], multiply(lower[row][i], variable(index)))
    for index, i in enumerate(reversed(word), 6):
        for row in range(4):
            upper[row][i] = plus(upper[row][i], multiply(upper[row][i-1], variable(index)))
    matrix = [[{} for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(4):
            for r in range(4):
                matrix[i][j] = plus(matrix[i][j], multiply(multiply(lower[i][r], variable(12+r)), upper[r][j]))
    monomials = []
    for g in nodes:
        poly = {degree_zero: Q(1)}
        for i in range(1, 5):
            poly = multiply(poly, matrix[i-1][g.index(i)-1])
        monomials.append(poly)
    expected_counts = [[164, 158, 187], [290, 224, 217], [126, 78, 61]]
    certificates = []
    for component in range(3):
        for m in range(3):
            poly = {}
            for c, term in zip(coefficients[component][m], monomials):
                if c:
                    poly = plus(poly, {e: c*v for e, v in term.items()})
            assert len(poly) == expected_counts[component][m]
            assert poly and all(value > 0 for value in poly.values())
            serialized = json.dumps([[list(e), str(c)] for e, c in sorted(poly.items())], separators=(",", ":"))
            certificates.append({"component": component+1, "ray": m+2,
                                 "positive_monomials": len(poly),
                                 "minimum_coefficient": str(min(poly.values())),
                                 "maximum_total_degree": max(sum(e) for e in poly),
                                 "coefficient_sha256": hashlib.sha256(serialized.encode()).hexdigest()})
    # Positive entries without total positivity do not suffice for the residual.
    bad = [[3, 5, 5, 3], [2, 6, 6, 5], [4, 4, 4, 1], [6, 3, 6, 6]]
    weights = [math.prod(bad[i-1][g.index(i)-1] for i in range(1, 5)) for g in nodes]
    assert bad[0][0]*bad[2][1]-bad[0][1]*bad[2][0] == -8
    assert sum(c*w for c, w in zip(coefficients[1][0], weights)) == -22825
    return {"k": 5, "partition": list(shape), "minus_dimension": 3,
            "Young_Coxeter_and_equivariance": "passed",
            "linear_field_and_metric": {"L": [[str(v) for v in row] for row in linear],
                                        "metric_diagonal": [str(v) for v in metric]},
            "formal_variables": variables, "lower_Jacobi_word": list(word),
            "upper_Jacobi_word": list(reversed(word)), "polynomials": certificates,
            "total_positive_coefficients": sum(r["positive_monomials"] for r in certificates),
            "positive_entries_without_TP_control": "passed",
            "scope": "Exact exhaustive expansion of nine fixed rational residual immanants under the 16-parameter positive LDU/Jacobi substitution, with all 1505 retained coefficients strictly positive. This is a formal polynomial identity certificate, not sampling. The written TP-factorization covering argument extends the positivity to every strictly TP 4x4 matrix. The analytic supersolution, Hilbert-space gap and ambient Hermite lift remain written proofs; other cell counts and the standard inequality are not certified."}


def standard_initial_and_antipodal_controls():
    """Exact bounded algebra for S77-S78; no covariance or analytic-limit proof."""
    dot = lambda x, y: sum(a*b for a, b in zip(x, y))
    initial = []
    for k in range(3, 13):
        p0 = [Q(1, 2)-Q(1, k)]*2+[Q(-1, k)]*(k-2)
        vectors = [[Q(1, 2)]*2+[Q(-int(i == l)) for i in range(2, k)]
                   for l in range(2, k)]
        assert dot(p0, p0) == Q(k-2, 2*k)
        assert [sum(w[i] for w in vectors) for i in range(k)] == [k*a for a in p0]
        assert all(dot(w, p0) == Q(1, 2) for w in vectors)
        transported = [k*p0[i]-sum(w[i]*dot(w, p0) for w in vectors) for i in range(k)]
        assert transported == [Q(k, 2)*a for a in p0]
        for j in range(2, k-1):
            value = [Q(0)]*k
            value[j], value[j+1] = Q(1), Q(-1)
            image = [k*value[i]-sum(w[i]*dot(w, value) for w in vectors) for i in range(k)]
            assert image == [(k-1)*a for a in value]
        initial.append({"k": k, "moment_direction_eigenvalue": str(Q(k, 2)),
                        "two_row_direction_eigenvalue": k-1})
    obstruction = []
    for k in range(7, 21):
        rays = {m: [Q(k-m)]*m+[Q(-m)]*(k-m) for m in range(2, k)}
        total = [sum(ray[i] for ray in rays.values()) for i in range(k)]
        x, y = rays[3], rays[k-1]
        swapped = list(y)
        swapped[2], swapped[k-1] = swapped[k-1], swapped[2]
        assert dot(x, y) == 3*k
        assert dot(x, swapped) == k*(3-k)
        epsilon = Q(1, k**4)
        interior_x = [a+epsilon*r for a, r in zip(x, total)]
        interior_y = [a+epsilon*r for a, r in zip(y, total)]
        assert sum(interior_x) == sum(interior_y) == 0
        assert interior_x[0] == interior_x[1] and interior_y[0] == interior_y[1]
        assert all(interior_x[i] > interior_x[i+1] and interior_y[i] > interior_y[i+1]
                   for i in range(1, k-1))
        swapped = list(interior_y)
        swapped[2], swapped[k-1] = swapped[k-1], swapped[2]
        gap = -dot(interior_x, swapped)-dot(interior_x, interior_y)
        assert gap > 0 and abs(gap-k*(k-6)) < 2
        assert 2*math.factorial(k-2)+(k-2)*math.factorial(k-2) == k*math.factorial(k-2)
        obstruction.append({"k": k, "epsilon": str(epsilon), "interior_exponent_gap": str(gap)})
    k = 6
    x, y = [3]*3+[-3]*3, [1]*5+[-5]
    swapped = list(y)
    swapped[2], swapped[5] = swapped[5], swapped[2]
    assert -dot(x, swapped) == dot(x, y) == 18
    return {"initial_rank_one_transport": initial, "antipodal_obstruction": obstruction,
            "six_cell_nonseparation_control": "passed",
            "scope": "Exact standard first-derivative transport algebra for 3<=k<=12 and ray/interior antipodal exponent separation for 7<=k<=20. The k=6 ray equality controls only this construction. These bounded calculations do not prove the convex-Gaussian covariance bound, Hilbert-Schmidt limit, unrestricted all-k exponential-scale obstruction, or standard stability at every correlation."}


def height_standard_reduction_controls():
    """Finite exact conditioning and coefficient controls for the written S84."""
    def product(a, b):
        return [[sum(a[i][r]*b[r][j] for r in range(2)) for j in range(2)] for i in range(2)]

    def determinant(a):
        return a[0][0]*a[1][1]-a[0][1]*a[1][0]

    def inverse(a):
        delta = determinant(a)
        return [[a[1][1]/delta, -a[0][1]/delta], [-a[1][0]/delta, a[0][0]/delta]]

    cases = []
    for rho in (Q(0), Q(1, 4), Q(1, 2), Q(3, 4)):
        covariance = [[Q(1), rho], [rho, Q(1)]]
        for adjacent in (False, True):
            ties = [[Q(2), rho if adjacent else 2*rho],
                    [rho if adjacent else 2*rho, Q(2)]]
            correction = product(product(covariance, inverse(ties)), covariance)
            conditional = [[covariance[i][j]-correction[i][j] for j in range(2)] for i in range(2)]
            if adjacent:
                expected = [[(2-rho*rho)/(4-rho*rho), rho/(4-rho*rho)],
                            [rho/(4-rho*rho), (2-rho*rho)/(4-rho*rho)]]
                assert determinant(conditional) == (1-rho*rho)/(4-rho*rho)
                assert inverse(conditional) == [[(2-rho*rho)/(1-rho*rho), -rho/(1-rho*rho)],
                                                [-rho/(1-rho*rho), (2-rho*rho)/(1-rho*rho)]]
                # Squared Gaussian density prefactors after factoring (2*pi)^-2.
                # The surface product is2 and the normalized tie covariance is ties/2.
                assert Q(4)/(determinant(ties)*determinant(conditional)) == Q(4)/(1-rho*rho)
            else:
                expected = [[v/2 for v in row] for row in covariance]
                assert determinant(conditional) == (1-rho*rho)/4
                assert Q(4)/(determinant(ties)*determinant(conditional)) == Q(4)/(1-rho*rho)**2
            assert conditional == expected
            cases.append({"rho": str(rho), "adjacent": adjacent,
                          "conditional_top_covariance": [[str(v) for v in row] for row in conditional]})

    coefficient_cases = []
    for k in range(3, 15):
        q = k-2
        joint = (4, 2*q)  # p^2*F^q and p*Fu*Fv*F^(q-1).
        differentiated = (2*(k-1), (k-1)*(k-2))
        assert tuple(Q(2, k-1)*v for v in differentiated) == joint
        assert 4+2*q == 2*k  # J_(k,0)=k*m_k*m_k.
        assert 4 != 2*k  # Dropping adjacent facets corrupts that normalization.
        coefficient_cases.append({"k": k, "joint_coefficients": list(joint),
                                  "noise_derivative_coefficients": list(differentiated)})
    # Treating the two adjacent tie normals as independent is incorrect.
    rho = Q(1, 2)
    covariance = [[Q(1), rho], [rho, Q(1)]]
    correction = product(product(covariance, [[Q(1, 2), Q(0)], [Q(0), Q(1, 2)]]), covariance)
    wrong = [[covariance[i][j]-correction[i][j] for j in range(2)] for i in range(2)]
    assert wrong != [[Q(7, 15), Q(2, 15)], [Q(2, 15), Q(7, 15)]]
    return {"conditional_gaussian_cases": cases, "coefficient_cases": coefficient_cases,
            "independent_adjacent_normals_rejection": "passed",
            "same_facet_only_rejection": "passed",
            "scope": "Exact conditional-top covariance, Gaussian-prefactor determinant and formal coefficient controls at four rational correlations and3<=k<=14. The all-k tied-height surface integral, ordinary mass constraint, compact operator compression, stability derivative identity and linear trace-rotation identification in S84 are written proofs. No sign of the restricted or full standard form is certified."}


def convex_nonlocal_mass_counterexample():
    """Exact directed rational margin for the actual Mehler counterexample S79."""
    # Algebra in Q(z)[d], where z represents sqrt(3), after multiplication by pi.
    # At rho=1/2, pi*alpha^2=81/32, pi*(alpha-mu0)^2=1/32,
    # and pi*(rho*D-C)=-9*z/16.
    for dimension in range(1, 21):
        constant = dimension*Q(81, 64)-dimension**2*Q(1, 64)
        root_coefficient = -dimension*Q(9, 16)
        assert (constant, root_coefficient) == (Q(dimension*(81-dimension), 64), -Q(36*dimension, 64))
    rho, dimension, root_lo, pi_hi = Q(1, 2), 19, Q(433, 250), Q(22, 7)
    assert root_lo**2 < 3
    assert 62**2 < 36**2*3 < 63**2
    assert 20**2 > 9**2*3
    numerator_hi = dimension*(81-dimension-36*root_lo)
    assert numerator_hi < 0
    upper = numerator_hi/(64*pi_hi)
    assert upper == -Q(133, 4000) < -Q(1, 32)
    return {"rho": str(rho), "dimension": dimension,
            "normalized_difference": "19*(62-36*sqrt(3))/(64*pi)",
            "normalized_difference_upper": str(upper), "negative_margin": "greater than 1/32",
            "dimension_polynomial": "d*(81-d-36*sqrt(3))/(64*pi)",
            "ordinary_mean_zero_linear_function": "sum_i x_i - d*sqrt(2/pi)",
            "first_negative_integer_in_this_linear_family": 19,
            "one_and_eighteen_dimensional_positive_controls": "passed",
            "scope": "Exact quadratic-in-dimension moment algebra and directed rational square/interval arithmetic for the linear counterexample on the positive orthant in dimension19 at rho1/2. Dimensions1 and18 are positive within this linear family; no lower-dimensional classification is implied. The Gaussian pair moment identities, ordinary-mean-zero realization and all-correlation finite-dimension family require the written analytic proof. This does not refute the simplex standard inequality S51."}


def standard_mass_weight_control():
    """Reject removing the translation eigenvector using the wrong mass weight."""
    rho, weights, a, h = Q(1, 2), [Q(1, 2)]*2, [Q(1), Q(4)], [Q(1), Q(-1)]
    kernel = [[Q(2), Q(2)], [Q(2), Q(14)]]
    root_a, b = [Q(1), Q(2)], [Q(1), Q(1, 2)]
    dot = lambda x, y: sum(w*xi*yi for w, xi, yi in zip(weights, x, y))
    apply = lambda matrix, x: [sum(v*y for v, y in zip(row, x)) for row in matrix]
    for i in range(2):
        assert sum(weights[j]*kernel[i][j] for j in range(2)) == a[i]/rho
    normalized = [[weights[j]*kernel[i][j]/(root_a[i]*root_a[j]) for j in range(2)]
                  for i in range(2)]
    assert apply(normalized, root_a) == [x/rho for x in root_a]
    trace = normalized[0][0]+normalized[1][1]
    determinant = normalized[0][0]*normalized[1][1]-normalized[0][1]*normalized[1][0]
    assert trace == Q(2)+Q(3, 4) and determinant == Q(2)*Q(3, 4)
    u = [root_a[i]*h[i] for i in range(2)]
    assert dot(h, [Q(1)]*2) == dot(u, b) == 0
    assert dot(u, root_a) != 0
    mass_term = sum(weights[i]*a[i]*h[i]**2 for i in range(2))
    interaction = sum(weights[i]*weights[j]*kernel[i][j]*h[i]*h[j]
                      for i in range(2) for j in range(2))
    assert mass_term-interaction == Q(-1, 2)
    inverse_form = [[Q(3), Q(-2)], [Q(-2), Q(0)]]
    form = [[Q(i == j)-normalized[i][j] for j in range(2)] for i in range(2)]
    for i in range(2):
        for j in range(2):
            assert sum(form[i][r]*inverse_form[r][j] for r in range(2)) == (i == j)
    assert dot(b, apply(inverse_form, b)) == Q(1, 2)
    return {"row_identity": "passed", "normalized_eigenvalues": ["2", "3/4"],
            "balanced_unnormalized_form": "-1/2", "mass_resolvent_form": "1/2",
            "scope": "Exact two-atom counterexample to using a second normalized eigenvalue bound alone with ordinary mass balance; it is not a simplex counterexample and does not certify either all-k test in S70."}


def iv(x):
    return (Q(x[0]), Q(x[1])) if isinstance(x, tuple) else (Q(x), Q(x))


def add(x, y):
    x, y = iv(x), iv(y)
    return x[0] + y[0], x[1] + y[1]


def neg(x):
    x = iv(x)
    return -x[1], -x[0]


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    values = [a * b for a in iv(x) for b in iv(y)]
    return min(values), max(values)


def div(x, y):
    y = iv(y)
    if y[0] * y[1] <= 0:
        raise ValueError("interval denominator must exclude zero")
    return mul(x, (1 / y[1], 1 / y[0]))


def round_out(x, scale=10**10):
    lo, hi = iv(x)
    return Q((lo * scale).__floor__(), scale), Q((hi * scale).__ceil__(), scale)


def sqrt_point(x, scale=10**12):
    x = Q(x)
    lo = Q(math.isqrt(x.numerator * scale**2 // x.denominator), scale)
    hi = lo + Q(1, scale)
    assert lo * lo <= x < hi * hi
    return lo, hi


def sqrt_iv(x):
    lo, hi = iv(x)
    return sqrt_point(lo)[0], sqrt_point(hi)[1]


def atan_point(x, count=20):
    x = Q(x)
    assert 0 < x < 1
    value = sum(((-1)**n * x**(2*n+1) / (2*n+1) for n in range(count)), Q(0))
    next_term = (-1)**count * x**(2*count+1) / (2*count+1)
    return min(value, value + next_term), max(value, value + next_term)


def pi_interval():
    # tan(4 atan(1/5)-atan(1/239))=1, with the angle in (0,pi/2).
    # Thus Machin's identity gives pi=16 atan(1/5)-4 atan(1/239).
    return round_out(sub(mul(16, atan_point(Q(1, 5))), mul(4, atan_point(Q(1, 239)))))


def asin_point(x, count=64):
    x = Q(x)
    if x < 0:
        return neg(asin_point(-x, count))
    assert 0 <= x < 1
    def coefficient(n):
        return Q(math.comb(2*n, n), 4**n * (2*n+1))
    lower = sum((coefficient(n) * x**(2*n+1) for n in range(count)), Q(0))
    # Positive decreasing coefficients bound the omitted tail geometrically.
    tail = coefficient(count) * x**(2*count+1) / (1 - x*x)
    return round_out((lower, lower + tail))


def asin_iv(x):
    lo, hi = iv(x)
    return asin_point(lo)[0], asin_point(hi)[1]


def operator_bounds(lower, upper):
    """Rational interval extensions of the analytic functions H and G."""
    rho, pi = (Q(lower), Q(upper)), pi_interval()
    rho2 = mul(rho, rho)
    sigma2 = sub(1, rho2)
    def quadrant(x):
        return add(Q(1, 4), div(asin_iv(x), mul(2, pi)))
    def J(u):
        base = add(1, rho2)
        delta = sub(mul(base, base), mul(rho2, u*u))
        return div(quadrant(div(mul(rho, u), base)), sqrt_iv(delta))
    j2, jhalf, jneg = J(Q(2)), J(Q(1, 2)), J(Q(-1))
    H = mul(div(Q(16, 9), rho2), add(sub(j2, mul(2, jhalf)), jneg))
    R = add(add(mul(4, j2), mul(4, jhalf)), jneg)
    L = add(div(mul(2, quadrant(rho)), sqrt_iv(sigma2)),
            div(quadrant(mul(Q(-1, 2), rho)), sqrt_iv(sub(1, mul(Q(1, 4), rho2)))))
    G = sub(mul(div(Q(4, 9), rho2), R),
            mul(div(mul(Q(16, 9), sigma2), rho2), mul(L, L)))
    return H, G


def interval_gap_check():
    # Check interval primitives on exact cases, including their rejection path.
    assert mul((-2, -1), (3, 4)) == (-8, -3)
    assert div((2, 4), (-2, -1)) == (-4, -1)
    try:
        div(1, (-1, 1))
    except ValueError:
        pass
    else:
        raise AssertionError("invalid denominator accepted")
    half_H, half_G = operator_bounds(Q(1, 2), Q(1, 2))
    assert half_H[1] < Q(7, 4)
    assert half_G[1] < Q(9, 16)
    maxima = [Q(0), Q(0)]
    for j in range(250, 600):
        H, G = operator_bounds(Q(j, 1000), Q(j+1, 1000))
        assert H[1] < Q(19, 10), (j, "H")
        assert G[1] < Q(4, 5), (j, "G")
        maxima = [max(maxima[0], H[1]), max(maxima[1], G[1])]
    rounded = [round_out(x)[1] for x in maxima]
    return {
        "interval": ["1/4", "3/5"],
        "closed_cover_mesh": "1/1000",
        "closed_cover_intervals": 350,
        "H_upper": str(rounded[0]),
        "G_upper": str(rounded[1]),
        "targets": {"H_less_than": "19/10", "G_less_than": "4/5"},
        "half_correlation_G_upper": str(round_out(half_G)[1]),
        "arithmetic": "Python Fraction, integer-square-root enclosures; positive arcsin series with geometric tail, alternating atan series and Machin pi identity; outward rational rounding.",
        "scope": "Uniform H/G bounds on every point of the stated interval, via the analytic operator reduction S24-S26. No sampling argument.",
        "full_three_cell_gap": "ell*(1-3/sqrt(10)) on [1/4,3/5]",
        "full_three_cell_gap_at_half": "ell*(2/sqrt(3)-1)",
    }


def quadrature():
    """Test S23 at four chosen correlations, without error certification."""
    import numpy as np
    phi0 = 1 / math.sqrt(2 * math.pi)
    cdf = np.vectorize(lambda x: (1 + math.erf(x / math.sqrt(2))) / 2)
    output = []
    for count in (180, 360):
        grid, weights = np.polynomial.legendre.leggauss(count)
        r, dr = 5 * (grid + 1), 5 * weights  # Finite [0,10] cutoff.
        w = dr * np.exp(-r * r / 2) / (2 * math.pi)
        rr, ss = r[:, None], r[None, :]
        for rho in (0.1, 0.5, 0.9, 0.99):
            sigma = math.sqrt(1 - rho * rho)
            u = rho * r / sigma
            a = rho / sigma * (2 * phi0 * cdf(u) + phi0 * np.exp(-3 * u * u / 8) * cdf(-u / 2))
            # Include quadrature weights before exponentiation to avoid overflow.
            base = (np.log(w)[:, None] + np.log(w)[None, :]) / 2
            base -= rho * rho * (rr * rr + ss * ss) / (2 * sigma * sigma)
            same = np.exp(base + rho * rr * ss / (sigma * sigma)) / (sigma * sigma)
            other = np.exp(base - rho * rr * ss / (2 * sigma * sigma)) / (sigma * sigma)
            scaled = (2 * same + other) / np.sqrt(a[:, None] * a[None, :])
            constraint = np.sqrt(w / a)
            # Euclidean coordinates are sqrt(a*w)*h, so the constraint is nu-mean.
            basis, _ = np.linalg.qr(constraint[:, None], mode="complete")
            reduced = basis[:, 1:].T @ scaled @ basis[:, 1:]
            maximum = float(np.linalg.eigvalsh(reduced)[-1])
            standard_basis, _ = np.linalg.qr(np.sqrt(w)[:, None], mode="complete")
            sign_basis, _ = np.linalg.qr((np.sqrt(w)*r)[:, None], mode="complete")
            def restricted_min(matrix, basis):
                return float(np.linalg.eigvalsh(basis[:, 1:].T @ matrix @ basis[:, 1:]) [0])
            standard_min = restricted_min(np.diag(a/rho)-(2*same+other)/rho, standard_basis)
            sign_min = restricted_min(np.diag(a/rho)-2*(same-other)/rho, sign_basis)
            mixed_rotation_min = restricted_min(np.diag(a/rho)-(2*same+other), standard_basis)
            ell = 1.5*phi0
            gap = ell*(1/sigma-math.sqrt(3)/2) if rho <= .5 else ell*(1-rho)*(3*rho-1)/(3*rho*sigma)
            assert standard_min >= gap-1e-8
            assert sign_min >= ell*sigma-1e-8
            assert mixed_rotation_min >= ell*sigma-1e-8
            output.append({"nodes": count, "cutoff": 10, "rho": rho,
                           "largest_constrained_generalized_eigenvalue": maximum,
                           "standard_L2_minimum": standard_min, "sign_L2_minimum_off_rotations": sign_min,
                           "extra_degree_one_standard_minimum_off_rotations": mixed_rotation_min,
                           "analytic_full_gap_lower": gap, "analytic_sign_gap_lower": ell*sigma})
    return {"numpy_version": np.__version__, "scope": "Uncertified Nystrom exploration on [0,10]; no tail, discretization or rounding error enclosure; not an all-rho or fixed-rho proof.", "cases": output}


def perimeter_controls():
    import numpy as np
    q = NormalDist().inv_cdf(2/3)
    beta = 3/(4*math.pi**2)-1/(4*math.pi*math.sqrt(3))
    c = (2*math.exp(-q*q/2)-1.5)/math.pi
    x, w = np.polynomial.legendre.leggauss(40)
    cases = []
    for theta in (.1, .05, .025):
        u = theta*(x+1)/2
        Cq = float(theta/2*np.sum(w*np.exp(-q*q/(1+np.cos(u))))/math.pi)
        NT = 2/3-3/(4*math.pi**2)*((math.pi-theta)**2-math.acos(math.cos(theta)/2)**2)
        J_upper = math.exp(-q*q/(1-math.cos(theta)))
        normalized_without_jump = (2*Cq-NT)/theta
        cases.append({"theta": theta, "simplex_beta_quotient": (3*theta/(2*math.pi)-NT)/theta**2,
                      "slab_normalized_deficit_without_jump": normalized_without_jump,
                      "slab_normalized_limit": c,
                      "jump_correction_upper": 2*J_upper/theta})
        assert normalized_without_jump-2*J_upper/theta > c
    for rho in (.1, .5, .9):
        def I(t):
            return (.25+math.asin(rho*t)/(2*math.pi))/(2*math.pi*math.sqrt(1-rho*rho*t*t))
        derivative = 3*math.acos(-rho)/(2*math.pi**2*math.sqrt(1-rho*rho))
        derivative += 3*math.acos(rho/2)/(4*math.pi**2*math.sqrt(1-rho*rho/4))
        assert abs(derivative-(6*I(1)+3*I(-.5))) < 1e-12
    return {"scope": "Uncertified floating normalization/expansion controls, not a proof of D0-D9 or monotonicity failure.",
            "q": q, "beta3": beta, "cases": cases}


def main():
    if not __debug__:
        raise RuntimeError("Run without -O: the exact certificate uses assertions.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quadrature", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    receipt = {
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "proof_source_sha256": hashlib.sha256((here / "SIMPLEX_INTERFACE_SYMMETRY.md").read_bytes()).hexdigest(),
        "perimeter_proof_source_sha256": hashlib.sha256((here / "THREE_CELL_DEFICIT_AND_PERIMETER.md").read_bytes()).hexdigest(),
        "exact_checks": {"characters": character_checks(), "half_correlation_bound": explicit_gap_check(),
                         "elementary_half_point_and_reflections": elementary_half_point_and_reflections(),
                         "higher_cell_reflections": higher_cell_reflection_checks(),
                         "hook_cone": hook_cone_checks(),
                         "sign_vandermonde": sign_vandermonde_controls(),
                         "universal_blocks": universal_block_controls(),
                         "triple_sign_reflections": triple_sign_reflection_controls(),
                         "antisymmetry_height": antisymmetry_height_controls(),
                         "standard_mass_weight": standard_mass_weight_control(),
                         "standard_hook_gaussian": standard_hook_gaussian_controls(),
                         "higher_exterior_hooks": higher_exterior_hook_checks(),
                         "representation_hermite_floors": representation_hermite_floor_checks(),
                         "height_standard_reduction": height_standard_reduction_controls(),
                         "five_cell_mixed_certificate": five_cell_mixed_certificate(),
                         "standard_initial_and_antipodal": standard_initial_and_antipodal_controls(),
                         "convex_nonlocal_mass_counterexample": convex_nonlocal_mass_counterexample(),
                         "interval_gap": interval_gap_check()},
        "analytic_claims_certified_by_this_script": False,
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
    }
    if args.quadrature:
        receipt["standard_sector_exploration"] = quadrature()
        receipt["perimeter_controls"] = perimeter_controls()
    rendered = json.dumps(receipt, indent=2) + "\n"
    if args.output:
        args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
