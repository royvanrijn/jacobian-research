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

    young_cases = []
    shapes = [shape for k in range(3, 6) for shape in partitions(k)]
    shapes += [(3, 2, 1), (2, 2, 2), (3, 1, 1, 1), (2, 2, 1, 1)]
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
        "orientation_and_non_TN_negative_controls": "passed",
        "scope": "Exact S68 coefficient identity for 3<=k<=7; exact rational seminormal Coxeter relations and nonnegative I+s_i matrices for every partition through k=5 and four selected k=6 partitions; selected nonsymmetric TP exponential-matrix sums and their (12)-minus corners. The all-k immanant formula and positivity depend on the written argument and retained primary Young-form/interval-decomposition statements; no spectral gap is certified.",
    }


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
                         "standard_mass_weight": standard_mass_weight_control(),
                         "standard_hook_gaussian": standard_hook_gaussian_controls(),
                         "higher_exterior_hooks": higher_exterior_hook_checks(),
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
