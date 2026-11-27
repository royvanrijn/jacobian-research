"""Exact corner, successor, and planar coefficient checks for Weyl degree 56."""
from collections import defaultdict
from fractions import Fraction
from math import lcm
from operator import add

from verify_weyl_bernstein_through_fifty_five import corners, successors

F = Fraction


def partitions(total, ceiling, excluded):
    if total == 0:
        yield ()
    else:
        for part in range(min(total, ceiling), 0, -1):
            if part == excluded:
                continue
            for tail in partitions(total - part, part, excluded):
                yield (part, *tail)


def candidate_scalars(denominator, multiplicity, deficit, previous):
    candidates = []
    for n in range(5, 2 * deficit + 1, 2):
        scaled_slope = denominator - F(deficit * n, n * multiplicity - 1)
        index = scaled_slope.denominator
        if index > 2 * multiplicity or (2 % index and (n - 2) % index):
            continue
        slope = scaled_slope / denominator
        if slope > previous:
            candidates.append((n - 2, slope))
    return candidates


def candidate_multiples(denominator, multiplicity, deficit, previous):
    _, lattice = successors(denominator, multiplicity, deficit, previous)
    candidates = set()
    for degree, orbit, zero, slope in lattice:
        if zero * degree == multiplicity:
            continue
        for parts in partitions((multiplicity - zero) // orbit,
                                multiplicity, F(multiplicity, degree)):
            num_roots = orbit * len(parts) + bool(zero)
            if not 2 <= num_roots <= degree:
                continue
            orders = (zero, *parts)
            if all(order * degree % multiplicity == 0 for order in orders):
                continue
            root = max(orders)
            new_denominator = denominator if root == zero else lcm(
                denominator, slope.denominator)
            corner_x = multiplicity - F(deficit, denominator)
            new_x = corner_x - slope * (multiplicity - root)
            new_deficit = new_denominator * (root - new_x)
            assert new_deficit.denominator == 1
            candidates.add((degree, orbit, slope, root,
                            new_denominator, int(new_deficit)))
    return candidates


def check_corners():
    rows = [
        ((7, 21), (3, 1), (2, 6), (5, -1)),
        ((7, 21), (2, 1), (3, 9), (4, -1)),
        ((7, 21), (1, 0), (5, 15), (7, -2)),
        ((8, 20), (3, 0), (2, 5), (4, -1)),
        ((10, 18), (1, 0), (5, 9), (2, -1)),
    ]
    assert corners(16) == [
        ((4, 12), (1, 0), (3, 9), (4, -1)),
        ((6, 10), (1, 0), (3, 5), (2, -1)),
    ]
    assert corners(18) == [((6, 12), (2, 0), (2, 4), (3, -1))]
    assert corners(21) == [((6, 15), (1, 0), (4, 10), (3, -1))]
    assert corners(27) == [((9, 18), (3, 0), (2, 4), (3, -1))]
    assert corners(28) == rows
    assert all(not corners(n) for n in (17, 19, 29))
    scalars = []
    multiples = []
    multiplicities = []
    for (u, v), (h, j), _, (rho, sigma) in rows:
        edge_degree = F(v - j, rho)
        assert edge_degree.denominator == 1
        weight = rho * h + sigma * j
        bracket_degree = ((rho + sigma) * edge_degree - (h - j)) / weight
        assert bracket_degree.denominator == 1
        roots = sorted({max(j, *part)
                        for part in partitions(int(edge_degree), int(edge_degree),
                                               F(weight, rho + sigma))
                        if len(part) <= bracket_degree})
        multiplicities.append(roots)
        for mu in roots:
            beta = (rho + sigma) * mu - weight
            data = (rho, mu, beta, F(-sigma, rho))
            scalars.extend((u, v, h, j, mu, *case)
                           for case in candidate_scalars(*data))
            for case in sorted(candidate_multiples(*data)):
                multiples.append((data, case))
                _, _, slope, root, new_denominator, new_deficit = case
                next_data = (new_denominator, root, new_deficit, slope)
                assert not candidate_multiples(*next_data)
                assert not candidate_scalars(*next_data)
    assert multiplicities == [[4], [3, 4, 5], [2, 3], [5], list(range(3, 10))]
    assert multiples == [
        ((5, 4, 2, F(1, 5)), (2, 1, F(4, 5), 3, 5, 1)),
        ((2, 6, 4, F(1, 2)), (9, 4, F(5, 8), 2, 2, 1)),
    ]
    assert scalars == [
        (7, 21, 2, 1, 4, 7, F(19, 28)), (7, 21, 2, 1, 5, 3, F(7, 12)),
        (7, 21, 1, 0, 2, 3, F(16, 21)), (7, 21, 1, 0, 3, 5, F(3, 5)),
    ]
    print('Five sum-28 crossings, two commuting outputs, four scalar endpoints.')


def add_bracket(left, right, output):
    for (i, j), a in left.items():
        for (k, m), b in right.items():
            determinant = i * m - j * k
            if determinant:
                output[i+k-1, j+m-1] += determinant * a * b


def build_layers(ring, name, degree, weight, left_weight, right_weight,
                 top, fixed):
    known = defaultdict(dict)
    unknown = defaultdict(list)
    for i in range(degree + 1):
        for j in range(degree + 1):
            total_degree = left_weight(i, j)
            layer = weight - right_weight(i, j)
            if layer < 0 or total_degree > degree:
                continue
            key = (i, j)
            coefficient = (top.get(key, 0) if total_degree == degree
                           else fixed.get((name, *key)))
            if coefficient is None:
                unknown[layer].append(key)
            elif coefficient:
                known[layer][key] = ring(coefficient)
    return known, unknown


def solve_layer(ring, left, right, previous, columns, free_parameters):
    exponents = sorted(set(previous).union(*(col.keys() for col in columns)))
    reduced = [[F(str(col[exp])) for col in columns] + [-previous[exp]]
               for exp in exponents]
    pivots = []
    for index in range(len(columns)):
        target = len(pivots)
        pivot = next((i for i in range(target, len(reduced))
                      if reduced[i][index]), None)
        if pivot is None:
            continue
        reduced[target], reduced[pivot] = reduced[pivot], reduced[target]
        scale = reduced[target][index]
        reduced[target] = [value / scale for value in reduced[target]]
        for i, row in enumerate(reduced):
            if i != target and row[index]:
                reduced[i] = [a - row[index] * b
                              for a, b in zip(row, reduced[target])]
        pivots.append(index)
    free_values = {k: free_parameters[key]
                   for k, key in enumerate(left + right) if k not in pivots}
    values = free_values.copy()
    for k, row in zip(pivots, reduced):
        values[k] = row[-1] - sum(
            (ring(row[j]) * value for j, value in free_values.items()), ring.zero)
    return values, [row[-1] for row in reduced[len(pivots):] if row[-1]]


def coefficient_system(ring, left_weight, right_weight, tops, caps,
                       fixed, rhs, free, last_layer, rhs_unknown=None):
    p, unknown_p = build_layers(ring, 'p', *caps[0], left_weight,
                                 right_weight, tops[0], fixed)
    q, unknown_q = build_layers(ring, 'q', *caps[1], left_weight,
                                 right_weight, tops[1], fixed)
    residuals = {}
    max_bracket_weight = caps[0][1] + caps[1][1] - right_weight(1, 1)
    for layer in range(last_layer + 1):
        previous = defaultdict(ring.zero.copy)
        for earlier in range(layer + 1):
            add_bracket(p[earlier], q[layer - earlier], previous)
        for key, coefficient in rhs.items():
            if right_weight(*key) == max_bracket_weight - layer:
                previous[key] -= ring(coefficient)
        # This order makes the free variable names stable across Python versions.
        keys_p = [('p', key) for key in sorted(unknown_p[layer], key=lambda e: f"{e[0]}_{e[1]}")]
        keys_q = [('q', key) for key in sorted(unknown_q[layer], key=lambda e: f"{e[0]}_{e[1]}")]
        columns = []
        for name, key in keys_p + keys_q:
            col = defaultdict(F)
            if name == 'p':
                add_bracket({key: 1}, q[0], col)
            else:
                add_bracket(p[0], {key: 1}, col)
            columns.append(col)
        if rhs_unknown and layer == rhs_unknown[0]:
            keys_q.append(('rhs', (0, 0)))
            columns.append(defaultdict(F, {(0, 0): rhs_unknown[1]}))
        values, residual = solve_layer(ring, keys_p, keys_q, previous, columns, free)
        for index, (name, key) in enumerate(keys_p + keys_q):
            if name != 'rhs' and values[index]:
                (p if name == 'p' else q)[layer][key] = values[index]
        if residual:
            residuals[layer] = residual
    return residuals


def check_five():
    from sympy.polys.domains import QQ
    from sympy.polys.rings import ring

    coefficients, *parameters = ring('A,B,C,D', QQ)
    fixed = {('p', 0, 0): 0, ('p', 1, 1): 1, ('q', 0, 0): 0,
             ('q', 2, 0): F(5, 16), ('q', 3, 5): F(15, 8), ('q', 4, 10): F(5, 2),
             ('q', 4, 12): 0, ('q', 2, 6): 0}
    free = dict(zip([('q', e) for e in [(4, 14), (4, 15), (3, 13), (3, 15)]],
                    parameters))
    residuals = coefficient_system(
        coefficients, lambda i, j: i + 2*j, lambda i, j: 5*i - j,
        ({(2, 6): 1}, {(5, 15): 1}), ((14, 4), (35, 10)),
        fixed, {(2, 0): -F(5, 8), (0, 1): F(5, 8)}, free, 11,
        rhs_unknown=(10, -F(5, 8)))
    assert set(residuals) == {11}
    assert len(residuals[11]) == 1
    assert residuals[11][0].is_ground and residuals[11][0]
    print('Ratio 2:5: a nonzero rational obstruction at layer 11.')


def seven_system(last_layer=77):
    from sympy.polys.domains import QQ
    from sympy.polys.rings import ring

    coefficients, *parameters = ring('A,B,C,D,E', QQ)
    _, s, v = ring('s,v', QQ)
    p_top = 16*s**8*(s+v)**4*(s-v/4)**2
    q_top = -16384*s**28*(s+v)**14*(s-v/4)**7
    fixed = {('p', 0, 0): 0, ('p', 1, 1): 1, ('q', 0, 0): 0,
             **{('q', i, i): 0 for i in (1, 2, 3)}}
    for k, value in enumerate((F(35, 128), F(35, 16), F(35, 8), F(7, 2), 1)):
        fixed['q', 7*k, 1+5*k] = value
    free = dict(zip([('q', e) for e in [(6, 5), (5, 4), (4, 3), (7, 5), (3, 2)]],
                    parameters))
    residuals = coefficient_system(
        coefficients, lambda i, j: i + j, lambda i, j: 7*j - 5*i,
        (p_top, q_top), ((14, 2), (49, 7)), fixed,
        {(0, 1): F(35, 128)}, free, last_layer)
    expected = (24, 29, 34, 36, 39, 41, 44, 46, 48, 49, 51, 53, 54, 56,
                58, 59, 60, 61, 63, 64, 65, 66, 68, 69, 70, 71, 72, 73,
                74, 75, 76, 77)
    assert tuple(residuals) == tuple(k for k in expected if k <= last_layer)
    assert all(len(row) == 1 for row in residuals.values())
    return coefficients, residuals


def weighted_monomials(weight, weights):
    if not weights:
        if weight == 0:
            yield ()
        return
    head, *tail = weights
    for exponent in range(weight // head + 1):
        for rest in weighted_monomials(weight - head * exponent, tail):
            yield (exponent, *rest)


def check_seven():
    """Find multipliers using modular pivots, then expand their exact identity."""
    from flint import fmpq, fmpq_mat, nmod_mat

    ring, rows = seven_system(56)
    weights = (2, 4, 6, 7, 8)
    certificate_weight = 72
    prime = 32003
    equations = [row[0] for row in rows.values()]
    for layer, equation in zip(rows, equations):
        for exponents in equation:
            degree = sum(a*b for a, b in zip(exponents, weights))
            assert layer >= degree and (layer - degree) % 12 == 0
    columns = [exponents for k in range(certificate_weight, -1, -12)
               for exponents in weighted_monomials(k, weights)]
    column_index = {exponents: index for index, exponents in enumerate(columns)}
    products = [(f, monomial)
                for f, layer in enumerate(rows)
                for k in range(certificate_weight - layer, -1, -12)
                for monomial in weighted_monomials(k, weights)]
    assert (len(products), len(columns)) == (1341, 1918)
    constant = column_index[(0, 0, 0, 0, 0)]

    def terms(index):
        equation, monomial = products[index]
        for exponents, coefficient in equations[equation].items():
            yield column_index[tuple(map(add, exponents, monomial))], coefficient

    def residue(coefficient):
        numerator = int(coefficient.numerator)
        denominator = int(coefficient.denominator)
        return numerator * pow(denominator, -1, prime) % prime

    modular = nmod_mat(len(columns), len(products)+1, prime)
    for k in range(len(products)):
        for j, coefficient in terms(k):
            modular[j, k] = residue(coefficient)
    modular[constant, len(products)] = 1
    echelon, rank = modular.rref()
    independent_products = [
        next(k for k in range(len(products)+1) if echelon[j, k])
        for j in range(rank)
    ]
    assert rank == 1258 and len(products) not in independent_products

    transpose = nmod_mat(rank, len(columns), prime)
    for k, index in enumerate(independent_products):
        for j, coefficient in terms(index):
            transpose[k, j] = residue(coefficient)
    echelon, transpose_rank = transpose.rref()
    assert transpose_rank == rank
    selected_columns = [next(k for k in range(len(columns)) if echelon[j, k])
                        for j in range(rank)]
    selected_index = {k: j for j, k in enumerate(selected_columns)}
    square = fmpq_mat(rank, rank)
    for k, index in enumerate(independent_products):
        for j, coefficient in terms(index):
            if j in selected_index:
                square[selected_index[j], k] = fmpq(
                    int(coefficient.numerator), int(coefficient.denominator))
    target = fmpq_mat([[int(k == constant)] for k in selected_columns])
    solution = square.solve(target, algorithm='dixon')
    multipliers = [ring.zero.copy() for _ in equations]
    for index, j in enumerate(independent_products):
        if solution[index, 0]:
            equation, monomial = products[j]
            multipliers[equation][monomial] = ring.domain(str(solution[index, 0]))
    assert sum((h*f for h, f in zip(multipliers, equations)), ring.zero) == 1
    print('Ratio 2:7: exact rational unit-ideal identity from layers through 56.')


if __name__ == '__main__':
    check_corners()
    check_five()
    check_seven()
