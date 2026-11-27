"""Exact rational supersolution certificate for the six-cell (3,2,1) block."""

import itertools
import json
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from math import lcm


def run():
    k, shape = 6, (3, 2, 1)
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
    dimension_minus = len(minus)
    assert dimension_minus == 8

    def restrict(matrix):
        return [[matrix[i][j] for j in minus] for i in minus]

    # The rational seminormal matrices are a positive diagonal rescaling
    # of the positive off-diagonal orthogonal Young basis.
    metric = {0: Q(1)}
    while len(metric) < dimension:
        before = len(metric)
        for matrix in generators:
            for i in list(metric):
                for j in range(dimension):
                    if i != j and matrix[i][j]:
                        candidate = metric[i] * matrix[i][j] / matrix[j][i]
                        assert candidate > 0
                        assert metric.setdefault(j, candidate) == candidate
        assert len(metric) > before
    for matrix in generators:
        assert all(metric[i]*matrix[i][j] == metric[j]*matrix[j][i]
                   for i in range(dimension) for j in range(dimension))

    seed = [Q(1)]*dimension_minus
    stabilizer_orbit = [apply(restrict(representation((0, 1, 2, *rest))), seed)
                        for rest in itertools.permutations(range(3,k))]
    seed_sum = [sum(row) for row in zip(*stabilizer_orbit)]
    linear = [[Q(0)]*k for _ in range(dimension_minus)]
    transports = {}
    for label in range(2,k):
        permutation = list(range(k))
        permutation[2], permutation[label] = permutation[label], permutation[2]
        vector = apply(restrict(representation(tuple(permutation))), seed_sum)
        for row,value in enumerate(vector):
            linear[row][label] = value*Q(16,169)
        permutation = list(range(k))
        permutation[1], permutation[label] = permutation[label], permutation[1]
        transports[label] = [[2*v for v in row]
                             for row in restrict(representation(tuple(permutation)))]
    assert all(not row[0] and not row[1] and not sum(row) for row in linear)
    for index in range(2,k-1):
        swapped = [row.copy() for row in linear]
        for row in swapped:
            row[index], row[index+1] = row[index+1], row[index]
        assert product(restrict(generators[index]), linear) == swapped
    rays = [[Q(k-m)]*m+[Q(-m)]*(k-m) for m in range(2,k)]
    images = [apply(linear,ray) for ray in rays]
    assert all(v >= 0 for image in images for v in image)
    assert all(sum(image[row] for image in images) > 0 for row in range(dimension_minus))

    nodes = [(0,*rest) for rest in itertools.permutations(range(1,k))]
    node_index = {g:i for i,g in enumerate(nodes)}
    coefficients = [[None]*len(nodes) for _ in rays]
    for rest in itertools.permutations(range(2,k)):
        pi = (0,1,*rest)
        inverse = [pi.index(i) for i in range(k)]
        for label in range(2,k):
            g = tuple(label if member==1 else 1 if member==label else member for member in pi)
            for m,ray in enumerate(rays):
                z = [ray[i] for i in inverse]
                field = apply(linear,z)
                moved = apply(transports[label],field)
                residual = [field[row]-moved[row]+2*linear[row][label]*(z[1]-z[label])
                            for row in range(dimension_minus)]
                coefficients[m][node_index[g]] = residual
    denominator = lcm(*(v.denominator for c in coefficients for row in c if row for v in row))
    integers = [[[int(v*denominator) for v in c] if c else [0]*dimension_minus
                 for c in row] for row in coefficients]
    n = k-1
    word = [i for last in range(n-1) for i in range(last,-1,-1)]
    permutation = list(range(n))
    for i in word:
        assert permutation[i] < permutation[i+1]
        permutation[i], permutation[i+1] = permutation[i+1], permutation[i]
    assert permutation == list(reversed(range(n)))

    # Each input variable has exponent at most five. Four-bit exponent
    # slots make a monomial product ordinary integer addition without carry.
    lower = [[({0:1} if i==j else {}) for j in range(n)] for i in range(n)]
    upper = [[({0:1} if i==j else {}) for j in range(n)] for i in range(n)]

    def add(target, source, scale=1, shift=0):
        for exponent,value in source.items():
            key = exponent+shift
            target[key] = target.get(key,0)+value*scale

    def multiply(left,right):
        result = {}
        for exponent,value in left.items():
            add(result,right,value,exponent)
        return result

    for h,i in enumerate(word):
        for row in range(n):
            add(lower[row][i],lower[row][i+1],shift=1<<(4*h))
    for h,i in enumerate(reversed(word),len(word)):
        for row in range(n):
            add(upper[row][i+1],upper[row][i],shift=1<<(4*h))
    matrix = [[{} for i in range(n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            for h in range(n):
                add(matrix[i][j],multiply(lower[i][h],upper[h][j]),shift=1<<(4*(2*len(word)+h)))

    totals = [[{} for _ in range(dimension_minus)] for _ in rays]
    for i,g in enumerate(nodes):
        if all(not any(integers[m][i]) for m in range(len(rays))):
            continue
        monomial = {0:1}
        factors = sorted((matrix[j-1][g.index(j)-1] for j in range(1,k)),key=len)
        for factor in factors:
            monomial = multiply(monomial,factor)
        for m in range(len(rays)):
            for r,value in enumerate(integers[m][i]):
                if value:
                    add(totals[m][r],monomial,value)

    results = []
    parameter_count = n*n
    for m,ray_totals in enumerate(totals,2):
        for r,row_totals in enumerate(ray_totals,1):
            terms = sorted((exponent,value) for exponent,value in row_totals.items() if value)
            assert terms and all(value > 0 for exponent,value in terms), (m,r)
            encoded = [([int((exponent>>(4*i))&15) for i in range(parameter_count)],value)
                       for exponent,value in terms]
            digest = sha256(json.dumps(encoded,separators=(",",":")).encode()).hexdigest()
            results.append({
                "ray":m,
                "component":r,
                "nonzero_coefficients":len(terms),
                "minimum_coefficient":str(Q(min(v for e,v in terms),denominator)),
                "coefficient_sha256":digest,
            })
    return {"k":k,"shape":shape,"linear_field":[[str(v) for v in row] for row in linear],
            "metric_diagonal":[str(metric[i]) for i in minus],
            "chamber_ray_values":[[str(v) for v in row] for row in images],
            "jacobi_lower_word":[i+1 for i in word],
            "coefficient_denominator":denominator,
            "total_nonzero_coefficients":sum(item["nonzero_coefficients"] for item in results),
            "residuals":results}

if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
