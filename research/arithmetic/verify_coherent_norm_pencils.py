"""Exact bounded examples for the coherent norm pencils (standard library)."""

from fractions import Fraction
from functools import cache
from itertools import combinations
from math import gcd, isqrt, prod


FIELDS = (23, 31, 59, 83)
PRIME_LIMIT = 151


@cache
def factor(n):
    factors = {}
    for p in range(2, isqrt(n) + 1):
        if p * p > n:
            break
        while n % p == 0:
            n //= p
            factors[p] = factors.get(p, 0) + 1
    if n > 1:
        factors[n] = 1
    return factors


def cube_free(n):
    result = 1
    for p, exponent in factor(abs(n)).items():
        result *= p ** (exponent % 3)
    return result


def prime_cubes(d):
    """Select the first signed generator (a+b sqrt(-d))/2 at each prime."""
    for p in range(3, PRIME_LIMIT + 1, 2):
        if factor(p) != {p: 1} or d % p == 0:
            continue
        for b in range(1, isqrt(4 * p**3 // d) + 1):
            a_squared = 4 * p**3 - d * b * b
            a = isqrt(a_squared)
            if a and a * a == a_squared and (a - b) % 2 == 0 and b % p:
                yield p, a, b
                break


def multiply(left, right, d):
    a, b = left
    c, e = right
    u, v = a * c - d * b * e, a * e + b * c
    assert (u % 2, v % 2) == (0, 0)
    return u // 2, v // 2


def partitions(edges, vertices):
    local = {}
    for a, b in combinations(vertices, 2):
        for p, exponent in factor(abs(edges[a, b][1])).items():
            for e in range(1, exponent + 1):
                local.setdefault((p, e), set()).add((a, b))
    for pairs in local.values():
        for a, b, c in combinations(vertices, 3):
            # An equivalence relation permits zero, one or three such edges.
            assert sum(edge in pairs for edge in ((a, b), (a, c), (b, c))) != 2


def normalized_traces(primes, coefficients):
    p1, p2, p3 = primes
    x, y, z = coefficients
    assert gcd(x, y) == gcd(x, z) == gcd(y, z) == 1
    assert gcd(x, p1 * p2) == gcd(y, p2 * p3) == gcd(z, p1 * p3) == 1
    a, b, c = p1**3 * y**2, p2**3 * z**2, p3**3 * x**2
    assert (a - b) % x == (c - b) % y == (a - c) % z == 0
    quotients = [divmod(n, d) for n, d in
                 ((a + b - c, y*z), (b + c - a, x*z), (b - a - c, x*y))]
    assert all(remainder == 0 for _, remainder in quotients)
    e, remainder = divmod(4*a*c - (a + c - b)**2, (x*y*z)**2)
    assert remainder == 0
    return tuple(trace for trace, _ in quotients), e


def odd_part(n):
    n = abs(n)
    return n // (n & -n)


def check_cayley(primes, coefficients, traces, defect, *, check_orbit=False):
    p1, p2, p3 = primes
    x, y, z = coefficients
    u12, u23, u13 = traces
    denominators = ((p2*p3)**3, (p1*p3)**3, (p1*p2)**3)
    xi = tuple(Fraction(u*u, denominator) - 2 for u, denominator in
               zip((u23, u13, u12), denominators))
    assert tuple(q.denominator for q in xi) == denominators
    assert sum(q*q for q in xi) - prod(xi) == 4
    w = prod(primes)**3
    a, b, c = xi
    q = (u12*u23*u13, defect*u12*y*z, defect*u23*x*z, -defect*u13*x*y)
    assert q == tuple(w * v for v in
                      (2+a+b+c, 2-a-b+c, 2+a-b-c, 2-a+b-c))
    assert sum(q) == 8*w and sum(Fraction(1, v) for v in q) == 0
    assert all(q[0] + v > 0 for v in q[1:])
    assert 8 % gcd(*q) == 0
    assert odd_part(gcd(*q[1:])) == odd_part(defect)
    assert odd_part(gcd(q[0], defect)) == 1
    for missing, expected in enumerate(
            (gcd(x, u23, u13), gcd(y, u12, u13), gcd(z, u12, u23)), 1):
        assert odd_part(gcd(*(v for i, v in enumerate(q) if i != missing))) == (
            odd_part(expected)**2
        )

    content = gcd(*q)
    assert content == (8 if defect % 2 == 0 else 1 if x*y*z % 2 else 2)
    node = tuple(gcd(*(v // content for i, v in enumerate(q) if i != missing))
                 for missing in range(4))
    assert node[0] == (defect if defect % 2 else defect // 4)
    divisor = 1 if defect % 2 else 2
    t12, t23, t13 = (value // divisor for value in traces)
    square_roots = (gcd(x, t23, t13), gcd(y, t12, t13), gcd(z, t12, t23))
    assert node[1:] == tuple(v*v for v in square_roots)
    a, b, c = p1**3*y*y, p2**3*z*z, p3**3*x*x
    for small, root, difference in zip(coefficients, square_roots,
                                       (a-b, c-b, a-c)):
        quotient, remainder = divmod(difference, small*root)
        assert remainder == 0
        assert odd_part(root) == odd_part(gcd(small, quotient*root))

    adjusted = [2 if defect % 2 and value % 2 == 0 else 1
                for value in coefficients]
    pair_factors = []
    reduced = tuple(value // content for value in q)
    residuals = [abs(reduced[v]) // prod(node[w] for w in range(4) if w != v)
                 for v in range(4)]
    for i, p in enumerate((p3, p1, p2), 1):
        j, k = (v for v in (1, 2, 3) if v != i)
        trace_scale = 2 if defect % 2 == 0 else adjusted[i - 1]
        a_i, a_remainder = divmod(abs(coefficients[i - 1]),
                                  adjusted[i - 1]*square_roots[i - 1])
        e_i, e_remainder = divmod(abs(traces[i - 1]),
                                  trace_scale*square_roots[j - 1]*square_roots[k - 1])
        assert a_remainder == e_remainder == 0
        assert gcd(residuals[j], residuals[k]) == a_i
        assert gcd(residuals[0], residuals[i]) == e_i
        factor, remainder = divmod(2*trace_scale**2, content)
        assert remainder == 0
        assert reduced[0] + reduced[i] == factor*p**3*node[j]*node[k]*e_i**2
        assert reduced[j] + reduced[k] == factor*p**3*node[0]*node[i]*a_i**2
        pair_factors.extend((a_i, e_i))
    assert all(gcd(a, b) == 1 for a, b in combinations(pair_factors, 2))

    p1_cube, p2_cube, p3_cube = (p**3 for p in primes)
    common = gcd(x, u23, u13)
    assert gcd(odd_part(common), defect*prod(primes)) == 1
    assert x*x*p3_cube == p2_cube*z*z - u12*z*y + p1_cube*y*y
    assert x*u23 == 2*p2_cube*z - u12*y
    assert (2*p2_cube*z - u12*y) % (x*common) == 0
    assert (u12*u12 - 4*p1_cube*p2_cube) % (x*x) == 0
    assert (x*x*p3_cube - p2_cube*z*z - p1_cube*y*y) % (y*z) == 0

    def walk(point, previous, prime, remaining_steps):
        if not remaining_steps:
            return
        for i in range(3):
            if i == previous:
                continue
            other = [point[j] for j in range(3) if j != i]
            changed = list(point)
            changed[i] = prod(other) - point[i]
            first_prime = primes[i] if prime is None else prime
            assert all(v.denominator % first_prime == 0 for v in changed)
            walk(changed, i, first_prime, remaining_steps - 1)

    if check_orbit:
        walk(xi, None, None, 6)


def run_field(d, orientations=None):
    generators = {p: (a, b) for p, a, b in prime_cubes(d)}
    if orientations is not None:
        generators = {p: (generators[p][0], sign * generators[p][1])
                      for p, sign in orientations.items()}
    roots = {}
    for p, (a, b) in generators.items():
        roots[p] = -a * pow(b, -1, p**3) % p**3
        assert (roots[p] ** 2 + d) % p**3 == 0
    edges = {}
    for first, second in combinations(generators, 2):
        a, b = generators[first]
        u, v = multiply((a, -b), generators[second], d)
        assert u and v and gcd(v, first * second) == 1
        assert u * u + d * v * v == 4 * (first * second) ** 3
        assert (u - v * roots[first]) % first**3 == 0
        assert (u + v * roots[second]) % second**3 == 0
        edges[first, second] = u, v

    def v(a, b):
        return edges[a, b][1]

    for index, (a, b, c) in enumerate(combinations(generators, 3)):
        u, anchor = edges[a, b]
        x, y = v(a, c), v(b, c)
        assert anchor % gcd(x, y) == 0
        assert anchor**2 * c**3 == b**3 * x**2 - u * x * y + a**3 * y**2
        up, vp = multiply(edges[a, b], edges[b, c], d)
        assert (up, vp) == tuple(b**3 * k for k in edges[a, c])
        original = (v(a, b), v(b, c), v(a, c))
        delta = gcd(*original)
        signed = tuple(value // delta for value in original)
        traces, e = normalized_traces((a, b, c), signed)
        assert e == d * delta**2 and gcd(e, a*b*c) == 1
        assert traces == (edges[a, b][0], edges[b, c][0], edges[a, c][0])

        check_cayley((a, b, c), signed, traces, e, check_orbit=(index == 0))

    counts = {2: 0, 3: 0}
    first_collision = None
    for a, b, c, e in combinations(generators, 4):
        products = v(a, b) * v(c, e), v(a, c) * v(b, e), v(a, e) * v(b, c)
        assert products[0] + products[2] == products[1]
        classes = tuple(map(cube_free, products))
        counts[len(set(classes))] += 1
        if len(set(classes)) == 2 and first_collision is None:
            first_collision = (a, b, c, e), products, classes

    five = tuple(generators)[:5]
    if orientations is None:
        assert len(five) == 5
    for four in combinations(five, 4):
        partitions(edges, four)
    partitions(edges, five)
    print(f"d={d}: {len(generators)} primes; opposite class counts {counts}")
    print(f"  signed generators: {[(p, generators[p]) for p in five]}")
    if first_collision:
        vertices, products, classes = first_collision
        print(f"  collision {vertices}: products {products}; classes {classes}")
        print(f"  signed edges: {[(pair, edges[pair]) for pair in combinations(vertices, 2)]}")
    return counts


def relaxed_rank_one_record():
    p1, p2, p3 = 11, 13, 19
    x, y, z = 389, 47, 299  # The full-modulus record uses j = 2*z.
    assert gcd(x, y) == gcd(x, z) == gcd(y, z) == 1
    assert gcd(x, p1 * p2) == gcd(y, p2 * p3) == gcd(z, p1 * p3) == 1
    assert cube_free(x) != cube_free(y)
    a, b, c = p1**3, p2**3, p3**3
    assert (y * y * a - z * z * b) % x == 0
    assert (x * x * c - z * z * b) % y == 0
    assert (y * y * a - x * x * c) % z == 0
    assert (4 * x * x * y * y * a * c - (y * y * a + x * x * c - z * z * b)**2
            == -700867140113231117)
    c1 = p1 * pow(p2, -1, x) % x
    c3 = p3 * pow(p2, -1, y) % y
    radius = 38
    small_vectors = {
        (i, j, k)
        for j in range(-radius, radius + 1)
        for i in range(-radius + (c1 * j + radius) % x, radius + 1, x)
        for k in range(-radius + (c3 * j + radius) % y, radius + 1, y)
        if i * i + j * j + k * k <= radius**2 and (p3 * i - p1 * k) % z == 0
    }
    assert small_vectors == {(0, 0, 0), (p1, p2, p3), (-p1, -p2, -p3)}
    assert normalized_traces((p1, p2, p3), (x, y, z))[1] < 0
    print("Relaxed rank-one record has no imaginary quadratic norm lift.")


if __name__ == "__main__":
    totals = {2: 0, 3: 0}
    for field in FIELDS:
        for key, value in run_field(field).items():
            totals[key] += value
    assert totals == {2: 0, 3: 5642}
    assert run_field(43, {17: 1, 59: -1, 127: 1, 139: 1}) == {2: 1, 3: 0}
    run_field(2)
    _, e = normalized_traces((3, 5, 7), (1, 2, 1))
    assert e == 10475 and gcd(e, 3*5*7) == 5
    print("Positive-Heron control fails the incident-prime unit test at 5.")
    relaxed_rank_one_record()
