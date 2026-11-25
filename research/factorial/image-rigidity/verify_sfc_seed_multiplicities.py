"""Exclude every ordinary split/merged multiplicity of the interior-collision seed."""
from math import comb, factorial, isqrt

SPLIT = ((5, 2, 3), (2, 5, 3))
MERGED = ((7, 7, 6),)
BASE = ((5, 4, 4), (4, 5, 4), (5, 5, 4))
PERIOD = 36


def homogeneous_coefficients(joint_moment, sum_degree, difference_degree):
    """Coefficient list in t of E[(U+V+tW)^sum_degree (U-V)^difference_degree]."""
    return [
        comb(sum_degree, k) * sum(
            (-1) ** (difference_degree - j) * comb(difference_degree, j)
            * comb(sum_degree - k, i)
            * joint_moment(i + j, sum_degree - k - i + difference_degree - j, k)
            for i in range(sum_degree - k + 1)
            for j in range(difference_degree + 1)
        )
        for k in range(sum_degree + 1)
    ]


def projectively_coprime(first, second, prime):
    """Test invertibility of the homogeneous Sylvester matrix over F_prime."""
    first_degree, second_degree = len(first) - 1, len(second) - 1
    size = first_degree + second_degree
    rows = [
        [0] * offset + polynomial + [0] * (shift_count - 1 - offset)
        for polynomial, shift_count in ((first, second_degree), (second, first_degree))
        for offset in range(shift_count)
    ]
    for column in range(size):
        pivot = next((i for i in range(column, size) if rows[i][column] % prime), None)
        if pivot is None:
            return False
        rows[column], rows[pivot] = rows[pivot], rows[column]
        pivot_inverse = pow(rows[column][column], -1, prime)
        for i in range(column + 1, size):
            multiplier = rows[i][column] * pivot_inverse % prime
            for j in range(column + 1, size):
                rows[i][j] = (rows[i][j] - multiplier * rows[column][j]) % prime
    return True


def excludes(joint_moment, prime):
    def moments(sum_degree, difference_degree):
        return homogeneous_coefficients(joint_moment, sum_degree, difference_degree)

    def multiply(first, second):
        product = [0] * (len(first) + len(second) - 1)
        for i, left in enumerate(first):
            for j, right in enumerate(second):
                product[i+j] = (product[i+j] + left*right) % prime
        return product

    v2, = moments(0, 2)
    if v2 % prime == 0:
        return False
    v4, = moments(0, 4)
    m2 = moments(2, 0)
    cubic = [v2 * value for value in moments(3, 0)]
    for i, value in enumerate(multiply(m2, moments(1, 2))):
        cubic[i] = (cubic[i] - 3*value) % prime
    quartic = [v2*v2 * value for value in moments(4, 0)]
    for i, value in enumerate(multiply(m2, moments(2, 2))):
        quartic[i] = (quartic[i] - 6*v2*value) % prime
    for i, value in enumerate(multiply(m2, m2)):
        quartic[i] = (quartic[i] + v4*value) % prime
    return projectively_coprime(cubic, quartic, prime)


def factors(prime):
    assert all(prime % divisor for divisor in range(2, isqrt(prime) + 1))
    assert (prime - 1) % PERIOD == 0
    assert all(4 * max(row) < prime for row in BASE + SPLIT + MERGED)
    factorials = [factorial(n) % prime for n in range(prime)]
    result = {}
    for n in (2, 3, 4):
        for i in range(n + 1):
            for j in range(n - i + 1):
                k = n - i - j
                factors_at_index = []
                for rows in (BASE, SPLIT, MERGED):
                    weight = 1
                    for a, b, c in rows:
                        weight = weight * factorials[a*i + b*j + c*k] % prime
                    factors_at_index.append(weight)
                result[i, j, k] = factors_at_index
    return result


def unresolved_lifts(unresolved, prime):
    table = factors(prime)

    def survives(u, v):
        def joint_moment(i, j, k):
            base, split, merged = table[i, j, k]
            return base * pow(split, u, prime) * pow(merged, v, prime) % prime

        return not excludes(joint_moment, prime)

    return {
        (u, v) for u, v in unresolved if any(
            survives(ul, vl)
            for ul in range(u, prime - 1, PERIOD)
            for vl in range(v, prime - 1, PERIOD)
        )
    }


def check_controls():
    assert not projectively_coprime([-1, 1], [-1, 0, 1], 37)
    assert not projectively_coprime([1, 1, 0, 0], [1, 0, 1, 0, 0], 37)
    assert projectively_coprime([1, 1], [1, 0, 1], 37)

    # The law below has the collision (1+i, 1-i, -1) by elementary moments:
    # R,B are independent and E R² = E R⁴ = E B² = E B⁴ = 1/3.
    weights = {-1: 1, 0: 4, 1: 1}

    def collision(i, j, k):
        return sum(
            wr * wb * (4+r+b)**i * (4+r-b)**j * 8**k
            for r, wr in weights.items()
            for b, wb in weights.items()
        )

    assert not excludes(collision, 37)


def main():
    check_controls()
    unresolved = {(u, v) for u in range(PERIOD) for v in range(PERIOD)}
    for prime, expected in ((37, 38), (73, 4), (109, 0)):
        unresolved = unresolved_lifts(unresolved, prime)
        assert len(unresolved) == expected
        print(f"p={prime}: {len(unresolved)} unresolved residue pairs")
        if len(unresolved) <= 4:
            print(sorted(unresolved))
    assert not unresolved


if __name__ == "__main__":
    main()
