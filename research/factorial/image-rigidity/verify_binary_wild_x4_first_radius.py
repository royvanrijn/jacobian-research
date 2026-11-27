"""Verify the first critical X^4 moment identities over Z[i,pi]/(i²+1, pi³-3).

The proof reduces all starting orders and all prime-power axes to the
162 representative systems generated here. Python 3.12, standard library.
"""

from functools import lru_cache
from itertools import product
from math import comb

PRECISION = 12
MODULUS = 81


def gaussian_multiply(left, right):
    a, b = left
    c, d = right
    return ((a * c - b * d) % MODULUS, (a * d + b * c) % MODULUS)


def gaussian_power(value, exponent):
    if exponent < 0:
        real, imag = value
        inverse_norm = pow(real**2 + imag**2, -1, MODULUS)
        return gaussian_power((real * inverse_norm, -imag * inverse_norm), -exponent)
    result = (1, 0)
    while exponent:
        if exponent % 2:
            result = gaussian_multiply(result, value)
        value = gaussian_multiply(value, value)
        exponent //= 2
    return result


def term(degree=0, a=0, b=0, c=0, coefficient=(1, 0)):
    if degree >= PRECISION:
        return {}
    shift, power = divmod(degree, 3)
    modulus = 3**((PRECISION - power + 2) // 3)
    value = tuple(3**shift * part % modulus for part in coefficient)
    return {(power, a, b, c): value} if any(value) else {}


ONE = term()
B = term(b=1)


def add(*polynomials):
    result = {}
    for poly in polynomials:
        for key, (x, y) in poly.items():
            a, b = result.get(key, (0, 0))
            value = ((a + x) % MODULUS, (b + y) % MODULUS)
            if any(value):
                result[key] = value
            else:
                result.pop(key, None)
    return result


def multiply(left, right):
    result = {}
    for (p, a, b, c), x in left.items():
        for (q, d, e, f), y in right.items():
            key = ((p + q) % 3, a + d, b + e, c + f)
            factor = 3 if p + q >= 3 else 1
            u, v = gaussian_multiply(x, y)
            r, s = result.get(key, (0, 0))
            value = ((r + factor * u) % MODULUS, (s + factor * v) % MODULUS)
            if any(value):
                result[key] = value
            else:
                result.pop(key, None)
    return result


def scale(poly, coefficient, degree=0):
    return multiply(term(degree, coefficient=coefficient), poly)


def residue_at(poly, order):
    factor = 3**(order // 3)
    result = {}
    for (power, a, b, c), (x, y) in poly.items():
        if power == order % 3:
            value = (x // factor % 3, y // factor % 3)
            if any(value):
                result[a, b, c] = value
    return result


def initial(poly):
    for order in range(PRECISION):
        part = residue_at(poly, order)
        if part:
            return order, part
    return PRECISION, {}


def evaluate(poly, values):
    result = (0, 0)
    for exponents, coefficient in poly.items():
        for value, exponent in zip(values, exponents):
            coefficient = gaussian_multiply(coefficient, gaussian_power(value, exponent))
        result = tuple((a + b) % 3 for a, b in zip(result, coefficient))
    return result


def lift_c(poly):
    """Substitute C = 1+i + pi*D. The last exponent now labels D."""
    result = {}
    for (power, a, b, c), value in poly.items():
        for j in range(c + 1):
            factor = gaussian_multiply((comb(c, j), 0), gaussian_power((1, 1), c-j))
            result = add(result, term(power+j, a, b, j, gaussian_multiply(value, factor)))
    return result


class Factorials:
    def __init__(self):
        self.free = [1]

    @lru_cache(None)
    def measure(self, n):
        """Return v_3(n!) and the 3-free factorial modulo 81."""
        if n == 0:
            return 0, 1
        exponent, unit = self.measure(n // 3)
        while len(self.free) <= n:
            j = len(self.free)
            self.free.append(self.free[-1] * (j if j % 3 else 1) % MODULUS)
        return exponent + n // 3, unit * self.free[n] % MODULUS

    def ratio(self, numerators, denominators):
        exponent, numerator, denominator = 0, 1, 1
        for n, power in numerators:
            value, unit = self.measure(n)
            exponent += power * value
            numerator = numerator * pow(unit, power, MODULUS) % MODULUS
        for n, power in denominators:
            value, unit = self.measure(n)
            exponent -= power * value
            denominator = denominator * pow(unit, power, MODULUS) % MODULUS
        return exponent, numerator * pow(denominator, -1, MODULUS) % MODULUS


class Moments:
    def __init__(self, qx, qy):
        self.factorials = Factorials()
        self.b = add(term(coefficient=(0, 1)), term(1, b=1))
        u = add(term(1, coefficient=(-1, -1)), term(2, a=1, coefficient=(-1, 0)), term(4, c=1))
        max_order = 8 + 27
        axis_powers = [ONE]
        foreign_powers = [ONE]
        for _ in range(max_order):
            axis_powers.append(multiply(axis_powers[-1], self.b))
        for _ in range(5):
            foreign_powers.append(multiply(foreign_powers[-1], u))
        self.native = []
        for n in range(max_order + 1):
            result = {}
            for m in range(min(n, 5) + 1):
                for y in range(n - m + 1):
                    x = n - m - y
                    value, unit = self.factorials.ratio(
                        [(qx*x + 4*m, 1), (qy*y, 1)],
                        [(x, 1), (y, 1), (m, 1), (qx, x), (qy, y), (4, m)],
                    )
                    if value < 0:
                        raise ArithmeticError("native coefficient lost integrality")
                    factor = unit * 3**value % MODULUS
                    if factor:
                        result = add(result, scale(multiply(axis_powers[y], foreign_powers[m]), (factor, 0)))
            self.native.append(result)

    @lru_cache(None)
    def unit_power(self, ell):
        """Expand (1+i+pi*A)^ell, including the negative powers."""
        result = {}
        for j in range(PRECISION):
            coefficient = comb(ell, j) if ell >= 0 else (-1)**j * comb(j-ell-1, j)
            value = gaussian_multiply((coefficient, 0), gaussian_power((1, 1), ell-j))
            result = add(result, term(j, a=j, coefficient=value))
        return result

    def forms(self, block):
        start = 9 * block
        scalar = {}
        for ell in range(-27, 9):
            if start + ell < 0:
                scalar[ell] = {}
                continue
            val, unit = self.factorials.ratio([(start, 1)], [(start + ell, 1)])
            order = ell + 3 * val
            if order < 0:
                raise ArithmeticError("nonintegral normalized constant factor")
            scalar[ell] = scale(self.unit_power(ell), (unit, 0), order)
        moments = [
            add(*(multiply(scalar[ell], self.native[i-ell]) for ell in range(-27, i+1)))
            for i in range(5, 9)
        ]
        e5, e6, e7, e8 = moments
        d = scale(add(ONE, self.b), (-1, 0))
        z7 = add(e7, multiply(d, e6), multiply(self.b, e5))
        z8 = add(e8, multiply(d, e7), multiply(self.b, e6))
        return e5, e6, z7, z8


def require(condition, context):
    if not condition:
        raise ArithmeticError(context)


def check(qx, qy, H, moments):
    e5, e6, z7, z8 = moments.forms(H)
    constant, cube_a, cube_b, c, c2 = (0, 0, 0), (3, 0, 0), (0, 3, 0), (0, 0, 1), (0, 0, 2)
    context = (qx, qy, H)
    require(initial(e5) == (3, {
        cube_a: (0, 2), cube_b: (0, 1),
        **({constant: (H % 3, H % 3)} if H % 3 else {}),
    }), (context, "E5 initial"))
    require(initial(e6) == (3, {cube_a: (1, 2), cube_b: (2, 0)}), (context, "E6 initial"))
    require(initial(z7) == (4, {constant: (2, 1), c: (0, 2)}), (context, "Z7 initial"))
    offset = (int(qx > 27) - H) % 3
    require(initial(z8) == (8, {
        c: (2, 1), c2: (0, 1), **({constant: (offset, 0)} if offset else {}),
    }), (context, "Z8 initial"))
    a0, b0, c0 = ((-H) % 3, H % 3), (H % 3, 0), (1, 1)
    if evaluate(initial(z8)[1], (a0, b0, c0)) != (0, 0):
        return
    e5, e6, z7, z8 = map(lift_c, (e5, e6, z7, z8))
    if qx == 27:
        order, quotient = initial(z8)
        require(order == 9 and all(m[2] == 0 for m in quotient), (context, "Qx=27 independence"))
        require(evaluate(quotient, (a0, b0, (0, 0))) == (2, 2), (context, "Qx=27 obstruction"))
    elif qy == 1:
        form = add(z8, scale(e5, (0, 1), 6), scale(e6, (2, 0), 6))
        order, quotient = initial(form)
        require(order == 10 and evaluate(quotient, ((0, 0),)*3) == (0, 2), (context, "Qy=1 obstruction"))
        require(initial(z7) == (5, {(1, 0, 0): (0, 1), (0, 0, 1): (0, 2)}), (context, "D residue"))
    else:
        form = add(z8, scale(e5, (2, 1), 6), scale(multiply(z7, z7), (0, 1)),
                   scale(multiply(B, e5), (1, 2), 7), scale(multiply(B, e6), (0, 1), 7),
                   scale(z7, (2, 2), 6))
        order, quotient = initial(form)
        require(order == 11, (context, "last order"))
        expected = {(0, 1, 0): (1, 0), (0, 5, 0): (2, 2), (3, 2, 0): (1, 0)}
        require({m: z for m, z in quotient.items() if m != constant} == expected,
                (context, "last polynomial"))
        fourth_e5, fourth_e6 = residue_at(e5, 4), residue_at(e6, 4)
        require(all(m[2] == 0 for m in (*fourth_e5, *fourth_e6)),
                (context, "degree-four independence"))
        require(constant not in fourth_e5 and fourth_e6[constant] == (0, 1),
                (context, "contradictory degree-four constants"))


def main():
    for qx, qy in product((27, 81, 243), (1, 3)):
        moments = Moments(qx, qy)
        for H in range(27):
            check(qx, qy, H, moments)
    print("Verified 162 linked X^4 systems at precision pi^12 (pi^3=3).")


if __name__ == "__main__":
    main()
