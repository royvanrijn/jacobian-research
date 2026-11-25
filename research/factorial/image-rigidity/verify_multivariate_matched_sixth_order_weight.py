#!/usr/bin/env python3
"""Certify the symbolic sixth-order weight and independent root-reversion controls."""

from fractions import Fraction as Q
from functools import lru_cache
from math import factorial

from sympy import Poly, QQ, bernoulli, binomial, expand, invert, symbols

from verify_multivariate_matched_second_order_angles import Z, direct_coefficients, require


def derive_sixth_model():
    degree = 6
    j, w, v, k = symbols("j w v k")
    unsigned = symbols("V2 V3 V4 V5 V6 V7")
    signed = symbols("D2 D3 D4 D5 D6 D7")
    variance = (1, *unsigned)
    log_parts = [0] * (degree + 1)
    for m, moment in enumerate(variance, 1):
        prefix = (bernoulli(m+1, k) - bernoulli(m+1)) / (m+1)
        centered = expand(
            Q((-1)**(m+1), m) * (prefix - k * prefix.subs(k, j)/j).subs(k, v+j/2)
        )
        for index in range(1, degree+1):
            power = 2*m-index
            coefficient = centered.coeff(v, power) if power >= 0 else 0
            if coefficient:
                log_parts[index] += coefficient * w**power * (
                    signed[m-2] if power % 2 else moment
                )
    exp_parts = [1]
    for index in range(1, degree+1):
        exp_parts.append(expand(sum(
            m*log_parts[m]*exp_parts[index-m] for m in range(1, index+1)
        ) / index))

    # Unit-variance Hermite recurrence and its normalized derivatives at a root.
    basis = [(1, 0), (Z, 1)]
    for index in range(1, degree):
        basis.append(tuple(expand(Z*a+(j+index)*b)
                           for a, b in zip(basis[index], basis[index-1], strict=True)))
    derivatives = [0, 1]
    for index in range(3*degree):
        derivatives.append(expand((j-index)*derivatives[index]-Z*derivatives[index+1]))
    diff_polys = tuple((Poly(a, Z), Poly(b, Z)) for a, b in basis)

    @lru_cache(maxsize=None)
    def shifted_derivative(offset, order):
        a, b = diff_polys[offset]
        return expand(sum(binomial(order, i)*(
            a.diff((Z, order-i)).as_expr()*derivatives[i]
            + b.diff((Z, order-i)).as_expr()*derivatives[i+1]
        ) for i in range(order+1)))

    @lru_cache(maxsize=None)
    def apply(index, order, offset):
        return expand(sum((-1)**power * coeff * shifted_derivative(offset, power+order)
                          for (power,), coeff in Poly(exp_parts[index], w).terms()))

    center = (1, j/24, j*(5*j-2)/5760,
              j/322560+j*(j-1)/46080+j*(j-1)*(j-2)/82944)

    @lru_cache(maxsize=None)
    def derivative(index, order):
        return expand(sum(center[p//2]*apply(index-p, order, p)
                          for p in range(0, index+1, 2))/factorial(order))

    shifts = {}

    @lru_cache(maxsize=None)
    def shift_power(power, index):
        if power == 0:
            return int(index == 0)
        if index < power:
            return 0
        return expand(sum(shifts[t]*shift_power(power-1, index-t)
                          for t in range(1, index-power+2)))

    for index in range(1, degree+1):
        shifts[index] = -expand(sum(
            derivative(p, order)*shift_power(order, index-p)
            for p in range(index+1)
            for order in range(index-p+1)
            if (p, order) != (0, 1)
        ))

    # Scaling back to variance R turns j^a z^(2b) into J^a Y^b R^(3-a-b).
    model = Poly(expand(shifts[degree]/Z), j, Z, *unsigned, *signed, domain=QQ)
    require(all(b % 2 == 0 and 2*a+b <= 6 for a, b, *_ in model.monoms()),
            "sixth-order weighted degree")
    require(max(a+b//2 for a, b, *_ in model.monoms()) == 3, "maximal weight")
    require(max(sum(powers) for _, _, *powers in model.monoms()) == 6, "moment degree")
    mass = sum(abs(coefficient) for coefficient in model.coeffs())
    require(mass == Q(57500100121, 195955200) and mass < 300, ("sixth l1", mass))
    return model


def root_reversion(polys):
    anchor = polys[0]
    inverse = Poly(invert(anchor.diff(), anchor), Z, domain=QQ)
    zero, one = Poly(0, Z, domain=QQ), Poly(1, Z, domain=QQ)
    shifts = [zero] * len(polys)
    for index in range(1, len(polys)):
        residual = zero
        powers = [one] + [zero] * index
        for m in range(index+1):
            for k in range(index+1-m):
                term = polys[k].diff((Z, m)) * Q(1, factorial(m))
                residual = (residual + term*powers[index-k]).rem(anchor)
            powers = [sum((powers[k]*shifts[h-k] for k in range(h+1)), zero).rem(anchor)
                      for h in range(index+1)]
        shifts[index] = (-inverse*residual).rem(anchor)
    return shifts[-1]


def main():
    model = derive_sixth_model()
    samples = (
        ((Q(1, 2), 1), (Q(1, 3), -1), (Q(1, 6), 1)),
        ((Q(1, 3), 1), (Q(1, 3), -1), (Q(1, 6), 1), (Q(1, 6), -1)),
        ((Q(1), -1),),
        ((Q(3, 12), 1), (Q(4, 12), 1), (Q(5, 12), -1)),
    )
    count = 0
    for order in range(2, 13):
        for factors in samples:
            require(sum(x for x, _ in factors) == 1, "unit variance control")
            unsigned = [sum(x**m for x, _ in factors) for m in range(2, 8)]
            signed = [sum(sign*x**m for x, sign in factors) for m in range(2, 8)]
            values = (order, Z, *unsigned, *signed)
            expected = Poly(Z * model.as_expr().subs(dict(zip(model.gens, values, strict=True))),
                            Z, domain=QQ)
            polys = direct_coefficients(order, factors, offset=6)
            actual = root_reversion(polys)
            require((actual-expected).rem(polys[0]).is_zero,
                    ("sixth-order root reversion", order, factors))
            count += 1
    print(f"Verified exact sixth-order degree and norm, and {count} independent controls.")


if __name__ == "__main__":
    main()
