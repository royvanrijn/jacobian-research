#!/usr/bin/env python3
"""Check matched angular coefficients through fourth order by exact algebra."""

from fractions import Fraction as Q
from math import comb, factorial

from sympy import Poly, QQ, Symbol, bernoulli, binomial, expand, invert, symbols


Z = Symbol("z")


def require(condition, context):
    if not condition:
        raise ArithmeticError(context)


def factor_logs(order, index, factors, max_power):
    logs = [Q(0)] * (max_power + 1)
    for reciprocal, orientation in factors:
        k = index if orientation == 1 else order - index
        for power in range(1, max_power + 1):
            chord = Q(sum(h**power for h in range(k))) - Q(
                k * sum(h**power for h in range(order)), order
            )
            logs[power] += Q((-1)**(power+1), power) * reciprocal**power * chord
    return logs


def exp_coefficients(logs):
    coeffs = [Q(1)]
    for degree in range(1, len(logs)):
        coeffs.append(sum(
            power * logs[power] * coeffs[degree-power]
            for power in range(1, degree+1)
        ) / degree)
    return coeffs


def direct_coefficients(order, factors, offset=2):
    limit = order + offset
    max_power = limit // 2
    terms = [dict() for _ in range(limit + 1)]
    for k in range(order + 1):
        center = Q(2*k-order, 2)
        coeffs = exp_coefficients(factor_logs(order, k, factors, max_power))
        signed_binomial = (-1)**k * comb(order, k)
        for power, coeff in enumerate(coeffs):
            for degree in range(limit + 1 - 2*power):
                total = terms[2*power + degree]
                contribution = signed_binomial * coeff * (-center)**degree / factorial(degree)
                total[(degree,)] = total.get((degree,), Q(0)) + contribution
    require(not any(any(row.values()) for row in terms[:order]), ("exponential-scale cancellation", order))
    return tuple(Poly.from_dict(row, Z, domain=QQ) for row in terms[order:])


def hermite(order, variance):
    return Poly.from_dict({
        (order - 2*power,): Q(factorial(order),
            2**power * factorial(power) * factorial(order - 2*power)) * variance**power
        for power in range(order // 2 + 1)
    }, Z, domain=QQ)


def verify_formula(order, factors):
    variance = sum(x for x, _ in factors)
    v2 = sum(x*x for x, _ in factors)
    v3 = sum(x**3 for x, _ in factors)
    signed2 = sum(sign*x*x for x, sign in factors)
    first, second, third = direct_coefficients(order, factors)
    anchor = hermite(order, variance)
    anchor_derivative = anchor.diff()
    d = signed2 / 6

    expected_second = d * anchor.diff((Z, 3))
    expected_third = (
        -variance * Q(order**2, 8) * anchor
        -Q(order - 1, 4) * v2 * anchor.diff((Z, 2))
        +v3 / 12 * anchor.diff((Z, 4))
        +d*d / 2 * anchor.diff((Z, 6))
        +Q(order, 24) * hermite(order + 2, variance)
    )
    require(first == anchor and second == expected_second and third == expected_third,
            ("formal coefficients", order, factors))

    # On an imaginary root, z**2 is minus the squared height.
    angular_model = Poly(
        Z * (
            -order * variance / 24 - Q(order-1, 4) * v2 / variance
            +v3 / (12 * variance**2) * (2*order-3 + Z**2 / variance)
            -d*d / variance**3 * (Q(7*order, 2)-5 + 2*Z**2 / variance)
        ), Z, domain=QQ
    )
    # Clear the denominators of the two root-displacement coefficients.
    # Division by the anchor tests the identity at all its roots simultaneously.
    numerator = (
        -2 * third * anchor_derivative**2
        +2 * second * second.diff() * anchor_derivative
        -second**2 * anchor_derivative.diff()
        -2 * angular_model * anchor_derivative**3
    )
    require(numerator.rem(anchor).is_zero, ("angular correction", order, factors))


def ordinary_factors(order, rows):
    return tuple(
        (Q(gap, order*lower + h), orientation)
        for lower, gap, orientation in rows
        for h in range(1, gap+1)
    )


def derive_fourth_model():
    j, r, v, w, k = symbols("j R v w k")
    v2, v3, v4, v5, d2, d3, d4 = symbols("V2 V3 V4 V5 D2 D3 D4")
    unsigned = (None, r, v2, v3, v4, v5)
    signed = {2: d2, 3: d3, 4: d4}
    # Centered Gauss logarithms, collected by powers of lambda after w=lambda*v.
    log_parts = [0] * 5
    for m in range(1, 6):
        prefix = (bernoulli(m+1, k) - bernoulli(m+1)) / (m+1)
        chord = prefix - k * prefix.subs(k, j) / j
        centered = expand(Q((-1)**(m+1), m) * chord.subs(k, j/2+v))
        for degree in range(1, 5):
            power = 2*m-degree
            coeff = centered.coeff(v, power) if power >= 0 else 0
            if coeff:
                mass = signed[m] if power % 2 else unsigned[m]
                log_parts[degree] += coeff * mass * w**power
    exp_parts = [1]
    for degree in range(1, 5):
        exp_parts.append(expand(sum(
            m * log_parts[m] * exp_parts[degree-m] for m in range(1, degree+1)
        ) / degree))

    # H_{j+p}=A_p H_j+B_p H_j'. At an anchor root, K_m=H_j^(m)/H_j'.
    basis = [(Poly(1, Z), Poly(0, Z)), (Poly(Z, Z), Poly(r, Z))]
    for p in range(1, 4):
        basis.append(tuple(Poly(Z*a.as_expr()+r*(j+p)*b.as_expr(), Z)
                           for a, b in zip(basis[p], basis[p-1], strict=True)))
    normalized_derivatives = [0, 1]
    for degree in range(12):
        normalized_derivatives.append(expand(
            ((j-degree)*normalized_derivatives[degree]
             -Z*normalized_derivatives[degree+1])/r
        ))

    def shifted_derivative(p, degree):
        a, b = basis[p]
        return sum(binomial(degree, i)*(
            a.diff((Z, degree-i)).as_expr()*normalized_derivatives[i]
            +b.diff((Z, degree-i)).as_expr()*normalized_derivatives[i+1]
        ) for i in range(degree+1))

    finite_differences = (1, j/24, j*(5*j-2)/5760)

    def transformed(degree, derivative):
        return expand(sum(
            finite_differences[p//2] * (-1)**m * coeff
            * shifted_derivative(p, m+derivative)
            for p in range(0, degree+1, 2)
            for (m,), coeff in Poly(exp_parts[degree-p], w).terms()
        ))

    u = {(degree, derivative): transformed(degree, derivative)
         for degree in range(1, 5) for derivative in range(5-degree)}
    k2, k3, k4 = normalized_derivatives[2:5]
    a1 = -u[1, 0]
    a2 = expand(-k2*a1*a1/2 - u[1, 1]*a1 - u[2, 0])
    a3 = expand(-k2*a1*a2 - k3*a1**3/6 - u[1, 1]*a2
                - u[1, 2]*a1*a1/2 - u[2, 1]*a1 - u[3, 0])
    a4 = expand(-(
        k2*(a1*a3+a2*a2/2) + k3*a1*a1*a2/2 + k4*a1**4/24
        + u[1, 1]*a3 + u[1, 2]*a1*a2 + u[1, 3]*a1**3/6
        + u[2, 1]*a2 + u[2, 2]*a1*a1/2 + u[3, 1]*a1 + u[4, 0]
    ) / Z)
    require(all(exponent % 2 == 0 for (exponent,), _ in Poly(a4, Z).terms()),
            "quartic coefficient is even in the anchor height")

    jr, y = symbols("J Y")
    moments = symbols("u2 u3 u4 u5 t2 t3 t4")
    variables = (jr, r, y, *moments)
    raw_moments = (v2, v3, v4, v5, d2, d3, d4)
    moment_orders = (2, 3, 4, 5, 2, 3, 4)
    moment_values = (r**m*a for m, a in zip(moment_orders, moments, strict=True))
    normalized = dict(zip(raw_moments, moment_values, strict=True))
    normalized.update({j: jr/r, Z**2: y})
    model = Poly(a4.subs(normalized), *variables, domain=QQ)
    require(all(sum(exponents[:3]) == 2 for exponents in model.monoms()),
            "quartic coefficient homogeneous in J, R, Y")
    require(max(sum(exponents[3:]) for exponents in model.monoms()) == 4,
            "quartic coefficient moment degree")
    mass = sum(abs(coefficient) for coefficient in model.coeffs())
    require(mass == Q(43859, 3240) and mass < 14, ("quartic coefficient l1", mass))
    return model


def fourth_root_reversion(polys):
    # Independent generic substitution in Q[z]/H, using direct Gauss coefficients.
    base = polys[0]
    inverse = Poly(invert(base.diff(), base), Z, domain=QQ)
    zero, one = Poly(0, Z, domain=QQ), Poly(1, Z, domain=QQ)
    shifts = [zero] * 5
    for degree in range(1, 5):
        residual = zero
        powers = [one] + [zero] * degree
        for m in range(degree+1):
            for k in range(degree+1-m):
                term = polys[k].diff((Z, m)) * (Q(1, factorial(m)))
                residual = (residual + term * powers[degree-k]).rem(base)
            powers = [sum((powers[k] * shifts[h-k] for k in range(h+1)), zero).rem(base)
                      for h in range(degree+1)]
        shifts[degree] = (-inverse * residual).rem(base)
    return shifts[4]


def verify_quartic_model(order, factors, model):
    r = sum(x for x, _ in factors)
    unsigned = [sum(x**m for x, _ in factors)/r**m for m in range(2, 6)]
    signed = [sum(sign*x**m for x, sign in factors)/r**m for m in range(2, 5)]
    values = (order*r, r, Z**2, *unsigned, *signed)
    predicted = Poly(Z * model.as_expr().subs(dict(zip(model.gens, values, strict=True))),
                     Z, domain=QQ)
    polys = direct_coefficients(order, factors, offset=4)
    actual = fourth_root_reversion(polys)
    require((actual-predicted).rem(polys[0]).is_zero,
            ("fourth-order angular correction", order, factors))


def main():
    samples = (
        ((Q(1, 3), 1), (Q(1, 4), 1), (Q(1, 7), -1)),
        ((Q(1, 2), 1), (Q(1, 2), -1)),
        ((Q(1, 2), 1), (Q(1, 3), -1), (Q(1, 6), -1)),
        ((Q(4, 5), 1),),
    )
    ordinary = (
        ((7, 1, 1), (11, 3, -1)),
        ((12, 2, 1), (12, 2, -1)),
        ((19, 1, 1), (37, 1, -1), (63, 2, -1)),
    )
    tested = 0
    for order in range(2, 12):
        for factors in samples:
            verify_formula(order, factors)
            tested += 1
        for rows in ordinary:
            lower = ordinary_factors(order, rows)
            upper = ordinary_factors(order+1, rows)
            total = sum(x for x, _ in lower)
            for degree, factors in ((order, lower), (order+1, upper)):
                verify_formula(degree, tuple((x/total, sign) for x, sign in factors))
                tested += 1
    print(f"Verified {tested} exact second-order coefficient and root-reduction checks.")
    model = derive_fourth_model()
    fourth_samples = (
        ((Q(1, 2), 1), (Q(1, 3), -1), (Q(1, 6), 1)),
        ((Q(1, 3), 1), (Q(1, 3), -1), (Q(1, 7), 1), (Q(1, 7), -1)),
        ((Q(4, 5), -1),),
        ((Q(3, 20), 1), (Q(4, 20), 1), (Q(5, 20), -1)),
    )
    rows = ((11, 2, 1), (16, 3, -1), (42, 1, 1))
    fourth_tested = 0
    for order in range(2, 15):
        for factors in fourth_samples:
            verify_quartic_model(order, factors, model)
            fourth_tested += 1
        lower = ordinary_factors(order, rows)
        upper = ordinary_factors(order+1, rows)
        total = sum(x for x, _ in lower)
        for degree, factors in ((order, lower), (order+1, upper)):
            verify_quartic_model(degree, tuple((x/total, sign) for x, sign in factors), model)
            fourth_tested += 1
    print(f"Verified exact symbolic quartic degree and {fourth_tested} independent controls.")


if __name__ == "__main__":
    main()
