"""Exact Newton residual certificates for constant-plus-prime-power moments.

All-index conclusions use PRIME_POWER_TWO_TERM_SFC.tex. The exporter builds
actual finite polynomials and checks their Newton edges; it never infers an
all-degree theorem from a scan. Only exact rational arithmetic is accepted.
"""
from fractions import Fraction
from math import factorial, comb, gcd
import sympy as s


def integer(value, minimum=0):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)) or value < minimum:
        raise ValueError('an exact integer >= '+str(minimum)+' is required')
    return int(value)


def prime_power(degree):
    degree = integer(degree, 2)
    factors = s.factorint(degree)
    if len(factors) != 1:
        raise ValueError('the exponent must be a prime power')
    prime, exponent = next(iter(factors.items()))
    return int(prime), int(exponent)


def valuation(value, prime):
    value = Fraction(value)
    if not value:
        raise ValueError('zero has no finite valuation')
    result = 0
    numerator, denominator = abs(value.numerator), value.denominator
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def residue(value, prime):
    value = Fraction(value)
    if value.denominator % prime == 0:
        raise ValueError('nonintegral residue input')
    return value.numerator*pow(value.denominator, -1, prime) % prime


def p_power(prime, exponent):
    return Fraction(prime)**exponent


def factorial_valuation(n, prime):
    result = 0
    while n:
        n //= prime
        result += n
    return result


def newton_residuals(coefficients, prime):
    """Literal lower hull and normalized edge residuals from rational coefficients."""
    prime = integer(prime, 2)
    if not s.isprime(prime):
        raise ValueError('a prime is required')
    coefficients = list(map(Fraction, coefficients))
    if not coefficients or not coefficients[0] or not coefficients[-1]:
        raise ValueError('nonzero constant and leading coefficients required')
    points = [(k, valuation(c, prime)) for k, c in enumerate(coefficients) if c]
    hull = []
    for point in points:
        while len(hull) >= 2:
            (x0, y0), (x1, y1) = hull[-2:]
            if (x1-x0)*(point[1]-y1)-(y1-y0)*(point[0]-x1) > 0:
                break
            hull.pop()
        hull.append(point)
    result = []
    T = s.Symbol('T')
    for left, right in zip(hull, hull[1:]):
        width, height = right[0]-left[0], right[1]-left[1]
        divisor = gcd(width, abs(height))
        e, h = width//divisor, height//divisor
        base_unit = residue(coefficients[left[0]]/p_power(prime, left[1]), prime)
        values = []
        for j in range(divisor+1):
            index = left[0]+j*e
            target = left[1]+j*h
            coefficient = coefficients[index]
            assert not coefficient or valuation(coefficient, prime) >= target
            values.append(residue(coefficient/p_power(prime, target), prime)*pow(base_unit, -1, prime) % prime)
        polynomial = s.Poly(sum(c*T**j for j, c in enumerate(values)), T, modulus=prime)
        result.append({'left': list(left), 'right': list(right), 'e': e, 'h': h,
                       'residual_coefficients': values,
                       'residual_squarefree': s.gcd(polynomial, polynomial.diff()).degree() == 0})
    return result


def expected_edges(n, prime, gamma):
    """The formula proved for falling-factorial coefficients with geometric units."""
    n = integer(n, 1)
    if not s.isprime(prime) or not gamma % prime:
        raise ValueError('prime and nonzero residue are required')
    remaining, place, left, height = n, 0, 0, 0
    result = []
    while remaining:
        digit = remaining % prime
        remaining //= prime
        if digit:
            e = prime**place
            h = (e-1)//(prime-1)
            epsilon = pow(gamma, e, prime)*((-1)**h) % prime
            falling = 1
            values = [1]
            for j in range(1, digit+1):
                falling *= digit-j+1
                values.append(falling*pow(epsilon, j, prime) % prime)
            result.append({'left': [left, height], 'right': [left+digit*e, height+digit*h],
                           'e': e, 'h': h, 'residual_coefficients': values,
                           'residual_squarefree': True})
            left += digit*e
            height += digit*h
        place += 1
    return result


def class_certificate(degree, shape=1):
    degree = integer(degree, 2)
    prime, exponent = prime_power(degree)
    shape = integer(shape, 1)
    if shape > degree:
        raise ValueError('this theorem requires integer Gamma shape <= exponent')
    return {'status': 'two_consecutive_moments_cannot_both_vanish',
            'exponent': degree, 'prime': prime, 'prime_exponent': exponent, 'gamma_shape': shape,
            'phase': 'a+b*Y^exponent', 'coefficients': 'arbitrary complex a,b; nonzero phase',
            'starts': 'every integer n>=1', 'window': 2,
            'all_index_mechanism': 'geometric-unit falling-factorial Newton residuals are reversed truncated exponentials',
            'theorem': 'FIRSFCMONOPRIME20261116', 'finite_search_used_for_theorem': False,
            'independent_replay': False, 'formal_verification': False, 'external_review': False}


def moment_polynomial(n, degree, shape=1, variable=None):
    n, degree, shape = integer(n), integer(degree, 1), integer(shape, 1)
    variable = s.Symbol('lambda') if variable is None else variable
    return s.Poly.from_list([s.Integer(comb(n, k))*s.rf(shape, degree*k)
                            for k in range(n, -1, -1)], variable, domain=s.QQ)


def polynomial_certificate(n, degree, shape=1):
    n = integer(n, 1)
    result = class_certificate(degree, shape)
    prime, exponent = result['prime'], result['prime_exponent']
    b = (degree-1)//(prime-1)
    polynomial = moment_polynomial(n, degree, shape)
    scaled = [Fraction(int(polynomial.nth(k)), prime**(b*k)) for k in range(n+1)]
    edges = newton_residuals(scaled, prime)
    gamma = (-1)**b % prime
    assert edges == expected_edges(n, prime, gamma)
    assert all(row['residual_squarefree'] for row in edges)
    result.update(n=n, moment_coefficients=[str(polynomial.nth(k)) for k in range(n+1)],
                  parameter_scaling='lambda=x/p^'+str(b), newton_edges=edges,
                  ordinary_pair_profile={'moment_indices': [n, n+1],
                      'complete_intersection_degrees': [n, n+1], 'length': n*(n+1),
                      'highest_nonzero_degree': 2*n-1, 'top_graded_dimension': 1,
                      'actual_coefficient_ideal_nilpotence_bound': 2*n,
                      'universal_bound_is_sharp': True,
                      'theorem': 'FIRPRIMEPAIRI20261116'},
                  certificate='all edge residuals squarefree; weighted-initial-product lemma excludes every square factor')
    return result


def marked_class_certificate(degree, divided_power_coefficients):
    """M(Y)=sum a_e*Y**e/e!, degree<degree; a_e/mass p-integral."""
    degree = integer(degree, 2)
    prime, exponent = prime_power(degree)
    values = [s.sympify(value) for value in divided_power_coefficients]
    if not values or any(value.has(s.Float) or value.is_Rational is not True for value in values):
        raise ValueError('exact rational divided-power marker coefficients required')
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    if len(values) > degree:
        raise ValueError('marker degree must be below the prime-power exponent')
    mass = sum(values)
    if mass == 0 or any(value and valuation(Fraction(value/mass), prime) < 0 for value in values):
        return {'status': 'outside_unit_mass_marker_criterion', 'period_is_zero': None,
                'reason': 'normalized divided-power coefficients must be p-integral and mass nonzero'}
    return {'status': 'two_consecutive_marked_moments_cannot_both_vanish',
            'exponent': degree, 'prime': prime, 'prime_exponent': exponent,
            'marker_divided_power_coefficients': [str(value) for value in values],
            'marker_mass': str(mass), 'window': 2, 'starts': 'every integer n>=1',
            'phase': 'a+b*Y^exponent, nonzero; Y is mean-one exponential',
            'theorem': 'FIRPRIMEMARKERS20261116',
            'all_index_mechanism': 'normalized block factorial multiplier reduces to constant1 modulo p',
            'scope': 'Outside this sufficient marker criterion no vanishing decision is made.'}
