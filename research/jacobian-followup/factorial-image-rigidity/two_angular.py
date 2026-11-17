"""Exact Stokes and arithmetic certificates for q=X**3-3*X+Y**2.

The all-exponent theorem is in TWO_ANGULAR_COORDINATES.tex. A zero
specialized remainder is not a zero-period certificate. No quadrature or
finite moment prefix is used to assert universal nonvanishing.
"""
from functools import lru_cache
import sympy as s

X, Y, h = s.symbols('X Y h')
q = X**3-3*X+Y**2
u, v, w = (2-X-2*Y)/4, (2-X+2*Y)/4, X/2


def _integer(value, minimum=0):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)) or value < minimum:
        raise ValueError('an exact integer >= '+str(minimum)+' is required')
    return int(value)


def _rational_poly(poly):
    poly = s.sympify(poly)
    if poly.has(s.Float):
        raise ValueError('inexact coefficients are unsupported')
    return s.Poly(poly, X, Y, domain=s.QQ)


@lru_cache(None)
def _x_certificate(k):
    if k < 2:
        return X**k, s.S.Zero
    r, p = _x_certificate(k-2)
    if k > 2:
        r0, p0 = _x_certificate(k-3)
        r, p = r-h*s.Rational(k-2, 3)*r0, p-h*s.Rational(k-2, 3)*p0
    return s.expand(r), s.expand(p+X**(k-2)/3)


@lru_cache(None)
def _y_certificate(k):
    if k == 0:
        return s.S.One, s.S.Zero
    r, p = s.S.Zero, s.S.Zero
    if k > 1:
        r0, p0 = _y_certificate(k-2)
        r, p = -h*s.Rational(k-1, 2)*r0, -h*s.Rational(k-1, 2)*p0
    return s.expand(r), s.expand(p+Y**(k-1)/2)


def normal_form(poly):
    """Return exact r=a+b*X+Dx(P)+Dy(Q), with Dx=h*dX+qX."""
    poly = _rational_poly(poly)
    rem = px = py = s.S.Zero
    for (i, j), coefficient in poly.terms():
        rx, cx = _x_certificate(i)
        ry, cy = _y_certificate(j)
        rem += coefficient*rx*ry
        px += coefficient*cx*ry
        py += coefficient*X**i*cy
    rem, px, py = map(s.expand, (rem, px, py))
    return {'a': rem.coeff(X, 0), 'b': rem.coeff(X, 1), 'P': px, 'Q': py}


def stokes_residual(poly, certificate):
    """An independent polynomial reconstruction, expected to equal zero."""
    a, b, px, py = (certificate[k] for k in ('a', 'b', 'P', 'Q'))
    return s.expand(poly-a-b*X-h*s.diff(px, X)-s.diff(q, X)*px
                    -h*s.diff(py, Y)-s.diff(q, Y)*py)


def barycentric_amplitude(exponents):
    if len(exponents) != 3:
        raise ValueError('three barycentric exponents are required')
    a, b, c = [_integer(k) for k in exponents]
    return s.expand(u**a*v**b*w**c)


def _exact_algebraic(value):
    value = s.sympify(value)
    if value.has(s.Float) or value.free_symbols or value.is_algebraic is not True:
        raise ValueError('an exact algebraic scale is required')
    return value


def _valuation5(value):
    value = s.Rational(value)
    if not value:
        return None
    numerator, denominator = map(abs, (int(value.p), int(value.q)))
    result = 0
    while numerator % 5 == 0:
        numerator //= 5
        result += 1
    while denominator % 5 == 0:
        denominator //= 5
        result -= 1
    return result


def integer_shapes_certificate(scale):
    """Uniform in all three positive integer block shapes; no amplitude scan."""
    scale = _exact_algebraic(scale)
    if scale == 0:
        return {'status': 'nonzero_constant_phase_limit', 'limit': '1',
                'all_positive_integer_block_shapes': True}
    t = s.Symbol('t')
    minimal = s.Poly(s.minimal_polynomial(scale, t), t, domain=s.QQ).monic()
    valuations = [_valuation5(c) for c in minimal.all_coeffs()]
    accepted = any(k is not None and k < 0 for k in valuations)
    return {
        'status': ('nonzero_for_all_integer_block_shapes' if accepted
                   else 'outside_uniform_arithmetic_criterion'),
        'scale': str(scale), 'monic_minimal_polynomial': str(minimal.as_expr()),
        'coefficient_valuations_at_5': valuations,
        'critical_point': ['1', '0'], 'barycentric_values': ['1/4', '1/4', '1/2'],
        'prime': 5, 'uniform_mechanism': 'h=1/scale reduces to zero at a place above five',
        'period_is_zero': False if accepted else None,
        'theorem': 'FIRTWOANGSHAPES20261116' if accepted else None,
        'scope': 'Written all-exponent theorem; failed test makes no zero-period or FC decision.'}


def amplitude_certificate(poly, scale):
    """A rational amplitude at an algebraic scale; retain the actual boundary."""
    scale = _exact_algebraic(scale)
    if scale == 0:
        raise ValueError('use direct polynomial integration at scale zero')
    certificate = normal_form(poly)
    assert stokes_residual(poly, certificate) == 0
    field = s.QQ.algebraic_field(scale)
    values = [field.from_sympy(certificate[k].subs(h, 1/scale)) for k in ('a', 'b')]
    nonzero = any(value != field.zero for value in values)
    return {
        'status': ('nonzero_transcendental_period' if nonzero
                   else 'zero_remainder_no_period_decision'),
        'scale': str(scale),
        'remainder': [str(certificate[k]) for k in ('a', 'b')],
        'specialized_remainder': [str(field.to_sympy(value)) for value in values],
        'stokes_P': str(certificate['P']), 'stokes_Q': str(certificate['Q']),
        'identity_verified': True, 'boundary_term_retained': True,
        'period_is_zero': False if nonzero else None,
        'theorem': 'FIRTWOANGLINEAR20261116' if nonzero else None}


def weighted_certificate(shapes, scale):
    if len(shapes) != 3:
        raise ValueError('three block shape sums are required')
    shapes = [_integer(a, 1) for a in shapes]
    result = amplitude_certificate(barycentric_amplitude([a-1 for a in shapes]), scale)
    result['block_shapes'] = shapes
    return result


def family_certificate(variables, blocks, shapes, degree, scale, lower=0):
    """Construct a certified Gamma family. Shapes here are positive integers."""
    variables = list(variables)
    if not variables or len(set(variables)) != len(variables):
        raise ValueError('distinct source variables are required')
    shapes = [_integer(a, 1) for a in shapes]
    if len(shapes) != len(variables):
        raise ValueError('one Gamma shape per variable is required')
    degree = _integer(degree, 4)
    if len(blocks) != 3 or any(not block for block in blocks):
        raise ValueError('three nonempty blocks are required')
    blocks = [[_integer(j) for j in block] for block in blocks]
    if sorted(j for block in blocks for j in block) != list(range(len(variables))):
        raise ValueError('blocks must partition all source variable indices')
    scale = _exact_algebraic(scale)
    lower = s.sympify(lower)
    if lower.has(s.Float):
        raise ValueError('inexact lower coefficients are unsupported')
    lower_poly = s.Poly(lower, *variables, domain='EX')
    if lower_poly.total_degree() > degree-2:
        raise ValueError('the lower polynomial must have degree <=D-2')
    block_shapes = [sum(shapes[j] for j in block) for block in blocks]
    if scale == 0:
        period = integer_shapes_certificate(scale)
        accepted = True
    else:
        period = integer_shapes_certificate(scale)
        accepted = period['status'] == 'nonzero_for_all_integer_block_shapes'
        if not accepted:
            period = weighted_certificate(block_shapes, scale)
            accepted = period['status'] == 'nonzero_transcendental_period'
    b1, b2, b3 = [sum(variables[j] for j in block) for block in blocks]
    radial = sum(variables)
    cubic = 8*b3**3-6*b3*radial**2+(b2-b1)**2*radial
    phase = s.expand(radial**degree+degree*scale*radial**(degree-4)*cubic+lower)
    return {'status': ('eventually_every_pure_moment_nonzero' if accepted
                       else 'outside_nonzero_limit_criterion'),
            'phase': phase, 'block_shapes': block_shapes, 'degree': degree,
            'normalizer': '(sum shapes)_(D*m)', 'period_certificate': period,
            'first_nonzero_index': None, 'strong_factorial_window': None,
            'theorem': 'FIRTWOANGGAMMA20261116' if accepted else None}
