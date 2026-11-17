"""Exact Gamma linear-power moment windows, retaining polynomial markers.

The rational identities support every exponent. Window size is a proved
bound, not a claimed optimal zero-gap length. Individual rational shapes
may aggregate to integer shapes at equal nonzero nodes. Inexact and
transcendental arithmetic inputs are outside this exporter.
"""
import sympy as s


def positive_integer(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)) or value < 1:
        raise ValueError('an exact positive integer is required')
    return int(value)


def linear_power_window(nodes, shapes, power=1, marker=1, variables=None):
    nodes, shapes = list(nodes), list(shapes)
    if not nodes or len(nodes) != len(shapes):
        raise ValueError('one node and positive integer shape per variable required')
    shapes = list(map(s.sympify, shapes))
    if any(a.has(s.Float) or a.is_Rational is not True or a <= 0 for a in shapes):
        raise ValueError('exact positive rational shapes are required')
    power = positive_integer(power)
    nodes = list(map(s.sympify, nodes))
    variables = list(s.symbols('Y:'+str(len(nodes)))) if variables is None else list(variables)
    if len(variables) != len(nodes) or len(set(variables)) != len(nodes):
        raise ValueError('distinct original marker variables required')
    marker = s.sympify(marker)
    if marker.has(s.Float) or any(c.has(s.Float) for c in nodes):
        raise ValueError('inexact coefficients are unsupported')
    marked = s.Poly(marker, *variables, domain='EX')
    coefficients = nodes+[c for _, c in marked.terms()]
    if any(c.free_symbols or c.is_algebraic is not True for c in coefficients):
        raise ValueError('exact algebraic nodes and marker coefficients required')
    extensions = [c for c in coefficients if c.is_Rational is not True]
    field = s.QQ.algebraic_field(*extensions) if extensions else s.QQ
    nodes = [field.to_sympy(field.from_sympy(c)) for c in nodes]
    z, t = s.Dummy('z'), s.Dummy('t')
    grouped = {}
    for i, node in enumerate(nodes):
        key = field.from_sympy(node)
        if key != field.zero:
            grouped.setdefault(key, []).append(i)
    groups = [(field.to_sympy(key), indices, sum(shapes[i] for i in indices))
              for key, indices in grouped.items()]
    if any(a.is_Integer is not True for _, _, a in groups):
        raise ValueError('shape sums at each distinct nonzero node must be integers')
    mass = s.expand(sum(c*s.prod(s.rf(a, k) for a, k in zip(shapes, exponents))
                        for exponents, c in marked.terms()))
    if not groups:
        return {'status': 'zero_phase', 'numerator': mass, 'denominator': s.S.One,
                'variable': t, 'window': None, 'all_positive_moments_zero': True,
                'constant_marked_moment': mass}
    if marked.is_zero:
        return {'status': 'identically_zero_marked_sequence', 'numerator': s.S.Zero,
                'denominator': s.S.One, 'variable': t, 'window': None,
                'all_marked_moments_zero_including_constant': True,
                'theorem': 'FIRMARKEDWINDOW20261116'}
    degrees = [int(marked.degree(variable)) for variable in variables]
    multiplicities = [(node, int(a)+sum(degrees[i] for i in indices)) for node, indices, a in groups]
    denominator = s.Poly(s.prod((1-node*z)**m for node, m in multiplicities), z, domain=field)
    sampled_denominator = s.Poly(s.prod((1-node**power*t)**m for node, m in multiplicities), t, domain=field)
    rational = sum(c*s.prod(s.rf(a, k) for a, k in zip(shapes, exponents)) /
                   s.prod((1-node*z)**(int(a)+sum(exponents[i] for i in indices))
                          for node, indices, a in groups)
                   for exponents, c in marked.terms())
    num, den = s.fraction(s.together(rational))
    raw_num, raw_den = s.Poly(num, z, domain=field), s.Poly(den, z, domain=field)
    # Multiply by the known denominator using exact polynomial division.
    common_num = raw_num*denominator.exquo(raw_den)
    quotient = s.Poly(sampled_denominator.as_expr().subs(t, z**power), z, domain=field).exquo(denominator)
    section_source = common_num*quotient
    sampled_num = s.Poly.from_dict({(exponents[0]//power,): c
                                   for exponents, c in section_source.terms()
                                   if exponents[0] % power == 0}, t, domain=field)
    if sampled_num.is_zero:
        return {'status': 'identically_zero_marked_sequence', 'numerator': s.S.Zero,
                'denominator': s.S.One, 'variable': t, 'window': None,
                'all_marked_moments_zero_including_constant': True,
                'theorem': 'FIRMARKEDWINDOW20261116'}
    gcd = s.gcd(sampled_num, sampled_denominator)
    sampled_num, sampled_denominator = sampled_num.exquo(gcd), sampled_denominator.exquo(gcd)
    normalization = sampled_denominator.eval(0)
    sampled_num, sampled_denominator = sampled_num.mul_ground(1/normalization), sampled_denominator.mul_ground(1/normalization)
    assert sampled_num.degree() < sampled_denominator.degree()
    return {'status': 'nonzero_in_every_certified_window',
            'numerator': sampled_num.as_expr(), 'denominator': sampled_denominator.as_expr(),
            'variable': t, 'window': sampled_denominator.degree(),
            'window_start': 'every integer n>=1',
            'normalized_coefficient': 'E[marker*L^(power*m)]/(power*m)!',
            'uncancelled_denominator_degree': denominator.degree(),
            'optimal_window_claim': False,
            'identity_scope': 'all exponents; exact rational section identity',
            'theorem': 'FIRMARKEDWINDOW20261116'}


def coefficients(certificate, through):
    if not isinstance(through, int) or through < 0:
        raise ValueError('nonnegative integer truncation required')
    t = certificate['variable']
    numerator = s.Poly(certificate['numerator'], t)
    denominator = s.Poly(certificate['denominator'], t)
    result = []
    for k in range(through+1):
        result.append(s.expand(numerator.nth(k)-sum(denominator.nth(j)*result[k-j]
                                                    for j in range(1, min(k, denominator.degree())+1))))
    return result


def export(certificate):
    """JSON-safe rational identity with a stable parameter name."""
    parameter = s.Symbol('t')
    return {k: (str(value.subs(certificate['variable'], parameter)) if isinstance(value, s.Basic) else value)
            for k, value in certificate.items() if k != 'variable'}
