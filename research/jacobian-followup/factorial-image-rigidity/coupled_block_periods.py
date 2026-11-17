"""Exact coupled cube periods and Gamma block-limit identities.

Arithmetic decisions require exact algebraic scales/markers. Block limits
accept rational positive shapes; their written proof allows all positive reals.
A zero leading period does not certify a whole Gamma marked zero sequence.
"""
from math import prod
import sympy as s
from sympy.functions.combinatorial.numbers import stirling


def _integer(value, minimum=0):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)) or value < minimum:
        raise ValueError('integer >= '+str(minimum)+' required')
    return int(value)


def _variables(variables):
    variables = tuple(variables)
    if not variables or len(set(variables)) != len(variables) or any(not isinstance(v, s.Symbol) for v in variables):
        raise ValueError('nonempty distinct variable symbols required')
    return variables


def _polynomial(value, variables):
    value = s.sympify(value)
    if value.has(s.Float) or value.free_symbols-set(variables):
        raise ValueError('exact algebraic polynomial coefficients required')
    polynomial = s.Poly(value, *variables, extension=True)
    if any(c.is_algebraic is not True for c in polynomial.coeffs()):
        raise ValueError('exact algebraic polynomial coefficients required')
    return polynomial


def _homogeneous(polynomial, degree):
    return s.expand(sum(c*prod(v**a for v, a in zip(polynomial.gens, powers))
                        for powers, c in polynomial.terms() if sum(powers) == degree))


def cube_reduction(marker, variables, parameter=None):
    """Return literal W(t), its verified partial fractions and the E-basis vector.

    Basis order is 1,exp(z),J_1(z),...,J_(r-1)(z). All coefficients
    are Laurent polynomials. Internal SymPy values are retained for replay.
    """
    variables = _variables(variables)
    polynomial = _polynomial(marker, variables)
    r = len(variables)
    t, z = s.Dummy('t'), s.Dummy('z') if parameter is None else parameter
    W = s.cancel(sum(c/prod(t+a+1 for a in powers) for powers, c in polynomial.terms()))
    capacities = {}
    for powers, c in polynomial.terms():
        if c == 0:
            continue
        for shift in set(a+1 for a in powers):
            capacities[shift] = max(capacities.get(shift, 0), sum(a+1 == shift for a in powers))
    partial, reconstructed = [], s.Integer(0)
    coefficients = [s.Integer(0)]*(r+1)
    for shift, capacity in sorted(capacities.items()):
        regular = s.cancel((t+shift)**capacity*W)
        for height in range(1, capacity+1):
            order = capacity-height
            amplitude = s.simplify(s.diff(regular, t, order).subs(t, -shift)/s.factorial(order))
            if amplitude == 0:
                continue
            partial.append(dict(shift=shift, height=height, amplitude=amplitude))
            reconstructed += amplitude/(t+shift)**height
            for ell in range(1, shift+1):
                factor = amplitude*stirling(shift, ell, kind=1, signed=True)/z**shift
                index = height-ell
                if index > 0:
                    assert index <= r-1
                    coefficients[index+1] += factor
                elif index == 0:
                    coefficients[0] -= factor
                    coefficients[1] += factor
                else:
                    touchard = sum(stirling(-index, j, kind=2)*z**j for j in range(-index+1))
                    coefficients[1] += factor*touchard
    assert s.cancel(W-reconstructed) == 0, 'partial fractions failed reconstruction'
    return dict(dimension=r, marker=polynomial.as_expr(), moment_parameter=t,
                rational_moment=W, parameter=z, partial_fractions=partial,
                basis=['1', 'exp(z)']+['J_'+str(j)+'(z)' for j in range(1, r)],
                coefficients=[s.cancel(c) for c in coefficients])


def period_certificate(marker, variables, scale):
    variables = _variables(variables)
    scale = s.simplify(s.sympify(scale))
    if scale.has(s.Float) or scale.is_algebraic is not True:
        raise ValueError('exact algebraic scale required')
    reduction = cube_reduction(marker, variables)
    z, t = reduction['parameter'], reduction['moment_parameter']
    marker_coefficients = _polynomial(marker, tuple(variables)).coeffs()
    generators = [a for a in marker_coefficients+[scale] if not a.is_Rational]
    field = s.QQ.algebraic_field(*generators) if generators else s.QQ
    canonical = lambda value: field.to_sympy(field.from_sympy(s.simplify(value)))
    if scale == 0:
        values = [canonical(reduction['rational_moment'].subs(t, 0))]+[s.Integer(0)]*reduction['dimension']
    else:
        values = [canonical(c.subs(z, scale)) for c in reduction['coefficients']]
    zero = all(c == 0 for c in values)
    transcendental = any(c != 0 for c in values[1:])
    return dict(status='zero_period_certified' if zero else 'nonzero_period_certified',
                scale=str(scale), dimension=reduction['dimension'], marker=str(reduction['marker']),
                rational_moment=str(reduction['rational_moment']),
                partial_fractions=[dict(shift=p['shift'], height=p['height'], amplitude=str(p['amplitude']))
                                   for p in reduction['partial_fractions']],
                basis=reduction['basis'], coefficients=[str(c) for c in reduction['coefficients']],
                specialized_coefficients=[str(c) for c in values],
                nonzero_ordinary_point_value_independence_used=scale != 0,
                transcendence_certified=transcendental,
                algebraic_value=None if transcendental else str(values[0]),
                cube_marked_power_sequence_identically_zero=reduction['rational_moment'] == 0,
                full_Gamma_marked_sequence_zero_claimed=False,
                zero_period_alone_decides_whole_sequence=False,
                theorem='FIRCOUPLEDPERIOD20261116')


def block_limit_data(phase, variables, blocks, exponents, shapes, marker=1):
    """Verify a block top on an active subset; integrate inactive marker factors."""
    variables = _variables(variables)
    blocks = tuple(tuple(_integer(i) for i in block) for block in blocks)
    exponents = tuple(_integer(d, 1) for d in exponents)
    active = tuple(i for block in blocks for i in block)
    if (not blocks or len(blocks) != len(exponents) or any(not block for block in blocks)
            or len(set(active)) != len(active) or any(i >= len(variables) for i in active)):
        raise ValueError('nonempty disjoint blocks of actual variable indices required')
    inactive = tuple(i for i in range(len(variables)) if i not in active)
    shapes = tuple(s.sympify(a) for a in shapes)
    if len(shapes) != len(variables) or any(a.is_Rational is not True or a <= 0 for a in shapes):
        raise ValueError('one positive rational Gamma shape per variable required')
    polynomial, H = _polynomial(phase, variables), _polynomial(marker, variables)
    D = sum(exponents)
    block_sums = [sum(variables[i] for i in block) for block in blocks]
    F0 = s.expand(prod(S**d for S, d in zip(block_sums, exponents)))
    top = _homogeneous(polynomial, D)
    corner = prod(variables[block[0]]**d for block, d in zip(blocks, exponents))
    c = polynomial.coeff_monomial(corner)
    if polynomial.total_degree() != D or c == 0 or s.expand(top-c*F0) != 0:
        return dict(status='outside_block_radial_top_no_decision', FC_decision=None)
    G = _homogeneous(polynomial, D-1)
    remainder = s.expand(polynomial.as_expr()-c*F0-G)
    degree = 0 if H.is_zero else max(sum(powers[i] for i in active) for powers, _ in H.terms())
    Htop = s.expand(sum(coeff*prod(v**a for v, a in zip(variables, powers))
                         for powers, coeff in H.terms() if sum(powers[i] for i in active) == degree))
    angular = tuple(s.Dummy('T'+str(i)) for i in range(len(variables)))
    substitution = {variables[i]: d*angular[i] for block, d in zip(blocks, exponents) for i in block}
    restricted_G = G.subs({variables[i]: 0 for i in inactive}, simultaneous=True)
    q = s.expand(restricted_G.subs(substitution, simultaneous=True)/(c*prod(d**d for d in exponents)))
    raw_marker = s.expand(Htop.subs(substitution, simultaneous=True))
    angular_marker = raw_marker
    if inactive:
        angular_marker = s.expand(sum(coeff*prod(s.rf(shapes[i], a) for i, a in zip(inactive, powers))
                                    for powers, coeff in s.Poly(raw_marker, *(variables[i] for i in inactive)).terms()))
    return dict(status='block_limit_identity_certified', variables=variables, blocks=blocks,
                exponents=exponents, shapes=shapes, active_indices=active, inactive_indices=inactive,
                block_shapes=tuple(sum(shapes[i] for i in block) for block in blocks),
                phase=polynomial.as_expr(), top=top, scalar=c, next_layer=G, remainder=remainder,
                marker_degree=degree, marker_top_active=Htop, marker_is_zero=H.is_zero, angular_variables=angular,
                angular_marker_with_inactive=raw_marker, restricted_next_layer=restricted_G,
                angular_phase=q, angular_marker=angular_marker,
                unit_nonvanishing_decision=None)


def block_limit_certificate(phase, variables, blocks, exponents, shapes, marker=1):
    data = block_limit_data(phase, variables, blocks, exponents, shapes, marker)
    if data['status'] != 'block_limit_identity_certified':
        return data
    return dict(status=data['status'], phase=str(data['phase']), scalar=str(data['scalar']),
                block_indices=[list(b) for b in data['blocks']], exponents=list(data['exponents']),
                active_indices=list(data['active_indices']), inactive_indices=list(data['inactive_indices']),
                shapes=[str(a) for a in data['shapes']], block_shapes=[str(a) for a in data['block_shapes']],
                next_layer=str(data['next_layer']), remainder=str(data['remainder']),
                angular_phase=str(data['angular_phase']), angular_marker=str(data['angular_marker']),
                marker_top_active=str(data['marker_top_active']),
                angular_marker_with_inactive=str(data['angular_marker_with_inactive']),
                inactive_expectation_taken=True, marker_degree_convention='active variables only',
                restricted_next_layer=str(data['restricted_next_layer']),
                angular_variables=[str(t) for t in data['angular_variables']],
                angular_constraints=[' + '.join(str(data['angular_variables'][i]) for i in block)+' = 1'
                                     for block in data['blocks']],
                marker_degree=data['marker_degree'], marker_is_zero=data['marker_is_zero'],
                normalization='m^marker_degree * scalar^m * product_j (block_shape_j)_(d_j*m)',
                theorem_shape_domain='all positive real', executable_shape_domain='positive rational',
                nonzero_leading_period_implies_eventual_nonzero=True,
                zero_leading_period_implies_zero_sequence=False, unit_nonvanishing_decision=None)


def monomial_top_certificate(phase, variables, exponents, shapes, marker=1):
    variables = _variables(variables)
    exponents = tuple(_integer(d) for d in exponents)
    if len(exponents) != len(variables) or not any(exponents):
        raise ValueError('one nonnegative exponent per variable, with at least one positive, required')
    active = tuple(i for i, d in enumerate(exponents) if d)
    data = block_limit_data(phase, variables, tuple((i,) for i in active),
                            tuple(exponents[i] for i in active), shapes, marker)
    if data['status'] != 'block_limit_identity_certified':
        return data
    substitution = {t: 1 for t in data['angular_variables']}
    q = s.simplify(data['angular_phase'].subs(substitution))
    amplitude = s.simplify(data['angular_marker'].subs(substitution))
    field = s.QQ if amplitude.is_Rational else s.QQ.algebraic_field(amplitude)
    nonzero = field.from_sympy(amplitude) != field.zero
    return dict(status='single_top_monomial_eventual_nonzero_certified', phase=str(data['phase']),
                exponents=list(exponents), shapes=[str(a) for a in data['shapes']],
                active_indices=list(active), inactive_indices=list(data['inactive_indices']),
                unit_limit_exponent=str(q), leading_marker_amplitude=str(amplitude),
                marker_degree=data['marker_degree'], marker_degree_convention='active variables only',
                marker_top_active=str(data['marker_top_active']), inactive_expectation_taken=True,
                unit_moments_eventually_all_nonzero=True,
                marked_eventual_nonzero_certified=nonzero,
                marked_decision='eventually_nonzero' if nonzero else 'zero_leading_amplitude_no_decision',
                zero_leading_amplitude_decides_whole_sequence=False,
                every_variable_occurs_in_top=not data['inactive_indices'],
                all_complex_lower_coefficients_in_theorem=True,
                all_positive_real_shapes_in_theorem=True,
                ordinary_universal_tail_moment_ideal_is_unit=True, SFC_claimed=False,
                theorem='FIRMONOMIALTOPFC20261116')


def coupled_family_certificate(exponents, scale, lower=0, marker=1, variables=None,
                               inactive_variables=(), inactive_shapes=None, coupling=0):
    exponents = tuple(_integer(d, 1) for d in exponents)
    if not exponents or exponents[0] < 2:
        raise ValueError('first block exponent >=2 and all other exponents >=1 required')
    r = len(exponents)
    variables = tuple(s.symbols('y:'+str(2*r))) if variables is None else _variables(variables)
    if len(variables) != 2*r:
        raise ValueError('two factorial variables per uniform block required')
    inactive_variables = tuple(inactive_variables)
    all_variables = _variables(variables+inactive_variables)
    inactive_shapes = ((s.Integer(1),)*len(inactive_variables) if inactive_shapes is None
                       else tuple(s.sympify(a) for a in inactive_shapes))
    if len(inactive_shapes) != len(inactive_variables):
        raise ValueError('one inactive shape per inactive variable required')
    scale = s.sympify(scale)
    if scale.has(s.Float) or scale.is_algebraic is not True:
        raise ValueError('the executable period certificate requires an exact algebraic scale')
    lower_poly = _polynomial(lower, all_variables)
    coupling_poly = _polynomial(coupling, all_variables)
    D = sum(exponents)
    if not lower_poly.is_zero and lower_poly.total_degree() > D-2:
        raise ValueError('lower remainder exceeds degree D-2')
    if (s.expand(_homogeneous(coupling_poly, D-1)-coupling_poly.as_expr()) != 0
            or s.expand(coupling_poly.as_expr().subs({v: 0 for v in inactive_variables})) != 0):
        raise ValueError('coupling must be homogeneous of degree D-1 and vanish at inactive=0')
    blocks = tuple((2*j, 2*j+1) for j in range(r))
    sums = [variables[2*j]+variables[2*j+1] for j in range(r)]
    F0 = prod(S**d for S, d in zip(sums, exponents))
    G = exponents[0]*scale*sums[0]**(exponents[0]-2)*prod(S**(d-1) for S, d in zip(sums[1:], exponents[1:]))*prod(variables[2*j] for j in range(r))
    phase = s.expand(F0+G+coupling_poly.as_expr()+lower_poly.as_expr())
    data = block_limit_data(phase, all_variables, blocks, exponents, (1,)*(2*r)+inactive_shapes, marker)
    cube = tuple(s.Dummy('U'+str(j)) for j in range(r))
    substitution = {data['angular_variables'][2*j+i]: cube[j] if i == 0 else 1-cube[j]
                    for j in range(r) for i in range(2)}
    q = s.expand(data['angular_phase'].subs(substitution, simultaneous=True))
    P = s.expand(data['angular_marker'].subs(substitution, simultaneous=True))
    assert s.expand(q-scale*prod(cube)) == 0
    unit = period_certificate(1, cube, scale)
    marked = period_certificate(P, cube, scale)
    assert unit['status'] == 'nonzero_period_certified'
    gaussian_inactive = [2*a for a in inactive_shapes]
    gaussian_count = 4*r+sum(int(a) for a in gaussian_inactive) if all(a.is_Integer for a in gaussian_inactive) else None
    return dict(status='coupled_factorial_family_certified', phase=str(phase), variables=[str(v) for v in all_variables],
                inactive_variables=[str(v) for v in inactive_variables],
                inactive_shapes=[str(a) for a in inactive_shapes], coupling=str(coupling_poly.as_expr()),
                marker_degree_convention='active variables only', inactive_expectation_taken=True,
                exponents=list(exponents), scale=str(scale), angular_phase=str(q),
                marker_degree=data['marker_degree'], leading_cube_marker=str(P),
                exact_zero_input_marker=data['marker_is_zero'],
                whole_marked_sequence_zero_certified=data['marker_is_zero'],
                unit_period=unit, leading_marker_period=marked,
                unit_moments_eventually_all_nonzero=True,
                marked_eventual_nonzero_certified=marked['status'] == 'nonzero_period_certified',
                marked_decision='eventually_nonzero' if marked['status'] == 'nonzero_period_certified'
                                else ('identically_zero_input_marker' if data['marker_is_zero'] else 'zero_leading_period_no_decision'),
                zero_leading_marker_period_decides_full_sequence=False,
                arbitrary_complex_lower_terms_in_theorem=True,
                all_complex_scale_FC_and_infinitely_many_moments_theorem=True,
                ordinary_universal_tail_moment_ideal_is_unit=True,
                all_Q_algebra_scalar_normalization_ideal_is_zero=True,
                effective_tail_cutoff_claimed=False,
                all_complex_scale_nonzero_limit_claimed=False, SFC_claimed=False,
                Gaussian_coordinate_count=gaussian_count, active_Gaussian_coordinate_count=4*r,
                theorem='FIRCOUPLEDBLOCKFC20261116')
