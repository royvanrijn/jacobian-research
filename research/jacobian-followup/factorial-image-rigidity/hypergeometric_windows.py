"""Exact monomial/marker windows. Universal proofs are in the companion TeX.

The executable domain is positive rational shapes, exact algebraic markers.
The mathematics allows every positive real shape. A degree-dependent window
is not a two-moment SFC certificate unless an actual support budget suffices.
"""
from dataclasses import dataclass
from math import comb, prod
import sympy as s


def _integer(value, minimum=0):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)) or value < minimum:
        raise ValueError('integer >= '+str(minimum)+' required')
    return int(value)


@dataclass(frozen=True)
class Model:
    law: str
    exponents: tuple
    shapes: tuple
    variables: tuple
    marker: s.Expr
    homogenized_marker: s.Expr
    homogenization_degree: int
    t: s.Symbol
    Q: s.Poly
    S: s.Poly
    R: s.Poly

    @property
    def degree(self):
        return sum(self.exponents)

    @property
    def marker_degree(self):
        return None if self.R.is_zero else int(self.R.degree())

    @property
    def beta(self):
        if self.law != 'dirichlet' or self.R.is_zero:
            return None
        inactive = sum(a for a, d in zip(self.shapes, self.exponents) if not d)
        active = sum(d > 0 for d in self.exponents)
        return inactive+self.homogenization_degree-self.marker_degree+s.Rational(active-1, 2)


def make_model(exponents, shapes, marker=1, variables=None, law='gamma'):
    exponents = tuple(_integer(d) for d in exponents)
    if not exponents or not sum(exponents):
        raise ValueError('a nonconstant monomial is required')
    shapes = tuple(s.sympify(a) for a in shapes)
    if len(shapes) != len(exponents) or any(a.is_Rational is not True or a <= 0 for a in shapes):
        raise ValueError('one positive exact rational shape per variable required')
    if law not in ('gamma', 'dirichlet'):
        raise ValueError('law must be gamma or dirichlet')
    if law == 'dirichlet' and len(shapes) < 2:
        raise ValueError('the one-coordinate Dirichlet point law is excluded')
    variables = tuple(s.symbols('y:'+str(len(shapes)))) if variables is None else tuple(variables)
    if len(variables) != len(shapes) or len(set(variables)) != len(variables) or any(not isinstance(v, s.Symbol) for v in variables):
        raise ValueError('distinct variable symbols required')
    marker = s.sympify(marker)
    if marker.has(s.Float) or marker.free_symbols-set(variables):
        raise ValueError('exact algebraic marker coefficients required')
    marker_poly = s.Poly(marker, *variables, extension=True)
    if any(c.is_algebraic is not True for c in marker_poly.coeffs()):
        raise ValueError('exact algebraic marker coefficients required')
    E = 0 if marker_poly.is_zero else int(marker_poly.total_degree())
    homogeneous = marker
    if law == 'dirichlet':
        homogeneous = s.expand(sum(c*prod(y**j for y, j in zip(variables, powers))
                                  *sum(variables)**(E-sum(powers))
                                  for powers, c in marker_poly.terms()))
    t = s.Dummy('t')
    R = s.Poly(s.expand(sum(c*prod(s.rf(a+d*t, j) for a, d, j in zip(shapes, exponents, powers))
                           for powers, c in s.Poly(homogeneous, *variables, extension=True).terms())), t, extension=True)
    Q = s.Poly(prod(a+d*t+j for a, d in zip(shapes, exponents) for j in range(d)), t)
    S = s.Poly(s.rf(sum(shapes)+E+sum(exponents)*t, sum(exponents)) if law == 'dirichlet' else 1, t)
    return Model(law, exponents, shapes, variables, marker, homogeneous, E, t, Q, S, R)


def coefficient(model, k):
    k = _integer(k)
    base = prod(s.rf(a, d*k) for a, d in zip(model.shapes, model.exponents))
    if model.law == 'dirichlet':
        base /= s.rf(sum(model.shapes), model.degree*k+model.homogenization_degree)
    return s.simplify(base*model.R.eval(k))


def moment(model, n, constant=1, coefficient_value=1):
    n = _integer(n)
    constant, coefficient_value = s.sympify(constant), s.sympify(coefficient_value)
    if constant.has(s.Float) or coefficient_value.has(s.Float):
        raise ValueError('exact phase coefficients required')
    return s.expand(sum(s.binomial(n, k)*constant**(n-k)*coefficient_value**k*coefficient(model, k)
                        for k in range(n+1)))


def moment_polynomial(model, n, variable=None):
    variable = s.Symbol('x') if variable is None else variable
    return s.Poly(moment(model, n, 1, variable), variable, extension=True)


def operator_polynomials(model, n):
    """Return P(t),B(t) for [P(theta)-x B(theta)]g_N=0; x acts last."""
    t = model.t
    R, Q, S = model.R.as_expr(), model.Q.as_expr(), model.S.as_expr()
    return (s.Poly(t*S.subs(t, t-1)*R.subs(t, t-1), t, extension=True),
            s.Poly((n-t)*Q*R.subs(t, t+1), t, extension=True))


def window_certificate(model, start=1):
    start = _integer(start, 1)
    base = dict(law=model.law, shapes=[str(a) for a in model.shapes], exponents=list(model.exponents),
                monomial_total_degree=model.degree, marker=str(model.marker),
                homogenization_degree=model.homogenization_degree if model.law == 'dirichlet' else None,
                literal_marker_polynomial=str(model.R.as_expr()),
                Q_coefficients=[str(c) for c in reversed(model.Q.all_coeffs())],
                S_coefficients=[str(c) for c in reversed(model.S.all_coeffs())],
                executable_shape_domain='positive rational', theorem_shape_domain='positive real',
                phase_requirement='a+b*Z with b!=0; all complex a,b', all_starting_indices=True,
                two_term_sfc_claimed=False, counterexample_claimed=False)
    if model.R.is_zero:
        return dict(base, status='marked_zero_sequence_identity', includes_index_zero=True,
                    window_length=None, marker_degree=None)
    L = model.degree+model.marker_degree+1
    N = s.Symbol('N', integer=True, nonnegative=True)
    P, B = operator_polynomials(model, N)
    result = dict(base, status='uniform_moment_window_certified', window_length=L,
                  marker_degree=model.marker_degree, start=start, terminal=start+L-1,
                  theta_operator_P=[str(c) for c in reversed(P.all_coeffs())],
                  theta_operator_B=[str(c) for c in reversed(B.all_coeffs())],
                  constant_phase_marker_mass=str(coefficient(model, 0)),
                  constant_phase_rule='b=0: moments equal a^n times the marker mass',
                  singular_parameters=['0'])
    if model.law == 'dirichlet':
        assert model.beta > 0
        result.update(singular_parameters=['0', str(-s.Rational(model.degree**model.degree, model.Q.LC()))],
                      extra_singular_indicial_exponent='N + '+str(model.beta),
                      extra_singular_beta=str(model.beta), extra_singular_beta_positive=True,
                      extra_singular_exponent_exceeds_polynomial_degree=True)
    return result


def ordinary_ideal_certificate(model, start=1):
    """Exhibit a regular SUBIDEAL, without asserting the window ideal is CI."""
    start = _integer(start, 1)
    if s.expand(model.marker-1) != 0:
        raise ValueError('this coefficient-ideal exporter is for unmarked moments')
    D, x = model.degree, s.Symbol('x')
    terminal = start+D
    f = moment_polynomial(model, start, x)
    candidates = [moment_polynomial(model, start+j, x) for j in range(1, D+1)]
    for value in range(start*(D-1)+1):
        G = s.Poly(sum(value**(j-1)*candidates[j-1].as_expr() for j in range(1, D+1)), x, domain=s.QQ)
        gcd = s.gcd(f, G)
        if gcd.degree() == 0:
            return dict(status='ordinary_window_ideal_bound_certified', start=start, terminal=terminal,
                        selected_integer=value, selection_integer_upper_bound=start*(D-1),
                        regular_subideal_degrees=[start, terminal],
                        dehomogenized_first=str(f.as_expr()), dehomogenized_second=str(G.as_expr()),
                        gcd_is_one=True, projective_infinity_excluded_by_first=True,
                        regular_subideal_length=start*terminal,
                        window_ideal_length_upper_bound=start*terminal,
                        window_nilpotence_exponent_upper_bound=2*start+D-1,
                        all_Q_algebras=True, ordinary_not_radical=True,
                        window_ideal_is_complete_intersection_claimed=False, sharpness_claimed=False)
    raise AssertionError('the proved finite selection bound failed')


def block_sum_sfc_certificate(block_sizes, exponents, extra_linear=0):
    sizes = tuple(_integer(r, 1) for r in block_sizes)
    exponents = tuple(_integer(d, 1) for d in exponents)
    extra_linear = _integer(extra_linear)
    if not sizes or len(sizes) != len(exponents):
        raise ValueError('one nonempty block per positive exponent required')
    support = 1+prod(comb(r+d-1, d) for r, d in zip(sizes, exponents))
    window = 1+sum(exponents)
    accepted = support >= window
    return dict(status='expanded_block_sum_sfc_certified' if accepted else 'support_budget_insufficient_no_decision',
                block_sizes=list(sizes), exponents=list(exponents), extra_independent_linear_terms=extra_linear,
                actual_expanded_support_when_a_b_nonzero=support+extra_linear,
                proved_window_length=window+extra_linear,
                all_starting_indices=True, applies_to_factorial_moments=True,
                base_Gaussian_moment_analogue=True,
                Gaussian_linear_adjunction_claimed=False,
                full_two_term_sfc_claimed=accepted and support == 2,
                powers_of_entire_phase_claimed=False, counterexample_claimed=False,
                degenerate_a_or_b_zero='base moments nonzero directly; factorial linear adjunctions still apply')
