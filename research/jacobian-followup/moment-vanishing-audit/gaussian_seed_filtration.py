"""Exact marker formulas for determinant compensators and pure power delays.

Universal assertions are proved in POLYNOMIAL_SEEDS_AND_PURE_MARKER_DELAY.md.
This module performs finite rational symbolic arithmetic, not a safety oracle.
"""
from functools import lru_cache
import sympy as s

t, u, z, w, y, v, a, T = s.symbols('t u z w y v a T')
D, gamma, theta = s.symbols('D gamma theta')


def nonnegative(value, name):
    if type(value) is not int or value < 0:
        raise ValueError(name + ' must be a nonnegative integer')


def compensator(M, variable=y):
    nonnegative(M, 'M')
    N = 2*M+4
    return s.Add(*(s.binomial(N, j+1)*variable**j for j in range(N)))


def phase(f, M, *, five_variables=False):
    """Return the polynomial in original variables; f is a polynomial in t."""
    nonnegative(M, 'M')
    f = s.Poly(f, t).as_expr()
    if f.subs(t, 0) != 1:
        raise ValueError('f(0) must equal one')
    B = s.cancel((1-f)/t)
    if five_variables:
        if t in B.free_symbols:
            raise ValueError('Five-variable construction requires constant B')
        return B*(w*z + v*(1+y) - (1+y)*compensator(M)*T**2/2)
    B = B.subs(t, u)
    return a + B*w*z + B*v*(1+y) - B*(1+y)*compensator(M)*T**2/2


@lru_cache(None)
def momentum_polynomial(M, c, d, j):
    """d-th gamma derivative at zero, with the half power expanded finitely."""
    N = 2*M+4
    R = s.expand(D**N + (1-D)*(1+gamma)*sum(
        (1+gamma)**(N-1-i)*D**i for i in range(N)))
    assert s.expand(R.subs(gamma, 0)) == 1
    # Work modulo gamma^(d+1); R-1 has no constant term.
    delta = s.Poly(R-1, gamma)
    power = s.Poly(1, gamma)
    total = s.Poly(0, gamma)
    for k in range(d+1):
        total += s.binomial(-s.Rational(2*j+1, 2), k)*power
        power = s.Poly(s.Add(*(coef*gamma**ex[0] for ex, coef in
            (power*delta).terms() if ex[0] <= d)), gamma)
    result = s.Poly((1-D+gamma)**c*total.as_expr(), gamma).nth(d)*s.factorial(d)
    return s.factor(result)


def marker_period(f, M, exponents):
    """Period of U^p A^q Z^az W^bw Y^c V^d T^h in the seven-variable phase.

    For the five-variable phase use p=q=0. Time coefficients are external.
    """
    nonnegative(M, 'M')
    if len(exponents) != 7:
        raise ValueError('Seven original exponents required')
    for e in exponents:
        nonnegative(e, 'exponent')
    p, q, az, bw, c, d, h = exponents
    f = s.Poly(f, t).as_expr()
    if f.subs(t, 0) != 1:
        raise ValueError('f(0) must equal one')
    if az != bw or h % 2:
        return s.Integer(0)
    j = h//2
    B = s.cancel((1-f)/t).subs(t, u)
    Du = 1-t*B
    inner = momentum_polynomial(M, c, d, j)
    value = u**p*Du**(M-az-c+(2*M+4)*j)*inner.subs(D, Du)
    return s.factor(s.factorial(az)*s.factorial2(2*j-1)*
                    s.diff(value, u, q).subs(u, t))


def filtration_generator(f, M, d):
    nonnegative(M, 'M')
    nonnegative(d, 'd')
    return f**(M-d)


def power_transfer(old_period, K):
    """Extract every K-th exponential coefficient: polynomial + EulerOp*h_K.

    old_period is rational in t with only possible pole at t=1.
    """
    if type(K) is not int or K < 1:
        raise ValueError('K must be positive')
    old_period = s.cancel(old_period)
    _, den = s.fraction(old_period)
    pd = s.Poly(den, t)
    if s.cancel(den/(1-t)**pd.degree()).has(t):
        raise ValueError('Only a pole at one is supported')
    num = s.fraction(old_period)[0]
    quotient, remainder = s.div(num, den, t)
    polynomial = sum(s.Poly(quotient, t).nth(K*n)*s.factorial(K*n)/s.factorial(n)*t**n
                     for n in range(s.degree(quotient, t)//K+1)) if quotient else s.Integer(0)
    pole_part = s.cancel(remainder/den)
    h_op = s.Integer(0)
    if pole_part:
        X = s.Symbol('X')
        laurent = s.expand(pole_part.subs(t, 1-X))
        for term in s.Add.make_args(laurent):
            coeff, exponent = term.as_coeff_exponent(X)
            if coeff.has(X) or not exponent.is_Integer or exponent >= 0:
                raise ValueError('Invalid pole decomposition')
            p = -int(exponent)
            h_op += coeff*s.prod(K*theta+j for j in range(1, p))/s.factorial(p-1)
    return s.expand(polynomial), s.expand(h_op)


def power_data(M, K):
    nonnegative(M, 'M')
    if type(K) is not int or K <= M:
        raise ValueError('Pure powers require K > M')
    return dict(M=M, K=K, gaussian_dimension=5, first_bad_marker_degree=M+1,
                bad_marker=(1+y)**(M+1), unit_period=s.Integer(1),
                actual_rank=K, quotient_operator=quotient_operator(K),
                phase=phase(1-t, M, five_variables=True)**K)


def h_coefficients(K, N):
    if type(K) is not int or K < 1:
        raise ValueError('K must be positive')
    nonnegative(N, 'N')
    return [s.factorial(K*n)/s.factorial(n) for n in range(N+1)]


def quotient_operator(K):
    if type(K) is not int or K < 1:
        raise ValueError('K must be positive')
    return 1-t*K*s.prod(K*theta+j for j in range(1, K))


def airy_cover(K):
    if type(K) is not int or K < 2:
        raise ValueError('Airy cover requires K >= 2')
    x = s.Symbol('x')
    return dict(order=K-1, t=(-1)**(K-1)/(K*x**K), gauge=x,
                homogeneous_operator='D_x**(K-1)-x',
                connected_quotient_group=('Gm' if K == 2 else
                    'Sp_'+str(K-1) if K % 2 else 'SL_'+str(K-1)),
                group_uses_classical_theorem=True, computed_Galois_group=False)
