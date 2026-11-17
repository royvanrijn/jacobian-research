"""Exact normalization/boundary controls, not a regularity/growth proof."""
from pathlib import Path
from math import prod
import hashlib
import json
import sympy as s


def exact_poly(expr, xs, max_degree=30):
    expr = s.sympify(expr)
    if expr.has(s.Float):
        raise ValueError('Exact coefficients required')
    p = s.Poly(expr, *xs)
    if p.total_degree() > max_degree:
        raise ValueError('Control degree cap exceeded')
    return p


def wick(expr, xs):
    p = exact_poly(expr, xs)
    return s.expand(sum(c * prod(s.factorial2(a - 1) if a else 1 for a in alpha)
                        for alpha, c in p.terms() if all(a % 2 == 0 for a in alpha)))


def sphere(expr, xs):
    p = exact_poly(expr, xs)
    return s.expand(sum(c * prod(s.rf(s.Rational(1, 2), a // 2) for a in alpha)
                        / s.rf(s.Rational(len(xs), 2), sum(alpha) // 2)
                        for alpha, c in p.terms() if all(a % 2 == 0 for a in alpha)))


def normalized(P, Q, xs, n):
    pp = exact_poly(P, xs)
    if pp.is_zero or len({sum(a) for a, c in pp.terms()}) != 1 or pp.total_degree() < 1:
        raise ValueError('Nonzero homogeneous positive-degree phase required')
    degree = pp.total_degree()
    k = (degree*n + 1)//2
    radial = 2**k * s.rf(s.Rational(len(xs), 2), k)
    qp = exact_poly(Q, xs)
    total = 0
    for alpha, c in qp.terms():
        ell = sum(alpha)
        if (degree*n + ell) % 2:
            continue
        epsilon = ell % 2 if degree % 2 else 0
        a = (ell-epsilon)//2
        multiplier = prod(len(xs)+degree*n+epsilon+2*j for j in range(a))
        monomial = prod(x**e for x, e in zip(xs, alpha))
        total += c*multiplier*sphere(monomial*P**n, xs)
    direct = s.simplify(wick(Q*P**n, xs)/radial)
    assert s.expand(total-direct) == 0
    return direct


def borel_class(degree, order, tail_eventually_zero=False):
    if not isinstance(degree, int) or degree < 1:
        raise ValueError('Positive integer homogeneous degree required')
    order = s.sympify(order)
    if order.has(s.Float) or order.is_Rational is not True:
        raise ValueError('Exact rational order required in the small calculator')
    if tail_eventually_zero:
        return dict(radius='infinite', function='polynomial')
    threshold = s.Rational(degree, 2)
    if order < threshold:
        return dict(radius='zero')
    if order == threshold:
        return dict(radius='finite positive')
    return dict(radius='infinite', function='entire', exact_order=str(1/(order-threshold)))


def regular_first_order_at_infinity(a, z):
    """For f'=a(z)f, regularity at infinity is deg numerator-denominator<=-1."""
    if a == 0:
        return True
    num, den = s.fraction(s.cancel(a))
    return s.degree(num, z)-s.degree(den, z) <= -1


def run():
    x, y, z = s.symbols('x y z')
    controls = []
    models = [([x],x), ([x,y],x+s.I*y), ([x,y],x*x+y*y),
              ([x,y],x**3+s.I*y**3), ([x,y],(x+s.I*y)**2*x),
              ([x,y],x**4+y**4)]
    for xs,P in models:
        Q = 1+xs[0]+xs[0]**2+xs[-1]**3+xs[0]**4
        values=[str(normalized(P,Q,xs,n)) for n in range(5)]
        controls.append(dict(dimension=len(xs), phase_degree=int(s.total_degree(P)), normalized=values))
    factorial=[]
    for ell in [1,2,3]:
        Z,W=(x+s.I*y)/s.sqrt(2),(x-s.I*y)/s.sqrt(2)
        for n in [1,2,3]:
            actual=wick((Z*W)**(ell*n),[x,y])
            assert s.simplify(actual-s.factorial(ell*n))==0
            factorial.append(dict(degree=2*ell,index=n,moment=str(actual)))
    odd=[]
    for n in range(7):
        v=wick(x**(3*n),[x]);expected=0 if n%2 else s.factorial2(3*n-1) if n else 1
        assert v==expected
        odd.append(str(v))
    # D=2,N=2,P=x^2+y^2, Q=1+x^2: b_n=n+2, so F=(2-z)/(1-z)^2.
    F=(2-z)/(1-z)**2
    assert s.simplify(2/(1-z)+z*s.diff(1/(1-z),z)-F)==0
    assert [normalized(x*x+y*y,1+x*x,[x,y],n) for n in range(5)]==list(range(2,7))
    # Odd degree with both marker parities retained: N=1,P=x^3,Q=1+x+x^2+x^3.
    # Sphere moments are 1 or0; eq(5) selects parity before Euler operations.
    H_even=1/(1-z*z);H_odd=z/(1-z*z)
    F_odd=H_even+H_odd+(H_even+3*z*s.diff(H_even,z))+(2*H_odd+3*z*s.diff(H_odd,z))
    coeff=[s.expand(s.series(F_odd,z,0,5).removeO()).coeff(z,n) for n in range(5)]
    assert coeff==[normalized(x**3,1+x+x*x+x**3,[x],n) for n in range(5)]
    regular=[]
    for name,a,expected in [('rational',1/(1-z),True),('algebraic',z/(1-z*z),True),
                            ('constant',0,True),('exponential',1,False),('exp_square',2*z,False)]:
        actual=bool(regular_first_order_at_infinity(a,z));assert actual==expected
        regular.append(dict(model=name,regular_at_infinity=actual))
    grid=[]
    for degree in [1,2,3,4,6]:
        t=s.Rational(degree,2)
        for order in [t-1,t,t+1,t+s.Rational(1,2)]:
            grid.append(dict(degree=degree,borel_order=str(order),classification=borel_class(degree,order)))
        assert borel_class(degree,t+1)['exact_order']=='1'
        assert borel_class(degree,t+s.Rational(1,2))['exact_order']=='2'
        assert borel_class(degree,0,True)['function']=='polynomial'
    # Homogeneity boundary: Z^4+T has exactly T's pure moments, in every order
    # by E[Z^j]=0 for all positive j. These are selected exact transcriptions.
    T=s.symbols('T');Z=(x+s.I*y)/s.sqrt(2)
    boundary=[]
    for n in range(5):
        actual=wick((Z**4+T)**n,[x,y,T]);assert actual==wick(T**n,[x,y,T]);boundary.append(str(actual))
    refused=0
    for action in [lambda:normalized(x+x*x,1,[x],1),lambda:exact_poly(s.Float('0.5')*x,[x]),
                   lambda:exact_poly(x**31,[x]),lambda:borel_class(0,1),lambda:borel_class(2,0.5)]:
        try:action()
        except ValueError:refused+=1
        else:raise AssertionError('Guard did not reject')
    return dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                normalization_controls=controls,factorial_examples=factorial,odd_degree_sequence=odd,
                even_euler_identity=str(F),odd_euler_coefficients=[str(v) for v in coeff],
                infinity_regularity_controls=regular,borel_classification_controls=grid,
                inhomogeneous_boundary=boundary,guard_refusals=refused,supporting_checker_hashes={},
                normalized_moment_calculator_implemented=True,borel_threshold_calculator_implemented=True,
                universal_regularity_proved_by_finite_controls=False,asymptotic_fit_used=False,
                complex_pointwise_lower_bound_claimed=False,every_nonzero_subsequence_claimed=False,
                full_exceptional_MZ_decided=False,general_period_solver_implemented=False,
                independent_replay=False,external_review=False,formal_verification=False,novelty_claimed=False)


if __name__=='__main__':
    output=Path(__file__).with_name('complex_homogeneous_growth_verification.json')
    receipt=run()
    if output.exists():
        assert json.loads(output.read_text())==receipt,'Existing receipt differs; preserve immutable evidence'
    else:
        output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(completed=receipt['completed'],normalizations=len(receipt['normalization_controls']),
                         factorials=len(receipt['factorial_examples']),borel_controls=len(receipt['borel_classification_controls']),
                         guard_refusals=receipt['guard_refusals'],checker_sha256=receipt['checker_sha256'])))
