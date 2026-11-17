"""Exact necessary budgets; a successful test never asserts realization."""
from math import gcd


def leading_budget(r, s, n, f, field_degree=None):
    if min(r, s, n, f) < 1 or n < 2:
        raise ValueError('positive orders/residue and nontrivial quotient required')
    g = gcd(r, s)
    a, b = r // g, s // g
    e = r + s - 1
    necessary = (r >= 2 and s >= 2 and g % (n*f) == 0
                 and e <= r*s//g)
    ell = r*s//n if r*s % n == 0 else None
    return dict(orders=[r, s], leading_gcd=g, primitive_orders=[a, b],
                transverse_index=e, upstairs_residue=n*f,
                projective_degree=r*s//g, marked_weights=[1, 1],
                necessary_local_conditions=necessary,
                local_degree=ell, minimum_field_degree=None if ell is None else ell+1,
                actual_field_budget=None if field_degree is None else
                    necessary and ell is not None and field_degree >= ell+1,
                realization_claimed=False)


def cusp_family(q, n):
    if q < 1 or n < 2 or q % n:
        raise ValueError('q must be positive and divisible by quotient order n>1')
    return dict(q=q, quotient_order=n, transverse_index=5*q-1,
                residue_degree=q//n, degree=6*q*q//n,
                finite_excess=q*(q+1)//n, marked_weights=[1, 1],
                non_Gorenstein=n>2, complement='C* times C',
                existence_input='Zannier passport theorem, Orevkov Theorem1',
                actual_Keller_source=False)


def explicit_pair(q):
    import sympy as S
    X, Y, t = S.symbols('X Y t')
    if q == 1:
        p=t*t+1; z=t**3+S.Rational(3,2)*t
    elif q == 2:
        p=t**4+t; z=t**6+S.Rational(3,2)*t**3+S.Rational(3,8)
    else:
        raise ValueError('only two explicit controls; unbounded existence is a theorem')
    P=S.expand(X**(2*q)*p.subs(t,Y/X))
    Q=S.expand(X**(3*q)*z.subs(t,Y/X))
    return dict(X=X,Y=Y,t=t,p=p,z=z,P=P,Q=Q)
