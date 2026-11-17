"""PF105 local models and necessary actual-Keller residue filters.

Written analytic proofs and declared primary boundary inputs supply the
universal theorem. Passing a filter does not realize any finite Keller map.
"""
import math
import sympy as S


def singular_triple_filter(quotient_order,residue_local_degree):
    if not isinstance(quotient_order,int) or not isinstance(residue_local_degree,int) or \
       quotient_order<2 or residue_local_degree<1:
        raise ValueError('nontrivial cyclic quotient order and positive local degree required')
    return dict(upstairs_residue_degree=quotient_order*residue_local_degree,
                necessary_bound_satisfied=quotient_order*residue_local_degree<=2,
                necessary_only=True)


def actual_keller_residue_filter(e,f,field_degree=None):
    if not isinstance(e,int) or not isinstance(f,int) or e<1 or f<1:
        raise ValueError('positive integral transverse and residue degrees required')
    if field_degree is not None and (not isinstance(field_degree,int) or field_degree<2):
        raise ValueError('nontrivial integral field degree required')
    return dict(low_index_residue_one=(f==1 or e>=4),
                generic_cost_budget=(field_degree is None or e*f<=field_degree-1),
                necessary_only=True,realization_asserted=False)


def rank_zero_model(b,perturb=0):
    if not isinstance(b,int) or b<2:
        raise ValueError('tangency order b>=2 required')
    X,Y=S.symbols('X Y')
    u=X**2/2+S.Rational(2,b-1)*X*Y**b-S.Rational(1,2*b)*Y**(2*b)
    u+=perturb*(X**2*Y+S.Rational(2,b)*X*Y**(b+1)-Y**(2*b+1)/S.Integer(2*b+1))
    v=X*Y;restricted=S.expand(u.subs(X,Y**b))
    residue_degree=math.gcd(b+1,*(int(monom[0]) for monom,coefficient in S.Poly(restricted,Y).terms()))
    return dict(X=X,Y=Y,u=u,v=v,critical=X-Y**b,unit=1+perturb*Y,
                restricted_u=restricted,upstairs_residue_degree=residue_degree,
                local_surface_degree=2*(b+1))


def quotient_counterfamily(s):
    """Genuine normal finite maps, but complement C* times C, not A2."""
    if not isinstance(s,int) or s<1:
        raise ValueError('positive integer s required')
    X,Y,g,A,B,C=S.symbols('X Y g A B C');p=2*s;m=4*s;b=4*s-1
    coefficients=[S.Rational((-1)**j*S.binomial(p,j),p-m*j) for j in range(p+1)]
    u=S.expand(sum(coefficients[j]*X**(p-j)*Y**(b*j) for j in range(p+1)));v=X*Y
    shifted=S.Poly(S.expand(u.subs(X,g+Y**b)),g,Y)
    downstairs=0
    for (i,j),coefficient in shifted.terms():
        assert (i+j)%2==0
        odd=i%2;assert j%2==odd
        downstairs+=coefficient*B**odd*A**((i-odd)//2)*C**((j-odd)//2)
    c=S.Rational(2**p*S.factorial(p-1),S.factorial2(2*p-1))
    return dict(s=s,p=p,m=m,b=b,X=X,Y=Y,g=g,A=A,B=B,C=C,
                u=u,v=v,coefficients=coefficients,downstairs_u=S.expand(downstairs),
                downstairs_v=B+C**(2*s),leading_boundary_coefficient=c,
                transverse_index=2*s+1,upstairs_residue_degree=2*s,residue_degree=s,
                upstairs_degree=8*s*s,quotient_degree=4*s*s,
                eta=2*s*s-s,complement='C* times C',finite_Keller_realization=False)
