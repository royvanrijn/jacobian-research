"""PF106 necessary quotient/scalar-Morse data, never a Keller realization test."""
import math
import sympy as S


def quotient_constraint(e,n):
    if not isinstance(e,int) or not isinstance(n,int) or e<1 or n<2:
        raise ValueError('positive transverse index and nontrivial quotient order required')
    return dict(necessary_bound_satisfied=(n<=e-1 and math.gcd(e,n)==1),
                marked_weights=[1,(-e)%n],necessary_only=True)


def morse_budget(e,n,f,a,field_degree=None):
    """Use ONLY if an actual scalar pullback has nondegenerate quadratic term."""
    if any(not isinstance(t,int) for t in [e,n,f,a]) or e<2 or n<2 or f<1 or a<1:
        raise ValueError('positive integral point data and nontrivial quotient required')
    if field_degree is not None and (not isinstance(field_degree,int) or field_degree<2):
        raise ValueError('nontrivial integral generic field degree required')
    k=e-1;m=n*f*a;ell=k*f*a;eta=ell-e*f
    allowed=(quotient_constraint(e,n)['necessary_bound_satisfied'] and k%(n*f)==0 and
             a>=2 and k%m!=0)
    return dict(necessary_local_conditions=allowed,local_degree=ell,eta=eta,
                minimum_field_degree=ell+1,quadratic_field_lower_bound=4*f*f+1,
                actual_field_budget=(field_degree is None or field_degree>=ell+1),
                scalar_Morse_hypothesis_required=True,necessary_only=True,realization_asserted=False)


def morse_model(k,m):
    """Genuine finite smooth germs; quotient models use n|gcd(k,m).

    A full polynomial model has exact upstairs residue gcd(k,m).
    Higher analytic terms can make that residue a proper divisor.
    """
    if not isinstance(k,int) or not isinstance(m,int) or k<2 or m<2:
        raise ValueError('k>=2 and m>=2 required')
    if k%m==0:
        raise ValueError('nonzero diagonal resonance obstructs an analytic primitive')
    X,Y=S.symbols('X Y');b=m-1
    coefficients=[S.Rational((-1)**j*S.binomial(k,j),k-m*j) for j in range(k+1)]
    u=S.expand(sum(coefficients[j]*X**(k-j)*Y**(b*j) for j in range(k+1)))
    q=math.gcd(k,m);a=m//q;beta=k*b//q
    c=-S.factorial(k)/(m*S.prod(S.Rational(j)-S.Rational(k,m) for j in range(k+1)))
    return dict(X=X,Y=Y,u=u,v=X*Y,k=k,m=m,b=b,critical=X-Y**b,
                coefficients=coefficients,coefficient=S.factor(c),upstairs_degree=k*m,
                upstairs_residue=q,image_pair=[a,beta],conductor=(a-1)*(beta-1))


def conductor_intersection(residue,conductor):
    if not isinstance(residue,int) or not isinstance(conductor,int) or residue<1 or conductor<0:
        raise ValueError('positive residue degree and nonnegative conductor required')
    return residue*(conductor-1)+1
