"""Exact original-marker coordinates for P=W(1+Z)^2Y.

Coordinates use (1,f,f',k), f=Phi[1], k=Phi[1-Z]. The companion proof
establishes independence, the full origin-saturated lattice and literal
original-image projection. No actual-period rank is inferred from jets.
"""
from __future__ import annotations
from functools import lru_cache
import sympy as sp

Z,W,Y,U,T=sp.symbols("Z W Y u t")
FRAME=("1","f","f_prime","k")
CONNECTION=sp.ImmutableMatrix([
    [0,0,0,0],
    [0,0,1,0],
    [0,-sp.Rational(3,4)/T**2,1/(16*T**3)-3/T,0],
    [0,1/T,0,-1/T],
])


def row(values):
    matrix=sp.ImmutableMatrix(values)
    if matrix.shape==(4,1):
        matrix=matrix.T
    if matrix.shape!=(1,4):
        raise ValueError("Expected four coordinates in the faithful frame")
    return sp.ImmutableMatrix([[sp.cancel(entry) for entry in matrix]])


@lru_cache(maxsize=None)
def laurent_period(n:int):
    """Coordinates of F_n=Phi[(1+Z)^(-n)] in the faithful frame."""
    n=int(n)
    if n==-1:
        return row([0,2,0,-1])
    if n==0:
        return row([0,1,0,0])
    if n==1:
        return row([sp.Rational(1,2),sp.Rational(1,2),0,0])
    if n==2:
        return row([sp.Rational(1,2),sp.Rational(1,2)-2*T*T,-4*T**3,0])
    if n>=3:
        return row((3*laurent_period(n-1)-laurent_period(n-2)
                    +(n-3)*T*T*laurent_period(n-4))/2)
    return row((2*laurent_period(n+4)-3*laurent_period(n+3)+laurent_period(n+2))
               /((n+1)*T*T))


def hermite_raise(b:int,j:int):
    """The polynomial (Y-partial_Y)^b Y^j."""
    value=Y**int(j)
    for _ in range(int(b)):
        value=sp.expand(Y*value-sp.diff(value,Y))
    return value


def laurent_amplitude(marker,allow_time=False):
    marker=sp.sympify(marker)
    if marker.has(U) or (marker.has(T) and not allow_time):
        raise ValueError("Original markers are independent of t and u")
    polynomial=sp.Poly(sp.expand(marker),Z,W,Y)
    answer=0
    for (a,b,j),coefficient in polynomial.terms():
        answer+=(coefficient*(U**-1-1)**a*(U*U/T)**b
                 *hermite_raise(b,j).subs(Y,(U-U*U)/T))
    return sp.expand(answer)


def reduce_laurent(amplitude):
    """Reduce a finite Laurent-u amplitude; coefficients may be rational in t."""
    answer=row([0,0,0,0])
    for term in sp.Add.make_args(sp.expand(amplitude)):
        exponent=term.as_powers_dict().get(U,sp.Integer(0))
        if not exponent.is_Integer:
            raise ValueError("The amplitude is not a Laurent polynomial in u")
        coefficient=sp.cancel(term/U**exponent)
        if coefficient.has(U):
            raise ValueError("The amplitude is not a Laurent polynomial in u")
        answer+=coefficient*laurent_period(int(exponent))
    return row(answer)


def reduce_marker(marker):
    """Exact actual-period coordinates; zero iff the original marker is in kerPhi."""
    return reduce_laurent(laurent_amplitude(marker))


def reduce_time_amplitude(amplitude):
    """The same linear map allowing external polynomial time coefficients."""
    for coefficient in sp.Poly(sp.expand(amplitude),Z,W,Y).coeffs():
        sp.Poly(coefficient,T)
    return reduce_laurent(laurent_amplitude(amplitude,allow_time=True))


def connection_derivative(coordinates):
    coordinates=row(coordinates)
    return row(coordinates.diff(T)+coordinates*CONNECTION)


def laurent_time_terms(expression):
    """A finite exponent->coefficient dictionary, or ValueError outside C[t,t^-1]."""
    numerator,denominator=sp.fraction(sp.cancel(expression))
    npoly,dpoly=sp.Poly(numerator,T),sp.Poly(denominator,T)
    if len(dpoly.terms())!=1:
        raise ValueError("A coordinate has a finite pole away from t=0")
    (power,),scale=dpoly.terms()[0]
    return {int(exponent-power):sp.cancel(coefficient/scale)
            for (exponent,),coefficient in npoly.terms() if coefficient}


def basis_series(order:int):
    """Exact Taylor coefficients through the requested order, including f'."""
    order=max(0,int(order))
    f=0
    for q in range((order+1)//2+2):
        coefficient=sp.factorial(4*q)/(2**q*sp.factorial(q)*sp.factorial(2*q))
        f+=coefficient*T**(2*q)
    fp=sp.diff(f,T)
    k=sum(sp.factorial(4*q)/(2**q*sp.factorial(q)*sp.factorial(2*q)*(2*q+1))
          *T**(2*q) for q in range(order//2+1))
    def truncate(value):
        return sum(coefficient*T**power for (power,),coefficient in sp.Poly(value,T).terms()
                   if power<=order)
    return [sp.Integer(1),truncate(f),truncate(fp),truncate(k)]


def truncated_series(coordinates,order:int):
    """Exact finite Laurent jet; this does not infer a series identity from jets."""
    coordinates=row(coordinates)
    terms=[laurent_time_terms(value) for value in coordinates]
    minimum=min([0]+[power for entries in terms for power in entries])
    basis=basis_series(int(order)-minimum)
    value=sp.expand(sum(coefficient*function for coefficient,function in zip(coordinates,basis)))
    return sp.expand(sum(value.coeff(T,power)*T**power
                         for power in range(minimum,int(order)+1)))


def is_in_full_image(coordinates):
    """Decide membership in M, NOT in the literal fixed-marker vector image V.

    The proved lattice description reduces this to a finite principal-part
    check whose length is fixed by the input Laurent denominators.
    """
    coordinates=row(coordinates)
    try:
        terms=[laurent_time_terms(value) for value in coordinates]
    except (ValueError,sp.PolynomialError):
        return False
    minimum=min([0]+[power for entries in terms for power in entries])
    if minimum>=0:
        return True
    return truncated_series(coordinates,-1)==0


def original_projection(coordinates):
    """Constant-line projection pi; pi(V)=t^-1 C[t^-1], pi(M)=C[t,t^-1]."""
    return row(coordinates)[0]


def original_projection_allows(coordinates):
    """Necessary condition for V membership; no sufficiency is asserted."""
    try:
        terms=laurent_time_terms(original_projection(coordinates))
    except (ValueError,sp.PolynomialError):
        return False
    return all(power<0 for power in terms)