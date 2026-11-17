"""Exact full-image coordinates for Phi_m[A]=E[A exp(tW(1+Z)^mY)].

The source frame is F_(1-m),...,F_m. Its faithfulness, actual rank 2m
and exact origin-saturated global lattice are proved in the companion
note, not inferred from the finite controls in the checker.
"""
from __future__ import annotations
from functools import lru_cache
import sympy as sp

Z,W,Y,U,T,Q=sp.symbols("Z W Y u t q")


def validate_m(m):
    value=sp.sympify(m)
    if not value.is_Integer or value<1:
        raise ValueError("m must be a positive integer")
    return int(value)


def frame_indices(m):
    m=validate_m(m)
    return tuple(range(1-m,m+1))


def row(values,m):
    m=validate_m(m)
    matrix=sp.ImmutableMatrix(values)
    if matrix.shape==(2*m,1):
        matrix=matrix.T
    if matrix.shape!=(1,2*m):
        raise ValueError("Expected 2m coordinates in the source window")
    return sp.ImmutableMatrix([[sp.cancel(value) for value in matrix]])


@lru_cache(maxsize=None)
def source_period(n,m):
    """The exact Laurent-recurrence coordinates of F_n."""
    m=validate_m(m)
    n=int(n)
    indices=frame_indices(m)
    if n in indices:
        return row([int(index==n) for index in indices],m)
    if n>m:
        s=n-2*m
        return row(((2*m-1)*source_period(n-1,m)-(m-1)*source_period(n-2,m)
                    +(s+m-1)*T*T*source_period(s,m))/m,m)
    if m==1:
        return row((source_period(n+2,m)-source_period(n+1,m))/(n*T*T),m)
    return row((m*source_period(n+2*m,m)-(2*m-1)*source_period(n+2*m-1,m)
                +(m-1)*source_period(n+2*m-2,m))/((n+m-1)*T*T),m)


@lru_cache(maxsize=None)
def connection_matrix(m):
    m=validate_m(m)
    rows=[]
    for n in frame_indices(m):
        rows.append(row(-source_period(n,m)/T
                    +(source_period(n+2*m-2,m)-2*source_period(n+2*m-1,m)
                      +source_period(n+2*m,m))/T**3,m))
    return sp.ImmutableMatrix.vstack(*rows)


def connection_derivative(coordinates,m):
    m=validate_m(m)
    coordinates=row(coordinates,m)
    return row(coordinates.diff(T)+coordinates*connection_matrix(m),m)


def reduce_marker(marker,m):
    """Exact coordinates of a fixed original polynomial marker."""
    m=validate_m(m)
    marker=sp.sympify(marker)
    if marker.has(T,U,Q):
        raise ValueError("Original markers cannot contain the reserved t,u,q symbols")
    amplitude=0
    for (a,b,j),coefficient in sp.Poly(sp.expand(marker),Z,W,Y).terms():
        hermite=Y**j
        for _ in range(b):
            hermite=sp.expand(Y*hermite-sp.diff(hermite,Y))
        amplitude+=(coefficient*(U**-1-1)**a*(U**m/T)**b
                    *hermite.subs(Y,U**(m-1)*(1-U)/T))
    answer=row([0]*(2*m),m)
    for term in sp.Add.make_args(sp.expand(amplitude)):
        power=term.as_powers_dict().get(U,sp.Integer(0))
        if not power.is_Integer:
            raise ValueError("Internal amplitude is not Laurent polynomial in u")
        coefficient=sp.cancel(term/U**power)
        if coefficient.has(U):
            raise ValueError("Internal amplitude is not Laurent polynomial in u")
        answer+=coefficient*source_period(int(power),m)
    return row(answer,m)


def primitive_marker(m):
    m=validate_m(m)
    if m==1:
        raise ValueError("The separate polynomial primitive exists for m>=2")
    r=m-1
    return sp.cancel((1-r*Z)*((1+Z)**r-1)/(r*Z)).expand()


def unit_operator(m):
    """Return C,N(q),D(q) for D(vartheta)-C t^2 N(vartheta+1)."""
    m=validate_m(m)
    if m==1:
        raise ValueError("Use the exceptional inhomogeneous equation for m=1")
    c=sp.Rational(2*m**(2*m),(m-1)**(2*m-2))
    n=sp.prod(Q-sp.Rational(j,2*m) for j in range(1,2*m) if j!=m)
    d=sp.prod(Q-sp.Rational(j,2*(m-1)) for j in range(2*m-2) if j!=m-1)
    return c,sp.Poly(n,Q).as_expr(),sp.Poly(d,Q).as_expr()


@lru_cache(maxsize=None)
def unit_primitive_frame(m):
    """Rows of the proved global R-frame in source-window coordinates."""
    m=validate_m(m)
    one=m*source_period(1,m)-(m-1)*source_period(0,m)
    f=source_period(0,m)
    if m==1:
        return sp.ImmutableMatrix.vstack(one,f)
    rows=[one,f]
    for _ in range(2*m-3):
        rows.append(row(T*connection_derivative(rows[-1],m)/2,m))
    rows.append(reduce_marker(primitive_marker(m),m))
    return sp.ImmutableMatrix.vstack(*rows)


def laurent_time_terms(value):
    numerator,denominator=sp.fraction(sp.cancel(value))
    npoly,dpoly=sp.Poly(numerator,T),sp.Poly(denominator,T)
    if len(dpoly.terms())!=1:
        raise ValueError("Coordinate is not Laurent polynomial in t")
    (power,),scale=dpoly.terms()[0]
    return {int(exponent-power):sp.cancel(coefficient/scale)
            for (exponent,),coefficient in npoly.terms() if coefficient}


def truncated_series(coordinates,m,order):
    """Finite exact Laurent jet; its input denominator fixes principal-part length."""
    m=validate_m(m)
    coordinates=row(coordinates,m)
    terms=[laurent_time_terms(value) for value in coordinates]
    minimum=min([0]+[power for entries in terms for power in entries])
    cutoff=int(order)-minimum
    bases=[]
    for n in frame_indices(m):
        bases.append(sum(sp.factorial2(2*q-1)*sp.binomial(2*m*q-n,2*q)*T**(2*q)
                         for q in range(max(0,cutoff//2+1))))
    value=sp.expand(sum(coefficient*basis for coefficient,basis in zip(coordinates,bases)))
    return sp.expand(sum(value.coeff(T,power)*T**power
                         for power in range(minimum,int(order)+1)))


def is_in_full_image(coordinates,m):
    """Decide membership in M_m, not in its literal fixed-marker subspace V_m."""
    m=validate_m(m)
    coordinates=row(coordinates,m)
    try:
        terms=[laurent_time_terms(value) for value in coordinates]
    except (ValueError,sp.PolynomialError):
        return False
    minimum=min([0]+[power for entries in terms for power in entries])
    return minimum>=0 or truncated_series(coordinates,m,-1)==0
