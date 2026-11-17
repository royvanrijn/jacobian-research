"""Finite polynomial recovery witnesses for amplitude coefficient content.

Symbols other than X and the supplied time variables are coefficient
parameters. No actual image kernel,Gaussian rank or safety is inferred.
"""
import itertools
import sympy as s

def multiindices(n,total):
    return [v for v in itertools.product(range(total+1),repeat=n) if sum(v)<=total]

def setup(phases,X,times,C):
    if not X or len(phases)!=len(times) or not times:
        raise ValueError('nonempty variables and matched phase/time lists required')
    C=s.Matrix(C)
    if C.shape!=(len(X),len(X)) or C!=C.T or C.det()==0 or any(not c.is_Rational for c in C):
        raise ValueError('fixed nonsingular rational symmetric covariance required')
    if len(set(X+times))!=len(X)+len(times):
        raise ValueError('spatial and time variables must be distinct')
    for P in phases:s.Poly(P,*X)
    return C.inv()

def adjoint(poly,i,phases,X,times,C):
    inv=setup(phases,X,times,C)
    drift=sum(inv[i,j]*X[j] for j in range(len(X)))
    drift-=sum(tj*s.diff(P,X[i]) for P,tj in zip(phases,times))
    return s.expand(drift*poly-s.diff(poly,X[i]))

def adjoint_multi(poly,powers,phases,X,times,C):
    for i,k in enumerate(powers):
        for _ in range(k):poly=adjoint(poly,i,phases,X,times,C)
    return s.expand(poly)

def coefficient_witness(weight,powers,amplitude,phases,X,times,C,top_only=False):
    setup(phases,X,times,C)
    if len(powers)!=len(X) or any(type(k) is not int or k<0 for k in powers):
        raise ValueError('nonnegative integer multiindex required')
    w=s.Poly(weight,*X)
    if w.is_zero:raise ValueError('nonzero weight required for a recovery witness')
    degree=int(w.total_degree())
    if sum(powers)>degree:raise ValueError('requested index exceeds weight degree')
    if top_only and sum(powers)!=degree:raise ValueError('top recovery needs a top-degree index')
    c=w.coeff_monomial(s.prod(x**k for x,k in zip(X,powers)))
    if top_only and c==0:raise ValueError('top recovery needs a nonzero coefficient')
    out=0;projection=0;fac=s.prod(s.factorial(k) for k in powers)
    betas=[(0,)*len(X)] if top_only else multiindices(len(X),degree-sum(powers))
    for beta in betas:
        gamma=tuple(a+b for a,b in zip(powers,beta));monomial=s.prod(x**k for x,k in zip(X,beta))
        scalar=s.Rational((-1)**sum(beta),fac*s.prod(s.factorial(k) for k in beta))
        out+=scalar*adjoint_multi(monomial*amplitude,gamma,phases,X,times,C)
        deriv=weight
        for x,k in zip(X,gamma):deriv=s.diff(deriv,x,k)
        projection+=scalar*monomial*deriv
    assert s.expand(projection-c)==0
    dp=max([int(s.Poly(P,*X).total_degree()) for P in phases]+[0]);delta=max(1,dp-1)
    bound=int(s.Poly(amplitude,*X).total_degree())+((delta*degree) if top_only else ((delta+1)*degree-sum(powers)))
    out=s.expand(out)
    assert s.Poly(out,*X).total_degree()<=bound
    assert s.Poly(out,*times).total_degree()<=degree
    return dict(coefficient=str(c),weight=str(weight),index=list(powers),amplitude=str(amplitude),phases=list(map(str,phases)),variables=list(map(str,X)),times=list(map(str,times)),covariance=[[str(z) for z in row] for row in s.Matrix(C).tolist()],recovery_amplitude=str(out),Taylor_projection=str(s.expand(projection)),top_only=top_only,spatial_degree_bound=bound,time_degree_bound=degree,requires_coefficient_unit_to_divide=bool(top_only),actual_image_kernel_computed=False)

def cutoff_transfer(C,E,weight_degree,phase_degree,direction):
    C=s.Rational(C);E=s.Rational(E)
    if C<0 or E<0 or type(weight_degree) is not int or weight_degree<0 or type(phase_degree) is not int or phase_degree<0:
        raise ValueError('nonnegative costs and integer degrees required')
    delta=max(1,phase_degree-1);D=weight_degree
    if direction=='original_to_weighted':intercept=E+C*D
    elif direction=='weighted_to_original':intercept=E+(C*delta+1)*D
    else:raise ValueError('unknown transfer direction')
    return dict(C=str(C),E=str(intercept),direction=direction,weight_degree=D,phase_degree=phase_degree,sufficient_only=True,requires_field_or_unit_top_coefficient=direction=='weighted_to_original',sharp_intercept_claimed=False)
