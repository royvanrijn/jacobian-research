#!/usr/bin/env python3
"""Exact terminal-seven paired replay using the linked sextic and septic skews."""
from fractions import Fraction as Q

from verify_multivariate_terminal_seven_light_cap_extension import (
    CAP,RAYS as PRECEDING_RAYS,tight_envelope,
)
from verify_multivariate_terminal_seven_paired_light import (
    ANCHOR_RATIO,check_anchor,exponents,nominal,radial,
)
from verify_multivariate_terminal_seven_private_linear import (
    multiply_polynomials,verify_skew_convexity,
)
from verify_multivariate_terminal_six_diffuse_bands import bernstein_coefficients
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID,I,ONE,ZERO,interval,log,rational,require,scalar,
)

RAYS = {
    480: PRECEDING_RAYS[480],
    500: (83,-148,-549,-703,-924),
    520: (43,-161,-569,-716,-928),
    540: (0,-195,-589,-730,-931),
    560: (-41,-230,-609,-743,-935),
    580: (-82,-264,-632,-747,-938),
    600: (-124,-298,-651,-761,-944),
    620: (-166,-333,-670,-774,-948),
}


def validate_transfer():
    eta=(0,-1,2,2,-1,-1)
    gamma=(0,-1,-1,6,-1,-1,-1)

    def derivative(weights,denominator,scale,weight):
        return [(Q(a,denominator)*n*scale*weight,n*scale)
                for n,a in enumerate(weights) if n and a]

    for entries in (
        derivative(gamma,14,Q(6,7),1)+derivative(eta,6,1,Q(-107,200)),
        derivative(eta,6,1,Q(29,50))+derivative(gamma,14,Q(48,55),-1),
    ):
        polynomial=[Q(0)]*len(entries)
        for i,(num,_) in enumerate(entries):
            basis=[num]
            for j,(_,den) in enumerate(entries):
                if i!=j:
                    basis=multiply_polynomials(basis,(1,den))
            polynomial=[a+b for a,b in zip(polynomial,basis,strict=True)]
        require(polynomial[0]==0,"skew transfer derivative at zero")
        require(min(bernstein_coefficients(polynomial[1:],Q(0),CAP))>0,
                ("skew transfer exact derivative",polynomial))


def model(order,u,skew,envelope,eta_cap):
    rho, variables = envelope
    mid=(u[0]+u[1])/2
    ref=nominal(order,mid)
    U=-log(interval(*u))
    eta=U*eta_cap*interval(*skew)
    mean=log(interval(u[0]/mid,u[1]/mid))
    if order==6:
        W=U
        primary=eta
    else:
        mean=ANCHOR_RATIO*mean+(ANCHOR_RATIO-rho)*U
        W=rho*U
        error=(scalar(1,40)*U*eta_cap).symmetric_error()
        gamma=scalar(14,25)*eta+error
        bound=variables[0][1]*W
        primary=I(max(gamma.lo,bound.lo),min(gamma.hi,bound.hi))
        if primary.lo>primary.hi:
            return None
    pieces=((primary,variables[0][0]),)
    pieces+=tuple((p*W,vec) for vec,p in variables[1:])
    pieces+=((mean,tuple(rational(e) for e in exponents(order))),)
    remainder=[]
    for k,a in enumerate(ref):
        radius=ZERO
        for parameter,vector in pieces:
            value=parameter*vector[k]
            radius=radius+I(0,max(-value.lo,value.hi))
        require(radius.hi<GRID,("exponential remainder",order,u))
        error=scalar(1,2)*radius.square()/(ONE-radius)
        remainder.append(a*(ONE+error))
    return ((ONE,tuple(remainder)),)+tuple(
        (p,tuple(a*v for a,v in zip(ref,vector,strict=True))) for p,vector in pieces
    )


def cosine_interval(left,right,index,u):
    a,b=RAYS[left][index],RAYS[right][index]
    values=tuple(Q(a,1000)+Q(b-a,right-left)*(e-Q(left,1000)) for e in u)
    return interval(min(values),max(values))


def verify(order,index,left,right,env,eta_cap):
    # Reflections allow the actual sextic inner skew to be nonnegative.
    u0,u1=Q(left,1000),Q(right,1000)
    pending=[((u0,u1),(Q(0),Q(1)),0,0)]
    count=0
    while pending:
        u,skew,du,ds=pending.pop()
        coefficients=model(order,u,skew,env,eta_cap)
        if coefficients is None:
            continue
        used=radial(order,cosine_interval(left,right,index,u),coefficients)
        if used is not None:
            count+=used
            continue
        require(max(du,ds)<15,("unresolved skew-conditioned ray",order,index,left,right,u,skew))
        if du<ds+1:
            mid=sum(u)/2
            pending.extend((((u[0],mid),skew,du+1,ds),((mid,u[1]),skew,du+1,ds)))
        else:
            mid=sum(skew)/2
            pending.extend(((u,(skew[0],mid),du,ds+1),(u,(mid,skew[1]),du,ds+1)))
    return count


def main():
    validate_transfer()
    verify_skew_convexity()
    check_anchor()
    require(RAYS[480]==PRECEDING_RAYS[480],"matching preceding Gaussian sectors")
    for row in RAYS.values():
        require(all(1000>a>b>-1000 for a,b in zip(row,row[1:])),"ordered transfer rays")
    bounds={j:tight_envelope(j) for j in (6,7)}
    eta_cap=I(bounds[6][1][0][1].hi,bounds[6][1][0][1].hi)
    knots=tuple(RAYS)
    for left,right in zip(knots,knots[1:]):
        counts=[verify(j,i,left,right,bounds[j],eta_cap) for j in (6,7) for i in range(5)]
        print(left,right,counts,flush=True)


if __name__=="__main__":
    main()
