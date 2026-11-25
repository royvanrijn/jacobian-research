#!/usr/bin/env python3
"""Exact paired Bernstein replay for shared terminal-eight profiles on .72 to .76."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from sys import argv

from verify_multivariate_terminal_eight_paired_bridge_72 import FINISH
from verify_multivariate_terminal_seven_paired_light import exponents, signed, halves
from verify_multivariate_terminal_seven_private_linear import multiply_polynomials
from verify_multivariate_terminal_six_diffuse_bands import bernstein_coefficients
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID,I,ONE,ZERO,interval,log,rational,require,scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import power_bounds


ANCHOR=Q(49,48)
INNER,OUTER=Q(1,12),Q(12)
CAPS=((720,730,Q(1,25)),(730,740,Q(1,25)),(740,750,Q(1,25)),
      (750,760,Q(1,25)))
RAYS={
    100:(828,684,132,-207,-607,-856),
    720:(-347,-447,-703,-770,-920,-946),
    730:(-369,-467,-714,-778,-923,-950),
    740:(-392,-486,-727,-787,-928,-950),
    750:(-414,-516,-738,-796,-931,-952),
    760:(-437,-535,-749,-804,-936,-952),
    770:(-460,-554,-760,-808,-938,-954),
    780:(-484,-573,-771,-817,-943,-955),
}
W7=(0,-1,-1,6,-1,-1,-1)
W8=(0,-1,-1,3,3,-1,-1,-1)


def validate_skew():
    def derivative(weights,denominator,scale,weight):
        return [(Q(a,denominator)*n*scale*weight,Q(n)*scale)
                for n,a in enumerate(weights) if n and a]
    for low,entries in (
        (True, derivative(W8,8,Q(7,8),1)+derivative(W7,14,1,Q(-189,100))),
        (False, derivative(W7,14,1,Q(1977,1000))+derivative(W8,8,Q(8,9),-1)),
    ):
        polynomial=[Q(0)]*len(entries)
        for i,(num,_) in enumerate(entries):
            basis=[num]
            for j,(_,den) in enumerate(entries):
                if i!=j:
                    basis=multiply_polynomials(basis,(1,den))
            polynomial=[a+b for a,b in zip(polynomial,basis,strict=True)]
        require(polynomial[0]==0,"skew transfer derivative at zero")
        if low:
            require(all(a>0 for a in polynomial[1:4]) and all(a<0 for a in polynomial[4:]),
                    "at most one turn in lower skew comparison")
            last=ZERO
            for num,den in entries:
                last=last+rational(num/den)*log(ONE+rational(den/8))
            require(last.lo>0,"lower skew comparison endpoint")
        else:
            require(min(bernstein_coefficients(polynomial[1:],Q(0),Q(1,8)))>0,
                    "upper skew comparison derivative")

    for weights in (W7,W8):
        m=len(weights)-1
        hs=tuple(k for k in range(1,m) if weights[k])
        polynomial=[Q(0)]*(2*len(hs)-1)
        for k in hs:
            factor=[Q(weights[k]*k*(m-k))]
            for j in hs:
                if k!=j:
                    factor=multiply_polynomials(factor,(1,2*j,j*j))
            polynomial=[a+b for a,b in zip(polynomial,factor,strict=True)]
        require(min(bernstein_coefficients(polynomial,Q(0),Q(1,8)))>0,"inner skew convexity")


@lru_cache(maxsize=None)
def nominal(order,mean):
    ratio=Q(1) if order==7 else ANCHOR
    return tuple(scalar(comb(order,k))*interval(*power_bounds(mean,ratio*e)) if e else ONE
                 for k,e in enumerate(exponents(order)))


def envelope(cap,order):
    t=rational(cap)
    ratio_lower=log(ONE+scalar(49,8)*t)/log(ONE+scalar(6)*t)
    ratio_upper=ONE-log(ONE-t*scalar(1,8))/log(ONE+scalar(6)*t)
    rho=I(ratio_lower.lo,ratio_upper.hi)
    require(rho.lo<=rational(ANCHOR).lo<=rational(ANCHOR).hi<=rho.hi,"Gaussian ratio")
    x=t if order==7 else scalar(7)*t/(scalar(8)-t)
    s=scalar(1,2)*log(ONE+scalar(order-1)*x)
    weights=W7 if order==7 else W8
    quotient=14 if order==7 else 8
    skew=ZERO
    for h,c in enumerate(weights):
        if h and c:
            skew=skew+scalar(c,quotient)*log(ONE+scalar(h)*x)
    bound=(skew/s).symmetric_error()
    if order==7:
        H=scalar(2,3)-scalar(1,2)*log((ONE+scalar(5)*x)/(ONE+x))/s
        L=scalar(3)*x**4/s
        K=scalar(60,7)*x**6/s
        E=s**4*scalar(1,25)
        pieces=(
            ((0,Q(5,2),Q(5,2),1,-1,Q(-5,2),Q(-5,2),0),bound),
            ((0,Q(3,2),Q(1,2),0,0,Q(-1,2),Q(-3,2),0),L.symmetric_error()),
            ((0,Q(-1,2),Q(-1,2),0,0,Q(1,2),Q(1,2),0),K.symmetric_error()),
            ((0,0,1,Q(9,5),Q(9,5),1,0,0),H.positive_error()),
            ((0,0,0,-1,-1,0,0,0),E.positive_error()),
        )
    else:
        H=scalar(5,7)-scalar(1,2)*log((ONE+scalar(6)*x)/(ONE+x))/s
        s2=s.square()
        a,b,c=(ONE-scalar(k,49)*s2 for k in (25,9,1))
        J3=scalar(120*16,343*49*5)*s**4/(a*b)
        J4=scalar(48*24*8,343*49**2*7)*s**6/(a*b*c)
        L=scalar(42)*x**4/s
        K=scalar(60)*x**6/s
        pieces=(
            ((0,Q(7,5),Q(8,5),1,0,-1,Q(-8,5),Q(-7,5),0),bound),
            ((0,Q(1,5),Q(3,35),0,0,0,Q(-3,35),Q(-1,5),0),L.symmetric_error()),
            ((0,0,Q(-1,7),0,0,0,Q(1,7),0,0),K.symmetric_error()),
            ((0,0,1,2,Q(12,5),2,1,0,0),H.positive_error()),
            ((0,0,0,-1,Q(-8,5),-1,0,0,0),J3.positive_error()),
            ((0,0,0,0,1,0,0,0,0),J4.positive_error()),
        )
    return rho,tuple((tuple(rational(v) for v in dx),scale) for dx,scale in pieces)


@lru_cache(maxsize=8192)
def coefficient_model(order,u,normalized,env,lower_bound):
    rho,variables=env
    mid=(u[0]+u[1])/2
    ref=nominal(order,mid)
    U=-log(interval(*u))
    log_shift=log(interval(u[0]/mid,u[1]/mid))
    primary=U*lower_bound*interval(*normalized)
    if order==7:
        strength=U
    else:
        log_shift=rational(ANCHOR)*log_shift+(rational(ANCHOR)-rho)*U
        strength=rho*U
        transferred=scalar(1933,1000)*primary+(scalar(44,1000)*lower_bound*U).symmetric_error()
        direct=variables[0][1]*strength
        primary=I(max(transferred.lo,direct.lo),min(transferred.hi,direct.hi))
        if primary.lo>primary.hi:
            return None
    pieces=((primary,variables[0][0]),)
    pieces+=tuple((p*strength,vector) for vector,p in variables[1:])
    pieces+=((log_shift,tuple(rational(e) for e in exponents(order))),)
    g,gv=pieces[0]
    m,mv=pieces[-1]
    cross=(
        (g.square(),tuple(v.square()*scalar(1,2) for v in gv)),
        (g*m,tuple(v*w for v,w in zip(gv,mv,strict=True))),
        (m.square(),tuple(v.square()*scalar(1,2) for v in mv)),
    )
    remainder=[]
    for k,base in enumerate(ref):
        shifted=g*gv[k]+m*mv[k]
        secondary=ZERO
        for parameter,vector in pieces[1:-1]:
            secondary=secondary+parameter*vector[k]
        signed_total=shifted+secondary
        radius=I(0,max(-signed_total.lo,signed_total.hi))
        require(radius.hi<GRID,"exponential remainder")
        cubic=(radius**3/(scalar(6)*(ONE-radius*scalar(1,4)))).symmetric_error()
        error=shifted*secondary+scalar(1,2)*secondary.square()+cubic
        remainder.append(base*(ONE+error))
    return ((ONE,tuple(remainder)),)+tuple(
        (p,tuple(a*b for a,b in zip(ref,v,strict=True))) for p,v in pieces+cross
    )


@lru_cache(maxsize=None)
def bernstein_transform(order):
    return tuple(
        tuple(rational(sum(
            Q(comb(k,l)*comb(i,l),comb(order,l))*INNER**(k-l)*(OUTER-INNER)**l
            for l in range(min(k,i)+1)
        )) for k in range(order+1))
        for i in range(order+1)
    )


@lru_cache(maxsize=None)
def sine_polynomials(order):
    result=[(0,),(1,)]
    for _ in range(1,order):
        a=(0,)+tuple(2*c for c in result[-1])
        b=result[-2]+(0,)*(len(a)-len(result[-2]))
        result.append(tuple(v-w for v,w in zip(a,b,strict=True)))
    return tuple(result)


@lru_cache(maxsize=16384)
def sine_intervals(order,cos_lo,cos_hi):
    lo,hi=Q(cos_lo,GRID),Q(cos_hi,GRID)
    out=[]
    for p in sine_polynomials(order):
        values=bernstein_coefficients(p,lo,hi)
        out.append(interval(min(values),max(values)))
    return tuple(out)


def projected(order,cosines,model):
    sines=sine_intervals(order,cosines.lo,cosines.hi)
    transform=bernstein_transform(order)
    result=[]
    for h in range(order+1):
        directions=tuple(sines[k-h] if k>=h else -sines[h-k] for k in range(order+1))
        groups=[]
        for _,coefficients in model:
            product=tuple(a*b for a,b in zip(coefficients,directions,strict=True))
            bounds=[]
            for row in transform:
                lo=hi=0
                for a,b in zip(row,product,strict=True):
                    vals=(a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi)
                    lo+=min(vals)
                    hi+=max(vals)
                bounds.append(I(lo//GRID,-(-hi//GRID)))
            groups.append(tuple(bounds))
        result.append(tuple(groups))
    return tuple(result)


@lru_cache(maxsize=None)
def binomial_weights(order):
    return tuple(comb(order,k) for k in range(order+1))


def midpoint_signed(projection,model,order):
    weights=binomial_weights(order)
    lo=hi=0
    for controls,(parameter,_) in zip(projection,model,strict=True):
        lower=sum(a.lo*b for a,b in zip(controls,weights,strict=True))
        upper=sum(a.hi*b for a,b in zip(controls,weights,strict=True))
        products=(parameter.lo*lower,parameter.lo*upper,parameter.hi*lower,parameter.hi*upper)
        lo+=min(products)
        hi+=max(products)
    return hi<0 or lo>0


def radial(order,cosines,model):
    stack=[(projected(order,cosines,model),0)]
    used=0
    while stack:
        polys,depth=stack.pop()
        if any(signed(candidate,model) for candidate in polys):
            used+=1
            continue
        if depth>=14 or not any(midpoint_signed(candidate,model,order) for candidate in polys):
            return None
        first,second=[],[]
        for polys_at_rotation in polys:
            children=tuple(halves(group) for group in polys_at_rotation)
            fl,fr=zip(*children,strict=True)
            first.append(fl)
            second.append(fr)
        stack.extend(((tuple(first),depth+1),(tuple(second),depth+1)))
    return used


def ray_cosines(left,right,index,u):
    a,b=RAYS[left][index],RAYS[right][index]
    values=tuple(Q(a,1000)+Q(b-a,right-left)*(v-Q(left,1000)) for v in u)
    return interval(min(values),max(values))


def verify(order,index,left,right,env,lower_bound):
    u0,u1=Q(left,1000),Q(right,1000)
    stack=[((u0,u1),(Q(0),Q(1)),0,0)]
    count=0
    while stack:
        u,normalized,du,ds=stack.pop()
        coefficients=coefficient_model(order,u,normalized,env,lower_bound)
        if coefficients is None:
            continue
        used=radial(order,ray_cosines(left,right,index,u),coefficients)
        if used is not None:
            count+=used
            continue
        require(max(du,ds)<15,("unresolved ray",order,index,left,right,u,normalized))
        if du<2*ds+3:
            mid=sum(u)/2
            stack.extend((((u[0],mid),normalized,du+1,ds),((mid,u[1]),normalized,du+1,ds)))
        else:
            mid=sum(normalized)/2
            stack.extend(((u,(normalized[0],mid),du,ds+1),(u,(mid,normalized[1]),du,ds+1)))
    return count


def check_anchor():
    require((1+INNER)**8<2,"annulus dominance")
    for row in RAYS.values():
        require(all(1000>a>b>-1000 for a,b in zip(row,row[1:])),"strict ray order")
    a,b,c=nominal(7,Q(1,10))[1:4]
    d,e,f,g=nominal(8,Q(1,10))[1:5]
    def seven(x):
        return x**3-x.square()-scalar(2)*x+ONE+a*(x.square()-x-ONE)+b*(x-ONE)+c
    def eight(x):
        return x**4-scalar(4)*x.square()+scalar(2)+d*(x**3-scalar(3)*x)+e*(x.square()-scalar(2))+f*x+g
    points=(scalar(2),)+tuple(scalar(x,500) for x in RAYS[100])+(scalar(-2),)
    for evaluator,signs in ((seven,(1,1,-1,-1,1,1,-1,-1)),(eight,(1,-1,-1,1,1,-1,-1,1))):
        for x,sign in zip(points,signs,strict=True):
            require((scalar(sign)*evaluator(x)).lo>0,("Gaussian anchor",evaluator.__name__,x.shown()))


def main():
    validate_skew()
    check_anchor()
    require(RAYS[720]==FINISH,"matching bridge rays at .72")
    selection=tuple(map(int,argv[1:])) if len(argv)>1 else ()
    for a,b,cap in CAPS:
        if selection and a not in selection:
            continue
        env={j:envelope(cap,j) for j in (7,8)}
        limit=env[7][1][0][1].hi
        lower_bound=I(limit,limit)
        knots=[k for k in RAYS if a<=k<=b]
        for left,right in zip(knots,knots[1:]):
            for j in (7,8):
                for i in range(6):
                    print(left,right,cap,j,i,verify(j,i,left,right,env[j],lower_bound),flush=True)


if __name__=="__main__":
    main()
