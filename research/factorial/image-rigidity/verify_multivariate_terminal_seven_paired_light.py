#!/usr/bin/env python3
"""Exact paired-ray separation for light shared rows at terminal seven."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from sys import argv

from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, interval, log, rational, require, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import power_bounds


RAY_KNOTS = {
    100: (789, 583, -59, -419, -813),
    120: (759, 540, -93, -428, -824),
    140: (727, 510, -121, -437, -831),
    160: (697, 465, -145, -446, -837),
    180: (666, 434, -172, -524, -847),
    200: (638, 388, -198, -541, -853),
    220: (608, 356, -224, -547, -858),
    240: (578, 310, -253, -554, -864),
    260: (548, 278, -278, -561, -872),
    280: (516, 246, -302, -568, -877),
    300: (483, 199, -326, -583, -882),
    320: (442, 126, -349, -590, -889),
    340: (401, 94, -372, -598, -894),
    360: (359, 75, -395, -613, -898),
    380: (315, 68, -417, -621, -902),
    400: (280, 74, -442, -635, -906),
    420: (248, -110, -464, -643, -910),
    440: (212, -92, -482, -671, -912),
    460: (164, -101, -506, -678, -916),
    480: (115, -122, -527, -686, -920),
    500: (75, -143, -544, -717, -921),
    520: (28, -177, -565, -730, -925),
    540: (-17, -211, -585, -743, -929),
    560: (-53, -234, -607, -750, -933),
    580: (-91, -268, -629, -763, -936),
    600: (-129, -293, -648, -776, -940),
    620: (-170, -346, -669, -788, -942),
    640: (-208, -388, -690, -800, -945),
    660: (-246, -430, -708, -809, -949),
    680: (-284, -455, -730, -821, -951),
    700: (-326, -475, -748, -814, -954),
    720: (-373, -557, -766, -824, -957),
    740: (-416, -550, -782, -856, -961),
    750: (-439, -561, -791, -862, -962),
    760: (-461, -567, -798, -866, -964),
    780: (-506, -597, -817, -877, -967),
    800: (-549, -633, -835, -889, -970),
    820: (-592, -670, -853, -901, -973),
    840: (-638, -698, -870, -912, -976),
    850: (-662, -736, -878, -918, -977),
    860: (-686, -753, -887, -923, -978),
    880: (-732, -788, -903, -933, -981),
}
BANDS = (
    (100, 420, Q(1, 10)),
    (420, 480, Q(3, 40)),
    (480, 620, Q(1, 20)),
    (620, 720, Q(1, 16)),
    (720, 750, Q(7, 200)),
    (750, 880, Q(1, 40)),
)
ANCHOR_RATIO = scalar(36,35)
INNER, OUTER = Q(1,10), Q(10)


@lru_cache(maxsize=None)
def exponents(order):
    return tuple(Q(k*(order-k),order-1) for k in range(order+1))


@lru_cache(maxsize=None)
def nominal(order, mean):
    multiplier = Q(1) if order==6 else Q(36,35)
    return tuple(
        scalar(comb(order,k))*interval(*power_bounds(mean,multiplier*e)) if e else ONE
        for k,e in enumerate(exponents(order))
    )


def envelope(t, order):
    t0=rational(t)
    lower=log(ONE+scalar(36,7)*t0)/log(ONE+scalar(5)*t0)
    upper=ONE-log(ONE-t0*scalar(1,7))/log(ONE+scalar(5)*t0)
    rho=I(lower.lo,upper.hi)
    require(rho.lo<ANCHOR_RATIO.lo<=ANCHOR_RATIO.hi<rho.hi,"outer ratio anchor")
    t=t0 if order==6 else scalar(6)*t0/(scalar(7)-t0)
    s=scalar(1,2)*log(ONE+scalar(order-1)*t)
    if order==6:
        a,b,c,d,e=(ONE+scalar(k)*t for k in range(1,6))
        outer=log(a*b*c*d/e.square())/(scalar(6)*s)
        linked=scalar(6)*t**4/s
        inner=(scalar(4)*outer+linked)*scalar(1,5)
        defect=scalar(3,5)-scalar(1,2)*log(d/a)/s
        J=scalar(6,5)*t**5/s
        groups=(
            ((0,Q(5,4),1,0,-1,Q(-5,4),0),inner.symmetric_error()),
            ((0,1,0,0,0,-1,0),(linked*scalar(1,4)).symmetric_error()),
            ((0,0,1,Q(3,2),1,0,0),defect.positive_error()),
            ((0,0,0,-1,0,0,0),J.positive_error()),
        )
    else:
        inner=t.square()/s
        L=scalar(3)*t**4/s
        K=scalar(60,7)*t**6/s
        defect=scalar(2,3)-scalar(1,2)*log((ONE+scalar(5)*t)/(ONE+t))/s
        E=s**4*scalar(1,25)
        groups=(
            ((0,Q(5,2),Q(5,2),1,-1,Q(-5,2),Q(-5,2),0),inner.symmetric_error()),
            ((0,Q(3,2),Q(1,2),0,0,Q(-1,2),Q(-3,2),0),L.symmetric_error()),
            ((0,Q(-1,2),Q(-1,2),0,0,Q(1,2),Q(1,2),0),K.symmetric_error()),
            ((0,0,1,Q(9,5),Q(9,5),1,0,0),defect.positive_error()),
            ((0,0,0,-1,-1,0,0,0),E.positive_error()),
        )
    return rho, tuple((tuple(rational(x) for x in vector),v) for vector,v in groups)


def coefficient_model(order,u,env):
    rho, variables=env
    mid=(u[0]+u[1])/2
    reference=nominal(order,mid)
    U=-log(interval(*u))
    mean=log(interval(u[0]/mid,u[1]/mid))
    if order==7:
        mean=ANCHOR_RATIO*mean+(ANCHOR_RATIO-rho)*U
        U=rho*U
    pieces=tuple((p*U,vector) for vector,p in variables)
    pieces+=((mean,tuple(rational(e) for e in exponents(order))),)
    data=[]
    for k,a in enumerate(reference):
        radius=ZERO
        for parameter,vector in pieces:
            term=parameter*vector[k]
            radius=radius+I(0,max(-term.lo,term.hi))
        require(radius.hi<GRID,("exponential remainder",order,u))
        error=scalar(1,2)*radius.square()/(ONE-radius)
        data.append(a*(ONE+error))
    return ((ONE,tuple(data)),)+tuple(
        (parameter,tuple(a*v for a,v in zip(reference,vector,strict=True)))
        for parameter,vector in pieces
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


def projected(order,cosines,model):
    sines=[ZERO,ONE]
    for k in range(1,order):
        sines.append(scalar(2)*cosines*sines[k]-sines[k-1])
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


def signed(projection,model):
    negative=positive=True
    for column in zip(*projection,strict=True):
        lo=hi=0
        for (a,_),b in zip(model,column,strict=True):
            vals=(a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi)
            lo+=min(vals)
            hi+=max(vals)
        negative=negative and hi<0
        positive=positive and lo>0
        if not (negative or positive):
            return False
    return True


def halves(values):
    left,right=[],[]
    while values:
        left.append(values[0])
        right.append(values[-1])
        values=tuple(
            I((a.lo+b.lo)//2,-(-(a.hi+b.hi)//2))
            for a,b in zip(values,values[1:])
        )
    return tuple(left),tuple(reversed(right))


def radial(order,cosines,model):
    stack=[(projected(order,cosines,model),0)]
    used=0
    while stack:
        polys,depth=stack.pop()
        if any(signed(candidate,model) for candidate in polys):
            used+=1
            continue
        if depth>=14:
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
    a,b=RAY_KNOTS[left][index],RAY_KNOTS[right][index]
    first=Q(a,1000)+Q(b-a,right-left)*(u[0]-Q(left,1000))
    last=Q(a,1000)+Q(b-a,right-left)*(u[1]-Q(left,1000))
    return interval(min(first,last),max(first,last))


def verify(order,index,left,right,t):
    env=envelope(t,order)
    pending=[((Q(left,1000),Q(right,1000)),0)]
    found=0
    while pending:
        u,depth=pending.pop()
        used=radial(order,ray_cosines(left,right,index,u),coefficient_model(order,u,env))
        if used is not None:
            found+=used
            continue
        require(depth<15,("unresolved ray interval",order,index,left,right,t,[float(x) for x in u]))
        middle=(u[0]+u[1])/2
        pending.extend((((u[0],middle),depth+1),((middle,u[1]),depth+1)))
    return found


def check_anchor():
    require((1+INNER)**7<2,"root annulus")
    for rays in RAY_KNOTS.values():
        require(all(1000>a>b>-1000 for a,b in zip(rays,rays[1:])),"ordered rays")
    u=Q(1,10)
    a,b,c=nominal(6,u)[1:4]
    d,e,f=nominal(7,u)[1:4]
    def even(x):
        return x**3-scalar(3)*x+a*(x.square()-scalar(2))+b*x+c
    def odd(x):
        return x**3-x.square()-scalar(2)*x+ONE+d*(x.square()-x-ONE)+e*(x-ONE)+f
    for test,signs in ((even,(1,1,-1,-1,1,1,-1)),(odd,(1,-1,-1,1,1,-1,-1))):
        abscissae=(scalar(2),)+tuple(scalar(x,500) for x in RAY_KNOTS[100])+(scalar(-2),)
        for x,sign in zip(abscissae,signs,strict=True):
            require((scalar(sign)*test(x)).lo>0,("Gaussian anchor",test.__name__,x.shown()))


def main():
    check_anchor()
    for a,b,t in BANDS:
        if len(argv)>1 and a not in map(int,argv[1:]):
            continue
        knots=[k for k in RAY_KNOTS if a<=k<=b]
        for left,right in zip(knots,knots[1:]):
            counts=[verify(order,ray,left,right,t) for order in (6,7) for ray in range(5)]
            print(left,right,t,counts,flush=True)


if __name__=="__main__":
    main()
