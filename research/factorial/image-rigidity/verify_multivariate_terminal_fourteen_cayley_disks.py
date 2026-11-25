"""Direct exact grouped Cayley root disks for terminal fourteen."""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial, prod
from sys import argv
from typing import NamedTuple

from flint import acb, acb_poly, arb, ctx

from verify_multivariate_terminal_twelve_cayley import (
    GRID, I, ONE, cayley_matrix, interval, log_bounds, log_shape, model_derivative,
    model_eval, rational, require, strength_factor, transform_family,
)
from verify_multivariate_terminal_nine_shared_matched_light import exp_small

ORDERS = LOWER, UPPER = (13, 14)
R = Q(169, 168)
DISK = Q(3)
WEIGHT = 11
EXPECTED_CELLS = {20: 350, 25: 318, 30: 318, 35: 339, 40: 364, 45: 417, 50: 484, 55: 498, 60: 352, 65: 312, 70: 265, 75: 216, 80: 337, 85: 318, 90: 264, 95: 240, 100: 240, 105: 240, 110: 241, 115: 297, 120: 369, 125: 405, 130: 503, 135: 668, 140: 732, 145: 821, 150: 886, 155: 1188, 160: 1543, 165: 1609, 170: 914, 175: 790, 180: 649, 185: 482, 190: 366, 195: 321, 200: 279, 205: 248, 210: 240, 215: 240, 220: 240, 225: 240, 230: 243, 235: 240}
ctx.prec = 112
RADII = tuple(arb(k)/2**e for e in range(22,5,-1) for k in (8,9,10,11,12,13,14,15))

class Limits(NamedTuple):
    cap: Q
    ratio: tuple
    qmax: Q
    disk: Q

def power(u,k):
    return exp_small(-rational(k)*interval(*log_bounds(1/u)))

def ceil_on_grid(value,scale):
    return Q(-(-value.numerator*scale//value.denominator),scale)

def lattice_cap(u,c):
    if not c:
        return c
    ell=max(1,-(-(c.denominator-c.numerator)//(LOWER*c.numerator)))
    return max(Q(1,LOWER*ell+1),min(c,(1-u)/(LOWER*u-1)))

def slab_limits(cap,left,right,refined=False):
    p_lo,p_hi=Q(left,1000),Q(right,1000)
    u=1-p_hi
    c=(u**-2-1)/(LOWER-1)
    if cap is not None:
        c=min(c,cap)
    if refined:
        c=lattice_cap(u,c)
    rho=1+log_bounds(1/(1-c/UPPER))[1]/log_bounds(1+(LOWER-1)*c)[0] if c else R
    require(R<=rho<2,("ratio upper",rho))
    r_lo=(1-Q(power(u,R).hi,GRID))/p_hi if p_hi else R
    r_hi=(1-Q(power(1-p_lo,rho).lo,GRID))/p_lo if p_lo else rho
    qmax=ceil_on_grid(p_hi*r_hi,10000)+Q(1,10000)
    require(Q((rational(p_hi).sqrt()*rational(r_hi).sqrt()).square().hi,GRID)<qmax,"deficit enclosure")
    return Limits(c,(r_lo,r_hi),qmax,DISK)

def upper_reciprocal(lower_cap):
    return LOWER*lower_cap/(UPPER-lower_cap)

def moment_ratio(n,left,right,limits):
    c=limits.cap
    p_lo=Q(left,1000)
    if n==UPPER:
        c=upper_reciprocal(c)
        p_lo*=limits.ratio[0]
    t=(n-1)*c/(2+(n-1)*c)
    c0=strength_factor(limits.qmax)
    return min(c0,t/p_lo) if p_lo else c0

def directions(n):
    return {d:log_shape(n,d) for d in range(2,WEIGHT+1)}

def high_exponential_weights(heights, minimum, first_unretained_power):
    """Sum the high-weight portion through the retained exponential powers."""
    coefficients = [Q(0)] * minimum
    coefficients[0] = Q(1)
    tail = total = Q(0)
    log_mass = sum(heights.values())
    positive = tuple((d,h) for d,h in heights.items() if h)
    for exponent in range(1, first_unretained_power):
        next_coefficients = [Q(0)] * minimum
        discarded = tail * log_mass
        for weight, coefficient in enumerate(coefficients):
            if coefficient:
                for d, h in positive:
                    product = coefficient * h
                    if weight+d<minimum:
                        next_coefficients[weight+d] += product
                    else:
                        discarded += product
        coefficients = [coefficient/exponent for coefficient in next_coefficients]
        tail = discarded/exponent
        total += tail
    return total


@lru_cache(maxsize=None)
def raw_error(n, a, qmax):
    c = strength_factor(qmax)
    bounds = []
    for k in range(n+1):
        r = min(k,n-k)
        omega = sum(abs((int(r<=h<n-r)-Q(n-2*r,n))/2) for h in range(n))
        m = WEIGHT+1
        log_tail = (omega*c*a**(m-1)/m+r*qmax*c*a**m/(m+1))/(1-(qmax*a)**2)
        heights = {d:abs(vec[k])*c*a**(d-1)*qmax**d/d for d,vec in directions(n).items()}
        retained_height = sum(heights.values())
        full_height = retained_height+log_tail*qmax**m
        require(full_height<1, ("small log", n, k))
        d0=WEIGHT//2+1
        exp_tail = retained_height**d0/(factorial(d0)*(1-retained_height))
        exp_tail += high_exponential_weights(heights,m,d0)
        bounds.append(log_tail/(1-full_height)+exp_tail/qmax**m)
    return tuple(bounds)

@lru_cache(maxsize=None)
def skew_coeff_terms(n,a,c,skew_lo,skew_hi):
    source=directions(n)
    bounds={2:interval(skew_lo,skew_hi)}
    signed_abs=min(abs(skew_lo),abs(skew_hi)) if skew_lo*skew_hi>0 else Q(0)
    positive=(c*c+2*skew_hi)/(2*c)
    negative=(c*c-2*skew_lo)/(2*c)
    whole=c//a
    rem=c-whole*a
    for d in range(3,WEIGHT+1):
        positive_bound=min(a**(d-2)*positive**2,positive**d)
        negative_bound=min(a**(d-2)*negative**2,negative**d)
        majorized=whole*a**d+rem**d
        if d%2:
            low=(2*signed_abs)**(d-1)/(c**(d-2)*d)
            high=min(majorized,positive_bound+negative_bound)/d
            bounds[d]=interval(low,high)
        else:
            bounds[d]=interval(-min(majorized,negative_bound)/d,min(majorized,positive_bound)/d)
    result=[]
    for weight,direction,counts in monomial_directions(n):
        parameter=prod((bounds[d].square()**(j//2)*(bounds[d] if j%2 else rational(Q(1)))
                        for d,j in counts),start=rational(Q(1)))
        result.append((weight,direction,parameter))
    return tuple(result)


@lru_cache(maxsize=None)
def monomial_directions(n):
    source=directions(n)
    output=[]
    for degree in range(1,WEIGHT//2+1):
        for indices in combinations_with_replacement(source,degree):
            weight=sum(indices)
            if weight>WEIGHT:
                continue
            counts=tuple(Counter(indices).items())
            denominator=prod(factorial(j) for _,j in counts)
            direction=tuple(prod(source[d][k] for d in indices)/Q(denominator) for k in range(n+1))
            output.append((weight,direction,counts))
    return tuple(output)

@lru_cache(maxsize=8192)
def skew_model(n,lo,hi,a,limits,linearized,skew_lo,skew_hi):
    lam=I(lo,hi)
    q=lam.square()
    family=transform_family(n,0,(Q(1),)*(n+1))
    if linearized:
        mid=I((lo+hi)//2,(lo+hi)//2)
        core=model_eval(family,mid,mid.square())
        variations=((lam-mid,model_derivative(family,lam,q)),)
    else:
        core=model_eval(family,lam,q)
        variations=()
    groups=variations+tuple((mu, model_eval(transform_family(n,m,vec),lam,q))
                     for m,vec,mu in skew_coeff_terms(n,a,strength_factor(limits.qmax),skew_lo,skew_hi))
    bounds=raw_error(n,a,limits.qmax)
    base=[]
    for (h,v),row in zip(core,cayley_matrix(n),strict=True):
        err=sum(Q(abs(x)*comb(n,k),2**n)*bounds[k] for k,x in enumerate(row))
        base.append((h,v+(rational(err)*lam**(2*(WEIGHT+1)+h-n)).symmetric_error()))
    box=[]
    for h,b in base:
        for mu,vec in groups:
            b+=mu*vec[h][1]
        box.append((h,b))
    r=limits.disk
    ok=box[-1][1].lo*r**n>sum(max(abs(v.lo),abs(v.hi))*r**h for h,v in box[:-1])
    return tuple(base),groups,tuple(box),ok


@lru_cache(maxsize=4096)
def ball(lo,hi):
    return arb(lo+hi,hi-lo)/(2*GRID)

def outward(i):
    return ball(i.lo,i.hi)

def mid(i):
    return arb(i.lo+i.hi)/(2*GRID)

def groups_for(basis,groups):
    data=[(arb(1),acb_poly([outward(v) for _,v in basis]))]
    reference=[mid(v) for _,v in basis]
    for weight,vec in groups:
        if weight.lo or weight.hi:
            data.append((outward(weight),acb_poly([outward(v) for _,v in vec])))
            for k,v in vec:
                reference[k]+=mid(weight)*mid(v)
    return tuple(reference),tuple(data)

def radius_for(order,center,groups):
    shift=acb_poly([center,1])
    coeff=sum((w*poly(shift) for w,poly in groups),acb_poly())
    magnitude=tuple(abs(coeff[m]).upper() for m in range(order+1))
    linear=abs(coeff[1]).lower()
    for radius in RADII:
        lhs=linear*radius
        if lhs>magnitude[0] and lhs>magnitude[0]+sum(radius**k*magnitude[k] for k in range(2,order+1)):
            return radius
    return None

def skew_cap(C,A):
    whole=C//A
    return (whole*A**2+(C-whole*A)**2)/2

def transport(limits,left,right):
    c=limits.cap
    T=12*c/(2+12*c)
    C=strength_factor(limits.qmax)
    q_ratio=limits.ratio[0]
    basic=13**4*C**3*Q(right,1000)/((168-T)**3*q_ratio**2)
    if not left or not c:
        return basic
    B=168*c/(26+168*c)
    E=12*(1+14*c)/(26+168*c)
    rho=1+log_bounds(1/(1-c/14))[1]/log_bounds(1+12*c)[0]
    U=log_bounds(1/(1-Q(right,1000)))[1]
    q0=Q(left,1000)*q_ratio
    excess=R**2*Q(56,12)*B*E**2*(rho-R)*U/q0**2
    return min(basic,excess)

def verify_cell(p0,p1,lim,skew):
    lam=interval(Q(p0,1000),Q(p1,1000)).sqrt()
    ratio=interval(*lim.ratio)
    T=strength_factor(lim.qmax)
    A=moment_ratio(LOWER,p0,p1,lim)
    lower_radius=skew_cap(T,A)
    P=interval(skew[0]*lower_radius,skew[1]*lower_radius)
    err=transport(lim,p0,p1)
    disks=[]
    for order in ORDERS:
        sr=ratio.sqrt() if order==UPPER else ONE
        scale=lam*sr
        h=scale.square()
        A=moment_ratio(order,p0,p1,lim)
        signed= P*(rational(R)/ratio).square()+interval(-err,err) if order==UPPER else P
        radius=skew_cap(T,A)
        blo,bhi=max(-radius,Q(signed.lo,GRID)),min(radius,Q(signed.hi,GRID))
        if blo>bhi:
            return True
        basis,groups,_,ok=skew_model(order,scale.lo,scale.hi,A,lim,True,blo,bhi)
        if not ok:
            return False
        reference,groups=groups_for(basis,groups)
        try:
            roots=acb_poly(list(reference)).roots(tol=arb(2)**-85,maxprec=512)
        except (ArithmeticError,ValueError):
            return False
        plus=[z for z in roots if z.imag>0]
        minus=[z for z in roots if z.imag<0]
        real=[z for z in roots if z.imag.contains(0)]
        if len(plus)!=order//2 or len(minus)!=order//2 or len(real)!=order%2:
            return False
        for z in real:
            if radius_for(order,acb(z.real.mid()),groups) is None:
                return False
        for z in plus:
            center=acb(z.real.mid(),z.imag.mid())
            radius=radius_for(order,center,groups)
            if radius is None or not center.imag>radius:
                return False
            radii=abs(center)+arb(0,radius)
            denom=1-outward(h)*radii**2
            if not denom>0:
                return False
            metric=outward(sr)*(center.imag+arb(0,radius))/denom
            disks.append(metric)
    disks.sort(key=lambda value:float(value.mid()))
    return all(a.upper()<b.lower() for a,b in zip(disks,disks[1:]))

def replay(cap,left,right):
    pending=[(Q(left),Q(right),Q(0),Q(1),0,0)]
    cells=0
    while pending:
        a,b,v,w,dp,ds=pending.pop()
        lim=slab_limits(cap,a,b,refined=True)
        if verify_cell(a,b,lim,(v,w)):
            cells+=1
        else:
            if dp+ds>=23:
                raise ArithmeticError(("unresolved",float(a),float(b),float(v),float(w),cells))
            if dp<ds+1:
                midp=(a+b)/2
                pending.extend(((a,midp,v,w,dp+1,ds),(midp,b,v,w,dp+1,ds)))
            else:
                mids=(v+w)/2
                pending.extend(((a,b,v,mids,dp,ds+1),(a,b,mids,w,dp,ds+1)))
    return cells


def main():
    left,right=map(int,argv[1:]) if len(argv)>1 else (20,240)
    require(left in EXPECTED_CELLS and right in (*EXPECTED_CELLS,240) and left<right,"slab range")
    total=0
    for slab in range(left,right,5):
        cap=None if slab<50 else Q(1,100) if slab<80 else Q(1,80) if slab<120 else Q(1,70) if slab<145 else Q(1,60)
        count=replay(cap,slab,slab+5)
        require(count==EXPECTED_CELLS[slab],("replay count",slab,count))
        total+=count
        print(slab,slab+5,count,flush=True)
        skew_model.cache_clear()
        skew_coeff_terms.cache_clear()
        raw_error.cache_clear()
    print(f"Verified {total} outward interval cells.",flush=True)

if __name__=="__main__":
    main()
