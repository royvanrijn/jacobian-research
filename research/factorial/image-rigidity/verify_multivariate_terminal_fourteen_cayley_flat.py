"""Exact grouped near-flat Cayley rays at terminal fourteen."""
from fractions import Fraction as Q
from sys import argv

from verify_multivariate_terminal_thirteen_cayley import (
    GRID, I, Limits, ONE, ZERO, accurate_projection, bisect, interval,
    log_bounds, model, model_eval, rational, require, ray_for, scalar,
    skew_model, strength_factor, transform_family,
)
from verify_multivariate_terminal_twelve_cayley_lower import skew_coeff_terms
from verify_multivariate_terminal_nine_shared_matched_light import exp_small

ORDERS=LOWER,UPPER=(13,14)
R=Q(169,168)
DISK=Q(3)
SCALE=100000
RAYS = {
    0: (119067, 107359, 95022, 85088, 74597, 65511, 55889, 47286, 38149, 29832, 20965, 12801),
    5: (119445, 107659, 95242, 85263, 74726, 65612, 55963, 47343, 38189, 29861, 20984, 12812),
    10: (119827, 107962, 95463, 85439, 74856, 65714, 56037, 47400, 38229, 29890, 21002, 12823),
    15: (120212, 108268, 95686, 85617, 74987, 65817, 56112, 47457, 38269, 29919, 21021, 12834),
    20: (120600, 108576, 95910, 85795, 75118, 65920, 56187, 47515, 38310, 29948, 21039, 12845),
}

EXPECTED_CELLS = {
    (13, 0): (866, 91, 1312, 171, 785, 156, 415, 97, 99, 22, 7, 1),
    (14, 0): (1391, 149, 5513, 455, 5165, 410, 586, 145, 140, 30, 18, 1),
    (13, 5): (1012, 67, 840, 246, 764, 217, 318, 66, 71, 26, 16, 5),
    (14, 5): (2095, 214, 1498, 445, 1177, 313, 716, 204, 196, 48, 27, 7),
    (13, 10): (1286, 146, 996, 185, 685, 172, 467, 110, 146, 42, 10, 10),
    (14, 10): (2378, 175, 2619, 688, 1944, 594, 664, 172, 172, 42, 42, 10),
    (13, 15): (1707, 128, 1183, 375, 921, 138, 419, 128, 128, 32, 29, 8),
    (14, 15): (5147, 244, 1391, 512, 1500, 512, 677, 155, 252, 43, 54, 13),
}


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
    return Limits(c,(r_lo,r_hi),qmax,DISK,True)

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

def radial_exponent(limits):
    c=limits.cap
    cu=upper_reciprocal(c)
    return 1+max((LOWER-1)*c/(2+(LOWER-1)*c),LOWER*cu/(2+LOWER*cu))

def anchor_for(limits,rays,first):
    p = rational(Q(first,1000))
    ratio = interval(*limits.ratio)
    for n in ORDERS:
        pending = [(ratio,0)]
        while pending:
            r,depth = pending.pop()
            root_r = r.sqrt() if n==UPPER else ONE
            h = p*r if n==UPPER else p
            coeff = tuple(v for _,v in model_eval(transform_family(n,0,(Q(1),)*(n+1)),h.sqrt(),h))
            levels = (rational(limits.disk),*(scalar(z,SCALE)/root_r for z in rays[first]),ZERO)
            polynomial = tuple(scalar((-1)**((h-n%2)//2))*coeff[h] for h in range(n%2,n+1,2))
            prior = 0
            for i,y in enumerate(levels):
                x = y.square()
                xm = (x.lo+x.hi)//2
                center = I(xm,xm)
                value = ZERO
                for term in reversed(polynomial):
                    value = value*center+term
                slope = ZERO
                for k in range(len(polynomial)-1,0,-1):
                    slope = slope*x+scalar(k)*polynomial[k]
                value += slope*(x-center)
                sign = 1 if value.lo>0 else -1 if value.hi<0 else 0
                if not sign or (i and (sign!=prior)!=(n==UPPER if i%2 else n==LOWER)):
                    break
                prior = sign
            else:
                continue
            require(depth<16 and n==UPPER, ("reference anchor signs", n,r.shown(),i,value.shown()))
            pending.extend((half,depth+1) for half in bisect(r))

def sharp_projection(order,index,left,right,lim,rays,window=None):
    ratio=interval(*lim.ratio)
    lo,hi=window or (left,right)
    lam=interval(Q(lo,1000),Q(hi,1000)).sqrt()
    A=moment_ratio(order,lo,hi,lim)
    C=strength_factor(lim.qmax)
    alpha=rational(radial_exponent(lim))
    skew_radius=C*A/2
    pending=[(lam,interval(Q(0),DISK),ratio,None,0)]
    count=0
    while pending:
        lam,a,r,sk,depth=pending.pop()
        sr=r.sqrt() if order==UPPER else ONE
        scale=lam*sr
        p,q=lam.square(),scale.square()
        cutoff=(alpha*scale/(scalar(2)-alpha*q)).hi
        if a.lo>cutoff:
            continue
        a=I(a.lo,min(a.hi,cutoff))
        b=ray_for(p,index,left,right,rays)/sr
        denom=ONE+q.square()*a.square()*b.square()
        x=a*(ONE+q*b.square())/denom
        y=b*(ONE-q*a.square())/denom
        if sk is None:
            base,groups,box,ok=model(order,scale.lo,scale.hi,A,lim,index<3)
        else:
            base,groups,box,ok=skew_model(order,scale.lo,scale.hi,A,lim,index<3,*sk)
        if not ok:
            require(depth<34,("scaled disk",order,index,left,right,lam.shown(),r.shown(),sk))
            pending.extend((half,a,r,sk,depth+1) for half in bisect(lam))
            continue
        if accurate_projection(order,x,y,base,groups,box):
            count+=1
            continue
        if sk is None and index%2==0 and depth>=11 and a.hi<=GRID//5:
            pending.extend((lam,a,r,(skew_radius*Q(k,4),skew_radius*Q(k+1,4)),depth)
                           for k in range(-4,4))
            continue
        require(depth<44,("unresolved",order,index,left,right,lam.shown(),a.shown(),r.shown(),sk,count))
        wl,wa,wr = 10*(lam.hi-lam.lo),a.hi-a.lo,10*(r.hi-r.lo) if order==UPPER else 0
        if wr>max(wl,wa):
            pending.extend((lam,a,half,sk,depth+1) for half in bisect(r))
        elif wl>wa:
            pending.extend((half,a,r,sk,depth+1) for half in bisect(lam))
        else:
            pending.extend((lam,half,r,sk,depth+1) for half in bisect(a))
    return count


def validate_rays(limits,rays,left,right):
    for endpoint in (left,right):
        a=(int(DISK*SCALE),*rays[endpoint],0)
        require(len(a)==UPPER and all(x>y for x,y in zip(a,a[1:])),"ray order")
        require(Q(a[1],SCALE)**2*limits.qmax<1,"physical ray")

def run(left, right, order, index):
    require(left in RAYS and right in RAYS and left<right, "slab range")
    require(order in (0, *ORDERS) and -1<=index<UPPER-2, "order or ray")
    anchor_for(slab_limits(None,0,0),RAYS,0)
    total=0
    for slab in range(left,right,5):
        end=slab+5
        lim=slab_limits(None,slab,end)
        validate_rays(lim,RAYS,slab,end)
        for j in ORDERS:
            if order and j!=order:
                continue
            for k in range(UPPER-2):
                if index>=0 and k!=index:
                    continue
                if k==0:
                    cells=0
                    for sub in range(slab,end):
                        inner=slab_limits(None,sub,sub+1,refined=True)
                        cells+=sharp_projection(j,k,slab,end,inner,RAYS,(sub,sub+1))
                else:
                    cells=sharp_projection(j,k,slab,end,lim,RAYS)
                require(cells==EXPECTED_CELLS[j,slab][k],("replay count",j,k,slab,cells))
                total+=cells
                print(j,k,slab,end,cells,flush=True)
        model.cache_clear()
        skew_model.cache_clear()
        skew_coeff_terms.cache_clear()
    print(f"Verified {total} projection cells.",flush=True)

if __name__=="__main__":
    run(*(tuple(map(int,argv[1:])) if len(argv)>1 else (0,20,0,-1)))
