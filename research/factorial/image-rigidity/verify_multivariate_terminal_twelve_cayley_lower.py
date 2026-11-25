"""Exact leading-skew Cayley separation on the XII band .70–.76."""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial, prod
from sys import argv

from verify_multivariate_terminal_twelve_cayley import (
    I, ONE, UPPER, WEIGHT, accurate_projection, anchor_for, bisect, cayley_matrix,
    check_ray, directions, height, interval, limits_for, model_derivative,
    model_eval, moment_ratio, rational, raw_error, ray_for, require, strength_factor,
    transform_family, validate_rays,
)


RAYS = {
    240: (138039, 121271, 100994, 87541, 71990, 63883, 47252, 37743, 23450, 12728),
    245: (138612, 121706, 101270, 87748, 72125, 63993, 47316, 37789, 23475, 12740),
    250: (139192, 122147, 101547, 87957, 72261, 64103, 47380, 37836, 23500, 12753),
    255: (139779, 122592, 101827, 88167, 72398, 64214, 47446, 37883, 23525, 12766),
    260: (140372, 123042, 102110, 88378, 72535, 64325, 47511, 37931, 23550, 12779),
    265: (140973, 123496, 102394, 88592, 72674, 64438, 47577, 37979, 23575, 12792),
    270: (141580, 123956, 102835, 88807, 72814, 64551, 47643, 38027, 23601, 12805),
    275: (142275, 124421, 103125, 89024, 72954, 64664, 47473, 38075, 23627, 12427),
    280: (142898, 124890, 103418, 89242, 73096, 64779, 47540, 38123, 23652, 12440),
    285: (143528, 125365, 103714, 88660, 73238, 64894, 47607, 38172, 23399, 12452),
    290: (144166, 125845, 104012, 88878, 73382, 65010, 47674, 38221, 23425, 12465),
    295: (144812, 126331, 104313, 89098, 73526, 65127, 47741, 38271, 23450, 12478),
    300: (145466, 126822, 104617, 89319, 73672, 65244, 47810, 38321, 23476, 12491),
}

EXPECTED_CELLS = {
    (11, 240): (1114, 49, 260, 62, 119, 51, 62, 30, 31, 18),
    (11, 245): (1060, 49, 252, 62, 108, 49, 62, 30, 31, 18),
    (11, 250): (1027, 44, 240, 58, 99, 48, 62, 29, 30, 17),
    (11, 255): (995, 43, 223, 58, 91, 46, 61, 29, 30, 17),
    (11, 260): (963, 43, 219, 58, 91, 44, 57, 30, 30, 17),
    (11, 265): (934, 42, 188, 55, 91, 42, 57, 30, 30, 16),
    (11, 270): (849, 42, 178, 54, 91, 40, 57, 28, 30, 16),
    (11, 275): (779, 42, 169, 50, 90, 40, 55, 26, 30, 15),
    (11, 280): (752, 41, 161, 50, 88, 40, 55, 25, 30, 15),
    (11, 285): (740, 39, 161, 50, 85, 40, 55, 24, 30, 15),
    (11, 290): (735, 38, 161, 50, 81, 40, 55, 24, 30, 15),
    (11, 295): (695, 36, 150, 50, 81, 40, 53, 24, 30, 15),
    (12, 240): (5284, 63, 1691, 175, 461, 111, 118, 41, 34, 21),
    (12, 245): (3988, 62, 1545, 175, 453, 101, 118, 41, 34, 21),
    (12, 250): (3210, 62, 1273, 176, 420, 77, 106, 41, 32, 21),
    (12, 255): (2635, 61, 811, 177, 420, 64, 102, 39, 31, 21),
    (12, 260): (2255, 63, 594, 178, 412, 64, 102, 39, 31, 21),
    (12, 265): (2146, 64, 766, 174, 375, 62, 91, 39, 30, 21),
    (12, 270): (2217, 64, 780, 169, 372, 62, 87, 38, 28, 21),
    (12, 275): (2157, 64, 575, 169, 334, 62, 74, 36, 27, 20),
    (12, 280): (1994, 64, 480, 170, 325, 62, 74, 34, 25, 20),
    (12, 285): (1905, 65, 439, 169, 325, 58, 69, 32, 24, 20),
    (12, 290): (1853, 65, 438, 169, 269, 58, 53, 28, 24, 20),
    (12, 295): (1789, 60, 434, 169, 227, 55, 46, 27, 24, 20),
}


@lru_cache(maxsize=None)
def skew_coeff_terms(n,a,c,skew_lo,skew_hi):
    source=directions(n)
    bounds={2:interval(skew_lo,skew_hi)}
    signed_abs=min(abs(skew_lo),abs(skew_hi)) if skew_lo*skew_hi>0 else Q(0)
    for d in range(3,WEIGHT+1):
        if d%2:
            bounds[d]=interval((2*signed_abs)**(d-1)/(c**(d-2)*d),c*a**(d-1)/d)
        else:
            bounds[d]=interval(-(c*a-2*skew_lo)*a**(d-2)/(2*d),
                                (c*a+2*skew_hi)*a**(d-2)/(2*d))
    out=[]
    for degree in range(1, WEIGHT//2+1):
        for indices in combinations_with_replacement(source,degree):
            weight=sum(indices)
            if weight>WEIGHT:
                continue
            count=Counter(indices)
            den=prod(factorial(j) for j in count.values())
            vec=tuple(prod(source[d][k] for d in indices)/Q(den) for k in range(n+1))
            param=prod((bounds[d].square()**(j//2)*(bounds[d] if j%2 else rational(Q(1)))
                        for d,j in count.items()),start=rational(Q(1)))
            out.append((weight,vec,param))
    return tuple(out)

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


def check_leading(order,left,right,rays,cap=Q(3,100)):
    disk=Q(3)
    limits=limits_for(cap,left,right,disk,True,refined=True)
    validate_rays(limits,rays,left,right)
    ratio=interval(*limits.ratio)
    lam=interval(Q(left,1000),Q(right,1000)).sqrt()
    a=moment_ratio(order,left,right,limits)
    radius=strength_factor(limits.qmax)*a/2
    pending=[(lam,interval(Q(0),disk),ratio,radius*Q(i-4,4),radius*Q(i-3,4),0)
             for i in range(8)]
    count=0
    while pending:
        lam,x,r,skew_lo,skew_hi,depth=pending.pop()
        root_ratio = r.sqrt() if order==UPPER else ONE
        scale=lam*root_ratio
        p,q=lam.square(),scale.square()
        b=ray_for(p,0,left,right,rays)/root_ratio
        base,groups,box,ok=skew_model(order,scale.lo,scale.hi,a,limits,True,skew_lo,skew_hi)
        if not ok:
            require(depth<32,("scaled disk",order,left,right,lam.shown(),skew_lo,skew_hi))
            pending.extend((half,x,r,skew_lo,skew_hi,depth+1) for half in bisect(lam))
            continue
        max_x=(ONE/q.sqrt()).hi
        if x.lo>max_x:
            continue
        x=I(x.lo,min(x.hi,max_x))
        y=height(q,x,b)
        if accurate_projection(order,x,y,base,groups,box):
            count+=1
            continue
        require(depth<44,("leading ray",order,left,right,lam.shown(),x.shown(),r.shown(),skew_lo,skew_hi,count))
        wl,wx,wr=10*(lam.hi-lam.lo),x.hi-x.lo,12*(r.hi-r.lo) if order==UPPER else 0
        if wr>max(wl,wx):
            pending.extend((lam,x,half,skew_lo,skew_hi,depth+1) for half in bisect(r))
        elif wl>wx:
            pending.extend((half,x,r,skew_lo,skew_hi,depth+1) for half in bisect(lam))
        else:
            pending.extend((lam,half,r,skew_lo,skew_hi,depth+1) for half in bisect(x))
    return count


def main():
    order, ray, low, high = map(int,argv[1:]) if len(argv)>1 else (0,-1,240,300)
    require(order in (0,11,12) and -1<=ray<10,"ray selection")
    require(low in RAYS and high in RAYS and low<high,"range")
    require(tuple(RAYS)==tuple(range(240,301,5)),"contiguous ray knots")
    anchor_for(limits_for(Q(3,100),240,240,Q(3),True),RAYS,240)
    total=0
    for left in range(low,high,5):
        right=left+5
        limits=limits_for(Q(3,100),left,right,Q(3),True)
        validate_rays(limits,RAYS,left,right)
        for n in (11,12):
            if order and order!=n:
                continue
            for index in range(10):
                if ray>=0 and ray!=index:
                    continue
                if index==0:
                    count=check_leading(n,left,right,RAYS)
                else:
                    linearized=index==1 or (n==12 and index==2)
                    count=check_ray(n,index,left,right,limits,RAYS,linearized)
                require(count==EXPECTED_CELLS[n,left][index],
                        ("replay cell count",n,left,index,count))
                print(n,left,right,index,count,flush=True)
                total+=count
    print(f"Verified {total} projection cells.",flush=True)


if __name__=="__main__":
    main()
