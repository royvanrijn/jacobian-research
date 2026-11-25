"""Exact stable Cayley certification at terminal eleven for outer mean at least .95."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial, prod
from sys import argv

from verify_multivariate_terminal_seven_upper_bridge import cayley_matrix, height
from verify_multivariate_terminal_six_intermediate_cayley import (
    I, ONE, ZERO, bisect, interval, rational, require, scalar,
)
from verify_multivariate_terminal_ten_near_flat import (
    accurate_projection, log_vectors, model_eval, sb, transform_family,
)

LOWER, UPPER = J = (10,11)
LOWER_DEFICIT=Q(1,20)
LAMBDA_TOP=Q(2236068,10**7)
QMAX=Q(51,1000)
RATIO=(Q(10093,10000),Q(5058,5000))
DISK=Q(2)
SCALE=100000
C=1/(1-QMAX)

RAYS = {
    0: (115496, 108473, 85985, 79764, 60586, 54742, 36907, 30849, 14350),
    5: (116165, 108779, 86471, 79915, 60671, 54187, 37236, 30879, 14362),
    10: (116521, 108080, 86648, 80067, 60981, 54258, 37275, 30909, 14375),
    15: (116799, 107372, 86984, 80219, 61068, 54330, 37314, 30415, 14387),
    20: (117160, 107673, 87165, 80373, 61156, 54401, 37353, 30445, 14037),
    25: (117442, 107977, 87346, 80527, 61244, 54473, 37392, 30474, 14049),
    30: (117809, 107252, 87529, 80683, 61332, 54546, 37431, 30504, 14061),
    35: (118095, 107554, 87713, 80839, 61421, 54619, 37470, 30534, 14073),
    40: (118468, 107859, 87899, 80996, 61511, 54692, 37510, 30564, 14448),
    45: (118759, 108167, 88085, 81155, 61600, 54766, 37550, 30595, 14461),
    50: (119137, 108477, 88273, 81314, 61691, 54840, 37590, 30097, 14473),
}
EXPECTED_CELLS = {
    (10, 0): (431, 258, 595, 377, 505, 352, 255, 156, 54),
    (10, 5): (300, 148, 321, 193, 286, 176, 140, 83, 34),
    (10, 10): (183, 91, 187, 107, 161, 90, 71, 46, 24),
    (10, 15): (238, 96, 191, 108, 170, 93, 71, 46, 28),
    (10, 20): (304, 103, 191, 108, 170, 99, 72, 46, 30),
    (10, 25): (365, 109, 202, 109, 170, 101, 80, 53, 32),
    (10, 30): (442, 118, 213, 109, 170, 102, 83, 63, 36),
    (10, 35): (509, 122, 250, 109, 170, 102, 83, 64, 39),
    (10, 40): (527, 65, 156, 62, 97, 57, 47, 36, 23),
    (10, 45): (719, 65, 161, 62, 103, 57, 62, 39, 23),
    (11, 0): (1490, 460, 1411, 557, 599, 428, 261, 184, 67),
    (11, 5): (2678, 248, 717, 287, 577, 221, 176, 92, 51),
    (11, 10): (1820, 128, 389, 148, 349, 117, 93, 52, 35),
    (11, 15): (3505, 130, 458, 151, 360, 121, 94, 52, 40),
    (11, 20): (6125, 133, 595, 152, 370, 124, 103, 56, 40),
    (11, 25): (10229, 136, 716, 154, 380, 124, 105, 56, 45),
    (11, 30): (14331, 143, 949, 168, 386, 124, 105, 56, 47),
    (11, 35): (17116, 144, 1124, 221, 413, 124, 105, 56, 47),
    (11, 40): (17224, 79, 699, 204, 220, 74, 65, 32, 27),
    (11, 45): (35446, 82, 739, 223, 243, 74, 81, 36, 30),
}

def parameter_bound(d):
    return C**d/d


@lru_cache(maxsize=None)
def terms(n):
    source=log_vectors(n)
    items=[]
    for degree in range(1,4):
        for indices in combinations_with_replacement(source,degree):
            weight=sum(indices)
            if weight>7:
                continue
            occurrences=Counter(indices)
            factor=prod(factorial(a) for a in occurrences.values())
            direction=tuple(prod(source[d][k] for d in indices)/Q(factor) for k in range(n+1))
            bound=prod(parameter_bound(d) for d in indices)
            signed=any(d%2==0 and a%2 for d,a in occurrences.items())
            parameter=sb(bound) if signed else interval(Q(0),bound)
            items.append((weight,direction,parameter))
    return tuple(items)


@lru_cache(maxsize=None)
def raw_error(n):
    source=log_vectors(n)
    bounds=[]
    for k in range(n+1):
        r=min(k,n-k)
        omega=sum(abs((int(r<=h<n-r)-Q(n-2*r,n))/2) for h in range(n))
        log_tail=(omega*C**8/8+r*QMAX*C**9/9)/(1-(QMAX*C)**2)
        h={d:abs(vec[k])*parameter_bound(d)*QMAX**d for d,vec in source.items()}
        height=sum(h.values())
        full_height=height+log_tail*QMAX**8
        require(full_height<1,("small log",n,k))
        exp_tail=height**4/(24*(1-height))
        for degree in range(1,4):
            for indices in combinations_with_replacement(h,degree):
                if sum(indices)<8:
                    continue
                exp_tail+=prod(h[d] for d in indices)/prod(factorial(a) for a in Counter(indices).values())
        bounds.append(log_tail/(1-full_height)+exp_tail/QMAX**8)
    return tuple(bounds)


@lru_cache(maxsize=8192)
def model(n,lo,hi):
    lam=I(lo,hi)
    q=lam.square()
    core=model_eval(transform_family(n,0,(Q(1),)*(n+1)),lam,q)
    groups=tuple((parameter,model_eval(transform_family(n,weight,direction),lam,q))
                 for weight,direction,parameter in terms(n))
    bounds=raw_error(n)
    matrices=cayley_matrix(n)
    base=[]
    for h,v in core:
        err=sum(Q(abs(a)*comb(n,k),2**n)*bounds[k] for k,a in enumerate(matrices[h]))
        base.append((h,v+(rational(err)*lam**(16+h-n)).symmetric_error()))
    box=[]
    for h,b in base:
        for variable,group in groups:
            b+=variable*group[h][1]
        box.append((h,b))
    dominates=box[-1][1].lo*DISK**n>sum(max(abs(a.lo),abs(a.hi))*DISK**h for h,a in box[:-1])
    return tuple(base),groups,tuple(box),dominates


def verify_constants():
    u=1-LOWER_DEFICIT
    lo,hi=RATIO
    # Linked outer ratios from product comparison and derivative of (1+s)^(-1/2).
    a=1+(1-u*u)/(LOWER**2-1)
    require((lo-1)**2*a**3 <=(u*(1+u)/(2*(LOWER**2-1)))**2,"outer ratio lower")
    require(hi>=Q(LOWER**2, LOWER-1)/(UPPER-(u**-2-1)/(LOWER-1)),"outer ratio upper")
    require(LAMBDA_TOP**2>=LOWER_DEFICIT and hi*LAMBDA_TOP**2<=QMAX,"deficit bounds")
    require(DISK**2*QMAX<1,"Cayley disk left half-plane")
    for slot in RAYS.values():
        a=(2*SCALE,*slot,0)
        require(all(x>y for x,y in zip(a,a[1:])), "strict ray order")


def ray_for(p,i,left,right):
    a,b=RAYS[left][i],RAYS[right][i]
    return rational(Q(a,SCALE)) + rational(Q(b-a,SCALE))*(p*scalar(1000)-scalar(left))/scalar(right-left)


def anchor():
    # N circle roots in the common reference scale, bounded over the entire ratio range.
    ratio=interval(*RATIO)
    for n in J:
        coeff=tuple(v for h,v in model_eval(transform_family(n,0,(Q(1),)*(n+1)),ZERO,ZERO))
        b=(scalar(2),*(scalar(u,SCALE)/(ratio.sqrt() if n==UPPER else ONE) for u in RAYS[0]),ZERO)
        prior=0
        for i,y in enumerate(b):
            polynomial=tuple(scalar((-1)**((h-n%2)//2))*coeff[h] for h in range(n%2,n+1,2))
            x=y.square()
            xm=(x.lo+x.hi)//2
            center=I(xm,xm)
            f=ZERO
            for coefficient in reversed(polynomial):
                f=f*center+coefficient
            slope=ZERO
            for k in range(len(polynomial)-1,0,-1):
                slope=slope*x+scalar(k)*polynomial[k]
            result=f+slope*(x-center)
            sign=1 if result.lo>0 else -1 if result.hi<0 else 0
            require(sign!=0,("Hermite anchor sign",n,i,result.shown()))
            if i:
                require((sign!=prior)==(n==UPPER if i%2 else n==LOWER),("Hermite anchor sectors",n,i))
            prior=sign


def check_ray(n,i,left,right):
    ratio=interval(*RATIO)
    lo=interval(Q(left,1000),Q(right,1000)).sqrt()
    pending=[(lo,interval(Q(0),DISK),ratio,0)]
    count=0
    while pending:
        lam,x,r,depth=pending.pop()
        root_ratio=r.sqrt() if n==UPPER else ONE
        scale=lam*root_ratio
        p=lam.square()
        q=scale.square()
        b=ray_for(p,i,left,right)
        if n==UPPER:
            b=b/root_ratio
        base,groups,box,ok=model(n,scale.lo,scale.hi)
        if not ok:
            require(depth<32,("scaled disk",n,i,left,right,lam.shown()))
            pending.extend((half,x,r,depth+1) for half in bisect(lam))
            continue
        y=height(q,x,b)
        if accurate_projection(n,x,y,base,groups,box):
            count+=1
            continue
        require(depth<44,("unresolved",n,i,left,right,lam.shown(),x.shown(),r.shown()))
        wl=10*(lam.hi-lam.lo)
        wx=x.hi-x.lo
        wr=12*(r.hi-r.lo) if n==UPPER else 0
        if wr>max(wl,wx):
            pending.extend((lam,x,half,depth+1) for half in bisect(r))
        elif wl>wx:
            pending.extend((half,x,r,depth+1) for half in bisect(lam))
        else:
            pending.extend((lam,half,r,depth+1) for half in bisect(x))
    return count


def main():
    verify_constants()
    anchor()
    n,i=map(int,argv[1:3]) if len(argv)>1 else (0,-1)
    low,high=map(int,argv[3:]) if len(argv)>3 else (0,max(RAYS))
    require(n in (0,*J) and -1<=i<UPPER-2,"replay ray")
    require(low in RAYS and high in RAYS and low<high,"replay range")
    total=0
    for a,b in zip(RAYS,tuple(RAYS)[1:]):
        if not low<=a<b<=high:
            continue
        for order in J:
            if n and n!=order:
                continue
            for index in range(UPPER-2):
                if i>=0 and i!=index:
                    continue
                count=check_ray(order,index,a,b)
                require(count==EXPECTED_CELLS[order,a][index],
                        ("replay cell count",order,a,b,index,count))
                print(order,a,b,index,count,flush=True)
                total+=count
    print("Certified",total,"projection cells",flush=True)

if __name__=='__main__':
    main()
