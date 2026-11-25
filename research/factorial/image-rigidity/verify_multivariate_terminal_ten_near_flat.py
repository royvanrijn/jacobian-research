"""Exact stable Cayley certification at terminal ten for outer mean at least .97."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial, prod
from sys import argv

from verify_multivariate_terminal_six_marked_row_tail import projected_interval
from verify_multivariate_terminal_seven_upper_bridge import cayley_matrix, height
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, bisect, interval, rational, require, scalar,
)

J=(9,10)
LOWER_DEFICIT=Q(3,100)
LAMBDA_TOP=Q(1732051,10**7)
QMAX=Q(31,1000)
RATIO=(Q(10119,10000),Q(5067,5000))
DISK=Q(2)
SCALE=100000
C=1/(1-QMAX)
TERMS=12

RAYS={
    0: (113946, 106334, 81994, 75199, 54354, 47275, 28422, 21767),
    5: (114755, 106626, 82336, 75333, 54424, 47331, 28449, 21786),
    10: (115099, 105845, 82681, 75466, 54494, 47388, 28475, 21806),
    15: (115447, 106135, 82842, 74784, 54565, 47444, 28502, 21286),
    20: (115701, 105340, 83191, 74916, 54906, 47501, 28529, 21305),
    25: (116054, 105629, 83355, 75050, 54978, 46898, 28556, 21324),
    30: (116312, 105920, 83519, 75184, 55050, 46954, 28583, 21343),
}
EXPECTED_CELLS = {
    (9, 0): (233, 134, 286, 214, 262, 163, 102, 73),
    (9, 5): (125, 77, 143, 115, 131, 83, 61, 52),
    (9, 10): (93, 45, 86, 62, 67, 44, 32, 30),
    (9, 15): (108, 50, 87, 50, 67, 44, 33, 32),
    (9, 20): (120, 50, 96, 54, 67, 44, 35, 34),
    (9, 25): (143, 50, 112, 59, 67, 44, 35, 34),
    (10, 0): (684, 190, 428, 355, 333, 200, 140, 96),
    (10, 5): (599, 103, 302, 188, 176, 105, 73, 49),
    (10, 10): (1067, 77, 230, 107, 89, 54, 43, 26),
    (10, 15): (1388, 90, 230, 110, 93, 54, 44, 26),
    (10, 20): (1724, 96, 233, 112, 101, 56, 45, 26),
    (10, 25): (1961, 102, 258, 112, 101, 54, 45, 26),
}



def sb(k):
    return interval(-k,k)


def theta_powers(n):
    return tuple(Q(1)-Q(2*h,n-1) for h in range(n))


@lru_cache(maxsize=None)
def log_shape(n,d):
    p=theta_powers(n)
    answer=[]
    for k in range(n+1):
        r=min(k,n-k)
        if d%2:
            value=sum(t-t**d for t in p[:r])
        else:
            omega=tuple((int(r<=h<n-r)-Q(n-2*r,n))/2 for h in range(n))
            value=sum(a*t**d for a,t in zip(omega,p,strict=True))
            value=-value if 2*k<=n else value
        answer.append(value)
    return tuple(answer)


def parameter_bound(d):
    return C**d/d


@lru_cache(maxsize=None)
def log_vectors(n):
    return {d:log_shape(n,d) for d in range(2,8)}


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


@lru_cache(maxsize=None)
def series(e):
    values=[Q(1)]
    for i in range(1,TERMS+1):
        values.append(values[-1]*(Q(i-1)-e)/i)
    return tuple(values)


@lru_cache(maxsize=None)
def transform_family(n,weight,direction):
    result=[]
    taylors=tuple(series(Q(k*(n-k),n-1)) for k in range(n+1))
    for row,h in zip(cayley_matrix(n),range(n+1)):
        weights=tuple(Q(a*comb(n,k),2**n)*direction[k] for k,a in enumerate(row))
        coefficients=tuple(sum(a*t[i] for a,t in zip(weights,taylors,strict=True)) for i in range(TERMS))
        tail=sum(abs(a*t[TERMS]) for a,t in zip(weights,taylors,strict=True))
        r=next((i for i,v in enumerate(coefficients) if v),TERMS)
        power=2*(weight+r)+h-n
        require(power>=0,("Cayley cancellation",n,weight,h,power))
        result.append((h,power,tuple(map(rational,coefficients[r:])),TERMS-r,rational(tail)))
    return tuple(result)


def model_eval(family,lam,q):
    result=[]
    for h,degree,polynomial,tail_degree,tail_coeff in family:
        v=ZERO
        for a in reversed(polynomial):
            v=v*q+a
        error=tail_coeff*q**tail_degree/(ONE-q)
        result.append((h,(v+error.symmetric_error())*lam**degree))
    return tuple(result)


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
    a=1+(1-u*u)/80
    require((lo-1)**2*a**3 <=(u*(1+u)/160)**2,"outer ratio lower")
    require(hi>=Q(81,8)/(10-(u**-2-1)/8),"outer ratio upper")
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
        b=(scalar(2),*(scalar(u,SCALE)/(ratio.sqrt() if n==10 else ONE) for u in RAYS[0]),ZERO)
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
                require((sign!=prior)==(n==10 if i%2 else n==9),("Hermite anchor sectors",n,i))
            prior=sign


def powers(n,x,y):
    re=[ONE]; im=[ZERO]
    for _ in range(n):
        re.append(re[-1]*x-im[-1]*y)
        im.append(re[-2]*y+im[-1]*x)
    return re,im


def accurate_projection(n,x,y,base,groups,coarse):
    def dot(terms,re,im):
        return projected_interval(terms,re), projected_interval(terms,im)

    re,im=powers(n,x,y)
    if any(a.lo>0 or a.hi<0 for a in dot(coarse,re,im)):
        return True

    cx=(x.lo+x.hi)//2
    cy=(y.lo+y.hi)//2
    dx=x-I(cx,cx)
    dy=y-I(cy,cy)
    re0,im0=powers(n,I(cx,cx),I(cy,cy))
    remainder_radius=dx.square()+dy.square()
    absolute_radius=(x.square()+y.square()).sqrt()
    radii=[ONE]
    for _ in range(n-2):
        radii.append(radii[-1]*absolute_radius)

    def taylor(terms):
        a,b=dot(terms,re0,im0)
        derivative=[(k-1,I(k*v.lo,k*v.hi)) for k,v in terms if k]
        c,d=dot(derivative,re0,im0)
        remainder=ZERO
        for k,v in terms:
            if k>=2:
                cap=max(abs(v.lo),abs(v.hi))*(k*(k-1)//2)
                remainder+=I(cap,cap)*radii[k-2]
        extra=(remainder*remainder_radius).symmetric_error()
        return a+c*dx-d*dy+extra,b+c*dy+d*dx+extra

    if any(a.lo>0 or a.hi<0 for a in taylor(coarse)):
        return True
    center=taylor(base)
    families=[(parameter,taylor(group)) for parameter,group in groups]

    def bounds(rect,a,b):
        re,im=rect
        return a*(re.lo if a>=0 else re.hi)+b*(im.lo if b>=0 else im.hi),a*(re.hi if a>=0 else re.lo)+b*(im.hi if b>=0 else im.lo)

    axes=[(1,0),(0,1)]
    axes.append(tuple(v.lo+v.hi for v in center))
    for par,(re1,im1) in families:
        if par.lo<par.hi:
            axes.append((im1.lo+im1.hi,-re1.lo-re1.hi))
    for r,i in zip(re0,im0,strict=True):
        axes.append((i.lo+i.hi,-r.lo-r.hi))
    for a,b in axes:
        if not (a or b):
            continue
        low,high=bounds(center,a,b)
        low*=GRID; high*=GRID
        for par,rect in families:
            u,v=bounds(rect,a,b)
            products=(u*par.lo,u*par.hi,v*par.lo,v*par.hi)
            low+=min(products)
            high+=max(products)
        if low>0 or high<0:
            return True
    return False


def check_ray(n,i,left,right):
    ratio=interval(*RATIO)
    lo=interval(Q(left,1000),Q(right,1000)).sqrt()
    pending=[(lo,interval(Q(0),DISK),ratio,0)]
    count=0
    while pending:
        lam,x,r,depth=pending.pop()
        root_ratio=r.sqrt() if n==10 else ONE
        scale=lam*root_ratio
        p=lam.square()
        q=scale.square()
        b=ray_for(p,i,left,right)
        if n==10:
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
        wr=12*(r.hi-r.lo) if n==10 else 0
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
    require(n in (0,*J) and -1<=i<8,"replay ray")
    require(low in RAYS and high in RAYS and low<high,"replay range")
    total=0
    for a,b in zip(RAYS,tuple(RAYS)[1:]):
        if not low<=a<b<=high:
            continue
        for order in J:
            if n and n!=order:
                continue
            for index in range(8):
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
