"""Exact terminal-ten stable Cayley exclusions on the upper marked-row complement."""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial, prod
from sys import argv

from verify_multivariate_terminal_ten_near_flat import accurate_projection
from verify_multivariate_terminal_seven_upper_bridge import cayley_matrix, height
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, bisect, interval, rational, require, scalar,
)

J=(9,10)
LOWER_DEFICIT=Q(3,25)
LAMBDA_TOP=Q(3464102,10**7)
QMAX=Q(61,500)
RATIO=(Q(5051,5000),Q(5072,5000))
DISK=Q(2)
SCALE=100000
C=1/(1-QMAX)
TERMS=12

RAYS={
    70: (118946, 106024, 84869, 76286, 55917, 47413, 28803, 21498),
    75: (119328, 107486, 84846, 76428, 55993, 47472, 28831, 21518),
    80: (119819, 107795, 85021, 76571, 56070, 47531, 28859, 20992),
    85: (120209, 108107, 85395, 76715, 56148, 47591, 28888, 21012),
    90: (120602, 108422, 85573, 76859, 56225, 47650, 28916, 21031),
    95: (120999, 108739, 85753, 77004, 56304, 47710, 28945, 21051),
    100: (121510, 109059, 85933, 77151, 56101, 47771, 28974, 21071),
    105: (121916, 109382, 86114, 77298, 56179, 47831, 29003, 21090),
    110: (122325, 109707, 86297, 77446, 56258, 47892, 29032, 21660),
    115: (122738, 110036, 86481, 77595, 56337, 47954, 29062, 21680),
    120: (123155, 110367, 86666, 77745, 56417, 48015, 29091, 21701),
}
EXPECTED_CELLS = {
    (9, 70): (483, 39, 79, 43, 60, 35, 26, 20),
    (9, 75): (324, 32, 79, 43, 50, 32, 24, 20),
    (9, 80): (251, 31, 79, 43, 46, 30, 24, 20),
    (9, 85): (210, 31, 75, 39, 45, 28, 22, 20),
    (9, 90): (193, 31, 67, 38, 45, 28, 21, 20),
    (9, 95): (175, 31, 67, 38, 41, 28, 21, 20),
    (9, 100): (148, 31, 66, 34, 41, 28, 21, 20),
    (9, 105): (121, 31, 66, 34, 39, 28, 21, 20),
    (9, 110): (108, 31, 62, 32, 39, 28, 21, 20),
    (9, 115): (108, 31, 61, 33, 39, 28, 21, 20),
    (10, 70): (60612, 111, 384, 105, 173, 36, 57, 25),
    (10, 75): (42160, 102, 294, 105, 119, 36, 53, 25),
    (10, 80): (26241, 101, 293, 106, 94, 36, 33, 25),
    (10, 85): (14894, 101, 291, 106, 88, 36, 29, 21),
    (10, 90): (11054, 101, 287, 106, 88, 36, 29, 19),
    (10, 95): (10514, 100, 288, 101, 88, 36, 29, 19),
    (10, 100): (9993, 89, 278, 98, 88, 36, 29, 19),
    (10, 105): (7288, 87, 269, 98, 88, 35, 26, 19),
    (10, 110): (5340, 88, 262, 98, 88, 35, 26, 19),
    (10, 115): (4229, 88, 262, 98, 88, 33, 26, 17),
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


def parameter_bound(d,a):
    return C*a**(d-1)/d


def reciprocal_cap(u):
    l=1
    while u*u>Q(9*l+1,9*l+9):
        l+=1
    return min(Q(1,55),(u**-2-1)/8,max(Q(1,9*l+1),(1-u)/(9*u-1)))


def moment_ratio(n,p_lo,p_hi):
    c=reciprocal_cap(1-p_hi)
    q=p_lo
    if n==10:
        c=9*c/(10-c)
        q*=RATIO[0]
    return min(C,(n-1)*c/(2+(n-1)*c)/q)


@lru_cache(maxsize=None)
def log_vectors(n):
    return {d:log_shape(n,d) for d in range(2,8)}


@lru_cache(maxsize=None)
def terms(n,a):
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
            bound=prod(parameter_bound(d,a) for d in indices)
            signed=any(d%2==0 and a%2 for d,a in occurrences.items())
            parameter=sb(bound) if signed else interval(Q(0),bound)
            items.append((weight,direction,parameter))
    return tuple(items)


@lru_cache(maxsize=None)
def raw_error(n,a):
    source=log_vectors(n)
    bounds=[]
    for k in range(n+1):
        r=min(k,n-k)
        omega=sum(abs((int(r<=h<n-r)-Q(n-2*r,n))/2) for h in range(n))
        log_tail=(omega*C*a**7/8+r*QMAX*C*a**8/9)/(1-(QMAX*a)**2)
        h={d:abs(vec[k])*parameter_bound(d,a)*QMAX**d for d,vec in source.items()}
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
def model(n,lo,hi,a):
    lam=I(lo,hi)
    q=lam.square()
    core=model_eval(transform_family(n,0,(Q(1),)*(n+1)),lam,q)
    groups=tuple((parameter,model_eval(transform_family(n,weight,direction),lam,q))
                 for weight,direction,parameter in terms(n,a))
    bounds=raw_error(n,a)
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
    require(hi>=Q(81,8)/(10-reciprocal_cap(u)),"outer ratio upper")
    require(LAMBDA_TOP**2>=LOWER_DEFICIT and hi*LAMBDA_TOP**2<=QMAX,"deficit bounds")
    require(DISK**2*QMAX<1,"Cayley disk left half-plane")
    for slot in RAYS.values():
        a=(2*SCALE,*slot,0)
        require(all(x>y for x,y in zip(a,a[1:])), "strict ray order")


def ray_for(p,i,left,right):
    a,b=RAYS[left][i],RAYS[right][i]
    return rational(Q(a,SCALE)) + rational(Q(b-a,SCALE))*(p*scalar(1000)-scalar(left))/scalar(right-left)


def anchor():
    p=rational(Q(min(RAYS),1000))
    ratio=interval(*RATIO)
    levels=RAYS[min(RAYS)]
    for n in J:
        pending=[(ratio,0)]
        while pending:
            r,depth=pending.pop()
            root_r=r.sqrt() if n==10 else ONE
            h=p*r if n==10 else p
            scale=h.sqrt()
            coeff=tuple(v for k,v in model_eval(transform_family(n,0,(Q(1),)*(n+1)),scale,h))
            rays=(scalar(2),*(scalar(z,SCALE)/root_r for z in levels),ZERO)
            polynomial=tuple(scalar((-1)**((h-n%2)//2))*coeff[h] for h in range(n%2,n+1,2))
            prior=0
            for i,y in enumerate(rays):
                x=y.square()
                xm=(x.lo+x.hi)//2
                center=I(xm,xm)
                f=ZERO
                for term in reversed(polynomial):
                    f=f*center+term
                slope=ZERO
                for k in range(len(polynomial)-1,0,-1):
                    slope=slope*x+scalar(k)*polynomial[k]
                value=f+slope*(x-center)
                sign=1 if value.lo>0 else -1 if value.hi<0 else 0
                if not sign or (i and (sign!=prior)!=(n==10 if i%2 else n==9)):
                    break
                prior=sign
            else:
                continue
            require(depth<16 and n==10,("reference anchor signs",n,r.shown(),i,value.shown()))
            pending.extend((half,depth+1) for half in bisect(r))


def check_ray(n,i,left,right):
    ratio=interval(*RATIO)
    lo=interval(Q(left,1000),Q(right,1000)).sqrt()
    aa=moment_ratio(n,Q(left,1000),Q(right,1000))
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
        base,groups,box,ok=model(n,scale.lo,scale.hi,aa)
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
    low,high=map(int,argv[3:]) if len(argv)>3 else (min(RAYS),max(RAYS))
    require(n in (0,*J) and -1<=i<8,"replay ray")
    require(low in RAYS and high in RAYS and low<high,"replay bounds")
    minimum_u=1-LOWER_DEFICIT
    require(Q(19,27)<minimum_u*minimum_u and (1-minimum_u)/(9*minimum_u-1)<Q(1,55),
            "unmarked unit and nonunit reciprocal cap")
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
