"""Exact cap-sensitive Cayley separation on the terminal-eleven middle band."""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial, prod
from sys import argv
from typing import NamedTuple

from verify_multivariate_terminal_seven_upper_bridge import cayley_matrix, height
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, bisect, interval, rational, require, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds
from verify_multivariate_terminal_ten_near_flat import (
    accurate_projection, log_vectors, model_eval, sb, transform_family,
)


ORDERS = (10, 11)
SCALE = 100000
REFERENCE_RATIO = Q(100, 99)


class Limits(NamedTuple):
    cap: Q
    ratio: tuple
    qmax: Q
    disk: Q
    clip: bool


LIMITS = Limits(Q(1,55), (Q(10084,10000),Q(10110,10000)), Q(243,1000), Q(13,5), True)
RAYS = {
    120: (124795, 111938, 91020, 82813, 62997, 55912, 38164, 30526, 14652),
    125: (125328, 112283, 91228, 82986, 63095, 55992, 38207, 31096, 14666),
    130: (125768, 112632, 91438, 83159, 63433, 56072, 38251, 31129, 14679),
    135: (126213, 112984, 91825, 83334, 63533, 56153, 38294, 31162, 14322),
    140: (126662, 113339, 92038, 83510, 63633, 56234, 38338, 31195, 14336),
    145: (127115, 113697, 92253, 83687, 63734, 56316, 38381, 31229, 14349),
    150: (127574, 114059, 92291, 83865, 63836, 56398, 38425, 31262, 14362),
    155: (128139, 114424, 92688, 84045, 63938, 56481, 38470, 31296, 14376),
    160: (128608, 114792, 92908, 84226, 64041, 56564, 38514, 31330, 14389),
    165: (129081, 115164, 93129, 84408, 64144, 56648, 38559, 31364, 14775),
    170: (129559, 115539, 93352, 84591, 64002, 56732, 38604, 31398, 14789),
    175: (130042, 115918, 93577, 84776, 64105, 56816, 38649, 31432, 14803),
    180: (130530, 116301, 93803, 84962, 64210, 56901, 38695, 31466, 14818),
    185: (131024, 116687, 94031, 85149, 64315, 56987, 38740, 31501, 14832),
    190: (131522, 117077, 94261, 85338, 64420, 57073, 38786, 31536, 14846),
    195: (132137, 117471, 94492, 85528, 64527, 57159, 38832, 31571, 14861),
    200: (132648, 117869, 94726, 85719, 64634, 57247, 38879, 31606, 14875),
    205: (133163, 118271, 94961, 85912, 64741, 57334, 38925, 31642, 14889),
    210: (133685, 118677, 95198, 86106, 64849, 57422, 38972, 31677, 14904),
    215: (134212, 119086, 95436, 86302, 64958, 57511, 39019, 31713, 14919),
    220: (134745, 119500, 95677, 86499, 65068, 57600, 39067, 31749, 14933),
    225: (135284, 119919, 95920, 86697, 65178, 57689, 39115, 31785, 14948),
    230: (135828, 120341, 96164, 86897, 65289, 57780, 39162, 31821, 14963),
    235: (136379, 120768, 96410, 87099, 65400, 57870, 39211, 31858, 14978),
    240: (136936, 121199, 96659, 87302, 65513, 57962, 39259, 31895, 14993),
}

EXPECTED_CELLS = {
    (10, 120): (472, 42, 104, 46, 79, 41, 45, 30, 16),
    (10, 125): (375, 42, 100, 46, 68, 41, 45, 27, 15),
    (10, 130): (332, 42, 96, 46, 57, 41, 44, 27, 17),
    (10, 135): (325, 42, 96, 46, 57, 37, 40, 28, 15),
    (10, 140): (299, 42, 95, 46, 57, 37, 36, 28, 15),
    (10, 145): (277, 42, 95, 46, 57, 37, 32, 28, 15),
    (10, 150): (237, 42, 95, 44, 57, 37, 31, 28, 16),
    (10, 155): (188, 42, 72, 42, 56, 37, 31, 28, 16),
    (10, 160): (182, 43, 71, 39, 56, 37, 30, 28, 16),
    (10, 165): (181, 43, 72, 39, 56, 37, 28, 28, 16),
    (10, 170): (171, 41, 72, 39, 55, 34, 27, 26, 16),
    (10, 175): (155, 41, 72, 38, 52, 34, 27, 25, 16),
    (10, 180): (156, 41, 70, 38, 52, 34, 27, 25, 16),
    (10, 185): (147, 41, 69, 38, 52, 33, 27, 25, 16),
    (10, 190): (134, 41, 60, 38, 51, 33, 27, 24, 15),
    (10, 195): (126, 42, 59, 38, 50, 33, 27, 22, 15),
    (10, 200): (126, 39, 59, 37, 50, 33, 27, 21, 15),
    (10, 205): (126, 39, 59, 37, 50, 33, 27, 20, 15),
    (10, 210): (118, 39, 55, 37, 50, 33, 27, 20, 14),
    (10, 215): (118, 37, 55, 38, 46, 32, 27, 20, 14),
    (10, 220): (113, 35, 55, 38, 46, 32, 27, 20, 14),
    (10, 225): (104, 32, 55, 38, 46, 32, 27, 20, 14),
    (10, 230): (104, 32, 55, 38, 46, 32, 27, 20, 14),
    (10, 235): (105, 32, 55, 38, 46, 32, 27, 19, 14),
    (11, 120): (44262, 105, 593, 133, 137, 45, 81, 34, 25),
    (11, 125): (38852, 97, 569, 133, 137, 44, 62, 31, 25),
    (11, 130): (30025, 93, 577, 132, 133, 44, 48, 29, 23),
    (11, 135): (3984, 88, 567, 133, 124, 44, 39, 28, 18),
    (11, 140): (3715, 73, 545, 134, 115, 44, 37, 28, 18),
    (11, 145): (3451, 74, 484, 134, 112, 44, 37, 26, 18),
    (11, 150): (3437, 66, 466, 125, 112, 44, 37, 25, 18),
    (11, 155): (3416, 65, 480, 122, 112, 42, 37, 25, 18),
    (11, 160): (3083, 79, 534, 124, 112, 42, 37, 25, 18),
    (11, 165): (2535, 78, 447, 123, 107, 42, 37, 25, 18),
    (11, 170): (1958, 72, 409, 107, 104, 42, 37, 25, 18),
    (11, 175): (1395, 67, 387, 107, 104, 41, 37, 25, 18),
    (11, 180): (979, 65, 386, 83, 105, 41, 36, 25, 18),
    (11, 185): (770, 57, 383, 83, 105, 41, 33, 25, 18),
    (11, 190): (921, 55, 369, 79, 105, 41, 33, 25, 18),
    (11, 195): (913, 53, 358, 79, 101, 41, 32, 25, 17),
    (11, 200): (735, 51, 353, 77, 99, 41, 32, 24, 17),
    (11, 205): (626, 50, 353, 70, 99, 41, 32, 24, 17),
    (11, 210): (590, 50, 353, 70, 99, 39, 32, 24, 17),
    (11, 215): (579, 50, 351, 64, 91, 39, 32, 24, 17),
    (11, 220): (569, 51, 351, 54, 89, 39, 32, 23, 17),
    (11, 225): (562, 52, 342, 52, 89, 39, 32, 23, 17),
    (11, 230): (541, 53, 332, 50, 88, 38, 32, 23, 17),
    (11, 235): (517, 51, 323, 49, 88, 36, 32, 22, 16),
}


def reciprocal_cap(u, cap):
    return min(cap, (u**-2-1)/9)


def strength_factor(qmax):
    return (1-qmax/2)/(1-qmax)


def moment_ratio(n, p_lo, p_hi, limits):
    cap = reciprocal_cap(1-p_hi, limits.cap)
    q_lo = p_lo
    if n == 11:
        cap = 10*cap/(11-cap)
        q_lo *= limits.ratio[0]
    tanh_cap = (n-1)*cap/(2+(n-1)*cap)
    return min(strength_factor(limits.qmax), tanh_cap/q_lo)


@lru_cache(maxsize=None)
def terms(n, a, c):
    source = log_vectors(n)
    items = []
    for degree in range(1, 4):
        for indices in combinations_with_replacement(source, degree):
            weight = sum(indices)
            if weight > 7:
                continue
            occurrences = Counter(indices)
            factor = prod(factorial(m) for m in occurrences.values())
            direction = tuple(prod(source[d][k] for d in indices)/Q(factor) for k in range(n+1))
            bound = prod(c*a**(d-1)/d for d in indices)
            signed = any(d%2 == 0 and m%2 for d,m in occurrences.items())
            parameter = sb(bound) if signed else interval(Q(0), bound)
            items.append((weight, direction, parameter))
    return tuple(items)


@lru_cache(maxsize=None)
def raw_error(n, a, qmax):
    c = strength_factor(qmax)
    bounds = []
    for k in range(n+1):
        r = min(k,n-k)
        omega = sum(abs((int(r<=h<n-r)-Q(n-2*r,n))/2) for h in range(n))
        log_tail = (omega*c*a**7/8+r*qmax*c*a**8/9)/(1-(qmax*a)**2)
        heights = {d:abs(vec[k])*c*a**(d-1)*qmax**d/d for d,vec in log_vectors(n).items()}
        retained_height = sum(heights.values())
        full_height = retained_height+log_tail*qmax**8
        require(full_height<1, ("small log", n, k))
        exp_tail = retained_height**4/(24*(1-retained_height))
        for degree in range(1,4):
            for indices in combinations_with_replacement(heights,degree):
                if sum(indices)>=8:
                    exp_tail += prod(heights[d] for d in indices)/prod(
                        factorial(m) for m in Counter(indices).values())
        bounds.append(log_tail/(1-full_height)+exp_tail/qmax**8)
    return tuple(bounds)


def model_derivative(family, lam, q):
    result = []
    for h, degree, polynomial, tail_degree, tail_coeff in family:
        value = derivative = ZERO
        for coefficient in reversed(polynomial):
            derivative = derivative*q+value
            value = value*q+coefficient
        tail = tail_coeff*q**tail_degree/(ONE-q)
        tail_slope = tail/(ONE-q)
        if tail_degree:
            tail_slope += tail_coeff*scalar(tail_degree)*q**(tail_degree-1)/(ONE-q)
        derivative = (derivative+tail_slope.symmetric_error())*scalar(2)*lam**(degree+1)
        if degree:
            derivative += (value+tail.symmetric_error())*scalar(degree)*lam**(degree-1)
        result.append((h, derivative))
    return tuple(result)


@lru_cache(maxsize=8192)
def model(n, lo, hi, a, limits, linearized):
    lam = I(lo, hi)
    q = lam.square()
    family = transform_family(n,0,(Q(1),)*(n+1))
    if linearized:
        mid = I((lo+hi)//2, (lo+hi)//2)
        core = model_eval(family, mid, mid.square())
        variations = ((lam-mid, model_derivative(family,lam,q)),)
    else:
        core = model_eval(family,lam,q)
        variations = ()
    groups = variations + tuple(
        (parameter, model_eval(transform_family(n,weight,direction),lam,q))
        for weight,direction,parameter in terms(n,a,strength_factor(limits.qmax)))
    bounds = raw_error(n,a,limits.qmax)
    base = []
    for (h,v), row in zip(core,cayley_matrix(n),strict=True):
        err = sum(Q(abs(entry)*comb(n,k),2**n)*bounds[k] for k,entry in enumerate(row))
        base.append((h,v+(rational(err)*lam**(16+h-n)).symmetric_error()))
    box = []
    for h,b in base:
        for variable,group in groups:
            b += variable*group[h][1]
        box.append((h,b))
    disk = limits.disk
    dominates = box[-1][1].lo*disk**n > sum(max(abs(v.lo),abs(v.hi))*disk**h for h,v in box[:-1])
    return tuple(base),groups,tuple(box),dominates


def verify_constants(limits, rays):
    knots = tuple(sorted(rays))
    require(all(b-a==5 for a,b in zip(knots,knots[1:])), "contiguous five-thousandth subbands")
    u = 1-Q(knots[-1],1000)
    lo,hi = limits.ratio
    cap = reciprocal_cap(u, limits.cap)
    require(lo <= REFERENCE_RATIO*(1-(REFERENCE_RATIO-1)*(1-u)/(2*u)), "outer ratio lower")
    require(hi >= 1+log_bounds(1/(1-cap/11))[1]/log_bounds(1+9*cap)[0], "outer ratio upper")
    max_scale = rational(Q(knots[-1],1000)).sqrt()*rational(hi).sqrt()
    require(Q(max_scale.square().hi,GRID)<limits.qmax, "deficit bound")
    if not limits.clip:
        require(limits.disk**2*limits.qmax<1,"Cayley disk left half-plane")
    for slot in rays.values():
        levels = (int(limits.disk*SCALE),*slot,0)
        require(len(slot)==9 and all(a>b for a,b in zip(levels,levels[1:])), "strict ray order")
        require(Q(max(slot),SCALE)**2*limits.qmax < 1, "left half-plane rays")


def ray_for(p,i,left,right,rays):
    a,b = rays[left][i],rays[right][i]
    return rational(Q(a,SCALE)) + rational(Q(b-a,SCALE))*(p*scalar(1000)-scalar(left))/scalar(right-left)


def verify_anchor(limits,rays):
    first = min(rays)
    p = rational(Q(first,1000))
    ratio = interval(*limits.ratio)
    for n in ORDERS:
        pending = [(ratio,0)]
        while pending:
            r,depth = pending.pop()
            root_r = r.sqrt() if n==11 else ONE
            h = p*r if n==11 else p
            scale = h.sqrt()
            coeff = tuple(v for _,v in model_eval(transform_family(n,0,(Q(1),)*(n+1)),scale,h))
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
                if not sign or (i and (sign!=prior)!=(n==11 if i%2 else n==10)):
                    break
                prior = sign
            else:
                continue
            require(depth<16 and n==11, ("reference anchor signs", n, r.shown(), i, value.shown()))
            pending.extend((half,depth+1) for half in bisect(r))


def check_ray(n,i,left,right,limits,rays,linearized=False):
    ratio = interval(*limits.ratio)
    lam = interval(Q(left,1000),Q(right,1000)).sqrt()
    a = moment_ratio(n,Q(left,1000),Q(right,1000),limits)
    pending = [(lam,interval(Q(0),limits.disk),ratio,0)]
    count = 0
    while pending:
        lam,x,r,depth = pending.pop()
        root_ratio = r.sqrt() if n==11 else ONE
        scale = lam*root_ratio
        p,q = lam.square(),scale.square()
        b = ray_for(p,i,left,right,rays)
        if n==11:
            b = b/root_ratio
        base,groups,box,ok = model(n,scale.lo,scale.hi,a,limits,linearized)
        if not ok:
            require(depth<32,("scaled disk",n,i,left,right,lam.shown()))
            pending.extend((half,x,r,depth+1) for half in bisect(lam))
            continue
        if limits.clip:
            max_x = (ONE/q.sqrt()).hi
            if x.lo>max_x:
                continue
            x = I(x.lo,min(x.hi,max_x))
        y = height(q,x,b)
        if accurate_projection(n,x,y,base,groups,box):
            count += 1
            continue
        require(depth<44,("unresolved",n,i,left,right,lam.shown(),x.shown(),r.shown()))
        wl,wx,wr = 10*(lam.hi-lam.lo),x.hi-x.lo,12*(r.hi-r.lo) if n==11 else 0
        if wr>max(wl,wx):
            pending.extend((lam,x,half,depth+1) for half in bisect(r))
        elif wl>wx:
            pending.extend((half,x,r,depth+1) for half in bisect(lam))
        else:
            pending.extend((lam,half,r,depth+1) for half in bisect(x))
    return count

def main():
    verify_constants(LIMITS,RAYS)
    verify_anchor(LIMITS,RAYS)
    n,i = map(int,argv[1:3]) if len(argv)>1 else (0,-1)
    lower,upper = map(int,argv[3:]) if len(argv)>3 else (120,240)
    require(n in (0,*ORDERS) and -1<=i<9,"replay ray")
    require(lower in RAYS and upper in RAYS and lower<upper,"replay bounds")
    total = 0
    for a,b in zip(tuple(RAYS),tuple(RAYS)[1:]):
        if not lower<=a<b<=upper:
            continue
        for order in ORDERS:
            if n and n!=order:
                continue
            for index in range(9):
                if i>=0 and i!=index:
                    continue
                linearized = order==11 and index==0 and a>=135
                count = check_ray(order,index,a,b,LIMITS,RAYS,linearized)
                require(count==EXPECTED_CELLS[order,a][index],
                        ("replay cell count",order,a,b,index,count))
                print(order,a,b,index,count,flush=True)
                total += count
    print(f"Verified {total} projection cells.",flush=True)


if __name__=="__main__":
    main()
