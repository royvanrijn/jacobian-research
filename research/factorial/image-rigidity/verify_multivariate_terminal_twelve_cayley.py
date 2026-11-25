"""Exact grouped Cayley separation for shared light rows at terminal twelve."""

from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial, prod
from sys import argv

from verify_multivariate_terminal_eleven_middle_cayley import (
    Limits, accurate_projection, model_derivative, strength_factor,
)
from verify_multivariate_terminal_seven_upper_bridge import cayley_matrix, height
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, bisect, interval, rational, require, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds
from verify_multivariate_terminal_ten_near_flat import log_shape, model_eval, sb, transform_family

ORDERS = (11, 12)
LOWER, UPPER = ORDERS
REFERENCE_RATIO = Q(121, 120)
SCALE = 100000

MODES = {
    "flat": (None, 0, 50, Q(13,5), False),
    "upper": (Q(1,70), 80, 120, Q(12,5), False),
    "seam": (Q(1,90), 50, 80, Q(12,5), False),
    "middle": (Q(1,55), 120, 240, Q(13,5), True),
}

def reciprocal_cap(u, cap):
    mass_cap = (u**-2 - 1)/(LOWER-1)
    return mass_cap if cap is None else min(cap, mass_cap)

def lattice_cap(u, cap):
    if not cap:
        return cap
    d1 = max(1,-(-(cap.denominator-cap.numerator)//(LOWER*cap.numerator)))
    return max(Q(1,LOWER*d1+1),min(cap,(1-u)/(LOWER*u-1)))

def limits_for(cap, left, right, disk, clip, refined=False):
    p_lo, p_hi = Q(left,1000), Q(right,1000)
    u_min = 1-p_hi
    cap = reciprocal_cap(u_min, cap)
    if refined:
        cap = lattice_cap(u_min, cap)
    lo = REFERENCE_RATIO*(1-(REFERENCE_RATIO-1)*p_hi/(2*u_min))
    if cap:
        limit = 1+log_bounds(1/(1-cap/UPPER))[1]/log_bounds(1+(LOWER-1)*cap)[0]
    else:
        limit = REFERENCE_RATIO
    require(1 < limit < 2 and REFERENCE_RATIO<=limit, ("ratio upper", limit))
    hi = limit*(1-(limit-1)*p_lo/2)
    qmax = Q(int(hi*p_hi*10000)+1,10000)
    requires = rational(p_hi).sqrt()*rational(hi).sqrt()
    require(Q(requires.square().hi,GRID)<qmax,"deficit bound")
    return Limits(cap,(lo,hi),qmax,disk,clip)

def moment_ratio(n,left,right,limits):
    cap = limits.cap
    q_lo = Q(left,1000)
    if n==UPPER:
        cap = LOWER*cap/(UPPER-cap)
        q_lo *= limits.ratio[0]
    tanh_cap = (n-1)*cap/(2+(n-1)*cap)
    c = strength_factor(limits.qmax)
    return min(c,tanh_cap/q_lo) if q_lo else c

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

def ray_for(p,i,left,right,rays):
    a,b = rays[left][i],rays[right][i]
    return rational(Q(a,SCALE)) + rational(Q(b-a,SCALE))*(p*scalar(1000)-scalar(left))/scalar(right-left)

def check_ray(n,i,left,right,limits,rays,linearized=False,window=None):
    ratio = interval(*limits.ratio)
    lo,hi = window or (left,right)
    lam = interval(Q(lo,1000),Q(hi,1000)).sqrt()
    a = moment_ratio(n,lo,hi,limits)
    pending = [(lam,interval(Q(0),limits.disk),ratio,0)]
    count = 0
    while pending:
        lam,x,r,depth = pending.pop()
        root_ratio = r.sqrt() if n==UPPER else ONE
        scale = lam*root_ratio
        p,q = lam.square(),scale.square()
        b = ray_for(p,i,left,right,rays)
        if n==UPPER:
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
        require(depth<44,("unresolved",n,i,left,right,lam.shown(),x.shown(),r.shown(),count))
        wl,wx,wr = 10*(lam.hi-lam.lo),x.hi-x.lo,12*(r.hi-r.lo) if n==UPPER else 0
        if wr>max(wl,wx):
            pending.extend((lam,x,half,depth+1) for half in bisect(r))
        elif wl>wx:
            pending.extend((half,x,r,depth+1) for half in bisect(lam))
        else:
            pending.extend((lam,half,r,depth+1) for half in bisect(x))
    return count

def validate_rays(limits,rays,left,right):
    if not limits.clip:
        require(limits.disk**2*limits.qmax<1,"Cayley disk left half-plane")
    for end in (left,right):
        slot = rays[end]
        levels = (int(limits.disk*SCALE),*slot,0)
        require(len(slot)==UPPER-2 and all(a>b for a,b in zip(levels,levels[1:])), "strict ray order")
        require(Q(max(slot),SCALE)**2*limits.qmax < 1, "left half-plane rays")


WEIGHT = 9

def directions(n):
    return {d:log_shape(n,d) for d in range(2,WEIGHT+1)}

@lru_cache(maxsize=None)
def terms(n, a, c):
    source = directions(n)
    items = []
    for degree in range(1, WEIGHT//2+1):
        for indices in combinations_with_replacement(source, degree):
            weight = sum(indices)
            if weight > WEIGHT:
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
        m = WEIGHT+1
        log_tail = (omega*c*a**(m-1)/m+r*qmax*c*a**m/(m+1))/(1-(qmax*a)**2)
        heights = {d:abs(vec[k])*c*a**(d-1)*qmax**d/d for d,vec in directions(n).items()}
        retained_height = sum(heights.values())
        full_height = retained_height+log_tail*qmax**m
        require(full_height<1, ("small log", n, k))
        d0=WEIGHT//2+1
        exp_tail = retained_height**d0/(factorial(d0)*(1-retained_height))
        for degree in range(1,d0):
            for indices in combinations_with_replacement(heights,degree):
                if sum(indices)>=m:
                    exp_tail += prod(heights[d] for d in indices)/prod(
                        factorial(j) for j in Counter(indices).values())
        bounds.append(log_tail/(1-full_height)+exp_tail/qmax**m)
    return tuple(bounds)

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
        base.append((h,v+(rational(err)*lam**(2*(WEIGHT+1)+h-n)).symmetric_error()))
    box = []
    for h,b in base:
        for variable,group in groups:
            b += variable*group[h][1]
        box.append((h,b))
    disk = limits.disk
    dominates = box[-1][1].lo*disk**n > sum(max(abs(v.lo),abs(v.hi))*disk**h for h,v in box[:-1])
    return tuple(base),groups,tuple(box),dominates


RAYS = {
    "flat": {
        0: (116829, 110185, 89396, 83510, 65890, 60344, 44106, 38371, 23170, 17731),
        5: (117715, 104820, 91063, 79289, 66480, 58880, 44367, 34923, 21888, 11823),
        10: (117965, 105105, 91149, 79440, 66747, 58963, 44417, 34959, 21907, 11463),
        15: (118217, 105392, 91120, 79591, 66850, 59046, 44468, 34995, 22188, 11472),
        20: (118529, 105681, 91321, 79744, 66955, 59129, 44519, 35031, 22207, 11482),
        25: (118844, 105973, 91407, 79897, 67059, 59214, 44570, 35067, 22227, 11492),
        30: (119162, 106267, 91493, 80690, 67165, 59298, 44621, 35543, 22247, 11502),
        35: (119481, 106564, 91580, 80848, 67271, 59383, 44886, 35580, 22267, 11884),
        40: (119865, 106863, 91784, 81007, 67378, 59468, 44938, 35618, 22550, 11894),
        45: (120190, 108018, 91872, 81167, 67485, 59554, 44991, 35655, 22571, 11904),
        50: (120518, 108327, 91959, 81327, 67592, 59641, 45044, 35693, 22591, 11914),
    },
    "seam": {
        50: (120518, 108327, 91959, 81327, 67592, 59641, 45044, 35693, 22591, 11914),
        55: (120974, 107774, 92286, 81489, 67701, 59727, 45096, 35731, 22612, 11924),
        60: (121371, 108083, 92495, 81651, 67809, 59814, 45149, 35769, 22633, 11935),
        65: (121771, 108394, 92826, 81815, 67918, 59902, 45203, 35807, 22654, 11945),
        70: (122241, 108708, 93038, 81980, 68028, 59990, 45256, 35846, 22410, 11955),
        75: (122649, 109024, 93252, 82145, 68138, 60078, 45310, 35884, 22431, 11966),
        80: (123128, 109343, 93590, 81647, 68249, 60167, 45364, 35923, 22451, 11601),
    },
    "upper": {
        80: (122931, 109345, 93224, 82313, 68251, 60168, 45365, 35924, 22717, 11977),
        85: (123346, 109667, 93563, 82481, 68362, 60257, 45419, 35962, 22738, 11987),
        90: (123833, 109991, 93780, 82650, 68475, 60347, 45473, 36002, 22759, 11998),
        95: (124256, 110319, 94000, 82820, 68587, 60437, 45528, 36041, 22781, 12008),
        100: (124684, 110649, 94346, 82991, 68701, 60528, 45583, 36080, 22802, 12019),
        105: (125185, 110983, 94569, 83163, 68815, 60619, 45639, 36120, 22824, 12030),
        110: (125622, 111319, 94794, 82655, 68930, 60711, 45694, 36160, 22578, 12040),
        115: (126062, 111658, 95020, 82826, 69046, 60804, 45750, 36200, 22600, 12051),
        120: (126507, 112001, 95377, 82999, 69162, 60897, 45806, 36240, 22622, 12062),
    },
    "middle": {
        120: (126368, 112003, 94994, 83689, 69163, 61444, 45807, 36692, 22890, 12062),
        125: (126816, 112348, 95352, 83866, 69280, 60991, 45864, 36281, 22912, 12073),
        130: (127269, 112697, 95583, 84044, 69398, 61085, 45920, 36322, 22934, 12084),
        135: (127800, 113048, 95815, 84224, 69516, 61179, 45978, 36363, 22956, 12095),
        140: (128262, 113403, 96049, 84404, 69635, 61274, 46035, 36404, 22978, 12106),
        145: (128729, 113761, 96417, 83883, 69755, 61370, 46093, 36445, 23001, 12117),
        150: (129201, 114123, 96655, 84063, 69876, 61466, 46151, 36487, 23023, 12128),
        155: (129678, 114488, 96896, 84245, 69997, 61563, 46209, 36528, 22776, 12140),
        160: (130160, 114856, 97138, 84427, 69938, 61660, 46267, 36570, 22798, 12151),
        165: (130723, 115228, 97381, 84611, 70060, 61758, 46326, 36613, 22821, 12162),
        170: (131216, 115604, 97627, 84796, 70183, 61856, 46385, 36655, 22843, 12174),
        175: (131713, 115982, 98011, 84983, 70306, 61955, 46445, 36697, 22866, 12185),
        180: (132216, 116365, 98262, 85171, 70431, 62054, 46278, 36740, 22889, 12196),
        185: (132724, 116751, 98514, 85360, 70556, 62155, 46338, 36783, 22912, 12208),
        190: (133238, 117141, 98768, 85551, 70682, 62255, 46398, 36826, 22935, 12220),
        195: (133757, 117535, 99025, 85743, 70808, 62357, 46458, 36870, 22959, 12231),
        200: (134364, 117933, 99283, 85936, 70936, 62459, 46519, 36914, 22982, 12243),
        205: (134896, 118335, 99543, 86131, 71064, 62561, 46580, 36958, 23005, 12255),
        210: (135433, 118741, 99806, 86327, 71193, 62664, 46641, 37002, 23029, 12267),
        215: (135976, 119150, 100213, 86525, 71323, 62768, 46703, 37046, 23053, 12278),
        220: (136525, 119564, 100481, 86724, 71453, 62873, 46765, 37091, 23077, 12290),
        225: (137080, 119983, 100751, 86925, 71585, 62978, 46827, 37136, 23100, 12302),
        230: (137641, 120405, 101023, 86367, 71717, 63083, 46890, 37181, 23125, 12314),
        235: (138209, 120832, 101297, 86567, 71851, 63190, 46953, 37226, 23149, 12327),
        240: (138783, 121263, 101574, 86768, 71985, 63297, 47016, 37272, 23173, 12339),
    },
}

EXPECTED_CELLS = {
    "flat": {
        (11, 0): (788, 491, 1232, 906, 1095, 781, 596, 359, 233, 151),
        (12, 0): (1778, 948, 2273, 1327, 1651, 1088, 772, 454, 224, 186),
        (11, 5): (268, 123, 325, 235, 305, 197, 176, 85, 71, 24),
        (12, 5): (1118, 277, 1035, 343, 569, 299, 232, 113, 54, 39),
        (11, 10): (342, 124, 364, 245, 315, 200, 184, 92, 75, 24),
        (12, 10): (1297, 278, 1003, 343, 573, 299, 235, 114, 69, 42),
        (11, 15): (438, 125, 390, 252, 346, 206, 186, 92, 77, 35),
        (12, 15): (1458, 285, 1002, 343, 591, 302, 247, 123, 79, 46),
        (11, 20): (477, 132, 405, 268, 379, 214, 186, 92, 77, 36),
        (12, 20): (1619, 291, 1005, 343, 600, 303, 261, 125, 95, 46),
        (11, 25): (398, 80, 250, 151, 209, 118, 100, 53, 47, 21),
        (12, 25): (1894, 158, 547, 206, 330, 166, 135, 68, 52, 27),
        (11, 30): (532, 86, 259, 149, 217, 118, 104, 53, 66, 23),
        (12, 30): (4098, 158, 559, 209, 330, 168, 139, 67, 54, 32),
        (11, 35): (659, 91, 274, 151, 229, 118, 111, 58, 66, 24),
        (12, 35): (6070, 164, 591, 213, 329, 171, 142, 68, 58, 30),
        (11, 40): (1012, 104, 326, 152, 237, 123, 111, 59, 59, 26),
        (12, 40): (8977, 174, 604, 216, 331, 174, 149, 71, 58, 34),
        (11, 45): (2107, 103, 369, 153, 244, 123, 111, 60, 71, 26),
        (12, 45): (12353, 176, 941, 227, 335, 174, 152, 76, 58, 39),
    },
    "seam": {
        (11, 50): (3638, 76, 413, 153, 238, 122, 112, 56, 71, 26),
        (12, 50): (13488, 168, 1672, 214, 336, 174, 152, 73, 59, 35),
        (11, 55): (1791, 76, 351, 149, 230, 118, 112, 56, 60, 26),
        (12, 55): (6307, 166, 981, 193, 334, 171, 152, 73, 59, 30),
        (11, 60): (1253, 76, 294, 141, 210, 118, 112, 56, 47, 26),
        (12, 60): (2563, 161, 915, 193, 330, 170, 143, 68, 58, 30),
        (11, 65): (744, 76, 255, 140, 206, 118, 111, 52, 56, 26),
        (12, 65): (14162, 156, 886, 185, 307, 162, 137, 66, 58, 29),
        (11, 70): (634, 75, 255, 140, 206, 109, 110, 50, 62, 22),
        (12, 70): (12016, 144, 797, 184, 306, 162, 134, 66, 46, 29),
        (11, 75): (525, 72, 233, 139, 198, 109, 101, 50, 43, 24),
        (12, 75): (9173, 137, 776, 184, 306, 162, 134, 66, 46, 29),
    },
    "upper": {
        (11, 80): (6725, 77, 398, 142, 229, 117, 111, 56, 69, 26),
        (12, 80): (5178, 165, 1368, 194, 335, 170, 151, 72, 58, 32),
        (11, 85): (2053, 76, 317, 141, 212, 117, 111, 56, 50, 26),
        (12, 85): (3077, 149, 943, 194, 332, 171, 144, 69, 58, 30),
        (11, 90): (1309, 76, 291, 141, 206, 117, 111, 52, 47, 26),
        (12, 90): (2077, 141, 878, 186, 307, 163, 139, 67, 58, 30),
        (11, 95): (993, 77, 259, 141, 205, 109, 111, 49, 47, 26),
        (12, 95): (1424, 141, 851, 186, 307, 162, 133, 67, 58, 29),
        (11, 100): (726, 77, 240, 140, 197, 109, 110, 49, 47, 25),
        (12, 100): (14421, 141, 818, 184, 307, 163, 133, 67, 58, 29),
        (11, 105): (625, 73, 232, 129, 189, 105, 101, 49, 50, 22),
        (12, 105): (13393, 141, 773, 184, 307, 156, 133, 67, 47, 29),
        (11, 110): (329, 45, 151, 79, 113, 63, 62, 32, 33, 16),
        (12, 110): (5867, 79, 403, 116, 173, 90, 73, 43, 30, 19),
        (11, 115): (309, 45, 150, 79, 113, 63, 57, 32, 28, 16),
        (12, 115): (3779, 79, 388, 109, 172, 89, 73, 43, 30, 18),
    },
    "middle": {
        (11, 120): (3559, 50, 261, 84, 132, 70, 66, 38, 43, 18),
        (12, 120): (1397, 95, 738, 139, 224, 95, 82, 46, 35, 26),
        (11, 125): (1921, 50, 223, 84, 131, 67, 62, 38, 42, 18),
        (12, 125): (1381, 93, 511, 137, 193, 95, 81, 42, 35, 25),
        (11, 130): (1053, 50, 207, 84, 127, 67, 62, 37, 43, 18),
        (12, 130): (1631, 89, 489, 133, 182, 95, 82, 41, 35, 23),
        (11, 135): (759, 50, 207, 84, 127, 67, 62, 36, 36, 17),
        (12, 135): (3281, 85, 473, 126, 174, 94, 77, 41, 35, 22),
        (11, 140): (623, 48, 191, 86, 121, 67, 62, 36, 34, 17),
        (12, 140): (3424, 85, 469, 125, 175, 91, 77, 41, 35, 20),
        (11, 145): (554, 49, 171, 80, 112, 67, 62, 34, 30, 17),
        (12, 145): (8142, 83, 467, 124, 173, 91, 76, 41, 35, 20),
        (11, 150): (450, 49, 163, 79, 112, 62, 62, 33, 38, 16),
        (12, 150): (7127, 83, 425, 117, 173, 87, 74, 41, 35, 19),
        (11, 155): (397, 50, 154, 80, 112, 62, 58, 33, 41, 16),
        (12, 155): (10495, 78, 730, 116, 253, 87, 74, 41, 33, 18),
        (11, 160): (356, 48, 152, 79, 113, 62, 58, 33, 36, 16),
        (12, 160): (9815, 74, 628, 109, 246, 87, 74, 41, 29, 18),
        (11, 165): (317, 48, 151, 79, 111, 62, 57, 33, 30, 16),
        (12, 165): (9638, 71, 588, 102, 246, 87, 73, 41, 29, 18),
        (11, 170): (294, 46, 144, 71, 111, 60, 57, 33, 28, 16),
        (12, 170): (5278, 70, 574, 102, 242, 87, 73, 41, 29, 18),
        (11, 175): (271, 45, 133, 67, 111, 52, 55, 33, 28, 16),
        (12, 175): (2953, 62, 595, 102, 238, 81, 71, 41, 29, 18),
        (11, 180): (241, 45, 133, 59, 95, 46, 53, 33, 27, 16),
        (12, 180): (2416, 62, 499, 102, 223, 81, 68, 41, 29, 18),
        (11, 185): (218, 44, 131, 57, 91, 46, 53, 33, 27, 16),
        (12, 185): (2268, 59, 284, 102, 223, 81, 68, 42, 26, 18),
        (11, 190): (189, 44, 128, 57, 82, 44, 53, 32, 27, 16),
        (12, 190): (2202, 59, 240, 102, 215, 81, 66, 42, 23, 18),
        (11, 195): (175, 45, 123, 57, 81, 42, 52, 30, 27, 15),
        (12, 195): (2205, 60, 236, 102, 215, 77, 66, 40, 22, 18),
        (11, 200): (157, 44, 119, 55, 81, 40, 52, 27, 27, 15),
        (12, 200): (2221, 60, 235, 102, 215, 75, 62, 40, 22, 18),
        (11, 205): (153, 43, 110, 56, 81, 40, 52, 26, 27, 15),
        (12, 205): (2141, 60, 235, 102, 214, 67, 62, 39, 23, 18),
        (11, 210): (146, 43, 110, 54, 81, 40, 52, 26, 27, 15),
        (12, 210): (2047, 60, 236, 101, 198, 59, 58, 39, 23, 18),
        (11, 215): (146, 42, 105, 52, 77, 40, 52, 26, 27, 15),
        (12, 215): (1949, 60, 234, 103, 181, 59, 56, 39, 23, 18),
        (11, 220): (146, 41, 105, 48, 69, 40, 52, 24, 27, 15),
        (12, 220): (1829, 60, 233, 103, 179, 57, 56, 38, 23, 18),
        (11, 225): (135, 42, 103, 49, 69, 40, 52, 23, 27, 14),
        (12, 225): (1720, 60, 225, 103, 159, 57, 46, 36, 23, 18),
        (11, 230): (136, 42, 95, 48, 65, 40, 48, 23, 27, 14),
        (12, 230): (1556, 61, 217, 101, 140, 57, 46, 36, 23, 18),
        (11, 235): (137, 42, 92, 49, 66, 40, 48, 23, 27, 14),
        (12, 235): (1267, 57, 209, 101, 124, 53, 44, 36, 23, 17),
    },
}


def run(mode, n, ray, lo, hi):
    cap, first, last, disk, clip = MODES[mode]
    rays = RAYS[mode]
    nodes = sorted(rays)
    require(nodes == list(range(first,last+1,5)) and lo in nodes and hi in nodes and lo<hi, "range")
    require(n in (0,*ORDERS) and -1<=ray<UPPER-2, "ray selection")
    anchor_for(limits_for(cap,first,first,disk,clip), rays, first)
    total = 0
    for left,right in zip(nodes,nodes[1:]):
        if not lo<=left<right<=hi:
            continue
        limits = limits_for(cap,left,right,disk,clip)
        validate_rays(limits,rays,left,right)
        for order in ORDERS:
            if n and n!=order:
                continue
            for index in range(UPPER-2):
                if ray>=0 and ray!=index:
                    continue
                linearized = index in (0,1) or (order==UPPER and index==2)
                refined = order==UPPER and index==0 and (
                    (mode=="flat" and left==45) or (mode=="seam" and left<65)
                    or (mode=="upper" and left<100) or (mode=="middle" and left<145))
                count = 0
                if refined:
                    local_cap = Q(1,60) if mode=="middle" else cap
                    for sub in range(left,right):
                        narrow = limits_for(local_cap,sub,sub+1,disk,clip,refined=True)
                        validate_rays(narrow,rays,left,right)
                        count += check_ray(order,index,left,right,narrow,rays,linearized,(sub,sub+1))
                else:
                    count = check_ray(order,index,left,right,limits,rays,linearized)
                require(count==EXPECTED_CELLS[mode][order,left][index],
                        ("replay cell count",mode,order,left,index,count))
                print(mode,order,left,right,index,count,flush=True)
                total += count
    return total


def main():
    if len(argv)==1:
        total = sum(run(mode,0,-1,MODES[mode][1],MODES[mode][2]) for mode in MODES)
    else:
        mode = argv[1]
        params = tuple(map(int,argv[2:]))
        require(mode in MODES and len(params) in (0,4), "select band or band/order/ray/left/right")
        total = run(mode,*(params or (0,-1,MODES[mode][1],MODES[mode][2])))
    print(f"Verified {total} projection cells.",flush=True)


if __name__=="__main__":
    main()
