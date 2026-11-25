#!/usr/bin/env python3
"""Exact terminal-eleven light-row separation on bounded outer-mean bands."""

from ctypes import CDLL
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_nine_high_curvature import anchor_value
from verify_multivariate_terminal_ten_shared_light import normalized_native, shape
from verify_multivariate_terminal_nine_shared_matched_upper import ANNULUS,native_replay,replay_python
from verify_multivariate_terminal_seven_private_rows import DISK_RADIUS,verify_disk,root_bounds
from verify_multivariate_terminal_six_joint_ray_bands import (
    INNER as SMALL_INNER, OUTER as SMALL_OUTER, radial_replay,
)


from fractions import Fraction as Q
from functools import lru_cache
from math import comb,prod

from verify_multivariate_terminal_nine_upper_marked_rows import precise_power
from verify_multivariate_terminal_nine_shared_matched_light import exp_small
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID,I,ONE,ZERO,interval,log,rational,scalar,require,
)
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds,power_bounds,rounded_interval
from verify_multivariate_terminal_six_private_angles import multiply_interval
from verify_multivariate_terminal_seven_private_linear import skew_weights

LOWER, UPPER = ORDERS = (10, 11)
REFERENCE_RATIO = Q(100, 99)
FLOOR=Q(1,40)

@lru_cache(maxsize=None)
def ratio_bounds(cap):
    lo,hi=log_bounds(1+(LOWER-1)*cap)
    lower=log_bounds(1+Q(LOWER**2, UPPER)*cap)[0]/hi
    upper=1+log_bounds(1/(1-cap/UPPER))[1]/lo
    require(lower<=REFERENCE_RATIO<=upper,("ratio",cap))
    return interval(lower,upper)

def nominal_exponent(order,k):
    return Q(k*(order-k),order-1)*(1 if order==LOWER else REFERENCE_RATIO)

@lru_cache(maxsize=None)
def independent_exponents(order,cap):
    d0,d1=log_bounds(1+(LOWER-1)*cap)
    if order==LOWER:
        x=cap
        ratio=(Q(1),Q(1))
    else:
        x=LOWER*cap/(UPPER-cap)
        ratio=(log_bounds(1+Q(LOWER**2, UPPER)*cap)[0]/d1,1+log_bounds(1/(1-cap/UPPER))[1]/d0)
    m=order-1
    dlo,dhi=log_bounds(1+m*x)
    assert (m*x/(2+m*x))**2<=Q(1,2)
    out=[]
    for k in range(1,order):
        r=min(k,order-k)
        lo=sum(log_bounds((1+(m-h)*x)/(1+h*x))[0] for h in range(r))/dhi
        skew=2*sum(multiply_interval((w,w),log_bounds(1+h*x))[1]
                   for h,w in enumerate(skew_weights(order,r)) if w and h)/dlo
        hi=Q(k*(order-k),m)
        box=rounded_interval(ratio[0]*(lo-skew),ratio[1]*(hi+skew))
        require(0<box[0]<=nominal_exponent(order,k)<=box[1],("independent exponents",order,k,cap,box))
        out.append(box)
    return tuple(out)

def independent_coefficients(order,cap,means):
    out=tuple((comb(order,k)*power_bounds(means[0],b)[0],
               comb(order,k)*power_bounds(means[1],a)[1])
              for k,(a,b) in enumerate(independent_exponents(order,cap),1))
    return ((Q(1),Q(1)),*out,(Q(1),Q(1)))


def exact_cap(declared, outer_lower):
    """Enumerate rows above FLOOR that can occur at this outer mean."""
    threshold = 1+(LOWER-1)*FLOOR/(1+FLOOR)
    cap = min(FLOOR, declared)
    gap = 1
    while threshold**(-gap) >= outer_lower**2:
        base = 1
        while LOWER*base+1 < gap/FLOOR:
            x = Q(gap, LOWER*base+1)
            if cap < x <= declared:
                row_square = prod(Q(LOWER*base+h, LOWER*base+(LOWER-1)*gap+h) for h in range(1, gap+1))
                if row_square >= outer_lower**2:
                    cap = x
            base += 1
        gap += 1
    return cap


@lru_cache(maxsize=None)
def constants(order, lower_cap):
    cap = lower_cap if order == LOWER else LOWER*lower_cap/(UPPER-lower_cap)
    m = order-1
    T = rational(m*cap/(2+m*cap))
    lo, hi = log_bounds(1+m*cap)
    s = interval(lo/2, hi/2)
    require(T.square().hi < GRID//2, "moment quotient monotonicity range")
    return T, s, T.square()/(scalar(2)*s)


@lru_cache(maxsize=None)
def leading_skew_direction(order, lower_cap, skew_terms, mean_terms):
    _, skew, _, _, _ = shape(order, skew_terms, mean_terms)
    T, _, _ = constants(order, lower_cap)
    weights = tuple(T**(2*l-2)*scalar(1, 2*l) for l in range(2, skew_terms+1))
    direction = tuple(first+sum((w*v for w,v in zip(weights,rest,strict=True)),ZERO)
                      for first,*rest in zip(*skew,strict=True))
    return direction, weights


@lru_cache(maxsize=8192)
def coefficient_model(order, lower_cap, means, signed, extended=False, *, precise=False, monotone=False):
    skew_terms,mean_terms=(4,3) if extended else (3,2)
    midpoint = sum(means)/2
    skew_midpoint = sum(signed)/2
    ref_ratio = Q(1) if order == LOWER else REFERENCE_RATIO
    T, s, skew_bound = constants(order, lower_cap)
    Pmax = constants(LOWER, lower_cap)[2]
    lower_strength = -log(interval(*means))
    leading = interval(*signed)*Pmax*lower_strength
    nominal = rational(skew_midpoint)*Pmax*(-log(rational(midpoint)))
    h = interval(means[0]/midpoint-1, means[1]/midpoint-1)
    magnitude = I(0, max(-h.lo,h.hi))
    log_error = rational(ref_ratio)*magnitude.square()/(scalar(2)*(ONE-magnitude))
    outer_mean = I(-log_error.hi,0)
    strength = lower_strength
    if order == UPPER:
        rho = ratio_bounds(lower_cap)
        if monotone:
            require(lower_cap <= Q(1, 8), "matched outer ratio scope")
            rho = I(rational(REFERENCE_RATIO).lo, rho.hi)
        outer_mean += (rational(ref_ratio)-rho)*lower_strength
        strength = rho*lower_strength
        eta = lower_cap/(2*UPPER+(LOWER**2-2)*lower_cap)
        transfer = rational(REFERENCE_RATIO**2*(1+eta**2))
        transfer_width = rational(2*REFERENCE_RATIO**2*eta)
        leading = transfer*leading+(transfer_width*Pmax*lower_strength).symmetric_error()
        nominal = transfer*nominal
    gaussian, skew, mean, skew_sizes, distances = shape(order,skew_terms,mean_terms)
    direction, weights = leading_skew_direction(order,lower_cap,skew_terms,mean_terms)
    reference = []
    for k, v in enumerate(direction):
        if k == 0 or k == order:
            reference.append(ONE)
            continue
        exponent = ref_ratio*Q(k*(order-k), order-1)
        power = precise_power(midpoint, exponent) if precise else interval(*power_bounds(midpoint, exponent))
        value = scalar(comb(order, k))*power
        reference.append(value*exp_small(nominal*v))
    pieces = [(h, tuple(rational(ref_ratio)*v for v in gaussian)),
              (outer_mean, gaussian), (leading-nominal, direction)]
    for vector, weight in zip(skew[1:], weights, strict=True):
        bound = weight*skew_bound*strength
        pieces.append((bound.symmetric_error(), vector))
    lower = I(leading.lo, leading.lo) if leading.lo > 0 else ZERO
    upper_strength = I(strength.hi, strength.hi)
    for j, vector in enumerate(mean, 1):
        lo = (scalar(2)*lower/upper_strength)**(2*j)*upper_strength*scalar(1, 2*j+1)
        hi = T**(2*j+1)*strength/(scalar(2*j+1)*s)
        if lo.lo > hi.hi:
            # This rectangle is empty by the lower and upper moment bounds.
            return None
        pieces.append((I(lo.lo, hi.hi), vector))
    skew_tail = T**(2*skew_terms+1)*strength/(scalar(2*skew_terms+2)*(ONE-T.square()))
    defect_tail = T**(2*mean_terms+2)*strength/(scalar(2*mean_terms+3)*(ONE-T.square()))
    errors = []
    for k, c in enumerate(reference):
        tail = skew_sizes[k]*skew_tail+scalar(distances[k])*defect_tail
        size = tail
        for parameter, vector in pieces:
            term = parameter*vector[k]
            size += I(0, max(-term.lo, term.hi))
        if size.hi >= GRID:
            # Subdivide before applying the exponential bound.
            return ()
        if extended:
            error = scalar(1,24)*size**4/(ONE-size)+tail*(ONE+size+scalar(1,2)*size.square())
        else:
            error = scalar(1,6)*size**3/(ONE-size)+tail*(ONE+size)
        errors.append(c*(ONE+error.symmetric_error()))
    products = []
    for i, (first, a) in enumerate(pieces):
        for j, (second, b) in enumerate(pieces[:i+1]):
            parameter = scalar(1, 2)*first.square() if i == j else first*second
            products.append((parameter, tuple(v*w for v, w in zip(a, b, strict=True))))
    cubic=[]
    if extended:
        for i,(p,a) in enumerate(pieces):
            for j,(q,b) in enumerate(pieces[:i+1]):
                for k,(r,c) in enumerate(pieces[:j+1]):
                    if i==k:
                        parameter=scalar(1,6)*p.square()*p
                    elif i==j:
                        parameter=scalar(1,2)*p.square()*r
                    elif j==k:
                        parameter=scalar(1,2)*p*q.square()
                    else:
                        parameter=p*q*r
                    cubic.append((parameter,tuple(v*w*x for v,w,x in zip(a,b,c,strict=True))))
    return ((ONE, tuple(errors)),) + tuple(
        (parameter, tuple(c*v for c, v in zip(reference, vector, strict=True)))
        for parameter, vector in (*pieces,*products,*cubic)
    )



SCALE = 10**6


def verify_anchor_order(order, core, centers, rays):
    centers=tuple((Q(x,100000),Q(y,100000)) for x,y in centers)
    upper=tuple((x,y) for x,y in centers if y>0)
    require(len(upper)==order//2 and len(centers)==order//2+(order%2),("disk count",order))
    require(all(y>DISK_RADIUS for _,y in upper),"upper disks")
    disks=upper+tuple((x,-y) for x,y in upper)
    if order%2:
        x,y=centers[-1]
        require(y==0 and x< -DISK_RADIUS,"negative real anchor")
        disks+=((x,y),)
    for i,(x,y) in enumerate(disks):
        verify_disk(core,(x,y))
        for a,b in disks[:i]:
            require((x-a)**2+(y-b)**2>(2*DISK_RADIUS)**2,"disjoint anchor disks")
    sectors=range(int(order==LOWER), UPPER-1, 2)
    for (x,y),sector in zip(upper,sectors,strict=True):
        for j,cosine in enumerate(rays):
            sine=root_bounds(1-cosine*cosine,2)
            lo,hi=multiply_interval(sine,(-x,-x))
            if j<sector:
                require(cosine*y+lo>DISK_RADIUS,("anchor ray side",order,sector,j))
            else:
                require(cosine*y+hi< -DISK_RADIUS,("anchor ray side",order,sector,j))


def gaussian_anchor(order, u, centers, rays):
    coeffs=[]
    for k in range(order+1):
        power=precise_power(u,nominal_exponent(order,k))
        coeffs.append((Q(power.lo*comb(order,k),GRID),Q(power.hi*comb(order,k),GRID)))
    verify_anchor_order(order,coeffs,centers,rays)


def verify_setup(rays, anchors, declared, small, *, knot_scale=100):
    lower, upper = Q(min(rays), knot_scale), Q(max(rays), knot_scale)
    for values in rays.values():
        bounds=(SCALE,*values,-SCALE)
        require(len(values)==UPPER-2 and all(a>b for a,b in zip(bounds,bounds[1:])),"ordered rays")
    annulus=(SMALL_INNER,SMALL_OUTER) if small else ANNULUS
    for order in ORDERS:
        box = independent_coefficients(order, declared, (Q(0), upper))
        inner, outer = annulus
        dominance = (sum(cap*inner**k for k, (_,cap) in enumerate(box) if k) < 1
                     and sum(cap*outer**k for k, (_,cap) in enumerate(box[:-1])) < outer**order)
        if not dominance:
            cap = exact_cap(declared, lower)
            bounds = independent_coefficients(order, cap, (lower, upper))
            require(all(a > inner*d and b < outer*c
                        for (a,b), (c,d) in zip(bounds, bounds[1:])),
                    ("Enestrom annulus", order, lower, upper, cap))
        cosine=tuple(Q(x,SCALE) for x in rays[min(rays)])
        if small:
            for i,x in enumerate((Q(1),*cosine,Q(-1))):
                expected=(i+(order==UPPER))//2
                require((-1)**expected*anchor_value(order,x)>0,("binomial anchor",order,i))
        else:
            gaussian_anchor(order,lower,anchors[order],cosine)


def ray_cosines(rays, left, right, index, means, knot_scale):
    a, b = rays[left][index], rays[right][index]
    return tuple(sorted(Q(a, SCALE) + Q(b-a, SCALE)*(knot_scale*v-left)/(right-left)
                        for v in means))


def verify_small_band(order,cap,rays,left,right, *, knot_scale=100):
    pending=[((Q(left,knot_scale),Q(right,knot_scale)),tuple(range(UPPER-2)),0)]
    total=0
    while pending:
        means,indices,depth=pending.pop()
        box=independent_coefficients(order,cap,means)
        bad=[]
        for index in indices:
            used=radial_replay(order,box,ray_cosines(rays,left,right,index,means,knot_scale))
            if used is None:
                bad.append(index)
            else:
                total+=used
        if bad:
            require(depth<16,("small band",order,cap,means,bad))
            mid=sum(means)/2
            pending.extend((((means[0],mid),bad,depth+1),((mid,means[1]),bad,depth+1)))
    return total


def verify_fused_band(order,cap,rays,left,right,replay,extended, *, precise=False, monotone=False, knot_scale=100):
    pending=[((Q(left,knot_scale),Q(right,knot_scale)),(Q(0),Q(1)),tuple(range(UPPER-2)),0,0)]
    total=0
    while pending:
        means,p,active,mean_depth,skew_depth=pending.pop()
        model=coefficient_model(order,cap,means,p,extended,precise=precise,monotone=monotone)
        if model is None:
            continue
        bad=[]
        for index in active:
            a,b=rays[left][index],rays[right][index]
            slope=Q(b-a,SCALE)*knot_scale*sum(means)/(2*(right-left))
            count=replay(order,ray_cosines(rays,left,right,index,means,knot_scale),model,slope) if model else None
            if count is None:
                bad.append(index)
            else:
                total+=count
        if not bad:
            continue
        require(max(mean_depth,skew_depth)<17,
                ("unresolved ray",order,cap,left,right,means,p,bad))
        if mean_depth<2*skew_depth+2:
            mid=sum(means)/2
            pending.extend((((means[0],mid),p,bad,mean_depth+1,skew_depth),
                            ((mid,means[1]),p,bad,mean_depth+1,skew_depth)))
        else:
            mid=sum(p)/2
            pending.extend(((means,(p[0],mid),bad,mean_depth,skew_depth+1),
                            (means,(mid,p[1]),bad,mean_depth,skew_depth+1)))
    return total




def scaled_native(library, *, compress=False):
    """Use a common positive scale for the base and its derivative group."""
    replay = normalized_native(library)

    def project(order, cosines, model, slope):
        maximum = max(max(abs(c.lo), abs(c.hi)) for _, vector in model[:2] for c in vector)
        shift = max(0, (maximum - 1).bit_length() - 57)
        if shift:
            model = tuple((parameter, tuple(I(c.lo >> shift, -((-c.hi) >> shift))
                                            for c in vector))
                          for parameter, vector in model)
        if compress and len(model) > 2 and model[0][0] == ONE:
            base = list(model[0][1])
            retained = []
            for parameter, vector in model[2:]:
                radius = max(abs(parameter.lo), abs(parameter.hi))
                size = max(max(abs(c.lo), abs(c.hi)) for c in vector)
                if radius*size <= GRID**2 // 2**15:
                    for k, coefficient in enumerate(vector):
                        base[k] += parameter*coefficient
                else:
                    retained.append((parameter, vector))
            coarser = ((ONE, tuple(base)), model[1], *retained)
            count = replay(order, cosines, coarser, slope)
            if count is not None:
                return count
        return replay(order, cosines, model, slope)

    return project


BANDS = (
    (Q('1/8'), 'small', {}, {
        0: (955275, 802959, 621323, 365050, 71157, -207708, -501600, -714519, -896155),
        1: (953140, 795685, 608028, 348397, 51895, -224889, -514864, -722902, -899448),
        2: (950193, 787118, 593182, 331004, 32977, -241128, -526897, -730323, -902324),
        3: (946413, 777577, 577421, 313366, 14557, -256586, -538080, -737134, -904940),
        4: (941798, 767229, 561049, 295666, -3361, -271380, -548601, -743487, -907365),
    }, {
        0: (1, Q('1/8'), (141, 213)),
        1: (2, Q('1/8'), (262, 470)),
        2: (3, Q('1/8'), (389, 852)),
        3: (4, Q('1/8'), (675, 1667)),
    }),
    (Q('1/8'), 'quad', {10: ((93605, 35187), (52369, 85191), (-7436, 99723), (-62881, 77756), (-95657, 29152)), 11: ((94755, 31962), (59841, 80119), (6764, 99771), (-46839, 88352), (-85817, 51337), (-100000, 0))}, {
        4: (941798, 767229, 561049, 295666, -3361, -271380, -548601, -743487, -907365),
        6: (930141, 744558, 527165, 260425, -37798, -299314, -568098, -755152, -911787),
        8: (915454, 719763, 492426, 225651, -70571, -325438, -585989, -765761, -915778),
        10: (898014, 693270, 457282, 191463, -101902, -350086, -602627, -775561, -919444),
        12: (878088, 665375, 421979, 157889, -131976, -373500, -618247, -784713, -922851),
        14: (855917, 636294, 386658, 124918, -160944, -395856, -633015, -793327, -926045),
        16: (831712, 606192, 351402, 92524, -188931, -417294, -647054, -801485, -929059),
    }, {
        4: (6, Q('5/41'), (54, 85)),
        6: (8, Q('5/41'), (55, 72)),
        8: (10, Q('5/41'), (54, 63)),
        10: (12, Q('5/41'), (58, 78)),
        12: (14, Q('5/41'), (59, 276)),
        14: (16, Q('5/41'), (59, 315)),
    }),
    (Q('1/8'), 'cubic-precise', {10: ((82025, 57201), (31067, 95052), (-25224, 96766), (-71176, 70242), (-96691, 25512)), 11: ((84318, 53764), (39214, 91991), (-12562, 99208), (-58235, 81294), (-89121, 45359), (-100000, 0))}, {
        16: (831712, 606192, 351402, 92524, -188931, -417294, -647054, -801485, -929059),
        18: (805657, 575197, 316262, 60675, -216040, -437924, -660460, -809249, -931919),
        20: (777912, 543413, 281267, 29340, -242358, -457833, -673309, -816668, -934643),
    }, {
        16: (18, Q('5/41'), (73, 252)),
        18: (20, Q('5/51'), (64, 81)),
    }),
    (Q('1/8'), 'quad', {10: ((76504, 64398), (24075, 97059), (-30265, 95310), (-73360, 67958), (-96955, 24491)), 11: ((79078, 61210), (32178, 94681), (-18207, 98329), (-61302, 79007), (-89974, 43643), (-100000, 0))}, {
        20: (777912, 543413, 281267, 29340, -242358, -457833, -673309, -816668, -934643),
        22: (748615, 510923, 246434, -1514, -267957, -477097, -685662, -823782, -937248),
        24: (715092, 477798, 211775, -31917, -292900, -495775, -697571, -830623, -939747),
        26: (682952, 444096, 177292, -61895, -317241, -513921, -709078, -837217, -942151),
        28: (649596, 409868, 142988, -91473, -341026, -531579, -720218, -843587, -944467),
        30: (612076, 375156, 108862, -120676, -364297, -548788, -731024, -849754, -946705),
        32: (576474, 339999, 74912, -149524, -387090, -565582, -741521, -855733, -948871),
        34: (536766, 304428, 41135, -178037, -409437, -581990, -751733, -861540, -950970),
        36: (499229, 268474, 7526, -206231, -431366, -598039, -761680, -867186, -953008),
    }, {
        20: (22, Q('5/51'), (60, 898)),
        22: (24, Q('4/41'), (85, 1877)),
        24: (26, Q('4/41'), (139, 2121)),
        26: (28, Q('4/41'), (124, 577)),
        28: (30, Q('4/41'), (146, 1077)),
        30: (32, Q('3/31'), (174, 717)),
        32: (34, Q('3/31'), (491, 503)),
        34: (36, Q('3/31'), (3540, 834)),
    }),
    (Q('1/8'), 'cubic-precise', {10: ((49292, 87007), (-2898, 99958), (-47924, 87768), (-80653, 59120), (-97817, 20780)), 11: ((52446, 85143), (4403, 99903), (-38349, 92355), (-71683, 69725), (-92784, 37297), (-100000, 0))}, {
        36: (499229, 268474, 7526, -206231, -431366, -598039, -761680, -867186, -953008),
        38: (460832, 232160, -25917, -234124, -452903, -613751, -771382, -872684, -954988),
    }, {
        36: (38, Q('3/31'), (3077, 478)),
    }),
    (Q('1/8'), 'cubic', {10: ((45448, 89075), (-6167, 99810), (-49923, 86647), (-81449, 58018), (-97910, 20339)), 11: ((48622, 87384), (984, 99995), (-40658, 91362), (-72827, 68529), (-93088, 36533), (-100000, 0))}, {
        38: (460832, 232160, -25917, -234124, -452903, -613751, -771382, -872684, -954988),
        39: (441327, 213876, -42578, -247962, -463531, -621488, -776145, -875380, -955958),
        40: (418441, 195511, -59199, -261730, -474070, -629148, -780853, -878042, -956916),
        41: (398545, 177067, -75781, -275429, -484522, -636735, -785507, -880672, -957861),
        42: (378465, 158547, -92325, -289062, -494889, -644249, -790108, -883271, -958793),
    }, {
        38: (39, Q('3/31'), (741, 12514)),
        39: (40, Q('3/31'), (1005, 6201)),
        40: (41, Q('2/21'), (3470, 305)),
        41: (42, Q('2/21'), (4399, 1795)),
    }),
    (Q('93/1000'), 'cubic-scaled-precise-monotone-compress', {10: ((37528, 92691), (-12646, 99197), (-53811, 84287), (-82983, 55801), (-98088, 19463)), 11: ((40710, 91338), (-5819, 99831), (-45166, 89219), (-75038, 66100), (-93671, 35011), (-100000, 0))}, {
        42: (378465, 158547, -92325, -289062, -494889, -644249, -790108, -883271, -958793),
        43: (358205, 139953, -108830, -302630, -505174, -651694, -794659, -885839, -959715),
        44: (335589, 121286, -125297, -316134, -515379, -659071, -799161, -888378, -960624),
    }, {
        42: (43, Q('1/11'), (4137, 256)),
        43: (44, Q('1/11'), (14287, 275)),
    }),
    (Q('9/100'), 'cubic', {10: ((33460, 94236), (-15857, 98735), (-55706, 83047), (-83724, 54684), (-98173, 19027)), 11: ((36631, 93049), (-9203, 99576), (-47370, 88069), (-76108, 64865), (-93952, 34250), (-100000, 0))}, {
        44: (344112, 121286, -125297, -316134, -515379, -659071, -799161, -888378, -960624),
        45: (323488, 102550, -141727, -329576, -525505, -666382, -803615, -890888, -961523),
        46: (302692, 83746, -158121, -342957, -535556, -673628, -808023, -893371, -962412),
        47: (281730, 64875, -174478, -356278, -545533, -680812, -812386, -895826, -963290),
    }, {
        44: (45, Q('3/41'), (197, 1640)),
        45: (46, Q('3/41'), (249, 5424)),
        46: (47, Q('3/41'), (247, 6135)),
    }),
    (Q('9/100'), 'cubic-scaled-precise-monotone', {10: ((27231, 96221), (-20639, 97847), (-58490, 81111), (-84805, 52992), (-98298, 18374)), 11: ((30370, 95277), (-14257, 98979), (-50617, 86244), (-77672, 62984), (-94360, 33107), (-100000, 0))}, {
        470: (281730, 64875, -174478, -356278, -545533, -680812, -812386, -895826, -963290),
        475: (271187, 55415, -182644, -362917, -550494, -684382, -814551, -897044, -963725),
        480: (263727, 45940, -190800, -369542, -555437, -687936, -816705, -898255, -964158),
    }, {
        470: (475, Q('3/41'), (180, 340)),
        475: (480, Q('2/31'), (162, 134)),
    }),
    (Q('9/100'), 'cubic', {10: ((25124, 96793), (-22225, 97499), (-59404, 80444), (-85157, 52423), (-98338, 18156)), 11: ((28247, 95928), (-15936, 98722), (-51684, 85608), (-78184, 62348), (-94494, 32726), (-100000, 0))}, {
        48: (263727, 45940, -190800, -369542, -555437, -687936, -816705, -898255, -964158),
        49: (242425, 26942, -207087, -382748, -565271, -695000, -820981, -900659, -965016),
        50: (220964, 7884, -223339, -395898, -575036, -702006, -825216, -903037, -965865),
        51: (196284, -11235, -239557, -408994, -584734, -708956, -829411, -905392, -966705),
        52: (174542, -30411, -255741, -422036, -594366, -715851, -833566, -907723, -967536),
        53: (152653, -49644, -271891, -435025, -603934, -722693, -837683, -910031, -968358),
        54: (130622, -68932, -288009, -447963, -613439, -729482, -841763, -912317, -969171),
        55: (105494, -88274, -304094, -460851, -622883, -736220, -845806, -914581, -969977),
        56: (83217, -107669, -320148, -473689, -632267, -742908, -849813, -916824, -970774),
    }, {
        48: (49, Q('2/31'), (237, 1081)),
        49: (50, Q('2/31'), (244, 3759)),
        50: (51, Q('2/31'), (283, 1413)),
        51: (52, Q('2/31'), (338, 507)),
        52: (53, Q('2/31'), (375, 762)),
        53: (54, Q('2/31'), (436, 1100)),
        54: (55, Q('2/31'), (679, 1104)),
        55: (56, Q('2/31'), (793, 2451)),
    }),
    (Q('9/100'), 'cubic-scaled-precise-monotone', {10: ((7737, 99700), (-34759, 93765), (-66475, 74707), (-87856, 47763), (-98646, 16400)), 11: ((10662, 99430), (-29270, 95620), (-59978, 80016), (-82106, 57084), (-95509, 29633), (-100000, 0))}, {
        56: (83217, -107669, -320148, -473689, -632267, -742908, -849813, -916824, -970774),
        57: (60809, -127115, -336170, -486478, -641593, -749547, -853785, -919046, -971563),
    }, {
        56: (57, Q('2/31'), (905, 3668)),
    }),
    (Q('3/50'), 'cubic-scaled-precise-monotone-compress', {10: ((5502, 99849), (-36309, 93176), (-67332, 73936), (-88179, 47164), (-98683, 16177)), 11: ((8395, 99647), (-30925, 95098), (-60987, 79250), (-82578, 56400), (-95630, 29239), (-100000, 0))}, {
        57: (60809, -127115, -336170, -486478, -641593, -749547, -853785, -919046, -971563),
        58: (38270, -146611, -352160, -499219, -650861, -756138, -857724, -921247, -972345),
        59: (12785, -166157, -368120, -511914, -660072, -762682, -861630, -923429, -973119),
        60: (-9966, -185751, -384050, -524563, -669229, -769181, -865503, -925592, -973886),
        61: (-32837, -205391, -399949, -537166, -678332, -775635, -869344, -927736, -974646),
        62: (-55827, -225078, -415819, -549725, -687382, -782046, -873154, -929861, -975398),
    }, {
        57: (58, Q('3/61'), (763, 693)),
        58: (59, Q('3/61'), (3487, 821)),
        59: (60, Q('2/41'), (17453, 969)),
        60: (61, Q('2/41'), (6342, 974)),
        61: (62, Q('2/41'), (4129, 966)),
    }),
    (Q('3/50'), 'cubic-scaled-precise-monotone-compress', {10: ((-5853, 99829), (-44001, 89799), (-71532, 69880), (-89754, 44094), (-98861, 15049)), 11: ((-3154, 99950), (-39163, 92012), (-65944, 75175), (-84877, 52876), (-96219, 27239), (-100000, 0))}, {
        62: (-55827, -225078, -415819, -549725, -687382, -782046, -873154, -929861, -975398),
        63: (-80757, -244810, -431660, -562240, -696380, -788413, -876934, -931969, -976144),
        64: (-96934, -264586, -447471, -574712, -705327, -794738, -880685, -934058, -976884),
    }, {
        62: (63, Q('2/41'), (5858, 986)),
        63: (64, Q('2/41'), (2730, 1222)),
    }),
    (Q('3/40'), 'cubic-scaled-precise-monotone-compress', {10: ((-10476, 99450), (-47053, 88238), (-73176, 68156), (-90365, 42827), (-98930, 14589)), 11: ((-7868, 99690), (-42441, 90547), (-67889, 73424), (-85772, 51412), (-96447, 26421), (-100000, 0))}, {
        64: (-96934, -264586, -447471, -574712, -705327, -794738, -880685, -934058, -976884),
        65: (-120360, -284405, -463255, -587142, -714224, -801023, -884406, -936131, -977616),
        66: (-143899, -304267, -479010, -599531, -723072, -807266, -888099, -938186, -978343),
        67: (-167548, -324170, -494737, -611879, -731873, -813471, -891763, -940225, -979063),
    }, {
        64: (65, Q('2/41'), (1384, 2183)),
        65: (66, Q('2/41'), (1621, 1849)),
        66: (67, Q('2/41'), (1770, 2243)),
    }),
    (Q('3/40'), 'cubic-scaled-precise-monotone-compress', {10: ((-17493, 98458), (-51607, 85655), (-75605, 65451), (-91264, 40877), (-99031, 13886)), 11: ((-15032, 98864), (-47341, 88084), (-70769, 70652), (-87089, 49148), (-96781, 25167), (-100000, 0))}, {
        67: (-167548, -324170, -494737, -611879, -731873, -813471, -891763, -940225, -979063),
        68: (-191307, -344114, -510437, -624187, -740626, -819636, -895401, -942248, -979777),
        69: (-217528, -364098, -526109, -636455, -749333, -825764, -899012, -944255, -980485),
        70: (-241443, -384120, -541755, -648684, -757994, -831854, -902596, -946246, -981187),
    }, {
        67: (68, Q('2/41'), (1906, 2917)),
        68: (69, Q('2/41'), (2621, 3176)),
        69: (70, Q('2/41'), (3521, 3582)),
    }),
    (Q('11/200'), 'cubic-scaled-precise-monotone', {10: ((-24604, 96926), (-56131, 82761), (-77993, 62586), (-92142, 38857), (-99130, 13164)), 11: ((-22305, 97481), (-52220, 85282), (-73606, 67692), (-88378, 46791), (-97108, 23877), (-100000, 0))}, {
        70: (-241443, -384120, -541755, -648684, -757994, -831854, -902596, -946246, -981187),
        71: (-265461, -404182, -557374, -660875, -766611, -837907, -906155, -948223, -981883),
        72: (-289579, -424281, -572967, -673028, -775184, -843924, -909689, -950184, -982574),
        73: (-315920, -444417, -588534, -685145, -783714, -849906, -913198, -952130, -983259),
    }, {
        70: (71, Q('1/21'), (3938, 3717)),
        71: (72, Q('1/21'), (5032, 5314)),
        72: (73, Q('1/21'), (9661, 6685)),
    }),
    (Q('9/200'), 'cubic-scaled-precise-monotone-compress', {10: ((-31804, 94808), (-60628, 79525), (-80342, 59542), (-93000, 36756), (-99226, 12419)), 11: ((-29680, 95494), (-57079, 82110), (-76401, 64520), (-89639, 44326), (-97426, 22543), (-100000, 0))}, {
        73: (-315920, -444417, -588534, -685145, -783714, -849906, -913198, -952130, -983259),
        74: (-340172, -464589, -604075, -697224, -792202, -855853, -916683, -954062, -983939),
        75: (-365893, -484797, -619591, -709268, -800647, -861767, -920143, -955980, -984613),
    }, {
        73: (74, Q('2/51'), (11007, 4498)),
        74: (75, Q('2/51'), (22025, 5857)),
    }),
    (Q('9/200'), 'cubic-scaled-precise-monotone-compress', {10: ((-36652, 93041), (-63611, 77160), (-81887, 57398), (-93562, 35301), (-99289, 11906)), 11: ((-34652, 93804), (-60308, 79768), (-78243, 62274), (-90467, 42612), (-97634, 21624), (-100000, 0))}, {
        75: (-356520, -484797, -619591, -709268, -800647, -861767, -920143, -955980, -984613),
        76: (-381216, -505039, -635082, -721276, -809052, -867646, -923581, -957885, -985283),
    }, {
        75: (76, Q('2/61'), (4508, 10856)),
    }),
)


def verify(library, selected):
    total = 0
    by_mode = {
        "plain": replay_python if library is None else normalized_native(library),
        "scaled": replay_python if library is None else scaled_native(library),
        "compress": replay_python if library is None else scaled_native(library, compress=True),
    }
    prior = Q(0)
    for declared, mode, anchors, rays, expected in BANDS:
        knots = tuple(sorted(rays))
        knot_scale = 1000 if knots[0] >= 100 else 100
        lower, upper = Q(knots[0], knot_scale), Q(knots[-1], knot_scale)
        require(lower == prior < upper, ("band coverage", prior, lower, upper))
        prior = upper
        if selected and lower not in selected:
            continue
        small = mode == "small"
        verify_setup(rays, anchors, declared, small, knot_scale=knot_scale)
        strategy = "compress" if "compress" in mode else "scaled" if "scaled" in mode else "plain"
        for left, right in zip(knots, knots[1:]):
            cap = declared if small else exact_cap(declared, Q(left, knot_scale))
            target, recorded_cap, receipt = expected[left]
            require(right == target and cap == recorded_cap, ("knot or cap mismatch", left, right, cap))
            cells = tuple(
                verify_small_band(j, cap, rays, left, right, knot_scale=knot_scale) if small else
                verify_fused_band(j, cap, rays, left, right, by_mode[strategy], mode.startswith("cubic"),
                                  precise="precise" in mode, monotone="monotone" in mode,
                                  knot_scale=knot_scale)
                for j in ORDERS)
            if library is not None:
                require(cells == receipt, ("radial cell counts", left, right, cells, receipt))
            total += sum(cells)
            print(left, right, cap, cells, flush=True)
    require(prior == Q(19, 25), "terminal-eleven light upper endpoint")
    print(f"Verified {total} exact radial cells.", flush=True)


def main():
    selected = {Q(item) for item in argv[1:] if item != "--native"}
    starts = {Q(min(rays), 1000 if min(rays) >= 100 else 100) for _, _, _, rays, _ in BANDS}
    require(selected <= starts, "select band lower endpoints")
    if "--native" not in argv:
        verify(None, selected)
        return
    with TemporaryDirectory(prefix="sfc-eleven-light-") as directory:
        library = Path(directory)/"bernstein.so"
        source = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source),
             "-o", str(library)), check=True)
        verify(CDLL(str(library)), selected)


if __name__ == "__main__":
    main()
