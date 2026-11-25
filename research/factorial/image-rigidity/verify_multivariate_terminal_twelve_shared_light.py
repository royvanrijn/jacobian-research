#!/usr/bin/env python3
"""Exact terminal-twelve light-row separation on bounded outer-mean bands."""

from ctypes import CDLL
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_nine_high_curvature import anchor_value
from verify_multivariate_terminal_ten_shared_light import normalized_native, shape
from verify_multivariate_terminal_eleven_shared_light import scaled_native
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

LOWER, UPPER = ORDERS = (11, 12)
REFERENCE_RATIO = Q(121, 120)
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


def rows_above_floor(declared, outer_lower, floor):
    """Yield rows above the chosen floor whose squared means reach the slab."""
    require(outer_lower > 0 and floor > 0, "positive row enumeration parameters")
    threshold = 1 + (LOWER-1)*floor/(1+floor)
    gap = 1
    while threshold**(-gap) >= outer_lower**2:
        base = 1
        while LOWER*base+1 < gap/floor:
            reciprocal = Q(gap, LOWER*base+1)
            if reciprocal <= declared:
                row_square = prod(Q(LOWER*base+h, LOWER*base+(LOWER-1)*gap+h)
                                  for h in range(1, gap+1))
                if row_square >= outer_lower**2:
                    yield base, gap, reciprocal
            base += 1
        gap += 1


def exact_cap(declared, outer_lower):
    """Enumerate rows above FLOOR that can occur at this outer mean."""
    if not outer_lower:
        return declared
    floor = min(FLOOR, declared)
    return max((reciprocal for _, _, reciprocal in
                rows_above_floor(declared, outer_lower, floor)), default=floor)


@lru_cache(maxsize=None)
def ordinary_ratio_upper(cap, outer_lower):
    """Bound the log ratio by exact row ratios above the reciprocal floor."""
    floor = min(FLOOR, cap)
    result = Q(ratio_bounds(floor).hi, GRID)
    for base, gap, _ in rows_above_floor(cap, outer_lower, floor):
        lower = sum(log_bounds(1+Q((LOWER-1)*gap, LOWER*base+h))[0]
                    for h in range(1, gap+1))
        upper = sum(log_bounds(1+Q(LOWER*gap, UPPER*base+h))[1]
                    for h in range(1, gap+1))
        result = max(result, upper/lower)
    return result


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


@lru_cache(maxsize=512)
def log_power(value, exponent):
    """Enclose a power without exponentiating large rational denominators."""
    low, high = log_bounds(1/value)
    argument = rational(exponent)*interval(-high, -low)
    magnitude = max(abs(argument.lo), abs(argument.hi))
    shift = max(0, magnitude.bit_length() - 46)
    power = exp_small(argument*scalar(1, 2**shift))
    for _ in range(shift):
        power = power.square()
    return power


@lru_cache(maxsize=128)
def coefficient_model(order, lower_cap, means, signed, extended=False, *, precise=False, monotone=False, coupled=False,
                      ratio_ceiling=None, loss_transport=False, center_ratio=False):
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
    outer_adjustment = ZERO
    strength = lower_strength
    if order == UPPER:
        rho = ratio_bounds(lower_cap)
        if ratio_ceiling is not None:
            rho = I(rho.lo, min(rho.hi, rational(ratio_ceiling).hi))
        if monotone or coupled or loss_transport:
            minimum = REFERENCE_RATIO
            if coupled and lower_strength.hi > 0:
                signed_floor = max(Q(0), signed[0], -signed[1])
                kappa = Q(LOWER-1, UPPER)
                skew_per_loss = signed_floor * Q(Pmax.lo, GRID)
                minimum += (4*REFERENCE_RATIO*kappa**2*skew_per_loss**2
                            / ((LOWER**2-1)*Q(lower_strength.hi, GRID)))
            floor = rational(minimum).lo
            if floor > rho.hi:
                return None
            rho = I(floor, rho.hi)
        if center_ratio:
            ref_ratio = Q(rho.lo+rho.hi, 2*GRID)
        outer_adjustment = (rational(ref_ratio)-rho)*lower_strength
        strength = rho*lower_strength
        eta = lower_cap/(2*UPPER+(LOWER**2-2)*lower_cap)
        if loss_transport:
            c = lower_cap
            scale = 2*LOWER + (LOWER**2-1)*c
            b = (LOWER**2-1)*c/scale
            e = (LOWER-1)*(1+UPPER*c)/scale
            constant = REFERENCE_RATIO**2*Q(4*UPPER, LOWER-1)*b*e**2
            surplus = I(0, max(0, rho.hi-rational(REFERENCE_RATIO).lo))
            excess = rational(constant)*surplus*lower_strength
            transfer = rational(REFERENCE_RATIO**2)
            individual = rational(REFERENCE_RATIO**2*(2*eta+eta**2))*Pmax*lower_strength
            error = min(excess.hi, individual.hi)
            leading = transfer*leading+I(-error, error)
        else:
            transfer = rational(REFERENCE_RATIO**2*(1+eta**2))
            transfer_width = rational(2*REFERENCE_RATIO**2*eta)
            leading = transfer*leading+(transfer_width*Pmax*lower_strength).symmetric_error()
        nominal = transfer*nominal
    log_error = rational(ref_ratio)*magnitude.square()/(scalar(2)*(ONE-magnitude))
    outer_mean = I(-log_error.hi,0)+outer_adjustment
    gaussian, skew, mean, skew_sizes, distances = shape(order,skew_terms,mean_terms)
    direction, weights = leading_skew_direction(order,lower_cap,skew_terms,mean_terms)
    reference = []
    for k, v in enumerate(direction):
        if k == 0 or k == order:
            reference.append(ONE)
            continue
        exponent = ref_ratio*Q(k*(order-k), order-1)
        if center_ratio and order == UPPER:
            power = log_power(midpoint, exponent)
        elif precise:
            power = precise_power(midpoint, exponent)
        else:
            power = interval(*power_bounds(midpoint, exponent))
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


def verify_fused_band(order,cap,rays,left,right,replay,extended, *, precise=False, monotone=False, coupled=False, ordinary=False,
                      loss_transport=False, center_ratio=False, balanced=False,
                      indices=None, knot_scale=100):
    if indices is None:
        indices = tuple(range(UPPER-2))
    pending=[((Q(left,knot_scale),Q(right,knot_scale)),(Q(0),Q(1)),indices,0,0)]
    total=0
    ratio_ceiling = ordinary_ratio_upper(cap, Q(left,knot_scale)) if ordinary and order==UPPER else None
    while pending:
        means,p,active,mean_depth,skew_depth=pending.pop()
        model=coefficient_model(order,cap,means,p,extended,precise=precise,monotone=monotone,
                                coupled=coupled,ratio_ceiling=ratio_ceiling,loss_transport=loss_transport,
                                center_ratio=center_ratio)
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
        split_mean = mean_depth+2 < skew_depth if balanced else mean_depth < 2*skew_depth+2
        if split_mean:
            mid=sum(means)/2
            pending.extend((((means[0],mid),p,bad,mean_depth+1,skew_depth),
                            ((mid,means[1]),p,bad,mean_depth+1,skew_depth)))
        else:
            mid=sum(p)/2
            pending.extend(((means,(p[0],mid),bad,mean_depth,skew_depth+1),
                            (means,(mid,p[1]),bad,mean_depth,skew_depth+1)))
    return total


BANDS = (
    (Q('1/8'), 'small', {}, {
        0: (962709, 833300, 680984, 456840, 200567, -58252, -337117, -561261, -774180, -903590),
        1: (960937, 826907, 669194, 441224, 181456, -76677, -353344, -573274, -781088, -906728),
        2: (958491, 819286, 655753, 424580, 162229, -94480, -368357, -584107, -787175, -909444),
        3: (955307, 810679, 641226, 407444, 143222, -111648, -382463, -594142, -792740, -911907),
        4: (950469, 801227, 625916, 390048, 124526, -128232, -395833, -603561, -797914, -914183),
    }, {
        0: (1, Q('1/8'), (269, 424)),
        1: (2, Q('1/8'), (639, 1178)),
        2: (3, Q('1/8'), (1033, 3750)),
        3: (4, Q('1/8'), (4882, 9606)),
    }),
    (Q('1/8'), 'quad-scaled-precise-monotone-compress', {11: ((94695, 32139), (59632, 80274), (6528, 99787), (-46992, 88271), (-85863, 51260), (-100000, 0)), 12: ((95575, 29417), (65551, 75519), (18377, 98297), (-32174, 94683), (-73720, 67567), (-96974, 24415))}, {
        4: (951350, 801227, 625916, 390048, 124526, -128232, -395833, -603561, -797914, -914183),
        6: (941094, 780190, 593703, 354957, 88168, -159854, -420804, -620967, -807378, -918324),
        8: (927820, 756801, 560134, 319884, 53166, -189710, -443893, -636895, -815949, -922053),
        10: (911734, 731493, 525766, 285072, 19419, -218080, -465486, -651677, -823840, -925472),
        12: (893070, 704579, 490919, 250629, -13191, -245177, -485848, -665529, -831187, -928646),
        14: (872052, 676292, 455791, 216595, -44776, -271171, -505170, -678606, -838085, -931618),
        16: (847021, 646815, 420505, 182983, -75434, -296196, -523596, -691022, -844604, -934420),
    }, {
        4: (6, Q('1/8'), (98, 106)),
        6: (8, Q('1/8'), (75, 104)),
        8: (10, Q('2/17'), (73, 98)),
        10: (12, Q('2/17'), (81, 305)),
        12: (14, Q('2/17'), (93, 430)),
        14: (16, Q('2/17'), (146, 592)),
    }),
    (Q('1/8'), 'quad-scaled-precise-monotone-compress', {11: ((83955, 54329), (38693, 92211), (-12990, 99153), (-58470, 81125), (-89187, 45229), (-100000, 0)), 12: ((85823, 51327), (45408, 89096), (-2096, 99978), (-46249, 88663), (-79734, 60353), (-97697, 21337))}, {
        16: (847021, 646815, 420505, 182983, -75434, -296196, -523596, -691022, -844604, -934420),
        18: (819762, 616291, 385145, 149787, -105252, -320360, -541241, -702865, -850795, -937076),
        20: (790475, 584838, 349764, 116996, -134305, -343753, -558195, -714204, -856701, -939606),
        22: (759341, 552549, 314400, 84594, -162658, -366449, -574533, -725096, -862354, -942023),
    }, {
        16: (18, Q('2/17'), (162, 4649)),
        18: (20, Q('2/17'), (417, 10726)),
        20: (22, Q('2/17'), (3911, 1698)),
    }),
    (Q('1/8'), 'cubic-scaled-precise-coupled-ordinary-loss-compress', {11: ((75711, 65329), (28081, 95976), (-21369, 97690), (-62986, 77671), (-90437, 42674), (-100000, 0)), 12: ((77940, 62653), (34799, 93750), (-11163, 99375), (-51921, 85465), (-82033, 57189), (-97967, 20061))}, {
        220: (759341, 552549, 314400, 84594, -162658, -366449, -574533, -725096, -862354, -942023),
        221: (760084, 550915, 312633, 82984, -164059, -367567, -575335, -725630, -862630, -942141),
        240: (735820, 519507, 279079, 52563, -190368, -388513, -590315, -735588, -867781, -944341),
        260: (701797, 485779, 243818, 20886, -217483, -409999, -605595, -745718, -873005, -946568),
    }, {
        220: (221, Q('2/17'), (1060, 424)),
        221: (240, Q('4/45'), (169, 125)),
        240: (260, Q('4/45'), (416, 180)),
    }),
    (Q('1/8'), 'cubic-scaled-precise-coupled-compress', {11: ((69457, 71942), (21065, 97756), (-26609, 96395), (-65728, 75365), (-91185, 41051), (-100000, 0)), 12: ((71865, 69537), (27699, 96087), (-16888, 98564), (-55391, 83258), (-83416, 55153), (-98128, 19258))}, {
        26: (701797, 485779, 243818, 20886, -217483, -409999, -605595, -745718, -873005, -946568),
        27: (685455, 468676, 226214, 5175, -230831, -420540, -613059, -750657, -875548, -947651),
        28: (668801, 451422, 208630, -10455, -244047, -430953, -620414, -755518, -878046, -948715),
        29: (651844, 434025, 191066, -26004, -257135, -441244, -627664, -760304, -880502, -949760),
        30: (634595, 416488, 173524, -41475, -270100, -451417, -634813, -765018, -882918, -950788),
        31: (614495, 398818, 156003, -56870, -282944, -461476, -641865, -769662, -885296, -951799),
        32: (596665, 381020, 138505, -72191, -295674, -471425, -648824, -774241, -887637, -952793),
        33: (578572, 363097, 121030, -87439, -308291, -481269, -655692, -778755, -889943, -953772),
        34: (560224, 345056, 103579, -102617, -320800, -491010, -662475, -783208, -892214, -954736),
    }, {
        26: (27, Q('4/45'), (299, 117)),
        27: (28, Q('4/45'), (133, 137)),
        28: (29, Q('4/45'), (144, 658)),
        29: (30, Q('3/34'), (150, 17096)),
        30: (31, Q('3/34'), (470, 10767)),
        31: (32, Q('3/34'), (659, 268)),
        32: (33, Q('3/34'), (469, 342)),
        33: (34, Q('3/34'), (439, 21893)),
    }),
    (Q('1/8'), 'cubic-scaled-precise-coupled-compress', {11: ((55496, 83188), (7200, 99740), (-36436, 93126), (-70729, 70692), (-92530, 37922), (-100000, 0)), 12: ((58129, 81370), (13515, 99082), (-27724, 96080), (-61766, 78645), (-85912, 51177), (-98417, 17724))}, {
        34: (560224, 345056, 103579, -102617, -320800, -491010, -662475, -783208, -892214, -954736),
        35: (538979, 326898, 86151, -117725, -333205, -500653, -669173, -787602, -894453, -955686),
        36: (520129, 308630, 68746, -132766, -345507, -510199, -675792, -791939, -896660, -956622),
    }, {
        34: (35, Q('3/34'), (1150, 24375)),
        35: (36, Q('3/34'), (3806, 31027)),
    }),
    (Q('3/50'), 'cubic-scaled-precise-coupled-ordinary-loss-centered-compress-balanced', {11: ((51746, 85570), (3770, 99929), (-38778, 92175), (-71897, 69504), (-92841, 37155), (-100000, 0)), 12: ((54410, 83902), (9979, 99501), (-30323, 95292), (-63261, 77447), (-86491, 50193), (-98483, 17350))}, {
        36: (530784, 308630, 68746, -132766, -345507, -510199, -675792, -791939, -896660, -956622),
        37: (511751, 290253, 51366, -147741, -357710, -519654, -682332, -796220, -898837, -957545),
    }, {
        36: (37, Q('4/67'), (158, 4097)),
    }),
    (Q('3/50'), 'cubic-scaled-precise-coupled-ordinary-loss-centered-compress-balanced', {11: ((49837, 86696), (2060, 99979), (-39934, 91680), (-72470, 68907), (-92993, 36773), (-100000, 0)), 12: ((52513, 85102), (8213, 99662), (-31608, 94873), (-63997, 76840), (-86774, 49701), (-98516, 17165))}, {
        37: (511751, 290253, 51366, -147741, -357710, -519654, -682332, -796220, -898837, -957545),
        38: (489799, 271772, 34010, -162651, -369818, -529018, -688797, -800449, -900984, -958454),
    }, {
        37: (38, Q('4/67'), (179, 2362)),
    }),
    (Q('3/50'), 'cubic-scaled-precise-coupled-ordinary-loss-centered-compress', {11: ((47906, 87778), (353, 99999), (-41080, 91172), (-73036, 68306), (-93143, 36392), (-100000, 0)), 12: ((50591, 86259), (6449, 99792), (-32883, 94439), (-64724, 76229), (-87054, 49210), (-98548, 16980))}, {
        38: (489799, 271772, 34010, -162651, -369818, -529018, -688797, -800449, -900984, -958454),
        39: (470297, 253190, 16677, -177498, -381832, -538296, -695189, -804626, -903103, -959352),
        40: (450576, 234510, -631, -192283, -393756, -547489, -701511, -808753, -905195, -960237),
    }, {
        38: (39, Q('4/67'), (190, 3088)),
        39: (40, Q('4/67'), (193, 2100)),
    }),
    (Q('3/40'), 'cubic-scaled-precise-coupled-compress', {11: ((43979, 89810), (-3050, 99953), (-43344, 90118), (-74148, 67097), (-93437, 35631), (-100000, 0)), 12: ((46676, 88438), (2923, 99957), (-35407, 93522), (-66154, 74991), (-87602, 48227), (-98611, 16611))}, {
        40: (447879, 234510, -631, -192283, -393756, -547489, -701511, -808753, -905195, -960237),
        41: (427942, 215735, -17915, -207009, -405591, -556601, -707764, -812832, -907260, -961111),
        42: (407799, 196868, -35174, -221675, -417341, -565632, -713952, -816865, -909300, -961973),
        43: (387456, 177911, -52410, -236283, -429006, -574587, -720075, -820852, -911314, -962825),
        44: (364221, 158868, -69622, -250835, -440591, -583466, -726135, -824795, -913305, -963666),
        45: (343498, 139739, -86810, -265331, -452096, -592272, -732135, -828696, -915272, -964497),
        46: (322590, 120529, -103974, -279773, -463523, -601007, -738076, -832555, -917217, -965318),
    }, {
        40: (41, Q('1/15'), (290, 1826)),
        41: (42, Q('1/15'), (317, 1179)),
        42: (43, Q('1/15'), (502, 964)),
        43: (44, Q('1/15'), (684, 1005)),
        44: (45, Q('1/15'), (923, 792)),
        45: (46, Q('1/15'), (1015, 3783)),
    }),
    (Q('3/40'), 'cubic-scaled-precise-coupled-ordinary-loss-compress', {11: ((31722, 94835), (-13178, 99128), (-49928, 86644), (-77342, 63389), (-94274, 33352), (-100000, 0)), 12: ((34405, 93895), (-7617, 99710), (-42777, 90389), (-70273, 71145), (-89169, 45265), (-98789, 15513))}, {
        460: (322590, 120529, -103974, -279773, -463523, -601007, -738076, -832555, -917217, -965318),
        470: (301502, 101238, -121114, -294162, -474875, -609672, -743960, -836374, -919139, -966129),
    }, {
        460: (470, Q('1/15'), (1083, 46321)),
    }),
    (Q('13/200'), 'cubic-scaled-precise-coupled-ordinary-loss-centered-compress', {11: ((29616, 95514), (-14855, 98891), (-50998, 86019), (-77855, 62758), (-94408, 32971), (-100000, 0)), 12: ((32289, 94644), (-9368, 99560), (-43978, 89811), (-70937, 70484), (-89420, 44767), (-98818, 15331))}, {
        47: (304175, 101238, -121114, -294162, -474875, -609672, -743960, -836374, -919139, -966129),
        48: (282901, 81869, -138231, -308499, -486154, -618270, -749788, -840155, -921041, -966932),
    }, {
        47: (48, Q('1/17'), (539, 3981)),
    }),
    (Q('3/40'), 'cubic-scaled-precise-coupled-compress', {11: ((27492, 96147), (-16529, 98625), (-52060, 85380), (-78363, 62123), (-94540, 32590), (-100000, 0)), 12: ((30153, 95346), (-11118, 99380), (-45171, 89216), (-71594, 69816), (-89668, 44268), (-98846, 15148))}, {
        48: (282901, 81869, -138231, -308499, -486154, -618270, -749788, -840155, -921041, -966932),
        49: (261453, 62424, -155325, -322785, -497360, -626802, -755562, -843897, -922921, -967725),
        50: (239835, 42906, -172395, -337020, -508496, -635269, -761283, -847602, -924781, -968509),
    }, {
        48: (49, Q('1/17'), (533, 1730)),
        49: (50, Q('1/17'), (661, 88969)),
    }),
    (Q('1/20'), 'cubic-scaled-precise-coupled-ordinary-loss-centered-compress', {11: ((23194, 97273), (-19867, 98007), (-54162, 84062), (-79365, 60838), (-94800, 31826), (-100000, 0)), 12: ((25827, 96607), (-14612, 98927), (-47537, 87978), (-72892, 68460), (-90156, 43266), (-98901, 14782))}, {
        50: (245102, 42906, -172395, -337020, -508496, -635269, -761283, -847602, -924781, -968509),
        51: (223284, 23315, -189443, -351206, -519563, -643674, -766952, -851271, -926622, -969285),
        52: (201301, 3653, -206467, -365345, -530564, -652017, -772572, -854905, -928443, -970052),
        53: (179155, -16077, -223469, -379435, -541498, -660301, -778142, -858505, -930246, -970811),
    }, {
        50: (51, Q('3/67'), (332, 2264)),
        51: (52, Q('3/67'), (365, 8074)),
        52: (53, Q('3/67'), (382, 27700)),
    }),
    (Q('1/20'), 'cubic-scaled-precise-coupled-ordinary-loss-centered-compress', {11: ((16627, 98608), (-24851, 96863), (-57264, 81981), (-80832, 58874), (-95180, 30671), (-100000, 0)), 12: ((19204, 98139), (-19842, 98012), (-51036, 85996), (-74796, 66374), (-90869, 41747), (-98982, 14231))}, {
        53: (179155, -16077, -223469, -379435, -541498, -660301, -778142, -858505, -930246, -970811),
        54: (154297, -35874, -240448, -393479, -552368, -668527, -783665, -862072, -932031, -971563),
        55: (131861, -55736, -257404, -407478, -563176, -676695, -789141, -865606, -933798, -972307),
        56: (109274, -75663, -274338, -421431, -573921, -684807, -794572, -869109, -935548, -973043),
    }, {
        53: (54, Q('3/67'), (469, 6782)),
        54: (55, Q('3/67'), (563, 1685)),
        55: (56, Q('3/67'), (598, 1997)),
    }),
    (Q('1/20'), 'cubic-scaled-precise-coupled-compress', {11: ((9926, 99506), (-29809, 95454), (-60308, 79768), (-82260, 56862), (-95548, 29505), (-100000, 0)), 12: ((12430, 99224), (-25058, 96810), (-54477, 83859), (-76654, 64220), (-91561, 40206), (-99060, 13677))}, {
        56: (109274, -75663, -274338, -421431, -573921, -684807, -794572, -869109, -935548, -973043),
        57: (86540, -95652, -291251, -435341, -584607, -692865, -799958, -872580, -937281, -973772),
        58: (63661, -115703, -308141, -449207, -595233, -700870, -805301, -876021, -938997, -974493),
        59: (40640, -135813, -325010, -463030, -605802, -708822, -810601, -879433, -940698, -975208),
        60: (15094, -155981, -341857, -476812, -616314, -716723, -815860, -882816, -942383, -975916),
        61: (-8169, -176207, -358683, -490552, -626770, -724574, -821079, -886171, -944052, -976618),
        62: (-31564, -196489, -375487, -504252, -637172, -732376, -826258, -889498, -945707, -977312),
    }, {
        56: (57, Q('3/67'), (702, 3109)),
        57: (58, Q('3/67'), (1013, 3973)),
        58: (59, Q('2/45'), (1416, 5247)),
        59: (60, Q('2/45'), (3532, 5959)),
        60: (61, Q('2/45'), (7488, 8695)),
        61: (62, Q('2/45'), (6635, 10042)),
    }),
    (Q('1/20'), 'cubic-scaled-precise-coupled-compress', {11: ((-3851, 99926), (-39651, 91803), (-66235, 74920), (-85011, 52661), (-96253, 27119), (-100000, 0)), 12: ((-1535, 99988), (-35447, 93507), (-61200, 79086), (-80241, 59678), (-92889, 37036), (-99210, 12547))}, {
        62: (-31564, -196489, -375487, -504252, -637172, -732376, -826258, -889498, -945707, -977312),
        63: (-55090, -216826, -392271, -517913, -647520, -740130, -831397, -892798, -947347, -978001),
        64: (-80983, -237217, -409034, -531534, -657816, -747837, -836499, -896072, -948972, -978683),
        65: (-104722, -257661, -425776, -545116, -668061, -755497, -841564, -899320, -950584, -979359),
        66: (-128583, -278156, -442497, -558660, -678254, -763113, -846592, -902542, -952182, -980030),
    }, {
        62: (63, Q('2/45'), (5127, 10322)),
        63: (64, Q('2/45'), (5018, 11942)),
        64: (65, Q('2/45'), (4675, 15195)),
        65: (66, Q('2/45'), (3365, 24333)),
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
                                  coupled="coupled" in mode, ordinary="ordinary" in mode,
                                  loss_transport="loss" in mode, center_ratio="centered" in mode,
                                  balanced="balanced" in mode, knot_scale=knot_scale)
                for j in ORDERS)
            if library is not None:
                require(cells == receipt, ("radial cell counts", left, right, cells, receipt))
            total += sum(cells)
            print(left, right, cap, cells, flush=True)
    require(prior == Q(33, 50), "terminal-twelve light upper endpoint")
    print(f"Verified {total} exact radial cells.", flush=True)


def main():
    selected = {Q(item) for item in argv[1:] if item != "--native"}
    starts = {Q(min(rays), 1000 if min(rays) >= 100 else 100) for _, _, _, rays, _ in BANDS}
    require(selected <= starts, "select band lower endpoints")
    if "--native" not in argv:
        verify(None, selected)
        return
    with TemporaryDirectory(prefix="sfc-twelve-light-") as directory:
        library = Path(directory)/"bernstein.so"
        source = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source),
             "-o", str(library)), check=True)
        verify(CDLL(str(library)), selected)


if __name__ == "__main__":
    main()
