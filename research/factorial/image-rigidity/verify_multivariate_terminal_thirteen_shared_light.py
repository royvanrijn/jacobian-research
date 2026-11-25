#!/usr/bin/env python3
"""Exact terminal-thirteen light-row separation on bounded outer-mean bands."""

from ctypes import CDLL
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_nine_high_curvature import anchor_value
from verify_multivariate_terminal_ten_shared_light import normalized_native, shape
from verify_multivariate_terminal_eleven_shared_light import scaled_native
from verify_multivariate_terminal_nine_shared_matched_upper import ANNULUS,replay_python
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

LOWER, UPPER = ORDERS = (12, 13)
REFERENCE_RATIO = Q(LOWER**2, LOWER**2 - 1)
FLOOR=Q(1,40)
SKEW_CENTERS = ((Q(3,4), Q(1,4)), (Q(2,3), Q(1,3)), (Q(5,8), Q(3,8)))

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
def leading_skew_direction(order, lower_cap, skew_terms, mean_terms, center_skew):
    _, skew, _, _, _ = shape(order, skew_terms, mean_terms)
    T, _, _ = constants(order, lower_cap)
    centers = (tuple(pair[0] for pair in SKEW_CENTERS) if center_skew
               else (Q(1,2),) * len(SKEW_CENTERS))
    weights = tuple(T**(2*l-2)*rational(centers[l-2]/l) for l in range(2, skew_terms+1))
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
                      ratio_ceiling=None, loss_transport=False, center_ratio=False, center_skew=False):
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
    direction, weights = leading_skew_direction(order,lower_cap,skew_terms,mean_terms,center_skew)
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
    for l, (vector, weight) in enumerate(zip(skew[1:], weights, strict=True), 2):
        bound = (rational(SKEW_CENTERS[l-2][1]/(2*l))*T**(2*l-1)*strength
                 if center_skew else weight*skew_bound*strength)
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
                      loss_transport=False, center_ratio=False, balanced=False, center_skew=False,
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
                                center_ratio=center_ratio,center_skew=center_skew)
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
        0: (968434, 857218, 727809, 530856, 306712, 69141, -189678, -413442, -637586, -796281, -925691),
        1: (966941, 851581, 717358, 516373, 288267, 50273, -207815, -428327, -648094, -802623, -928126),
        2: (964886, 844792, 705235, 500657, 269295, 31660, -224939, -442001, -657480, -808186, -930243),
        3: (961554, 837029, 691917, 484246, 250260, 13480, -241213, -454795, -666121, -813260, -932163),
    }, {
        0: (1, Q('1/8'), (488, 946)),
        1: (2, Q('1/8'), (1355, 3668)),
        2: (3, Q('1/8'), (5817, 16383)),
    }),
    (Q('1/8'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((95905, 28325), (66882, 74342), (20085, 97962), (-30854, 95121), (-73119, 68218), (-96899, 24708)), 13: ((96532, 26108), (71501, 69911), (29967, 95404), (-17389, 98476), (-60105, 79921), (-89533, 44540), (-100000, 0))}, {
        3: (962181, 837029, 691917, 484246, 250260, 13480, -241213, -454795, -666121, -813260, -932163),
        4: (958767, 828403, 677689, 467397, 231332, -4246, -256764, -466885, -674194, -817970, -933937),
        5: (954611, 819005, 662755, 450271, 212592, -21528, -271687, -478388, -681806, -822390, -935596),
        6: (949704, 808910, 647266, 432972, 194079, -38387, -286056, -489388, -689033, -826571, -937161),
        7: (944050, 798185, 631337, 415572, 175812, -54848, -299933, -499949, -695930, -830548, -938645),
        8: (937665, 786887, 615054, 398122, 157797, -70936, -313366, -510123, -702537, -834347, -940061),
        9: (930568, 775065, 598486, 380658, 140033, -86676, -326397, -519949, -708890, -837992, -941416),
        10: (922784, 762761, 581684, 363205, 122516, -102088, -339062, -529462, -715015, -841498, -942717),
        11: (913133, 750011, 564692, 345781, 105239, -117194, -351391, -538690, -720933, -844879, -943970),
        12: (903979, 736850, 547541, 328400, 88195, -132013, -363409, -547658, -726665, -848149, -945180),
        13: (892869, 723305, 530260, 311071, 71374, -146561, -375141, -556386, -732226, -851315, -946350),
        14: (882457, 709401, 512869, 293801, 54769, -160854, -386606, -564893, -737629, -854388, -947484),
        15: (870006, 695163, 495387, 276594, 38371, -174906, -397821, -573194, -742887, -857374, -948585),
        16: (856897, 680610, 477829, 259455, 22172, -188731, -408803, -581304, -748010, -860279, -949655),
        17: (844767, 665761, 460207, 242386, 6164, -202339, -419566, -589234, -753008, -863110, -950696),
        18: (832155, 650632, 442530, 225387, -9661, -215742, -430123, -596996, -757888, -865872, -951711),
    }, {
        3: (4, Q('6/49'), (57, 92)),
        4: (5, Q('6/49'), (54, 55)),
        5: (6, Q('6/49'), (58, 75)),
        6: (7, Q('6/49'), (57, 71)),
        7: (8, Q('6/49'), (58, 67)),
        8: (9, Q('6/49'), (59, 102)),
        9: (10, Q('3/25'), (62, 979)),
        10: (11, Q('3/25'), (67, 861)),
        11: (12, Q('3/25'), (85, 137)),
        12: (13, Q('3/25'), (122, 141)),
        13: (14, Q('3/25'), (123, 1118)),
        14: (15, Q('3/25'), (198, 3229)),
        15: (16, Q('3/25'), (631, 1525)),
        16: (17, Q('3/25'), (698, 684)),
        17: (18, Q('3/25'), (2316, 1597)),
    }),
    (Q('3/40'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((83049, 55704), (41428, 91015), (-5581, 99844), (-48457, 87475), (-80636, 59143), (-97803, 20845)), 13: ((84716, 53133), (47078, 88225), (3649, 99933), (-37567, 92675), (-70942, 70479), (-92539, 37902), (-100000, 0))}, {
        18: (838825, 650632, 442530, 225387, -9661, -215742, -430123, -596996, -757888, -865872, -951711),
        19: (825974, 635238, 424809, 208459, -25310, -228951, -440485, -604600, -762658, -868568, -952700),
        20: (812667, 619594, 407051, 191602, -40789, -241974, -450663, -612055, -767325, -871203, -953667),
        21: (798924, 603712, 389261, 174817, -56106, -254820, -460667, -619369, -771894, -873780, -954611),
        22: (784761, 587604, 371446, 158103, -71266, -267496, -470504, -626549, -776371, -876303, -955535),
        23: (770193, 571279, 353609, 141459, -86275, -280011, -480185, -633603, -780761, -878775, -956439),
        24: (753273, 554749, 335756, 124883, -101139, -292371, -489714, -640537, -785068, -881198, -957324),
    }, {
        18: (19, Q('6/85'), (69, 109)),
        19: (20, Q('5/73'), (95, 86)),
        20: (21, Q('5/73'), (127, 109)),
        21: (22, Q('5/73'), (99, 172)),
        22: (23, Q('5/73'), (80, 4994)),
        23: (24, Q('5/73'), (89, 3024)),
    }),
    (Q('3/40'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((74542, 66659), (30744, 95157), (-14461, 98949), (-53929, 84212), (-82836, 56020), (-98061, 19598)), 13: ((76505, 64398), (36407, 93137), (-5767, 99834), (-44014, 89793), (-74178, 67064), (-93404, 35717), (-100000, 0))}, {
        24: (753273, 554749, 335756, 124883, -101139, -292371, -489714, -640537, -785068, -881198, -957324),
        25: (737900, 538022, 317888, 108377, -115862, -304583, -499101, -647356, -789297, -883575, -958192),
        26: (722167, 521106, 300011, 91937, -130451, -316652, -508350, -654066, -793450, -885908, -959044),
        27: (706086, 504010, 282125, 75564, -144908, -328584, -517469, -660671, -797533, -888199, -959879),
        28: (689670, 486740, 264234, 59255, -159239, -340384, -526461, -667176, -801547, -890450, -960700),
        29: (672929, 469304, 246340, 43011, -173448, -352057, -535333, -673585, -805496, -892663, -961505),
        30: (655874, 451708, 228444, 26829, -187538, -363608, -544089, -679903, -809382, -894840, -962297),
        31: (636324, 433958, 210548, 10708, -201514, -375041, -552734, -686133, -813209, -896982, -963076),
        32: (618647, 416059, 192653, -5352, -215380, -386360, -561271, -692277, -816979, -899090, -963842),
        33: (600686, 398017, 174761, -21353, -229137, -397569, -569705, -698341, -820693, -901166, -964596),
        34: (580199, 379837, 156873, -37296, -242790, -408671, -578039, -704325, -824354, -903211, -965339),
    }, {
        24: (25, Q('5/73'), (123, 162)),
        25: (26, Q('5/73'), (138, 1861)),
        26: (27, Q('4/61'), (165, 6390)),
        27: (28, Q('4/61'), (172, 11885)),
        28: (29, Q('4/61'), (174, 7584)),
        29: (30, Q('4/61'), (179, 3994)),
        30: (31, Q('4/61'), (196, 943)),
        31: (32, Q('4/61'), (245, 1650)),
        32: (33, Q('4/61'), (246, 6229)),
        33: (34, Q('4/61'), (504, 1716)),
    }),
    (Q('3/50'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((57570, 81766), (12977, 99154), (-28122, 95964), (-61996, 78464), (-86001, 51027), (-98427, 17667)), 13: ((59821, 80134), (18398, 98293), (-20436, 97890), (-53612, 84414), (-78869, 61479), (-94641, 32298), (-100000, 0))}, {
        34: (584702, 379837, 156873, -37296, -242790, -408671, -578039, -704325, -824354, -903211, -965339),
        35: (566216, 361523, 138990, -53182, -256342, -419670, -586277, -710234, -827965, -905226, -966070),
        36: (547468, 343080, 121112, -69012, -269795, -430570, -594422, -716070, -831526, -907213, -966790),
        37: (528466, 324512, 103240, -84788, -283152, -441372, -602478, -721836, -835039, -909171, -967500),
        38: (509216, 305823, 85376, -100510, -296416, -452081, -610446, -727533, -838507, -911104, -968200),
    }, {
        34: (35, Q('4/73'), (235, 342)),
        35: (36, Q('4/73'), (244, 420)),
        36: (37, Q('4/73'), (280, 565)),
        37: (38, Q('4/73'), (737, 1110)),
    }),
    (Q('3/50'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((50000, 86603), (5911, 99825), (-33271, 94303), (-64944, 76042), (-87138, 49060), (-98558, 16924)), 13: ((52304, 85231), (11165, 99375), (-26013, 96557), (-57146, 82063), (-80563, 59242), (-95082, 30973), (-100000, 0))}, {
        38: (509216, 305823, 85376, -100510, -296416, -452081, -610446, -727533, -838507, -911104, -968200),
        39: (487414, 287016, 67520, -116179, -309590, -462699, -618330, -733165, -841930, -913010, -968890),
        40: (467684, 268096, 49672, -131797, -322675, -473228, -626133, -738732, -845311, -914892, -969570),
        41: (447728, 249065, 31832, -147364, -335675, -483671, -633856, -744238, -848650, -916750, -970242),
        42: (427551, 229927, 14002, -162881, -348591, -494031, -641503, -749684, -851950, -918584, -970905),
        43: (407160, 210685, -3818, -178349, -361425, -504309, -649075, -755072, -855210, -920396, -971559),
    }, {
        38: (39, Q('4/73'), (971, 1062)),
        39: (40, Q('4/73'), (1000, 549)),
        40: (41, Q('2/37'), (887, 631)),
        41: (42, Q('2/37'), (660, 969)),
        42: (43, Q('2/37'), (725, 13987)),
    }),
    (Q('1/20'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((40018, 91643), (-2882, 99958), (-39497, 91869), (-68450, 72901), (-88478, 46601), (-98711, 16007)), 13: ((42344, 90593), (2118, 99978), (-32788, 94472), (-61365, 78958), (-82564, 56420), (-95601, 29333), (-100000, 0))}, {
        43: (409485, 210685, -3818, -178349, -361425, -504309, -649075, -755072, -855210, -920396, -971559),
        44: (386560, 191341, -21629, -193769, -374179, -514509, -656575, -760403, -858433, -922186, -972205),
        45: (365755, 171900, -39429, -209143, -386856, -524631, -664004, -765680, -861619, -923955, -972844),
        46: (344750, 152362, -57218, -224470, -399457, -534678, -671365, -770904, -864770, -925703, -973474),
        47: (323551, 132731, -74997, -239751, -411984, -544652, -678659, -776075, -867886, -927432, -974097),
    }, {
        43: (44, Q('3/61'), (598, 5964)),
        44: (45, Q('3/61'), (753, 1089)),
        45: (46, Q('3/61'), (810, 1253)),
        46: (47, Q('3/61'), (1153, 4542)),
    }),
    (Q('1/22'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew-balanced', {12: ((31663, 94855), (-9883, 99510), (-44329, 89638), (-71131, 70288), (-89493, 44621), (-98826, 15277)), 13: ((33971, 94053), (-5116, 99869), (-38067, 92471), (-64601, 76333), (-84084, 54128), (-95993, 28023), (-100000, 0))}, {
        47: (323551, 132731, -74997, -239751, -411984, -544652, -678659, -776075, -867886, -927432, -974097),
        48: (302161, 113009, -92764, -254989, -424439, -554555, -685889, -781197, -870969, -929141, -974713),
        49: (278295, 93199, -110520, -270182, -436823, -564388, -693056, -786269, -874019, -930831, -975322),
        50: (256550, 73302, -128265, -285332, -449138, -574153, -700161, -791294, -877038, -932503, -975923),
    }, {
        47: (48, Q('4/97'), (3756, 6391)),
        48: (49, Q('3/73'), (5733, 7737)),
        49: (50, Q('3/73'), (6592, 8827)),
    }),
    (Q('1/20'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((25199, 96773), (-15114, 98851), (-47875, 87795), (-73076, 68263), (-90225, 43121), (-98909, 14730)), 13: ((27477, 96151), (-10539, 99443), (-41953, 90774), (-66956, 74276), (-85182, 52383), (-96275, 27038), (-100000, 0))}, {
        50: (256550, 73302, -128265, -285332, -449138, -574153, -700161, -791294, -877038, -932503, -975923),
        51: (234628, 53321, -145998, -300440, -461385, -583852, -707206, -796273, -880025, -934157, -976519),
        52: (212533, 33257, -163719, -315506, -473566, -593487, -714193, -801206, -882983, -935793, -977107),
    }, {
        50: (51, Q('3/61'), (2061, 2483)),
        51: (52, Q('3/61'), (2576, 5317)),
    }),
    (Q('1/22'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew-balanced', {12: ((20804, 97812), (-18592, 98257), (-50204, 86485), (-74345, 66879), (-90700, 42112), (-98963, 14364)), 13: ((23052, 97307), (-14152, 98994), (-44510, 89548), (-68494, 72860), (-85896, 51204), (-96458, 26378), (-100000, 0))}, {
        52: (217030, 33257, -163719, -315506, -473566, -593487, -714193, -801206, -882983, -935793, -977107),
        53: (192502, 13114, -181428, -330532, -485682, -603057, -721123, -806096, -885912, -937413, -977689),
        54: (170055, -7108, -199125, -345517, -497735, -612567, -727997, -810942, -888812, -939017, -978266),
    }, {
        52: (53, Q('3/73'), (11442, 23004)),
        53: (54, Q('3/73'), (14122, 26082)),
    }),
    (Q('1/30'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((16341, 98656), (-22062, 97536), (-52506, 85107), (-75592, 65467), (-91166, 41095), (-99016, 13996)), 13: ((18555, 98264), (-17763, 98410), (-47041, 88245), (-70008, 71407), (-86597, 50010), (-96637, 25714), (-100000, 0))}, {
        54: (170055, -7108, -199125, -345517, -497735, -612567, -727997, -810942, -888812, -939017, -978266),
        55: (147444, -27406, -216811, -360463, -509726, -622015, -734817, -815746, -891685, -940604, -978836),
    }, {
        54: (55, Q('4/121'), (1699, 2185)),
    }),
    (Q('1/20'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((14086, 99003), (-23795, 97128), (-53647, 84392), (-76208, 64749), (-91395, 40582), (-99042, 13812)), 13: ((16280, 98666), (-19568, 98067), (-48298, 87563), (-70756, 70666), (-86942, 49408), (-96726, 25380), (-100000, 0))}, {
        55: (147444, -27406, -216811, -360463, -509726, -622015, -734817, -815746, -891685, -940604, -978836),
        56: (124674, -47779, -234484, -375369, -521656, -631405, -741583, -820510, -894530, -942176, -979400),
    }, {
        55: (56, Q('3/73'), (3328, 5245)),
    }),
    (Q('1/25'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew-balanced', {12: ((11816, 99299), (-25525, 96687), (-54783, 83659), (-76818, 64023), (-91622, 40067), (-99067, 13627)), 13: ((13988, 99017), (-21372, 97690), (-49549, 86862), (-71498, 69914), (-87284, 48801), (-96813, 25046), (-100000, 0))}, {
        56: (124674, -47779, -234484, -375369, -521656, -631405, -741583, -820510, -894530, -942176, -979400),
        57: (101746, -68226, -252144, -390237, -533527, -640736, -748298, -825233, -897349, -943733, -979959),
    }, {
        56: (57, Q('1/25'), (28554, 61986)),
    }),
    (Q('1/40'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((9530, 99545), (-27254, 96214), (-55912, 82909), (-77424, 63289), (-91847, 39548), (-99093, 13441)), 13: ((11679, 99316), (-23175, 97278), (-50794, 86139), (-72235, 69152), (-87622, 48190), (-96899, 24709), (-100000, 0))}, {
        57: (101746, -68226, -252144, -390237, -533527, -640736, -748298, -825233, -897349, -943733, -979959),
        58: (76540, -88743, -269793, -405068, -545339, -650011, -754962, -829917, -900142, -945274, -980512),
        59: (53334, -109331, -287429, -419861, -557093, -659231, -761577, -834563, -902910, -946802, -981059),
        60: (29981, -129988, -305053, -434618, -568792, -668396, -768143, -839172, -905654, -948315, -981602),
    }, {
        57: (58, Q('1/40'), (2286, 3251)),
        58: (59, Q('1/40'), (3056, 3758)),
        59: (60, Q('1/40'), (3797, 5029)),
    }),
    (Q('1/36'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew-balanced', {12: ((2584, 99967), (-32429, 94596), (-59264, 80547), (-79214, 61035), (-92510, 37973), (-99167, 12879)), 13: ((4656, 99892), (-28581, 95829), (-54494, 83847), (-74415, 66801), (-88621, 46328), (-97153, 23691), (-100000, 0))}, {
        60: (34125, -129988, -305053, -434618, -568792, -668396, -768143, -839172, -905654, -948315, -981602),
        61: (10571, -150711, -322664, -449339, -580435, -677508, -774662, -843745, -908373, -949814, -982139),
    }, {
        60: (61, Q('3/109'), (44151, 73175)),
    }),
    (Q('1/40'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((240, 100000), (-34151, 93988), (-60370, 79721), (-79801, 60265), (-92727, 37440), (-99192, 12690)), 13: ((2283, 99974), (-30382, 95273), (-55717, 83040), (-75131, 65994), (-88948, 45697), (-97236, 23348), (-100000, 0))}, {
        61: (10571, -150711, -322664, -449339, -580435, -677508, -774662, -843745, -908373, -949814, -982139),
        62: (-13126, -171500, -340263, -464025, -592025, -686567, -781134, -848281, -911069, -951300, -982671),
        63: (-36965, -192353, -357850, -478675, -603561, -695576, -787561, -852783, -913742, -952772, -983198),
        64: (-62891, -213270, -375424, -493291, -615045, -704534, -793943, -857251, -916393, -954232, -983720),
        65: (-86971, -234248, -392985, -507873, -626478, -713444, -800282, -861685, -919021, -955678, -984238),
        66: (-111183, -255287, -410534, -522422, -637861, -722305, -806578, -866087, -921628, -957113, -984751),
    }, {
        61: (62, Q('1/40'), (5256, 10856)),
        62: (63, Q('1/40'), (5921, 10770)),
        63: (64, Q('1/40'), (6896, 10182)),
        64: (65, Q('1/40'), (9328, 10854)),
        65: (66, Q('1/40'), (10281, 13455)),
    }),
    (Q('1/50'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((-11682, 99315), (-42732, 90410), (-65819, 75285), (-82674, 56258), (-93782, 34711), (-99310, 11726)), 13: ((-9803, 99518), (-39375, 91922), (-61753, 78655), (-78641, 61770), (-90543, 42449), (-97640, 21597), (-100000, 0))}, {
        66: (-111183, -255287, -410534, -522422, -637861, -722305, -806578, -866087, -921628, -957113, -984751),
        67: (-135525, -276385, -428071, -536937, -649194, -731118, -812832, -870456, -924214, -958535, -985259),
        68: (-161800, -297541, -445595, -551421, -660479, -739886, -819045, -874794, -926780, -959946, -985763),
        69: (-186357, -318754, -463107, -565872, -671716, -748607, -825218, -879101, -929325, -961345, -986263),
        70: (-211036, -340024, -480606, -580291, -682906, -757284, -831352, -883378, -931850, -962733, -986759),
    }, {
        66: (67, Q('1/50'), (10898, 13015)),
        67: (68, Q('1/50'), (11305, 15068)),
        68: (69, Q('1/50'), (13512, 20505)),
        69: (70, Q('1/50'), (18286, 26431)),
    }),
    (Q('1/50'), 'cubic-scaled-precise-coupled-ordinary-loss-compress-skew', {12: ((-23920, 97097), (-51269, 85857), (-71143, 70276), (-85449, 51947), (-94795, 31841), (-99423, 10724)), 13: ((-22237, 97496), (-48350, 87535), (-67667, 73629), (-82040, 57179), (-92076, 39013), (-98027, 19768), (-100000, 0))}, {
        71: (-237519, -361349, -498093, -594680, -694050, -765917, -837446, -887625, -934356, -964109, -987250),
        72: (-259113, -382727, -515568, -609037, -705148, -774506, -843503, -891843, -936843, -965475, -987737),
    }, {
        71: (72, Q('1/50'), (30813, 44944)),
    }),
)


def verify(library, selected):
    total = 0
    prior = Q(0)
    for exponent, (center, radius) in enumerate(SKEW_CENTERS, 1):
        q = 2*exponent
        require(radius >= 1-center, ("positive skew extremum", q))
        require((q*center)**q*center <= (q+1)**(q+1)*radius**q,
                ("negative skew extremum", q))
    by_mode = {
        "plain": replay_python if library is None else normalized_native(library),
        "scaled": replay_python if library is None else scaled_native(library),
        "compress": replay_python if library is None else scaled_native(library, compress=True),
    }
    for declared, mode, anchors, rays, expected in BANDS:
        knots = tuple(sorted(rays))
        lower, upper = Q(knots[0], 100), Q(knots[-1], 100)
        require(prior <= lower < upper, ("ordered bands", prior, lower, upper))
        prior = upper
        if selected and lower not in selected:
            continue
        small = mode == "small"
        verify_setup(rays, anchors, declared, small)
        strategy = "compress" if "compress" in mode else "scaled" if "scaled" in mode else "plain"
        for left, right in zip(knots, knots[1:]):
            cap = declared if small else exact_cap(declared, Q(left, 100))
            target, recorded_cap, receipt = expected[left]
            require(right == target and cap == recorded_cap, ("knot or cap mismatch", left, right, cap))
            cells = tuple(
                verify_small_band(j, cap, rays, left, right) if small else
                verify_fused_band(j, cap, rays, left, right, by_mode[strategy], mode.startswith("cubic"),
                                  precise="precise" in mode, monotone="monotone" in mode,
                                  coupled="coupled" in mode, ordinary="ordinary" in mode,
                                  loss_transport="loss" in mode, center_ratio="centered" in mode,
                                  center_skew="skew" in mode, balanced="balanced" in mode)
                for j in ORDERS)
            if library is not None:
                require(cells == receipt, ("radial cell counts", left, right, cells, receipt))
            total += sum(cells)
            print(left, right, cap, cells, flush=True)
    print(f"Verified {total} exact radial cells.", flush=True)


def main():
    selected = {Q(item) for item in argv[1:] if item != "--native"}
    require(selected <= {Q(min(rays), 100) for _, _, _, rays, _ in BANDS},
            "select band lower endpoints")
    if "--native" not in argv:
        verify(None, selected)
        return
    with TemporaryDirectory(prefix="sfc-thirteen-light-") as directory:
        library = Path(directory)/"bernstein.so"
        source = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source),
             "-o", str(library)), check=True)
        verify(CDLL(str(library)), selected)


if __name__ == "__main__":
    main()
