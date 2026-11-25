#!/usr/bin/env python3
"""Exact terminal-ten light-row separation on bounded outer-mean bands."""

from ctypes import CDLL
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_nine_high_curvature import anchor_value
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

ORDERS=(9,10)
REFERENCE_RATIO=Q(81,80)
FLOOR=Q(1,40)

@lru_cache(maxsize=None)
def shape(order, skew_terms, mean_terms):
    m=order-1
    distances=tuple(min(k,order-k) for k in range(order+1))
    gaussian=tuple(rational(Q(k*(order-k),m)) for k in range(order+1))
    weights=tuple(skew_weights(order,k) for k in distances)
    points=tuple(Q(1)-Q(2*h,m) for h in range(order))
    skew=[]
    for l in range(1,skew_terms+1):
        moment=tuple(sum(w*p**(2*l) for w,p in zip(row,points,strict=True))
                     for row in weights)
        skew.append(tuple(rational(-a if 2*k<=order else a) for k,a in enumerate(moment)))
    defect=tuple(tuple(rational(sum(p-p**(2*l+1) for p in points[:r]))
                       for r in distances) for l in range(1,mean_terms+1))
    sizes=tuple(rational(sum(abs(w) for w in row)) for row in weights)
    return gaussian,tuple(skew),defect,sizes,distances

@lru_cache(maxsize=None)
def ratio_bounds(cap):
    lo,hi=log_bounds(1+8*cap)
    lower=log_bounds(1+Q(81,10)*cap)[0]/hi
    upper=1+log_bounds(1/(1-cap/10))[1]/lo
    require(lower<=REFERENCE_RATIO<=upper,("ratio",cap))
    return interval(lower,upper)

def nominal_exponent(order,k):
    return Q(k*(order-k),order-1)*(1 if order==9 else REFERENCE_RATIO)

@lru_cache(maxsize=None)
def independent_exponents(order,cap):
    d0,d1=log_bounds(1+8*cap)
    if order==9:
        x=cap
        ratio=(Q(1),Q(1))
    else:
        x=9*cap/(10-cap)
        ratio=(log_bounds(1+Q(81,10)*cap)[0]/d1,1+log_bounds(1/(1-cap/10))[1]/d0)
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
    threshold = 1+8*FLOOR/(1+FLOOR)
    cap = min(FLOOR, declared)
    gap = 1
    while threshold**(-gap) >= outer_lower**2:
        base = 1
        while 9*base+1 < gap/FLOOR:
            x = Q(gap, 9*base+1)
            if cap < x <= declared:
                row_square = prod(Q(9*base+h, 9*base+8*gap+h) for h in range(1, gap+1))
                if row_square >= outer_lower**2:
                    cap = x
            base += 1
        gap += 1
    return cap


@lru_cache(maxsize=None)
def constants(order, cap9):
    cap = cap9 if order == 9 else 9*cap9/(10-cap9)
    m = order-1
    T = rational(m*cap/(2+m*cap))
    lo, hi = log_bounds(1+m*cap)
    s = interval(lo/2, hi/2)
    require(T.square().hi < GRID//2, "moment quotient monotonicity range")
    return T, s, T.square()/(scalar(2)*s)


@lru_cache(maxsize=None)
def leading_skew_direction(order, cap9, skew_terms, mean_terms):
    _, skew, _, _, _ = shape(order, skew_terms, mean_terms)
    T, _, _ = constants(order, cap9)
    weights = tuple(T**(2*l-2)*scalar(1, 2*l) for l in range(2, skew_terms+1))
    direction = tuple(first+sum((w*v for w,v in zip(weights,rest,strict=True)),ZERO)
                      for first,*rest in zip(*skew,strict=True))
    return direction, weights


@lru_cache(maxsize=8192)
def coefficient_model(order, cap9, means, signed, extended=False):
    skew_terms,mean_terms=(4,3) if extended else (3,2)
    midpoint = sum(means)/2
    skew_midpoint = sum(signed)/2
    ref_ratio = Q(1) if order == 9 else REFERENCE_RATIO
    T, s, skew_bound = constants(order, cap9)
    Pmax = constants(9, cap9)[2]
    U9 = -log(interval(*means))
    leading = interval(*signed)*Pmax*U9
    nominal = rational(skew_midpoint)*Pmax*(-log(rational(midpoint)))
    h = interval(means[0]/midpoint-1, means[1]/midpoint-1)
    magnitude = I(0, max(-h.lo,h.hi))
    log_error = rational(ref_ratio)*magnitude.square()/(scalar(2)*(ONE-magnitude))
    outer_mean = I(-log_error.hi,0)
    strength = U9
    if order == 10:
        rho = ratio_bounds(cap9)
        outer_mean += (rational(ref_ratio)-rho)*U9
        strength = rho*U9
        eta = cap9/(20+79*cap9)
        transfer = rational(REFERENCE_RATIO**2*(1+eta**2))
        transfer_width = rational(2*REFERENCE_RATIO**2*eta)
        leading = transfer*leading+(transfer_width*Pmax*U9).symmetric_error()
        nominal = transfer*nominal
    gaussian, skew, mean, skew_sizes, distances = shape(order,skew_terms,mean_terms)
    direction, weights = leading_skew_direction(order,cap9,skew_terms,mean_terms)
    reference = []
    for k, v in enumerate(direction):
        if k == 0 or k == order:
            reference.append(ONE)
            continue
        exponent = ref_ratio*Q(k*(order-k), order-1)
        value = scalar(comb(order, k))*interval(*power_bounds(midpoint, exponent))
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



SCALE=10**6

SMALL = {
    0: (945375, 763739, 543893, 250000, -86824, -380717, -676915, -858550),
    1: (942757, 755448, 528967, 232577, -105363, -395569, -686150, -862963),
    2: (939162, 745826, 512711, 214794, -123109, -409283, -694382, -866800),
    3: (934637, 735273, 495784, 197042, -140124, -422170, -701963, -870293),
    4: (929223, 723978, 478463, 179438, -156500, -434394, -709050, -873533),
}

EARLY = {
    4: (929223, 723978, 478463, 179438, -156500, -434394, -709050, -873533),
    6: (915916, 699628, 443215, 144845, -187610, -457259, -722096, -879450),
    8: (899610, 673438, 407670, 111148, -216886, -478449, -733992, -884803),
    10: (880648, 645809, 372147, 78333, -244644, -498308, -745003, -889731),
    12: (859325, 617011, 336800, 46349, -271115, -517072, -755302, -894319),
    14: (835890, 587235, 301707, 15134, -296475, -534910, -765009, -898627),
    16: (810558, 556625, 266905, -15375, -320864, -551951, -774211, -902699),
    18: (783508, 525290, 232410, -45233, -344393, -568295, -782978, -906567),
    20: (754895, 493317, 198226, -74491, -367152, -584021, -791363, -910257),
    22: (724852, 460778, 164349, -103195, -389217, -599197, -799409, -913791),
    24: (693493, 427731, 130772, -131385, -410654, -613875, -807151, -917184),
    26: (657392, 394223, 97487, -159097, -431515, -628102, -814619, -920450),
    28: (623609, 360298, 64482, -186363, -451847, -641917, -821838, -923603),
    30: (585113, 325989, 31748, -213213, -471692, -655353, -828830, -926651),
    32: (549276, 291328, -725, -239672, -491086, -668440, -835612, -929603),
    34: (508774, 256340, -32949, -265762, -510059, -681204, -842201, -932467),
    36: (471180, 221051, -64932, -291506, -528639, -693666, -848612, -935250),
    38: (428994, 185479, -96684, -316922, -546851, -705847, -863289, -937956),
}

EARLY_ANCHOR = {9: ((92154, 38829), (43051, 90259), (-24137, 97043), (-79068, 61223), (-100000, 0)), 10: ((93691, 34957), (52642, 85023), (-7163, 99743), (-62742, 77868), (-95639, 29210))}

BANDS = (
    (Q(1,8), {
        38: (436625, 185479, -96684, -316922, -546851, -705847, -854856, -937956),
        40: (397525, 149644, -128214, -342027, -564719, -717766, -860945, -940592),
        42: (357720, 113562, -159530, -366838, -582261, -729437, -866888, -943162),
        44: (313454, 77248, -190640, -391369, -599496, -740876, -872695, -945670),
        46: (272375, 40715, -221552, -415633, -616441, -752095, -878373, -948120),
        48: (226977, 3976, -252273, -439642, -633110, -763107, -883930, -950515),
        50: (184785, -32959, -282809, -463409, -649518, -773922, -889372, -952858),
        52: (139576, -70079, -313167, -486942, -665677, -784550, -894706, -955152),
    }, {9: ((42518, 90511), (-13915, 99027), (-59901, 80074), (-89702, 44198), (-100000, 0)), 10: ((46334, 88618), (-5422, 99853), (-49470, 86907), (-81269, 58270), (-97889, 20439))}),
    (Q(99,1000), {
        52: (156528, -70079, -313167, -486942, -665677, -784550, -894706, -955152),
        54: (109477, -107375, -343353, -510252, -681599, -795000, -899937, -957399),
        56: (65529, -144839, -373372, -533349, -697294, -805282, -905069, -959603),
        58: (21089, -182462, -403230, -556239, -712773, -815402, -910109, -961764),
        60: (-27119, -220238, -432931, -578932, -728045, -825367, -915060, -963885),
        62: (-72386, -258160, -462481, -601435, -743118, -835185, -919925, -965968),
    }, {9: ((13845, 99037), (-34773, 93759), (-70520, 70901), (-92551, 37872), (-100000, 0)), 10: ((17461, 98464), (-27860, 96041), (-62615, 77970), (-86390, 50367), (-98479, 17375))}),
    (Q(3,40), {
        62: (-72386, -258160, -462481, -601435, -743118, -835185, -919925, -965968),
        64: (-121168, -296221, -491884, -623754, -758000, -844862, -924710, -968014),
        Q(641, 10): (-123171, -298131, -493347, -624861, -758735, -845339, -924945, -968115),
        66: (-161228, -334415, -521143, -645896, -772700, -854403, -929416, -970026),
        68: (-207864, -372738, -550265, -667867, -787223, -863813, -934049, -972004),
        70: (-254892, -411184, -579251, -689673, -801577, -873099, -938609, -973950),
    }, {9: ((-8197, 99664), (-49061, 87138), (-77397, 63322), (-94345, 33150), (-100000, 0)), 10: ((-5004, 99875), (-43435, 90074), (-71226, 70191), (-89640, 44325), (-98848, 15134))}),
    (Q(11,200), {
        70: (-254892, -411184, -579251, -689673, -801577, -873099, -938609, -973950),
        Q(701, 10): (-257391, -413112, -580694, -690755, -802286, -873557, -938834, -974046),
        72: (-304875, -449749, -608106, -711319, -815767, -882264, -943101, -975865),
        74: (-352506, -488427, -636833, -732810, -829800, -891312, -947527, -977751),
        76: (-402764, -527215, -665436, -754151, -843680, -900249, -951889, -979608),
    }, {9: ((-26574, 96404), (-60187, 79859), (-82568, 56414), (-95670, 29107), (-100000, 0)), 10: ((-23862, 97111), (-55663, 83076), (-77747, 62892), (-92052, 39071), (-99120, 13240))}),
)

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
    sectors=range(order%2,9,2)
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


def verify_setup(rays, anchors, declared, small):
    lower,upper=Q(min(rays),100),Q(max(rays),100)
    for values in rays.values():
        bounds=(SCALE,*values,-SCALE)
        require(len(values)==8 and all(a>b for a,b in zip(bounds,bounds[1:])),"ordered rays")
    annulus=(SMALL_INNER,SMALL_OUTER) if small else ANNULUS
    for order in ORDERS:
        box=independent_coefficients(order,declared,(Q(0),upper))
        require(sum(cap*annulus[0]**k for k,(_,cap) in enumerate(box) if k)<1,("inner radius",order))
        require(sum(cap*annulus[1]**k for k,(_,cap) in enumerate(box[:-1]))<annulus[1]**order,
                ("outer radius",order))
        cosine=tuple(Q(x,SCALE) for x in rays[min(rays)])
        if small:
            for i,x in enumerate((Q(1),*cosine,Q(-1))):
                expected=(i+(order==10))//2
                require((-1)**expected*anchor_value(order,x)>0,("binomial anchor",order,i))
        else:
            gaussian_anchor(order,lower,anchors[order],cosine)


def cosines(rays,left,right,index,means):
    a,b=rays[left][index],rays[right][index]
    return tuple(sorted(Q(a,SCALE)+Q(b-a,SCALE)*(100*v-left)/(right-left) for v in means))


def verify_small_band(order,cap,rays,left,right):
    pending=[((Q(left,100),Q(right,100)),tuple(range(8)),0)]
    total=0
    while pending:
        means,indices,depth=pending.pop()
        box=independent_coefficients(order,cap,means)
        bad=[]
        for index in indices:
            used=radial_replay(order,box,cosines(rays,left,right,index,means))
            if used is None:
                bad.append(index)
            else:
                total+=used
        if bad:
            require(depth<16,("small band",order,cap,means,bad))
            mid=sum(means)/2
            pending.extend((((means[0],mid),bad,depth+1),((mid,means[1]),bad,depth+1)))
    return total


def verify_fused_band(order,cap,rays,left,right,replay,extended):
    pending=[((Q(left,100),Q(right,100)),(Q(0),Q(1)),tuple(range(8)),0,0)]
    total=0
    while pending:
        means,p,active,mean_depth,skew_depth=pending.pop()
        model=coefficient_model(order,cap,means,p,extended)
        if model is None:
            continue
        bad=[]
        for index in active:
            a,b=rays[left][index],rays[right][index]
            slope=Q(b-a,SCALE)*100*sum(means)/(2*(right-left))
            count=replay(order,cosines(rays,left,right,index,means),model,slope) if model else None
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



def normalized_native(library):
    fast=native_replay(library)

    def replay(order, cosines, model, slope):
        normalized=[]
        for index,(parameter,vector) in enumerate(model):
            largest=max(max(abs(c.lo),abs(c.hi)) for c in vector)
            # The first two vectors are fused. Their relative scale
            # must be preserved, so normalize only the other groups.
            shift=max(0,(largest-1).bit_length()-57) if index>=2 else 0
            if shift:
                parameter=I(parameter.lo<<shift,parameter.hi<<shift)
                vector=tuple(I(c.lo>>shift,-((-c.hi)>>shift)) for c in vector)
            normalized.append((parameter,vector))
        return fast(order,cosines,tuple(normalized),slope)

    return replay


def verify(replay, selected):
    definitions=((Q(1,8),SMALL,{},True,False),
                 (Q(1,8),EARLY,EARLY_ANCHOR,False,False),
                 *((cap,rays,anchors,False,True) for cap,rays,anchors in BANDS))
    for declared,rays,anchors,small,extended in definitions:
        if selected and min(rays) not in selected:
            continue
        verify_setup(rays,anchors,declared,small)
        knots=tuple(sorted(rays))
        total=0
        for left,right in zip(knots,knots[1:]):
            cap=declared if small else exact_cap(declared,Q(left,100))
            counts=tuple(verify_small_band(j,cap,rays,left,right) if small else
                         verify_fused_band(j,cap,rays,left,right,replay,extended)
                         for j in ORDERS)
            total+=sum(counts)
            print(left,right,cap,counts,flush=True)
        print(f"Verified {total} radial cells for [{min(rays)}/100,{max(rays)}/100].",flush=True)


def main():
    selected=tuple(map(int,(arg for arg in argv[1:] if arg!="--native")))
    require(set(selected)<={0,4,*(min(rays) for _,rays,_ in BANDS)},"invalid band selection")
    if "--native" not in argv:
        verify(replay_python,selected)
        return
    with TemporaryDirectory(prefix="sfc-ten-") as directory:
        library=Path(directory)/"bernstein.so"
        source=Path(__file__).with_suffix(".cpp")
        run(("c++","-std=c++17","-O3","-fPIC","-shared",str(source),"-o",str(library)),check=True)
        verify(normalized_native(CDLL(str(library))),selected)


if __name__=="__main__":
    main()
