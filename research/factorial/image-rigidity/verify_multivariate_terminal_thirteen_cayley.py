"""Exact grouped Cayley separation for light shared rows at terminal thirteen."""

from fractions import Fraction as Q
from sys import argv

from verify_multivariate_terminal_twelve_cayley import (
    GRID, I, ONE, ZERO, Limits, accurate_projection, bisect, height, interval,
    log_bounds, model, model_eval, rational, require, scalar, strength_factor,
    transform_family,
)
from verify_multivariate_terminal_twelve_cayley_lower import skew_model

ORDERS = (12, 13)
LOWER, UPPER = ORDERS
REFERENCE_RATIO = Q(144, 143)
SCALE = 100000

MODES = {
    "flat": (None, 0, 50, Q(13,5), False),
    "upper": (Q(1,70), 80, 120, Q(13,5), False),
    "seam": (Q(1,90), 50, 80, Q(13,5), False),
    "middle": (Q(1,60), 120, 145, Q(13,5), False),
    "later": (Q(1,55), 145, 240, Q(13,5), True),
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

def check_ray(n,i,left,right,limits,rays,linearized=False,window=None,conditioned=False):
    ratio = interval(*limits.ratio)
    lo,hi = window or (left,right)
    lam = interval(Q(lo,1000),Q(hi,1000)).sqrt()
    a = moment_ratio(n,lo,hi,limits)
    radius = strength_factor(limits.qmax)*a/2
    skew = [(radius*Q(k-4,4),radius*Q(k-3,4)) for k in range(8)]
    pending = [(lam,interval(Q(0),limits.disk),ratio,None,0)]
    count = 0
    while pending:
        lam,x,r,sk,depth = pending.pop()
        root_ratio = r.sqrt() if n==UPPER else ONE
        scale = lam*root_ratio
        p,q = lam.square(),scale.square()
        b = ray_for(p,i,left,right,rays)
        if n==UPPER:
            b = b/root_ratio
        if sk is None:
            base,groups,box,ok = model(n,scale.lo,scale.hi,a,limits,linearized)
        else:
            base,groups,box,ok = skew_model(n,scale.lo,scale.hi,a,limits,linearized,*sk)
        if not ok:
            require(depth<32,("scaled disk",n,i,left,right,lam.shown()))
            pending.extend((half,x,r,sk,depth+1) for half in bisect(lam))
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
        if conditioned and sk is None and depth>=11 and x.hi<=GRID//4:
            pending.extend((lam,x,r,skew_bin,depth) for skew_bin in skew)
            continue
        require(depth<44,("unresolved",n,i,left,right,lam.shown(),x.shown(),r.shown(),sk,count))
        wl,wx,wr = 10*(lam.hi-lam.lo),x.hi-x.lo,12*(r.hi-r.lo) if n==UPPER else 0
        if wr>max(wl,wx):
            pending.extend((lam,x,half,sk,depth+1) for half in bisect(r))
        elif wl>wx:
            pending.extend((half,x,r,sk,depth+1) for half in bisect(lam))
        else:
            pending.extend((lam,half,r,sk,depth+1) for half in bisect(x))
    return count

def validate_rays(limits,rays,left,right):
    if not limits.clip:
        require(limits.disk**2*limits.qmax<1,"Cayley disk left half-plane")
    for end in (left,right):
        slot = rays[end]
        levels = (int(limits.disk*SCALE),*slot,0)
        require(len(slot)==UPPER-2 and all(a>b for a,b in zip(levels,levels[1:])), "strict ray order")
        require(Q(max(slot),SCALE)**2*limits.qmax < 1, "left half-plane rays")


RAYS = {
    "flat": {
        0: (118044, 105664, 92562, 82034, 70861, 61206, 50927, 41757, 31960, 23062, 13511),
        5: (118414, 105954, 92769, 82196, 70977, 61295, 50990, 41804, 31991, 23083, 13523),
        10: (118787, 106246, 92977, 82359, 71094, 61385, 51053, 41851, 32022, 23104, 13534),
        15: (119163, 106541, 93186, 82522, 71211, 61475, 51116, 41898, 32054, 23126, 13545),
        20: (119543, 106837, 93397, 82687, 71330, 61566, 51180, 41945, 32085, 23147, 13557),
        25: (119926, 107137, 93609, 82853, 71448, 61657, 51244, 41993, 32117, 23168, 13568),
        30: (120312, 107439, 93823, 83020, 71568, 61749, 51308, 42040, 32149, 23189, 13580),
        35: (120632, 107743, 94038, 83188, 71688, 61841, 51373, 42088, 32181, 23211, 13591),
        40: (120952, 108050, 94255, 83357, 71809, 61934, 51438, 42137, 32213, 23233, 13603),
        45: (121274, 108359, 94473, 83528, 71930, 62027, 51503, 42185, 32245, 23254, 13615),
        50: (121597, 108671, 94693, 83699, 72052, 62121, 51569, 42234, 32277, 23276, 13626),
    },
    "seam": {
        50: (121597, 108671, 94693, 83699, 72052, 62121, 51569, 42234, 32277, 23276, 13626),
        55: (121922, 108986, 94914, 83872, 72175, 62215, 51635, 42283, 32310, 23298, 13638),
        60: (122366, 109302, 95137, 84045, 72298, 62309, 51701, 42331, 32342, 23320, 13650),
        65: (122835, 109622, 95362, 84219, 72423, 62404, 51768, 42380, 32375, 23342, 13662),
        70: (123301, 109943, 95588, 84394, 72548, 62499, 51835, 42430, 32408, 23364, 13674),
        75: (123764, 110268, 95815, 84570, 72673, 62595, 51902, 42479, 32441, 23386, 13686),
        80: (124226, 110595, 96045, 84748, 72800, 62691, 51970, 42529, 32475, 23408, 13698),
    },
    "upper": {
        80: (123924, 110598, 96045, 84750, 72800, 62693, 51970, 42530, 32475, 23409, 13698),
        85: (124348, 110929, 96276, 84929, 72927, 62790, 52038, 42580, 32508, 23431, 13710),
        90: (124826, 111262, 96508, 85109, 73055, 62887, 52106, 42630, 32542, 23453, 13722),
        95: (125304, 111598, 96743, 85290, 73184, 62985, 52175, 42681, 32575, 23476, 13734),
        100: (125783, 111937, 96979, 85472, 73313, 63083, 52244, 42731, 32609, 23499, 13747),
        105: (126263, 112279, 97217, 85656, 73443, 63183, 52313, 42782, 32643, 23521, 13759),
        110: (126744, 112624, 97456, 85841, 73574, 63282, 52383, 42834, 32678, 23544, 13771),
        115: (127227, 112972, 97697, 86027, 73706, 63382, 52453, 42885, 32712, 23567, 13784),
        120: (127713, 113324, 97941, 86214, 73839, 63483, 52524, 42937, 32747, 23590, 13796),
    },
    "middle": {
        120: (127488, 113326, 97941, 86216, 73839, 63484, 52524, 42938, 32747, 23591, 13796),
        125: (127983, 113681, 98186, 86405, 73972, 63586, 52595, 42990, 32781, 23614, 13809),
        130: (128481, 114039, 98433, 86595, 74106, 63687, 52666, 43042, 32816, 23637, 13822),
        135: (128982, 114400, 98681, 86786, 74241, 63790, 52737, 43095, 32851, 23660, 13834),
        140: (129485, 114764, 98932, 86979, 74377, 63893, 52809, 43147, 32887, 23684, 13847),
        145: (129992, 115132, 99185, 87173, 74514, 63997, 52882, 43200, 32922, 23708, 13860),
    },
    "later": {
        145: (129857, 115134, 99185, 87174, 74514, 63998, 52882, 43201, 32922, 23708, 13860),
        150: (130370, 115505, 99439, 87370, 74652, 64102, 52954, 43254, 32958, 23732, 13873),
        155: (130887, 115880, 99696, 87567, 74790, 64207, 53028, 43308, 32993, 23756, 13885),
        160: (131408, 116259, 99954, 87765, 74929, 64312, 53101, 43362, 33029, 23779, 13898),
        165: (131933, 116641, 100215, 87965, 75069, 64419, 53175, 43416, 33065, 23804, 13911),
        170: (132462, 117027, 100477, 88167, 75211, 64525, 53249, 43471, 33102, 23828, 13925),
        175: (132996, 117416, 100742, 88369, 75352, 64633, 53324, 43525, 33138, 23852, 13938),
        180: (133534, 117809, 101009, 88574, 75495, 64741, 53399, 43581, 33175, 23877, 13951),
        185: (134077, 118207, 101277, 88780, 75639, 64850, 53475, 43636, 33212, 23901, 13964),
        190: (134624, 118608, 101549, 88987, 75784, 64959, 53551, 43692, 33249, 23926, 13978),
        195: (135177, 119013, 101822, 89196, 75930, 65069, 53627, 43748, 33286, 23951, 13991),
        200: (135735, 119422, 102097, 89406, 76076, 65180, 53704, 43804, 33323, 23976, 14004),
        205: (136298, 119835, 102375, 89618, 76224, 65292, 53781, 43861, 33361, 24001, 14018),
        210: (136867, 120253, 102655, 89832, 76372, 65404, 53859, 43918, 33399, 24026, 14032),
        215: (137441, 120674, 102937, 90047, 76522, 65517, 53937, 43975, 33437, 24051, 14045),
        220: (138022, 121100, 103222, 90264, 76672, 65630, 54016, 44032, 33475, 24077, 14059),
        225: (138608, 121531, 103509, 90483, 76824, 65745, 54095, 44090, 33514, 24103, 14073),
        230: (139200, 121965, 103798, 90703, 76976, 65860, 54174, 44149, 33552, 24128, 14087),
        235: (139799, 122405, 104090, 90925, 77130, 65976, 54254, 44207, 33591, 24154, 14101),
        240: (140404, 122849, 104385, 91149, 77284, 66092, 54335, 44266, 33630, 24180, 14115),
    },
}


EXPECTED_CELLS = {
    "flat": {
        (12, 0): (1044, 1099, 2067, 1305, 2035, 1097, 957, 476, 397, 229, 161),
        (13, 0): (1653, 1355, 3444, 2908, 2951, 1402, 1942, 1009, 615, 363, 196),
        (12, 5): (799, 278, 642, 359, 561, 291, 281, 120, 102, 82, 46),
        (13, 5): (1535, 372, 1191, 821, 1004, 600, 597, 282, 212, 105, 63),
        (12, 10): (1133, 282, 730, 367, 566, 291, 319, 125, 111, 98, 51),
        (13, 10): (1855, 396, 1254, 846, 1138, 646, 613, 282, 223, 109, 74),
        (12, 15): (796, 291, 820, 370, 579, 297, 401, 137, 116, 99, 52),
        (13, 15): (1541, 432, 1335, 860, 1230, 675, 625, 295, 234, 109, 74),
        (12, 20): (1895, 299, 3469, 373, 598, 307, 432, 159, 116, 99, 53),
        (13, 20): (3721, 436, 6923, 874, 1322, 711, 631, 298, 234, 115, 87),
        (12, 25): (2227, 158, 1651, 201, 366, 163, 248, 101, 78, 57, 30),
        (13, 25): (5741, 230, 3340, 464, 787, 402, 357, 157, 126, 63, 52),
        (12, 30): (2538, 158, 1728, 218, 377, 169, 258, 103, 88, 58, 30),
        (13, 30): (8710, 230, 3313, 472, 816, 402, 360, 160, 132, 63, 52),
        (12, 35): (3011, 162, 1754, 224, 400, 170, 272, 110, 101, 58, 31),
        (13, 35): (10942, 237, 3364, 487, 849, 407, 362, 164, 137, 65, 54),
        (12, 40): (3380, 166, 1810, 230, 458, 172, 276, 119, 106, 58, 33),
        (13, 40): (11058, 242, 3469, 494, 864, 424, 372, 166, 143, 71, 62),
        (12, 45): (4475, 174, 1863, 240, 491, 174, 275, 119, 112, 60, 34),
        (13, 45): (12236, 249, 3470, 498, 893, 434, 408, 173, 195, 73, 77),
    },
    "upper": {
        (12, 80): (8266, 177, 2029, 250, 585, 175, 320, 130, 128, 67, 35),
        (13, 80): (21959, 243, 5520, 506, 2109, 451, 1204, 184, 1475, 77, 195),
        (12, 85): (5984, 177, 1887, 230, 514, 174, 300, 122, 120, 67, 34),
        (13, 85): (13152, 240, 4729, 504, 1478, 435, 928, 183, 336, 77, 112),
        (12, 90): (4411, 169, 1812, 208, 470, 171, 276, 117, 120, 66, 34),
        (13, 90): (12253, 233, 4411, 496, 1391, 435, 628, 174, 257, 77, 97),
        (12, 95): (4003, 94, 1058, 145, 312, 96, 161, 69, 67, 38, 22),
        (13, 95): (11227, 135, 2448, 283, 978, 230, 290, 97, 124, 45, 52),
        (12, 100): (3515, 94, 1009, 142, 290, 95, 146, 63, 62, 36, 22),
        (13, 100): (10285, 126, 2281, 267, 851, 226, 232, 94, 113, 43, 46),
        (12, 105): (3206, 92, 995, 141, 270, 95, 146, 63, 60, 36, 20),
        (13, 105): (8677, 126, 2158, 266, 729, 225, 217, 94, 93, 41, 41),
        (12, 110): (2940, 92, 935, 141, 254, 95, 146, 59, 59, 37, 20),
        (13, 110): (7283, 126, 1893, 266, 725, 214, 192, 91, 82, 41, 35),
        (12, 115): (2825, 88, 908, 133, 230, 89, 143, 57, 58, 36, 20),
        (13, 115): (6164, 128, 1785, 262, 691, 213, 192, 91, 79, 41, 32),
    },
    "seam": {
        (12, 50): (5664, 176, 1906, 248, 537, 174, 292, 123, 120, 66, 34),
        (13, 50): (15497, 253, 5091, 504, 1418, 435, 547, 181, 240, 76, 97),
        (12, 55): (5786, 177, 2057, 250, 580, 174, 303, 127, 123, 67, 34),
        (13, 55): (11648, 265, 5312, 504, 1492, 446, 912, 186, 305, 78, 112),
        (12, 60): (3963, 177, 1891, 242, 513, 174, 289, 121, 119, 67, 34),
        (13, 60): (7332, 249, 4793, 503, 1411, 436, 577, 180, 247, 77, 97),
        (12, 65): (3159, 169, 1835, 223, 465, 171, 275, 118, 110, 60, 34),
        (13, 65): (5891, 233, 4418, 496, 1303, 424, 399, 169, 206, 73, 83),
        (12, 70): (2744, 165, 1661, 194, 377, 171, 264, 110, 100, 58, 33),
        (13, 70): (5199, 229, 4218, 472, 1241, 407, 362, 165, 147, 69, 72),
        (12, 75): (2590, 163, 1581, 194, 361, 164, 260, 98, 99, 58, 30),
        (13, 75): (4478, 229, 3877, 472, 1155, 403, 361, 165, 138, 66, 54),
    },
    "middle": {
        (12, 120): (4503, 95, 1010, 141, 299, 95, 151, 64, 67, 38, 22),
        (13, 120): (13067, 133, 2285, 269, 987, 229, 311, 94, 127, 45, 54),
        (12, 125): (4227, 94, 977, 142, 278, 95, 146, 64, 62, 36, 22),
        (13, 125): (11941, 130, 2012, 268, 882, 225, 235, 94, 120, 42, 48),
        (12, 130): (3813, 85, 966, 142, 261, 95, 146, 64, 60, 35, 20),
        (13, 130): (11458, 130, 1876, 267, 760, 217, 220, 92, 107, 42, 44),
        (12, 135): (3482, 85, 933, 134, 248, 89, 146, 58, 59, 35, 20),
        (13, 135): (10469, 130, 1827, 266, 726, 213, 197, 91, 90, 41, 41),
        (12, 140): (3282, 83, 911, 134, 212, 89, 144, 58, 58, 35, 20),
        (13, 140): (9236, 131, 1766, 265, 692, 214, 192, 90, 80, 41, 35),
    },
    "later": {
        (12, 145): (3674, 84, 956, 139, 259, 94, 147, 62, 60, 35, 20),
        (13, 145): (4682, 130, 1950, 268, 899, 219, 266, 92, 115, 42, 48),
        (12, 150): (3493, 85, 939, 134, 248, 89, 147, 58, 60, 35, 20),
        (13, 150): (4905, 131, 1825, 268, 781, 216, 221, 90, 99, 41, 41),
        (12, 155): (3477, 82, 891, 132, 212, 89, 145, 54, 58, 35, 20),
        (13, 155): (5928, 130, 1811, 270, 741, 214, 192, 90, 85, 41, 40),
        (12, 160): (3466, 81, 869, 126, 211, 89, 131, 48, 49, 35, 20),
        (13, 160): (8582, 131, 1613, 270, 723, 202, 190, 86, 77, 39, 33),
        (12, 165): (3379, 73, 818, 110, 209, 85, 131, 46, 40, 36, 20),
        (13, 165): (8510, 129, 1549, 270, 714, 202, 190, 86, 77, 38, 32),
        (12, 170): (3346, 72, 774, 110, 204, 83, 131, 44, 38, 37, 20),
        (13, 170): (7318, 132, 1402, 268, 706, 186, 189, 87, 77, 38, 32),
        (12, 175): (3304, 64, 725, 110, 195, 83, 131, 44, 38, 37, 20),
        (13, 175): (6315, 125, 2113, 413, 1120, 234, 284, 85, 77, 36, 32),
        (12, 180): (3288, 63, 725, 110, 192, 77, 121, 43, 38, 35, 20),
        (13, 180): (7689, 122, 2038, 406, 1096, 166, 272, 81, 74, 35, 31),
        (12, 185): (3264, 59, 717, 110, 183, 77, 109, 40, 38, 34, 20),
        (13, 185): (7519, 112, 1889, 406, 1043, 142, 270, 81, 73, 35, 29),
        (12, 190): (3236, 60, 703, 110, 174, 77, 93, 40, 38, 32, 20),
        (13, 190): (7400, 113, 1812, 406, 1009, 142, 267, 81, 70, 35, 29),
        (12, 195): (3078, 60, 688, 109, 162, 77, 93, 40, 38, 31, 20),
        (13, 195): (7306, 101, 1707, 391, 975, 135, 251, 77, 68, 35, 29),
        (12, 200): (2900, 60, 600, 109, 160, 73, 84, 39, 38, 31, 20),
        (13, 200): (7186, 100, 1593, 372, 888, 135, 244, 75, 64, 31, 29),
        (12, 205): (2720, 60, 576, 110, 160, 71, 84, 39, 38, 31, 20),
        (13, 205): (6998, 93, 1505, 371, 866, 135, 228, 71, 62, 29, 29),
        (12, 210): (2523, 60, 563, 110, 156, 72, 76, 39, 37, 30, 20),
        (13, 210): (6882, 93, 1464, 339, 743, 127, 220, 64, 62, 27, 28),
        (12, 215): (2343, 60, 544, 104, 156, 73, 77, 39, 37, 30, 20),
        (13, 215): (6736, 89, 1428, 319, 594, 119, 208, 56, 61, 26, 26),
        (12, 220): (2183, 58, 475, 105, 148, 69, 73, 39, 37, 30, 19),
        (13, 220): (6554, 87, 1365, 319, 568, 112, 188, 54, 61, 25, 26),
        (12, 225): (2006, 59, 456, 105, 148, 68, 71, 39, 37, 28, 19),
        (13, 225): (6396, 85, 1327, 287, 528, 112, 185, 50, 61, 25, 25),
        (12, 230): (1832, 55, 421, 107, 149, 60, 71, 39, 35, 24, 19),
        (13, 230): (6189, 80, 1298, 223, 480, 112, 177, 48, 57, 25, 25),
        (12, 235): (1620, 55, 400, 108, 149, 58, 71, 39, 34, 24, 19),
        (13, 235): (5998, 76, 1229, 218, 480, 112, 143, 48, 53, 25, 25),
    },
}


def run(mode, order, ray, lo, hi):
    cap, first, last, disk, clip = MODES[mode]
    rays = RAYS[mode]
    nodes = sorted(rays)
    require(nodes == list(range(first,last+1,5)) and lo in nodes and hi in nodes and lo<hi, "range")
    require(order in (0,*ORDERS) and -2<=ray<UPPER-2, "ray selection")
    anchor_for(limits_for(cap,first,first,disk,clip), rays, first)
    total = 0
    for left,right in zip(nodes,nodes[1:]):
        if not lo<=left<right<=hi:
            continue
        limits = limits_for(cap,left,right,disk,clip)
        validate_rays(limits,rays,left,right)
        for n in ORDERS:
            if order and order != n:
                continue
            for index in range(UPPER-2):
                if (ray>=0 and ray!=index) or (ray==-2 and index==0):
                    continue
                linearized = index in (0,1,2)
                conditioned = index in (0,2) and left>=20
                count = 0
                if index==0:
                    for sub in range(left,right):
                        narrow = limits_for(cap,sub,sub+1,disk,clip,refined=True)
                        validate_rays(narrow,rays,left,right)
                        count += check_ray(n,index,left,right,narrow,rays,linearized,
                                           (sub,sub+1),conditioned)
                else:
                    count = check_ray(n,index,left,right,limits,rays,linearized,
                                      conditioned=conditioned)
                require(count == EXPECTED_CELLS[mode][n,left][index],
                        ("replay cell count", mode, n, left, index, count))
                print(mode,n,left,right,index,count,flush=True)
                total += count
        model.cache_clear()
        skew_model.cache_clear()
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
