#!/usr/bin/env python3
"""Exact terminal-ten angular separation for four marked shared unit rows."""

from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from sys import argv

from verify_multivariate_terminal_eight_late_bridge import marked_disk_check
from verify_multivariate_terminal_nine_private_linear import shape
from verify_multivariate_terminal_nine_upper_marked_rows import (
    marked_core, matrix, precise_power, sector_check, table_value,
)
from verify_multivariate_terminal_seven_upper_bridge import bisect_u, height
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, bisect, interval, log, rational, require, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds
from verify_multivariate_terminal_six_marked_row_tail import projected_interval
from verify_multivariate_terminal_ten_near_flat import accurate_projection


ORDERS = (9, 10)
REFERENCE_RATIO = Q(81, 80)
SCALE = 100000
COORDINATE_P = Q(1, 8)
P = rational(COORDINATE_P)
ROOT_P = P.sqrt()
DISK = Q(11, 4)
RADII = tuple(Q(d, 2000) for d in range(4, 81, 2))
OUTER_CUTOFFS = {2: Q(19, 25), 3: Q(21, 25), 4: Q(22, 25), 5: Q(22, 25)}
RESIDUAL_CUTOFFS = {2: Q(113, 125), 3: Q(119, 125), 4: Q(97, 100), 5: Q(119, 125)}
RESIDUAL_CAPS = {2: Q(1, 37), 3: Q(1, 82), 4: Q(1, 136), 5: Q(1, 82)}

RAYS = {
    2: {
        904: (186116, 157802, 127819, 106215, 83269, 64598, 44581, 27283),
        908: (184323, 156361, 126740, 105348, 82618, 64104, 44249, 27081),
        912: (182531, 154919, 125658, 104479, 81966, 63608, 43914, 26879),
        916: (180740, 153477, 124574, 103607, 81311, 63110, 43578, 26675),
        920: (178952, 152034, 123488, 102732, 80653, 62609, 43241, 26470),
        924: (177165, 150591, 122399, 101855, 79993, 62106, 42901, 26264),
        928: (175379, 149147, 121308, 100975, 79330, 61601, 42560, 26057),
        932: (173594, 147701, 120213, 100092, 78664, 61094, 42217, 25849),
        936: (171810, 146255, 119116, 99206, 77995, 60584, 41873, 25640),
        940: (170027, 144808, 118015, 98317, 77323, 60071, 41526, 25429),
        944: (168244, 143359, 116912, 97424, 76647, 59556, 41177, 25217),
        948: (166462, 141908, 115804, 96528, 75968, 59037, 40826, 25004),
        952: (164680, 140455, 114693, 95628, 75286, 58516, 40473, 24789),
        956: (162899, 139001, 113578, 94724, 74599, 57992, 40117, 24573),
        960: (161117, 137545, 112460, 93816, 73909, 57464, 39759, 24355),
        964: (159335, 136086, 111336, 92903, 73215, 56932, 39398, 24136),
        968: (157553, 134625, 110209, 91986, 72516, 56398, 39035, 23915),
        972: (155771, 133161, 109076, 91065, 71813, 55859, 38669, 23692),
        976: (153989, 131694, 107939, 90138, 71105, 55316, 38300, 23467),
        980: (152206, 130225, 106796, 89206, 70392, 54769, 37927, 23241),
        984: (150422, 128752, 105648, 88268, 69674, 54218, 37552, 23012),
        988: (148638, 127276, 104495, 87325, 68950, 53662, 37173, 22781),
        992: (146854, 125797, 103335, 86376, 68220, 53101, 36790, 22548),
        996: (145069, 124314, 102169, 85420, 67484, 52535, 36403, 22312),
        1000: (143284, 122828, 100996, 84457, 66742, 51964, 36013, 22073),
    },
    3: {
        952: (142863, 122261, 100285, 83790, 66131, 51466, 35642, 21846),
        956: (140929, 120663, 99038, 82770, 65348, 50866, 35233, 21596),
        960: (138989, 119057, 97782, 81743, 64560, 50260, 34819, 21344),
        964: (137043, 117444, 96518, 80708, 63764, 49648, 34402, 21090),
        968: (135089, 115823, 95245, 79665, 62961, 49031, 33980, 20833),
        972: (133129, 114194, 93963, 78614, 62151, 48407, 33554, 20573),
        976: (131161, 112556, 92671, 77553, 61333, 47778, 33123, 20310),
        980: (129185, 110909, 91369, 76483, 60507, 47141, 32688, 20044),
        984: (127201, 109252, 90056, 75403, 59672, 46497, 32247, 19775),
        988: (125209, 107585, 88732, 74313, 58827, 45846, 31801, 19503),
        992: (123207, 105908, 87396, 73212, 57973, 45186, 31348, 19227),
        996: (121197, 104219, 86047, 72098, 57108, 44518, 30890, 18946),
        1000: (119177, 102519, 84684, 70972, 56232, 43840, 30424, 18662),
    },
    4: {
        970: (120498, 103524, 85360, 71488, 56584, 44099, 30586, 18760),
        974: (118382, 101752, 83950, 70326, 55684, 43404, 30110, 18469),
        978: (116250, 99964, 82526, 69151, 54772, 42700, 29626, 18174),
        982: (114102, 98160, 81085, 67961, 53847, 41985, 29136, 17874),
        986: (111938, 96339, 79626, 66756, 52910, 41261, 28638, 17570),
        990: (109755, 94500, 78150, 65535, 51958, 40525, 28132, 17260),
        994: (107554, 92641, 76655, 64296, 50992, 39776, 27617, 16946),
        998: (105333, 90763, 75139, 63039, 50009, 39015, 27093, 16625),
        1000: (104215, 89815, 74372, 62403, 49512, 38629, 26827, 16462),
    },
    5: {
        952: (121098, 103943, 85594, 71647, 56670, 44153, 30612, 18774),
        956: (118939, 102137, 84160, 70465, 55754, 43447, 30128, 18478),
        960: (116764, 100314, 82709, 69269, 54827, 42731, 29637, 18179),
        964: (114570, 98475, 81242, 68059, 53888, 42006, 29139, 17875),
        968: (112356, 96616, 79758, 66834, 52936, 41270, 28634, 17566),
        972: (110122, 94737, 78255, 65593, 51971, 40524, 28122, 17253),
        976: (107865, 92838, 76733, 64334, 50991, 39766, 27601, 16935),
        980: (105584, 90915, 75189, 63057, 49995, 38996, 27071, 16611),
        984: (103278, 88969, 73623, 61759, 48983, 38213, 26532, 16281),
        988: (100946, 86996, 72033, 60441, 47953, 37415, 25983, 15945),
        992: (98585, 84997, 70417, 59100, 46904, 36602, 25422, 15602),
        996: (96193, 82968, 68773, 57734, 45834, 35772, 24850, 15252),
        1000: (93771, 80908, 67098, 56341, 44741, 34923, 24264, 14893),
    },
}
ANCHOR_CENTRES = {
    2: (((14996, 141139), (2523, 97787), (-3745, 62352), (-6860, 30443), (-7813, 0)), ((16157, 144736), (3867, 104181), (-2615, 71115), (-6154, 41546), (-7743, 13676))),
    3: (((10264, 117601), (1774, 82061), (-2694, 52579), (-4980, 25743), (-5689, 0)), ((11016, 120466), (2709, 87297), (-1874, 59877), (-4452, 35089), (-5630, 11568))),
    4: (((7798, 102919), (1368, 72100), (-2105, 46316), (-3909, 22711), (-4473, 0)), ((8353, 105362), (2084, 76639), (-1461, 52703), (-3488, 30935), (-4423, 10207))),
    5: (((6287, 92645), (1112, 65066), (-1728, 41865), (-3218, 20547), (-3686, 0)), ((6725, 94807), (1693, 69127), (-1198, 47614), (-2868, 27976), (-3643, 9235))),
}
EXPECTED_CELLS = {
    (2, 904, 9): (99, 53, 82, 56, 57, 43, 39, 22),
    (2, 904, 10): (143, 69, 106, 71, 89, 58, 49, 27),
    (2, 908, 9): (100, 53, 82, 56, 57, 43, 40, 22),
    (2, 908, 10): (143, 69, 106, 75, 89, 58, 49, 27),
    (2, 912, 9): (99, 53, 82, 56, 59, 43, 39, 23),
    (2, 912, 10): (142, 69, 107, 78, 89, 58, 49, 28),
    (2, 916, 9): (91, 53, 82, 56, 59, 43, 39, 23),
    (2, 916, 10): (138, 69, 109, 79, 89, 58, 49, 28),
    (2, 920, 9): (90, 53, 82, 56, 63, 43, 39, 23),
    (2, 920, 10): (127, 68, 110, 79, 89, 56, 49, 28),
    (2, 924, 9): (89, 53, 82, 56, 63, 43, 39, 23),
    (2, 924, 10): (124, 69, 117, 79, 89, 56, 49, 28),
    (2, 928, 9): (87, 53, 81, 56, 63, 43, 39, 23),
    (2, 928, 10): (124, 69, 120, 80, 89, 55, 50, 28),
    (2, 932, 9): (81, 53, 81, 56, 63, 43, 38, 23),
    (2, 932, 10): (124, 69, 128, 81, 89, 55, 50, 28),
    (2, 936, 9): (81, 53, 81, 56, 63, 43, 38, 23),
    (2, 936, 10): (124, 68, 129, 83, 89, 55, 48, 28),
    (2, 940, 9): (81, 53, 81, 56, 64, 43, 38, 23),
    (2, 940, 10): (115, 68, 125, 87, 89, 55, 48, 28),
    (2, 944, 9): (81, 53, 81, 56, 64, 44, 38, 23),
    (2, 944, 10): (113, 68, 126, 87, 89, 55, 48, 28),
    (2, 948, 9): (82, 53, 81, 56, 64, 44, 38, 23),
    (2, 948, 10): (113, 68, 127, 87, 89, 55, 48, 28),
    (2, 952, 9): (81, 53, 81, 56, 64, 44, 39, 23),
    (2, 952, 10): (113, 68, 128, 87, 90, 55, 48, 29),
    (2, 956, 9): (82, 52, 81, 56, 62, 46, 39, 23),
    (2, 956, 10): (113, 68, 129, 87, 90, 55, 48, 29),
    (2, 960, 9): (81, 52, 82, 56, 61, 46, 39, 23),
    (2, 960, 10): (110, 68, 124, 88, 90, 55, 49, 29),
    (2, 964, 9): (80, 52, 82, 56, 62, 46, 40, 23),
    (2, 964, 10): (111, 68, 123, 90, 90, 55, 51, 29),
    (2, 968, 9): (78, 53, 83, 56, 61, 47, 41, 23),
    (2, 968, 10): (113, 69, 121, 93, 90, 55, 51, 29),
    (2, 972, 9): (76, 53, 83, 56, 64, 47, 41, 23),
    (2, 972, 10): (113, 70, 122, 93, 91, 55, 51, 31),
    (2, 976, 9): (79, 53, 84, 56, 65, 48, 43, 23),
    (2, 976, 10): (114, 73, 129, 93, 94, 55, 51, 33),
    (2, 980, 9): (83, 53, 86, 56, 65, 49, 43, 23),
    (2, 980, 10): (112, 75, 129, 92, 94, 55, 51, 33),
    (2, 984, 9): (84, 55, 89, 56, 64, 50, 43, 23),
    (2, 984, 10): (107, 75, 131, 91, 94, 55, 51, 35),
    (2, 988, 9): (84, 55, 93, 56, 65, 50, 43, 23),
    (2, 988, 10): (108, 75, 131, 92, 94, 56, 51, 37),
    (2, 992, 9): (84, 57, 94, 56, 66, 50, 43, 23),
    (2, 992, 10): (108, 75, 132, 92, 94, 59, 51, 37),
    (2, 996, 9): (84, 58, 94, 56, 66, 50, 43, 23),
    (2, 996, 10): (110, 75, 133, 94, 94, 63, 51, 38),
    (3, 952, 9): (105, 59, 88, 62, 73, 50, 42, 22),
    (3, 952, 10): (128, 86, 156, 104, 100, 57, 53, 33),
    (3, 956, 9): (105, 59, 89, 62, 73, 50, 42, 23),
    (3, 956, 10): (125, 86, 156, 105, 100, 57, 54, 35),
    (3, 960, 9): (105, 59, 89, 62, 73, 50, 42, 23),
    (3, 960, 10): (118, 85, 156, 106, 100, 58, 54, 37),
    (3, 964, 9): (105, 61, 92, 62, 74, 50, 43, 23),
    (3, 964, 10): (122, 86, 159, 107, 100, 58, 54, 40),
    (3, 968, 9): (104, 61, 92, 62, 75, 50, 44, 24),
    (3, 968, 10): (123, 87, 160, 108, 100, 59, 54, 40),
    (3, 972, 9): (96, 61, 92, 62, 77, 49, 44, 26),
    (3, 972, 10): (133, 88, 168, 108, 100, 62, 54, 40),
    (3, 976, 9): (96, 62, 92, 62, 78, 49, 44, 26),
    (3, 976, 10): (134, 91, 168, 109, 101, 68, 54, 40),
    (3, 980, 9): (93, 62, 100, 64, 79, 48, 44, 26),
    (3, 980, 10): (134, 96, 168, 109, 102, 72, 55, 41),
    (3, 984, 9): (93, 62, 100, 72, 79, 48, 44, 27),
    (3, 984, 10): (132, 96, 168, 109, 104, 74, 56, 42),
    (3, 988, 9): (89, 63, 104, 80, 79, 48, 43, 27),
    (3, 988, 10): (131, 96, 168, 112, 112, 76, 59, 42),
    (3, 992, 9): (89, 64, 116, 81, 79, 48, 43, 28),
    (3, 992, 10): (128, 100, 170, 113, 118, 76, 61, 42),
    (3, 996, 9): (90, 64, 117, 87, 81, 47, 43, 29),
    (3, 996, 10): (130, 100, 171, 116, 121, 79, 66, 42),
    (4, 970, 9): (99, 63, 108, 77, 80, 48, 46, 27),
    (4, 970, 10): (146, 106, 172, 115, 114, 84, 65, 43),
    (4, 974, 9): (99, 66, 112, 81, 81, 48, 46, 27),
    (4, 974, 10): (146, 109, 170, 116, 118, 86, 65, 43),
    (4, 978, 9): (102, 66, 116, 83, 86, 50, 46, 28),
    (4, 978, 10): (148, 109, 170, 116, 129, 86, 66, 43),
    (4, 982, 9): (100, 69, 117, 87, 86, 50, 46, 31),
    (4, 982, 10): (147, 112, 172, 116, 134, 86, 71, 43),
    (4, 986, 9): (96, 71, 117, 88, 84, 50, 46, 32),
    (4, 986, 10): (141, 114, 169, 120, 135, 91, 71, 43),
    (4, 990, 9): (96, 73, 117, 92, 86, 50, 46, 33),
    (4, 990, 10): (152, 118, 173, 120, 136, 94, 72, 43),
    (4, 994, 9): (104, 74, 119, 92, 88, 51, 45, 34),
    (4, 994, 10): (148, 120, 179, 122, 139, 97, 75, 44),
    (4, 998, 9): (66, 52, 75, 60, 56, 38, 34, 26),
    (4, 998, 10): (88, 78, 104, 77, 83, 63, 50, 32),
    (5, 952, 9): (115, 63, 108, 67, 84, 50, 46, 30),
    (5, 952, 10): (161, 110, 180, 118, 114, 88, 67, 44),
    (5, 956, 9): (115, 64, 111, 76, 84, 50, 46, 31),
    (5, 956, 10): (171, 113, 178, 120, 125, 88, 67, 44),
    (5, 960, 9): (115, 65, 112, 80, 84, 50, 46, 33),
    (5, 960, 10): (170, 113, 178, 120, 127, 88, 70, 44),
    (5, 964, 9): (115, 64, 123, 82, 84, 50, 46, 33),
    (5, 964, 10): (167, 113, 178, 120, 134, 88, 75, 44),
    (5, 968, 9): (112, 65, 123, 86, 86, 50, 45, 34),
    (5, 968, 10): (160, 117, 178, 120, 139, 93, 76, 44),
    (5, 972, 9): (114, 72, 124, 87, 86, 50, 45, 35),
    (5, 972, 10): (165, 117, 181, 120, 139, 96, 78, 44),
    (5, 976, 9): (110, 75, 124, 88, 88, 50, 45, 35),
    (5, 976, 10): (173, 117, 182, 120, 143, 96, 77, 46),
    (5, 980, 9): (106, 78, 132, 92, 90, 51, 44, 35),
    (5, 980, 10): (174, 120, 183, 120, 149, 97, 77, 52),
    (5, 984, 9): (119, 79, 133, 92, 89, 53, 44, 35),
    (5, 984, 10): (169, 121, 185, 120, 151, 97, 80, 54),
    (5, 988, 9): (131, 80, 133, 92, 89, 55, 46, 35),
    (5, 988, 10): (170, 123, 188, 120, 153, 96, 82, 52),
    (5, 992, 9): (120, 86, 137, 96, 88, 61, 46, 36),
    (5, 992, 10): (173, 128, 188, 122, 159, 100, 82, 53),
    (5, 996, 9): (113, 87, 143, 97, 95, 67, 48, 38),
    (5, 996, 10): (185, 129, 197, 133, 162, 103, 82, 55),
}



@lru_cache(maxsize=None)
def ratio_bounds(cap):
    lo,hi=log_bounds(1+8*cap)
    lower=log_bounds(1+Q(81,10)*cap)[0]/hi
    upper=1+log_bounds(1/(1-cap/10))[1]/lo
    require(lower<=REFERENCE_RATIO<=upper,"linked outer ratio")
    return interval(lower,upper)


@lru_cache(maxsize=8192)
def model(order,base,v):
    c9 = min(RESIDUAL_CAPS[base], (v[0]**-2-1)/8)
    cap=c9 if order==9 else 9*c9/(10-c9)
    mid=sum(v)/2
    a=tuple(interval(*q) for q in marked_core(order,base))
    # Marked core already includes all binomial coefficients.
    e0=1 if order==9 else REFERENCE_RATIO
    reference=tuple(coefficient*precise_power(mid,e0*Q(k*(order-k),order-1))
                    if 0<k<order else ONE for k,coefficient in enumerate(a))
    W=-log(interval(*v))
    mean=log(interval(v[0]/mid,v[1]/mid))
    if order==10:
        flat=rational(REFERENCE_RATIO)
        rho=ratio_bounds(c9)
        mean=flat*mean+(flat-rho)*W
        W=rho*W
    m=order-1
    T=rational(m*cap/(2+m*cap))
    gauss,skew,defect,sizes,distances=shape(order)
    vectors=(gauss,*skew,*defect)
    parameters=(mean,
                *((T**(2*k-1)*W*scalar(1,2*k)).symmetric_error() for k in range(1,4)),
                *((T**(2*k)*W*scalar(1,2*k+1)).positive_error() for k in range(1,3)))
    skew_tail=T**7*W/(scalar(8)*(ONE-T.square()))
    defect_tail=T**6*W/(scalar(7)*(ONE-T.square()))
    errors=[]
    for k,c in enumerate(reference):
        tail=sizes[k]*skew_tail+scalar(distances[k])*defect_tail
        radius=tail
        for parameter,vector in zip(parameters,vectors,strict=True):
            z=parameter*vector[k]
            radius+=I(0,max(-z.lo,z.hi))
        require(radius.hi<GRID,("marked exponent radius",v,order))
        err=radius**3*scalar(1,6)/(ONE-radius)+tail*(ONE+radius)
        errors.append(c*(ONE+err.symmetric_error()))

    mat=matrix(order)
    lam = ROOT_P
    def transform(vec):
        pairs=tuple(enumerate(vec))
        return tuple(projected_interval(pairs,row)*lam**k for k,row in enumerate(mat))
    groups=[(ONE,transform(errors))]
    groups.extend((p,transform(tuple(a*b for a,b in zip(reference,d,strict=True))))
                  for p,d in zip(parameters,vectors,strict=True))
    for i,j in combinations_with_replacement(range(len(vectors)),2):
        p=parameters[i]*parameters[j] if i!=j else scalar(1,2)*parameters[i].square()
        d=tuple(c*a*b for c,a,b in zip(reference,vectors[i],vectors[j],strict=True))
        groups.append((p,transform(d)))
    return tuple(groups)


@lru_cache(maxsize=8192)
def taylor_form(order,base,v):
    data=model(order,base,v)
    core=tuple(enumerate(data[0][1]))
    groups=tuple((h,tuple(enumerate(a))) for h,a in data[1:])
    full=[]
    for k,b in core:
        for h,a in groups:
            b+=h*a[k][1]
        full.append((k,b))
    safe=full[-1][1].lo*DISK**order>sum(max(abs(a.lo),abs(a.hi))*DISK**k for k,a in full[:-1])
    return core,groups,tuple(full),safe


@lru_cache(maxsize=8192)
def circle_safe(base,order,v):
    core,groups,full,safe=taylor_form(order,base,v)
    if safe:
        return True
    radius_square=rational(DISK**2)
    pending=[(interval(-DISK,DISK),0)]
    while pending:
        x,depth=pending.pop()
        delta=radius_square-x.square()
        y=I(max(0,delta.lo),max(0,delta.hi)).sqrt()
        if accurate_projection(order,x,y,core,groups,full):
            continue
        if depth>=12:
            return False
        pending.extend((h,depth+1) for h in bisect(x))
    return True


def verify(base,order,left,right,i):
    table=RAYS[base]
    a,b=table[left][i],table[right][i]
    pending=[((Q(left,1000),Q(right,1000)),interval(-DISK,DISK),0)]
    count=0
    while pending:
        v,x,depth=pending.pop()
        core,groups,full,safe=taylor_form(order,base,v)
        if not safe and not circle_safe(base,order,v):
            require(depth<25,("wide disk",base,order,v))
            pending.extend((h,x,depth+1) for h in bisect_u(v))
            continue
        levels=sorted(table_value(end,left,right,a,b) for end in v)
        y=height(P,x,interval(*levels))
        if accurate_projection(order,x,y,core,groups,full):
            count+=1
            continue
        require(depth<36,("unresolved wide ray",base,order,left,right,i,v,x.shown()))
        if 35*(v[1]-v[0])>Q(x.hi-x.lo,GRID):
            pending.extend((h,x,depth+1) for h in bisect_u(v))
        else:
            pending.extend((v,h,depth+1) for h in bisect(x))
    return count



def verify_anchor(base, order):
    core = tuple(enumerate(interval(*c) for c in marked_core(order, base)))
    coeffs = tuple(projected_interval(core, row) * ROOT_P**k
                   for k,row in enumerate(matrix(order)))
    box = ((ONE, coeffs),)
    rays = tuple(rational(Q(x,SCALE)) for x in RAYS[base][1000])
    sectors = (1,3,5,7) if order == 9 else (0,2,4,6,8)
    centers = tuple((Q(x,SCALE),Q(y,SCALE)) for x,y in ANCHOR_CENTRES[base][order-9])
    located = []
    for i,(x,y) in enumerate(centers):
        accept_sector = lambda r: sector_check((x,y),rays,sectors[i],r) if i<len(sectors) else y==0
        radius = next((r for r in RADII if x*x+y*y<(DISK-r)**2 and accept_sector(r)
                       and marked_disk_check(order,(x,y),box,r)),None)
        require(radius is not None,("marked anchor disk",base,order,x,y))
        located.append((x,y,radius))
    full=located+[(x,-y,r) for x,y,r in located if y>0]
    require(len(full)==order,"anchor root count")
    require(all((x-a)**2+(y-b)**2>(r+s)**2
                for i,(x,y,r) in enumerate(full) for a,b,s in full[:i]),
            "anchor disks disjoint")


def verify_thresholds(base):
    u=OUTER_CUTOFFS[base]
    v=RESIDUAL_CUTOFFS[base]
    one_row=Q(9*base+1,9*base+9)
    require(u>one_row and u*u>one_row*v*v,"mark uniqueness and residual cutoff")
    cap=RESIDUAL_CAPS[base]
    first=(cap.denominator-1)//9
    require(Q(9*first+1)==1/cap,"unit cap form")
    require(Q(9*first-8,9*first)<v*v and (1-v)/(9*v-1)<cap,"ordinary residual cap")
    require(DISK**2*COORDINATE_P<1,"Cayley disk in left half-plane")
    require(Q(min(RAYS[base]),1000)==v and max(RAYS[base])==1000,"residual interval")
    for rays in RAYS[base].values():
        levels=(2*SCALE,*rays,0)
        require(all(a>b for a,b in zip(levels,levels[1:])),"strict ray ordering")


def main():
    bases=tuple(RAYS) if len(argv)<2 else (int(argv[1]),)
    orders=ORDERS if len(argv)<3 else (int(argv[2]),)
    require(all(b in RAYS for b in bases) and all(j in ORDERS for j in orders),"replay selection")
    cells=0
    for base in bases:
        verify_thresholds(base)
        knots=tuple(RAYS[base])
        for order in orders:
            verify_anchor(base,order)
            for left,right in zip(knots,knots[1:]):
                counts=tuple(verify(base,order,left,right,i) for i in range(8))
                require(counts==EXPECTED_CELLS[base,left,order],
                        ("replay counts",base,left,right,order,counts))
                cells+=sum(counts)
                print(base,order,left,right,counts,flush=True)
    print(f"Certified {cells} marked ray cells.",flush=True)


if __name__=="__main__":
    main()
