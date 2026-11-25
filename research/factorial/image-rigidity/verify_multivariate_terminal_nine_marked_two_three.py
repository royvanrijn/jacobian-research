#!/usr/bin/env python3
"""Exact terminal-nine angular separation with a marked (2,3) row."""

from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from sys import argv

from verify_multivariate_terminal_eight_late_bridge import marked_disk_check
from verify_multivariate_terminal_nine_private_linear import ratio_bounds, shape
from verify_multivariate_terminal_nine_upper_bridge import projected_sign
from verify_multivariate_terminal_nine_upper_marked_rows import (
    marked_core, matrix, precise_power, sector_check, table_value,
)
from verify_multivariate_terminal_seven_upper_bridge import bisect_u, height
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, bisect, interval, log, rational, require, scalar,
)
from verify_multivariate_terminal_six_marked_row_tail import projected_interval


SCALE = 100000
ORDERS = (8, 9)
REFERENCE_RATIO = Q(64, 63)
P = rational(Q(1, 8))
SQUARE_ROOT_P = P.sqrt()
DISK = Q(11, 4)

RAYS = {
    903: (181624, 151417, 119224, 96144, 71418, 51302, 29450),
    907: (179877, 150039, 118225, 95365, 70864, 50912, 29230),
    911: (178131, 148659, 117223, 94583, 70308, 50519, 29009),
    915: (176387, 147279, 116219, 93799, 69749, 50125, 28787),
    919: (174643, 145898, 115212, 93012, 69188, 49729, 28564),
    923: (172900, 144516, 114202, 92222, 68625, 49331, 28340),
    924: (172464, 144170, 113949, 92024, 68483, 49231, 28284),
    928: (170722, 142786, 112936, 91231, 67916, 48830, 28058),
    932: (168980, 141401, 111919, 90434, 67347, 48427, 27831),
    936: (167238, 140014, 110899, 89635, 66775, 48022, 27603),
    940: (165497, 138625, 109876, 88832, 66199, 47615, 27373),
    944: (163755, 137234, 108849, 88025, 65621, 47206, 27141),
    948: (162013, 135841, 107818, 87215, 65039, 46794, 26909),
    952: (160271, 134445, 106784, 86401, 64454, 46379, 26674),
    956: (158528, 133047, 105745, 85583, 63865, 45962, 26438),
    960: (156784, 131647, 104701, 84761, 63273, 45541, 26201),
    964: (155040, 130243, 103653, 83934, 62676, 45118, 25961),
    968: (153294, 128837, 102601, 83103, 62076, 44692, 25720),
    972: (151547, 127427, 101543, 82267, 61472, 44262, 25476),
    976: (149800, 126014, 100480, 81426, 60863, 43829, 25231),
    980: (148051, 124597, 99411, 80580, 60249, 43393, 24983),
    984: (146300, 123176, 98336, 79727, 59630, 42953, 24733),
    988: (144548, 121752, 97255, 78870, 59006, 42508, 24481),
    992: (142795, 120323, 96168, 78005, 58377, 42060, 24226),
    996: (141040, 118890, 95074, 77135, 57741, 41607, 23968),
    1000: (139284, 117452, 93972, 76257, 57100, 41149, 23708),
}
ANCHOR_CENTRES = {
    8: ((13677, 136851), (1066, 90150), (-4867, 51822), (-7367, 16943)),
    9: ((14996, 141139), (2523, 97787), (-3745, 62352), (-6860, 30443), (-7813, 0)),
}
EXPECTED_CELLS = {
    903: ((2050, 231, 1403, 205, 238, 72, 22),
          (10352, 821, 4239, 664, 1281, 131, 50)),
    907: ((2046, 231, 1371, 212, 238, 72, 22),
          (10632, 834, 4308, 670, 1289, 135, 53)),
    911: ((1795, 229, 1329, 213, 237, 74, 22),
          (10259, 826, 4225, 670, 1305, 133, 55)),
    915: ((1705, 229, 1271, 212, 232, 74, 23),
          (9783, 831, 4176, 640, 1267, 132, 62)),
    919: ((1571, 228, 1206, 212, 229, 71, 24),
          (9193, 825, 4188, 636, 1228, 130, 62)),
    923: ((395, 70, 306, 65, 70, 29, 18),
          (2221, 213, 1051, 172, 307, 47, 28)),
    924: ((1363, 229, 1142, 212, 229, 66, 24),
          (8403, 804, 4059, 621, 1147, 132, 62)),
    928: ((1327, 231, 1092, 211, 222, 69, 25),
          (7823, 804, 3996, 614, 1083, 132, 62)),
    932: ((1188, 238, 1039, 211, 215, 66, 26),
          (7272, 783, 4028, 612, 1006, 131, 62)),
    936: ((1116, 240, 991, 206, 206, 66, 26),
          (6813, 762, 3826, 643, 961, 131, 62)),
    940: ((1008, 236, 976, 205, 210, 68, 26),
          (6336, 767, 3538, 648, 917, 131, 62)),
    944: ((941, 236, 924, 213, 211, 74, 29),
          (5701, 758, 3305, 655, 854, 131, 62)),
    948: ((939, 235, 839, 213, 208, 78, 30),
          (5140, 757, 3010, 651, 790, 131, 62)),
    952: ((834, 235, 759, 219, 206, 78, 33),
          (4553, 749, 2751, 651, 745, 143, 62)),
    956: ((751, 235, 713, 219, 206, 79, 34),
          (4013, 726, 2499, 651, 708, 143, 62)),
    960: ((729, 236, 715, 219, 204, 79, 34),
          (3499, 723, 2227, 663, 656, 143, 62)),
    964: ((717, 239, 681, 214, 198, 80, 34),
          (3080, 716, 1963, 674, 618, 145, 64)),
    968: ((703, 243, 642, 209, 192, 81, 34),
          (2741, 717, 1768, 681, 581, 151, 68)),
    972: ((700, 242, 602, 207, 189, 81, 34),
          (2436, 709, 1616, 689, 558, 151, 78)),
    976: ((620, 245, 592, 207, 188, 81, 35),
          (2270, 696, 1533, 700, 525, 151, 83)),
    980: ((590, 245, 554, 209, 188, 81, 35),
          (2180, 693, 1463, 714, 510, 151, 83)),
    984: ((579, 253, 528, 210, 189, 81, 35),
          (2056, 686, 1431, 744, 500, 151, 83)),
    988: ((620, 275, 541, 222, 202, 81, 37),
          (2048, 685, 1394, 740, 482, 151, 83)),
    992: ((617, 292, 563, 223, 206, 85, 41),
          (2031, 689, 1398, 743, 491, 153, 86)),
    996: ((568, 297, 575, 225, 210, 88, 42),
          (2048, 707, 1446, 769, 530, 159, 95)),
}


@lru_cache(maxsize=8192)
def coefficient_model(order, v):
    c8 = min(Q(1, 33), (v[0]**-2 - 1)/7)
    cap=c8 if order==8 else 8*c8/(9-c8)
    mid=sum(v)/2
    a=tuple(interval(*q) for q in marked_core(order, 2))
    # Marked core already includes all binomial coefficients.
    e0=1 if order==8 else REFERENCE_RATIO
    reference=tuple(coefficient*precise_power(mid,e0*Q(k*(order-k),order-1))
                    if 0<k<order else ONE for k,coefficient in enumerate(a))
    W=-log(interval(*v))
    mean=log(interval(v[0]/mid,v[1]/mid))
    if order==9:
        flat=rational(REFERENCE_RATIO)
        rho=ratio_bounds(c8)
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

    mat = matrix(order)
    lam = SQUARE_ROOT_P
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
def ray_model(order, v):
    contributions = coefficient_model(order, v)
    base = tuple(enumerate(contributions[0][1]))
    groups = tuple((parameter, tuple(enumerate(coefficients)))
                   for parameter, coefficients in contributions[1:])
    full = []
    for k, bound in base:
        for parameter, coefficients in groups:
            bound += parameter * coefficients[k][1]
        full.append((k, bound))
    dominates = full[-1][1].lo * DISK**order > sum(
        max(abs(bound.lo), abs(bound.hi)) * DISK**k for k, bound in full[:-1]
    )
    return base, groups, tuple(full), dominates


def verify_ray(order, left, right, index):
    ray_left, ray_right = RAYS[left][index], RAYS[right][index]
    pending = [((Q(left, 1000), Q(right, 1000)), interval(-DISK, DISK), 0)]
    accepted = 0
    while pending:
        v, x, depth = pending.pop()
        base, groups, full, dominates = ray_model(order, v)
        if not dominates:
            require(depth < 18, ("root disk", order, v))
            pending.extend((half, x, depth+1) for half in bisect_u(v))
            continue
        levels = sorted(table_value(t, left, right, ray_left, ray_right) for t in v)
        y = height(P, x, interval(*levels))
        if projected_sign(order, x, y, base, groups, full, False):
            accepted += 1
            continue
        require(depth < 34, ("ray", order, left, right, index, v, x.shown()))
        if 80 * (v[1]-v[0]) > Q(x.hi-x.lo, GRID):
            pending.extend((half, x, depth+1) for half in bisect_u(v))
        else:
            pending.extend((v, half, depth+1) for half in bisect(x))
    return accepted


def verify_anchor(order):
    core = tuple(interval(*bounds) for bounds in marked_core(order, 2))
    terms = tuple(enumerate(core))
    coefficients = tuple(projected_interval(terms, row) * SQUARE_ROOT_P**degree
                         for degree, row in enumerate(matrix(order)))
    model = ((ONE, coefficients),)
    rays = tuple(rational(Q(b, SCALE)) for b in RAYS[1000])
    sectors = (1, 3, 5, 7) if order == 8 else (0, 2, 4, 6)
    radius = Q(13, 1000)
    centers = tuple((Q(x, SCALE), Q(y, SCALE)) for x, y in ANCHOR_CENTRES[order])
    full = centers + tuple((x, -y) for x, y in centers if y > 0)
    require(len(full) == order, "anchor root count")
    require(order == 8 or centers[-1][1] == 0, "real nonic anchor")
    for i, (x, y) in enumerate(full):
        require(x*x + y*y < (DISK-radius)**2, "anchor disk inside the Cayley disk")
        require(all((x-a)**2+(y-b)**2 > (2*radius)**2 for a,b in full[:i]),
                "disjoint anchor disks")
    for center, sector in zip(centers[:4], sectors, strict=True):
        require(sector_check(center, rays, sector, radius), ("anchor sector", order, center))
    for center in centers:
        require(marked_disk_check(order, center, model, radius),
                ("anchor root disk", order, center))


def verify_thresholds():
    u = Q(19, 25)
    v = Q(903, 1000)
    unit_square = Q(17, 24)
    require(u > unit_square and u*u > unit_square * v*v, "unique mark and residual cutoff")
    require(Q(25, 32) < v*v and (1-v)/(8*v-1) < Q(1, 33), "ordinary residual cap")
    require(DISK**2/8 < 1, "Cayley disk in the left half-plane")
    for rays in RAYS.values():
        levels = (int(DISK*SCALE), *rays, 0)
        require(all(a > b for a,b in zip(levels, levels[1:])), "strict ray order")


def main():
    verify_thresholds()
    knots = tuple(RAYS)
    first, last = (knots[0], knots[-1]) if len(argv) == 1 else tuple(map(int, argv[1:3]))
    require(first in knots and last in knots and first < last, "replay band")
    orders = ORDERS if len(argv) < 4 else (int(argv[3]),)
    require(all(order in ORDERS for order in orders), "replay order")
    total = 0
    for order in orders:
        verify_anchor(order)
        for left,right in zip(knots, knots[1:]):
            if first <= left < right <= last:
                counts = tuple(verify_ray(order, left, right, i) for i in range(7))
                require(counts == EXPECTED_CELLS[left][order-8],
                        ("replay cell count", order, left, right, counts))
                total += sum(counts)
                print(order, left, right, counts, flush=True)
    print(f"Certified {total} Cayley projection cells.", flush=True)


if __name__ == "__main__":
    main()
