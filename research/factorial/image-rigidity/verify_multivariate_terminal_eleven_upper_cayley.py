#!/usr/bin/env python3
"""Exact cap-sensitive Cayley separation from outer mean .88 to .94."""

from fractions import Fraction as Q
from sys import argv

from verify_multivariate_terminal_eleven_middle_cayley import (
    ORDERS, REFERENCE_RATIO, Limits, check_ray, reciprocal_cap, verify_anchor, verify_constants,
)
from verify_multivariate_terminal_nine_upper_marked_rows import precise_power
from verify_multivariate_terminal_six_intermediate_cayley import GRID, rational, require
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds


LIMITS = Limits(Q(1,70), (Q(10094,10000),Q(4043,4000)), Q(1213,10000), Q(12,5), False)
TIGHT = Limits(Q(1,70), (Q(50487,50000),Q(20209,20000)), Q(66,1000), Q(12,5), False)

RAYS = {
    60: (119893, 108027, 88643, 81633, 61868, 55628, 37667, 30155, 14497),
    65: (120281, 108338, 88834, 81795, 61959, 55061, 37708, 30186, 14510),
    70: (120674, 108652, 89027, 81959, 62051, 55137, 37749, 30747, 14523),
    75: (121159, 110060, 89221, 82123, 62144, 55213, 37790, 30778, 14536),
    80: (121559, 110387, 89416, 81478, 62237, 55289, 37831, 30809, 14548),
    85: (121964, 110716, 89613, 81641, 62565, 55366, 37872, 30841, 14195),
    90: (122463, 111048, 89811, 81806, 62660, 55443, 37914, 30872, 14207),
    95: (122876, 111383, 90180, 81972, 62755, 55520, 37955, 30904, 14220),
    100: (123292, 111720, 90382, 82138, 62851, 55598, 37997, 30936, 14233),
    105: (123713, 112061, 90585, 82306, 62948, 55676, 38039, 30968, 14245),
    110: (124137, 112405, 90789, 82475, 63045, 55755, 38082, 31000, 14258),
    115: (124662, 112752, 90995, 82645, 62904, 55834, 38124, 31033, 14640),
    120: (125096, 113102, 91202, 82817, 63001, 55914, 38167, 31065, 14653),
}
EXPECTED_CELLS = {
    (10, 60): (1636, 67, 216, 96, 186, 86, 100, 48, 26),
    (10, 65): (1099, 67, 190, 85, 162, 66, 84, 49, 25),
    (10, 70): (666, 70, 171, 74, 124, 58, 76, 42, 24),
    (10, 75): (521, 69, 170, 71, 100, 57, 70, 42, 23),
    (10, 80): (476, 69, 168, 67, 95, 56, 52, 42, 21),
    (10, 85): (357, 68, 151, 67, 95, 56, 44, 42, 21),
    (10, 90): (273, 67, 142, 62, 95, 56, 44, 41, 21),
    (10, 95): (258, 67, 118, 62, 95, 56, 44, 41, 21),
    (10, 100): (241, 67, 113, 61, 95, 56, 44, 41, 21),
    (10, 105): (241, 67, 105, 60, 95, 56, 40, 41, 21),
    (10, 110): (155, 42, 71, 39, 55, 36, 27, 27, 15),
    (10, 115): (126, 42, 71, 39, 54, 36, 27, 25, 15),
    (11, 60): (39411, 122, 1413, 157, 470, 88, 207, 56, 49),
    (11, 65): (80580, 110, 824, 153, 327, 73, 114, 56, 47),
    (11, 70): (46303, 109, 483, 152, 177, 73, 106, 48, 43),
    (11, 75): (25894, 95, 397, 150, 133, 73, 78, 41, 35),
    (11, 80): (15985, 89, 397, 141, 130, 73, 57, 36, 29),
    (11, 85): (14134, 86, 395, 141, 129, 73, 56, 35, 26),
    (11, 90): (12030, 84, 383, 141, 129, 69, 56, 35, 26),
    (11, 95): (7515, 77, 382, 140, 129, 69, 56, 35, 26),
    (11, 100): (4452, 77, 346, 129, 129, 68, 56, 34, 26),
    (11, 105): (3829, 77, 338, 121, 121, 67, 56, 34, 26),
    (11, 110): (3900, 58, 203, 79, 72, 43, 35, 24, 18),
    (11, 115): (3649, 56, 204, 79, 72, 43, 35, 24, 18),
}



def verify_tight_ratio():
    u_min,u_max = Q(187,200),Q(47,50)
    lo,hi = TIGHT.ratio
    cap = reciprocal_cap(u_min,TIGHT.cap)
    r = REFERENCE_RATIO
    require(lo <= r*(1-(r-1)*(1-u_min)/(2*u_min)), "tight ratio lower")
    log_ratio = Q(4043,4000)
    require(log_ratio >= 1+log_bounds(1/(1-cap/11))[1]/log_bounds(1+9*cap)[0],
            "rational logarithmic upper bound")
    power_lower = Q(precise_power(u_max,log_ratio).lo,GRID)
    require((1-power_lower)/(1-u_max) <= hi, "tight ratio upper")
    scale = rational(1-u_min).sqrt()*rational(hi).sqrt()
    require(Q(scale.square().hi,GRID)<TIGHT.qmax, "tight deficit")
    require(TIGHT.disk*TIGHT.disk*TIGHT.qmax<1, "tight left half-plane")
    common = Q(101,100)
    require(LIMITS.ratio[0] <= lo <= common <= hi <= LIMITS.ratio[1],
            "common comparison ratio")


def main():
    verify_constants(LIMITS,RAYS)
    verify_tight_ratio()
    verify_anchor(LIMITS,RAYS)
    n,i = map(int,argv[1:3]) if len(argv)>1 else (0,-1)
    lower,upper = map(int,argv[3:]) if len(argv)>3 else (60,120)
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
                tight = (order,index,a)==(11,0,60)
                count = check_ray(order,index,a,b,TIGHT if tight else LIMITS,RAYS,tight)
                require(count==EXPECTED_CELLS[order,a][index],
                        ("replay cell count",order,a,b,index,count))
                print(order,a,b,index,count,flush=True)
                total += count
    print(f"Verified {total} projection cells.",flush=True)


if __name__=="__main__":
    main()
