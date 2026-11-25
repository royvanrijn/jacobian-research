"""Exact conditioned Cayley separation on the XII band .66–.70."""

from fractions import Fraction as Q
from sys import argv

from verify_multivariate_terminal_twelve_cayley import (
    anchor_for, check_ray, limits_for, require, validate_rays,
)
from verify_multivariate_terminal_twelve_cayley_lower import check_leading


CAP = Q(1, 28)
DISK = Q(3)

RAYS = {
    300: (145301, 126826, 104303, 90136, 73474, 65245, 47811, 38321, 23757, 13671),
    305: (145961, 127322, 104607, 89545, 73620, 65364, 47879, 38372, 23784, 13292),
    310: (146630, 127824, 105076, 89770, 73766, 65483, 47948, 38422, 23811, 13306),
    315: (147307, 128332, 105387, 89997, 73914, 65603, 48018, 38472, 23837, 13320),
    320: (148081, 128846, 105700, 90226, 74063, 65724, 48087, 38523, 23864, 13334),
    325: (148776, 129366, 106017, 90457, 74213, 65845, 48157, 38575, 23609, 13348),
    330: (149480, 129892, 106336, 90690, 74364, 65968, 48228, 38626, 23636, 13362),
    335: (150194, 130425, 106659, 90924, 74516, 66091, 48299, 38678, 23663, 12980),
    340: (150916, 130963, 106984, 91161, 74670, 66216, 48371, 38730, 23690, 12994),
}

EXPECTED_CELLS = {
    (11, 300): (1423, 40, 318, 55, 103, 41, 61, 26, 30, 17),
    (11, 305): (1335, 38, 263, 50, 94, 41, 61, 26, 30, 16),
    (11, 310): (1282, 37, 236, 50, 90, 41, 61, 24, 30, 16),
    (11, 315): (1141, 37, 205, 50, 90, 41, 57, 24, 30, 16),
    (11, 320): (1001, 38, 204, 50, 90, 41, 57, 24, 30, 16),
    (11, 325): (986, 36, 201, 50, 90, 41, 57, 24, 30, 16),
    (11, 330): (980, 35, 196, 49, 90, 41, 55, 23, 30, 16),
    (11, 335): (975, 35, 195, 49, 90, 41, 55, 23, 30, 16),
    (12, 300): (4255, 66, 887, 183, 396, 66, 104, 39, 32, 22),
    (12, 305): (3376, 66, 1157, 180, 361, 66, 98, 39, 32, 22),
    (12, 310): (2923, 67, 1250, 181, 327, 64, 80, 34, 30, 22),
    (12, 315): (3614, 68, 1107, 189, 361, 70, 83, 33, 28, 22),
    (12, 320): (3387, 66, 836, 189, 331, 65, 83, 29, 27, 22),
    (12, 325): (2972, 61, 647, 181, 327, 60, 75, 28, 25, 21),
    (12, 330): (3253, 59, 551, 181, 315, 57, 61, 27, 24, 21),
    (12, 335): (10140, 59, 517, 181, 310, 56, 60, 25, 24, 20),
}


def verify_cap():
    require(min(Q(1,23),Q(2,45))>Q(1,25),"one- and two-factor row cutoffs")
    require(Q(1,34)<CAP and Q(2,56)==CAP,"one- and two-factor lattice cap")
    require((1+10*CAP/(1+CAP))**3>Q(50,33)**2,"three-factor mean cutoff")


def main():
    order, ray, low, high = map(int,argv[1:]) if len(argv)>1 else (0,-1,300,340)
    require(order in (0,11,12) and -1<=ray<10,"ray selection")
    require(low in RAYS and high in RAYS and low<high,"range")
    require(tuple(RAYS)==tuple(range(300,341,5)),"contiguous ray knots")
    verify_cap()
    anchor_for(limits_for(CAP,300,300,DISK,True),RAYS,300)
    total=0
    for left in range(low,high,5):
        right=left+5
        limits=limits_for(CAP,left,right,DISK,True)
        validate_rays(limits,RAYS,left,right)
        for n in (11,12):
            if order and order!=n:
                continue
            for index in range(10):
                if ray>=0 and ray!=index:
                    continue
                if index==0:
                    count=check_leading(n,left,right,RAYS,CAP)
                else:
                    linearized=index==1 or (n==12 and index==2)
                    count=check_ray(n,index,left,right,limits,RAYS,linearized)
                require(count==EXPECTED_CELLS[n,left][index],
                        ("replay cell count",n,left,index,count))
                print(n,left,right,index,count,flush=True)
                total+=count
    print(f"Verified {total} projection cells.",flush=True)


if __name__=="__main__":
    main()
