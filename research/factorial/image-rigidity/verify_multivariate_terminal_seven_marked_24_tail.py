#!/usr/bin/env python3
"""Exact terminal-seven separation with one marked (2,4) shared row."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial
from sys import argv

from verify_multivariate_terminal_seven_marked_light_tail import (
    BANDS, coefficient_model, core_weights, verify_constants,
)
from verify_multivariate_terminal_seven_paired_light import envelope, radial
from verify_multivariate_terminal_six_intermediate_cayley import (
    ONE, ZERO, complex_mul, interval, modulus, nth_root, require, scalar,
)



SINGLE_RAYS={
    149: (807, 600, -38, -402, -805),
    160: (799, 591, -47, -408, -808),
    180: (784, 557, -64, -420, -817),
    200: (767, 541, -87, -432, -821),
    220: (750, 541, -103, -430, -825),
    240: (730, 524, -114, -441, -829),
    260: (712, 489, -130, -452, -833),
    280: (695, 454, -146, -463, -836),
    300: (678, 437, -161, -521, -840),
    320: (662, 419, -177, -519, -843),
    340: (643, 401, -192, -528, -847),
    360: (627, 383, -207, -526, -850),
    380: (610, 348, -221, -536, -853),
    400: (593, 329, -236, -534, -856),
    420: (576, 311, -246, -543, -860),
    440: (558, 293, -260, -542, -863),
    460: (540, 274, -278, -551, -866),
    480: (522, 189, -292, -560, -872),
    500: (502, 171, -306, -569, -875),
    520: (482, 153, -320, -568, -878),
    540: (459, 151, -329, -577, -880),
    560: (436, 132, -342, -586, -883),
    580: (414, 114, -355, -595, -886),
    600: (392, 95, -372, -603, -888),
    620: (369, 61, -385, -612, -891),
    640: (346, 43, -398, -629, -893),
    660: (324, 40, -411, -629, -896),
    662: (322, 38, -412, -629, -896),
    675: (308, 26, -420, -643, -895),
    680: (302, 21, -423, -645, -895),
    685: (299, 16, -426, -647, -896),
    695: (288, 22, -433, -651, -900),
    700: (282, 17, -436, -653, -901),
    705: (276, -17, -439, -655, -901),
    715: (267, -26, -445, -659, -902),
    720: (262, -31, -448, -661, -903),
    725: (256, -50, -451, -663, -904),
    735: (245, -60, -457, -667, -905),
    740: (239, -50, -460, -669, -905),
    750: (228, -59, -467, -673, -907),
    760: (217, -68, -476, -676, -908),
    780: (195, -87, -488, -684, -910),
    800: (170, -105, -500, -692, -912),
    820: (147, -124, -512, -692, -915),
    840: (123, -143, -524, -700, -917),
    860: (100, -161, -535, -707, -919),
    880: (77, -180, -547, -715, -921),
    900: (52, -186, -559, -722, -925),
    920: (30, -205, -570, -730, -927),
    940: (7, -224, -582, -737, -930),
    960: (-15, -242, -596, -744, -932),
    980: (-38, -261, -607, -752, -934),
    1000: (-61, -269, -618, -759, -936),
}
SINGLE_DISKS={
    6: ((-87, 911), (-678, 763), (-1029, 292)),
    7: ((-25, 903), (-559, 835), (-941, 490), (-1078, 0)),
}

@lru_cache(maxsize=None)
def factors(order):
    return tuple(nth_root(Q(
        factorial(2*order+2*k)**order,
        factorial(2*order)**(order-k)*factorial(4*order)**k,
    ), order) for k in range(order+1))


@lru_cache(maxsize=None)
def relative_factors(order):
    return tuple(a/b for a,b in zip(factors(order),core_weights(order,(1,0)),strict=True))


def retained_model(order,u,env):
    relative=relative_factors(order)
    return tuple((parameter,tuple(a*b for a,b in zip(coefficients,relative,strict=True)))
                 for parameter,coefficients in coefficient_model(order,u,env,(1,0)))


def ray_cosines(left,right,index,u):
    a,b=SINGLE_RAYS[left][index],SINGLE_RAYS[right][index]
    cosines=tuple(Q(a,1000)+Q(b-a,right-left)*(endpoint-Q(left,1000)) for endpoint in u)
    return interval(min(cosines),max(cosines))


def verify(order,index,left,right,t):
    env=envelope(t,order)
    pending=[((Q(left,1000),Q(right,1000)),0)]
    count=0
    while pending:
        u,depth=pending.pop()
        used=radial(order,ray_cosines(left,right,index,u),retained_model(order,u,env))
        if used is not None:
            count+=used
            continue
        require(depth<15,("unresolved marked ray",order,index,left,right,t,[float(x) for x in u]))
        middle=(u[0]+u[1])/2
        pending.extend((((u[0],middle),depth+1),((middle,u[1]),depth+1)))
    return count


def check_anchor():
    for row in SINGLE_RAYS.values():
        require(all(1000>a>b>-1000 for a,b in zip(row,row[1:])),"strict ray order")
    radius=scalar(1,500)
    sector_order={6:(1,3,5),7:(0,2,4,5)}
    for order in (6,7):
        coefficients=tuple(scalar(comb(order,k))*f for k,f in enumerate(factors(order)))
        points=SINGLE_DISKS[order]
        for (a,b),sector in zip(points,sector_order[order],strict=True):
            x,y=scalar(a,1000),scalar(b,1000)
            powers=[(ONE,ZERO)]
            for _ in range(order):
                powers.append(complex_mul(powers[-1],(x,y)))
            terms=[]
            for k in range(order+1):
                real=imag=ZERO
                for h in range(k,order+1):
                    factor=scalar(comb(h,k))*coefficients[h]
                    p,q=powers[h-k]
                    real,imag=real+factor*p,imag+factor*q
                terms.append(modulus((real,imag)))
            remainder=terms[0]
            for k in range(2,order+1):
                remainder=remainder+terms[k]*radius**k
            require((radius*terms[1]).lo>remainder.hi,("Rouche disk",order,a,b))
            for k,c in enumerate(SINGLE_RAYS[1000]):
                cos=scalar(c,1000)
                sine=(ONE-cos.square()).sqrt()
                side=cos*y-sine*x
                require(side.lo>radius.hi if k<sector else side.hi<(-radius).lo,
                        ("sector disk",order,a,b,k))
            if b:
                require(y.lo>radius.hi,"upper disk")
        full=points+tuple((a,-b) for a,b in points if b)
        require(len(full)==order,"all roots isolated")
        for i,(a,b) in enumerate(full):
            for x,y in full[:i]:
                require((a-x)**2+(b-y)**2>16,"disjoint disks")


def main():
    verify_constants()
    check_anchor()
    for low,high,cap in BANDS:
        if len(argv)>1 and low!=int(argv[1]):
            continue
        knots=sorted(k for k in SINGLE_RAYS if max(149,low)<=k<=high)
        for left,right in zip(knots,knots[1:]):
            cells=[verify(j,ray,left,right,cap) for j in (6,7) for ray in range(5)]
            print(left,right,cap,cells,flush=True)


if __name__=="__main__":
    main()
