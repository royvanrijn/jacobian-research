#!/usr/bin/env python3
"""Exact terminal-seven separation with one or two marked shared rows."""
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from sys import argv

from verify_multivariate_terminal_seven_paired_light import (
    ANCHOR_RATIO, envelope, exponents, nominal, radial,
)
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, complex_mul, interval, log, marked_row, modulus, rational,
    require, scalar,
)



CORE_RAYS = {
    (1,0): {
        110: (808, 602, -41, -400, -810),
        120: (799, 592, -52, -408, -812),
        140: (778, 552, -75, -424, -818),
        150: (766, 541, -86, -418, -821),
        160: (754, 548, -97, -426, -828),
        180: (730, 525, -119, -441, -833),
        200: (704, 467, -140, -443, -838),
        220: (683, 444, -161, -446, -843),
        240: (659, 420, -181, -530, -847),
        260: (637, 396, -196, -531, -852),
        280: (616, 355, -216, -543, -856),
        300: (593, 331, -235, -545, -860),
        320: (569, 306, -255, -557, -864),
        340: (547, 282, -278, -558, -872),
        360: (524, 240, -296, -570, -876),
        380: (499, 183, -314, -572, -879),
        400: (470, 158, -332, -574, -883),
        420: (439, 134, -346, -586, -886),
        440: (407, 109, -363, -588, -890),
        460: (375, 85, -381, -600, -896),
        480: (345, 60, -402, -612, -899),
        500: (310, 66, -418, -623, -902),
        520: (285, 55, -435, -634, -905),
        540: (259, -88, -451, -645, -908),
        560: (233, -83, -471, -656, -911),
        580: (201, -93, -487, -659, -912),
        600: (166, -103, -503, -670, -915),
        620: (128, -101, -519, -680, -918),
        640: (96, -126, -531, -698, -920),
        660: (58, -152, -546, -701, -923),
        662: (53, -155, -548, -702, -923),
        675: (38, -171, -558, -715, -923),
        680: (29, -178, -562, -718, -924),
        685: (21, -171, -565, -720, -925),
        695: (4, -184, -573, -725, -926),
        700: (-4, -191, -577, -728, -927),
        705: (-12, -197, -580, -730, -927),
        715: (-25, -210, -588, -735, -929),
        720: (-33, -217, -592, -738, -929),
        725: (-39, -223, -595, -740, -930),
        735: (-53, -225, -603, -745, -931),
        740: (-61, -231, -607, -748, -932),
        750: (-75, -245, -617, -747, -933),
        760: (-89, -280, -624, -752, -935),
        780: (-118, -306, -639, -761, -937),
        800: (-147, -332, -653, -771, -940),
        820: (-177, -347, -667, -781, -942),
        840: (-208, -374, -681, -790, -945),
        860: (-241, -400, -695, -800, -948),
        880: (-274, -426, -709, -809, -950),
        900: (-308, -444, -723, -818, -952),
        920: (-340, -470, -739, -824, -955),
        940: (-372, -489, -753, -833, -957),
        960: (-405, -516, -768, -852, -961),
        980: (-437, -543, -781, -861, -963),
        1000: (-470, -596, -796, -867, -967),
    },
    (2,0): {
        150: (805, 597, -40, -404, -806),
        160: (798, 589, -49, -410, -808),
        180: (782, 555, -67, -422, -817),
        200: (763, 556, -89, -433, -822),
        220: (746, 538, -106, -432, -826),
        240: (728, 521, -117, -443, -829),
        260: (709, 468, -133, -454, -833),
        280: (692, 451, -149, -465, -837),
        300: (674, 433, -165, -523, -841),
        320: (657, 415, -180, -521, -844),
        340: (641, 396, -196, -531, -848),
        360: (622, 361, -211, -529, -851),
        380: (605, 342, -226, -538, -854),
        400: (587, 324, -240, -537, -857),
        420: (570, 305, -255, -546, -861),
        440: (552, 287, -269, -545, -864),
        460: (534, 251, -283, -555, -867),
        480: (515, 182, -297, -564, -873),
        500: (495, 164, -311, -563, -876),
        520: (474, 146, -325, -572, -879),
        540: (451, 143, -334, -581, -881),
        560: (428, 124, -348, -590, -884),
        580: (405, 106, -361, -599, -887),
        600: (381, 87, -378, -607, -889),
        620: (357, 53, -391, -616, -892),
        640: (335, 49, -404, -633, -895),
        660: (311, 30, -417, -633, -897),
        662: (309, 28, -418, -634, -897),
        675: (297, 16, -427, -647, -896),
        680: (293, 11, -430, -649, -899),
        685: (288, 21, -433, -651, -900),
        695: (276, -18, -439, -655, -901),
        700: (271, -22, -442, -657, -902),
        705: (267, -27, -445, -659, -903),
        715: (256, -36, -452, -663, -904),
        720: (250, -56, -455, -665, -904),
        725: (244, -60, -458, -667, -905),
        735: (233, -55, -464, -671, -906),
        740: (227, -60, -467, -673, -907),
        750: (216, -69, -473, -677, -908),
        760: (205, -79, -483, -681, -909),
        780: (182, -98, -495, -689, -911),
        800: (157, -117, -507, -696, -914),
        820: (132, -136, -519, -704, -916),
        840: (108, -154, -531, -705, -918),
        860: (83, -173, -543, -712, -920),
        880: (60, -192, -555, -720, -923),
        900: (36, -199, -566, -728, -925),
        920: (14, -218, -578, -735, -929),
        940: (-9, -237, -590, -743, -931),
        960: (-33, -256, -604, -750, -933),
        980: (-55, -264, -615, -757, -935),
        1000: (-79, -284, -627, -765, -937),
    },
    (1,1): {
        150: (805, 597, -40, -404, -806),
        160: (798, 589, -49, -410, -813),
        180: (782, 573, -66, -422, -817),
        200: (765, 556, -84, -433, -821),
        220: (746, 538, -100, -432, -825),
        240: (727, 521, -117, -443, -829),
        260: (709, 486, -133, -454, -833),
        280: (691, 468, -149, -453, -837),
        300: (674, 433, -164, -523, -845),
        320: (657, 415, -180, -532, -848),
        340: (640, 396, -190, -530, -851),
        360: (623, 378, -205, -540, -855),
        380: (604, 360, -225, -538, -858),
        400: (587, 324, -240, -547, -861),
        420: (569, 305, -254, -556, -864),
        440: (551, 286, -268, -555, -867),
        460: (533, 268, -282, -564, -870),
        480: (515, 249, -296, -563, -873),
        500: (494, 213, -310, -572, -876),
        520: (473, 194, -324, -571, -878),
        540: (450, 175, -337, -580, -881),
        560: (427, 156, -350, -589, -884),
        580: (403, 121, -359, -597, -886),
        600: (380, 102, -373, -606, -892),
        620: (356, 83, -385, -615, -894),
        640: (334, 64, -398, -623, -894),
        660: (309, 60, -411, -631, -896),
        662: (307, 58, -412, -632, -897),
        675: (296, 60, -421, -645, -898),
        680: (291, 55, -424, -647, -899),
        685: (285, 50, -427, -649, -899),
        695: (274, 55, -433, -653, -900),
        700: (270, 50, -436, -655, -901),
        705: (264, 45, -439, -657, -902),
        715: (253, -66, -445, -661, -903),
        720: (247, -56, -449, -663, -903),
        725: (241, -61, -452, -665, -904),
        735: (230, -70, -458, -669, -905),
        740: (224, -75, -461, -671, -906),
        750: (213, -70, -470, -675, -907),
        760: (201, -79, -476, -678, -908),
        780: (176, -85, -488, -686, -910),
        800: (150, -104, -500, -693, -912),
        820: (125, -123, -512, -701, -915),
        840: (99, -130, -524, -708, -917),
        860: (73, -149, -535, -715, -919),
        880: (49, -169, -547, -716, -921),
        900: (23, -188, -558, -723, -923),
        920: (3, -208, -572, -730, -925),
        940: (-22, -228, -583, -737, -927),
        960: (-47, -236, -594, -739, -929),
        980: (-70, -256, -605, -746, -931),
        1000: (-95, -277, -618, -753, -933),
    },
}
DISKS={
    (1,0): {
        6: ((-443, 792), (-837, 588), (-1056, 216)),
        7: ((-397, 801), (-761, 657), (-1002, 366), (-1088, 0)),
    },
    (2,0): {
        6: ((-104, 911), (-686, 756), (-1029, 288)),
        7: ((-40, 905), (-568, 829), (-942, 484), (-1077, 0)),
    },
    (1,1): {
        6: ((-125, 992), (-664, 748), (-961, 275)),
        7: ((-59, 998), (-557, 831), (-885, 466), (-1000, 0)),
    },
}

CORES=((1,0),(2,0),(1,1))
for core in CORES[1:]:
    CORE_RAYS[core][149]=CORE_RAYS[core][150]
BANDS=((110,620,Q(1,8)),(620,662,Q(2,19)),(662,900,Q(2,25)),(900,1000,Q(1,20)))

@lru_cache(maxsize=None)
def core_weights(order, core):
    one=tuple(f/scalar(comb(order,k)) for k,f in enumerate(marked_row(order)))
    forward,backward=core
    return tuple(a**forward*b**backward for a,b in zip(one,one[::-1],strict=True))


def coefficient_model(order,u,env,core):
    rho, variables=env
    mid=(u[0]+u[1])/2
    reference=tuple(a*b for a,b in zip(nominal(order,mid),core_weights(order,core),strict=True))
    U=-log(interval(*u))
    mean=log(interval(u[0]/mid,u[1]/mid))
    if order==7:
        mean=ANCHOR_RATIO*mean+(ANCHOR_RATIO-rho)*U
        U=rho*U
    pieces=tuple((p*U,vector) for vector,p in variables)
    pieces+=((mean,tuple(rational(e) for e in exponents(order))),)
    remainder=[]
    for k,a in enumerate(reference):
        radius=ZERO
        for parameter,vector in pieces:
            term=parameter*vector[k]
            radius=radius+I(0,max(-term.lo,term.hi))
        require(radius.hi<GRID,("exponential remainder",order,u))
        error=scalar(1,2)*radius.square()/(ONE-radius)
        remainder.append(a*(ONE+error))
    return ((ONE,tuple(remainder)),)+tuple(
        (parameter,tuple(a*v for a,v in zip(reference,vector,strict=True)))
        for parameter,vector in pieces
    )


def ray_cosines(core,left,right,index,u):
    rays=CORE_RAYS[core]
    a,b=rays[left][index],rays[right][index]
    first=Q(a,1000)+Q(b-a,right-left)*(u[0]-Q(left,1000))
    last=Q(a,1000)+Q(b-a,right-left)*(u[1]-Q(left,1000))
    return interval(min(first,last),max(first,last))


def verify(core,order,index,left,right,t):
    env=envelope(t,order)
    pending=[((Q(left,1000),Q(right,1000)),0)]
    count=0
    while pending:
        u,depth=pending.pop()
        used=radial(order,ray_cosines(core,left,right,index,u),coefficient_model(order,u,env,core))
        if used is not None:
            count+=used
            continue
        require(depth<15,("unresolved marked ray",core,order,index,left,right,t,[float(x) for x in u]))
        middle=(u[0]+u[1])/2
        pending.extend((((u[0],middle),depth+1),((middle,u[1]),depth+1)))
    return count


def check_anchor(core):
    rays=CORE_RAYS[core]
    for row in rays.values():
        require(all(1000>a>b>-1000 for a,b in zip(row,row[1:])),"strict ray order")
    radius=scalar(1,500)
    sector_order={6:(1,3,5),7:(0,2,4,5)}
    for order in (6,7):
        weights=core_weights(order,core)
        coefficients=tuple(scalar(comb(order,k))*f for k,f in enumerate(weights))
        points=DISKS[core][order]
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
            require((radius*terms[1]).lo>remainder.hi,("Rouche disk",core,order,a,b))
            for k,c in enumerate(rays[1000]):
                cos=scalar(c,1000)
                sine=(ONE-cos.square()).sqrt()
                side=cos*y-sine*x
                require(side.lo>radius.hi if k<sector else side.hi<(-radius).lo,
                        ("sector disk",core,order,a,b,k))
            if b:
                require(y.lo>radius.hi,"upper disk")
        full=points+tuple((a,-b) for a,b in points if b)
        require(len(full)==order,"all roots isolated")
        for i,(a,b) in enumerate(full):
            for x,y in full[:i]:
                require((a-x)**2+(b-y)**2>16,"disjoint disks")


def verify_constants():
    require(Q(11,100)**16<Q(7,12)**65,"singleton curvature endpoint")
    require(Q(149,1000)**48<Q(7,12)**169,"double curvature endpoint")
    require(Q(11,10)**7<2,"annulus")
    require(Q(62,100)**2>Q(13*14,23*24), "first residual gap-two cap")
    require(Q(662,1000)**2>Q(19*20,29*30), "second residual gap-two cap")
    for u,t in ((Q(62,100),Q(2,19)),(Q(662,1000),Q(2,25))):
        require(((1+t)/(1+6*t))**3<u*u,"larger residual gap caps")
    require((Q(9,10)**-2-1)/5<Q(1,20),"outer loss cap")


def main():
    verify_constants()
    for number,core in enumerate(CORES):
        if len(argv)>1 and number!=int(argv[1]):
            continue
        check_anchor(core)
        start=110 if number==0 else 149
        for low,high,cap in BANDS:
            left0=max(start,low)
            if len(argv)>2 and low!=int(argv[2]):
                continue
            knots=sorted(k for k in CORE_RAYS[core] if left0<=k<=high)
            for left,right in zip(knots,knots[1:]):
                cells=[verify(core,j,ray,left,right,cap) for j in (6,7) for ray in range(5)]
                print(core,left,right,cap,cells,flush=True)


if __name__=="__main__":
    main()
