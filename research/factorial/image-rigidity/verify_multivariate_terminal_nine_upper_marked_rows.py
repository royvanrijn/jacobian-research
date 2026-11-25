#!/usr/bin/env python3
"""Exact terminal-nine marked unit rows in the upper moderate band."""

from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial
from sys import argv

from verify_multivariate_terminal_seven_upper_bridge import bisect_u, cayley_matrix
from verify_multivariate_terminal_eight_late_bridge import marked_disk_check
from verify_multivariate_terminal_nine_private_linear import ratio_bounds, shape
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, interval, log, rational, require, scalar,
)
from verify_multivariate_terminal_six_marked_row_tail import projected_interval


ORDERS = (8, 9)
REFERENCE_RATIO = Q(64, 63)
SCALE = 100000
COORDINATE_P = Q(1,8)
RADII = tuple(Q(d,1000) for d in (13,14,15,16,18,20,24,28,32))
OUTER_CUTOFFS = {3: Q(21,25), 4: Q(22,25)}
RESIDUAL_CUTOFFS = {3: Q(19,20), 4: Q(121,125)}

MARKED = {
    3: {
        950: ((140242, 118022, 94143, 76320, 57060, 41104, 23662), (((6927, 138025), (552, 90334), (-2458, 51811), (-3731, 16927)), ((7590, 142310), (1290, 97951), (-1886, 62303), (-3471, 30389), (-3957, 0)))),
        952: ((139296, 117255, 93563, 75860, 56726, 40866, 23526), (((7000, 137096), (558, 89780), (-2489, 51509), (-3780, 16830)), ((7668, 141345), (1304, 97344), (-1910, 61936), (-3515, 30215), (-4008, 0)))),
        954: ((138349, 116486, 92981, 75398, 56390, 40627, 23391), (((7074, 136165), (564, 89224), (-2521, 51205), (-3829, 16733)), ((7749, 140378), (1319, 96735), (-1934, 61568), (-3561, 30040), (-4060, 0)))),
        956: ((137400, 115715, 92396, 74934, 56053, 40387, 23254), (((7150, 135233), (571, 88666), (-2553, 50900), (-3879, 16636)), ((7831, 139410), (1334, 96124), (-1959, 61198), (-3608, 29865), (-4114, 0)))),
        958: ((136449, 114942, 91810, 74468, 55714, 40146, 23117), (((7228, 134299), (577, 88106), (-2587, 50594), (-3931, 16538)), ((7915, 138440), (1349, 95511), (-1985, 60827), (-3655, 29688), (-4169, 0)))),
        960: ((135497, 114167, 91221, 74001, 55374, 39903, 22979), (((7308, 133363), (583, 87544), (-2621, 50287), (-3984, 16439)), ((8001, 137469), (1365, 94896), (-2011, 60454), (-3704, 29511), (-4225, 0)))),
        962: ((134543, 113390, 90631, 73532, 55032, 39660, 22841), (((7389, 132425), (590, 86980), (-2656, 49977), (-4037, 16340)), ((8089, 136496), (1381, 94279), (-2037, 60080), (-3754, 29333), (-4282, 0)))),
        964: ((133587, 112611, 90038, 73061, 54689, 39415, 22702), (((7473, 131485), (597, 86414), (-2692, 49667), (-4093, 16241)), ((8180, 135521), (1397, 93660), (-2065, 59703), (-3805, 29154), (-4340, 0)))),
        966: ((132629, 111830, 89443, 72587, 54344, 39169, 22562), (((7559, 130544), (603, 85846), (-2729, 49355), (-4149, 16141)), ((8272, 134544), (1414, 93038), (-2093, 59325), (-3857, 28974), (-4400, 0)))),
        968: ((131669, 111047, 88846, 72112, 53997, 38922, 22421), (((7646, 129600), (610, 85276), (-2766, 49041), (-4207, 16041)), ((8367, 133565), (1431, 92414), (-2122, 58946), (-3911, 28793), (-4462, 0)))),
        970: ((130708, 110262, 88247, 71635, 53649, 38674, 22280), (((7736, 128654), (618, 84703), (-2805, 48726), (-4266, 15940)), ((8464, 132584), (1449, 91788), (-2151, 58564), (-3965, 28611), (-4524, 0)))),
        972: ((129744, 109474, 87645, 71156, 53299, 38424, 22138), (((7828, 127706), (625, 84127), (-2844, 48409), (-4327, 15838)), ((8563, 131601), (1467, 91159), (-2181, 58180), (-4021, 28428), (-4589, 0)))),
        974: ((128778, 108684, 87040, 70674, 52947, 38173, 21995), (((7922, 126756), (632, 83550), (-2885, 48091), (-4389, 15736)), ((8665, 130616), (1485, 90528), (-2212, 57795), (-4079, 28244), (-4655, 0)))),
        976: ((127810, 107892, 86433, 70191, 52593, 37920, 21851), (((8019, 125804), (640, 82970), (-2926, 47771), (-4453, 15633)), ((8770, 129629), (1504, 89894), (-2244, 57408), (-4138, 28059), (-4722, 0)))),
        978: ((126840, 107098, 85824, 69704, 52238, 37666, 21706), (((8118, 124849), (648, 82387), (-2968, 47449), (-4518, 15529)), ((8877, 128640), (1523, 89258), (-2276, 57018), (-4198, 27873), (-4792, 0)))),
        980: ((125868, 106300, 85211, 69216, 51880, 37411, 21561), (((8220, 123892), (656, 81801), (-3012, 47125), (-4585, 15425)), ((8987, 127648), (1543, 88619), (-2309, 56626), (-4260, 27686), (-4863, 0)))),
        982: ((124893, 105501, 84596, 68725, 51520, 37154, 21414), (((8325, 122932), (664, 81213), (-3057, 46799), (-4654, 15321)), ((9100, 126655), (1563, 87977), (-2343, 56233), (-4324, 27498), (-4936, 0)))),
        984: ((123916, 104698, 83979, 68232, 51159, 36896, 21267), (((8432, 121970), (673, 80622), (-3103, 46471), (-4725, 15215)), ((9216, 125659), (1584, 87332), (-2378, 55837), (-4389, 27308), (-5011, 0)))),
        986: ((122937, 103893, 83358, 67736, 50795, 36636, 21119), (((8542, 121006), (682, 80029), (-3150, 46142), (-4797, 15109)), ((9335, 124660), (1606, 86684), (-2414, 55438), (-4456, 27117), (-5087, 0)))),
        988: ((121955, 103085, 82734, 67237, 50429, 36374, 20969), (((8655, 120039), (691, 79432), (-3198, 45810), (-4872, 15002)), ((9457, 123659), (1628, 86034), (-2451, 55038), (-4525, 26925), (-5166, 0)))),
        990: ((120971, 102275, 82108, 66736, 50060, 36110, 20819), (((8771, 119069), (700, 78832), (-3248, 45476), (-4948, 14895)), ((9582, 122656), (1650, 85380), (-2489, 54634), (-4596, 26732), (-5247, 0)))),
        992: ((119984, 101461, 81478, 66231, 49690, 35845, 20667), (((8891, 118097), (709, 78229), (-3299, 45140), (-5027, 14786)), ((9711, 121650), (1674, 84723), (-2528, 54229), (-4668, 26537), (-5331, 0)))),
        994: ((118995, 100645, 80845, 65724, 49316, 35578, 20515), (((9013, 117122), (719, 77623), (-3351, 44802), (-5108, 14677)), ((9844, 120642), (1697, 84063), (-2568, 53820), (-4743, 26341), (-5416, 0)))),
        996: ((118003, 99825, 80208, 65214, 48941, 35309, 20361), (((9139, 116145), (729, 77014), (-3405, 44461), (-5191, 14567)), ((9980, 119631), (1722, 83399), (-2609, 53409), (-4819, 26143), (-5505, 0)))),
        998: ((117008, 99002, 79568, 64701, 48562, 35038, 20206), (((9269, 115164), (740, 76401), (-3460, 44118), (-5277, 14456)), ((10120, 118617), (1748, 82732), (-2651, 52996), (-4898, 25944), (-5595, 0)))),
        1000: ((116011, 98177, 78925, 64184, 48181, 34765, 20050), (((9403, 114181), (751, 75785), (-3517, 43772), (-5365, 14344)), ((10264, 117601), (1774, 82061), (-2694, 52579), (-4980, 25743), (-5689, 0)))),
    },
    4: {
        968: ((118600, 100210, 80375, 65309, 48967, 35319, 20356), (((5498, 116816), (445, 77182), (-2045, 44495), (-3123, 14570)), ((6001, 120301), (1042, 83566), (-1565, 53435), (-2898, 26138), (-3312, 0)))),
        970: ((117565, 99359, 79718, 64784, 48581, 35043, 20199), (((5578, 115798), (451, 76554), (-2080, 44146), (-3176, 14458)), ((6088, 119247), (1058, 82881), (-1591, 53013), (-2947, 25936), (-3368, 0)))),
        972: ((116527, 98505, 79058, 64256, 48194, 34766, 20041), (((5661, 114776), (458, 75922), (-2115, 43795), (-3231, 14345)), ((6178, 118190), (1074, 82193), (-1618, 52589), (-2998, 25732), (-3426, 0)))),
        974: ((115484, 97646, 78394, 63725, 47803, 34487, 19882), (((5747, 113750), (465, 75287), (-2152, 43441), (-3287, 14231)), ((6270, 117128), (1091, 81500), (-1646, 52162), (-3050, 25527), (-3486, 0)))),
        976: ((114437, 96783, 77727, 63191, 47411, 34206, 19721), (((5835, 112720), (472, 74648), (-2189, 43085), (-3346, 14116)), ((6366, 116063), (1108, 80804), (-1674, 51732), (-3103, 25321), (-3548, 0)))),
        978: ((113387, 95917, 77055, 62653, 47015, 33923, 19559), (((5926, 111686), (479, 74006), (-2228, 42727), (-3406, 14000)), ((6464, 114993), (1126, 80104), (-1704, 51299), (-3159, 25113), (-3611, 0)))),
        980: ((112331, 95046, 76380, 62112, 46617, 33638, 19397), (((6020, 110648), (487, 73359), (-2268, 42366), (-3467, 13884)), ((6566, 113919), (1144, 79399), (-1735, 50863), (-3216, 24904), (-3677, 0)))),
        982: ((111272, 94171, 75700, 61568, 46216, 33351, 19232), (((6118, 109605), (495, 72709), (-2310, 42002), (-3531, 13767)), ((6671, 112840), (1163, 78691), (-1766, 50425), (-3275, 24693), (-3745, 0)))),
        984: ((110208, 93291, 75017, 61019, 45812, 33061, 19067), (((6218, 108557), (503, 72054), (-2353, 41636), (-3598, 13648)), ((6780, 111758), (1183, 77978), (-1799, 49983), (-3336, 24480), (-3815, 0)))),
        986: ((109140, 92407, 74328, 60467, 45405, 32770, 18900), (((6322, 107505), (511, 71395), (-2397, 41267), (-3666, 13529)), ((6892, 110670), (1203, 77260), (-1833, 49537), (-3399, 24266), (-3887, 0)))),
        988: ((108066, 91518, 73636, 59911, 44994, 32475, 18732), (((6430, 106448), (520, 70731), (-2443, 40895), (-3737, 13409)), ((7008, 109578), (1224, 76538), (-1868, 49089), (-3465, 24050), (-3962, 0)))),
        990: ((106988, 90624, 72938, 59351, 44581, 32179, 18562), (((6541, 105386), (529, 70063), (-2490, 40520), (-3810, 13287)), ((7128, 108481), (1246, 75812), (-1904, 48636), (-3532, 23832), (-4040, 0)))),
        992: ((105905, 89726, 72236, 58787, 44163, 31880, 18391), (((6656, 104320), (538, 69391), (-2539, 40141), (-3886, 13165)), ((7253, 107379), (1268, 75080), (-1941, 48180), (-3602, 23612), (-4120, 0)))),
        994: ((104817, 88822, 71529, 58218, 43743, 31578, 18219), (((6776, 103247), (548, 68713), (-2590, 39760), (-3965, 13041)), ((7382, 106272), (1292, 74343), (-1980, 47720), (-3675, 23390), (-4204, 0)))),
        996: ((103724, 87913, 70816, 57645, 43319, 31274, 18044), (((6900, 102170), (558, 68030), (-2643, 39375), (-4046, 12916)), ((7516, 105159), (1316, 73601), (-2020, 47257), (-3750, 23165), (-4290, 0)))),
        998: ((102626, 86999, 70099, 57067, 42890, 30966, 17868), (((7028, 101087), (568, 67342), (-2698, 38986), (-4131, 12790)), ((7654, 104042), (1341, 72853), (-2062, 46789), (-3828, 22939), (-4380, 0)))),
        1000: ((101522, 86079, 69375, 56484, 42458, 30656, 17690), (((7161, 99999), (579, 66648), (-2755, 38594), (-4219, 12663)), ((7798, 102919), (1368, 72100), (-2105, 46316), (-3909, 22711), (-4473, 0)))),
    },
}



@lru_cache(maxsize=8192)
def precise_power(value,exponent):
    if value in (0,1):
        return rational(value)
    d=2**52
    a,b=exponent.numerator,exponent.denominator
    right=d**b*value.numerator**a
    factor=value.denominator**a
    middle=int(float(value)**float(exponent)*d)
    lo,hi=max(0,middle-4),min(d,middle+5)
    if not lo**b*factor<=right<=hi**b*factor:
        lo,hi=0,d
    while hi-lo>1:
        middle=(lo+hi)//2
        if middle**b*factor<right:
            lo=middle
        else:
            hi=middle
    return interval(Q(lo,d),Q(hi,d))


@lru_cache(maxsize=8)
def marked_core(order,base):
    denominator=10**22
    result=[]
    for k in range(order+1):
        if k in (0,order):
            result.append((Q(1),Q(1)))
            continue
        target=Q((comb(order,k)*factorial(order*base+k)*denominator)**order,
                 factorial(order*base)**(order-k)*factorial(order*(base+1))**k)
        lo,hi=0,(1+comb(order,order//2))*denominator
        require(hi**order>=target,"unit core root bracket")
        while hi-lo>1:
            mid=(lo+hi)//2
            if mid**order<target:
                lo=mid
            else:
                hi=mid
        result.append((Q(lo,denominator),Q(hi,denominator)))
    return tuple(result)


@lru_cache(maxsize=8192)
def model(order,base,v):
    c8=(v[0]**-2-1)/7
    cap=c8 if order==8 else 8*c8/(9-c8)
    mid=sum(v)/2
    a=tuple(interval(*q) for q in marked_core(order,base))
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

    mat=matrix(order)
    lam=rational(COORDINATE_P).sqrt()
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


@lru_cache(maxsize=2)
def matrix(order):
    return tuple(tuple(map(scalar,row)) for row in cayley_matrix(order))


def table_value(v,left,right,a,b):
    return Q(a,SCALE)+Q(b-a,SCALE)*(1000*v-left)/(right-left)


def sector_check(centre,rays,sector,radius):
    x,y=(interval(q-radius,q+radius) for q in centre)
    if y.lo<=0:
        return False
    p=rational(COORDINATE_P)
    for i,b in enumerate(rays):
        val=(y-b)*(ONE+p*b*y)+p*b*x.square()
        if i<sector:
            if val.hi>=0:
                return False
        elif val.lo<=0:
            return False
    return True


def verify_band(base,order,left,right):
    table=MARKED[base]
    rays1,centres1=table[left]
    rays2,centres2=table[right]
    pending=[((Q(left,1000),Q(right,1000)),0)]
    sectors=(1,3,5,7) if order==8 else (0,2,4,6)
    count=0
    while pending:
        v,depth=pending.pop()
        midpoint=sum(v)/2
        ray_box=tuple(interval(*sorted(table_value(end,left,right,a,b) for end in v))
                      for a,b in zip(rays1,rays2,strict=True))
        centres=tuple(tuple(table_value(midpoint,left,right,a,b) for a,b in zip(w,z,strict=True))
                      for w,z in zip(centres1[order-8],centres2[order-8],strict=True))
        box=model(order,base,v)
        chosen=[]
        for i,(x,y) in enumerate(centres):
            sector=sectors[i] if i<len(sectors) else None
            radius=next((r for r in RADII if x*x+y*y<(2-r)**2
                         and (sector_check((x,y),ray_box,sector,r) if sector is not None else y==0)
                         and marked_disk_check(order,(x,y),box,r)),None)
            if radius is None:
                break
            chosen.append((x,y,radius))
        else:
            full=chosen+[(x,-y,r) for x,y,r in chosen if y>0]
            require(len(full)==order,"marked candidate count")
            if all((a-x)**2+(b-y)**2>(t+r)**2
                   for i,(x,y,r) in enumerate(full) for a,b,t in full[:i]):
                count+=1
                continue
        require(depth<14,("unresolved marked disk",base,order,v,len(chosen)))
        pending.extend((half,depth+1) for half in bisect_u(v))
    return count


def main():
    for base,table in MARKED.items():
        u0_square=Q(8*base+1,8*base+8)
        require(OUTER_CUTOFFS[base]>u0_square,"unique marked row")
        require(OUTER_CUTOFFS[base]**2>u0_square*RESIDUAL_CUTOFFS[base]**2,"residual lower cutoff")
        if len(argv)>1 and int(argv[1])!=base:
            continue
        knots=tuple(table)
        require(Q(knots[0],1000)==RESIDUAL_CUTOFFS[base] and knots[-1]==1000,"marked coverage")
        for rays,_ in table.values():
            z=(2*SCALE,*rays,0)
            require(all(a>b for a,b in zip(z,z[1:])),"marked ray order")
        for order in ORDERS:
            total=0
            for a,b in zip(knots,knots[1:]):
                if len(argv)>2 and int(argv[2])!=order:
                    continue
                count=verify_band(base,order,a,b)
                total+=count
                print(base,order,a,b,count,flush=True)
            print("Marked",base,order,total,flush=True)
    require(Q(25,32)*Q(33,40)<Q(22,25)**2,"marked alternatives do not coexist")

if __name__=="__main__":
    main()
