#!/usr/bin/env python3
"""Terminal-ten retained complements for the fully shared .52–.76 bands."""

from ctypes import CDLL
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, prod
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_nine_upper_marked_rows import precise_power
from verify_multivariate_terminal_six_intermediate_cayley import GRID, ONE, ZERO, interval, require
from verify_multivariate_terminal_six_private_angles import multiply_interval
from verify_multivariate_terminal_six_retained_rows import core_powers
from verify_multivariate_terminal_ten_shared_light import (
    ANNULUS, coefficient_model, exact_cap, independent_coefficients, nominal_exponent,
    normalized_native, replay_python, verify_anchor_order,
)

SCALE = 10**6
KNOT_SCALE = 10000

BANDS=((Q(52,100),Q(62,100),Q(99,1000)),
       (Q(62,100),Q(70,100),Q(3,40)),
       (Q(70,100),Q(76,100),Q(11,200)))
LIMIT=Q(1,8)


@lru_cache(maxsize=None)
def row_square(base,gap):
    return prod(Q(9*base+h,9*base+8*gap+h) for h in range(1,gap+1))

def candidates(lower,cap):
    gap=1
    threshold=1+8*cap/(1+cap)
    found=[]
    while threshold**(-gap)>=lower**2:
        base=1
        while 9*base+1<gap/cap:
            if Q(gap,9*base+1)<=LIMIT:
                value=row_square(base,gap)
                if value>=lower**2:
                    found.extend(((base,gap,sign),value) for sign in (-1,1))
            base+=1
        gap+=1
    return tuple(sorted(found))

def enumerate_cores(lower,cap):
    rows=candidates(lower,cap)
    result={}
    def visit(i,selected,squared):
        if i==len(rows):
            if selected:
                reversal=tuple(sorted((a,d,-o) for a,d,o in selected))
                if selected<=reversal:
                    result[selected]=squared
            return
        row,v2=rows[i]
        require(Q(0)<v2<Q(1),"retained row mean")
        while squared>=lower**2:
            visit(i+1,selected,squared)
            selected+=(row,)
            squared*=v2
    visit(0,(),Q(1))
    return tuple(sorted(result.items()))

def bracket_tail(lower,upper,core_square):
    l,h=0,KNOT_SCALE
    while l+1<h:
        m=(l+h)//2
        if Q(m,KNOT_SCALE)**2*core_square<=lower**2:
            l=m
        else:
            h=m
    first=Q(l,KNOT_SCALE)
    if core_square<=upper**2:
        return first,Q(1)
    l,h=0,KNOT_SCALE
    while l+1<h:
        m=(l+h)//2
        if Q(m,KNOT_SCALE)**2*core_square<upper**2:
            l=m
        else:
            h=m
    return first,Q(h,KNOT_SCALE)


CASES = (
    (0, ((1, 1, -1),), True, {
        6976: (149988, -72064, -310535, -483590, -661874, -781757, -892974, -954411),
        7056: (136931, -83193, -319468, -490475, -666560, -784831, -894510, -955071),
        7136: (123826, -94337, -328384, -497338, -671224, -787888, -896036, -955727),
        7216: (110674, -105495, -337283, -504180, -675867, -790930, -897554, -956378),
        7296: (97476, -116666, -346165, -511001, -680488, -793956, -899063, -957026),
        7376: (84233, -127850, -355029, -517802, -685089, -796967, -900562, -957670),
        7456: (70947, -139046, -363877, -524583, -689669, -799963, -902054, -958310),
        7536: (57616, -150255, -372709, -531344, -694230, -802943, -903536, -958946),
        7616: (44244, -161476, -381523, -538085, -698770, -805909, -905011, -959579),
        7696: (30830, -172708, -390321, -544807, -703290, -808861, -906477, -960207),
        7776: (17377, -183951, -399103, -551510, -707792, -811799, -907935, -960833),
        7856: (3884, -195204, -407868, -558193, -712274, -814722, -909385, -961454),
        7936: (-9647, -206468, -416617, -564858, -716738, -817632, -910828, -962072),
        8016: (-23215, -217740, -425350, -571504, -721184, -820529, -912263, -962687),
        8096: (-36819, -229022, -434068, -578132, -725611, -823413, -913690, -963299),
        8176: (-50456, -240312, -442769, -584742, -730022, -826284, -915111, -963907),
        8256: (-64127, -251610, -451454, -591335, -734415, -829142, -916524, -964512),
        8320: (-75086, -260653, -458392, -596597, -737917, -831420, -917650, -964994),
    }, {9: ((13860, 104387), (-34897, 94840), (-69130, 70153), (-89724, 37042), (-96612, 0)), 10: ((17791, 104167), (-28008, 97634), (-61552, 77515), (-83920, 49434), (-95141, 16951))}),
    (0, ((1, 1, -1), (1, 1, -1)), False, {
        9360: (151005, -71573, -310349, -483696, -662003, -781961, -893057, -954483),
        9440: (141478, -79777, -317022, -488869, -665553, -784297, -894230, -954988),
        9520: (131937, -87984, -323688, -494033, -669094, -786627, -895398, -955491),
        9600: (122382, -96193, -330346, -499189, -672626, -788951, -896564, -955992),
        9680: (112814, -104404, -336996, -504336, -676150, -791268, -897726, -956492),
        9760: (103235, -112618, -343639, -509475, -679666, -793580, -898884, -956990),
        9840: (93644, -120832, -350274, -514605, -683175, -795886, -900040, -957487),
        9920: (84044, -129047, -356902, -519728, -686676, -798187, -901193, -957983),
        10000: (74435, -137263, -363523, -524843, -690170, -800483, -902342, -958478),
    }, {9: ((14638, 110015), (-35250, 95695), (-68084, 69019), (-87158, 35946), (-93437, 0)), 10: ((19012, 110134), (-28332, 99035), (-60861, 76670), (-81770, 48161), (-92055, 16397))}),
    (0, ((1, 1, -1), (1, 1, 1)), False, {
        9360: (135226, -76874, -305488, -476847, -654231, -776046, -889489, -952882),
        9440: (125062, -85305, -311998, -481801, -657540, -778204, -890556, -953340),
        9520: (114857, -93747, -318494, -486737, -660832, -780349, -891617, -953795),
        9600: (104611, -102201, -324974, -491657, -664107, -782481, -892670, -954247),
        9680: (94325, -110668, -331440, -496559, -667364, -784601, -893717, -954695),
        9760: (83998, -119146, -337890, -501443, -670605, -786709, -894757, -955141),
        9840: (73629, -127636, -344325, -506310, -673828, -788804, -895789, -955583),
        9920: (63218, -136137, -350744, -511160, -677034, -790887, -896815, -956023),
        10000: (52765, -144651, -357147, -515991, -680223, -792957, -897833, -956459),
    }, {9: ((11703, 99313), (-34020, 94035), (-69497, 71904), (-92185, 38754), (-100000, 0)), 10: ((15343, 98816), (-27077, 96264), (-61349, 78970), (-85712, 51511), (-98391, 17867))}),
    (0, ((2, 2, -1),), False, {
        9468: (150389, -71983, -310717, -483873, -662189, -782009, -893127, -954481),
        9548: (140982, -80087, -317318, -488992, -665704, -784322, -894289, -954981),
        9628: (131561, -88194, -323912, -494101, -669210, -786630, -895447, -955480),
        9708: (122129, -96304, -330498, -499203, -672708, -788931, -896602, -955977),
        9788: (112685, -104414, -337076, -504296, -676199, -791227, -897754, -956472),
        9868: (103231, -112526, -343647, -509381, -679682, -793517, -898902, -956966),
        9948: (93767, -120638, -350211, -514459, -683158, -795802, -900048, -957459),
        10000: (87611, -125912, -354474, -517755, -685414, -797285, -900791, -957779),
    }, {9: ((14699, 110473), (-35269, 95771), (-67993, 68940), (-86953, 35867), (-93187, 0)), 10: ((18955, 110626), (-28447, 99121), (-60854, 76573), (-81623, 48040), (-91823, 16346))}),
    (1, ((1, 1, -1),), True, {
        8318: (-74744, -260370, -458175, -596432, -737808, -831349, -917615, -964979),
        8398: (-88469, -271680, -466833, -602994, -742172, -834186, -919016, -965579),
        8478: (-102223, -282996, -475476, -609540, -746520, -837012, -920411, -966176),
        8558: (-116004, -294316, -484104, -616069, -750852, -839826, -921800, -966770),
        8638: (-129810, -305641, -492717, -622582, -755169, -842629, -923183, -967361),
        8718: (-143638, -316969, -501315, -629080, -759472, -845422, -924559, -967950),
        8798: (-157487, -328299, -509899, -635562, -763761, -848205, -925931, -968536),
        8878: (-171355, -339631, -518469, -642030, -768037, -850979, -927297, -969120),
        8958: (-185238, -350963, -527025, -648484, -772300, -853743, -928658, -969702),
        9038: (-199133, -362294, -535568, -654925, -776551, -856499, -930014, -970282),
        9118: (-213038, -373622, -544098, -661353, -780791, -859247, -931366, -970859),
        9198: (-226950, -384947, -552614, -667768, -785020, -861987, -932714, -971435),
        9278: (-240864, -396266, -561119, -674172, -789240, -864721, -934058, -972009),
        9358: (-254776, -407579, -569612, -680565, -793451, -867448, -935400, -972582),
        9392: (-260688, -412384, -573218, -683279, -795238, -868606, -935969, -972825),
    }, {9: ((-9823, 107531), (-49406, 88708), (-75209, 62453), (-90137, 32156), (-95044, 0)), 10: ((-6352, 108362), (-43968, 92375), (-69551, 69712), (-85918, 43213), (-93952, 14627))}),
    (2, ((1, 1, -1),), True, {
        9391: (-260514, -412243, -573112, -683200, -795186, -868572, -935952, -972818),
        9471: (-274417, -423542, -581589, -689579, -799386, -871292, -937289, -973389),
        9551: (-288307, -434830, -590055, -695951, -803581, -874008, -938624, -973960),
        9631: (-302180, -446103, -598512, -702314, -807771, -876721, -939958, -974529),
        9711: (-316029, -457360, -606959, -708672, -811957, -879432, -941291, -975098),
        9791: (-329848, -468598, -615397, -715024, -816141, -882142, -942624, -975667),
        9871: (-343631, -479813, -623828, -721372, -820326, -884853, -943957, -976237),
        9951: (-357370, -491002, -632250, -727718, -824511, -887564, -945291, -976807),
        10000: (-365760, -497842, -637405, -731603, -827076, -889227, -946110, -977156),
    }, {9: ((-30597, 107223), (-60874, 81984), (-79405, 55401), (-89763, 27920), (-93115, 0)), 10: ((-27706, 108870), (-56768, 86183), (-75339, 62377), (-86797, 37748), (-92312, 12636))}),
)

@lru_cache(maxsize=None)
def core_factor(order,rows):
    squared=core_powers(order,rows)
    out=[]
    denominator=2**55
    for val in squared:
        lo,hi=0,denominator
        rhs=val.numerator*denominator**order
        norm=val.denominator
        require(hi**order*norm>=rhs,"marked root bracket")
        while lo+1<hi:
            m=(lo+hi)//2
            if m**order*norm<rhs:
                lo=m
            else:
                hi=m
        out.append((Q(lo,denominator),Q(hi,denominator)))
    out[0]=out[-1]=(Q(1),Q(1))
    return tuple(out)

@lru_cache(maxsize=8192)
def model(order,rows,cap,means,signed,extended):
    factors=core_factor(order,rows)
    if means==(Q(1),Q(1)):
        core=tuple(interval(lo*comb(order,k),hi*comb(order,k))
                   for k,(lo,hi) in enumerate(factors))
        return ((ONE,core),(ZERO,tuple(ZERO for _ in core)))
    contents=coefficient_model(order,cap,means,signed,extended)
    if not contents:
        return contents
    c=tuple(interval(lo,hi) for lo,hi in factors)
    return tuple((p,tuple(a*x for a,x in zip(vector,c,strict=True))) for p,vector in contents)

def anchor(order,rows,cap,lower,upper,rays,centers):
    factors=core_factor(order,rows)
    independent=independent_coefficients(order,cap,(Q(0),upper))
    box=tuple(multiply_interval(a,b) for a,b in zip(factors,independent,strict=True))
    small,large=ANNULUS
    require(sum(cap*small**k for k,(_,cap) in enumerate(box) if k)<1,("inner radius",rows,order))
    require(sum(cap*large**k for k,(_,cap) in enumerate(box[:-1]))<large**order,
            ("outer radius",rows,order))
    rows=[]
    for k,(lo,hi) in enumerate(factors):
        coef=interval(comb(order,k)*lo,comb(order,k)*hi)*precise_power(lower,nominal_exponent(order,k))
        rows.append((Q(coef.lo,GRID),Q(coef.hi,GRID)))
    verify_anchor_order(order,rows,centers,tuple(Q(x,SCALE) for x in rays))

def verify_slab(order,rows,cap,table,left,right,extended,replay):
    pending=[((Q(left,KNOT_SCALE),Q(right,KNOT_SCALE)),(Q(-1),Q(1)),tuple(range(8)),0,0)]
    total=0
    while pending:
        means,p,rays,md,sd=pending.pop()
        polynomial=model(order,rows,cap,means,p,extended)
        if polynomial is None:
            continue
        missing=[]
        for index in rays:
            a,b=table[left][index],table[right][index]
            cosines=tuple(sorted(Q(a,SCALE)+Q(b-a,SCALE)*(KNOT_SCALE*v-left)/(right-left) for v in means))
            slope=Q(b-a,SCALE)*KNOT_SCALE*sum(means)/(2*(right-left))
            used=replay(order,cosines,polynomial,slope) if polynomial else None
            if used is None:
                missing.append(index)
            else:
                total+=used
        if not missing:
            continue
        require(max(md,sd)<17,("unresolved retained ray",order,rows,means,p,missing))
        if md<2*sd+2:
            m=sum(means)/2
            pending.extend((((means[0],m),p,missing,md+1,sd),((m,means[1]),p,missing,md+1,sd)))
        else:
            m=sum(p)/2
            pending.extend(((means,(p[0],m),missing,md,sd+1),(means,(m,p[1]),missing,md,sd+1)))
    return total

def verify(replay):
    boxes = tuple(dict(enumerate_cores(lower, cap)) for lower, _, cap in BANDS)
    keys = {(band, rows) for band, box in enumerate(boxes) for rows in box}
    require(len(CASES) == 6 and keys == {(band, rows) for band, rows, _, _, _ in CASES},
            ("retained list enumeration", boxes))
    for band, rows, extended, rays, centers in CASES:
        outer_lower, outer_upper, declared = BANDS[band]
        actual_lower, actual_upper = bracket_tail(outer_lower, outer_upper, boxes[band][rows])
        knots = tuple(sorted(rays))
        lower, upper = Q(knots[0], KNOT_SCALE), Q(knots[-1], KNOT_SCALE)
        require(lower <= actual_lower <= actual_upper <= upper,
                ("residual coverage", band, rows, actual_lower, actual_upper))
        for values in rays.values():
            edges = (SCALE, *values, -SCALE)
            require(len(values) == 8 and all(a > b for a, b in zip(edges, edges[1:])),
                    "ordered rays")
        for order in (9, 10):
            anchor(order, rows, exact_cap(declared, lower), lower, upper, rays[knots[0]],
                   centers[order])
        total = 0
        for left, right in zip(knots, knots[1:]):
            cap = exact_cap(declared, Q(left, KNOT_SCALE))
            total += sum(verify_slab(order, rows, cap, rays, left, right, extended, replay)
                         for order in (9, 10))
        print(f"Band {outer_lower}–{outer_upper}, retained {rows}: {total} radial cells.",
              flush=True)


def main():
    require(argv[1:] in ([], ["--native"]), "usage: [--native]")
    if "--native" not in argv:
        verify(replay_python)
        return
    from verify_multivariate_terminal_ten_shared_light import __file__ as light
    with TemporaryDirectory(prefix="sfc-ten-retained-") as directory:
        library = Path(directory)/"bernstein.so"
        source = Path(light).with_suffix(".cpp")
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source), "-o", str(library)),
            check=True)
        verify(normalized_native(CDLL(str(library))))


if __name__ == "__main__":
    main()
