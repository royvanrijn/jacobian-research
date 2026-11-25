#!/usr/bin/env python3
"""Exact terminal-nine separation on the shared .55–.62 band."""

from ctypes import CDLL
from fractions import Fraction as Q
from functools import lru_cache
from math import comb
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

from flint import ctx

from verify_multivariate_terminal_nine_private_rows import (
    ORDERS, REFERENCE_RATIO, outer_coefficient_box, outer_exponents,
    verify_private_anchor, verify_skew_convexity,
)
from verify_multivariate_terminal_nine_shared_matched_upper import (
    ANNULUS, RAY_SCALE, coefficient_model, exact_cap, native_replay, replay_python,
    verify as verify_light_bands,
)
from verify_multivariate_terminal_nine_retained_bands import enumerate_band, verify_band_core
from verify_multivariate_terminal_nine_upper_marked_rows import marked_core, precise_power
from verify_multivariate_terminal_six_intermediate_cayley import GRID, interval, require


KNOT_SCALE = 3000
DECLARED_CAP = Q(1, 10)
TAIL_CAP = Q(1, 17)
LOWER, UPPER = Q(11, 20), Q(31, 50)
MARKED = ((1, 1, -1),)
EXPECTED = {
    MARKED,
    ((2, 2, -1),),
    ((1, 1, -1), (1, 1, -1)),
    ((1, 1, -1), (1, 1, 1)),
}

RAYS = {
    2200: (57651, -171810, -418671, -588122, -759528, -864609, -953133),
    2210: (52188, -176300, -422056, -590580, -761012, -865456, -953430),
    2220: (46717, -180792, -425438, -593035, -762493, -866301, -953726),
    2230: (41240, -185285, -428816, -595486, -763972, -867145, -954021),
    2240: (35756, -189779, -432191, -597934, -765447, -867986, -954316),
    2250: (30265, -194276, -435562, -600378, -766919, -868826, -954610),
    2260: (24768, -198773, -438931, -602819, -768389, -869664, -954903),
    2270: (19264, -203272, -442296, -605257, -769855, -870499, -955196),
    2280: (13754, -207772, -445657, -607691, -771319, -871334, -955488),
    2290: (8237, -212273, -449016, -610121, -772780, -872166, -955779),
    2300: (2714, -216776, -452371, -612549, -774237, -872996, -956070),
    2310: (-2816, -221280, -455723, -614973, -775692, -873825, -956359),
    2320: (-8351, -225786, -459071, -617394, -777145, -874652, -956648),
    2330: (-13893, -230292, -462417, -619811, -778594, -875477, -956937),
    2340: (-19440, -234800, -465759, -622226, -780041, -876301, -957225),
    2350: (-24994, -239309, -469098, -624637, -781485, -877122, -957512),
    2360: (-30554, -243819, -472434, -627044, -782926, -877942, -957798),
    2370: (-36119, -248330, -475767, -629449, -784364, -878761, -958084),
    2380: (-41691, -252842, -479097, -631851, -785800, -879577, -958369),
    2390: (-47268, -257355, -482423, -634249, -787233, -880392, -958653),
    2400: (-52850, -261869, -485746, -636644, -788664, -881206, -958937),
    2410: (-58438, -266385, -489067, -639036, -790092, -882018, -959220),
    2420: (-64032, -270901, -492384, -641425, -791517, -882828, -959503),
    2430: (-69631, -275418, -495698, -643811, -792940, -883636, -959785),
    2440: (-75235, -279936, -499009, -646194, -794360, -884443, -960066),
    2450: (-80845, -284455, -502317, -648574, -795778, -885248, -960347),
    2460: (-86459, -288974, -505622, -650950, -797193, -886052, -960627),
    2470: (-92079, -293495, -508924, -653324, -798605, -886854, -960907),
    2480: (-97703, -298016, -512223, -655695, -800016, -887655, -961186),
}

ANCHORS = {8: ((342, 9485), (-4553, 8844), (-8185, 6127), (-10126, 2181)), 9: ((749, 9419), (-3753, 9146), (-7306, 7072), (-9594, 3832), (-10385, 0))}


@lru_cache(maxsize=None)
def multipliers(order):
    return tuple(interval(lo/comb(order,k),hi/comb(order,k))
                 for k,(lo,hi) in enumerate(marked_core(order,1)))


@lru_cache(maxsize=8192)
def coefficients(order, means, signed):
    model = coefficient_model(order, TAIL_CAP, means, signed)
    if not model:
        return model
    factors = multipliers(order)
    return tuple((parameter, tuple(a*b for a,b in zip(vector,factors,strict=True)))
                 for parameter,vector in model)


def verify_setup():
    verify_skew_convexity()
    lower, upper = Q(min(RAYS),KNOT_SCALE), Q(max(RAYS),KNOT_SCALE)
    require(Q(3,4)*lower == LOWER and Q(3,4)*upper == UPPER, "tail endpoints")
    require(exact_cap(DECLARED_CAP,lower) == TAIL_CAP, "residual reciprocal cap")
    for values in RAYS.values():
        sequence=(RAY_SCALE,*values,-RAY_SCALE)
        require(len(values)==7 and all(a>b for a,b in zip(sequence,sequence[1:])),
                "ordered rays")
    small,large=ANNULUS
    for order in ORDERS:
        core=marked_core(order,1)
        exponents=outer_exponents(order,TAIL_CAP,8*TAIL_CAP/(9-TAIL_CAP))
        box=outer_coefficient_box(core,(Q(0),upper),exponents)
        require(sum(cap*small**k for k,(_,cap) in enumerate(box) if k)<1, ("inner radius",order))
        require(sum(cap*large**k for k,(_,cap) in enumerate(box[:-1]))<large**order,
                ("outer radius",order))
        ratio=Q(1) if order==8 else REFERENCE_RATIO
        anchor=[]
        for k,coefficient in enumerate(core):
            e=ratio*Q(k*(order-k),order-1)
            value=interval(*coefficient)*precise_power(lower,e)
            anchor.append((Q(value.lo,GRID),Q(value.hi,GRID)))
        verify_private_anchor(order,anchor,ANCHORS[order],
                              tuple(Q(x,RAY_SCALE) for x in RAYS[min(RAYS)]))


def verify_band(order,left,right,replay):
    means=Q(left,KNOT_SCALE),Q(right,KNOT_SCALE)
    pending=[(means,(Q(-1),Q(1)),tuple(range(7)),0,0)]
    total=0
    while pending:
        u,p,active,mean_depth,skew_depth=pending.pop()
        model=coefficients(order,u,p)
        if model is None:
            continue
        missing=[]
        for ray in active:
            a,b=RAYS[left][ray],RAYS[right][ray]
            cosines=tuple(sorted(Q(a,RAY_SCALE)+Q(b-a,RAY_SCALE)*(KNOT_SCALE*v-left)/(right-left)
                                 for v in u))
            slope=Q(b-a,RAY_SCALE)*KNOT_SCALE*sum(u)/(2*(right-left))
            count=replay(order,cosines,model,slope) if model else None
            if count is None:
                missing.append(ray)
            else:
                total+=count
        if not missing:
            continue
        require(max(mean_depth,skew_depth)<16,("unresolved marked ray",order,u,p,missing))
        if mean_depth<2*skew_depth+2:
            center=sum(u)/2
            pending.extend((((u[0],center),p,missing,mean_depth+1,skew_depth),
                            ((center,u[1]),p,missing,mean_depth+1,skew_depth)))
        else:
            center=sum(p)/2
            pending.extend(((u,(p[0],center),missing,mean_depth,skew_depth+1),
                            (u,(center,p[1]),missing,mean_depth,skew_depth+1)))
    return total


def verify_marked(replay):
    verify_setup()
    total=0
    knots=tuple(sorted(RAYS))
    for left,right in zip(knots,knots[1:]):
        counts=tuple(verify_band(j,left,right,replay) for j in ORDERS)
        total+=sum(counts)
        print(f"{left}–{right} / {KNOT_SCALE}: {counts}",flush=True)
    print(f"Verified {total} radial cells for the marked unit row.")


def verify(replay):
    ctx.prec = 128
    verify_light_bands(replay, (55,))
    budget, lists = enumerate_band(DECLARED_CAP, LOWER)
    require(len(lists) == len(EXPECTED) and {rows for rows, _ in lists} == EXPECTED,
            ("retained list enumeration", lists))
    leaves = cells = 0
    for rows, remaining in lists:
        if rows == MARKED:
            continue
        count, radial = verify_band_core(rows, remaining, budget-remaining,
                                         DECLARED_CAP, LOWER, UPPER, True)
        leaves += count
        cells += radial
        print(f"Retained list {rows}: {count} outer cells, {radial} radial cells.", flush=True)
    print(f"Verified {len(lists)} lists, including the marked singleton.", flush=True)
    print(f"Other retained lists: {leaves} outer cells and {cells} radial cells.", flush=True)
    verify_marked(replay)


def main():
    require(argv[1:] in ([],["--native"]), "usage: [--native]")
    if "--native" not in argv:
        verify(replay_python)
        return
    from verify_multivariate_terminal_nine_shared_matched_upper import __file__ as existing
    with TemporaryDirectory(prefix="sfc-nine-marked-") as directory:
        library=Path(directory)/"bernstein.so"
        source=Path(existing).with_suffix(".cpp")
        run(("c++","-std=c++17","-O3","-fPIC","-shared",str(source),"-o",str(library)),check=True)
        verify(native_replay(CDLL(str(library))))


if __name__=="__main__":
    main()
