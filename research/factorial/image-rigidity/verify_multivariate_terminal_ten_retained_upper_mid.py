#!/usr/bin/env python3
"""Exact retained-core separation between outer means .76 and .88."""

from ctypes import CDLL
from fractions import Fraction as Q
from math import prod
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

from flint import arb, ctx

from verify_multivariate_terminal_nine_retained_rows import power_bounds
from verify_multivariate_terminal_nine_shared_matched_upper import ANNULUS, replay_python
from verify_multivariate_terminal_six_intermediate_cayley import ONE, interval, require
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds
from verify_multivariate_terminal_six_retained_rows import ball
from verify_multivariate_terminal_ten_retained_profiles import (
    POSITIONS, RAY_SCALE, anchor_angles, bounds, cap_at, core_coefficients,
)
from verify_multivariate_terminal_ten_shared_light import (
    ORDERS, coefficient_model, normalized_native,
)


LOWER, UPPER = Q(19, 25), Q(22, 25)
RETAINED_FLOOR = Q(1, 45)
REFINED_ANCHORS = ((), ((3, 1, -1),), ((4, 1, -1),))
EXPECTED = (
    (50, 11449), (16, 1943), (8, 872), (8, 888), (32, 4201),
    (8, 900), (8, 891), (8, 930), (8, 910), (32, 4973), (16, 1871),
    (13, 1823), (8, 899), (8, 934), (17, 1989), (16, 1914), (8, 904),
)


def row_mean_square(lower, gap):
    return prod(Q(9*lower+h, 9*lower+8*gap+h) for h in range(1,gap+1))


def cores():
    strength_limit = Q(9, 4)*log_bounds(1/LOWER)[1]
    gap_limit = int(strength_limit*(1+1/(9*RETAINED_FLOOR)))
    eligible = []
    for gap in range(1, gap_limit+1):
        lower_limit = int((gap/RETAINED_FLOOR-1)/9)
        for lower in range(1, lower_limit+1):
            strength = Q(gap**2, lower+gap)
            ratio = Q(gap, 9*lower+1)
            if RETAINED_FLOOR < ratio <= Q(1,8) and strength <= strength_limit:
                square = row_mean_square(lower,gap)
                if square >= LOWER**2:
                    eligible.extend((lower,gap,sign,strength,square) for sign in (-1,1))

    def collect(index, rows, budget, square):
        if index == len(eligible):
            if rows <= tuple(sorted((a,d,-sign) for a,d,sign in rows)):
                yield rows,budget,square
            return
        lower,gap,sign,strength,row_square = eligible[index]
        for copies in range(1+int(budget/strength)):
            combined = square*row_square**copies
            if combined >= LOWER**2:
                yield from collect(index+1,rows+((lower,gap,sign),)*copies,
                                   budget-copies*strength,combined)

    return sorted(collect(0,(),strength_limit,Q(1)))


def grouped_projection(core, order, means, cap, cosine, replay):
    factors = tuple(interval(*factor) for factor in core[order])
    pending = [((Q(-1),Q(1)),0)]
    count = 0
    while pending:
        signed,depth = pending.pop()
        model = coefficient_model(order,cap,means,signed,extended=True)
        if model is None:
            continue
        cells = None
        if model:
            groups = tuple((parameter,tuple(a*b for a,b in zip(vector,factors,strict=True)))
                           for parameter,vector in model)
            cells = replay(order,(cosine,cosine),groups,Q(0))
        if cells is not None:
            count += cells
        elif depth < 5:
            middle = sum(signed)/2
            pending.extend((((signed[0],middle),depth+1),((middle,signed[1]),depth+1)))
        else:
            return None
    return count


def verify_rays(core, boxes, means, cap, anchor, grouped, replay):
    independent = {j: ((ONE,tuple(interval(*entry) for entry in boxes[j])),) for j in ORDERS}
    count = 0
    for left,right in zip(anchor,anchor[1:]):
        for position in POSITIONS:
            proposal = float(left.mid())*float(position)+float(right.mid())*float(1-position)
            cosine = Q(round(RAY_SCALE*proposal),RAY_SCALE)
            if not left > ball(cosine) > right:
                continue
            cells = []
            for j in ORDERS:
                result = replay(j,(cosine,cosine),independent[j],Q(0))
                if result is None and grouped:
                    result = grouped_projection(core,j,means,cap,cosine,replay)
                if result is None:
                    break
                cells.append(result)
            if len(cells) == 2:
                count += sum(cells)
                break
        else:
            return None
    return count


def verify_core(rows, budget, square, replay):
    core = core_coefficients(rows,False)
    core_lower,core_upper = power_bounds(square,Q(1,2))
    lower = max(power_bounds(Q(1,2),Q(17,18)*budget)[0],LOWER/core_upper)
    upper = min(Q(1),UPPER/core_lower)
    if upper < lower:
        return 0,0
    tolerance = arb(2)**-50 if rows in REFINED_ANCHORS else None
    pending = [((lower,upper),0)]
    outer_count = radial_count = 0
    while pending:
        means,depth = pending.pop()
        cap = min(RETAINED_FLOOR,cap_at(means[0]))
        midpoint = sum(means)/2
        anchor = anchor_angles(rows,False,midpoint,tolerance)
        count = None
        if anchor is not None:
            boxes = bounds(core,means,cap)
            # Strict Enestrom bounds for every admissible array on the cell.
            inner,outer = ANNULUS
            if all(a > inner*d and b < outer*c for j in ORDERS
                   for (a,b),(c,d) in zip(boxes[j],boxes[j][1:])):
                count = verify_rays(core,boxes,means,cap,anchor,depth>=3,replay)
        if count is not None:
            outer_count += 1
            radial_count += count
            continue
        require(depth < 11, ("unresolved residual interval", rows, means))
        pending.extend((((means[0],midpoint),depth+1),((midpoint,means[1]),depth+1)))
    return outer_count,radial_count


def verify(replay, indices, check_counts):
    ctx.prec = 160
    items = cores()
    require(len(items) == len(EXPECTED), "exact retained list count")
    require(all(0 <= index < len(items) for index in indices), "valid retained list index")
    for index in indices or range(len(items)):
        rows,budget,square = items[index]
        counts = verify_core(rows,budget,square,replay)
        if check_counts:
            require(counts == EXPECTED[index], ("replay counts",rows,counts))
        print(index,rows,*counts,flush=True)


def main():
    use_native = "--native" in argv
    selected = tuple(int(item) for item in argv[1:] if item != "--native")
    if not use_native:
        verify(replay_python,selected,False)
        return
    with TemporaryDirectory(prefix="sfc-ten-upper-mid-") as folder:
        library = Path(folder)/"bernstein.so"
        source = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
        run(("c++","-std=c++17","-O3","-fPIC","-shared",str(source),"-o",str(library)),check=True)
        verify(normalized_native(CDLL(str(library))),selected,True)


if __name__ == "__main__":
    main()
