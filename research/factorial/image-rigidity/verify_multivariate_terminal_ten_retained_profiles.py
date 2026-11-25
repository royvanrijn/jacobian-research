#!/usr/bin/env python3
"""Exact terminal-ten paired rays for retained ordinary exponent rows."""

from ctypes import CDLL
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from sys import argv

from flint import acb_poly, arb, ctx
from verify_multivariate_terminal_nine_shared_matched_upper import native_replay, replay_python
from verify_multivariate_terminal_six_retained_rows import ball, core_powers
from verify_multivariate_terminal_six_joint_ray_bands import log_bounds, radial_replay
from verify_multivariate_terminal_six_private_angles import multiply_interval, require
from verify_multivariate_terminal_seven_private_rows import verify_annulus
from verify_multivariate_terminal_six_intermediate_cayley import interval
from verify_multivariate_terminal_nine_retained_rows import power_bounds as powq
from verify_multivariate_terminal_ten_shared_light import (
    ORDERS, coefficient_model, independent_exponents, nominal_exponent,
)

CAP=Q(1,8)
SHARED_CUTOFF = Q(77, 10)
PRIVATE_CUTOFF = Q(57, 10)
RAY_SCALE=1000000
POSITIONS=(Q(1,2),Q(2,5),Q(3,5),Q(3,10),Q(7,10),Q(1,5),Q(4,5),Q(1,10),Q(9,10),Q(1,16),Q(15,16))
def cores(private):
    cutoff = PRIVATE_CUTOFF if private else SHARED_CUTOFF
    eligible = []
    gap_bound = int(cutoff*(1+1/(9*CAP)))
    for gap in range(1, gap_bound+1):
        lower_bound = int((gap/CAP-1)/9)
        for lower in range(1, lower_bound+1):
            strength = Q(gap**2, lower+gap)
            if Q(gap, 9*lower+1) > CAP and strength <= cutoff:
                eligible.extend((lower, gap, sign, strength) for sign in (-1, 1))

    # Traverse in the same lexicographic order used after column reversal.
    if not private:
        eligible.sort()

    def collect(index, rows, budget):
        if index == len(eligible):
            if private or rows and rows <= tuple(sorted((a, d, -sign) for a, d, sign in rows)):
                yield rows, budget
            return
        lower, gap, sign, strength = eligible[index]
        for copies in range(1+int(budget/strength)):
            yield from collect(index+1, rows+((lower, gap, sign),)*copies, budget-copies*strength)

    return sorted(collect(0, (), cutoff))


def fixed_powers(order, rows, private):
    retained = core_powers(order, rows) if rows else (Q(1),)*(order+1)
    if private:
        return tuple(p*Q(factorial(k)**order, factorial(order)**k) for k, p in enumerate(retained))
    return retained


@lru_cache(maxsize=32)
def core_coefficients(rows, private):
    return {order: tuple(powq(value, Q(1, order)) for value in fixed_powers(order, rows, private))
            for order in ORDERS}


def anchor_angles(rows, private, mean, tolerance=None):
    upper_roots = []
    for order in ORDERS:
        factors = fixed_powers(order, rows, private)
        coefficients = [arb(1)]+[
            comb(order,k)*ball(factors[k])**ball(Q(1,order))*ball(mean)**ball(nominal_exponent(order,k))
            for k in range(1,order)
        ]+[arb(1)]
        try:
            roots = acb_poly(coefficients).roots(tol=tolerance, maxprec=512)
        except (ArithmeticError, ValueError):
            return None
        upper = [z.real/abs(z) for z in roots if z.imag > 0]
        if len(upper) != order//2 or sum(z.imag < 0 for z in roots) != order//2:
            return None
        upper_roots.extend((cosine, order) for cosine in upper)
    upper_roots.sort(key=lambda item: float(item[0].mid()), reverse=True)
    if (tuple(order for _, order in upper_roots) != (10,9)*4+(10,)
        or any(not a>b for (a,_),(b,_) in zip(upper_roots, upper_roots[1:]))):
        return None
    return tuple(cosine for cosine, order in upper_roots)


def bounds(core,v,cap):
    tail={}
    for j in ORDERS:
        tail[j]=[(Q(1),Q(1))]+[(powq(v[0],b)[0],powq(v[1],a)[1])
                  for a,b in independent_exponents(j,cap)]+[(Q(1),Q(1))]
    out={}
    for j in ORDERS:
        out[j]=[(a*comb(j,k),b*comb(j,k)) for k,(a,b) in enumerate(
            multiply_interval(c,t) for c,t in zip(core[j],tail[j],strict=True))]
    for k in range(1,10):
        p=tail[9][k][1]*powq(v[1],Q(k*k,126))[1]
        q=tail[9][k-1][1]*powq(v[1],Q((10-k)**2,126))[1]
        a,b=out[10][k]
        out[10][k]=a,min(b,comb(10,k)*core[10][k][1]*min(p,q))
    return out

def cap_at(v):
    cap=min(CAP,(v**-2-1)/8)
    # Enlarge to a small rational grid for cheap power and log enclosures.
    u=cap*10000
    return min(CAP,Q(max(1,-(-u.numerator//u.denominator)),10000))

def verify_grouped(core, j, v, cap, c, replay):
    leaves=0
    factors=tuple(interval(*y) for y in core[j])
    todo=[((Q(-1),Q(1)),0)]
    while todo:
        p,depth=todo.pop()
        model=coefficient_model(j,cap,v,p)
        if model is None:
            continue
        if model:
            model=tuple((q, tuple(a*b for a,b in zip(r,factors,strict=True))) for q,r in model)
            cells=replay(j,(c,c),model,Q(0))
        else:
            cells=None
        if cells is not None:
            leaves+=cells
        elif depth<5:
            middle=sum(p)/2
            todo.extend((((p[0],middle),depth+1),((middle,p[1]),depth+1)))
        else:
            return None
    return leaves

def verify_rays(core, coefficients, v, cap, targets, grouped, replay):
    count=0
    for left,right in zip(targets,targets[1:]):
        for position in POSITIONS:
            a=float(left.mid())
            b=float(right.mid())
            c=Q(round(RAY_SCALE*(a*float(position)+b*float(1-position))),RAY_SCALE)
            if not left>ball(c)>right:
                continue
            cells=[]
            for j in ORDERS:
                answer=radial_replay(j,coefficients[j],(c,c))
                if answer is None and grouped:
                    answer=verify_grouped(core,j,v,cap,c,replay)
                if answer is None:
                    break
                cells.append(answer)
            if len(cells)==2:
                count+=sum(cells)
                break
        else:
            return None
    return count

def verify_core(rows, budget, private, replay):
    core=core_coefficients(rows, private)
    for j in ORDERS:
        verify_annulus([(a*comb(j,k),b*comb(j,k)) for k,(a,b) in enumerate(core[j])], j)
    lower=powq(Q(1,2),Q(17,18)*budget)[0]
    todo=[((lower,Q(1)),0)]
    leaves=cells=0
    while todo:
        v,depth=todo.pop()
        cap=cap_at(v[0])
        midpoint=sum(v)/2
        targets=anchor_angles(rows,private,midpoint)
        count=None
        if targets is not None:
            coeffs=bounds(core,v,cap)
            count=verify_rays(core,coeffs,v,cap,targets,depth>=5,replay)
        if count is not None:
            leaves+=1
            cells+=count
            continue
        if depth>=11:
            raise ArithmeticError((rows,tuple(map(float,v)),cap,targets))
        todo.extend((((v[0],midpoint),depth+1),((midpoint,v[1]),depth+1)))
    return leaves,cells

def verify(private, shard, shards, replay):
    ctx.prec = 160
    _, log_two_upper = log_bounds(2)
    require(log_bounds(1+Q(81,80))[0]/log_two_upper > 1+Q(1,126), "Gaussian attenuation")
    require(Q(190*17,18)*log_two_upper < 126, "tail attenuation conversion")
    items = cores(private)
    require(len(items) == (280 if private else 1025), "retained-list count")
    require(0 <= shard < shards, "valid shard index")
    total_leaves = total_cells = 0
    selected = items[shard::shards]
    for rows, budget in selected:
        leaves, cells = verify_core(rows, budget, private, replay)
        total_leaves += leaves
        total_cells += cells
    family = "one private row" if private else "fully shared"
    print(f"Verified {family}: {len(selected)} lists, {total_leaves} outer cells and {total_cells} radial cells.")


def main():
    args = [item for item in argv[1:] if item != "--native"]
    require(len(args) in (1, 3) and args[0] in ("private", "shared"),
            "expected private or shared, and optional shard index and total shard count")
    private = args[0] == "private"
    shard, shards = map(int, args[1:]) if len(args) == 3 else (0, 1)
    if "--native" not in argv:
        verify(private, shard, shards, replay_python)
        return
    with TemporaryDirectory(prefix="sfc-ten-retained-") as folder:
        library = Path(folder)/"bernstein.so"
        source = Path(__file__).with_name("verify_multivariate_terminal_ten_shared_light.cpp")
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source), "-o", str(library)),
            check=True)
        verify(private, shard, shards, native_replay(CDLL(str(library))))


if __name__ == "__main__":
    main()
