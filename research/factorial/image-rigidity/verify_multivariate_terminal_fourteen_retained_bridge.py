#!/usr/bin/env python3
"""Exact terminal-fourteen retained complements of the high-u bridge."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import comb, factorial, prod
from tempfile import TemporaryDirectory

from flint import acb_poly, arb

from verify_multivariate_terminal_fourteen_retained_profiles import (
    CAP, ORDERS, REFERENCE_RATIO, build_projector, configure_replay, retained,
)
from verify_multivariate_terminal_fourteen_retained_upper import (
    DISK_RADII, choose_centers, verify_disks,
)
from verify_multivariate_terminal_thirteen_retained_disks import (
    EXP_TERMS, MEAN_TERMS, SKEW_TERMS, log_shapes, symmetric,
)
from verify_multivariate_terminal_six_retained_rows import ball


LOWER, UPPER = Q(19, 25), Q(171, 200)
RESIDUAL_CAP = Q(1, 60)
SINGLETON = ((4, 1, -1),)
SINGLETON_ENDPOINTS = (Q(526, 625), Q(87, 100), Q(91, 100), Q(9469, 10000))
SKEW_CENTERS = (Q(3, 4), Q(2, 3), Q(5, 8))
PRIVATE = ((), ())

# Outer cells, radial cells and disk cells, in exact retained-list order.
COUNTS = (
    (43, 11735, 0), (96, 21459, 40), (4, 851, 0), (4, 835, 0),
    (16, 3443, 0), (16, 3401, 0), (160, 0, 160),
    (16, 4172, 0), (16, 4555, 0), (4, 846, 0), (16, 3622, 0),
    (16, 3997, 0), (42, 9996, 11),
)


def centered_model(order, fixed_powers, means, cap, signed):
    midpoint = sum(means)/2
    ref_ratio = Q(1) if order == ORDERS[0] else REFERENCE_RATIO
    reference = tuple(comb(order,k) * ball(value)**ball(Q(1,order))
                      * ball(midpoint)**ball(ref_ratio*Q(k*(order-k),order-1))
                      for k,value in enumerate(fixed_powers))
    loss = -ball(means[0]).union(ball(means[1])).log()
    if order == ORDERS[0] or cap == 0:
        rho = ball(ref_ratio)
        reciprocal = cap
    else:
        upper = 1 + (-ball(1-cap/ORDERS[1]).log()/ball(1+(ORDERS[0]-1)*cap).log())
        rho = ball(ref_ratio).union(upper)
        reciprocal = ORDERS[0]*cap/(ORDERS[1]-cap)
    alpha = -rho*loss-ball(ref_ratio)*ball(midpoint).log()
    W = arb((rho*loss).upper())
    T = ball(Q((order-1)*reciprocal, 2+(order-1)*reciprocal))
    n,N = ORDERS
    lower_T = ball(Q((n-1)*cap, 2+(n-1)*cap))
    pmax = lower_T*loss/2
    first = ball(signed[0]).union(ball(signed[1]))*pmax
    nominal = ball(sum(signed)/2)*lower_T*(-ball(midpoint).log())/2
    if order == N:
        eta = cap/(2*N+(n*n-2)*cap)
        transport = ball(REFERENCE_RATIO**2*(1+eta**2))
        first = transport*first+symmetric(ball(2*REFERENCE_RATIO**2*eta)*pmax)
        nominal *= transport
    shape,sizes,distances = log_shapes(order)
    vectors = [tuple(ball(x) for x in row) for row in shape]
    vectors[1] = tuple(vectors[1][k] + sum(ball(SKEW_CENTERS[l-2]/l)*T**(2*l-2)*vectors[l][k]
                         for l in range(2,SKEW_TERMS+1)) for k in range(order+1))
    skew = tuple(symmetric(ball((1-SKEW_CENTERS[l-2])/(2*l))*T**(2*l-1)*W)
                 for l in range(2,SKEW_TERMS+1))
    lower_abs = max(arb(0),first.lower(),(-first).lower())
    means_moments = []
    for l in range(1,MEAN_TERMS+1):
        hi = T**(2*l)*W/(2*l+1)
        low = ((2*lower_abs)**(2*l)/W**(2*l-1)/(2*l+1) if W>0 else arb(0))
        means_moments.append(low.union(hi))
    params = (alpha,first-nominal,*skew,*means_moments)
    reference = tuple(c*(nominal*direction).exp()
                      for c,direction in zip(reference,vectors[1],strict=True))
    skew_tail = T**(2*SKEW_TERMS+1)*W/((2*SKEW_TERMS+2)*(1-T*T))
    mean_tail = T**(2*MEAN_TERMS+2)*W/((2*MEAN_TERMS+3)*(1-T*T))
    base_errors = []
    for k in range(order+1):
        tail = ball(sizes[k])*skew_tail + distances[k]*mean_tail
        length = tail + sum(abs(param*row[k]) for param,row in zip(params,vectors,strict=True))
        if not length < 1:
            return None
        error = (length**(EXP_TERMS+1)/factorial(EXP_TERMS+1)+tail)/(1-length)
        base_errors.append(error)

    groups = []
    negligible = ball(Q(1,10**23))
    magnitudes = tuple(max(abs(entry).upper() for entry in row) for row in vectors)
    for degree in range(1,EXP_TERMS+1):
        for indices in combinations_with_replacement(range(len(params)),degree):
            denominator = prod(factorial(indices.count(i)) for i in set(indices))
            param = prod(params[i] for i in indices)/denominator
            crude = abs(param)*prod(magnitudes[i] for i in indices)
            if crude < negligible:
                for k in range(order+1):
                    base_errors[k] += crude
                continue
            polynomial = acb_poly([reference[k]*prod(vectors[i][k] for i in indices)
                                   for k in range(order+1)])
            groups.append((param,polynomial))
    base = [value*(1+symmetric(error)) for value,error in zip(reference,base_errors,strict=True)]
    return reference, ((arb(1),acb_poly(base)), *groups)


def verify_centered_signed_disks(powers, means, cap, signed):
    models = {order:centered_model(order, powers[order], means, cap, signed) for order in ORDERS}
    if any(model is None for model in models.values()):
        return False
    centers = choose_centers({order:model[0] for order,model in models.items()})
    if centers is None:
        return False
    upper, (odd,negative) = centers
    disks = []
    for order, center in (*upper, (odd,negative)):
        shift = acb_poly([center,1])
        coefficients = sum((param*polynomial(shift) for param,polynomial in models[order][1]),
                           acb_poly())
        left = abs(coefficients[1]).lower()
        magnitudes = tuple(abs(coefficients[k]).upper() for k in range(order+1))
        for radius in map(ball,DISK_RADII):
            if radius*left > magnitudes[0] + sum(radius**k*magnitudes[k] for k in range(2,order+1)):
                break
        else:
            return False
        disks.append((center,radius))
    upper_disks, (negative, radius) = disks[:-1], disks[-1]
    return (negative.real < -radius
            and all(z.imag > r for z,r in upper_disks)
            and all((a.conjugate()*b).imag > r*abs(b)+s*abs(a)+r*s
                    for (a,r),(b,s) in zip(upper_disks,upper_disks[1:])))


def verify_centered_disks(powers, means, cap):
    pending = [((Q(-1),Q(1)),0)]
    while pending:
        signed,depth = pending.pop()
        if verify_centered_signed_disks(powers,means,cap,signed):
            continue
        if depth>=6:
            return False
        midpoint = sum(signed)/2
        pending.extend((((signed[0],midpoint),depth+1),((midpoint,signed[1]),depth+1)))
    return True


def verify_singleton(means):
    retained.require((SINGLETON_ENDPOINTS[0], SINGLETON_ENDPOINTS[-1]) == means,
                     ("singleton outer coverage", means))
    for l, center in enumerate(SKEW_CENTERS, start=2):
        q = 2*l - 2
        retained.require(q**q * center**(q+1) <= (q+1)**(q+1) * (1-center)**q,
                         ("centered skew inequality", l, center))
    powers = {order: retained.fixed_powers(order, SINGLETON, PRIVATE) for order in ORDERS}
    expected_counts = (32, 64, 64)
    shards = zip(SINGLETON_ENDPOINTS[:-1], SINGLETON_ENDPOINTS[1:], strict=True)
    for expected, shard in zip(expected_counts, shards, strict=True):
        pending = [(shard, 0)]
        accepted = 0
        while pending:
            cell, depth = pending.pop()
            cap = retained.cap_at(cell[0], RESIDUAL_CAP)
            if cell[1]-cell[0] <= Q(1, 1000) and verify_centered_disks(powers, cell, cap):
                accepted += 1
                continue
            retained.require(depth < 17, ("unresolved centered disk cell", cell, cap))
            midpoint = sum(cell)/2
            pending.extend((((cell[0], midpoint), depth+1), ((midpoint, cell[1]), depth+1)))
        retained.require(accepted == expected, ("centered disk count", shard, accepted, expected))
        print("Centered singleton shard", shard, accepted, flush=True)
    total = sum(expected_counts)
    return total, 0, total


def verify_interval(rows, means, replay):
    """Replay shared rays, with signed root disks above residual mean 17/20."""
    core = retained.core_coefficients(rows, PRIVATE)
    powers = {order: retained.fixed_powers(order, rows, PRIVATE) for order in ORDERS}
    coefficients = retained.bounds(core, means, retained.cap_at(means[0], RESIDUAL_CAP))
    global_annulus = retained.enclosing_annulus(coefficients)
    pending = [(means, 0)]
    outer = radial = disks = 0
    while pending:
        cell, depth = pending.pop()
        midpoint = sum(cell)/2
        cap = retained.cap_at(cell[0], RESIDUAL_CAP)
        targets = retained.anchor_angles(rows, PRIVATE, midpoint)
        if targets is not None:
            coefficients = retained.bounds(core, cell, cap)
            if global_annulus or retained.enclosing_annulus(coefficients):
                count = retained.verify_rays(
                    core, coefficients, cell, cap, targets, cell[1]-cell[0] <= Q(1, 200),
                    replay, model_builder=retained.refined_model)
                if count is not None:
                    outer += 1
                    radial += count
                    continue
        if cell[0] >= Q(17, 20) and cell[1]-cell[0] <= Q(1, 1000):
            if verify_disks(powers, cell, cap):
                outer += 1
                disks += 1
                continue
        retained.require(depth < 17, ("unresolved retained bridge cell", rows, cell, cap))
        pending.extend((((cell[0], midpoint), depth+1), ((midpoint, cell[1]), depth+1)))
    return outer, radial, disks


def verify(replay, selection):
    # Import after binding the common enumerator to orders thirteen and fourteen.
    from verify_multivariate_terminal_thirteen_retained_bands import bracket_tail, retained_cores
    import verify_multivariate_terminal_thirteen_shared_light as light

    cores = retained_cores(LOWER, RESIDUAL_CAP, CAP)
    retained.require(len(cores) == len(COUNTS), ("retained bridge list count", len(cores)))
    retained.require(cores[6][0] == SINGLETON, ("exact singleton orientation", cores[6]))
    totals = [0, 0, 0]
    for index in selection:
        rows, squared = cores[index]
        means = bracket_tail(LOWER, UPPER, squared)
        floor = Q(1, 40) if index < 6 else Q(1, 150)
        if light.FLOOR != floor:
            light.FLOOR = floor
            light.ordinary_ratio_upper.cache_clear()
        result = verify_singleton(means) if rows == SINGLETON else verify_interval(rows, means, replay)
        retained.require(result == COUNTS[index], ("retained bridge cell counts", index, result))
        totals = [total + count for total, count in zip(totals, result, strict=True)]
        print(index, rows, means, *result, flush=True)
    print(f"Verified {len(selection)} retained lists, {totals[0]} outer cells, "
          f"{totals[1]} radial cells, and {totals[2]} outer disk covers.", flush=True)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("indices", type=int, choices=range(len(COUNTS)), nargs="*")
    args = parser.parse_args()
    configure_replay()
    with TemporaryDirectory(prefix="sfc-fourteen-retained-bridge-") as directory:
        verify(build_projector(directory), args.indices or tuple(range(len(COUNTS))))


if __name__ == "__main__":
    main()
