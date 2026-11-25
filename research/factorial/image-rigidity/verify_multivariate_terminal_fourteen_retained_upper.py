#!/usr/bin/env python3
"""Signed root disks for the terminal-fourteen upper retained complements."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from math import comb, factorial, prod

from flint import acb, acb_poly, arb

from verify_multivariate_terminal_fourteen_retained_profiles import configure_replay, retained
from verify_multivariate_terminal_thirteen_retained_disks import (
    EXP_TERMS, MEAN_TERMS, SKEW_TERMS, log_shapes, symmetric,
)
from verify_multivariate_terminal_six_retained_rows import ball


ORDERS = (13, 14)
REFERENCE_RATIO = Q(169, 168)
DISK_RADII = tuple(Q(k,2**p) for p in range(24,3,-1) for k in range(8,16))


def residual_model(order, fixed_powers, means, cap, signed):
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
    params = (alpha,first-nominal,
              *(symmetric(T**(2*l-1)*W/(2*l)) for l in range(2,SKEW_TERMS+1)),
              *(arb(0).union(T**(2*l)*W/(2*l+1)) for l in range(1,MEAN_TERMS+1)))
    shape, sizes, distances = log_shapes(order)
    vectors = tuple(tuple(ball(x) for x in row) for row in shape)
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


def choose_centers(references):
    upper = []
    real = []
    for order, reference in references.items():
        try:
            roots = acb_poly(list(reference)).roots(tol=arb(2)**-91, maxprec=512)
        except (ArithmeticError, ValueError):
            return None
        left = [z for z in roots if z.imag < 0]
        right = [z for z in roots if z.imag > 0]
        if len(left) != order//2 or len(right) != order//2:
            return None
        upper.extend((order, acb(z.real.mid(),z.imag.mid())) for z in right)
        real.extend((order,acb(z.real.mid())) for z in roots if z.imag.contains(0))
    if len(real) != 1 or real[0][0] != ORDERS[0]:
        return None
    upper.sort(key=lambda entry: float((entry[1].real/abs(entry[1])).mid()), reverse=True)
    if tuple(order for order,_ in upper) != (ORDERS[1],ORDERS[0])*(ORDERS[0]//2)+(ORDERS[1],):
        return None
    return upper, real[0]


def verify_signed_disks(powers, means, cap, signed):
    models = {order:residual_model(order, powers[order], means, cap, signed) for order in ORDERS}
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


def verify_disks(powers, means, cap):
    pending = [((Q(-1),Q(1)),0)]
    while pending:
        signed,depth = pending.pop()
        if verify_signed_disks(powers,means,cap,signed):
            continue
        if depth>=3:
            return False
        midpoint = sum(signed)/2
        pending.extend((((signed[0],midpoint),depth+1),((midpoint,signed[1]),depth+1)))
    return True


# Outer endpoints, reciprocal cap, expected core lists, and disk cells.
BANDS = {
    "855-88": (Q(171, 200), Q(22, 25), Q(1, 70), 3, 145),
    "88-92": (Q(22, 25), Q(23, 25), Q(1, 80), 3, 178),
    "92-95": (Q(23, 25), Q(19, 20), Q(1, 100), 2, 48),
}


def residual_cells(bands, rows, means, cap):
    if rows != ((5, 1, -1),) or cap != Q(1, 70):
        yield rows, means, cap
        return
    tighter = Q(1, 85)
    yield rows, means, tighter
    for extra, squared in bands.retained_cores(means[0], tighter, cap, canonical=False):
        yield tuple(sorted(rows + extra)), bands.bracket_tail(*means, squared), tighter


def verify_interval(rows, means, declared):
    powers = {order: retained.fixed_powers(order, rows, ((), ())) for order in ORDERS}
    pending = [(means, 0)]
    count = 0
    while pending:
        cell, depth = pending.pop()
        cap = retained.cap_at(cell[0], declared)
        if cell[1] - cell[0] <= Q(1, 1000) and verify_disks(powers, cell, cap):
            count += 1
            continue
        retained.require(depth < 17, ("unresolved retained disk cell", rows, cell, cap))
        midpoint = sum(cell) / 2
        pending.extend((((cell[0], midpoint), depth + 1), ((midpoint, cell[1]), depth + 1)))
    return count


def verify_band(name):
    # Import after the order parameters have been configured.
    import verify_multivariate_terminal_thirteen_retained_bands as bands

    lower, upper, cap, expected_lists, expected_cells = BANDS[name]
    import verify_multivariate_terminal_thirteen_shared_light as light
    light.FLOOR = Q(1, 150) if name == "855-88" else Q(1, 40)
    lists = bands.retained_cores(lower, cap, Q(1, 8))
    retained.require(len(lists) == expected_lists, ("retained list count", name, len(lists)))
    total = 0
    for rows, squared in lists:
        means = bands.bracket_tail(lower, upper, squared)
        count = sum(verify_interval(core, cell, residual_cap)
                    for core, cell, residual_cap in residual_cells(bands, rows, means, cap))
        total += count
        print(name, rows, count, flush=True)
    retained.require(total == expected_cells, ("disk count", name, total, expected_cells))
    print(f"Verified {name}: {len(lists)} retained lists and {total} disk cells.", flush=True)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("band", choices=BANDS, nargs="?")
    args = parser.parse_args()
    configure_replay()
    for name in (args.band,) if args.band else BANDS:
        verify_band(name)


if __name__ == "__main__":
    main()
