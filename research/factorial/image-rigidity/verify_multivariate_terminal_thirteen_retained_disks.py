#!/usr/bin/env python3
"""Outward root-disk tests for retained profiles with a small residual loss."""

from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations_with_replacement
from math import comb, factorial, prod

from flint import acb, acb_poly, arb

from verify_multivariate_terminal_seven_private_linear import skew_weights
from verify_multivariate_terminal_six_retained_rows import ball


ORDERS = (12, 13)
REFERENCE_RATIO = Q(144, 143)
SKEW_TERMS = 4
MEAN_TERMS = 4
EXP_TERMS = 4
DISK_RADII = tuple(Q(1, 2**p) for p in range(3, 22))


@lru_cache(maxsize=None)
def log_shapes(order):
    points = tuple(1-Q(2*h, order-1) for h in range(order))
    distances = tuple(min(k, order-k) for k in range(order+1))
    weights = tuple(skew_weights(order, r) for r in distances)
    gaussian = tuple(Q(k*(order-k), order-1) for k in range(order+1))
    skew = tuple(tuple((1 if 2*k>order else -1)
                       * sum(w*p**(2*l) for w,p in zip(row, points, strict=True))
                       for k,row in enumerate(weights)) for l in range(1,SKEW_TERMS+1))
    mean = tuple(tuple(sum(p-p**(2*l+1) for p in points[:r]) for r in distances)
                 for l in range(1,MEAN_TERMS+1))
    return (gaussian, *skew, *mean), tuple(sum(map(abs,row)) for row in weights), distances


def symmetric(value):
    return arb(0, value.upper())


def residual_model(order, fixed_powers, means, cap):
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
    params = (alpha,
              *(symmetric(T**(2*l-1)*W/(2*l)) for l in range(1,SKEW_TERMS+1)),
              *(arb(0).union(T**(2*l)*W/(2*l+1)) for l in range(1,MEAN_TERMS+1)))
    shape, sizes, distances = log_shapes(order)
    vectors = tuple(tuple(ball(x) for x in row) for row in shape)
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
    if len(real) != 1 or real[0][0] != ORDERS[1]:
        return None
    upper.sort(key=lambda entry: float((entry[1].real/abs(entry[1])).mid()), reverse=True)
    if tuple(order for order,_ in upper) != (ORDERS[1],ORDERS[0])*(ORDERS[0]//2):
        return None
    return upper, real[0]


def verify_disks(powers, means, cap):
    models = {order:residual_model(order, powers[order], means, cap) for order in ORDERS}
    if any(model is None for model in models.values()):
        return False
    centers = choose_centers({order:model[0] for order,model in models.items()})
    if centers is None:
        return False
    upper, (odd,negative) = centers
    for radius in map(ball,DISK_RADII):
        if not negative.real < -radius or any(not z.imag > radius for _,z in upper):
            continue
        if any(not (a.conjugate()*b).imag > radius*(abs(a)+abs(b)+radius)
               for (_,a),(_,b) in zip(upper,upper[1:])):
            continue
        break
    else:
        return False
    for order, center in (*upper, (odd,negative)):
        shift = acb_poly([center,1])
        coefficients = sum((param*polynomial(shift) for param,polynomial in models[order][1]),
                           acb_poly())
        if not radius*abs(coefficients[1]).lower() > (
                abs(coefficients[0]).upper()
                + sum(radius**k*abs(coefficients[k]).upper() for k in range(2,order+1))):
            return False
    return True
