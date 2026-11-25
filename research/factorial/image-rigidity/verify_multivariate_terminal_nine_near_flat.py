#!/usr/bin/env python3
"""Exact Cayley ray separation at the fully shared terminal-nine flat endpoint."""

from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial
from sys import argv

from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID,
    I,
    ONE,
    ZERO,
    bisect,
    interval,
    projections,
    rational,
    require,
)
from verify_multivariate_terminal_six_marked_row_tail import projected_interval
from verify_multivariate_terminal_seven_near_flat import positive, uniform
from verify_multivariate_terminal_eight_near_flat import (
    coefficients as octic_coefficients,
    validate_constants as validate_octic,
)


NONIC = 9
MAX_LAM = Q(22361, 100000)
UPPER_LAM = Q(113, 500)
R0, R1 = Q(507, 500), Q(10177, 10000)
Q0 = Q(13, 250)
EXPONENTS = (Q(0), Q(1), Q(7, 4), Q(9, 4), Q(5, 2))
ALL_ONES = (1, 1, 1, 1, 1)
MEAN_DIRECTION = (0, 0, Q(1), Q(15, 7), Q(20, 7))

SUM_S_CAP = Q0 / (1 - Q0)
RECIPROCAL_CAP = (2 - Q0) / (8 * (1 - Q0)**2)
BASE_MEAN = Q(3, 4)
BASE_DERIVATIVE = BASE_MEAN * (1 - BASE_MEAN**2)
H_CAP = BASE_DERIVATIVE / (
    3 * (1 - BASE_MEAN**2 * SUM_S_CAP**2) * (1 - Q0)**3
)


def mean_difference_cap(a):
    return a * (1 - a*a) * (BASE_MEAN**2 - a*a) / (
        5 * (1 - BASE_MEAN**2 * SUM_S_CAP**2)
        * (1 - a*a * SUM_S_CAP**2) * (1 - Q0)**5
    )


J3_CAP = mean_difference_cap(Q(1, 2))
J4_CAP = J3_CAP + mean_difference_cap(Q(1, 4))


def skew_weight(power, i):
    moments = tuple((m - 4)**(2 * power) for m in range(NONIC))
    return (Q(NONIC - 2*i, NONIC) * sum(moments) - sum(moments[i:NONIC-i])) / (
        4 * power
    )


FIRST_SKEW = tuple(skew_weight(1, i) for i in range(5))
SECOND_SKEW = tuple(skew_weight(2, i) for i in range(5))
SKEW_TAIL = NONIC * (4 * RECIPROCAL_CAP)**6 / (
    12 * (1 - (4 * RECIPROCAL_CAP * Q0)**2)
)


def error_caps():
    mean_errors = [Q(0)]
    skew_errors = [Q(0)]
    differences = (0, 0, 0, J3_CAP, J4_CAP)
    for i in range(1, 5):
        c = MEAN_DIRECTION[i] * H_CAP
        a = FIRST_SKEW[i] * RECIPROCAL_CAP**2
        b = SECOND_SKEW[i] * RECIPROCAL_CAP**4 + SKEW_TAIL * Q0**2
        full_skew = a + b * Q0**2
        difference = differences[i]
        mean_errors.append(
            c*c + (full_skew + a) * b / 2 + full_skew**4 * Q0**2 / 20
            + Q(501, 500) * c * full_skew**2 * Q0
        )
        skew_errors.append(
            SKEW_TAIL + full_skew**3 / 5 + Q0 * (difference*a + c*b)
            + Q0**2 * c*c * (full_skew + full_skew**3 * Q0**4 / 5)
            + Q0**3 * c * full_skew**3 / 5
        )
    return tuple(mean_errors), tuple(skew_errors)


MEAN_ERRORS, SKEW_ERRORS = error_caps()


def validate_constants():
    validate_octic()
    u = 1 - MAX_LAM**2
    require(Q(1, 20) <= MAX_LAM**2 and 2 * UPPER_LAM < 1, "left half-plane")
    require(R1 * MAX_LAM**2 < UPPER_LAM**2 < Q0, "nonic scale upper")
    require(SUM_S_CAP < 1 and 4 * RECIPROCAL_CAP * Q0 < 1, "series domains")
    lower_reciprocal_cap = Q(2, 125)
    require((u**-2 - 1) / 7 < lower_reciprocal_cap
            and Q(64, 7) / (9 - lower_reciprocal_cap) < R1, "upper ratio")
    require((u * (1 + u) / 126)**2 > (R0 - 1)**2 * (1 + (1 - u*u) / 63)**3,
            "lower ratio")
    require(min(FIRST_SKEW + SECOND_SKEW) >= 0, "skew vector signs")
    require(max(MEAN_DIRECTION) * H_CAP * Q0**3 < Q(1, 100), "mean amplitude")
    skew_limit = (max(FIRST_SKEW) * RECIPROCAL_CAP**2 * Q0**2
                  + max(SECOND_SKEW) * RECIPROCAL_CAP**4 * Q0**4
                  + SKEW_TAIL * Q0**6)
    require(skew_limit < Q(1, 100), "skew amplitude")


@lru_cache(maxsize=None)
def cayley(k, i):
    return Q(comb(NONIC, i), 2**NONIC) * sum(
        (-1)**(i-l) * comb(i, l) * comb(NONIC-i, k-l)
        for l in range(max(0, k-NONIC+i), min(k, i)+1)
    )


@lru_cache(maxsize=None)
def weights(k):
    return tuple(2 * cayley(k, i) if k % 2 or i else Q(0) for i in range(5))


@lru_cache(maxsize=None)
def power_series(exponent):
    coefficients = [Q(1)]
    for m in range(1, 10):
        coefficients.append(coefficients[-1] * (m - 1 - exponent) / m)
    return tuple(coefficients)


@lru_cache(maxsize=None)
def scalar_series(k, vector):
    terms = tuple(a*b for a, b in zip(weights(k), vector, strict=True))
    series = tuple(power_series(exponent) for exponent in EXPONENTS)
    coefficients = [
        sum(a*b[m] for a, b in zip(terms, series, strict=True)) for m in range(9)
    ]
    first = next((m for m, coefficient in enumerate(coefficients) if coefficient), None)
    if first is None:
        require(not any(terms), "zero family")
        return None
    tail = sum(abs(a*b[9]) for a, b in zip(terms, series, strict=True))
    return first, coefficients[first], tuple(map(rational, coefficients[first:])), rational(tail)


def scaled_family(q, lam, power, vector, mean):
    terms = []
    min_degree = 0 if power == 0 and vector == ALL_ONES else 1
    for k in range(1 if mean else 0, NONIC+1, 2):
        row = scalar_series(k, vector)
        if row is None:
            continue
        offset, _, coefficients, tail = row
        value = coefficients[-1]
        for coefficient in reversed(coefficients[:-1]):
            value = coefficient + q*value
        degree = 2 * (power + offset) - (NONIC - k)
        require(degree >= min_degree, ("scaled Cayley divisibility", k, power, vector))
        remainder = (tail * q**(9-offset) / (ONE-q)).symmetric_error()
        terms.append((k, (value + remainder) * lam**degree))
    return tuple(terms)


def coefficients(order, lam):
    if order == 8:
        return octic_coefficients(8, lam, ONE)
    q = lam.square()
    central = dict(scaled_family(q, lam, 0, ALL_ONES, True))
    for k in range(NONIC + 1):
        caps = MEAN_ERRORS if k % 2 else SKEW_ERRORS
        error = sum(abs(a)*c for a, c in zip(weights(k), caps, strict=True))
        central[k] = central.get(k, ZERO) + uniform(error) * lam**(12-NONIC+k)
    return tuple(central.items()), (
        (positive(H_CAP), scaled_family(q, lam, 3, MEAN_DIRECTION, True)),
        (interval(-J3_CAP, Q(0)), scaled_family(q, lam, 5, (0, 0, 0, 1, 0), True)),
        (interval(-J4_CAP, Q(0)), scaled_family(q, lam, 5, (0, 0, 0, 0, 1), True)),
        (positive(RECIPROCAL_CAP**4),
         scaled_family(q, lam, 4, tuple(a*a/2 for a in FIRST_SKEW), True)),
        (uniform(RECIPROCAL_CAP**2), scaled_family(q, lam, 2, FIRST_SKEW, False)),
        (uniform(RECIPROCAL_CAP**4), scaled_family(q, lam, 4, SECOND_SKEW, False)),
        (uniform(H_CAP * RECIPROCAL_CAP**2),
         scaled_family(q, lam, 5, tuple(a*h for a, h in zip(FIRST_SKEW, MEAN_DIRECTION, strict=True)), False)),
    )


def chosen_projections(order,x,y,families):
    real=[ONE]
    imag=[ZERO]
    for k in range(order):
        real.append(real[-1]*x-imag[-1]*y)
        imag.append(real[-2]*y+imag[-1]*x)
    yield imag
    yield real
    # The largest signed term is the leading skew. Use a
    # projection orthogonal to its value at a provisional centre.
    for _, terms in (families[4:5] if order==9 else families[-3:-2]):
        ar=projected_interval(terms,real)
        ai=projected_interval(terms,imag)
        mr=(ar.lo+ar.hi)//2
        mi=(ai.lo+ai.hi)//2
        scale=max(abs(mr),abs(mi))
        if scale:
            a=rational(Q(mr,scale))
            b=rational(Q(mi,scale))
            yield [a*v-b*u for u,v in zip(real,imag,strict=True)]
    yield from projections(order,x,y)



def quadratic_range(alo, ahi, blo, bhi):
    # Lower and upper endpoint quadratics on each half of [-1, 1].
    low = min(0, alo+blo, -ahi+blo)
    high = max(0, ahi+bhi, -alo+bhi)
    if blo>0:
        for a in (alo,-ahi):
            if -2*blo<a<0:
                low=min(low,(-a*a)//(4*blo))
    if bhi<0:
        for a in (ahi,-alo):
            if 0<a<-2*bhi:
                high=max(high,-((a*a)//(4*bhi)))
    return I(low,high)


MEAN_LINK = BASE_DERIVATIVE * (1 - SUM_S_CAP**2/3)**2 * (4 * RECIPROCAL_CAP)**3 / 3
RANGE_SKEW_2 = rational(RECIPROCAL_CAP**2)
RANGE_SKEW_4 = rational(RECIPROCAL_CAP**4)
RANGE_SKEW_ERROR = uniform(Q(2, 5) * RECIPROCAL_CAP**4)
RANGE_MEAN_LINK = rational(MEAN_LINK)
RANGE_MEAN = rational(H_CAP)


def projected_nonic_value(nominal, h, j, k, b, a, c, mixed, fams):
    # Combine the two skew directions before bounding their shared signed moment.
    nominal=nominal+fams[1][0]*j+fams[2][0]*k
    nominal=nominal+RANGE_SKEW_ERROR*c+fams[6][0]*mixed
    a=RANGE_SKEW_2*a+RANGE_SKEW_4*c
    b=RANGE_SKEW_4*b
    low_h,high_h=I(h.lo,h.lo),I(h.hi,h.hi)
    low_quad=b.lo+(RANGE_MEAN_LINK*low_h).lo if h.lo>=0 else b.lo
    high_quad=b.hi+(RANGE_MEAN_LINK*high_h).hi if h.hi<=0 else b.hi
    quad=quadratic_range(a.lo,a.hi,low_quad,high_quad)
    remainder_lo=(RANGE_MEAN*low_h).lo if h.lo<0 else 0
    remainder_hi=(RANGE_MEAN*high_h).hi if h.hi>0 else 0
    return nominal+quad+I(remainder_lo,remainder_hi)


def nonic_projection(base,fams,direction):
    nominal=projected_interval(base,direction)
    values=(projected_interval(terms,direction) for _,terms in fams)
    return projected_nonic_value(nominal,*values,fams)


def enclosing_axis(rectangle, a, b):
    left,right=rectangle
    low=a*(left.lo if a>=0 else left.hi)+b*(right.lo if b>=0 else right.hi)
    high=a*(left.hi if a>=0 else left.lo)+b*(right.hi if b>=0 else right.lo)
    return I(low//GRID,-(-high//GRID))


def ray_height(p,x,b):
    num=p*b*x.square()
    y=I(0,b.hi)
    for _ in range(3):
        update=b-num/(ONE+p*b*y)
        y=I(max(y.lo,update.lo),min(y.hi,update.hi))
    return y


def powers(x,y):
    real=[ONE]
    imag=[ZERO]
    for k in range(9):
        real.append(real[-1]*x-imag[-1]*y)
        imag.append(real[-2]*y+imag[-1]*x)
    return real,imag


def nonic_rectangles_sign(x, y, base, fams, coarse, taylor):
    # Centering the complex Taylor enclosure reduces cancellation loss over a wide ray tube.
    real,imag=powers(x,y)
    if taylor:
        xm=(x.lo+x.hi)//2
        ym=(y.lo+y.hi)//2
        real0,imag0=powers(I(xm,xm),I(ym,ym))
        dx=x-I(xm,xm)
        dy=y-I(ym,ym)
        radius_square=dx.square()+dy.square()
        radius=(x.square()+y.square()).sqrt()
        modulus=[ONE]
        for k in range(1,8):
            modulus.append(modulus[-1]*radius)
        def box(terms):
            re=projected_interval(terms,real0)
            im=projected_interval(terms,imag0)
            derivative=[(k-1,I(k*a.lo,k*a.hi)) for k,a in terms if k>=1]
            pr=projected_interval(derivative,real0)
            pi=projected_interval(derivative,imag0)
            second=ZERO
            for k,a in terms:
                if k>=2:
                    amax=max(abs(a.lo),abs(a.hi))
                    second+=I(amax*(k*(k-1)//2),amax*(k*(k-1)//2))*modulus[k-2]
            rem=(second*radius_square).symmetric_error()
            return re+pr*dx-pi*dy+rem,im+pr*dy+pi*dx+rem
    else:
        def box(terms):
            return projected_interval(terms,real),projected_interval(terms,imag)
    if any(z.lo>0 or z.hi<0 for z in box(coarse)):
        return True
    center=box(base)
    grouped=tuple(box(terms) for _,terms in fams)
    candidates=[(GRID,0),(0,GRID)]
    def propose(rectangle, perpendicular=True):
        # The center proposes the axis. The interval test below certifies its sign.
        xr,yi=rectangle
        u,v=xr.lo+xr.hi,yi.lo+yi.hi
        scale=max(abs(u),abs(v))
        if scale:
            u,v=(u*GRID)//scale,(v*GRID)//scale
            candidates.append((v,-u) if perpendicular else (u,v))
    propose(center,False)
    for i in (4,5,0,3,6,1,2):
        propose(grouped[i])
    skx=RANGE_SKEW_2*grouped[4][0]+RANGE_SKEW_4*grouped[5][0]
    sky=RANGE_SKEW_2*grouped[4][1]+RANGE_SKEW_4*grouped[5][1]
    propose((skx,sky))
    for a,b in candidates:
        nominal=enclosing_axis(center,a,b)
        values=(enclosing_axis(pair,a,b) for pair in grouped)
        z=projected_nonic_value(nominal,*values,fams)
        if z.lo>0 or z.hi<0:
            return True
    return False


def enclose_families(base,fams):
    result=dict(base)
    for parameter,terms in fams:
        for k,a in terms:
            result[k]=result[k]+parameter*a
    return tuple(result.items())


RATIO=interval(R0,R1)
ROOT_RATIO=RATIO.sqrt()

def verify_ray(order, index):
    limit = MAX_LAM if order == 8 else UPPER_LAM
    initial = (interval(Q(0), limit), interval(Q(0), Q(2)))
    widths = tuple(max(1, value.hi - value.lo) for value in initial)

    @lru_cache(maxsize=16000)
    def shape(lo, hi):
        lam = I(lo, hi)
        base, families = coefficients(order, lam)
        return base, families, enclose_families(base, families)

    _, _, coarse = shape(initial[0].lo, initial[0].hi)
    circle = dict(coarse)
    require(circle[order].lo * 2**order > sum(
        max(abs(value.lo), abs(value.hi)) * 2**k
        for k, value in circle.items() if k != order
    ), ("scaled root disk", order))
    ray0, ray1 = map(rational, RAYS[index])
    total = 0
    stack = [(initial, 0)]
    while stack:
        (lam, x), depth = stack.pop()
        total += 1
        base, families, coarse = shape(lam.lo, lam.hi)
        deficit = lam.square()
        if order == 8:
            level = ray0 + ray1 * deficit
        else:
            level = (ray0 + ray1 * deficit / RATIO) / ROOT_RATIO
        y = ray_height(deficit, x, level)
        if order == 9 and nonic_rectangles_sign(x, y, base, families, coarse, False):
            continue
        if order == 9 and nonic_rectangles_sign(x, y, base, families, coarse, True):
            continue
        for direction in chosen_projections(order, x, y, families):
            if order == 9:
                value = nonic_projection(base, families, direction)
            else:
                value = projected_interval(base, direction)
                for parameter, terms in families:
                    value = value + parameter * projected_interval(terms, direction)
            if value.lo > 0 or value.hi < 0:
                break
        else:
            box = (lam, x)
            require(depth < 42, ("unresolved ray", order, index,
                                 [value.shown() for value in box]))
            weights = (100, 60) if order == 9 else (100, 50)
            k = max(range(2), key=lambda i: weights[i] * (box[i].hi - box[i].lo)
                    * GRID // widths[i])
            for half in bisect(box[k]):
                stack.append((box[:k] + (half,) + box[k+1:], depth + 1))
    return total


RAYS=(
    (Q(2245,2000),Q(27,40)),(Q(9577,10000),Q(241,500)),
    (Q(973,1250),Q(36,125)),(Q(3181,5000),Q(199,1000)),
    (Q(1201,2500),Q(29,250)),(Q(869,2500),Q(73,1000)),
    (Q(201,1000),Q(7,200)),
)

def anchor():
    for k in range(1, 10, 2):
        offset, leading, _, _ = scalar_series(k, ALL_ONES)
        require(2 * offset == 9 - k, ("nonic Hermite degree", k))
        expected = Q(factorial(9), 32**offset * factorial(offset) * factorial(k))
        require(leading == expected, ("nonic Hermite coefficient", k))
    def f(j,v):
        a=Q(1,2*(j-1))
        return sum(Q((-1)**k*factorial(j),2**k*factorial(k)*factorial(j-2*k))*a**k*v**(j//2-k)
                   for k in range(j//2+1))
    brackets=[]
    for order,values in ((8,(20,191,560,1226)),(9,(65,269,642,1272))):
        enclosures=[]
        for number in values:
            lo,hi=Q(number,1000),Q(number+1,1000)
            require(f(order,lo)*f(order,hi)<0,("Hermite bracket",order,lo,hi))
            enclosures.append((lo,hi) if order==8 else (R0*lo,R1*hi))
        brackets.append(enclosures)
    a,b=brackets
    sequence=tuple(pair for i in range(3,-1,-1) for pair in (b[i],a[i]))
    for (high,low),(ray,_) in zip(zip(sequence,sequence[1:]),RAYS,strict=True):
        require(high[0]>ray**2>low[1],"anchor positions")
    for p in (0,MAX_LAM**2):
        heights=[r+s*p for r,s in RAYS]
        require(all(a>b>0 for a,b in zip(heights,heights[1:])),"ray order")

def main():
    validate_constants()
    anchor()
    cases = [(order, ray) for order in (8, 9) for ray in range(7)]
    if len(argv) > 1:
        require(len(argv) == 3, "expected the order and ray index")
        case = tuple(map(int, argv[1:]))
        require(case in cases, ("order and ray", case))
        cases = [case]
    total = 0
    for order, ray in cases:
        cells = verify_ray(order, ray)
        total += cells
        print(f"Order {order}, ray {ray}: {cells} cells.", flush=True)
    print(f"Verified {total} Cayley cells.", flush=True)


if __name__ == "__main__":
    main()
