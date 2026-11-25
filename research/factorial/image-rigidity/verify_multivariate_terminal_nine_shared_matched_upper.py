#!/usr/bin/env python3
"""Exact terminal-nine matched rays between outer means .48 and .76."""

from ctypes import CDLL, POINTER, Structure, c_int, c_int64, c_uint64
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, prod
from pathlib import Path
from subprocess import run
from sys import argv
from tempfile import TemporaryDirectory

from verify_multivariate_terminal_eight_paired_bridge import (
    midpoint_signed, sine_intervals, sine_polynomials,
)
from verify_multivariate_terminal_eight_paired_light import bernstein_transform, halves, signed
from verify_multivariate_terminal_six_diffuse_bands import bernstein_coefficients
from verify_multivariate_terminal_nine_private_linear import ratio_bounds, shape
from verify_multivariate_terminal_nine_private_rows import (
    ORDERS, REFERENCE_RATIO, outer_coefficient_box, outer_exponents,
    verify_private_anchor, verify_skew_convexity,
)
from verify_multivariate_terminal_nine_shared_matched_light import ANNULUS, exp_small
from verify_multivariate_terminal_six_intermediate_cayley import (
    GRID, I, ONE, ZERO, interval, log, rational, scalar,
)
from verify_multivariate_terminal_six_joint_ray_bands import power_bounds, log_bounds
from verify_multivariate_terminal_six_private_angles import require


RAY_SCALE = 10**6
FLOOR = Q(1, 40)
BANDS = (
    (Q(1, 8), {
        48: (196334, 81204, -329959, -433921, -732385, -812464, -955756),
        50: (153983, -54668, -360416, -457613, -744124, -815016, -957755),
        52: (109612, -76651, -387531, -481049, -755658, -823444, -957932),
        54: (66325, -99673, -414433, -524412, -766996, -831716, -959931),
        55: (44512, -118994, -427807, -535487, -772594, -835797, -959226),
    },
     {8: ((1896, 9819), (-3691, 9294), (-7667, 6419), (-9737, 2277)), 9: ((2349, 9720), (-2811, 9597), (-6721, 7405), (-9167, 3995), (-10000, 0))}),
    (Q(1, 10), {
        55: (47927, -114315, -429985, -541292, -774108, -837197, -963291),
        56: (26088, -123360, -443334, -545543, -779701, -841267, -964223),
        57: (4141, -142943, -456634, -556610, -787147, -845303, -965146),
        58: (-17911, -152574, -469885, -567620, -792605, -849304, -966059),
        59: (-40066, -172408, -483088, -572275, -798021, -850001, -966964),
        60: (-62324, -182649, -496245, -583330, -803396, -854022, -966865),
        61: (-83735, -202744, -509355, -594332, -808730, -861109, -967780),
        62: (-106206, -232115, -522420, -599406, -814024, -864980, -968687),
    },
     {8: ((423, 9991), (-4609, 8874), (-8034, 5955), (-9780, 2086)), 9: ((848, 9964), (-3836, 9235), (-7227, 6911), (-9301, 3672), (-10000, 0))}),
    (Q(3, 40), {
        62: (-102271, -211324, -524192, -601371, -815852, -863417, -970942),
        63: (-124764, -240863, -537247, -606692, -821138, -867366, -971793),
        64: (-147351, -261036, -550259, -623504, -826386, -871284, -972637),
        65: (-170030, -281244, -563228, -634497, -831597, -872407, -973473),
        66: (-192800, -301486, -576154, -640135, -838234, -876350, -974302),
        67: (-215660, -321761, -590469, -651188, -843327, -880263, -975123),
        68: (-238607, -333906, -603271, -662196, -848385, -886662, -975938),
        69: (-261641, -354465, -616033, -668285, -853409, -890436, -976013),
        70: (-284759, -375059, -628755, -679359, -857116, -894183, -976839),
    },
     {8: ((-1099, 9939), (-5502, 8351), (-8380, 5457), (-9820, 1890)), 9: ((-714, 9974), (-4840, 8751), (-7709, 6370), (-9426, 3338), (-10000, 0))}),
    (Q(11, 200), {
        70: (-282391, -373279, -631528, -680935, -861078, -892976, -978529),
        71: (-305518, -393883, -644194, -691992, -866012, -894523, -979303),
        72: (-328726, -407110, -656820, -720505, -870915, -900587, -980071),
        73: (-352015, -428012, -669408, -726678, -875785, -904353, -980832),
        74: (-375382, -448949, -681959, -737143, -880625, -908094, -981588),
        75: (-398827, -469921, -694473, -743671, -885434, -911809, -982337),
        76: (-422348, -510232, -707999, -754209, -890213, -915500, -983081),
    },
     {8: ((-2891, 9573), (-6493, 7605), (-8753, 4835), (-9862, 1655)), 9: ((-2565, 9665), (-5964, 8027), (-8232, 5678), (-9561, 2931), (-10000, 0))}),
)


def exact_cap(declared, outer_lower):
    """Enumerate rows above FLOOR that can occur at this outer mean."""
    threshold = 1+7*FLOOR/(1+FLOOR)
    cap = min(FLOOR, declared)
    gap = 1
    while threshold**(-gap) >= outer_lower**2:
        base = 1
        while 8*base+1 < gap/FLOOR:
            x = Q(gap, 8*base+1)
            if cap < x <= declared:
                row_square = prod(Q(8*base+h, 8*base+7*gap+h) for h in range(1, gap+1))
                if row_square >= outer_lower**2:
                    cap = x
            base += 1
        gap += 1
    return cap


@lru_cache(maxsize=None)
def constants(order, cap8):
    cap = cap8 if order == 8 else 8*cap8/(9-cap8)
    m = order-1
    T = rational(m*cap/(2+m*cap))
    lo, hi = log_bounds(1+m*cap)
    s = interval(lo/2, hi/2)
    require(T.square().hi < GRID//2, "moment quotient monotonicity range")
    return T, s, T.square()/(scalar(2)*s)


@lru_cache(maxsize=None)
def leading_skew_direction(order, cap8):
    _, skew, _, _, _ = shape(order)
    T, _, _ = constants(order, cap8)
    weights = tuple(T**(2*l-2)*scalar(1, 2*l) for l in (2, 3))
    second, third = weights
    direction = tuple(a+second*b+third*c for a, b, c in zip(*skew, strict=True))
    return direction, weights


@lru_cache(maxsize=8192)
def coefficient_model(order, cap8, means, signed):
    midpoint = sum(means)/2
    skew_midpoint = sum(signed)/2
    ref_ratio = Q(1) if order == 8 else REFERENCE_RATIO
    T, s, skew_bound = constants(order, cap8)
    Pmax = constants(8, cap8)[2]
    U8 = -log(interval(*means))
    leading = interval(*signed)*Pmax*U8
    nominal = rational(skew_midpoint)*Pmax*(-log(rational(midpoint)))
    h = interval(means[0]/midpoint-1, means[1]/midpoint-1)
    magnitude = I(0, max(-h.lo,h.hi))
    log_error = rational(ref_ratio)*magnitude.square()/(scalar(2)*(ONE-magnitude))
    outer_mean = I(-log_error.hi,0)
    strength = U8
    if order == 9:
        rho = ratio_bounds(cap8)
        outer_mean += (rational(ref_ratio)-rho)*U8
        strength = rho*U8
        eta = cap8/(18+62*cap8)
        transfer = rational(REFERENCE_RATIO**2*(1+eta**2))
        transfer_width = rational(2*REFERENCE_RATIO**2*eta)
        leading = transfer*leading+(transfer_width*Pmax*U8).symmetric_error()
        nominal = transfer*nominal
    gaussian, skew, mean, skew_sizes, distances = shape(order)
    direction, weights = leading_skew_direction(order, cap8)
    reference = []
    for k, v in enumerate(direction):
        if k == 0 or k == order:
            reference.append(ONE)
            continue
        exponent = ref_ratio*Q(k*(order-k), order-1)
        value = scalar(comb(order, k))*interval(*power_bounds(midpoint, exponent))
        reference.append(value*exp_small(nominal*v))
    pieces = [(h, tuple(rational(ref_ratio)*v for v in gaussian)),
              (outer_mean, gaussian), (leading-nominal, direction)]
    for vector, weight in zip(skew[1:], weights, strict=True):
        bound = weight*skew_bound*strength
        pieces.append((bound.symmetric_error(), vector))
    lower = I(leading.lo, leading.lo) if leading.lo > 0 else ZERO
    upper_strength = I(strength.hi, strength.hi)
    for j, vector in enumerate(mean, 1):
        lo = (scalar(2)*lower/upper_strength)**(2*j)*upper_strength*scalar(1, 2*j+1)
        hi = T**(2*j+1)*strength/(scalar(2*j+1)*s)
        if lo.lo > hi.hi:
            # This rectangle is empty by the lower and upper moment bounds.
            return None
        pieces.append((I(lo.lo, hi.hi), vector))
    skew_tail = T**7*strength/(scalar(8)*(ONE-T.square()))
    defect_tail = T**6*strength/(scalar(7)*(ONE-T.square()))
    errors = []
    for k, c in enumerate(reference):
        tail = skew_sizes[k]*skew_tail+scalar(distances[k])*defect_tail
        size = tail
        for parameter, vector in pieces:
            term = parameter*vector[k]
            size += I(0, max(-term.lo, term.hi))
        if size.hi >= GRID:
            # Subdivide before applying the exponential bound.
            return ()
        error = scalar(1, 6)*size**3/(ONE-size)+tail*(ONE+size)
        errors.append(c*(ONE+error.symmetric_error()))
    products = []
    for i, (first, a) in enumerate(pieces):
        for j, (second, b) in enumerate(pieces[:i+1]):
            parameter = scalar(1, 2)*first.square() if i == j else first*second
            products.append((parameter, tuple(v*w for v, w in zip(a, b, strict=True))))
    return ((ONE, tuple(errors)),) + tuple(
        (parameter, tuple(c*v for c, v in zip(reference, vector, strict=True)))
        for parameter, vector in (*pieces, *products)
    )


def verify_band(order, cap8, rays, left, right, replay):
    means = Q(left, 100), Q(right, 100)
    pending = [(means, (Q(0), Q(1)), tuple(range(7)), 0, 0)]
    total = 0
    while pending:
        u, p, active, mean_depth, skew_depth = pending.pop()
        model = coefficient_model(order, cap8, u, p)
        if model is None:
            continue
        bad = []
        for ray in active:
            a, b = rays[left][ray], rays[right][ray]
            cosines = tuple(sorted(Q(a, RAY_SCALE)+Q(b-a, RAY_SCALE)*(100*v-left)/(right-left)
                                   for v in u))
            slope = Q(b-a,RAY_SCALE)*100*sum(u)/(2*(right-left))
            count = replay(order, cosines, model, slope) if model else None
            if count is None:
                bad.append(ray)
            else:
                total += count
        if not bad:
            continue
        require(max(mean_depth, skew_depth) < 16,
                ("unresolved ray", order, cap8, left, right, u, p, bad))
        if mean_depth < 2*skew_depth+2:
            center = sum(u)/2
            pending.extend((((u[0], center), p, bad, mean_depth+1, skew_depth),
                            ((center, u[1]), p, bad, mean_depth+1, skew_depth)))
        else:
            center = sum(p)/2
            pending.extend(((u, (p[0], center), bad, mean_depth, skew_depth+1),
                            (u, (center, p[1]), bad, mean_depth, skew_depth+1)))
    return total


def verify_setup(cap8, rays, anchors):
    left, right = Q(min(rays), 100), Q(max(rays), 100)
    for values in rays.values():
        bounds = (RAY_SCALE, *values, -RAY_SCALE)
        require(len(values) == 7 and all(a > b for a,b in zip(bounds, bounds[1:])), "ordered rays")
    small, large = ANNULUS
    for order in ORDERS:
        raw = tuple((Q(comb(order,k)),)*2 for k in range(order+1))
        exponents = outer_exponents(order, cap8, 8*cap8/(9-cap8))
        box = outer_coefficient_box(raw, (Q(0), right), exponents)
        require(sum(cap*small**k for k,(_,cap) in enumerate(box) if k) < 1, ("inner radius",order))
        require(sum(cap*large**k for k,(_,cap) in enumerate(box[:-1])) < large**order,
                ("outer radius",order))
        ratio = Q(1) if order == 8 else REFERENCE_RATIO
        anchor = [(Q(1),Q(1))]
        for k in range(1,order):
            lo,hi = power_bounds(left, ratio*Q(k*(order-k),order-1))
            anchor.append((comb(order,k)*lo,comb(order,k)*hi))
        anchor.append((Q(1),Q(1)))
        verify_private_anchor(order, anchor, anchors[order],
                              tuple(Q(x,RAY_SCALE) for x in rays[min(rays)]))


def polynomial_at(coefficients, value):
    result = Q(0)
    for coefficient in reversed(coefficients):
        result = result*value+coefficient
    return result


def derivative(coefficients):
    return tuple(k*coefficient for k, coefficient in enumerate(coefficients) if k) or (0,)


@lru_cache(maxsize=4096)
def directional_parts(order, left, right, slope):
    center = (left+right)/2
    square_radius_half = (right-left)**2/8
    base, linear = [], []
    for polynomial in sine_polynomials(order):
        first = derivative(polynomial)
        second = derivative(first)
        curvature = bernstein_coefficients(second, left, right)
        value = polynomial_at(polynomial, center)
        base.append(interval(value+min(0,min(curvature))*square_radius_half,
                             value+max(0,max(curvature))*square_radius_half))
        linear.append(rational(polynomial_at(first, center)*slope))

    bounds = interval(left, right)
    full = sine_intervals(order, bounds.lo, bounds.hi)
    return full, tuple(base), tuple(linear)


def signed_direction(directions, index):
    return directions[index] if index >= 0 else -directions[-index]


def projected(order, directions, model):
    full, base, linear = directions
    transform = bernstein_transform(order)
    result = []
    for pivot in range(order+1):
        groups = []
        for group, (_, coefficients) in enumerate(model):
            direction = base if group == 0 else full
            products = []
            for k, coefficient in enumerate(coefficients):
                product = coefficient*signed_direction(direction, k-pivot)
                if group == 1:
                    product += model[0][1][k]*signed_direction(linear, k-pivot)
                products.append(product)
            bounds = []
            for row in transform:
                lo = hi = 0
                for a, b in zip(row, products, strict=True):
                    values = (a.lo*b.lo, a.lo*b.hi, a.hi*b.lo, a.hi*b.hi)
                    lo += min(values)
                    hi += max(values)
                bounds.append(I(lo//GRID, -(-hi//GRID)))
            groups.append(tuple(bounds))
        result.append(tuple(groups))
    return tuple(result)


def replay_python(order, cosines, model, slope):
    directions = directional_parts(order, *cosines, slope)
    pending = [(projected(order, directions, model), 0)]
    count = 0
    while pending:
        projections, depth = pending.pop()
        if any(signed(p, model) for p in projections):
            count += 1
            continue
        if depth >= 14 or not any(midpoint_signed(p, model, order) for p in projections):
            return None
        left, right = [], []
        for projection in projections:
            children = tuple(halves(group) for group in projection)
            a, b = zip(*children, strict=True)
            left.append(a)
            right.append(b)
        pending.extend(((tuple(left), depth+1), (tuple(right), depth+1)))
    return count


class Words(Structure):
    _fields_ = [("low", c_uint64), ("high", c_int64)]


@lru_cache(maxsize=None)
def native_transform(order):
    values = tuple(Words(x & (2**64-1), x>>64)
                   for row in bernstein_transform(order) for v in row for x in (v.lo, v.hi))
    return (Words*len(values))(*values)


@lru_cache(maxsize=4096)
def native_directions(order, left, right, slope):
    result = []
    for directions in directional_parts(order, left, right, slope):
        values = tuple(x for v in directions for x in (v.lo, v.hi))
        if any(x < -2**63 or x >= 2**63 for x in values):
            return None
        result.append((c_int64*len(values))(*values))
    return tuple(result)


def native_replay(library):
    function = library.sfc_fused_radial
    function.argtypes = (c_int, c_int, POINTER(c_int64), POINTER(Words),
                         POINTER(c_int64), POINTER(c_int64), POINTER(c_int64))
    function.restype = c_int

    def replay(order, cosines, model, slope):
        values = []
        for parameter, coefficients in model:
            require(len(coefficients) == order+1, ("coefficient count", order))
            values.extend((parameter.lo, parameter.hi))
            for coefficient in coefficients:
                values.extend((coefficient.lo, coefficient.hi))
        directions = native_directions(order, *cosines, slope)
        if directions is None or any(x < -2**63 or x >= 2**63 for x in values):
            return replay_python(order, cosines, model, slope)
        inputs = (c_int64*len(values))(*values)
        full, base, linear = directions
        result = function(order+1, len(model), full, native_transform(order),
                          inputs, base, linear)
        if result == -2:
            return replay_python(order, cosines, model, slope)
        return None if result < 0 else result

    return replay


def verify(replay, selection):
    verify_skew_convexity()
    for declared, rays, anchors in BANDS:
        knots = tuple(sorted(rays))
        if selection and knots[0] not in selection:
            continue
        cap8 = exact_cap(declared, Q(knots[0],100))
        verify_setup(cap8, rays, anchors)
        total = 0
        for left, right in zip(knots,knots[1:]):
            counts = tuple(verify_band(order, cap8, rays, left, right, replay) for order in ORDERS)
            total += sum(counts)
            print(left, right, counts, flush=True)
        print(f"Verified {total} radial cells for {knots[0]}–{knots[-1]} at cap {declared}.", flush=True)


def main():
    arguments = tuple(x for x in argv[1:] if x != "--native")
    selection = tuple(map(int, arguments))
    require(set(selection) <= {min(rays) for _, rays, _ in BANDS}, "invalid band selection")
    if "--native" not in argv:
        verify(replay_python, selection)
        return
    with TemporaryDirectory(prefix="sfc-nine-") as directory:
        library = Path(directory)/"bernstein.so"
        source = Path(__file__).with_suffix(".cpp")
        run(("c++", "-std=c++17", "-O3", "-fPIC", "-shared", str(source), "-o", str(library)),
            check=True)
        verify(native_replay(CDLL(str(library))), selection)


if __name__ == "__main__":
    main()
