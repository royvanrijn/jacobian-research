#!/usr/bin/env python3
"""Exact terminal-nine certificates for private factorial profiles."""

from fractions import Fraction as Q
from functools import lru_cache
from math import factorial, prod

from verify_multivariate_terminal_nine_high_curvature import anchor_value
from verify_multivariate_terminal_seven_private_rows import (
    DISK_RADIUS, coefficient_box, private_core, root_bounds, verify_annulus, verify_disk,
)
from verify_multivariate_terminal_seven_private_linear import (
    multiply_polynomials, skew_weights, weighted_log_upper,
)
from verify_multivariate_terminal_six_joint_ray_bands import (
    log_bounds, power_bounds, radial_replay, rounded_interval,
)
from verify_multivariate_terminal_six_private_angles import (
    coefficient_bounds, exponential_lower, multiply_interval, require,
)


ORDERS = (8, 9)
RAY_SCALE = 100000
RECIPROCAL_FLOOR = Q(1, 35)
REFERENCE_RATIO = Q(64, 63)
OUTER_DIVISIONS = 40
LOSS = (
    (1, Q(11, 10)), (2, Q(6, 5)), (3, Q(5, 4)), (7, Q(29, 20)),
    (10, Q(3, 2)), (17, Q(9, 5)), (25, Q(21, 10)),
)

HIGH = (
    ((1, 1, 1, 1), (), 0, (92832, 74840, 42638, 10439, -29476, -62457, -86603)),
    ((3,), (), 1, (92467, 76202, 41318, 11334, -30467, -63050, -87556)),
    ((1, 1, 1), (1,), 0, (92935, 74147, 43388, 10439, -28481, -62457, -85866)),
    ((1, 1), (1, 1), 0, (92986, 74147, 43575, 10439, -28149, -62457, -85617)),
    ((1, 1, 1), (), 8, (92441, 74840, 41507, 10439, -30137, -63640, -87321)),
    ((1, 2), (), 4, (92414, 75525, 41318, 10439, -30467, -63050, -87556)),
    ((1, 1), (), 25, (92441, 74840, 41696, 9543, -29806, -63050, -87083)),
    ((2,), (), 21, (92414, 75525, 41696, 10439, -30137, -63050, -87083)),
    ((1,), (1,), 24, (92441, 74147, 42073, 9543, -29476, -63640, -86603)),
    ((1, 1), (1,), 6, (92441, 73446, 42262, 9543, -29144, -63640, -86603)),
    ((2,), (1,), 2, (92414, 74147, 42262, 9543, -29144, -63050, -86603)),
)
LOW = (
    ((3,), (), (
        (91420, 69948, 37217, 6559, -34115, -63513, -86884),
        (91548, 70762, 37643, 6987, -33774, -63299, -86801),
    ), (((8371, 3993), (2993, 9372), (-4573, 9266), (-9894, 3819)), ((8566, 3545), (4203, 8791), (-2474, 9903), (-8324, 6432), (-10628, 0)))),
    ((1, 1, 1), (), (
        (90042, 66278, 33884, 2722, -36485, -64986, -87581),
        (90291, 66763, 34526, 3344, -36012, -64694, -87471),
        (90531, 67783, 35157, 4556, -35544, -64405, -87361),
        (90827, 68859, 35884, 5269, -35000, -64446, -87234),
        (91037, 69279, 36465, 5841, -34562, -64176, -87130),
        (91228, 70197, 37150, 6378, -34148, -63921, -87033),
        (91402, 71075, 37656, 6881, -33536, -63681, -86940),
        (91582, 71446, 38193, 7418, -33118, -63809, -87000),
        (91727, 71750, 38639, 7866, -32768, -63594, -86758),
    ), (((8072, 4224), (2635, 9476), (-4842, 9217), (-10033, 3773)), ((8305, 3722), (3896, 8906), (-2742, 9894), (-8491, 6382), (-10744, 0)))),
    ((1, 2), (), (
        (90812, 68233, 35649, 4420, -35249, -64220, -87301),
        (90996, 69156, 36188, 5555, -34836, -63963, -87204),
        (91174, 69527, 36717, 6079, -34428, -63709, -87107),
        (91398, 70479, 37327, 6687, -33954, -63799, -86993),
        (91554, 71333, 37813, 7173, -33572, -63562, -86902),
    ), (((8223, 4092), (2824, 9421), (-4704, 9249), (-9969, 3801)), ((8435, 3620), (4057, 8843), (-2603, 9903), (-8410, 6413), (-10692, 0)))),
    ((2,), (), (
        (84403, 60922, 23375, -8226, -43718, -67589, -89580),
        (85133, 61958, 24602, -7691, -42720, -67091, -89277),
        (85802, 62972, 25664, -6607, -41932, -66599, -89106),
        (86560, 64137, 27062, -4774, -41015, -65672, -89045),
        (87119, 64501, 28180, -3759, -40275, -65204, -88886),
        (87664, 64242, 29228, -1639, -39368, -65123, -88736),
        (88129, 65050, 30353, -150, -38709, -64709, -88593),
        (88646, 65903, 31540, 1408, -38000, -64263, -88584),
        (89034, 66603, 32407, 2812, -37406, -64259, -88456),
        (89414, 67250, 33217, 3577, -36845, -64282, -88335),
        (90051, 69041, 34816, 5697, -35503, -63961, -88388),
        (90364, 70123, 35543, 6396, -34983, -63639, -88277),
        (90619, 70625, 36216, 7048, -34495, -63720, -88172),
        (90896, 71601, 36839, 7654, -34039, -63823, -87923),
        (91121, 72015, 37415, 8219, -33611, -63561, -87829),
        (91322, 72392, 37949, 8744, -33210, -63316, -87741),
        (91503, 72736, 38443, 9233, -32835, -63475, -87658),
        (91784, 74047, 39370, 9826, -32126, -63297, -87607),
        (91917, 74309, 39774, 9623, -31806, -63101, -87536),
        (92035, 74055, 40009, 9997, -31509, -63315, -87310),
        (92140, 74271, 40349, 9730, -31233, -63147, -87247),
        (92254, 74466, 40660, 10048, -30977, -62991, -87189),
    ), (((6940, 4759), (1622, 9652), (-5629, 9179), (-10641, 3730)), ((7181, 4239), (2867, 9144), (-3597, 9911), (-9168, 6333), (-11346, 0)))),
    ((1, 1), (1,), (
        (90242, 69239, 36874, 6021, -33004, -64182, -86455),
        (90513, 69654, 37323, 6429, -32727, -64021, -86394),
        (90769, 70056, 37766, 6835, -32449, -63859, -86333),
        (91051, 71018, 38422, 6714, -32120, -63666, -86424),
        (91267, 71373, 38834, 7097, -31850, -63508, -86365),
        (91459, 71698, 39220, 7459, -31594, -63745, -86307),
        (91632, 71996, 39582, 7801, -31349, -63602, -86252),
    ), (((8612, 4469), (2977, 9524), (-4440, 9105), (-9502, 3701)), ((8872, 3929), (4277, 8978), (-2357, 9813), (-7999, 6280), (-10194, 0)))),
    ((2,), (1,), (
        (91090, 71175, 38534, 7605, -31873, -63503, -86195),
        (91280, 71484, 38892, 7335, -31644, -63369, -86144),
        (91451, 71769, 39229, 7648, -31424, -63240, -86095),
    ), (((8669, 4246), (3153, 9450), (-4341, 9178), (-9528, 3756)), ((8905, 3747), (4419, 8887), (-2242, 9852), (-7987, 6349), (-10235, 0)))),
)


OUTER_CASES = (
    ((1, 1), (), 25, (
        (93250, 74917, 44735, 11235, -27077, -60791, -85196),
        (93128, 74620, 44027, 11283, -27930, -61183, -85556),
        (93005, 74829, 43284, 10626, -28507, -61553, -85899),
        (92864, 75017, 42672, 10661, -29325, -61907, -86040),
        (92707, 74646, 42049, 10005, -29871, -62248, -86363),
        (92535, 74803, 41419, 10041, -30404, -62578, -86492),
        (92348, 74401, 40784, 9391, -30926, -62899, -86800),
        (92148, 73984, 40145, 9430, -31438, -62768, -86919),
        (91970, 73562, 39498, 8774, -31963, -63092, -86862),
        (91748, 73121, 38854, 8815, -32460, -63395, -87156),
        (91513, 72669, 38209, 8177, -32949, -63257, -87268),
        (91268, 72207, 37563, 8220, -33430, -63552, -87376),
        (91011, 71735, 36917, 7588, -34155, -63842, -87656),
        (90745, 70673, 36270, 6960, -34623, -64127, -87760),
        (90438, 70180, 35623, 6336, -35085, -63982, -87861),
        (90151, 69086, 34977, 5715, -35541, -64262, -87961),
        (89856, 68576, 34331, 5098, -35993, -64117, -88226),
        (89575, 67474, 33681, 4469, -36470, -64415, -88000),
        (89264, 66950, 33037, 3859, -36913, -64688, -88097),
        (88945, 66420, 32393, 3252, -37351, -64957, -88193),
        (88618, 65886, 31750, 2648, -37786, -64813, -88450),
        (88285, 65347, 31107, 1390, -38216, -65079, -88542),
        (87944, 64804, 30466, 795, -38643, -65342, -88634),
        (87598, 64259, 29825, 201, -39068, -65201, -88724),
        (87245, 63709, 29186, -1038, -39489, -65463, -88814),
        (86885, 63156, 28381, -2270, -40142, -65722, -88902),
        (86520, 63233, 27743, -3495, -40555, -65979, -88990),
        (86150, 62677, 27106, -4071, -40966, -66234, -89077),
        (85774, 62119, 26470, -5283, -41374, -66488, -89163),
        (85394, 62195, 25667, -5851, -41780, -66351, -89249),
        (85009, 61635, 25032, -6418, -42183, -66603, -89334),
        (84630, 60438, 24566, -6987, -42594, -67245, -89271),
        (84236, 59871, 23933, -7548, -42993, -67491, -89356),
        (83839, 59950, 23301, -8108, -43390, -67736, -89440),
        (83439, 59384, 22670, -8667, -43788, -67982, -89525),
        (83035, 58815, 22040, -9224, -44406, -68224, -89608),
        (82577, 58896, 21411, -9779, -44797, -68466, -89691),
        (82165, 58327, 20782, -10332, -45187, -68706, -89774),
        (81751, 58408, 19985, -10883, -45575, -68946, -89856),
        (81334, 57840, 19357, -10821, -45962, -69184, -89938),
        (80915, 57272, 18730, -11373, -46348, -69422, -90020),
    )),
    ((1,), (1,), 24, (
        (93250, 74917, 44735, 11235, -27077, -60791, -85196),
        (93143, 74668, 44117, 10691, -27833, -61585, -85529),
        (93037, 74385, 43616, 10834, -28312, -61890, -85653),
        (92912, 74622, 42951, 10283, -28775, -62180, -85962),
        (92770, 74296, 42426, 9736, -29223, -62456, -86071),
        (92612, 73951, 41895, 9193, -29658, -62721, -86175),
        (92439, 73588, 41361, 9343, -30080, -62976, -86461),
        (92249, 73764, 40823, 8810, -30491, -63222, -86556),
        (92080, 73382, 40278, 8271, -30910, -63472, -86468),
        (91863, 72980, 39738, 8429, -31302, -63263, -86740),
        (91632, 72565, 39198, 7910, -31684, -63491, -86825),
        (91387, 72700, 38659, 7398, -32057, -63711, -86907),
        (91128, 72265, 38121, 7566, -32422, -63926, -86987),
        (90855, 71819, 37585, 7064, -32778, -64135, -87064),
        (90569, 71928, 37051, 7237, -33127, -63909, -87139),
        (90270, 71466, 36519, 6745, -33468, -64109, -87388),
        (89958, 71558, 35990, 6259, -34053, -64304, -87457),
        (89658, 69962, 35450, 5761, -34152, -64510, -87357),
        (89321, 70046, 34925, 5945, -34722, -64695, -87424),
        (88937, 69551, 34404, 5475, -35034, -64458, -87489),
        (88575, 69049, 33887, 5664, -35341, -64637, -87552),
        (88201, 67963, 33373, 5202, -35640, -64812, -87613),
        (87776, 67446, 32863, 4747, -35934, -64983, -87673),
        (87379, 66925, 32355, 4943, -36223, -65151, -87732),
        (86967, 66396, 31853, 4497, -36505, -64903, -87789),
        (86503, 65862, 31354, 4056, -37026, -65065, -87845),
        (86110, 65323, 30860, 4261, -37296, -65223, -87899),
        (85663, 64779, 30370, 3829, -37561, -65377, -87952),
        (85205, 64230, 29718, 3404, -37820, -65124, -88003),
        (84736, 63676, 29237, 3615, -38074, -65273, -88054),
        (84256, 62536, 28761, 3199, -38323, -65420, -88103),
        (83775, 61981, 28448, 2781, -38574, -65567, -87987),
        (83271, 61415, 27981, 3000, -38813, -65707, -88034),
        (82755, 60846, 27353, 2599, -39286, -65843, -88245),
        (82280, 59693, 26895, -5831, -39516, -65978, -88126),
        (81742, 59119, 26445, -6187, -39740, -66109, -88170),
        (81193, 58542, 25999, -6537, -39960, -66237, -88213),
        (80632, 57961, 25559, -6272, -40174, -65971, -88255),
        (80060, 57952, 25125, -6612, -40384, -66095, -88296),
        (79424, 57364, 24698, -6947, -40590, -66216, -88336),
        (78828, 56774, 24276, -7277, -40792, -66335, -88376),
    )),
)


def verify_row_loss():
    count = 0
    for last, factor in LOSS:
        limit = Q(last, 8)
        base_cap = limit * factor / (factor - 1)**2
        base_cap = base_cap.numerator // base_cap.denominator
        gap_stop = 1
        while Q(gap_stop**2, base_cap + gap_stop) <= limit:
            gap_stop += 1
        for base in range(1, base_cap + 1):
            for gap in range(1, gap_stop):
                strength = Q(gap**2, base + gap)
                if strength > limit or Q(gap, base) <= factor - 1:
                    continue
                for order in ORDERS:
                    first = factorial(order * base)
                    endpoint = factorial(order * (base + gap))
                    for k in range(1, order):
                        power = Q(first**(order-k) * endpoint**k, factorial(order*base+k*gap)**order)
                        exponent = factor * strength * k * (order-k) / 2
                        require(power <= exponential_lower(exponent),
                                ("row loss", limit, base, gap, order, k))
                        count += 1
    return count


def verify_high():
    count = 0
    for increasing, decreasing, last, ray_values in HIGH:
        cosines = [Q(x, RAY_SCALE) for x in ray_values]
        knots = [Q(1), *cosines, Q(-1)]
        require(len(cosines) == 7, "high ray count")
        require(all(a > b for a, b in zip(knots, knots[1:])), "high ray order")
        for order in ORDERS:
            for k, x in enumerate(knots):
                previous_roots = (k + (order == 9)) // 2
                require((-1)**previous_roots * anchor_value(order, x) > 0,
                        ("binomial anchor signs", order, k))
            caps = coefficient_bounds(order, increasing, Q(last, 8), decreasing)
            box = [(Q(1), Q(1))] + [(Q(0), cap) for cap in caps[1:-1]] + [(Q(1), Q(1))]
            verify_annulus(box, order)
            for cosine in cosines:
                cells = radial_replay(order, box, (cosine, cosine))
                require(cells is not None, ("high ray", increasing, decreasing, last, order, cosine))
                count += cells
    return count


def verify_private_anchor(order, core, centers, cosines):
    upper = tuple((Q(x, 10000), Q(y, 10000)) for x, y in centers if y > 0)
    require(len(upper) == 4 and len(centers) == 4 + (order == 9), "private anchor counts")
    disks = upper + tuple((x, -y) for x, y in upper)
    require(all(y > DISK_RADIUS for _, y in upper), ("upper disks", order))
    if order == 9:
        x, y = centers[-1]
        require(y == 0 and Q(x, 10000) < -DISK_RADIUS, "negative anchor disk")
        disks += ((Q(x, 10000), Q(0)),)
    for i, (x, y) in enumerate(disks):
        verify_disk(core, (x, y))
        for a, b in disks[:i]:
            require((x-a)**2 + (y-b)**2 > (2*DISK_RADIUS)**2, "disjoint anchor disks")

    sectors = range(1 - (order == 9), 8, 2)
    for (x, y), sector in zip(upper, sectors, strict=True):
        for j, c in enumerate(cosines):
            sine = root_bounds(1-c*c, 2)
            lo, hi = multiply_interval(sine, (-x, -x))
            if j < sector:
                require(c*y + lo > DISK_RADIUS, ("anchor ray side", sector, j, order))
            else:
                require(c*y + hi < -DISK_RADIUS, ("anchor ray side", sector, j, order))


def cosine_at(knots, index, ray, strength):
    a, b = knots[index][ray], knots[index+1][ray]
    return Q(a, RAY_SCALE) + Q(b-a, RAY_SCALE)*(8*strength-index)


def verify_low_band(order, core, knots, index, ray):
    factor = next(f for last, f in LOSS if index + 1 <= last)
    pending = [((Q(index, 8), Q(index+1, 8)), 0)]
    count = 0
    while pending:
        strength, depth = pending.pop()
        cosines = tuple(sorted(cosine_at(knots, index, ray, s) for s in strength))
        cells = radial_replay(order, coefficient_box(order, core, factor, strength), cosines)
        if cells is not None:
            count += cells
            continue
        require(depth < 12, ("low ray", order, index, ray, strength))
        middle = sum(strength) / 2
        pending.extend((((strength[0], middle), depth+1), ((middle, strength[1]), depth+1)))
    return count


def verify_low_case(increasing, decreasing, knots, centers):
    require(any((increasing, decreasing, len(knots)-1) == row[:3] for row in HIGH),
            ("upper shared curvature cutoff", increasing, decreasing))
    for values in knots:
        require(len(values) == 7, "low ray count")
        sequence = (RAY_SCALE, *values, -RAY_SCALE)
        require(all(a > b for a, b in zip(sequence, sequence[1:])), "low ray order")
    total = 0
    for order, disks in zip(ORDERS, centers, strict=True):
        core = private_core(order, increasing, decreasing)
        verify_annulus(core, order)
        verify_private_anchor(order, core, disks, tuple(Q(x, RAY_SCALE) for x in knots[0]))
        total += sum(verify_low_band(order, core, knots, index, ray)
                     for index in range(len(knots)-1) for ray in range(7))
    return total


def possible_rows(limit):
    largest_gap = int(limit * (1 + 1/(8*RECIPROCAL_FLOOR)))
    for gap in range(1, largest_gap+1):
        largest_base = int((gap/RECIPROCAL_FLOOR-1)/8)
        for base in range(1, largest_base+1):
            x8 = Q(gap, 8*base+1)
            if x8 <= RECIPROCAL_FLOOR or Q(gap**2, base+gap) > limit:
                continue
            mean_square = prod(Q(8*base+h, 8*base+7*gap+h) for h in range(1, gap+1))
            yield x8, Q(gap, 9*base+1), mean_square


def reciprocal_caps(rows, lower_mean):
    cap8 = cap9 = RECIPROCAL_FLOOR
    for x8, x9, mean_square in rows:
        if lower_mean**2 <= mean_square:
            cap8, cap9 = max(cap8, x8), max(cap9, x9)
    return cap8, cap9


def verify_skew_convexity():
    count = 0
    for order in ORDERS:
        m = order - 1
        for k in range(1, (order+1)//2):
            weights = skew_weights(order, k)
            require(sum(h*w for h, w in enumerate(weights)) == 0, "skew derivative")
            numerator = [Q(0)] * (2*m-3)
            for h in range(1, m):
                term = [weights[h]*h*(m-h)]
                for a in range(1, m):
                    if a != h:
                        term = multiply_polynomials(term, (1, 2*a, a*a))
                numerator = [a+b for a, b in zip(numerator, term, strict=True)]
            for i in range(len(numerator)-1, 1, -1):
                if numerator[i] < 0:
                    numerator[i-2] += numerator[i]
                    numerator[i] = 0
            require(min(numerator) >= 0, ("skew convexity", order, k))
            count += 1
    return count


@lru_cache(maxsize=None)
def outer_exponents(order, cap8, cap9):
    if order == 8:
        cap, ratio, reference = cap8, (Q(1), Q(1)), Q(1)
    else:
        cap, reference = cap9, REFERENCE_RATIO
        low, high = log_bounds(1+7*cap8)
        ratio = (log_bounds(1+Q(64, 9)*cap8)[0]/high,
                 1+log_bounds(1/(1-cap8/9))[1]/low)
        require(ratio[0] <= reference <= ratio[1], "reference ratio enclosure")
    m = order-1
    low, high = log_bounds(1+m*cap)
    bounds = []
    for k in range(1, order):
        r = min(k, order-k)
        mean = sum(log_bounds((1+(m-h)*cap)/(1+h*cap))[0] for h in range(r))/high
        gaussian = Q(k*(order-k), m)
        skew = 2*weighted_log_upper(skew_weights(order, r), cap)/low
        a, b = rounded_interval(ratio[0]*(mean-skew), ratio[1]*(gaussian+skew))
        require(0 < a <= reference*gaussian <= b, ("coefficient exponent", order, k, cap))
        bounds.append((a, b))
    return bounds


def outer_coefficient_box(core, means, exponents):
    result = [(Q(1), Q(1))]
    for c, (lo, hi) in zip(core[1:-1], exponents, strict=True):
        bounds = (power_bounds(means[0], hi)[0], power_bounds(means[1], lo)[1])
        result.append(multiply_interval(c, bounds))
    return result + [(Q(1), Q(1))]


def verify_outer_case(increasing, decreasing, cutoff, knots):
    require(any((increasing, decreasing, cutoff) == row[:3] for row in HIGH),
            "paired high shared-curvature cutoff")
    require(len(knots) == OUTER_DIVISIONS+1, "outer ray knot count")
    for values in knots:
        require(len(values) == 7, "outer ray count")
        sequence = (RAY_SCALE, *values, -RAY_SCALE)
        require(all(a > b for a, b in zip(sequence, sequence[1:])), "outer ray order")
    anchor = (Q(1), *(Q(x, RAY_SCALE) for x in knots[0]), Q(-1))
    rows = tuple(possible_rows(Q(cutoff, 8)))
    count = 0
    for order in ORDERS:
        for k, x in enumerate(anchor):
            previous_roots = (k + (order == 9)) // 2
            require((-1)**previous_roots * anchor_value(order, x) > 0, ("outer binomial signs", order, k))
        core = private_core(order, increasing, decreasing)
        verify_annulus(core, order)
        for index in range(OUTER_DIVISIONS):
            first = Q(index, OUTER_DIVISIONS)
            last = Q(index+1, OUTER_DIVISIONS)
            pending = [((first, last), tuple(range(7)), 0)]
            while pending:
                means, rays, depth = pending.pop()
                exponents = outer_exponents(order, *reciprocal_caps(rows, means[0]))
                box = outer_coefficient_box(core, means, exponents)
                failed = []
                for ray in rays:
                    a, b = knots[index][ray], knots[index+1][ray]
                    cosines = tuple(sorted(
                        (a + (b-a)*(OUTER_DIVISIONS*mu-index))/RAY_SCALE for mu in means))
                    cells = radial_replay(order, box, cosines)
                    if cells is None:
                        failed.append(ray)
                    else:
                        count += cells
                if failed:
                    require(depth < 12, ("outer ray", increasing, decreasing, order, means, failed))
                    middle = sum(means)/2
                    pending.extend((((means[0], middle), failed, depth+1),
                                    ((middle, means[1]), failed, depth+1)))
    return count, len(rows)


def main():
    print(f"Verified {verify_row_loss()} shared row-index inequalities.", flush=True)
    print(f"Verified {verify_high()} high-band radial cells.", flush=True)
    total = 0
    for increasing, decreasing, knots, centers in LOW:
        cells = verify_low_case(increasing, decreasing, knots, centers)
        total += cells
        print(f"{increasing!r}, {decreasing!r}, through S0={Q(len(knots)-1, 8)}: "
              f"{cells} radial cells.", flush=True)
    print(f"Verified {total} low-band radial cells.", flush=True)
    print(f"Verified {verify_skew_convexity()} skew inequalities.", flush=True)
    for increasing, decreasing, cap, knots in OUTER_CASES:
        cells, rows = verify_outer_case(increasing, decreasing, cap, knots)
        print(f"{increasing!r}, {decreasing!r}: {cells} outer-band radial cells,"
              f" {rows} retained row types.", flush=True)


if __name__ == "__main__":
    main()
