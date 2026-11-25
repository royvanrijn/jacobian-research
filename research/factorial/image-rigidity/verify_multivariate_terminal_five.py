"""Outward Arb replay for arbitrary exponent profiles at terminal five.

Run with Python and python-flint. The middle-range rays are fixed exact
rational inputs. All sign decisions use outward enclosures, and a box
that cannot be excluded before the search limit causes failure.
"""
from fractions import Fraction

from flint import acb, arb, ctx


# Upper endpoints, followed by the three ray cosines scaled by 10**12.
# The initial lower endpoint is 13/20.
HIGH_RAYS = """53/80 -342784287741 -581233353482 -880863663225
27/40 -365641971116 -597234868492 -885584617232
11/16 -388583496615 -613201312765 -890274035140
7/10 -411608264135 -629129172499 -894931982859
57/80 -434715662702 -645015059011 -899558565258
29/40 -457905072483 -660855711117 -904153924562
3/4 -492841585594 -684525058787 -910988818777
61/80 -527959076867 -708075589681 -917754611625
31/40 -551470219914 -723705282639 -922227189713
63/80 -575060204511 -739275378942 -926669757909
4/5 -598728393604 -754783392558 -931082647681
13/16 -622474152295 -770226963203 -935466216375
33/40 -646296849027 -785603854272 -939820845714
67/80 -670195856691 -800911950471 -944146940354
17/20 -694170553649 -816149255214 -948444926511
69/80 -718220324687 -831313887809 -952715250663
7/8 -742344561895 -846404080483 -956958378317
71/80 -766542665480 -861418175268 -961174792872
9/10 -790814044519 -876354620788 -965364994547
91/100 -812720455467 -889729832711 -969114190256
23/25 -832241915553 -901564323910 -972429646981
93/100 -851809282314 -913346761120 -975729240636
187/200 -866514762936 -922149059080 -978193695954
47/50 -876332615623 -928000678098 -979831857358
189/200 -886161787120 -933838951996 -981466214148
19/20 -896002243196 -939663813528 -983096803052
191/200 -905853949936 -945475197665 -984723661115
24/25 -915716873739 -951273041561 -986346825690
193/200 -925590981327 -957057284527 -987966334435
97/100 -935476239747 -962827867984 -989582225309"""


# The initial lower endpoint for this table is 1/500.
LOW_RAYS = """9/4000 748050451449 221878721019 -544732322914
1/400 747334991222 221170065593 -545786275723
3/1000 746314019074 220164028542 -547254393876
1/250 744415025018 218310306525 -549879678089
3/500 740994755114 215031482456 -554297888587
1/100 734981187048 209447570670 -561264048642
1/20 699569621102 179713580243 -590634288211
3/40 652414907696 143942940587 -616567333101
1/10 614313517700 117804459531 -632224651725
1/8 573650446875 91823179147 -646164115584
3/20 532209068787 65638663602 -659097103118
27/160 495999069122 41301994915 -669885410447
3/16 466116361608 19056363131 -678861734006
33/160 436672137018 -3238176502 -687654200530
9/40 407081356271 -25607111268 -696305955277
39/160 377331780566 -48071425026 -704847458102
21/80 347413914240 -70647841823 -713300182036
87/320 324858529284 -87661514725 -719590663203
9/32 309757624209 -99045018193 -723763578099
93/320 294612001990 -110462442440 -727921132351
3/10 279420987469 -121914469143 -732064149351
5/16 261639946934 -135319525313 -736880152607
13/40 241240944730 -150698567174 -742361983306
27/80 220757813784 -166140574934 -747820903612
7/20 200189511231 -181645074425 -753257371441
29/80 179535121545 -197211078348 -758671579201
3/8 158793846924 -212837124558 -764063505520
31/80 137964999212 -228521314458 -769432958868
2/5 117047992953 -244261351275 -774779614426
33/80 96042339276 -260054577988 -780103045313
17/40 74947640369 -275898014709 -785402749033
7/16 53763584393 -291788395268 -790678169877
9/20 32489940705 -307722202813 -795928717867
37/80 11126555322 -323695704219 -801153784727
19/40 -10326653456 -339704983149 -806352757318
39/80 -31869699266 -355745971602 -811525028854
1/2 -53502531807 -371814479837 -816670008228
41/80 -75225040848 -387906224568 -821787127676
21/40 -97037060148 -404016855363 -826875849020
43/80 -118938371276 -420141979199 -831935668663
11/20 -140928707335 -436277183153 -836966121518
9/16 -163007756577 -452418055224 -841966783996
91/160 -179625050590 -464524785614 -845697500701
183/320 -187952285821 -470577844129 -847556449211
23/40 -193510646089 -474612947424 -848793362728
19/32 -213008084937 -488733052990 -853107458993
49/80 -246587833910 -512922330068 -860447996614
101/160 -280362706752 -537078452460 -867718374808
13/20 -314330979501 -561187691434 -874918088347"""


def real(value):
    if isinstance(value, Fraction):
        return arb(value.numerator) / value.denominator
    return arb(value)


def span(lower, upper):
    return real(lower).union(real(upper))


def symmetric(radius):
    return span(-radius, radius)


def nonnegative_square(value):
    lower = value.lower().max(0)
    upper = value.upper()
    return span(lower * lower, upper * upper)


def cosh_minus_one(value):
    sine = value.sinh()
    return sine * sine / (value.cosh() + 1)


def remainder(t, parameter):
    scale = 1 - t
    original = t * ((1 + parameter * t * t) / (parameter + t * t)).sqrt()
    first_derivative = Fraction(4, 3) if parameter == 2 else Fraction(3, 2)
    return (original - 1 + real(first_derivative) * scale) / (scale * scale)


def quartic_ray_excluded(u, radius, ray):
    mu = 1 - u
    nu = 1 + ray
    lower = u.lower()
    g = (1 / (lower * lower) - 1) / 3
    skew_cap = ((1 + g) * (1 + 2 * g) / (1 + 3 * g)).root(4).log()
    sine = symmetric(skew_cap.sinh())
    cosine_correction = span(0, cosh_minus_one(skew_cap))

    radius_squared = nonnegative_square(radius)
    cosine = (1 + radius_squared).sqrt()
    radial_correction = radius_squared / (cosine + 1)
    imaginary_coefficient = nu - mu + ray * radial_correction + u * cosine_correction
    if not (radius * imaginary_coefficient).overlaps(u * cosine * sine):
        return True

    shape = span(Fraction(2, 9), remainder(lower, 2))
    symmetric_value = 4 * nu**2 - 8 * mu * nu + 6 * shape * mu**2
    real_value = symmetric_value + 8 * u * ray * cosine_correction
    real_value -= radius_squared * (4 + 8 * u * ray * (1 + cosine_correction) / (cosine * (cosine + 1)))
    return not real_value.contains(0)


def quintic_ray_excluded(w, radius, ray):
    mu = 1 - w
    nu = 1 + ray
    lower = w.lower()
    g = (1 / (lower * lower) - 1) / 4

    b_base = ((1 + g) * (1 + 2 * g) * (1 + 3 * g))**2 / (1 + 4 * g)**3
    c_base = (1 + 2 * g)**4 / ((1 + g) * (1 + 3 * g) * (1 + 4 * g))
    link_base = ((1 + 3 * g) * (1 + g))**4 / ((1 + 2 * g)**6 * (1 + 4 * g))
    b_cap = b_base.log() / 10
    link_cap = link_base.log() / 10
    c_cap = (c_base.log() / 10).min((b_cap + link_cap) / 2)

    c_sine = c_cap.sinh()
    c_cosine_correction = cosh_minus_one(c_cap)
    d = 5 * w * symmetric(b_cap.sinh())
    sinh_difference_cap = link_cap * (2 * c_cap).max(b_cap).cosh() + 2 * c_sine * c_cosine_correction
    delta = symmetric((5 * mu * c_sine + 5 * w * sinh_difference_cap).upper())

    p_correction = 5 * w * span(0, cosh_minus_one(b_cap))
    q_correction = 10 * w * span(0, c_cosine_correction)
    p = 5 * w - 1 + p_correction
    radius_squared = nonnegative_square(radius)
    cosine = (1 + radius_squared).sqrt()
    radial_correction = radius_squared / (cosine + 1)
    imaginary_coefficient = 4 * nu - 5 * mu + 4 * ray * radial_correction + p_correction
    imaginary_value = (radius * imaginary_coefficient - d * cosine) * (2 * (nu + radial_correction))
    if not imaginary_value.overlaps(delta):
        return True

    shape = span(Fraction(3, 8), remainder(lower, 3))
    symmetric_value = 4 * nu**2 - 10 * mu * nu + 10 * shape * mu**2
    real_value = symmetric_value + p_correction * (2 * ray - 1) + q_correction
    k_zero = 8 * nu**2 - 20 * nu + 15 * mu - 5 * mu * nu
    k = k_zero - 8 * ray * radial_correction - p * ray * radial_correction / (cosine + 1)
    k += p_correction * (ray - 2)
    real_value += radius_squared * k + 2 * d * radius * (cosine - ray)
    return not real_value.contains(0)


def check_angle_anchor(lower, upper, rays):
    u = real((lower + upper) / 2)
    w = u.root(15)**16
    v = u.root(3)**4
    y = w.root(2)**3

    def quartic(x):
        return x**2 + 4 * u * x + 6 * v - 2

    def quintic(x):
        return x**2 + (5 * w - 1) * x + 10 * y - 5 * w - 1

    first, second, third = [2 * real(ray) for ray in rays]
    if not all((
        2 > first > second > third > -2,
        quartic(first) > 0,
        first + 2 * u > 0,
        quartic(second) < 0,
        quartic(third) < 0,
        quintic(first) < 0,
        quintic(second) < 0,
        quintic(third) > 0,
        third + (5 * w - 1) / 2 < 0,
        quartic(real(-2)) > 0,
        quintic(real(-2)) > 0,
        quartic(real(2)) > 0,
        quintic(real(2)) > 0,
    )):
        raise ArithmeticError(("angle anchor", lower, upper, rays))


def visit(coefficient_range, radius_cap, ray, exclude):
    stack = [(coefficient_range, (Fraction(0), Fraction(radius_cap)))]
    count = 0
    while stack:
        box = stack.pop()
        count += 1
        if exclude(span(*box[0]), span(*box[1]), ray):
            continue
        if count >= 10000:
            raise ArithmeticError(("unresolved ray box", box, ray))
        coefficient_width = box[0][1] - box[0][0]
        radius_width = box[1][1] - box[1][0]
        index = 0 if 10 * coefficient_width >= radius_width else 1
        lower, upper = box[index]
        midpoint = (lower + upper) / 2
        stack.append(box[:index] + ((lower, midpoint),) + box[index + 1:])
        stack.append(box[:index] + ((midpoint, upper),) + box[index + 1:])
    return count


def verify_high_angle_cover():
    lower = Fraction(13, 20)
    visited = 0
    slabs = 0
    for row in HIGH_RAYS.splitlines():
        upper_text, *ray_texts = row.split()
        upper = Fraction(upper_text)
        if not lower < upper <= Fraction(97, 100):
            raise ArithmeticError(("angle cover order", lower, upper))
        rays = tuple(Fraction(int(item), 10**12) for item in ray_texts)
        check_angle_anchor(lower, upper, rays)

        wl_value = real(lower).root(3)**4
        uq = real(upper)**2
        wh_value = (15 * uq / (16 - uq)).sqrt()
        # Round outward so the two-dimensional bisection uses fractions only.
        wl = Fraction(int((10**12 * wl_value.lower()).floor().unique_fmpz()), 10**12)
        wh = Fraction(int((10**12 * wh_value.upper()).ceil().unique_fmpz()), 10**12)
        if not Fraction(1, 2) <= wl <= wh < 1:
            raise ArithmeticError(("Taylor remainder domain", wl, wh))
        for ray in map(real, rays):
            visited += visit((lower, upper), 8, ray, quartic_ray_excluded)
            visited += visit((wl, wh), 16, ray, quintic_ray_excluded)
        slabs += 1
        lower = upper
    if lower != Fraction(97, 100):
        raise ArithmeticError(("uncovered angle range", lower))
    return slabs, visited



def outer_skew_cap(u):
    x = (1 / u**2 - 1) / 3
    return ((1 + x) * (1 + 2*x) / (1 + 3*x)).root(4)


def quintic_outer_skew_cap(w):
    x = (1 / w**2 - 1) / 4
    return (((1 + x) * (1 + 2*x) * (1 + 3*x))**2 / (1 + 4*x)**3).root(10)


def low_quartic_ray_excluded(u, radius, ray, cap):
    lower, upper = u.lower(), u.upper()
    center = span(
        lower.root(3)**4,
        upper * ((1 + 2*upper**2) / (2 + upper**2)).sqrt(),
    )
    square = radius**2
    inverse = 1 / square
    difference = -ray * (square - inverse) / (2*u)
    radical = (4 + difference**2).sqrt()
    skew = radius * (radical - difference) / 2
    if not skew.overlaps(cap):
        return True
    value = (square + inverse) * (2*ray**2 - 1) + 4*u*ray*radical + 6*center
    return not value.contains(0)


def low_quintic_ray_excluded(w, radius, outer_skew, ray):
    w_lower, w_upper = w.lower(), w.upper()
    b_lower, b_upper = outer_skew.lower(), outer_skew.upper()
    if quintic_outer_skew_cap(w_lower) < b_lower:
        return True
    x = (1 / w_lower**2 - 1) / 4
    link_cap = (((1 + 3*x)*(1 + x))**4 / ((1 + 2*x)**6*(1 + 4*x))).root(10)
    inner_cap = ((1 + 2*x)**4 / ((1 + x)*(1 + 3*x)*(1 + 4*x))).root(10)
    c_squared_lower = (b_lower / link_cap).max(1 / inner_cap**2)
    c_squared_upper = (b_upper * link_cap).min(inner_cap**2)
    if c_squared_lower > c_squared_upper:
        return True
    c_squared = span(c_squared_lower, c_squared_upper)
    inner_skew = c_squared.sqrt()
    center = span(
        w_lower.root(2)**3,
        w_upper * ((1 + 3*w_upper**2) / (3 + w_upper**2)).sqrt(),
    )
    if not all(value.is_finite() for value in (link_cap, inner_skew, center)):
        return False

    radius2, radius4, radius5 = (radius**exponent for exponent in (2, 4, 5))
    even = 8*ray**2 - 2
    linear = 4*ray
    forward = w * outer_skew
    backward = w / outer_skew
    plus = even*radius5 + 20*ray*backward*radius4 - 10*forward*radius - linear
    minus = -linear*radius5 - 10*backward*radius4 + 20*ray*forward*radius + even
    if plus > 0 or minus > 0:
        return True
    if not (minus + 20*center*radius2*inner_skew).contains(0):
        return True
    if not (plus + 20*center*radius2*radius/inner_skew).contains(0):
        return True
    if not (radius*minus - c_squared*plus).contains(0):
        return True
    return not (plus*minus).overlaps((20*center)**2 * radius5)


def outward_fraction(value, lower):
    denominator = 10**12
    if lower:
        endpoint = (denominator * value.lower()).floor()
    else:
        endpoint = (denominator * value.upper()).ceil()
    return Fraction(int(endpoint.unique_fmpz()), denominator)


def visit_boxes(box, exclude, weights, relative, limit):
    stack = [box]
    count = 0
    while stack:
        box = stack.pop()
        count += 1
        if exclude(*(span(*pair) for pair in box)):
            continue
        if count >= limit:
            raise ArithmeticError(("unresolved ray box", box))
        widths = [(upper - lower) * weight for (lower, upper), weight in zip(box, weights)]
        if relative:
            widths = [width / (lower + upper) for width, (lower, upper) in zip(widths, box)]
        index = max(range(len(box)), key=widths.__getitem__)
        lower, upper = box[index]
        midpoint = (lower + upper) / 2
        for half in ((lower, midpoint), (midpoint, upper)):
            stack.append(box[:index] + (half,) + box[index + 1:])
    return count


def verify_low_angle_cover():
    lower = Fraction(1, 500)
    visited = 0
    endpoint = Fraction(13, 20)
    rows = LOW_RAYS.splitlines()
    for row in rows:
        upper_text, *ray_texts = row.split()
        upper = Fraction(upper_text)
        if not lower < upper <= endpoint:
            raise ArithmeticError(("low cover order", lower, upper))
        rays = tuple(Fraction(int(item), 10**12) for item in ray_texts)
        check_angle_anchor(lower, upper, rays)

        quartic_cap = span(1, outer_skew_cap(real(lower)))
        w_lower = outward_fraction(real(lower).root(3)**4, lower=True)
        squared = real(upper)**2
        w_upper = outward_fraction((15*squared/(16-squared)).sqrt(), lower=False)
        b_upper = outward_fraction(quintic_outer_skew_cap(real(w_lower)), lower=False)
        quartic_box = ((lower, upper), (Fraction(1, 16), Fraction(16)))
        quintic_box = ((w_lower, w_upper), (Fraction(1, 32), Fraction(32)), (Fraction(1), b_upper))
        for ray in map(real, rays):
            visited += visit_boxes(
                quartic_box, lambda u, radius: low_quartic_ray_excluded(u, radius, ray, quartic_cap),
                (1, 3), relative=True, limit=200000,
            )
            visited += visit_boxes(
                quintic_box, lambda w, radius, b: low_quintic_ray_excluded(w, radius, b, ray),
                (1, 3, 1), relative=True, limit=200000,
            )
        lower = upper
    if lower != endpoint:
        raise ArithmeticError(("uncovered low angle range", lower))
    return len(rows), visited


def near_quartic_ray_excluded(scale, abscissa, alpha_scaled, ray):
    scale2 = scale**2
    scale4 = scale2**2
    u = 1 - scale2
    lower = u.lower()
    alpha_cap = ((1 + lower) / (3*lower**2))**2 / 2
    if alpha_scaled.lower() > alpha_cap:
        return True
    alpha = span(alpha_scaled.lower(), alpha_scaled.upper().min(alpha_cap))
    shape = span(Fraction(2, 9), remainder(lower, 2) if lower < real(".999999")
                 else 3*lower*(3 + 4*lower**2 - lower**4)
                 / ((2 + lower**2)**2 * (1 + 2*lower**2)
                    * ((2 + lower**2)*(1 + 2*lower**2)).sqrt()))
    sinh1 = span(1, "1.001")
    cosh2 = span(".5", ".501")
    scaled_root = acb(abscissa, ray * (1 - scale*abscissa))
    z = -1 + scale*scaled_root
    value = scaled_root**4 - 4*z*scaled_root**2 + 6*shape*z**2
    value += 4*u*alpha*sinh1*scale*z*(1-z)*scaled_root
    value += 4*u*alpha**2*cosh2*scale4*(z + z**3)
    return not value.contains(0)


def near_quintic_ray_excluded(scale, abscissa, curvature_ratio, beta_scaled, ray):
    scale2 = scale**2
    scale4 = scale2**2
    scale6 = scale2*scale4
    w = 1 - curvature_ratio*scale2
    lower = w.lower()
    shift = (1 + lower) / (4*lower**2)
    beta_cap = (curvature_ratio.upper()*shift)**2
    if beta_scaled.lower() > beta_cap:
        return True
    beta = span(beta_scaled.lower(), beta_scaled.upper().min(beta_cap))
    link_cap = real(".6")*(curvature_ratio.upper()*shift)**4
    link = span(-link_cap, link_cap)
    chi = (beta - link*scale4)/2
    shape = span(Fraction(3, 8), remainder(lower, 3) if lower < real(".999999")
                 else 12*lower*(3 + 6*lower**2 - lower**4)
                 / ((3 + lower**2)**2 * (1 + 3*lower**2)
                    * ((3 + lower**2)*(1 + 3*lower**2)).sqrt()))
    y = 1 - real("1.5")*curvature_ratio*scale2 + shape*curvature_ratio**2*scale4
    sinh1 = span(1, "1.001")
    cosh2 = span(".5", ".501")
    sinh3 = span(Fraction(1, 6), ".167")
    scaled_root = acb(abscissa, ray * (1 - scale*abscissa))
    z = -1 + scale*scaled_root

    correction = 10*(-curvature_ratio/2 + shape*curvature_ratio**2*scale2)*chi*sinh1
    correction += 5*w*(-link*scale2 + (2*chi**3*sinh3 - beta**3*sinh3)*scale6)
    value = scaled_root**5 - 5*curvature_ratio*z*scaled_root**3
    value += 10*curvature_ratio**2*shape*z**2*scaled_root
    value += 5*w*beta*sinh1*scale*z*(1-z)*scaled_root**2 + scale*correction*z**2*(1-z)
    value += scale4*(5*w*beta**2*cosh2*z*(1-z+z**2) + 10*y*chi**2*cosh2*z**2)*scaled_root
    return not value.contains(0)


def verify_near_flat_guards():
    maximum_scale = real(Fraction(1, 4))
    u = 1 - maximum_scale**2
    ratio = real(Fraction(4, 3))
    w = 1 - ratio*maximum_scale**2
    link_cap = real(".6")*(ratio*(1 + w)/(4*w**2))**4
    lower_u = 1 - real(".18")**2
    minimum_ratio = 16*real(15).sqrt() / ((16 - lower_u**2)*(16 - lower_u**2).sqrt())
    check = [
        ((1 + u)/(3*u**2))**2 / 2 < real(".3"),
        (ratio*(1 + w)/(4*w**2))**2 < real(Fraction(2, 3)),
        link_cap < real(".21"),
        (real(Fraction(2, 3)) + link_cap*maximum_scale**4)/2 < real(".334"),
        remainder(u, 2) < real(".3"),
        remainder(w, 3) < real(".5"),
        minimum_ratio > real("1.05"),
    ]
    a, b, c, link, fv, fy = map(real, (".3", Fraction(2, 3), ".334", ".21", ".3", ".5"))
    l = maximum_scale
    correction = 10*(ratio/2 + fy*ratio**2*l**2)*c*real("1.001")
    correction += 5*(link*l**2 + (2*c**3+b**3)*real(".167")*l**6)
    z, z_minus, z_quad, inverse_radius = map(real, (".30", ".35", ".1075", ".05"))
    fourth = 4*z*inverse_radius + 6*fv*z**2*inverse_radius**2
    fourth += 4*a*real("1.001")*l*z*z_minus*inverse_radius
    fourth += 4*a*a*real(".501")*l**4*(z*inverse_radius**2+z**3)*inverse_radius
    fifth = 5*ratio*z*inverse_radius + 10*ratio**2*fy*z**2*inverse_radius**2
    fifth += 5*b*real("1.001")*l*z*z_minus*inverse_radius
    fifth += l*correction*z**2*z_minus*inverse_radius**2
    fifth += l**4*real(".501")*(5*b*b*z*z_quad*inverse_radius
                              +10*c*c*z**2*inverse_radius**2)
    check += [correction < real("2.5"), fourth < real(".07"), fifth < real(".11")]
    check += [cap*l**4 < real(".01") for cap in (a, b, c)]
    # Power series tails with common geometric ratio at most 0.01^2.
    tail = real(".01")**2/(1-real(".01")**2)
    check += [1 + tail/6 < real("1.001"), real(".5") + tail/24 < real(".501"),
              real(Fraction(1, 6)) + tail/120 < real(".167")]
    p_min = 2 - real(Fraction(8, 3)).sqrt()
    p_max = 2 + real(Fraction(8, 3)).sqrt()
    q_min = (5 - real(10).sqrt()) / 2
    q_max = (5 + real(10).sqrt()) / 2
    check += [p_min < real(".8")**2 < real("1.05")*q_min,
              ratio*q_min < real("1.5")**2 < p_max,
              p_max < 4 < real("1.05")*q_max]
    if not all(check):
        raise ArithmeticError(("near-flat guards", check))


def verify_near_flat_cover():
    verify_near_flat_guards()
    scale = (Fraction(0), Fraction(9, 50))
    abscissa = (Fraction(-20), Fraction(20))
    alpha_scaled = (Fraction(0), Fraction(3, 10))
    beta_scaled = (Fraction(0), Fraction(2, 3))
    curvature_ratio = (Fraction(21, 20), Fraction(4, 3))
    visited = 0
    for ray in map(real, (2, Fraction(3, 2), Fraction(4, 5))):
        visited += visit_boxes(
            (scale, abscissa, alpha_scaled),
            lambda l, x, a: near_quartic_ray_excluded(l, x, a, ray),
            (1, 1, scale[1]), relative=False, limit=200000,
        )
        visited += visit_boxes(
            (scale, abscissa, curvature_ratio, beta_scaled),
            lambda l, x, r, b: near_quintic_ray_excluded(l, x, r, b, ray),
            (1, 1, 1, scale[1]), relative=False, limit=200000,
        )
    return visited


def verify_small_guards():
    a, b = map(Fraction, (".193", ".497"))
    u = Fraction(1, 500)
    root = Fraction(".045")
    outer = Fraction(".084")
    angle_bound = a/4 + a**3/(24*(1-a*a)) + b/5 + b**3/(30*(1-b*b))
    if not (
        root**2 > u and outer**5 > u**2
        and 4*(root+root**3)+6*u < a
        and 5*(outer+outer**4)+10*(outer**2+outer**3) < b
        and angle_bound < Fraction(".154") and real(".154") < arb.pi()/20
    ):
        raise ArithmeticError("small-mean bounds")


if __name__ == "__main__":
    ctx.prec = 110
    verify_small_guards()
    print("Low-angle slabs and boxes", verify_low_angle_cover(), flush=True)
    print("High-angle slabs and boxes", verify_high_angle_cover(), flush=True)
    print("Near-flat boxes", verify_near_flat_cover(), flush=True)
