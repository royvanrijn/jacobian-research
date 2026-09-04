"""Exact, outcome-independent model and binary-quartic size transformations.

Run under Sage Python. No factorization of quartic discriminants, point search,
rank labels, or known exceptional points enter coordinate selection.
"""
from __future__ import annotations

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, gcd, lcm

R = PolynomialRing(QQ, "z")
z = R.gen()


def qstr(x):
    return str(QQ(x))


def bits(x):
    x = QQ(x)
    return max(abs(x.numerator()).nbits(), x.denominator().nbits())


def size(values):
    heights = [int(bits(v)) for v in values]
    return {"maximum_bits": max(heights), "total_bits": sum(heights)}


def point_record(p):
    return {"infinity": True} if not p else dict(zip(("x", "y"), map(qstr, p.xy())))


def read_point(E, p):
    return E(0) if p.get("infinity") else E(QQ(p["x"]), QQ(p["y"]))


def map_record(phi):
    # Sage's convention: old x = u^2 new x+r;
    # old y = u^3 new y+s*u^2 new x+t.
    return {"u_r_s_t": list(map(qstr, phi.tuple())),
            "convention": "source_x=u^2*target_x+r; source_y=u^3*target_y+s*u^2*target_x+t"}


def transport(E, F, points):
    phi = E.isomorphism_to(F)
    inverse = ~phi
    transported = [phi(p) for p in points]
    if any(inverse(q) != p or q not in F for p, q in zip(points, transported)):
        raise ArithmeticError("model/section round trip failed")
    return phi, transported


def pointed(E, p):
    if not p or E.a1() or E.a3():
        raise ArithmeticError("expected a finite point on a completed-square model")
    x, y = p.xy()
    a2, a4 = E.a2(), E.a4()
    return R([a2*a2-2*a2*x-3*x*x-4*a4, -8*y, -6*x-2*a2, 0, 1])


def normalize(f):
    """Return integral g = ordinate_scale^2*f, removing exact square content.

    Only valuations at fixed small primes and an exact square-root test are
    used. Non-square residual content is retained; the ordinate scale is
    rational, not incorrectly required to be integral.
    """
    den = ZZ(lcm([v.denominator() for v in f]))
    g = R(f * den**2)
    content = abs(gcd([ZZ(v) for v in g]))
    if not content:
        raise ArithmeticError("zero quartic")
    residual, root = content, ZZ(1)
    # Testing only whether the WHOLE content is square loses an arbitrarily
    # large square factor in content = 2*k^2. Peel fixed small-prime valuations
    # first, preserving the nonsquare factor in the resulting polynomial.
    for prime in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        exponent = residual.valuation(prime)
        root *= ZZ(prime)**(exponent//2)
        residual //= ZZ(prime)**exponent
    rest_root = residual.isqrt()
    if rest_root**2 == residual:
        root *= rest_root
    scale = QQ(den)/root
    g /= root**2
    if any(v.denominator() != 1 for v in g) or g != f*scale**2:
        raise ArithmeticError("integral normalization identity failed")
    return g, scale


def binary_transform(f, matrix):
    a, b, c, d = map(QQ, matrix)
    if a*d == b*c:
        raise ArithmeticError("singular binary coordinate change")
    return sum((f[i]*(a*z+b)**i*(c*z+d)**(4-i) for i in range(5)), R(0))


def quartic_j(f):
    e, d, c, b, a = [f[i] for i in range(5)]
    I = 12*a*e-3*b*d+c*c
    J = 72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
    return 6912*I**3/(4*I**3-J**2)


def integral_quartic_bit_lower_bound(j):
    """Necessary coefficient bound for ANY integral binary quartic with this j.

    For |a,...,e| < 2^b, |I| <= 16 H^2 and |J| <= 137 H^3.
    Reduced j numerator < 2^(25+6b), denominator < 2^(16+6b).
    This bounds integral coefficients, not rational-coordinate height boxes.
    """
    j = QQ(j)
    return max(1, int((abs(j.numerator()).nbits()-25+5)//6),
               int((j.denominator().nbits()-16+5)//6))


def transform_menu(E, p, points):
    """Fixed finite menu: scale, translate, invert, and three-section PGL2.

    All slopes come from the supplied sixteen generic sections. Their bit
    sizes, never search responses or exceptional points, determine the menu.
    """
    x, y = p.xy()
    exponent = int((abs(x.numerator()).nbits()-x.denominator().nbits())//2)
    scales = [QQ(1), QQ(1)/x.denominator().isqrt()]
    scales += [QQ(2)**(exponent+i) for i in (-1, 0, 1)]
    slopes = sorted({(q[1]+y)/(q[0]-x) for q in points if q and q[0] != x},
                    key=lambda v: (bits(v), v))
    shifts = [QQ(0)] + slopes[:2]
    for k, scale in enumerate(dict.fromkeys(scales)):
        for h, shift in enumerate(shifts):
            yield f"affine-scale{k}-shift{h}", (scale, shift, QQ(0), QQ(1))
            yield f"reciprocal-scale{k}-shift{h}", (shift, scale, QQ(1), QQ(0))
    if len(slopes) >= 3:
        v0, v1, vi = slopes[:3]
        c = (v1-v0)/(vi-v1)
        yield "three-generic-section-slopes", (vi*c, v0, c, QQ(1))


def select_chart(E, p, points):
    f = pointed(E, p)
    trials = []
    selected = None
    for name, matrix in transform_menu(E, p, points):
        g, scale = normalize(binary_transform(f, matrix))
        metric = size([g[i] for i in range(5)])
        trial = {"name": name, **metric}
        trials.append(trial)
        score = (metric["maximum_bits"], metric["total_bits"], name)
        if selected is None or score < selected[0]:
            selected = (score, {"name": name, "matrix_a_b_c_d": list(map(qstr, matrix)),
                "ordinate_scale": qstr(scale),
                "integral_coefficients_ascending": [qstr(g[i]) for i in range(5)],
                **metric})
    record = selected[1]
    g = R(record["integral_coefficients_ascending"])
    if quartic_j(g) != E.j_invariant():
        raise ArithmeticError("selected binary quartic changed j")
    record["trials"] = trials
    return record


def quartic_point_to_source(record, coordinate, ordinate, E, base_point, phi):
    """Replay GL2, ordinate scaling, pointed map, then Weierstrass inverse.

    A GL2 pole maps to one of the two original points at infinity, whose
    images are the elliptic origin and Q. Return None there; callers record
    this explicitly and do not claim an affine point or new direction.
    """
    t, v = QQ(coordinate), QQ(ordinate)
    g = R(record["integral_coefficients_ascending"])
    if v*v != g(t):
        raise ArithmeticError("search point is off selected quartic")
    a, b, c, d = map(QQ, record["matrix_a_b_c_d"])
    denominator = c*t+d
    if not denominator:
        return None
    raw_t = (a*t+b)/denominator
    raw_v = v/QQ(record["ordinate_scale"])/denominator**2
    qx, qy = base_point.xy()
    x = (raw_t**2-qx-E.a2()+raw_v)/2
    y = raw_t*(x-qx)-qy
    p = E(x, y)
    answer = (~phi)(p)
    if phi(answer) != p:
        raise ArithmeticError("point transport failed")
    return answer
