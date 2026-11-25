"""Certify the numerical thresholds in the uniform 4/7 proof."""

from fractions import Fraction as Q
from itertools import combinations
from math import isqrt


GRID = 10**12


def root_lower(value):
    value = Q(value)
    return Q(isqrt(value.numerator * GRID**2 // value.denominator), GRID)


def squared_difference(a, b):
    a, b = Q(a), Q(b)
    if a <= b:
        return Q(0)
    return (root_lower(a) + Q(1, GRID) - root_lower(b))**2


def check(name, bound, cap):
    if bound >= Q(cap):
        raise ValueError(f"{name}: expected {bound} < {cap}")
    print(f"{name}: < {Q(cap)}")


def polygon_vertices(i, j, k, u, h, x):
    """Vertices on both sides of the affine split in the mass bound."""
    limit = 1 - k * u
    split = 1 - (k - 1) * u - x
    lines = ((1, 0, u), (1, 0, h), (0, 1, u), (0, 1, h),
             (i, j, limit), (i - 1, j - 1, split))
    for (a, b, r), (c, d, s) in combinations(lines, 2):
        determinant = a * d - b * c
        if determinant:
            first = (r * d - b * s) / determinant
            second = (a * s - c * r) / determinant
            if u <= first <= h and u <= second <= h and i * first + j * second <= limit:
                yield first, second


def main():
    d0, x, z = Q(4, 7), Q(3, 7), Q(1, 7)
    fourier = 4 - 6 * d0
    c0, h, u = Q(2, 7), Q(".152"), Q(".09")
    half = Q(1, 2)

    def p(y, w):
        return squared_difference(2 * y, w)

    def g(m, w):
        return squared_difference(m + 1, m * w / x)

    def G(b):
        return b + p(x, 4 * b)

    if not ((3+x)/4-c0 == 1-z-c0 == x+c0/2 == d0):
        raise ValueError("rational endpoint identities")
    check("tall degree one",
          max((1+x)/2-h+p(x,Q(".75")), (1+x)/2-Q(".25")+p(x,1-x)),
          ".56587")
    check("short degree one with short height >= 1/2", max(
        (1+mass)/2-a+(p(x,half-a)+p(1-mass-x,half-a))/2
        for a in (c0,x) for mass in (a,x)
    ), ".539")
    check("short degree one and early tall coordinate",
          max(1-c0-h, half-c0+2*Q(".175")), d0)
    if not (1-3*Q(".175") >= x and 1-3*Q(".161") >= 6*x/5 and
            x+(1-3*Q(".175"))/2 > 2*c0):
        raise ValueError("early moment interpolation thresholds")
    for m, lower, upper in ((3,h,Q(".161")), (2,Q(".161"),Q(".175"))):
        bounds = []
        for b in (lower, upper):
            w = 1-3*b
            penalty = max(x*g(2,w),
                          2*max(Q(0),c0-w/2)*(1-c0/x)+(x-w/2)*g(m,w))
            bounds.append((1+x)/2-b+penalty/2)
        check(f"early coordinate, restricted moments m={m}", max(bounds), ".56705")
    check("early coordinate, unrestricted moments", max(
        (1+x)/2-b+p(x,1-3*b)/2 for b in (Q(".175"),Q(11,42))
    ), ".568")

    if not (h > 8*x/25 and Q(".192") > h and 3*x+h > 2-d0 and 9*h > 1):
        raise ValueError("tall-coordinate thresholds")
    check("two tall G terms when mass sum <= .42", G(Q(".192"))+G(Q(".228")), ".42245")
    check("two tall masses if one chosen mass exceeds .228",
          1-3*Q(".228")+2*G(Q(".192"))-z, ".56204")
    check("narrow two-tall case half the mass sum", G(Q(".25"))+G(x-Q(".25"))-z, ".29223")
    tau = (4*d0-2)/3
    check("short degree at least four singleton",
          max(Q(".29223")+a+p(x,4*a) for a in (tau,Q(".25"))), ".54224")
    check("short degree three singleton",
          max(Q(".29223")+a+p(x,3*a) for a in (Q(".08"),fourier-Q(".42"))),
          ".56226")
    check("short excess is insufficient", x+fourier+3-d0-7*Q(".42"), ".489")

    check("three large coordinates, highest at least six",
          half+h-3*u/2+p(x,6*u), ".55348")
    for i,j,k in combinations(range(1,6),3):
        values = []
        for a,b in polygon_vertices(i,j,k,u,h,x):
            mass = min(x, 1-(i-1)*a-(j-1)*b-(k-1)*u)
            values.append(half+mass/2-a-b+p(x,1-i*a-j*b))
        if values:
            check(f"three large coordinates of degrees {i}, {j}, {k}", max(values), ".56822")

    g2,g3 = g(2,half),g(3,half)
    short_penalty = max(
        x*g2,
        Q(".06")*(1-Q(".28")/x)+(x-Q(".25"))*g2,
        (x-Q(".25"))**2/(2*x)+(x-Q(".28"))*g3,
    )
    tall_penalty = max(x*g3, x*(3-1/x)**2/24+(x-Q(1,6))*g3)
    check("balanced short penalty", short_penalty, ".03969")
    check("balanced tall penalty", tall_penalty, ".01231")
    check("balanced two bins", half+u/2+(short_penalty+tall_penalty)/2, d0)
    check("unbalanced two bins", Q(".75")-(d0-2*h+u)/2, ".57129")


if __name__ == "__main__":
    main()

