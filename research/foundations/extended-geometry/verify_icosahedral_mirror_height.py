"""Exact checks of the icosahedral mirror-height formulas."""

from math import comb
import sympy as sp


q, p, r, w, u = sp.symbols("q p r w u")
h = r - 7*q**2*p
f = sp.expand(
    -64*q**9*p**2 + 2240*q**6*p**3 - 21760*q**3*p**4 + 8640*p**5
    + 16*q**7*p*h + 80*q**4*p**2*h - 3600*q*p**3*h
    + 4*q**5*h**2 - 200*q**2*p*h**2 - 25*h**3
)
a = 45*p**2 + 13*p*q**3 + q**6 - 5*q*r
b = 19*p**2*q**2 + 3*p*q**5 + (15*p - q**3)*r
f1, a1, b1 = (poly.subs(q, 1) for poly in (f, a, b))


def check_fibre():
    metric = sp.Matrix([
        [4*q, 12*p, 20*r],
        [12*p, h, 3*q**2*r - 5*q**4*p - 40*q*p**2],
        [20*r, 3*q**2*r - 5*q**4*p - 40*q*p**2,
         (25*q**4 + 240*q*p)*r - 15*p*(q**3+8*p)**2],
    ])
    gradient = sp.Matrix([f.diff(x) for x in (q, p, r)])
    assert sp.expand((gradient.T*metric*gradient)[0] + 320*a*b*f) == 0

    z = 5*u + 15*p - 1
    g = 225*f1 + 80*w*a1*b1
    eq = (2880*p**2 + 240*u*(w-3)*p - 225*u**3
          + (80*w-180)*u**2 - 16*w*u)
    identity = sp.cancel(z*g - 5*a1**2*eq)
    assert sp.cancel(identity.subs(u, b1/a1)) == 0
    d = w**2 - 22*w + 45
    assert sp.expand((24*p+(w-3)*u)**2
                     - u*(45*u**2+d*u+16*w/5) - eq/5) == 0
    r_new = 3*p*(1+u)-u**2-u+16*w*u/45
    assert sp.rem(z*(45*p**2+13*p+1-5*r_new)
                  - (5*p+1)**2*(27*p-1), eq, p) == 0


def check_sheet_poles():
    s = sp.symbols("s")
    old_r = (45*p**2+13*p+1-s)/5
    x = -b1.subs(r, old_r)/s
    wx = 45*f1.subs(r, old_r)/(16*s**2)
    new_p = sp.expand((wx-3*x)/12-p)
    new_r = sp.expand(3*new_p*(1-x)-x**2+x-16*wx/45)
    c = (5*p+1)**2*(27*p-1)
    assert sp.cancel(s**2*new_p).subs(s, 0).expand() == (-3*c**2/320).expand()
    assert sp.cancel(s**3*new_r).subs(s, 0).expand() == (-9*c**3/1600).expand()
    convergent = sp.cancel(s**5*(5*new_r**2+192*new_p**3))
    assert convergent.subs(s, 0).expand() == (
        27*c**5*(45*p+31)/512000
    ).expand()


def conditional_moment(i, j, k):
    x, g = sp.symbols("x g")
    p_chart = ((w-3)*x+g)/24
    r_chart = 3*p_chart*(1-x)-x**2+x-16*w*x/45
    a_chart = 45*p_chart**2+13*p_chart+1-5*r_chart
    f_chart = 16*w*a_chart**2*x/45
    d = w**2 - 22*w + 45
    norm = d**2 - 576*w
    phase = sp.Poly(p_chart**i*r_chart**j*f_chart**k, g, x)
    square = x*(-45*x**2+d*x-16*w/5)
    even = sp.Poly(sum(
        coefficient * square**(g_exp//2)*x**x_exp
        for (g_exp, x_exp), coefficient in phase.terms()
        if g_exp % 2 == 0
    ), x)
    result = sp.Poly(0, w)
    for (degree,), coefficient in even.terms():
        arcsine = sum(
            sp.Rational(comb(degree, 2*l)*comb(2*l, l), 4**l*90**degree)
            * d**(degree-2*l)*norm**l
            for l in range(degree//2+1)
        )
        result += sp.Poly(coefficient*arcsine, w)
    return sum(value/sp.Integer(2*degree+1)
               for (degree,), value in result.terms())


def check_moments():
    coords = sp.symbols("s1 s2 s3")
    s1, s2, s3 = coords
    root5 = sp.sqrt(5)
    field = sp.QQ.algebraic_field(root5)
    kappa = (3+root5)/2
    form = (s2-kappa*s3)*(s1-kappa*s2)*(kappa*s1-s3)/(2+root5)
    p0 = sp.Poly(form, *coords, domain=field)
    q0 = sp.Poly(s1+s2+s3, *coords, domain=field)
    h0 = 4*sum(sp.Poly(coord, *coords, domain=field)*p0.diff(coord)**2
               for coord in coords)
    r0 = h0+7*q0**2*p0
    f0 = sp.Poly(0, *coords, domain=field)
    for (q_exp,p_exp,r_exp), coefficient in sp.Poly(f,q,p,r).terms():
        f0 += coefficient * q0**q_exp * p0**p_exp * r0**r_exp

    for i, j, k in [(1,0,0), (0,1,0), (2,0,0), (1,1,0), (0,2,0),
                    (3,0,0), (0,0,1), (1,0,1)]:
        form = p0**i*r0**j*f0**k
        direct = field.zero
        for exponents, value in form.terms():
            factor = (sp.prod(sp.rf(sp.S.Half, exponent) for exponent in exponents)
                      /sp.rf(sp.Rational(3,2), sum(exponents)))
            direct += field.convert(factor*value)
        conditional = conditional_moment(i, j, k)
        assert direct == field.convert(conditional), (i,j,k, direct, conditional)
        print((i, j, k), conditional)


if __name__ == "__main__":
    check_fibre()
    check_sheet_poles()
    check_moments()
    print("PASS: mirror chart, sheet poles and independent sphere moments.")
