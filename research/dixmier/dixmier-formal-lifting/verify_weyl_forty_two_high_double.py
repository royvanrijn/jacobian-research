"""Check the two leading double-point computations at a high root."""

import sympy as sp


def first_face():
    z, w = sp.symbols("z w")
    q = (5 + sp.sqrt(15) * sp.I) / 4
    field = sp.QQ.algebraic_field(q)
    d = (5 - q) / 3
    r = z * (z - 1)**4 * (z - q)
    a = z * (z - 1) * (z - q) * (z - d)
    identity = 3*r*sp.diff(a, z) - 2*a*sp.diff(r, z) + q*d*r
    assert sp.Poly(sp.expand(identity), z, domain=field).is_zero

    dr = sp.Poly(sp.cancel((r.subs(z, w)-r)/(w-z)), w, z, domain=field)
    da = sp.Poly(sp.cancel((a.subs(z, w)-a)/(w-z)), w, z, domain=field)
    resultant = sp.Poly(dr.resultant(da).as_expr(), z, domain=field)
    sextic = sp.Poly(
        z**6 - 6*z**5 + (151-q)*z**4/10 + (10*q-302)*z**3/15
        + (1127-93*q)*z**2/80 + (139*q-953)*z/240 + 7*(q-3)/240,
        z, domain=field,
    )
    zero_factor = sp.Poly(z*(z-1)*(z-q), z, domain=field)
    assert resultant.monic() == zero_factor**2 * sextic
    assert sextic.gcd(sextic.diff()).degree() == 0
    assert sextic.gcd(zero_factor).degree() == 0
    determinant = dr.diff(z)*da.diff(w)-dr.diff(w)*da.diff(z)
    basis = sp.groebner(
        [dr.as_expr(), da.as_expr(), determinant.as_expr()], w, z, domain=field
    )
    assert list(basis) == [1]
    assert field.from_sympy(2*q**2-5*q+5) == field.zero
    assert field.from_sympy((q-1)**4-1) != field.zero


def scalar_face():
    rho, S = sp.symbols("rho S")
    G = 35+280*S+560*S**2+448*S**3+128*S**4
    assert sp.expand((1+8*S)*G-2*S*(1+S)*sp.diff(G,S)) == 35
    H = sum(rho**j for j in range(8))
    seventh = sum(rho**j for j in range(7))
    reciprocal = 5*rho**6-4*rho**5+rho**4+6*rho**3+rho**2-4*rho+5
    numerator = sum(
        coefficient*(-1)**degree*(rho**(7*degree)-1)*H**(4-degree)
        for degree, coefficient in enumerate((35,280,560,448,128))
    )
    factors = -8*(rho-1)*(rho+1)**5*(rho**2+1)**5*seventh*reciprocal
    assert sp.expand(numerator-factors) == 0
    assert sp.gcd(reciprocal, sp.diff(reciprocal, rho)) == 1
    assert sp.gcd(reciprocal, rho*H*(rho**7-1)) == 1
    assert sp.gcd(seventh, H) == 1

    cubic = sp.expand(G**2-128**2*S*(1+S)**7)
    assert cubic == 896*S**3+2912*S**2+3216*S+1225
    assert sp.gcd(cubic, sp.diff(cubic,S)) == 1
    assert sp.gcd(cubic, S*(1+S)*G) == 1


def main():
    first_face()
    scalar_face()
    assert 6+84+84+12+12 == 198
    assert 6*(-3)+84*sp.Rational(-16,7)+84*(-3)+12*(-3) == -498
    print("Verified the first-face sextic, the twelve scalar ratios and the three scalar matching values.")


if __name__=="__main__":
    main()
