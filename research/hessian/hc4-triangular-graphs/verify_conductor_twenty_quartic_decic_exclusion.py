"""Verify the parameter identities and collision equations for (10,4,13)."""

from itertools import product
from sympy import Poly, Rational, binomial, diff, expand, factorial, resultant, symbols


def monic_binomial_coefficients(parameters, power, maximum):
    weights = (2, 3, 4)
    coefficients = [0] * (maximum+1)
    for orders in product(*(range(maximum//weight+1) for weight in weights)):
        degree = sum(weight*order for weight, order in zip(weights, orders))
        if degree > maximum:
            continue
        total = sum(orders)
        term = binomial(power, total)*factorial(total)
        for parameter, order in zip(parameters, orders):
            term *= parameter**order/factorial(order)
        coefficients[degree] += term
    return coefficients


def check_quartic_decic_identities():
    z, a, b, delta, p, q, r, Q, R = symbols("z a b Delta p q r Q R")
    a_y, b_y, delta_y, p_y, q_y, r_y = symbols("a_y b_y Delta_y p_y q_y r_y")
    c = (a**2-delta)/4
    quartic = z**4+a*z**2+b*z+c
    coefficients = monic_binomial_coefficients((a, b, c), Rational(5, 2), 10)
    decic_base = sum(value*z**(10-index) for index, value in enumerate(coefficients))
    decic = decic_base+p*z**3+q*z**2+r*z

    def dy(f):
        return (diff(f, a)*a_y+diff(f, b)*b_y+diff(f, delta)*delta_y
                + diff(f, p)*p_y+diff(f, q)*q_y+diff(f, r)*r_y)

    jacobian = Poly(diff(quartic, z)*dy(decic)-diff(decic, z)*dy(quartic), z)
    assert jacobian.degree() == 6
    assert expand(jacobian.nth(6)-4*p_y) == 0
    assert expand(jacobian.nth(5).subs(p_y, 0)-4*q_y) == 0
    assert expand(jacobian.nth(4).subs({p_y: 0,q_y: 0})-4*r_y+3*p*a_y) == 0

    first = -(15*delta**2*b-48*delta*p+128*Q*b+64*R*a+24*a**2*p-40*a*b**3)/64
    second = (5*delta**3+128*delta*Q-120*delta*a*b**2-256*R*b+192*a*b*p+40*b**4)/256
    primitive = (8*delta*R+6*delta*a*p-5*delta*b**3-12*R*a**2-5*a**3*p+12*b**2*p)/32-a*first/4
    scaled = {p_y: 0, q_y: 0, r_y: 3*p*a_y/4, q: Q, r: R+3*p*a/4}
    adjusted = Poly(jacobian.as_expr().subs(scaled)
                    - jacobian.nth(3).subs(scaled)*diff(quartic, z)/4,z)
    assert adjusted.degree() == 2
    assert expand(adjusted.nth(2)-dy(first).subs(p_y, 0)) == 0
    assert expand(adjusted.nth(1)-dy(second).subs(p_y, 0)) == 0
    assert expand(adjusted.nth(0)-dy(primitive).subs(p_y, 0)-a*dy(first).subs(p_y, 0)/4) == 0

    s, g, Y = symbols("s g Y")
    upper, lower = (s+g)/2, (s-g)/2
    decic = decic.subs({q: Q,r: R+3*p*a/4})
    quartic_difference = expand((quartic.subs(z, upper)-quartic.subs(z, lower))/g)
    decic_difference = expand((decic.subs(z, upper)-decic.subs(z, lower))/g)
    squared_gap = -s**2-2*a-2*b/s
    assert expand(quartic_difference.subs(g**2, squared_gap)) == 0
    replaced_powers = {g**i: squared_gap**(i//2) for i in (8, 6, 4, 2)}
    sum_poly = Poly(expand(128*s*decic_difference.subs(replaced_powers, simultaneous=True)),s)
    last_terms = (15*delta**2+128*Q-40*a*b**2)*s**2+(128*R+32*a*p)*s-64*b*p
    factored = (8*s**6*(s**2+2*a)**2+20*delta*s**4*(s**2+2*a)-8*b*s**5*(s**2+2*a)
                -8*b**2*s**4+(64*p-20*delta*b)*s**3+last_terms)
    assert expand(sum_poly.as_expr()-factored) == 0
    assert sum_poly.degree() == 10 and sum_poly.LC() == 8 and sum_poly.nth(0) == -64*b*p
    gap_equation = (8*s**6*Y**2+(20*delta*s**4-40*b*s**5)*Y+40*b**2*s**4
                    +(64*p-60*delta*b)*s**3+last_terms)
    assert expand(gap_equation.subs(Y,-squared_gap)-factored) == 0

    leading = {a: Rational(3, 8), b: 1, delta: 1}
    assert resultant(diff(quartic, z).subs(leading),diff(decic_base,z).subs(leading),z) == 208000


if __name__ == "__main__":
    check_quartic_decic_identities()
    print("Verified the conductor-twenty quartic–decic parameter and collision identities.")
