"""Check Laurent, parameter and collision identities for sequence (10,4,11)."""

from functools import cache
from itertools import product

from sympy import Poly, Rational, binomial, diff, expand, factorial, resultant, symbols


def coefficient_factory(a, b, c):
    @cache
    def coefficient(power, degree):
        if degree < 0:
            return 0
        value = 0
        ranges = [range(degree // weight + 1) for weight in (2, 3, 4)]
        for exponents in product(*ranges):
            if sum(index*weight for index, weight in zip(exponents, (2, 3, 4))) != degree:
                continue
            total = sum(exponents)
            term = binomial(power, total)*factorial(total)
            for parameter, index in zip((a, b, c), exponents):
                term *= parameter**index/factorial(index)
            value += term
        return value
    return coefficient


def check_quartic_decic():
    z, a, b, delta, alpha, beta, p, kappa = symbols("z a b Delta alpha beta p kappa")
    a_y, b_y, delta_y, alpha_y, beta_y, p_y = symbols("a_y b_y Delta_y alpha_y beta_y p_y")
    c = (a**2-delta)/4
    quartic = z**4+a*z**2+b*z+c
    coefficient = coefficient_factory(a, b, c)

    def square_root_coefficient(k):
        return (coefficient(Rational(5, 2), k)
                + alpha*coefficient(Rational(1, 2), k-8)/2
                + beta*coefficient(-Rational(1, 2), k-12)/2)

    decic = expand(sum(square_root_coefficient(k)*z**(10-k) for k in range(11)) + p*z)
    e1, e2, e3 = (-square_root_coefficient(k) for k in (11, 12, 13))

    def dy(f):
        return (diff(f, a)*a_y + diff(f, b)*b_y + diff(f, delta)*delta_y
                + diff(f, alpha)*alpha_y + diff(f, beta)*beta_y + diff(f, p)*p_y)

    jacobian = Poly(diff(quartic, z)*dy(decic)-diff(decic, z)*dy(quartic), z)
    constant_parameters = {alpha_y: 0, p_y: 0}
    assert jacobian.degree() == 5
    assert expand(jacobian.nth(5)-2*alpha_y) == 0
    assert expand(jacobian.nth(4).subs(alpha_y, 0)-4*p_y) == 0
    assert expand(jacobian.nth(3).subs(constant_parameters)) == 0
    assert expand((jacobian.nth(2)-4*dy(e1)+p*a_y).subs(constant_parameters)) == 0
    assert expand((jacobian.nth(1)-4*dy(e2)-2*beta_y+p*b_y).subs(constant_parameters)) == 0

    first = b*(15*delta**2-40*a*b**2+64*alpha)/256+p*a/4
    second = -delta**3+24*delta*a*b**2-8*b**4-64*alpha*delta/5+256*p*b/5
    primitive = a*first-3*p*a**2/8+p*delta/4-5*delta*b**3/32
    assert expand(first+e1-p*a/4) == 0
    assert expand(5*second/1024+e2+beta/2-p*b/4) == 0
    assert expand((jacobian.nth(0)-dy(primitive)+a*dy(first)).subs(constant_parameters)) == 0
    first_residual = a*(40*b**3-64*p)-b*(15*delta**2+64*alpha)+256*kappa
    assert expand(first_residual+256*(first-kappa)) == 0

    total, gap, Y = symbols("s g Y")
    upper, lower = (total+gap)/2, (total-gap)/2
    quartic_difference = expand((quartic.subs(z, upper)-quartic.subs(z, lower))/gap)
    decic_difference = expand((decic.subs(z, upper)-decic.subs(z, lower))/gap)
    squared_gap = -total**2-2*a-2*b/total
    assert expand(quartic_difference.subs(gap**2, squared_gap)) == 0
    powers = {gap**i: squared_gap**(i//2) for i in (8, 6, 4, 2)}
    sum_poly = Poly(expand(128*decic_difference.subs(powers, simultaneous=True)), total)
    assert sum_poly.degree() == 9
    assert sum_poly.LC() == 8
    assert sum_poly.nth(0) == 128*p

    alpha_value = (-15*delta**2+40*a*b**2+(256*kappa-64*p*a)/b)/64
    sum_equation = (2*total**4*(total**2+2*a)**2
                    + 5*delta*total**2*(total**2+2*a)-2*b*total**3*(total**2+2*a)
                    - 2*b**2*total**2-5*delta*b*total+(64*kappa-16*p*a)/b+32*p/total)
    assert expand(sum_poly.as_expr().subs(alpha, alpha_value)-4*total*sum_equation) == 0
    gap_equation = (2*total**4*Y**2+(5*delta*total**2-10*b*total**3)*Y
                    + 10*b**2*total**2-15*delta*b*total+(64*kappa-16*p*a)/b+32*p/total)
    assert expand(gap_equation.subs(Y, -squared_gap)-sum_equation) == 0

    limiting_shape = {a: Rational(3, 8), b: 1, delta: 1, alpha: 0, p: 0}
    quartic_derivative = diff(quartic, z).subs(limiting_shape)
    decic_derivative = diff(decic, z).subs(limiting_shape)
    assert resultant(quartic_derivative, decic_derivative, z) == 208000


if __name__ == "__main__":
    check_quartic_decic()
    print("Verified the conductor-eighteen Laurent, parameter and collision identities.")
