"""Check the Laurent, parameter and collision identities for (10,4,9)."""

from functools import cache
from itertools import product

from sympy import Poly, Rational, binomial, diff, expand, factorial, resultant, symbols, together


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
    z, a, b, delta, alpha, beta, kappa = symbols("z a b Delta alpha beta kappa")
    a_y, b_y, delta_y, alpha_y, beta_y = symbols("a_y b_y Delta_y alpha_y beta_y")
    c = (a**2-delta)/4
    quartic = z**4+a*z**2+b*z+c
    coefficient = coefficient_factory(a, b, c)

    def square_root_coefficient(k):
        return (coefficient(Rational(5, 2), k)
                + alpha*coefficient(Rational(1, 2), k-8)/2
                + beta*coefficient(-Rational(1, 2), k-12)/2)

    decic = expand(sum(square_root_coefficient(k)*z**(10-k) for k in range(11)))
    e1, e2, e3 = (-square_root_coefficient(k) for k in (11, 12, 13))

    def dy(f):
        return (diff(f, a)*a_y + diff(f, b)*b_y + diff(f, delta)*delta_y
                + diff(f, alpha)*alpha_y + diff(f, beta)*beta_y)

    relative_jacobian = Poly(diff(quartic, z)*dy(decic)-diff(decic, z)*dy(quartic), z)
    fixed_alpha = {alpha_y: 0}
    assert relative_jacobian.degree() == 5
    assert expand(relative_jacobian.nth(5)-2*alpha_y) == 0
    assert expand(relative_jacobian.nth(3).subs(fixed_alpha)) == 0
    assert expand((relative_jacobian.nth(2)-4*dy(e1)).subs(fixed_alpha)) == 0
    assert expand((relative_jacobian.nth(1)-4*dy(e2)-2*beta_y).subs(fixed_alpha)) == 0
    constant = 4*dy(e3)+2*a*dy(e1)+a_y*e1
    assert expand((relative_jacobian.nth(0)-constant).subs(fixed_alpha)) == 0

    assert expand(-e1-b*(15*delta**2-40*a*b**2+64*alpha)/256) == 0
    expected_second = 5*(-delta**3+24*delta*a*b**2-8*b**4)/1024-alpha*delta/16
    assert expand(-e2-beta/2-expected_second) == 0
    assert expand(4*e3+2*a*e1+5*delta*b**3/32) == 0
    a_value = (15*delta**2+64*alpha-256*kappa/b)/(40*b**2)
    invariant = delta**3+16*alpha*delta/5-96*kappa*delta/(5*b)-b**4
    assert expand(expected_second.subs(a, a_value)-5*invariant/128) == 0

    total, gap, Y = symbols("s g Y")
    upper, lower = (total+gap)/2, (total-gap)/2
    quartic_difference = expand((quartic.subs(z, upper)-quartic.subs(z, lower))/gap)
    decic_difference = expand((decic.subs(z, upper)-decic.subs(z, lower))/gap)
    squared_gap = -total**2-2*a-2*b/total
    assert expand(quartic_difference.subs(gap**2, squared_gap)) == 0
    alpha_value = (-15*delta**2+40*a*b**2+256*kappa/b)/64
    sum_equation = (2*total**4*(total**2+2*a)**2+5*delta*total**2*(total**2+2*a)
                    -2*b*total**3*(total**2+2*a)-2*b**2*total**2-5*delta*b*total+64*kappa/b)
    replaced = {gap**power: squared_gap**(power//2) for power in (8, 6, 4, 2)}
    decic_reduced = decic_difference.subs(replaced, simultaneous=True).subs(alpha, alpha_value)
    assert expand(together(32*decic_reduced-total*sum_equation).as_numer_denom()[0]) == 0
    gap_equation = (2*total**4*Y**2+(5*delta*total**2-10*b*total**3)*Y
                    +10*b**2*total**2-15*delta*b*total+64*kappa/b)
    assert expand(gap_equation.subs(Y, -squared_gap)-sum_equation) == 0

    leading = {a: Rational(3, 8), b: 1, delta: 1, alpha: 0}
    assert resultant(diff(quartic, z).subs(leading), diff(decic, z).subs(leading), z) == 208000


if __name__ == "__main__":
    check_quartic_decic()
    print("Verified the conductor-sixteen Laurent, parameter and collision identities.")
