"""Exact identities used for the global octic conductor exclusion."""

from sympy import Poly, Rational, diff, expand, symbols


def check_parameter_equations():
    z, a, b, c, a_y, b_y, c_y = symbols("z a b c a_y b_y c_y")
    quartic = z**4 + a*z**2 + b*z + c
    delta = a**2 - 4*c
    quadratic = z**2 + a/2
    sextic = 3*quartic*quadratic/2 - quadratic**3/2 + 3*b**2/8
    square_difference = Poly(sextic**2 - quartic**3, z)
    assert square_difference.degree() == 5
    assert expand(square_difference.nth(5) - 3*b*delta/8) == 0

    def dy(polynomial):
        return diff(polynomial, a)*a_y + diff(polynomial, b)*b_y + diff(polynomial, c)*c_y

    first = b*delta
    second = delta**2-8*a*b**2
    primitive = b**3/4-3*first*a/16
    jacobian = Poly(diff(quartic, z)*dy(sextic) - dy(quartic)*diff(sextic, z), z)
    assert jacobian.degree() == 2
    assert expand(jacobian.nth(2) - 3*dy(first)/4) == 0
    assert expand(jacobian.nth(1) + 3*dy(second)/32) == 0
    assert expand(jacobian.nth(0) - dy(primitive) - 3*a*dy(first)/16) == 0
    relative = Poly(3*diff(quartic, z)*sextic - 2*quartic*diff(sextic, z), z)
    assert relative.degree() <= 2


def check_collision_elimination():
    total, gap, a, b, delta = symbols("s g a b Delta")
    lower = (total-gap)/2
    higher = (total+gap)/2
    lower_d = lower**2 + a/2
    higher_d = higher**2 + a/2
    lower_u = lower_d**2 + b*lower-delta/4
    higher_u = higher_d**2 + b*higher-delta/4
    middle_u = (lower_u+higher_u)/2

    quartic_condition = total*(higher_d+lower_d)+b
    assert expand((higher_u-lower_u) - gap*quartic_condition) == 0
    sextic_condition = 3*middle_u - (higher_d**2+higher_d*lower_d+lower_d**2)
    gap_residual = total**2*gap**2 - Rational(3, 2)*(delta-2*b*total)
    assert expand(sextic_condition - gap_residual/2) == 0
    sum_equation = total**4 + 2*a*total**2 - b*total + 3*delta/2
    assert expand(sum_equation - (2*total*quartic_condition-gap_residual)) == 0


if __name__ == "__main__":
    check_parameter_equations()
    check_collision_elimination()
    print("Verified the octic approximate-root, parameter and collision identities.")
