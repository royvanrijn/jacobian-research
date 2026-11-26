"""Check the coefficient, collision and leading-pole identities in the cubic–quintic proof."""

from sympy import Poly, Rational, diff, expand, resultant, symbols


def check_coefficient_identities():
    z, a, b, d, e = symbols("z a b d e")
    a_y, b_y, d_y, e_y = symbols("a_y b_y d_y e_y")
    c, j, alpha, beta, gamma, delta = symbols("c j alpha beta gamma delta")
    cubic = z**3 + a*z
    quintic = z**5 + b*z**4 + d*z**2 + e*z
    jacobian = diff(cubic, z) * (b_y*z**4 + d_y*z**2 + e_y*z)
    jacobian -= a_y * z * diff(quintic, z)
    primitive = alpha*cubic + beta*quintic + gamma*cubic**2 + delta*cubic*quintic
    residual = Poly(jacobian - j - c*diff(cubic, z)*quintic - diff(primitive, z), z)
    scale = {delta: -3*c/8, b_y: b*c/8}
    assert residual.nth(7).subs(scale) == 0
    assert expand(residual.nth(6).subs(scale)) == 0

    normalized = {
        c: 0,
        delta: 0,
        b_y: 0,
        gamma: -5*a_y/6,
        d_y: (5*beta + 4*b*a_y)/3,
        e_y: 4*b*beta/3 - 20*a*a_y/9,
        alpha: 5*a*beta/9 + 4*a*a_y*b/9 - 2*a_y*d/3,
    }
    for degree in range(2, 8):
        assert expand(residual.nth(degree).subs(normalized)) == 0

    m, m_y, constant = symbols("m m_y C")
    move_origin = {d: m + 4*b*a/3, e: constant - 10*a**2/9 + 4*b*m/5}
    constant_y = e_y + 20*a*a_y/9 - Rational(4, 5)*b*(d_y-4*b*a_y/3)
    assert expand(constant_y.subs(normalized)) == 0

    first_integral = 5*a**3/3 - 9*constant*a - 36*b*a*m/5 - 27*m**2/5
    first_derivative = diff(first_integral, a)*a_y + diff(first_integral, m)*m_y
    coefficient_one = residual.nth(1).subs(normalized).subs(move_origin)
    assert expand(first_derivative.subs(m_y, 5*beta/3) - 9*coefficient_one) == 0

    primitive_time = ((a**2/3 - 3*constant/5)*m - 6*b*m**2/25 + 4*b*a**3/27)
    primitive_derivative = diff(primitive_time, a)*a_y + diff(primitive_time, m)*m_y
    coefficient_zero = residual.nth(0).subs(normalized).subs(move_origin)
    assert expand(primitive_derivative.subs(m_y, 5*beta/3) - j - coefficient_zero) == 0


def check_collision_equations():
    total, product, a, b, d, e = symbols("s p a b d e")
    differences = [0, 1]
    for _ in range(2, 6):
        differences.append(expand(total*differences[-1] - product*differences[-2]))
    assert differences[3] + a == total**2 - product + a
    quintic = differences[5] + b*differences[4] + d*differences[2] + e
    elimination = total**4 + b*total**3 + a*total**2 + (2*a*b-d)*total - (a**2+e)
    assert expand(quintic.subs(product, total**2+a) + elimination) == 0
    assert expand((total**2 - 4*product).subs(product, total**2+a)) == -3*total**2 - 4*a

    lead_a, lead_m = symbols("A M")
    leading_sum = total**4 + lead_a*total**2 - lead_m*total + lead_a**2/9
    collision = resultant(leading_sum, 3*total**2 + 4*lead_a, total)
    collision = expand(collision).subs(lead_m**2, 25*lead_a**3/81)
    assert expand(collision - Rational(175, 3)*lead_a**4) == 0


if __name__ == "__main__":
    check_coefficient_identities()
    check_collision_equations()
    print("Verified the cubic–quintic coefficient, first-integral and gap identities.")
