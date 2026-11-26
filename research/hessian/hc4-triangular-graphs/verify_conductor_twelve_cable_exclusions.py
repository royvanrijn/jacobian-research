"""Replay exact algebra used by the conductor-twelve cable exclusions."""

from sympy import Poly, Rational, diff, discriminant, expand, integrate, rem, symbols


def sextic_parameter_and_gaps():
    z, a, b, c, p, R, a_y, b_y, c_y = symbols("z a b c p R a_y b_y c_y")
    coefficient, p_y, coefficient_y = symbols("r p_y r_y")
    delta = a**2 - 4*c
    quartic = z**4 + a*z**2 + b*z + c
    quadratic = z**2 + a/2
    sextic = (3*quartic*quadratic/2 - quadratic**3/2 + 3*b**2/8
              + p*z**3 + (R + 3*p*a/4)*z)
    first = (b+p)*delta - p*a**2/2 - 4*R*a/3
    second = delta**2 - 8*a*b*(b+p) + 32*R*b/3
    primitive = (b**3/4 + 3*p*b**2/8 - 3*first*a/16 - 5*p*a**3/32
                 + 3*p*a*delta/16 - 3*R*a**2/8 + R*delta/4)

    def dy(f):
        return diff(f, a)*a_y + diff(f, b)*b_y + diff(f, c)*c_y

    uncorrected = sextic + (coefficient - R - 3*p*a/4)*z
    uncorrected_y = diff(uncorrected, p)*p_y + diff(uncorrected, coefficient)*coefficient_y
    uncorrected_y += dy(uncorrected)
    initial = Poly(diff(quartic, z)*uncorrected_y - dy(quartic)*diff(uncorrected, z), z)
    assert initial.nth(6) == 4*p_y
    assert expand(initial.nth(4).subs(p_y, 0) - 4*coefficient_y + 3*p*a_y) == 0

    residual = Poly(diff(quartic, z)*dy(sextic) - dy(quartic)*diff(sextic, z)
                    + 3*p*b_y*diff(quartic, z)/4, z)
    assert residual.degree() == 2
    assert expand(residual.nth(2) - 3*dy(first)/4) == 0
    assert expand(residual.nth(1) + 3*dy(second)/32) == 0
    assert expand(residual.nth(0) - dy(primitive) - 3*a*dy(first)/16) == 0

    total, gap = symbols("s g")
    upper, lower = (total+gap)/2, (total-gap)/2
    quartic_diff = expand((quartic.subs(z, upper)-quartic.subs(z, lower))/gap)
    sextic_diff = expand((sextic.subs(z, upper)-sextic.subs(z, lower))/gap)
    squared_gap = -total**2-2*a-2*b/total
    assert expand(quartic_diff.subs(gap**2, squared_gap)) == 0
    sum_polynomial = (2*total**6+4*a*total**4-(2*b+4*p)*total**3
                      + 3*delta*total**2-(2*p*a+8*R)*total+4*p*b)
    reduced = sextic_diff.subs(gap**4, squared_gap**2).subs(gap**2, squared_gap)
    assert expand(8*total*reduced + sum_polynomial) == 0


def degree_nine_model():
    z, k, r, T = symbols("z k r T")
    cubic = z**3+k*z
    quadratic = r*(z**2+k/2)
    sextic = expand(cubic**2+quadratic)
    ninth = expand(cubic**3+3*cubic*quadratic/2+3*r**2*z/8)
    alpha = r**2*(2*r-3*k**2)/16
    root = expand(ninth**2-sextic**3-alpha*sextic)
    assert expand(
        root+r**3*(2*k**3+8*k**2*z**2+4*k*r+4*k*z**4-r*z**2)/64
    ) == 0
    assert Poly(root, z).nth(4) == -k*r**3/16

    special = {k: 1, r: Rational(3, 2)}
    sextic = sextic.subs(special)
    ninth = ninth.subs(special)
    root = root.subs(special)
    assert expand(diff(root, z)+27*z*(16*z**2+13)/512) == 0
    assert diff(ninth, z).subs(z, 0) == Rational(63, 32)
    assert expand(diff(sextic, z)/(2*z)).subs(z**2, -Rational(13, 16)) == Rational(315, 256)

    g = expand(ninth/z).subs(z**2, T)
    assert discriminant(g, T) == Rational(35721, 4096)

    def action_mod_g(P):
        primitive = integrate(-diff(P, z)*ninth, z)
        action = expand((primitive-primitive.subs(z, -z))/z).subs(z**2, T)
        return rem(action, g, T)

    assert action_mod_g(sextic) == Rational(189, 160)
    expected = (-Rational(1161, 36608)*T**3-Rational(8883, 56320)*T**2
                -Rational(4473, 18304)*T-Rational(200151, 585728))
    assert expand(action_mod_g(root)-expected) == 0


if __name__ == "__main__":
    sextic_parameter_and_gaps()
    degree_nine_model()
    print("Verified cable parameter coefficients, collision equations and reflected actions.")
