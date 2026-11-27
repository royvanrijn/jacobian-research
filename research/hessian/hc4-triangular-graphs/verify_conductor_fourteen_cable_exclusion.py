"""Check exact generator, parameter, and collision identities in the cable proof."""

from sympy import Poly, Rational, diff, expand, solve, symbols


def check_cable_identities():
    z, a, b, c, p, r, t, R, T = symbols("z a b c p r t R T")
    a_y, b_y, c_y, p_y, r_y, t_y = symbols("a_y b_y c_y p_y r_y t_y")
    h_u, h_v, h_u2 = symbols("h_u h_v h_u2")
    quartic = z**4 + a*z**2 + b*z + c
    quadratic = z**2 + a/2
    delta = a**2 - 4*c
    sextic = (3*quartic*quadratic/2 - quadratic**3/2 + 3*b**2/8
              + p*z**5 + r*z**3 + t*z)

    def dy(f):
        return (
            diff(f, a)*a_y + diff(f, b)*b_y + diff(f, c)*c_y
            + diff(f, p)*p_y + diff(f, r)*r_y + diff(f, t)*t_y
        )

    hamiltonian = h_u*quartic + h_v*sextic + h_u2*quartic**2
    comparison = Poly(
        diff(quartic, z)*dy(sextic) - dy(quartic)*diff(sextic, z)
        - diff(hamiltonian, z), z
    )
    r_value = R + 5*p*a/4
    t_value = T + 5*p*a**2/32 + 3*R*a/4 + 5*p*c/4 - 25*p**2*b/24
    integrations = {r: r_value, t: t_value}
    expected = {
        p_y: 0,
        h_u2: 0,
        r_y: dy(r_value).subs(p_y, 0),
        h_v: -5*p*b_y/6,
        t_y: dy(t_value).subs(p_y, 0),
        h_u: 5*p*dy(a*b)/16 - 3*R*b_y/4,
    }
    substitutions = {}
    for degree, variable in zip(range(8, 2, -1), expected, strict=True):
        coefficient = comparison.nth(degree).subs(integrations).subs(substitutions)
        solution = solve(coefficient, variable)
        assert len(solution) == 1
        assert expand(solution[0] - expected[variable]) == 0
        substitutions[variable] = expected[variable]

    first = (3*delta*R/4 + 5*delta*a*p/16 - 3*R*a**2/8 - T*a
             - 5*a**3*p/32 + 15*b**2*p/8 + b*(3*delta/4 + 5*R*p/2 + 25*a*p**2/24))
    second = -3*delta**2/32 + 3*a*b**2/4 + b*(-5*delta*p/16 + 3*R*a/4 - T + 5*a**2*p/32)
    primitive = (-5*delta**2*p/128 + delta*T/4 + 5*delta*a**2*p/128
                 - R*a**3/16 - T*a**2/8 - 15*a**4*p/512 + b**3/4
                 + b**2*(3*R/8 - 5*a*p/16 - 125*p**3/288)
                 + b*(-3*delta*a/16 - 25*delta*p**2/96 + 5*T*p/6 + 25*a**2*p**2/192))
    residual = Poly(comparison.as_expr().subs(integrations).subs(substitutions), z)
    assert residual.degree() == 2
    assert expand(residual.nth(2) - dy(first).subs(p_y, 0)) == 0
    assert expand(residual.nth(1) - dy(second).subs(p_y, 0)) == 0
    assert expand(residual.nth(0) - dy(primitive + a*first/4).subs(p_y, 0) + a_y*first/4) == 0

    total, gap = symbols("s g")
    upper, lower = (total+gap)/2, (total-gap)/2
    quartic_diff = expand((quartic.subs(z, upper)-quartic.subs(z, lower))/gap)
    normalized_sextic = sextic.subs(integrations)
    sextic_diff = expand((normalized_sextic.subs(z, upper)-normalized_sextic.subs(z, lower))/gap)
    squared_gap = -total**2-2*a-2*b/total
    assert expand(quartic_diff.subs(gap**2, squared_gap)) == 0
    sum_equation = (2*total**7 + 2*p*total**6 + 4*a*total**5 + (-4*R+3*a*p-2*b)*total**4
                    + (3*delta+8*b*p)*total**3
                    + (5*delta*p/2-2*R*a-8*T-3*a**2*p/4+25*b*p**2/3)*total**2
                    + (4*R*b+a*b*p)*total - 2*p*b**2)
    reduced = sextic_diff.subs(gap**4, squared_gap**2).subs(gap**2, squared_gap)
    assert expand(8*total**2*reduced + sum_equation) == 0

    # The smallest-sum face in the a-leading pole regime must not have s=-b/a.
    lead_delta = 5*p*a**3/(24*b)
    small_face = 3*lead_delta*total**3 - 3*a**2*p*total**2/4 + a*b*p*total - 2*p*b**2
    assert expand(small_face.subs(total, -b/a) + 35*p*b**2/8) == 0


if __name__ == "__main__":
    check_cable_identities()
    print("Verified the conductor-fourteen cable parameter equations and collision identities.")
