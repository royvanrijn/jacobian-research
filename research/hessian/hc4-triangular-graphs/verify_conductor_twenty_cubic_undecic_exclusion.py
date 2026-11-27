"""Replay the cubic–undecic parameter equations and leading collision tests."""

from sympy import Poly, Rational, diff, expand, integrate, rem, resultant, symbols


z, a, b, c, d, e, f, g, h = symbols("z a b c d e f g h")
B, D, E, F, G, H = symbols("B D E F G H")
da, db, dc, dd, de, df, dg, dh = symbols("da db dc dd de df dg dh")
coefficients = dict(zip((a, b, c, d, e, f, g, h), (da, db, dc, dd, de, df, dg, dh)))
lower = z**3 + a*z
upper = z**11 + b*z**10 + c*z**8 + d*z**7 + e*z**5 + f*z**4 + g*z**2 + h*z


def derivative(expression):
    return sum(diff(expression, coefficient)*rate
               for coefficient, rate in coefficients.items())


def primitive(form, variables):
    exact = 0
    for variable in variables:
        residual = expand(form - derivative(exact))
        exact += integrate(residual.coeff(coefficients[variable]), variable)
    assert expand(form - derivative(exact)) == 0
    return expand(exact)


def parameter_equations():
    jacobian = Poly(diff(lower, z)*derivative(upper) - derivative(lower)*diff(upper, z), z)
    assert jacobian.degree() == 12 and jacobian.nth(12) == 3*db
    replacements = {b: B, db: 0}
    residual = Poly(jacobian.as_expr().subs(replacements, simultaneous=True), z)
    for potential in (lower**4, upper.subs(b, B), lower**3, lower**2, lower):
        term = Poly(diff(potential, z), z)
        residual -= term*residual.nth(term.degree())/term.LC()

    variables = [a, c, d, e, f, g, h]
    weights = {a: 2, c: 3, B: 1, D: 4, E: 6, F: 7, G: 9, H: 10}
    scale = symbols("scale")
    scaling = {variable: variable*scale**weight for variable, weight in weights.items()}
    for degree, variable, constant in ((9, d, D), (7, e, E), (6, f, F), (4, g, G), (3, h, H)):
        assert residual.degree() == degree
        rate = coefficients[variable]
        form = expand(residual.nth(degree))
        assert form.coeff(rate) == 3
        variables.remove(variable)
        integral = constant + primitive(rate - form/3, variables)
        assert expand(integral.subs(scaling, simultaneous=True)
                      - integral*scale**weights[constant]) == 0
        replacement = {variable: integral, rate: derivative(integral)}
        residual = Poly(residual.as_expr().subs(replacement, simultaneous=True), z)
        replacements.update(replacement)

    assert residual.degree() == 1
    first_integral = primitive(residual.nth(1), variables)
    scale_primitive = primitive(residual.nth(0), variables)
    for integral, weight in ((first_integral, 12), (scale_primitive, 13)):
        assert expand(integral.subs(scaling, simultaneous=True) - integral*scale**weight) == 0
    return upper.subs(replacements, simultaneous=True), first_integral, scale_primitive


def leading_collision_tests(upper_polynomial, first_integral, scale_primitive):
    S, p, A, C, t = symbols("S p A C t")
    monomials = [0, 1]
    for degree in range(2, 12):
        monomials.append(expand(S*monomials[-1] - p*monomials[-2]))

    def divided_difference(polynomial):
        coefficients_in_z = Poly(polynomial, z)
        return sum(coefficients_in_z.nth(i)*monomials[i] for i in range(12))

    assert expand(divided_difference(lower) - (S**2+a-p)) == 0
    collision = -divided_difference(upper_polynomial).subs(p, S**2+a)
    collision = Poly(expand(collision), S)
    assert collision.degree() == 10 and collision.LC() == 1
    scaling = {a: A/t**2, c: C/t**3, S: S/t}

    def initial(expression, weight):
        return Poly(expand(t**weight*expression.subs(scaling, simultaneous=True)), t).nth(0)

    relation = 29282*A**6 - 98010*A**3*C**2 + 10935*C**4
    assert expand(initial(first_integral, 12) + 2*relation/2910897) == 0
    assert expand(initial(scale_primitive, 13) + 4*A**2*C*(242*A**3-135*C**2)/29403) == 0
    polynomial = initial(collision.as_expr(), 10).subs(A, 1)
    gap_resultant = resultant(polynomial, 3*S**2+4, S)
    assert rem(gap_resultant, relation.subs(A, 1), C) == 4576*(45*C**2-11)/9
    u = symbols("u")
    relation_in_u = relation.subs(A, 1).subs(C**2, u)
    assert relation_in_u.subs(u, Rational(242, 135)) == Rational(-556358, 5)
    assert relation_in_u.subs(u, Rational(11, 45)) == Rational(29887, 5)


if __name__ == "__main__":
    leading_collision_tests(*parameter_equations())
    print("Verified the cubic–undecic parameter integrations and leading gap obstruction.")
