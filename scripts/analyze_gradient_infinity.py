#!/usr/bin/env python3
"""Exact first chart calculation for the gradient flow at infinity."""
import sympy as sp

x, y, z = sp.symbols("x y z")
r, Y, Z = sp.symbols("r Y Z")
U, gamma = sp.symbols("U gamma")

u = 1 + x*y
F = sp.Matrix([
    u**3*z + y**2*u*(4 + 3*x*y),
    y + 3*x*u**2*z + 3*x*y**2*(4 + 3*x*y),
    2*x - 3*x**2*y - x**3*z,
])
q = sp.Matrix([-sp.Rational(1, 4), 0, 0])
L = 16*(F - q).dot(F - q)
negative_gradient = -sp.Matrix([sp.diff(L, variable) for variable in (x, y, z)])

# Weighted infinity chart: x=1/r, Y=xy, Z=x^2*z. Its natural construction
# coordinates are U=1+Y and gamma=1-3Y/2-Z/2.
chart = {x: 1/r, y: Y*r, z: Z*r**2}
vx, vy, vz = [sp.factor(component.subs(chart)) for component in negative_gradient]
chart_field = sp.Matrix([
    -r**2*vx,
    Y*r*vx + vy/r,
    2*Z*r*vx + vz/r**2,
])

natural_coordinates = {Y: U - 1, Z: 5 - 3*U - 2*gamma}
natural_field = sp.Matrix([
    chart_field[0],
    chart_field[1],
    -sp.Rational(3, 2)*chart_field[1] - sp.Rational(1, 2)*chart_field[2],
]).subs(natural_coordinates).applyfunc(sp.cancel)

# r^6 clears all chart poles. This changes only the time parametrization in
# the interior r != 0.
desingularized = (r**6*natural_field).applyfunc(sp.cancel)
boundary = sp.simplify(desingularized.subs(r, 0))
assert boundary == sp.Matrix([0, 0, -32*gamma])

jacobian_on_equilibrium_line = sp.simplify(
    desingularized.jacobian((r, U, gamma)).subs({r: 0, gamma: 0})
)
assert jacobian_on_equilibrium_line == sp.diag(0, 0, -32)

# The Palais--Smale curve is U=gamma=0, but it is not an orbit of the
# negative-gradient flow.
field_on_ps_curve = natural_field.subs({U: 0, gamma: 0}).applyfunc(sp.factor)
expected_ps_field = sp.Matrix([
    -392*r**5,
    -8*(49*r**4 + 33),
    -4*(343*r**4 - 99),
])
assert (field_on_ps_curve - expected_ps_field).applyfunc(sp.simplify) == sp.zeros(3, 1)


def terms_through_total_degree(expression, variables, degree):
    """Return the terms whose total degree is at most ``degree``."""
    polynomial = sp.Poly(sp.expand(expression), *variables)
    return sp.Add(*[
        coefficient*sp.prod(variable**power for variable, power in zip(variables, powers))
        for powers, coefficient in polynomial.terms()
        if sum(powers) <= degree
    ])


# First formal center-manifold jet gamma=h(r,U). The invariance defect has no
# terms through total degree seven. Higher terms of h are not fixed here.
h = (
    3*U**2*r**4
    + sp.Rational(49, 8)*U**3*r**4
    + sp.Rational(99, 8)*r**6
    + sp.Rational(99, 4)*U*r**6
)
reduced_r = sp.expand(desingularized[0].subs(gamma, h))
reduced_U = sp.expand(desingularized[1].subs(gamma, h))
invariance_defect = sp.expand(
    desingularized[2].subs(gamma, h)
    - sp.diff(h, r)*reduced_r
    - sp.diff(h, U)*reduced_U
)
assert terms_through_total_degree(invariance_defect, (r, U), 7) == 0
assert terms_through_total_degree(reduced_r, (r, U), 11) == -392*r**11
assert terms_through_total_degree(reduced_U, (r, U), 7) == -264*r**6 - 528*U*r**6

print("PASS: r^6 desingularizes the weighted infinity chart")
print("PASS: boundary field is", tuple(boundary))
print("PASS: boundary equilibria are the line r=gamma=0")
print("PASS: transverse eigenvalues along that line are (0, 0, -32)")
print("PASS: the explicit Palais--Smale curve is not a gradient trajectory")
print("PASS: formal center-manifold and reduced-flow jets verified")
