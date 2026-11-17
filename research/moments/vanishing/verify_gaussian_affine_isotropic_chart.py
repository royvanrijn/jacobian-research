"""Exact identities, a fixed-family certificate and bounded Wick controls.

The written proof supplies completeness and unrestricted degree coverage.
This checker supplies symbolic identity checks, explicitly bounded controls,
and a complete coefficient-elimination certificate for the stated resonant
six-parameter quartic family. That certificate is distinct from a vanishing
moment prefix: it forces a form whose all-order safety is proved by weight.
Run from the repository root with .venv/bin/python and --output PATH to retain
a receipt. No receipt is written unless that option is supplied.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from fractions import Fraction
from pathlib import Path

import sympy as sp


Z, W, Y, T = sp.symbols("z w y t")


def zero(expr):
    assert sp.cancel(sp.expand(expr)) == 0, expr


def expectation(expr):
    result = 0
    for (p, q, r), coefficient in sp.Poly(sp.expand(expr), Z, W, Y).terms():
        if p == q and r % 2 == 0:
            result += coefficient * math.factorial(p) * math.factorial(r) / (
                2 ** (r // 2) * math.factorial(r // 2)
            )
    return sp.expand(result)


def sparse(expr):
    return {powers: Fraction(coefficient) for powers, coefficient in
            sp.Poly(expr, Z, W, Y).terms()}


def multiply(left, right):
    result = {}
    for alpha, a in left.items():
        for beta, b in right.items():
            gamma = tuple(x + y for x, y in zip(alpha, beta))
            result[gamma] = result.get(gamma, Fraction(0)) + a * b
    return {gamma: value for gamma, value in result.items() if value}


def wick(poly, marker=(0, 0, 0)):
    result = Fraction(0)
    for powers, coefficient in poly.items():
        p, q, r = (a + b for a, b in zip(powers, marker))
        if p == q and r % 2 == 0:
            result += coefficient * math.factorial(p) * math.factorial(r) / (
                2 ** (r // 2) * math.factorial(r // 2)
            )
    return result


def radial_controls():
    """Compare the composed divergent-series formula with direct factorial Wick."""
    radial, y, t = sp.symbols("R y t")
    cases = [(1, 2, 3, -1, 2, -4), (0, 1, 0, 0, 0, 0),
             (1, 0, 0, 0, 0, -2), (0, 0, 2, -2, 1, 0)]
    order = 8
    for params in cases:
        alpha, beta, a, gamma, b, k = map(sp.Rational, params)
        phase = alpha * radial**2 + beta * radial * y + a * y**2 + gamma * radial + b * y + k
        d = 1 - 2 * a * t
        quadratic = alpha * t + beta**2 * t**2 / (2 * d)
        linear = gamma * t + beta * b * t**2 / d
        scalar = k * t + b**2 * t**2 / (2 * d)
        argument = sp.series(quadratic / (1 - linear)**2, t, 0, order + 1).removeO()
        composition = sp.Integer(1)
        argument_power = sp.Integer(1)
        for n in range(1, order + 1):
            argument_power = sp.series(argument_power * argument, t, 0, order + 1).removeO()
            composition += sp.Rational(math.factorial(2*n), math.factorial(n)) * argument_power
        prefactor = sp.series(sp.exp(scalar) / (sp.sqrt(d) * (1-linear)), t, 0, order+1).removeO()
        formula = sp.series(prefactor * composition, t, 0, order+1).removeO().expand()
        power = sp.Integer(1)
        for m in range(order+1):
            contraction = sp.Integer(0)
            for (i, j), coefficient in sp.Poly(power, radial, y).terms():
                if j % 2 == 0:
                    contraction += coefficient * math.factorial(i) * sp.Rational(
                        math.factorial(j), 2**(j//2) * math.factorial(j//2))
            zero(formula.coeff(t, m) * math.factorial(m) - contraction)
            power = sp.expand(power * phase)
    return {"parameter_cases": cases, "moment_orders": [0, order]}


def relative_endpoint_controls():
    """Bounded checks of the exact relative relations, including the endpoint."""
    radial, t = sp.symbols("R t")
    d = 1-2*t
    hnorm = 1/sp.sqrt(d)
    potential = t*(radial**3+2*radial) + t**2*(radial**2+1)**2/(2*d)
    phi = -radial+potential
    order = 4
    exponential = sp.series(sp.exp(potential), t, 0, order+1).removeO()
    normalization = sp.series(hnorm, t, 0, order+1).removeO()
    def functional(q):
        expr = sp.series(q*exponential, t, 0, order+1).removeO().expand()
        contracted = sum(coefficient*math.factorial(power[0])
                         for power, coefficient in sp.Poly(expr, radial).terms())
        return sp.series(normalization*contracted, t, 0, order+1).removeO().expand()
    for j in range(1, 4):
        f = radial**j
        zero(functional(sp.diff(f, radial)+sp.diff(phi, radial)*f))
    endpoint_period = -sp.series(hnorm*sp.exp(potential.subs(radial, 0)), t, 0, order+1).removeO()
    zero(functional(sp.diff(phi, radial))-endpoint_period)
    # Universal boundary-class connection, not just the selected test phase.
    abstract_phi = sp.Function("phi")(t, radial)
    px, pt = sp.diff(abstract_phi, radial), sp.diff(abstract_phi, t)
    endpoint_pt = pt.subs(radial, 0)
    boundary_action = sp.diff(px, t)+(pt-endpoint_pt)*px
    zero(boundary_action-(sp.diff(pt-endpoint_pt, radial)+px*(pt-endpoint_pt)))
    return {"g": "R^3+2R", "B": "R^2+1", "a": 1,
            "parameter_orders": [0, order], "relative_relations": ["D(R)", "D(R^2)", "D(R^3)"],
            "endpoint_relation": "L(D(1))=-H exp(phi(0))"}


def resonant_quartic_certificate():
    """Verify complete elimination for a fixed mixed quartic coefficient family.

    P = alpha*Z^3 + beta*W^3 + gamma*Z^2*W*Y + delta*Z*W^2*Y
        + epsilon*((ZW)^2-2) + eta*(Y^2-1).

    Universal Wick identities force epsilon=eta=0 and all four products
    between (alpha,gamma) and (beta,delta) to vanish. Over C this is exactly
    strict one-sided support. The written weight argument then gives all
    marked-period cutoffs; no finite prefix is treated as sufficient alone.
    """
    alpha, beta, gamma, delta, eps, eta = sp.symbols(
        "alpha beta gamma delta epsilon eta")
    inv_a, inv_d = sp.symbols("A D")
    phase = (alpha*Z**3 + beta*W**3 + gamma*Z**2*W*Y + delta*Z*W**2*Y
             + eps*(Z**2*W**2-2) + eta*(Y**2-1))
    moments = {n: expectation(phase**n) for n in (2, 3, 4, 5)}
    zero(moments[2]-2*(6*alpha*beta+6*gamma*delta+10*eps**2+eta**2))
    zero(moments[3]-8*(81*eps*(alpha*beta+gamma*delta)
                       + 9*eta*gamma*delta+74*eps**3+eta**3))

    # This establishes invariant dependence before using a convenient section.
    # Consequently alpha=gamma=1 below loses no coefficient-zero faces.
    for n, value in moments.items():
        for powers, coefficient in sp.Poly(
                value, alpha, beta, gamma, delta, eps, eta).terms():
            i, j, k, ell, _, _ = powers
            assert i == j and k == ell, (n, powers, coefficient)

    invariant_moments = {
        n: sp.expand(value.subs({alpha: 1, beta: inv_a, gamma: 1, delta: inv_d}))
        for n, value in moments.items()
    }
    f_coefficients = [24400, 2364, 6864, -652, 468, -63]
    g_coefficients = [4196800, 537276, 1192320, -13760, 63780, -1215, -722]
    f = sp.Poly.from_list(f_coefficients, gens=eps)
    g = sp.Poly.from_list(g_coefficients, gens=eps)

    # If eta is nonzero, scale P by 1/eta. Moments two and three then give
    # these substitutions, and moments four and five give coprime eliminants.
    a_eta_one = -inv_d-sp.Rational(5, 3)*eps**2-sp.Rational(1, 6)
    d_eta_one = (122*eps**3+27*eps-2)/18
    zero(invariant_moments[2].subs({eta: 1, inv_a: a_eta_one}))
    zero(invariant_moments[3].subs(eta, 1).subs(inv_a, a_eta_one)
         .subs(inv_d, d_eta_one))
    zero(invariant_moments[4].subs(eta, 1).subs(inv_a, a_eta_one)
         .subs(inv_d, d_eta_one)+4*f.as_expr())
    zero(invariant_moments[5].subs(eta, 1).subs(inv_a, a_eta_one)
         .subs(inv_d, d_eta_one)+sp.Rational(16, 3)*g.as_expr())

    # Construct the Sylvester matrix explicitly and evaluate its determinant
    # over F_7 by independent modular row elimination. A nonzero reduction
    # certifies a nonzero integer resultant, hence no common complex root.
    prime = 7
    assert sp.Poly(f.as_expr(), eps, modulus=prime).degree() == 5
    assert sp.Poly(g.as_expr(), eps, modulus=prime).degree() == 6
    sylvester = []
    for shift in range(6):
        sylvester.append([0]*shift+f_coefficients+[0]*(5-shift))
    for shift in range(5):
        sylvester.append([0]*shift+g_coefficients+[0]*(4-shift))
    assert len(sylvester) == 11 and all(len(row) == 11 for row in sylvester)
    matrix = [[coefficient % prime for coefficient in row] for row in sylvester]
    determinant = 1
    for column in range(11):
        pivot = next((row for row in range(column, 11)
                      if matrix[row][column] % prime), None)
        assert pivot is not None, (column, matrix)
        if pivot != column:
            matrix[pivot], matrix[column] = matrix[column], matrix[pivot]
            determinant = -determinant
        pivot_value = matrix[column][column] % prime
        determinant = determinant*pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        for row in range(column+1, 11):
            factor = matrix[row][column]*inverse % prime
            matrix[row] = [(left-factor*right) % prime
                           for left, right in zip(matrix[row], matrix[column])]
    assert determinant % prime == 5
    assert int(sp.resultant(f.as_expr(), g.as_expr(), eps)) % prime == 5

    # Eta=0: the first two remaining relations force epsilon=0.
    a_eta_zero = -inv_d-sp.Rational(5, 3)*eps**2
    zero(invariant_moments[2].subs(eta, 0).subs(inv_a, a_eta_zero))
    zero(invariant_moments[3].subs(eta, 0).subs(inv_a, a_eta_zero)
         + 488*eps**3)

    core = phase.subs({eps: 0, eta: 0})
    moment_six = expectation(core**6)
    expected_six = math.factorial(9)*(
        20*(alpha*beta)**3 + 180*(alpha*beta)**2*(gamma*delta)
        + 540*(alpha*beta)*(gamma*delta)**2 + 300*(gamma*delta)**3)
    zero(moment_six-expected_six)
    zero(expected_six.subs({alpha: 1, beta: inv_a, gamma: 1, delta: inv_d})
         .subs(inv_d, -inv_a)-80*math.factorial(9)*inv_a**3)

    # Reduce the universal eighth moment modulo A=alpha*beta=D=gamma*delta=0.
    reduced_eight = sp.Integer(0)
    for powers, coefficient in sp.Poly(
            expectation(core**8), alpha, beta, gamma, delta).terms():
        i, j, k, ell = powers
        if (i and j) or (k and ell):
            continue
        reduced_eight += coefficient*alpha**i*beta**j*gamma**k*delta**ell
    left_cross, right_cross = alpha**2*delta**6, beta**2*gamma**6
    zero(reduced_eight-420*math.factorial(12)*(left_cross+right_cross))
    zero(left_cross*right_cross-(alpha*beta)**2*(gamma*delta)**6)

    # All four cross products vanish exactly when one angular side is empty.
    # The eight-moment identity gives both residual cross products zero over C:
    # the two summands sum to zero and their product is already zero.
    for mask in range(16):
        ap, bp, gp, dp = (bool(mask & (1 << bit)) for bit in range(4))
        cross_products_zero = not ((ap and bp) or (gp and dp)
                                   or (ap and dp) or (bp and gp))
        one_side_empty = (not ap and not gp) or (not bp and not dp)
        assert cross_products_zero == one_side_empty

    return {
        "kind": "complete coefficient-elimination certificate for this fixed family",
        "coefficient_domain": "C; universal polynomial computations over Q",
        "phase": "alpha Z^3+beta W^3+gamma Z^2 W Y+delta Z W^2 Y+epsilon((ZW)^2-2)+eta(Y^2-1)",
        "moment_orders_used": [2, 3, 4, 5, 6, 8],
        "invariant_dependence_orders": [2, 3, 4, 5],
        "invariants": {"A": "alpha beta", "D": "gamma delta"},
        "eliminants": {"variable": "epsilon",
                       "F_coefficients_descending": f_coefficients,
                       "G_coefficients_descending": g_coefficients},
        "modular_resultant": {"prime": prime, "degrees_preserved": [5, 6],
                              "sylvester_dimension": 11, "determinant_mod_prime": 5,
                              "independent_sympy_resultant_mod_prime": 5},
        "necessity": [
            "eta!=0: scale to eta=1; moments2,3 eliminate A,D; coprime F,G exclude this case",
            "eta=0: moments2,3 give -488 epsilon^3=0, hence epsilon=0",
            "moment6 with D=-A gives 80*9!*A^3=0, hence A=D=0",
            "moment8 modulo A=D=0 gives 420*12!*(alpha^2 delta^6+beta^2 gamma^6)=0; summands have product0",
            "all cross products between (alpha,gamma) and (beta,delta) vanish",
        ],
        "pure_iff": "epsilon=eta=0 and either alpha=gamma=0 or beta=delta=0",
        "all_order_sufficiency": "strict one-sided angular support; written weight proof gives m>deg Q for every marker",
        "finite_prefix_scope": "The verified finite equations force the explicit all-order safe form; no general finite-prefix sufficiency is asserted.",
    }


def common_root_cubic_barrier():
    """Exact zero-order local chart; the original phase is explicitly nonpure."""
    r = sp.symbols("barrier_r", nonzero=True)
    generic_c, generic_h = sp.symbols("barrier_c barrier_h", nonzero=True)
    generic_u = -generic_h**2/(6*Z*generic_c)
    generic_k = -generic_h**3/(54*Z*generic_c**2)
    source = (r*Z-generic_h)/generic_c
    phase = -source**2/2+(generic_u*source+generic_k)/r
    cube = -(r*Z-sp.Rational(2, 3)*generic_h)**3/(2*r*Z*generic_c**2)
    zero(phase-cube)

    parameter_h = sp.symbols("barrier_H", nonzero=True)
    c = (Z-1)**9
    h = parameter_h*Z**2*(Z-1)**6
    u = -parameter_h**2*Z**3*(Z-1)**3/6
    k = -parameter_h**3*Z**5/54
    zero(u+ h**2/(6*Z*c))
    zero(k+ h**3/(54*Z*c**2))
    polynomial_phase = W*(h+c*Y)+u*Y+k
    zero(expectation(polynomial_phase))
    zero(expectation(polynomial_phase**2)-306)

    # Critical curve and derivative are exact rational identities.
    zeta = sp.symbols("zeta")
    factor = r*Z-sp.Rational(2, 3)*h
    critical_r = sp.Rational(2, 3)*parameter_h*zeta*(zeta-1)**6
    zero(factor.subs({Z: zeta, r: critical_r}))
    critical_fz = sp.diff(factor, Z).subs({Z: zeta, r: critical_r})
    expected_fz = -sp.Rational(2, 3)*parameter_h*zeta*(zeta-1)**5*(7*zeta-1)
    zero(critical_fz-expected_fz)

    # A unit reparametrization avoids an approximate Puiseux centre:
    # H=3/(2A^6), zeta=1+A R, and r=zeta R^6 satisfy the critical curve exactly.
    chart_a = sp.symbols("barrier_A", nonzero=True)
    chart_r, chart_x = sp.symbols("barrier_R barrier_X")
    exact_h = sp.Rational(3, 2)/chart_a**6
    exact_zeta = 1+chart_a*chart_r
    exact_r = exact_zeta*chart_r**6
    chart_z = exact_zeta+chart_r**3*chart_x
    zero(exact_r-critical_r.subs({zeta: exact_zeta, parameter_h: exact_h}))
    exact_critical_fz = critical_fz.subs({zeta: exact_zeta, parameter_h: exact_h})
    zero(exact_critical_fz+exact_zeta*chart_r**5*(7*exact_zeta-1)/chart_a)
    leading_fz = sp.cancel(exact_critical_fz/chart_r**5).subs(chart_r, 0)
    zero(leading_fz+6/chart_a)

    chart_c = c.subs(Z, chart_z)
    zero(chart_c-chart_r**9*(chart_a+chart_r**2*chart_x)**9)
    chart_factor = factor.subs({r: exact_r, parameter_h: exact_h, Z: chart_z})
    scaled_factor = sp.cancel(chart_factor/chart_r**8)
    assert sp.Poly(scaled_factor, chart_r, chart_x).degree(chart_r) >= 0
    zero(scaled_factor.subs(chart_r, 0)+6*chart_x/chart_a)
    # Cancel the exact R^24 factors before taking the leading phase.
    chart_phase = -scaled_factor**3/(
        2*exact_zeta*chart_z*(chart_a+chart_r**2*chart_x)**18)
    leading_phase = sp.cancel(chart_phase.subs(chart_r, 0))
    zero(leading_phase-108*chart_x**3/chart_a**21)
    chart_volume = sp.cancel(exact_r*chart_r**3/chart_c)
    leading_volume = sp.cancel(chart_volume.subs(chart_r, 0))
    zero(leading_volume-chart_a**-9)

    return {
        "kind": "exact limitation of a structural strict-positive-volume proof",
        "coefficient_domain": "C, H!=0",
        "coefficients": {"c": "(Z-1)^9", "h": "H Z^2 (Z-1)^6",
                         "u": "-H^2 Z^3 (Z-1)^3/6", "k": "-H^3 Z^5/54", "a": 0},
        "cube_identity": "phi=-(r Z-2h/3)^3/(2 r Z c^2)",
        "exact_chart": {"H": "3/(2 A^6)", "zeta": "1+A R",
                        "r": "zeta R^6", "Z": "zeta+R^3 X",
                        "parameter_note": "unit reparametrization of the sixth-root parameter; no approximate critical centre"},
        "critical_F_z_leading": "-6 A^-1 R^5",
        "phase_leading": "108 A^-21 X^3",
        "normalized_volume_leading": "A^-9 dX",
        "critical_action": "identically zero",
        "normalized_volume_order": 0,
        "original_moments": {"1": 0, "2": 306},
        "pure": False,
        "not_an_H21_counterexample": True,
        "does_not_assert": ["unit period", "rational seed period", "a pure unsafe phase"],
        "scope": "The W-affine shape alone does not force strict positive local volume. Extra consequences of purity may still exclude this chart.",
    }


def varying_quadratic_controls():
    """Algebraic controls for clean localization and a varying W^2 coefficient."""
    z, w, y, yp, t, r, a = sp.symbols("z w y yp t r a")
    ell, c, h, u, k = (sp.Function(name)(z) for name in ("ell", "c", "h", "u", "k"))
    d = 1-2*a*t
    precision = sp.Matrix([[-2*t*ell, -t*c], [-t*c, d]])
    determinant = -2*t*ell+(4*a*ell-c**2)*t**2
    pencil = -2*ell+t*(4*a*ell-c**2)
    zero(precision.det()-determinant)
    zero(determinant-t*pencil)
    source = sp.Matrix([t*h-z, t*u])
    nullvector = sp.Matrix([d, t*c])
    nullsource = d*(t*h-z)+t**2*c*u
    for entry in precision*nullvector-sp.Matrix([determinant, 0]):
        zero(entry)
    zero((source.T*nullvector)[0]-nullsource)
    adjugate_numerator = (source.T*precision.adjugate()*source)[0]
    zero(d*adjugate_numerator-nullsource**2-determinant*(t*u)**2)

    # Completing the nonnull square leaves the actual linear null source.
    fibre = (-(sp.Matrix([w, y]).T*precision*sp.Matrix([w, y]))[0]/2
             + (source.T*sp.Matrix([w, y]))[0]+t*k)
    completed = fibre.subs(y, yp+t*c*w/d)
    expected = -d*yp**2/2+t*u*yp+(nullsource/d)*w-determinant*w**2/(2*d)+t*k
    zero(completed-expected)

    # Primitive-pencil source obstruction, with no functional assumptions.
    ell_value, c_value, h_value, u_value = sp.symbols("E C H U")
    field_f = 4*a*ell_value-c_value**2
    field_s = (1-2*a*t)*(t*h_value-z)+t**2*c_value*u_value
    bracket = z*c_value*field_f-2*ell_value*c_value*h_value+4*ell_value**2*u_value
    zero(field_f**2*field_s.subs(t, 2*ell_value/field_f)-c_value*bracket)
    zero(bracket.subs(ell_value, 0)+z*c_value**3)

    curve_q = c**2-4*a*ell+2*ell/t
    zero(determinant+t**2*curve_q)
    zero(-sp.trace(precision.inv()*sp.diff(precision, z))/2
         + sp.diff(curve_q, z)/(2*curve_q))
    zero(-sp.trace(precision.inv()*sp.diff(precision, t))/2
         + 1/t+sp.diff(curve_q, t)/(2*curve_q))
    effective = t*k+(source.T*precision.inv()*source)[0]/2
    universal_numerator = (
        z**2+t*(-2*z*h-2*a*z**2)
        +t**2*(h**2+4*a*z*h-2*c*z*u-4*ell*k)
        +t**3*(-2*a*h**2+2*c*h*u-2*ell*u**2+2*(4*a*ell-c**2)*k))
    zero(2*determinant*effective-universal_numerator)

    # This finite-degree symbolic control illustrates the general inverse formula.
    b = sp.symbols("b", nonzero=True)
    coefficients = sp.symbols("f0:7")
    polynomial = sum(coefficient*w**power for power, coefficient in enumerate(coefficients))
    inverse = sum((-1)**j*sp.diff(polynomial, w, j)/b**(j+1) for j in range(7))
    zero(sp.diff(inverse, w)+b*inverse-polynomial)
    m = sp.symbols("m", nonzero=True)
    normal_d = sp.Matrix([[0, 0], [-m, 0]])
    normal_h = sp.Matrix([[0, -1/m], [0, 0]])
    assert normal_d*normal_h+normal_h*normal_d == sp.eye(2)

    # Exact function-field examples, including repeated ell roots.
    examples = []
    for epoly, cpoly, aval, hpoly, upoly in [
        (z**2, sp.Integer(1), 0, sp.Integer(0), sp.Integer(0)),
        ((z-1)**5, z, 1, z**3+1, z**2),
        (z**2+1, z+1, 0, z**2, z),
    ]:
        assert sp.gcd(epoly, cpoly) == 1
        jpoly = -2*epoly+t*(4*aval*epoly-cpoly**2)
        spoly = (1-2*aval*t)*(t*hpoly-z)+t**2*cpoly*upoly
        field = sp.QQ.frac_field(t)
        jpol = sp.Poly(jpoly, z, domain=field)
        spol = sp.Poly(spoly, z, domain=field)
        assert sp.gcd(jpol, spol).degree() == 0
        assert sp.gcd(jpol, jpol.diff()).degree() == 0
        examples.append({"ell": str(epoly), "c": str(cpoly), "a": aval})
    exceptional_phase = effective.subs({ell: z**2, c: 1, h: 0, u: 0, k: 0, a: 0})
    zero(sp.limit(t*exceptional_phase, t, 0)+sp.Rational(1, 4))
    assert sp.cancel(sp.diff(exceptional_phase, z)) != 0

    fixed_delta = determinant.subs(c, 0)
    zero(fixed_delta+t**2*(2/t-4*a)*ell)
    normalizer_square = r**2/(2*r-4*a)
    zero((1/fixed_delta).subs(t, 1/r)+normalizer_square/ell)
    fixed_phase = effective.subs(c, 0)
    zero(fixed_phase-(t*k-(t*h-z)**2/(4*t*ell)+t**2*u**2/(2*d)))
    fixed_edge = fixed_phase.subs({ell: z**2, h: 1, u: 0, k: 0, a: 0})
    edge_constant = sp.limit(fixed_edge+1/(4*t), t, 0)
    zero(edge_constant-1/(2*z))
    assert sp.diff(edge_constant, z) != 0

    return {
        "kind": "algebraic identity controls, not a moment-prefix search",
        "symbolic": [
            "varying precision determinant and nullvector",
            "Schur numerator d*N=S^2+det*(tu)^2",
            "rank-one fibre square completion and actual linear null source",
            "F^2 S substitution and modulo-ell obstruction",
            "both spatial and parameter half-density connections",
            "universal effective numerator",
            "finite inverse of partial_w+b on a symbolic degree-six polynomial",
            "higher-pole normal complex contraction",
            "fixed-cover scalar normalizer and phase",
        ],
        "exact_function_field_examples": examples,
        "edge_controls": ["ell=z^2,c=1: constant leading phase retains an actual spatial pole",
                          "c=0,ell=z^2,h=1: next spatial pole retained"],
        "scope": "The written clean-localization and curve-period arguments supply arbitrary-degree theorem coverage. These computations verify algebraic identities and the stated exact examples.",
    }


def shared_quadratic_source_controls():
    """Exact shared-factor source criterion and its genuine exception."""
    z, t, a, g, e, C, c, h, u = sp.symbols("z t a g e C c h u")
    d = 1-2*a*t
    f = 4*a*e-C
    residual = -2*e+t*f
    source = d*(t*h-z)+t**2*c*u
    obstruction = C*(z*f-2*e*h)+4*e**2*c*u

    # Relation g*C=c^2 is exactly the common-factor identity.
    original_delta = -2*t*g*e+t**2*(4*a*g*e-c**2)
    zero(original_delta-t*g*residual-t**2*(g*C-c**2))
    zero(original_delta.subs(c**2, g*C)-t*g*residual)

    # This gives an exact iff because the primitive pencil is irreducible.
    zero(f**2*source.subs(t, 2*e/f)-obstruction)
    zero(obstruction.subs(e, 0)+z*C**2)
    fixed_source = source.subs({e: 0, c: 0})
    zero(fixed_source-d*(t*h-z))

    # Explicit shared-factor example outside the previous B/C statements:
    # ell has a fourth-order root, while its gcd with c^2 has order two.
    ell_poly = z*(z-1)**4
    c_poly = z-1
    g_poly = sp.gcd(ell_poly, c_poly**2)
    e_poly = sp.cancel(ell_poly/g_poly)
    C_poly = sp.cancel(c_poly**2/g_poly)
    assert g_poly == (z-1)**2 or sp.expand(g_poly-(z-1)**2) == 0
    assert sp.degree(e_poly, z) >= 2
    j_poly = -2*e_poly-t*C_poly
    s_poly = -z
    r_poly = obstruction.subs({e: e_poly, C: C_poly, c: c_poly, a: 0, h: 0, u: 0})
    zero(r_poly+z)
    field = sp.QQ.frac_field(t)
    j = sp.Poly(j_poly, z, domain=field)
    s = sp.Poly(s_poly, z, domain=field)
    assert sp.gcd(j, s).degree() == 0
    assert sp.gcd(j, j.diff()).degree() == 0

    # The exception is genuine despite clean fixed roots. At alpha=1,
    # t*h(alpha)-alpha=-t/2-1 is nonzero, but the moving nullsource vanishes.
    exceptional_ell = z*(z-1)**2
    exceptional_c = z-1
    exceptional_g = sp.gcd(exceptional_ell, exceptional_c**2)
    exceptional_e = sp.cancel(exceptional_ell/exceptional_g)
    exceptional_C = sp.cancel(exceptional_c**2/exceptional_g)
    exceptional_subs = {e: exceptional_e, C: exceptional_C,
                        c: exceptional_c, a: 0, h: -sp.Rational(1, 2), u: 0}
    zero(obstruction.subs(exceptional_subs))
    exceptional_j = residual.subs(exceptional_subs)
    exceptional_s = source.subs(exceptional_subs)
    zero(exceptional_j-2*exceptional_s)
    assert (t*(-sp.Rational(1, 2))-z).subs(z, 1) != 0

    # Endpoint simplification in the finite valuation argument, for l<=2:
    # endpoint-1/2=(2-l)/(2(q-l)) is nonnegative whenever q>l.
    q, l = sp.symbols("q l")
    endpoint = (q/sp.Integer(2)+1-l)/(q-l)
    zero(endpoint-sp.Rational(1, 2)-(2-l)/(2*(q-l)))

    return {
        "status": "passed",
        "symbolic": ["determinant shared-factor identity under g*C=c^2",
                     "exact moving nullsource obstruction R",
                     "R modulo e forces e|z",
                     "fixed-root source",
                     "finite valuation endpoint bound"],
        "new_clean_example": {"ell": str(ell_poly), "c": str(c_poly),
                              "g": str(g_poly), "R": str(sp.cancel(r_poly))},
        "excluded_source_kernel_example": {
            "ell": str(exceptional_ell), "c": str(exceptional_c),
            "h": "-1/2", "a": 0, "u": 0,
            "J0": str(exceptional_j), "S": str(exceptional_s), "R": 0},
        "scope": "Exact identities and source examples, with no moment prefixes.",
    }


def moving_nullsource_zero_controls():
    """Actual-seed identities behind the written arbitrary-polynomial exclusion."""
    n = sp.symbols("n", integer=True, nonnegative=True)
    family_n = sp.symbols("family_n", integer=True, positive=True)
    # The general-j proof is coefficient extraction from
    # (1+x)^(n+2*family_n*j) (1-(1+x))^j. Both parameters stay symbolic.
    for j in range(8):
        def choose_polynomial(top, bottom):
            return sp.prod(top-index for index in range(bottom))/sp.factorial(bottom)
        difference = sum(
            (-1)**a*sp.binomial(j, a)*choose_polynomial(n+2*family_n*j+a, 2*j)
            for a in range(j+1))
        zero(difference-(-1)**j*choose_polynomial(n+2*family_n*j, j))

    def base_coefficient(marker_degree, order, nu=1):
        if order < marker_degree or (order-marker_degree) % (2*nu):
            return sp.Integer(0)
        j = (order-marker_degree)//(2*nu)
        return ((-sp.Rational(1, 2))**marker_degree*(-1)**j
                * sp.factorial2(2*j-1)*sp.Rational(1, 2)**(2*(nu-1)*j)
                * sp.binomial(marker_degree+2*nu*j, j))

    # Independent direct Wick contractions of the original three-variable P0,
    # rather than contractions of its eliminated or localized source.
    for nu in range(1, 5):
        base = sparse(Z**(2*nu+1)*W**2+Z**nu*W*Y-W/2)
        power = {(0, 0, 0): Fraction(1)}
        for order in range(11):
            for marker_degree in range(5):
                assert (wick(power, (marker_degree, 0, 0))/math.factorial(order)
                        == base_coefficient(marker_degree, order, nu))
            power = multiply(power, base)

    # The algebraic Gaussian kernel reproduces every B_n by the written
    # binomial-series identity. These controls retain several nonzero markers.
    x = sp.symbols("x")
    for marker_degree in range(5):
        kernel = 2**marker_degree/(sp.sqrt(1+4*x)
                 * (1+sp.sqrt(1+4*x))**marker_degree)
        expansion = sp.series(kernel, x, 0, 8).removeO().expand()
        for j in range(8):
            assert expansion.coeff(x, j) == (-1)**j*sp.binomial(marker_degree+2*j, j)

    # Explicit R=0 example, including the actual precision and residual source.
    z, t, w, y = sp.symbols("z t w y")
    precision = sp.Matrix([[-2*t*z**3, -t*z], [-t*z, 1]])
    source = sp.Matrix([-z-t/2, 0])
    residual = -2*z-t
    zero(precision.det()-t*z**2*residual)
    zero(source[0]-residual/2)
    zero(source[0].subs(z, -t/2))
    effective = (source.T*precision.inv()*source)[0]/2
    zero(effective+1/(4*t*z)+1/(8*z**2))

    # Rational parametrization checks the cover normal form and its density;
    # the actual-seed identity itself was derived above by Wick contraction.
    cover_u = sp.symbols("cover_u")
    cover_z = t*(cover_u**2-1)/2
    gaussian_s = cover_u/(t*(1-cover_u**2))
    cover_v = (1+cover_u**2)/(1-cover_u**2)
    zero(cover_v**2-1-4*t**2*gaussian_s**2)
    zero(sp.together(cover_z+t/(1+cover_v)))
    zero(effective.subs(z, cover_z)+gaussian_s**2/2)
    volume_ratio = (sp.diff(cover_z, cover_u)/(t*cover_z*cover_u)
                    / sp.diff(gaussian_s, cover_u))
    zero(volume_ratio+2/cover_v)

    k0, k1 = sp.symbols("k0 k1")
    p0 = Z**3*W**2+Z*W*Y-W/2
    zero(expectation(p0+k0+k1*Z)-k0)
    zero(expectation((p0+k1*Z)**2)+4+k1)
    p = sp.symbols("p", integer=True, positive=True)
    leading = sp.factorial(2*p)**2/(2**p*sp.factorial(p)**3)
    previous = sp.factorial(2*p-2)**2/(2**(p-1)*sp.factorial(p-1)**3)
    zero(sp.combsimp(previous/leading)-p/(2*(2*p-1)**2))
    for nu in range(1, 5):
        exact_ratio = (sp.Rational(1, 2)**(2*nu-1)*(2*p)*(2*p-1)/p**2
                       * sp.prod(2*nu*p-index for index in range(2*nu))
                       / sp.prod((2*nu-1)*p-index for index in range(2*nu-1)))
        leading_nu = (sp.factorial(2*p)*sp.factorial(2*nu*p)
                      / (2**((2*nu-1)*p)*sp.factorial(p)**2
                         * sp.factorial((2*nu-1)*p)))
        zero(sp.combsimp(leading_nu/leading_nu.subs(p, p-1))-exact_ratio)
        ratio_constant = (sp.Rational(2)**(3-2*nu)*(2*nu)**(2*nu)
                          / sp.Integer(2*nu-1)**(2*nu-1))
        zero(sp.limit(exact_ratio/p, p, sp.oo)-ratio_constant)
        assert ratio_constant > 0
        for marker_degree in range(5):
            for j in range(6):
                marker_ratio = (sp.Rational(1, 2)**marker_degree
                                * sp.binomial(marker_degree+2*nu*j, j)
                                / sp.binomial(2*nu*j, j))
                product_ratio = sp.prod(
                    sp.Rational(2*nu*j+index, 2*((2*nu-1)*j+index))
                    for index in range(1, marker_degree+1))
                assert marker_ratio == product_ratio and 0 < marker_ratio <= 1
    return {
        "kind": "actual-seed identities supporting a written all-order obstruction",
        "phase": "Z^(2nu+1) W^2+Z^nu WY-W/2+k(Z), integer nu>=1, arbitrary k in C[Z]",
        "conclusion": "no member of this all-degree family is pure",
        "marked_base_formula": "B_N=(-t/2)^N sum_j (-1)^j (2j-1)!! 2^(-2(nu-1)j) binom(N+2nu*j,j)t^(2nu*j)",
        "direct_Wick_controls": {"nu": [1, 4], "marker_degrees": [0, 4], "moment_orders": [0, 10]},
        "finite_difference_controls": {"N": "symbolic", "nu": "symbolic", "j": [0, 7]},
        "Gaussian_kernel_controls": {"nu": 1, "marker_degrees": [0, 4], "series_orders": [0, 7]},
        "source_exception": {"ell": "Z^3", "c": "Z", "h": "-1/2",
                             "a": 0, "u": 0, "R": 0, "J0": "-2Z-t"},
        "relative_even_coefficient_bound": "(exp(A)-1)*p/[2(2p-1)^2], A=sum_l |k_l|, k0=0",
        "general_ratio_limit": "(L_j/L_(j-1))/j -> 2^(3-2nu)(2nu)^(2nu)/(2nu-1)^(2nu-1)>0",
        "ratio_controls": {"nu": [1, 4], "j": "symbolic"},
        "scope": "The written finite-difference identity and uniform complex-coefficient estimate prove the arbitrary-degree exclusion for every nu>=1. These computations check exact algebraic identities and finite Wick controls; no prefix sufficiency is used. The displayed explicit relative bound is for nu=1.",
    }


def divisible_mixed_source_controls():
    """Phase-group and scale identities for regular h/c at common roots."""
    z, r, c, v, u, k = sp.symbols("z r c v u k")
    h = c*v
    source = (r*z-h)/c
    phase = -source**2/2+(u*source+k)/r
    groups = (-r**2*z**2/(2*c**2)+r*z*v/c
              -v**2/2+z*u/c+(k-u*v)/r)
    zero(phase-groups)
    for power, coefficient in {
        2: -z**2/(2*c**2), 1: z*v/c,
        0: -v**2/2+z*u/c, -1: k-u*v,
    }.items():
        zero(sp.expand(phase).coeff(r, power)-coefficient)

    m, N, epsilon, L, mu, A, j1 = sp.symbols("m N epsilon L mu A j1")
    # H1 and H3 combine without assumptions on the Puiseux-centre shifts.
    volume_order = 1-(m-1)*N+epsilon
    E = 1-(m-1)*N-L/mu
    zero(volume_order.subs(epsilon, -L/mu)-E)
    level_two = 2-2*m*N
    zero(E.subs({L: level_two, mu: 1})-(-1+(m+1)*N))
    zero(E.subs({L: level_two, mu: 2})-N)
    zero((-1+(m+1)*N).subs(N, 1/(m-1))-2/(m-1))

    level_zero_E = E.subs({L: -A*N, mu: 1})
    zero(level_zero_E-(1+(A-m+1)*N))
    zero(level_zero_E.subs(N, 2/(2*m-A))-(A+2)/(2*m-A))

    # For j1<=m, this remainder is nonnegative; on H2 the lower bound
    # m*N-1 is strictly positive, as its endpoint is 1/(m-1).
    level_one = 1-j1*N
    zero((level_one-level_two)-(m*N-1)-(m-j1)*N)
    zero((m*N-1).subs(N, 1/(m-1))-1/(m-1))

    # Ideal-membership form of the zero-root conditions for c|h:
    # c=z^2*C implies h=z^2*C*V, while u=z*U and k=z*K.
    # These are arbitrary functions, not selected polynomial examples.
    Cfun, Vfun, Ufun, Kfun = (sp.Function(name)(z)
                             for name in ("C", "V", "U", "K"))
    c_zero = z**2*Cfun
    h_zero = c_zero*Vfun
    u_zero = z*Ufun
    k_zero = z*Kfun
    zero(h_zero/z**2-Cfun*Vfun)
    zero(c_zero.subs(z, 0))
    zero(sp.diff(c_zero, z).subs(z, 0))
    zero(h_zero.subs(z, 0))
    zero(sp.diff(h_zero, z).subs(z, 0))
    zero(u_zero.subs(z, 0))
    zero(k_zero.subs(z, 0))

    return {
        "status": "passed",
        "symbolic": [
            "h=c*v phase with r-levels 2,1,0,-1",
            "H1/H3 volume lower bound",
            "level-two simple and double-root bounds",
            "level-zero endpoint (A+2)/(2m-A)",
            "level-one weight separation and positive H2 endpoint",
            "c in (z^2), h in (c) imply h in (z^2)",
            "zero-root conditions u,k in (z)",
        ],
        "scope": "Exact identities and ideal conditions; no moment prefixes or sampled phases.",
    }


def common_root_sharp_range_controls():
    """Exact identities behind the written all-multiplicity local range."""
    z, r = sp.symbols("sharp_z sharp_r", nonzero=True)
    c, h, u, k = sp.symbols("sharp_c sharp_h sharp_u sharp_k", nonzero=True)
    source = (r*z-h)/c
    phase = -source**2/2+(u*source+k)/r
    groups = (-r**2*z**2/(2*c**2)+r*z*h/c**2
              -h**2/(2*c**2)+z*u/c+(k-u*h/c)/r)
    zero(phase-groups)

    ell, p = sp.symbols("sharp_ell sharp_p", integer=True, positive=True)
    m = ell+p
    N = sp.symbols("sharp_N", positive=True)
    A, B, L, mu = sp.symbols("sharp_A sharp_B sharp_L sharp_mu")
    E = 1-(m-1)*N-L/mu
    level_two = 2-2*m*N
    level_one = 1-(m+p)*N
    level_zero = -A*N
    level_minus_one = -1-B*N

    # Highest-level-one faces satisfy ell*N<=1 by comparison with level two.
    # The displayed remainders prove the exact endpoint lower bounds used
    # in the written case split (ell<=2 versus ell>2).
    zero((level_one-level_two)-(ell*N-1))
    one_simple = E.subs({L: level_one, mu: 1})
    one_double = E.subs({L: level_one, mu: 2})
    zero(one_simple-(p+1)*N)
    zero(one_simple.subs(N, 1/(m-1))-(p+1)/(m-1))
    zero(one_double-(sp.Rational(1, 2)+(2-ell)*N/2))
    zero(one_double.subs(N, 1/ell)-1/ell)
    zero(one_double-1/ell-(ell-2)*(1/ell-N)/2)

    # A four-level face forces these weights exactly; no multiplicity
    # samples or finite search over root orders is used.
    forced = {N: 1/ell, A: 2*p, B: 2*p-ell}
    common_weight = -2*p/ell
    for weight in (level_two, level_one, level_zero, level_minus_one):
        zero(weight.subs(forced)-common_weight)
    zero((level_zero-level_two).subs(N, 1/ell)-(2*p-A)/ell)
    zero((level_minus_one-level_two).subs(N, 1/ell)-(2*p-ell-B)/ell)
    triple_volume = E.subs({N: 1/ell, L: common_weight, mu: 3})
    zero(triple_volume-(3-p)/(3*ell))
    zero((m-3*p)-(ell-2*p))
    zero((3*ell-2*m)-(ell-2*p))

    # When B=2p-ell>0, regularity of k forces ord(u)=ell-p.
    # The correlated four coefficients then have this exact factorization.
    alpha, lead_c, lead_h, lead_u = sp.symbols(
        "sharp_alpha sharp_C sharp_H sharp_U", nonzero=True)
    tau_power = sp.symbols("sharp_T")
    face = (-alpha**2*tau_power**3/(2*lead_c**2)
            +alpha*lead_h*tau_power**2/lead_c**2
            +(-lead_h**2/(2*lead_c**2)+alpha*lead_u/lead_c)*tau_power
            -lead_u*lead_h/lead_c)
    bracket = -tau_power*(alpha*tau_power-lead_h)/(2*lead_c**2)+lead_u/lead_c
    zero(face-(alpha*tau_power-lead_h)*bracket)
    zero(bracket.subs(tau_power, lead_h/alpha)-lead_u/lead_c)
    assert (lead_u/lead_c).is_zero is False

    # The generic cube identity is already checked by common_root_cubic_barrier.
    # For h=H*z^2*(z-1)^ell, retain ell symbolically and compute through the
    # logarithmic derivative to avoid replacing an arbitrary order by samples.
    parameter_h = sp.symbols("sharp_parameter_H", nonzero=True)
    model_h = parameter_h*z**2*(z-1)**ell
    log_derivative = sp.cancel(sp.diff(model_h, z)/model_h)
    zero(log_derivative-2/z-ell/(z-1))
    critical_r = sp.Rational(2, 3)*parameter_h*z*(z-1)**ell
    critical_fz = critical_r-sp.Rational(2, 3)*model_h*log_derivative
    zero(critical_fz/model_h
         +sp.Rational(2, 3)*((ell+1)*z-1)/(z*(z-1)))
    # This is -(2H/3) z (z-1)^(ell-1) ((ell+1)z-1),
    # expressed with a common symbolic power for exact rational reduction.
    expected_fz = (-sp.Rational(2, 3)*model_h
                   *((ell+1)*z-1)/(z*(z-1)))
    zero(critical_fz-expected_fz)

    # With r=R^(3ell), zeta-1~A*R^3 and displacement R^(2p+3)X:
    # numerator and denominator orders agree, while eta has order 3-p.
    fz_order = 3*ell-3
    displacement_order = 2*p+3
    factor_order = fz_order+displacement_order
    denominator_order = 3*ell+6*m
    zero(factor_order-(3*ell+2*p))
    zero(3*factor_order-(9*ell+6*p))
    zero(denominator_order-(9*ell+6*p))
    volume_R_order = 3*ell+displacement_order-3*m
    zero(volume_R_order-(3-p))
    zero(volume_R_order/(3*ell)-(3-p)/(3*ell))
    chart_a = sp.symbols("sharp_chart_A", nonzero=True)
    leading_cubic = -(-ell/chart_a)**3/(2*chart_a**(2*m))
    zero(sp.powsimp(leading_cubic-ell**3/(2*chart_a**(2*m+3))))

    return {
        "status": "passed",
        "symbolic": [
            "generic four-level W-affine phase",
            "highest-level-one simple/double-root volume and endpoint bounds",
            "forced four-level weights N=1/ell,A=2p,B=2p-ell",
            "triple-root volume (3-p)/(3ell)",
            "four-level factorization and nonzero bracket at the linear root",
            "symbolic-ell critical F_z identity through logarithmic derivative",
            "cubic numerator/denominator order 9ell+6p and volume order 3-p",
        ],
        "written_local_range": {
            "orders": "m=ord c, ell=ord h, p=m-ell>0 at a nonzero common root",
            "uniform_positive_if": "p<=2 or ell<2p",
            "sharp_exceptional_strip": "p>=3 and ell>=2p; equivalently ceil(2m/3)<=ell<=m-3",
            "all_multiplicities": "The exhaustive face proof establishes this range for all integer orders; it is not finite test coverage.",
            "multiplicity_eight": "m<=8 admits no exceptional pair because the strip requires m>=3p>=9",
        },
        "generic_critical_F_z": "-(2H/3) z (z-1)^(ell-1) ((ell+1)z-1)",
        "generic_chart": {"r": "R^(3ell)", "zeta": "1+A R^3+O(R^6)",
                          "displacement": "R^(2p+3) X",
                          "phase_leading": "ell^3 X^3/(2 A^(2m+3))",
                          "normalized_volume": "R^(3-p)(A^-m+O(R)) dX"},
        "representative_actual_chart": "common_root_cubic_barrier checks m=9,ell=6,p=3 exactly, including M2=306",
        "scope": "Exact symbolic identities only. The written proof supplies exhaustive faces, the all-order sharp range and complete period comparison. No finite root-order prefix is used or claimed.",
    }




def tilted_nullsource_zero_controls():
    """Actual Wick controls for the written analytic-amplitude exclusion."""
    z, w, y, t = Z, W, Y, T
    a = sp.symbols("tilted_a")
    d = 1-2*a*t

    def tilted_coefficients(nu, marker_degree, order):
        value = 0
        for j in range(max(0, (order-marker_degree)//(2*nu)+1)):
            coefficient = ((-sp.Rational(1, 2))**marker_degree*(-1)**j
                           *sp.factorial2(2*j-1)*sp.Rational(1, 2)**(2*(nu-1)*j)
                           *sp.binomial(marker_degree+2*nu*j, j))
            value += (coefficient*t**(marker_degree+2*nu*j)
                      *d**(-marker_degree-sp.Rational(3, 2)-(2*nu+1)*j))
        return sp.series(value, t, 0, order+1).removeO().expand()

    # The Gaussian tilt is a polynomial-marker identity. Its two factors
    # are j!*d^(-j-1) and (2r-1)!!*d^(-r-1/2); this verifies their symbolic
    # density and covariance exponents for arbitrary balanced monomials.
    pair_order, gaussian_order = sp.symbols(
        "tilted_pair_order tilted_gaussian_order", integer=True, nonnegative=True)
    tilt_left = (sp.factorial(pair_order)*d**(-pair_order-1)
                 *sp.factorial2(2*gaussian_order-1)
                 *d**(-gaussian_order-sp.Rational(1, 2)))
    tilt_right = (d**(-sp.Rational(3, 2))*sp.factorial(pair_order)
                  *sp.factorial2(2*gaussian_order-1)*d**(-pair_order-gaussian_order))
    assert sp.powsimp(tilt_left/tilt_right, force=True) == 1

    # Independent original three-variable Wick contractions with symbolic a.
    # These are bounded controls, separate from the unrestricted kernel proof.
    for nu in range(1, 4):
        p0 = z**(2*nu+1)*w**2+z**nu*w*y-w/2
        base = p0+a*(2*z*w+y*y)
        power = sp.Integer(1)
        marked = [sp.Integer(0)]*4
        for order in range(9):
            if order:
                power = sp.expand(power*base)
            for marker_degree in range(4):
                marked[marker_degree] += (expectation(z**marker_degree*power)
                                          *t**order/math.factorial(order))
        for marker_degree in range(4):
            zero(marked[marker_degree]-tilted_coefficients(nu, marker_degree, 8))

        lam = d**sp.Rational(1, 2*nu)
        tau = t*d**(-sp.Rational(2*nu+1, 2*nu))
        for marker_degree in range(4):
            for j in range(4):
                lhs = (d**(-marker_degree-sp.Rational(3, 2)-(2*nu+1)*j)
                       *t**(marker_degree+2*nu*j))
                rhs = d**(-sp.Rational(3, 2))*lam**marker_degree*tau**(marker_degree+2*nu*j)
                # All powers are analytic germs of the same unit d, with
                # constant coefficient one; their exponents therefore combine.
                assert sp.powsimp(lhs/rhs, force=True) == 1

    # The reparametrization exponent balance also remains symbolic in nu,N,j.
    nu = sp.symbols("tilted_nu", integer=True, positive=True)
    marker_degree, j = sp.symbols("tilted_N tilted_j", integer=True, nonnegative=True)
    left_exponent = -marker_degree-sp.Rational(3, 2)-(2*nu+1)*j
    right_exponent = (-sp.Rational(3, 2)+marker_degree/(2*nu)
                      -(2*nu+1)*(marker_degree+2*nu*j)/(2*nu))
    zero(left_exponent-right_exponent)
    reparametrization = t*d**(-(2*nu+1)/(2*nu))
    zero(reparametrization.subs(t, 0))
    zero(sp.diff(reparametrization, t).subs(t, 0)-1)

    # Exact algebraic phase identity for arbitrary polynomial v and c=z^nu.
    # Treating c and v as symbols already verifies this identity without
    # a degree bound; the polynomial Wick shift is checked separately below.
    c_symbol, v_symbol, k_symbol = sp.symbols("tilted_c tilted_v tilted_k")
    u_symbol = c_symbol*v_symbol
    residual_k = k_symbol-z*u_symbol**2
    generic_base = z*c_symbol**2*w**2+c_symbol*w*y-w/2+a*(2*z*w+y*y)
    original = (a*(2*z*w+y*y)+z*(c_symbol*w+u_symbol)**2
                +(c_symbol*w+u_symbol)*y-w/2+residual_k)
    translated = generic_base+v_symbol/2+residual_k-2*a*z*v_symbol
    zero(original-translated.subs(w, w+v_symbol))

    v0, v1, k0, k1 = sp.symbols("tilted_v0 tilted_v1 tilted_k0 tilted_k1")
    v = v0+v1*z
    k = k0+k1*z
    for nu_value in range(1, 3):
        c = z**nu_value
        u = c*v
        residual_k = k-z*u*u
        p0 = z**(2*nu_value+1)*w*w+z**nu_value*w*y-w/2
        original = a*(2*z*w+y*y)+z*(c*w+u)**2+(c*w+u)*y-w/2+residual_k
        translated = p0+a*(2*z*w+y*y)+v/2+residual_k-2*a*z*v
        zero(original-translated.subs(w, w+v))
        for marker in (sp.Integer(1), z, z*z, y*y, z*y*y):
            for w_power in range(5):
                left = expectation(marker*(w+v)**w_power)
                # Only finitely many Z powers contract this fixed W degree.
                # Thus the verifier uses polynomial markers exclusively.
                truncation = sp.series(sp.exp(z*v), z, 0, w_power+1).removeO()
                right = expectation(marker*w**w_power*truncation)
                zero(left-right)
        amplitude = sp.exp(d*z*v+t*(v/2+residual_k))
        zero(amplitude.subs({t: 0, z: 0})-1)

    # Full original residual-locus source, independently compared with its
    # translated tilted kernel, for symbolic a,v0,k0,k1 through order four.
    u = z*v0
    residual_k = k0+k1*z-z*u*u
    original = a*(2*z*w+y*y)+z*(z*w+u)**2+(z*w+u)*y-w/2+residual_k
    direct = sp.Integer(0)
    power = sp.Integer(1)
    for order in range(5):
        if order:
            power = sp.expand(power*original)
        direct += expectation(power)*t**order/math.factorial(order)
    amplitude_time = -2*a*z*v0+v0/2+residual_k
    amplitude = (sum((z*v0)**N/math.factorial(N) for N in range(5))
                 *sum((t*amplitude_time)**m/math.factorial(m) for m in range(5)))
    amplitude = sp.expand(amplitude)
    expected = sum(amplitude.coeff(z, N)*tilted_coefficients(1, N, 4)
                   for N in range(5))
    expected = sp.series(expected, t, 0, 5).removeO().expand()
    zero(direct-expected)

    return {
        "status": "passed",
        "written_theorem": {
            "integer_parameter": "every nu>=1",
            "coefficients": "a in C arbitrary, k in C[Z] arbitrary, u in Z^nu C[Z]",
            "phase": "a(2ZW+Y^2)+Z(Z^nu W+u)^2+(Z^nu W+u)Y-W/2+k-Zu^2",
            "conclusion": "no pure member for the actual ordinary-Gaussian Wick functional",
        },
        "analytic_amplitude_lemma": "For every jointly analytic A(T,Z) with A(0,0)=1, the written proof gives [T^(2nu*q)] sum_N [Z^N]A B_N=(-1)^q L_q(1+O(1/q)); this is unrestricted proof scope, not bounded test coverage.",
        "symbolic": [
            "Gaussian tilt density and covariance exponents for arbitrary balanced monomials",
            "analytic reparametrization exponent identity for symbolic nu,N,j and derivative one at zero",
            "generic polynomial phase identity under W -> W+v(Z)",
        ],
        "bounded_controls": {
            "tilted_kernel": {"nu": [1, 3], "Z_marker_degrees": [0, 3],
                              "moment_orders": [0, 8], "a": "symbolic"},
            "analytic_reparametrization_germs": {"nu": [1, 3], "N_and_j": [0, 3]},
            "polynomial_translation": {"nu": [1, 2], "v": "v0+v1 Z", "W_powers": [0, 4],
                                       "markers": ["1", "Z", "Z^2", "Y^2", "Z Y^2"]},
            "full_original_source": {"nu": 1, "v": "v0", "k": "k0+k1 Z",
                                     "a": "symbolic", "moment_orders": [0, 4]},
        },
        "actual_source_equivalence": "The polynomial Wick translation identities and the written locally finite B_N=O(T^N) pairing identify the translated actual source. No localization or evaluation of nonpolynomial markers is assumed.",
        "scope": "Exact symbolic identities and explicitly bounded original-source Wick controls. The all-nu, all-a, arbitrary-polynomial-k and u-divisible-by-Z^nu exclusion follows from the written analytic-amplitude argument, not from finite prefixes.",
        "external_review": False,
        "formal_verification": False,
    }


def checks():
    c, d = sp.symbols("c d", nonzero=True)
    ell = sp.Function("ell")(Z)
    h = c + d * Z
    v = -d - 3 * d**2 * Z / (2 * c) - d**3 * Z**2 / (2 * c**2)
    b = h**2 * ell
    k = -c**2 * Z * h * ell**2 / 2
    denominator = h - 2 * Z * v
    zero(denominator - h**3 / c**2)
    zero((h - Z * sp.diff(h, Z))**2 * denominator - h**3)
    zero(2 * denominator * k + Z * b**2)
    root = c * T / (1 - d * T)
    zero(root - T * h.subs(Z, root))
    j = 1 - d * T
    precision = 1 - 2 * T * v.subs(Z, root)
    zero(j**2 * precision - 1)

    alpha, beta, gamma = sp.symbols("alpha beta gamma")
    phase = W * h + v * Y**2 + b * Y + k
    def direction(expr):
        return (alpha * sp.diff(expr, Z) + beta * sp.diff(expr, W)
                + gamma * sp.diff(expr, Y))
    second = sp.expand(direction(direction(phase)))
    zero(second.coeff(Y, 2) + alpha**2 * d**3 / c**2)
    zero(second.subs(alpha, 0) - 2 * gamma**2 * v)

    abstract_b = sp.Function("b")(Z)
    forward_y = Y - Z * abstract_b / c
    forward_w = W + abstract_b * Y / c - Z * abstract_b**2 / (2 * c**2)
    inverse_y = Y + Z * abstract_b / c
    inverse_w = W - abstract_b * Y / c - Z * abstract_b**2 / (2 * c**2)
    zero(2 * Z * forward_w + forward_y**2 - 2 * Z * W - Y**2)
    zero(forward_y.subs({Y: inverse_y, W: inverse_w}, simultaneous=True) - Y)
    zero(forward_w.subs({Y: inverse_y, W: inverse_w}, simultaneous=True) - W)

    # Mixed constant source: phase expansion and the exact moving-coordinate bridge.
    lam = sp.symbols("lambda", nonzero=True)
    a = sp.symbols("a")
    abstract_h, abstract_u, abstract_k = (
        sp.Function(name)(Z) for name in ("h", "u", "k")
    )
    s = (Z / T - abstract_h) / lam
    phi = -s**2 / 2 + T * (a * s**2 + abstract_u * s + abstract_k)
    expanded_phi = (
        -Z**2 / (2 * lam**2 * T**2)
        + (Z * abstract_h + a * Z**2) / (lam**2 * T)
        - abstract_h**2 / (2 * lam**2)
        - 2 * a * Z * abstract_h / lam**2 + Z * abstract_u / lam
        + T * (a * abstract_h**2 / lam**2 - abstract_u * abstract_h / lam + abstract_k)
    )
    zero(phi - expanded_phi)
    j_z = 1 - T * sp.diff(abstract_h, Z)
    r_y = T * lam / j_z
    r_t = Z / (T * j_z)  # root equation: h(r)+lambda Y=r/t
    # F_Y = -Y + phi_z r_Y + Y; check the spatial Stein identity at s=Y.
    f_z = T * (sp.diff(abstract_u, Z) * Y + sp.diff(abstract_k, Z))
    f_y = T * (2 * a * Y + abstract_u) + f_z * r_y
    zero(sp.diff(phi, Z).subs(Z / T - abstract_h, lam * Y) * r_y - (f_y - Y))
    # Check r_Yt=r_tY along the root and H'/H=-1/t.
    total_t_r_y = sp.diff(r_y, T) + sp.diff(r_y, Z) * r_t
    total_y_r_t = sp.diff(r_t, Z) * r_y
    zero(total_t_r_y - total_y_r_t)

    # General mixed source: spatial relation before logarithmic specialization.
    abstract_c = sp.Function("c")(Z)
    general_s = (Z / T - abstract_h) / abstract_c
    general_phi = -general_s**2 / 2 + T * (
        a * general_s**2 + abstract_u * general_s + abstract_k
    )
    general_j = 1 - T * (sp.diff(abstract_h, Z) + sp.diff(abstract_c, Z) * Y)
    general_r_y = T * abstract_c / general_j
    general_f_y = T * (2 * a * Y + abstract_u) + f_z * general_r_y
    zero(sp.diff(general_phi, Z).subs(
        Z / T - abstract_h, abstract_c * Y
    ) * general_r_y - (general_f_y - Y))
    # Explicit time-volume identity for a general amplitude q(z).
    amp = sp.Function("Q")(Z)
    general_r_t = Z / (T * general_j)
    volume = 1 / (T * abstract_c)
    # With F_t along the root = phi_t + phi_z*r_t, cancellation reduces to this.
    lhs_remainder = (
        (sp.diff(volume, T) + sp.diff(volume, Z) * general_r_t)
        * general_r_y * amp
        + volume * (sp.diff(general_r_y, T)
                    + sp.diff(general_r_y, Z) * general_r_t) * amp
        + volume * general_r_y * sp.diff(amp, Z) * general_r_t
        + volume * general_r_y * amp / T
    )
    rhs_remainder = sp.diff(volume * general_r_t * amp, Z) * general_r_y + (
        sp.diff(volume * general_r_t * amp, Y)
    )
    zero(lhs_remainder - rhs_remainder)

    # Root-zero tilt: exact reduction of the transformed denominator.
    eta, gamma, u0 = sp.symbols("eta gamma u0")
    delta_root = sp.symbols("sqrt_delta", nonzero=True)
    mu = T * u0 / delta_root**2
    denom = 1 - T * eta - T * gamma * mu
    qtilt = T * gamma / (delta_root * denom)
    zero(1 - T * (eta + gamma * (mu + Y / delta_root)) - denom * (1 - qtilt * Y))

    # Fixed logarithmic volume: partial fractions at several simple punctures.
    squarefree_cases = [Z*(Z-1), (Z+2)*(Z-1)*(Z-3),
                        (Z+3)*(Z+1)*(Z-2)*(Z-4), Z*(Z+1)*(Z+2)*(Z-2)*(Z-3)]
    for polynomial in squarefree_cases:
        roots = sp.polys.polytools.ground_roots(polynomial, Z)
        assert all(multiplicity == 1 for multiplicity in roots.values())
        zero(1/polynomial - sum(1/(sp.diff(polynomial, Z).subs(Z, root) * (Z-root))
                               for root in roots))
    # A repeated-root rejection control: this nonzero phase is pure by weight.
    assert sp.discriminant(Z**2, Z) == 0
    repeated_phase = sparse(Z**2 * W * Y)
    repeated_power = {(0, 0, 0): Fraction(1)}
    for m in range(1, 9):
        repeated_power = multiply(repeated_power, repeated_phase)
        assert wick(repeated_power) == 0

    # Repeated-square reciprocal phase: exact top coefficient and t-infinity powers.
    xr, reciprocal_r = sp.symbols("x r", nonzero=True)
    kappa = sp.symbols("kappa", nonzero=True)
    root_alpha, h0, h1, h2, h3, h4 = sp.symbols("root_alpha h0 h1 h2 h3 h4")
    u0x, u1x, u2x, k0x, k1x = sp.symbols("u0x u1x u2x k0x k1x")
    reciprocal_s = ((root_alpha*reciprocal_r-h0)*xr**2
                    + (reciprocal_r-h1)*xr - h2 - h3/xr - h4/xr**2) / kappa
    reciprocal_phi = sp.expand(-reciprocal_s**2/2 + (
        a*reciprocal_s**2 + (u0x+u1x/xr+u2x/xr**2)*reciprocal_s + k0x+k1x/xr
    )/reciprocal_r)
    top_reciprocal = -(1-2*a/reciprocal_r)*(root_alpha*reciprocal_r-h0)**2/(2*kappa**2)
    zero(reciprocal_phi.coeff(xr, 4) - top_reciprocal)
    # Multiplication by r clears every parameter pole, irrespective of x powers.
    sp.Poly(sp.expand(reciprocal_r*xr**4*reciprocal_phi), reciprocal_r, xr)
    # An additional polynomial factor is a unit at this selected double root.
    q0 = sp.symbols("q0", nonzero=True)
    q1, q2 = sp.symbols("q1 q2")
    extra_factor = q0 + q1/xr + q2/xr**2
    general_double_s = reciprocal_s / extra_factor
    general_double_phi = -general_double_s**2/2 + (
        a*general_double_s**2 + (u0x+u1x/xr+u2x/xr**2)*general_double_s + k0x+k1x/xr
    )/reciprocal_r
    zero(sp.limit(general_double_phi/xr**4, xr, sp.oo) - top_reciprocal/q0**2)
    nu, vv = sp.symbols("u v", nonzero=True)
    moving_centre = sp.Function("xi")(nu*vv)
    zero(nu*sp.diff(moving_centre, nu) - vv*sp.diff(moving_centre, vv))

    # Algebraic coefficients used in the polynomial quadratic-form orbit exclusion.
    fcoord, gcoord, leading, constant = sp.symbols("F G A B")
    quadform = Y**2 + 2*Z*W
    transformed = (1+fcoord) * (quadform-gcoord**2)/(2*fcoord) - (1+fcoord)*(2+fcoord)*gcoord**2/2
    zero(transformed - (1+fcoord)/(2*fcoord) * (quadform-(1+fcoord)**2*gcoord**2))
    coeff_y2 = sp.expand(transformed.subs(gcoord, leading*Y+constant)).coeff(Y, 2)
    zero(coeff_y2 - (1+fcoord)*(1-(1+fcoord)**2*leading**2)/(2*fcoord))
    orbit_f = leading*Y + constant
    orbit_h = Y/(2*leading) + sp.symbols("E")
    orbit_polynomial = (1+orbit_f)*orbit_h - (1+orbit_f)*(2+orbit_f)*gcoord**2/2
    zero(sp.expand(orbit_polynomial).coeff(Y, 2) - (1-leading**2*gcoord**2)/2)

    alpha, beta, radial_a = sp.symbols("alpha beta radial_a")
    a_numerator = sp.expand(alpha*T*(1-2*radial_a*T) + beta**2*T**2/2)
    zero(a_numerator.coeff(T, 1) - alpha)
    zero(a_numerator.subs(alpha, 0).coeff(T, 2) - beta**2/2)
    radial_result = radial_controls()
    endpoint_result = relative_endpoint_controls()
    resonant_result = resonant_quartic_certificate()
    barrier_result = common_root_cubic_barrier()
    varying_result = varying_quadratic_controls()
    shared_result = shared_quadratic_source_controls()
    divisible_result = divisible_mixed_source_controls()
    sharp_range_result = common_root_sharp_range_controls()
    moving_nullsource_result = moving_nullsource_zero_controls()
    tilted_nullsource_result = tilted_nullsource_zero_controls()

    # Generic elimination of W,Y, with constant W^2 coefficient nonzero.
    lam_quad = sp.symbols("lambda_quad", nonzero=True)
    c_quad, a_quad = sp.symbols("c_quad a_quad")
    precision_wy = sp.Matrix([[-2*lam_quad*T, -c_quad*T], [-c_quad*T, 1-2*a_quad*T]])
    delta_wy = -2*lam_quad*T+(4*lam_quad*a_quad-c_quad**2)*T**2
    zero(precision_wy.det()-delta_wy)
    wy_source = sp.Matrix([T*abstract_h-Z, T*abstract_u])
    effective = T*abstract_k+(wy_source.T*precision_wy.inv()*wy_source)[0]/2
    numerator_wy = (
        Z**2+T*(-2*Z*abstract_h-2*a_quad*Z**2)
        + T**2*(abstract_h**2+4*a_quad*Z*abstract_h-2*c_quad*Z*abstract_u-4*lam_quad*abstract_k)
        + T**3*(-2*a_quad*abstract_h**2+2*c_quad*abstract_h*abstract_u
                -2*lam_quad*abstract_u**2+2*(4*lam_quad*a_quad-c_quad**2)*abstract_k)
    )
    zero(2*delta_wy*effective-numerator_wy)
    # Spatially varying WY coefficient: half-density connection on the double cover.
    variable_precision = precision_wy.subs(c_quad, abstract_c)
    variable_delta = delta_wy.subs(c_quad, abstract_c)
    curve_q = abstract_c**2-4*lam_quad*a_quad+2*lam_quad/T
    zero(variable_delta+T**2*curve_q)
    zero(-sp.trace(variable_precision.inv()*sp.diff(variable_precision, Z))/2
         + sp.diff(curve_q, Z)/(2*curve_q))
    zero(-sp.trace(variable_precision.inv()*sp.diff(variable_precision, T))/2
         + 1/T + sp.diff(curve_q, T)/(2*curve_q))
    variable_effective = effective.subs(c_quad, abstract_c)
    zero(sp.limit(T*variable_effective, T, 0)+Z**2/(4*lam_quad))
    curve_r = sp.symbols("curve_r")
    for cpoly in (Z, Z**2, Z**3+Z, Z**4-Z):
        assert sp.discriminant(cpoly**2+2*curve_r, Z) != 0
    # Pure W^2 is the threshold-equality control for the generic elimination.
    zero(numerator_wy.subs({abstract_h: 0, abstract_u: 0, abstract_k: 0,
                           a_quad: 0, c_quad: 0})-Z**2)

    # Exact all-order polynomials from the shear; direct Wick controls are bounded.
    shear_phase = W + Z**2 * Y - Z**5 / 2
    transformed_w = W - Z**2 * Y - Z**5 / 2
    marker_periods = {
        q: expectation((transformed_w**q).subs(Z, Z + T))
        for q in range(1, 5)
    }
    expected = {
        1: -T**5 / 2,
        2: T**10 / 4 - 4 * T**4,
        3: -T**15 / 8 + 6 * T**9 - 18 * T**3,
        4: T**20 / 16 - 6 * T**14 + 84 * T**8 - 48 * T**2,
    }
    for q in marker_periods:
        zero(marker_periods[q] - expected[q])
        assert sp.degree(marker_periods[q], T) == 5 * q
    power = {(0, 0, 0): Fraction(1)}
    base = sparse(shear_phase)
    for m in range(33):
        if m:
            assert wick(power) == 0
        for q, period in marker_periods.items():
            assert wick(power, (0, q, 0)) == (
                sp.expand(period).coeff(T, m) * math.factorial(m)
            )
        power = multiply(power, base)

    # Long's phase: pure through the stated control range, one fixed unsafe marker.
    long_phase = (1 + Z) * (W - (2 + Z) * Y**2 / 2)
    power = {(0, 0, 0): Fraction(1)}
    base = sparse(long_phase)
    for m in range(9):
        assert wick(power) == (1 if m == 0 else 0)
        assert wick(power, (1, 0, 0)) == (0 if m == 0 else math.factorial(m))
        power = multiply(power, base)

    # A positive-weight example with nonconstant v.
    weight_phase = W * (Z**2 + Z**4) + Z**3 * Y**2 + (Z + Z**2) * Y + Z**5
    power = {(0, 0, 0): Fraction(1)}
    base = sparse(weight_phase)
    for m in range(8):
        assert wick(power) == (1 if m == 0 else 0)
        if m > 2:
            assert wick(power, (0, 2, 0)) == 0
        power = multiply(power, base)

    # Multiple double roots are allowed by the zero-root positive-weight branch.
    mixed_positive = sparse(W*(Z**2 + Z**2*(Z-1)**2*Y) + Z*Y + Z)
    mixed_power = {(0, 0, 0): Fraction(1)}
    for m in range(1, 9):
        mixed_power = multiply(mixed_power, mixed_positive)
        assert wick(mixed_power) == 0
        if m > 2:
            assert wick(mixed_power, (0, 2, 0)) == 0

    # Rejection controls: changing the shear correction or adding constant Y².
    assert expectation((W + Z**2 * Y - Z**5 / 3)**6) != 0
    assert expectation(Y**2 + Z * W) != 0
    return {
        "symbolic": [
            "arbitrary-ell unsafe-branch rational identities",
            "root equation and determinant cancellation",
            "directional second derivative and isotropic-axis obstruction",
            "shear inverse and Gaussian quadratic-form preservation",
            "mixed constant-source phase expansion and moving-coordinate identities",
            "general mixed-source spatial and time-volume identities",
            "root-zero Gaussian tilt denominator identity",
            "multipuncture squarefree partial-fraction identities",
            "repeated-square reciprocal quartic and moving-centre derivative identities",
            "selected double-root leading coefficient with additional root factors",
            "Long polynomial-orbit coefficient identities",
            "radial quartic quadratic-coefficient numerator",
            "relative endpoint boundary-class parameter action",
            "constant W^2 generic precision, determinant and universal phase numerator",
            "hyperelliptic spatial and parameter half-density connections and nonconstant phase",
            "four exact all-order marker polynomials from the shear",
        ],
        "bounded_controls": {
            "shear_B_z_squared": {"moment_orders": [0, 32], "W_marker_powers": [1, 4]},
            "long_phase": {"moment_orders": [0, 8], "markers": ["1", "Z"]},
            "positive_weight_example": {"moment_orders": [0, 7], "markers": ["1", "W^2"]},
            "rejection_controls": 2,
            "squarefree_partial_fractions": {"degrees": [2, 3, 4, 5]},
            "repeated_root_pure_rejection_control": {"phase": "Z^2 W Y", "moment_orders": [1, 8]},
            "multiple_double_root_positive_control": {"c": "Z^2 (Z-1)^2", "moment_orders": [1, 8], "markers": ["1", "W^2"]},
            "radial_quartic_formula": radial_result,
            "all_degree_radial_relative_bridge": endpoint_result,
            "generic_hyperelliptic_smoothness_controls": {"c": ["Z", "Z^2", "Z^3+Z", "Z^4-Z"]},
        },
        "marker_periods_B_z_squared": {str(q): str(p) for q, p in marker_periods.items()},
        "complete_certificates": {"resonant_quartic_family": resonant_result},
        "method_barriers": {"common_root_cubic_chart": barrier_result},
        "symbolic_identity_controls": {"varying_quadratic_source": varying_result,
                                       "shared_quadratic_source": shared_result,
                                       "regular_h_over_c": divisible_result,
                                       "common_root_sharp_range": sharp_range_result,
                                       "moving_nullsource_zero": moving_nullsource_result,
                                       "tilted_nullsource_zero": tilted_nullsource_result},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = checks()
    here = Path(__file__).resolve()
    note = here.with_name("GAUSSIAN_AFFINE_ISOTROPIC_CHART.md")
    receipt = {
        "status": "passed",
        "scope": "Exact symbolic identities, bounded Wick controls, the stated complete resonant-quartic coefficient certificate, actual-seed controls for the written all-degree nullsource-zero exclusion, and an exact nonpure local-method barrier. The certificate forces an all-order safe form in its fixed family; arbitrary-degree exclusions and analytic classifications use the written proofs rather than finite controls. The barrier is not a unit-period or H(2,1) counterexample.",
        "checked_on": "2026-11-17",
        "command": ".venv/bin/python research/moments/vanishing/verify_gaussian_affine_isotropic_chart.py --output research/moments/vanishing/gaussian_affine_isotropic_chart_verification.json",
        "environment": {"python": platform.python_version(), "sympy": sp.__version__},
        "checker_sha256": hashlib.sha256(here.read_bytes()).hexdigest(),
        "proof_sha256": hashlib.sha256(note.read_bytes()).hexdigest(),
        "checks": result,
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
    }
    if args.output:
        args.output.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
