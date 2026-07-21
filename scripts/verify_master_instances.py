#!/usr/bin/env python3
"""Exact map, polynomiality, Jacobian, and collision checks for displayed cases."""
import sympy as sp

from master_cancellation import DISPLAYED_INSTANCES, hensel_jet, parameter_polynomial

x, y, z, A, q = sp.symbols("x y z A q")

displayed_moduli = {
    (instance.m, instance.r): parameter_polynomial(instance.m, instance.r, q)
    for instance in DISPLAYED_INSTANCES
}
assert sp.discriminant(displayed_moduli[(2, 1)], q) == -8
assert sp.discriminant(displayed_moduli[(1, 2)], q) == -15
quartic = displayed_moduli[(2, 2)]
assert sp.Poly(quartic, q).is_irreducible
quartic_group, _ = sp.polys.numberfields.galois_group(quartic, q)
assert quartic_group.order() == 8

for instance in DISPLAYED_INSTANCES:
    m, r = instance.m, instance.r
    modulus = parameter_polynomial(m, r, q)
    hA = hensel_jet(m, r, A, q)
    actual_A = 1 + x*y**m
    actual_h = hA.subs(A, actual_A)
    B = actual_A**(r+1)*z + y**(m+1)*actual_h
    P = sp.expand(actual_A*B)
    Q = sp.expand(y + x*B)

    # Expand the integral first, then substitute its polynomial upper limit.
    t = sp.symbols("t")
    integrand = sp.expand((1 - t*(Q-P*t)**m)**r)
    anti = sp.integrate(integrand, t)
    rational_R = sp.cancel(anti.subs(t, x/actual_A) - anti.subs(t, 0))
    numerator, denominator = sp.fraction(rational_R)
    reduced_numerator = sp.rem(
        sp.Poly(numerator, q), sp.Poly(modulus, q)
    ).as_expr()
    reduced_denominator = denominator
    polynomial_R = sp.cancel(reduced_numerator/reduced_denominator)
    assert sp.denom(polynomial_R) == 1

    # Reduce the exact Jacobian modulo the defining parameter polynomial.
    determinant = sp.cancel(
        sp.Matrix([P, Q, polynomial_R]).jacobian((x, y, z)).det() + 1
    )
    det_num = sp.fraction(determinant)[0]
    assert sp.rem(sp.Poly(det_num, q), sp.Poly(modulus, q)).is_zero

    # A full n-point collision at (P,Q,R)=(1,0,0), over the splitting field
    # of E_{m,r}(U).  It is certified through the exact factorization in T.
    T, U = sp.symbols("T U")
    collision = sp.integrate((1 - (-1)**m*T**(m+1))**r, T)
    collision = sp.expand(collision - collision.subs(T, 0))
    E = sp.cancel(collision/T).subs(T**(m+1), U)
    assert sp.Poly(E, U).degree() == r
    assert sp.factor(collision - T*E.subs(U, T**(m+1))) == 0
    assert sp.gcd(sp.Poly(collision, T), sp.Poly(sp.diff(collision, T), T)) == 1

    print(
        f"PASS {instance.label}: (m,r)=({m},{r}), "
        f"deg M={m*r}, deg fiber={r*(m+1)+1}"
    )
