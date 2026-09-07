#!/usr/bin/env sage
"""Search the second q=4 pencil on the exact CM24 first child.

The first q=4 neighbor is ``U=(x-T)/T^2``.  Its CM24 Jacobian has fibers
I5* at infinity and I6 at U=-3/2.  The second lattice divisor is again an
old-degree-two, zero-MW-projection class in ``L(2O+4F)``.  This bounded exact
diagnostic therefore tests

    W = (X+f(v))/v^2,       v=U+3/2,

where ``deg(f)<=4`` and the ``v^2`` coefficient is fixed to zero by the
translation gauge.  It solves every square-clearing chart compatible with a
binary quartic.  A hit is evidence only after its fiber signature and local
A4+D9 component branches match the transported divisor.
"""

from sage.all import PolynomialRing, QQ, ZZ


# Exact CM24 first-child Jacobian, derived from U=(x-T)/T^2.
old_base = PolynomialRing(QQ, "U")
U = old_base.gen()
A1 = -27*(
    U**6-QQ(9)/2*U**4-QQ(135)/8*U**3-QQ(81)/8*U**2
    -QQ(1701)/32*U+QQ(5103)/64
)
B1 = 54*(
    U**9-QQ(27)/4*U**7-QQ(405)/16*U**6-QQ(243)/32*U**5
    -QQ(8019)/64*U**4-QQ(12393)/32*U**3
    -QQ(229635)/256*U**2+QQ(19683)/512*U-QQ(925101)/1024
)
Delta1 = 4*A1**3+27*B1**2
assert tuple((factor.degree(), ZZ(exponent)) for factor, exponent in Delta1.factor()) == (
    (1, 6), (2, 2), (3, 1)
)


def outside_coefficients(polynomial, low, high):
    return tuple(
        coefficient
        for degree, coefficient in enumerate(polynomial.list())
        if coefficient and not low <= degree <= high
    )


def fiber_signature(A, B, Delta):
    signature = []
    root_rank = ZZ(0)
    for factor, delta_order in Delta.factor():
        a_order = A.valuation(factor)
        b_order = B.valuation(factor)
        degree = factor.degree()
        if a_order == 0 or b_order == 0:
            fiber = f"I{delta_order}"
            local_rank = max(ZZ(0), ZZ(delta_order)-1)
        elif a_order == 2 and b_order == 3 and delta_order >= 6:
            fiber = f"I{delta_order-6}*"
            local_rank = ZZ(delta_order)-2
        elif a_order >= 3 and b_order == 4 and delta_order == 8:
            fiber = "IV*"
            local_rank = ZZ(6)
        elif a_order == 3 and b_order >= 5 and delta_order == 9:
            fiber = "III*"
            local_rank = ZZ(7)
        elif a_order >= 4 and b_order == 5 and delta_order == 10:
            fiber = "II*"
            local_rank = ZZ(8)
        else:
            raise RuntimeError(
                f"unclassified finite valuations {(a_order,b_order,delta_order)}"
            )
        signature.append((str(factor.monic()), degree, a_order, b_order, ZZ(delta_order), fiber))
        root_rank += degree*local_rank
    infinity = (ZZ(8-A.degree()), ZZ(12-B.degree()), ZZ(24-Delta.degree()))
    a_order, b_order, delta_order = infinity
    if a_order == 0 or b_order == 0:
        fiber = f"I{delta_order}"
        local_rank = max(ZZ(0), delta_order-1)
    elif a_order == 2 and b_order == 3 and delta_order >= 6:
        fiber = f"I{delta_order-6}*"
        local_rank = delta_order-2
    elif a_order >= 3 and b_order == 4 and delta_order == 8:
        fiber = "IV*"
        local_rank = ZZ(6)
    elif a_order == 3 and b_order >= 5 and delta_order == 9:
        fiber = "III*"
        local_rank = ZZ(7)
    elif a_order >= 4 and b_order == 5 and delta_order == 10:
        fiber = "II*"
        local_rank = ZZ(8)
    else:
        raise RuntimeError(f"unclassified infinity valuations {infinity}")
    signature.append(("infinity", 1, a_order, b_order, delta_order, fiber))
    root_rank += local_rank
    return root_rank, tuple(signature)


hits = []
for clearing_order in range(1, 4):
    coefficient_ring = PolynomialRing(
        QQ, names=("c0", "c1", "c3", "c4"), order="degrevlex"
    )
    c0, c1, c3, c4 = coefficient_ring.gens()
    local = PolynomialRing(coefficient_ring, "v")
    v = local.gen()
    base_change = v-QQ(3)/2
    A = local(A1(base_change))
    B = local(B1(base_change))
    f = c0+c1*v+c3*v**3+c4*v**4
    g = v**2
    low = 2*clearing_order
    high = low+4
    polynomials = (
        g**3,
        -3*g**2*f,
        3*g*f**2+A*g,
        -f**3-A*f+B,
    )
    equations = tuple(
        coefficient
        for polynomial in polynomials
        for coefficient in outside_coefficients(polynomial, low, high)
    )
    ideal = coefficient_ring.ideal(equations)
    solutions = ()
    dimension = -1 if ideal.is_one() else ideal.dimension()
    if not ideal.is_one() and dimension == 0:
        solutions = tuple(ideal.variety(ring=QQ))
    print(
        f"Q80SECONDQ4RR|clearing={clearing_order}|variables=4|"
        f"equations={len(equations)}|dimension={dimension}|"
        f"rational_solutions={len(solutions)}",
        flush=True,
    )
    for solution in solutions:
        f_solution = local(
            QQ(solution[c0])+QQ(solution[c1])*v
            +QQ(solution[c3])*v**3+QQ(solution[c4])*v**4
        )
        quartic_coefficients = PolynomialRing(QQ, "W")
        W = quartic_coefficients.gen()
        quartic_ring = PolynomialRing(quartic_coefficients, "v")
        vv = quartic_ring.gen()
        quartic, remainder = (
            (vv**2*W-quartic_ring(f_solution))**3
            +quartic_ring(A1(vv-QQ(3)/2))*(vv**2*W-quartic_ring(f_solution))
            +quartic_ring(B1(vv-QQ(3)/2))
        ).quo_rem(vv**(2*clearing_order))
        assert remainder == 0 and quartic.degree() <= 4
        q0, q1, q2, q3, q4 = [quartic[index] for index in range(5)]
        invariant_i = 12*q4*q0-3*q3*q1+q2**2
        invariant_j = (
            72*q4*q2*q0+9*q3*q2*q1-27*q4*q1**2
            -27*q3**2*q0-2*q2**3
        )
        A2 = -27*invariant_i
        B2 = -27*invariant_j
        Delta2 = 4*A2**3+27*B2**2
        root_rank, signature = fiber_signature(A2, B2, Delta2)
        factors = tuple(
            (str(factor.monic()), factor.degree(), ZZ(exponent))
            for factor, exponent in Delta2.factor()
        )
        hits.append((clearing_order, f_solution, quartic, factors))
        print(
            f"Q80SECONDQ4RR|hit|clearing={clearing_order}|f={f_solution}|"
            f"quartic={quartic}|Delta_factors={factors}|root_rank={root_rank}|"
            f"fibers={signature}",
            flush=True,
        )

print(
    f"Q80SECONDQ4RR|hits={len(hits)}|status=PASS_BOUNDED_ANSATZ",
    flush=True,
)
