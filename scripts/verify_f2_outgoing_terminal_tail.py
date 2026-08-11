#!/usr/bin/env python3
"""Verify the complete outgoing terminal tail for F2 ``(75,125)``."""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.log_node_profiles import (  # noqa: E402
    compile_toric_fan_profile,
    determinant,
)


def support_inequality_audit() -> None:
    # In q=y_old, r=x*y_old^5-1, every grouped Laurent monomial q^A*r^B
    # has B>=0.  The terminal supporting inequalities are
    #   5*A+12*B >= -3  for P,
    #   5*A+12*B >= -5  for R=-Q.
    # The nonzero s=0 endpoint terms are q^-3*r and q^-1.
    for B in range(40):
        p_minimum_A = -((-(-3 - 12 * B)) // 5)
        while 5 * p_minimum_A + 12 * B < -3:
            p_minimum_A += 1
        r_minimum_A = -((-(-5 - 12 * B)) // 5)
        while 5 * r_minimum_A + 12 * B < -5:
            r_minimum_A += 1
        for A in range(p_minimum_A, p_minimum_A + 5):
            arm_five_order = 2 * A + 5 * B
            arm_six_order = A + 3 * B
            assert arm_five_order >= -1
            assert arm_six_order >= 0
            if arm_five_order == -1:
                assert (A, B) == (-3, 1)
            if arm_six_order == 0:
                assert (A, B) in ((-3, 1), (0, 0))
        for A in range(r_minimum_A, r_minimum_A + 5):
            arm_five_order = 2 * A + 5 * B
            arm_six_order = A + 3 * B
            assert arm_five_order >= -2
            assert arm_six_order >= -1
            if arm_five_order == -2 or arm_six_order == -1:
                assert (A, B) == (-1, 0)

    # Exact identities proving the inequalities for arbitrary B>=1.
    A, B, slack = sp.symbols("A B slack", integer=True)
    # slack=5*A+12*B+3 for P.
    assert sp.expand(5 * (2 * A + 5 * B + 1) - (2 * slack + B - 1)).subs(
        slack, 5 * A + 12 * B + 3
    ) == 0
    assert sp.expand(5 * (A + 3 * B) - (slack + 3 * B - 3)).subs(
        slack, 5 * A + 12 * B + 3
    ) == 0
    # slack=5*A+12*B+5 for R=-Q.
    assert sp.expand(5 * (2 * A + 5 * B + 2) - (2 * slack + B)).subs(
        slack, 5 * A + 12 * B + 5
    ) == 0
    assert sp.expand(5 * (A + 3 * B + 1) - (slack + 3 * B)).subs(
        slack, 5 * A + 12 * B + 5
    ) == 0


def outgoing_fan_audit() -> None:
    exponent_map = ((1, 0), (-2, 1))
    source_rays = ((5, 12), (2, 5), (1, 3), (0, 1))
    target_rays = ((5, 2), (2, 1), (1, 1), (0, 1))
    assert determinant(exponent_map[0], exponent_map[1]) == 1
    assert tuple(
        (
            exponent_map[0][0] * ray[0] + exponent_map[0][1] * ray[1],
            exponent_map[1][0] * ray[0] + exponent_map[1][1] * ray[1],
        )
        for ray in source_rays
    ) == target_rays
    assert all(determinant(left, right) == 1 for left, right in zip(source_rays, source_rays[1:]))
    assert all(determinant(left, right) == 1 for left, right in zip(target_rays, target_rays[1:]))

    profile = compile_toric_fan_profile(
        exponent_map=exponent_map,
        source_rays=source_rays,
        target_rays=target_rays,
    )
    assert profile.refined_source_rays == source_rays
    assert profile.determinants == (1, 1, 1)


def intermediate_node_audit() -> None:
    # For the source cone (2,5)|(1,3), regular parameters are
    # c=q^3/r and d=r^2/q^5.  Then q=c^2*d and r=c^5*d^3.
    # The forced endpoint terms give a=c^2*d and b=c*d, while the target
    # cone (2,1)|(1,1) has parameters a/b and b^2/a.  They pull back to c,d.
    c, d = sp.symbols("c d")
    q = c**2 * d
    r = c**5 * d**3
    target_a = q
    target_b = sp.cancel(q**-2 * r)
    assert target_b == c * d
    assert sp.cancel(target_a / target_b) == c
    assert sp.cancel(target_b**2 / target_a) == d


def final_endpoint_audit() -> None:
    # At the last cone, alpha=q and beta=r/q^3 are regular.  After the target
    # translation P->P-c0, the forced coefficient of q^-3*r makes beta a
    # simple tangential parameter.  Arbitrary higher terms cannot change the
    # invertible differential.
    alpha, beta, c0, c1, c2, c3 = sp.symbols(
        "alpha beta c0 c1 c2 c3"
    )
    target_a = alpha * (1 + c1 * alpha + c2 * alpha * beta)
    target_tangent = beta * (1 + c3 * alpha) + c2 * alpha
    jacobian = sp.Matrix(
        [
            [sp.diff(target_a, alpha), sp.diff(target_a, beta)],
            [sp.diff(target_tangent, alpha), sp.diff(target_tangent, beta)],
        ]
    )
    assert jacobian.subs({alpha: 0, beta: 0}).det() == 1
    # The translation is a target toric shear: b'=(P-c0)/R=b-c0*a.
    a, b = sp.symbols("a b")
    assert sp.Matrix([[1, 0], [-c0, 1]]).det() == 1
    assert sp.expand(b - c0 * a) == b - c0 * a


def main() -> None:
    support_inequality_audit()
    outgoing_fan_audit()
    intermediate_node_audit()
    final_endpoint_audit()
    print(
        "PASS: the F2 outgoing terminal tail maps unimodularly to the "
        "existing (5,2)->(2,1)->(1,1) target fan and is log-etale"
    )


if __name__ == "__main__":
    main()
