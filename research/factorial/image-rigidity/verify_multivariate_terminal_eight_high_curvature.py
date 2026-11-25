#!/usr/bin/env python3
"""Prove terminal-eight angular separation for total curvature at least 15/2."""

from fractions import Fraction as Q

from verify_multivariate_terminal_six_joint_ray_bands import INNER, OUTER, radial_replay
from verify_multivariate_terminal_six_private_angles import coefficient_bounds, require


CURVATURE = Q(15, 2)
COSINES = (
    Q(4511, 5000), Q(687, 1000), Q(3037, 10000), Q(-841, 10000),
    Q(-5081, 10000), Q(-799, 1000),
)


def verify_anchor():
    knots = (Q(1), *COSINES, Q(-1))
    require(all(left > right for left, right in zip(knots, knots[1:])), "ray order")
    expected = {
        7: (1, 1, -1, -1, 1, 1, -1, -1),
        8: (1, -1, -1, 1, 1, -1, -1, 1),
    }
    expressions = {
        7: (8*x**3 - 4*x**2 - 4*x + 1 for x in knots),
        8: (8*x**4 - 8*x**2 + 1 for x in knots),
    }
    for order in (7, 8):
        require(all(sign * value > 0 for sign, value
                    in zip(expected[order], expressions[order], strict=True)),
                ("anchor signs", order))


def main():
    verify_anchor()
    total = 0
    for order in (7, 8):
        caps = coefficient_bounds(order, (), CURVATURE)
        require(sum(a*INNER**k for k, a in enumerate(caps) if k) < 1,
                ("inner radius", order))
        require(sum(a*OUTER**k for k, a in enumerate(caps[:-1])) < OUTER**order,
                ("outer radius", order))
        boxes = [(Q(1), Q(1)) if k in (0, order) else (Q(0), a)
                 for k, a in enumerate(caps)]
        for cosine in COSINES:
            cells = radial_replay(order, boxes, (cosine, cosine))
            require(cells is not None, ("unresolved ray", order, cosine))
            total += cells
    print(f"Verified {total} terminal-eight radial cells.")


if __name__ == "__main__":
    main()
