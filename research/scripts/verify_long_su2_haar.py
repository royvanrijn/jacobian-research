#!/usr/bin/env python3
"""Exact symbolic certificate for the local SU(2) Haar proof.

The proof identifies SU(2) with the unit 3-sphere, checks that left
multiplication is orthogonal, and computes the normalized volume density in
Hopf coordinates.  Together with the monomial calculation in
verify_long_xz_mathieu.py, this removes the former conditional dependence on
the Mueger--Tuset integration formula.
"""

import sympy as sp


def main() -> None:
    # Unit-quaternion left multiplication preserves the Euclidean metric.
    a, b, c, d = sp.symbols("a b c d", real=True)
    left = sp.Matrix(
        [
            [a, -b, -c, -d],
            [b, a, -d, c],
            [c, d, a, -b],
            [d, -c, b, a],
        ]
    )
    norm_squared = a**2 + b**2 + c**2 + d**2
    assert sp.simplify(left.T * left - norm_squared * sp.eye(4)) == sp.zeros(4)

    # Hopf coordinates on the full-measure chart 0<x<1:
    # alpha=sqrt(1-x)e^(i theta_2), beta=sqrt(x)e^(i theta_1).
    x, theta_1, theta_2 = sp.symbols("x theta_1 theta_2", positive=True, real=True)
    embedding = sp.Matrix(
        [
            sp.sqrt(1 - x) * sp.cos(theta_2),
            sp.sqrt(1 - x) * sp.sin(theta_2),
            sp.sqrt(x) * sp.cos(theta_1),
            sp.sqrt(x) * sp.sin(theta_1),
        ]
    )
    assert sp.trigsimp(embedding.dot(embedding)) == 1
    derivative = embedding.jacobian((x, theta_1, theta_2))
    gram = sp.simplify(sp.trigsimp(derivative.T * derivative))
    expected_gram = sp.diag(1 / (4 * x * (1 - x)), x, 1 - x)
    assert sp.simplify(gram - expected_gram) == sp.zeros(3)
    assert sp.simplify(gram.det()) == sp.Rational(1, 4)

    # Surface density is (1/2) dx dtheta_1 dtheta_2.  Dividing by
    # area(S^3)=2*pi^2 makes x uniform and both angles normalized Haar.
    normalized_density = sp.simplify(sp.Rational(1, 2) / (2 * sp.pi**2))
    product_density = sp.simplify(1 / (2 * sp.pi) ** 2)
    assert normalized_density == product_density

    print("PASS Long SU(2): unit-quaternion left action is orthogonal")
    print("PASS Long SU(2): Hopf Gram determinant is exactly 1/4")
    print("PASS Long SU(2): normalized Haar measure is dx dtheta1 dtheta2/(2pi)^2")


if __name__ == "__main__":
    main()
