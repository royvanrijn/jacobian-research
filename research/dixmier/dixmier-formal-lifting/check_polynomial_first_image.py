#!/usr/bin/env python3
"""Exact controls for the global polynomial-first-image conjugate criterion."""
from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
import tempfile

import sympy as s

HERE = Path(__file__).resolve().parent


def zero(expression):
    assert s.cancel(expression) == 0, s.factor(expression)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    X, Q, P, Z, r = s.symbols("X Q P Z r")
    h, alpha = s.symbols("h alpha", nonzero=True)
    beta = s.symbols("beta")
    R = X * (2 - 3 * X * Q)
    w = 1 - 3 * X * Q
    d0 = ((1 + 3 * X * Q) / 2, -3 * Q**2)
    tangent = (3 * X**2, 2 * w)

    def derivative(vector, function):
        return s.cancel(vector[0] * s.diff(function, X)
                        + vector[1] * s.diff(function, Q))

    matrix = s.Matrix([d0, tangent])
    zero(matrix.det() - 1)
    Dsymbol = d0[0] * P + d0[1] * Z
    Bsymbol = tangent[0] * P + tangent[1] * Z
    zero(2 * w * Dsymbol + 3 * Q**2 * Bsymbol - P)
    zero(-3 * X**2 * Dsymbol + (1 + 3 * X * Q) * Bsymbol / 2 - Z)
    zero(derivative(d0, R) - 1)
    zero(derivative(tangent, R))
    for i in range(2):
        zero(derivative(d0, tangent[i]) - derivative(tangent, d0[i])
             - s.Rational(9, 2) * Q * tangent[i])

    # Independent recursive Ore multiplication in [D0,r]=h. These finite
    # controls support, rather than prove, the written arbitrary-order formula.
    def left_D(operator):
        out = {}
        for order, coefficient in operator.items():
            out[order + 1] = out.get(order + 1, 0) + coefficient
            out[order] = out.get(order, 0) + h * s.diff(coefficient, r)
        return {j: s.expand(c) for j, c in out.items() if c != 0}

    controls = []
    for degree in range(1, 6):
        f = r**degree
        for order in range(1, 7):
            product = {0: f}
            for _ in range(order):
                product = left_D(product)
            product[order] = s.expand(product.get(order, 0) - f)
            commutator = {j: c for j, c in product.items() if c != 0}
            assert max(commutator) == order - 1
            zero(commutator[order - 1] - order * h * s.diff(f, r))
            for j in range(order):
                expected = s.binomial(order, order - j) * h**(order - j) * s.diff(f, r, order - j)
                zero(commutator.get(j, 0) - expected)
        controls.append({"degree": degree, "composition_degree": s.Poly(f.subs(r, R), X, Q).total_degree()})
        assert controls[-1]["composition_degree"] == 3 * degree

    affine = alpha * r + beta
    zero(h * s.diff(affine, r) / alpha - h)
    nonlinear = r + h * r**4
    derivative_f = s.diff(nonlinear, r)
    zero(h * derivative_f / derivative_f - h)
    assert s.Poly(derivative_f.subs(r, R), X, Q).total_degree() == 9
    # Its formal inverse is geometric; a finite prefix has the expected
    # nonzero remainder, so no polynomial-in-r inverse is being asserted.
    for n in range(1, 5):
        prefix = sum((-4 * h * r**3)**j for j in range(n))
        zero(derivative_f * prefix - 1 + (-4 * h * r**3)**n)

    files = [
        "POLYNOMIAL_FIRST_IMAGE_REPARAMETRIZATION_OBSTRUCTION.md",
        "check_polynomial_first_image.py",
        "EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md",
    ]
    receipt = {
        "claim": "A polynomial f(r) has a global finite conjugate exactly when f is affine",
        "status": "PASS",
        "environment": {"python": platform.python_version(), "sympy": s.__version__},
        "coefficient_scope": "Characteristic zero; L=K((h)); globally finite PBW source algebra; arbitrary finite partner order and scalar parameter poles.",
        "D0_B_coefficient_determinant": "1",
        "global_polynomial_inverse_frame": True,
        "D0_r_and_B_r_commutators": True,
        "D0_B_commutator": "(9/2)*h*Q*B",
        "recursive_Ore_controls": "orders1..6 and polynomial degrees1..5; not an all-order finite test",
        "composition_degree_controls": controls,
        "affine_positive_control": "D=alpha^-1*D0 for f=alpha*r+beta",
        "nonlinear_rational_localized_positive_control": "D=(1+4*h*r^3)^-1*D0 for f=r+h*r^4",
        "nonlinear_derivative_is_nonconstant_global_polynomial": True,
        "finite_geometric_prefixes_reject_polynomial_inverse": True,
        "hashes": {name: sha256((HERE / name).read_bytes()).hexdigest() for name in files},
        "coverage_boundary": "Exact frame, vector-field and finite Ore identities only. Global PBW equivalence, centralizer, units and arbitrary-order conjugate classification are written proofs. Rational-position conjugates are allowed, not excluded; no general odd-jet or arbitrary first-image theorem.",
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
    }
    output = args.output or Path(tempfile.mkdtemp(prefix="polynomial-first-image-")) / "polynomial_first_image_controls.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS: determinant-one global frame; exact vector-field relations; Ore and affine/localized controls")
    print("Receipt:", output)


if __name__ == "__main__":
    main()
