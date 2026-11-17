#!/usr/bin/env python3
"""Exact rational controls for the deformed-first-image descent obstruction.

The universal finite centralizer and h-pole arguments are written proofs;
this checker verifies their explicit position, vector-field and Long inputs.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from math import comb
from pathlib import Path
import platform
import tempfile

import sympy as s

HERE = Path(__file__).resolve().parent
INPUT = HERE.parent / "dixmier-audit" / "verify_long_pbw.py"
X, Q, P, Z = s.symbols("X Q P Z")
u, w, r, piu, pir = s.symbols("u w r pi_u pi_r")
h, epsilon = s.symbols("h epsilon", nonzero=True)


def zero(expression):
    assert s.cancel(expression) == 0, s.factor(expression)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    wx = 1 - 3 * X * Q
    ux = wx**2 / Q
    rx = X * (2 - 3 * X * Q)
    inverse = {X: u * (1 - w) / (3 * w**2), Q: w**2 / u}
    zero(wx.subs(inverse, simultaneous=True) - w)
    zero(ux.subs(inverse, simultaneous=True) - u)
    zero(rx.subs(inverse, simultaneous=True) - u * (1 - w**2) / (3 * w**2))
    zero(Q - 1 / (3 * rx + ux))
    zero(wx**2 - ux / (3 * rx + ux))

    v0 = ((1 + 3 * X * Q) / 2, -3 * Q**2)
    vu = (-3 * X**2 * Q**2 / (2 * wx), -Q**2)

    def derivative(vector, function):
        return s.cancel(vector[0] * s.diff(function, X)
                        + vector[1] * s.diff(function, Q))

    zero(derivative(v0, rx) - 1)
    zero(derivative(v0, ux))
    zero(derivative(v0, wx) + 3 * wx**3 / (2 * ux))
    zero(derivative(vu, rx))
    zero(derivative(vu, ux) - 1)
    for i in range(2):
        zero(derivative(v0, vu[i]) - derivative(vu, v0[i]))
    determinant = s.cancel(v0[0] * vu[1] - v0[1] * vu[0])
    zero(determinant + Q**2 / (2 * wx))
    assert determinant != 0

    d0_symbol = v0[0] * P + v0[1] * Z
    piu_symbol = vu[0] * P + vu[1] * Z
    momentum_inverse = {
        P: 2 * wx * (pir - 3 * piu),
        Z: -3 * X**2 * pir - wx * (1 + 3 * X * Q) * piu / Q**2,
    }
    zero(d0_symbol.subs(momentum_inverse, simultaneous=True) - pir)
    zero(piu_symbol.subs(momentum_inverse, simultaneous=True) - piu)

    specification = importlib.util.spec_from_file_location("canonical_long_input", INPUT)
    assert specification is not None and specification.loader is not None
    producer = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(producer)
    raw = producer.long_images()  # Polynomial INPUT only.

    def convert(polynomial):
        return s.expand(sum(s.Rational(c.numerator, c.denominator)
                            * X**e[0] * Q**e[1] * P**e[2] * Z**e[3]
                            for e, c in polynomial.items()))

    classical_R, _, _, classical_S = map(convert, raw)
    zero(classical_R - rx)
    B = 3 * X**2 * P + 2 * wx * Z
    beta = B - 9 * Q**2
    y = Q - X * beta / 3
    compact_S = y + 3 * X * (1 + X * y)**2 * beta + 3 * X * y**2 * (4 + 3 * X * y)
    zero(compact_S - classical_S)
    transformed_S = s.Poly(s.expand(classical_S.subs(momentum_inverse, simultaneous=True)), pir, piu)
    assert transformed_S.degree(pir) == 0
    assert transformed_S.degree(piu) == 2
    quadratic = s.factor(transformed_S.coeff_monomial(piu**2))
    expected_quadratic = -4 * rx * ux * X**2 / (3 * Q**3)
    zero(quadratic - expected_quadratic)
    zero(B.subs(momentum_inverse, simultaneous=True) + 2 * wx * piu / Q**2)

    deck_X = 2 / (3 * Q) - X
    for invariant in (rx, ux, Q):
        zero(invariant.subs(X, deck_X) - invariant)
    zero(wx.subs(X, deck_X) + wx)
    difference = s.factor(quadratic.subs(X, deck_X) - quadratic)
    expected_difference = -16 * rx * ux * wx / (27 * Q**5)
    zero(difference - expected_difference)
    assert difference != 0
    assert s.cancel(difference + expected_difference) != 0  # Wrong sign rejected.

    # Positive control: the undeformed exact-r second image really is finite
    # and has the desired classical reduction. Its exact r-commutation follows
    # from [B,r]=0 and its polynomial expression in X,Q,B.
    vb = (3 * X**2, 2 * wx)
    zero(derivative(vb, rx))
    A0 = Q + 3 * X * Q**2 - X * B / 3
    Sb0 = -3 * rx * A0**2 + (18 * X * Q - 8) * A0 + 9 * Q
    zero(Sb0 - classical_S)

    # One central epsilon symbol verifies identities for every substitution
    # epsilon(h) in h*K[[h]]; no truncation/order bound on that series is used.
    zero(derivative(v0, derivative(v0, rx)))

    def ore_add(*operators):
        result = {}
        for operator in operators:
            for degree, coefficient in operator.items():
                result[degree] = result.get(degree, 0) + coefficient
        return {degree: s.factor(coefficient) for degree, coefficient in result.items()
                if s.factor(coefficient) != 0}

    def ore_scale(operator, scalar):
        return {degree: s.factor(scalar * coefficient)
                for degree, coefficient in operator.items()
                if s.factor(scalar * coefficient) != 0}

    def ore_mul(left, right):
        result = {}
        for i, ai in left.items():
            for j, bj in right.items():
                for k in range(i + 1):
                    degree = i + j - k
                    result[degree] = result.get(degree, 0) + ai * comb(i, k) * h**k * s.diff(bj, r, k)
        return ore_add(result)

    def ore_comm(left, right):
        return ore_add(ore_mul(left, right), ore_scale(ore_mul(right, left), -1))

    d_op = {1: s.Integer(1)}
    r_op = {0: r}
    deformed = {0: r, 1: epsilon}
    generator = {2: epsilon / (2 * h)}
    assert ore_comm(d_op, deformed) == {0: h}
    first_adjoint = ore_comm(generator, r_op)
    assert first_adjoint == {1: epsilon}
    assert ore_comm(generator, first_adjoint) == {}
    assert ore_comm(generator, d_op) == {}

    files = [
        "DEFORMED_FIRST_IMAGE_AND_QUADRATIC_DESCENT_OBSTRUCTION.md",
        "check_deformed_first_image.py",
        "GLOBAL_ADMISSION_AND_QUARTIC_CRITICAL_FACE.md",
    ]
    hashes = {name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
              for name in files}
    hashes["../dixmier-audit/verify_long_pbw.py"] = hashlib.sha256(INPUT.read_bytes()).hexdigest()
    receipt = {
        "claim": "Every nonzero scalar-series deformation R_epsilon=r+epsilon(h)*D0, epsilon in h*K[[h]], fails Long quadratic descent",
        "status": "PASS",
        "environment": {"python": platform.python_version(), "sympy": s.__version__},
        "hypotheses": "Characteristic zero; arbitrary nonzero epsilon(h) in h*K[[h]], of any positive order; finite rational-position operators; h-regular prescribed canonical Long second residue.",
        "rational_chart_and_quadratic_deck": True,
        "commuting_derivations_and_inverse_momentum_frame": True,
        "D0_r_bracket_and_formal_exponential_input": True,
        "arbitrary_scalar_series_identities": {
            "first_image": "R_epsilon=r+epsilon(h)*D0",
            "epsilon_scope": "Every nonzero element of h*K[[h]]; positive order unbounded",
            "formal_generator": "(epsilon/h)*D0^2/2",
            "first_adjoint_r": "epsilon*D0",
            "second_adjoint_r_zero": True,
            "D0_fixed": True,
            "convergence": "Written proof: adjoint raises h valuation by ord_h(epsilon)>=1",
            "epsilon_zero_boundary": "Excluded hypothesis; the undeformed finite second image is a positive control"
        },
        "canonical_Long_input_and_compact_S_match": True,
        "quadratic_Piu_coefficient": "-4*r*u*X^2/(3*Q^3)",
        "deck_difference": "-16*r*u*w/(27*Q^5)",
        "nonzero_deck_difference_and_wrong_sign_rejected": True,
        "undeformed_finite_second_image_positive_control": True,
        "hashes": hashes,
        "coverage_boundary": "Exact rational/vector-field/Long coefficient identities only. The all-order finite centralizer, constants, parameter-pole noncancellation, formal convergence and exclusion are written proofs. No arbitrary deformed-first-image, matrix, algebraic-cover or division-ring exclusion.",
        "independent_replay": False,
        "formal_verification": False,
        "external_review": False,
    }
    output = args.output or Path(tempfile.mkdtemp(prefix="deformed-first-image-")) / "deformed_first_image_controls.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(receipt, indent=2) + "\n")
    print("PASS: rational frame, exact D0 first-pair input, Long quadratic coefficient and nonzero deck difference")
    print("Receipt:", output)


if __name__ == "__main__":
    main()
