#!/usr/bin/env python3
"""Verify that the two Case-1 identities compose to one unit certificate.

The large archived identity is written for four degree-five-field generators
obtained after eliminating ``u3`` on ``h != 0``. This verifier checks the
missing algebraic bridges exactly:

* every eliminated generator is already in the pre-division seven-generator
  ideal, with no extra power of ``h``;
* the weighted degree-five descent is an invertible scaling and field
  embedding;
* the specialized ``h=0`` identity lifts to ``1 = sum(A_i G_i) + B h``.

Together with the independently replayed identity ``h = sum(T_j F_j)``, these
checks prove a single identity ``1 = sum(U_i G_i)``. The combined multipliers
are deliberately kept as a factored straight-line certificate instead of
expanding the 89 MB hard multiplier list again.
"""

from __future__ import annotations

import argparse
import hashlib
import pickle
import sys
from pathlib import Path

from flint import fmpq, fmpq_poly


BRANCH1_SHA256 = "368a1dafdb6d26708b85d652a437c848a7676ba1419e4575ae74186f022621b9"
H0_SHA256 = "664de005e99bc6a0e61ba479ba64ed57ddbfa9c5399b955cbb12237ac70f8186"
DEG35_SHA256 = "082471d05a2a7ceebca9fd3a615d8fb6fddaee8ff80afae161a043f71edcb575"
DEG5_SHA256 = "5a6e423d74ef09fc9c7a7282c500bda566018d7e56a93124665796bbe417cedf"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def default_replay_root() -> Path:
    external = Path(__file__).resolve().parents[1] / "external" / "zenodo-21479814"
    matches = sorted(external.glob("*/release_bundle/exact_replay"))
    if len(matches) != 1:
        raise RuntimeError("pass the exact_replay directory explicitly")
    return matches[0]


def add(left, right, zero):
    out = dict(left)
    for monomial, coefficient in right.items():
        value = out.get(monomial, zero) + coefficient
        if value:
            out[monomial] = value
        else:
            out.pop(monomial, None)
    return out


def scale(poly, scalar):
    return {m: c * scalar for m, c in poly.items() if c * scalar}


def multiply(left, right, zero):
    out = {}
    for m, c in left.items():
        for n, d in right.items():
            monomial = tuple(a + b for a, b in zip(m, n))
            out[monomial] = out.get(monomial, zero) + c * d
    return {m: c for m, c in out.items() if c}


def multiply_variable(poly, variable):
    out = {}
    for monomial, coefficient in poly.items():
        target = list(monomial)
        target[variable] += 1
        out[tuple(target)] = coefficient
    return out


def evaluate(poly, value, modulus):
    result = fmpq_poly()
    for degree in range(len(poly) - 1, -1, -1):
        result = (result * value + fmpq_poly([poly[degree]])) % modulus
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("replay_root", nargs="?", type=Path)
    args = parser.parse_args()
    root = (args.replay_root or default_replay_root()).resolve()
    sys.path.insert(0, str(root))

    import exact_core as ec  # type: ignore
    from degree5_core import L, decode_poly as decode_l_poly  # type: ignore

    branch_path = root / "case1_branch1_after_w.pkl"
    h0_path = root / "h0_branch1_exact_certificate.pkl"
    deg35_path = root / "hne0_deg35.pkl"
    deg5_path = root / "hne0_polred.pkl"
    expected_hashes = {
        branch_path: BRANCH1_SHA256,
        h0_path: H0_SHA256,
        deg35_path: DEG35_SHA256,
        deg5_path: DEG5_SHA256,
    }
    for path, expected in expected_hashes.items():
        actual = sha256(path)
        if actual != expected:
            raise RuntimeError(f"unexpected SHA-256 for {path.name}: {actual}")

    def decode_k(serialized):
        coefficients = [fmpq(0) for _ in range(35)]
        for degree, (numerator, denominator) in serialized.items():
            coefficients[int(degree)] = fmpq(int(numerator), int(denominator))
        return ec.K(fmpq_poly(coefficients))

    def decode_k_system(path):
        raw = pickle.loads(path.read_bytes())
        return [
            {tuple(m): decode_k(c) for m, c in equation.items()}
            for equation in raw
        ]

    zero = ec.K(0)
    generators = decode_k_system(branch_path)
    eliminated = decode_k_system(deg35_path)
    if len(generators) != 7 or len(eliminated) != 6:
        raise RuntimeError("unexpected Case-1 generator count")

    # G0 = b_h (h*u3-N). The coefficient of u3 must be exactly b_h*h.
    g0 = generators[0]
    if max(m[3] for m in g0) != 1:
        raise RuntimeError("distinguished generator is not linear in u3")
    u3_coefficient = {m[:3]: c for m, c in g0.items() if m[3] == 1}
    if set(u3_coefficient) != {(1, 0, 0)}:
        raise RuntimeError("u3 coefficient is not a nonzero scalar times h")
    bh = u3_coefficient[(1, 0, 0)]

    # If Gj=a_j*u3+b_j, then Fj=h*Gj-(a_j/b_h)*G0.
    for index, (generator, expected) in enumerate(zip(generators[1:], eliminated), 1):
        if max(m[3] for m in generator) > 1:
            raise RuntimeError(f"generator {index} is nonlinear in u3")
        coefficient = {m[:3] + (0,): c for m, c in generator.items() if m[3] == 1}
        correction = multiply(coefficient, g0, zero)
        lifted = add(
            multiply_variable(generator, 0),
            scale(correction, -(bh.inv())),
            zero,
        )
        if any(m[3] for m in lifted):
            raise RuntimeError(f"u3 did not cancel in generator {index}")
        if {m[:3]: c for m, c in lifted.items()} != expected:
            raise RuntimeError(f"elimination lift failed for generator {index}")
    print("CASE1_ALL_SIX_ELIMINATION_LIFTS_PASS")

    # Recheck the weighted variable scaling and degree-five subfield embedding.
    weights = (5, 6, 5)
    row_weights = (6, 5, 1, 0, 6, 5)
    u = ec.K(fmpq_poly([0, 1]))
    minpoly5 = fmpq_poly([26, 0, 3, 3, -1, 1])
    phi = fmpq_poly(
        [
            fmpq(-9725570295901, 12623962),
            fmpq(-1170753213563, 971074),
            fmpq(-387111042229, 12623962),
            fmpq(1578225240619, 12623962),
            fmpq(-469713794365, 6311981),
        ]
    )
    degree5_factor = fmpq_poly([ec.MOD[7 * k] for k in range(6)])
    if not evaluate(degree5_factor, phi, minpoly5).is_zero():
        raise RuntimeError("degree-five field embedding failed")

    def descend(value):
        if any(value.p[i] and i % 7 for i in range(len(value.p))):
            raise RuntimeError("coefficient is outside the degree-five subfield")
        compressed = fmpq_poly(
            [value.p[7 * j] if 7 * j < len(value.p) else 0 for j in range(5)]
        )
        return evaluate(compressed, phi, minpoly5)

    transformed = [decode_l_poly(q) for q in pickle.loads(deg5_path.read_bytes())]
    for row, (old, new) in enumerate(zip(eliminated, transformed)):
        expected_terms = {}
        for monomial, coefficient in old.items():
            exponent = sum(a * b for a, b in zip(weights, monomial)) - row_weights[row]
            field_value = L(descend(coefficient * u**exponent))
            if field_value:
                expected_terms[monomial] = field_value
        if expected_terms != new.terms:
            raise RuntimeError(f"degree-five descent failed on row {row + 1}")
    print("CASE1_DEGREE_FIVE_DESCENT_ISOMORPHISM_PASS")

    # Lift 1=sum Ai*Gi|h=0 to 1=sum Ai*Gi+B*h.
    payload = pickle.loads(h0_path.read_bytes())
    if payload.get("format") != "jc2-h0-nullstellensatz-v1" or payload.get("branch") != "s=c":
        raise RuntimeError("unexpected h=0 certificate payload")
    multipliers3 = [
        {tuple(m): decode_k(c) for m, c in poly.items()}
        for poly in payload["multipliers"]
    ]
    multipliers4 = [
        {(0,) + m: c for m, c in multiplier.items()}
        for multiplier in multipliers3
    ]
    b_poly = {}
    lifted_sum = {}
    for generator, multiplier in zip(generators, multipliers4):
        specialized = {m: c for m, c in generator.items() if m[0] == 0}
        difference = add(generator, scale(specialized, ec.K(-1)), zero)
        if any(m[0] == 0 for m in difference):
            raise RuntimeError("Gi-Gi|h=0 is not divisible by h")
        quotient = {(m[0] - 1,) + m[1:]: c for m, c in difference.items()}
        b_poly = add(b_poly, scale(multiply(multiplier, quotient, zero), ec.K(-1)), zero)
        lifted_sum = add(lifted_sum, multiply(multiplier, generator, zero), zero)
    lifted_sum = add(lifted_sum, multiply_variable(b_poly, 0), zero)
    if lifted_sum != {(0, 0, 0, 0): ec.K(1)}:
        raise RuntimeError("h=0 identity did not lift to 1=sum(Ai*Gi)+B*h")
    print("CASE1_H0_IDENTITY_POLYNOMIAL_LIFT_PASS")
    print("CASE1_FACTORED_DIRECT_UNIT_COMPOSITION_PASS")


if __name__ == "__main__":
    main()
