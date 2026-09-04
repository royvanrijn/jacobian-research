#!/usr/bin/env sage-python
"""Exhaust the finite trace-zero genus-five R17 normalization chart.

The trace-zero ``P.O=0`` bisections are obtained from

    x(u)=x0+x1*u+...+x4*u^4.

For every one of the ``p^5`` polynomials in this subchart, this script factors
``x^3+A*x+B`` over ``GF(p)`` and retains it exactly when its polynomial
squareclass is a nonconstant squarefree quadratic. Cover keys keep both the
monic quadratic and the leading scalar character. This exhausts polynomial x
only. General trace-zero bisections can have x=X/Q, so it is not the complete
``P.O=0`` twist-section chart and is not a characteristic-zero rank result.
"""

# status: UNPROMOTED_RESULT
# claim: exhaustive finite-field census of the polynomial-x trace-zero subchart
# inputs: artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json
# outputs: artifacts/generated-results/elkies-k3-r17-074d9-tracezero-genus5-polynomial-x-p19-v1.json

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from sage.all import GF, PolynomialRing, QQ, ZZ


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-norm12-wgxli-lineage-fibres-v1.json"
)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if not ZZ(args.prime).is_prime() or args.prime < 5:
        parser.error("--prime must be a prime at least five")

    model_path = args.model.resolve()
    document = json.loads(model_path.read_text())
    if document.get("status") == "PROVED_EXACT_LINEAGE_REALIZATION_AND_DISPLAYED_QUOTIENTS":
        equation = document["representative"]
        if equation.get("chart") != "norm12-orbit-074d9":
            raise ValueError("the exact lineage representative changed")
        coordinate = "norm12-orbit-074d9 native u"
    elif document.get("status") == "PASS_TRANSCRIBED_PUBLISHED_R17_MODEL":
        equation = document
        coordinate = document["coordinate"]
    else:
        raise ValueError("the input is not a supported certified R17 model")

    prime = int(args.prime)
    field = GF(prime)
    ring = PolynomialRing(field, "u")
    u = ring.gen()

    def coefficient_mod(value):
        value = QQ(value)
        denominator = int(value.denominator()) % prime
        if denominator == 0:
            raise ZeroDivisionError(f"model coefficient denominator vanishes at p={prime}")
        return field(int(value.numerator()) % prime) / field(denominator)

    A = ring([coefficient_mod(value) for value in equation["A_coefficients_low_to_high"]])
    B = ring([coefficient_mod(value) for value in equation["B_coefficients_low_to_high"]])
    if A.degree() != 8 or B.degree() != 12:
        raise ValueError("the model degree drops at the selected prime")

    survivors = []
    buckets = defaultdict(list)
    degree_drop_count = 0
    repeated_quadratic_count = 0
    for coefficients in product(field, repeat=5):
        x = ring(coefficients)
        branch = x**3 + A * x + B
        if branch.degree() != 12:
            degree_drop_count += 1
            continue
        factorization = branch.factor()
        square_factor = ring.one()
        reduced = ring(factorization.unit())
        for factor, exponent in factorization:
            square_factor *= factor ** (int(exponent) // 2)
            if int(exponent) % 2:
                reduced *= factor
        if square_factor**2 * reduced != branch:
            raise ArithmeticError("polynomial squareclass reconstruction failed")
        if reduced.degree() != 2:
            continue
        q0, q1, q2 = (reduced[index] for index in range(3))
        if q1**2 - 4 * q0 * q2 == 0:
            repeated_quadratic_count += 1
            continue
        inverse = q2**-1
        scalar_character = int(q2 ** ((prime - 1) // 2))
        key = (int(q0 * inverse), int(q1 * inverse), scalar_character)
        index = len(survivors)
        survivors.append(
            {
                "x_coefficients_low_to_high": list(map(int, coefficients)),
                "removed_square_factor_coefficients_low_to_high": [
                    int(value) for value in square_factor
                ],
                "reduced_quadratic_coefficients_low_to_high": [
                    int(q0), int(q1), int(q2)
                ],
                "cover_key_monic_q0_q1_and_scalar_character": list(key),
            }
        )
        buckets[key].append(index)

    collisions = []
    for key, indices in sorted(buckets.items()):
        if len(indices) < 2:
            continue
        collisions.append(
            {
                "cover_key_monic_q0_q1_and_scalar_character": list(key),
                "multiplicity": len(indices),
                "survivor_indices": indices,
            }
        )
    payload = {
        "schema": "elkies-k3.r17-tracezero-genus5-normalization-modp-search.v1",
        "status": "PASS_COMPLETE_MODP_TRACEZERO_GENUS5_POLYNOMIAL_X_CENSUS",
        "prime": prime,
        "coordinate": coordinate,
        "search": {
            "x_degree_bound": 4,
            "total_x_polynomial_count": prime**5,
            "degree_drop_count": degree_drop_count,
            "repeated_quadratic_count": repeated_quadratic_count,
        },
        "survivor_count": len(survivors),
        "distinct_cover_count": len(buckets),
        "collision_cover_count": len(collisions),
        "triple_collision_cover_count": sum(
            record["multiplicity"] >= 3 for record in collisions
        ),
        "maximum_cover_multiplicity": max(map(len, buckets.values()), default=0),
        "survivors": survivors,
        "cover_collisions": collisions,
        "inputs": {relative(model_path): digest(model_path)},
        "proof_boundary": (
            "The p^5 polynomial-x trace-zero subchart is exhaustive at this prime "
            "and retains the full scalar squareclass. It omits the general x=X/Q "
            "trace-zero form, degree drops, repeated covers, denominator sections, "
            "parameter-at-infinity charts, and solutions with bad reduction. Modular "
            "cover collisions are necessary conditions only; characteristic-zero "
            "reconstruction and height-pairing independence are separate gates."
        ),
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "search_r17_tracezero_genus5_normalizations_modp.sage "
            f"--model {relative(model_path)} --prime {prime} "
            f"--output {relative(args.output if args.output.is_absolute() else ROOT / args.output)}"
        ),
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if args.check:
        if output.read_text() != encoded:
            raise SystemExit("stored artifact differs from replay")
        print(f"PASS check {relative(output)}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(
        "R17TRACEZEROGENUS5"
        f"|p={prime}|x={prime**5}|survivors={len(survivors)}"
        f"|covers={len(buckets)}|collisions={len(collisions)}"
        f"|triples={payload['triple_collision_cover_count']}"
        f"|max_multiplicity={payload['maximum_cover_multiplicity']}"
        f"|output={relative(output)}|status={payload['status']}"
    )


if __name__ == "__main__":
    main()
