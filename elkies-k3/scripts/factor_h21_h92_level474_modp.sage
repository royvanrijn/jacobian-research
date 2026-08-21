#!/usr/bin/env sage -python
"""Extract a normalized modular image of the H21/H92 level-474 factor.

This is the modular-image worker for characteristic-zero reconstruction.  It
reuses the pinned parsing and sparse-pullback implementation from
``verify_h21_h92_level474_branch.sage``, factors at a caller-selected prime,
selects the unique degree-21 factor, and verifies its incidence with the
sixth CM-24 tangent before exporting portable JSON coefficients.
"""

from sage.all import GF, PolynomialRing, factor

import argparse
import hashlib
import json
from pathlib import Path
import runpy
import time


NORMALIZING_EXPONENT = (13, 8)


def stage(name, **values):
    payload = "|".join(f"{key}={value}" for key, value in values.items())
    print(f"H21H92MODP|stage={name}" + (f"|{payload}" if payload else ""), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime", required=True, type=int)
    parser.add_argument("--h21", required=True, type=Path)
    parser.add_argument("--h92", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()

    prime = arguments.prime
    if prime in (2, 3, 37) or not GF(prime).is_field():
        raise ValueError("prime must be a prime other than 2, 3, or 37")

    verifier_path = Path(__file__).with_name("verify_h21_h92_level474_branch.sage")
    namespace = runpy.run_path(str(verifier_path))
    globals_dictionary = namespace["build_pullback"].__globals__
    globals_dictionary["PRIME"] = prime

    h21_hash = namespace["verify_input"](
        arguments.h21, namespace["H21_SHA256"], "H21"
    )
    h92_hash = namespace["verify_input"](
        arguments.h92, namespace["H92_SHA256"], "H92"
    )
    h21, unused_content = namespace["extract_h21"](arguments.h21)
    unused_ring, h92 = namespace["extract_h92"](arguments.h92)

    pullback = namespace["build_pullback"](h21, h92)
    exponents = tuple(tuple(exponent) for exponent in pullback.dict())
    valuations = tuple(min(exponent[index] for exponent in exponents) for index in range(2))
    r, s = pullback.parent().gens()
    residual = pullback // (r**valuations[0] * s**valuations[1])
    stage(
        "pullback",
        prime=prime,
        base_valuations=valuations,
        residual_degree=residual.total_degree(),
        residual_terms=len(residual.dict()),
    )

    start = time.monotonic()
    factorization = residual.factor()
    degrees = tuple(
        (int(component.total_degree()), len(component.dict()), int(multiplicity))
        for component, multiplicity in factorization
    )
    stage("factored", seconds=f"{time.monotonic() - start:.2f}", degrees=degrees)

    degree_21 = [
        component for component, multiplicity in factorization
        if component.total_degree() == 21 and multiplicity == 1
    ]
    if len(degree_21) != 1:
        raise AssertionError(f"expected one squarefree degree-21 factor, found {len(degree_21)}")
    target = degree_21[0]

    base = GF(prime)
    univariate = PolynomialRing(base, "Z")
    Z = univariate.gen()
    cm_polynomial = Z**2 - Z + 1
    if cm_polynomial.is_irreducible():
        extension = GF(prime**2, "w", modulus=cm_polynomial)
        cm_roots = (extension.gen(),)
    else:
        extension = base
        cm_roots = tuple(root for root, multiplicity in cm_polynomial.roots())
    local_ring = PolynomialRing(extension, names=("x", "y"))
    x, y = local_ring.gens()
    slope_ring = PolynomialRing(extension, "z")
    z = slope_ring.gen()
    target_extension = target.change_ring(extension)
    incidences = []
    for w in cm_roots:
        if target_extension(w, -w):
            continue
        translated = local_ring(target_extension(w + x, -w + y))
        order = min(sum(exponent) for exponent in translated.dict())
        initial = local_ring(
            {
                tuple(exponent): coefficient
                for exponent, coefficient in translated.dict().items()
                if sum(exponent) == order
            }
        )
        tangent = slope_ring(initial(z, 1))
        branch_6_slope = (-7 * w - 33) / 37
        if tangent(branch_6_slope) == 0:
            incidences.append(
                {
                    "cm_root": str(w),
                    "intersection_multiplicity": int(order),
                    "slope": str(branch_6_slope),
                    "tangent_polynomial": str(factor(tangent)),
                }
            )
    if not incidences:
        raise AssertionError("degree-21 factor does not contain the sixth CM-24 tangent")

    coefficient = target[NORMALIZING_EXPONENT]
    if not coefficient:
        raise AssertionError(f"normalizing coefficient {NORMALIZING_EXPONENT} vanishes")
    normalized = target / coefficient
    coefficient_records = [
        {
            "r": int(exponent[0]),
            "s": int(exponent[1]),
            "coefficient": int(normalized[exponent]),
        }
        for exponent in sorted(
            (tuple(exponent) for exponent in normalized.dict()), reverse=True
        )
    ]

    output = {
        "schema": "elkies-k3.h21-h92-level474-modular-factor.v1",
        "status": "PASS_MODULAR_FACTOR",
        "prime": prime,
        "inputs": {
            "h21": {"path": str(arguments.h21), "sha256": h21_hash},
            "h92": {"path": str(arguments.h92), "sha256": h92_hash},
        },
        "pullback": {
            "degree": int(pullback.total_degree()),
            "terms": len(pullback.dict()),
            "base_valuations": list(valuations),
            "residual_degree": int(residual.total_degree()),
            "residual_terms": len(residual.dict()),
            "factor_degrees_terms_multiplicities": [list(record) for record in degrees],
        },
        "target": {
            "degree": int(target.total_degree()),
            "terms": len(target.dict()),
            "normalizing_exponent": list(NORMALIZING_EXPONENT),
            "normalizing_coefficient_before_scaling": int(coefficient),
            "incidences": incidences,
            "coefficients": coefficient_records,
        },
        "proof_boundary": (
            "One modular image for CRT/rational reconstruction; no "
            "characteristic-zero factor is asserted."
        ),
    }
    encoded = json.dumps(output, indent=2, sort_keys=True) + "\n"
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded)
    stage(
        "complete",
        status=output["status"],
        output=arguments.output,
        sha256=hashlib.sha256(encoded.encode()).hexdigest(),
    )


if __name__ == "__main__":
    main()
