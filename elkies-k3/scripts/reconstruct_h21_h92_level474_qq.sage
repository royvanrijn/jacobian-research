#!/usr/bin/env sage -python
"""Reconstruct and exactly certify the H21/H92 level-474 factor over QQ.

The input files are modular degree-21 factors emitted by
``factor_h21_h92_level474_modp.sage``.  Their normalized coefficients are
combined by CRT and centered.  The resulting integer polynomial is then
certified as an exact characteristic-zero factor without expanding the full
degree-700 pullback over QQ.

The exact certificate uses the change r=x/(s-1).  It makes the candidate
monic of x-degree 13.  With weights wt(x)=3 and wt(s)=2, the transformed
factor has weight at most 39 and the transformed H21 pullback has weight at
most 1920.  Polynomial reduction therefore produces a remainder of x-degree
at most 12 whose coefficients have s-degree at most 960.  Exact vanishing at
961 distinct rational s-values proves that remainder is identically zero.
"""

from sage.all import Integer, PolynomialRing, QQ, ZZ, prod
from sage.arith.misc import CRT

import argparse
import hashlib
import json
from pathlib import Path
import runpy
import time


SCHEMA = "elkies-k3.h21-h92-level474-qq-factor.v1"
MODULAR_SCHEMA = "elkies-k3.h21-h92-level474-modular-factor.v1"
HORNER_ORDER = (2, 1, 3, 0)
EXPECTED_TERMS = 133
EXPECTED_MAX_ABS_COEFFICIENT = 26378808832


def stage(name, **values):
    payload = "|".join(f"{key}={value}" for key, value in values.items())
    print(f"H21H92QQ|stage={name}" + (f"|{payload}" if payload else ""), flush=True)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficient_map(record):
    return {
        (int(item["r"]), int(item["s"])): ZZ(item["coefficient"])
        for item in record["target"]["coefficients"]
    }


def compile_horner(terms, order, level=0):
    variable = order[level]
    if level == len(order) - 1:
        return (
            "leaf",
            variable,
            tuple((QQ(coefficient), tuple(exponent)[variable]) for exponent, coefficient in terms.items()),
        )
    groups = {}
    for exponent, coefficient in terms.items():
        exponent = tuple(exponent)
        groups.setdefault(exponent[variable], {})[exponent] = coefficient
    return (
        "node",
        variable,
        tuple(
            (exponent, compile_horner(groups[exponent], order, level + 1))
            for exponent in sorted(groups, reverse=True)
        ),
    )


def evaluate_horner(circuit, power_tables, modulus, polynomial_ring):
    kind, variable, payload = circuit
    if kind == "leaf":
        return sum(
            (coefficient * power_tables[variable][exponent] for coefficient, exponent in payload),
            polynomial_ring.zero(),
        )
    first_exponent, first_child = payload[0]
    previous = first_exponent
    result = evaluate_horner(first_child, power_tables, modulus, polynomial_ring)
    for exponent, child in payload[1:]:
        result = (
            result * power_tables[variable][previous - exponent]
            + evaluate_horner(child, power_tables, modulus, polynomial_ring)
        ).mod(modulus)
        previous = exponent
    if previous:
        result = (result * power_tables[variable][previous]).mod(modulus)
    return result


def polynomial_weight(polynomial, x_weight=3, s_weight=2):
    maximum = -1
    for x_exponent, coefficient in enumerate(polynomial.list()):
        for s_exponent, scalar in coefficient.dict().items():
            if not scalar:
                continue
            if not isinstance(s_exponent, (int, Integer)):
                s_exponent = s_exponent[0]
            maximum = max(maximum, x_weight * x_exponent + s_weight * s_exponent)
    return maximum


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h21", required=True, type=Path)
    parser.add_argument("--h92", required=True, type=Path)
    parser.add_argument(
        "--modular-factor",
        required=True,
        action="append",
        type=Path,
        help="repeat for each modular JSON image",
    )
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()

    verifier_path = Path(__file__).with_name("verify_h21_h92_level474_branch.sage")
    namespace = runpy.run_path(str(verifier_path))
    h21_hash = namespace["verify_input"](
        arguments.h21, namespace["H21_SHA256"], "H21"
    )
    h92_hash = namespace["verify_input"](
        arguments.h92, namespace["H92_SHA256"], "H92"
    )

    modular_records = []
    primes = []
    maps = []
    for path in arguments.modular_factor:
        record = json.loads(path.read_text())
        if record.get("schema") != MODULAR_SCHEMA or record.get("status") != "PASS_MODULAR_FACTOR":
            raise ValueError(f"not a passing modular-factor artifact: {path}")
        prime = int(record["prime"])
        if prime in primes:
            raise ValueError(f"duplicate prime {prime}")
        if record["inputs"]["h21"]["sha256"] != h21_hash:
            raise ValueError(f"H21 hash mismatch in {path}")
        if record["inputs"]["h92"]["sha256"] != h92_hash:
            raise ValueError(f"H92 hash mismatch in {path}")
        if int(record["target"]["degree"]) != 21 or not record["target"]["incidences"]:
            raise ValueError(f"artifact is not a branch-6 degree-21 factor: {path}")
        primes.append(prime)
        maps.append(coefficient_map(record))
        modular_records.append(
            {"path": str(path), "sha256": sha256(path), "prime": prime}
        )
    if len(primes) < 2:
        raise ValueError("at least two modular images are required")

    support = set().union(*(set(coefficient_map0) for coefficient_map0 in maps))
    modulus = ZZ(prod(primes))
    coefficients = {}
    for exponent in support:
        residue = ZZ(CRT([mapping.get(exponent, 0) for mapping in maps], primes))
        coefficients[exponent] = residue if residue <= modulus // 2 else residue - modulus
    coefficients = {exponent: value for exponent, value in coefficients.items() if value}
    for prime, mapping in zip(primes, maps):
        for exponent in support:
            if coefficients.get(exponent, 0) % prime != mapping.get(exponent, 0) % prime:
                raise AssertionError(f"CRT reduction mismatch at p={prime}, exponent={exponent}")
    maximum_coefficient = max(abs(value) for value in coefficients.values())
    stage(
        "crt",
        primes=",".join(map(str, primes)),
        modulus=modulus,
        terms=len(coefficients),
        max_abs=maximum_coefficient,
    )
    if len(coefficients) != EXPECTED_TERMS:
        raise AssertionError(f"expected {EXPECTED_TERMS} terms, got {len(coefficients)}")
    if maximum_coefficient != EXPECTED_MAX_ABS_COEFFICIENT:
        raise AssertionError(
            f"expected maximum coefficient {EXPECTED_MAX_ABS_COEFFICIENT}, got {maximum_coefficient}"
        )
    if coefficients.get((13, 8)) != 1:
        raise AssertionError("candidate normalization coefficient is not one")

    integer_ring = PolynomialRing(ZZ, names=("r", "s"))
    candidate = integer_ring(coefficients)
    if candidate.content() != 1 or candidate.total_degree() != 21:
        raise AssertionError("candidate is not primitive of total degree 21")

    h21_integer, unused_content = namespace["extract_h21"](arguments.h21)
    h92_ring, h92_coefficients = namespace["extract_h92"](arguments.h92)
    A1, A, B1, B, B2 = h92_coefficients
    satake = (
        -36 * A,
        162 * B,
        -2430 * (A * B + 2 * A1 * B2),
        -2916 * A**3 + 4374 * B**2 + 17496 * B1 * B2,
    )
    satake_weights = (2, 3, 5, 6)

    coefficient_ring = PolynomialRing(QQ, "s")
    s = coefficient_ring.gen()
    difference = s - 1
    x_ring = PolynomialRing(coefficient_ring, "x")
    x = x_ring.gen()

    transformed_factor = x_ring.zero()
    row_divisibilities = []
    for r_exponent in range(candidate.degree(integer_ring.gen(0)) + 1):
        row = coefficient_ring(
            sum(
                QQ(value) * s**s_exponent
                for (candidate_r_exponent, s_exponent), value in coefficients.items()
                if candidate_r_exponent == r_exponent
            )
        )
        if r_exponent <= 5:
            transformed_row = row * difference ** (5 - r_exponent)
            removed = 0
        else:
            divisor = difference ** (r_exponent - 5)
            transformed_row, remainder = row.quo_rem(divisor)
            if remainder:
                raise AssertionError(
                    f"r^{r_exponent} coefficient lacks the required (s-1) divisibility"
                )
            removed = r_exponent - 5
        row_divisibilities.append(removed)
        transformed_factor += transformed_row * x**r_exponent
    if transformed_factor.degree() != 13 or not transformed_factor.is_monic():
        raise AssertionError("transformed factor is not monic of x-degree 13")
    factor_weight = polynomial_weight(transformed_factor)
    if factor_weight > 39:
        raise AssertionError(f"transformed factor has weight {factor_weight}, not at most 39")

    def transform_satake(polynomial, weight):
        if polynomial.degree(h92_ring.gen(0)) > 4 * weight:
            raise AssertionError("Satake r-degree exceeds its denominator-clearing exponent")
        return x_ring(
            sum(
                QQ(coefficient)
                * s**exponent[1]
                * difference ** (4 * weight - exponent[0])
                * x**exponent[0]
                for exponent, coefficient in polynomial.dict().items()
            )
        )

    transformed_satake = tuple(
        transform_satake(polynomial, weight)
        for polynomial, weight in zip(satake, satake_weights)
    )
    satake_transformed_weights = tuple(polynomial_weight(value) for value in transformed_satake)
    if satake_transformed_weights != tuple(16 * weight for weight in satake_weights):
        raise AssertionError(
            f"unexpected transformed Satake weights {satake_transformed_weights}"
        )

    satake_ring = PolynomialRing(QQ, names=("s2", "s3", "s5", "s6"))
    h21 = satake_ring(h21_integer)
    h21_weights = set(
        sum(weight * exponent for weight, exponent in zip(satake_weights, tuple(monomial)))
        for monomial in h21.dict()
    )
    if h21_weights != {120}:
        raise AssertionError(f"H21 is not weight-120 homogeneous: {h21_weights}")
    pullback_weight = max(
        sum(
            exponent * transformed_weight
            for exponent, transformed_weight in zip(tuple(monomial), satake_transformed_weights)
        )
        for monomial in h21.dict()
    )
    if pullback_weight != 1920:
        raise AssertionError(f"unexpected transformed pullback weight {pullback_weight}")
    remainder_degree_bound = max(
        (pullback_weight - 3 * x_exponent) // 2 for x_exponent in range(13)
    )
    if remainder_degree_bound != 960:
        raise AssertionError(f"unexpected remainder degree bound {remainder_degree_bound}")

    transformed_factor_coefficients = tuple(transformed_factor.list())
    transformed_satake_coefficients = tuple(tuple(value.list()) for value in transformed_satake)

    def specialize(coefficient_list, value, target_ring):
        return target_ring([QQ(coefficient(value)) for coefficient in coefficient_list])

    maxima = tuple(
        max(tuple(monomial)[index] for monomial in h21.dict()) for index in range(4)
    )
    circuit = compile_horner(h21.dict(), HORNER_ORDER)
    points = tuple(range(-remainder_degree_bound // 2, remainder_degree_bound // 2 + 1))
    if len(points) != remainder_degree_bound + 1 or len(set(points)) != len(points):
        raise AssertionError("specialization point set has the wrong cardinality")

    exact_start = time.monotonic()
    for point_index, point in enumerate(points, 1):
        univariate_ring = PolynomialRing(QQ, "x")
        factor_at_point = specialize(
            transformed_factor_coefficients, point, univariate_ring
        )
        if factor_at_point.degree() != 13 or not factor_at_point.is_monic():
            raise AssertionError(f"factor specialization is not monic at s={point}")
        values_at_point = tuple(
            specialize(coefficients0, point, univariate_ring).mod(factor_at_point)
            for coefficients0 in transformed_satake_coefficients
        )
        power_tables = []
        for maximum, value in zip(maxima, values_at_point):
            powers = [univariate_ring.one()]
            for unused_exponent in range(maximum):
                powers.append((powers[-1] * value).mod(factor_at_point))
            power_tables.append(tuple(powers))
        image = evaluate_horner(
            circuit, tuple(power_tables), factor_at_point, univariate_ring
        )
        if image:
            raise AssertionError(f"nonzero exact remainder at s={point}: {image}")
        if point_index % 100 == 0 or point_index == len(points):
            stage(
                "exact_points",
                checked=point_index,
                total=len(points),
                last_s=point,
                seconds=f"{time.monotonic() - exact_start:.2f}",
            )

    coefficient_records = [
        {"r": exponent[0], "s": exponent[1], "coefficient": int(value)}
        for exponent, value in sorted(coefficients.items(), reverse=True)
    ]
    output = {
        "schema": SCHEMA,
        "status": "PASS_CHARACTERISTIC_ZERO_FACTOR",
        "inputs": {
            "h21": {"path": str(arguments.h21), "sha256": h21_hash},
            "h92": {"path": str(arguments.h92), "sha256": h92_hash},
            "modular_factors": modular_records,
        },
        "crt": {
            "primes": primes,
            "modulus": int(modulus),
            "centered_lift": True,
        },
        "factor": {
            "variables": ["r", "s"],
            "total_degree": int(candidate.total_degree()),
            "r_degree": int(candidate.degree(integer_ring.gen(0))),
            "s_degree": int(candidate.degree(integer_ring.gen(1))),
            "terms": len(coefficients),
            "content": int(candidate.content()),
            "maximum_absolute_coefficient": int(maximum_coefficient),
            "normalizing_exponent": [13, 8],
            "coefficients": coefficient_records,
        },
        "exact_divisibility_certificate": {
            "coordinate_change": "r=x/(s-1)",
            "transformed_factor_multiplier": "(s-1)^5",
            "transformed_factor_x_degree": int(transformed_factor.degree()),
            "transformed_factor_weight_bound": int(factor_weight),
            "weights": {"x": 3, "s": 2},
            "row_divisibilities_removed": row_divisibilities,
            "satake_denominator_exponents": [4 * weight for weight in satake_weights],
            "satake_transformed_weights": list(satake_transformed_weights),
            "h21_weight": 120,
            "pullback_weight_bound": int(pullback_weight),
            "remainder_x_degree_bound": 12,
            "remainder_coefficient_s_degree_bound": int(remainder_degree_bound),
            "exact_rational_specializations": {
                "first": points[0],
                "last": points[-1],
                "count": len(points),
                "all_zero": True,
            },
        },
        "conclusion": (
            "The primitive degree-21 polynomial is an exact factor over QQ "
            "of the H21 equation pulled back to the H92 chart, and its modular "
            "images are the branch-6 CM24 components."
        ),
        "proof_boundary": (
            "This proves the characteristic-zero component equation.  It does "
            "not yet give its normalization or a birational map to the "
            "published level-474 sextic."
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
