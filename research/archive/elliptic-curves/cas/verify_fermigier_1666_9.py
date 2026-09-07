#!/usr/bin/env python3
"""Verify the low-conductor Fermigier specialization ``T=1666/9``.

The exact Python part specializes the binary-quartic Jacobian, evaluates the
degree-20 discriminant factor, and checks twelve rational points.  PARI/GP is
then used for the global minimal model, conductor, certified factorization,
local reduction data, root number, and a numerical height-pairing determinant.

No ``ellrank`` call is made.  In particular, the nonzero numerical determinant
is recorded only as numerical evidence and is not an exact Mordell--Weil
independence or rank certificate.
"""

from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import json
from pathlib import Path
import platform
import shutil
import subprocess
from typing import Any, Sequence

from ek_k3 import rational_to_string, valuation
from fermigier_mestre import FermigierMestreFamily, ROOTS
from pari_bridge import minimal_curve_data, pari_version


PARAMETER = Fraction(1666, 9)
TARGET_LOG_CONDUCTOR = Decimal("182.72")

EXPECTED_SPECIALIZED_MODEL = (
    Fraction(0),
    Fraction(0),
    Fraction(0),
    Fraction(-4130694945279060548750364892489, 1594323),
    Fraction(
        16797885092800609712297229925634166953636434826,
        10460353203,
    ),
)
EXPECTED_BINARY_INVARIANTS = (
    Fraction(4130694945279060548750364892489, 43046721),
    Fraction(
        -16797885092800609712297229925634166953636434826,
        282429536481,
    ),
)
EXPECTED_DISCRIMINANT_FACTOR = Fraction(
    -39796602799841062365899369832967190529119013443184458671488468796760736707934400,
    5559060566555523,
)
EXPECTED_H_VALUATIONS = {7: 18, 17: 4, 37: 3}

EXPECTED_MINIMAL_MODEL = (
    1,
    0,
    1,
    -35841792875182741121324144,
    82627122352018241034203681114802552242,
)
EXPECTED_CONDUCTOR = (
    101523255017246417712694892860237179024105368632978033830
)
EXPECTED_CONDUCTOR_FACTORIZATION = (
    (2, 1),
    (3, 1),
    (5, 1),
    (7, 2),
    (13, 1),
    (17, 2),
    (37, 1),
    (43, 2),
    (177298878520409, 1),
    (1515522470494692043926333382781, 1),
)
BAD_PRIMES = tuple(prime for prime, _ in EXPECTED_CONDUCTOR_FACTORIZATION)
EXPECTED_MINIMAL_DISCRIMINANT = (
    -2578507259649968205788921034197730786223615384053775609005819111033566871300
)
EXPECTED_ROOT_NUMBER = 1
EXPECTED_LOCAL_REDUCTION = {
    2: (1, 6, 2, 0, 2, -1),
    3: (1, 19, 15, 0, 15, 1),
    5: (1, 6, 2, 0, 2, -1),
    7: (2, -1, 4, 3, 6, 0),
    13: (1, 6, 2, 0, 2, -1),
    17: (2, 4, 3, 2, 4, 0),
    37: (1, 7, 3, 0, 3, 1),
    43: (2, 3, 2, 1, 3, 0),
    177298878520409: (1, 5, 1, 0, 1, 1),
    1515522470494692043926333382781: (1, 5, 1, 0, 1, -1),
}

REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/verify_fermigier_1666_9.py"
)


def factorization_product(factorization: Sequence[tuple[int, int]]) -> int:
    """Return the exact product represented by ``(prime, exponent)`` pairs."""

    answer = 1
    for prime, exponent in factorization:
        if prime < 2 or exponent < 1:
            raise ValueError("factorization entries must have p >= 2 and e >= 1")
        answer *= prime**exponent
    return answer


def integral_weierstrass_invariants(
    coefficients: Sequence[int],
) -> dict[str, int]:
    """Compute ``c4``, ``c6``, and the discriminant of an integral model."""

    if len(coefficients) != 5:
        raise ValueError("an extended Weierstrass vector has five coefficients")
    a1, a2, a3, a4, a6 = coefficients
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    b8 = (
        a1**2 * a6
        + 4 * a2 * a6
        - a1 * a3 * a4
        + a2 * a3**2
        - a4**2
    )
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    discriminant = -b2**2 * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    return {"c4": c4, "c6": c6, "discriminant": discriminant}


def kodaira_symbol(code: int) -> str:
    """Translate PARI's integer Kodaira code to its conventional symbol."""

    exceptional = {
        1: "I0",
        2: "II",
        3: "III",
        4: "IV",
        -1: "IV*",
        -2: "III*",
        -3: "II*",
        -4: "I0*",
    }
    if code in exceptional:
        return exceptional[code]
    if code >= 5:
        return f"I{code - 4}"
    if code <= -5:
        return f"I{-code - 4}*"
    raise ValueError(f"unrecognized PARI Kodaira code {code}")


def selected_known_points() -> tuple[tuple[Fraction, Fraction], ...]:
    """Return exact family images 2--13, the scanner's twelve-point seed."""

    all_points = FermigierMestreFamily.known_jacobian_points(PARAMETER)
    if len(all_points) != 13:
        raise AssertionError("the family no longer returned thirteen known images")
    points = all_points[1:]
    if len(points) != 12 or len(set(points)) != 12:
        raise AssertionError("the selected Jacobian seed must contain 12 distinct points")
    return points


def exact_specialization_data() -> dict[str, Any]:
    """Run all dependency-free exact checks and return serializable data."""

    coefficients = FermigierMestreFamily.coefficients(PARAMETER)
    if coefficients != EXPECTED_SPECIALIZED_MODEL:
        raise AssertionError("the pinned specialized Jacobian model changed")
    binary_invariants = FermigierMestreFamily.binary_invariants(PARAMETER)
    if binary_invariants != EXPECTED_BINARY_INVARIANTS:
        raise AssertionError("the pinned binary-quartic invariants changed")

    discriminant_factor = FermigierMestreFamily.discriminant_factor(PARAMETER)
    if discriminant_factor != EXPECTED_DISCRIMINANT_FACTOR:
        raise AssertionError("the pinned value H(1666/9) changed")
    h_valuations = {
        prime: valuation(discriminant_factor, prime)
        for prime in EXPECTED_H_VALUATIONS
    }
    if h_valuations != EXPECTED_H_VALUATIONS:
        raise AssertionError("the pinned H-valuations changed")

    _, _, _, coefficient_a, coefficient_b = coefficients
    serialized_points = []
    for selection_index, (x_value, y_value) in enumerate(
        selected_known_points(), 1
    ):
        residual = y_value**2 - (
            x_value**3 + coefficient_a * x_value + coefficient_b
        )
        if residual != 0:
            raise AssertionError("a selected family point left the specialized curve")
        serialized_points.append(
            {
                "selection_index": selection_index,
                "family_known_point_index": selection_index + 1,
                "x": rational_to_string(x_value),
                "y": rational_to_string(y_value),
                "exact_curve_residual": "0",
                "on_specialized_model": True,
            }
        )

    return {
        "parameter": rational_to_string(PARAMETER),
        "binary_quartic_invariants": {
            "I": rational_to_string(binary_invariants[0]),
            "J": rational_to_string(binary_invariants[1]),
        },
        "specialized_short_weierstrass_model": {
            "equation": "y^2 = x^3 + a4*x + a6",
            "coefficients_a1_a2_a3_a4_a6": [
                rational_to_string(value) for value in coefficients
            ],
        },
        "discriminant_factor": {
            "name": "H(T) = disc_X(R_T)/16",
            "value": rational_to_string(discriminant_factor),
            "exact_valuations": {
                str(prime): exponent for prime, exponent in h_valuations.items()
            },
        },
        "known_jacobian_points": {
            "selection": (
                "family known images 2--13; image 1 is omitted to match the "
                "twelve-point search seed"
            ),
            "coordinate_model": "the specialized short Weierstrass model above",
            "count": len(serialized_points),
            "all_checked_exactly_on_curve": True,
            "points": serialized_points,
        },
    }


def pari_factorization(integer: int, *, timeout: float) -> list[dict[str, Any]]:
    """Factor an integer with PARI and certify each returned base with ``isprime``."""

    if integer < 1:
        raise ValueError("only positive integers are factored")
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    program = "\n".join(
        (
            f"N={integer};",
            "F=factor(N);",
            'print("FACTOR_BEGIN");',
            (
                "for(i=1,matsize(F)[1],"
                'print(F[i,1]," ",F[i,2]," ",isprime(F[i,1])));'
            ),
            'print("FACTOR_END");',
            "quit",
        )
    )
    result = subprocess.run(
        [executable, "-q"],
        input=program + "\n",
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI/GP factorization failed: {result.stderr.strip()}")
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    try:
        start = lines.index("FACTOR_BEGIN") + 1
        end = lines.index("FACTOR_END")
    except ValueError as error:
        raise RuntimeError("PARI/GP did not emit a complete factor block") from error
    answer = []
    for line in lines[start:end]:
        prime_text, exponent_text, isprime_text = line.split()
        answer.append(
            {
                "prime": int(prime_text),
                "exponent": int(exponent_text),
                "pari_isprime": bool(int(isprime_text)),
            }
        )
    if not answer:
        raise RuntimeError("PARI/GP emitted an empty factorization")
    return answer


def _normalize_local_reduction(curve: dict[str, Any]) -> dict[str, dict[str, Any]]:
    normalized: dict[str, dict[str, Any]] = {}
    field_order = (
        "conductor_exponent",
        "kodaira_code",
        "tamagawa_number",
        "minimal_c4_valuation",
        "minimal_discriminant_valuation",
        "ellap",
    )
    for prime in BAD_PRIMES:
        entry = dict(curve["local_reduction"][str(prime)])
        observed = tuple(entry[field] for field in field_order)
        if observed != EXPECTED_LOCAL_REDUCTION[prime]:
            raise AssertionError(f"the pinned local reduction at p={prime} changed")
        entry["kodaira_symbol"] = kodaira_symbol(entry["kodaira_code"])
        normalized[str(prime)] = entry
    return normalized


def build_artifact(*, timeout: float, stack_bytes: int) -> dict[str, Any]:
    """Run the exact Python and PARI checks and assemble the result artifact."""

    exact = exact_specialization_data()
    points = selected_known_points()
    curve = minimal_curve_data(
        EXPECTED_SPECIALIZED_MODEL,
        timeout=timeout,
        known_points=points,
        local_primes=BAD_PRIMES,
        stack_bytes=stack_bytes,
    )
    if tuple(curve["minimal_model"]) != EXPECTED_MINIMAL_MODEL:
        raise AssertionError("the pinned global minimal model changed")
    if curve["conductor"] != EXPECTED_CONDUCTOR:
        raise AssertionError("the pinned conductor changed")
    if curve["minimal_discriminant"] != EXPECTED_MINIMAL_DISCRIMINANT:
        raise AssertionError("the pinned minimal discriminant changed")
    if curve["root_number"] != EXPECTED_ROOT_NUMBER:
        raise AssertionError("the pinned global root number changed")
    if "pari_ellrank" in curve:
        raise AssertionError("this verifier must not run ellrank")
    if curve["supplied_points"] != {
        "count": 12,
        "on_curve_count": 12,
        "all_on_curve": True,
        "height_pairing_determinant_approx": curve["supplied_points"][
            "height_pairing_determinant_approx"
        ],
    }:
        raise AssertionError("PARI did not accept all twelve selected points")
    if Decimal(
        curve["supplied_points"]["height_pairing_determinant_approx"]
    ) <= 0:
        raise AssertionError("the numerical height-pairing determinant is not positive")

    factor_rows = pari_factorization(EXPECTED_CONDUCTOR, timeout=timeout)
    observed_factorization = tuple(
        (row["prime"], row["exponent"]) for row in factor_rows
    )
    if observed_factorization != EXPECTED_CONDUCTOR_FACTORIZATION:
        raise AssertionError("the pinned conductor factorization changed")
    if not all(row["pari_isprime"] for row in factor_rows):
        raise AssertionError("PARI did not certify every conductor factor as prime")
    if factorization_product(observed_factorization) != EXPECTED_CONDUCTOR:
        raise AssertionError("the conductor factors do not multiply back exactly")

    minimal_invariants = integral_weierstrass_invariants(EXPECTED_MINIMAL_MODEL)
    if minimal_invariants["discriminant"] != EXPECTED_MINIMAL_DISCRIMINANT:
        raise AssertionError("the exact Python minimal-model discriminant changed")
    local_reduction = _normalize_local_reduction(curve)
    factor_exponents = dict(EXPECTED_CONDUCTOR_FACTORIZATION)
    for prime in BAD_PRIMES:
        local = local_reduction[str(prime)]
        if local["conductor_exponent"] != factor_exponents[prime]:
            raise AssertionError(f"local/global conductor mismatch at p={prime}")
        if valuation(Fraction(minimal_invariants["c4"]), prime) != local[
            "minimal_c4_valuation"
        ]:
            raise AssertionError(f"minimal c4 valuation mismatch at p={prime}")
        if valuation(Fraction(EXPECTED_MINIMAL_DISCRIMINANT), prime) != local[
            "minimal_discriminant_valuation"
        ]:
            raise AssertionError(f"minimal discriminant valuation mismatch at p={prime}")

    with localcontext() as context:
        context.prec = 80
        log_conductor = Decimal(curve["log_conductor"])
        threshold_margin = TARGET_LOG_CONDUCTOR - log_conductor
    if threshold_margin <= 0:
        raise AssertionError("the candidate no longer lies below the log-conductor bound")

    curve_record = {
        "global_minimal_model": {
            "equation": (
                "y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6"
            ),
            "coefficients_a1_a2_a3_a4_a6": list(curve["minimal_model"]),
        },
        "minimal_model_invariants": {
            key: str(value) for key, value in minimal_invariants.items()
        },
        "conductor": str(curve["conductor"]),
        "conductor_factorization": factor_rows,
        "factorization_product_checked_exactly": True,
        "natural_log_conductor_approx": curve["log_conductor"],
        "strict_threshold": str(TARGET_LOG_CONDUCTOR),
        "threshold_minus_log_conductor_approx": str(threshold_margin),
        "below_strict_threshold": True,
        "minimal_discriminant": str(curve["minimal_discriminant"]),
        "global_root_number": curve["root_number"],
        "local_reduction_at_every_bad_prime": local_reduction,
    }

    numerical_pairing = dict(curve["supplied_points"])
    numerical_pairing["status"] = (
        "numerical evidence only; a nonzero floating determinant is not an "
        "exact independence certificate"
    )
    return {
        "schema_version": 1,
        "status": (
            "verified exact specialization and point equations plus PARI global/local "
            "arithmetic; rank remains unknown and no target hit is certified"
        ),
        "family": {
            "name": "normalized Fermigier--Mestre family",
            "root_tuple": list(ROOTS),
        },
        "exact_specialization": exact,
        "curve": curve_record,
        "numerical_height_pairing": numerical_pairing,
        "rank_status": {
            "status": "unknown",
            "ellrank_invoked": False,
            "exact_independence_certificate": None,
            "exact_rank_lower_bound_claimed_here": None,
            "explanation": (
                "Twelve exact points are checked on the curve, but this artifact "
                "contains no exact proof that they are independent."
            ),
        },
        "targets": {
            "small_conductor_target": {
                "rank_at_least": 21,
                "log_conductor_strict_upper_bound": str(TARGET_LOG_CONDUCTOR),
                "conductor_condition_met": True,
                "rank_condition_certified": False,
                "certified_hit": False,
            },
            "alternative_rank_target": {
                "rank_at_least": 30,
                "rank_condition_certified": False,
                "certified_hit": False,
            },
        },
        "method": {
            "exact_python": (
                "Fraction specialization, H-valuations, twelve curve equations, "
                "factor-product identity, and minimal-model invariant/valuation checks"
            ),
            "pari_gp": (
                "ellminimalmodel, ellglobalred, elllocalred at all bad primes, "
                "ellrootno, factor plus isprime, and ellheightmatrix"
            ),
            "excluded": "ellrank is intentionally not called",
        },
        "software": {"python": platform.python_version(), "pari_gp": pari_version()},
        "reproducing_command": REPRODUCING_COMMAND,
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_fermigier_1666_9.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.timeout <= 0:
        raise SystemExit("--timeout must be positive")
    if args.stack_bytes < 8_000_000:
        raise SystemExit("--stack-bytes must be at least 8000000")
    artifact = build_artifact(timeout=args.timeout, stack_bytes=args.stack_bytes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    curve = artifact["curve"]
    print(f"wrote {args.output}")
    print(f"conductor={curve['conductor']}")
    print(f"log(N)={curve['natural_log_conductor_approx']}")
    print(f"threshold margin={curve['threshold_minus_log_conductor_approx']}")
    print("rank status: unknown; no target hit is certified")


if __name__ == "__main__":
    main()
