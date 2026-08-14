#!/usr/bin/env python3
"""Exact finite-reduction certificate for Kihara's 14 sections at ``t=2``.

The binary-quartic covariant map sends the fifteen printed quartic points to
the short Jacobian.  Subtracting the image of ``P15`` gives fourteen exact
points.  These differences are divisible by two; PARI is used only to find
halves, which are then checked by the dependency-free exact group law.
One further index-two saturation is found and checked.  The final basis is
certified independent using reductions modulo good primes and a separate
modular certificate excluding rational 2-torsion.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
from typing import Sequence

from alternate_quartic_covers import point_on_short_curve, short_add
from kihara_rank14 import (
    kihara_specialization,
    known_quartic_points,
    short_jacobian_coefficients,
)
from mod2_reduction_independence import (
    Mod2ReductionSignature,
    combined_mod2_rank,
    find_mod2_reduction_certificate,
    find_two_torsion_certificate_prime,
)
from nagao_1994 import quartic_covariants_at


Q = Fraction
CERTIFICATE_PARAMETER = Q(2)
SATURATED_SUBSET_INDICES = tuple(range(12))
REPRODUCE_COMMAND = (
    ".venv/bin/python elliptic-curves/cas/certify_kihara_rank14_basis.py "
    "--output artifacts/generated-results/elliptic_kihara_rank14_basis.json"
)


def _q_string(value: Fraction) -> str:
    value = Q(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _point_record(point: tuple[Fraction, Fraction]) -> list[str]:
    return [_q_string(point[0]), _q_string(point[1])]


def _negate(point: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    return point[0], -point[1]


def covariant_images(
    parameter_t: Fraction = CERTIFICATE_PARAMETER,
) -> tuple[
    tuple[Fraction, ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
]:
    """Return the short model, fifteen images, and fourteen ``Pi-P15`` images."""

    specialization = kihara_specialization(parameter_t)
    coefficients = short_jacobian_coefficients(parameter_t)
    images: list[tuple[Fraction, Fraction]] = []
    for x_value, y_value in known_quartic_points(parameter_t):
        g_value, h_value = quartic_covariants_at(
            specialization.quartic_coefficients, x_value
        )
        image = (36 * g_value / y_value**2, 108 * h_value / y_value**3)
        if not point_on_short_curve(coefficients, image):
            raise AssertionError("a covariant image missed the short Jacobian")
        images.append(image)
    if len(set(images)) != 15:
        raise AssertionError("the covariant map collapsed two printed points")
    negative_origin = _negate(images[14])
    differences = tuple(
        short_add(coefficients, image, negative_origin) for image in images[:14]
    )
    if any(point is None for point in differences):
        raise AssertionError("a printed section specialized to the chosen origin")
    return coefficients, tuple(images), tuple(differences)  # type: ignore[arg-type]


def _terminate_process_group(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=2)
    except ProcessLookupError:
        return
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=2)


def _run_gp(program: str, *, timeout: float) -> str:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    if timeout <= 0 or timeout > 60:
        raise ValueError("GP timeout must lie in (0,60]")
    process = subprocess.Popen(
        [executable, "-q", "-s", "512000000"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(program, timeout=timeout)
    except BaseException:
        _terminate_process_group(process)
        raise
    if process.returncode != 0 or "***" in stderr:
        raise RuntimeError(f"PARI/GP failed: {stderr.strip()}")
    return stdout


def _pari_halves(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
    *,
    timeout: float,
) -> tuple[tuple[Fraction, Fraction], ...]:
    vector = ",".join(f"({_q_string(Q(value))})" for value in coefficients)
    point_text = ",".join(
        f"[({_q_string(x_value)}),({_q_string(y_value)})]"
        for x_value, y_value in points
    )
    program = "\n".join(
        (
            f"E=ellinit([{vector}]);",
            f"P=[{point_text}];",
            'print("HALVES_BEGIN");',
            (
                "for(i=1,#P,Qh=0;"
                'if(!ellisdivisible(E,P[i],2,&Qh),error("not divisible by 2"));'
                "print(Qh));"
            ),
            'print("HALVES_END");',
            "quit",
        )
    )
    lines = [line.strip() for line in _run_gp(program + "\n", timeout=timeout).splitlines()]
    start = lines.index("HALVES_BEGIN") + 1
    end = lines.index("HALVES_END")
    answer: list[tuple[Fraction, Fraction]] = []
    pattern = re.compile(r"\[([^,]+), ([^\]]+)\]")
    for line in lines[start:end]:
        match = pattern.fullmatch(line)
        if match is None:
            raise ValueError(f"could not parse PARI point {line!r}")
        answer.append((Q(match.group(1)), Q(match.group(2))))
    if len(answer) != len(points):
        raise AssertionError("PARI returned the wrong number of halves")
    for point, half in zip(points, answer):
        if not point_on_short_curve(coefficients, half):
            raise AssertionError("a PARI half is not on the exact curve")
        if short_add(coefficients, half, half) != point:
            raise AssertionError("a PARI half failed exact doubling")
    return tuple(answer)


def _point_sum(
    coefficients: Sequence[Fraction],
    points: Sequence[tuple[Fraction, Fraction]],
) -> tuple[Fraction, Fraction]:
    answer = None
    for point in points:
        answer = short_add(coefficients, answer, point)
    if answer is None:
        raise ValueError("the exact point sum was the identity")
    return answer


def derive_saturated_basis(
    *, gp_timeout: float = 20.0
) -> tuple[
    tuple[Fraction, ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
    tuple[tuple[Fraction, Fraction], ...],
]:
    coefficients, images, differences = covariant_images()
    first_halves = _pari_halves(
        coefficients, differences, timeout=gp_timeout
    )
    subset_sum = _point_sum(
        coefficients,
        tuple(first_halves[index] for index in SATURATED_SUBSET_INDICES),
    )
    saturation_point = _pari_halves(
        coefficients, (subset_sum,), timeout=gp_timeout
    )[0]
    basis = (saturation_point,) + first_halves[1:]
    if short_add(coefficients, saturation_point, saturation_point) != subset_sum:
        raise AssertionError("the saturation relation failed exact doubling")
    return coefficients, images, differences, first_halves, basis


def build_certificate(*, gp_timeout: float = 20.0) -> dict[str, object]:
    coefficients, images, differences, first_halves, basis = derive_saturated_basis(
        gp_timeout=gp_timeout
    )
    signatures = find_mod2_reduction_certificate(
        coefficients, basis, prime_bound=500
    )
    rank = combined_mod2_rank(signatures, len(basis))
    two_torsion_prime = find_two_torsion_certificate_prime(
        coefficients, prime_bound=200
    )
    if rank != 14:
        raise AssertionError("the saturated Kihara basis lacked full mod-2 rank")
    return {
        "schema": "elliptic-kihara-rank14-basis-v1",
        "reproduce_command": REPRODUCE_COMMAND,
        "parameter_t": "2",
        "short_jacobian_coefficients": [_q_string(value) for value in coefficients],
        "covariant_image_count": len(images),
        "covariant_images": [_point_record(point) for point in images],
        "origin": "covariant image of P15",
        "difference_images": [_point_record(point) for point in differences],
        "first_halves": [_point_record(point) for point in first_halves],
        "saturation_relation": {
            "twice_new_basis_point_1_equals_sum_first_halves": [
                index + 1 for index in SATURATED_SUBSET_INDICES
            ]
        },
        "saturated_basis": [_point_record(point) for point in basis],
        "no_rational_2_torsion_prime": two_torsion_prime,
        "mod2_reduction_signatures": [
            {
                "prime": signature.prime,
                "group_order": signature.group_order,
                "doubled_subgroup_order": signature.doubled_subgroup_order,
                "quotient_dimension": signature.quotient_dimension,
                "rows": [list(row) for row in signature.rows],
            }
            for signature in signatures
        ],
        "combined_mod2_rank": rank,
        "exact_specialized_rank_lower_bound": 14,
        "generic_section_consequence": (
            "P1-P15,...,P14-P15 are generically independent: a generic relation "
            "would specialize at smooth t=2 and then contradict independence of "
            "their differences under the nonconstant covariant map"
        ),
        "status": "exact finite-reduction certificate",
    }


def signature_from_record(record: dict[str, object]) -> Mod2ReductionSignature:
    return Mod2ReductionSignature(
        prime=int(record["prime"]),
        group_order=int(record["group_order"]),
        doubled_subgroup_order=int(record["doubled_subgroup_order"]),
        quotient_dimension=int(record["quotient_dimension"]),
        rows=tuple(tuple(int(value) for value in row) for row in record["rows"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gp-timeout", type=float, default=20.0)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    certificate = build_certificate(gp_timeout=arguments.gp_timeout)
    encoded = json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    if arguments.output:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
