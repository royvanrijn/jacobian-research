#!/usr/bin/env python3
"""Exclude the maximal Case-2 derivative-gcd stratum exactly.

The certified Case-2 residue has ``deg(C),deg(G)=(8,12)``.  If
``deg(gcd(C',G'))=7``, then ``C'`` divides ``G'``.  After solving
``(J3),(J2)`` and the determined part of ``(J1)``, divide ``G'`` by ``C'``.
Three of the seven remainder coefficients, together with the coefficient of
``t^19`` in ``(J0)=B*G'-C'*F``, generate the unit ideal over the exact
degree-35 first-block field.

No residual ``(J1)`` compatibility equation is used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

from audit_case2_residue_strata import (
    default_replay_root,
    load_exact_core,
    singular_minpoly,
    solve_case2_through_j1,
)


EXPECTED_INPUT_SHA256 = (
    "1ac0b4db7ddd0b50fcbef6c93d49c28f"
    "7f80cb4133a73b5f2158af6c78f3b069"
)


@dataclass(frozen=True)
class MaximalGcdAudit:
    gcd_degree: int
    total_remainder_coefficients: int
    selected_remainder_degrees: tuple[int, ...]
    j0_coefficient_degree: int
    selected_constraint_count: int
    selected_constraint_term_counts: tuple[int, ...]
    selected_constraint_parameter_degrees: tuple[int, ...]
    uses_j1_compatibility: bool
    uses_j0: bool
    singular_input_sha256: str
    unit_ideal: bool


def maximal_gcd_remainders(ec, C, G):
    divisor = ec.tder(C)
    remainder = list(ec.tder(G))
    divisor_degree = len(divisor) - 1
    leading = divisor[divisor_degree].d[(0,) * ec.NP]
    leading_inverse = leading.inv()
    for degree in range(len(remainder) - 1, divisor_degree - 1, -1):
        coefficient = remainder[degree] * leading_inverse
        shift = degree - divisor_degree
        for index, value in enumerate(divisor):
            remainder[index + shift] -= coefficient * value
    if any(remainder[divisor_degree:]):
        raise AssertionError("exact division of G' by C' failed")
    result = tuple(remainder[:divisor_degree])
    if len(result) != 7 or any(not polynomial for polynomial in result):
        raise AssertionError("unexpected maximal-gcd remainder support")
    return result


def selected_constraints(ec, B, C, F, G):
    remainders = maximal_gcd_remainders(ec, C, G)
    j0 = ec.tadd(
        ec.tmul(B, ec.tder(G)),
        ec.tscale(ec.tmul(ec.tder(C), F), ec.K(-1)),
    )
    if len(j0) != 20 or not j0[19]:
        raise AssertionError("the terminal J0 coefficient is missing")
    return remainders, (
        remainders[0],
        remainders[1],
        remainders[2],
        j0[19],
    )


def singular_source(ec, constraints) -> str:
    return "\n".join(
        (
            'LIB "resources.lib";',
            "Resources::setcores(1);",
            'LIB "nfmodstd.lib";',
            "ring R=(0,u),(r,s,h),dp;",
            f"minpoly={singular_minpoly(ec)};",
            "option(redSB);",
            "ideal I=",
            ",\n".join(polynomial.sing(3) for polynomial in constraints)
            + ";",
            "ideal S=nfmodStd(I);",
            (
                'if(size(S)>0 && S[1]==1)'
                '{print("CASE2_MAXIMAL_GCD_UNIT");}'
            ),
            "quit;",
            "",
        )
    )


def run_audit(ec, B, C, F, G, singular: str, timeout: int):
    remainders, constraints = selected_constraints(ec, B, C, F, G)
    term_counts = tuple(len(polynomial.d) for polynomial in constraints)
    parameter_degrees = tuple(
        max(map(sum, polynomial.d)) for polynomial in constraints
    )
    if term_counts != (155, 155, 155, 9):
        raise AssertionError("maximal-gcd constraint term counts drifted")
    if parameter_degrees != (13, 13, 13, 4):
        raise AssertionError("maximal-gcd constraint degrees drifted")

    source = singular_source(ec, constraints)
    digest = hashlib.sha256(source.encode()).hexdigest()
    if EXPECTED_INPUT_SHA256 and digest != EXPECTED_INPUT_SHA256:
        raise AssertionError("maximal-gcd Singular input hash drifted")

    with tempfile.TemporaryDirectory(prefix="jc2-max-gcd-") as directory:
        path = Path(directory) / "maximal_gcd.sing"
        path.write_text(source)
        completed = subprocess.run(
            [singular, "-q", str(path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    marker = "CASE2_MAXIMAL_GCD_UNIT"
    if completed.returncode != 0 or marker not in completed.stdout:
        raise RuntimeError(
            "maximal-gcd exact elimination failed\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )

    return MaximalGcdAudit(
        gcd_degree=7,
        total_remainder_coefficients=len(remainders),
        selected_remainder_degrees=(0, 1, 2),
        j0_coefficient_degree=19,
        selected_constraint_count=len(constraints),
        selected_constraint_term_counts=term_counts,
        selected_constraint_parameter_degrees=parameter_degrees,
        uses_j1_compatibility=False,
        uses_j0=True,
        singular_input_sha256=digest,
        unit_ideal=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replay-root", type=Path)
    parser.add_argument("--singular", default=shutil.which("Singular"))
    parser.add_argument("--timeout", type=int, default=300)
    args = parser.parse_args()
    if not args.singular:
        raise RuntimeError("Singular is required")

    replay_root = (args.replay_root or default_replay_root()).resolve()
    ec = load_exact_core(replay_root)
    B, _, C, F, G, compatibility = solve_case2_through_j1(ec)
    if not compatibility:
        raise AssertionError("the unused J1 compatibility system is missing")
    audit = run_audit(ec, B, C, F, G, args.singular, args.timeout)
    print(json.dumps(asdict(audit), sort_keys=True))
    print("CASE2_MAXIMAL_GCD_STRATUM_PASS")


if __name__ == "__main__":
    main()
