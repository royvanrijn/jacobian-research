#!/usr/bin/env python3
"""Exclude the Case-2 derivative-gcd degree-six stratum exactly.

The lower-jet audit forces ``t|gcd(C',G')``.  If the gcd has degree six,
write

    C' = H(t) (t+v),  deg(H)=6,  H(0)=0.

Synthetic division reconstructs ``H`` in the three exact high-layer
parameters and ``v``.  Every such common divisor satisfies ``C'(0)=0``,
``H(0)=0``, and all six coefficients of ``G' mod H``.  The first two
conditions, the last two remainder coefficients, and the coefficient of
``t^19`` in ``(J0)=B*G'-C'*F`` generate the unit ideal over the exact
degree-35 first-block field.  No residual ``(J1)`` compatibility equation
is used.
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
from probe_case2_gcd6 import (
    constraint_term_count,
    gcd6_constraints,
    j0_constraints,
    polynomial_singular,
)


EXPECTED_INPUT_SHA256 = (
    "08bca008eae9dfe29fd383406a1d5b80"
    "a5e75b934a453dff7560cc05cdd30033"
)


@dataclass(frozen=True)
class GcdSixAudit:
    gcd_degree: int
    linear_cofactor_parameter: str
    selected_factor_conditions: tuple[str, ...]
    total_g_remainder_coefficients: int
    selected_g_remainder_degrees: tuple[int, ...]
    j0_coefficient_degree: int
    selected_constraint_count: int
    selected_constraint_term_counts: tuple[int, ...]
    selected_constraint_parameter_degrees: tuple[int, ...]
    uses_j1_compatibility: bool
    uses_j0: bool
    singular_input_sha256: str
    unit_ideal: bool


def selected_constraints(ec, B, C, F, G):
    gcd_constraints = gcd6_constraints(ec, C, G)
    if len(gcd_constraints) != 8:
        raise AssertionError("unexpected gcd-six constraint count")
    j0 = dict(j0_constraints(ec, B, C, F, G))
    if 19 not in j0:
        raise AssertionError("the terminal J0 coefficient is missing")
    # Indices 0,1 are C'(0),H(0); indices 6,7 are remainder degrees 4,5.
    return (
        gcd_constraints[0],
        gcd_constraints[1],
        gcd_constraints[6],
        gcd_constraints[7],
        j0[19],
    )


def singular_source(ec, constraints) -> str:
    return "\n".join(
        (
            'LIB "resources.lib";',
            "Resources::setcores(1);",
            'LIB "nfmodstd.lib";',
            "ring R=(0,u),(r,s,h,v),dp;",
            f"minpoly={singular_minpoly(ec)};",
            "option(redSB);",
            "ideal I="
            + ",".join(
                polynomial_singular(polynomial)
                for polynomial in constraints
            )
            + ";",
            "ideal S=nfmodStd(I);",
            (
                'if(size(S)>0 && S[1]==1)'
                '{print("CASE2_GCD6_UNIT");}'
            ),
            "quit;",
            "",
        )
    )


def run_audit(ec, B, C, F, G, singular: str, timeout: int):
    constraints = selected_constraints(ec, B, C, F, G)
    term_counts = tuple(
        constraint_term_count(polynomial) for polynomial in constraints
    )
    parameter_degrees = tuple(
        max(
            v_degree + sum(monomial)
            for v_degree, parameter_polynomial
            in polynomial.coefficients.items()
            for monomial in parameter_polynomial.d
        )
        for polynomial in constraints
    )
    if term_counts != (5, 30, 656, 352, 9):
        raise AssertionError("gcd-six constraint term counts drifted")
    if parameter_degrees != (2, 7, 16, 15, 4):
        raise AssertionError("gcd-six constraint degrees drifted")

    source = singular_source(ec, constraints)
    digest = hashlib.sha256(source.encode()).hexdigest()
    if EXPECTED_INPUT_SHA256 and digest != EXPECTED_INPUT_SHA256:
        raise AssertionError("gcd-six Singular input hash drifted")

    with tempfile.TemporaryDirectory(prefix="jc2-gcd6-") as directory:
        path = Path(directory) / "gcd6.sing"
        path.write_text(source)
        completed = subprocess.run(
            [singular, "-q", str(path)],
            check=False,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    marker = "CASE2_GCD6_UNIT"
    if completed.returncode != 0 or marker not in completed.stdout:
        raise RuntimeError(
            "gcd-six exact elimination failed\n"
            f"stdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )

    return GcdSixAudit(
        gcd_degree=6,
        linear_cofactor_parameter="v in C'=H*(t+v)",
        selected_factor_conditions=("C'(0)=0", "H(0)=0"),
        total_g_remainder_coefficients=6,
        selected_g_remainder_degrees=(4, 5),
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
    print("CASE2_GCD6_STRATUM_PASS")


if __name__ == "__main__":
    main()
