#!/usr/bin/env python3
"""Compute an ambient fixed-fiber Picard--Fuchs ideal with Macaulay2.

This is a long research computation, not a verifier.  It rewrites the
moment generating integrand as

    u^2 / (u^3 - z Q(u,t))

and asks Macaulay2 for a first-order holonomic annihilator of that
specific rational function and its D-module pushforward in u,t.  The
pushforward concerns closed de-Rham classes.  An operator returned here
must still be audited against
the relative interval endpoints t=0,1 before it can certify the original
period.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_two_pair_sic_bidegree33_rank_two_relative_jacobian import (  # noqa: E402
    q_expression,
)


OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_dmodule_picard_fuchs_research.json"
)
LOG = OUTPUT.with_suffix(".log")
V_CHECKPOINT = OUTPUT.with_name(
    "two_pair_sic_bidegree33_rank_two_dmodule_picard_fuchs_v_checkpoint.m2"
)


def macaulay2_code(
    include_pushforward: bool,
    v_checkpoint: Path,
    strategy: str,
) -> str:
    q = q_expression(0).replace("t", "v")
    common = f"""
needsPackage "BernsteinSato";
R = QQ[z,u,v];
Q = {q};
H = u^3-z*Q;
print "FIRST_ORDER_ANNIHILATOR_START";
A1 = kOrderAnnFa(1,H,-1);
print("FIRST_ORDER_ANNIHILATOR_GENERATORS " | toString numgens A1);
print("FIRST_ORDER_HOLONOMIC " | toString isHolonomic A1);
print("FIRST_ORDER_RANK " | toString holonomicRank A1);
W = ring A1;
u = W_1;
M = W^1/A1;
print "CYCLIC_ANNIHILATOR_START";
K = ker map(M,W^1,matrix{{{{u^2}}}});
IratRaw = ideal flatten entries generators K;
print("CYCLIC_RAW_GENERATORS " | toString numgens IratRaw);
Irat = trim IratRaw;
print("CYCLIC_ANNIHILATOR_GENERATORS " | toString numgens Irat);
print("CYCLIC_HOLONOMIC " | toString isHolonomic Irat);
print("CYCLIC_RANK " | toString holonomicRank Irat);
"""
    if not include_pushforward:
        return common + """
print "ANNIHILATOR_ONLY_END";
exit 0;
"""
    m2_strategy = {
        "schreyer": "Schreyer",
        "vhomogenize": "Vhomogenize",
    }[strategy]
    return common + f"""
print "PUSHFORWARD_V_START";
Jv = DintegrationIdeal(
    Irat,{{0,0,1}},Strategy=>{m2_strategy}
    );
print("PUSHFORWARD_V_GENERATORS " | toString numgens Jv);
checkpoint = openOut "{v_checkpoint.as_posix()}";
checkpoint << "Wv = " << toExternalString ring Jv << ";" << endl;
checkpoint << "Jv = " << toExternalString Jv << ";" << endl;
checkpoint << close;
print "PUSHFORWARD_V_CHECKPOINT_WRITTEN";
print "PUSHFORWARD_U_START";
J = DintegrationIdeal(Jv,{{0,1}},Strategy=>{m2_strategy});
print "TELESCOPER_BEGIN";
print J;
print "TELESCOPER_END";
exit 0;
"""


def marked(output: str, start: str, end: str) -> str:
    begin = output.index(start) + len(start)
    finish = output.index(end, begin)
    return output[begin:finish].strip()


def marker_integer(output: str, marker: str) -> int:
    match = re.search(
        rf"(?m)^{re.escape(marker)}\s+(-?\d+)\s*$",
        output,
    )
    if match is None:
        raise RuntimeError(f"missing integer marker {marker}")
    return int(match.group(1))


def has_exact_marker(output: str, marker: str) -> bool:
    return any(line.strip() == marker for line in output.splitlines())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--timeout",
        type=int,
        default=1800,
        help="wall-clock limit in seconds; 0 disables the limit",
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--log", type=Path, default=LOG)
    parser.add_argument("--v-checkpoint", type=Path, default=V_CHECKPOINT)
    parser.add_argument(
        "--strategy",
        choices=("schreyer", "vhomogenize"),
        default="schreyer",
    )
    parser.add_argument("--annihilator-only", action="store_true")
    arguments = parser.parse_args()

    macaulay2 = shutil.which("M2")
    normaliz = shutil.which("normaliz")
    if macaulay2 is None:
        raise RuntimeError("Macaulay2 is required")
    if normaliz is None:
        raise RuntimeError(
            "Normaliz is required by Macaulay2's BernsteinSato package"
        )

    arguments.log.parent.mkdir(parents=True, exist_ok=True)
    arguments.v_checkpoint.parent.mkdir(parents=True, exist_ok=True)
    timeout = None if arguments.timeout == 0 else arguments.timeout
    try:
        with arguments.log.open("w") as log:
            completed = subprocess.run(
                [macaulay2, "--silent", "--stop"],
                input=macaulay2_code(
                    include_pushforward=not arguments.annihilator_only,
                    v_checkpoint=arguments.v_checkpoint.resolve(),
                    strategy=arguments.strategy,
                ),
                text=True,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=timeout,
                check=True,
            )
    except subprocess.TimeoutExpired as error:
        raise RuntimeError(
            f"Macaulay2 timed out after {arguments.timeout} seconds; "
            f"partial output is retained in {arguments.log}"
        ) from error
    output = arguments.log.read_text()
    required_markers = [
        "FIRST_ORDER_ANNIHILATOR_START",
        "FIRST_ORDER_HOLONOMIC true",
        "FIRST_ORDER_RANK 1",
        "CYCLIC_ANNIHILATOR_START",
        "CYCLIC_HOLONOMIC true",
        "CYCLIC_RANK 1",
    ]
    if arguments.annihilator_only:
        required_markers.append("ANNIHILATOR_ONLY_END")
    else:
        required_markers.extend(
            [
                "PUSHFORWARD_V_START",
                "PUSHFORWARD_V_CHECKPOINT_WRITTEN",
                "PUSHFORWARD_U_START",
                "TELESCOPER_BEGIN",
                "TELESCOPER_END",
            ]
        )
    for marker in required_markers:
        if not has_exact_marker(output, marker):
            raise RuntimeError(
                f"missing Macaulay2 marker {marker}\n{output}\n"
                f"{completed.stderr}"
            )

    payload = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "dmodule-picard-fuchs-research-v1"
        ),
        "status": (
            (
                "exact characteristic-zero rank-one rational-function "
                "annihilator at one rank-two fiber; pushforward pending"
            )
            if arguments.annihilator_only
            else (
                "exact characteristic-zero ambient D-module pushforward "
                "at one rank-two fiber; the relative endpoint audit "
                "remains"
            )
        ),
        "integrand": "u^2/(u^3-z*Q(u,t))",
        "point": 0,
        "software": {
            "Macaulay2": subprocess.run(
                [macaulay2, "--version"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip(),
            "normaliz": subprocess.run(
                [normaliz, "--version"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip(),
        },
        "denominator_annihilator": {
            "algorithm": "kOrderAnnFa(1,H,-1)",
            "generators": marker_integer(
                output,
                "FIRST_ORDER_ANNIHILATOR_GENERATORS",
            ),
            "holonomic": True,
            "holonomic_rank": 1,
        },
        "specific_rational_function_annihilator": {
            "numerator": "u^2",
            "raw_generators": marker_integer(
                output,
                "CYCLIC_RAW_GENERATORS",
            ),
            "generators": marker_integer(
                output,
                "CYCLIC_ANNIHILATOR_GENERATORS",
            ),
            "presentation_optimization": "exact ideal trim",
            "holonomic": True,
            "holonomic_rank": 1,
        },
        "retained_weyl_variables": ["z", "Dz"],
        "remaining_gate": (
            "extract an operator, compare it with the sampled order-14 "
            "shift factor, and verify the t=0,1 certificate boundaries "
            "for the interval period"
        ),
    }
    if arguments.annihilator_only:
        payload["pushforward_ideal"] = None
    else:
        payload["integration"] = {
            "method": "sequential",
            "resolution_strategy": arguments.strategy,
            "first_variable": "v=t",
            "first_weight": [0, 0, 1],
            "intermediate_generators": marker_integer(
                output,
                "PUSHFORWARD_V_GENERATORS",
            ),
            "second_variable": "u",
            "second_weight": [0, 1],
        }
        payload["pushforward_ideal"] = marked(
            output,
            "TELESCOPER_BEGIN",
            "TELESCOPER_END",
        )
    arguments.output.write_text(json.dumps(payload, indent=2) + "\n")
    print("PASS exact rational-function annihilator")
    if not arguments.annihilator_only:
        print("PASS exact ambient D-module pushforward")
    try:
        display_output = arguments.output.relative_to(ROOT)
    except ValueError:
        display_output = arguments.output
    print(f"PASS wrote {display_output}")


if __name__ == "__main__":
    main()
