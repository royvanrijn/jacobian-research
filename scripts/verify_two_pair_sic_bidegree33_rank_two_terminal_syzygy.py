#!/usr/bin/env python3
"""Independently verify the modular terminal-syzygy correction."""

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

from research_two_pair_sic_bidegree33_rank_two_relative_divergence import (  # noqa: E402
    OPERATOR_INPUT,
    singular_string,
)
from verify_two_pair_sic_bidegree33_rank_two_relative_jacobian import (  # noqa: E402
    q_expression,
)


TERMINAL = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_interior_divergence_checkpoint_m0.poly"
)
R_CERTIFICATE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_terminal_syzygy_R.sing"
)


def singular_code(
    prime: int,
    terminal: str,
    certificate: Path,
) -> str:
    return f"""
ring r={prime},(u,U,t),dp;
ideal relation=std(ideal(u*U-1));
poly Q={q_expression(0)};
poly A=u*diff(Q,u)-3*Q;
poly C=diff(Q,t);
poly terminal={terminal};
execute(read({singular_string(certificate)}));
poly H=C*(u*diff(R,u)-U*diff(R,U)+3*R)-A*diff(R,t);
H=reduce(H,relation);
poly error=reduce(terminal-Q*H,relation);
print("VERIFY_BEGIN");
print(size(R));
print(size(H));
print(size(error));
print("VERIFY_END");
if(error!=0){{ERROR("terminal syzygy identity");}}
print("PASS terminal=Q*H(R)");
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--operator", type=Path, default=OPERATOR_INPUT)
    parser.add_argument("--terminal", type=Path, default=TERMINAL)
    parser.add_argument(
        "--certificate",
        type=Path,
        default=R_CERTIFICATE,
    )
    parser.add_argument("--timeout", type=int, default=300)
    arguments = parser.parse_args()

    operator = json.loads(arguments.operator.read_text())
    prime = int(operator["modular_operator"]["prime"])
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    completed = subprocess.run(
        [singular, "-q"],
        input=singular_code(
            prime,
            arguments.terminal.read_text().strip(),
            arguments.certificate,
        ),
        text=True,
        capture_output=True,
        timeout=arguments.timeout,
        check=True,
    )
    combined = completed.stdout + completed.stderr
    if "?" in combined or "error occurred" in combined:
        raise RuntimeError(combined)
    if "PASS terminal=Q*H(R)" not in combined:
        raise RuntimeError(combined)
    match = re.search(
        r"VERIFY_BEGIN\s+(\d+)\s+(\d+)\s+(\d+)\s+VERIFY_END",
        combined,
    )
    if match is None:
        raise RuntimeError(combined)
    r_terms, h_terms, error_terms = (
        int(value) for value in match.groups()
    )
    if error_terms != 0:
        raise RuntimeError("nonzero terminal-syzygy error")
    print(
        f"PASS exact modular terminal identity with {r_terms}-term R "
        f"and {h_terms}-term H(R)"
    )
    print("PASS terminal residual is removed in the Laurent quotient")


if __name__ == "__main__":
    main()
