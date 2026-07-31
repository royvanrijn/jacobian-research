#!/usr/bin/env python3
"""Independently replay an interior-divergence certificate chunk.

The research producer writes exact polynomials X_r,Y_r and a restart
residual.  This verifier reconstructs the sampled Ore target directly from
the operator artifact, checks every coefficient identity modulo u*U-1, and
checks that the recomputed final residual is exactly the checkpoint value.

The input may be a finite-field operator or the reconstructed primitive
characteristic-zero operator.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from research_two_pair_sic_bidegree33_rank_two_relative_divergence import (  # noqa: E402
    EXPECTED_M_DEGREE,
    EXPECTED_ORDER,
    OPERATOR_INPUT,
    coefficient_polynomials,
    singular_string,
)
from verify_two_pair_sic_bidegree33_rank_two_relative_jacobian import (  # noqa: E402
    q_expression,
)


def checkpoint(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if payload["format"] != (
        "two-pair-sic-bidegree33-rank-two-"
        "relative-divergence-checkpoint-v1"
    ):
        raise ValueError(f"unsupported checkpoint format in {path}")
    return payload


def singular_code(
    prime: int,
    coefficients: list[list[object]],
    certificate_path: Path,
    start_degree: int,
    stop_degree: int,
    initial_residual: str | None,
    claimed_residual: str | None,
    recovered_residual: Path | None,
    order: int,
    m_degree: int,
) -> str:
    expressions = coefficient_polynomials(coefficients, order, m_degree)
    q = q_expression(0)
    residual = (
        initial_residual
        if initial_residual is not None
        else f"F[{m_degree + 1}]"
    )
    steps = []
    for degree in range(start_degree, stop_degree - 1, -1):
        certificate_degree = degree - 1
        steps.append(
            f"""
X=X{certificate_degree};
Y=Y{certificate_degree};
poly lift_error=reduce(residual-X*A-Y*C,relation);
if(lift_error!=0){{ERROR("certificate lift at m^{degree}");}}
divergence=u*diff(X,u)-U*diff(X,U)+diff(Y,t);
print("VERIFIED_STEP_BEGIN");
print({degree});
print(size(subst(Y,t,0)));
print(size(subst(Y,t,1)));
print("VERIFIED_STEP_END");
residual=reduce(F[{degree}]-Q*divergence,relation);
"""
        )
    terminal = ""
    if stop_degree == 1:
        terminal = """
poly terminal=reduce(residual,relation);
print("TERMINAL_BEGIN");
print(size(terminal));
print("TERMINAL_END");
"""
    if claimed_residual is not None:
        checkpoint_check = f"""
poly claimed_residual={claimed_residual};
if(residual-claimed_residual!=0)
{{
  ERROR("checkpoint residual mismatch");
}}
"""
    elif recovered_residual is not None:
        checkpoint_check = f"""
write(":w "+{singular_string(recovered_residual)},string(residual));
"""
    else:
        raise ValueError("a claimed or recovered residual is required")
    return f"""
ring r={prime},(u,U,t),dp;
poly Q={q};
poly P=Q*U^3;
poly A=u*diff(Q,u)-3*Q;
poly C=diff(Q,t);
poly rel=u*U-1;
ideal relation=std(ideal(rel));
list F={",".join(expressions)};
execute(read({singular_string(certificate_path)}));
poly residual={residual};
residual=reduce(residual,relation);
poly X;
poly Y;
poly divergence;
{"".join(steps)}
{checkpoint_check}
print("CHECKPOINT_MATCH_BEGIN");
print(size(residual));
print("CHECKPOINT_MATCH_END");
{terminal}
print("PASS exact certificate chunk");
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--operator", type=Path, default=OPERATOR_INPUT)
    parser.add_argument("--certificate", type=Path, required=True)
    parser.add_argument(
        "--input-checkpoint",
        type=Path,
        help="checkpoint immediately before this chunk; omit at m^58",
    )
    parser.add_argument(
        "--output-checkpoint",
        type=Path,
        help="checkpoint written immediately after this chunk",
    )
    parser.add_argument(
        "--recover-checkpoint",
        type=Path,
        help="independently reconstruct a missing checkpoint from the chunk",
    )
    parser.add_argument("--timeout", type=int, default=900)
    arguments = parser.parse_args()
    if (
        (arguments.output_checkpoint is None)
        == (arguments.recover_checkpoint is None)
    ):
        raise ValueError(
            "give exactly one of --output-checkpoint or "
            "--recover-checkpoint"
        )

    operator_bytes = arguments.operator.read_bytes()
    operator_sha256 = hashlib.sha256(operator_bytes).hexdigest()
    operator = json.loads(operator_bytes)
    modular = operator.get("modular_operator")
    if modular is not None:
        prime = int(modular["prime"])
        coefficients = modular["coefficients"]
    else:
        prime = 0
        coefficients = operator.get(
            "primitive_integer_coefficients",
            operator["primitive_coefficients"],
        )
    operator_shape = operator.get("operator", {})
    order = int(operator_shape.get("order", EXPECTED_ORDER))
    m_degree = int(operator_shape.get("m_degree", EXPECTED_M_DEGREE))
    assert len(coefficients) == order + 1
    assert all(
        len(polynomial) == m_degree + 1
        for polynomial in coefficients
    )

    expected = {
        "prime": prime,
        "point": operator["point"],
        "mode": "interior",
        "operator_sha256": operator_sha256,
    }

    initial_residual = None
    input_checkpoint = None
    if arguments.input_checkpoint is not None:
        input_checkpoint = checkpoint(arguments.input_checkpoint)
        for key, value in expected.items():
            if input_checkpoint.get(key) != value:
                raise ValueError(f"input checkpoint {key} mismatch")
        input_residual_path = (
            arguments.input_checkpoint.parent
            / input_checkpoint["residual_file"]
        )
        initial_residual = input_residual_path.read_text().strip()
    start_degree = (
        int(input_checkpoint["next_m_degree"])
        if input_checkpoint is not None
        else m_degree
    )
    output_checkpoint = None
    claimed_residual = None
    recovered_residual = None
    if arguments.output_checkpoint is not None:
        output_checkpoint = checkpoint(arguments.output_checkpoint)
        for key, value in expected.items():
            if output_checkpoint.get(key) != value:
                raise ValueError(f"output checkpoint {key} mismatch")
        recorded_start, stop_degree = (
            int(value)
            for value in output_checkpoint["completed_m_degrees"]
        )
        if recorded_start != start_degree:
            raise ValueError("output checkpoint does not follow input")
        if output_checkpoint["next_m_degree"] != stop_degree - 1:
            raise ValueError(
                "output checkpoint degree range is inconsistent"
            )
        output_residual_path = (
            arguments.output_checkpoint.parent
            / output_checkpoint["residual_file"]
        )
        claimed_residual = output_residual_path.read_text().strip()
    else:
        certificate_degrees = []
        with arguments.certificate.open() as source:
            for line in source:
                match = re.match(r"poly Y(\d+)=", line)
                if match is not None:
                    certificate_degrees.append(int(match.group(1)))
        if not certificate_degrees:
            raise ValueError("certificate contains no Y polynomials")
        stop_degree = min(certificate_degrees) + 1
        expected_degrees = list(
            range(start_degree - 1, stop_degree - 2, -1)
        )
        if certificate_degrees != expected_degrees:
            raise ValueError(
                "certificate degrees are not a contiguous descending chunk"
            )
        arguments.recover_checkpoint.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        recovered_residual = arguments.recover_checkpoint.with_suffix(
            ".poly"
        )
    if arguments.input_checkpoint is None and start_degree != m_degree:
        raise ValueError(
            f"only the m^{m_degree} chunk can omit --input-checkpoint"
        )
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    completed = subprocess.run(
        [singular, "-q"],
        input=singular_code(
            prime,
            coefficients,
            arguments.certificate,
            start_degree,
            stop_degree,
            initial_residual,
            claimed_residual,
            recovered_residual,
            order,
            m_degree,
        ),
        text=True,
        capture_output=True,
        timeout=arguments.timeout,
        check=True,
    )
    combined = completed.stdout + completed.stderr
    if "?" in combined or "error occurred" in combined:
        raise RuntimeError(combined)
    if "PASS exact certificate chunk" not in combined:
        raise RuntimeError(combined)
    verified = [
        {
            "m_degree": int(values[0]),
            "Y_terms_at_t0": int(values[1]),
            "Y_terms_at_t1": int(values[2]),
        }
        for values in re.findall(
            r"VERIFIED_STEP_BEGIN\s+(\d+)\s+(\d+)\s+(\d+)\s+"
            r"VERIFIED_STEP_END",
            combined,
        )
    ]
    expected_steps = start_degree - stop_degree + 1
    if len(verified) != expected_steps:
        raise RuntimeError(combined)
    residual_match = re.search(
        r"CHECKPOINT_MATCH_BEGIN\s+(\d+)\s+CHECKPOINT_MATCH_END",
        combined,
    )
    if residual_match is None:
        raise RuntimeError(combined)
    residual_terms = int(residual_match.group(1))
    if (
        output_checkpoint is not None
        and residual_terms != output_checkpoint["residual_terms"]
    ):
        raise RuntimeError("checkpoint residual term count mismatch")
    terminal_match = re.search(
        r"TERMINAL_BEGIN\s+(\d+)\s+TERMINAL_END",
        combined,
    )
    terminal_terms = (
        int(terminal_match.group(1))
        if terminal_match is not None
        else None
    )
    if arguments.recover_checkpoint is not None:
        if (
            recovered_residual is None
            or not recovered_residual.is_file()
            or recovered_residual.stat().st_size == 0
        ):
            raise RuntimeError("Singular did not recover the residual")
        recovery_payload = {
            "format": (
                "two-pair-sic-bidegree33-rank-two-"
                "relative-divergence-checkpoint-v1"
            ),
            **expected,
            "recovered_by_independent_replay": True,
            "completed_m_degrees": [start_degree, stop_degree],
            "next_m_degree": stop_degree - 1,
            "residual_terms": residual_terms,
            "residual_file": recovered_residual.name,
        }
        arguments.recover_checkpoint.write_text(
            json.dumps(recovery_payload, indent=2) + "\n"
        )

    print(
        f"PASS independently verified m^{start_degree} through "
        f"m^{stop_degree} over "
        f"{'QQ' if prime == 0 else f'F_{prime}'}"
    )
    print(
        f"PASS exact checkpoint residual match ({residual_terms} terms)"
    )
    if terminal_terms is not None:
        print(
            "PASS terminal residual modulo u*U-1 has "
            f"{terminal_terms} terms"
        )
    print(
        "PASS retained endpoint traces at every verified "
        "certificate degree"
    )
    if arguments.recover_checkpoint is not None:
        print(
            f"PASS recovered checkpoint {arguments.recover_checkpoint}"
        )


if __name__ == "__main__":
    main()
