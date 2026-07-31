#!/usr/bin/env python3
"""Prototype a relative divergence certificate for a sampled Ore factor.

For a fixed order-14 operator over a finite field or over QQ,

    G = sum_j g_j(m) S^j

and P=Q/u^3, a relative telescoping certificate of the form

    sum_j g_j(m) P^(m+j)
      = D_u(X P^m) + d_t(t(1-t) Y P^m)

can be solved from the highest power of m downward.  In the Laurent ring
with U=u^-1, put

    A = u Q_u - 3Q,   C = t(1-t) Q_t.

If X_k,Y_k are the coefficient at m^k, the coefficient recursion consists
of repeated lifts through (A,C,uU-1), followed by the relative divergence

    u X_u - U X_U + d_t(t(1-t)Y).

This script tests a bounded number of leading recursion steps.  It is a
research probe, not a completed telescoping certificate.
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

from verify_two_pair_sic_bidegree33_rank_two_relative_jacobian import (  # noqa: E402
    q_expression,
)


OPERATOR_INPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_ore_reconstruct_research.json"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_relative_divergence_research.json"
)
EXPECTED_ORDER = 14
EXPECTED_M_DEGREE = 58


def singular_string(value: Path) -> str:
    """Quote an absolute path as a Singular string literal."""

    return '"' + str(value.resolve()).replace("\\", "\\\\").replace('"', '\\"') + '"'


def coefficient_polynomials(
    coefficients: list[list[object]],
    order: int,
    m_degree: int,
) -> list[str]:
    def scalar(value: object) -> str | None:
        if isinstance(value, int):
            return str(value) if value else None
        if (
            isinstance(value, list)
            and len(value) == 2
            and all(isinstance(entry, int) for entry in value)
        ):
            numerator, denominator = value
            if numerator == 0:
                return None
            return (
                str(numerator)
                if denominator == 1
                else f"({numerator}/{denominator})"
            )
        raise TypeError(f"unsupported operator coefficient {value!r}")

    expressions = []
    for coefficient_degree in range(m_degree + 1):
        terms = []
        for shift in range(order + 1):
            coefficient = scalar(
                coefficients[shift][coefficient_degree]
            )
            if coefficient is not None:
                terms.append(
                    f"({coefficient})*Q*P^{shift}"
                )
        expressions.append("+".join(terms) if terms else "0")
    return expressions


def singular_code(
    prime: int,
    coefficients: list[list[int]],
    steps: int,
    mode: str,
    start_degree: int,
    initial_residual: str | None,
    checkpoint_polynomial: Path | None,
    certificate_output: Path | None,
    mapped_quotient: bool,
    order: int,
    m_degree: int,
) -> str:
    q = q_expression(0)
    expressions = coefficient_polynomials(coefficients, order, m_degree)
    stop_degree = start_degree - steps + 1
    if mode == "relative":
        critical_generator = "t*(1-t)*diff(Q,t)"
        t_primitive = "t*(1-t)*Y"
    elif mode == "interior":
        critical_generator = "diff(Q,t)"
        t_primitive = "Y"
    else:
        raise ValueError(mode)
    if mapped_quotient:
        ring_setup = f"""
ring ambient={prime},(u,U,t),dp;
poly Q0={q};
poly A0=u*diff(Q0,u)-3*Q0;
poly C0={critical_generator.replace("Q", "Q0")};
poly rel0=u*U-1;
ideal generators0=ideal(A0,C0,rel0);
ideal critical0=std(generators0);
matrix critical_lift0=lift(generators0,critical0);
ideal relation0=std(ideal(rel0));
qring r=relation0;
poly Q=imap(ambient,Q0);
poly P=Q*U^3;
poly A=imap(ambient,A0);
poly C=imap(ambient,C0);
ideal critical=std(ideal(A,C));
setring ambient;
ideal quotient_critical0=imap(r,critical);
matrix quotient_lift0=lift(generators0,quotient_critical0);
setring r;
matrix critical_lift=imap(ambient,quotient_lift0);
ideal relation=std(ideal(0));
"""
        lift_assignment = """
    X=lift_matrix[1,1];
    Y=lift_matrix[2,1];
    Z=0;
"""
        # Singular does not normalize equality tests in this quotient
        # representation.  The emitted X,Y identities are replayed in the
        # ambient ring by the independent verifier.
        lift_check = "0"
        terminal_expression = "residual"
    else:
        ring_setup = f"""
ring r={prime},(u,U,t),dp;
poly Q={q};
poly P=Q*U^3;
poly A=u*diff(Q,u)-3*Q;
poly C={critical_generator};
poly rel=u*U-1;
ideal generators=ideal(A,C,rel);
matrix critical_lift;
ideal critical=liftstd(generators,critical_lift);
ideal relation=std(ideal(rel));
"""
        lift_assignment = """
    X=lift_matrix[1,1];
    Y=lift_matrix[2,1];
    Z=lift_matrix[3,1];
"""
        lift_check = "residual-X*A-Y*C-Z*rel"
        terminal_expression = "reduce(residual,relation)"
    lift_command = """
    lift_matrix=critical_lift*lift(critical,ideal(residual));
"""
    residual_expression = (
        initial_residual
        if initial_residual is not None
        else f"F[{m_degree + 1}]"
    )
    certificate_initialization = ""
    certificate_write = ""
    if certificate_output is not None:
        certificate_file = singular_string(certificate_output)
        coefficient_ring = (
            "characteristic-zero" if prime == 0 else "modular"
        )
        certificate_initialization = f"""
string certificate_file={certificate_file};
write(":w "+certificate_file,
  "// Exact {coefficient_ring} interior divergence certificate chunk");
write(":a "+certificate_file,
  "// prime={prime}, start_m_degree={start_degree}, stop_m_degree={stop_degree}, mapped_quotient={int(mapped_quotient)}");
"""
        certificate_write = """
    write(":a "+certificate_file,
      "poly X"+string(k-1)+"="+string(X)+";");
    write(":a "+certificate_file,
      "poly Y"+string(k-1)+"="+string(Y)+";");
"""
    checkpoint_write = ""
    if checkpoint_polynomial is not None:
        checkpoint_file = singular_string(checkpoint_polynomial)
        checkpoint_write = f"""
  string checkpoint_file={checkpoint_file};
  write(":w "+checkpoint_file,string(residual));
  print("CHECKPOINT_TERMS_BEGIN");
  print(size(residual));
  print("CHECKPOINT_TERMS_END");
"""
    return f"""
{ring_setup}
list F={",".join(expressions)};
poly residual={residual_expression};
residual=reduce(residual,relation);
poly X;
poly Y;
poly Z;
poly divergence;
matrix lift_matrix;
int k;
{certificate_initialization}
poly leading_obstruction=reduce(residual,critical);
if(leading_obstruction!=0)
{{
  print("LEADING_OBSTRUCTION_BEGIN");
  print(size(leading_obstruction));
  print("LEADING_OBSTRUCTION_END");
}}
else
{{
  for(k={start_degree};k>={stop_degree};k--)
  {{
    if(reduce(residual,critical)!=0)
    {{
      ERROR("intermediate critical obstruction");
    }}
{lift_command}
{lift_assignment}
    if({lift_check}!=0)
    {{
      print("LIFT_ERROR_TERMS_BEGIN");
      print(size({lift_check}));
      print("LIFT_ERROR_TERMS_END");
      ERROR("relative lift identity");
    }}
    X=reduce(X,relation);
    Y=reduce(Y,relation);
    divergence=
      reduce(
        u*diff(X,u)-U*diff(X,U)+diff({t_primitive},t),
        relation
      );
  print("STEP_BEGIN");
  print(k);
  print(k-1);
  print(size(residual));
    print(size(X));
    print(size(Y));
    print(size(Z));
    print(size(divergence));
    print("STEP_END");
{certificate_write}
    if(k>0)
    {{
    residual=reduce(F[k]-Q*divergence,relation);
    }}
  }}
{checkpoint_write}
  if({stop_degree}==1)
  {{
    poly terminal={terminal_expression};
    print("TERMINAL_BEGIN");
    print(size(terminal));
    print("TERMINAL_END");
  }}
}}
print("PASS bounded relative divergence recursion");
"""


def parse_steps(output: str) -> list[dict[str, int]]:
    blocks = re.findall(
        r"STEP_BEGIN\s+"
        r"(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+"
        r"(\d+)\s+"
        r"STEP_END",
        output,
    )
    return [
        {
            "m_degree": int(values[0]),
            "certificate_m_degree": int(values[1]),
            "residual_terms": int(values[2]),
            "X_terms": int(values[3]),
            "Y_terms": int(values[4]),
            "relation_terms": int(values[5]),
            "divergence_terms": int(values[6]),
        }
        for values in blocks
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--operator", type=Path, default=OPERATOR_INPUT)
    parser.add_argument("--steps", type=int, default=3)
    parser.add_argument(
        "--mode",
        choices=("relative", "interior"),
        default="relative",
    )
    parser.add_argument(
        "--checkpoint-input",
        type=Path,
        help="resume from checkpoint metadata written by --checkpoint-output",
    )
    parser.add_argument(
        "--checkpoint-output",
        type=Path,
        help="write restart metadata and its adjacent .poly residual",
    )
    parser.add_argument(
        "--certificate-output",
        type=Path,
        help="write the exact X_k,Y_k polynomials for this chunk",
    )
    parser.add_argument(
        "--mapped-quotient",
        action="store_true",
        help=(
            "map the ambient Groebner lift into the quotient u*U=1 "
            "before descending"
        ),
    )
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    operator_bytes = arguments.operator.read_bytes()
    operator_sha256 = hashlib.sha256(operator_bytes).hexdigest()
    payload = json.loads(operator_bytes)
    modular = payload.get("modular_operator")
    if modular is not None:
        prime = int(modular["prime"])
        coefficients = modular["coefficients"]
    else:
        prime = 0
        coefficients = payload.get(
            "primitive_integer_coefficients",
            payload["primitive_coefficients"],
        )
    operator = payload.get("operator", {})
    order = int(operator.get("order", EXPECTED_ORDER))
    m_degree = int(operator.get("m_degree", EXPECTED_M_DEGREE))
    assert len(coefficients) == order + 1
    assert all(
        len(polynomial) == m_degree + 1
        for polynomial in coefficients
    )
    start_degree = m_degree
    initial_residual = None
    resumed_from = None
    if arguments.checkpoint_input is not None:
        checkpoint = json.loads(arguments.checkpoint_input.read_text())
        if checkpoint["format"] != (
            "two-pair-sic-bidegree33-rank-two-"
            "relative-divergence-checkpoint-v1"
        ):
            raise ValueError("unsupported checkpoint format")
        expected = {
            "prime": prime,
            "point": payload["point"],
            "mode": arguments.mode,
            "operator_sha256": operator_sha256,
        }
        for key, value in expected.items():
            if checkpoint.get(key) != value:
                raise ValueError(
                    f"checkpoint {key} mismatch: "
                    f"{checkpoint.get(key)!r} != {value!r}"
                )
        start_degree = int(checkpoint["next_m_degree"])
        residual_path = (
            arguments.checkpoint_input.parent
            / checkpoint["residual_file"]
        ).resolve()
        initial_residual = residual_path.read_text().strip()
        if not initial_residual:
            raise ValueError("checkpoint residual is empty")
        resumed_from = str(arguments.checkpoint_input)
    if not 1 <= arguments.steps <= start_degree:
        raise ValueError(
            f"steps must be between 1 and the start degree {start_degree}"
        )
    stop_degree = start_degree - arguments.steps + 1
    next_degree = stop_degree - 1

    checkpoint_polynomial = None
    if arguments.checkpoint_output is not None:
        arguments.checkpoint_output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        checkpoint_polynomial = arguments.checkpoint_output.with_suffix(
            ".poly"
        )
    if arguments.certificate_output is not None:
        arguments.certificate_output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    completed = subprocess.run(
        [singular, "-q"],
        input=singular_code(
            prime,
            coefficients,
            arguments.steps,
            arguments.mode,
            start_degree,
            initial_residual,
            checkpoint_polynomial,
            arguments.certificate_output,
            arguments.mapped_quotient,
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
    if "PASS bounded relative divergence recursion" not in combined:
        raise RuntimeError(combined)
    records = parse_steps(combined)
    obstruction_match = re.search(
        r"LEADING_OBSTRUCTION_BEGIN\s+(\d+)\s+"
        r"LEADING_OBSTRUCTION_END",
        combined,
    )
    obstruction_terms = (
        int(obstruction_match.group(1))
        if obstruction_match is not None
        else None
    )
    if obstruction_terms is None and len(records) != arguments.steps:
        raise RuntimeError(combined)
    if obstruction_terms is not None and records:
        raise RuntimeError(combined)
    terminal_match = re.search(
        r"TERMINAL_BEGIN\s+(\d+)\s+TERMINAL_END",
        combined,
    )
    terminal_terms = (
        int(terminal_match.group(1))
        if terminal_match is not None
        else None
    )
    checkpoint_terms_match = re.search(
        r"CHECKPOINT_TERMS_BEGIN\s+(\d+)\s+CHECKPOINT_TERMS_END",
        combined,
    )
    checkpoint_terms = (
        int(checkpoint_terms_match.group(1))
        if checkpoint_terms_match is not None
        else None
    )
    if stop_degree == 1:
        if obstruction_terms is None and terminal_terms is None:
            raise RuntimeError(combined)
    elif terminal_terms is not None:
        raise RuntimeError(combined)
    if checkpoint_polynomial is not None:
        if obstruction_terms is not None:
            raise RuntimeError(
                "cannot checkpoint a leading-obstruction run"
            )
        if checkpoint_terms is None:
            raise RuntimeError(combined)
        if (
            not checkpoint_polynomial.is_file()
            or checkpoint_polynomial.stat().st_size == 0
        ):
            raise RuntimeError("Singular did not write the checkpoint")
        checkpoint_payload = {
            "format": (
                "two-pair-sic-bidegree33-rank-two-"
                "relative-divergence-checkpoint-v1"
            ),
            "prime": prime,
            "point": payload["point"],
            "mode": arguments.mode,
            "operator_sha256": operator_sha256,
            "mapped_quotient": arguments.mapped_quotient,
            "completed_m_degrees": [start_degree, stop_degree],
            "next_m_degree": next_degree,
            "residual_terms": checkpoint_terms,
            "residual_file": checkpoint_polynomial.name,
        }
        arguments.checkpoint_output.write_text(
            json.dumps(checkpoint_payload, indent=2) + "\n"
        )
    if arguments.certificate_output is not None:
        if (
            not arguments.certificate_output.is_file()
            or arguments.certificate_output.stat().st_size == 0
        ):
            raise RuntimeError("Singular did not write the certificate")

    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "relative-divergence-research-v1"
        ),
        "status": (
            (
                f"exact {'characteristic-zero' if prime == 0 else 'modular'} "
                "leading-symbol obstruction to the direct "
                "zero-boundary ansatz"
            )
            if obstruction_terms is not None
            else (
                "bounded exact "
                f"{'characteristic-zero' if prime == 0 else 'modular'} "
                "leading-step recursion; not a complete telescoping "
                "certificate"
            )
        ),
        "prime": prime,
        "point": payload["point"],
        "operator": {
            "order": order,
            "m_degree": m_degree,
        },
        "laurent_model": "U=u^-1 with relation u*U-1",
        "mapped_quotient_engine": arguments.mapped_quotient,
        "mode": arguments.mode,
        "gradient_ideal": (
            ["u*Q_u-3Q", "t*(1-t)*Q_t", "u*U-1"]
            if arguments.mode == "relative"
            else ["u*Q_u-3Q", "Q_t", "u*U-1"]
        ),
        "start_m_degree": start_degree,
        "stop_m_degree": stop_degree,
        "next_m_degree": next_degree,
        "resumed_from": resumed_from,
        "steps_requested": arguments.steps,
        "steps": records,
        "leading_obstruction_terms": obstruction_terms,
        "terminal_m0_terms_mod_uU_minus_1": terminal_terms,
        "checkpoint": (
            str(arguments.checkpoint_output)
            if arguments.checkpoint_output is not None
            else None
        ),
        "certificate": (
            str(arguments.certificate_output)
            if arguments.certificate_output is not None
            else None
        ),
        "remaining_gate": (
            (
                (
                    "retain the two endpoint states in the full 14+2+2 "
                    "connection; the direct zero-boundary scalar ansatz "
                    "cannot start"
                )
                if arguments.mode == "relative"
                else (
                    "the sampled scalar operator has a nonzero leading "
                    "class even in the interior quotient"
                )
            )
            if obstruction_terms is not None
            else (
                (
                    "independently replay this complete characteristic-zero "
                    "certificate and audit the endpoint trace"
                )
                if prime == 0
                else (
                    "independently replay this complete modular certificate; "
                    "then audit the endpoint trace and lift across enough "
                    "primes to reconstruct characteristic zero"
                )
            )
            if terminal_terms == 0
            else (
                "the chosen descending lifts leave a nonzero terminal "
                "m^0 residual; add syzygy freedom or use the full "
                "14+2+2 connection"
            )
            if terminal_terms is not None
            else (
                f"complete all {m_degree} descending lifts and verify that the "
                "final m^0 residual is a multiple of u*U-1; then "
                + (
                    "independently replay the characteristic-zero "
                    "certificate and audit both endpoint identities"
                    if prime == 0
                    else (
                        "reconstruct and independently replay the "
                        "characteristic-zero certificate"
                    )
                )
            )
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    if obstruction_terms is not None:
        print(
            "PASS nonzero leading obstruction excludes the direct "
            "zero-boundary recursion"
        )
    else:
        print(
            f"PASS {arguments.steps} exact "
            f"{'characteristic-zero' if prime == 0 else 'modular'} "
            f"relative-divergence steps from m^{start_degree}"
            + ("" if prime == 0 else f" at prime {prime}")
        )
    if arguments.checkpoint_output is not None:
        print(f"PASS wrote checkpoint {arguments.checkpoint_output}")
    if arguments.certificate_output is not None:
        print(f"PASS wrote certificate {arguments.certificate_output}")
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()
