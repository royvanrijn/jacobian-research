#!/usr/bin/env python3
"""Solve the leading syzygy correction for the sampled order-14 operator.

The naive descending divergence lift fails because the m^58 coefficient is
not in the relative gradient ideal.  A degree-58 certificate may nevertheless
start with the Koszul syzygy

    (X_58,Y_58) = R * (C,-A),

where A=u*Q_u-3Q and C=t*(1-t)*Q_t.  Its critical contribution cancels at
degree 59, while its divergence contributes at degree 58.

This bounded modular probe constructs the resulting linear map on the exact
18-dimensional saturated relative critical algebra and solves for R in the
known standard-monomial basis.  It is not a full telescoping certificate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from flint import nmod_mat


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
    / "two_pair_sic_bidegree33_rank_two_leading_syzygy_research.json"
)
BASIS = (
    (0, 0),
    (1, 0),
    (2, 0),
    (3, 0),
    (4, 0),
    (0, 1),
    (1, 1),
    (2, 1),
    (3, 1),
    (0, 2),
    (1, 2),
    (2, 2),
    (0, 3),
    (1, 3),
    (2, 3),
    (0, 4),
    (1, 4),
    (0, 5),
)
EXPECTED_ORDER = 14
EXPECTED_M_DEGREE = 58


def basis_polynomial(exponents: tuple[int, int]) -> str:
    u_degree, t_degree = exponents
    factors = []
    if u_degree:
        factors.append(f"u^{u_degree}")
    if t_degree:
        factors.append(f"t^{t_degree}")
    return "*".join(factors) if factors else "1"


def singular_code(
    prime: int,
    coefficients: list[list[int]],
) -> str:
    q = q_expression(0)
    leading = [
        coefficients[shift][EXPECTED_M_DEGREE]
        for shift in range(EXPECTED_ORDER + 1)
    ]
    target = "+".join(
        f"({coefficient})*Q*P^{shift}"
        for shift, coefficient in enumerate(leading)
        if coefficient
    )
    basis = ",".join(
        basis_polynomial(exponents) for exponents in BASIS
    )
    return f"""
LIB "elim.lib";
ring r={prime},(u,t),dp;
poly Q={q};
poly A=u*diff(Q,u)-3*Q;
poly T=t*(1-t);
poly C=T*diff(Q,t);
list saturation=sat_with_exp(ideal(A,C),ideal(u));
ideal I=std(saturation[1]);
if(vdim(I)!=18){{ERROR("relative quotient length");}}
ideal unit_system=ideal(u)+I;
matrix inverse_lift=lift(unit_system,ideal(1));
poly inverse_u=inverse_lift[1,1];
if(reduce(u*inverse_u-1,I)!=0){{ERROR("u inverse");}}
poly P=reduce(Q*inverse_u^3,I);
list basis={basis};

proc dump_poly(poly f)
{{
  intvec exponent;
  print("POLY_BEGIN");
  while(f!=0)
  {{
    exponent=leadexp(f);
    print(leadcoef(f));
    print(exponent[1]);
    print(exponent[2]);
    f=f-lead(f);
  }}
  print("POLY_END");
}}

poly target=reduce({target},I);
dump_poly(target);
int index;
poly R;
poly image;
for(index=1;index<=18;index++)
{{
  R=basis[index];
  image=reduce(
    Q*(u*diff(C*R,u)+diff(-T*A*R,t)),
    I
  );
  dump_poly(image);
}}
print("PASS leading syzygy quotient matrix");
"""


def parse_polynomials(output: str, prime: int) -> list[list[int]]:
    blocks = re.findall(
        r"POLY_BEGIN\s*(.*?)\s*POLY_END",
        output,
        flags=re.DOTALL,
    )
    index = {exponents: position for position, exponents in enumerate(BASIS)}
    vectors = []
    for block in blocks:
        values = [int(value) for value in block.split()]
        if len(values) % 3:
            raise RuntimeError(block)
        vector = [0] * len(BASIS)
        for offset in range(0, len(values), 3):
            coefficient, u_degree, t_degree = values[offset : offset + 3]
            exponents = (u_degree, t_degree)
            if exponents not in index:
                raise RuntimeError(
                    f"unexpected quotient monomial {exponents}"
                )
            vector[index[exponents]] = coefficient % prime
        vectors.append(vector)
    return vectors


def expression(coefficients: list[int], prime: int) -> str:
    centered = [
        value if value <= prime // 2 else value - prime
        for value in coefficients
    ]
    terms = []
    for coefficient, exponents in zip(centered, BASIS, strict=True):
        if coefficient:
            terms.append(
                f"({coefficient})*{basis_polynomial(exponents)}"
            )
    return "+".join(terms) if terms else "0"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--operator", type=Path, default=OPERATOR_INPUT)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    operator_payload = json.loads(arguments.operator.read_text())
    modular = operator_payload["modular_operator"]
    prime = int(modular["prime"])
    coefficients = modular["coefficients"]

    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    completed = subprocess.run(
        [singular, "-q"],
        input=singular_code(prime, coefficients),
        text=True,
        capture_output=True,
        timeout=arguments.timeout,
        check=True,
    )
    combined = completed.stdout + completed.stderr
    if "?" in combined or "error occurred" in combined:
        raise RuntimeError(combined)
    if "PASS leading syzygy quotient matrix" not in combined:
        raise RuntimeError(combined)
    vectors = parse_polynomials(combined, prime)
    if len(vectors) != len(BASIS) + 1:
        raise RuntimeError(combined)
    target = vectors[0]
    images = vectors[1:]
    rows = [
        [images[column][row] for column in range(len(BASIS))]
        for row in range(len(BASIS))
    ]
    matrix = nmod_mat(rows, prime)
    target_column = nmod_mat([[value] for value in target], prime)
    rank = matrix.rank()
    target_nonzero = any(target)
    if rank != 0:
        raise RuntimeError(f"unexpected leading syzygy rank {rank}")
    if not target_nonzero:
        raise RuntimeError("unexpected zero leading obstruction")

    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-"
            "leading-syzygy-research-v1"
        ),
        "status": "exact bounded modular no-go for the leading syzygy repair",
        "prime": prime,
        "point": operator_payload["point"],
        "operator": {
            "order": EXPECTED_ORDER,
            "m_degree": EXPECTED_M_DEGREE,
        },
        "relative_critical_length": len(BASIS),
        "leading_syzygy_map_rank": rank,
        "leading_obstruction_nonzero": target_nonzero,
        "basis_exponents_u_t": [list(exponents) for exponents in BASIS],
        "conclusion": (
            "the leading m^58 obstruction is nonzero, while every "
            "Koszul correction R*(C,-A) has zero divergence class; "
            "there is no direct zero-boundary polynomial certificate "
            "of the tested form"
        ),
        "remaining_gate": (
            "construct the full 14+2+2 m-dependent connection with "
            "endpoint states retained, rather than prescribing the "
            "sampled scalar factor in a zero-boundary ansatz"
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS leading m^58 relative obstruction is nonzero")
    print("PASS Koszul leading-syzygy divergence has quotient rank zero")
    print("PASS direct zero-boundary polynomial ansatz is excluded")
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()
