#!/usr/bin/env python3
"""Relative-Jacobian and order-18 recurrence probe in bidegree (3,3).

The exact characteristic-zero part computes the logarithmic Jacobian
quotient of the beta-period polynomial at an integral exact-rank-two point.
The modular part checks the natural order-18 recurrence shape at two points
and three primes and audits its forward coefficient.

The recurrence remains a fitted modular identity, not a universal
creative-telescoping certificate.
"""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from verify_two_pair_sic_bidegree33_rank_two_holonomic_probe import (
    POINTS,
    matrix_product,
    substituted_polynomial,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "scripts"
    / "explore_two_pair_sic_bidegree33_rank_two_recurrence.cpp"
)
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_rank_two_relative_jacobian.json"
)
PRIMES = (1_000_003, 1_000_033, 1_000_037)
ORDER = 18
DEGREE = 18
MAXIMUM_MOMENT = 460
COMMON_FORWARD_SHIFTS = (44, 46, 47, 49, 50, 52, 53, 55)
EXPECTED_BASIS = (
    "1",
    "u",
    "u2",
    "u3",
    "u4",
    "t",
    "ut",
    "u2t",
    "u3t",
    "t2",
    "ut2",
    "u2t2",
    "t3",
    "ut3",
    "u2t3",
    "t4",
    "ut4",
    "t5",
)


def q_expression(point: int, prime: int | None = None) -> str:
    raw_u, raw_w = POINTS[point]
    u = [[Fraction(value) for value in row] for row in raw_u]
    w = [[Fraction(value) for value in row] for row in raw_w]
    polynomial = substituted_polynomial(matrix_product(u, w))
    terms = []
    for (u_exponent, t_exponent), coefficient in sorted(
        polynomial.items()
    ):
        u_exponent += 3
        if prime is None:
            scalar = (
                str(coefficient.numerator)
                if coefficient.denominator == 1
                else (
                    f"({coefficient.numerator}/"
                    f"{coefficient.denominator})"
                )
            )
        else:
            scalar = str(
                coefficient.numerator
                * pow(coefficient.denominator % prime, -1, prime)
                % prime
            )
        factors = []
        if u_exponent:
            factors.append(f"u^{u_exponent}")
        if t_exponent:
            factors.append(f"t^{t_exponent}")
        monomial = "*".join(factors) or "1"
        terms.append(f"{scalar}*{monomial}")
    return "+".join(terms)


def run_singular(code: str) -> str:
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    completed = subprocess.run(
        [singular, "-q"],
        input=code,
        text=True,
        capture_output=True,
        timeout=60,
        check=True,
    )
    if "?" in completed.stdout or completed.stderr:
        raise RuntimeError(completed.stdout + completed.stderr)
    return completed.stdout


def exact_relative_jacobian() -> dict[str, object]:
    output = run_singular(
        f"""
LIB "elim.lib";
ring r=0,(u,t),dp;
poly Q={q_expression(0)};
poly A=u*diff(Q,u)-3*Q;
poly C=t*(1-t)*diff(Q,t);
ideal raw=std(A,C);
if(vdim(raw)!=30){{ERROR("raw length");}}
list saturation=sat_with_exp(raw,ideal(u));
ideal relative=saturation[1];
if(saturation[2]!=6){{ERROR("saturation exponent");}}
if(vdim(relative)!=18){{ERROR("relative length");}}
ideal endpoint0=sat(ideal(A,t),ideal(u));
ideal endpoint1=sat(ideal(A,t-1),ideal(u));
ideal interior=sat(
  sat(ideal(A,diff(Q,t)),ideal(u)),
  ideal(t*(1-t))
);
if(vdim(endpoint0)!=2){{ERROR("endpoint zero length");}}
if(vdim(endpoint1)!=2){{ERROR("endpoint one length");}}
if(vdim(interior)!=14){{ERROR("interior length");}}
print("PASS exact logarithmic Jacobian lengths 30 -> 18");
print("PASS saturation exponent 6");
print("PASS decomposition 2 + 2 + 14");
print(kbase(relative));
"""
    )
    for marker in (
        "PASS exact logarithmic Jacobian lengths 30 -> 18",
        "PASS saturation exponent 6",
        "PASS decomposition 2 + 2 + 14",
    ):
        assert marker in output
    observed_basis = {
        line.strip()
        .rstrip(",")
        .replace("^", "")
        .replace("*", "")
        .replace(" ", "")
        for line in output.splitlines()
        if line.strip().rstrip(",") in EXPECTED_BASIS
    }
    assert observed_basis == set(EXPECTED_BASIS), observed_basis
    return {
        "raw_length": 30,
        "u_saturation_exponent": 6,
        "relative_length": 18,
        "endpoint_t0_length": 2,
        "endpoint_t1_length": 2,
        "interior_length": 14,
        "standard_monomial_basis": list(EXPECTED_BASIS),
    }


def parse_probe(output: str) -> list[list[int]]:
    lines = output.strip().splitlines()
    assert lines[0] == "FOUND order=18 degree=18"
    coefficients = [[0] * (DEGREE + 1) for _ in range(ORDER + 1)]
    for line in lines[1:]:
        values = [int(value) for value in line.split()]
        assert len(values) == DEGREE + 2
        coefficients[values[0]] = values[1:]
    return coefficients


def polynomial_multiply(
    left: list[int],
    right: list[int],
    prime: int,
) -> list[int]:
    answer = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            answer[i + j] = (
                answer[i + j] + left_value * right_value
            ) % prime
    return answer


def polynomial_divmod(
    numerator: list[int],
    denominator: list[int],
    prime: int,
) -> tuple[list[int], list[int]]:
    remainder = [value % prime for value in numerator]
    quotient = [0] * (len(numerator) - len(denominator) + 1)
    inverse = pow(denominator[-1], -1, prime)
    for shift in range(
        len(numerator) - len(denominator),
        -1,
        -1,
    ):
        coefficient = remainder[shift + len(denominator) - 1] * inverse
        coefficient %= prime
        quotient[shift] = coefficient
        for index, value in enumerate(denominator):
            remainder[shift + index] = (
                remainder[shift + index] - coefficient * value
            ) % prime
    while remainder and remainder[-1] == 0:
        remainder.pop()
    return quotient, remainder


def common_forward_polynomial(prime: int) -> list[int]:
    answer = [1]
    for shift in COMMON_FORWARD_SHIFTS:
        answer = polynomial_multiply(
            answer,
            [shift % prime, 3],
            prime,
        )
    return answer


def leading_remainder_check(
    forward_coefficients: list[int],
    prime: int,
) -> int:
    q = q_expression(0, prime)
    terms = [
        (
            f"{coefficient}*u^{54 - 3 * shift}"
            f"*Q^{shift}"
        )
        for shift, coefficient in enumerate(forward_coefficients)
        if coefficient
    ]
    output = run_singular(
        f"""
LIB "elim.lib";
ring r={prime},(u,t),dp;
poly Q={q};
poly A=u*diff(Q,u)-3*Q;
poly C=t*(1-t)*diff(Q,t);
ideal relative=sat(ideal(A,C),ideal(u));
poly leading={"+".join(terms)};
poly remainder=reduce(leading,relative);
if(remainder==0){{ERROR("unexpected degree-17 certificate");}}
if(size(remainder)!=18){{ERROR("leading remainder support");}}
print("PASS leading coefficient has full 18-coordinate remainder");
"""
    )
    assert "PASS leading coefficient has full 18-coordinate remainder" in output
    return 18


def modular_recurrences() -> tuple[list[dict[str, object]], int]:
    compiler = shutil.which("g++")
    if compiler is None:
        raise RuntimeError("g++ is required")
    records = []
    point_zero_forward: list[int] | None = None
    with tempfile.TemporaryDirectory(prefix="sic33-relative-jacobian-") as path:
        executable = Path(path) / "recurrence-probe"
        subprocess.run(
            [
                compiler,
                "-O3",
                "-std=c++17",
                str(SOURCE),
                "-o",
                str(executable),
            ],
            check=True,
            timeout=30,
        )
        for point in range(2):
            for prime in PRIMES:
                completed = subprocess.run(
                    [
                        str(executable),
                        str(prime),
                        str(MAXIMUM_MOMENT),
                        str(ORDER),
                        str(DEGREE),
                        str(point),
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=120,
                )
                coefficients = parse_probe(completed.stdout)
                forward = coefficients[-1]
                fixed = common_forward_polynomial(prime)
                quotient, remainder = polynomial_divmod(
                    forward,
                    fixed,
                    prime,
                )
                assert not remainder
                assert len(quotient) == 11
                assert quotient[-1] != 0
                if point == 0 and prime == PRIMES[0]:
                    point_zero_forward = [
                        value % prime for value in forward
                    ]
                records.append(
                    {
                        "point": point,
                        "prime": prime,
                        "order": ORDER,
                        "m_degree": DEGREE,
                        "moments_computed": MAXIMUM_MOMENT + 1,
                        "fit_equations": (
                            (ORDER + 1) * (DEGREE + 1) - 1
                        ),
                        "unused_verification_equations": (
                            MAXIMUM_MOMENT
                            - ORDER
                            + 1
                            - (
                                (ORDER + 1)
                                * (DEGREE + 1)
                                - 1
                            )
                        ),
                        "forward_fixed_factor_degree": 8,
                        "forward_parameter_factor_degree": 10,
                    }
                )

        for point in range(2):
            completed = subprocess.run(
                [
                    str(executable),
                    str(PRIMES[0]),
                    str(MAXIMUM_MOMENT),
                    str(ORDER),
                    str(DEGREE - 1),
                    str(point),
                ],
                check=False,
                capture_output=True,
                text=True,
                timeout=120,
            )
            assert completed.stdout.strip() == "NONE"

    assert point_zero_forward is not None
    remainder_support = leading_remainder_check(
        point_zero_forward,
        PRIMES[0],
    )
    return records, remainder_support


def main() -> None:
    jacobian = exact_relative_jacobian()
    recurrence_records, remainder_support = modular_recurrences()
    artifact = {
        "format": "two-pair-sic-bidegree33-relative-jacobian-v1",
        "status": (
            "exact relative-Jacobian calculation plus modular recurrence "
            "evidence; not a universal telescoping certificate"
        ),
        "relative_jacobian": jacobian,
        "modular_order_18_probe": {
            "records": recurrence_records,
            "degree_17_status": (
                "no order-18 degree-17 fit at either point modulo 1000003"
            ),
            "common_forward_factor": (
                "product_(k in {44,46,47,49,50,52,53,55})(3m+k)"
            ),
            "remaining_forward_factor": (
                "a point-dependent monic decic after normalization"
            ),
            "leading_m_coefficient_remainder_support": remainder_support,
        },
        "interpretation": {
            "natural_relative_rank": 18,
            "order_27_role": (
                "higher-order low-m-degree desingularized tradeoff"
            ),
            "certificate_gate": (
                "a degree-17 polynomial certificate is impossible at the "
                "sample because its leading coefficient has nonzero image "
                "in every relative-Jacobian basis coordinate"
            ),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS exact logarithmic Jacobian length 18")
    print("PASS endpoint/interior decomposition 2+2+14")
    print("PASS six modular order-18 degree-18 recurrence fits")
    print("PASS 83 unused recurrence equations in every modular fit")
    print("PASS common degree-eight forward factor and variable decic")
    print("PASS degree-17 ansatz fails at both sampled points")
    print("PASS leading coefficient has full relative-Jacobian remainder")
    print("PASS result remains modular evidence, not a recurrence certificate")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
