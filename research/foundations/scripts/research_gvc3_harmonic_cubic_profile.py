#!/usr/bin/env python3
"""Exact obstruction and modular replay for the harmonic-cubic GVC(3) repair.

On the sphere rho=1 the degree-six family is

    F = alpha*E + H3*O,

where E=xy-2*t^2-x^2*t^2, O=y-3*x*t^2, and H3 is an arbitrary
Delta-harmonic cubic.  The first Reynolds moment removes the weight -1
coordinate of H3.  This script compiles every later moment from the seven
remaining weight-channel Laurent polynomials, without expanding F^m in
x,y,t, and sends selected finite-field ideals to Singular.

Modular Singular elimination supplies the discovery replay.  Exact msolve
unit bases over QQ exclude the nonzero-alpha projective cover, while exact
radical containments and phase weight handle the alpha=0 terminal planes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path
from typing import Iterator


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "gvc3_harmonic_cubic_profile_modular.json"
)

VARIABLES = ("p3", "p2", "p1", "p0", "n2", "n3")
N3_LOCALIZED_VARIABLES = ("z", "p3", "p2", "p1", "p0", "n2")
D_K_LOCALIZED_VARIABLES = ("p1", "n2", "z", "p3", "p2", "p0", "n3")
D_K_BOUNDARY_VARIABLES = ("p1", "p0", "p3", "p2", "n2", "n3")
ALPHA0_RADICAL_VARIABLES = ("z",) + VARIABLES
ALPHA0_RADICAL_GENERATORS = (
    "p1",
    "p0",
    "p3*n2",
    "p3*n3",
    "p2*n2",
    "p2*n3",
)

# Coefficients of alpha, h_3, h_2, h_1, h_0, h_-2, h_-3 after putting
# y=(1-t^2)/x on rho=1.  Keys are (phase exponent of x, exponent of t).
CHANNELS: tuple[dict[tuple[int, int], int], ...] = (
    {(2, 2): -1, (0, 2): -3, (0, 0): 1},
    {(4, 2): -3, (2, 0): 1, (2, 2): -1},
    {(3, 3): -3, (1, 3): -1, (1, 1): 1},
    {(0, 4): 5, (0, 2): -6, (2, 4): 15, (2, 2): -3, (0, 0): 1},
    {(1, 5): 15, (1, 3): -9, (-1, 5): 5, (-1, 3): -8, (-1, 1): 3},
    {
        (-1, 7): -3,
        (-1, 5): 6,
        (-1, 3): -3,
        (-3, 7): -1,
        (-3, 5): 3,
        (-3, 3): -3,
        (-3, 1): 1,
    },
    {
        (-2, 8): 3,
        (-2, 6): -9,
        (-2, 4): 9,
        (-2, 2): -3,
        (-4, 8): 1,
        (-4, 6): -4,
        (-4, 4): 6,
        (-4, 2): -4,
        (-4, 0): 1,
    },
)

ScalarPolynomial = dict[tuple[int, int], Fraction]
MomentPolynomial = dict[tuple[int, ...], Fraction]
ModularPolynomial = dict[tuple[int, ...], int]


def primitive_integer_polynomial(polynomial: MomentPolynomial) -> dict[tuple[int, ...], int]:
    denominator_lcm = 1
    for coefficient in polynomial.values():
        denominator_lcm = math.lcm(denominator_lcm, coefficient.denominator)
    integral = {
        exponent: int(coefficient * denominator_lcm)
        for exponent, coefficient in polynomial.items()
    }
    content = 0
    for coefficient in integral.values():
        content = math.gcd(content, abs(coefficient))
    if not content:
        return {}
    primitive = {
        exponent: coefficient // content
        for exponent, coefficient in integral.items()
        if coefficient
    }
    leading_coefficient = next(iter(sorted(primitive.items(), reverse=True)))[1]
    if leading_coefficient < 0:
        primitive = {exponent: -coefficient for exponent, coefficient in primitive.items()}
    return primitive


def multiply_scalar(left: ScalarPolynomial, right: ScalarPolynomial) -> ScalarPolynomial:
    answer: ScalarPolynomial = {}
    for (left_phase, left_height), left_coefficient in left.items():
        for (right_phase, right_height), right_coefficient in right.items():
            exponent = (left_phase + right_phase, left_height + right_height)
            answer[exponent] = (
                answer.get(exponent, Fraction(0))
                + left_coefficient * right_coefficient
            )
    return {exponent: coefficient for exponent, coefficient in answer.items() if coefficient}


def compositions(total: int, length: int, prefix: tuple[int, ...] = ()) -> Iterator[tuple[int, ...]]:
    if length == 1:
        yield prefix + (total,)
        return
    for entry in range(total + 1):
        yield from compositions(total - entry, length - 1, prefix + (entry,))


def reynolds_value(polynomial: ScalarPolynomial) -> Fraction:
    """Extract phase zero and apply normalized height integration."""
    return sum(
        (
            coefficient * Fraction(1, height + 1)
            for (phase, height), coefficient in polynomial.items()
            if phase == 0 and height % 2 == 0
        ),
        Fraction(0),
    )


def moment(order: int) -> MomentPolynomial:
    powers: list[list[ScalarPolynomial]] = []
    for channel in CHANNELS:
        channel_powers = [{(0, 0): Fraction(1)}]
        integral_channel = {
            exponent: Fraction(coefficient) for exponent, coefficient in channel.items()
        }
        for _ in range(order):
            channel_powers.append(
                multiply_scalar(channel_powers[-1], integral_channel)
            )
        powers.append(channel_powers)

    answer: MomentPolynomial = {}
    for occupation in compositions(order, len(CHANNELS)):
        scalar = {(0, 0): Fraction(1)}
        multinomial = math.factorial(order)
        for channel_index, count in enumerate(occupation):
            multinomial //= math.factorial(count)
            if count:
                scalar = multiply_scalar(scalar, powers[channel_index][count])
        value = reynolds_value(scalar)
        if value:
            answer[occupation] = multinomial * value
    return answer


def validate_compiler(compiled_moments: dict[int, MomentPolynomial]) -> dict[str, object]:
    """Exact structural checks independent of the finite-field eliminations."""
    second = primitive_integer_polynomial(specialize(compiled_moments[2], "alpha1"))
    expected_second = {
        (0, 0, 0, 1, 1, 0): 720,
        (0, 1, 0, 0, 1, 0): -1872,
        (0, 0, 1, 0, 0, 1): -2592,
        (1, 0, 0, 0, 0, 1): -4992,
        (0, 0, 0, 0, 0, 1): 936,
        (0, 0, 0, 2, 0, 0): -4524,
        (0, 1, 0, 1, 0, 0): 2288,
        (0, 0, 2, 0, 0, 0): 4576,
        (0, 0, 1, 0, 0, 0): 10296,
        (0, 0, 0, 0, 0, 0): 9009,
    }
    negative_second = {exponent: -coefficient for exponent, coefficient in second.items()}
    if second != expected_second and negative_second != expected_second:
        raise AssertionError("unexpected primitive second moment")

    terminal_checks: dict[str, bool] = {}
    # Positive terminal plane: only p3,p2 can be nonzero.  Negative terminal
    # plane: only n2,n3 can be nonzero.
    for name, active_indices in (
        ("positive", {0, 1}),
        ("negative", {4, 5}),
    ):
        vanishes = True
        for polynomial in compiled_moments.values():
            alpha0 = specialize(polynomial, "alpha0")
            if any(
                coefficient
                and all(
                    exponent[index] == 0
                    for index in range(len(VARIABLES))
                    if index not in active_indices
                )
                for exponent, coefficient in alpha0.items()
            ):
                vanishes = False
                break
        terminal_checks[name] = vanishes
    if not all(terminal_checks.values()):
        raise AssertionError("terminal plane does not annihilate every compiled moment")
    return {
        "second_moment_primitive": {
            "720*n2*p0-1872*n2*p2-2592*n3*p1-4992*n3*p3+936*n3"
            "-4524*p0^2+2288*p0*p2+4576*p1^2+10296*p1+9009": True
        },
        "terminal_planes_exact_through_compiled_cutoff": terminal_checks,
    }


def specialize(moment_polynomial: MomentPolynomial, chart: str) -> MomentPolynomial:
    if chart == "alpha1":
        return {
            exponent[1:]: coefficient
            for exponent, coefficient in moment_polynomial.items()
        }
    if chart == "alpha0":
        return {
            exponent[1:]: coefficient
            for exponent, coefficient in moment_polynomial.items()
            if exponent[0] == 0
        }
    raise ValueError(f"unknown chart: {chart}")


def singular_expression(polynomial: MomentPolynomial, prime: int) -> str:
    terms: list[str] = []
    for exponent, coefficient in sorted(polynomial.items(), reverse=True):
        modular_coefficient = (
            coefficient.numerator
            * pow(coefficient.denominator, -1, prime)
        ) % prime
        if not modular_coefficient:
            continue
        factors = []
        for variable, power in zip(VARIABLES, exponent, strict=True):
            if power == 1:
                factors.append(variable)
            elif power > 1:
                factors.append(f"{variable}^{power}")
        monomial = "*".join(factors)
        if not monomial:
            terms.append(str(modular_coefficient))
        elif modular_coefficient == 1:
            terms.append(monomial)
        else:
            terms.append(f"{modular_coefficient}*{monomial}")
    return "+".join(terms) or "0"


def add_modular(
    left: ModularPolynomial,
    right: ModularPolynomial,
    prime: int,
) -> ModularPolynomial:
    answer = dict(left)
    for exponent, coefficient in right.items():
        answer[exponent] = (answer.get(exponent, 0) + coefficient) % prime
        if answer[exponent] == 0:
            del answer[exponent]
    return answer


def multiply_modular(
    left: ModularPolynomial,
    right: ModularPolynomial,
    prime: int,
) -> ModularPolynomial:
    answer: ModularPolynomial = {}
    for left_exponent, left_coefficient in left.items():
        for right_exponent, right_coefficient in right.items():
            exponent = tuple(
                left_entry + right_entry
                for left_entry, right_entry in zip(
                    left_exponent, right_exponent, strict=True
                )
            )
            answer[exponent] = (
                answer.get(exponent, 0)
                + left_coefficient * right_coefficient
            ) % prime
    return {exponent: coefficient for exponent, coefficient in answer.items() if coefficient}


def power_modular(
    polynomial: ModularPolynomial,
    exponent: int,
    prime: int,
) -> ModularPolynomial:
    variables = len(next(iter(polynomial))) if polynomial else 5
    answer: ModularPolynomial = {(0,) * variables: 1}
    base = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = multiply_modular(answer, base, prime)
        remaining >>= 1
        if remaining:
            base = multiply_modular(base, base, prime)
    return answer


def modularize(polynomial: MomentPolynomial, prime: int) -> ModularPolynomial:
    return {
        exponent: (
            coefficient.numerator * pow(coefficient.denominator, -1, prime)
        ) % prime
        for exponent, coefficient in polynomial.items()
        if coefficient
    }


def n3_localized_equations(
    compiled_moments: dict[int, MomentPolynomial],
    cutoff: int,
    prime: int,
) -> list[ModularPolynomial]:
    """Use moment two to solve n3 on its nonzero-coefficient chart."""
    second = modularize(specialize(compiled_moments[2], "alpha1"), prime)
    coefficient: ModularPolynomial = {}
    remainder: ModularPolynomial = {}
    for exponent, value in second.items():
        assert exponent[-1] <= 1
        target = coefficient if exponent[-1] else remainder
        reduced_exponent = exponent[:-1]
        target[reduced_exponent] = value

    negative_remainder = {
        exponent: (-value) % prime for exponent, value in remainder.items()
    }
    maximum_power = cutoff
    coefficient_powers = [power_modular(coefficient, power, prime) for power in range(maximum_power + 1)]
    remainder_powers = [power_modular(negative_remainder, power, prime) for power in range(maximum_power + 1)]

    localized: list[ModularPolynomial] = []
    for order in range(3, cutoff + 1):
        original = modularize(specialize(compiled_moments[order], "alpha1"), prime)
        highest_n3 = max((exponent[-1] for exponent in original), default=0)
        pieces: list[ModularPolynomial] = [dict() for _ in range(highest_n3 + 1)]
        for exponent, value in original.items():
            pieces[exponent[-1]][exponent[:-1]] = value
        transformed: ModularPolynomial = {}
        for n3_power, piece in enumerate(pieces):
            term = multiply_modular(piece, remainder_powers[n3_power], prime)
            term = multiply_modular(
                term,
                coefficient_powers[highest_n3 - n3_power],
                prime,
            )
            transformed = add_modular(transformed, term, prime)
        localized.append({(0,) + exponent: value for exponent, value in transformed.items()})

    inverse_equation: ModularPolynomial = {
        (1,) + exponent: value for exponent, value in coefficient.items()
    }
    constant = (0,) * len(N3_LOCALIZED_VARIABLES)
    inverse_equation[constant] = (inverse_equation.get(constant, 0) - 1) % prime
    return [inverse_equation] + localized


def singular_expression_modular(
    polynomial: ModularPolynomial,
    variables: tuple[str, ...],
) -> str:
    terms: list[str] = []
    for exponent, coefficient in sorted(polynomial.items(), reverse=True):
        if not coefficient:
            continue
        factors = []
        for variable, power in zip(variables, exponent, strict=True):
            if power == 1:
                factors.append(variable)
            elif power > 1:
                factors.append(f"{variable}^{power}")
        monomial = "*".join(factors)
        if not monomial:
            terms.append(str(coefficient))
        elif coefficient == 1:
            terms.append(monomial)
        else:
            terms.append(f"{coefficient}*{monomial}")
    return "+".join(terms) or "0"


def singular_program(
    equations: list[MomentPolynomial],
    prime: int,
    extra_equations: tuple[str, ...] = (),
) -> str:
    expressions = [singular_expression(equation, prime) for equation in equations]
    expressions.extend(extra_equations)
    return f"""option(redSB);
ring R={prime},({','.join(VARIABLES)}),dp;
ideal I={','.join(expressions)};
timer=1;
ideal G=slimgb(I);
int elapsed=timer;
print(\"RESULT_BEGIN\");
print(\"basis_size=\"+string(size(G)));
print(\"dimension=\"+string(dim(G)));
print(\"elapsed_ticks=\"+string(elapsed));
if ((size(G)==1)&&(G[1]==1)) {{ print(\"unit=1\"); }}
else {{ print(\"unit=0\"); }}
print(\"BASIS_BEGIN\");
print(G);
print(\"BASIS_END\");
print(\"RESULT_END\");
quit;
"""


def parse_singular(stdout: str) -> dict[str, object]:
    result: dict[str, object] = {}
    for line in stdout.splitlines():
        stripped = line.strip()
        for key in ("basis_size", "dimension", "elapsed_ticks", "unit"):
            if stripped.startswith(key + "="):
                result[key] = int(stripped.split("=", 1)[1])
    if "BASIS_BEGIN" in stdout and "BASIS_END" in stdout:
        basis = stdout.split("BASIS_BEGIN", 1)[1].split("BASIS_END", 1)[0].strip()
        result["basis_sha256"] = hashlib.sha256(basis.encode()).hexdigest()
        result["basis_preview"] = basis[:4000]
    return result


def parse_radical_audit(stdout: str) -> dict[str, object]:
    results: list[dict[str, object]] = []
    for index, generator in enumerate(ALPHA0_RADICAL_GENERATORS):
        prefix = f"radical_{index}="
        value = None
        for line in stdout.splitlines():
            stripped = line.strip()
            if stripped.startswith(prefix):
                value = int(stripped.split("=", 1)[1])
                break
        results.append({"generator": generator, "unit": value})
    return {
        "radical_tests": results,
        "all_radical_tests_unit": int(all(item["unit"] == 1 for item in results)),
    }


def run_case(
    singular: str,
    chart: str,
    prime: int,
    cutoff: int,
    timeout: int,
    compiled_moments: dict[int, MomentPolynomial],
) -> dict[str, object]:
    radical_audit = chart == "alpha0_radical"
    if chart == "alpha1_n3":
        modular_equations = n3_localized_equations(compiled_moments, cutoff, prime)
        expressions = [
            singular_expression_modular(equation, N3_LOCALIZED_VARIABLES)
            for equation in modular_equations
        ]
        source = f"""option(redSB);
ring R={prime},({','.join(N3_LOCALIZED_VARIABLES)}),dp;
ideal I={','.join(expressions)};
timer=1;
ideal G=slimgb(I);
int elapsed=timer;
print(\"RESULT_BEGIN\");
print(\"basis_size=\"+string(size(G)));
print(\"dimension=\"+string(dim(G)));
print(\"elapsed_ticks=\"+string(elapsed));
if ((size(G)==1)&&(G[1]==1)) {{ print(\"unit=1\"); }}
else {{ print(\"unit=0\"); }}
print(\"BASIS_BEGIN\"); print(G); print(\"BASIS_END\");
print(\"RESULT_END\"); quit;
"""
        term_counts = [len(equation) for equation in modular_equations]
    elif chart in {"alpha1_d_k", "alpha1_dk"}:
        equations = [
            specialize(compiled_moments[order], "alpha1")
            for order in range(2, cutoff + 1)
        ]
        expressions = [singular_expression(equation, prime) for equation in equations]
        d_equation = "(-2592)*p1+(-4992)*p3+936"
        k_equation = "720*p0+(-1872)*p2"
        if chart == "alpha1_d_k":
            variables = D_K_LOCALIZED_VARIABLES
            expressions.extend((d_equation, f"z*({k_equation})-1"))
        else:
            variables = D_K_BOUNDARY_VARIABLES
            expressions.extend((d_equation, k_equation))
        source = f"""option(redSB);
ring R={prime},({','.join(variables)}),dp;
ideal I={','.join(expressions)};
timer=1; ideal G=slimgb(I); int elapsed=timer;
print(\"RESULT_BEGIN\");
print(\"basis_size=\"+string(size(G)));
print(\"dimension=\"+string(dim(G)));
print(\"elapsed_ticks=\"+string(elapsed));
if ((size(G)==1)&&(G[1]==1)) {{ print(\"unit=1\"); }}
else {{ print(\"unit=0\"); }}
print(\"BASIS_BEGIN\"); print(G); print(\"BASIS_END\");
print(\"RESULT_END\"); quit;
"""
        term_counts = [len(equation) for equation in equations] + [1, 2]
    elif radical_audit:
        equations = [
            specialize(compiled_moments[order], "alpha0")
            for order in range(2, cutoff + 1)
        ]
        expressions = [singular_expression(equation, prime) for equation in equations]
        tests: list[str] = []
        for index, generator in enumerate(ALPHA0_RADICAL_GENERATORS):
            tests.append(
                f"ideal J{index}=I,z*({generator})-1; "
                f"ideal G{index}=slimgb(J{index}); "
                f"if ((size(G{index})==1)&&(G{index}[1]==1)) "
                f'{{ print("radical_{index}=1"); }} '
                f'else {{ print("radical_{index}=0"); }}'
            )
        source = f"""option(redSB);
ring R={prime},({','.join(ALPHA0_RADICAL_VARIABLES)}),dp;
ideal I={','.join(expressions)};
timer=1;
{chr(10).join(tests)}
int elapsed=timer;
print("RESULT_BEGIN");
print("elapsed_ticks="+string(elapsed));
print("RESULT_END"); quit;
"""
        term_counts = [len(equation) for equation in equations]
    else:
        equations = [
            specialize(compiled_moments[order], chart)
            for order in range(2, cutoff + 1)
        ]
        source = singular_program(equations, prime)
        term_counts = [len(equation) for equation in equations]
    record: dict[str, object] = {
        "chart": chart,
        "prime": prime,
        "cutoff": cutoff,
        "term_counts": term_counts,
        "input_sha256": hashlib.sha256(source.encode()).hexdigest(),
    }
    try:
        started = time.monotonic()
        completed = subprocess.run(
            [singular, "-q"],
            input=source,
            text=True,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
        record["seconds"] = time.monotonic() - started
        record["returncode"] = completed.returncode
        record["status"] = "completed" if completed.returncode == 0 else "failed"
        if radical_audit:
            record.update(parse_radical_audit(completed.stdout))
        else:
            record.update(parse_singular(completed.stdout))
        record["stdout_tail"] = completed.stdout[-4000:]
        record["stderr"] = completed.stderr[-4000:]
    except subprocess.TimeoutExpired as error:
        record["seconds"] = time.monotonic() - started
        record["status"] = "timeout"
        record["stdout_tail"] = (error.stdout or "")[-4000:]
        record["stderr"] = (error.stderr or "")[-4000:]
    return record


def exact_alpha0_radical_audit(
    compiled_moments: dict[int, MomentPolynomial],
    timeout: int,
    threads: int,
) -> dict[str, object]:
    """Prove the six alpha=0 radical containments over QQ with msolve."""
    import sympy as sp

    from jcsearch import msolve

    z, p3, p2, p1, p0, n2, n3 = sp.symbols("z p3 p2 p1 p0 n2 n3")
    variables = (p3, p2, p1, p0, n2, n3)
    solver_variables = (z,) + variables

    def expression(polynomial: MomentPolynomial) -> sp.Expr:
        return sum(
            (
                sp.Rational(coefficient.numerator, coefficient.denominator)
                * sp.prod(
                    variable**power
                    for variable, power in zip(variables, exponent, strict=True)
                )
                for exponent, coefficient in polynomial.items()
            ),
            sp.Integer(0),
        )

    equations = [
        expression(specialize(compiled_moments[order], "alpha0"))
        for order in range(2, 8)
    ]
    generators = (p1, p0, p3 * n2, p3 * n3, p2 * n2, p2 * n3)
    tests: list[dict[str, object]] = []
    for generator in generators:
        started = time.monotonic()
        result = msolve.run(
            equations + [z * generator - 1],
            solver_variables,
            prime=0,
            threads=threads,
            groebner=True,
            timeout=timeout,
        )
        elapsed = time.monotonic() - started
        exact_unit = (
            result.returncode == 0
            and result.empty
            and "#field characteristic: 0" in result.output
            and "#length of basis:      1 element" in result.output
            and result.output.rstrip().endswith("[1]:")
        )
        tests.append(
            {
                "generator": str(generator).replace("**", "^"),
                "unit": int(exact_unit),
                "seconds": elapsed,
                "returncode": result.returncode,
                "output_sha256": hashlib.sha256(result.output.encode()).hexdigest(),
                "output": result.output[-1000:],
                "stderr": result.stderr[-1000:],
            }
        )
    return {
        "backend": "msolve exact characteristic-zero F4",
        "moments": [2, 3, 4, 5, 6, 7],
        "radical_generators": list(ALPHA0_RADICAL_GENERATORS),
        "tests": tests,
        "all_unit": int(all(test["unit"] == 1 for test in tests)),
        "conclusion": (
            "radical(I_2,...,I_7)=(p1,p0,p3*n2,p3*n3,p2*n2,p2*n3)"
        ),
    }


def exact_alpha1_projective_audit(
    compiled_moments: dict[int, MomentPolynomial],
    timeout: int,
    threads: int,
) -> dict[str, object]:
    """Certify the exhaustive alpha=1 pivot cover over QQ with msolve."""
    import sympy as sp

    from jcsearch import msolve

    z, p3, p2, p1, p0, n2, n3 = sp.symbols("z p3 p2 p1 p0 n2 n3")
    variables = (p3, p2, p1, p0, n2, n3)

    def expression(polynomial: MomentPolynomial) -> sp.Expr:
        return sum(
            (
                sp.Rational(coefficient.numerator, coefficient.denominator)
                * sp.prod(
                    variable**power
                    for variable, power in zip(variables, exponent, strict=True)
                )
                for exponent, coefficient in polynomial.items()
            ),
            sp.Integer(0),
        )

    moments = {
        order: expression(specialize(compiled_moments[order], "alpha1"))
        for order in range(2, 9)
    }
    d_pivot = 936 - 2592 * p1 - 4992 * p3
    k_pivot = 720 * p0 - 1872 * p2
    cases = (
        (
            "D_nonzero",
            [moments[order] for order in range(2, 9)] + [z * d_pivot - 1],
            (z,) + variables,
            [2, 3, 4, 5, 6, 7, 8],
        ),
        (
            "D_zero_K_nonzero",
            [moments[order] for order in range(2, 8)]
            + [d_pivot, z * k_pivot - 1],
            (z,) + variables,
            [2, 3, 4, 5, 6, 7],
        ),
        (
            "D_zero_K_zero",
            [moments[order] for order in range(2, 8)] + [d_pivot, k_pivot],
            variables,
            [2, 3, 4, 5, 6, 7],
        ),
    )
    tests: list[dict[str, object]] = []
    for name, equations, solver_variables, orders in cases:
        started = time.monotonic()
        result = msolve.run(
            equations,
            solver_variables,
            prime=0,
            threads=threads,
            groebner=True,
            timeout=timeout,
        )
        elapsed = time.monotonic() - started
        exact_unit = (
            result.returncode == 0
            and result.empty
            and "#field characteristic: 0" in result.output
            and "#length of basis:      1 element" in result.output
            and result.output.rstrip().endswith("[1]:")
        )
        tests.append(
            {
                "chart": name,
                "moments": orders,
                "unit": int(exact_unit),
                "seconds": elapsed,
                "returncode": result.returncode,
                "output_sha256": hashlib.sha256(result.output.encode()).hexdigest(),
                "output": result.output[-1000:],
                "stderr": result.stderr[-1000:],
            }
        )
    return {
        "backend": "msolve exact characteristic-zero F4",
        "pivots": {
            "D": "936-2592*p1-4992*p3",
            "K": "720*p0-1872*p2",
        },
        "tests": tests,
        "all_unit": int(all(test["unit"] == 1 for test in tests)),
        "conclusion": "the alpha=1 pure-moment scheme through moment 8 is empty",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--singular", default="Singular")
    parser.add_argument(
        "--charts",
        nargs="+",
        choices=(
            "alpha1",
            "alpha1_n3",
            "alpha1_d_k",
            "alpha1_dk",
            "alpha0",
            "alpha0_radical",
        ),
        default=["alpha1_n3"],
    )
    parser.add_argument(
        "--cases",
        nargs="+",
        metavar="CHART:CUTOFF",
        help="run chart-specific cutoffs instead of the charts/cutoffs product",
    )
    parser.add_argument("--primes", nargs="+", type=int, default=[101, 103, 107])
    parser.add_argument("--cutoffs", nargs="+", type=int, default=[6, 7])
    parser.add_argument("--timeout", type=int, default=1200)
    parser.add_argument(
        "--exact-boundary",
        action="store_true",
        help="also certify the alpha=0 radical over QQ with msolve",
    )
    parser.add_argument(
        "--exact-all",
        action="store_true",
        help="certify both alpha=0 and the full alpha=1 projective cover over QQ",
    )
    parser.add_argument(
        "--augment-exact-boundary",
        action="store_true",
        help="add the exact alpha=0 audit to an existing output artifact",
    )
    parser.add_argument(
        "--augment-exact-all",
        action="store_true",
        help="add every exact QQ audit to an existing output artifact",
    )
    parser.add_argument("--msolve-threads", type=int, default=4)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()

    if arguments.augment_exact_boundary or arguments.augment_exact_all:
        if not arguments.output.is_file():
            raise SystemExit(f"missing artifact to augment: {arguments.output}")
        artifact = json.loads(arguments.output.read_text())
        if not str(artifact.get("format", "")).startswith(
            "gvc3-harmonic-cubic-profile-"
        ):
            raise SystemExit("output artifact has the wrong format")
        compiled_moments = {order: moment(order) for order in range(2, 8)}
        if arguments.augment_exact_all:
            compiled_moments[8] = moment(8)
        validate_compiler(compiled_moments)
        exact_boundary = exact_alpha0_radical_audit(
            compiled_moments,
            arguments.timeout,
            arguments.msolve_threads,
        )
        artifact["format"] = "gvc3-harmonic-cubic-profile-modular-v2"
        artifact["status"] = (
            "bounded modular discovery on alpha!=0; exact characteristic-zero "
            "radical and terminal proof on alpha=0"
        )
        artifact["exact_characteristic_zero_boundary"] = exact_boundary
        if arguments.augment_exact_all:
            exact_projective = exact_alpha1_projective_audit(
                compiled_moments,
                arguments.timeout,
                arguments.msolve_threads,
            )
            artifact["exact_characteristic_zero_projective_cover"] = exact_projective
            artifact["format"] = "gvc3-harmonic-cubic-profile-v3"
            artifact["status"] = (
                "exact characteristic-zero obstruction; modular cases are "
                "discovery replay"
            )
        else:
            exact_projective = None
        arguments.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
        print("DONE exact alpha0 radical unit=", exact_boundary["all_unit"], flush=True)
        if exact_boundary["all_unit"] != 1:
            raise RuntimeError("exact alpha=0 radical audit did not certify every generator")
        if exact_projective is not None and exact_projective["all_unit"] != 1:
            raise RuntimeError("exact alpha=1 projective audit did not certify every chart")
        return

    selected_cases: list[tuple[str, int]]
    if arguments.cases:
        valid_charts = {
            "alpha1",
            "alpha1_n3",
            "alpha1_d_k",
            "alpha1_dk",
            "alpha0",
            "alpha0_radical",
        }
        selected_cases = []
        for specification in arguments.cases:
            try:
                chart, cutoff_text = specification.rsplit(":", 1)
                cutoff = int(cutoff_text)
            except ValueError as error:
                raise SystemExit(f"invalid case {specification!r}; expected CHART:CUTOFF") from error
            if chart not in valid_charts or cutoff < 2:
                raise SystemExit(f"invalid case {specification!r}")
            selected_cases.append((chart, cutoff))
    else:
        selected_cases = [
            (chart, cutoff)
            for chart in arguments.charts
            for cutoff in arguments.cutoffs
        ]

    maximum_cutoff = max(cutoff for _, cutoff in selected_cases)
    if arguments.exact_boundary or arguments.exact_all:
        maximum_cutoff = max(maximum_cutoff, 7)
    if arguments.exact_all:
        maximum_cutoff = max(maximum_cutoff, 8)
    compiled_moments = {order: moment(order) for order in range(2, maximum_cutoff + 1)}
    exact_checks = validate_compiler(compiled_moments)
    artifact: dict[str, object] = {
        "format": "gvc3-harmonic-cubic-profile-modular-v2",
        "status": (
            "bounded modular discovery on alpha!=0; alpha=0 is certified exactly "
            "when exact_characteristic_zero_boundary is present"
        ),
        "family": "F=alpha*E+H3*O with H3 an arbitrary harmonic cubic",
        "first_moment": "16*h_-1/35",
        "elimination": "h_-1=0",
        "harmonic_basis": [
            "x^3",
            "x^2*t",
            "x^2*y-4*x*t^2",
            "3*x*y*t-2*t^3",
            "x*y^2-4*y*t^2",
            "y^2*t",
            "y^3",
        ],
        "moment_term_counts": {
            str(order): len(compiled_moments[order])
            for order in compiled_moments
        },
        "exact_structural_checks": exact_checks,
        "cases": [],
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    for chart, cutoff in selected_cases:
        for prime in arguments.primes:
            print(f"RUN chart={chart} p={prime} cutoff={cutoff}", flush=True)
            result = run_case(
                arguments.singular,
                chart,
                prime,
                cutoff,
                arguments.timeout,
                compiled_moments,
            )
            artifact["cases"].append(result)
            arguments.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
            print(
                "DONE",
                result.get("status"),
                "unit=", result.get("unit"),
                "dim=", result.get("dimension"),
                "basis=", result.get("basis_size"),
                "radical=", result.get("all_radical_tests_unit"),
                "seconds=", round(float(result.get("seconds", 0)), 3),
                flush=True,
            )

    if arguments.exact_boundary or arguments.exact_all:
        print("RUN exact alpha0 radical over QQ", flush=True)
        exact_boundary = exact_alpha0_radical_audit(
            compiled_moments,
            arguments.timeout,
            arguments.msolve_threads,
        )
        artifact["exact_characteristic_zero_boundary"] = exact_boundary
        arguments.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
        print(
            "DONE exact alpha0 radical unit=",
            exact_boundary["all_unit"],
            flush=True,
        )
        if exact_boundary["all_unit"] != 1:
            raise RuntimeError("exact alpha=0 radical audit did not certify every generator")
    if arguments.exact_all:
        print("RUN exact alpha1 projective cover over QQ", flush=True)
        exact_projective = exact_alpha1_projective_audit(
            compiled_moments,
            arguments.timeout,
            arguments.msolve_threads,
        )
        artifact["exact_characteristic_zero_projective_cover"] = exact_projective
        artifact["format"] = "gvc3-harmonic-cubic-profile-v3"
        artifact["status"] = (
            "exact characteristic-zero obstruction; modular cases are discovery replay"
        )
        arguments.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
        print(
            "DONE exact alpha1 projective cover unit=",
            exact_projective["all_unit"],
            flush=True,
        )
        if exact_projective["all_unit"] != 1:
            raise RuntimeError("exact alpha=1 projective audit did not certify every chart")


if __name__ == "__main__":
    main()
