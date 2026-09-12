#!/usr/bin/env python3
"""Exact and bounded audits for a three-variable GVC lift.

The script has four logically separate parts.

1. It replays the coordinate-only multiplier z2 for the rank-five SIC(2)
   witness through the declared cutoff.
2. It checks the canonical five-channel auxiliary-exponent lift over
   GF(101); its full-support chart is empty through moment four.
3. It proves by an exact rational Groebner basis that the minimal tagged
   rank-two lift

       Lambda = d_t d_z + B_3(d_t,d_y),
       P      = (t-y) z + t^2 (t-y)

   cannot have its first five pure moments zero, for arbitrary binary cubic
   B_3.
4. Over GF(101), it checks that the complete factor-compatible profile

       P = (t-y) (z + t^2 + q1*t*y + q2*y^2)

   has a one-dimensional pure-zero fiber through moment five and empty
   normalized affine fiber through moment six.
5. It performs deterministic exact GF(101) fiber solves for 200 general
   normalized cubic profiles.  This last part is a bounded discovery search,
   not a characteristic-zero theorem.

Run with the repository virtual environment because SymPy is required.
Singular is required for the modular Groebner calculations.
"""

from __future__ import annotations

import importlib.util
import json
import random
import subprocess
from math import comb, factorial
from pathlib import Path

import sympy as sp
from sympy.polys.domains import GF
from sympy.polys.rings import ring


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "three_variable_gvc_tagged_lift.json"
)
SIC_CHECKER = ROOT / "scripts" / "verify_two_pair_image_mathieu_counterexample.py"
COORDINATE_MULTIPLIER_CUTOFF = 6
GENERAL_PROFILE_SAMPLES = 200
MODULUS = 101
RANDOM_SEED = 20260730


def load_sic_checker():
    spec = importlib.util.spec_from_file_location("sic2_checker", SIC_CHECKER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {SIC_CHECKER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def coordinate_multiplier_replay() -> list[int]:
    """Check E_2(z2 F^m)=known_m*z1 through the declared cutoff."""
    module = load_sic_checker()
    f, _, _ = module.witness()
    z2 = module.monomial((0, 0, 0, 1))
    checked: list[int] = []
    for order in range(1, COORDINATE_MULTIPLIER_CUTOFF + 1):
        contracted = module.contraction(
            module.multiply(z2, module.power(f, order))
        )
        expected_coefficient = (
            factorial(4 * order + 2)
            * factorial(order)
            // module.double_factorial_odd(2 * order + 1)
        )
        expected = {(0, 0, 1, 0): expected_coefficient}
        assert contracted == expected
        checked.append(order)
    return checked


def rank_five_auxiliary_ladder() -> dict[str, object]:
    """Check the canonical rank-five exponent ladder through moment four.

    Decompose the SIC tensor by rows M[k,*].  Put the binary dual monomial
    in channel i at auxiliary exponent 4-i, and put row k of M at the same
    coordinate exponent.  All ten channel scalars are free; the two global
    scalings normalize a0=b0=1.  A unit Groebner basis over GF(101) excludes
    the complete full-support chart through moment four.
    """
    a_variables = sp.symbols("a1:5")
    b_variables = sp.symbols("b1:5")
    parameters = (*a_variables, *b_variables)
    a_coefficients = (sp.Integer(1), *a_variables)
    b_coefficients = (sp.Integer(1), *b_variables)
    matrix = (
        (-1, 2, 0, 0, 0),
        (sp.Rational(-3, 2), 2, 6, 0, 0),
        (sp.Rational(-1, 2), sp.Rational(3, 2), 6, 6, 0),
        (0, 1, sp.Rational(3, 2), 2, 2),
        (0, 0, sp.Rational(-1, 2), sp.Rational(-3, 2), -1),
    )
    auxiliary_exponents = (4, 3, 2, 1, 0)
    operator = {
        (index, 4 - index, auxiliary_exponents[index]):
        a_coefficients[index]
        for index in range(5)
    }
    polynomial: dict[tuple[int, int, int], object] = {}
    for channel in range(5):
        for coordinate_t in range(5):
            coefficient = matrix[channel][coordinate_t]
            if coefficient:
                exponent = (
                    coordinate_t,
                    4 - coordinate_t,
                    auxiliary_exponents[channel],
                )
                polynomial[exponent] = (
                    polynomial.get(exponent, sp.Integer(0))
                    + b_coefficients[channel] * coefficient
                )

    def multiply_ternary(left, right):
        result = {}
        for left_exponent, left_coefficient in left.items():
            for right_exponent, right_coefficient in right.items():
                exponent = tuple(
                    left_exponent[index] + right_exponent[index]
                    for index in range(3)
                )
                result[exponent] = (
                    result.get(exponent, sp.Integer(0))
                    + left_coefficient * right_coefficient
                )
        return result

    def contract_ternary(operator_power, polynomial_power):
        result = {}
        for derivative, operator_coefficient in operator_power.items():
            for exponent, polynomial_coefficient in polynomial_power.items():
                if any(
                    derivative[index] > exponent[index]
                    for index in range(3)
                ):
                    continue
                residual = tuple(
                    exponent[index] - derivative[index]
                    for index in range(3)
                )
                falling_factorial = 1
                for index in range(3):
                    falling_factorial *= (
                        factorial(exponent[index])
                        // factorial(residual[index])
                    )
                result[residual] = (
                    result.get(residual, sp.Integer(0))
                    + operator_coefficient
                    * polynomial_coefficient
                    * falling_factorial
                )
        return result

    operator_power = {(0, 0, 0): sp.Integer(1)}
    polynomial_power = {(0, 0, 0): sp.Integer(1)}
    equations = []
    output_term_counts = []
    for _order in range(1, 5):
        operator_power = multiply_ternary(operator_power, operator)
        polynomial_power = multiply_ternary(polynomial_power, polynomial)
        contraction = contract_ternary(operator_power, polynomial_power)
        output_term_counts.append(len(contraction))
        for coefficient in contraction.values():
            numerator = sp.together(coefficient).as_numer_denom()[0]
            equations.append(
                singular_expression(
                    sp.Poly(
                        numerator,
                        *parameters,
                        modulus=MODULUS,
                    ).as_expr()
                )
            )
    summary = singular_summary(
        tuple(str(parameter) for parameter in parameters),
        equations,
    )
    assert summary == (1, -1, "0")
    return {
        "field": f"GF({MODULUS})",
        "operator_auxiliary_exponents": list(auxiliary_exponents),
        "polynomial_auxiliary_exponents": list(auxiliary_exponents),
        "normalization": "a0=b0=1",
        "orders": [1, 2, 3, 4],
        "output_term_counts": output_term_counts,
        "groebner_summary": {
            "basis_size": summary[0],
            "dimension": summary[1],
            "normal_form_of_1": summary[2],
        },
        "status": "bounded exact modular full-support obstruction",
    }


BinaryPolynomial = dict[int, object]


def binary_multiply(
    left: BinaryPolynomial,
    right: BinaryPolynomial,
    zero,
) -> BinaryPolynomial:
    """Multiply binary homogeneous forms, indexed by the t-exponent."""
    result: BinaryPolynomial = {}
    for left_t, left_coefficient in left.items():
        for right_t, right_coefficient in right.items():
            exponent = left_t + right_t
            result[exponent] = (
                result.get(exponent, zero)
                + left_coefficient * right_coefficient
            )
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def tagged_moment(
    order: int,
    b_symbol: BinaryPolynomial,
    cubic_profile: BinaryPolynomial,
    linear_profile: BinaryPolynomial,
    one,
    zero,
) -> object:
    """Return Lambda^m(P^m) using the exact diagonal channel formula.

    Here Lambda=d_t*d_z+B(d_t,d_y), P=z*L+C, deg(B)=deg(C)=3,
    and deg(L)=1.  The z-degree and binary total-degree inequalities force
    the same channel count k in the two binomial expansions.
    """
    linear_powers: list[BinaryPolynomial] = [{0: one}]
    for _ in range(order):
        linear_powers.append(
            binary_multiply(
                linear_powers[-1],
                linear_profile,
                zero,
            )
        )

    b_power: BinaryPolynomial = {0: one}
    c_power: BinaryPolynomial = {0: one}
    result = zero
    for channel_count in range(order + 1):
        if channel_count:
            b_power = binary_multiply(b_power, b_symbol, zero)
            c_power = binary_multiply(c_power, cubic_profile, zero)
        polynomial = binary_multiply(
            linear_powers[order - channel_count],
            c_power,
            zero,
        )
        pairing = zero
        for derivative_t, operator_coefficient in b_power.items():
            target_t = derivative_t + order - channel_count
            target_y = 3 * channel_count - derivative_t
            polynomial_coefficient = polynomial.get(target_t, zero)
            if polynomial_coefficient:
                pairing += (
                    operator_coefficient
                    * polynomial_coefficient
                    * factorial(target_t)
                    * factorial(target_y)
                )
        result += (
            comb(order, channel_count) ** 2
            * factorial(order - channel_count)
            * pairing
        )
    return result


def exact_minimal_tag_obstruction() -> dict[str, object]:
    """Exact Q-Groebner obstruction for C=t^2(t-y), moments 1..5."""
    a0, a1, a2, a3 = sp.symbols("a0 a1 a2 a3")
    parameters = (a0, a1, a2, a3)
    b_symbol = {3: a0, 2: a1, 1: a2, 0: a3}
    linear_profile = {1: sp.Integer(1), 0: sp.Integer(-1)}
    # C=t^2(t-y)=t^3-t^2*y.
    cubic_profile = {3: sp.Integer(1), 2: sp.Integer(-1)}
    equations = []
    groebner_sizes = []
    for order in range(1, 6):
        equation = sp.Poly(
            sp.expand(
                tagged_moment(
                    order,
                    b_symbol,
                    cubic_profile,
                    linear_profile,
                    sp.Integer(1),
                    sp.Integer(0),
                )
            ),
            *parameters,
            domain=sp.QQ,
        ).primitive()[1].as_expr()
        equations.append(equation)
        basis = sp.groebner(equations, *parameters, order="grevlex")
        groebner_sizes.append(len(basis.polys))
        if order <= 4:
            assert not (
                len(basis.polys) == 1
                and basis.polys[0].as_expr() == 1
            )
        else:
            assert len(basis.polys) == 1
            assert basis.polys[0].as_expr() == 1
    return {
        "field": "Q",
        "orders": list(range(1, 6)),
        "groebner_basis_sizes": groebner_sizes,
        "terminal_order": 5,
        "terminal_ideal": "(1)",
    }


def singular_summary(
    variable_names: tuple[str, ...],
    equations: list[str],
) -> tuple[int, int, str]:
    """Return size, dimension, and normal form of 1 for a modular ideal."""
    script = (
        f"ring r={MODULUS},({','.join(variable_names)}),dp;\n"
        "option(redSB);\n"
        f"ideal I={','.join(equations) if equations else '0'};\n"
        "ideal G=std(I);\n"
        "print(size(G));\n"
        "print(dim(G));\n"
        "reduce(1,G);\n"
        "quit;\n"
    )
    completed = subprocess.run(
        ["Singular", "-q"],
        input=script,
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    values = [
        line.strip()
        for line in completed.stdout.splitlines()
        if line.strip()
    ]
    if len(values) < 3:
        raise RuntimeError(
            "unexpected Singular output:\n"
            + completed.stdout
            + completed.stderr
        )
    return int(values[-3]), int(values[-2]), values[-1]


def singular_expression(polynomial) -> str:
    return str(polynomial).replace("**", "^")


def factor_compatible_modular_obstruction() -> dict[str, object]:
    """Check the normalized C=(t-y)(t^2+q1*t*y+q2*y^2) chart."""
    (
        coefficient_ring,
        a0,
        a1,
        a2,
        a3,
        q1,
        q2,
    ) = ring("a0,a1,a2,a3,q1,q2", GF(MODULUS))
    b_symbol = {3: a0, 2: a1, 1: a2, 0: a3}
    quadratic_profile = {2: coefficient_ring.one, 1: q1, 0: q2}
    linear_profile = {1: coefficient_ring.one, 0: -coefficient_ring.one}
    equations = []
    summaries = {}
    linear_power: BinaryPolynomial = {0: coefficient_ring.one}
    for order in range(1, 7):
        linear_power = binary_multiply(
            linear_power,
            linear_profile,
            coefficient_ring.zero,
        )
        # L^(m-k) C^k = L^m Q^k when C=LQ.
        b_power: BinaryPolynomial = {0: coefficient_ring.one}
        q_power: BinaryPolynomial = {0: coefficient_ring.one}
        value = coefficient_ring.zero
        for channel_count in range(order + 1):
            if channel_count:
                b_power = binary_multiply(
                    b_power,
                    b_symbol,
                    coefficient_ring.zero,
                )
                q_power = binary_multiply(
                    q_power,
                    quadratic_profile,
                    coefficient_ring.zero,
                )
            polynomial = binary_multiply(
                linear_power,
                q_power,
                coefficient_ring.zero,
            )
            pairing = coefficient_ring.zero
            for derivative_t, operator_coefficient in b_power.items():
                target_t = derivative_t + order - channel_count
                target_y = 3 * channel_count - derivative_t
                polynomial_coefficient = polynomial.get(
                    target_t,
                    coefficient_ring.zero,
                )
                if polynomial_coefficient:
                    pairing += (
                        operator_coefficient
                        * polynomial_coefficient
                        * (
                            factorial(target_t)
                            * factorial(target_y)
                            % MODULUS
                        )
                    )
            value += (
                comb(order, channel_count) ** 2
                * factorial(order - channel_count)
                * pairing
            )
        equations.append(singular_expression(value))
        if order in (5, 6):
            summaries[str(order)] = singular_summary(
                ("a0", "a1", "a2", "a3", "q1", "q2"),
                equations,
            )
    assert summaries["5"][1] == 1
    assert summaries["5"][2] == "1"
    assert summaries["6"] == (1, -1, "0")
    return {
        "field": f"GF({MODULUS})",
        "normalization": "coefficient of t^2 in Q is 1",
        "through_moment_5": {
            "basis_size": summaries["5"][0],
            "dimension": summaries["5"][1],
            "normal_form_of_1": summaries["5"][2],
        },
        "through_moment_6": {
            "basis_size": summaries["6"][0],
            "dimension": summaries["6"][1],
            "normal_form_of_1": summaries["6"][2],
        },
    }


def specialize_parameter_polynomial(polynomial, values: tuple[int, int, int]) -> str:
    """Specialize c1,c2,c3 and return a Singular polynomial in a0..a3."""
    terms: list[str] = []
    for monomial, coefficient in polynomial.terms():
        scalar = int(coefficient) % MODULUS
        for index, value in zip((4, 5, 6), values):
            scalar = (
                scalar
                * pow(value, monomial[index], MODULUS)
                % MODULUS
            )
        if scalar == 0:
            continue
        term = str(scalar)
        for index, name in enumerate(("a0", "a1", "a2", "a3")):
            exponent = monomial[index]
            if exponent:
                if exponent == 1:
                    term += f"*{name}"
                else:
                    term += f"*{name}^{exponent}"
        terms.append(term)
    return "+".join(terms) or "0"


def general_cubic_profile_fibers() -> dict[str, object]:
    """Solve 200 normalized cubic-profile fibers exactly over GF(101)."""
    (
        coefficient_ring,
        a0,
        a1,
        a2,
        a3,
        c1,
        c2,
        c3,
    ) = ring("a0,a1,a2,a3,c1,c2,c3", GF(MODULUS))
    b_symbol = {3: a0, 2: a1, 1: a2, 0: a3}
    cubic_profile = {3: coefficient_ring.one, 2: c1, 1: c2, 0: c3}
    linear_profile = {1: coefficient_ring.one, 0: -coefficient_ring.one}
    equations = [
        tagged_moment(
            order,
            b_symbol,
            cubic_profile,
            linear_profile,
            coefficient_ring.one,
            coefficient_ring.zero,
        )
        for order in range(1, 7)
    ]

    random_generator = random.Random(RANDOM_SEED)
    profiles = [
        (0, 0, 0),
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 1, 1),
        (MODULUS - 1, MODULUS - 1, MODULUS - 1),
    ]
    profiles.extend(
        tuple(random_generator.randrange(MODULUS) for _ in range(3))
        for _ in range(GENERAL_PROFILE_SAMPLES - len(profiles))
    )

    first_empty_histogram: dict[int, int] = {}
    exceptional_profiles: dict[int, list[list[int]]] = {}
    survivors = []
    for profile in profiles:
        specialized = [
            specialize_parameter_polynomial(equation, profile)
            for equation in equations
        ]
        first_empty = None
        for cutoff in (4, 5, 6):
            summary = singular_summary(
                ("a0", "a1", "a2", "a3"),
                specialized[:cutoff],
            )
            if summary == (1, -1, "0"):
                first_empty = cutoff
                break
        if first_empty is None:
            survivors.append(list(profile))
        else:
            first_empty_histogram[first_empty] = (
                first_empty_histogram.get(first_empty, 0) + 1
            )
            if first_empty > 5:
                exceptional_profiles.setdefault(first_empty, []).append(
                    list(profile)
                )

    assert first_empty_histogram == {5: 198, 6: 2}
    assert survivors == []
    return {
        "field": f"GF({MODULUS})",
        "normalization": "coefficient of t^3 in C is 1",
        "sample_count": len(profiles),
        "random_seed": RANDOM_SEED,
        "first_empty_moment_histogram": {
            str(key): value
            for key, value in sorted(first_empty_histogram.items())
        },
        "profiles_first_empty_at_6": exceptional_profiles.get(6, []),
        "survivors_through_6": survivors,
        "status": "bounded exact fiber search, not a global theorem",
    }


def main() -> None:
    coordinate_orders = coordinate_multiplier_replay()
    auxiliary_ladder = rank_five_auxiliary_ladder()
    exact_tag = exact_minimal_tag_obstruction()
    factor_profile = factor_compatible_modular_obstruction()
    general_profiles = general_cubic_profile_fibers()

    artifact = {
        "format": "three-variable-gvc-tagged-lift-v1",
        "coordinate_only_sic_multiplier": {
            "identity": (
                "E_2(z2*F^m)=((4m+2)!*m!/(2m+1)!!)*z1"
            ),
            "bounded_replay_orders": coordinate_orders,
            "all_order_proof": (
                "E(xi1*H)=d_z1 E(H), while the xi2*z2 derivative "
                "vanishes by odd Hopf height"
            ),
        },
        "rank_five_auxiliary_ladder": auxiliary_ladder,
        "tagged_channel_formula": (
            "sum_k binom(m,k)^2 (m-k)! "
            "d_t^(m-k) B^k(L^(m-k) C^k)"
        ),
        "exact_minimal_tag_obstruction": exact_tag,
        "factor_compatible_modular_obstruction": factor_profile,
        "general_cubic_profile_fibers": general_profiles,
        "conclusion": {
            "gvc3_counterexample_found": False,
            "closed_exactly": (
                "the fixed Long-tag profile C=t^2(t-y), arbitrary "
                "binary cubic B, through pure moment five over Q"
            ),
            "remaining": (
                "general cubic profile in characteristic zero; higher "
                "order/degree tagged channels; non-factor-compatible lifts"
            ),
        },
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS coordinate-only SIC(2) multiplier replay")
    print("PASS canonical rank-five auxiliary ladder obstruction")
    print("PASS exact Q obstruction for the minimal tagged GVC(3) lift")
    print("PASS GF(101) complete quadratic-profile obstruction through m=6")
    print("PASS 200 exact general-cubic GF(101) fibers; no m=6 survivor")
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
