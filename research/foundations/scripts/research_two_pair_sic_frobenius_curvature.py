#!/usr/bin/env python3
"""Research the Frobenius/curvature bridge for the SIC2C4 moments.

The script keeps three arithmetic objects separate:

* the integral mixed moment of the denominator-cleared radial family;
* its first-order recurrence and recurrence-operator p-curvature; and
* the Picard--Fuchs operator of the normalized angular period.

The exact derivations are recorded in
``extended-geometry/TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md``.
The differential p-curvature calculation is a bounded experiment; the
recurrences, their good-prime recurrence curvatures, and the valuation
identities are checked exactly.
"""

from __future__ import annotations

from fractions import Fraction
from math import factorial
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import GF
from sympy.polys.fields import field


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_frobenius_curvature.json"
)

M = sp.symbols("m")
X = sp.symbols("x")
DIFFERENTIAL_PRIME_BOUND = 101
DIFFERENCE_DIRECT_PRIME_BOUND = 43
VALUATION_PRIME_BOUND = 101
VALUATION_ORDER_MULTIPLIER = 2
CASES = (
    (4, 1),
    (5, 1),
    (6, 1),
    (7, 1),
    (8, 1),
    (8, 2),
    (9, 2),
    (12, 3),
    (13, 3),
)


def primes_through(bound: int) -> list[int]:
    return list(sp.primerange(2, bound + 1))


def odd_double_factorial(number: int) -> int:
    answer = 1
    for value in range(1, number + 1, 2):
        answer *= value
    return answer


def cleared_moment(degree: int, seed_power: int, order: int) -> int:
    """Return E_2(Z*(R^k*(2F)^r)^m), where d=4r+k."""
    assert degree >= 4 * seed_power >= 4
    return (
        2 ** (seed_power * order)
        * factorial(degree * order + 2)
        * factorial(seed_power * order)
        // odd_double_factorial(2 * seed_power * order + 1)
    )


def full_ratio(degree: int, seed_power: int) -> sp.Expr:
    """Return M_(m+1)/M_m for the integral denominator-cleared moments."""
    radial = sp.prod(
        degree * M + shift for shift in range(3, degree + 3)
    )
    numerator_factorials = sp.prod(
        seed_power * M + shift for shift in range(1, seed_power + 1)
    ) ** 2
    denominator_factorial = sp.prod(
        2 * seed_power * M + shift
        for shift in range(2, 2 * seed_power + 2)
    )
    return sp.cancel(
        4**seed_power
        * radial
        * numerator_factorials
        / denominator_factorial
    )


def angular_ratio(seed_power: int) -> sp.Expr:
    """Return b_r(m+1)/b_r(m) for b_r(m)=(rm)!/(2rm+1)!!."""
    numerator = sp.prod(
        seed_power * M + shift for shift in range(1, seed_power + 1)
    )
    denominator = sp.prod(
        2 * seed_power * M + 2 * shift + 1
        for shift in range(1, seed_power + 1)
    )
    return sp.cancel(numerator / denominator)


def recurrence_pair(ratio: sp.Expr) -> tuple[sp.Poly, sp.Poly]:
    """Return coprime primitive A,B with a_(m+1)/a_m=B/A."""
    numerator, denominator = sp.fraction(sp.cancel(ratio))
    numerator_poly = sp.Poly(numerator, M, domain=sp.QQ)
    denominator_poly = sp.Poly(denominator, M, domain=sp.QQ)
    common_denominator = sp.ilcm(
        *[
            coefficient.q
            for coefficient in (
                list(numerator_poly.all_coeffs())
                + list(denominator_poly.all_coeffs())
            )
        ]
    )
    numerator_poly = sp.Poly(
        sp.expand(common_denominator * numerator), M, domain=sp.ZZ
    )
    denominator_poly = sp.Poly(
        sp.expand(common_denominator * denominator), M, domain=sp.ZZ
    )
    common_content = sp.igcd(
        int(numerator_poly.content()), int(denominator_poly.content())
    )
    numerator_poly = sp.Poly(
        numerator_poly.as_expr() / common_content, M, domain=sp.ZZ
    )
    denominator_poly = sp.Poly(
        denominator_poly.as_expr() / common_content, M, domain=sp.ZZ
    )
    if denominator_poly.LC() < 0:
        numerator_poly = -numerator_poly
        denominator_poly = -denominator_poly
    assert sp.gcd(numerator_poly, denominator_poly).degree() == 0
    return denominator_poly, numerator_poly


def factor_string(poly: sp.Poly) -> str:
    return sp.sstr(sp.factor(poly.as_expr()))


def factorial_valuation(number: int, prime: int) -> int:
    answer = 0
    while number:
        number //= prime
        answer += number
    return answer


def integer_valuation(number: int, prime: int) -> int:
    answer = 0
    while number and number % prime == 0:
        answer += 1
        number //= prime
    return answer


def moment_valuation(
    degree: int, seed_power: int, order: int, prime: int
) -> int:
    power_of_two = 2 * seed_power * order if prime == 2 else 0
    return (
        power_of_two
        + factorial_valuation(degree * order + 2, prime)
        + 2 * factorial_valuation(seed_power * order, prime)
        - factorial_valuation(2 * seed_power * order + 1, prime)
    )


def ratio_valuation(
    degree: int, seed_power: int, order: int, prime: int
) -> int:
    """Return v_p(M_(m+1)/M_m) from the coprime recurrence."""
    denominator, numerator = recurrence_pair(
        full_ratio(degree, seed_power)
    )
    return recurrence_pair_valuation(
        denominator, numerator, order, prime
    )


def recurrence_pair_valuation(
    denominator: sp.Poly,
    numerator: sp.Poly,
    order: int,
    prime: int,
) -> int:
    numerator_value = abs(int(numerator.eval(order)))
    denominator_value = abs(int(denominator.eval(order)))
    return integer_valuation(numerator_value, prime) - integer_valuation(
        denominator_value, prime
    )


def matrix_multiply(left: list[list], right: list[list]) -> list[list]:
    zero = left[0][0] * 0
    return [
        [
            sum(
                (
                    left[row][index] * right[index][column]
                    for index in range(len(right))
                ),
                zero,
            )
            for column in range(len(right[0]))
        ]
        for row in range(len(left))
    ]


def matrix_add(left: list[list], right: list[list]) -> list[list]:
    return [
        [
            left[row][column] + right[row][column]
            for column in range(len(left[0]))
        ]
        for row in range(len(left))
    ]


def matrix_derivative(matrix: list[list], variable) -> list[list]:
    return [
        [entry.diff(variable) for entry in row]
        for row in matrix
    ]


def differential_p_curvature(prime: int) -> tuple[object, list[list]]:
    """Compute p-curvature for theta(2theta+1)-x(theta+1)^2.

    This is the minimal Picard--Fuchs operator for
    sum_m m!/(2m+1)!! x^m.  The denominator-cleared angular series differs
    only by the dilation x -> 2x.
    """
    rational_field, variable = field("x", GF(prime))
    zero = rational_field.zero
    one = rational_field.one
    denominator = variable * (2 - variable)

    # For L=a2 D^2+a1 D+a0, horizontal jets satisfy Y'=B Y.
    # We use the connection D-B, so A=-B below.
    connection = [
        [zero, -one],
        [
            -one / denominator,
            3 * (1 - variable) / denominator,
        ],
    ]
    curvature = connection
    for _ in range(1, prime):
        curvature = matrix_add(
            matrix_derivative(curvature, variable),
            matrix_multiply(connection, curvature),
        )
    return variable, curvature


def zero_matrix(matrix: list[list]) -> bool:
    return all(entry == 0 for row in matrix for entry in row)


def differential_curvature_audit() -> dict:
    audited_primes: list[int] = []
    pole_exponents: dict[str, dict[str, int]] = {}
    for prime in primes_through(DIFFERENTIAL_PRIME_BOUND):
        if prime == 2:
            continue
        variable, curvature = differential_p_curvature(prime)
        trace = curvature[0][0] + curvature[1][1]
        determinant = (
            curvature[0][0] * curvature[1][1]
            - curvature[0][1] * curvature[1][0]
        )
        square = matrix_multiply(curvature, curvature)
        assert trace == 0
        assert determinant == 0
        assert not zero_matrix(curvature)
        assert zero_matrix(square)

        supports: set[str] = set()
        maximum_exponents: dict[str, int] = {}
        for row in curvature:
            for entry in row:
                if entry == 0:
                    continue
                _, factors = entry.denom.factor_list()
                for factor, exponent in factors:
                    if factor == variable.numer:
                        label = "x"
                    elif factor == (variable - 2).numer:
                        label = "x-2"
                    else:
                        raise AssertionError(
                            f"unexpected p-curvature pole {factor} at p={prime}"
                        )
                    supports.add(label)
                    maximum_exponents[label] = max(
                        exponent, maximum_exponents.get(label, 0)
                    )
        assert supports == {"x", "x-2"}
        audited_primes.append(prime)
        pole_exponents[str(prime)] = maximum_exponents

    return {
        "operator": "theta*(2*theta+1)-x*(theta+1)^2",
        "period": "sum_(m>=0) m!/(2m+1)!! * x^m",
        "cleared_period_relation": "x -> 2*x",
        "primes": audited_primes,
        "result": "nonzero rank-one nilpotent p-curvature",
        "nilpotency_index": 2,
        "pole_support": ["x", "x-2"],
        "maximum_pole_exponents": pole_exponents,
        "status": "bounded exact computation, not an all-prime proof",
    }


def direct_difference_norm(
    ratio: sp.Expr, prime: int
) -> tuple[object, object]:
    rational_field, variable = field("x", GF(prime))
    product = rational_field.one
    for shift in range(prime):
        shifted = ratio.subs(M, X + shift)
        product *= rational_field.from_expr(shifted)
    return variable, product


def fraction_field_equal(left, right) -> bool:
    """Compare FracElements without relying on unit normalization."""
    return left.numer * right.denom == right.numer * left.denom


def recurrence_case(degree: int, seed_power: int) -> dict:
    ratio = full_ratio(degree, seed_power)
    denominator, numerator = recurrence_pair(ratio)
    assert numerator.degree() - denominator.degree() == degree

    for order in range(0, 12):
        left = denominator.eval(order) * cleared_moment(
            degree, seed_power, order + 1
        )
        right = numerator.eval(order) * cleared_moment(
            degree, seed_power, order
        )
        assert left == right

    good_primes = [
        prime
        for prime in primes_through(DIFFERENCE_DIRECT_PRIME_BOUND)
        if prime > degree + 2
    ]
    direct_primes = good_primes
    if (degree, seed_power) not in {(4, 1), (5, 1)}:
        direct_primes = good_primes[:2]
    for prime in direct_primes:
        variable, norm = direct_difference_norm(ratio, prime)
        expected = (
            (degree**degree % prime)
            * (variable**prime - variable) ** degree
        )
        assert fraction_field_equal(norm, expected)

    denominator_roots = [
        str(root)
        for root in sp.roots(denominator.as_expr(), M)
    ]
    return {
        "degree": degree,
        "seed_power": seed_power,
        "radial_padding": degree - 4 * seed_power,
        "minimal_recurrence": (
            f"({factor_string(denominator)})*M_(m+1)"
            f" - ({factor_string(numerator)})*M_m = 0"
        ),
        "A_forward": factor_string(denominator),
        "B_backward": factor_string(numerator),
        "degree_A": denominator.degree(),
        "degree_B": numerator.degree(),
        "degree_difference": numerator.degree() - denominator.degree(),
        "forward_singular_roots": denominator_roots,
        "good_prime_p_curvature_characteristic_polynomial": (
            f"T - ({degree}^{degree})*(m^p-m)^{degree}"
        ),
        "direct_product_primes": direct_primes,
    }


def angular_recurrence_case(seed_power: int) -> dict:
    ratio = angular_ratio(seed_power)
    denominator, numerator = recurrence_pair(ratio)
    for order in range(0, 12):
        current = Fraction(
            factorial(seed_power * order),
            odd_double_factorial(2 * seed_power * order + 1),
        )
        following = Fraction(
            factorial(seed_power * (order + 1)),
            odd_double_factorial(2 * seed_power * (order + 1) + 1),
        )
        assert denominator.eval(order) * following == numerator.eval(
            order
        ) * current
    return {
        "seed_power": seed_power,
        "A_forward": factor_string(denominator),
        "B_backward": factor_string(numerator),
        "minimal_recurrence": (
            f"({factor_string(denominator)})*b_(m+1)"
            f" - ({factor_string(numerator)})*b_m = 0"
        ),
        "good_prime_difference_p_curvature": f"2^(-{seed_power})",
    }


def valuation_audit() -> dict:
    summaries: dict[str, dict] = {}
    first_reentries: list[dict] = []
    for degree, seed_power in CASES:
        key = f"d={degree},r={seed_power}"
        denominator, numerator = recurrence_pair(
            full_ratio(degree, seed_power)
        )
        negative_jumps = 0
        reentries = 0
        first_reentry = None
        checked_primes = [
            prime
            for prime in primes_through(VALUATION_PRIME_BOUND)
            if prime > 2
        ]
        for prime in checked_primes:
            cutoff = VALUATION_ORDER_MULTIPLIER * prime
            valuations = [
                moment_valuation(degree, seed_power, order, prime)
                for order in range(1, cutoff + 1)
            ]
            for order in range(1, cutoff):
                jump = valuations[order] - valuations[order - 1]
                assert jump == recurrence_pair_valuation(
                    denominator, numerator, order, prime
                )
                if jump < 0:
                    negative_jumps += 1
            for exponent in (1, 2, 3):
                survival = [value < exponent for value in valuations]
                for index in range(1, len(survival)):
                    if not survival[index - 1] and survival[index]:
                        reentries += 1
                        candidate = {
                            "degree": degree,
                            "seed_power": seed_power,
                            "prime_power": f"{prime}^{exponent}",
                            "vanishing_order": index,
                            "reentry_order": index + 1,
                            "valuations": [
                                valuations[index - 1],
                                valuations[index],
                            ],
                        }
                        if first_reentry is None:
                            first_reentry = candidate
                        first_reentries.append(candidate)
        summaries[key] = {
            "primes_through": VALUATION_PRIME_BOUND,
            "orders_through": f"{VALUATION_ORDER_MULTIPLIER}*p",
            "prime_power_exponents": [1, 2, 3],
            "negative_consecutive_jumps": negative_jumps,
            "prime_power_reentries": reentries,
            "first_reentry": first_reentry,
        }

    # The exact re-entry already isolated in the characteristic-p note.
    assert moment_valuation(5, 1, 4, 11) == 2
    assert moment_valuation(5, 1, 5, 11) == 1
    assert ratio_valuation(5, 1, 4, 11) == -1
    assert summaries["d=4,r=1"]["negative_consecutive_jumps"] == 0
    assert summaries["d=4,r=1"]["prime_power_reentries"] == 0

    return {
        "summaries": summaries,
        "first_reentries": first_reentries[:20],
        "exact_control": {
            "family": "R*(2F), d=5, r=1",
            "prime_power": "11^2",
            "v_11(M_4)": 2,
            "v_11(M_5)": 1,
            "ratio_jump_at_m=4": -1,
            "responsible_forward_factor": "2*m+3",
            "factor_value_at_m=4": 11,
        },
    }


def main() -> None:
    recurrence_data = [
        recurrence_case(degree, seed_power)
        for degree, seed_power in CASES
    ]
    angular_data = [
        angular_recurrence_case(seed_power)
        for seed_power in range(1, 4)
    ]
    differential_data = differential_curvature_audit()
    valuations = valuation_audit()

    artifact = {
        "format": "two-pair-sic-frobenius-curvature-v1",
        "integral_moment": (
            "M_(d,r)(m)=4^(r*m)*(d*m+2)!*((r*m)!)^2/(2*r*m+1)!"
        ),
        "normalization": {
            "angular_period": "b_r(m)=(r*m)!/(2*r*m+1)!!",
            "cleared_moment_relation": (
                "M_(d,r)(m)=(d*m+2)!*2^(r*m)*b_r(m)"
            ),
            "warning": (
                "the factorial inverse-Borel factor depends on d, while "
                "the angular Picard-Fuchs period depends only on r"
            ),
        },
        "minimal_recurrences": recurrence_data,
        "angular_recurrences": angular_data,
        "differential_picard_fuchs_p_curvature": differential_data,
        "recurrence_p_curvature": {
            "definition": (
                "tau^p multiplier=product_(i=0)^(p-1) "
                "B(m+i)/A(m+i)"
            ),
            "constant_field_coordinate": "Z=m^p-m",
            "good_prime_full_moment": "d^d*Z^d",
            "good_prime_angular_period": "2^(-r)",
            "nilpotency": (
                "both rank-one recurrence curvatures are generically "
                "nonzero, hence not nilpotent"
            ),
            "singularity_loss": (
                "the norm sends every good linear factor to a scalar "
                "multiple of Z, so numerator and denominator factors "
                "cancel down to their degree difference"
            ),
        },
        "valuation_correlation": valuations,
        "conclusion": {
            "classical_differential_p_curvature": (
                "identical for all first radial lifts R^k*(2F), so it "
                "cannot recover their d-dependent valuation phases"
            ),
            "recurrence_p_curvature": (
                "its generic norm d^d*Z^d also forgets the separate "
                "zero/pole ledger needed for prime-power re-entry"
            ),
            "effective_bridge": (
                "the coprime first-order recurrence plus p-adic "
                "valuations of its individual linear factors gives the "
                "exact local transition rule"
            ),
        },
        "downstream_use": {
            "holonomic_pipeline": (
                "after creative telescoping, store primitive integral "
                "period/raw-moment recurrences, then the local factor or "
                "step-matrix Smith ledger, and only then p-curvature"
            ),
            "bidegree_33_rank_two": (
                "apply after the order-18/27 characteristic-zero "
                "recurrence and parameter denominator are certified"
            ),
            "bidegree_44_rank_two": (
                "apply after an explicit exact-rank-two closed point or "
                "component makes scalar recurrence derivation possible"
            ),
        },
        "written_analysis": (
            "extended-geometry/"
            "TWO_PAIR_SIC_FROBENIUS_CURVATURE_BRIDGE.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")

    print("PASS SIC2C4 curvature: coprime minimal recurrences")
    print("PASS SIC2C4 curvature: good-prime recurrence norm d^d*(m^p-m)^d")
    print(
        "PASS SIC2C4 curvature: normalized differential p-curvature is "
        f"nonzero nilpotent through p={DIFFERENTIAL_PRIME_BOUND}"
    )
    print("PASS SIC2C4 curvature: differential poles stay at x=0,2")
    print("PASS SIC2C4 curvature: recurrence singularity explains 11^2 re-entry")
    print(f"PASS SIC2C4 curvature: wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
