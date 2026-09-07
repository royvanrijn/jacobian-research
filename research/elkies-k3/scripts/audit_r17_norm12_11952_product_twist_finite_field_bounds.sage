#!/usr/bin/env sage -python
"""Audit the finite-field rank-bound route for the 17 product twists.

This script deliberately stops short of a Mordell--Weil rank claim.  It
certifies the good-reduction geometry of every selected product twist and
computes the first two Frobenius power sums of its elliptic L-polynomial at
the two small good primes already native to the direct alternate-Q80 model.

For a squarefree quartic ``d`` coprime to the 24 nodal fibres, the twist

    y^2 = x^3 + d^2*A*x + d^3*B

has ``chi=4``, four ``I0*`` fibres, and the original 24 ``I1`` fibres.  Its
conductor has degree 32, hence its nonconstant elliptic L-polynomial has
degree 28.  The first two power sums do not decide whether ``1-p*T`` divides
that polynomial.  The output therefore remains UNKNOWN unless a future
implementation supplies the full global Frobenius calculation (or an
equivalent Selmer bound).
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shlex
import sys

import numpy as np
from sage.all import GF, PolynomialRing, QQ, ZZ, prime_range


ROOT = Path(__file__).resolve().parents[2]
DIRECT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
SHORTLIST = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json"
)
RANK_SCREEN = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json"
)
OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-product-twist-finite-field-bound-audit-v1.json"
)
SCHEMA = "elkies-k3.r17-norm12-11952-product-twist-finite-field-bound-audit.v1"
GOOD_PRIMES = (131, 137)


def digest(path: Path) -> str:
    result = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def legendre_table(prime: int) -> np.ndarray:
    table = np.empty(prime, dtype=np.int8)
    table[0] = 0
    for value in range(1, prime):
        symbol = pow(value, (prime - 1) // 2, prime)
        table[value] = 1 if symbol == 1 else -1
    return table


def reduce_rational_coefficients(values, prime: int) -> list[int]:
    result = []
    for entry in values:
        value = QQ(entry)
        denominator = int(value.denominator() % prime)
        if denominator == 0:
            raise ZeroDivisionError(f"coefficient denominator vanishes modulo {prime}")
        result.append(
            int(value.numerator() % prime) * pow(denominator, -1, prime) % prime
        )
    return result


def eval_fp(coefficients, value: int, prime: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * value + coefficient) % prime
    return answer


def mul_fp2(left, right, prime: int, nonsquare: int):
    a, b = left
    c, d = right
    return (
        (a * c + nonsquare * b * d) % prime,
        (a * d + b * c) % prime,
    )


def eval_fp2(coefficients, value, prime: int, nonsquare: int):
    answer = (0, 0)
    for coefficient in reversed(coefficients):
        answer = mul_fp2(answer, value, prime, nonsquare)
        answer = ((answer[0] + coefficient) % prime, answer[1])
    return answer


def first_nonsquare(prime: int) -> int:
    for value in range(2, prime):
        if pow(value, (prime - 1) // 2, prime) == prime - 1:
            return value
    raise ArithmeticError(prime)


def factor_degrees(polynomial) -> list[int]:
    return sorted(
        int(factor.degree())
        for factor, multiplicity in polynomial.factor()
        for _unused in range(int(multiplicity))
    )


def reduce_model(a_values, b_values, prime: int):
    field = GF(prime)
    ring = PolynomialRing(field, "u")
    a = ring(reduce_rational_coefficients(a_values, prime))
    b = ring(reduce_rational_coefficients(b_values, prime))
    discriminant_core = 4 * a**3 + 27 * b**2
    if a.degree() != 8 or b.degree() != 12 or discriminant_core.degree() != 24:
        raise ArithmeticError("degree_drop")
    if not discriminant_core.is_squarefree():
        raise ArithmeticError("discriminant_not_squarefree")
    if discriminant_core.gcd(a) != 1:
        raise ArithmeticError("non_nodal_bad_fibre")
    return ring, a, b, discriminant_core


def reduce_twist(ring, values, discriminant_core, prime: int):
    d = ring(reduce_rational_coefficients(values, prime))
    if d.degree() != 4:
        raise ArithmeticError("twist_degree_drop")
    if not d.is_squarefree():
        raise ArithmeticError("twist_not_squarefree")
    if d.gcd(discriminant_core) != 1:
        raise ArithmeticError("twist_meets_original_bad_fibre")
    return d


def reduction_failure(a_values, b_values, targets, prime: int):
    try:
        ring, _a, _b, discriminant_core = reduce_model(a_values, b_values, prime)
    except ZeroDivisionError:
        return "coefficient_denominator_bad"
    except ArithmeticError as error:
        return str(error)
    for target in targets:
        try:
            reduce_twist(
                ring,
                target["product_quartic_coefficients_low_to_high"],
                discriminant_core,
                prime,
            )
        except ZeroDivisionError:
            return "twist_denominator_bad"
        except ArithmeticError as error:
            return str(error)
    return None


def fp_power_sums(a_coefficients, b_coefficients, twists, prime: int):
    character = legendre_table(prime)
    xs = np.arange(prime, dtype=np.int64)
    x3 = xs * xs % prime * xs % prime
    sums = np.zeros(len(twists), dtype=np.int64)
    local_character_sums = []
    for parameter in range(prime):
        av = eval_fp(a_coefficients, parameter, prime)
        bv = eval_fp(b_coefficients, parameter, prime)
        fibre_sum = int(character[(x3 + av * xs + bv) % prime].sum())
        local_character_sums.append(fibre_sum)
        for index, twist in enumerate(twists):
            sums[index] += fibre_sum * int(
                character[eval_fp(twist, parameter, prime)]
            )
    infinity_sum = int(
        character[
            (x3 + a_coefficients[8] * xs + b_coefficients[12]) % prime
        ].sum()
    )
    local_character_sums.append(infinity_sum)
    for index, twist in enumerate(twists):
        sums[index] += infinity_sum * int(character[twist[4]])
    return [int(value) for value in sums], local_character_sums


def fp2_power_sums(a_coefficients, b_coefficients, twists, prime: int):
    nonsquare = first_nonsquare(prime)
    character = legendre_table(prime)
    coordinates = np.arange(prime, dtype=np.int64)
    xa = np.repeat(coordinates, prime)
    xb = np.tile(coordinates, prime)
    x2a = (xa * xa + nonsquare * xb * xb) % prime
    x2b = 2 * xa * xb % prime
    x3a = (x2a * xa + nonsquare * x2b * xb) % prime
    x3b = (x2a * xb + x2b * xa) % prime
    sums = np.zeros(len(twists), dtype=np.int64)

    def fibre_character_sum(av, bv) -> int:
        aa, ab = av
        ba, bb = bv
        rhs_a = (x3a + aa * xa + nonsquare * ab * xb + ba) % prime
        rhs_b = (x3b + aa * xb + ab * xa + bb) % prime
        norm = (rhs_a * rhs_a - nonsquare * rhs_b * rhs_b) % prime
        return int(character[norm].sum())

    def fp2_character(value) -> int:
        a, b = value
        norm = (a * a - nonsquare * b * b) % prime
        return int(character[norm])

    local_character_sum_sha = sha256()
    for ta in range(prime):
        for tb in range(prime):
            parameter = (ta, tb)
            fibre_sum = fibre_character_sum(
                eval_fp2(a_coefficients, parameter, prime, nonsquare),
                eval_fp2(b_coefficients, parameter, prime, nonsquare),
            )
            local_character_sum_sha.update(f"{fibre_sum},".encode("ascii"))
            for index, twist in enumerate(twists):
                sums[index] += fibre_sum * fp2_character(
                    eval_fp2(twist, parameter, prime, nonsquare)
                )
    infinity_sum = fibre_character_sum(
        (a_coefficients[8], 0), (b_coefficients[12], 0)
    )
    local_character_sum_sha.update(f"{infinity_sum}".encode("ascii"))
    for index, twist in enumerate(twists):
        sums[index] += infinity_sum * fp2_character((twist[4], 0))
    return [int(value) for value in sums], nonsquare, local_character_sum_sha.hexdigest()


def load_inputs():
    direct = json.loads(DIRECT.read_text())
    shortlist = json.loads(SHORTLIST.read_text())
    rank_screen = json.loads(RANK_SCREEN.read_text())
    if direct.get("status") != "PASS_EXACT_DIRECT_TWO_NEIGHBOR_EQUATION_FRAME_AND_SECTIONS":
        raise ValueError("unexpected direct-fibration status")
    if shortlist.get("schema") != "elkies-k3.r17-norm12-11952-v4-pair-shortlist.v1":
        raise ValueError("unexpected V4 shortlist schema")
    if rank_screen.get("schema") != "elkies-k3.r17-norm12-11952-v4-base-rank-screen.v1":
        raise ValueError("unexpected base-rank screen schema")
    exact_ranks = {
        int(row["shortlist_rank"])
        for row in rank_screen["results"]
        if (row.get("rank_lower_bound"), row.get("rank_upper_bound")) == (1, 1)
    }
    targets = [
        row for row in shortlist["pairs"]
        if int(row["shortlist_rank"]) in exact_ranks
    ]
    targets.sort(key=lambda row: int(row["shortlist_rank"]))
    if len(targets) != len(exact_ranks) or len(targets) != 17:
        raise ArithmeticError("the exact-rank-one target set is no longer the pinned 17")
    model = direct["weierstrass_model"]
    if model.get("degrees_A_B_Delta") != [8, 12, 24]:
        raise ArithmeticError("direct model degree pattern changed")
    return model, targets


def build_payload():
    model, targets = load_inputs()
    a_values = model["A_coefficients_low_to_high"]
    b_values = model["B_coefficients_low_to_high"]

    prime_scan = []
    for prime in prime_range(5, max(GOOD_PRIMES) + 1):
        failure = reduction_failure(a_values, b_values, targets, int(prime))
        prime_scan.append(
            {
                "prime": int(prime),
                "all_17_good": failure is None,
                "first_failure": failure,
            }
        )
    common_good = [row["prime"] for row in prime_scan if row["all_17_good"]]
    if common_good != list(GOOD_PRIMES):
        raise ArithmeticError(f"unexpected common good-prime scan: {common_good}")

    target_records = {
        row["pair_key"]: {
            "shortlist_rank": int(row["shortlist_rank"]),
            "pair_key": row["pair_key"],
            "labels": row["labels"],
            "product_quartic_coefficients_low_to_high": row[
                "product_quartic_coefficients_low_to_high"
            ],
            "reductions": [],
        }
        for row in targets
    }
    reduction_summaries = []
    for prime in GOOD_PRIMES:
        ring, a, b, discriminant_core = reduce_model(a_values, b_values, prime)
        twists = [
            reduce_twist(
                ring,
                row["product_quartic_coefficients_low_to_high"],
                discriminant_core,
                prime,
            )
            for row in targets
        ]
        a_coefficients = [int(a[index]) for index in range(9)]
        b_coefficients = [int(b[index]) for index in range(13)]
        twist_coefficients = [
            [int(d[index]) for index in range(5)] for d in twists
        ]
        first_sums, fp_local_sums = fp_power_sums(
            a_coefficients, b_coefficients, twist_coefficients, prime
        )
        second_sums, nonsquare, fp2_local_sha = fp2_power_sums(
            a_coefficients, b_coefficients, twist_coefficients, prime
        )
        reduction_summaries.append(
            {
                "prime": prime,
                "fp2_nonsquare": nonsquare,
                "original_discriminant_factor_degrees": factor_degrees(
                    discriminant_core
                ),
                "original_nodal_fibre_count_geometric": 24,
                "fp_local_character_sums_sha256": sha256(
                    ",".join(map(str, fp_local_sums)).encode("ascii")
                ).hexdigest(),
                "fp2_local_character_sums_sha256": fp2_local_sha,
            }
        )
        for index, (target, twist) in enumerate(zip(targets, twists)):
            s1 = ZZ(first_sums[index])
            s2 = ZZ(second_sums[index])
            c2_numerator = s1**2 - s2
            if c2_numerator % 2:
                raise ArithmeticError("Newton coefficient is not integral")
            record = {
                "prime": prime,
                "good_reduction": True,
                "twist_factor_degrees": factor_degrees(twist),
                "twist_discriminant_gcd_degree": int(
                    twist.gcd(discriminant_core).degree()
                ),
                "fibre_configuration_geometric": "4I0*+24I1",
                "elliptic_L_frobenius_power_sums_n1_n2": [int(s1), int(s2)],
                "elliptic_L_coefficients_through_T2": [
                    1,
                    int(-s1),
                    int(c2_numerator // 2),
                ],
                "moments_computed": 2,
                "moments_needed_with_functional_equation": 14,
                "rank_zero_decided": False,
                "status": "UNKNOWN_PARTIAL_FROBENIUS_DATA",
            }
            target_records[target["pair_key"]]["reductions"].append(record)

    inputs = (Path(__file__).resolve(), DIRECT, SHORTLIST, RANK_SCREEN)
    return {
        "schema": SCHEMA,
        "status": "UNKNOWN_NO_FINITE_FIELD_MORDELL_WEIL_UPPER_BOUND",
        "source_label": "norm12-orbit-11952",
        "target_selection": {
            "criterion": "exact V4 base-Jacobian rank interval [1,1]",
            "target_count": 17,
            "shortlist_ranks": [int(row["shortlist_rank"]) for row in targets],
        },
        "surface_geometry": {
            "twist_equation": "Y^2=X^3+d^2*A*X+d^3*B",
            "degree_d": 4,
            "arithmetic_genus_chi": 4,
            "h11": 40,
            "geometric_fibre_configuration": "4I0*+24I1",
            "trivial_lattice_rank": 18,
            "tame_conductor_degree": 32,
            "nonconstant_elliptic_L_polynomial_degree": 28,
            "geometric_mordell_weil_hodge_upper_bound": 22,
        },
        "displayed_model_prime_scan": {
            "range": "odd primes 5 through 137",
            "common_good_primes": common_good,
            "first_common_good_prime": common_good[0],
            "rows": prime_scan,
            "qualification": (
                "This scans the displayed rational Weierstrass model. A separate "
                "integral re-minimization could be needed before interpreting a "
                "coefficient-denominator failure as intrinsic bad reduction."
            ),
        },
        "reduction_summaries": reduction_summaries,
        "targets": list(target_records.values()),
        "rank_bound_gate": {
            "characteristic_zero_specialization": (
                "A rank-zero good reduction over F_p(u) would exclude every "
                "non-torsion characteristic-zero product-twist section."
            ),
            "sufficient_frobenius_test": (
                "Compute the complete degree-28 elliptic L-polynomial and verify "
                "that 1-p*T is not a factor."
            ),
            "functional_equation_workload": (
                "With the functional-equation sign and exact reconstruction "
                "bounds, power sums through n=14 determine a degree-28 polynomial."
            ),
            "completed_power_sums_per_reduction": 2,
            "missing_power_sums_per_reduction": 12,
            "existing_k3_shortcut_inapplicable": (
                "The p,p^2 endpoint argument had only a three-dimensional "
                "transcendental remainder after 19 known K3 classes. Here the "
                "chi=4 twist has a degree-28 nonconstant L-polynomial."
            ),
            "function_field_descent_boundary": (
                "The 2-torsion cubic defines a genus-ten cover; SageMath and PARI "
                "in this repository do not provide the required general "
                "function-field 2-Selmer computation."
            ),
        },
        "method": {
            "local_trace_identity": (
                "a_(p^n)(E^d_u)=chi_(p^n)(d(u))*a_(p^n)(E_u); "
                "the I0* branch fibres have local trace zero"
            ),
            "frobenius_power_sum": (
                "Tr(Frob^n on H1(P1bar,j_*R1))=-sum_(u in P1(F_(p^n))) "
                "a_(p^n)(E^d_u)"
            ),
            "fp_counts": "exact vectorized quadratic-character sums",
            "fp2_counts": "exact vectorized norm-character sums",
            "numpy_version": np.__version__,
        },
        "proof_boundary": (
            "The good-reduction geometry, conductor/L-degree calculation, and "
            "the first two Frobenius power sums at p=131 and p=137 are exact for "
            "all seventeen selected product twists. No complete L-polynomial, "
            "finite-field Mordell--Weil upper bound, or characteristic-zero "
            "nonexistence theorem is obtained. Every target remains UNKNOWN."
        ),
        "inputs": {display_path(path): digest(path) for path in inputs},
        "reproducing_command": shlex.join(sys.argv),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build_payload()
    if args.check:
        stored = json.loads(args.output.read_text())
        # argv is intentionally allowed to differ between production and check.
        stored.pop("reproducing_command", None)
        payload.pop("reproducing_command", None)
        if stored != payload:
            raise ArithmeticError("stored finite-field audit does not replay")
        print(
            f"ALTPRODUCTFFCHECK|targets=17|primes=131,137|moments=2/14|"
            f"status={payload['status']}|output={display_path(args.output)}",
            flush=True,
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"ALTPRODUCTFF|targets=17|primes=131,137|moments=2/14|"
        f"status={payload['status']}|output={display_path(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
