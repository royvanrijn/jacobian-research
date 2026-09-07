#!/usr/bin/env python3
"""Bounded CRT/Gauss neighborhood search around Nagao's rank-21 curve.

Nagao denotes the published specialization by ``E_{14721/376}``, whereas the
six-root constructor ``q(X-T)q(X+T)`` in :mod:`nagao_1994` uses
``T=14721/188``.  This script keeps both parameters visible and searches the
constructor parameter.

The search profile copies five exact p-adic features of the published curve.
It combines them by CRT, Gauss-reduces the lattice

    {(a,b): a-r*b == 0 (mod M)},

and exhausts a declared square of lattice coefficients.  Good-prime scores
are applied in leakage-free stages and conductor computations are restricted
to the final survivors.  No point search or rank computation is performed:
matching a rank-21 curve's local profile is not evidence of rank 21.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from decimal import Decimal
from fractions import Fraction
import hashlib
import json
from math import gcd, log
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
from typing import Any, Sequence

from crt_lattice import crt_pair, gauss_reduce
from ek_k3 import legendre_symbol, primes_up_to, rational_to_string
from multiple_root_lifting import (
    affine_variable_coefficients,
    fixed_divisor_valuation,
)
from nagao_1994 import (
    PRIMARY_SOURCE,
    RANK21_CONSTRUCTION,
    RANK21_CONSTRUCTOR_PARAMETER,
    RANK21_CONDUCTOR_FACTORIZATION,
    RANK21_PUBLISHED_CONDUCTOR,
    RANK21_PUBLISHED_MODEL,
    RANK21_PUBLISHED_PARAMETER,
    RANK21_PUBLISHED_POINTS,
    RANK21_ROOTS,
    factorization_product,
    point_on_extended_weierstrass,
    rank21_short_jacobian_coefficients,
    short_jacobian_coefficients,
)
from pari_bridge import minimal_curve_data, pari_version


Q = Fraction
TARGET_LOG_CONDUCTOR = Decimal("182.72")
LOCAL_PRIMES = (5, 7, 13, 17, 23)
REPRODUCING_COMMAND = (
    "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
    "elliptic-curves/cas/search_nagao_rank21_neighborhood.py"
)


def _integer_discriminant_polynomial() -> tuple[int, ...]:
    coefficients = RANK21_CONSTRUCTION.primitive_discriminant_polynomial
    if any(coefficient.denominator != 1 for coefficient in coefficients):
        raise AssertionError("the rank-21 primitive discriminant is not integral")
    answer = tuple(int(coefficient) for coefficient in coefficients)
    if len(answer) != 21:
        raise AssertionError("the rank-21 discriminant did not have degree 20")
    return answer


DISCRIMINANT_POLYNOMIAL = _integer_discriminant_polynomial()


@dataclass(frozen=True)
class LocalCondition:
    label: str
    prime: int
    exponent: int
    residue: int
    forced_discriminant_valuation: int
    use_in_search: bool

    @property
    def modulus(self) -> int:
        return self.prime**self.exponent


# Both p=5 statements are retained: the first verifies the broad ball and the
# second is the record's narrower lift.  The p=23 fixed divisor is verified
# separately from the residue ball that buys the record one additional power.
LOCAL_CONDITIONS = (
    LocalCondition("p5-broad", 5, 1, 2, 3, False),
    LocalCondition("p5-record-lift", 5, 2, 17, 4, True),
    LocalCondition("p7-record", 7, 1, 0, 4, True),
    LocalCondition("p13-record", 13, 1, 3, 4, True),
    LocalCondition("p17-record", 17, 1, 16, 3, True),
    LocalCondition("p23-record", 23, 1, 6, 3, True),
)


@dataclass(frozen=True)
class NeighborhoodCandidate:
    identifier: str
    numerator: int
    denominator: int
    height: int

    @property
    def parameter(self) -> Fraction:
        return Q(self.numerator, self.denominator)

    @property
    def coefficients(self) -> tuple[Fraction, ...]:
        return short_jacobian_coefficients(RANK21_CONSTRUCTION, self.parameter)


@dataclass(frozen=True)
class ScoreCandidate:
    identifier: str
    coefficients: tuple[Fraction, ...]


def p_adic_valuation(value: int, prime: int) -> int:
    if value == 0:
        raise ValueError("the p-adic valuation of zero is not finite")
    answer = 0
    value = abs(value)
    while value % prime == 0:
        answer += 1
        value //= prime
    return answer


def forced_valuation(condition: LocalCondition) -> int:
    """Prove a lower bound on D(T) throughout one integral p-adic ball."""

    transformed = affine_variable_coefficients(
        DISCRIMINANT_POLYNOMIAL,
        condition.residue,
        condition.modulus,
    )
    return fixed_divisor_valuation(transformed, condition.prime)


def fixed_discriminant_valuation(prime: int) -> int:
    return fixed_divisor_valuation(DISCRIMINANT_POLYNOMIAL, prime)


def _fraction_mod(value: Fraction, prime: int) -> int:
    if value.denominator % prime == 0:
        raise ValueError("a coefficient denominator is not a p-adic unit")
    return value.numerator * pow(value.denominator, -1, prime) % prime


def classify_condition(condition: LocalCondition) -> dict[str, Any]:
    """Classify the uniform reduction on a bad p-adic ball for p>=5."""

    # T=0 is a removable specialization of the polynomial formulas but the
    # direct constructor divides by T^2.  Evaluating at T=p has the same
    # residue and gives the polynomial continuation modulo p.
    representative = condition.residue % condition.prime
    if representative == 0:
        representative = condition.prime
    coefficients = short_jacobian_coefficients(
        RANK21_CONSTRUCTION, Q(representative)
    )
    coefficient_a = _fraction_mod(coefficients[3], condition.prime)
    coefficient_b = _fraction_mod(coefficients[4], condition.prime)
    c4 = (-48 * coefficient_a) % condition.prime
    if c4 == 0:
        return {
            "reduction": "additive or unresolved",
            "conductor_exponent": None,
            "split_multiplicative": None,
            "c4_mod_prime": 0,
            "proof": (
                "the presented c4 vanishes modulo p, so this coarse test does "
                "not resolve the minimal reduction"
            ),
        }
    double_root = (
        -3
        * coefficient_b
        * pow(2 * coefficient_a, -1, condition.prime)
        % condition.prime
    )
    tangent_value = 3 * double_root % condition.prime
    tangent_symbol = legendre_symbol(tangent_value, condition.prime)
    if tangent_symbol == 0:
        raise AssertionError("a multiplicative nodal fiber became cuspidal")
    split = tangent_symbol == 1
    return {
        "reduction": "split multiplicative" if split else "nonsplit multiplicative",
        "conductor_exponent": 1,
        "split_multiplicative": split,
        "c4_mod_prime": c4,
        "double_root_mod_prime": double_root,
        "tangent_value_mod_prime": tangent_value,
        "tangent_legendre_symbol": tangent_symbol,
        "proof": (
            "D vanishes on the ball, c4 is a unit modulo p, and the nodal "
            "tangent symbol depends only on T modulo p"
        ),
    }


def discriminant_valuation_at_rational(
    numerator: int, denominator: int, prime: int
) -> int:
    """Return v_p(b^20 D(a/b)) when b is a p-adic unit."""

    if denominator % prime == 0:
        raise ValueError("the rational denominator must be a p-adic unit")
    degree = len(DISCRIMINANT_POLYNOMIAL) - 1
    value = sum(
        coefficient * numerator**index * denominator ** (degree - index)
        for index, coefficient in enumerate(DISCRIMINANT_POLYNOMIAL)
    )
    return p_adic_valuation(value, prime)


def search_conditions() -> tuple[LocalCondition, ...]:
    answer = tuple(condition for condition in LOCAL_CONDITIONS if condition.use_in_search)
    primes = [condition.prime for condition in answer]
    if len(set(primes)) != len(primes):
        raise AssertionError("the CRT search profile repeats a prime")
    return answer


def profile_crt() -> tuple[int, int]:
    residue, modulus = 0, 1
    for condition in search_conditions():
        residue, modulus = crt_pair(
            residue, modulus, condition.residue, condition.modulus
        )
    return residue, modulus


def enumerate_neighborhood(coefficient_radius: int) -> tuple[NeighborhoodCandidate, ...]:
    """Exhaust the sign-quotiented square in the Gauss-reduced CRT basis."""

    if coefficient_radius < 1:
        raise ValueError("the coefficient radius must be positive")
    residue, modulus = profile_crt()
    basis = gauss_reduce((modulus, 0), (residue, 1))
    primitive_vectors: set[tuple[int, int]] = set()
    for left in range(-coefficient_radius, coefficient_radius + 1):
        for right in range(-coefficient_radius, coefficient_radius + 1):
            if left == 0 and right == 0:
                continue
            numerator = left * basis[0][0] + right * basis[1][0]
            denominator = left * basis[0][1] + right * basis[1][1]
            if denominator == 0 or gcd(denominator, modulus) != 1:
                continue
            common = gcd(abs(numerator), abs(denominator))
            numerator //= common
            denominator //= common
            if denominator < 0:
                numerator, denominator = -numerator, -denominator
            if numerator == 0 or gcd(denominator, modulus) != 1:
                continue
            if (numerator - residue * denominator) % modulus:
                raise AssertionError("primitive normalization lost the CRT class")
            degree = len(DISCRIMINANT_POLYNOMIAL) - 1
            discriminant = sum(
                coefficient * numerator**index * denominator ** (degree - index)
                for index, coefficient in enumerate(DISCRIMINANT_POLYNOMIAL)
            )
            if discriminant:
                primitive_vectors.add((numerator, denominator))
    ordered = sorted(
        primitive_vectors,
        key=lambda pair: (max(abs(pair[0]), pair[1]), abs(pair[0]), pair[1]),
    )
    return tuple(
        NeighborhoodCandidate(
            identifier=f"nagao-rank21-neighbor-{index:04d}",
            numerator=numerator,
            denominator=denominator,
            height=max(abs(numerator), denominator),
        )
        for index, (numerator, denominator) in enumerate(ordered, start=1)
    )


def _gp_rational(value: Fraction) -> str:
    return f"({value.numerator}/{value.denominator})"


def score_candidates_with_pari(
    candidates: Sequence[ScoreCandidate],
    *,
    cutoff: int,
    timeout: float,
    stack_bytes: int,
) -> dict[str, dict[str, Any]]:
    """Compute the declared Nagao-style good-prime score in one bounded batch."""

    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    if not candidates or cutoff < 5 or timeout <= 0 or stack_bytes < 8_000_000:
        raise ValueError("invalid PARI score bounds")
    last_prime = primes_up_to(cutoff)[-1]
    commands = ["default(realprecision,80);"]
    for index, candidate in enumerate(candidates):
        vector = ",".join(_gp_rational(value) for value in candidate.coefficients)
        commands.extend(
            (
                f"E=ellminimalmodel(ellinit([{vector}]));",
                "S=0;USED=0;BAD=0;",
                (
                    f"forprime(p=5,{cutoff},"
                    "if(valuation(E.disc,p)>0,BAD++,"
                    "A=ellap(E,p);S+=(2-A)/(p+1-A)*log(p);USED++));"
                ),
                f'print("ROW|{index}|",S,"|",USED,"|",BAD,"|{last_prime}");',
            )
        )
    commands.append("quit")
    result = subprocess.run(
        [executable, "-q", "-s", str(stack_bytes)],
        input="\n".join(commands) + "\n",
        text=True,
        capture_output=True,
        timeout=timeout,
    )
    if result.returncode != 0 or "***" in result.stderr:
        raise RuntimeError(f"PARI/GP score failed: {result.stderr.strip()}")
    answer: dict[str, dict[str, Any]] = {}
    for line in result.stdout.splitlines():
        if not line.startswith("ROW|"):
            continue
        _, index_text, score, used, bad, observed_last = line.split("|")
        candidate = candidates[int(index_text)]
        answer[candidate.identifier] = {
            "score": score,
            "good_primes_used": int(used),
            "bad_primes_skipped": int(bad),
            "last_numerical_prime": int(observed_last),
        }
    if len(answer) != len(candidates):
        raise RuntimeError("PARI omitted one or more score records")
    return answer


def parse_integer_tuple(value: str) -> tuple[int, ...]:
    try:
        answer = tuple(int(part) for part in value.split(",") if part)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected comma-separated integers") from error
    if not answer or any(integer < 1 for integer in answer):
        raise argparse.ArgumentTypeError("all values must be positive")
    return answer


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficient-radius", type=int, default=12)
    parser.add_argument("--stages", type=parse_integer_tuple, default=(200, 2_000, 10_000))
    parser.add_argument("--keep-counts", type=parse_integer_tuple, default=(32, 12, 6))
    parser.add_argument("--score-timeout", type=float, default=180.0)
    parser.add_argument("--conductor-timeout", type=float, default=45.0)
    parser.add_argument("--calibration-timeout", type=float, default=60.0)
    parser.add_argument("--stack-bytes", type=int, default=256_000_000)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            root
            / "artifacts"
            / "generated-results"
            / "elliptic_nagao_rank21_neighborhood.json"
        ),
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.coefficient_radius < 1:
        raise SystemExit("--coefficient-radius must be positive")
    if len(args.stages) != len(args.keep_counts):
        raise SystemExit("--stages and --keep-counts must have equal length")
    if any(left >= right for left, right in zip(args.stages, args.stages[1:])):
        raise SystemExit("--stages must be strictly increasing")
    if any(left < right for left, right in zip(args.keep_counts, args.keep_counts[1:])):
        raise SystemExit("--keep-counts must be nonincreasing")
    if min(args.score_timeout, args.conductor_timeout, args.calibration_timeout) <= 0:
        raise SystemExit("all timeouts must be positive")
    if args.stack_bytes < 8_000_000:
        raise SystemExit("the PARI stack bound is too small")

    local_verification = []
    for condition in LOCAL_CONDITIONS:
        observed = forced_valuation(condition)
        if observed != condition.forced_discriminant_valuation:
            raise AssertionError(
                f"{condition.label}: expected forced valuation "
                f"{condition.forced_discriminant_valuation}, observed {observed}"
            )
        local_verification.append(
            {
                **asdict(condition),
                "modulus": condition.modulus,
                "observed_exact_forced_discriminant_valuation": observed,
                "classification": classify_condition(condition),
            }
        )
    fixed_at_23 = fixed_discriminant_valuation(23)
    if fixed_at_23 != 2:
        raise AssertionError("the degree-20 discriminant lost its fixed 23^2 divisor")

    residue, modulus = profile_crt()
    reduced_basis = gauss_reduce((modulus, 0), (residue, 1))
    published_numerator = RANK21_CONSTRUCTOR_PARAMETER.numerator
    published_denominator = RANK21_CONSTRUCTOR_PARAMETER.denominator
    if (published_numerator - residue * published_denominator) % modulus:
        raise AssertionError("Nagao's published constructor parameter missed the profile")

    candidates = enumerate_neighborhood(args.coefficient_radius)
    if args.keep_counts[0] > len(candidates):
        raise SystemExit("the first keep count exceeds the candidate population")
    by_identifier = {candidate.identifier: candidate for candidate in candidates}

    survivors = list(candidates)
    stages: list[dict[str, Any]] = []
    latest_scores: dict[str, dict[str, Any]] = {}
    for cutoff, keep_count in zip(args.stages, args.keep_counts):
        score_inputs = tuple(
            ScoreCandidate(candidate.identifier, candidate.coefficients)
            for candidate in survivors
        )
        scores = score_candidates_with_pari(
            score_inputs,
            cutoff=cutoff,
            timeout=args.score_timeout,
            stack_bytes=args.stack_bytes,
        )
        ranked = sorted(
            (
                {
                    "candidate_id": candidate.identifier,
                    "constructor_parameter": rational_to_string(candidate.parameter),
                    "height": candidate.height,
                    **scores[candidate.identifier],
                }
                for candidate in survivors
            ),
            key=lambda record: (
                -Decimal(record["score"]),
                record["height"],
                record["candidate_id"],
            ),
        )
        retained = ranked[: min(keep_count, len(ranked))]
        for record in retained:
            latest_scores[record["candidate_id"]] = record
        stages.append(
            {
                "cutoff": cutoff,
                "population_scored": len(ranked),
                "keep_count": len(retained),
                "ranked_population": ranked,
                "retained_candidate_ids": [
                    record["candidate_id"] for record in retained
                ],
            }
        )
        survivors = [by_identifier[record["candidate_id"]] for record in retained]

    conductor_records: list[dict[str, Any]] = []
    conductor_errors: list[dict[str, str]] = []
    for candidate in survivors:
        local_checks = []
        for condition in search_conditions():
            if candidate.denominator % condition.prime == 0:
                raise AssertionError("a CRT candidate acquired a bad local denominator")
            if (
                candidate.numerator - condition.residue * candidate.denominator
            ) % condition.modulus:
                raise AssertionError("a CRT candidate missed a profile condition")
            valuation = discriminant_valuation_at_rational(
                candidate.numerator, candidate.denominator, condition.prime
            )
            if valuation < condition.forced_discriminant_valuation:
                raise AssertionError("a candidate lost a forced discriminant power")
            local_checks.append(
                {
                    "label": condition.label,
                    "prime": condition.prime,
                    "observed_discriminant_valuation": valuation,
                    "forced_lower_bound": condition.forced_discriminant_valuation,
                    "proved_reduction_on_ball": classify_condition(condition)["reduction"],
                }
            )
        try:
            curve = minimal_curve_data(
                candidate.coefficients,
                timeout=args.conductor_timeout,
                local_primes=LOCAL_PRIMES,
                stack_bytes=args.stack_bytes,
            )
            for check in local_checks:
                pari_local = curve["local_reduction"][str(check["prime"])]
                if (
                    pari_local["conductor_exponent"] != 1
                    or pari_local["minimal_c4_valuation"] != 0
                    or pari_local["ellap"] != 1
                    or pari_local["minimal_discriminant_valuation"]
                    != check["observed_discriminant_valuation"]
                ):
                    raise AssertionError("PARI disagreed with an exact split local profile")
            conductor_records.append(
                {
                    **latest_scores[candidate.identifier],
                    "local_checks": local_checks,
                    "curve": curve,
                    "below_strict_log_conductor_target": (
                        Decimal(curve["log_conductor"]) < TARGET_LOG_CONDUCTOR
                    ),
                    "rank_claim": None,
                }
            )
        except Exception as error:
            conductor_errors.append(
                {"candidate_id": candidate.identifier, "error": str(error)}
            )

    published_score = {}
    published_score_input = ScoreCandidate(
        "nagao-published-rank21",
        rank21_short_jacobian_coefficients(RANK21_PUBLISHED_PARAMETER),
    )
    for cutoff in args.stages:
        published_score[str(cutoff)] = score_candidates_with_pari(
            (published_score_input,),
            cutoff=cutoff,
            timeout=args.score_timeout,
            stack_bytes=args.stack_bytes,
        )[published_score_input.identifier]

    final_record_score = Decimal(published_score[str(args.stages[-1])]["score"])
    for record in conductor_records:
        record["score_exceeds_published_rank21_at_final_cutoff"] = (
            Decimal(record["score"]) > final_record_score
        )
    conductor_records.sort(
        key=lambda record: (
            not record["below_strict_log_conductor_target"],
            -Decimal(record["score"]),
            Decimal(record["curve"]["log_conductor"]),
        )
    )

    exact_point_checks = [
        point_on_extended_weierstrass(RANK21_PUBLISHED_MODEL, point)
        for point in RANK21_PUBLISHED_POINTS
    ]
    calibration: dict[str, Any] = {
        "published_parameter_notation": rational_to_string(
            RANK21_PUBLISHED_PARAMETER
        ),
        "constructor_parameter": rational_to_string(RANK21_CONSTRUCTOR_PARAMETER),
        "factor_two_identity_checked": (
            RANK21_CONSTRUCTOR_PARAMETER == 2 * RANK21_PUBLISHED_PARAMETER
        ),
        "printed_model": list(RANK21_PUBLISHED_MODEL),
        "printed_conductor": RANK21_PUBLISHED_CONDUCTOR,
        "printed_conductor_factorization": [
            list(factor) for factor in RANK21_CONDUCTOR_FACTORIZATION
        ],
        "factorization_product_checked": (
            factorization_product(RANK21_CONDUCTOR_FACTORIZATION)
            == RANK21_PUBLISHED_CONDUCTOR
        ),
        "printed_points_checked_exactly": len(exact_point_checks),
        "all_printed_points_on_printed_model": all(exact_point_checks),
        "published_independence_status": (
            "Nagao's independence result is cited; this neighborhood search "
            "does not provide a new independence proof"
        ),
        "scores": published_score,
    }
    try:
        pari_calibration = minimal_curve_data(
            rank21_short_jacobian_coefficients(RANK21_PUBLISHED_PARAMETER),
            timeout=args.calibration_timeout,
            local_primes=LOCAL_PRIMES,
            stack_bytes=args.stack_bytes,
        )
        if tuple(pari_calibration["minimal_model"]) != RANK21_PUBLISHED_MODEL:
            raise AssertionError("PARI did not reproduce Nagao's printed model")
        if pari_calibration["conductor"] != RANK21_PUBLISHED_CONDUCTOR:
            raise AssertionError("PARI did not reproduce Nagao's conductor")
        calibration["pari_replay"] = pari_calibration
        calibration["pari_replay_matches_printed_model_and_conductor"] = True
    except Exception as error:
        calibration["pari_replay_error"] = str(error)
        calibration["pari_replay_matches_printed_model_and_conductor"] = False

    script_path = Path(__file__).resolve()
    artifact = {
        "schema_version": 1,
        "status": (
            "bounded exact CRT/Gauss and staged score/conductor experiment; "
            "published rank 21 is calibration only, no searched specialization "
            "has a rank claim, and no target hit is claimed"
        ),
        "primary_source": PRIMARY_SOURCE,
        "strict_target": {
            "rank_at_least": 21,
            "log_conductor_less_than": str(TARGET_LOG_CONDUCTOR),
            "alternative_rank_at_least": 30,
        },
        "target_hits": [],
        "family": {
            "root_tuple": list(RANK21_ROOTS),
            "published_parameter_convention": "Nagao E_t uses t=14721/376",
            "constructor_parameter_convention": (
                "q(X-T)q(X+T) uses T=2t=14721/188"
            ),
            "discriminant_polynomial_degree": len(DISCRIMINANT_POLYNOMIAL) - 1,
            "discriminant_polynomial_fixed_23_valuation": fixed_at_23,
        },
        "local_verification": local_verification,
        "search_profile": {
            "conditions": [asdict(condition) for condition in search_conditions()],
            "crt_residue": residue,
            "crt_modulus": modulus,
            "gauss_reduced_basis": [list(vector) for vector in reduced_basis],
            "published_constructor_parameter_satisfies_profile": True,
            "local_power_savings_proxy": sum(
                (condition.forced_discriminant_valuation - 1)
                * log(condition.prime)
                for condition in search_conditions()
            ),
        },
        "neighborhood": {
            "definition": (
                "every nonzero integer coefficient pair (i,j) in [-R,R]^2 "
                "in the exact Gauss-reduced basis, quotiented by primitive "
                "normalization and simultaneous sign; zero-denominator, "
                "nonunit-denominator, T=0, and singular specializations removed"
            ),
            "coefficient_radius": args.coefficient_radius,
            "raw_nonzero_coefficient_pairs": (2 * args.coefficient_radius + 1) ** 2 - 1,
            "distinct_primitive_nonsingular_parameters": len(candidates),
            "minimum_height": min(candidate.height for candidate in candidates),
            "maximum_height": max(candidate.height for candidate in candidates),
            "published_parameter_in_neighborhood": any(
                candidate.parameter == RANK21_CONSTRUCTOR_PARAMETER
                for candidate in candidates
            ),
        },
        "selection_protocol": {
            "score": (
                "sum over good numerical primes 5<=p<=cutoff of "
                "((2-a_p)/(p+1-a_p))*log(p)"
            ),
            "stages": list(args.stages),
            "keep_counts": list(args.keep_counts),
            "leakage_control": (
                "each score cutoff sees only survivors retained at the prior "
                "cutoff; conductor and root-number work is final-stage only"
            ),
            "rank_work": "none",
        },
        "stages": stages,
        "final_conductor_records": conductor_records,
        "conductor_errors": conductor_errors,
        "published_record_calibration": calibration,
        "summary": {
            "population": len(candidates),
            "finalists": len(survivors),
            "conductor_calls_completed": len(conductor_records),
            "conductor_timeouts_or_errors": len(conductor_errors),
            "below_strict_log_conductor_target_without_rank_claim": sum(
                record["below_strict_log_conductor_target"]
                for record in conductor_records
            ),
            "finalists_beating_published_rank21_score": sum(
                record["score_exceeds_published_rank21_at_final_cutoff"]
                for record in conductor_records
            ),
            "target_hits": 0,
        },
        "bounds": {
            "score_timeout_seconds_per_batch": args.score_timeout,
            "conductor_timeout_seconds_per_candidate": args.conductor_timeout,
            "calibration_timeout_seconds": args.calibration_timeout,
            "pari_stack_bytes": args.stack_bytes,
        },
        "software": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "pari_gp": pari_version(),
        },
        "reproducing_command": REPRODUCING_COMMAND,
        "invocation": " ".join(shlex.quote(part) for part in [sys.executable, *sys.argv]),
        "script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    print(
        f"profile T={residue} mod {modulus}, candidates={len(candidates)}, "
        f"conductor completions={len(conductor_records)}, "
        f"timeouts/errors={len(conductor_errors)}"
    )
    for record in conductor_records:
        print(
            f"T={record['constructor_parameter']} score={record['score']} "
            f"logN={record['curve']['log_conductor']} "
            f"root={record['curve']['root_number']}"
        )


if __name__ == "__main__":
    main()
