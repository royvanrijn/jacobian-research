#!/usr/bin/env sage -python
"""Blind half-lattice / fake-2-descent replay on the published R17 fibre.

The search half of this program deliberately imports neither the public rank-28
point fixture nor an artifact derived from those points.  It reconstructs only
the published R17 equation, the seventeen published sections, and their
specialization at ``t=-9529/5471``.  The finite search is:

    parity class -> specialized shortest representative -> quartic chart
      -> integral model -> PARI minimization -> Cremona--Stoll reduction
      -> modular square-density audit -> hyperellratpoints -> E(Q).

The quartic attached to ``P=(x_P,y_P)`` is

    w^2 = m^4 - 6*x_P*m^2 - 8*y_P*m - 3*x_P^2 - 4*A.

It is a pointed genus-one model birational to ``E``, not an everywhere
unpointed Selmer torsor.  Thus local solubility is automatic; the recorded
finite-prime data measure the square-sieve density for nontrivial affine
search coordinates instead.

This is a bounded point search and not a complete 2-descent or rank upper
bound.  Run ``verify_half_lattice_fake_descent_replay.sage`` afterwards to
compare its frozen output with the verification-only exceptional-point fixture.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal, getcontext
from fractions import Fraction
from hashlib import sha256
import json
from math import lcm
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time
from typing import Any, Iterable, Sequence

import numpy as np
from fpylll import Enumeration, GSO, IntegerMatrix
from sage.all import EllipticCurve, QQ, pari, prime_range


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
CAS = ELLIPTIC / "cas"
sys.path[:0] = [str(ELLIPTIC), str(CAS)]

from alternate_quartic_covers import alternate_cover, point_on_short_curve  # noqa: E402
from ecsearch.q12o5867_specialization import (  # noqa: E402
    evaluate_projective_specialization,
    global_minimal_model_with_change,
    load_q12o5867_data,
    short_certificate_model,
)
from ek_k3 import rational_to_string  # noqa: E402
from elliptic_candidate_record import source_point_to_target  # noqa: E402
from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from search_nagao_u135_alternate_covers import relation_proposals  # noqa: E402
from search_nagao_u42_skew_height import exact_linear_combination  # noqa: E402


MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "half_lattice_fake_descent_rank28_blind_v1.json"
)
PARAMETER = (-9529, 5471)
DIMENSION = 17
GENERIC_DEEPEST_COUNT = 43
SPECIALIZED_DEEPEST_COUNT = 43
ROUNDING_SCALES = (10_000, 100_000, 1_000_000)
SQUARE_SIEVE_PRIMES = (7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)

# Theorem 4 of arXiv:2608.25406v1, in the published P1,...,P17 basis.
GENERIC_GRAM = (
    (4, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, 1),
    (-2, 4, 2, 1, 0, 1, 0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1),
    (-2, 2, 4, 2, 1, 0, 1, 2, 1, 1, 1, 0, 1, 1, 2, 1, 0),
    (-2, 1, 2, 4, 1, 0, 1, 1, 2, 1, 1, 1, 0, 1, 0, 1, 0),
    (-2, 0, 1, 1, 4, 2, 0, 1, 1, 0, 1, 1, 1, 2, 1, 0, -1),
    (-2, 1, 0, 0, 2, 4, 1, 0, 0, 1, 2, 1, 1, 1, 1, 0, -1),
    (-2, 0, 1, 1, 0, 1, 4, 0, 0, 1, 1, 0, 1, 1, 2, 2, -2),
    (-2, 1, 2, 1, 1, 0, 0, 4, 1, 2, 1, 0, 1, 1, 1, 1, 0),
    (-2, 1, 1, 2, 1, 0, 0, 1, 4, 2, 2, 1, 2, 2, 1, 1, -1),
    (-2, 1, 1, 1, 0, 1, 1, 2, 2, 4, 2, 0, 2, 1, 2, 1, -1),
    (-2, 1, 1, 1, 1, 2, 1, 1, 2, 2, 4, 0, 1, 2, 2, 1, -2),
    (-2, 2, 0, 1, 1, 1, 0, 0, 1, 0, 0, 4, 0, 1, 0, 1, 1),
    (-2, 1, 1, 0, 1, 1, 1, 1, 2, 2, 1, 0, 4, 2, 2, 1, -2),
    (-2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 1, 2, 4, 2, 1, -2),
    (-2, 1, 2, 0, 1, 1, 2, 1, 1, 2, 2, 0, 2, 2, 4, 1, -2),
    (-2, 1, 1, 1, 0, 0, 2, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0),
    (1, 1, 0, 0, -1, -1, -2, 0, -1, -1, -2, 1, -2, -2, -2, 0, 4),
)
EXPECTED_GENERIC_MINIMUM_HISTOGRAM = {0: 1, 4: 1311, 6: 26672, 8: 63925, 10: 39120, 12: 43}


def file_digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def point_record(point: tuple[Fraction, Fraction]) -> dict[str, str]:
    return {"x": rational_to_string(point[0]), "y": rational_to_string(point[1])}


def bit_height(value: Fraction | int) -> int:
    value = Fraction(value)
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def binary_rank(rows: Iterable[int]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        value = int(row)
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def restricted_signature_rank(signatures, columns: Sequence[int]) -> int:
    packed = []
    for signature in signatures:
        for row in signature.rows:
            packed.append(sum((int(row[column]) & 1) << offset for offset, column in enumerate(columns)))
    return binary_rank(packed)


class CosetOracle:
    """fplll floating CVP decisions with exact integral norm recomputation."""

    def __init__(self, gram: Sequence[Sequence[int]], degree: int = 2) -> None:
        self.gram = tuple(tuple(int(value) for value in row) for row in gram)
        self.degree = int(degree)
        self.gso = GSO.Mat(
            IntegerMatrix.from_matrix(self.gram),
            gram=True,
            float_type="dd",
            update=True,
        )
        self.mu = tuple(
            tuple(self.gso.get_mu(i, j) if i > j else 0.0 for j in range(DIMENSION))
            for i in range(DIMENSION)
        )
        self.distance_bound = (
            (degree - 1) ** 2
            * sum(abs(value) for row in self.gram for value in row)
            / (degree * degree)
            + 1.0
        )

    def solve(self, mask: int) -> tuple[int, tuple[int, ...], float]:
        residue = tuple((mask >> index) & 1 for index in range(DIMENSION))
        target = [
            -(
                residue[i]
                + sum(residue[j] * self.mu[j][i] for j in range(i + 1, DIMENSION))
            )
            / self.degree
            for i in range(DIMENSION)
        ]
        solutions = Enumeration(self.gso).enumerate(
            0,
            DIMENSION,
            self.distance_bound,
            0,
            target=target,
        )
        if not solutions:
            raise ArithmeticError("CVP enumeration returned no solution")
        reported_distance, coordinates = solutions[0]
        closest = tuple(int(round(value)) for value in coordinates)
        if any(abs(value - integer) > 1.0e-7 for value, integer in zip(coordinates, closest)):
            raise ArithmeticError("CVP enumeration returned nonintegral coordinates")
        representative = tuple(
            residue[index] + self.degree * closest[index] for index in range(DIMENSION)
        )
        norm = sum(
            representative[i] * self.gram[i][j] * representative[j]
            for i in range(DIMENSION)
            for j in range(DIMENSION)
        )
        error = abs(self.degree * self.degree * float(reported_distance) - norm)
        if error > 1.0e-6 or norm < 0:
            raise ArithmeticError(f"invalid CVP norm={norm}, error={error}")
        return norm, representative, error


def canonical_height_gram(short_model, generic_points) -> tuple[tuple[Decimal, ...], ...]:
    getcontext().prec = 110
    pari.default("realprecision", 110)
    curve = pari(EllipticCurve(QQ, list(short_model)))
    raw = curve.ellheightmatrix([list(point) for point in generic_points])
    gram = tuple(
        tuple(Decimal(str(raw[i, j])) for j in range(DIMENSION))
        for i in range(DIMENSION)
    )
    maximum_asymmetry = max(
        abs(gram[i][j] - gram[j][i]) for i in range(DIMENSION) for j in range(DIMENSION)
    )
    if maximum_asymmetry > Decimal("1e-90"):
        raise ArithmeticError(f"height Gram is unexpectedly asymmetric: {maximum_asymmetry}")
    return gram


def quadratic_decimal(gram, vector: Sequence[int]) -> Decimal:
    return sum(
        Decimal(vector[i]) * gram[i][j] * Decimal(vector[j])
        for i in range(DIMENSION)
        for j in range(DIMENSION)
    )


def rank_half_classes(specialized_gram) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    started = time.monotonic()
    generic_oracle = CosetOracle(GENERIC_GRAM)
    generic = []
    generic_histogram: Counter[int] = Counter()
    maximum_generic_error = 0.0
    for mask in range(1 << DIMENSION):
        norm, representative, error = generic_oracle.solve(mask)
        generic.append((norm, representative))
        generic_histogram[norm] += 1
        maximum_generic_error = max(maximum_generic_error, error)
    if dict(sorted(generic_histogram.items())) != EXPECTED_GENERIC_MINIMUM_HISTOGRAM:
        raise ArithmeticError(f"generic R17 CVP histogram changed: {generic_histogram}")

    rounded_runs = []
    for scale in ROUNDING_SCALES:
        rounded = tuple(
            tuple(int((value * Decimal(scale)).to_integral_value()) for value in row)
            for row in specialized_gram
        )
        oracle = CosetOracle(rounded)
        run = []
        maximum_error = 0.0
        for mask in range(1 << DIMENSION):
            norm, representative, error = oracle.solve(mask)
            actual = quadratic_decimal(specialized_gram, representative) / 4
            run.append((actual, representative, norm))
            maximum_error = max(maximum_error, error)
        rounded_runs.append((scale, run, maximum_error))

    final = rounded_runs[-1][1]
    classes = []
    for mask in range(1 << DIMENSION):
        generic_norm, generic_representative = generic[mask]
        actual, specialized_representative, unused = final[mask]
        classes.append(
            {
                "mask": mask,
                "hex": f"0x{mask:05x}",
                "generic_minimum_norm": generic_norm,
                "generic_depth": Decimal(generic_norm) / 4,
                "generic_representative": generic_representative,
                "specialized_depth": actual,
                "specialized_representative": specialized_representative,
            }
        )
    generic_order = sorted(classes, key=lambda row: (-row["generic_depth"], row["mask"]))
    specialized_order = sorted(classes, key=lambda row: (-row["specialized_depth"], row["mask"]))
    for rank, row in enumerate(generic_order, 1):
        row["generic_rank_with_mask_tiebreak"] = rank
    for rank, row in enumerate(specialized_order, 1):
        row["specialized_rank"] = rank

    stability = {}
    for size in (10, 20, 43, 50, 100, 500, 1000):
        sets = []
        for unused_scale, run, unused_error in rounded_runs:
            order = sorted(range(1 << DIMENSION), key=lambda mask: (-run[mask][0], mask))
            sets.append(set(order[:size]))
        stability[str(size)] = {
            "all_scale_intersection": len(set.intersection(*sets)),
            "pairwise_with_final": [len(current & sets[-1]) for current in sets[:-1]],
        }
    representative_disagreements = []
    for left in range(len(rounded_runs) - 1):
        earlier = rounded_runs[left][1]
        representative_disagreements.append(
            sum(earlier[mask][1] != final[mask][1] for mask in range(1 << DIMENSION))
        )

    generic_top = {row["mask"] for row in generic_order[:GENERIC_DEEPEST_COUNT]}
    specialized_top = {
        row["mask"] for row in specialized_order[:SPECIALIZED_DEEPEST_COUNT]
    }
    selected_masks = sorted(
        generic_top | specialized_top,
        key=lambda mask: (
            classes[mask]["specialized_rank"],
            classes[mask]["generic_rank_with_mask_tiebreak"],
        ),
    )
    summary = {
        "elapsed_seconds": time.monotonic() - started,
        "generic_minimum_norm_histogram": {
            str(key): value for key, value in sorted(generic_histogram.items())
        },
        "generic_maximum_cvp_distance_error": maximum_generic_error,
        "specialized_rounding_scales": list(ROUNDING_SCALES),
        "specialized_maximum_cvp_distance_errors": [run[2] for run in rounded_runs],
        "specialized_representative_disagreements_with_final": representative_disagreements,
        "top_set_stability": stability,
        "generic_deepest_count": len(generic_top),
        "specialized_deepest_count": len(specialized_top),
        "deepest_intersection_count": len(generic_top & specialized_top),
        "selected_union_count": len(selected_masks),
        "selected_masks": selected_masks,
        "specialized_depth_quantiles": {
            name: str(value)
            for name, value in zip(
                ("minimum", "q25", "median", "q75", "q90", "q99", "maximum"),
                np.quantile(
                    np.asarray([float(row["specialized_depth"]) for row in classes]),
                    (0, 0.25, 0.5, 0.75, 0.9, 0.99, 1),
                ),
            )
        },
    }
    return classes, summary


def gp_rational(value: Fraction) -> str:
    value = Fraction(value)
    return f"({value.numerator}/{value.denominator})"


def gp_polynomial(coefficients: Sequence[Fraction]) -> str:
    return "+".join(
        f"({Fraction(value).numerator})*x^{index}" for index, value in enumerate(coefficients)
    )


@dataclass(frozen=True)
class QuarticSearchResult:
    record: dict[str, Any]
    curve_points: tuple[tuple[Fraction, Fraction], ...]


def parse_gp_vector(text: str) -> list[str]:
    body = text.strip()
    if not body.startswith("[") or not body.endswith("]"):
        raise ValueError(f"not a GP vector: {text!r}")
    body = body[1:-1].strip()
    return [] if not body else [item.strip() for item in body.split(",")]


def modular_square_density(p_coefficients: Sequence[int], q_coefficients: Sequence[int], prime: int) -> dict[str, Any]:
    if prime == 2:
        raise ValueError("the simple square-density audit uses odd primes")
    count = 0
    for x_value in range(prime):
        p_value = sum((coefficient % prime) * pow(x_value, index, prime) for index, coefficient in enumerate(p_coefficients)) % prime
        q_value = sum((coefficient % prime) * pow(x_value, index, prime) for index, coefficient in enumerate(q_coefficients)) % prime
        discriminant = (q_value * q_value + 4 * p_value) % prime
        if discriminant == 0 or pow(discriminant, (prime - 1) // 2, prime) == 1:
            count += 1
    return {"prime": prime, "affine_x_survivors": count, "fraction": f"{count}/{prime}"}


def run_quartic_search(
    *,
    mask: int,
    representative: Sequence[int],
    short_model,
    generic_points,
    height_bound: int,
    timeout_seconds: float,
    stack_bytes: int,
) -> QuarticSearchResult:
    base_point = exact_linear_combination(
        Fraction(short_model[3]), generic_points, representative
    )
    if base_point is None:
        raise ArithmeticError("a nonzero selected class produced the point at infinity")
    cover = alternate_cover(short_model, base_point)
    denominator = 1
    for coefficient in cover.coefficients:
        denominator = lcm(denominator, Fraction(coefficient).denominator)
    integral_coefficients = tuple(
        Fraction(coefficient) * denominator * denominator for coefficient in cover.coefficients
    )
    if any(value.denominator != 1 for value in integral_coefficients):
        raise ArithmeticError("quartic denominator clearing failed")
    integral_coefficients = tuple(int(value) for value in integral_coefficients)
    polynomial = gp_polynomial(tuple(Fraction(value) for value in integral_coefficients))
    x_base, y_base = base_point
    program = f"""
C0=[{polynomial},0];
gettime(); C1=hyperellminimalmodel(C0,&m1); minms=gettime();
gettime(); C2=hyperellred(C1,&m2); redms=gettime();
print("MINMS|",minms);
print("REDMS|",redms);
print("REDP|",Vec(C2[1]));
print("REDQ|",Vec(C2[2]));
print("REDDISC|",hyperelldisc(C2));
gettime(); R=hyperellratpoints(C2,{height_bound}); searchms=gettime();
print("SEARCHMS|",searchms);
print("SEARCHCOUNT|",#R);
for(i=1,#R,p=R[i];z=m2[2][2,1]*p[1]+m2[2][2,2];if(z==0,print("TRIVIAL_INFINITY|stage=red"),p1=[(m2[2][1,1]*p[1]+m2[2][1,2])/z,(m2[1]*p[2]+subst(m2[3],x,p[1]))/z^2];z=m1[2][2,1]*p1[1]+m1[2][2,2];if(z==0,print("TRIVIAL_INFINITY|stage=min"),p0=[(m1[2][1,1]*p1[1]+m1[2][1,2])/z,(m1[1]*p1[2]+subst(m1[3],x,p1[1]))/z^2];if(!hyperellisoncurve(C0,p0),error("inverse model map left C0"));ex=(p0[1]^2-{gp_rational(x_base)}+p0[2]/{denominator})/2;ey=p0[1]*(ex-{gp_rational(x_base)})-{gp_rational(y_base)};print("POINT|",p0[1],"|",p0[2]/{denominator},"|",ex,"|",ey))));
quit
"""
    started = time.monotonic()
    try:
        completed = subprocess.run(
            ["gp", "-q", "-s", str(stack_bytes)],
            input=program,
            text=True,
            capture_output=True,
            timeout=timeout_seconds,
            check=False,
        )
        wall_seconds = time.monotonic() - started
    except subprocess.TimeoutExpired:
        return QuarticSearchResult(
            {
                "mask": mask,
                "hex": f"0x{mask:05x}",
                "status": "bounded_search_timeout",
                "height_bound": height_bound,
                "timeout_seconds": timeout_seconds,
                "wall_seconds": time.monotonic() - started,
                "representative": list(map(int, representative)),
                "base_point": point_record(base_point),
                "raw_quartic_coefficients_ascending": [rational_to_string(value) for value in cover.coefficients],
                "raw_rational_coefficient_maximum_bits": max(bit_height(value) for value in cover.coefficients),
                "denominator_clearing_factor_bits": denominator.bit_length(),
                "integral_model_maximum_coefficient_bits": max(abs(value).bit_length() for value in integral_coefficients),
            },
            (),
        )
    if completed.returncode != 0 or "***" in completed.stderr:
        return QuarticSearchResult(
            {
                "mask": mask,
                "hex": f"0x{mask:05x}",
                "status": "pari_failure",
                "height_bound": height_bound,
                "timeout_seconds": timeout_seconds,
                "wall_seconds": wall_seconds,
                "error": completed.stderr.strip()[-2000:],
                "representative": list(map(int, representative)),
                "base_point": point_record(base_point),
                "raw_quartic_coefficients_ascending": [rational_to_string(value) for value in cover.coefficients],
                "raw_rational_coefficient_maximum_bits": max(bit_height(value) for value in cover.coefficients),
                "denominator_clearing_factor_bits": denominator.bit_length(),
                "integral_model_maximum_coefficient_bits": max(abs(value).bit_length() for value in integral_coefficients),
            },
            (),
        )

    markers: dict[str, str] = {}
    curve_points = []
    raw_points = []
    infinity_count = 0
    for line in completed.stdout.splitlines():
        if line.startswith("POINT|"):
            unused, raw_x, raw_y, curve_x, curve_y = line.split("|", 4)
            raw_point = (Fraction(raw_x), Fraction(raw_y))
            curve_point = (Fraction(curve_x), Fraction(curve_y))
            if raw_point[1] ** 2 != cover.value(raw_point[0]):
                raise ArithmeticError("mapped PARI point left the raw quartic")
            if cover.cover_point_to_curve(raw_point) != curve_point:
                raise ArithmeticError("PARI/Python quartic maps disagree")
            if not point_on_short_curve(short_model, curve_point):
                raise ArithmeticError("mapped quartic point left E")
            raw_points.append(raw_point)
            curve_points.append(curve_point)
        elif line.startswith("TRIVIAL_INFINITY|"):
            infinity_count += 1
        elif "|" in line:
            key, value = line.split("|", 1)
            markers[key] = value.strip()
    for required in ("MINMS", "REDMS", "REDP", "REDQ", "REDDISC", "SEARCHMS", "SEARCHCOUNT"):
        if required not in markers:
            raise ArithmeticError(f"PARI omitted marker {required} for mask {mask:#x}")
    red_p_high = [int(value) for value in parse_gp_vector(markers["REDP"])]
    red_q_high = [int(value) for value in parse_gp_vector(markers["REDQ"])]
    red_p = list(reversed(red_p_high))
    red_q = list(reversed(red_q_high))
    local_profile = [
        modular_square_density(red_p, red_q, prime)
        for prime in SQUARE_SIEVE_PRIMES
    ]
    record = {
        "mask": mask,
        "hex": f"0x{mask:05x}",
        "status": "bounded_search_complete",
        "height_bound": height_bound,
        "timeout_seconds": timeout_seconds,
        "wall_seconds": wall_seconds,
        "representative": list(map(int, representative)),
        "base_point": point_record(base_point),
        "raw_quartic_coefficients_ascending": [rational_to_string(value) for value in cover.coefficients],
        "raw_rational_coefficient_maximum_bits": max(bit_height(value) for value in cover.coefficients),
        "denominator_clearing_factor_bits": denominator.bit_length(),
        "integral_model_maximum_coefficient_bits": max(abs(value).bit_length() for value in integral_coefficients),
        "minimalization_milliseconds": int(markers["MINMS"]),
        "reduction_milliseconds": int(markers["REDMS"]),
        "reduced_model": {
            "P_coefficients_ascending": red_p,
            "Q_coefficients_ascending": red_q,
            "maximum_coefficient_bits": max(abs(value).bit_length() for value in red_p + red_q),
            "discriminant": markers["REDDISC"],
            "discriminant_bits": abs(int(markers["REDDISC"])).bit_length(),
        },
        "local_stage": {
            "solubility_filter": "not_applicable_pointed_model_birational_to_E",
            "reason": "the monic raw quartic has rational points at infinity",
            "affine_modular_square_sieve_profile": local_profile,
            "joint_independent_density_product": str(
                np.prod([row["affine_x_survivors"] / row["prime"] for row in local_profile])
            ),
        },
        "search_milliseconds": int(markers["SEARCHMS"]),
        "signed_affine_points_reported": int(markers["SEARCHCOUNT"]),
        "trivial_points_mapping_to_infinity": infinity_count,
        "finite_raw_points": [point_record(point) for point in raw_points],
        "finite_curve_points": [point_record(point) for point in curve_points],
    }
    return QuarticSearchResult(record, tuple(curve_points))


def load_blind_r17_fibre():
    data = load_q12o5867_data(MODEL, SECTIONS)
    specialization = evaluate_projective_specialization(data, *PARAMETER)
    minimal_model, minimal_change, _ = global_minimal_model_with_change(specialization.model)
    minimal_points = tuple(
        source_point_to_target(point, minimal_change) for point in specialization.points
    )
    short_model, short_change = short_certificate_model(minimal_model)
    short_points = tuple(
        source_point_to_target(point, short_change) for point in minimal_points
    )
    if len(short_points) != DIMENSION or any(
        not point_on_short_curve(short_model, point) for point in short_points
    ):
        raise ArithmeticError("published R17 specialization did not yield seventeen points")
    signatures = find_mod2_reduction_certificate(short_model, short_points, prime_bound=500)
    if combined_mod2_rank(signatures, DIMENSION) != DIMENSION:
        raise ArithmeticError("the specialized generic subgroup failed mod-2 independence")
    return minimal_model, short_model, short_points, signatures


def normalized_certificate(payload: Any) -> Any:
    """Drop deliberately machine-dependent effort fields for ``--check``.

    The artifact records wall-clock effort because it is experimentally useful,
    but replay validity must not depend on scheduler noise or PARI millisecond
    counters.  Everything else, including all models, points, depths, ranks, and
    declared bounds, is compared exactly.
    """

    if isinstance(payload, dict):
        return {
            key: normalized_certificate(value)
            for key, value in payload.items()
            if key
            not in {
                "elapsed_seconds",
                "wall_seconds",
                "minimalization_milliseconds",
                "reduction_milliseconds",
                "search_milliseconds",
            }
        }
    if isinstance(payload, list):
        return [normalized_certificate(value) for value in payload]
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--rank-only", action="store_true")
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    args = parser.parse_args()
    if args.height_bound <= 0 or not 0 < args.timeout_seconds <= 60:
        raise SystemExit("height and timeout bounds must be positive; timeout is capped at 60 seconds")
    if shutil.which("gp") is None:
        raise SystemExit("PARI/GP executable 'gp' was not found")

    minimal_model, short_model, generic_points, generic_signatures = load_blind_r17_fibre()
    specialized_gram = canonical_height_gram(short_model, generic_points)
    classes, ranking_summary = rank_half_classes(specialized_gram)
    class_by_mask = {row["mask"]: row for row in classes}
    selected_masks = ranking_summary["selected_masks"]

    search_records = []
    discoveries: dict[tuple[Fraction, Fraction], set[int]] = {}
    if not args.rank_only:
        for index, mask in enumerate(selected_masks, 1):
            row = class_by_mask[mask]
            outcome = run_quartic_search(
                mask=mask,
                representative=row["specialized_representative"],
                short_model=short_model,
                generic_points=generic_points,
                height_bound=args.height_bound,
                timeout_seconds=args.timeout_seconds,
                stack_bytes=args.stack_bytes,
            )
            outcome.record.update(
                {
                    "generic_depth": str(row["generic_depth"]),
                    "generic_rank_with_mask_tiebreak": row["generic_rank_with_mask_tiebreak"],
                    "specialized_depth": str(row["specialized_depth"]),
                    "specialized_rank": row["specialized_rank"],
                }
            )
            search_records.append(outcome.record)
            for point in outcome.curve_points:
                discoveries.setdefault(point, set()).add(mask)
            print(
                f"HALFLATTICE|class={index}/{len(selected_masks)}|mask={mask:#07x}|"
                f"status={outcome.record['status']}|points={len(outcome.curve_points)}",
                flush=True,
            )

    basis_with_signs = {
        signed
        for point in generic_points
        for signed in (point, (point[0], -point[1]))
    }
    candidates = tuple(
        sorted(
            (point for point in discoveries if point not in basis_with_signs),
            key=lambda point: (bit_height(point[0]), bit_height(point[1]), point),
        )
    )
    proposals = relation_proposals(
        short_model,
        generic_points,
        candidates,
        timeout=min(60.0, max(15.0, args.timeout_seconds * 2)),
        stack_bytes=args.stack_bytes,
    ) if candidates else ()
    unexplained = tuple(
        point for point, (unused_relation, exact) in zip(candidates, proposals) if not exact
    )
    all_signatures = find_mod2_reduction_certificate(
        short_model, generic_points + unexplained, prime_bound=1000
    ) if unexplained else generic_signatures
    total_rank = combined_mod2_rank(all_signatures, DIMENSION + len(unexplained))
    selected_candidate_offsets = []
    selected_columns = list(range(DIMENSION))
    current_rank = restricted_signature_rank(all_signatures, selected_columns)
    for offset in range(len(unexplained)):
        trial = selected_columns + [DIMENSION + offset]
        rank = restricted_signature_rank(all_signatures, trial)
        if rank > current_rank:
            selected_candidate_offsets.append(offset)
            selected_columns = trial
            current_rank = rank
    if current_rank != total_rank:
        raise ArithmeticError("greedy finite-code basis did not reach the combined rank")

    candidate_rows = []
    unexplained_index = {point: index for index, point in enumerate(unexplained)}
    for point, (relation, exact) in zip(candidates, proposals):
        offset = unexplained_index.get(point)
        candidate_rows.append(
            {
                "point": point_record(point),
                "source_masks": sorted(discoveries[point]),
                "source_hex": [f"0x{mask:05x}" for mask in sorted(discoveries[point])],
                "exact_relation_in_generic_subgroup": exact,
                "generic_relation_if_exact": list(relation) if exact else None,
                "selected_for_blind_independent_quotient_basis": (
                    offset in selected_candidate_offsets if offset is not None else False
                ),
            }
        )

    selected_classes = []
    for mask in selected_masks:
        row = class_by_mask[mask]
        selected_classes.append(
            {
                "mask": mask,
                "hex": row["hex"],
                "generic_minimum_norm": row["generic_minimum_norm"],
                "generic_depth": str(row["generic_depth"]),
                "generic_rank_with_mask_tiebreak": row["generic_rank_with_mask_tiebreak"],
                "generic_representative": list(row["generic_representative"]),
                "specialized_depth": str(row["specialized_depth"]),
                "specialized_rank": row["specialized_rank"],
                "specialized_representative": list(row["specialized_representative"]),
            }
        )

    script_path = Path(__file__).resolve()
    payload = {
        "schema": "elliptic-curves.half-lattice-fake-descent-rank28-blind.v1",
        "status": (
            "PASS_BLIND_BOUNDED_RANK28_HALF_LATTICE_REPLAY"
            if not args.rank_only
            else "PASS_BLIND_RANK28_HALF_LATTICE_RANKING_ONLY"
        ),
        "blindness_boundary": {
            "search_loaded_exceptional_point_fixture": False,
            "forbidden_fixture": "elliptic-curves/cas/elkies_rank28.py",
            "forbidden_positive_control_artifact": "elkies_2026_high_rank_positive_controls_v2.json",
            "permitted_inputs": [str(MODEL.relative_to(ROOT)), str(SECTIONS.relative_to(ROOT))],
            "verification_is_separate": "elliptic-curves/cas/verify_half_lattice_fake_descent_replay.sage",
        },
        "input_hashes": {
            str(MODEL.relative_to(ROOT)): file_digest(MODEL),
            str(SECTIONS.relative_to(ROOT)): file_digest(SECTIONS),
            str(script_path.relative_to(ROOT)): file_digest(script_path),
        },
        "fibre": {
            "parameter": "-9529/5471",
            "minimal_model": [rational_to_string(Fraction(value)) for value in minimal_model],
            "short_model": [rational_to_string(Fraction(value)) for value in short_model],
            "generic_point_count": len(generic_points),
            "generic_points": [point_record(point) for point in generic_points],
            "generic_mod2_independence_rank": combined_mod2_rank(generic_signatures, DIMENSION),
            "generic_mod2_certificate_primes": [signature.prime for signature in generic_signatures],
            "specialized_canonical_height_gram": [
                [str(value) for value in row] for row in specialized_gram
            ],
        },
        "ranking": {
            **ranking_summary,
            "selected_classes": selected_classes,
            "generic_depth_status": "exact integral CVP; histogram matches the complete R17 M/2M census",
            "specialized_depth_status": (
                "complete CVP on three rounded canonical-height Grams; representatives and top sets "
                "are cross-scale numerical evidence, not interval-certified canonical-height CVP"
            ),
        },
        "pipeline": {
            "representative_policy": "shortest representative on the scale-10^6 rounded specialized canonical-height Gram",
            "quartic": "w^2=m^4-6*x_P*m^2-8*y_P*m-3*x_P^2-4*A",
            "integralization": "W=D*w with D=lcm of rational coefficient denominators",
            "minimization": "PARI hyperellminimalmodel",
            "reduction": "PARI hyperellred (Cremona--Stoll)",
            "local_filter": "pointed-model solubility is automatic; record affine modular square-sieve densities",
            "search": f"PARI hyperellratpoints on every selected reduced model through height {args.height_bound}",
            "map_back": "compose exact inverse hyperelliptic model changes, undo W=D*w, then apply the exact quartic-to-E map",
        },
        "declared_search_budget": {
            "rank_only": args.rank_only,
            "selected_class_count": len(selected_masks),
            "height_bound_each": args.height_bound,
            "timeout_seconds_each": args.timeout_seconds,
            "stack_bytes_each": args.stack_bytes,
            "single_pass_no_retry": True,
        },
        "search_records": search_records,
        "blind_results": {
            "distinct_nonbasis_candidates": len(candidates),
            "candidates_unexplained_by_exact_generic_group_law": len(unexplained),
            "combined_finite_mod2_rank": total_rank,
            "certified_quotient_gain": total_rank - DIMENSION,
            "certified_rank_lower_bound_from_blind_search": total_rank,
            "selected_independent_candidate_count": len(selected_candidate_offsets),
            "candidate_points": candidate_rows,
        },
        "software": {
            "python": platform.python_version(),
            "sage": subprocess.run([sys.executable, "-c", "import sage.env;print(sage.env.SAGE_VERSION)"], text=True, capture_output=True, check=False).stdout.strip(),
            "pari_gp": subprocess.run(["gp", "-fq"], input="print(version());quit\n", text=True, capture_output=True, check=False).stdout.strip(),
            "platform": platform.platform(),
        },
        "claim_boundary": [
            "All curve equations, group-law maps, returned quartic points, and finite-reduction independence checks are exact.",
            "The generic half-lattice depths are exact integral CVP results.",
            "The specialized canonical heights and their rounded-form CVP ordering are high-precision numerical evidence, not interval certificates.",
            "A hyperellratpoints miss is bounded by the declared reduced-coordinate height and per-model timeout.",
            "These pointed quartics are birational coordinate models of E, not nontrivial locally soluble 2-covering torsors.",
            "No rank upper bound or full Mordell--Weil saturation claim is made.",
        ],
        "reproducing_command": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
        ),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.exists():
            raise SystemExit(f"stale or missing output: {display_path(args.output)}")
        frozen = json.loads(args.output.read_text())
        if normalized_certificate(frozen) != normalized_certificate(payload):
            raise SystemExit(f"stale or invalid output: {display_path(args.output)}")
        print(f"HALFLATTICE|status=PASS_CHECK|output={display_path(args.output)}")
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(
        f"HALFLATTICE|status={payload['status']}|selected={len(selected_masks)}|"
        f"gain={payload['blind_results']['certified_quotient_gain']}|"
        f"output={display_path(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
