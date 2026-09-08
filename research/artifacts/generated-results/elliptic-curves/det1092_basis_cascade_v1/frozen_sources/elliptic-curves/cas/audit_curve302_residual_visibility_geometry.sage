#!/usr/bin/env sage-python
"""Retrospective half-lattice visibility audit for the 302 M17 plateau.

This is deliberately a known-point diagnostic.  It fixes a direct-summand
basis of the recovered and residual quotients in the public rank-31 group,
then studies how those points look in degree-two charts from the original
specialized M17 lattice.

For a target R and a parity p in M17/2M17, write q_p for the reduced centre
in p and choose z so that p+2z is closest to 2R in the rounded height metric.
The translated target R+(u-z), where q_p=p+2u, has

    2(R+(u-z))-q_p = 2R-(p+2z).

Thus it is the correct target-relative translate for the pointed chart C_qp.
The all-parity scan is only Babai reduction; a fixed 512-parity shortlist is
then re-solved by DD and MPFR CVP, and the 32 best such charts per target are
constructed and evaluated exactly.  Hence the output gives an exact minimum
coordinate *in that declared vetted set*, never a global coordinate minimum.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import importlib.machinery
import json
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path

from fpylll import Enumeration, FPLLL, GSO, IntegerMatrix
from sage.all import EllipticCurve, GF, QQ, ZZ, matrix, pari, vector
from sage.version import version as SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
ART = ROOT / "artifacts" / "generated-results" / "elliptic-curves"
LOCAL = ROOT / "artifacts" / "local" / "elliptic-curves" / "curve302-residual-visibility-geometry-v1"

PARENT = ART / "curve302_recovered_mw17_parent_v1.json"
SPAN = ART / "curve302_recovered_public_span_v1" / "result.json"
FILTRATION = ART / "curve302_recovered_quotient_local_filtration_v1.json"
DEGREE_TWO = ART / "curve302_parent_degree2_multisection_lattice_v1.json"
PUBLIC_SOURCE = CAS / "icarm_curve302.py"
MAPPER_SOURCE = CAS / "factor_free_pari_mapping.sage"
STORE_SOURCE = CAS / "research_runtime" / "store.py"
OUTPUT = ART / "curve302_residual_visibility_geometry_v1.json"

DIMENSION = 17
METRIC_SCALE = 1_000_000
PARITY_COUNT = (1 << DIMENSION) - 1
CVP_SHORTLIST = 512
EXACT_CHARTS = 32
CVP_PRECISION = 192

sys.path.insert(0, str(CAS))
import icarm_curve302 as curve  # noqa: E402
from alternate_quartic_covers import alternate_cover  # noqa: E402
from half_lattice_pointed_sieve import linear_combination  # noqa: E402
from research_runtime.store import checkpoint  # noqa: E402
from search_observability import primitive  # noqa: E402

mapper = importlib.machinery.SourceFileLoader("curve302_visibility_mapper", str(MAPPER_SOURCE)).load_module()


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.relative_to(ROOT))


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def put_new(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def parity(mask: int, dimension: int = DIMENSION):
    return vector(ZZ, [(mask >> index) & 1 for index in range(dimension)])


def parity_mask(row) -> int:
    return sum((int(value) & 1) << index for index, value in enumerate(row))


def as_ints(row) -> list[int]:
    return [int(value) for value in row]


def fraction_text(value) -> str:
    value = Fraction(value)
    return str(value)


def fraction_square_root(value: Fraction) -> Fraction:
    value = Fraction(value)
    if value < 0:
        raise ArithmeticError("a quartic value was negative")
    numerator, denominator = isqrt(value.numerator), isqrt(value.denominator)
    if numerator * numerator != value.numerator or denominator * denominator != value.denominator:
        raise ArithmeticError("an asserted quartic value was not a rational square")
    return Fraction(numerator, denominator)


class RoundedMetricCVP:
    """Cross-precision CVP for one integral rounded Gram matrix.

    The returned vectors are checked in the integral Gram form.  This proves
    consistency of the numerical CVP calls, not an exact canonical-height
    theorem: the input height pairing was rounded first by design.
    """

    def __init__(self, gram, float_type: str, precision: int) -> None:
        self.gram = matrix(ZZ, gram)
        self.dimension = self.gram.nrows()
        if self.gram.ncols() != self.dimension or not self.gram.is_positive_definite():
            raise ArithmeticError("the rounded height Gram is not positive definite")
        if float_type == "mpfr":
            FPLLL.set_precision(precision)
        self.gso = GSO.Mat(
            IntegerMatrix.from_matrix([as_ints(row) for row in self.gram.rows()]),
            gram=True,
            float_type=float_type,
            update=True,
        )
        self.mu = [
            [self.gso.get_mu(i, j) if i > j else 0.0 for j in range(self.dimension)]
            for i in range(self.dimension)
        ]

    def gso_target(self, target) -> list[float]:
        return [
            float(target[i])
            + sum(float(target[j]) * self.mu[j][i] for j in range(i + 1, self.dimension))
            for i in range(self.dimension)
        ]

    def babai(self, target):
        return vector(ZZ, [int(value) for value in self.gso.babai(self.gso_target(target), gso=True)])

    def closest(self, target):
        target = vector(QQ, target)
        trial = vector(ZZ, [int(round(float(value))) for value in target])
        difference = vector(QQ, trial) - target
        bound = difference * self.gram * difference + 1
        if bound <= 0:
            raise ArithmeticError("CVP trial bound was nonpositive")
        # Enumeration retains its previous target internally, so each CVP call
        # must get a fresh enumerator.  Reusing one silently turns later calls
        # into a stale-target query.
        reported, coordinates = Enumeration(self.gso).enumerate(
            0, self.dimension, float(bound), 0, target=self.gso_target(target)
        )[0]
        answer = vector(ZZ, [int(round(value)) for value in coordinates])
        if any(abs(value - integer) > 1.0e-7 for value, integer in zip(coordinates, answer)):
            raise ArithmeticError("CVP returned a nonintegral coefficient")
        exact_distance = (vector(QQ, answer) - target) * self.gram * (vector(QQ, answer) - target)
        tolerance = 1.0e-6 * max(1.0, abs(float(exact_distance)))
        if abs(float(reported) - float(exact_distance)) > tolerance:
            raise ArithmeticError(
                "CVP floating distance disagrees with the integral metric: "
                "reported={!r}, exact={!r}".format(reported, exact_distance)
            )
        return answer, exact_distance


def cross_precision_closest(dd: RoundedMetricCVP, mpfr: RoundedMetricCVP, target):
    left, left_distance = dd.closest(target)
    right, right_distance = mpfr.closest(target)
    if left_distance != right_distance:
        raise ArithmeticError("DD and MPFR CVP minima differ in the rounded metric")
    return min(left, right), left_distance


def build_targets(generic, recovered, filtration):
    """Make a deterministic 4-local + 3-strict + 7-residual basis.

    The residual complement is not canonical mathematically.  We first use
    the deterministic Smith complement and then add M24 vectors so every
    selected residual lift is strict modulo two.  This retains a direct
    integral complement while exposing the strict filtration honestly.
    """

    smith, left, right = recovered.smith_form()
    if list(map(int, smith.diagonal())) != [1] * 24 or left * recovered * right != smith:
        raise ArithmeticError("M24 was not a primitive rank-24 public sublattice")
    residual_raw = left.inverse()[:, 24:]
    if abs(recovered.augment(residual_raw).det()) != 1:
        raise ArithmeticError("Smith residual complement was not unimodular")

    strict = matrix(GF(2), filtration["local_filtration_mod_2"]["strict_kernel_public_words"])
    strict_recovered = matrix(
        GF(2), filtration["local_filtration_mod_2"]["M24_strict_kernel_public_words"]
    )
    recovered_mod2 = recovered.change_ring(GF(2))
    residual = []
    residual_adjustments = []
    for index in range(7):
        raw = residual_raw.column(index)
        choices = []
        for mask in range(1 << strict.nrows()):
            kernel_word = sum(
                (((mask >> row) & 1) * strict.row(row) for row in range(strict.nrows())),
                vector(GF(2), 31),
            )
            difference = kernel_word - raw.change_ring(GF(2))
            if difference in recovered_mod2.column_space():
                coefficient = recovered_mod2.solve_right(difference)
                choices.append((tuple(map(int, kernel_word)), tuple(map(int, coefficient))))
        if len(choices) != 8:
            raise ArithmeticError("a residual mod-two class did not have eight strict M24 lifts")
        kernel_word, adjustment = min(choices)
        lift = raw + recovered * vector(ZZ, adjustment)
        if lift.change_ring(GF(2)) != vector(GF(2), kernel_word):
            raise ArithmeticError("strict residual lift has wrong parity")
        residual.append(lift)
        residual_adjustments.append(list(adjustment))
    residual_matrix = matrix(ZZ, [as_ints(item) for item in residual]).transpose()
    if abs(recovered.augment(residual_matrix).det()) != 1:
        raise ArithmeticError("strict residual adjustment changed the integral complement")
    if any(item.change_ring(GF(2)) not in strict.row_space() for item in residual):
        raise ArithmeticError("a residual lift escaped the strict kernel")

    strict_quotient = []
    strict_full_coefficients = []
    for word in strict_recovered.rows():
        coefficient = recovered_mod2.solve_right(word)
        strict_quotient.append(vector(GF(2), coefficient[17:]))
        strict_full_coefficients.append(vector(ZZ, as_ints(coefficient)))
    if matrix(GF(2), strict_quotient).rank() != 3:
        raise ArithmeticError("recovered strict kernel did not give three quotient directions")
    quotient_columns = list(strict_quotient)
    span = matrix(GF(2), 7, 0)
    for item in quotient_columns:
        span = span.augment(matrix(GF(2), 7, 1, item))
    for index in range(7):
        candidate = vector(GF(2), [item == index for item in range(7)])
        if span.augment(matrix(GF(2), 7, 1, candidate)).rank() > span.rank():
            quotient_columns.append(candidate)
            span = span.augment(matrix(GF(2), 7, 1, candidate))
    if len(quotient_columns) != 7:
        raise ArithmeticError("could not complete the recovered quotient basis")
    # Place the four locally non-strict directions first, matching the stated
    # 4-local + 3-strict recovered filtration.
    adapted_columns = quotient_columns[3:] + quotient_columns[:3]
    adapted = matrix(ZZ, [as_ints(item) for item in adapted_columns]).transpose()
    if abs(adapted.det()) != 1:
        raise ArithmeticError("the recovered local/strict basis was not unimodular")
    recovered_targets = []
    for index in range(7):
        if index < 4:
            coefficient = vector(ZZ, [0] * 17 + as_ints(adapted.column(index)))
        else:
            # The strict public word can have a nonzero generic M17 component.
            # Retaining it is essential: deleting it changes the local signature.
            coefficient = strict_full_coefficients[index - 4]
        word = recovered * coefficient
        recovered_targets.append(word)
    if any(item.change_ring(GF(2)) in strict.row_space() for item in recovered_targets[:4]):
        raise ArithmeticError("a declared recovered local direction was strict")
    if any(item.change_ring(GF(2)) not in strict.row_space() for item in recovered_targets[4:]):
        raise ArithmeticError("a declared recovered strict direction was non-strict")

    total = generic.augment(matrix(ZZ, [as_ints(item) for item in recovered_targets]).transpose()).augment(residual_matrix)
    if abs(total.det()) != 1:
        raise ArithmeticError("the fourteen displayed quotient targets do not complement M17")
    targets = []
    for index, word in enumerate(recovered_targets):
        targets.append(
            {
                "id": "recovered-local-{:02d}".format(index + 1)
                if index < 4
                else "recovered-strict-{:02d}".format(index - 3),
                "cohort": "recovered_local" if index < 4 else "recovered_strict",
                "public_word": as_ints(word),
                "strict_mod_2": index >= 4,
            }
        )
    for index, word in enumerate(residual):
        targets.append(
            {
                "id": "residual-strict-{:02d}".format(index + 1),
                "cohort": "residual_strict",
                "public_word": as_ints(word),
                "strict_mod_2": True,
                "M24_adjustment_mod_2": residual_adjustments[index],
            }
        )
    return targets, {
        "M17_to_M24_to_D_integral_direct_sum": True,
        "recovered_quotient_basis": "deterministic unimodular 4-local plus 3-strict mod-two adaptation",
        "residual_quotient_basis": "Smith complement adjusted by M24 to strict mod-two lifts",
        "residual_strict_lift_count_per_mod_two_class": 8,
    }


def rounded_public_height_gram():
    getcontext().prec = 110
    pari.default("realprecision", 110)
    pari.allocatemem(256_000_000, 1_073_741_824, silent=True)
    E = EllipticCurve(QQ, list(curve.short_coefficients()))
    height = pari(E).ellheightmatrix([list(point) for point in curve.SHORT_POINTS], precision=384)
    decimals = [[Decimal(str(height[i, j])) for j in range(31)] for i in range(31)]
    asymmetry = max(abs(decimals[i][j] - decimals[j][i]) for i in range(31) for j in range(31))
    if asymmetry > Decimal("1e-90"):
        raise ArithmeticError("PARI height Gram was asymmetric")
    rounded = matrix(
        ZZ,
        31,
        31,
        [int((entry * Decimal(METRIC_SCALE)).to_integral_value()) for row in decimals for entry in row],
    )
    if not rounded.is_positive_definite():
        raise ArithmeticError("rounded public height Gram was not positive definite")
    return rounded, str(asymmetry)


def metric_setup(public_gram, generic):
    gram = generic.transpose() * public_gram * generic
    change = matrix(ZZ, pari(gram).qflllgram()).transpose()
    reduced = change * gram * change.transpose()
    if abs(change.det()) != 1 or not reduced.is_positive_definite():
        raise ArithmeticError("M17 rounded-height LLL reduction failed")
    return gram, change, reduced


def reduced_parity(original, inverse_change):
    answer = vector(ZZ, (matrix(ZZ, 1, DIMENSION, original) * inverse_change).row(0))
    return vector(ZZ, [int(value) % 2 for value in answer])


def target_projection(target, generic, public_gram, original_gram, change):
    cross = generic.transpose() * public_gram * (2 * target)
    continuous = original_gram.solve_right(cross)
    return vector(QQ, (matrix(QQ, 1, DIMENSION, continuous).row(0) * change.inverse()))


def approximate_scan(target, generic, public_gram, change, reduced_gram, inverse_change):
    """Babai-reduce every nonzero M17 parity; retain only a fixed shortlist."""

    projection = target_projection(target, generic, public_gram, generic.transpose() * public_gram * generic, change)
    oracle = RoundedMetricCVP(reduced_gram, "dd", CVP_PRECISION)
    entries = []
    for mask in range(1, 1 << DIMENSION):
        original = parity(mask)
        reduced_parity_word = reduced_parity(original, inverse_change)
        lattice_target = vector(QQ, [(projection[index] - reduced_parity_word[index]) / 2 for index in range(DIMENSION)])
        closest = oracle.babai(lattice_target)
        centre_reduced = reduced_parity_word + 2 * closest
        centre_original = vector(ZZ, (matrix(ZZ, 1, DIMENSION, centre_reduced) * change).row(0))
        delta = 2 * target - generic * centre_original
        value = int(delta * public_gram * delta)
        entries.append((value, mask))
    entries.sort()
    if len(entries) != PARITY_COUNT or len({mask for _, mask in entries}) != PARITY_COUNT:
        raise ArithmeticError("the exhaustive Babai parity scan lost a class")
    return projection, entries[:CVP_SHORTLIST], entries[0][0]


def exact_cvp_rows(target, generic, public_gram, original_gram, change, reduced_gram, inverse_change, shortlist):
    """Re-solve the declared shortlist, including its fixed chart centres."""

    projection = target_projection(target, generic, public_gram, original_gram, change)
    dd = RoundedMetricCVP(reduced_gram, "dd", CVP_PRECISION)
    mpfr = RoundedMetricCVP(reduced_gram, "mpfr", CVP_PRECISION)
    rows = []
    for _, mask in shortlist:
        original = parity(mask)
        p = reduced_parity(original, inverse_change)
        # q_p=p+2u is the fixed reduced chart centre for this parity.
        u, centre_distance = cross_precision_closest(dd, mpfr, -vector(QQ, p) / 2)
        centre_reduced = p + 2 * u
        centre_original = vector(ZZ, (matrix(ZZ, 1, DIMENSION, centre_reduced) * change).row(0))
        # p+2z is closest to 2R.  The chart target is R+(u-z).
        z, target_distance = cross_precision_closest(
            dd, mpfr, vector(QQ, [(projection[index] - p[index]) / 2 for index in range(DIMENSION)])
        )
        closest_original = vector(ZZ, (matrix(ZZ, 1, DIMENSION, p + 2 * z) * change).row(0))
        translation = vector(ZZ, (matrix(ZZ, 1, DIMENSION, u - z) * change).row(0))
        delta = 2 * target - generic * closest_original
        value = int(delta * public_gram * delta)
        if value <= 0:
            raise ArithmeticError("a quotient target unexpectedly fell in M17")
        rows.append(
            {
                "parity": mask,
                "centre_m17_word": as_ints(centre_original),
                "closest_to_double_target_m17_word": as_ints(closest_original),
                "target_translation_m17_word": as_ints(translation),
                "delta_public_word": as_ints(delta),
                "centred_metric_numerator": value,
                "centred_height_in_rounded_metric": fraction_text(Fraction(value, 4 * METRIC_SCALE)),
                "centre_distance_in_reduced_metric": fraction_text(centre_distance),
                "target_distance_in_reduced_metric": fraction_text(target_distance),
            }
        )
    rows.sort(key=lambda item: (item["centred_metric_numerator"], item["parity"]))
    if len(rows) != CVP_SHORTLIST or len({row["parity"] for row in rows}) != CVP_SHORTLIST:
        raise ArithmeticError("the exact CVP shortlist was malformed")
    return rows


def degree_two_label(parity_word, degree_change, degree_oracle_dd, degree_oracle_mpfr):
    reduced = reduced_parity(parity_word, degree_change.inverse())
    u, distance = cross_precision_closest(degree_oracle_dd, degree_oracle_mpfr, -vector(QQ, reduced) / 2)
    representative = reduced + 2 * u
    norm = int(representative * degree_oracle_dd.gram * representative)
    mask = parity_mask(reduced)
    if norm == 10:
        category = "rational_bisection"
    elif mask and norm >= 8 and norm % 4 == 0:
        category = "genus_one_bisection_lattice_candidate"
    else:
        category = "outside_section_nonnegative_degree_two_survivors"
    return {
        "M_over_2M_orbit_mask_in_degree_two_LLL_coordinates": mask,
        "minimum_generic_MW17_norm": norm,
        "category": category,
        "minimum_distance_in_generic_reduced_metric": fraction_text(distance),
    }


def exact_chart_row(target, row, generic, basis_points, degree_change, degree_oracle_dd, degree_oracle_mpfr):
    """Construct C_q and evaluate the prescribed translated target exactly."""

    centre_word = vector(ZZ, row["centre_m17_word"])
    translation = vector(ZZ, row["target_translation_m17_word"])
    target_word = target + generic * translation
    short_model = curve.short_coefficients()
    centre = linear_combination(short_model, basis_points, centre_word)
    target_point = linear_combination(short_model, curve.SHORT_POINTS, target_word)
    if centre is None or target_point is None:
        raise ArithmeticError("a nontrivial quotient chart acquired an infinite endpoint")
    raw_cover = alternate_cover(short_model, centre)
    raw_parameter, raw_ordinate = raw_cover.curve_point_to_cover(target_point)
    chart = mapper.mapping(
        short_model,
        basis_points,
        {"representative": row["centre_m17_word"]},
    )
    a, b, c, d = [Fraction(value) for value in chart["matrix"]]
    numerator, denominator = primitive(d * raw_parameter - b, -c * raw_parameter + a)
    coefficients = [Fraction(value) for value in chart["discriminant_quartic"]]
    value = sum(coefficient * numerator**i * denominator ** (4 - i) for i, coefficient in enumerate(coefficients))
    transformed_denominator = c * numerator + d * denominator
    ratio = Fraction(chart["square_ratio"])
    expected_root = abs(raw_ordinate * transformed_denominator**2 / fraction_square_root(ratio))
    if value != expected_root * expected_root:
        raise ArithmeticError("the reduced pointed-quartic coordinate failed its exact square check")
    at_infinity = denominator == 0
    height = None if at_infinity else max(abs(numerator), denominator)
    label = degree_two_label(parity(row["parity"]), degree_change, degree_oracle_dd, degree_oracle_mpfr)
    return {
        **row,
        "translated_target_public_word": as_ints(target_word),
        "exact_target_point_short_model": [fraction_text(value) for value in target_point],
        "exact_centre_point_short_model": [fraction_text(value) for value in centre],
        "degree_two_translation_orbit": label,
        "raw_pointed_quartic_coordinate": [fraction_text(raw_parameter), fraction_text(raw_ordinate)],
        "reduced_quartic_coefficients_ascending": chart["discriminant_quartic"],
        "reduced_to_raw_horizontal_matrix": chart["matrix"],
        "reduced_coordinate": [str(numerator), str(denominator)],
        "reduced_quartic_square_root_absolute": fraction_text(expected_root),
        "at_reduced_coordinate_infinity": at_infinity,
        "reduced_coordinate_height": None if height is None else str(height),
        "reduced_coordinate_decimal_digits": None if height is None else len(str(height)),
        "inside_historical_height_125000": True if at_infinity else height <= 125000,
    }


def build_target(target_spec, generic, basis_points, public_gram, original_gram, change, reduced_gram, inverse_change, degree_change, degree_dd, degree_mpfr):
    target = vector(ZZ, target_spec["public_word"])
    projection, shortlist, babai_minimum = approximate_scan(
        target, generic, public_gram, change, reduced_gram, inverse_change
    )
    exact_rows = exact_cvp_rows(
        target, generic, public_gram, original_gram, change, reduced_gram, inverse_change, shortlist
    )
    charts = [
        exact_chart_row(target, row, generic, basis_points, degree_change, degree_dd, degree_mpfr)
        for row in exact_rows[:EXACT_CHARTS]
    ]
    finite_charts = [row for row in charts if not row["at_reduced_coordinate_infinity"]]
    best_coordinate = min(finite_charts, key=lambda row: (int(row["reduced_coordinate_height"]), row["parity"])) if finite_charts else None
    best_metric = exact_rows[0]
    return {
        **target_spec,
        "all_nonzero_M17_parities_babai_scanned": PARITY_COUNT,
        "Babai_scan_best_centred_metric_numerator": babai_minimum,
        "CVP_vetted_parities": CVP_SHORTLIST,
        "exact_pointed_quartics_constructed": EXACT_CHARTS,
        "best_CVP_centred_height_in_rounded_metric": best_metric["centred_height_in_rounded_metric"],
        "best_CVP_parity": best_metric["parity"],
        "projective_infinity_rows_in_vetted_set": sum(row["at_reduced_coordinate_infinity"] for row in charts),
        "minimum_finite_reduced_coordinate_height_in_vetted_set": None if best_coordinate is None else best_coordinate["reduced_coordinate_height"],
        "minimum_finite_reduced_coordinate_decimal_digits_in_vetted_set": None if best_coordinate is None else best_coordinate["reduced_coordinate_decimal_digits"],
        "minimum_finite_coordinate_parity_in_vetted_set": None if best_coordinate is None else best_coordinate["parity"],
        "CVP_rows": exact_rows,
        "exact_pointed_quartic_rows": charts,
    }


def summarise(rows):
    groups = {}
    for row in rows:
        groups.setdefault(row["cohort"], []).append(row)
    answer = {}
    for cohort, entries in sorted(groups.items()):
        finite = [item for item in entries if item["minimum_finite_reduced_coordinate_height_in_vetted_set"] is not None]
        digits = sorted(item["minimum_finite_reduced_coordinate_decimal_digits_in_vetted_set"] for item in finite)
        heights = sorted(int(item["minimum_finite_reduced_coordinate_height_in_vetted_set"]) for item in finite)
        answer[cohort] = {
            "directions": len(entries),
            "directions_with_projective_infinity_in_vetted_set": sum(item["projective_infinity_rows_in_vetted_set"] > 0 for item in entries),
            "finite_coordinate_directions": len(finite),
            "minimum_finite_coordinate_decimal_digits_range_in_vetted_sets": None if not digits else [digits[0], digits[-1]],
            "smallest_finite_coordinate_height_in_vetted_sets": None if not heights else str(heights[0]),
            "largest_of_directionwise_minimum_finite_coordinate_heights": None if not heights else str(heights[-1]),
        }
    return answer


def protocol(bindings):
    return {
        "schema": "elliptic-curves.curve302-residual-visibility-geometry-protocol.v1",
        "bindings": bindings,
        "metric": {
            "source": "PARI ellheightmatrix at 384 bits, rounded entrywise",
            "scale": METRIC_SCALE,
            "height_identity": "hhat(R+m-q/2)=hhat(2(R+m)-q)/4",
            "warning": "The rounded metric is an exact finite decision metric, not an exact canonical-height pairing.",
        },
        "parity_schedule": {
            "M17_over_2M17_nonzero_classes": PARITY_COUNT,
            "all_class_stage": "DD Babai reduction only; its order is a deterministic shortlist screen",
            "cross_precision_CVP_shortlist_per_target": CVP_SHORTLIST,
            "exact_pointed_quartic_rows_per_target": EXACT_CHARTS,
            "CVP_precision_bits": CVP_PRECISION,
        },
        "boundary": [
            "Every target is a known public rank-31 point combination, used only retrospectively.",
            "The output does not claim a global minimum reduced coordinate: unvetted parity classes can have a better coordinate because a PGL2 reduction has chart-dependent distortion.",
            "The degree-two M/2M bijection labels divisor translation orbits. It does not construct all bisection equations or identify a split fibre.",
            "The residual seven are a deterministic strict integral complement, not a canonical arithmetic splitting of D/M24.",
        ],
    }


def bindings():
    return {
        relative(path): digest(path)
        for path in (PARENT, SPAN, FILTRATION, DEGREE_TWO, PUBLIC_SOURCE, MAPPER_SOURCE, STORE_SOURCE, Path(__file__))
    }


def build(*, reuse_checkpoint: bool):
    parent, span, filtration, degree_two = load(PARENT), load(SPAN), load(FILTRATION), load(DEGREE_TWO)
    generic = matrix(ZZ, parent["basis_embedding_in_public_D"])
    recovered = matrix(ZZ, [relation["word"] for relation in span["relations"]]).transpose()
    if generic.dimensions() != (31, 17) or recovered.dimensions() != (31, 24) or recovered[:, :17] != generic:
        raise ArithmeticError("the M17/M24 public-word chain changed")
    if span["status"] != filtration["status"] != "PASS":
        raise ArithmeticError("the public-span or local-filtration certificate is unavailable")
    if degree_two["status"] != "PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT":
        raise ArithmeticError("the complete degree-two M/2M enumeration is unavailable")
    current_bindings = bindings()
    frozen_protocol = protocol(current_bindings)
    LOCAL.mkdir(parents=True, exist_ok=True)
    if (LOCAL / "protocol.json").exists():
        if load(LOCAL / "protocol.json") != frozen_protocol:
            raise ArithmeticError("the local checkpoint belongs to a different protocol")
    else:
        checkpoint(LOCAL / "protocol.json", frozen_protocol)

    public_gram, asymmetry = rounded_public_height_gram()
    original_gram, change, reduced_gram = metric_setup(public_gram, generic)
    inverse_change = change.inverse()
    degree_gram = matrix(ZZ, parent["generic_height_gram"])
    degree_change = matrix(ZZ, degree_two["enumeration"]["row_LLL_change_to_certified_parent_basis"])
    if degree_change * degree_gram * degree_change.transpose() != matrix(ZZ, degree_two["enumeration"]["reduced_gram"]):
        raise ArithmeticError("degree-two parent lattice change of basis changed")
    degree_reduced = degree_change * degree_gram * degree_change.transpose()
    degree_dd = RoundedMetricCVP(degree_reduced, "dd", CVP_PRECISION)
    degree_mpfr = RoundedMetricCVP(degree_reduced, "mpfr", CVP_PRECISION)
    basis_points = tuple(
        linear_combination(curve.short_coefficients(), curve.SHORT_POINTS, generic.column(index))
        for index in range(DIMENSION)
    )
    if any(point is None for point in basis_points):
        raise ArithmeticError("a displayed M17 specialization basis point is infinite")
    targets, target_boundary = build_targets(generic, recovered, filtration)

    metric_record = {
        "public_rounded_height_gram": [as_ints(row) for row in public_gram.rows()],
        "maximum_unrounded_height_gram_asymmetry": asymmetry,
        "M17_rounded_height_gram": [as_ints(row) for row in original_gram.rows()],
        "M17_LLL_change": [as_ints(row) for row in change.rows()],
        "M17_reduced_rounded_height_gram": [as_ints(row) for row in reduced_gram.rows()],
    }
    rows_path = LOCAL / "rows.json"
    completed = load(rows_path) if reuse_checkpoint and rows_path.exists() else []
    if any(row["id"] != targets[index]["id"] for index, row in enumerate(completed)):
        raise ArithmeticError("the local row checkpoint is out of protocol order")
    for index in range(len(completed), len(targets)):
        completed.append(
            build_target(
                targets[index], generic, basis_points, public_gram, original_gram, change, reduced_gram, inverse_change,
                degree_change, degree_dd, degree_mpfr,
            )
        )
        if reuse_checkpoint:
            checkpoint(rows_path, completed)
        print("CURVE302VISIBILITY|completed_target={}/{}|{}".format(index + 1, len(targets), targets[index]["id"]), flush=True)
    if len(completed) != 14:
        raise ArithmeticError("the retrospective target panel was incomplete")
    return {
        "schema": "elliptic-curves.curve302-residual-visibility-geometry.v1",
        "status": "PASS_RETROSPECTIVE_VETTED_VISIBILITY_DIAGNOSTIC",
        "bindings": current_bindings,
        "protocol": frozen_protocol,
        "target_basis_boundary": target_boundary,
        "metric": metric_record,
        "degree_two_half_lattice_correspondence": {
            "same_finite_set": "The degree-two divisor translation action w -> w+2x and the pointed-chart centre action Q -> Q+2x are both M17/2M17.",
            "basis_transport": "For a chart parity p in the certified parent M17 basis, its degree-two orbit mask is mask(p * B^{-1} mod 2), where B is the retained degree-two row LLL change to the parent basis.",
            "cardinality": 1 << DIMENSION,
            "rational_bisection_orbits": degree_two["rational_bisections"]["translation_orbits"],
            "rational_bisection_shell": "minimum generic M17 norm 10",
            "interpretation": "This is an exact label correspondence of torsors. The short rational-bisection shell is not thereby a deep-centre or visibility theorem.",
        },
        "directions": completed,
        "cohort_summary": summarise(completed),
        "reproducing_command": "sage -python elliptic-curves/cas/audit_curve302_residual_visibility_geometry.sage --check",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.build == args.check:
        parser.error("choose exactly one of --build or --check")
    result = build(reuse_checkpoint=args.build)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.build:
        if OUTPUT.exists():
            raise FileExistsError("refusing to overwrite an immutable visibility certificate")
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        OUTPUT.write_text(rendered)
    elif OUTPUT.read_text() != rendered:
        raise ArithmeticError("the stored visibility diagnostic did not replay")
    print("CURVE302VISIBILITY|targets=14|status=PASS_RETROSPECTIVE_VETTED_VISIBILITY_DIAGNOSTIC", flush=True)


if __name__ == "__main__":
    main()
