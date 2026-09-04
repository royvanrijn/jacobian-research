#!/usr/bin/env sage -python
"""Blind iterative half-lattice search on ICARM curve 385.

The input boundary is the already-frozen rank-29 holdout search: this program
reads only its curve-385 short model, its seventeen generic points, and points
that the generic-deepest 43 arm found without the public complement.  It never
imports the 29-point public fixture.

Starting from the three directions independently certified among those blind
discoveries, each round:

1. proves that the current basis generates the complete discovered subgroup,
   enlarging it by an exact primitive relation if necessary;
2. computes the canonical-height form of that basis;
3. lifts the old 43 deep classes through every nonzero quotient-bit word;
4. chooses the shortest representative in each lifted parity class, orders
   the charts by decreasing new canonical depth, minimizes and reduces the
   associated pointed quartic, and runs the fixed bounded point search;
5. incorporates every newly found independent or saturating point.

The run stops at the first stable round or at a declared pre-search limit.
Search misses are bounded negative results, never rank upper bounds.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from decimal import Decimal, getcontext
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from math import gcd, lcm
from pathlib import Path
import platform
import resource
import shutil
import sys
import time
from typing import Any, Iterable, Sequence

from fpylll import Enumeration, GSO, IntegerMatrix
from sage.all import EllipticCurve, GF, Matrix, QQ, ZZ, pari


ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"
BLIND = ART / "half_lattice_search_ablation_rank29_holdout_blind_v1.json"
ENGINE_SOURCE = ROOT / "elliptic-curves/cas/half_lattice_fake_descent_replay.sage"
OUTPUT = ART / "curve385_iterated_half_lattice_blind_v1.json"

EXPECTED_BLIND_SHA256 = "1ee832ce6ecebc0550c008f8a10ccc2d75e727dfe9d5625802624c160e7969e6"
LABEL = "curve385-rank29"
GENERIC_DIMENSION = 17
OLD_CLASS_COUNT = 43
OPERATIVE_SCALE = 1_000_000
AUDIT_SCALE = 100_000
CERTIFICATE_PRIME_BOUND = 1_000
RATIONAL_RELATION_DENOMINATOR_BOUNDS = (2, 4, 8, 16, 32, 64, 128, 256, 1024, 4096)

sys.path[:0] = [str(ROOT / "elliptic-curves"), str(ROOT / "elliptic-curves/cas")]

from mod2_reduction_independence import (  # noqa: E402
    combined_mod2_rank,
    find_mod2_reduction_certificate,
)
from search_nagao_u135_alternate_covers import relation_proposals  # noqa: E402
from search_nagao_u42_skew_height import (  # noqa: E402
    exact_linear_combination,
    short_multiply,
)


engine = SourceFileLoader("curve385_iterated_half_lattice_engine", str(ENGINE_SOURCE)).load_module()


Point = tuple[Fraction, Fraction]


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def canonical_text(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def canonical_hash(value: Any) -> str:
    return sha256(canonical_text(value).encode()).hexdigest()


def cpu_clock() -> float:
    own = resource.getrusage(resource.RUSAGE_SELF)
    children = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime + own.ru_stime + children.ru_utime + children.ru_stime


def point_record(point: Point) -> dict[str, str]:
    return {"x": str(point[0]), "y": str(point[1])}


def read_point(record: dict[str, str]) -> Point:
    return Fraction(record["x"]), Fraction(record["y"])


def canonical_point(point: Point) -> Point:
    opposite = point[0], -point[1]
    return min(point, opposite)


def point_identifier(point: Point) -> str:
    point = canonical_point(point)
    return f"{point[0]}|{point[1]}"


def point_sort_key(point: Point):
    return (
        max(abs(point[0].numerator).bit_length(), point[0].denominator.bit_length()),
        max(abs(point[1].numerator).bit_length(), point[1].denominator.bit_length()),
        point,
    )


def signature_record(signatures, column_count: int) -> dict[str, Any]:
    return {
        "column_count": column_count,
        "combined_rank": combined_mod2_rank(signatures, column_count),
        "prime_bound": CERTIFICATE_PRIME_BOUND,
        "signatures": [
            {
                "prime": row.prime,
                "group_order": row.group_order,
                "doubled_subgroup_order": row.doubled_subgroup_order,
                "quotient_dimension": row.quotient_dimension,
                "rows": [list(vector) for vector in row.rows],
            }
            for row in signatures
        ],
    }


class CosetOracle:
    """Floating fplll CVP decisions with exact rounded-form recomputation."""

    def __init__(self, gram: Sequence[Sequence[int]], degree: int = 2) -> None:
        self.gram = tuple(tuple(int(value) for value in row) for row in gram)
        self.dimension = len(self.gram)
        if not self.dimension or any(len(row) != self.dimension for row in self.gram):
            raise ValueError("the CVP Gram matrix must be nonempty and square")
        self.degree = int(degree)
        self.gso = GSO.Mat(
            IntegerMatrix.from_matrix(self.gram),
            gram=True,
            float_type="dd",
            update=True,
        )
        self.mu = tuple(
            tuple(self.gso.get_mu(i, j) if i > j else 0.0 for j in range(self.dimension))
            for i in range(self.dimension)
        )
        self.distance_bound = (
            (degree - 1) ** 2
            * sum(abs(value) for row in self.gram for value in row)
            / (degree * degree)
            + 1.0
        )

    def solve(self, residue: Sequence[int]) -> tuple[int, tuple[int, ...], float]:
        residue = tuple(int(value) for value in residue)
        if len(residue) != self.dimension or any(value not in (0, 1) for value in residue):
            raise ValueError("a parity residue has the wrong shape")
        target = [
            -(
                residue[i]
                + sum(residue[j] * self.mu[j][i] for j in range(i + 1, self.dimension))
            )
            / self.degree
            for i in range(self.dimension)
        ]
        solutions = Enumeration(self.gso).enumerate(
            0,
            self.dimension,
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
            residue[index] + self.degree * closest[index]
            for index in range(self.dimension)
        )
        norm = sum(
            representative[i] * self.gram[i][j] * representative[j]
            for i in range(self.dimension)
            for j in range(self.dimension)
        )
        error = abs(self.degree * self.degree * float(reported_distance) - norm)
        if error > 1.0e-6 or norm < 0:
            raise ArithmeticError(f"invalid CVP norm={norm}, error={error}")
        return norm, representative, error


def canonical_height_gram(model, basis: Sequence[Point]):
    getcontext().prec = 110
    pari.default("realprecision", 110)
    curve = pari(EllipticCurve(QQ, list(model)))
    raw = curve.ellheightmatrix([list(point) for point in basis])
    dimension = len(basis)
    gram = tuple(
        tuple(Decimal(str(raw[i, j])) for j in range(dimension))
        for i in range(dimension)
    )
    maximum_asymmetry = max(
        abs(gram[i][j] - gram[j][i])
        for i in range(dimension)
        for j in range(dimension)
    )
    if maximum_asymmetry > Decimal("1e-90"):
        raise ArithmeticError(f"canonical-height Gram is asymmetric by {maximum_asymmetry}")
    return gram, maximum_asymmetry


def quadratic_decimal(gram, vector: Sequence[int]) -> Decimal:
    return sum(
        Decimal(vector[i]) * gram[i][j] * Decimal(vector[j])
        for i in range(len(vector))
        for j in range(len(vector))
    )


def rounded_gram(gram, scale: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(int((value * Decimal(scale)).to_integral_value()) for value in row)
        for row in gram
    )


def exact_relations_in_chunks(
    model,
    basis: Sequence[Point],
    points: Sequence[Point],
    *,
    chunk_size: int,
    timeout_seconds: float,
    stack_bytes: int,
):
    answer = []
    for start in range(0, len(points), chunk_size):
        chunk = points[start : start + chunk_size]
        answer.extend(
            relation_proposals(
                model,
                basis,
                chunk,
                timeout=timeout_seconds,
                stack_bytes=stack_bytes,
            )
        )
    return tuple(answer)


def rational_relation(model, basis: Sequence[Point], point: Point):
    """Recover and exactly verify a small-denominator relation, if present."""

    pari.default("realprecision", 150)
    curve = pari(EllipticCurve(QQ, list(model)))
    raw_basis = [list(row) for row in basis]
    height = curve.ellheightmatrix(raw_basis)
    pairings = pari([curve.ellheight(list(row), list(point)) for row in basis]).Col()
    coordinates = height.matsolve(pairings)
    decimals = tuple(str(coordinates[index]) for index in range(len(basis)))
    for bound in RATIONAL_RELATION_DENOMINATOR_BOUNDS:
        fractions = tuple(Fraction(value).limit_denominator(bound) for value in decimals)
        denominator = lcm(*(value.denominator for value in fractions))
        numerators = tuple(int(value * denominator) for value in fractions)
        common = denominator
        for value in numerators:
            common = gcd(common, abs(value))
        if common > 1:
            denominator //= common
            numerators = tuple(value // common for value in numerators)
        replay = exact_linear_combination(model[3], basis, numerators)
        multiplied = short_multiply(model[3], point, denominator)
        if replay == multiplied:
            relation = numerators + (-denominator,)
            if gcd(*(abs(value) for value in relation)) != 1:
                raise ArithmeticError("an exact rational relation was not primitive")
            return {
                "basis_coefficients": list(numerators),
                "point_coefficient": -denominator,
                "denominator": denominator,
                "reconstruction_bound": bound,
                "height_coordinate_decimals": list(decimals),
                "exact_group_law_replay": True,
            }
    return None


def basis_from_primitive_relation(
    model,
    basis: Sequence[Point],
    point: Point,
    relation: Sequence[int],
):
    relation_matrix = Matrix(ZZ, [list(map(int, relation))])
    smith, left, right = relation_matrix.smith_form(transformation=True)
    if smith[0, 0] != 1 or left * relation_matrix * right != smith:
        raise ArithmeticError("primitive relation failed Smith completion")
    completion = right.inverse()
    if tuple(int(value) for value in completion.row(0)) != tuple(map(int, relation)):
        if tuple(int(value) for value in completion.row(0)) == tuple(-int(x) for x in relation):
            completion.rescale_row(0, -1)
        else:
            raise ArithmeticError("Smith completion did not retain the primitive relation")
    generators = tuple(basis) + (point,)
    new_basis = []
    for row in completion.rows()[1:]:
        new_point = exact_linear_combination(model[3], generators, tuple(map(int, row)))
        if new_point is None:
            raise ArithmeticError("a saturated basis vector became torsion")
        new_basis.append(new_point)
    if len(new_basis) != len(basis):
        raise ArithmeticError("saturation changed the rational rank")
    return tuple(new_basis), [[int(value) for value in row] for row in completion.rows()]


def classify_discovered_group(
    *,
    model,
    basis: Sequence[Point],
    discoveries: dict[Point, set[str]],
    relation_chunk_size: int,
    relation_timeout_seconds: float,
    stack_bytes: int,
):
    """Return a basis of the group generated by basis and all discoveries."""

    basis = tuple(basis)
    events = []
    while True:
        basis_signs = {
            canonical_point(point)
            for point in basis
        }
        pending = tuple(
            sorted(
                (point for point in discoveries if canonical_point(point) not in basis_signs),
                key=point_sort_key,
            )
        )
        proposals = exact_relations_in_chunks(
            model,
            basis,
            pending,
            chunk_size=relation_chunk_size,
            timeout_seconds=relation_timeout_seconds,
            stack_bytes=stack_bytes,
        )
        integral = []
        nonintegral = []
        for point, (relation, exact) in zip(pending, proposals):
            if exact:
                integral.append((point, relation))
            else:
                nonintegral.append(point)

        changed = False
        for point in nonintegral:
            trial = basis + (point,)
            signatures = find_mod2_reduction_certificate(
                model, trial, prime_bound=CERTIFICATE_PRIME_BOUND
            )
            rank = combined_mod2_rank(signatures, len(trial))
            if rank == len(trial):
                basis = trial
                events.append(
                    {
                        "type": "NEW_Q_INDEPENDENT_DIRECTION",
                        "point": point_record(point),
                        "sources": sorted(discoveries[point]),
                        "basis_rank_after": len(basis),
                        "finite_reduction_certificate": signature_record(
                            signatures, len(trial)
                        ),
                    }
                )
                changed = True
                break
        if changed:
            continue

        for point in nonintegral:
            recovered = rational_relation(model, basis, point)
            if recovered is None:
                continue
            relation = tuple(recovered["basis_coefficients"]) + (
                int(recovered["point_coefficient"]),
            )
            if recovered["denominator"] == 1:
                raise ArithmeticError("PARI missed an exactly replayed integral relation")
            old_basis = basis
            basis, completion = basis_from_primitive_relation(
                model, basis, point, relation
            )
            events.append(
                {
                    "type": "FINITE_INDEX_SATURATION_INSIDE_DISCOVERED_GROUP",
                    "point": point_record(point),
                    "sources": sorted(discoveries[point]),
                    "rank": len(basis),
                    "index_factor": recovered["denominator"],
                    "primitive_relation": list(relation),
                    "relation_recovery": recovered,
                    "unimodular_completion_rows": completion,
                    "old_basis_sha256": canonical_hash([point_record(row) for row in old_basis]),
                    "new_basis_sha256": canonical_hash([point_record(row) for row in basis]),
                }
            )
            changed = True
            break
        if changed:
            continue

        unresolved = [
            {
                "point": point_record(point),
                "sources": sorted(discoveries[point]),
                "reason": (
                    "no exact integral or bounded-denominator relation and no full "
                    "finite-reduction independence certificate"
                ),
            }
            for point in nonintegral
        ]
        integral_rows = [
            {
                "point": point_record(point),
                "sources": sorted(discoveries[point]),
                "coordinates": list(map(int, relation)),
            }
            for point, relation in integral
        ]
        return basis, {
            "status": (
                "PASS_BASIS_EQUALS_DISCOVERED_GROUP" if not unresolved else "UNKNOWN_UNCLASSIFIED_POINTS"
            ),
            "discovered_nonbasis_point_count": len(pending),
            "exact_integral_relation_count": len(integral_rows),
            "exact_integral_relations": integral_rows,
            "unresolved": unresolved,
            "events": events,
            "basis_rank": len(basis),
            "basis_sha256": canonical_hash([point_record(row) for row in basis]),
        }


def gf2_lift_data(model, basis: Sequence[Point], generic: Sequence[Point], args):
    proposals = exact_relations_in_chunks(
        model,
        basis,
        generic,
        chunk_size=args.relation_chunk_size,
        timeout_seconds=args.relation_timeout_seconds,
        stack_bytes=args.stack_bytes,
    )
    if not all(exact for unused_relation, exact in proposals):
        raise ArithmeticError("a generic point is not integral in the discovered-group basis")
    coordinates = tuple(tuple(map(int, relation)) for relation, unused_exact in proposals)
    field = GF(2)
    rows = Matrix(field, coordinates)
    if rows.rank() != GENERIC_DIMENSION:
        raise ArithmeticError("the old generic subgroup lost mod-2 dimension inside the new group")
    running = rows
    complement = []
    for index in range(len(basis)):
        unit = Matrix(field, 1, len(basis), lambda unused_row, column: int(column == index))
        trial = running.stack(unit)
        if trial.rank() > running.rank():
            complement.append(tuple(int(value) for value in unit.row(0)))
            running = trial
    if running.rank() != len(basis) or len(complement) != len(basis) - GENERIC_DIMENSION:
        raise ArithmeticError("failed to construct the quotient-bit complement")
    return coordinates, tuple(complement)


def add_mod2(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    return tuple((int(a) + int(b)) & 1 for a, b in zip(left, right))


def scale_mod2(bit: int, vector: Sequence[int]) -> tuple[int, ...]:
    return tuple((int(bit) & 1) * int(value) for value in vector)


def lifted_residues(
    old_masks: Sequence[int],
    generic_coordinates: Sequence[Sequence[int]],
    complement: Sequence[Sequence[int]],
):
    dimension = len(generic_coordinates[0])
    rows = []
    for old_mask in old_masks:
        base = (0,) * dimension
        for index, vector in enumerate(generic_coordinates):
            if (old_mask >> index) & 1:
                base = add_mod2(base, vector)
        for quotient_word in range(1, 1 << len(complement)):
            residue = base
            for index, vector in enumerate(complement):
                if (quotient_word >> index) & 1:
                    residue = add_mod2(residue, vector)
            rows.append((old_mask, quotient_word, residue))
    if len({row[2] for row in rows}) != len(rows):
        raise ArithmeticError("lifted half-classes collided")
    return rows


def rank_lifts(model, basis: Sequence[Point], old_masks, args):
    height_started = time.monotonic()
    gram, asymmetry = canonical_height_gram(model, basis)
    height_seconds = time.monotonic() - height_started
    generic_coordinates, complement = gf2_lift_data(model, basis, args.generic_points, args)
    raw = lifted_residues(old_masks, generic_coordinates, complement)
    runs = {}
    for scale in (AUDIT_SCALE, OPERATIVE_SCALE):
        oracle = CosetOracle(rounded_gram(gram, scale))
        rows = []
        maximum_error = 0.0
        for old_mask, quotient_word, residue in raw:
            unused_norm, representative, error = oracle.solve(residue)
            depth = quadratic_decimal(gram, representative) / 4
            rows.append((depth, old_mask, quotient_word, residue, representative))
            maximum_error = max(maximum_error, error)
        rows.sort(key=lambda row: (-row[0], row[1], row[2]))
        runs[scale] = (rows, maximum_error)
    operative = runs[OPERATIVE_SCALE][0]
    audit = runs[AUDIT_SCALE][0]
    audit_map = {(row[1], row[2]): row[4] for row in audit}
    return operative, {
        "canonical_height_seconds": height_seconds,
        "canonical_height_maximum_asymmetry": str(asymmetry),
        "operative_rounding_scale": OPERATIVE_SCALE,
        "audit_rounding_scale": AUDIT_SCALE,
        "maximum_cvp_distance_error": {
            str(scale): error for scale, (unused_rows, error) in runs.items()
        },
        "representative_disagreement_count": sum(
            audit_map[(row[1], row[2])] != row[4] for row in operative
        ),
        "priority_order_identical_between_scales": [
            (row[1], row[2]) for row in audit
        ] == [(row[1], row[2]) for row in operative],
        "generic_coordinate_rows_in_current_basis": [list(row) for row in generic_coordinates],
        "quotient_complement_rows_mod2": [list(row) for row in complement],
        "quotient_bit_count": len(complement),
        "all_lifts_including_zero_word": OLD_CLASS_COUNT * (1 << len(complement)),
        "nonzero_quotient_word_lifts": len(operative),
        "ranked_lifts_sha256": canonical_hash(
            [
                {
                    "old_mask": row[1],
                    "quotient_word": row[2],
                    "residue": list(row[3]),
                    "representative": list(row[4]),
                    "depth": str(row[0]),
                }
                for row in operative
            ]
        ),
    }


def load_blind_input():
    if digest(BLIND) != EXPECTED_BLIND_SHA256:
        raise ArithmeticError("the frozen blind holdout artifact changed")
    blind = json.loads(BLIND.read_text())
    if blind.get("status") != "PASS_BLIND_ABLATION_SEARCH":
        raise ArithmeticError("the blind holdout search is not complete")
    if blind["blindness_boundary"]["exceptional_point_fixture_loaded"] is not False:
        raise ArithmeticError("the source search crossed its blindness boundary")
    case = next(row for row in blind["results"] if row["label"] == LABEL)
    model = tuple(Fraction(value) for value in case["short_model"])
    generic = tuple(read_point(row) for row in case["generic_points"])
    arm = next(row for row in case["arms"] if row["id"] == "generic-deepest43")
    old_masks = tuple(map(int, arm["masks"]))
    if len(generic) != GENERIC_DIMENSION or len(old_masks) != OLD_CLASS_COUNT:
        raise ArithmeticError("the curve-385 frozen input dimensions changed")
    candidate_indices = tuple(map(int, arm["candidate_point_indices"]))
    candidates = tuple(
        canonical_point(read_point(case["candidate_points"][index]["point"]))
        for index in candidate_indices
    )
    sources = {
        point: {f"initial-deep43-candidate:{index}"}
        for point, index in zip(candidates, candidate_indices)
    }
    old_mask_set = set(old_masks)
    initial_chart_keys = set()
    for row in case["cover_records"]:
        if row["mask"] not in old_mask_set:
            continue
        base_point = exact_linear_combination(
            model[3], generic, tuple(map(int, row["specialized_representative"]))
        )
        if base_point is None:
            raise ArithmeticError("an old deep-43 representative became torsion")
        initial_chart_keys.add(point_identifier(base_point))
    return case, model, generic, old_masks, candidates, sources, initial_chart_keys


def choose_initial_basis(model, generic: Sequence[Point], candidates: Sequence[Point]):
    basis = tuple(generic)
    selections = []
    for point in candidates:
        trial = basis + (point,)
        signatures = find_mod2_reduction_certificate(
            model, trial, prime_bound=CERTIFICATE_PRIME_BOUND
        )
        rank = combined_mod2_rank(signatures, len(trial))
        if rank == len(trial):
            basis = trial
            selections.append(
                {
                    "point": point_record(point),
                    "basis_rank_after": len(basis),
                    "finite_reduction_certificate": signature_record(signatures, len(trial)),
                }
            )
    return basis, selections


def discovery_records(discoveries: dict[Point, set[str]]):
    return [
        {"point": point_record(point), "sources": sorted(sources)}
        for point, sources in sorted(discoveries.items(), key=lambda row: point_sort_key(row[0]))
    ]


def write_payload(output: Path, payload: dict[str, Any]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def build_payload(args, *, case, model, generic, old_masks, initial_selections, initial_saturation):
    return {
        "schema": "elliptic-curves.curve385-iterated-half-lattice-blind.v1",
        "status": "PARTIAL_CHECKPOINT",
        "blindness_boundary": {
            "public_rank29_fixture_loaded": False,
            "public_exceptional_point_coordinates_loaded": False,
            "input_is_prior_frozen_blind_search": True,
            "initial_points_restricted_to_generic_deepest43_arm": True,
        },
        "curve": {
            "label": LABEL,
            "short_model": [str(value) for value in model],
            "generic_points": [point_record(point) for point in generic],
            "generic_rank": GENERIC_DIMENSION,
        },
        "old_deep43": {
            "masks": list(old_masks),
            "masks_sha256": canonical_hash(list(old_masks)),
            "source_arm_cover_count": len(old_masks),
        },
        "initial_transition": {
            "source": "generic-deepest43 blind discoveries only",
            "selected_new_direction_count": len(initial_selections),
            "rank_before": GENERIC_DIMENSION,
            "rank_after": GENERIC_DIMENSION + len(initial_selections),
            "selections": initial_selections,
            "discovered_group_saturation": initial_saturation,
        },
        "declared_budget": {
            "height_bound_each_quartic": args.height_bound,
            "wall_timeout_seconds_each_quartic_including_minimize_reduce_search": args.timeout_seconds,
            "gp_stack_bytes_each_quartic": args.stack_bytes,
            "finite_reduction_prime_bound": CERTIFICATE_PRIME_BOUND,
            "relation_chunk_size": args.relation_chunk_size,
            "relation_timeout_seconds_each_chunk": args.relation_timeout_seconds,
            "maximum_quotient_bits": args.max_quotient_bits,
            "maximum_nonzero_lifts_per_iteration": args.max_planned_lifts,
            "maximum_iterations": args.max_iterations,
            "checkpoint_every_completed_searches": args.checkpoint_every,
            "retries": 0,
        },
        "selection_policy": {
            "old_class_count": OLD_CLASS_COUNT,
            "lift_words": "every nonzero word in a deterministic mod-2 complement to old M17",
            "representative": "shortest representative at rounded canonical-height scale 10^6",
            "priority": "decreasing actual decimal canonical depth, then old mask, then quotient word",
            "audit_scale": AUDIT_SCALE,
            "operative_scale": OPERATIVE_SCALE,
            "unchanged_exact_base_point_charts_are_not_rerun": True,
            "stable": "a complete iteration changes neither rank nor finite index in the discovered group",
        },
        "current_basis": [],
        "discoveries": [],
        "searched_base_point_keys": [],
        "iterations": [],
        "input_hashes": {
            relative(BLIND): digest(BLIND),
            relative(ENGINE_SOURCE): digest(ENGINE_SOURCE),
            relative(Path(__file__)): digest(Path(__file__)),
        },
        "generation": {
            "python": platform.python_version(),
            "sage": str(sys.modules.get("sage")),
            "pari": str(pari("default(parisizemax)")),
        },
        "claim_boundary": [
            "Every accepted point and group relation is checked by exact rational group law.",
            "Full finite-reduction rank certifies Q-linear independence; a deficient bounded rank does not certify dependence.",
            "Discovered-group saturation means equality with the group generated by returned blind points, not saturation in the unknown full E(Q).",
            "Canonical-height ordering and CVP decisions use high-precision numerical heights with two rounding scales.",
            "Every point-search miss is bounded by the declared height and wall limits.",
            "No exact-rank upper bound or completeness of the public rank-29 subgroup is claimed.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--height-bound", type=int, default=100_000)
    parser.add_argument("--timeout-seconds", type=float, default=15.0)
    parser.add_argument("--stack-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--relation-chunk-size", type=int, default=64)
    parser.add_argument("--relation-timeout-seconds", type=float, default=180.0)
    parser.add_argument("--max-quotient-bits", type=int, default=4)
    parser.add_argument("--max-planned-lifts", type=int, default=688)
    parser.add_argument("--max-iterations", type=int, default=8)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    args = parser.parse_args()
    if args.height_bound <= 0 or not 0 < args.timeout_seconds <= 60:
        raise SystemExit("invalid quartic search budget")
    if args.relation_chunk_size <= 0 or not 0 < args.relation_timeout_seconds <= 300:
        raise SystemExit("invalid relation budget")
    if args.max_quotient_bits < 1 or args.max_planned_lifts < OLD_CLASS_COUNT:
        raise SystemExit("invalid lift limits")
    if args.max_iterations < 1 or args.checkpoint_every < 1:
        raise SystemExit("invalid iteration/checkpoint limits")
    if shutil.which("gp") is None:
        raise SystemExit("PARI/GP executable 'gp' was not found")

    (
        case,
        model,
        generic,
        old_masks,
        candidates,
        discoveries,
        initial_chart_keys,
    ) = load_blind_input()
    args.generic_points = generic
    basis, initial_selections = choose_initial_basis(model, generic, candidates)
    basis, initial_saturation = classify_discovered_group(
        model=model,
        basis=basis,
        discoveries=discoveries,
        relation_chunk_size=args.relation_chunk_size,
        relation_timeout_seconds=args.relation_timeout_seconds,
        stack_bytes=args.stack_bytes,
    )
    if initial_saturation["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
        raise ArithmeticError("initial discovered group could not be saturated exactly")
    if len(basis) != 20:
        raise ArithmeticError(f"the blind deep-43 transition reached rank {len(basis)}, not 20")

    payload = build_payload(
        args,
        case=case,
        model=model,
        generic=generic,
        old_masks=old_masks,
        initial_selections=initial_selections,
        initial_saturation=initial_saturation,
    )
    searched_keys = set(initial_chart_keys)
    payload["current_basis"] = [point_record(point) for point in basis]
    payload["discoveries"] = discovery_records(discoveries)
    payload["searched_base_point_keys"] = sorted(searched_keys)
    write_payload(args.output, payload)
    print(
        f"C385ITER|transition=M17->M{len(basis)}|"
        f"initial_directions={len(initial_selections)}|saturation=PASS",
        flush=True,
    )

    for iteration_index in range(1, args.max_iterations + 1):
        rank_before = len(basis)
        quotient_bits = rank_before - GENERIC_DIMENSION
        planned_count = OLD_CLASS_COUNT * ((1 << quotient_bits) - 1)
        if quotient_bits > args.max_quotient_bits or planned_count > args.max_planned_lifts:
            payload["status"] = "STOPPED_AT_DECLARED_LIFT_LIMIT"
            payload["stop"] = {
                "reason": "next iteration exceeds a predeclared quotient-bit or lift-count limit",
                "basis_rank": rank_before,
                "quotient_bits": quotient_bits,
                "next_nonzero_lift_count": planned_count,
            }
            write_payload(args.output, payload)
            print(
                f"C385ITER|status=LIMIT|rank={rank_before}|bits={quotient_bits}|"
                f"planned={planned_count}",
                flush=True,
            )
            return

        iteration_started_wall = time.monotonic()
        iteration_started_cpu = cpu_clock()
        ranked, ranking = rank_lifts(model, basis, old_masks, args)
        if ranking["quotient_bit_count"] != quotient_bits or len(ranked) != planned_count:
            raise ArithmeticError("lift census has the wrong size")
        iteration = {
            "iteration": iteration_index,
            "status": "SEARCHING",
            "basis_rank_before": rank_before,
            "basis_before": [point_record(point) for point in basis],
            "basis_before_sha256": canonical_hash([point_record(point) for point in basis]),
            "ranking": ranking,
            "cover_records": [],
            "unchanged_previously_searched_chart_count": 0,
            "unchanged_previously_searched_chart_keys": [],
        }
        payload["iterations"].append(iteration)
        write_payload(args.output, payload)
        print(
            f"C385ITER|iteration={iteration_index}|rank={rank_before}|bits={quotient_bits}|"
            f"nonzero_lifts={planned_count}|status=START",
            flush=True,
        )

        searched_this_iteration = 0
        for priority, (depth, old_mask, quotient_word, residue, representative) in enumerate(ranked, 1):
            base_point = exact_linear_combination(model[3], basis, representative)
            if base_point is None:
                raise ArithmeticError("a lifted nonzero class produced the point at infinity")
            base_key = point_identifier(base_point)
            if base_key in searched_keys:
                iteration["unchanged_previously_searched_chart_count"] += 1
                iteration["unchanged_previously_searched_chart_keys"].append(base_key)
                continue
            started_cpu = cpu_clock()
            outcome = engine.run_quartic_search(
                mask=sum(int(bit) << index for index, bit in enumerate(residue)),
                representative=representative,
                short_model=model,
                generic_points=basis,
                height_bound=args.height_bound,
                timeout_seconds=args.timeout_seconds,
                stack_bytes=args.stack_bytes,
            )
            cpu_seconds = cpu_clock() - started_cpu
            searched_keys.add(base_key)
            source = f"iteration:{iteration_index}:priority:{priority}"
            for point in outcome.curve_points:
                point = canonical_point(point)
                discoveries.setdefault(point, set()).add(source)
            iteration["cover_records"].append(
                {
                    "priority": priority,
                    "old_mask": old_mask,
                    "old_hex": f"0x{old_mask:05x}",
                    "quotient_word": quotient_word,
                    "quotient_word_binary": f"{quotient_word:0{quotient_bits}b}",
                    "current_basis_residue": list(residue),
                    "canonical_depth": str(depth),
                    "representative": list(representative),
                    "base_point_key": base_key,
                    "cpu_seconds": cpu_seconds,
                    "search": outcome.record,
                }
            )
            searched_this_iteration += 1
            if searched_this_iteration % args.checkpoint_every == 0:
                payload["discoveries"] = discovery_records(discoveries)
                payload["searched_base_point_keys"] = sorted(searched_keys)
                write_payload(args.output, payload)
            print(
                f"C385ITER|iteration={iteration_index}|priority={priority}/{planned_count}|"
                f"old={old_mask:#07x}|qword={quotient_word}|"
                f"status={outcome.record['status']}|points={len(outcome.curve_points)}",
                flush=True,
            )

        iteration["status"] = "SEARCH_COMPLETE"
        iteration["searched_new_chart_count"] = len(iteration["cover_records"])
        iteration["bounded_complete_count"] = sum(
            row["search"]["status"] == "bounded_search_complete"
            for row in iteration["cover_records"]
        )
        iteration["timeout_count"] = sum(
            row["search"]["status"] == "bounded_search_timeout"
            for row in iteration["cover_records"]
        )
        iteration["pari_failure_count"] = sum(
            row["search"]["status"] == "pari_failure"
            for row in iteration["cover_records"]
        )
        old_basis_hash = canonical_hash([point_record(point) for point in basis])
        basis, saturation = classify_discovered_group(
            model=model,
            basis=basis,
            discoveries=discoveries,
            relation_chunk_size=args.relation_chunk_size,
            relation_timeout_seconds=args.relation_timeout_seconds,
            stack_bytes=args.stack_bytes,
        )
        if saturation["status"] != "PASS_BASIS_EQUALS_DISCOVERED_GROUP":
            iteration["status"] = "UNKNOWN_UNCLASSIFIED_DISCOVERIES"
            iteration["discovered_group_saturation"] = saturation
            payload["status"] = "STOPPED_FAIL_CLOSED_UNCLASSIFIED_DISCOVERIES"
            payload["current_basis"] = [point_record(point) for point in basis]
            payload["discoveries"] = discovery_records(discoveries)
            payload["searched_base_point_keys"] = sorted(searched_keys)
            write_payload(args.output, payload)
            return
        new_basis_hash = canonical_hash([point_record(point) for point in basis])
        changed = old_basis_hash != new_basis_hash
        iteration.update(
            {
                "status": "CLASSIFIED",
                "discovered_group_saturation": saturation,
                "basis_rank_after": len(basis),
                "basis_after": [point_record(point) for point in basis],
                "basis_after_sha256": new_basis_hash,
                "group_changed": changed,
                "new_independent_direction_count": sum(
                    event["type"] == "NEW_Q_INDEPENDENT_DIRECTION"
                    for event in saturation["events"]
                ),
                "finite_index_saturation_event_count": sum(
                    event["type"] == "FINITE_INDEX_SATURATION_INSIDE_DISCOVERED_GROUP"
                    for event in saturation["events"]
                ),
                "wall_seconds": time.monotonic() - iteration_started_wall,
                "cpu_seconds": cpu_clock() - iteration_started_cpu,
            }
        )
        payload["current_basis"] = [point_record(point) for point in basis]
        payload["discoveries"] = discovery_records(discoveries)
        payload["searched_base_point_keys"] = sorted(searched_keys)
        if not changed:
            payload["status"] = "PASS_STABLE_DISCOVERED_SUBGROUP"
            payload["stable_rank"] = len(basis)
            payload["stable_after_iteration"] = iteration_index
            payload["reproducing_command"] = (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elliptic-curves/cas/run_curve385_iterated_half_lattice_search.sage "
                f"--height-bound {args.height_bound} --timeout-seconds {args.timeout_seconds} "
                f"--stack-bytes {args.stack_bytes} --max-quotient-bits {args.max_quotient_bits} "
                f"--max-planned-lifts {args.max_planned_lifts}"
            )
            write_payload(args.output, payload)
            print(
                f"C385ITER|iteration={iteration_index}|status=STABLE|rank={len(basis)}|"
                f"output={relative(args.output)}",
                flush=True,
            )
            return
        write_payload(args.output, payload)
        print(
            f"C385ITER|iteration={iteration_index}|status=GROW|"
            f"rank={rank_before}->{len(basis)}|events={len(saturation['events'])}",
            flush=True,
        )

    payload["status"] = "STOPPED_AT_DECLARED_ITERATION_LIMIT"
    payload["stop"] = {
        "reason": "the predeclared iteration limit was reached before stability",
        "basis_rank": len(basis),
    }
    write_payload(args.output, payload)


if __name__ == "__main__":
    main()
