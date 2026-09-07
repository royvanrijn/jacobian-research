#!/usr/bin/env python3
"""Measure multisection diversity on the published rootless R17 frame.

Status: ACTIVE_SEARCH diagnostic with complete degree-two subcertificate.
Claim: computes the exact/sampled invariants and boundaries recorded below.
Inputs: pinned R17 Gram, complete bisection atlas, and complete d=3 census.
Output: artifacts/generated-results/elkies-k3-r17-multisection-diversity-v1.json.
Supersedes: no certificate; refines the count-only discovery coordinate.

The exact degree-two calculation treats a low-genus translation class as a
vertex in ``M/2M`` and uses the quotient metric

    mu_d(c) = min {(w,w) : w == c mod dM}.

For vertices of genera ``g`` and ``h``, the minimum intersection after
independent section translations is

    mu_d(c-c')/2 + g + h - 2.

Degree two is exhausted.  Degree three combines the existing complete
minimum-norm census with a deterministic Aut(M)-closed graph sample.  Degree
four is sampled, but the embedded two-torsion overlap and its norm scaling are
exact.  Equation-complexity and finite-place squareclass reductions are added
for the complete rational-bisection atlas.

With ``--comparison-only``, the same implementation instead exhausts degree
two on R17 and the requested rootless foundry frames, computes the full
automorphism action on ``M/2M``, and writes the exact comparative profile.

This is a computational lattice/equation profile.  Except for the already
certified rational bisections and their quadratic extensions, low-genus coset
vertices do not by themselves prove nefness, irreducibility, descent, or a
specialization rank jump.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter, deque
from fractions import Fraction
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SHORT_GRAM = ROOT / "elkies-k3/data/lattice/short_vector_basis_gram.txt"
PRIORITY = (
    ROOT
    / "artifacts/generated-results/elkies-2026-bisection-equation-priority-full.tsv"
)
BISECTIONS = (
    ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
)
COLLISIONS = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-2026-equation-bisection-collisions-full-compact.json"
)
DEGREE3 = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-degree3-complete-current-source-top5-v1.json"
)
FOUNDRY = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
SPECTRUM = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-multisection-spectrum-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-r17-multisection-diversity-v1.json"
)
DEFAULT_COMPARISON_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-r17-foundry-degree2-diversity-comparison-v1.json"
)
DEFAULT_COMPARISON_FRAMES = ("NS0032-F011", "NS0028-F005")
BISECTION_CONTROLS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "elkies_2026_bisection_specialization_controls_v1.json"
)
D2_COMPARISON = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-r17-foundry-degree2-diversity-comparison-v1.json"
)
DEFAULT_CONTROL_CALIBRATION_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-r17-bisection-control-diversity-calibration-v1.json"
)
DIMENSION = 17


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def gram_digest(gram: list[list[int]]) -> str:
    text = "\n".join(" ".join(map(str, row)) for row in gram) + "\n"
    return hashlib.sha256(text.encode()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load_gram(path: Path) -> list[list[int]]:
    return [
        [int(entry) for entry in line.split()]
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def exact_norm(vector: tuple[int, ...] | list[int], gram: list[list[int]]) -> int:
    return sum(
        vector[i] * gram[i][j] * vector[j]
        for i in range(DIMENSION)
        for j in range(DIMENSION)
    )


def residue_id(residue: tuple[int, ...], degree: int) -> int:
    result = 0
    power = 1
    for entry in residue:
        result += entry * power
        power *= degree
    return result


def decode_residue(identifier: int, degree: int) -> tuple[int, ...]:
    entries = []
    for unused in range(DIMENSION):
        entries.append(identifier % degree)
        identifier //= degree
    assert identifier == 0
    return tuple(entries)


def negate_residue(residue: tuple[int, ...], degree: int) -> tuple[int, ...]:
    return tuple((-entry) % degree for entry in residue)


class CosetOracle:
    """Floating CVP branch decisions with exact integral output checks."""

    def __init__(
        self,
        gram: list[list[int]],
        degree: int,
        *,
        float_type: str = "dd",
        precision: int = 160,
    ) -> None:
        from fpylll import FPLLL, GSO, IntegerMatrix

        if float_type == "mpfr":
            FPLLL.set_precision(precision)
        self.gram = gram
        self.degree = degree
        self.gso = GSO.Mat(
            IntegerMatrix.from_matrix(gram),
            gram=True,
            float_type=float_type,
            update=True,
        )
        self.mu = [
            [self.gso.get_mu(i, j) if i > j else 0.0 for j in range(DIMENSION)]
            for i in range(DIMENSION)
        ]
        # The digit vector itself is a representative.  This deliberately
        # generous bound prevents an incorrect empty CVP ball.
        self.distance_bound = (
            (degree - 1) ** 2
            * sum(abs(entry) for row in gram for entry in row)
            / (degree * degree)
            + 1.0
        )

    def solve(self, residue: tuple[int, ...]) -> tuple[int, tuple[int, ...], float]:
        from fpylll import Enumeration

        degree = self.degree
        target = [
            -(
                residue[i]
                + sum(
                    residue[j] * self.mu[j][i]
                    for j in range(i + 1, DIMENSION)
                )
            )
            / degree
            for i in range(DIMENSION)
        ]
        solutions = Enumeration(self.gso).enumerate(
            0, DIMENSION, self.distance_bound, 0, target=target
        )
        if not solutions:
            raise RuntimeError("CVP enumeration returned no solution")
        reported_distance, coordinates = solutions[0]
        closest = tuple(int(round(entry)) for entry in coordinates)
        if any(
            abs(entry - integer) > 1.0e-7
            for entry, integer in zip(coordinates, closest)
        ):
            raise RuntimeError("nonintegral CVP output coordinates")
        representative = tuple(
            residue[index] + degree * closest[index]
            for index in range(DIMENSION)
        )
        if any(
            (representative[index] - residue[index]) % degree
            for index in range(DIMENSION)
        ):
            raise RuntimeError("CVP output lies in the wrong coset")
        norm = exact_norm(representative, self.gram)
        error = abs(degree * degree * reported_distance - norm)
        if error > 1.0e-7 or norm < 0 or norm % 2:
            raise RuntimeError(f"invalid CVP output norm={norm}, error={error}")
        return norm, representative, error


def fwht(values: np.ndarray, *, inverse: bool = False) -> np.ndarray:
    result = values.astype(np.int64, copy=True)
    width = 1
    while width < len(result):
        blocks = result.reshape(-1, 2 * width)
        left = blocks[:, :width].copy()
        right = blocks[:, width:].copy()
        blocks[:, :width] = left + right
        blocks[:, width:] = left - right
        width *= 2
    if inverse:
        if np.any(result % len(result)):
            raise RuntimeError("nonintegral inverse Walsh transform")
        result //= len(result)
    return result


def xor_correlation(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return fwht(fwht(left) * fwht(right), inverse=True)


def binary_rank(masks: list[int] | np.ndarray) -> int:
    pivots: dict[int, int] = {}
    for raw in masks:
        value = int(raw)
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def prime_field_rank(rows: list[tuple[int, ...]], prime: int) -> int:
    values = [[entry % prime for entry in row] for row in rows]
    rank = 0
    for column in range(DIMENSION):
        pivot = next(
            (index for index in range(rank, len(values)) if values[index][column]),
            None,
        )
        if pivot is None:
            continue
        values[rank], values[pivot] = values[pivot], values[rank]
        inverse = pow(values[rank][column], -1, prime)
        values[rank] = [(entry * inverse) % prime for entry in values[rank]]
        for index in range(len(values)):
            if index == rank or not values[index][column]:
                continue
            factor = values[index][column]
            values[index] = [
                (left - factor * right) % prime
                for left, right in zip(values[index], values[rank])
            ]
        rank += 1
        if rank == DIMENSION:
            break
    return rank


def entropy_bits(counts: Counter | dict) -> float:
    total = sum(counts.values())
    if not total:
        return 0.0
    return -sum(
        (count / total) * math.log2(count / total)
        for count in counts.values()
        if count
    )


def quantiles(values: list[int] | np.ndarray) -> dict[str, float | int]:
    array = np.asarray(values)
    return {
        "minimum": int(array.min()),
        "q10": float(np.quantile(array, 0.10)),
        "median": float(np.quantile(array, 0.50)),
        "q90": float(np.quantile(array, 0.90)),
        "q99": float(np.quantile(array, 0.99)),
        "maximum": int(array.max()),
    }


def parse_priority() -> tuple[list[dict[str, str]], np.ndarray, np.ndarray]:
    with PRIORITY.open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    masks = np.asarray([int(row["orbit_mask"], 0) for row in rows], dtype=np.int32)
    representatives = np.asarray(
        [[int(entry) for entry in row["short_basis_w"].split()] for row in rows],
        dtype=np.int64,
    )
    if len(set(map(int, masks))) != len(rows) != 0:
        raise RuntimeError("priority table has duplicate masks")
    return rows, masks, representatives


def integral_automorphism_group(gram: list[list[int]]) -> list[list[list[int]]]:
    """Enumerate the full integral automorphism group using PARI qfauto."""

    from sage.all import ZZ, identity_matrix, matrix, pari

    form = matrix(ZZ, gram)
    data = pari(form).qfauto()
    claimed_order = int(data[0])
    generators = [matrix(ZZ, raw) for raw in data[1]]
    identity = identity_matrix(ZZ, DIMENSION)
    group = {tuple(map(int, identity.list())): identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = current * generator
            if candidate.transpose() * form * candidate != form:
                raise RuntimeError("PARI qfauto generator convention mismatch")
            key = tuple(map(int, candidate.list()))
            if key not in group:
                group[key] = candidate
                frontier.append(candidate)
    if len(group) != claimed_order:
        raise RuntimeError(
            f"qfauto enumeration found {len(group)} elements, expected {claimed_order}"
        )
    return [
        [[int(item[row, column]) for column in range(DIMENSION)] for row in range(DIMENSION)]
        for item in group.values()
    ]


def mod_two_action(mask: int, automorphism: list[list[int]]) -> int:
    """Apply a column-convention integral automorphism to an F2 mask."""

    result = 0
    for row in range(DIMENSION):
        parity = 0
        for column in range(DIMENSION):
            if (mask >> column) & 1:
                parity ^= automorphism[row][column] & 1
        result |= parity << row
    return result


def exact_mod_two_orbits(
    masks: np.ndarray, automorphisms: list[list[list[int]]]
) -> dict:
    """Return exact Aut(M)-orbit statistics for a selected subset of M/2M."""

    selected = set(map(int, masks))
    unseen = set(selected)
    orbit_sizes = Counter()
    while unseen:
        initial = min(unseen)
        orbit = {mod_two_action(initial, automorphism) for automorphism in automorphisms}
        if not orbit <= selected:
            raise RuntimeError("qualifying d=2 subset is not Aut(M)-stable")
        unseen.difference_update(orbit)
        orbit_sizes[len(orbit)] += 1
    orbit_count = sum(orbit_sizes.values())
    probabilities = [
        size / len(masks)
        for size, count in orbit_sizes.items()
        for unused in range(count)
    ]
    shannon = -sum(probability * math.log2(probability) for probability in probabilities)
    return {
        "orbits": orbit_count,
        "orbit_size_histogram": {
            str(size): count for size, count in sorted(orbit_sizes.items())
        },
        "vertex_weighted_orbit_entropy_bits": shannon,
        "effective_orbit_count": 2**shannon,
        "uniform_orbit_entropy_bits": math.log2(orbit_count),
        "induced_action_order": len(
            {
                tuple(mod_two_action(1 << column, automorphism) for column in range(DIMENSION))
                for automorphism in automorphisms
            }
        ),
    }


def complete_degree_two_structure(
    gram: list[list[int]],
    automorphisms: list[list[list[int]]],
    *,
    expected_rational_count: int | None = None,
    audit_stride: int = 4096,
) -> tuple[dict, np.ndarray]:
    """Compute a complete invariant d=2 profile for a rootless rank-17 frame."""

    size = 1 << DIMENSION
    oracle = CosetOracle(gram, 2)
    minima = np.empty(size, dtype=np.int16)
    maximum_error = 0.0
    for identifier in range(size):
        norm, unused, error = oracle.solve(decode_residue(identifier, 2))
        minima[identifier] = norm
        maximum_error = max(maximum_error, error)

    histogram = Counter(map(int, minima))
    rational_masks = np.flatnonzero(
        (minima >= 10) & ((minima % 4) == 2)
    ).astype(np.int32)
    genus_one_masks = np.flatnonzero(
        (minima >= 8) & ((minima % 4) == 0)
    ).astype(np.int32)
    if expected_rational_count is not None and len(rational_masks) != expected_rational_count:
        raise RuntimeError(
            f"rational d=2 count {len(rational_masks)} != expected {expected_rational_count}"
        )

    rational = np.zeros(size, dtype=np.int64)
    rational[rational_masks] = 1
    genus_one = np.zeros(size, dtype=np.int64)
    genus_one[genus_one_masks] = 1
    rr = xor_correlation(rational, rational)
    gg = xor_correlation(genus_one, genus_one)
    rg = xor_correlation(rational, genus_one)
    identifiers = np.arange(size)

    def same_set_histogram(correlation: np.ndarray) -> dict[str, int]:
        return {
            str(norm): int(
                correlation[(minima == norm) & (identifiers != 0)].sum() // 2
            )
            for norm in sorted(histogram)
            if norm
        }

    def cross_set_histogram(correlation: np.ndarray) -> dict[str, int]:
        return {
            str(norm): int(correlation[minima == norm].sum())
            for norm in sorted(histogram)
        }

    rr_histogram = same_set_histogram(rr)
    gg_histogram = same_set_histogram(gg)
    rg_histogram = cross_set_histogram(rg)
    if sum(rr_histogram.values()) != len(rational_masks) * (len(rational_masks) - 1) // 2:
        raise RuntimeError("rational pair accounting failed")
    if sum(gg_histogram.values()) != len(genus_one_masks) * (len(genus_one_masks) - 1) // 2:
        raise RuntimeError("genus-one pair accounting failed")
    if sum(rg_histogram.values()) != len(rational_masks) * len(genus_one_masks):
        raise RuntimeError("mixed pair accounting failed")

    norm_four = np.flatnonzero(minima == 4).astype(np.int32)
    present = rational.astype(bool)
    degrees = np.zeros(len(rational_masks), dtype=np.int32)
    for delta in norm_four:
        degrees += present[rational_masks ^ delta]

    unseen = set(map(int, rational_masks))
    components = []
    while unseen:
        initial = unseen.pop()
        queue = deque([initial])
        count = 0
        while queue:
            vertex = queue.popleft()
            count += 1
            for delta in norm_four:
                neighbour = vertex ^ int(delta)
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    queue.append(neighbour)
        components.append(count)

    delta_set = set(map(int, norm_four))
    compatible_delta_pairs = []
    for index, left in enumerate(norm_four):
        compatible_delta_pairs.extend(
            (int(left), int(right))
            for right in norm_four[index + 1 :]
            if int(left ^ right) in delta_set
        )
    triangle_anchor_count = sum(
        int(np.count_nonzero(present[rational_masks ^ left] & present[rational_masks ^ right]))
        for left, right in compatible_delta_pairs
    )
    if triangle_anchor_count % 3:
        raise RuntimeError("triangle anchor count is not divisible by three")
    triangles = triangle_anchor_count // 3
    edge_count = int(rr_histogram.get("4", 0))
    wedges = int(sum(int(degree) * (int(degree) - 1) // 2 for degree in degrees))
    pair_count = len(rational_masks) * (len(rational_masks) - 1) // 2

    maximum_minimum = int(minima.max())
    audit = CosetOracle(gram, 2, float_type="mpfr", precision=256)
    audit_ids = sorted(
        set(range(0, size, audit_stride))
        | set(map(int, np.flatnonzero(minima == maximum_minimum)))
    )
    maximum_audit_error = 0.0
    for identifier in audit_ids:
        norm, unused, error = audit.solve(decode_residue(identifier, 2))
        if norm != int(minima[identifier]):
            raise RuntimeError("degree-two cross-precision minimum mismatch")
        maximum_audit_error = max(maximum_audit_error, error)

    rational_orbits = exact_mod_two_orbits(rational_masks, automorphisms)
    genus_one_orbits = exact_mod_two_orbits(genus_one_masks, automorphisms)
    inherited = {
        str(4 * norm): count
        for norm, count in sorted(histogram.items())
        if norm >= 8
    }
    return (
        {
            "status": "PASS_COMPLETE_COMPUTATIONAL_COSET_AND_GRAPH_PROFILE",
            "translation_cosets": size,
            "minimum_norm_histogram": {
                str(norm): count for norm, count in sorted(histogram.items())
            },
            "discrete_covering_radius": {
                "maximum_minimum_norm": maximum_minimum,
                "squared_radius": f"{maximum_minimum}/4",
                "radius": f"sqrt({maximum_minimum})/2",
            },
            "vertices": {
                "rational": {
                    "genus": 0,
                    "threshold": 10,
                    "count": len(rational_masks),
                    "minimum_norm_histogram": {
                        str(norm): int(np.count_nonzero(minima[rational_masks] == norm))
                        for norm in sorted(set(map(int, minima[rational_masks])))
                    },
                    "span_dimension_over_F2": binary_rank(rational_masks),
                    "aut_M": rational_orbits,
                },
                "genus_one": {
                    "genus": 1,
                    "threshold": 8,
                    "count": len(genus_one_masks),
                    "minimum_norm_histogram": {
                        str(norm): int(np.count_nonzero(minima[genus_one_masks] == norm))
                        for norm in sorted(set(map(int, minima[genus_one_masks])))
                    },
                    "span_dimension_over_F2": binary_rank(genus_one_masks),
                    "aut_M": genus_one_orbits,
                    "geometric_boundary": (
                        "All-section nonnegativity is exact; nefness, irreducibility, "
                        "descent, and rank gain are not inferred."
                    ),
                },
            },
            "quotient_pairing": {
                "formula": "min intersection=mu_2(c-c')/2+g+h-2",
                "rational_rational_minimum_norm_distribution": rr_histogram,
                "rational_rational_minimum_intersection_distribution": {
                    str(int(norm) // 2 - 2): count
                    for norm, count in rr_histogram.items()
                },
                "genus_one_genus_one_minimum_norm_distribution": gg_histogram,
                "rational_genus_one_minimum_norm_distribution": rg_histogram,
                "separation_entropy_bits_rational_pairs": entropy_bits(rr_histogram),
                "total_rational_pairs": pair_count,
            },
            "rational_zero_intersection_graph": {
                "vertices": len(rational_masks),
                "edges": edge_count,
                "edge_density": edge_count / pair_count,
                "connected_components": len(components),
                "component_sizes": sorted(components, reverse=True),
                "degree_histogram": {
                    str(key): int(value)
                    for key, value in zip(*np.unique(degrees, return_counts=True))
                },
                "degree_quantiles": quantiles(degrees),
                "wedges": wedges,
                "triangles": triangles,
                "global_transitivity": 3 * triangles / wedges if wedges else 0.0,
                "compatible_norm_four_delta_pairs": len(compatible_delta_pairs),
                "clique_counts": {
                    "k=1": len(rational_masks),
                    "k=2": edge_count,
                    "k=3": triangles,
                    "k>=4": "not counted",
                },
                "clique_boundary": (
                    "Pairwise quotient minima permit independent translations; a clique "
                    "does not certify one simultaneous representative choice."
                ),
            },
            "exact_d2_into_d4": {
                "embedding": "c mod 2M maps to 2c mod 4M",
                "minimum_norm_scaling": "mu_4(2c)=4*mu_2(c)",
                "inherited_genus_one_vertices": sum(inherited.values()),
                "minimum_norm_histogram": inherited,
                "inherited_rational_bisection_mechanisms": len(rational_masks),
            },
            "numerical_certificate": {
                "every_returned_norm_recomputed_integrally": True,
                "maximum_dd_distance_to_integral_norm_error": maximum_error,
                "mpfr_precision_bits": 256,
                "mpfr_audited_residues": len(audit_ids),
                "mpfr_audit_includes_every_deepest_coset": True,
                "maximum_mpfr_distance_to_integral_norm_error": maximum_audit_error,
            },
        },
        minima,
    )
def degree_two_profile(
    gram: list[list[int]],
    priority_rows: list[dict[str, str]],
    rational_masks: np.ndarray,
    trace_vectors: np.ndarray,
    *,
    seed: int,
    angle_pairs: int,
) -> tuple[dict, np.ndarray]:
    size = 1 << DIMENSION
    oracle = CosetOracle(gram, 2)
    minima = np.empty(size, dtype=np.int16)
    max_error = 0.0
    for identifier in range(size):
        norm, unused, error = oracle.solve(decode_residue(identifier, 2))
        minima[identifier] = norm
        max_error = max(max_error, error)

    histogram = Counter(map(int, minima))
    expected = {0: 1, 4: 1311, 6: 26672, 8: 63925, 10: 39120, 12: 43}
    if dict(sorted(histogram.items())) != expected:
        raise RuntimeError(f"unexpected complete d=2 histogram {histogram}")
    if set(map(int, minima[rational_masks])) != {10}:
        raise RuntimeError("priority masks do not equal the norm-ten frontier")

    # Cross-precision audit of a deterministic residue subset, including all
    # 43 deepest holes.
    audit = CosetOracle(gram, 2, float_type="mpfr", precision=256)
    audit_ids = sorted(
        set(range(0, size, 4096)) | set(map(int, np.flatnonzero(minima == 12)))
    )
    audit_error = 0.0
    for identifier in audit_ids:
        norm, unused, error = audit.solve(decode_residue(identifier, 2))
        if norm != int(minima[identifier]):
            raise RuntimeError("degree-two cross-precision minimum mismatch")
        audit_error = max(audit_error, error)

    rational = np.zeros(size, dtype=np.int64)
    rational[rational_masks] = 1
    genus_one_masks = np.flatnonzero(
        (minima >= 8) & ((minima % 4) == 0)
    ).astype(np.int32)
    genus_one = np.zeros(size, dtype=np.int64)
    genus_one[genus_one_masks] = 1

    rr = xor_correlation(rational, rational)
    gg = xor_correlation(genus_one, genus_one)
    rg = xor_correlation(rational, genus_one)
    identifiers = np.arange(size)

    def same_set_histogram(correlation: np.ndarray) -> dict[str, int]:
        return {
            str(norm): int(
                correlation[(minima == norm) & (identifiers != 0)].sum() // 2
            )
            for norm in sorted(histogram)
            if norm
        }

    def cross_set_histogram(correlation: np.ndarray) -> dict[str, int]:
        return {
            str(norm): int(correlation[minima == norm].sum())
            for norm in sorted(histogram)
        }

    rr_histogram = same_set_histogram(rr)
    gg_histogram = same_set_histogram(gg)
    rg_histogram = cross_set_histogram(rg)
    if sum(rr_histogram.values()) != len(rational_masks) * (len(rational_masks) - 1) // 2:
        raise RuntimeError("rational pair accounting failed")
    if sum(gg_histogram.values()) != len(genus_one_masks) * (len(genus_one_masks) - 1) // 2:
        raise RuntimeError("genus-one pair accounting failed")
    if sum(rg_histogram.values()) != len(rational_masks) * len(genus_one_masks):
        raise RuntimeError("mixed pair accounting failed")

    norm_four = np.flatnonzero(minima == 4).astype(np.int32)
    present = rational.astype(bool)
    degrees = np.zeros(len(rational_masks), dtype=np.int32)
    for delta in norm_four:
        degrees += present[rational_masks ^ delta]

    unseen = set(map(int, rational_masks))
    components = []
    while unseen:
        initial = unseen.pop()
        queue = deque([initial])
        count = 0
        while queue:
            vertex = queue.popleft()
            count += 1
            for delta in norm_four:
                neighbour = vertex ^ int(delta)
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    queue.append(neighbour)
        components.append(count)

    delta_set = set(map(int, norm_four))
    compatible_delta_pairs = []
    for index, left in enumerate(norm_four):
        compatible_delta_pairs.extend(
            (int(left), int(right))
            for right in norm_four[index + 1 :]
            if int(left ^ right) in delta_set
        )
    triangle_anchor_count = sum(
        int(np.count_nonzero(present[rational_masks ^ left] & present[rational_masks ^ right]))
        for left, right in compatible_delta_pairs
    )
    if triangle_anchor_count % 3:
        raise RuntimeError("triangle anchor count is not divisible by three")
    triangles = triangle_anchor_count // 3

    gram_array = np.asarray(gram, dtype=np.int64)
    if np.any(np.einsum("ij,jk,ik->i", trace_vectors, gram_array, trace_vectors) != 10):
        raise RuntimeError("priority trace representatives do not all have norm ten")
    index_by_mask = np.full(size, -1, dtype=np.int32)
    index_by_mask[rational_masks] = np.arange(len(rational_masks), dtype=np.int32)
    weighted = trace_vectors @ gram_array
    edge_inner_products = Counter()
    for delta in norm_four:
        right_masks = rational_masks ^ delta
        valid = (rational_masks < right_masks) & present[right_masks]
        left_indices = np.flatnonzero(valid)
        right_indices = index_by_mask[right_masks[valid]]
        products = np.einsum(
            "ij,ij->i", weighted[left_indices], trace_vectors[right_indices]
        )
        values, counts = np.unique(products, return_counts=True)
        edge_inner_products.update(
            {int(value): int(count) for value, count in zip(values, counts)}
        )
    if sum(edge_inner_products.values()) != int(rr_histogram["4"]):
        raise RuntimeError("edge angle accounting failed")

    rng = np.random.default_rng(seed)
    sampled_inner_products = Counter()
    remaining = angle_pairs
    while remaining:
        batch = min(remaining, 250_000)
        left = rng.integers(0, len(rational_masks), size=batch)
        right = rng.integers(0, len(rational_masks), size=batch)
        equal = left == right
        right[equal] = (right[equal] + 1) % len(rational_masks)
        products = np.einsum("ij,ij->i", weighted[left], trace_vectors[right])
        values, counts = np.unique(products, return_counts=True)
        sampled_inner_products.update(
            {int(value): int(count) for value, count in zip(values, counts)}
        )
        remaining -= batch

    complexity = {}
    for key in (
        "group_addition_upper_bound",
        "support_count",
        "dependency_count",
        "coordinate_input_bits",
        "maximum_absolute_coefficient",
        "coefficient_l1",
        "minimal_unoriented_count",
    ):
        complexity[key] = quantiles([int(row[key]) for row in priority_rows])

    rational_pairs = len(rational_masks) * (len(rational_masks) - 1) // 2
    return (
        {
            "degree": 2,
            "status": "PASS_COMPLETE_COMPUTATIONAL_COSET_AND_GRAPH_PROFILE",
            "translation_cosets": size,
            "minimum_norm_histogram": {
                str(key): value for key, value in sorted(histogram.items())
            },
            "discrete_covering_radius": {
                "maximum_minimum_norm": int(minima.max()),
                "squared_radius": "12/4=3",
                "radius": "sqrt(3)",
            },
            "vertices": {
                "rational": {
                    "genus": 0,
                    "threshold": 10,
                    "count": len(rational_masks),
                    "frontier_count": int(np.count_nonzero(minima == 10)),
                    "deep_hole_count": 0,
                    "aut_M_orbits": len(rational_masks),
                    "orbit_entropy_bits": math.log2(len(rational_masks)),
                    "span_dimension_over_F2": binary_rank(rational_masks),
                },
                "genus_one": {
                    "genus": 1,
                    "threshold": 8,
                    "count": len(genus_one_masks),
                    "frontier_count": int(np.count_nonzero(minima == 8)),
                    "deep_hole_count": int(np.count_nonzero(minima == 12)),
                    "deep_hole_norms": {"12": int(np.count_nonzero(minima == 12))},
                    "aut_M_orbits": len(genus_one_masks),
                    "orbit_entropy_bits": math.log2(len(genus_one_masks)),
                    "span_dimension_over_F2": binary_rank(genus_one_masks),
                    "geometric_boundary": (
                        "All-section nonnegativity is exact; nefness, irreducibility, "
                        "descent, and rank gain are not inferred."
                    ),
                },
            },
            "quotient_pairing": {
                "formula": "min intersection=mu_2(c-c')/2+g+h-2",
                "rational_rational_minimum_norm_distribution": rr_histogram,
                "rational_rational_minimum_intersection_distribution": {
                    str(int(norm) // 2 - 2): count
                    for norm, count in rr_histogram.items()
                },
                "genus_one_genus_one_minimum_norm_distribution": gg_histogram,
                "rational_genus_one_minimum_norm_distribution": rg_histogram,
                "separation_entropy_bits_rational_pairs": entropy_bits(rr_histogram),
                "total_rational_pairs": rational_pairs,
            },
            "rational_zero_intersection_graph": {
                "vertices": len(rational_masks),
                "edges": int(rr_histogram["4"]),
                "connected_components": len(components),
                "component_sizes": sorted(components, reverse=True),
                "degree_histogram": {
                    str(key): int(value)
                    for key, value in zip(*np.unique(degrees, return_counts=True))
                },
                "degree_quantiles": quantiles(degrees),
                "zero_intersection_graph_cliques": {
                    "k=1": len(rational_masks),
                    "k=2": int(rr_histogram["4"]),
                    "k=3": triangles,
                    "k>=4": "not counted in this first profile",
                },
                "clique_boundary": (
                    "A clique says every pair has minimum intersection zero after "
                    "independent translations. It does not assert that one simultaneous "
                    "choice of representatives realizes all zero intersections."
                ),
                "compatible_norm_four_delta_pairs": len(compatible_delta_pairs),
            },
            "priority_trace_angle_gauge": {
                "warning": (
                    "Angles use the equation-priority norm-ten representative of each "
                    "coset. They are reproducible but not invariant under independent "
                    "section translations; the quotient separation spectrum above is "
                    "the invariant replacement."
                ),
                "cosine_formula": "cos(theta)=inner_product/10",
                "zero_intersection_edges_exact_inner_product_histogram": {
                    str(key): value for key, value in sorted(edge_inner_products.items())
                },
                "all_pairs_sample_with_replacement": angle_pairs,
                "all_pairs_sample_inner_product_histogram": {
                    str(key): value
                    for key, value in sorted(sampled_inner_products.items())
                },
            },
            "equation_complexity_weights": complexity,
            "numerical_certificate": {
                "every_returned_norm_recomputed_integrally": True,
                "maximum_dd_distance_to_integral_norm_error": max_error,
                "mpfr_precision_bits": 256,
                "mpfr_audited_residues": len(audit_ids),
                "mpfr_audit_includes_all_norm_12_deep_holes": True,
                "maximum_mpfr_distance_to_integral_norm_error": audit_error,
            },
        },
        minima,
    )


def aut_closed_sample(degree: int, sample_size: int, seed: int) -> list[tuple[int, ...]]:
    rng = random.Random(seed ^ (degree << 24))
    residues = {(0,) * DIMENSION}
    while len(residues) < sample_size:
        residue = tuple(rng.randrange(degree) for unused in range(DIMENSION))
        residues.add(residue)
        residues.add(negate_residue(residue, degree))
    return sorted(residues)


class DisjointSet:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.size = [1] * size

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        left = self.find(left)
        right = self.find(right)
        if left == right:
            return
        if self.size[left] < self.size[right]:
            left, right = right, left
        self.parent[right] = left
        self.size[left] += self.size[right]


def sampled_graph_profile(
    gram: list[list[int]],
    degree: int,
    genera: tuple[int, ...],
    sample_size: int,
    seed: int,
) -> dict:
    residues = aut_closed_sample(degree, sample_size, seed)
    oracle = CosetOracle(gram, degree)
    minima: dict[int, int] = {}
    representatives: dict[int, tuple[int, ...]] = {}
    max_error = 0.0
    vertices = []
    categories: dict[int, list[tuple[int, ...]]] = {genus: [] for genus in genera}
    for residue in residues:
        identifier = residue_id(residue, degree)
        norm, representative, error = oracle.solve(residue)
        minima[identifier] = norm
        representatives[identifier] = representative
        max_error = max(max_error, error)
        residue_norm = exact_norm(residue, gram)
        for genus in genera:
            threshold = 2 * degree * degree - 2 * genus + 2
            if (residue_norm + 2 * genus - 2) % (2 * degree):
                continue
            if norm >= threshold:
                categories[genus].append(residue)
                vertices.append((residue, genus))

    category_rows = []
    for genus in genera:
        selected = categories[genus]
        identifiers = [residue_id(row, degree) for row in selected]
        orbits = {
            min(identifier, residue_id(negate_residue(row, degree), degree))
            for identifier, row in zip(identifiers, selected)
        }
        norm_histogram = Counter(minima[identifier] for identifier in identifiers)
        category_rows.append(
            {
                "genus": genus,
                "threshold": 2 * degree * degree - 2 * genus + 2,
                "sampled_vertices": len(selected),
                "sampled_aut_M_orbits": len(orbits),
                "sampled_orbit_entropy_bits": entropy_bits(
                    Counter(
                        min(
                            identifier,
                            residue_id(negate_residue(row, degree), degree),
                        )
                        for identifier, row in zip(identifiers, selected)
                    )
                ),
                "sampled_minimum_norm_histogram": {
                    str(key): value for key, value in sorted(norm_histogram.items())
                },
                "sampled_span_dimension_mod_prime_divisor": prime_field_rank(
                    selected, 3 if degree == 3 else 2
                ),
            }
        )

    difference_cache: dict[int, tuple[int, tuple[int, ...]]] = {}
    intersection_histogram = Counter()
    label_pair_histogram: dict[str, Counter] = {}
    edges = []
    for left in range(len(vertices)):
        left_residue, left_genus = vertices[left]
        for right in range(left + 1, len(vertices)):
            right_residue, right_genus = vertices[right]
            difference = tuple(
                (left_residue[index] - right_residue[index]) % degree
                for index in range(DIMENSION)
            )
            negative = negate_residue(difference, degree)
            difference_id = min(
                residue_id(difference, degree), residue_id(negative, degree)
            )
            if difference_id not in difference_cache:
                canonical = decode_residue(difference_id, degree)
                norm, representative, unused = oracle.solve(canonical)
                difference_cache[difference_id] = (norm, representative)
            difference_norm = difference_cache[difference_id][0]
            intersection = difference_norm // 2 + left_genus + right_genus - 2
            intersection_histogram[intersection] += 1
            label = f"g{min(left_genus, right_genus)}-g{max(left_genus, right_genus)}"
            label_pair_histogram.setdefault(label, Counter())[intersection] += 1
            if intersection <= 1:
                edges.append((left, right))

    disjoint = DisjointSet(len(vertices))
    for left, right in edges:
        disjoint.union(left, right)
    component_histogram = Counter(
        disjoint.size[index]
        for index in range(len(vertices))
        if disjoint.find(index) == index
    )
    sample_minimum_histogram = Counter(minima.values())
    return {
        "sample_status": "PASS_DETERMINISTIC_AUT_CLOSED_COMPUTATIONAL_SAMPLE",
        "degree": degree,
        "sample_seed": seed,
        "sampled_translation_cosets": len(residues),
        "aut_M_closed": True,
        "sampled_minimum_norm_histogram": {
            str(key): value for key, value in sorted(sample_minimum_histogram.items())
        },
        "sampled_maximum_minimum_norm": max(minima.values()),
        "sampled_discrete_covering_radius_lower_bound": (
            f"sqrt({max(minima.values())})/{degree}"
        ),
        "categories": category_rows,
        "small_intersection_graph": {
            "vertices": len(vertices),
            "edge_rule": "minimum intersection at most one",
            "edges": len(edges),
            "component_size_histogram": {
                str(key): value for key, value in sorted(component_histogram.items())
            },
            "minimum_intersection_histogram": {
                str(key): value for key, value in sorted(intersection_histogram.items())
            },
            "label_pair_minimum_intersection_histograms": {
                label: {str(key): value for key, value in sorted(histogram.items())}
                for label, histogram in sorted(label_pair_histogram.items())
            },
        },
        "numerical_certificate": {
            "every_returned_norm_recomputed_integrally": True,
            "maximum_dd_distance_to_integral_norm_error": max_error,
            "difference_cosets_solved": len(difference_cache),
        },
        "boundary": (
            "The graph is induced by a deterministic Aut(M)-closed sample; its "
            "edge, component, orbit, span, and covering data are not a full census."
        ),
    }


def complete_degree_three_profile() -> dict:
    payload = json.loads(DEGREE3.read_text())
    row = next(
        item for item in payload["spectra"] if item["frame_id"] == "NS0001-F001"
    )
    histogram = {int(key): value for key, value in row[
        "minimum_norm_histogram_all_translation_cosets"
    ].items()}

    def category(genus: int) -> dict:
        threshold = 2 * 3 * 3 - 2 * genus + 2
        admissible = {
            norm: count
            for norm, count in histogram.items()
            if norm >= threshold and (norm + 2 * genus - 2) % 6 == 0
        }
        count = sum(admissible.values())
        if count % 2:
            raise RuntimeError("nonzero degree-three vertices must pair under inversion")
        return {
            "genus": genus,
            "threshold": threshold,
            "count": count,
            "minimum_norm_histogram": {
                str(key): value for key, value in sorted(admissible.items())
            },
            "frontier_count": admissible.get(threshold, 0),
            "deep_hole_count": count - admissible.get(threshold, 0),
            "aut_M_orbits": count // 2,
            "orbit_entropy_bits": math.log2(count // 2),
        }

    return {
        "degree": 3,
        "status": "PASS_COMPLETE_COMPUTATIONAL_COSET_SPECTRUM",
        "translation_cosets": row["translation_cosets"],
        "minimum_norm_histogram": {
            str(key): value for key, value in sorted(histogram.items())
        },
        "discrete_covering_radius": {
            "maximum_minimum_norm": row["maximum_coset_minimum_norm"],
            "squared_radius": "26/9",
            "radius": "sqrt(26)/3",
        },
        "vertices": {
            "rational": category(0),
            "genus_one": category(1),
        },
        "aut_M_action": (
            "Aut(M)={+I,-I}; every qualifying nonzero coset forms a two-element orbit."
        ),
        "boundary": row["geometric_boundary"],
        "numerical_certificate": row["numerical_certificate"],
    }


def mod_fraction(fraction: Fraction, prime: int) -> int | None:
    denominator = fraction.denominator % prime
    if denominator == 0:
        return None
    return (fraction.numerator % prime) * pow(denominator, -1, prime) % prime


def p_valuation(value: Fraction, prime: int) -> int:
    if not value:
        return 10**9
    numerator = abs(value.numerator)
    denominator = value.denominator
    result = 0
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def polynomial_squareclass(coefficients: list[str], prime: int) -> tuple | None:
    fractions = [Fraction(value) for value in coefficients]
    minimum_valuation = min(p_valuation(value, prime) for value in fractions)
    if minimum_valuation == 10**9:
        return None
    scale = Fraction(prime) ** (-minimum_valuation)
    values = [mod_fraction(value * scale, prime) for value in fractions]
    if any(value is None for value in values):
        return None
    while values and values[-1] == 0:
        values.pop()
    if not values:
        return None
    degree = len(values) - 1
    leading = values[-1]
    legendre = 1 if pow(leading, (prime - 1) // 2, prime) == 1 else -1
    if degree == 0:
        return (minimum_valuation % 2, 0, legendre)
    inverse = pow(leading, -1, prime)
    monic = tuple((value * inverse) % prime for value in values[:-1])
    if degree == 1:
        return (minimum_valuation % 2, 1, legendre, *monic)
    if degree != 2:
        raise RuntimeError(f"expected quadratic squareclass, got degree {degree}")
    discriminant = (values[1] * values[1] - 4 * values[0] * values[2]) % prime
    if discriminant == 0:
        return (minimum_valuation % 2, 0, legendre)
    return (minimum_valuation % 2, 2, legendre, *monic)


def local_squareclass_profile(primes: list[int]) -> dict:
    collision_payload = json.loads(COLLISIONS.read_text())
    if collision_payload["distinct_quadratic_extensions"] != 39120:
        raise RuntimeError("global squareclass certificate is not the complete atlas")
    bisections = json.loads(BISECTIONS.read_text())["bisections"]
    if len(bisections) != 39120:
        raise RuntimeError("equation bisection input is incomplete")
    counters = {prime: Counter() for prime in primes}
    bad = Counter()
    joint = Counter()
    joint_good = 0
    for row in bisections:
        coefficients = row["residual_chord"]["q_coefficients"]
        keys = []
        good = True
        for prime in primes:
            key = polynomial_squareclass(coefficients, prime)
            if key is None:
                bad[prime] += 1
                good = False
            else:
                counters[prime][key] += 1
            keys.append(key)
        if good:
            joint[tuple(keys)] += 1
            joint_good += 1
    result = {}
    for prime in primes:
        counts = counters[prime]
        good = sum(counts.values())
        result[str(prime)] = {
            "good_reductions": good,
            "bad_or_undefined_reductions": bad[prime],
            "distinct_Gauss_local_squareclass_signatures": len(counts),
            "entropy_bits": entropy_bits(counts),
            "effective_number_of_classes": 2 ** entropy_bits(counts),
            "largest_collision_bucket": max(counts.values()),
            "singleton_buckets": sum(count == 1 for count in counts.values()),
        }
    return {
        "global_Q_t_squareclasses": 39120,
        "global_collisions": 0,
        "finite_place_reductions": result,
        "joint_signature": {
            "primes": primes,
            "records_good_at_every_prime": joint_good,
            "distinct_joint_signatures": len(joint),
            "entropy_bits": entropy_bits(joint),
            "largest_collision_bucket": max(joint.values()),
            "singleton_buckets": sum(count == 1 for count in joint.values()),
        },
        "convention": (
            "Each key records the Gauss valuation parity and the squarefree residual "
            "polynomial class in F_p(t), including the square/nonsquare leading scalar. "
            "Thus it is invariant under rational-square rescaling of q and represents "
            "the corresponding odd-p Gauss-local squareclass signature."
        ),
    }


def cross_degree_overlap(d2_minima: np.ndarray) -> dict:
    counts = Counter(map(int, d2_minima))
    embedded = {
        4 * norm: count
        for norm, count in counts.items()
        if norm >= 8 and norm % 4 == 0
    }
    # The norm-ten rational classes also satisfy the genus-one degree-four
    # congruence after multiplication by two.
    embedded[40] = counts[10]
    return {
        "common_torsion_point_model": "view M/dM as the d-torsion subgroup (1/d)M/M",
        "literal_intersection_formula": "T_d intersect T_e = T_gcd(d,e)",
        "d2_d3_literal_low_genus_overlap": 0,
        "d3_d4_literal_low_genus_overlap": 0,
        "reason_for_coprime_degrees": (
            "The torsion subgroups meet only at zero, and zero is not a qualifying "
            "low-genus vertex in the displayed categories."
        ),
        "d2_into_d4": {
            "embedding": "c mod 2M maps to 2c mod 4M",
            "minimum_norm_scaling": "mu_4(2c)=4*mu_2(c)",
            "genus_one_quadrisection_vertices_from_d2": sum(embedded.values()),
            "minimum_norm_histogram": {
                str(key): value for key, value in sorted(embedded.items())
            },
            "rational_bisections_reappearing_at_d4": counts[10],
            "rational_bisection_d4_label": (
                "genus one, minimum norm 40 (deep, not frontier)"
            ),
            "d2_genus_one_frontier_reappearing_at_d4": counts[8],
            "d2_genus_one_frontier_d4_label": (
                "genus one, minimum norm 32 (frontier)"
            ),
            "d2_norm12_deep_holes_reappearing_at_d4": counts[12],
            "d2_norm12_deep_hole_d4_label": "genus one, minimum norm 48",
        },
        "warning": (
            "There is no natural reduction map between coprime degrees. Any stronger "
            "d=2 versus d=3 comparison needs a declared M/6M CRT compatibility metric, "
            "not a representative-dependent overlap count."
        ),
    }


def graph_profile(adjacency: list[int]) -> dict:
    """Exact component, degree, and all-clique counts for a small bit graph."""

    vertex_count = len(adjacency)
    degrees = [row.bit_count() for row in adjacency]
    unseen = (1 << vertex_count) - 1
    components = []
    while unseen:
        first = unseen & -unseen
        unseen ^= first
        frontier = first
        size = 0
        while frontier:
            bit = frontier & -frontier
            frontier ^= bit
            vertex = bit.bit_length() - 1
            size += 1
            neighbours = adjacency[vertex] & unseen
            unseen ^= neighbours
            frontier |= neighbours
        components.append(size)

    clique_counts = Counter()

    def extend(candidates: int, size: int) -> None:
        while candidates:
            bit = candidates & -candidates
            candidates ^= bit
            vertex = bit.bit_length() - 1
            clique_counts[size + 1] += 1
            extend(candidates & adjacency[vertex], size + 1)

    extend((1 << vertex_count) - 1, 0)
    edges = sum(degrees) // 2
    pairs = vertex_count * (vertex_count - 1) // 2
    return {
        "vertices": vertex_count,
        "edges": edges,
        "edge_density": edges / pairs if pairs else None,
        "connected_components": len(components),
        "component_sizes": sorted(components, reverse=True),
        "degree_histogram": {
            str(key): value for key, value in sorted(Counter(degrees).items())
        },
        "degree_quantiles": quantiles(degrees),
        "clique_counts_by_size": {
            str(key): value for key, value in sorted(clique_counts.items())
        },
        "maximum_clique_size": max(clique_counts, default=0),
    }


def exceptional_pair_type(left: int, right: int) -> str:
    if not left and not right:
        return "both_zero"
    if not left or not right:
        return "one_zero"
    if left == right:
        return "same_nonzero_direction"
    return "independent_directions"


def mutual_information_bits(joint: Counter) -> float:
    total = sum(joint.values())
    if not total:
        return 0.0
    left = Counter()
    right = Counter()
    for (left_key, right_key), count in joint.items():
        left[left_key] += count
        right[right_key] += count
    result = 0.0
    for (left_key, right_key), count in joint.items():
        result += (count / total) * math.log2(
            count * total / (left[left_key] * right[right_key])
        )
    return result


def average_ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start + 1
        while stop < len(order) and values[order[stop]] == values[order[start]]:
            stop += 1
        average = (start + 1 + stop) / 2
        for position in range(start, stop):
            result[order[position]] = average
        start = stop
    return result


def pearson_correlation(left: list[float], right: list[float]) -> float | None:
    if len(left) != len(right) or len(left) < 2:
        return None
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    numerator = sum(
        (x - left_mean) * (y - right_mean) for x, y in zip(left, right)
    )
    left_norm = sum((x - left_mean) ** 2 for x in left)
    right_norm = sum((y - right_mean) ** 2 for y in right)
    if not left_norm or not right_norm:
        return None
    return numerator / math.sqrt(left_norm * right_norm)


def spearman_correlation(left: list[float], right: list[float]) -> float | None:
    return pearson_correlation(average_ranks(left), average_ranks(right))


def run_control_diversity_calibration(output: Path) -> None:
    """Profile the exact R17 bisections splitting at five certified controls."""

    if not D2_COMPARISON.exists():
        raise RuntimeError(
            "run --comparison-only first to create the pinned global d=2 profile"
        )
    gram = load_gram(SHORT_GRAM)
    priority_rows, priority_masks, unused = parse_priority()
    row_by_label = {
        f"orbit-{int(row['orbit_mask'], 0):05x}": row for row in priority_rows
    }
    if len(row_by_label) != 39120:
        raise RuntimeError("priority atlas labels are not unique and complete")
    global_profile = json.loads(D2_COMPARISON.read_text())["profiles"]["NS0001-F001"]
    global_pairs = {
        int(norm): count
        for norm, count in global_profile["quotient_pairing"][
            "rational_rational_minimum_norm_distribution"
        ].items()
    }
    global_pair_total = sum(global_pairs.values())
    global_zero_probability = global_pairs[4] / global_pair_total
    global_near_probability = (global_pairs[4] + global_pairs[6]) / global_pair_total

    controls = json.loads(BISECTION_CONTROLS.read_text())
    oracle = CosetOracle(gram, 2)
    difference_norm_cache: dict[int, int] = {}
    profiles = {}
    ordered_labels = []
    for fibre in controls["fibres"]:
        label = fibre["label"]
        ordered_labels.append(label)
        hits = sorted(fibre["hits"], key=lambda row: row["label"])
        masks = [int(hit["lattice_orbit_mask"]) for hit in hits]
        if len(masks) != len(set(masks)) != 0:
            raise RuntimeError(f"duplicate or empty split mask set for {label}")
        for hit, mask in zip(hits, masks):
            if hit["label"] != f"orbit-{mask:05x}" or hit["label"] not in row_by_label:
                raise RuntimeError(f"split label/mask mismatch for {label}")
            norm, representative, unused_error = oracle.solve(decode_residue(mask, 2))
            if norm != 10 or exact_norm(representative, gram) != 10:
                raise RuntimeError(f"split mask is not a rational bisection for {label}")

        quotient_masks = []
        for hit in hits:
            coordinates = hit["finite_quotient_class_modulo_generic_17"][
                "coordinates_over_f2"
            ]
            quotient_masks.append(
                sum((int(entry) & 1) << index for index, entry in enumerate(coordinates))
            )

        pair_histogram = Counter()
        relation_histogram = Counter()
        joint_histogram = Counter()
        pair_rows = []
        zero_adjacency = [0] * len(hits)
        near_adjacency = [0] * len(hits)
        for left in range(len(hits)):
            for right in range(left + 1, len(hits)):
                difference = masks[left] ^ masks[right]
                if difference not in difference_norm_cache:
                    norm, representative, unused_error = oracle.solve(
                        decode_residue(difference, 2)
                    )
                    if exact_norm(representative, gram) != norm:
                        raise RuntimeError("pair-difference norm recomputation failed")
                    difference_norm_cache[difference] = norm
                norm = difference_norm_cache[difference]
                relation = exceptional_pair_type(
                    quotient_masks[left], quotient_masks[right]
                )
                pair_histogram[norm] += 1
                relation_histogram[relation] += 1
                joint_histogram[(norm, relation)] += 1
                pair_rows.append(
                    {
                        "left": hits[left]["label"],
                        "right": hits[right]["label"],
                        "quotient_minimum_norm": norm,
                        "minimum_intersection": norm // 2 - 2,
                        "exceptional_pair_type": relation,
                    }
                )
                if norm == 4:
                    zero_adjacency[left] |= 1 << right
                    zero_adjacency[right] |= 1 << left
                if norm <= 6:
                    near_adjacency[left] |= 1 << right
                    near_adjacency[right] |= 1 << left

        zero_graph = graph_profile(zero_adjacency)
        near_graph = graph_profile(near_adjacency)
        pair_count = len(hits) * (len(hits) - 1) // 2
        for graph, probability in (
            (zero_graph, global_zero_probability),
            (near_graph, global_near_probability),
        ):
            expected = pair_count * probability
            variance = pair_count * probability * (1 - probability)
            graph["global_random_subset_expected_edges"] = expected
            graph["binomial_reference_z_score"] = (
                (graph["edges"] - expected) / math.sqrt(variance)
                if variance
                else None
            )
            graph["reference_boundary"] = (
                "The binomial z-score is a descriptive random-pair baseline; graph "
                "edges are dependent, so it is not a calibrated p-value."
            )

        complexity_keys = (
            "group_addition_upper_bound",
            "coordinate_input_bits",
            "support_count",
            "dependency_count",
            "disjoint_degree_in_pool",
        )
        complexities = {
            key: quantiles([int(row_by_label[hit["label"]][key]) for hit in hits])
            for key in complexity_keys
        }
        rank_lower_bound = int(fibre["rank_result"]["existing_unconditional_rank_lower_bound"])
        quotient_gain = int(
            fibre["split_class_span"]["known_public_complement_dimension"]
        )
        exceptional_span = int(
            fibre["split_class_span"]["dimension_modulo_generic_17"]
        )
        profiles[label] = {
            "parameter": fibre["parameter"],
            "certified_rank_lower_bound": rank_lower_bound,
            "displayed_quotient_rank_beyond_R17": quotient_gain,
            "split_bisection_vertices": len(hits),
            "split_bisection_labels": [hit["label"] for hit in hits],
            "lattice_mask_span_dimension_F2": binary_rank(masks),
            "exceptional_quotient_span_dimension_F2": exceptional_span,
            "exceptional_visibility_fraction": exceptional_span / quotient_gain,
            "independent_direction_efficiency_per_split_cover": exceptional_span / len(hits),
            "zero_exceptional_class_covers": sum(mask == 0 for mask in quotient_masks),
            "distinct_exceptional_quotient_classes_including_zero": len(set(quotient_masks)),
            "pair_separation": {
                "minimum_norm_histogram": {
                    str(key): value for key, value in sorted(pair_histogram.items())
                },
                "entropy_bits": entropy_bits(pair_histogram),
                "exceptional_pair_type_histogram": dict(sorted(relation_histogram.items())),
                "minimum_norm_by_exceptional_pair_type": {
                    f"{norm}|{relation}": count
                    for (norm, relation), count in sorted(joint_histogram.items())
                },
                "mutual_information_bits": mutual_information_bits(joint_histogram),
                "pairs": pair_rows,
            },
            "zero_intersection_graph": zero_graph,
            "intersection_at_most_one_graph": near_graph,
            "equation_and_global_graph_weights": complexities,
        }

    overlap = {}
    for left_index, left in enumerate(ordered_labels):
        left_set = set(profiles[left]["split_bisection_labels"])
        for right in ordered_labels[left_index + 1 :]:
            common = sorted(left_set & set(profiles[right]["split_bisection_labels"]))
            overlap[f"{left}|{right}"] = {
                "common_split_bisections": len(common),
                "labels": common,
            }

    metric_names = (
        "split_bisection_vertices",
        "lattice_mask_span_dimension_F2",
        "exceptional_quotient_span_dimension_F2",
        "exceptional_visibility_fraction",
        "independent_direction_efficiency_per_split_cover",
    )

    def correlations(labels: list[str]) -> dict:
        gains = [profiles[label]["displayed_quotient_rank_beyond_R17"] for label in labels]
        result = {}
        for metric in metric_names:
            values = [float(profiles[label][metric]) for label in labels]
            result[metric] = {
                "pearson_with_displayed_quotient_rank": pearson_correlation(values, gains),
                "spearman_with_displayed_quotient_rank": spearman_correlation(values, gains),
            }
        return result

    high_rank_labels = [label for label in ordered_labels if label != "icarm_curve394_rank_at_least_21"]
    all_split_labels = [
        split_label
        for label in ordered_labels
        for split_label in profiles[label]["split_bisection_labels"]
    ]
    payload = {
        "schema": "elkies-k3.r17-bisection-control-diversity-calibration.v1",
        "status": "EXACT_INPUT_PROFILES_DESCRIPTIVE_FIVE_CONTROL_CALIBRATION",
        "scope": (
            "The split-cover sets, lattice quotient minima, induced graph/clique "
            "profiles, equation weights, and displayed exceptional quotient classes "
            "are exact for the five certified fibres. Correlations are descriptive."
        ),
        "inputs": {
            relative(SHORT_GRAM): digest(SHORT_GRAM),
            relative(PRIORITY): digest(PRIORITY),
            relative(BISECTION_CONTROLS): digest(BISECTION_CONTROLS),
            relative(D2_COMPARISON): digest(D2_COMPARISON),
        },
        "global_R17_pair_background": {
            "minimum_norm_distribution": {
                str(key): value for key, value in sorted(global_pairs.items())
            },
            "zero_intersection_probability": global_zero_probability,
            "intersection_at_most_one_probability": global_near_probability,
        },
        "profiles": profiles,
        "cross_control_overlap": overlap,
        "cross_control_novelty": {
            "split_bisection_occurrences": len(all_split_labels),
            "distinct_split_bisections": len(set(all_split_labels)),
            "repeated_split_bisections": len(all_split_labels) - len(set(all_split_labels)),
            "novelty_fraction": len(set(all_split_labels)) / len(all_split_labels),
        },
        "descriptive_correlations": {
            "all_five_controls": correlations(ordered_labels),
            "four_high_rank_controls_only": correlations(high_rank_labels),
            "warning": (
                "n=5 and n=4, with historically selected positive controls. These "
                "coefficients diagnose the current sample and have no predictive or "
                "significance interpretation."
            ),
        },
        "proof_boundary": {
            "exact": (
                "Every selected bisection, pair minimum, graph edge, component, clique, "
                "lattice span, and exceptional quotient class in the pinned controls."
            ),
            "not_claimed": (
                "The controls do not estimate a population relationship, omitted "
                "degree-three/four visibility is not zero, and no diversity coordinate "
                "is asserted to cause or predict a rank jump."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/analyze_r17_multisection_diversity.py "
            "--control-calibration-only"
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "D2CONTROLCAL|{}|status={}".format(
            "|".join(
                f"{label}={profiles[label]['split_bisection_vertices']}"
                for label in ordered_labels
            ),
            payload["status"],
        ),
        flush=True,
    )
    print(f"D2CONTROLCAL|output={output}|sha256={digest(output)}")


def foundry_frame_records(payload: dict) -> dict[str, dict]:
    result = {}
    for ns_class in payload["ns_classes"]:
        for frame in ns_class["frames"]:
            result[frame["frame_id"]] = {
                **frame,
                "ns_id": ns_class["ns_id"],
            }
    return result


def run_degree_two_comparison(
    output: Path, comparison_frames: list[str], audit_stride: int
) -> None:
    """Build the exact R17/foundry d=2 structural comparison artifact."""

    foundry = json.loads(FOUNDRY.read_text())
    spectrum = json.loads(SPECTRUM.read_text())
    records = foundry_frame_records(foundry)
    pilot_rows = {row["frame_id"]: row for row in spectrum["targets"]}
    frame_ids = ["NS0001-F001"] + [
        frame_id for frame_id in comparison_frames if frame_id != "NS0001-F001"
    ]
    if len(frame_ids) != len(set(frame_ids)):
        raise RuntimeError("comparison frame ids must be distinct")

    profiles = {}
    lattice_rows = {}
    for frame_id in frame_ids:
        if frame_id not in records or frame_id not in pilot_rows:
            raise RuntimeError(f"frame {frame_id} is absent from a pinned input")
        frame = records[frame_id]
        if frame["root_type"] != "0" or int(frame["mw_rank_for_rho_19"]) != DIMENSION:
            raise RuntimeError(f"frame {frame_id} is not a rootless rank-17 target")
        gram = [[int(entry) for entry in row] for row in frame["reduced_gram"]]
        automorphisms = integral_automorphism_group(gram)
        recorded_order = int(frame["rootless_intrinsics"]["automorphism_group_order"])
        if len(automorphisms) != recorded_order:
            raise RuntimeError(f"automorphism order mismatch for {frame_id}")
        expected = int(
            pilot_rows[frame_id]["degree_two"]["rational_bisections"][
                "translation_orbits_with_minimum_norm_ten"
            ]
        )
        profile, unused = complete_degree_two_structure(
            gram,
            automorphisms,
            expected_rational_count=expected,
            audit_stride=audit_stride,
        )
        profiles[frame_id] = profile
        lattice_rows[frame_id] = {
            "name": "published R17" if frame_id == "NS0001-F001" else frame_id,
            "ns_id": frame["ns_id"],
            "determinant": int(frame["determinant"]),
            "reduced_gram_sha256": gram_digest(gram),
            "automorphism_group_order": recorded_order,
            "induced_mod_two_action_order": profile["vertices"]["rational"][
                "aut_M"
            ]["induced_action_order"],
        }

    comparison_vectors = {}
    for frame_id, profile in profiles.items():
        rational = profile["vertices"]["rational"]
        genus_one = profile["vertices"]["genus_one"]
        graph = profile["rational_zero_intersection_graph"]
        comparison_vectors[frame_id] = {
            "rational_vertices": rational["count"],
            "rational_Aut_orbits": rational["aut_M"]["orbits"],
            "rational_vertex_weighted_orbit_entropy_bits": rational["aut_M"][
                "vertex_weighted_orbit_entropy_bits"
            ],
            "genus_one_vertices": genus_one["count"],
            "maximum_coset_minimum_norm": profile["discrete_covering_radius"][
                "maximum_minimum_norm"
            ],
            "rational_span_dimension_F2": rational["span_dimension_over_F2"],
            "rational_pair_separation_entropy_bits": profile["quotient_pairing"][
                "separation_entropy_bits_rational_pairs"
            ],
            "zero_intersection_edges": graph["edges"],
            "zero_intersection_edge_density": graph["edge_density"],
            "zero_intersection_degree_median": graph["degree_quantiles"]["median"],
            "zero_intersection_triangles": graph["triangles"],
            "zero_intersection_global_transitivity": graph["global_transitivity"],
            "inherited_genus_one_d4_vertices": profile["exact_d2_into_d4"][
                "inherited_genus_one_vertices"
            ],
        }

    r17 = comparison_vectors["NS0001-F001"]
    deltas_from_r17 = {}
    for frame_id, vector in comparison_vectors.items():
        if frame_id == "NS0001-F001":
            continue
        deltas_from_r17[frame_id] = {
            key: value - r17[key]
            for key, value in vector.items()
            if isinstance(value, (int, float)) and not isinstance(value, bool)
        }

    payload = {
        "schema": "elkies-k3.r17-foundry-degree2-diversity-comparison.v1",
        "status": "PASS_COMPLETE_R17_FOUNDRY_DEGREE2_DIVERSITY_COMPARISON",
        "scope": (
            "Every M/2M coset, rational/genus-one vertex, pair separation, rational "
            "zero-intersection edge, graph component, and triangle is counted exactly. "
            "The d=2 to d=4 inherited mass is exact; no nonembedded d=4 census is made."
        ),
        "inputs": {
            relative(FOUNDRY): digest(FOUNDRY),
            relative(SPECTRUM): digest(SPECTRUM),
        },
        "lattices": lattice_rows,
        "profiles": profiles,
        "comparison_vectors": comparison_vectors,
        "deltas_from_R17": deltas_from_r17,
        "interpretation_boundary": {
            "what_is_discriminated": (
                "The vector separates raw abundance from symmetry reduction, quotient "
                "separation, local graph density, triangle abundance, covering depth, "
                "and mechanisms inherited at degree four."
            ),
            "what_is_not_claimed": (
                "No scalar diversity score is fitted, and no profile coordinate is "
                "claimed to predict specialization rank jumps without calibration."
            ),
            "clique_warning": (
                "A quotient-graph clique permits a separate minimizing translation for "
                "each pair; it need not lift to one simultaneously pairwise-disjoint set."
            ),
            "geometry_warning": (
                "Genus-one vertices are lattice candidates. Nefness, irreducibility, "
                "field of definition, independence, and rank gain remain separate gates."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/analyze_r17_multisection_diversity.py "
            "--comparison-only"
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "D2DIVERSITY|{}|status={}".format(
            "|".join(
                f"{frame_id}={profiles[frame_id]['vertices']['rational']['count']}"
                for frame_id in frame_ids
            ),
            payload["status"],
        ),
        flush=True,
    )
    print(f"D2DIVERSITY|output={output}|sha256={digest(output)}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--comparison-output", type=Path, default=DEFAULT_COMPARISON_OUTPUT
    )
    parser.add_argument(
        "--control-calibration-output",
        type=Path,
        default=DEFAULT_CONTROL_CALIBRATION_OUTPUT,
    )
    parser.add_argument(
        "--comparison-frame-id", action="append", dest="comparison_frames"
    )
    parser.add_argument("--comparison-only", action="store_true")
    parser.add_argument("--control-calibration-only", action="store_true")
    parser.add_argument("--comparison-audit-stride", type=int, default=4096)
    parser.add_argument("--sample-size", type=int, default=1025)
    parser.add_argument("--seed", type=int, default=20260902)
    parser.add_argument("--angle-pairs", type=int, default=1_000_000)
    parser.add_argument(
        "--local-primes",
        default="97,101,103,107,109,113,127,131,137,139,149,151,157,163",
        help="comma-separated odd primes (default: every prime from 97 through 163)",
    )
    parser.add_argument("--skip-local-squareclasses", action="store_true")
    arguments = parser.parse_args()
    arguments.output = arguments.output.resolve()
    arguments.comparison_output = arguments.comparison_output.resolve()
    arguments.control_calibration_output = arguments.control_calibration_output.resolve()
    primes = [int(value) for value in arguments.local_primes.split(",") if value]
    if arguments.sample_size < 3 or arguments.angle_pairs < 1:
        parser.error("sample-size must be at least 3 and angle-pairs must be positive")
    if arguments.comparison_audit_stride < 1:
        parser.error("comparison-audit-stride must be positive")
    if arguments.comparison_only and arguments.control_calibration_only:
        parser.error("comparison-only and control-calibration-only are mutually exclusive")
    if arguments.comparison_only:
        run_degree_two_comparison(
            arguments.comparison_output,
            arguments.comparison_frames or list(DEFAULT_COMPARISON_FRAMES),
            arguments.comparison_audit_stride,
        )
        return
    if arguments.control_calibration_only:
        run_control_diversity_calibration(arguments.control_calibration_output)
        return

    gram = load_gram(SHORT_GRAM)
    if len(gram) != DIMENSION or any(len(row) != DIMENSION for row in gram):
        raise RuntimeError("R17 Gram matrix must be 17 by 17")
    priority_rows, rational_masks, trace_vectors = parse_priority()
    degree_two, d2_minima = degree_two_profile(
        gram,
        priority_rows,
        rational_masks,
        trace_vectors,
        seed=arguments.seed,
        angle_pairs=arguments.angle_pairs,
    )
    degree_three = complete_degree_three_profile()
    degree_three["graph_sample"] = sampled_graph_profile(
        gram, 3, (0, 1), arguments.sample_size, arguments.seed
    )
    degree_four = sampled_graph_profile(
        gram, 4, (0, 1, 2), arguments.sample_size, arguments.seed
    )
    degree_four["status"] = "PASS_DETERMINISTIC_AUT_CLOSED_COMPUTATIONAL_SAMPLE"
    degree_four["exact_covering_radius_lower_bound_from_d2_overlap"] = {
        "maximum_minimum_norm_at_least": 48,
        "squared_radius_at_least": "48/16=3",
        "radius_at_least": "sqrt(3)",
    }

    foundry = json.loads(FOUNDRY.read_text())
    r17 = next(
        row for row in foundry["rootless_targets"] if row["frame_id"] == "NS0001-F001"
    )
    if int(r17["invariants"]["automorphism_group_order"]) != 2:
        raise RuntimeError("the recorded R17 automorphism group is no longer order two")

    local_squareclasses = (
        {"status": "SKIPPED_BY_COMMAND_LINE"}
        if arguments.skip_local_squareclasses
        else local_squareclass_profile(primes)
    )
    inputs = [SHORT_GRAM, PRIORITY, COLLISIONS, DEGREE3, FOUNDRY]
    if not arguments.skip_local_squareclasses:
        inputs.append(BISECTIONS)
    payload = {
        "schema": "elkies-k3.r17-multisection-diversity.v1",
        "status": "PASS_R17_MULTISECTION_DIVERSITY_FIRST_PROFILE",
        "scope": (
            "Degree two is complete. Degree three has a complete one-vertex spectrum "
            "and a sampled graph. Degree four is sampled apart from its exact embedded "
            "two-torsion overlap. Only the rational d=2 equation atlas carries certified "
            "curve and squareclass data."
        ),
        "inputs": {relative(path): digest(path) for path in inputs},
        "lattice": {
            "name": "published R17 / NS0001-F001",
            "rank": DIMENSION,
            "determinant": 948,
            "automorphism_group_order": 2,
            "automorphism_group": "{+I,-I}",
            "consequence": (
                "Orbit entropy mostly measures the small automorphism group here; it "
                "does not compress the d=2 atlas."
            ),
        },
        "definition": {
            "vertex": "a qualifying low-genus coset c in M/dM",
            "weight": (
                "coset minimum mu_d(c), arithmetic genus, depth above threshold, and "
                "equation-complexity proxy when an equation lift exists"
            ),
            "pair_metric": "mu_d(c-c')",
            "minimum_intersection": "mu_d(c-c')/2+g+h-2",
            "small_intersection_edge": "minimum intersection at most one",
            "aut_action": "the natural action of Aut(M) on M/dM",
        },
        "degree_two": degree_two,
        "degree_three": degree_three,
        "degree_four": degree_four,
        "cross_degree_overlap": cross_degree_overlap(d2_minima),
        "squareclass_diversity": local_squareclasses,
        "multisection_diversity_profile": {
            "recommendation": (
                "Keep a vector profile rather than collapsing immediately to one score: "
                "Aut-orbit entropy, covering radius/depth, mod-l span, quotient-separation "
                "entropy, small-intersection clique counts, equation cost, local "
                "squareclass entropy, and cross-degree novelty."
            ),
            "r17_findings": [
                "Aut(M)={+I,-I}, so d=2 orbit counts equal raw counts.",
                "The rational d=2 zero-intersection graph is connected despite broad degree variation.",
                "There are 43 d=2 norm-12 deep holes invisible to threshold-only counts.",
                "The complete d=3 rational and genus-one sets also contain deep vertices.",
                "All rational bisection squareclasses are globally distinct, but finite-place reductions quantify local collisions.",
                "The natural d=2 to d=4 embedding changes the label to genus one and moves rational bisections to norm 40.",
            ],
        },
        "proof_boundary": {
            "exact": (
                "d=2 coset minima and pair distributions, d=2 graph components/degrees/"
                "triangles, d=3 complete minimum spectrum and Aut orbit counts, the "
                "d=2 to d=4 norm scaling, equation-complexity inputs, and finite-place "
                "squareclass reductions (subject to the pinned equation artifact)."
            ),
            "computational": (
                "CVP branch choices use fplll double-double arithmetic with exact norm "
                "recomputation; a deterministic d=2 subset including every deepest hole "
                "is independently MPFR-audited."
            ),
            "sampled": "d=3 graph structure, and d=4 structure outside embedded two-torsion.",
            "geometric": (
                "For g>=1 or d>=3, coset qualification proves integrality and "
                "all-section nonnegativity only. It does not prove a curve, cover, "
                "descent, independence, or rank jump."
            ),
        },
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/analyze_r17_multisection_diversity.py"
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        "R17MULTIDIVERSITY|d2_rational=39120|d2_edges={}|d2_triangles={}|"
        "d3_rational={}|d4_sample={}|status={}".format(
            degree_two["rational_zero_intersection_graph"]["edges"],
            degree_two["rational_zero_intersection_graph"][
                "zero_intersection_graph_cliques"
            ]["k=3"],
            degree_three["vertices"]["rational"]["count"],
            degree_four["sampled_translation_cosets"],
            payload["status"],
        ),
        flush=True,
    )
    print(f"R17MULTIDIVERSITY|output={arguments.output}|sha256={digest(arguments.output)}")


if __name__ == "__main__":
    main()
