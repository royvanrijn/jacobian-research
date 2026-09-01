"""Basis-independent restrictions of finite codes to candidate sublattices.

The raw maps in :mod:`latent_lattice.finite` and :mod:`latent_lattice.local`
depend on choices of displayed Mordell--Weil generators and quotient bases.
This module forgets those choices.  It retains only invariants of the
restriction to a primitive candidate sublattice and of the induced finite
relation complex.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from hashlib import sha256
from itertools import permutations
import json
from math import gcd, log
from typing import Iterable, Sequence

from .finite import FiniteQuotientBlock
from .integer import rational_nullspace, rational_rank
from .local import ComponentBlock
from .relations import RelationComplex


def _modular_rank(rows: Sequence[Sequence[int]], prime: int) -> int:
    """Return exact row rank over ``F_prime`` for a prime modulus."""

    matrix = [
        [int(value) % prime for value in row]
        for row in rows
        if any(int(value) % prime for value in row)
    ]
    if not matrix:
        return 0
    width = len(matrix[0])
    if any(len(row) != width for row in matrix):
        raise ValueError("matrix row widths differ")
    rank = 0
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [(value * inverse) % prime for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank:
                continue
            coefficient = matrix[index][column]
            if coefficient:
                matrix[index] = [
                    (value - coefficient * pivot_value) % prime
                    for value, pivot_value in zip(matrix[index], matrix[rank])
                ]
        rank += 1
        if rank == len(matrix):
            break
    return rank


def _histogram(counter: Counter[tuple[int, ...]]) -> tuple[tuple[tuple[int, ...], int], ...]:
    return tuple(sorted((key, value) for key, value in counter.items()))


def _unoriented_class(values: Sequence[int], modulus: int) -> tuple[int, ...]:
    positive = tuple(int(value) % modulus for value in values)
    negative = tuple((-value) % modulus for value in positive)
    return min(positive, negative)


def _element_order(value: int, modulus: int) -> int:
    return modulus // gcd(modulus, int(value) % modulus)


def _candidate_mask(
    vertices: Sequence[Sequence[int]], basis_rows: Sequence[Sequence[int]]
) -> tuple[bool, ...]:
    if not basis_rows:
        raise ValueError("candidate basis must be nonempty")
    width = len(basis_rows[0])
    if any(len(row) != width for row in basis_rows):
        raise ValueError("candidate basis row widths differ")
    if any(len(vertex) != width for vertex in vertices):
        raise ValueError("candidate and relation-complex ambient widths differ")
    if rational_rank(basis_rows) != len(basis_rows):
        raise ValueError("candidate rows are not independent over Q")
    if len(basis_rows) == width:
        return tuple(True for _ in vertices)
    kernel = rational_nullspace(basis_rows)
    return tuple(
        all(sum(int(a) * int(b) for a, b in zip(vertex, normal)) == 0 for normal in kernel)
        for vertex in vertices
    )


def _restricted_prime_matrix(
    basis_rows: Sequence[Sequence[int]], block: FiniteQuotientBlock
) -> tuple[tuple[int, ...], ...]:
    if block.rows and len(block.rows[0]) != len(basis_rows[0]):
        raise ValueError("finite block and candidate ambient widths differ")
    ell = block.relation_prime
    return tuple(
        tuple(
            sum(int(code) * int(basis) for code, basis in zip(code_row, basis_row))
            % ell
            for basis_row in basis_rows
        )
        for code_row in block.rows
    )


@dataclass(frozen=True)
class PrimeCodeSignature:
    """Restriction invariant for one ``F_ell`` quotient block."""

    relation_prime: int
    ambient_quotient_dimension: int
    candidate_image_dimension: int
    candidate_kernel_dimension: int
    retained_ray_count: int
    zero_ray_count: int
    nonzero_unoriented_class_multiplicities: tuple[int, ...]
    ternary_type_histogram: tuple[tuple[tuple[int, ...], int], ...]
    scaled_type_histogram: tuple[tuple[tuple[int, ...], int], ...]

    @property
    def block_type(self) -> tuple[str, int, int]:
        return ("prime", self.relation_prime, self.ambient_quotient_dimension)

    def to_record(self) -> dict[str, object]:
        return {
            "relation_prime": self.relation_prime,
            "ambient_quotient_dimension": self.ambient_quotient_dimension,
            "candidate_image_dimension": self.candidate_image_dimension,
            "candidate_kernel_dimension": self.candidate_kernel_dimension,
            "retained_ray_count": self.retained_ray_count,
            "zero_ray_count": self.zero_ray_count,
            "nonzero_unoriented_class_multiplicities": list(
                self.nonzero_unoriented_class_multiplicities
            ),
            "ternary_type_histogram": [
                {"type": list(key), "count": count}
                for key, count in self.ternary_type_histogram
            ],
            "scaled_type_histogram": [
                {"type": list(key), "count": count}
                for key, count in self.scaled_type_histogram
            ],
        }


@dataclass(frozen=True)
class ComponentCodeSignature:
    """Restriction invariant for one cyclic component-group code."""

    modulus: int
    candidate_image_order: int
    candidate_kernel_index: int
    retained_ray_count: int
    element_order_histogram: tuple[tuple[int, int], ...]
    ternary_type_histogram: tuple[tuple[tuple[int, ...], int], ...]
    scaled_type_histogram: tuple[tuple[tuple[int, ...], int], ...]

    @property
    def block_type(self) -> tuple[str, int]:
        return ("component", self.modulus)

    def to_record(self) -> dict[str, object]:
        return {
            "modulus": self.modulus,
            "candidate_image_order": self.candidate_image_order,
            "candidate_kernel_index": self.candidate_kernel_index,
            "retained_ray_count": self.retained_ray_count,
            "element_order_histogram": [
                {"order": order, "count": count}
                for order, count in self.element_order_histogram
            ],
            "ternary_type_histogram": [
                {"type": list(key), "count": count}
                for key, count in self.ternary_type_histogram
            ],
            "scaled_type_histogram": [
                {"type": list(key), "count": count}
                for key, count in self.scaled_type_histogram
            ],
        }


@dataclass(frozen=True)
class CandidateFiniteSignature:
    """Canonical finite-aware invariant of one candidate relation subcomplex."""

    candidate_dimension: int
    retained_ray_count: int
    retained_ternary_relation_count: int
    retained_scaled_relation_count: int
    prime_blocks: tuple[PrimeCodeSignature, ...]
    component_blocks: tuple[ComponentCodeSignature, ...]
    joint_vertex_type_histogram: tuple[tuple[tuple[int, ...], int], ...]
    joint_ternary_type_histogram: tuple[tuple[tuple[int, ...], int], ...]
    canonical_digest: str

    def to_record(self) -> dict[str, object]:
        return {
            "candidate_dimension": self.candidate_dimension,
            "retained_ray_count": self.retained_ray_count,
            "retained_ternary_relation_count": self.retained_ternary_relation_count,
            "retained_scaled_relation_count": self.retained_scaled_relation_count,
            "prime_blocks": [block.to_record() for block in self.prime_blocks],
            "component_blocks": [block.to_record() for block in self.component_blocks],
            "joint_vertex_type_histogram": [
                {"type": list(key), "count": count}
                for key, count in self.joint_vertex_type_histogram
            ],
            "joint_ternary_type_histogram": [
                {"type": list(key), "count": count}
                for key, count in self.joint_ternary_type_histogram
            ],
            "canonical_digest": self.canonical_digest,
        }


def candidate_finite_signature_from_record(
    record: dict[str, object],
) -> CandidateFiniteSignature:
    """Rehydrate a signature record without any source-specific code maps."""

    def histogram(items) -> tuple[tuple[tuple[int, ...], int], ...]:
        return tuple(
            (tuple(map(int, item["type"])), int(item["count"])) for item in items
        )

    prime_blocks = tuple(
        PrimeCodeSignature(
            relation_prime=int(block["relation_prime"]),
            ambient_quotient_dimension=int(block["ambient_quotient_dimension"]),
            candidate_image_dimension=int(block["candidate_image_dimension"]),
            candidate_kernel_dimension=int(block["candidate_kernel_dimension"]),
            retained_ray_count=int(block["retained_ray_count"]),
            zero_ray_count=int(block["zero_ray_count"]),
            nonzero_unoriented_class_multiplicities=tuple(
                map(int, block["nonzero_unoriented_class_multiplicities"])
            ),
            ternary_type_histogram=histogram(block["ternary_type_histogram"]),
            scaled_type_histogram=histogram(block["scaled_type_histogram"]),
        )
        for block in record["prime_blocks"]
    )
    component_blocks = tuple(
        ComponentCodeSignature(
            modulus=int(block["modulus"]),
            candidate_image_order=int(block["candidate_image_order"]),
            candidate_kernel_index=int(block["candidate_kernel_index"]),
            retained_ray_count=int(block["retained_ray_count"]),
            element_order_histogram=tuple(
                (int(item["order"]), int(item["count"]))
                for item in block["element_order_histogram"]
            ),
            ternary_type_histogram=histogram(block["ternary_type_histogram"]),
            scaled_type_histogram=histogram(block["scaled_type_histogram"]),
        )
        for block in record["component_blocks"]
    )
    return CandidateFiniteSignature(
        candidate_dimension=int(record["candidate_dimension"]),
        retained_ray_count=int(record["retained_ray_count"]),
        retained_ternary_relation_count=int(record["retained_ternary_relation_count"]),
        retained_scaled_relation_count=int(record["retained_scaled_relation_count"]),
        prime_blocks=prime_blocks,
        component_blocks=component_blocks,
        joint_vertex_type_histogram=histogram(record["joint_vertex_type_histogram"]),
        joint_ternary_type_histogram=histogram(record["joint_ternary_type_histogram"]),
        canonical_digest=str(record["canonical_digest"]),
    )


def _normalized_counter(
    histogram: Sequence[tuple[tuple[int, ...], int]],
) -> dict[tuple[int, ...], float]:
    total = sum(count for _key, count in histogram)
    if not total:
        return {}
    return {key: count / total for key, count in histogram}


def _total_variation(
    left: Sequence[tuple[tuple[int, ...], int]],
    right: Sequence[tuple[tuple[int, ...], int]],
) -> float:
    a = _normalized_counter(left)
    b = _normalized_counter(right)
    return 0.5 * sum(abs(a.get(key, 0.0) - b.get(key, 0.0)) for key in a.keys() | b.keys())


def _multiplicity_profile(values: Sequence[int], total: int) -> tuple[float, ...]:
    if total <= 0:
        return ()
    return tuple(sorted((value / total for value in values), reverse=True))


def _profile_l1(left: Sequence[float], right: Sequence[float]) -> float:
    width = max(len(left), len(right))
    return sum(
        abs((left[index] if index < len(left) else 0.0) - (right[index] if index < len(right) else 0.0))
        for index in range(width)
    )


def _prime_block_distance(left: PrimeCodeSignature, right: PrimeCodeSignature) -> float:
    if left.block_type != right.block_type:
        return float("inf")
    left_total = max(1, left.retained_ray_count)
    right_total = max(1, right.retained_ray_count)
    return (
        abs(left.candidate_image_dimension - right.candidate_image_dimension)
        + abs(left.zero_ray_count / left_total - right.zero_ray_count / right_total)
        + 0.5
        * _profile_l1(
            _multiplicity_profile(left.nonzero_unoriented_class_multiplicities, left_total),
            _multiplicity_profile(right.nonzero_unoriented_class_multiplicities, right_total),
        )
        + _total_variation(left.ternary_type_histogram, right.ternary_type_histogram)
        + 0.5 * _total_variation(left.scaled_type_histogram, right.scaled_type_histogram)
    )


def _component_block_distance(
    left: ComponentCodeSignature, right: ComponentCodeSignature
) -> float:
    if left.block_type != right.block_type:
        return float("inf")
    left_orders = tuple(((order,), count) for order, count in left.element_order_histogram)
    right_orders = tuple(((order,), count) for order, count in right.element_order_histogram)
    return (
        abs(left.candidate_image_order / left.modulus - right.candidate_image_order / right.modulus)
        + _total_variation(left_orders, right_orders)
        + _total_variation(left.ternary_type_histogram, right.ternary_type_histogram)
        + 0.5 * _total_variation(left.scaled_type_histogram, right.scaled_type_histogram)
    )


def _minimum_type_matching(
    left: Sequence[object],
    right: Sequence[object],
    distance,
    *,
    allow_unmatched: bool = False,
) -> float:
    if len(left) != len(right) and not allow_unmatched:
        return float("inf")
    if not left or not right:
        return float("inf")
    smaller, larger = (left, right) if len(left) <= len(right) else (right, left)
    if len(larger) > 8:
        raise ValueError("exact block matching is bounded to at most eight equal-type blocks")
    return min(
        sum(distance(a, b) for a, b in zip(smaller, ordering))
        for ordering in permutations(larger, len(smaller))
    ) / len(smaller)


def finite_signature_distance(
    left: CandidateFiniteSignature,
    right: CandidateFiniteSignature,
    *,
    include_components: bool = True,
    active_prime_blocks_only: bool = False,
    allow_unmatched_blocks: bool = False,
) -> float:
    """Compare normalized source-free code profiles with exact block matching.

    Equal-type blocks are matched by a bounded exhaustive assignment, so their
    construction order and reduction primes play no role.  The distance is a
    calibration/search score, not a metric theorem and not an exact isometry
    certificate.  Conditioning on active blocks and allowing unmatched blocks
    are explicit statistical options for separating fibre-specific saturation
    from active-code shape; they must not be used to erase the strict profile.
    """

    if left.candidate_dimension != right.candidate_dimension:
        return float("inf")
    left_prime: dict[tuple[object, ...], list[PrimeCodeSignature]] = defaultdict(list)
    right_prime: dict[tuple[object, ...], list[PrimeCodeSignature]] = defaultdict(list)
    for block in left.prime_blocks:
        if not active_prime_blocks_only or block.candidate_image_dimension:
            left_prime[block.block_type].append(block)
    for block in right.prime_blocks:
        if not active_prime_blocks_only or block.candidate_image_dimension:
            right_prime[block.block_type].append(block)
    if not allow_unmatched_blocks and set(left_prime) != set(right_prime):
        return float("inf")
    common_prime_types = set(left_prime) & set(right_prime)
    pieces = [
        _minimum_type_matching(
            left_prime[key],
            right_prime[key],
            _prime_block_distance,
            allow_unmatched=allow_unmatched_blocks,
        )
        for key in sorted(common_prime_types, key=repr)
    ]
    if include_components:
        left_component: dict[tuple[object, ...], list[ComponentCodeSignature]] = defaultdict(list)
        right_component: dict[tuple[object, ...], list[ComponentCodeSignature]] = defaultdict(list)
        for block in left.component_blocks:
            left_component[block.block_type].append(block)
        for block in right.component_blocks:
            right_component[block.block_type].append(block)
        if not allow_unmatched_blocks and set(left_component) != set(right_component):
            return float("inf")
        common_component_types = set(left_component) & set(right_component)
        pieces.extend(
            _minimum_type_matching(
                left_component[key],
                right_component[key],
                _component_block_distance,
                allow_unmatched=allow_unmatched_blocks,
            )
            for key in sorted(common_component_types, key=repr)
        )
    if not pieces:
        return float("inf")
    return sum(pieces) / len(pieces)


def _prime_signature(
    basis_rows: Sequence[Sequence[int]],
    complex_: RelationComplex,
    retained: Sequence[bool],
    block: FiniteQuotientBlock,
    classes: tuple[tuple[int, ...], ...] | None = None,
) -> tuple[PrimeCodeSignature, tuple[tuple[int, ...], ...]]:
    ell = block.relation_prime
    if classes is None:
        classes = tuple(block.vector_class(vertex) for vertex in complex_.vertices)
    restricted = _restricted_prime_matrix(basis_rows, block)
    image_dimension = _modular_rank(restricted, ell)
    orbit_counts = Counter(
        _unoriented_class(value, ell)
        for value, keep in zip(classes, retained)
        if keep and any(value)
    )
    zero_count = sum(keep and not any(value) for value, keep in zip(classes, retained))
    ternary = Counter()
    for edge in complex_.ternary_relations:
        if not all(retained[index] for index in edge):
            continue
        values = tuple(classes[index] for index in edge)
        zeroes = sum(not any(value) for value in values)
        rank = _modular_rank(tuple(zip(*values)), ell) if values and values[0] else 0
        ternary[(zeroes, rank)] += 1
    scaled = Counter()
    for left, right, target, multiplier in complex_.scaled_relations:
        if not (retained[left] and retained[right] and retained[target]):
            continue
        values = (classes[left], classes[right], classes[target])
        source_zeroes = int(not any(values[0])) + int(not any(values[1]))
        target_zero = int(not any(values[2]))
        rank = _modular_rank(tuple(zip(*values)), ell) if values[0] else 0
        scaled[(multiplier % ell, source_zeroes, target_zero, rank)] += 1
    signature = PrimeCodeSignature(
        relation_prime=ell,
        ambient_quotient_dimension=block.quotient_dimension,
        candidate_image_dimension=image_dimension,
        candidate_kernel_dimension=len(basis_rows) - image_dimension,
        retained_ray_count=sum(retained),
        zero_ray_count=zero_count,
        nonzero_unoriented_class_multiplicities=tuple(sorted(orbit_counts.values())),
        ternary_type_histogram=_histogram(ternary),
        scaled_type_histogram=_histogram(scaled),
    )
    return signature, classes


def _component_signature(
    basis_rows: Sequence[Sequence[int]],
    complex_: RelationComplex,
    retained: Sequence[bool],
    block: ComponentBlock,
    classes: tuple[int, ...] | None = None,
) -> tuple[ComponentCodeSignature, tuple[int, ...]]:
    modulus = block.modulus
    if modulus < 1:
        raise ValueError("component modulus must be positive")
    if len(block.classes) != len(basis_rows[0]):
        raise ValueError("component block and candidate ambient widths differ")
    if classes is None:
        classes = tuple(block.vector_class(vertex) for vertex in complex_.vertices)
    restricted_values = tuple(block.vector_class(row) for row in basis_rows)
    divisor = modulus
    for value in restricted_values:
        divisor = gcd(divisor, value)
    image_order = modulus // divisor
    orders = Counter(
        _element_order(value, modulus)
        for value, keep in zip(classes, retained)
        if keep
    )
    ternary = Counter()
    for edge in complex_.ternary_relations:
        if not all(retained[index] for index in edge):
            continue
        values = tuple(classes[index] for index in edge)
        generated_order = modulus // gcd(modulus, *values)
        ternary[tuple(sorted(_element_order(value, modulus) for value in values)) + (generated_order,)] += 1
    scaled = Counter()
    for left, right, target, multiplier in complex_.scaled_relations:
        if not (retained[left] and retained[right] and retained[target]):
            continue
        source_orders = sorted(
            (_element_order(classes[left], modulus), _element_order(classes[right], modulus))
        )
        target_order = _element_order(classes[target], modulus)
        generated_order = modulus // gcd(
            modulus, classes[left], classes[right], classes[target]
        )
        scaled[(gcd(multiplier, modulus), *source_orders, target_order, generated_order)] += 1
    signature = ComponentCodeSignature(
        modulus=modulus,
        candidate_image_order=image_order,
        candidate_kernel_index=image_order,
        retained_ray_count=sum(retained),
        element_order_histogram=tuple(sorted(orders.items())),
        ternary_type_histogram=_histogram(ternary),
        scaled_type_histogram=_histogram(scaled),
    )
    return signature, classes


def _flatten_grouped_features(
    features: Iterable[tuple[tuple[object, ...], tuple[int, ...]]]
) -> tuple[int, ...]:
    groups: dict[tuple[object, ...], list[tuple[int, ...]]] = defaultdict(list)
    for block_type, feature in features:
        groups[block_type].append(feature)
    answer: list[int] = []
    for block_type, values in sorted(groups.items(), key=lambda item: repr(item[0])):
        encoded_type = json.dumps(block_type, separators=(",", ":"))
        answer.extend((len(encoded_type), *encoded_type.encode(), len(values)))
        for value in sorted(values):
            answer.extend((len(value), *value))
    return tuple(answer)


def candidate_finite_signature(
    basis_rows: Sequence[Sequence[int]],
    complex_: RelationComplex,
    *,
    finite_blocks: Sequence[FiniteQuotientBlock] = (),
    component_blocks: Sequence[ComponentBlock] = (),
    _finite_classes: Sequence[tuple[tuple[int, ...], ...]] | None = None,
    _component_classes: Sequence[tuple[int, ...]] | None = None,
) -> CandidateFiniteSignature:
    """Restrict exact finite codes to a candidate rational subspace.

    The digest is invariant under integral rebasing of the candidate, changes
    of basis in each good-reduction quotient, sign choices for the unoriented
    rays, and automorphisms of each cyclic component group.  Blocks of the
    same source-free type are treated as an unordered multiset.
    """

    retained = _candidate_mask(complex_.vertices, basis_rows)
    if _finite_classes is not None and len(_finite_classes) != len(finite_blocks):
        raise ValueError("precomputed finite-class block count differs")
    if _component_classes is not None and len(_component_classes) != len(component_blocks):
        raise ValueError("precomputed component-class block count differs")
    prime_results = [
        _prime_signature(
            basis_rows,
            complex_,
            retained,
            block,
            None if _finite_classes is None else _finite_classes[index],
        )
        for index, block in enumerate(finite_blocks)
    ]
    component_results = [
        _component_signature(
            basis_rows,
            complex_,
            retained,
            block,
            None if _component_classes is None else _component_classes[index],
        )
        for index, block in enumerate(component_blocks)
    ]
    prime_signatures = tuple(sorted((item[0] for item in prime_results), key=repr))
    component_signatures = tuple(
        sorted((item[0] for item in component_results), key=repr)
    )

    joint_vertices = Counter()
    for vertex, keep in enumerate(retained):
        if not keep:
            continue
        features = []
        for (signature, classes), _block in zip(prime_results, finite_blocks):
            features.append((signature.block_type, (int(not any(classes[vertex])),)))
        for (signature, classes), _block in zip(component_results, component_blocks):
            features.append(
                (signature.block_type, (_element_order(classes[vertex], signature.modulus),))
            )
        joint_vertices[_flatten_grouped_features(features)] += 1

    joint_edges = Counter()
    for edge in complex_.ternary_relations:
        if not all(retained[index] for index in edge):
            continue
        features = []
        for (signature, classes), _block in zip(prime_results, finite_blocks):
            values = tuple(classes[index] for index in edge)
            rank = _modular_rank(tuple(zip(*values)), signature.relation_prime) if values[0] else 0
            features.append(
                (signature.block_type, (sum(not any(value) for value in values), rank))
            )
        for (signature, classes), _block in zip(component_results, component_blocks):
            values = tuple(classes[index] for index in edge)
            features.append(
                (
                    signature.block_type,
                    tuple(sorted(_element_order(value, signature.modulus) for value in values))
                    + (signature.modulus // gcd(signature.modulus, *values),),
                )
            )
        joint_edges[_flatten_grouped_features(features)] += 1

    retained_ternary = sum(
        all(retained[index] for index in edge) for edge in complex_.ternary_relations
    )
    retained_scaled = sum(
        retained[left] and retained[right] and retained[target]
        for left, right, target, _multiplier in complex_.scaled_relations
    )
    provisional = {
        "candidate_dimension": len(basis_rows),
        "retained_ray_count": sum(retained),
        "retained_ternary_relation_count": retained_ternary,
        "retained_scaled_relation_count": retained_scaled,
        "prime_blocks": [block.to_record() for block in prime_signatures],
        "component_blocks": [block.to_record() for block in component_signatures],
        "joint_vertex_type_histogram": [
            [list(key), count] for key, count in sorted(joint_vertices.items())
        ],
        "joint_ternary_type_histogram": [
            [list(key), count] for key, count in sorted(joint_edges.items())
        ],
    }
    digest = sha256(
        json.dumps(provisional, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return CandidateFiniteSignature(
        candidate_dimension=len(basis_rows),
        retained_ray_count=sum(retained),
        retained_ternary_relation_count=retained_ternary,
        retained_scaled_relation_count=retained_scaled,
        prime_blocks=prime_signatures,
        component_blocks=component_signatures,
        joint_vertex_type_histogram=_histogram(joint_vertices),
        joint_ternary_type_histogram=_histogram(joint_edges),
        canonical_digest=digest,
    )


def candidate_finite_signatures(
    basis_matrices: Sequence[Sequence[Sequence[int]]],
    complex_: RelationComplex,
    *,
    finite_blocks: Sequence[FiniteQuotientBlock] = (),
    component_blocks: Sequence[ComponentBlock] = (),
) -> tuple[CandidateFiniteSignature, ...]:
    """Batch candidate restrictions with exact finite classes cached once.

    This is mathematically identical to repeated
    :func:`candidate_finite_signature` calls.  It avoids recomputing every
    ambient ray class for every candidate, which is important in a frozen
    dimension scan with hundreds of proper subspaces.
    """

    finite_classes = tuple(
        tuple(block.vector_class(vertex) for vertex in complex_.vertices)
        for block in finite_blocks
    )
    component_classes = tuple(
        tuple(block.vector_class(vertex) for vertex in complex_.vertices)
        for block in component_blocks
    )
    return tuple(
        candidate_finite_signature(
            basis_rows,
            complex_,
            finite_blocks=finite_blocks,
            component_blocks=component_blocks,
            _finite_classes=finite_classes,
            _component_classes=component_classes,
        )
        for basis_rows in basis_matrices
    )


def finite_joint_class_key(
    coordinates: Sequence[int],
    *,
    finite_blocks: Sequence[FiniteQuotientBlock] = (),
    component_blocks: Sequence[ComponentBlock] = (),
) -> tuple[int, ...]:
    """Return a source-local unoriented key for rare-code seed selection.

    Unlike :func:`candidate_finite_signature`, this key deliberately retains
    the chosen quotient coordinates.  It is suitable for ranking rare seeds
    inside one fibre, but must never be compared literally across fibres.
    """

    positive: list[int] = []
    negative: list[int] = []
    for block in finite_blocks:
        values = block.vector_class(coordinates)
        positive.extend((block.relation_prime, len(values), *values))
        negative.extend(
            (block.relation_prime, len(values), *((-value) % block.relation_prime for value in values))
        )
    for block in component_blocks:
        value = block.vector_class(coordinates)
        positive.extend((block.modulus, value))
        negative.extend((block.modulus, (-value) % block.modulus))
    return min(tuple(positive), tuple(negative))


def finite_rarity_weights(
    vectors: Sequence[Sequence[int]],
    *,
    finite_blocks: Sequence[FiniteQuotientBlock] = (),
    component_blocks: Sequence[ComponentBlock] = (),
) -> tuple[int, ...]:
    """Return exact inverse-frequency denominators for source-local seed keys.

    The returned integer for a vector is the population count of its joint
    finite-code orbit; smaller values are rarer and should be branched first.
    """

    keys = tuple(
        finite_joint_class_key(
            vector, finite_blocks=finite_blocks, component_blocks=component_blocks
        )
        for vector in vectors
    )
    counts = Counter(keys)
    return tuple(counts[key] for key in keys)


def finite_rarity_scores(
    vectors: Sequence[Sequence[int]],
    *,
    finite_blocks: Sequence[FiniteQuotientBlock] = (),
    component_blocks: Sequence[ComponentBlock] = (),
) -> tuple[float, ...]:
    """Return additive source-local rarity scores from separate code blocks.

    Scoring blocks separately avoids the degenerate product-code regime in
    which almost every joint key is unique.  These scores only order proposals
    inside one fibre; they are not literal cross-fibre fingerprints.
    """

    if not vectors:
        return ()
    populations: list[tuple[tuple[int, ...], ...]] = []
    for block in finite_blocks:
        populations.append(
            tuple(
                _unoriented_class(block.vector_class(vector), block.relation_prime)
                for vector in vectors
            )
        )
    for block in component_blocks:
        populations.append(
            tuple(
                (_element_order(block.vector_class(vector), block.modulus),)
                for vector in vectors
            )
        )
    if not populations:
        return tuple(0.0 for _ in vectors)
    counters = [Counter(classes) for classes in populations]
    total = len(vectors)
    return tuple(
        sum(
            -log(counter[classes[index]] / total)
            for counter, classes in zip(counters, populations)
        )
        for index in range(total)
    )
