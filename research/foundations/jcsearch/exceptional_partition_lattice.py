"""Weighted partition-lattice data for exceptional-component intersections.

An atomic component is indexed by nonnegative integers ``(a,b)`` with
``N=2*a+3*b``.  Distinct components differ by transfers

    3 double atoms <-> 2 triple atoms.

At the minimal common boundary, every transfer occupies its own sixfold
collision cluster.  The universal transverse cluster algebra is one dual
number, so ``t`` independent transfers have length ``2**t`` and binomial
tangent-cone Hilbert vector.

The same data also describes the first off-diagonal normalization
correspondence inside one component: it exchanges the allocations at two
sixfold roots and therefore has two independent transfer clusters.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import comb


@dataclass(frozen=True, order=True)
class AtomicComponent:
    """The component ``C_(a,b)`` associated with ``2**a 3**b``."""

    a: int
    b: int

    def __post_init__(self) -> None:
        if self.a < 0 or self.b < 0:
            raise ValueError("atomic multiplicities must be nonnegative")
        if self.degree < 3:
            raise ValueError("the exceptional theory starts in degree three")

    @property
    def degree(self) -> int:
        return 2 * self.a + 3 * self.b

    @property
    def dimension(self) -> int:
        return self.a + self.b - 1

    @property
    def partition(self) -> tuple[int, ...]:
        return (3,) * self.b + (2,) * self.a

    @property
    def label(self) -> str:
        return f"C_({self.a},{self.b})"


@dataclass(frozen=True)
class IntersectionRecord:
    """One generic minimal-boundary intersection calculation."""

    degree: int
    kind: str
    left: AtomicComponent
    right: AtomicComponent
    minimal_common_coarsening: tuple[int, ...]
    support_dimension: int
    transfer_cluster_count: int
    transverse_intersection_length: int
    transverse_embedding_dimension: int
    intersection_tangent_dimension: int
    tangent_cone_hilbert_vector: tuple[int, ...]
    factors_over_independent_clusters: bool

    def to_dict(self) -> dict:
        result = asdict(self)
        result["left"]["partition"] = list(self.left.partition)
        result["left"]["label"] = self.left.label
        result["right"]["partition"] = list(self.right.partition)
        result["right"]["label"] = self.right.label
        result["minimal_common_coarsening"] = list(
            self.minimal_common_coarsening
        )
        result["tangent_cone_hilbert_vector"] = list(
            self.tangent_cone_hilbert_vector
        )
        result["universal_collision_factors"] = [
            {"type": "Z_1", "algebra": "k[e]/(e^2)", "length": 2}
            for _ in range(self.transfer_cluster_count)
        ]
        return result


def atomic_components(degree: int) -> tuple[AtomicComponent, ...]:
    """List the solutions of ``degree=2*a+3*b`` in increasing ``b``."""
    degree = int(degree)
    if degree < 3:
        return ()
    result = []
    for b in range(degree // 3 + 1):
        remainder = degree - 3 * b
        if remainder % 2 == 0:
            result.append(AtomicComponent(remainder // 2, b))
    return tuple(result)


def _record(
    kind: str,
    left: AtomicComponent,
    right: AtomicComponent,
    coarsening: tuple[int, ...],
    transfers: int,
) -> IntersectionRecord:
    if transfers <= 0:
        raise ValueError("an intersection record needs a nontrivial transfer")
    coarsening = tuple(sorted(coarsening, reverse=True))
    support_dimension = len(coarsening) - 1
    hilbert = tuple(comb(transfers, degree) for degree in range(transfers + 1))
    return IntersectionRecord(
        degree=left.degree,
        kind=kind,
        left=left,
        right=right,
        minimal_common_coarsening=coarsening,
        support_dimension=support_dimension,
        transfer_cluster_count=transfers,
        transverse_intersection_length=sum(hilbert),
        transverse_embedding_dimension=transfers,
        intersection_tangent_dimension=support_dimension + transfers,
        tangent_cone_hilbert_vector=hilbert,
        factors_over_independent_clusters=True,
    )


def component_pair_intersection(
    left: AtomicComponent, right: AtomicComponent
) -> IntersectionRecord:
    """Return the generic intersection on the minimal common boundary.

    Distinct solutions of ``2*a+3*b=N`` satisfy
    ``a-c=3*k`` and ``b-d=-2*k``.  Cancelling common atoms leaves ``abs(k)``
    primitive identities ``3*2=2*3=6``.  Keeping those primitive identities
    in separate clusters gives the unique minimal integer-partition
    coarsening.
    """
    if left.degree != right.degree:
        raise ValueError("components must have the same degree")
    if left == right:
        raise ValueError("use first_self_intersection for one component")
    if (left.a - right.a) % 3 or (left.b - right.b) % 2:
        raise AssertionError("equal-degree atomic types have integral transfers")
    transfers = abs(left.a - right.a) // 3
    assert transfers == abs(left.b - right.b) // 2
    coarsening = (
        (6,) * transfers
        + (3,) * min(left.b, right.b)
        + (2,) * min(left.a, right.a)
    )
    return _record("component_pair", left, right, coarsening, transfers)


def first_self_intersection(
    component: AtomicComponent,
) -> IntersectionRecord | None:
    """Return the first off-diagonal normalization pair, when it exists.

    It uses three double and two triple atoms twice, with the two allocations
    exchanged between distinct sixfold roots.  This is a branch-pair
    statement, not the scheme-theoretic intersection of a component with
    itself.
    """
    if component.a < 3 or component.b < 2:
        return None
    coarsening = (
        (6, 6)
        + (3,) * (component.b - 2)
        + (2,) * (component.a - 3)
    )
    return _record(
        "normalization_self_pair", component, component, coarsening, 2
    )


def intersection_census(
    minimum_degree: int = 3, maximum_degree: int = 20
) -> tuple[IntersectionRecord, ...]:
    """Enumerate distinct component pairs and first self-pairs by degree."""
    if minimum_degree < 3 or maximum_degree < minimum_degree:
        raise ValueError("use 3 <= minimum_degree <= maximum_degree")
    records = []
    for degree in range(minimum_degree, maximum_degree + 1):
        components = atomic_components(degree)
        for left_index, left in enumerate(components):
            self_record = first_self_intersection(left)
            if self_record is not None:
                records.append(self_record)
            for right in components[left_index + 1 :]:
                records.append(component_pair_intersection(left, right))
    return tuple(
        sorted(
            records,
            key=lambda record: (
                record.degree,
                record.kind != "component_pair",
                record.left.b,
                record.right.b,
            ),
        )
    )

