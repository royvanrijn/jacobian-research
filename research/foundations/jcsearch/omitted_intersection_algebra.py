"""Combinatorics of primitive omitted-component merger faces.

At an ordered collection of distinct sixfold collision roots, a normalization
sheet chooses one of the two allocations ``(3, 0)`` and ``(0, 2)`` at every
root.  This module records the partition hypergraph of a finite collection of
such sheets and its branch-merger graph.

The algebra theorem is proved in the accompanying mathematical note.  The
combinatorics deliberately counts *active Hensel clusters*, not the rank of
an incidence matrix: distinct root clusters remain algebraically independent
even when the corresponding merger paths contain cycles.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from math import comb


BinaryBranch = tuple[int, ...]


@dataclass(frozen=True, order=True)
class MergerEdge:
    """An elementary compensated swap between two fixed-weight branches."""

    left: int
    right: int
    removed_root: int
    added_root: int


@dataclass(frozen=True)
class PrimitiveMergerHypergraph:
    """The binary cut hypergraph of separated primitive collision blocks."""

    branches: tuple[BinaryBranch, ...]
    active_roots: tuple[int, ...]
    root_cuts: tuple[tuple[int, ...], ...]
    merger_edges: tuple[MergerEdge, ...]
    component_count: int

    @property
    def transfer_cluster_count(self) -> int:
        return len(self.active_roots)

    @property
    def expected_tangent_excess(self) -> int:
        return self.transfer_cluster_count

    @property
    def predicted_transverse_length(self) -> int:
        return 2**self.transfer_cluster_count

    @property
    def predicted_hilbert_vector(self) -> tuple[int, ...]:
        t = self.transfer_cluster_count
        return tuple(comb(t, degree) for degree in range(t + 1))

    @property
    def cycle_rank(self) -> int:
        """Dimension of the cycle space of the branch-merger graph."""

        return len(self.merger_edges) - len(self.branches) + self.component_count


def fixed_weight_branches(root_count: int, chosen_count: int) -> tuple[BinaryBranch, ...]:
    """Return all binary allocation words of a fixed weight."""

    if root_count < 1 or not 0 <= chosen_count <= root_count:
        raise ValueError("use root_count >= 1 and 0 <= chosen_count <= root_count")
    result = []
    for chosen in combinations(range(root_count), chosen_count):
        chosen_set = set(chosen)
        result.append(tuple(int(root in chosen_set) for root in range(root_count)))
    return tuple(result)


def _components(vertex_count: int, edges: tuple[MergerEdge, ...]) -> int:
    adjacency = [set() for _ in range(vertex_count)]
    for edge in edges:
        adjacency[edge.left].add(edge.right)
        adjacency[edge.right].add(edge.left)
    unseen = set(range(vertex_count))
    count = 0
    while unseen:
        count += 1
        stack = [unseen.pop()]
        while stack:
            vertex = stack.pop()
            neighbours = adjacency[vertex] & unseen
            unseen.difference_update(neighbours)
            stack.extend(neighbours)
    return count


def primitive_merger_hypergraph(
    branches: tuple[BinaryBranch, ...] | list[BinaryBranch],
) -> PrimitiveMergerHypergraph:
    """Build root cuts and the Johnson-type merger graph.

    Each branch must have the same Hamming weight.  A merger edge joins two
    branches that differ by exchanging precisely one chosen and one unchosen
    root.  The root cuts are the active hyperedges: for each active root they
    list the sheets using allocation ``(3, 0)`` there.
    """

    branches = tuple(tuple(int(value) for value in branch) for branch in branches)
    if not branches:
        raise ValueError("at least one branch is required")
    root_count = len(branches[0])
    if root_count == 0 or any(len(branch) != root_count for branch in branches):
        raise ValueError("branches must be nonempty binary words of equal length")
    if any(value not in (0, 1) for branch in branches for value in branch):
        raise ValueError("primitive branches are binary")
    if len(set(branches)) != len(branches):
        raise ValueError("branches must be distinct")
    weights = {sum(branch) for branch in branches}
    if len(weights) != 1:
        raise ValueError("the branch-merger graph is defined within one component")

    root_cuts = tuple(
        tuple(index for index, branch in enumerate(branches) if branch[root])
        for root in range(root_count)
    )
    active_roots = tuple(
        root for root, cut in enumerate(root_cuts) if 0 < len(cut) < len(branches)
    )

    edges = []
    for left, right in combinations(range(len(branches)), 2):
        removed = [
            root
            for root in range(root_count)
            if branches[left][root] == 1 and branches[right][root] == 0
        ]
        added = [
            root
            for root in range(root_count)
            if branches[left][root] == 0 and branches[right][root] == 1
        ]
        if len(removed) == len(added) == 1:
            edges.append(MergerEdge(left, right, removed[0], added[0]))
    edge_tuple = tuple(edges)
    return PrimitiveMergerHypergraph(
        branches=branches,
        active_roots=active_roots,
        root_cuts=root_cuts,
        merger_edges=edge_tuple,
        component_count=_components(len(branches), edge_tuple),
    )
