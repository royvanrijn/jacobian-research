"""Exact finite-permutation-group utilities for normal-cover certificates.

The implementation is deliberately dependency-free.  It is intended for
small explicit certificates and for independent replay of GAP-produced data,
not as a replacement for GAP's subgroup libraries in large searches.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
from itertools import combinations
from typing import Iterable, Mapping, Sequence


Permutation = tuple[int, ...]
Subgroup = frozenset[Permutation]


def identity(degree: int) -> Permutation:
    return tuple(range(degree))


def validate_permutation(
    values: Sequence[int],
    degree: int,
    name: str = "permutation",
) -> Permutation:
    result = tuple(values)
    if len(result) != degree:
        raise ValueError(f"{name}: expected degree {degree}, got {len(result)}")
    if set(result) != set(range(degree)):
        raise ValueError(f"{name}: values are not a permutation")
    return result


def multiply(left: Permutation, right: Permutation) -> Permutation:
    """Return ``left * right``, applying ``right`` first."""

    if len(left) != len(right):
        raise ValueError("cannot multiply permutations of different degrees")
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(value: Permutation) -> Permutation:
    result = [0] * len(value)
    for source, target in enumerate(value):
        result[target] = source
    return tuple(result)


def generated_group(
    generators: Iterable[Permutation],
    *,
    max_order: int = 10_000,
) -> tuple[Permutation, ...]:
    generators = tuple(generators)
    if not generators:
        raise ValueError("at least one generator is required")
    degree = len(generators[0])
    if any(len(generator) != degree for generator in generators):
        raise ValueError("all generators must have the same degree")
    one = identity(degree)
    steps = tuple(
        dict.fromkeys(generators + tuple(inverse(value) for value in generators))
    )
    seen = {one}
    queue: deque[Permutation] = deque([one])
    while queue:
        current = queue.popleft()
        for step in steps:
            candidate = multiply(current, step)
            if candidate in seen:
                continue
            seen.add(candidate)
            if len(seen) > max_order:
                raise ValueError(f"generated group exceeds max_order={max_order}")
            queue.append(candidate)
    return tuple(sorted(seen))


def subgroup_generated_by(
    group: Sequence[Permutation],
    generators: Iterable[Permutation],
) -> Subgroup:
    group_set = set(group)
    generators = tuple(generators)
    if not generators:
        return frozenset({identity(len(group[0]))})
    subgroup = frozenset(generated_group(generators, max_order=len(group)))
    if not subgroup <= group_set:
        raise ValueError("subgroup generators do not generate a subgroup of the group")
    return subgroup


def conjugate_subgroup(subgroup: Subgroup, element: Permutation) -> Subgroup:
    element_inverse = inverse(element)
    return frozenset(
        multiply(multiply(element_inverse, member), element)
        for member in subgroup
    )


def subgroup_core(group: Sequence[Permutation], subgroup: Subgroup) -> Subgroup:
    conjugates = [conjugate_subgroup(subgroup, element) for element in group]
    return frozenset.intersection(*conjugates)


def stabilizer(group: Sequence[Permutation], point: int) -> Subgroup:
    return frozenset(element for element in group if element[point] == point)


def orbit(group: Sequence[Permutation], point: int) -> tuple[int, ...]:
    return tuple(sorted({element[point] for element in group}))


def conjugate_union(group: Sequence[Permutation], subgroup: Subgroup) -> Subgroup:
    return frozenset(
        member
        for element in group
        for member in conjugate_subgroup(subgroup, element)
    )


def enumerate_subgroups(
    group: Sequence[Permutation],
    *,
    max_subgroups: int = 10_000,
) -> tuple[Subgroup, ...]:
    """Enumerate all subgroups by adjoining one group element at a time."""

    one = identity(len(group[0]))
    trivial = frozenset({one})
    known = {trivial}
    queue: deque[Subgroup] = deque([trivial])
    while queue:
        subgroup = queue.popleft()
        for element in group:
            if element in subgroup:
                continue
            candidate = subgroup_generated_by(group, tuple(subgroup) + (element,))
            if candidate in known:
                continue
            known.add(candidate)
            if len(known) > max_subgroups:
                raise ValueError(
                    f"subgroup enumeration exceeds max_subgroups={max_subgroups}"
                )
            queue.append(candidate)
    return tuple(sorted(known, key=lambda item: (len(item), tuple(sorted(item)))))


def subgroup_conjugacy_classes(
    group: Sequence[Permutation],
    subgroups: Sequence[Subgroup],
) -> tuple[tuple[Subgroup, ...], ...]:
    remaining = set(subgroups)
    classes: list[tuple[Subgroup, ...]] = []
    while remaining:
        representative = min(
            remaining,
            key=lambda item: (len(item), tuple(sorted(item))),
        )
        conjugates = {
            conjugate_subgroup(representative, element)
            for element in group
        }
        current = tuple(
            sorted(
                conjugates & remaining,
                key=lambda item: (len(item), tuple(sorted(item))),
            )
        )
        classes.append(current)
        remaining -= conjugates
    return tuple(classes)


def normal_covering_number(
    group: Sequence[Permutation],
    *,
    max_subgroups: int = 10_000,
) -> int | None:
    """Return the exact normal covering number, or ``None`` for cyclic groups."""

    full = frozenset(group)
    subgroups = tuple(
        subgroup
        for subgroup in enumerate_subgroups(group, max_subgroups=max_subgroups)
        if subgroup != full
    )
    classes = subgroup_conjugacy_classes(group, subgroups)
    covers = tuple(
        conjugate_union(group, conjugacy_class[0])
        for conjugacy_class in classes
    )
    if frozenset().union(*covers) != full:
        return None
    for size in range(1, len(covers) + 1):
        for selected in combinations(covers, size):
            if frozenset().union(*selected) == full:
                return size
    raise AssertionError("finite noncyclic group has no computed normal covering")


@dataclass(frozen=True)
class ComponentCover:
    name: str
    basepoint: int
    orbit: tuple[int, ...]
    stabilizer: Subgroup
    core: Subgroup

    @property
    def degree(self) -> int:
        return len(self.orbit)

    @property
    def index(self) -> int:
        return self.degree


@dataclass(frozen=True)
class NormalCoverAnalysis:
    degree: int
    group: tuple[Permutation, ...]
    components: tuple[ComponentCover, ...]
    covered_elements: Subgroup
    common_core: Subgroup
    normal_covering_number: int | None

    @property
    def is_normal_cover(self) -> bool:
        return self.covered_elements == frozenset(self.group)

    @property
    def is_faithful(self) -> bool:
        return len(self.common_core) == 1

    @property
    def factorization_shape(self) -> tuple[int, ...]:
        return tuple(sorted(component.degree for component in self.components))

    @property
    def index_sum(self) -> int:
        return sum(component.index for component in self.components)


def analyze_component_action(
    *,
    degree: int,
    generators: Mapping[str, Sequence[int]],
    components: Sequence[Mapping[str, object]],
    max_group_order: int = 10_000,
    compute_gamma: bool = True,
) -> NormalCoverAnalysis:
    parsed_generators = tuple(
        validate_permutation(values, degree, name)
        for name, values in generators.items()
    )
    group = generated_group(parsed_generators, max_order=max_group_order)
    component_records: list[ComponentCover] = []
    seen_points: set[int] = set()
    for entry in components:
        name = str(entry["name"])
        basepoint = int(entry["basepoint"])
        component_orbit = orbit(group, basepoint)
        overlap = seen_points.intersection(component_orbit)
        if overlap:
            raise ValueError(f"{name}: component orbit overlaps at {sorted(overlap)}")
        seen_points.update(component_orbit)
        component_stabilizer = stabilizer(group, basepoint)
        component_records.append(
            ComponentCover(
                name=name,
                basepoint=basepoint,
                orbit=component_orbit,
                stabilizer=component_stabilizer,
                core=subgroup_core(group, component_stabilizer),
            )
        )
    if seen_points != set(range(degree)):
        raise ValueError("component orbits do not partition the action")
    covered = frozenset().union(
        *(
            conjugate_union(group, component.stabilizer)
            for component in component_records
        )
    )
    common_core = frozenset.intersection(
        *(component.core for component in component_records)
    )
    gamma = normal_covering_number(group) if compute_gamma else None
    return NormalCoverAnalysis(
        degree=degree,
        group=group,
        components=tuple(component_records),
        covered_elements=covered,
        common_core=common_core,
        normal_covering_number=gamma,
    )
