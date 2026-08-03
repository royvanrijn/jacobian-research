"""Support-first census machinery for normalized Keller collisions in A^3.

The module deliberately separates finite combinatorics from coefficient
algebra.  A collision frame has

    F(0) = 0,  JF(0) = I,  F(e_1) = 0,

so the three linear monomials are fixed and every other displayed monomial
has its own nonzero coefficient.  For fixed supports, the coefficient of a
Jacobian-determinant monomial is grouped without expanding a CAS expression:

    exponent = alpha + beta + gamma - (1, 1, 1),
    multiplier = det(alpha, beta, gamma).

A nonconstant bucket with one contribution is impossible on the exact
coefficient torus.  ``SupportEnumerator`` uses that observation as a closure
rule; unlike a coefficient-box search, every branch adds a monomial needed
to balance a currently singleton bucket.

Integer weights are infinite.  ``enumerate_valuation_faces`` instead returns
the finite equivalence classes determined by exposed Newton faces, including
all seven coordinate strata.  This is the finite object meant by a valuation
census.

The code is specific to three variables only where the determinant triples
are compiled.  It makes no weighted or equivariance assumption.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from typing import Any, Iterable, Mapping, Sequence

import sympy as sp
import z3


Exponent = tuple[int, int, int]
Atom = tuple[int, Exponent]

ZERO: Exponent = (0, 0, 0)
ONES: Exponent = (1, 1, 1)
LINEAR: tuple[Exponent, Exponent, Exponent] = (
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)


def _canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def sha256_json(payload: Any) -> str:
    """Return a stable SHA-256 digest of a JSON-compatible object."""

    return hashlib.sha256(_canonical_json(payload).encode()).hexdigest()


def total_degree(exponent: Exponent) -> int:
    return sum(exponent)


def determinant3(alpha: Exponent, beta: Exponent, gamma: Exponent) -> int:
    """Determinant of the matrix with the three exponents as rows."""

    return (
        alpha[0] * (beta[1] * gamma[2] - beta[2] * gamma[1])
        - alpha[1] * (beta[0] * gamma[2] - beta[2] * gamma[0])
        + alpha[2] * (beta[0] * gamma[1] - beta[1] * gamma[0])
    )


def add_exponents(*exponents: Exponent) -> Exponent:
    return tuple(sum(exponent[index] for exponent in exponents) for index in range(3))  # type: ignore[return-value]


def nonlinear_monomials(max_degree: int) -> tuple[Exponent, ...]:
    """All ordinary monomials of degrees two through ``max_degree``."""

    if max_degree < 2:
        return ()
    result = [
        (x_degree, y_degree, z_degree)
        for x_degree in range(max_degree + 1)
        for y_degree in range(max_degree + 1 - x_degree)
        for z_degree in range(max_degree + 1 - x_degree - y_degree)
        if 2 <= x_degree + y_degree + z_degree <= max_degree
    ]
    return tuple(sorted(result, key=lambda exponent: (sum(exponent), exponent)))


def degree_profiles_below(target: Sequence[int]) -> tuple[tuple[int, int, int], ...]:
    """Enumerate nonincreasing positive degree flags lexicographically below target."""

    if len(target) != 3:
        raise ValueError("a three-dimensional degree profile has three entries")
    target_tuple = tuple(int(value) for value in target)
    if tuple(sorted(target_tuple, reverse=True)) != target_tuple:
        raise ValueError("the target degree profile must be nonincreasing")
    profiles = []
    for first in range(1, target_tuple[0] + 1):
        for second in range(1, first + 1):
            for third in range(1, second + 1):
                profile = (first, second, third)
                if profile < target_tuple:
                    profiles.append(profile)
    return tuple(profiles)


def profile_rank_function(profile: Sequence[int], threshold: int) -> int:
    """Rank of the output coefficient matrix in degrees above threshold.

    If ``V`` is the span of the three coordinate polynomials and its invariant
    degree flag is ``d1 >= d2 >= d3``, this rank is

        dim(V / (V intersect polynomials_of_degree_at_most_threshold)).
    """

    return sum(int(degree) > threshold for degree in profile)


@dataclass(frozen=True, order=True)
class Support:
    """The exact nonlinear support of a normalized three-coordinate map."""

    rows: tuple[tuple[Exponent, ...], tuple[Exponent, ...], tuple[Exponent, ...]]

    def __post_init__(self) -> None:
        if len(self.rows) != 3:
            raise ValueError("a support needs three coordinate rows")
        normalized = tuple(
            tuple(sorted(set(row), key=lambda exponent: (sum(exponent), exponent)))
            for row in self.rows
        )
        if normalized != self.rows:
            raise ValueError("support rows must be sorted and duplicate-free")
        for row in self.rows:
            for exponent in row:
                if len(exponent) != 3 or min(exponent) < 0:
                    raise ValueError(f"invalid exponent {exponent}")
                if total_degree(exponent) < 2:
                    raise ValueError("Support stores nonlinear exponents only")

    @classmethod
    def from_atoms(cls, atoms: Iterable[Atom]) -> "Support":
        rows: list[list[Exponent]] = [[], [], []]
        for component, exponent in atoms:
            if component not in range(3):
                raise ValueError(f"invalid component {component}")
            rows[component].append(tuple(exponent))  # type: ignore[arg-type]
        ordered = tuple(
            tuple(sorted(set(row), key=lambda exponent: (sum(exponent), exponent)))
            for row in rows
        )
        return cls(ordered)  # type: ignore[arg-type]

    @property
    def atoms(self) -> tuple[Atom, ...]:
        return tuple(
            (component, exponent)
            for component, row in enumerate(self.rows)
            for exponent in row
        )

    @property
    def nonlinear_size(self) -> int:
        return sum(map(len, self.rows))

    @property
    def max_degree(self) -> int:
        return max((total_degree(exponent) for _, exponent in self.atoms), default=1)

    @property
    def full_rows(self) -> tuple[tuple[Exponent, ...], ...]:
        return tuple(
            (LINEAR[component], *row) for component, row in enumerate(self.rows)
        )

    @property
    def identifier(self) -> str:
        return "s-" + sha256_json(self.to_json())[:16]

    def to_json(self) -> list[list[list[int]]]:
        return [
            [list(exponent) for exponent in row]
            for row in self.rows
        ]

    def collision_axis_counts(self, axis: int = 0) -> tuple[int, int, int]:
        if axis not in range(3):
            raise ValueError("the collision axis must be 0, 1, or 2")
        return tuple(
            sum(
                all(exponent[index] == 0 for index in range(3) if index != axis)
                for exponent in row
            )
            for row in self.rows
        )  # type: ignore[return-value]

    def collision_support_possible(self, axis: int = 0) -> bool:
        """Coefficient-free necessity for ``F(e_axis)=0`` on an exact support.

        The axis component needs at least one nonlinear pure-axis term to
        cancel its fixed linear value one.  Every other component may have no
        pure-axis terms, or at least two; exactly one would force its declared
        nonzero coefficient to vanish.
        """

        counts = self.collision_axis_counts(axis)
        return counts[axis] >= 1 and all(
            counts[component] != 1 for component in range(3) if component != axis
        )

    def swapped_23(self) -> "Support":
        """Residual simultaneous source/target swap fixing the e1 collision."""

        def swap_exponent(exponent: Exponent) -> Exponent:
            return (exponent[0], exponent[2], exponent[1])

        return Support(
            (
                tuple(sorted(map(swap_exponent, self.rows[0]), key=lambda e: (sum(e), e))),
                tuple(sorted(map(swap_exponent, self.rows[2]), key=lambda e: (sum(e), e))),
                tuple(sorted(map(swap_exponent, self.rows[1]), key=lambda e: (sum(e), e))),
            )
        )

    def canonical_under_collision_stabilizer(self) -> "Support":
        swapped = self.swapped_23()
        return min(self, swapped)


@dataclass(frozen=True)
class BucketTerm:
    alpha: Exponent
    beta: Exponent
    gamma: Exponent
    multiplier: int

    @property
    def coefficient_atoms(self) -> tuple[Atom, ...]:
        return tuple(
            (component, exponent)
            for component, exponent in enumerate((self.alpha, self.beta, self.gamma))
            if exponent != LINEAR[component]
        )

    def to_json(self) -> dict[str, Any]:
        return {
            "alpha": list(self.alpha),
            "beta": list(self.beta),
            "gamma": list(self.gamma),
            "multiplier": self.multiplier,
        }


@dataclass(frozen=True)
class DeterminantBucket:
    exponent: Exponent
    terms: tuple[BucketTerm, ...]

    def to_json(self) -> dict[str, Any]:
        return {
            "exponent": list(self.exponent),
            "terms": [term.to_json() for term in self.terms],
        }


def determinant_buckets(support: Support) -> tuple[DeterminantBucket, ...]:
    """Compile the exact Jacobian-determinant buckets for ``support``."""

    grouped: dict[Exponent, list[BucketTerm]] = defaultdict(list)
    for alpha, beta, gamma in itertools.product(*support.full_rows):
        multiplier = determinant3(alpha, beta, gamma)
        if not multiplier:
            continue
        exponent_sum = add_exponents(alpha, beta, gamma)
        bucket_exponent = tuple(
            exponent_sum[index] - ONES[index] for index in range(3)
        )
        if min(bucket_exponent) < 0:
            raise AssertionError("a nonzero determinant triple gave a Laurent exponent")
        grouped[bucket_exponent].append(
            BucketTerm(alpha, beta, gamma, multiplier)
        )
    return tuple(
        DeterminantBucket(exponent, tuple(grouped[exponent]))
        for exponent in sorted(grouped, key=lambda item: (sum(item), item))
    )


def singleton_buckets(support: Support) -> tuple[DeterminantBucket, ...]:
    return tuple(
        bucket
        for bucket in determinant_buckets(support)
        if bucket.exponent != ZERO and len(bucket.terms) == 1
    )


class SupportEnumerator:
    """Exhaust exact sparse supports using determinant-bucket closure."""

    def __init__(self, max_degree: int, max_nonlinear_support: int) -> None:
        if max_degree < 2:
            raise ValueError("the sparse compiler needs max_degree at least two")
        if max_nonlinear_support < 1:
            raise ValueError("the support bound must be positive")
        self.max_degree = max_degree
        self.max_nonlinear_support = max_nonlinear_support
        self.exponents = nonlinear_monomials(max_degree)
        self.atoms: tuple[Atom, ...] = tuple(
            (component, exponent)
            for component in range(3)
            for exponent in self.exponents
        )
        self.atom_index = {atom: index for index, atom in enumerate(self.atoms)}
        self.axis_atoms = tuple(
            tuple(
                self.atom_index[(component, exponent)]
                for exponent in self.exponents
                if exponent[1] == 0 and exponent[2] == 0
            )
            for component in range(3)
        )
        self.completions, self.universe_statistics = self._compile_completions()
        self._bucket_cache: dict[frozenset[int], tuple[DeterminantBucket, ...]] = {}

    def _compile_completions(
        self,
    ) -> tuple[dict[Exponent, tuple[tuple[tuple[int, ...], int], ...]], dict[str, int]]:
        rows: list[list[tuple[Exponent, int | None]]] = []
        for component in range(3):
            row: list[tuple[Exponent, int | None]] = [(LINEAR[component], None)]
            row.extend(
                (exponent, self.atom_index[(component, exponent)])
                for exponent in self.exponents
            )
            rows.append(row)

        grouped: dict[Exponent, list[tuple[tuple[int, ...], int]]] = defaultdict(list)
        triple_count = 0
        for alpha, beta, gamma in itertools.product(*rows):
            multiplier = determinant3(alpha[0], beta[0], gamma[0])
            if not multiplier:
                continue
            exponent_sum = add_exponents(alpha[0], beta[0], gamma[0])
            bucket_exponent = tuple(exponent_sum[index] - 1 for index in range(3))
            required = tuple(
                atom
                for atom in (alpha[1], beta[1], gamma[1])
                if atom is not None
            )
            grouped[bucket_exponent].append((required, multiplier))
            triple_count += 1

        frozen = {key: tuple(value) for key, value in grouped.items()}
        statistics = {
            "nonzero_determinant_triples": triple_count,
            "determinant_bucket_count": len(frozen),
            "maximum_bucket_size": max(map(len, frozen.values())),
        }
        return frozen, statistics

    def _support_from_indices(self, selected: frozenset[int]) -> Support:
        return Support.from_atoms(self.atoms[index] for index in sorted(selected))

    def _buckets(self, selected: frozenset[int]) -> tuple[DeterminantBucket, ...]:
        cached = self._bucket_cache.get(selected)
        if cached is None:
            cached = determinant_buckets(self._support_from_indices(selected))
            self._bucket_cache[selected] = cached
        return cached

    def _singleton_exponents(self, selected: frozenset[int]) -> tuple[Exponent, ...]:
        return tuple(
            bucket.exponent
            for bucket in self._buckets(selected)
            if bucket.exponent != ZERO and len(bucket.terms) == 1
        )

    def _collision_offenders(self, selected: frozenset[int]) -> tuple[int, ...]:
        return tuple(
            component
            for component in (1, 2)
            if sum(atom in selected for atom in self.axis_atoms[component]) == 1
        )

    def enumerate(self) -> tuple[dict[int, tuple[Support, ...]], dict[str, int]]:
        """Return every no-singleton, collision-admissible support through the bound.

        Completeness is by induction inside any hypothetical final support.
        A currently active singleton must acquire another triple in its bucket;
        every such triple is a branch below.  Once balanced, an arbitrary next
        monomial is allowed, so balanced supersets are not lost.  A temporary
        single pure-axis term in component two or three is repaired by branching
        over every possible second axis term.
        """

        bound = self.max_nonlinear_support
        seen: set[frozenset[int]] = set()
        balanced: dict[int, set[Support]] = defaultdict(set)

        def visit(selected: frozenset[int]) -> None:
            if selected in seen or len(selected) > bound:
                return
            seen.add(selected)
            singleton_exponents = self._singleton_exponents(selected)
            if singleton_exponents:
                best_options: set[frozenset[int]] | None = None
                best_key: tuple[int, int, Exponent] | None = None
                for exponent in singleton_exponents:
                    options = {
                        frozenset(required).difference(selected)
                        for required, _multiplier in self.completions[exponent]
                        if frozenset(required).difference(selected)
                        and len(selected.union(required)) <= bound
                    }
                    if not options:
                        return
                    key = (len(options), sum(map(len, options)), exponent)
                    if best_key is None or key < best_key:
                        best_key = key
                        best_options = options
                assert best_options is not None
                for missing in sorted(best_options, key=lambda item: (len(item), tuple(sorted(item)))):
                    visit(selected.union(missing))
                return

            offenders = self._collision_offenders(selected)
            if offenders:
                component = offenders[0]
                for atom in self.axis_atoms[component]:
                    if atom not in selected:
                        visit(selected.union((atom,)))
                return

            support = self._support_from_indices(selected)
            if not support.collision_support_possible():
                raise AssertionError("collision repair terminated in an invalid support")
            balanced[len(selected)].add(support)
            if len(selected) == bound:
                return
            for atom in range(len(self.atoms)):
                if atom not in selected:
                    visit(selected.union((atom,)))

        for atom in self.axis_atoms[0]:
            visit(frozenset((atom,)))

        result = {
            size: tuple(sorted(supports))
            for size, supports in sorted(balanced.items())
        }
        return result, {
            "visited_partial_supports": len(seen),
            "cached_bucket_compilations": len(self._bucket_cache),
        }


def support_orbits(supports: Iterable[Support]) -> tuple[dict[str, Any], ...]:
    grouped: dict[Support, set[Support]] = defaultdict(set)
    for support in supports:
        grouped[support.canonical_under_collision_stabilizer()].add(support)
    return tuple(
        {
            "support_id": representative.identifier,
            "support": representative.to_json(),
            "orbit_size": len(members),
            "member_ids": sorted(member.identifier for member in members),
        }
        for representative, members in sorted(grouped.items())
    )


def collision_support_space_counts(max_degree: int) -> dict[str, int]:
    """Exact Boolean support-space and residual S2-orbit counts."""

    exponents = nonlinear_monomials(max_degree)
    axis_count = sum(exponent[1] == exponent[2] == 0 for exponent in exponents)
    off_axis_count = len(exponents) - axis_count
    first_count = (2**axis_count - 1) * 2**off_axis_count
    other_count = (2**axis_count - axis_count) * 2**off_axis_count
    raw = first_count * other_count**2

    fixed_exponents = sum(exponent[1] == exponent[2] for exponent in exponents)
    exponent_orbits = fixed_exponents + (len(exponents) - fixed_exponents) // 2
    fixed_first = (2**axis_count - 1) * 2 ** (exponent_orbits - axis_count)
    fixed = fixed_first * other_count
    return {
        "nonlinear_monomials_per_component": len(exponents),
        "optional_support_atoms": 3 * len(exponents),
        "pure_collision_axis_monomials_per_component": axis_count,
        "collision_admissible_labelled_supports": raw,
        "supports_fixed_by_residual_swap": fixed,
        "collision_stabilizer_orbits": (raw + fixed) // 2,
    }


def collision_support_counts_by_size(
    max_degree: int,
    maximum_size: int,
) -> dict[int, int]:
    """Count labelled collision-admissible supports by nonlinear size."""

    exponents = nonlinear_monomials(max_degree)
    axis_count = sum(exponent[1] == exponent[2] == 0 for exponent in exponents)
    off_axis_count = len(exponents) - axis_count

    def component_polynomial(axis_rule: str) -> list[int]:
        values = [0] * (maximum_size + 1)
        for axis_terms in range(min(axis_count, maximum_size) + 1):
            if axis_rule == "positive" and axis_terms == 0:
                continue
            if axis_rule == "not_one" and axis_terms == 1:
                continue
            for off_axis_terms in range(
                min(off_axis_count, maximum_size - axis_terms) + 1
            ):
                values[axis_terms + off_axis_terms] += (
                    math.comb(axis_count, axis_terms)
                    * math.comb(off_axis_count, off_axis_terms)
                )
        return values

    def multiply(left: Sequence[int], right: Sequence[int]) -> list[int]:
        result = [0] * (maximum_size + 1)
        for left_degree, left_value in enumerate(left):
            for right_degree, right_value in enumerate(right):
                if left_degree + right_degree <= maximum_size:
                    result[left_degree + right_degree] += left_value * right_value
        return result

    first = component_polynomial("positive")
    other = component_polynomial("not_one")
    product = multiply(multiply(first, other), other)
    return {size: product[size] for size in range(1, maximum_size + 1)}


def bucket_summary(support: Support) -> dict[str, Any]:
    buckets = determinant_buckets(support)
    nonconstant = tuple(bucket for bucket in buckets if bucket.exponent != ZERO)
    payload = [bucket.to_json() for bucket in buckets]
    histogram: dict[str, int] = defaultdict(int)
    for bucket in nonconstant:
        histogram[str(len(bucket.terms))] += 1
    return {
        "support_id": support.identifier,
        "bucket_count": len(buckets),
        "nonconstant_bucket_count": len(nonconstant),
        "singleton_bucket_count": sum(len(bucket.terms) == 1 for bucket in nonconstant),
        "contribution_count_histogram": dict(sorted(histogram.items(), key=lambda item: int(item[0]))),
        "bucket_sha256": sha256_json(payload),
        "buckets": payload,
    }


def _active_coordinate_masks() -> tuple[tuple[bool, bool, bool], ...]:
    return tuple(
        tuple(bool(mask & (1 << index)) for index in range(3))  # type: ignore[misc]
        for mask in range(1, 8)
    )


def _restricted_row(row: Sequence[Exponent], active: Sequence[bool]) -> tuple[Exponent, ...]:
    return tuple(
        exponent
        for exponent in row
        if all(active[index] or exponent[index] == 0 for index in range(3))
    )


def _canonical_weight_for_faces(
    rows: Sequence[Sequence[Exponent]],
    active: Sequence[bool],
    faces: Sequence[Sequence[int]],
) -> tuple[int, int, int]:
    weights = tuple(z3.Int(f"cw{index}") for index in range(3))
    maxima = tuple(z3.Int(f"cm{index}") for index in range(3))
    optimizer = z3.Optimize()
    for index, is_active in enumerate(active):
        if not is_active:
            optimizer.add(weights[index] == 0)
    optimizer.add(z3.Or(*(weights[index] > 0 for index in range(3) if active[index])))

    for component, row in enumerate(rows):
        if not row:
            optimizer.add(maxima[component] == 0)
            continue
        face = set(faces[component])
        if not face:
            raise AssertionError("a nonempty restricted support needs an exposed face")
        for monomial_index, exponent in enumerate(row):
            dot = sum(weights[index] * exponent[index] for index in range(3))
            if monomial_index in face:
                optimizer.add(dot == maxima[component])
            else:
                optimizer.add(dot <= maxima[component] - 1)

    absolute = tuple(z3.Int(f"ca{index}") for index in range(3))
    for index in range(3):
        optimizer.add(absolute[index] >= weights[index])
        optimizer.add(absolute[index] >= -weights[index])
        optimizer.add(absolute[index] >= 0)
    optimizer.set(priority="lex")
    optimizer.minimize(sum(absolute))
    for weight in weights:
        optimizer.minimize(weight)
    if optimizer.check() != z3.sat:
        raise AssertionError("a face pattern lost its valuation representative")
    model = optimizer.model()
    result = tuple(model.eval(weight).as_long() for weight in weights)
    common = math.gcd(*map(abs, result))
    if common > 1:
        result = tuple(value // common for value in result)
    return result  # type: ignore[return-value]


def enumerate_valuation_faces(support: Support) -> tuple[dict[str, Any], ...]:
    """Enumerate all infinity weights modulo equality of exposed faces.

    Coordinates that vanish identically on the escaping curve are represented
    by a coordinate stratum.  On every nonempty restricted component support,
    a positive maximum must be attained at least twice; otherwise that output
    cannot remain bounded along a torus-leading Puiseux arc.
    """

    records: list[dict[str, Any]] = []
    for active in _active_coordinate_masks():
        rows = tuple(_restricted_row(row, active) for row in support.full_rows)
        weights = tuple(z3.Int(f"w{index}") for index in range(3))
        maxima = tuple(z3.Int(f"m{index}") for index in range(3))
        solver = z3.Solver()
        for index, is_active in enumerate(active):
            if not is_active:
                solver.add(weights[index] == 0)
        solver.add(z3.Or(*(weights[index] > 0 for index in range(3) if active[index])))

        selectors: list[z3.BoolRef] = []
        row_selectors: list[list[z3.BoolRef]] = []
        for component, row in enumerate(rows):
            selected: list[z3.BoolRef] = []
            for monomial_index, exponent in enumerate(row):
                selector = z3.Bool(f"face_{component}_{monomial_index}")
                dot = sum(weights[index] * exponent[index] for index in range(3))
                solver.add(maxima[component] >= dot)
                solver.add(selector == (maxima[component] == dot))
                selected.append(selector)
                selectors.append(selector)
            if selected:
                count = z3.Sum(*(z3.If(selector, 1, 0) for selector in selected))
                solver.add(z3.Or(*selected))
                solver.add(z3.Implies(maxima[component] > 0, count >= 2))
            else:
                solver.add(maxima[component] == 0)
            row_selectors.append(selected)

        while solver.check() == z3.sat:
            model = solver.model()
            pattern = tuple(
                z3.is_true(model.eval(selector, model_completion=True))
                for selector in selectors
            )
            faces: list[tuple[int, ...]] = []
            offset = 0
            for selected in row_selectors:
                faces.append(
                    tuple(
                        index
                        for index in range(len(selected))
                        if pattern[offset + index]
                    )
                )
                offset += len(selected)

            weight = _canonical_weight_for_faces(rows, active, faces)
            face_exponents = [
                [list(row[index]) for index in face]
                for row, face in zip(rows, faces)
            ]
            maxima_values = [
                max(
                    sum(weight[index] * exponent[index] for index in range(3))
                    for exponent in row
                )
                if row
                else None
                for row in rows
            ]
            records.append(
                {
                    "active_coordinates": [
                        index + 1 for index, value in enumerate(active) if value
                    ],
                    "weight": [weight[index] if active[index] else None for index in range(3)],
                    "maxima": maxima_values,
                    "faces": face_exponents,
                }
            )
            if selectors:
                solver.add(
                    z3.Or(
                        *(
                            selector != z3.BoolVal(value)
                            for selector, value in zip(selectors, pattern)
                        )
                    )
                )
            else:
                break

    unique = {_canonical_json(record): record for record in records}
    return tuple(unique[key] for key in sorted(unique))


def sign_smt(support: Support) -> dict[str, Any]:
    """Necessary real/rational sign feasibility without polynomial arithmetic."""

    variables = {
        atom: z3.Bool(
            "positive_" + str(atom[0] + 1) + "_" + "_".join(map(str, atom[1]))
        )
        for atom in support.atoms
    }
    solver = z3.Solver()

    def term_positive(term: BucketTerm) -> z3.BoolRef:
        negative_bits = [z3.Not(variables[atom]) for atom in term.coefficient_atoms]
        if len(negative_bits) > 1:
            negative: z3.BoolRef = z3.Xor(*negative_bits)
        elif negative_bits:
            negative = negative_bits[0]
        else:
            negative = z3.BoolVal(False)
        if term.multiplier < 0:
            negative = z3.Not(negative)
        return z3.Not(negative)

    for bucket in determinant_buckets(support):
        if bucket.exponent == ZERO:
            continue
        signs = tuple(term_positive(term) for term in bucket.terms)
        solver.add(z3.Or(*signs))
        solver.add(z3.Or(*(z3.Not(sign) for sign in signs)))

    for component, row in enumerate(support.rows):
        axis_atoms = tuple(
            (component, exponent)
            for exponent in row
            if exponent[1] == exponent[2] == 0
        )
        if component == 0:
            solver.add(z3.Or(*(z3.Not(variables[atom]) for atom in axis_atoms)))
        elif axis_atoms:
            solver.add(z3.Or(*(variables[atom] for atom in axis_atoms)))
            solver.add(z3.Or(*(z3.Not(variables[atom]) for atom in axis_atoms)))

    status = solver.check()
    result: dict[str, Any] = {
        "support_id": support.identifier,
        "status": str(status),
        "scope": "necessary over ordered coefficient fields; not a complex-field exclusion",
    }
    if status == z3.sat:
        model = solver.model()
        result["sign_model"] = [
            {
                "component": atom[0] + 1,
                "exponent": list(atom[1]),
                "sign": 1 if z3.is_true(model.eval(variable, model_completion=True)) else -1,
            }
            for atom, variable in sorted(variables.items())
        ]
    return result


def coefficient_system(
    support: Support,
) -> tuple[tuple[sp.Symbol, ...], tuple[sp.Expr, ...]]:
    """Construct determinant, collision, and exact-support equations."""

    atoms = support.atoms
    coefficients = sp.symbols("c0:" + str(len(atoms)))
    coefficient_by_atom = dict(zip(atoms, coefficients))
    rho = sp.Symbol("rho")
    equations: list[sp.Expr] = []

    for bucket in determinant_buckets(support):
        if bucket.exponent == ZERO:
            continue
        equations.append(
            sp.Add(
                *(
                    term.multiplier
                    * sp.prod(
                        coefficient_by_atom[atom]
                        for atom in term.coefficient_atoms
                    )
                    for term in bucket.terms
                )
            )
        )

    for component, row in enumerate(support.rows):
        axis_coefficients = [
            coefficient_by_atom[(component, exponent)]
            for exponent in row
            if exponent[1] == exponent[2] == 0
        ]
        if component == 0:
            equations.append(1 + sp.Add(*axis_coefficients))
        elif axis_coefficients:
            equations.append(sp.Add(*axis_coefficients))

    equations.append(rho * sp.prod(coefficients) - 1)
    return (rho, *coefficients), tuple(map(sp.expand, equations))


def groebner_status(support: Support, modulus: int | None = None) -> dict[str, Any]:
    variables, equations = coefficient_system(support)
    kwargs: dict[str, Any] = {"order": "grevlex", "method": "f5b"}
    if modulus is not None:
        kwargs["modulus"] = modulus
    basis = sp.groebner(equations, *variables, **kwargs)
    unit = len(basis.polys) == 1 and basis.polys[0].as_expr() == 1
    result: dict[str, Any] = {
        "support_id": support.identifier,
        "field": "QQ" if modulus is None else f"F_{modulus}",
        "unit_ideal": unit,
        "basis_length": len(basis.polys),
    }
    if not unit:
        result["zero_dimensional"] = bool(basis.is_zero_dimensional)
        result["basis_sha256"] = sha256_json(
            [str(polynomial.as_expr()) for polynomial in basis.polys]
        )
    return result


def singular_batch_status(
    supports: Sequence[Support],
    characteristic: int = 0,
    *,
    timeout_seconds: int = 300,
    singular_binary: str | None = None,
) -> tuple[dict[str, Any], ...]:
    """Run exact coefficient-torus Gröbner decisions in one Singular process."""

    if characteristic < 0:
        raise ValueError("the characteristic must be nonnegative")
    binary = singular_binary or shutil.which("Singular")
    if binary is None:
        raise RuntimeError("Singular is required for the batched algebra stage")

    lines = ["option(redSB);"]
    for index, support in enumerate(supports):
        variables, equations = coefficient_system(support)
        ring_name = f"census_ring_{index}"
        variable_list = ",".join(map(str, variables))
        equation_list = ",".join(
            str(equation).replace("**", "^") for equation in equations
        )
        lines.extend(
            [
                f"ring {ring_name}={characteristic},({variable_list}),dp;",
                f"ideal census_ideal={equation_list};",
                "ideal census_basis=std(census_ideal);",
                "int census_unit=(reduce(1,census_basis)==0);",
                "int census_dimension=-1;",
                "if (census_unit==0) { census_dimension=dim(census_basis); }",
                (
                    f'print("RESULT|{support.identifier}|"'
                    '+string(census_unit)+"|"+string(size(census_basis))'
                    '+"|"+string(census_dimension));'
                ),
                f"kill {ring_name};",
            ]
        )
    lines.append("exit;")
    completed = subprocess.run(
        [binary, "-q"],
        input="\n".join(lines) + "\n",
        text=True,
        capture_output=True,
        timeout=timeout_seconds,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "Singular batch failed:\n"
            + completed.stdout[-4000:]
            + completed.stderr[-4000:]
        )

    parsed: dict[str, dict[str, Any]] = {}
    for line in completed.stdout.splitlines():
        if not line.startswith("RESULT|"):
            continue
        _marker, support_id, unit, basis_length, dimension = line.split("|")
        parsed[support_id] = {
            "support_id": support_id,
            "field": "QQ" if characteristic == 0 else f"F_{characteristic}",
            "unit_ideal": unit == "1",
            "basis_length": int(basis_length),
            "dimension": None if dimension == "-1" else int(dimension),
        }
    expected = {support.identifier for support in supports}
    if set(parsed) != expected:
        missing = sorted(expected - set(parsed))
        raise RuntimeError(f"Singular omitted {len(missing)} support results")
    return tuple(parsed[support.identifier] for support in supports)


def dense_quadratic_collision_status() -> dict[str, Any]:
    """Exact dense proof that normalized degree-at-most-two collisions are empty."""

    x, y, z = sp.symbols("x y z")
    source = (x, y, z)
    quadratic = (x**2, x * y, x * z, y**2, y * z, z**2)
    coefficients = sp.symbols("q0:18")
    outputs = sp.Matrix(
        [
            source[component]
            + sum(
                coefficients[6 * component + index] * monomial
                for index, monomial in enumerate(quadratic)
            )
            for component in range(3)
        ]
    )
    determinant = sp.Poly(sp.expand(outputs.jacobian(source).det() - 1), source)
    equations = [coefficient for _monomial, coefficient in determinant.terms()]
    substitutions = {
        coefficients[0]: -1,
        coefficients[6]: 0,
        coefficients[12]: 0,
    }
    equations = [sp.expand(equation.subs(substitutions)) for equation in equations]
    remaining = tuple(
        coefficient for coefficient in coefficients if coefficient not in substitutions
    )
    basis = sp.groebner(equations, *remaining, order="grevlex", method="f5b")
    unit = len(basis.polys) == 1 and basis.polys[0].as_expr() == 1
    return {
        "coordinate_degree_bound": 2,
        "coefficient_count_before_collision": 18,
        "coefficient_count_after_collision": len(remaining),
        "determinant_equation_count": len(equations),
        "unit_ideal_over_QQ": unit,
        "basis": [str(polynomial.as_expr()) for polynomial in basis.polys],
        "eliminated_degree_profiles": [
            list(profile)
            for profile in degree_profiles_below((3, 1, 1))
            if profile[0] <= 2
        ],
    }


def profile_structural_metadata(support: Support, profile: Sequence[int]) -> dict[str, Any]:
    """Describe the later coefficient-rank gates for an invariant degree flag.

    This is metadata rather than a support-only decision: coefficients in
    different output rows may cancel at top degree.  Rank equations and a
    Rabinowitsch equation for the required nonzero minors belong in the exact
    algebra stage.
    """

    rows = support.rows
    thresholds = []
    for threshold in range(1, max(profile) + 1):
        columns = sorted(
            {
                exponent
                for row in rows
                for exponent in row
                if total_degree(exponent) > threshold
            },
            key=lambda exponent: (sum(exponent), exponent),
        )
        zero_pattern = [
            [exponent in row for exponent in columns]
            for row in rows
        ]
        thresholds.append(
            {
                "threshold": threshold,
                "required_rank": profile_rank_function(profile, threshold),
                "columns": [list(exponent) for exponent in columns],
                "zero_pattern": zero_pattern,
            }
        )
    return {
        "support_id": support.identifier,
        "profile": list(profile),
        "threshold_rank_conditions": thresholds,
    }
