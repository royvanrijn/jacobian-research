#!/usr/bin/env python3
"""Bounded complete-chain obstruction for the (8,40) triple-root branch.

After translating the only residual triple root in the repeated-tail
(96,144) row, the vertical edge is

    (8,40) -> (8,12).

This module implements the relevant exact rational part of the published
GetGeneratedCorners, GetCornerChildrenList, and GetCompleteChains algorithms.
For positive integral lower endpoints it deliberately uses only the coarse
necessary bound

    b < a and b <= (a-b-1)^2,

plus b=0.  This is a superset of the possible-last-lower-corner list, so an
empty result here is stronger than an empty result after the full filter.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import json
import math


@dataclass(frozen=True)
class Corner:
    """The paper's integral encoding ``(a/l,b)``."""

    a: int
    l: int
    b: int

    def __post_init__(self) -> None:
        if self.a < 0 or self.l <= 0 or self.b < 0:
            raise ValueError("corner coordinates must be nonnegative and l positive")

    @property
    def x(self) -> Fraction:
        return Fraction(self.a, self.l)

    @property
    def diagonal_value(self) -> Fraction:
        return self.x - self.b

    def as_pair(self) -> tuple[str, int]:
        x = str(self.x.numerator)
        if self.x.denominator != 1:
            x += f"/{self.x.denominator}"
        return x, self.b


@dataclass(frozen=True)
class Edge:
    end: Corner
    start: Corner


Chain = tuple[Edge, ...]


def _primitive_direction(
    first_x: Fraction,
    first_y: int | Fraction,
    second_x: Fraction,
    second_y: int | Fraction,
) -> tuple[int, int]:
    """Primitive weight ``(rho,sigma)`` constant on the two points."""

    dx = first_x - second_x
    dy = Fraction(first_y) - Fraction(second_y)
    if dx == 0 and dy == 0:
        raise ValueError("a direction needs two distinct points")
    scale = math.lcm(dx.denominator, dy.denominator)
    integral_dx = int(dx * scale)
    integral_dy = int(dy * scale)
    divisor = math.gcd(abs(integral_dx), abs(integral_dy))
    rho = integral_dy // divisor
    sigma = -integral_dx // divisor
    if rho < 0:
        rho, sigma = -rho, -sigma
    if rho <= 0:
        raise ValueError("complete-chain directions must have rho positive")
    return rho, sigma


def edge_direction(edge: Edge) -> tuple[int, int]:
    return _primitive_direction(
        edge.end.x,
        edge.end.b,
        edge.start.x,
        edge.start.b,
    )


def _valuation(weight: tuple[int, int], corner: Corner) -> Fraction:
    return weight[0] * corner.x + weight[1] * corner.b


def _coarse_possible_lower_corner(corner: Corner) -> bool:
    """A proven superset of integral possible last lower corners."""

    if corner.l != 1 or corner.diagonal_value <= 0:
        return False
    if corner.b == 0:
        return True
    return (
        corner.b < corner.a
        and corner.b <= (corner.a - corner.b - 1) ** 2
    )


def _is_simple(edge: Edge) -> bool:
    rho, sigma = edge_direction(edge)
    gap = rho // math.gcd(rho, edge.end.l)
    ratio = Fraction(rho + sigma, 1) / _valuation((rho, sigma), edge.end)
    f2 = ratio * edge.end.b
    return f2 - 1 == gap and (
        gap > 1 or edge.start.diagonal_value > 0
    )


def generated_corners(edge: Edge) -> tuple[Corner, ...]:
    """Published GetGeneratedCorners with exact rational arithmetic."""

    end, start = edge.end, edge.start
    if start.diagonal_value < 0:
        return (start,)
    if start.diagonal_value == 0:
        return ()

    rho, sigma = edge_direction(edge)
    next_l = math.lcm(end.l, rho)
    gap = rho // math.gcd(rho, end.l)
    gamma_max = min(
        Fraction(end.b - start.b, gap),
        Fraction(end.b - 1),
    )
    if gamma_max.denominator != 1:
        raise ArithmeticError("gamma_max must be integral")
    maximum = int(gamma_max)
    gammas = (
        (maximum,)
        if _is_simple(edge)
        else tuple(range(start.b + 1, maximum + 1))
    )

    result = []
    for gamma in gammas:
        x = end.x + (gamma - end.b) * Fraction(-sigma, rho)
        numerator = x * next_l
        if numerator.denominator != 1:
            raise ArithmeticError("generated corner missed the next lattice")
        corner = Corner(int(numerator), next_l, gamma)
        admissible = (
            corner.diagonal_value < 0
            and (
                Fraction(next_l) - Fraction(corner.a, corner.b) > 1
                or math.gcd(corner.a, corner.b) > 1
            )
        )
        if admissible:
            result.append(corner)
    return tuple(result)


def corner_children(parent: Edge, corner: Corner) -> tuple[Edge, ...]:
    """Published GetCornerChildrenList with a permissive lower-corner gate."""

    parent_weight = edge_direction(parent)
    divisor = math.gcd(corner.a, corner.b)
    lower = math.floor(
        1
        + Fraction(
            divisor * (parent_weight[0] + parent_weight[1]),
            1,
        )
        / _valuation(parent_weight, corner)
    )
    if corner.l == 1:
        upper = divisor
    else:
        upper = math.floor(
            corner.l * (corner.b * corner.l - corner.a)
            + Fraction(divisor, corner.b)
        )

    result = []
    for mu in range(lower, upper + 1):
        f_x = Fraction(mu, divisor) * corner.x
        f_y = Fraction(mu, divisor) * corner.b
        rho, sigma = _primitive_direction(f_x, f_y, Fraction(1), 1)
        gap = rho // math.gcd(rho, corner.l)
        if gap > corner.b or mu % divisor == 0:
            continue
        for index in range(1, corner.b // gap + 1):
            start_x = corner.x + Fraction(index * gap * sigma, rho)
            start_b = corner.b - index * gap
            numerator = start_x * corner.l
            if numerator.denominator != 1:
                raise ArithmeticError("child endpoint missed the current lattice")
            start = Corner(int(numerator), corner.l, start_b)
            if corner.l > 1:
                accepted = start.diagonal_value != 0
            else:
                accepted = (
                    start.diagonal_value < 0
                    or _coarse_possible_lower_corner(start)
                )
            if accepted:
                result.append(Edge(corner, start))
    return tuple(result)


def _is_final(corner: Corner) -> bool:
    return (
        corner.b > 0
        and Fraction(corner.l) - Fraction(corner.a, corner.b) > 1
    )


def children_and_finals(edge: Edge) -> tuple[tuple[Edge, ...], tuple[Corner, ...]]:
    children = []
    finals = []
    for corner in generated_corners(edge):
        if _is_final(corner):
            finals.append(corner)
        else:
            children.extend(corner_children(edge, corner))
    return tuple(children), tuple(finals)


def _omega(value: int) -> int:
    factors = 0
    prime = 2
    while prime * prime <= value:
        while value % prime == 0:
            factors += 1
            value //= prime
        prime += 1
    if value > 1:
        factors += 1
    return factors


def maximum_chain_length(initial: Edge) -> int:
    rho, _ = edge_direction(initial)
    quotient = Fraction(initial.end.b - initial.start.b, rho)
    if quotient.denominator != 1:
        raise ArithmeticError("the complete-chain length input must be integral")
    return _omega(math.gcd(initial.end.b, int(quotient))) + 1


def complete_chains(
    initial: Edge,
) -> tuple[tuple[Chain, ...], tuple[int, ...]]:
    """Enumerate through the theorem's complete maximum length."""

    open_chains: tuple[Chain, ...] = ((initial,),)
    complete: list[Chain] = []
    open_counts = [1]
    for _ in range(maximum_chain_length(initial)):
        next_open = []
        for chain in open_chains:
            children, finals = children_and_finals(chain[-1])
            complete.extend(chain + (Edge(final, final),) for final in finals)
            next_open.extend(chain + (child,) for child in children)
        open_chains = tuple(next_open)
        open_counts.append(len(open_chains))
    return tuple(complete), tuple(open_counts)


@dataclass(frozen=True)
class NoEscapeAudit:
    translated_initial_edge: tuple[tuple[str, int], tuple[str, int]]
    maximum_length: int
    open_chain_counts: tuple[int, ...]
    first_child_endpoints: tuple[tuple[str, int], ...]
    complete_chain_count: int
    permissive_lower_corner_gate: str
    published_fixture_reproduced: bool
    excluded: bool

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def triple_root_no_escape_audit() -> NoEscapeAudit:
    """Certify that the translated triple-root edge has no complete chain."""

    translated = Edge(Corner(8, 1, 40), Corner(8, 1, 12))
    complete, counts = complete_chains(translated)
    first_corner = generated_corners(translated)
    assert first_corner == (Corner(8, 1, 12),)
    first_children = corner_children(translated, first_corner[0])

    # Regression against the published companion row.  The algorithm must
    # recover its displayed final corner before its empty result is trusted.
    fixture = Edge(Corner(8, 1, 32), Corner(8, 1, 28))
    fixture_complete, _ = complete_chains(fixture)
    published_final = Corner(11, 4, 7)
    fixture_reproduced = any(
        chain[-1].end == published_final for chain in fixture_complete
    )

    if (
        maximum_chain_length(translated) != 3
        or counts != (1, 6, 3, 0)
        or complete
        or not fixture_reproduced
    ):
        raise AssertionError("unexpected complete-chain continuation")

    return NoEscapeAudit(
        translated_initial_edge=(
            translated.end.as_pair(),
            translated.start.as_pair(),
        ),
        maximum_length=3,
        open_chain_counts=counts,
        first_child_endpoints=tuple(
            child.start.as_pair() for child in first_children
        ),
        complete_chain_count=0,
        permissive_lower_corner_gate=(
            "b=0 or (b<a and b<=(a-b-1)^2); this contains every true "
            "possible last lower corner"
        ),
        published_fixture_reproduced=fixture_reproduced,
        excluded=True,
    )


if __name__ == "__main__":
    print(json.dumps(triple_root_no_escape_audit().as_dict(), indent=2, sort_keys=True))
