#!/usr/bin/env python3
"""Test a p-typical Witt recursion on exposed binary quotient packets.

For a color-homogeneous scroll packet, fix color counts ``d_j``, a
total level ``w``, and color polynomials

    C_j(z) = sum_i c_(j,i) z^i.

At scale ``N`` its ordinary packet, before removal of the common radial
factorial, is

    M_N = (N*w)! (N*(capacity-w))!
          [z^(N*w)] product_j C_j(z)^(N*d_j).

Choose a prime p larger than ``w`` and ``capacity-w``.  The common
factorial valuation at scale p^k is

    capacity * (1+p+...+p^(k-1)).

After removing it and its Wilson sign, the proposed first two
p-typical ghost coordinates are

    X_1 = (G_1-G_0^p)/p,
    X_2 = (G_2-G_0^(p^2)-p*X_1^p)/p^2.

The first divisibility is the known Frobenius/one-carry congruence.  The
second is the first falsifiable Witt--Rees strictness test.  We run it
on the support-five quotient obstruction, the minimal S(6) and S(5,4)
exceptional scroll fibers, and the first larger reversal-symmetric
multi-state fiber.

These are exact integer computations at the displayed coefficient
specializations.  Passing them is a bounded regression, not a proof of
the universal Witt--Rees realization or of GVC(2).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import factorial


@dataclass(frozen=True)
class Color:
    levels: tuple[int, ...]
    count: int
    span: int


@dataclass(frozen=True)
class Packet:
    name: str
    colors: tuple[Color, ...]
    weight: int
    assignments: tuple[tuple[tuple[int, ...], ...], ...]

    @property
    def capacity(self) -> int:
        return sum(color.count * color.span for color in self.colors)


def multiply_by_sparse(
    polynomial: list[int],
    levels: tuple[int, ...],
    coefficients: tuple[int, ...],
) -> list[int]:
    assert len(levels) == len(coefficients)
    result = [0] * (len(polynomial) + max(levels))
    for degree, value in enumerate(polynomial):
        if not value:
            continue
        for level, coefficient in zip(levels, coefficients):
            result[degree + level] += value * coefficient
    return result


def polynomial_power(
    levels: tuple[int, ...],
    coefficients: tuple[int, ...],
    exponent: int,
) -> list[int]:
    result = [1]
    for _ in range(exponent):
        result = multiply_by_sparse(result, levels, coefficients)
    return result


def packet_coefficient(
    packet: Packet,
    scale: int,
    assignment: tuple[tuple[int, ...], ...],
) -> int:
    assert len(assignment) == len(packet.colors)
    color_powers = [
        polynomial_power(
            color.levels,
            coefficients,
            scale * color.count,
        )
        for color, coefficients in zip(packet.colors, assignment)
    ]
    target = scale * packet.weight
    combined = [1]
    for color_power in color_powers:
        new_length = min(
            target + 1,
            len(combined) + len(color_power) - 1,
        )
        new = [0] * new_length
        for left_degree, left in enumerate(combined):
            if not left:
                continue
            largest_right = min(
                len(color_power) - 1,
                target - left_degree,
            )
            for right_degree in range(largest_right + 1):
                right = color_power[right_degree]
                if right:
                    new[left_degree + right_degree] += left * right
        combined = new
    return combined[target] if target < len(combined) else 0


def packet_moment(
    packet: Packet,
    scale: int,
    assignment: tuple[tuple[int, ...], ...],
) -> int:
    complement = packet.capacity - packet.weight
    assert complement >= 0
    coefficient = packet_coefficient(packet, scale, assignment)
    return (
        factorial(scale * packet.weight)
        * factorial(scale * complement)
        * coefficient
    )


def valuation(number: int, prime: int) -> int:
    if number == 0:
        return 10**9
    number = abs(number)
    answer = 0
    while number % prime == 0:
        answer += 1
        number //= prime
    return answer


def geometric_sum(prime: int, height: int) -> int:
    return sum(prime**index for index in range(height))


def normalized_ghost(
    packet: Packet,
    assignment: tuple[tuple[int, ...], ...],
    prime: int,
    height: int,
) -> int:
    scale = prime**height
    moment = packet_moment(packet, scale, assignment)
    baseline = packet.capacity * geometric_sum(prime, height)
    divisor = prime**baseline
    assert moment % divisor == 0
    sign_exponent = baseline
    return (-1) ** sign_exponent * (moment // divisor)


def test_packet(packet: Packet, prime: int) -> tuple[int, int]:
    complement = packet.capacity - packet.weight
    assert prime > max(
        packet.weight,
        complement,
        *(color.count for color in packet.colors),
    )
    minimum_first_residual = 10**9
    minimum_second_residual = 10**9

    for assignment_index, assignment in enumerate(packet.assignments, 1):
        assert all(
            len(coefficients) == len(color.levels)
            for color, coefficients in zip(packet.colors, assignment)
        )
        assert all(
            coefficient % prime
            for coefficients in assignment
            for coefficient in coefficients
        )

        ghost_0 = normalized_ghost(packet, assignment, prime, 0)
        ghost_1 = normalized_ghost(packet, assignment, prime, 1)
        assert (ghost_1 - ghost_0) % prime == 0
        first_residual = ghost_1 - ghost_0**prime
        first_order = valuation(first_residual, prime)
        assert first_order >= 1
        coordinate_1 = first_residual // prime

        ghost_2 = normalized_ghost(packet, assignment, prime, 2)
        gauss_residual = ghost_2 - ghost_1
        gauss_order = valuation(gauss_residual, prime)
        assert gauss_order >= 2
        second_residual = (
            ghost_2
            - ghost_0 ** (prime**2)
            - prime * coordinate_1**prime
        )
        second_order = valuation(second_residual, prime)

        minimum_first_residual = min(
            minimum_first_residual,
            first_order,
        )
        minimum_second_residual = min(
            minimum_second_residual,
            second_order,
        )
        print(
            f"{packet.name} assignment {assignment_index}: "
            f"v_p(first residual)={first_order}, "
            f"v_p(G2-G1)={gauss_order}, "
            f"v_p(second residual)={second_order}"
        )

    return minimum_first_residual, minimum_second_residual


def packets() -> tuple[Packet, ...]:
    return (
        Packet(
            name="support-five",
            colors=(
                Color((0, 3), 1, 3),
                Color((1, 2, 3), 2, 3),
            ),
            weight=6,
            assignments=(
                ((1, 2), (3, 4, 5)),
                ((-1, 3), (2, -4, 5)),
            ),
        ),
        Packet(
            name="S(6)",
            colors=(Color(tuple(range(7)), 3, 6),),
            weight=8,
            assignments=(
                ((1, 2, 3, 4, 5, 6, 7),),
                ((-1, 2, -3, 4, -5, 6, -7),),
            ),
        ),
        Packet(
            name="S(5,4)",
            colors=(
                Color((0, 1, 4, 5), 2, 5),
                Color((0, 1, 4), 2, 4),
            ),
            weight=8,
            assignments=(
                ((1, 2, 3, 4), (5, 6, 7)),
                ((-1, 2, -3, 4), (5, -6, 7)),
            ),
        ),
        Packet(
            name="symmetric-r5",
            colors=(
                Color((0, 5), 1, 5),
                Color((1, 2, 3, 4), 3, 5),
            ),
            weight=8,
            assignments=(
                ((1, 2), (3, 4, 5, 6)),
                ((-1, 2), (3, -4, 5, -6)),
            ),
        ),
    )


def main() -> None:
    prime = 13
    results = {}
    for packet in packets():
        results[packet.name] = test_packet(packet, prime)

    print(f"prime: {prime}")
    print(f"minimum residual valuations: {results}")
    if all(second >= 2 for _, second in results.values()):
        print(
            "STATUS: the first p-squared Witt recursion passes on every "
            "tested exposed packet; universal Rees strictness remains "
            "to be proved"
        )
    else:
        failing = {
            name: orders
            for name, orders in results.items()
            if orders[1] < 2
        }
        print(
            "STATUS: the naive p-typical recursion fails at the second "
            f"ghost coordinate: {failing}"
        )


if __name__ == "__main__":
    main()
