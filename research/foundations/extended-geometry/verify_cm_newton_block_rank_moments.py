#!/usr/bin/env python3
"""Replay two signature bounds and a marked rank-two example."""

from math import comb


def coherent_dimension(index, conductor):
    n = conductor // 6
    return (index + 2) // 3 - index // 6 - index * (n + 1) // conductor


def check_rank_bound(conductor, multiplier, genus, expected_rank, expected_maximum):
    dimensions = [coherent_dimension(i, conductor) for i in range(conductor)]
    assert sum(dimensions) == genus
    rank_bound = 0
    rank_maximum = 0
    for index in range(1, conductor):
        target = multiplier * index % conductor
        de_rham_dimension = dimensions[index] + dimensions[-index]
        assert de_rham_dimension == dimensions[target] + dimensions[-target]
        rank_bound += max(0, dimensions[index] + dimensions[target] - de_rham_dimension)
        rank_maximum += min(dimensions[index], dimensions[target])
        assert (dimensions[index], dimensions[target]) != (2, 2)
    assert rank_bound == expected_rank and rank_maximum == expected_maximum
    print(f"Conductor {conductor}, genus {genus}: rank at least {rank_bound}, generic rank {rank_maximum}")


def coefficient(index, positive_factor, negative_factor, prime):
    lo = max(0, index - negative_factor)
    hi = min(positive_factor, index)
    return sum(
        (-1)**(positive_factor-j) * comb(positive_factor, j) * comb(negative_factor, index-j)
        for j in range(lo, hi+1)
    ) % prime


def check_two_dimensional_target():
    p = 229
    assert (5*p-11)//162 == 7 and (5*p-11)%162 == 0
    assert [coherent_dimension(i, 162) for i in (5, 11)] == [2, 2]
    entries = [
        tuple(coefficient(8-i+p*j, 190, 197, p) for j in range(2))
        for i in (1, 2)
    ]
    assert entries == [(148, 141), (16, 27)]
    (a, b), (c, d) = entries
    assert (a*d-b*c) % p == 137
    pivot = 2 * pow(167*165, -1, p) % p
    assert (8*237 * (a*d-b*c) * pow(5, -1, p) * pivot) % p == 26
    assert all((16*first * pow(5, -1, p)) % p for first, _ in entries)
    marker_counts = []
    for count in range(9):
        top = 1296 - 152*count
        bottom = top - 2*(8-count)
        if any((height+1) % p == 0 for height in range(bottom, top+1)):
            marker_counts.append(count)
    assert marker_counts == [1, 4]
    assert (8+p)*162+1 < p*p
    print("At 229: the reconstructed coherent block has determinant 137")


if __name__ == "__main__":
    check_rank_bound(54, 11, 40, 8, 30)
    check_rank_bound(162, 65, 121, 30, 90)
    check_two_dimensional_target()
