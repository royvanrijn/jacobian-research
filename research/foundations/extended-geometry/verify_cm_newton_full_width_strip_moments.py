#!/usr/bin/env python3
"""Check a full-strip moment and the extra upper-half pole."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_cm_newton_two_input_strip_moments import normalized_moment_residue, power


def build_phase(n, root_multiplicity, prime, lower_terms):
    angular = power({(2, 0): 1, (0, 0): -1}, root_multiplicity, prime)
    phase = {
        (index - 1, n - 2 * (index - 1)): coefficient
        for (index, _), coefficient in angular.items()
    }
    phase[-1, n] = -angular[0, 0]
    phase.update(lower_terms)
    for angle, height in phase:
        assert height >= 0 and height + 2 * angle <= n
        assert (n - 1) * height + 2 * n * angle >= n * (n - 3)
    return phase


def moment_pair(phase, order, prime):
    return tuple(
        normalized_moment_residue(power(phase, exponent, prime), prime)
        for exponent in (order, order + prime)
    )


def check_wider_strip():
    n, a, prime, order = 10, 5, 17, 5
    lower_terms = {(0, 8): 4, (4, 1): 9}
    for angle, height in lower_terms:
        assert (a + 1) * height + 2 * a * angle < (a + 1) * n - 2 * a
    phase = build_phase(n, 3, prime, lower_terms)
    assert moment_pair(phase, order, prime) == (0, 16)
    print("Full-width control: residues (0, 16) modulo 17")


def check_upper_half_pole():
    n, s, prime = 18, 13, 43
    order = (s * prime - 1) // n
    phase = build_phase(n, 4, prime, {(0, n - 1): 6})
    powers = power(phase, order, prime)
    top_height = s * prime - 1
    inputs = (
        powers.get((0, top_height), 0),
        powers.get((prime, top_height - 2 * prime), 0),
        powers.get((0, top_height - prime), 0),
    )
    assert inputs == (0, 5, 13)
    assert moment_pair(phase, order, prime) == (0, 0)
    print("Upper-half correction: inputs (0, 5, 13), residues (0, 0) modulo 43")


if __name__ == "__main__":
    check_wider_strip()
    check_upper_half_pole()
