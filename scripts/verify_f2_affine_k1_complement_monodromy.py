#!/usr/bin/env python3
"""Verify cyclic-complement certificates on the F2 ``k=1`` chart.

Run this file with ``sage -python``.  SageMath's Zariski--van Kampen
implementation requires the optional ``sirocco`` package.
"""

from __future__ import annotations

from dataclasses import dataclass

try:
    from sage.all import PolynomialRing, QQ
    from sage.schemes.curves.zariski_vankampen import (
        braid_monodromy,
        fundamental_group,
    )
except ImportError as exc:  # pragma: no cover - environment diagnostic
    raise SystemExit(
        "run with `sage -python`; SageMath and its optional sirocco package "
        "are required"
    ) from exc


@dataclass(frozen=True)
class Certificate:
    """One exact parametrized-curve and van Kampen certificate."""

    label: str
    parameters: tuple[int, int, int, int]
    polynomial: object
    braids: tuple[tuple[int, ...], ...]
    relators: tuple[tuple[int, ...], ...]
    resultant_sign: int = 1


def certificates():
    """Return witnesses for all immersed partitions and two cusp strata."""

    ring = PolynomialRing(QQ, names=("P", "Q"))
    P, Q = ring.gens()

    return (
        Certificate(
            "4A1",
            (1, 0, 0, 0),
            -P**5 + 5 * P**2 * Q + 5 * P * Q**2 + Q**3 + Q,
            (
                (2, 1, 2, 1, 2, 1, 2, 2, -1, -2, -1, -2, -1, -2),
                (2, 1, 2, 1, 2, 1, -2, -1, -1, -2),
                (2, 1, 1, 1, 2, 1, 2, -1, -2, -1, -1, -2),
                (2, 1, -2),
                (2, 1, 2, -1, -2),
                (2, 1, 1, -2),
            ),
            (
                (2, 1, -2, -1),
                (-3, 1),
                (2, -1),
                (-3, 1, 3, -1),
                (3, 2, -3, -2),
                (-3, 1, 3, -1),
            ),
            -1,
        ),
        Certificate(
            "A3+2A1",
            (0, 0, 4, -3),
            P**5
            + 12 * P**4
            + 48 * P**3
            - 9 * P**2 * Q
            - Q**3
            + 64 * P**2
            - 36 * P * Q
            - 27 * P,
            (
                (2, 1, 2, -1),
                (2, 1, 2, 1, 2, 1, 2, 1, 2, -1, -2, -2, -1, -2, -1, -2),
                (2, 1, 2, 1, 2, 2, -1, -2, -1, -2),
                (2, 1, 2, 1, 2, 2, 1, 1, 1, 1, -2, -2, -1, -2, -1, -2),
            ),
            (
                (-3, 1, 3, -1),
                (-2, -1, 3, 1, 2, -1, 3, 1, 2, -1, -3, 1, -2, -1, -3, 1),
                (1, -2),
                (2, 1, -2, -1),
                (3, -2),
            ),
        ),
        Certificate(
            "2A3",
            (2, 0, 0, -5),
            P**5
            - 40 * P**3
            - 35 * P**2 * Q
            - 10 * P * Q**2
            - Q**3
            - 5 * P
            - 2 * Q,
            (
                (1, 1, 2, -1, -1),
                (1, 2, 1, 2, 1, 2, -1, -2, 1, 1, -2, -1),
                (1, 2, 1, 2, 1, 2, -1, 2, 1, -2, -1, -2, -1, -2, -1),
                (1, 2, 2, 1, 2, 1, 1, -2, -1, -1),
            ),
            (
                (-3, -2, -1, 2, 1, 2),
                (-3, -2, 1, 2, 3, -2, 1, 2, 3, -2, -1, 2, -3, -2, -1, 2),
                (-3, 1),
                (3, 2, 3, 2, -3, -2, -3, -2),
            ),
        ),
        Certificate(
            "A5+A1",
            (3, -3, -17, -9),
            P**5
            - 87 * P**4
            - 9 * P**3 * Q
            + 363 * P**3
            - 27 * P**2 * Q
            - 15 * P * Q**2
            - Q**3
            - 224 * P**2
            + 432 * P * Q
            + 48 * Q**2
            - 1728 * P
            - 576 * Q,
            (
                (-2, 1, 2, 1, -2),
                (1, 2, 1, -2, -2),
                (1, 2, 1, 2, 2, 2, 2, 2, 2, -1, -2, -1),
                (2, 2),
            ),
            (
                (-3, 1),
                (2, -1),
                (2, 1, 2, 1, 2, 1, -2, -1, -2, -1, -2, -1),
                (3, 2, -3, -2),
            ),
        ),
        Certificate(
            "A7",
            (6, -4, -44, -37),
            P**5
            - 220 * P**4
            - 12 * P**3 * Q
            + 2088 * P**3
            - 147 * P**2 * Q
            - 30 * P * Q**2
            - Q**3
            - 20 * P**2
            + 4332 * P * Q
            + 240 * Q**2
            - 88837 * P
            - 14406 * Q,
            (
                (-2, -2, 1, 2, 1),
                (2, 1, 2, 1, -2, -2, -2),
                (2, 2, 2, 2, 2, 2, 2, 2),
            ),
            (
                (-3, -2, 1, 2),
                (-3, -2, 3, 1),
                (3, 2, 3, 2, 3, 2, 3, 2, -3, -2, -3, -2, -3, -2, -3, -2),
            ),
        ),
        Certificate(
            "A2+3A1",
            (-3, 0, 1, -7),
            P**5
            + 3 * P**4
            + 87 * P**3
            - 66 * P**2 * Q
            + 15 * P * Q**2
            - Q**3
            + 79 * P**2
            - 48 * P * Q
            + 6 * Q**2
            - 7 * P
            + 3 * Q,
            (
                (1, 2, 1, 2, 1, 2, 1, 2, 1, 2, -1, -2, -2, -1, -2, -1, -2, -1),
                (1, 2, 1, 2, 1, 2, 2, 2, -1, -2, -1, -2, -1),
                (1, 2, 1, 2, 1, -2, 1, 1, 2, -1, -2, -1, -2, -1),
                (1, 2, 1, 2, 1, -2, -1, -1),
                (1,),
            ),
            (
                (3, 2, 3, -2, -3, -2),
                (2, -1),
                (2, 1, -2, -1),
                (-3, -2, 1, 2, 3, -2, -1, 2),
                (-3, -2, -1, 2, 3, 2, -3, -2, 1, 2, 3, -2),
            ),
        ),
        Certificate(
            "2A2+2A1",
            (-3, 0, 0, -5),
            P**5
            + 60 * P**3
            - 60 * P**2 * Q
            + 15 * P * Q**2
            - Q**3
            - 80 * P
            + 48 * Q,
            (
                (1, 2, 1, -2, -1, 2),
                (1, 2, 1, -2, 1, 1, 1, 2, -1, -2, -1),
                (1, 2, 1, 2, 1, 2, 1, 2, -1, -2, -2, -1, -2, -1),
                (1, 2, 1, 2, 2, 2, -1, -2, -1),
            ),
            (
                (-3, -2, 1, 2, 3, -2, -1, 2),
                (-3, 1, 3, 1, -3, -1),
                (2, 1, 2, -1, -2, -1),
                (3, 2, -3, -2),
            ),
        ),
    )


def audit_collision_partitions() -> None:
    """Check the five off-diagonal root partitions and the bicusp factor."""

    ring = PolynomialRing(QQ, names=("u",))
    (u,) = ring.gens()
    expected = {
        "4A1": u**4 + u**2 - 1,
        "A3+2A1": (u - 1) ** 2 * (u**2 + 2 * u + 3),
        "2A3": (u**2 + 1) ** 2,
        "A5+A1": u * (u - 1) ** 3,
        "A7": (u - 1) ** 4,
        "2A2+2A1": (u**2 - 4) * (u**2 + 1),
    }
    for certificate in certificates():
        if certificate.label not in expected:
            continue
        a, b, c, d = certificate.parameters
        collision = (
            u**4 + b * u**3 + a * u**2 + (2 * a * b - c) * u - (a**2 + d)
        )
        assert collision == expected[certificate.label]
        diagonal = 3 * u**2 + 4 * a
        if certificate.label == "2A2+2A1":
            assert collision.gcd(diagonal).monic() == u**2 - 4
        else:
            assert collision.gcd(diagonal) == 1

        # Distinct roots of the collision quartic have distinct target
        # values at every witness.  This removes hidden triple-image mergers.
        image_ring = PolynomialRing(QQ, names=("left", "right", "inverse"))
        left, right, inverse_variable = image_ring.gens()
        collision_left = image_ring(collision(left))
        collision_right = image_ring(collision(right))

        def x_value(value):
            return -value * (value**2 + a)

        def y_value(value):
            return (value**2 + a) * (
                value**3 + 2 * a * value + a * b - c
            )

        distinct_images = image_ring.ideal(
            collision_left,
            collision_right,
            x_value(left) - x_value(right),
            y_value(left) - y_value(right),
            inverse_variable * (left - right) - 1,
        )
        assert distinct_images.is_one()

        if certificate.label != "2A2+2A1":
            # At every off-diagonal witness, p itself is a local coordinate
            # on both normalization branches, so root multiplicity equals
            # their intersection multiplicity.
            pair_ring = PolynomialRing(QQ, names=("z", "u"))
            z, pair_u = pair_ring.gens()
            pair = z**2 - pair_u * z + pair_u**2 + a
            critical = 3 * z**2 + a
            critical_u = pair.resultant(critical, z)
            assert ring(critical_u).gcd(collision) == 1


def implicit_equation(certificate: Certificate):
    """Recompute one exact resultant in ``QQ[P,Q]``."""

    ambient = PolynomialRing(QQ, names=("P", "Q", "t"))
    P, Q, t = ambient.gens()
    a, b, c, d = certificate.parameters
    p = t**3 + a * t
    q = t**5 + b * t**4 + c * t**2 + d * t
    resultant = (P - p).resultant(Q - q, t)
    assert resultant == certificate.resultant_sign * ambient(certificate.polynomial)
    return certificate.polynomial


def audit_group(certificate: Certificate) -> None:
    """Check certified braids and both van Kampen presentations."""

    polynomial = implicit_equation(certificate)
    braids, component_map, vertical_map, degree = braid_monodromy(polynomial)
    assert degree == 3
    assert component_map == {0: 0, 1: 0, 2: 0}
    assert vertical_map == {}
    assert sorted(tuple(braid.Tietze()) for braid in braids) == sorted(
        certificate.braids
    )

    raw_group = fundamental_group(
        polynomial,
        simplified=False,
        projective=False,
        puiseux=True,
    )
    assert raw_group.ngens() == 3
    actual_relators = sorted(
        tuple(relation.Tietze()) for relation in raw_group.relations()
    )
    assert actual_relators == sorted(certificate.relators)

    simplified_group = fundamental_group(
        polynomial,
        simplified=True,
        projective=False,
        puiseux=True,
    )
    assert simplified_group.ngens() == 1
    assert tuple(simplified_group.relations()) == ()
    assert tuple(simplified_group.abelian_invariants()) == (0,)


def compose(left, right):
    """Compose zero-based tuple permutations as ``left after right``."""

    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    """Invert a zero-based tuple permutation."""

    result = [0] * len(permutation)
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def evaluate_word(images, word):
    """Evaluate a signed Tietze word in tuple permutations."""

    result = tuple(range(len(images[0])))
    for letter in word:
        image = images[letter - 1] if letter > 0 else inverse(images[-letter - 1])
        result = compose(result, image)
    return result


def audit_first_escape() -> None:
    """Certify the first noncyclic cusp stratum and its degree-six action."""

    collision_ring = PolynomialRing(QQ, names=("u",))
    (u,) = collision_ring.gens()
    assert u**4 + u**3 == u**3 * (u + 1)

    ring = PolynomialRing(QQ, names=("P", "Q"))
    P, Q = ring.gens()
    polynomial = P**5 + P**4 + 3 * P**3 * Q - Q**3
    expected_braids = (
        (1, -2, -1, 2, 1, 2),
        (1, 2, 1, 2, 1, 2, 1, 1, 2, -1),
    )
    expected_relators = (
        (-2, -1, -3, -2, 3, 2, 3, 1),
        (2, 1, -2, -1),
        (2, 3, 1, 3, -1, -3, -2, -1),
    )

    ambient = PolynomialRing(QQ, names=("P", "Q", "t"))
    P3, Q3, t = ambient.gens()
    resultant = (P3 - t**3).resultant(Q3 - (t**5 + t**4), t)
    assert resultant == ambient(polynomial)

    braids, component_map, vertical_map, degree = braid_monodromy(polynomial)
    assert degree == 3
    assert component_map == {0: 0, 1: 0, 2: 0}
    assert vertical_map == {}
    assert sorted(tuple(braid.Tietze()) for braid in braids) == sorted(
        expected_braids
    )
    raw_group = fundamental_group(
        polynomial,
        simplified=False,
        projective=False,
        puiseux=True,
    )
    assert sorted(
        tuple(relation.Tietze()) for relation in raw_group.relations()
    ) == sorted(expected_relators)

    # Each geometric meridian has cycle type 2+2+1+1.  The relations hold,
    # the generated action is transitive, and every meridian fixes two
    # sheets.  Thus complement monodromy alone does not exclude degree six.
    images = (
        (0, 1, 3, 2, 5, 4),
        (0, 1, 4, 5, 2, 3),
        (2, 5, 0, 3, 4, 1),
    )
    identity = tuple(range(6))
    assert all(
        evaluate_word(images, relation) == identity
        for relation in expected_relators
    )
    assert all(
        sum(image[index] == index for index in range(6)) == 2
        for image in images
    )
    assert compose(images[0], images[2]) != compose(images[2], images[0])
    orbit = {0}
    frontier = [0]
    while frontier:
        point = frontier.pop()
        for image in images:
            next_point = image[point]
            if next_point not in orbit:
                orbit.add(next_point)
                frontier.append(next_point)
    assert orbit == set(range(6))


def main() -> None:
    audit_collision_partitions()
    checked = certificates()
    for certificate in checked:
        audit_group(certificate)
    audit_first_escape()
    print(
        "verified F2 k=1 complement monodromy: all five immersed collision "
        "partitions and the generic one-cusp and two-cusp strata have affine "
        "complement group Z; the E6+A1 escape has an exact transitive "
        "degree-six meridian action with fixed sheets"
    )


if __name__ == "__main__":
    main()
