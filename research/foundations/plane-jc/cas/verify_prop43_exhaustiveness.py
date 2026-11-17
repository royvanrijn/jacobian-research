#!/usr/bin/env python3
"""Exact finite audit of GGHV 2022, Proposition 4.3, for the (8,28) chain.

This checker deliberately does not re-prove the general Newton-polygon
results cited by Proposition 4.3.  It verifies every finite arithmetic and
lattice step that the proposition abbreviates by “as in Proposition 4.1” and
“the edges ... would have no way of being parallel”, together with the two
Laurent polygons produced by the final monomial map.

Dependencies: Python standard library only.
"""

from __future__ import annotations

from math import gcd
from typing import Iterable

Point = tuple[int, int]
Monomial = tuple[int, int, int]  # powers of x, y, alpha
Polynomial = dict[Monomial, int]


def weight(direction: Point, point: Point) -> int:
    return direction[0] * point[0] + direction[1] * point[1]


def determinant(left: Point, right: Point) -> int:
    return left[0] * right[1] - left[1] * right[0]


def subtract(left: Point, right: Point) -> Point:
    return left[0] - right[0], left[1] - right[1]


def scale(factor: int, point: Point) -> Point:
    return factor * point[0], factor * point[1]


def primitive_unsigned(vector: Point) -> Point:
    divisor = gcd(abs(vector[0]), abs(vector[1]))
    if divisor == 0:
        raise ValueError("zero vector has no primitive direction")
    return abs(vector[0]) // divisor, abs(vector[1]) // divisor


def prime_defect_data(current: Point, candidate: Point) -> tuple[int, Point, int]:
    """Return scale, primitive incoming step, and scale-free defect.

    If current=(A,B), candidate=(a,b), and
        current-candidate = g*(p,q),
    then
        det(current,candidate) = g*(B*p-A*q).
    The auxiliary-F condition therefore descends from
        det(current,candidate) | (A-B)*g
    to
        B*p-A*q | A-B.
    """
    difference = subtract(current, candidate)
    divisor = gcd(abs(difference[0]), abs(difference[1]))
    if divisor == 0:
        raise ValueError("current and candidate must differ")
    step = difference[0] // divisor, difference[1] // divisor
    defect = current[1] * step[0] - current[0] * step[1]
    assert determinant(current, candidate) == divisor * defect
    return divisor, step, defect


def poly_add(*polynomials: Polynomial) -> Polynomial:
    out: Polynomial = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            out[monomial] = out.get(monomial, 0) + coefficient
            if out[monomial] == 0:
                del out[monomial]
    return out


def poly_scale(polynomial: Polynomial, coefficient: int) -> Polynomial:
    return {
        monomial: coefficient * value
        for monomial, value in polynomial.items()
        if coefficient * value
    }


def poly_multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    out: Polynomial = {}
    for (lx, ly, la), left_coefficient in left.items():
        for (rx, ry, ra), right_coefficient in right.items():
            monomial = lx + rx, ly + ry, la + ra
            out[monomial] = (
                out.get(monomial, 0) + left_coefficient * right_coefficient
            )
            if out[monomial] == 0:
                del out[monomial]
    return out


def poly_power(polynomial: Polynomial, exponent: int) -> Polynomial:
    if exponent < 0:
        raise ValueError("negative polynomial powers are unsupported")
    out: Polynomial = {(0, 0, 0): 1}
    base = polynomial
    while exponent:
        if exponent & 1:
            out = poly_multiply(out, base)
        base = poly_multiply(base, base)
        exponent >>= 1
    return out


def monomial(
    x: int = 0,
    y: int = 0,
    alpha: int = 0,
    coefficient: int = 1,
) -> Polynomial:
    return {(x, y, alpha): coefficient} if coefficient else {}


def cross(origin: Point, first: Point, second: Point) -> int:
    return determinant(subtract(first, origin), subtract(second, origin))


def convex_hull(points: Iterable[Point]) -> tuple[Point, ...]:
    """Strict monotone-chain hull; collinear interior points are removed."""
    ordered = sorted(set(points))
    if len(ordered) <= 1:
        return tuple(ordered)

    lower: list[Point] = []
    for point in ordered:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)

    upper: list[Point] = []
    for point in reversed(ordered):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)

    return tuple(lower[:-1] + upper[:-1])


def admissible_k(point: Point) -> tuple[int, ...]:
    a, b = point
    if b <= 0:
        raise ValueError("Proposition 8.2 requires a positive second coordinate")
    return tuple(k for k in range(1, a // b + 1) if (k + 1) * b < a)


def endpoint_determinants(point: Point, k: int) -> tuple[int, int]:
    """Parallelism determinants for the two Proposition-8.2 assignments."""
    a, b = point

    # Assignment 1: P -> (-k,0), Q -> (k+1,1).
    direct = determinant(
        subtract((-k, 0), scale(2, point)),
        subtract((k + 1, 1), scale(3, point)),
    )

    # Assignment 2: P -> (k+1,1), Q -> (-k,0).
    swapped = determinant(
        subtract((k + 1, 1), scale(2, point)),
        subtract((-k, 0), scale(3, point)),
    )

    # Closed forms used in the accompanying proof note.
    assert direct == -2 * a + (5 * k + 2) * b - k
    assert swapped == 3 * a - (5 * k + 3) * b + k
    if 2 * a == 7 * b - 1:
        assert direct == (k - 1) * (5 * b - 1)
        assert 2 * swapped == (3 - 2 * k) * (5 * b - 1)
    return direct, swapped


def final_map(point: Point) -> Point:
    """Exponent map induced by x -> x^-1, y -> x^4 y."""
    i, j = point
    return -i + 4 * j, j


def main() -> None:
    # ------------------------------------------------------------------
    # A. The finite Proposition-2.5 check inside the predecessor step.
    # ------------------------------------------------------------------
    a, ell = 7, 1
    deltas = [
        delta
        for delta in range(ell + 1, (a - 1) // 2 + 1)
        if (delta - ell) % (a - 2 * delta) == 0
    ]
    assert deltas == [3]
    initial_directions = {(1, -2), (ell, -deltas[0])}
    assert initial_directions == {(1, -2), (1, -3)}
    print("PASS predecessor directions: (1,-2), (1,-3)")

    # ------------------------------------------------------------------
    # B. The upper-edge shear is an exact two-term identity.
    # ------------------------------------------------------------------
    x4 = monomial(x=4)
    y_var = monomial(y=1)
    alpha_var = monomial(alpha=1)
    shifted_y = poly_add(y_var, monomial(x=-4, alpha=1))
    shifted_inner = poly_add(
        poly_multiply(x4, shifted_y),
        poly_scale(alpha_var, -1),
    )
    assert shifted_inner == monomial(x=4, y=1)
    shifted_edge = poly_multiply(shifted_y, poly_power(shifted_inner, 7))
    assert shifted_edge == {
        (28, 8, 0): 1,
        (24, 7, 1): 1,
    }
    print("PASS edge shear: y(x^4 y-a)^7 -> x^28 y^8 + a x^24 y^7")

    # ------------------------------------------------------------------
    # C. Verify the uniform scale-free prime-defect lemma first on the
    #    printed Proposition-4.1 corner and then on the live corner.
    # ------------------------------------------------------------------
    printed_current = (21, 8)
    printed_survivors = [(13, 5), (5, 2), (1, 1)]
    printed_data = [prime_defect_data(printed_current, point) for point in printed_survivors]
    assert printed_data == [
        (1, (8, 3), 1),
        (2, (8, 3), 1),
        (1, (20, 7), 13),
    ]
    assert all(13 % defect == 0 for _, _, defect in printed_data)

    live_current = (24, 7)
    live_survivors = [(17, 5), (10, 3), (3, 1), (1, 1)]
    live_data = [prime_defect_data(live_current, point) for point in live_survivors]
    assert live_data == [
        (1, (7, 2), 1),
        (2, (7, 2), 1),
        (3, (7, 2), 1),
        (1, (23, 6), 17),
    ]
    assert all(17 % defect == 0 for _, _, defect in live_data)
    print("PASS prime-defect lemma: Proposition 4.1 and live (24,7) corner")

    # ------------------------------------------------------------------
    # D. Reconstruct the omitted 'as in Proposition 4.1' lattice census.
    # ------------------------------------------------------------------
    current = live_current
    old_direction = (-1, 4)
    assert weight(old_direction, current) == 4

    geometric_candidates = [
        (x, y)
        for x in range(-3, current[0] + 1)
        for y in range(0, current[1] + 1)
        if weight(old_direction, (x, y)) < weight(old_direction, current)
        and determinant(current, (x, y)) > 0
    ]
    expected_geometric = [
        (-3, 0),
        (-2, 0),
        (-1, 0),
        (1, 1),
        (2, 1),
        (3, 1),
        (5, 2),
        (6, 2),
        (9, 3),
        (10, 3),
        (13, 4),
        (17, 5),
    ]
    assert geometric_candidates == expected_geometric
    print(f"PASS geometric census: {len(geometric_candidates)} candidates")

    # If F has its other corner on the ray Q_{>0}(24,7), the same argument
    # printed in Proposition 4.1 gives
    #       24*y - 7*x  |  17*gcd(24-x, 7-y).
    f_compatible: list[tuple[Point, int]] = []
    print("\nF-compatibility table")
    print(" point       delta  gcd  17*gcd  c")
    for point in geometric_candidates:
        x, y = point
        step_gcd = gcd(current[0] - x, current[1] - y)
        transverse_distance = determinant(current, point)
        compatible = 17 * step_gcd % transverse_distance == 0
        c = 17 * step_gcd // transverse_distance if compatible else None
        if c is not None:
            primitive_step = (
                (current[0] - x) // step_gcd,
                (current[1] - y) // step_gcd,
            )
            f_corner = (1 + c * primitive_step[0], 1 + c * primitive_step[1])
            assert determinant(current, f_corner) == 0
            assert f_corner[0] > 0 and f_corner[1] > 0
        print(
            f" {str(point):10s} {transverse_distance:5d}"
            f" {step_gcd:4d} {17 * step_gcd:7d}"
            f" {'-' if c is None else c}"
        )
        if compatible:
            f_compatible.append((point, c))

    assert f_compatible == [
        ((1, 1), 1),
        ((3, 1), 17),
        ((10, 3), 17),
        ((17, 5), 17),
    ]
    proper_proportional = [
        point for point, _ in f_compatible if point != (1, 1)
    ]
    assert proper_proportional == [(3, 1), (10, 3), (17, 5)]
    print("PASS F census: only (3,1), (10,3), (17,5) survive off diagonal")

    # The immediate nonproportional branch is represented by (24,7).
    recursive_candidates = [current, (17, 5), (10, 3), (3, 1)]
    assert all(2 * x == 7 * y - 1 for x, y in recursive_candidates)
    assert all(
        primitive_unsigned(subtract(current, point)) == (7, 2)
        for point in proper_proportional
    )
    successor_normal = (-2, 7)
    assert all(weight(successor_normal, point) == 1 for point in recursive_candidates)
    print("PASS common line: 2a=7b-1, primitive edge (7,2), normal (-2,7)")

    # ------------------------------------------------------------------
    # E. Replace 'no way of being parallel' by exact determinants.
    # ------------------------------------------------------------------
    expected_k = {
        (24, 7): (1, 2),
        (17, 5): (1, 2),
        (10, 3): (1, 2),
        (3, 1): (1,),
    }
    assert {point: admissible_k(point) for point in recursive_candidates} == expected_k

    parallelism_table: dict[tuple[Point, int], tuple[int, int]] = {}
    print("\nParallelism determinant table")
    print(" (a,b)      k   P->(-k,0)   swapped")
    for point in recursive_candidates:
        for k in admissible_k(point):
            values = endpoint_determinants(point, k)
            parallelism_table[(point, k)] = values
            print(f" {str(point):10s} {k:2d} {values[0]:11d} {values[1]:10d}")

    assert parallelism_table == {
        ((24, 7), 1): (0, 17),
        ((24, 7), 2): (34, -17),
        ((17, 5), 1): (0, 12),
        ((17, 5), 2): (24, -12),
        ((10, 3), 1): (0, 7),
        ((10, 3), 2): (14, -7),
        ((3, 1), 1): (0, 2),
    }

    for point in recursive_candidates:
        # k=1 has exactly one parallel assignment:
        #       en(P)=(-1,0), en(Q)=(2,1).
        direct, swapped = endpoint_determinants(point, 1)
        assert direct == 0 and swapped != 0
        p_edge = subtract((-1, 0), scale(2, point))
        q_edge = subtract((2, 1), scale(3, point))
        assert determinant(p_edge, q_edge) == 0
        assert primitive_unsigned(p_edge) == primitive_unsigned(q_edge) == (7, 2)

    # A proper proportional candidate would be followed by an edge on the
    # same ray on which it was entered, hence it is not a distinct vertex.
    for point in proper_proportional:
        previous_p_edge = subtract(scale(2, point), scale(2, current))
        next_p_edge = subtract((-1, 0), scale(2, point))
        previous_q_edge = subtract(scale(3, point), scale(3, current))
        next_q_edge = subtract((2, 1), scale(3, point))
        assert determinant(previous_p_edge, next_p_edge) == 0
        assert determinant(previous_q_edge, next_q_edge) == 0
        assert previous_p_edge[0] < 0 and next_p_edge[0] < 0
        assert previous_q_edge[0] < 0 and next_q_edge[0] < 0

    print("PASS k=2 excluded exactly; k=1 has the unique endpoint assignment")
    print("PASS proper proportional candidates collapse as collinear nonvertices")

    # ------------------------------------------------------------------
    # F. Verify the final Laurent map and both output polygons.
    # ------------------------------------------------------------------
    case2_p_before = {(-1, 0), (0, 0), (56, 16), (48, 14)}
    case2_q_before = {(2, 1), (0, 0), (84, 24), (72, 21)}
    case1_p_before = case2_p_before | {(32, 8)}
    case1_q_before = case2_q_before | {(48, 12)}

    case2_p_after = {final_map(point) for point in case2_p_before}
    case2_q_after = {final_map(point) for point in case2_q_before}
    case1_p_after = {final_map(point) for point in case1_p_before}
    case1_q_after = {final_map(point) for point in case1_q_before}

    assert case2_p_after == {(0, 0), (1, 0), (8, 14), (8, 16)}
    assert case2_q_after == {(0, 0), (2, 1), (12, 21), (12, 24)}
    assert case1_p_after == case2_p_after | {(0, 8)}
    assert case1_q_after == case2_q_after | {(0, 12)}

    # The PDF proposition statement prints (0,1) in Case 1.  It is collinear
    # on the x=0 boundary and is not a vertex.  The proof's final computation
    # and the exponent map above force (8,14), as used by the repository.
    printed_case1_p = {(0, 0), (1, 0), (0, 1), (8, 16), (0, 8)}
    assert (0, 1) not in convex_hull(printed_case1_p)
    assert set(convex_hull(case1_p_after)) == case1_p_after
    assert (8, 14) in case1_p_after

    # For X=x^a y^b and Y=x^c y^d,
    # [X,Y]=(ad-bc)x^(a+c-1)y^(b+d-1).
    map_a, map_b, map_c, map_d = -1, 0, 4, 1
    jacobian_coefficient = map_a * map_d - map_b * map_c
    jacobian_exponent = (
        map_a + map_c - 1,
        map_b + map_d - 1,
    )
    assert jacobian_coefficient == -1
    assert jacobian_exponent == (2, 0)

    print("PASS final map: exactly the two corrected Proposition-4.3 polygons")
    print("PASS chain rule: coordinate Jacobian is -x^2")
    print("PASS printed typo: Case-1 (0,1) must read (8,14)")
    print("PROP43_EXHAUSTIVENESS_AUDIT_PASS")


if __name__ == "__main__":
    main()
