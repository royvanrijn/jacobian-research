#!/usr/bin/env python3
"""Exact rational/group/count controls for the joint 2,3,17 elliptic sieve."""

from __future__ import annotations

import argparse
from fractions import Fraction as F
from hashlib import sha256
import itertools
import json
from pathlib import Path
import platform


A2, A4 = -10656, 27091730
P = (F(373248, 49), F(59459616, 343))
S = (F(92179201, 6400), F(553221611649, 512000))
T = (F(0), F(0))


def on_curve(point):
    if point is None:
        return True
    x, y = point
    return y*y == x**3+A2*x*x+A4*x


def add(left, right):
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u:
        if y+v == 0:
            return None
        slope = (3*x*x+2*A2*x+A4)/(2*y)
    else:
        slope = (v-y)/(u-x)
    new_x = slope*slope-A2-x-u
    return new_x, -y+slope*(x-new_x)


def multiply(point, integer):
    if integer < 0:
        return multiply((point[0], -point[1]), -integer)
    result = None
    while integer:
        if integer & 1:
            result = add(result, point)
        point = add(point, point)
        integer >>= 1
    return result


def valuation(value: F, prime: int):
    assert value != 0
    numerator, denominator = value.numerator, value.denominator
    result = 0
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def parameter(point):
    return -(point[0]-3552)/point[1]


def reduce_fraction(value, prime):
    return value.numerator*pow(value.denominator, -1, prime) % prime


def reduce_point(point, prime):
    return tuple(reduce_fraction(value, prime) for value in point)


def add_mod(left, right, prime):
    if left is None:
        return right
    if right is None:
        return left
    x, y = left
    u, v = right
    if x == u:
        if (y+v) % prime == 0:
            return None
        slope = (3*x*x+2*A2*x+A4)*pow(2*y, -1, prime) % prime
    else:
        slope = (v-y)*pow(u-x, -1, prime) % prime
    new_x = (slope*slope-A2-x-u) % prime
    return new_x, (-y+slope*(x-new_x)) % prime


def multiples_mod(point, prime):
    result, current = [None], point
    while current is not None:
        assert current not in result
        result.append(current)
        current = add_mod(current, point, prime)
    return result


def run():
    checks = []

    def require(name, condition):
        assert condition, name
        checks.append(name)

    require("three_exact_rational_points_on_E0", all(on_curve(point) for point in [P, S, T]))
    R = add(P, T)
    require("exact_R_equals_P_plus_T", R == (F(663747385, 186624), -F(6525490184045, 80621568)))
    require("exact_R_on_curve", on_curve(R))
    require("two_adic_formal_generators_have_depth_four",
            valuation(parameter(R), 2) == valuation(parameter(S), 2) == 4)
    require("normalized_two_adic_generator_coefficients_are_odd",
            reduce_fraction(parameter(R)/16, 2) == reduce_fraction(parameter(S)/16, 2) == 1)
    require("three_adic_unit_determinant_depth_witnesses",
            valuation(parameter(R), 3) == 3 and valuation(parameter(multiply(S, 2)), 3) == 2)
    require("formal_and_C9_matrix_determinant_mod_three", (0*1-1*1) % 3 != 0)
    prime_q = 2709173
    require("global_descent_class_of_R", R[0] == 5*prime_q*F(7, 432)**2)
    require("global_descent_class_of_S", S[0] == F(9601, 80)**2)
    require("global_descent_class_of_T", A4 == 10*prime_q)

    prime = 17
    cubic_discriminant = A4*A4*(A2*A2-4*A4)
    require("good_reduction_at_seventeen", cubic_discriminant % prime != 0)
    all_points = [None]+[(x, y) for x in range(prime) for y in range(prime)
                         if (y*y-x**3-A2*x*x-A4*x) % prime == 0]
    require("complete_finite_field_count_eighteen", len(all_points) == 18)
    pp, ss, tt, rr = [reduce_point(point, prime) for point in [P, S, T, R]]
    require("finite_field_reductions", (pp, ss, tt, rr) == ((2, 14), (2, 14), (0, 0), (10, 15)))
    g_multiples = multiples_mod(pp, prime)
    r_multiples = multiples_mod(rr, prime)
    require("P_and_S_order_nine", len(g_multiples) == 9 and pp == ss)
    require("T_order_two", len(multiples_mod(tt, prime)) == 2)
    require("R_generates_cyclic_group_eighteen", len(r_multiples) == 18 and set(r_multiples) == set(all_points))
    branch_points = {None, (16, 4), (16, 13)}
    require("finite_branch_values_are_exact_C9_subgroup", branch_points == {g_multiples[j] for j in [0, 3, 6]})
    require("finite_branch_points_are_all_F_zero_points",
            {point for point in all_points if point is not None and (128*point[0]-25865) % prime == 0}
            == branch_points-{None})
    require("finite_branch_points_have_unit_y", all(point[1] % prime for point in branch_points if point is not None))
    p9, r18 = multiply(P, 9), multiply(R, 18)
    require("nine_P_formal_seventeen_depth_one", valuation(parameter(p9), prime) == 1)
    require("nine_P_formal_seventeen_unit_ten", reduce_fraction(parameter(p9)/prime, prime) == 10)
    require("eighteen_R_equals_twice_nine_P", r18 == multiply(p9, 2))
    require("eighteen_R_formal_seventeen_depth_one", valuation(parameter(r18), prime) == 1)
    require("eighteen_R_formal_seventeen_unit_three", reduce_fraction(parameter(r18)/prime, prime) == 3)

    # Exhaust every two-primary parity, formal two-class, formal three-class,
    # and C9 residue in the chosen finite joint quotient. Each branch-centred
    # seventeen-adic class has one surviving formal class out of 17^2;
    # its actual translated centre need not be the rational identity.
    total = before17 = branch = allowed_weight = 0
    row_counts = {"00": [0, 0, 0], "10": [0, 0, 0],
                  "01": [0, 0, 0], "11": [0, 0, 0]}
    field_counts = {"1": 0, "5q": 0, "2": 0, "10q": 0}
    for c, d, u, j3, k in itertools.product(range(2), range(2), range(4), range(3), range(9)):
        total += 1
        soluble2_necessary = bool(c) or u == 0
        soluble3_necessary = ((c, d) == (1, 1)
                               or ((c, d) in [(0, 0), (0, 1)] and j3 == 0))
        if not (soluble2_necessary and soluble3_necessary):
            continue
        before17 += 1
        row = row_counts[f"{c}{d}"]
        row[0] += 1
        is_branch = (u+d+c) % 2 == 0 and k % 3 == 0
        if is_branch:
            branch += 1
            row[1] += 1
        weight = 1 if is_branch else 17**2
        allowed_weight += weight
        row[2] += weight
        m_parity = (u+d) % 2
        field_label = ({(0, 0): "1", (1, 0): "5q",
                        (1, 1): "2", (0, 1): "10q"})[(m_parity, c)]
        field_counts[field_label] += weight
    require("finite_joint_quotient_has_four_hundred_thirty_two_base_classes", total == 432)
    require("finite_two_three_retained_fraction", F(before17, total) == F(7, 24))
    require("finite_joint_branch_mass", F(branch, total) == F(7, 144))
    finite_ratio = F(allowed_weight, total*17**2)
    require("finite_three_prime_ratio", finite_ratio == F(1687, 6936))
    require("finite_ratio_strictly_below_quarter", F(1, 4)-finite_ratio == F(47, 6936))
    require("independent_finite_weight_formula", F(7, 24)-F(7, 144)*F(288, 289) == finite_ratio)
    require("finite_fixed_field_weight_partition", field_counts == {"1": 1737, "5q": 2601,
                                                                  "2": 10422, "10q": 15606})

    a2, aO, aB, aAB = F(1, 7), F(1, 13), F(1, 3), F(1)
    ratio23 = (a2*aO+a2*aB+aAB)/4
    require("exact_joint_two_three_baseline", ratio23 == F(289, 1092))
    exact_branch_mass = F(1, 24)+a2*aO/12
    require("exact_branch_mass", exact_branch_mass == F(31, 728))
    require("formal_seventeen_annulus_geometric_sum", F(1, 307) == F(16, 17)*F(17, 17**3)/(1-F(1, 17**3)))
    limit_ratio = ratio23-exact_branch_mass*F(306, 307)
    require("exact_local_only_three_prime_ratio", limit_ratio == F(5321, 23946))
    require("alternate_unreduced_ratio_agrees", limit_ratio == F(37247, 167622))
    before17_fields = {"1": a2*aO/4, "5q": a2*aB/4, "2": F(1, 8), "10q": F(1, 8)}
    require("fixed_field_two_three_baselines", before17_fields == {"1": F(1, 364), "5q": F(1, 84),
                                                                   "2": F(1, 8), "10q": F(1, 8)})
    seventeen_factor = 1-F(1, 3)*F(306, 307)
    require("fixed_field_seventeen_branch_factor", seventeen_factor == F(205, 307))
    field_coefficients = {label: value*(seventeen_factor if label in ["1", "2"] else 1)
                          for label, value in before17_fields.items()}
    require("four_exact_fixed_field_coefficients", field_coefficients == {
        "1": F(205, 111748), "5q": F(1, 84), "2": F(205, 2456), "10q": F(1, 8)})
    require("fixed_field_coefficients_sum_to_whole_slice", sum(field_coefficients.values()) == limit_ratio)

    here = Path(__file__).resolve().parent
    primary_path = here/"sources"/"milne-elliptic-curves-2006.pdf"
    primary_data = primary_path.read_bytes()
    require("retained_primary_formal_group_source_exists", primary_data.startswith(b"%PDF-"))
    return {
        "status": "pass", "python": platform.python_version(),
        "command": ".venv/bin/python research/arithmetic/verify_elliptic_exception_joint_2_3_17.py",
        "checker_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "proof_sha256": sha256((here/"ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md").read_bytes()).hexdigest(),
        "exact_controls": checks,
        "finite_field_17": {"point_count": len(all_points), "group": "C18", "P_order": 9,
                            "R_order": 18, "branch_points": [None, [16, 4], [16, 13]]},
        "formal_depths": {"R_at_2": 4, "S_at_2": 4, "R_at_3": 3, "two_S_at_3": 2,
                          "nine_P_at_17": 1, "eighteen_R_at_17": 1},
        "finite_weight_count": {"base_classes": total, "two_three_accepted": before17,
                                "branch_classes": branch, "seventeen_classes_per_base": 289,
                                "retained_classes": allowed_weight, "row_counts": row_counts,
                                "fixed_field_counts": field_counts,
                                "ratio": str(finite_ratio)},
        "local_only_ratio": str(limit_ratio),
        "fixed_field_coefficients": {label: str(value) for label, value in field_coefficients.items()},
        "external_primary_source": {"path": "sources/milne-elliptic-curves-2006.pdf",
                                    "sha256": sha256(primary_data).hexdigest(),
                                    "statement": "Chapter II Theorem 4.1 and Lemma 5.9: formal filtration and addition/depth control",
                                    "proof_imported_not_replayed": True},
        "scope": ("Exact finite-field enumeration, rational elliptic arithmetic, formal-depth "
                  "witnesses and an exhaustive finite parity/weight count. Global saturation, "
                  "the 2/3 local classifications, ELV1, formal groups and the height lattice "
                  "remain separate written/imported premises. No complete Mordell--Weil basis, "
                  "uniform annulus error or all-primes Hasse density is asserted."),
        "formal_verification": False, "independent_whole_proof_replay": False,
        "external_review": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-receipt", action="store_true")
    args = parser.parse_args()
    result = run()
    output = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if args.write_receipt:
        Path(__file__).resolve().with_name("elliptic_exception_joint_2_3_17_checks.json").write_text(output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
