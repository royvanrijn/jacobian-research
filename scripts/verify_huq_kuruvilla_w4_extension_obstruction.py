#!/usr/bin/env python3
"""Exact W_4 extension obstructions for two stabilized W_3 lifts.

The preferred degree-19 W_3 certificate has a nonzero top de Rham class in
its next determinant digit, so that exact representative has no W_4 lift in
any polynomial degree.  The canonical degree-25 W_3 representative has zero
top class.  Exact GF(2) elimination proves that its complete necessary z^0
correction equation has no solution in degree at most 51, while the explicit
third-coordinate correction in degree 52 gives determinant one modulo 16.

These are fixed-representative results.  Together with d_3=19 they imply only
19 <= d_4 <= 52 for the unrestricted stabilized lifting problem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from affine_support_decoder import (
    affine_inconsistency_certificate,
    solve_affine_rows,
)


Monomial2 = tuple[int, int]
Monomial3 = tuple[int, int, int]
Polynomial3 = dict[Monomial3, int]

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument(
    "--search-degree19-extension",
    action="store_true",
    help="find the exact W_4 extension degree of the pinned class-zero W_3 lift",
)
parser.add_argument(
    "--write-extension-certificate",
    type=Path,
    help="write the returned W_4 correction support under generated-results",
)
parser.add_argument(
    "--replay-extension-certificate",
    type=Path,
    help="replay a pinned degree-52 W_4 correction without elimination",
)
args = parser.parse_args()
assert args.write_extension_certificate is None or args.search_degree19_extension

PREFERRED_CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/"
    "huq_kuruvilla_w3_degree19_witness_first_280_second_160.json"
)
PREFERRED_SHA256 = (
    "a79984550854ce01d903783156baa6f7d4720f56ec2824bd5823cd13088a5d7f"
)
CLASS_ZERO_CERTIFICATE = (
    ROOT
    / "artifacts/generated-results/"
    "huq_kuruvilla_w3_degree19_w4_class_zero.json"
)
CLASS_ZERO_SHA256 = (
    "9ad15068593af7cca87169c25eed2ff53068cc466d183ce320c2fc7d0e2c1aaa"
)
CLASS_ZERO_ERROR_SHA256 = (
    "846526f5c8b999ca35dcc838054eb67a6215bd8501429964f0056afb0b42fbdd"
)
FROZEN_EXTENSION_SHA256 = {
    "huq_kuruvilla_w4_degree52_from_degree19.json": (
        "438a189da33fdd081f61f9410186ca7d1b22c454dfb2cae5fbc02060f1b838ae"
    ),
}
AFFINE_DECODER_SHA256 = (
    "cd829ec30aca757c052fadf936c5cfb55abbf6622e09917e7bae841f263e834f"
)
assert hashlib.sha256(
    Path(__file__).with_name("affine_support_decoder.py").read_bytes()
).hexdigest() == AFFINE_DECODER_SHA256


def add_term(
    polynomial: Polynomial3,
    monomial: Monomial3,
    coefficient: int,
    modulus: int | None = None,
) -> None:
    value = polynomial.get(monomial, 0) + coefficient
    if modulus is not None:
        value %= modulus
    if value:
        polynomial[monomial] = value
    else:
        polynomial.pop(monomial, None)


def derivative(
    polynomial: Polynomial3,
    axis: int,
    modulus: int | None = None,
) -> Polynomial3:
    output: Polynomial3 = {}
    for monomial, coefficient in polynomial.items():
        exponent = monomial[axis]
        if not exponent:
            continue
        result = list(monomial)
        result[axis] -= 1
        add_term(output, tuple(result), coefficient * exponent, modulus)
    return output


def product(
    left: Polynomial3,
    right: Polynomial3,
    modulus: int | None = None,
) -> Polynomial3:
    output: Polynomial3 = {}
    for (i, j, k), left_coefficient in left.items():
        for (u, v, w), right_coefficient in right.items():
            add_term(
                output,
                (i + u, j + v, k + w),
                left_coefficient * right_coefficient,
                modulus,
            )
    return output


def determinant(coordinates: list[Polynomial3], modulus: int) -> Polynomial3:
    jacobian = [
        [derivative(coordinate, axis, modulus) for axis in range(3)]
        for coordinate in coordinates
    ]
    output: Polynomial3 = {}
    for sign, indices in (
        (1, (0, 1, 2)),
        (1, (1, 2, 0)),
        (1, (2, 0, 1)),
        (-1, (2, 1, 0)),
        (-1, (1, 0, 2)),
        (-1, (0, 2, 1)),
    ):
        term = product(
            product(jacobian[0][indices[0]], jacobian[1][indices[1]], modulus),
            jacobian[2][indices[2]],
            modulus,
        )
        for monomial, coefficient in term.items():
            add_term(output, monomial, sign * coefficient, modulus)
    return output


def determinant_error_digit(
    coordinates: list[Polynomial3], constant_mod_8: int = 1
) -> Polynomial3:
    """Return (det-c)/8 modulo two, asserting det=c modulo eight."""

    jacobian = determinant(coordinates, 16)
    error: Polynomial3 = {}
    for monomial in set(jacobian) | {(0, 0, 0)}:
        residue = (
            jacobian.get(monomial, 0)
            - constant_mod_8 * int(monomial == (0, 0, 0))
        ) % 16
        assert residue in (0, 8), (monomial, residue)
        if residue == 8:
            error[monomial] = 1
    return error


def maximum_degree(polynomial: Polynomial3) -> int:
    return max(map(sum, polynomial), default=-1)


P: Polynomial3 = {
    (1, 0, 0): 1,
    (2, 1, 0): 1,
    (4, 0, 0): 1,
    (6, 2, 0): 1,
}
Q: Polynomial3 = {
    (0, 1, 0): 1,
    (5, 0, 0): 1,
    (6, 1, 0): 1,
    (7, 2, 0): 1,
    (8, 3, 0): 1,
}
Z: Polynomial3 = {(0, 0, 1): 1}

SUPPORT_LAYOUT = {
    "p": (2, 0, 0),
    "q": (2, 1, 0),
    "c": (2, 2, 0),
    "a": (2, 0, 1),
    "b": (2, 1, 1),
    "r": (2, 2, 1),
    "u0": (4, 0, 0),
    "v0": (4, 1, 0),
    "w1": (4, 2, 1),
    "u1": (4, 0, 1),
    "v1": (4, 1, 1),
    "u2": (4, 0, 2),
    "v2": (4, 1, 2),
    "w3": (4, 2, 3),
}


def degree_19_lift_from_certificate(
    path: Path, expected_sha256: str
) -> tuple[list[Polynomial3], dict]:
    raw = path.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == expected_sha256
    data = json.loads(raw)
    assert data["claim"] == "HKM2W3"
    assert data["degree"] == 19
    assert data["modulus"] == 8
    assert data["total_support"] == sum(data["support_counts"].values())
    coordinates = [dict(P), dict(Q), dict(Z)]
    assert set(data["supports"]) == set(SUPPORT_LAYOUT)
    for name, support in data["supports"].items():
        scale, coordinate, z_degree = SUPPORT_LAYOUT[name]
        for i, j in support:
            add_term(coordinates[coordinate], (i, j, z_degree), scale, 16)
    assert max(maximum_degree(coordinate) for coordinate in coordinates) == 19
    return coordinates, data


preferred_degree_19, preferred_data = degree_19_lift_from_certificate(
    PREFERRED_CERTIFICATE, PREFERRED_SHA256
)
assert preferred_data["jacobian_constant_mod_8"] == 1
assert preferred_data["total_support"] == 440
preferred_error = determinant_error_digit(preferred_degree_19)
preferred_top_class = sorted(
    monomial
    for monomial in preferred_error
    if all(exponent % 2 for exponent in monomial)
)
preferred_class_encoding = json.dumps(
    preferred_top_class, separators=(",", ":")
).encode()
PREFERRED_TOP_CLASS_SHA256 = (
    "e219eeab3de8badeaf76c9cb393a1b1f0d8a791ae794ef66e97fa4b70c77b9fb"
)
assert len(preferred_error) == 1250
assert maximum_degree(preferred_error) == 54
assert len(preferred_top_class) == 48
assert maximum_degree({monomial: 1 for monomial in preferred_top_class}) == 35
assert {monomial[2] for monomial in preferred_top_class} == {1}
assert (1, 1, 1) in preferred_top_class
assert (
    hashlib.sha256(preferred_class_encoding).hexdigest()
    == PREFERRED_TOP_CLASS_SHA256
)
print("PASS: preferred degree-19 W_3 determinant error has support 1250")
print("PASS: its H^3_dR class has 48 odd-odd-odd monomials and is nonzero")
print(
    "AUDIT: preferred H^3 class SHA256 "
    f"{hashlib.sha256(preferred_class_encoding).hexdigest()}"
)
print("PASS: the preferred degree-19 W_3 representative has no W_4 extension")


# Integral half-Jacobian error K=(det D(P,Q)-1)/2.
K_INTEGER: dict[Monomial2, int] = {
    (1, 1): 1,
    (3, 0): 2,
    (5, 2): 3,
    (6, 0): -2,
    (7, 1): -1,
    (9, 0): 2,
    (9, 3): -1,
    (10, 1): -1,
    (11, 2): 3,
    (12, 3): -1,
    (13, 4): 1,
}
K_SUPPORT = {monomial for monomial, coefficient in K_INTEGER.items() if coefficient % 2}
DELTA = {(2, 2), (10, 4)}
C3_SUPPORT = K_SUPPORT.symmetric_difference(DELTA)
A3_SUPPORT = {(15, 2), (19, 6)}


def canonical_degree_25_lift() -> list[Polynomial3]:
    coordinates = [dict(P), dict(Q), dict(Z)]
    for i, j in A3_SUPPORT:
        add_term(coordinates[0], (i, j, 0), 4, 16)
    for (i, j), coefficient in K_INTEGER.items():
        add_term(coordinates[2], (i, j, 1), 2 * coefficient, 16)
    for i, j in C3_SUPPORT:
        add_term(coordinates[2], (i, j, 1), 4, 16)
    assert max(maximum_degree(coordinate) for coordinate in coordinates) == 25
    return coordinates


canonical_w3 = canonical_degree_25_lift()
canonical_error3 = determinant_error_digit(canonical_w3)
assert len(canonical_error3) == 35
assert maximum_degree(canonical_error3) == 51
assert {monomial[2] for monomial in canonical_error3} == {0}
canonical_top_class = {
    monomial
    for monomial in canonical_error3
    if all(exponent % 2 for exponent in monomial)
}
assert not canonical_top_class
L: dict[Monomial2, int] = {
    (i, j): coefficient for (i, j, _), coefficient in canonical_error3.items()
}
L_encoding = json.dumps(sorted(L), separators=(",", ":")).encode()
CANONICAL_ERROR_SHA256 = (
    "c9e4b2139db532cac4af47e976f0b2373da24996c66a3557ad00d77f61e6655d"
)
assert hashlib.sha256(L_encoding).hexdigest() == CANONICAL_ERROR_SHA256
print("PASS: canonical degree-25 W_3 error has support 35 and degree 51")
print("PASS: its H^3_dR obstruction vanishes")
print(
    "AUDIT: canonical next-error SHA256 "
    f"{hashlib.sha256(L_encoding).hexdigest()}"
)


PX = ((0, 0),)
PY = ((2, 0),)
QX = ((4, 0), (6, 2))
QY = ((0, 0), (6, 0), (8, 2))


def monomials(maximum: int) -> list[Monomial2]:
    return [
        (i, total - i)
        for total in range(maximum + 1)
        for i in range(total + 1)
    ]


def z0_completion_rows(
    degree: int,
) -> tuple[list[tuple[int, int]], int, list[Monomial2]]:
    """Encode D_F(R,S)+t=L for degree-d W_4 corrections."""

    variable_keys = (
        [("R", monomial) for monomial in monomials(degree)]
        + [("S", monomial) for monomial in monomials(degree)]
        + [("t", monomial) for monomial in monomials(degree - 1)]
    )
    contributions: dict[Monomial2, int] = {}

    def contribute(target: Monomial2, variable_index: int) -> None:
        contributions[target] = contributions.get(target, 0) ^ (
            1 << variable_index
        )

    for variable_index, (name, (i, j)) in enumerate(variable_keys):
        if name == "R":
            if i % 2:
                for u, v in QY:
                    contribute((i - 1 + u, j + v), variable_index)
            if j % 2:
                for u, v in QX:
                    contribute((i + u, j - 1 + v), variable_index)
        elif name == "S":
            if j % 2:
                for u, v in PX:
                    contribute((i + u, j - 1 + v), variable_index)
            if i % 2:
                for u, v in PY:
                    contribute((i - 1 + u, j + v), variable_index)
        else:
            contribute((i, j), variable_index)

    targets = sorted(set(contributions) | set(L))
    rows = [
        (contributions.get(monomial, 0), int(monomial in L))
        for monomial in targets
    ]
    return rows, len(variable_keys), targets


rows_51, variables_51, targets_51 = z0_completion_rows(51)
assert solve_affine_rows(rows_51) is None
dual_indices = affine_inconsistency_certificate(rows_51)
assert dual_indices is not None
dual_targets = sorted(targets_51[index] for index in dual_indices)
dual_encoding = json.dumps(dual_targets, separators=(",", ":")).encode()
DEGREE_51_DUAL_SHA256 = (
    "6b5bb61d2fc7e79fce3d1623eef931c4a398dfa70376caf67afdc9ab50acd280"
)
assert hashlib.sha256(dual_encoding).hexdigest() == DEGREE_51_DUAL_SHA256
assert dual_targets == [(39, 12), (40, 13)]
# The eliminated dual certificate is the transparent all-degree identity
#   Lambda(D_F(R,S)) = [x^41 y^13]R,
# where Lambda extracts the sum of the two displayed target coefficients.
# Its right side vanishes for deg R<=51, as does Lambda(t) for deg t<=50,
# while Lambda(L)=1.  Thus the exclusion is independent of the elimination
# ordering and also allows an arbitrary constant determinant.
functional_profile: dict[str, set[Monomial2]] = {"R": set(), "S": set()}


def toggle_profile(name: str, monomial: Monomial2) -> None:
    functional_profile[name].symmetric_difference_update({monomial})


for target_i, target_j in dual_targets:
    for u, v in QY:
        source = (target_i - u + 1, target_j - v)
        if min(source) >= 0 and source[0] % 2:
            toggle_profile("R", source)
    for u, v in QX:
        source = (target_i - u, target_j - v + 1)
        if min(source) >= 0 and source[1] % 2:
            toggle_profile("R", source)
    for u, v in PX:
        source = (target_i - u, target_j - v + 1)
        if min(source) >= 0 and source[1] % 2:
            toggle_profile("S", source)
    for u, v in PY:
        source = (target_i - u + 1, target_j - v)
        if min(source) >= 0 and source[0] % 2:
            toggle_profile("S", source)
assert functional_profile == {"R": {(41, 13)}, "S": set()}
assert sum(L.get(target, 0) for target in dual_targets) % 2 == 1
dual_left = 0
dual_right = 0
for index in dual_indices:
    left, right = rows_51[index]
    dual_left ^= left
    dual_right ^= right
assert dual_left == 0 and dual_right == 1
rows_51_any_constant = [
    row
    for target, row in zip(targets_51, rows_51, strict=True)
    if target != (0, 0)
]
assert solve_affine_rows(rows_51_any_constant) is None
print(
    "PASS: exact elimination excludes every degree-at-most-51 extension "
    f"({variables_51} variables, {len(rows_51)} coefficient equations)"
)
print(
    f"AUDIT: degree-51 dual certificate has {len(dual_targets)} targets and "
    f"SHA256 {hashlib.sha256(dual_encoding).hexdigest()}: {dual_targets}"
)
print(
    "PASS: Lambda(D_F(R,S))=[x^41*y^13]R gives a two-coefficient "
    "degree-51 obstruction"
)

rows_52, variables_52, _ = z0_completion_rows(52)
assert solve_affine_rows(rows_52) is not None
print(
    "PASS: the degree-52 affine completion system is consistent "
    f"({variables_52} variables, {len(rows_52)} coefficient equations)"
)

# The transparent solution is T=zL.  Adding 8T changes the next determinant
# digit by L because det D(P,Q)=1 modulo two.
canonical_w4 = [dict(coordinate) for coordinate in canonical_w3]
for (i, j), coefficient in L.items():
    add_term(canonical_w4[2], (i, j, 1), 8 * coefficient, 16)
assert max(maximum_degree(coordinate) for coordinate in canonical_w4) == 52
assert determinant(canonical_w4, 16) == {(0, 0, 0): 1}
print("PASS: adding 8*z*L gives a degree-52 determinant-one W_4 lift")
print("PASS: 52 is exact for extensions of the fixed degree-25 W_3 lift")
print("PASS: unrestricted stabilized degree satisfies 19 <= d_4 <= 52")


def full_w4_completion_rows(
    error: Polynomial3,
    degree: int,
) -> tuple[list[tuple[int, int]], list[tuple[str, Monomial3]], list[Monomial3]]:
    """Encode the full first-variation equation for a degree-d correction."""

    contributions: dict[Monomial3, int] = {}
    variable_keys: list[tuple[str, Monomial3]] = []

    def add_variable(
        name: str,
        source: Monomial3,
        targets: set[Monomial3],
    ) -> None:
        if not targets:
            return
        index = len(variable_keys)
        variable_keys.append((name, source))
        for target in targets:
            contributions[target] = contributions.get(target, 0) ^ (1 << index)

    z_layers = sorted({monomial[2] for monomial in error})
    for z_degree in z_layers:
        plane_bound = degree - z_degree
        if plane_bound >= 0:
            for i, j in monomials(plane_bound):
                r_targets: set[Monomial3] = set()
                if i % 2:
                    r_targets.update(
                        (i - 1 + u, j + v, z_degree) for u, v in QY
                    )
                if j % 2:
                    r_targets.update(
                        (i + u, j - 1 + v, z_degree) for u, v in QX
                    )
                add_variable("R", (i, j, z_degree), r_targets)

                s_targets: set[Monomial3] = set()
                if j % 2:
                    s_targets.update(
                        (i + u, j - 1 + v, z_degree) for u, v in PX
                    )
                if i % 2:
                    s_targets.update(
                        (i - 1 + u, j + v, z_degree) for u, v in PY
                    )
                add_variable("S", (i, j, z_degree), s_targets)

        # d/dz(z^(z_degree+1)) is nonzero in characteristic two exactly
        # when the output z-degree is even.
        t_plane_bound = degree - z_degree - 1
        if z_degree % 2 == 0 and t_plane_bound >= 0:
            for i, j in monomials(t_plane_bound):
                add_variable(
                    "T",
                    (i, j, z_degree + 1),
                    {(i, j, z_degree)},
                )

    targets = sorted(set(contributions) | set(error))
    rows = [
        (contributions.get(target, 0), int(target in error))
        for target in targets
        if target != (0, 0, 0)
    ]
    nonconstant_targets = [
        target for target in targets if target != (0, 0, 0)
    ]
    return rows, variable_keys, nonconstant_targets


if args.search_degree19_extension:
    class_zero_w3, class_zero_data = degree_19_lift_from_certificate(
        CLASS_ZERO_CERTIFICATE,
        CLASS_ZERO_SHA256,
    )
    class_zero_constant = class_zero_data["jacobian_constant_mod_8"]
    assert class_zero_constant == 5
    class_zero_error = determinant_error_digit(
        class_zero_w3,
        class_zero_constant,
    )
    assert not {
        monomial
        for monomial in class_zero_error
        if all(exponent % 2 for exponent in monomial)
    }
    class_zero_error_encoding = json.dumps(
        sorted(class_zero_error), separators=(",", ":")
    ).encode()
    assert (
        hashlib.sha256(class_zero_error_encoding).hexdigest()
        == CLASS_ZERO_ERROR_SHA256
    )
    print(
        "PASS: pinned degree-19 class-zero witness has next-error support "
        f"{len(class_zero_error)} and degree {maximum_degree(class_zero_error)}"
    )

    cache: dict[int, tuple[set[int] | None, list[tuple[str, Monomial3]], int]] = {}

    def solve_degree(
        degree: int,
    ) -> tuple[set[int] | None, list[tuple[str, Monomial3]], int]:
        if degree not in cache:
            rows, variable_keys, _ = full_w4_completion_rows(
                class_zero_error,
                degree,
            )
            cache[degree] = (
                solve_affine_rows(rows),
                variable_keys,
                len(rows),
            )
        return cache[degree]

    lower = 18
    upper = max(19, maximum_degree(class_zero_error) + 10)
    while solve_degree(upper)[0] is None:
        lower = upper
        upper *= 2
        assert upper <= 160, "unexpectedly large W_4 primitive degree"
    while lower + 1 < upper:
        middle = (lower + upper) // 2
        solution, variable_keys, equation_count = solve_degree(middle)
        status = "SAT" if solution is not None else "UNSAT"
        print(
            f"SEARCH: fixed class-zero W_3 extension degree {middle}: "
            f"{status} ({len(variable_keys)} variables, "
            f"{equation_count} equations)"
        )
        if solution is None:
            lower = middle
        else:
            upper = middle

    solution, variable_keys, equation_count = solve_degree(upper)
    assert solution is not None
    assert solve_degree(upper - 1)[0] is None
    correction_supports: dict[str, list[Monomial3]] = {
        "R": [],
        "S": [],
        "T": [],
    }
    class_zero_w4 = [dict(coordinate) for coordinate in class_zero_w3]
    coordinate_for_name = {"R": 0, "S": 1, "T": 2}
    for index in sorted(solution):
        name, monomial = variable_keys[index]
        correction_supports[name].append(monomial)
        add_term(
            class_zero_w4[coordinate_for_name[name]],
            monomial,
            8,
            16,
        )
    actual_degree = max(
        maximum_degree(coordinate) for coordinate in class_zero_w4
    )
    assert actual_degree == upper
    class_zero_w4_determinant = determinant(class_zero_w4, 16)
    assert set(class_zero_w4_determinant) == {(0, 0, 0)}
    determinant_constant = class_zero_w4_determinant[(0, 0, 0)]
    assert determinant_constant % 2 == 1
    support_counts = {
        name: len(support) for name, support in correction_supports.items()
    }
    print(
        f"PASS: exact extension degree of the pinned class-zero W_3 lift is {upper}"
    )
    print(
        "PASS: W_4 correction support is "
        f"R:{support_counts['R']}, S:{support_counts['S']}, "
        f"T:{support_counts['T']}"
    )
    print(
        "PASS: direct W_4 replay gives constant determinant "
        f"{determinant_constant} modulo 16"
    )

    if args.write_extension_certificate is not None:
        generated_directory = (
            ROOT / "artifacts" / "generated-results"
        ).resolve()
        target = args.write_extension_certificate.resolve()
        assert target.parent == generated_directory
        assert not target.exists(), f"refusing to overwrite {target}"
        certificate = {
            "schema_version": 1,
            "claim": "HKM2W4",
            "degree": upper,
            "modulus": 16,
            "fixed_w3_certificate_sha256": CLASS_ZERO_SHA256,
            "fixed_w3_degree": 19,
            "fixed_w3_jacobian_constant_mod_8": class_zero_constant,
            "minimum_extension_degree": upper,
            "degree_lower_bound": upper,
            "proved_fixed_representative_minimum": True,
            "system_at_minimum": {
                "variables": len(variable_keys),
                "equations": equation_count,
            },
            "support_counts": support_counts,
            "total_support": sum(support_counts.values()),
            "jacobian_constant_mod_16": determinant_constant,
            "supports": {
                name: [list(monomial) for monomial in support]
                for name, support in correction_supports.items()
            },
        }
        target.write_text(
            json.dumps(certificate, indent=2, sort_keys=True) + "\n"
        )
        print(f"PASS: wrote frozen W_4 extension certificate {target}")

if args.replay_extension_certificate is not None:
    path = args.replay_extension_certificate
    raw = path.read_bytes()
    if path.name in FROZEN_EXTENSION_SHA256:
        assert hashlib.sha256(raw).hexdigest() == FROZEN_EXTENSION_SHA256[path.name]
    data = json.loads(raw)
    assert data["schema_version"] == 1
    assert data["claim"] == "HKM2W4"
    assert data["degree"] == 52
    assert data["modulus"] == 16
    assert data["fixed_w3_certificate_sha256"] == CLASS_ZERO_SHA256
    assert data["proved_fixed_representative_minimum"] is True
    assert data["minimum_extension_degree"] == 52
    assert data["degree_lower_bound"] == 52
    assert data["support_counts"] == {"R": 314, "S": 98, "T": 674}
    assert data["total_support"] == 1086
    coordinates, base_data = degree_19_lift_from_certificate(
        CLASS_ZERO_CERTIFICATE,
        CLASS_ZERO_SHA256,
    )
    assert base_data["jacobian_constant_mod_8"] == 5
    coordinate_for_name = {"R": 0, "S": 1, "T": 2}
    for name, support in data["supports"].items():
        assert len(support) == data["support_counts"][name]
        for monomial in support:
            add_term(
                coordinates[coordinate_for_name[name]],
                tuple(monomial),
                8,
                16,
            )
    assert max(maximum_degree(coordinate) for coordinate in coordinates) == 52
    replayed_determinant = determinant(coordinates, 16)
    assert replayed_determinant == {
        (0, 0, 0): data["jacobian_constant_mod_16"]
    }
    assert data["jacobian_constant_mod_16"] == 13
    print("PASS: frozen degree-52 W_4 extension certificate replays directly")
