#!/usr/bin/env python3
"""Compare the first seventeen public points with the selected rank-17 cores.

The integral-coordinate, finite-quotient, and Neron-component calculations are
exact.  Canonical heights are recomputed by PARI/GP at the declared decimal
precision and are deliberately labelled numerical.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess

from compare_record_height_lattices import CURVES, CurveData, gp_vector
from pari_bridge import pari_version


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "record_first17_subgroups_v1.json"
)
CORE_PATH = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "record_rank17_core_candidates_v1.json"
)
RANK_CERTIFICATES = {
    "curve273": ROOT
    / "artifacts/generated-results/elliptic-curves/icarm_curve273_rank30_v1.json",
    "curve302": ROOT
    / "artifacts/generated-results/elliptic-curves/icarm_curve302_rank31_v1.json.gz",
}

# (prime, Kodaira symbol, component-group order, component rule).  Primes with
# trivial component group are omitted.  For curve 273 the exact factorization
# is independently checked below; curve 302's pinned certificate already
# contains and checks the complete discriminant factorization and Tate data.
COMPONENT_PLACES = {
    "curve273": (
        (2, "I16", 16, "split_multiplicative"),
        (3, "I12", 12, "split_multiplicative"),
        (5, "I8", 8, "split_multiplicative"),
        (7, "I5", 5, "split_multiplicative"),
        (13, "I5", 5, "split_multiplicative"),
        (31, "I2", 2, "split_multiplicative"),
        (41, "I2", 2, "split_multiplicative"),
        (47, "I4", 4, "split_multiplicative"),
        (53, "I3", 3, "split_multiplicative"),
        (67, "I3", 3, "split_multiplicative"),
        (379, "I2", 2, "split_multiplicative"),
    ),
    "curve302": (
        (2, "I15", 15, "split_multiplicative"),
        (3, "I4", 2, "nonsplit_multiplicative"),
        (5, "IV", 3, "additive_IV"),
        (7, "I6", 6, "split_multiplicative"),
        (11, "I4", 4, "split_multiplicative"),
        (13, "I5", 5, "split_multiplicative"),
        (19, "I2", 2, "split_multiplicative"),
        (23, "I2", 2, "split_multiplicative"),
        (29, "I3", 3, "split_multiplicative"),
        (37, "I2", 2, "split_multiplicative"),
        (41, "I2", 2, "split_multiplicative"),
        (73, "I2", 2, "split_multiplicative"),
        (131, "I2", 2, "split_multiplicative"),
        (167, "I2", 2, "split_multiplicative"),
    ),
}
CURVE273_DISCRIMINANT_FACTORIZATION = (
    (-1, 1),
    (2, 16),
    (3, 12),
    (5, 8),
    (7, 5),
    (13, 5),
    (31, 2),
    (41, 2),
    (47, 4),
    (53, 3),
    (67, 3),
    (379, 2),
    (4349, 1),
    (25721454817, 1),
    (
        97018222656318846556561979214040553412450110580812087282349817173780902099339117104673990259247421230916714670243202937,
        1,
    ),
)
THETA_THRESHOLDS = (3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5)
SHELL_RANKS = (1, 132, 328, 656, 984, 1180, 1311, 1312)


def rational_rank(rows: list[list[int]]) -> int:
    matrix = [[Fraction(value) for value in row] for row in rows]
    if not matrix:
        return 0
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pivot_value = matrix[rank][column]
        matrix[rank] = [value / pivot_value for value in matrix[rank]]
        for index in range(len(matrix)):
            if index == rank or not matrix[index][column]:
                continue
            multiplier = matrix[index][column]
            matrix[index] = [
                value - multiplier * pivot_entry
                for value, pivot_entry in zip(matrix[index], matrix[rank])
            ]
        rank += 1
    return rank


def f2_rref(rows: list[list[int]], width: int) -> tuple[list[list[int]], list[int]]:
    matrix = [[value & 1 for value in row] for row in rows if any(row)]
    rank = 0
    pivots: list[int] = []
    for column in range(width):
        pivot = next(
            (index for index in range(rank, len(matrix)) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        for index in range(len(matrix)):
            if index != rank and matrix[index][column]:
                matrix[index] = [
                    left ^ right for left, right in zip(matrix[index], matrix[rank])
                ]
        pivots.append(column)
        rank += 1
    return matrix[:rank], pivots


def f2_rank(rows: list[list[int]], width: int) -> int:
    return len(f2_rref(rows, width)[1])


def f2_nullspace(rows: list[list[int]], width: int) -> list[list[int]]:
    rref, pivots = f2_rref(rows, width)
    free = [column for column in range(width) if column not in pivots]
    basis = []
    for column in free:
        vector = [0] * width
        vector[column] = 1
        for row_index, pivot in enumerate(pivots):
            vector[pivot] = rref[row_index][column]
        basis.append(vector)
    return basis


def transpose(rows: list[list[int]]) -> list[list[int]]:
    return [list(column) for column in zip(*rows)] if rows else []


def multiply_f2(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    right_columns = transpose(right)
    return [
        [sum(a * b for a, b in zip(row, column)) & 1 for column in right_columns]
        for row in left
    ]


def read_certificate(label: str) -> dict[str, object]:
    path = RANK_CERTIFICATES[label]
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return json.load(handle)
    return json.loads(path.read_text(encoding="utf-8"))


def finite_kummer_code(label: str, point_count: int) -> dict[str, object]:
    certificate = read_certificate(label)["independence_certificate"]
    blocks = certificate.get("rows", certificate.get("signatures"))
    rows: list[list[int]] = []
    row_labels: list[dict[str, int]] = []
    local = []
    for block in blocks:
        block_rows = [[int(value) for value in row] for row in block["matrix_rows"]]
        if any(len(row) != point_count for row in block_rows):
            raise ArithmeticError("finite-quotient row width changed")
        offset = len(rows)
        rows.extend(block_rows)
        row_labels.extend(
            {"prime": int(block["prime"]), "local_coordinate": index + 1}
            for index in range(len(block_rows))
        )
        g_columns = transpose([row[:17] for row in block_rows])
        full_columns = transpose(block_rows)
        local.append(
            {
                "prime": int(block["prime"]),
                "row_offset": offset,
                "local_dimension": len(block_rows),
                "g17_image_dimension": f2_rank(g_columns, len(block_rows)),
                "displayed_image_dimension": f2_rank(full_columns, len(block_rows)),
            }
        )

    columns = transpose(rows)
    g_columns = columns[:17]
    quotient_columns = columns[17:]
    g_rank = f2_rank(g_columns, len(rows))
    full_rank = f2_rank(columns, len(rows))
    annihilator = f2_nullspace(g_columns, len(rows))
    quotient_signatures = multiply_f2(annihilator, transpose(quotient_columns))
    quotient_signatures_by_point = transpose(quotient_signatures)
    quotient_rank = f2_rank(quotient_signatures_by_point, len(annihilator))
    if g_rank != 17 or full_rank != point_count or quotient_rank != point_count - 17:
        raise ArithmeticError("the pinned mod-2 Kummer dimensions changed")
    return {
        "method": "faithful product of exact good-reduction E(F_p)/2E(F_p) quotients",
        "relation_prime": 2,
        "source_combined_row_count": len(rows),
        "source_row_labels": row_labels,
        "g17_image_dimension": g_rank,
        "displayed_image_dimension": full_rank,
        "quotient_image_dimension": quotient_rank,
        "g17_annihilator_rows": annihilator,
        "remaining_point_quotient_signatures": quotient_signatures_by_point,
        "local_blocks": local,
    }


def valuation(value: Fraction, prime: int) -> int:
    value = Fraction(value)
    if not value:
        return 10**9
    numerator = abs(value.numerator)
    denominator = value.denominator
    result = 0
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def add_points(
    left: tuple[Fraction, Fraction] | None,
    right: tuple[Fraction, Fraction] | None,
    coefficients: tuple[Fraction, ...],
) -> tuple[Fraction, Fraction] | None:
    if left is None:
        return right
    if right is None:
        return left
    x1, y1 = map(Fraction, left)
    x2, y2 = map(Fraction, right)
    a1, a2, a3, a4, a6 = map(Fraction, coefficients)
    if x1 == x2:
        if y2 == -y1 - a1 * x1 - a3:
            return None
        denominator = 2 * y1 + a1 * x1 + a3
        slope = (3 * x1**2 + 2 * a2 * x1 + a4 - a1 * y1) / denominator
        intercept = (-x1**3 + a4 * x1 + 2 * a6 - a3 * y1) / denominator
    else:
        slope = (y2 - y1) / (x2 - x1)
        intercept = (y1 * x2 - y2 * x1) / (x2 - x1)
    x3 = slope**2 + a1 * slope - a2 - x1 - x2
    y3 = -(slope + a1) * x3 - intercept - a3
    return x3, y3


def singular_reduction(
    point: tuple[Fraction, Fraction] | None,
    coefficients: tuple[Fraction, ...],
    prime: int,
) -> bool:
    if point is None:
        return False
    x_value, y_value = map(Fraction, point)
    if valuation(x_value, prime) < 0 or valuation(y_value, prime) < 0:
        return False
    a1, a2, a3, a4, _a6 = map(Fraction, coefficients)
    derivative_x = a1 * y_value - 3 * x_value**2 - 2 * a2 * x_value - a4
    derivative_y = 2 * y_value + a1 * x_value + a3
    return valuation(derivative_x, prime) > 0 and valuation(derivative_y, prime) > 0


def multiplicative_depth(
    point: tuple[Fraction, Fraction] | None,
    coefficients: tuple[Fraction, ...],
    prime: int,
    fibre_order: int,
) -> int:
    if not singular_reduction(point, coefficients, prime):
        return 0
    assert point is not None
    x_value, y_value = map(Fraction, point)
    a1, _a2, a3, _a4, _a6 = map(Fraction, coefficients)
    return min(valuation(2 * y_value + a1 * x_value + a3, prime), fibre_order // 2)


def orient_component_classes(
    curve: CurveData,
    prime: int,
    modulus: int,
    rule: str,
) -> list[int]:
    points = list(curve.points)
    if rule == "nonsplit_multiplicative":
        values = [int(singular_reduction(point, curve.coefficients, prime)) for point in points]
        status = lambda point: int(singular_reduction(point, curve.coefficients, prime))
        if any(
            (values[i] + values[j]) % 2
            != status(add_points(points[i], points[j], curve.coefficients))
            for i in range(len(points))
            for j in range(len(points))
        ):
            raise ArithmeticError("nonsplit component code failed the pair-sum replay")
        return values

    if rule == "additive_IV":
        depths = [int(singular_reduction(point, curve.coefficients, prime)) for point in points]
        options = [{0} if not depth else {1, 2} for depth in depths]
        observed_status = lambda point: int(
            singular_reduction(point, curve.coefficients, prime)
        )
        observed_depth = lambda point: observed_status(point)
    else:
        fibre_order = int(rule.removeprefix("split_multiplicative_I"))
        depths = [
            multiplicative_depth(point, curve.coefficients, prime, fibre_order)
            for point in points
        ]
        options = [{depth % modulus, (-depth) % modulus} for depth in depths]
        observed_depth = lambda point: multiplicative_depth(
            point, curve.coefficients, prime, fibre_order
        )

    values: list[int | None] = [
        next(iter(option)) if len(option) == 1 else None for option in options
    ]
    anchor = next((index for index, option in enumerate(options) if len(option) == 2), None)
    if anchor is not None:
        values[anchor] = depths[anchor] if rule != "additive_IV" else 1
    changed = True
    while changed:
        changed = False
        for index, option in enumerate(options):
            if values[index] is not None:
                continue
            candidates = []
            for candidate in option:
                valid = True
                for other, other_value in enumerate(values):
                    if other_value is None:
                        continue
                    sum_point = add_points(points[index], points[other], curve.coefficients)
                    residue = (candidate + other_value) % modulus
                    expected = int(residue != 0) if rule == "additive_IV" else min(residue, (-residue) % modulus)
                    if expected != observed_depth(sum_point):
                        valid = False
                        break
                if valid:
                    candidates.append(candidate)
            if len(candidates) == 1:
                values[index] = candidates[0]
                changed = True
    if any(value is None for value in values):
        raise ArithmeticError("component orientations were not resolved")
    oriented = [int(value) for value in values]
    for left in range(len(points)):
        for right in range(len(points)):
            residue = (oriented[left] + oriented[right]) % modulus
            expected = int(residue != 0) if rule == "additive_IV" else min(residue, (-residue) % modulus)
            observed_value = observed_depth(
                add_points(points[left], points[right], curve.coefficients)
            )
            if expected != observed_value:
                raise ArithmeticError("component code failed the pair-sum replay")
    return oriented


def cokernel_order(moduli: list[int], generators: list[list[int]]) -> tuple[list[int], int]:
    from sympy import Matrix, ZZ
    from sympy.matrices.normalforms import smith_normal_form

    size = len(moduli)
    columns = [[moduli[row] if row == column else 0 for row in range(size)] for column in range(size)]
    columns.extend(generators)
    presentation = Matrix(size, len(columns), lambda i, j: columns[j][i])
    smith = smith_normal_form(presentation, domain=ZZ)
    diagonal = [abs(int(smith[index, index])) for index in range(size)]
    invariants = [value for value in diagonal if value > 1]
    return invariants, math.prod(invariants)


def component_code(curve: CurveData) -> dict[str, object]:
    places = []
    point_columns: list[list[int]] = [[] for _ in curve.points]
    moduli = []
    for prime, symbol, component_order, rule in COMPONENT_PLACES[curve.label]:
        oriented_rule = rule
        if rule == "split_multiplicative":
            fibre_order = int(symbol[1:])
            oriented_rule = f"split_multiplicative_I{fibre_order}"
        values = orient_component_classes(
            curve, prime, component_order, oriented_rule
        )
        for column, value in zip(point_columns, values):
            column.append(value)
        moduli.append(component_order)
        places.append(
            {
                "prime": prime,
                "kodaira_symbol": symbol,
                "component_group_order": component_order,
                "rule": rule,
                "orientation": (
                    "first ambiguous submitted point is positive; all pair sums replayed"
                    if component_order > 2 and rule != "nonsplit_multiplicative"
                    else "canonical"
                ),
                "g17_classes": values[:17],
                "remaining_classes": values[17:],
            }
        )
    g_generators = point_columns[:17]
    invariants, quotient_order = cokernel_order(moduli, g_generators)
    remaining = []
    current_generators = list(g_generators)
    current_order = quotient_order
    for point_index, column in enumerate(point_columns[17:], 18):
        _new_invariants, new_order = cokernel_order(moduli, current_generators + [column])
        class_order = current_order // new_order
        remaining.append(
            {
                "point_index": point_index,
                "component_vector": column,
                "order_mod_previous_classes": class_order,
                "cumulative_component_quotient_order": new_order,
            }
        )
        current_generators.append(column)
        current_order = new_order
    return {
        "ambient_component_group_moduli": moduli,
        "g17_component_cokernel_invariants": invariants,
        "g17_component_cokernel_order": quotient_order,
        "places": places,
        "remaining_point_classes": remaining,
        "pair_sum_replay": "PASS",
    }


def parse_matrix(lines: list[str], begin: str, end: str) -> list[list[str]]:
    start = lines.index(begin) + 1
    stop = lines.index(end, start)
    return [line.split("|") for line in lines[start:stop]]


def gp_matrix(rows: list[list[str]]) -> str:
    return "[" + ";".join(",".join(row) for row in rows) + "]"


def height_and_theta(
    curve: CurveData,
    core_gram: list[list[str]],
    digits: int,
) -> dict[str, object]:
    points = "[" + ",".join(gp_vector(point) for point in curve.points[:17]) + "]"
    threshold_text = ",".join(str(value) for value in THETA_THRESHOLDS)
    rank_text = ",".join(str(value) for value in SHELL_RANKS)
    program = f"""
default(parisizemax,4000000000);
default(realprecision,{digits});
E=ellinit({gp_vector(curve.coefficients)});P={points};HG=ellheightmatrix(E,P);
HC={gp_matrix(core_gram)};
T=[{threshold_text}];K=[{rank_text}];
print("G17_BEGIN_GRAM");for(i=1,17,for(j=1,17,if(j>1,print1("|"));print1(HG[i,j]));print());print("G17_END_GRAM");
analyze(tag,G)={{my(D=matdet(G),L=exp((log(D)-log(948))/17),GW=matrix(17,17,i,j,G[i,j]),U=qflllgram(GW),R=U~*GW*U,RO=matrix(17,17,i,j,R[i,j]),Q=qfminim(R,5*L,100000,2),V);V=vecsort(vector(matsize(Q[3])[2],i,Q[3][,i]~*RO*Q[3][,i]));print(tag,"_DET|",D);print(tag,"_LAMBDA|",L);print(tag,"_LINES_TO_5LAMBDA|",#V);print1(tag,"_THETA");for(i=1,#T,print1("|",sum(j=1,#V,V[j]<=T[i]*L)));print();print1(tag,"_SHELL");for(i=1,#K,print1("|",V[K[i]]/L));print();print(tag,"_MINIMUM|",V[1]);print(tag,"_BEGIN_LLL_TRANSFORM");for(i=1,17,for(j=1,17,if(j>1,print1("|"));print1(U[i,j]));print());print(tag,"_END_LLL_TRANSFORM");print(tag,"_BEGIN_LLL_GRAM");for(i=1,17,for(j=1,17,if(j>1,print1("|"));print1(RO[i,j]));print());print(tag,"_END_LLL_GRAM");}};
analyze("G17",HG);analyze("CORE",HC);
"""
    completed = subprocess.run(
        ["gp", "-q"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    stderr = "\n".join(
        line
        for line in completed.stderr.splitlines()
        if "Warning: new" not in line and "increasing stack size" not in line
    ).strip()
    if stderr:
        raise RuntimeError(stderr)
    lines = [line.strip() for line in completed.stdout.splitlines() if line.strip()]

    def scalar(tag: str) -> str:
        return next(line.split("|", 1)[1] for line in lines if line.startswith(tag + "|"))

    def profile(tag: str) -> dict[str, object]:
        theta = next(line for line in lines if line.startswith(tag + "_THETA|"))
        shell = next(line for line in lines if line.startswith(tag + "_SHELL|"))
        theta_values = [int(value) for value in theta.split("|")[1:]]
        shell_values = shell.split("|")[1:]
        return {
            "determinant": scalar(tag + "_DET"),
            "determinant_forced_lambda_from_det_eq_948_lambda_pow_17": scalar(
                tag + "_LAMBDA"
            ),
            "minimum": scalar(tag + "_MINIMUM"),
            "unoriented_lines_through_5_lambda": int(
                scalar(tag + "_LINES_TO_5LAMBDA")
            ),
            "theta_cumulative_unoriented_lines": {
                str(bound): count for bound, count in zip(THETA_THRESHOLDS, theta_values)
            },
            "normalized_shell_height_by_nearest_rank": {
                str(rank): value for rank, value in zip(SHELL_RANKS, shell_values)
            },
            "lll_transform_columns": [
                [int(value) for value in row]
                for row in parse_matrix(
                    lines,
                    tag + "_BEGIN_LLL_TRANSFORM",
                    tag + "_END_LLL_TRANSFORM",
                )
            ],
            "lll_gram": parse_matrix(
                lines, tag + "_BEGIN_LLL_GRAM", tag + "_END_LLL_GRAM"
            ),
        }

    g17 = profile("G17")
    core = profile("CORE")
    g17_lambda = float(g17["determinant_forced_lambda_from_det_eq_948_lambda_pow_17"])
    core_lambda = float(core["determinant_forced_lambda_from_det_eq_948_lambda_pow_17"])
    return {
        "status": f"PARI canonical heights at {digits} decimal digits; not exact algebraic data",
        "g17_gram": parse_matrix(lines, "G17_BEGIN_GRAM", "G17_END_GRAM"),
        "g17": g17,
        "candidate_core": core,
        "g17_lambda_over_core_lambda": f"{g17_lambda / core_lambda:.17g}",
    }


def analyze_curve(
    curve: CurveData,
    core: dict[str, object],
    digits: int,
) -> dict[str, object]:
    point_count = len(curve.points)
    saturated_columns = [
        [int(value) for value in row]
        for row in core["saturated_basis_columns_in_public_point_coordinates"]
    ]
    core_rows = transpose(saturated_columns)
    g17_rows = [
        [int(row == column) for column in range(point_count)] for row in range(17)
    ]
    union_rank = rational_rank(g17_rows + core_rows)
    intersection_rank = 34 - union_rank
    if rational_rank(core_rows) != 17 or rational_rank(g17_rows) != 17:
        raise ArithmeticError("a rank-17 input ceased to have rank 17")

    if curve.label == "curve273":
        certificate = read_certificate(curve.label)
        discriminant = int(certificate["curve"]["discriminant"])
        product = 1
        for prime, exponent in CURVE273_DISCRIMINANT_FACTORIZATION:
            product *= prime**exponent
        if product != discriminant:
            raise ArithmeticError("curve 273 discriminant factorization changed")

    return {
        "label": curve.label,
        "displayed_rank": point_count,
        "g17_point_indices": list(range(1, 18)),
        "exact_coordinate_lattice": {
            "g17_rank": 17,
            "relative_saturation_index_in_displayed_subgroup": 1,
            "smith_invariants": [1] * 17,
            "reason": "the displayed points are an exact independent ordered basis, and G17 is its first-coordinate direct summand",
            "candidate_core_rank": 17,
            "union_rank": union_rank,
            "intersection_rank": intersection_rank,
            "remaining_quotient_rank": point_count - 17,
            "remaining_quotient_basis_point_indices": list(range(18, point_count + 1)),
            "remaining_quotient_classes": [
                {
                    "point_index": index,
                    "coordinates_in_free_quotient_basis": [
                        int(index == other) for other in range(18, point_count + 1)
                    ],
                    "order": "infinite",
                }
                for index in range(18, point_count + 1)
            ],
        },
        "canonical_height_and_theta": height_and_theta(
            curve, core["core_gram"], digits
        ),
        "finite_kummer_code": finite_kummer_code(curve.label, point_count),
        "bad_fibre_component_code": component_code(curve),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--digits", type=int, default=100)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.digits < 50:
        raise SystemExit("--digits must be at least 50")

    cores = {
        item["label"]: item
        for item in json.loads(CORE_PATH.read_text(encoding="utf-8"))["curves"]
    }
    curves = [curve for curve in CURVES if curve.label in ("curve273", "curve302")]
    import sympy

    payload = {
        "schema": "record-first17-subgroups-v1",
        "claim_scope": (
            "Exact inside the independently certified displayed free subgroups; "
            "no saturation or exact-rank claim for the full Mordell-Weil groups."
        ),
        "height_scope": (
            "Canonical-height matrices and theta profiles are high-precision "
            "PARI computations, since canonical heights are real-valued."
        ),
        "theta_normalization": "det(L)=948*lambda^17; cumulative counts use height/lambda",
        "decimal_precision_digits": args.digits,
        "software": {
            "python": platform.python_version(),
            "pari_gp": pari_version(),
            "sympy": sympy.__version__,
        },
        "generation": {
            "command": (
                "PYTHONPATH=elliptic-curves/cas .venv/bin/python "
                "elliptic-curves/cas/analyze_record_first17_subgroups.py"
            ),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "input_sha256": {
                str(CORE_PATH.relative_to(ROOT)): hashlib.sha256(
                    CORE_PATH.read_bytes()
                ).hexdigest(),
                **{
                    str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                    for path in RANK_CERTIFICATES.values()
                },
            },
        },
        "curves": [analyze_curve(curve, cores[curve.label], args.digits) for curve in curves],
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"FAIL: {args.output} differs from recomputation")
        print(
            f"PASS|{args.output}|sha256={hashlib.sha256(rendered.encode()).hexdigest()}"
        )
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered, encoding="utf-8")
    print(f"WROTE|{args.output}|sha256={hashlib.sha256(rendered.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
