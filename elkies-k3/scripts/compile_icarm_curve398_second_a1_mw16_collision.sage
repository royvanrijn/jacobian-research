#!/usr/bin/env sage-python
"""Compile curve 398's second A1/MW16 parent and compare both parents.

The complete norm-eight screen has two exact rational curve-398 survivors.  The
lower-complexity survivor (priority 16875) is already compiled.  This program
compiles the other survivor (priority 63669), specializes a saturated generic
MW16 basis, recovers its exact coordinates in the displayed public M30, and
computes the exact intersection rank, sum rank, and Smith quotient of the two
specialized generic subgroups inside M30.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, block_diagonal_matrix, matrix, pari, prod, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
TABLE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
SURVIVORS = ROOT / "artifacts/generated-results/elkies-k3-curve398-11952-norm8-a1-exact-survivors-v1.json"
TARGET = ROOT / "elliptic-curves/cas/icarm_curve398.py"
CHORD = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
SCREEN = ROOT / "elkies-k3/scripts/screen_icarm_curve398_norm8_a1_fibrations.sage"
FIRST_COMPILER = ROOT / "elkies-k3/scripts/compile_icarm_curve398_hidden_a1_mw16.sage"
FIRST_PARENT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve398_hidden_a1_mw16_v1.json"
PUBLIC_RANK = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve398_rank30_and_construction_v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elliptic-curves/icarm_curve398_two_parent_collision_v1.json"

PRIORITY_RANK = 63669
EXPECTED_PARAMETER = QQ(
    -541266381922712529166100960678122326542295329017811351186978386511278040283284966392829974955759690589708833207806994323443840
) / QQ(
    1966455527134683136777607542029510829585376789066249361045523577208160221833556912096256713936098199933472678271
)

# The first row is the new zero.  The other sixteen rows are the first
# lexicographically enumerated independent old degree-one sections.  Their
# Shioda Gram has determinant 474, so they are already a saturated MW16 basis.
SECTION_VECTORS = (
    (0, 0, 1, -1, 0, -1, -1, 0, 1, -1, 1, 1, 1, 0, 0, -1, 0),
    (-1, 0, 2, 0, 0, 0, -2, 0, -1, -1, 2, 1, -1, -1, 1, -2, 2),
    (-1, 1, 2, -1, 0, -1, -3, 0, 1, -2, 1, 1, 1, 0, 0, -1, 0),
    (-1, 1, 2, -1, 1, -1, -3, -1, 0, -1, 2, 1, 0, -1, 1, -2, 1),
    (-1, 1, 1, -1, 1, -1, -3, 0, 1, -1, 1, 1, 0, 0, 0, -1, 0),
    (0, 1, 1, -2, 0, -1, -2, 1, 1, -1, 0, 1, 1, 1, -1, 0, -1),
    (0, 0, 2, -1, 0, -1, -2, 1, 0, -1, 1, 1, 0, 0, 0, -1, 0),
    (0, 0, 2, -1, -1, 0, -2, 1, 0, -1, 1, 1, 0, 0, 0, -1, 0),
    (-1, 1, 2, -1, 0, -1, -3, 1, 0, -1, 1, 1, 0, 0, 0, -1, 0),
    (-1, 1, 2, -1, 0, 0, -3, 1, 0, -2, 1, 1, 0, 0, 0, -1, 0),
    (-1, 1, 2, -1, 1, -1, -3, 0, 0, -1, 1, 1, 0, 0, 0, -1, 0),
    (-1, 0, 2, 0, 0, 0, -2, 0, -1, -1, 2, 1, -1, -1, 1, -2, 1),
    (-1, 1, 3, -2, 0, -1, -4, 1, 1, -2, 1, 1, 1, 0, 0, -1, -1),
    (-1, 1, 2, -1, 1, 0, -3, 0, -1, -1, 1, 1, -1, 0, 1, -2, 1),
    (-1, 1, 2, -1, 0, -1, -3, 0, 1, -2, 2, 1, 1, 0, 0, -2, 0),
    (0, -1, 1, 0, -1, 1, 0, 1, -1, 0, 0, 0, -1, 0, 1, 0, 0),
    (-1, 0, 1, 1, 0, 1, 0, 0, -2, 0, 1, 0, -2, -1, 1, -1, 2),
)


def load(name: str, path: Path):
    return SourceFileLoader(name, str(path)).load_module()


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    shared = load("curve398_first_parent_compiler_helpers", FIRST_COMPILER)
    screen = load("curve398_second_parent_screen", SCREEN)
    chord = load("curve398_second_parent_chord", CHORD)
    target = load("curve398_second_parent_target", TARGET)
    model = json.loads(MODEL.read_text())
    survivor_document = json.loads(SURVIVORS.read_text())
    first_parent = json.loads(FIRST_PARENT.read_text())
    public_rank = json.loads(PUBLIC_RANK.read_text())
    if survivor_document["status"] != "PASS_EXACT_RATIONAL_PARAMETER_CANDIDATES":
        raise ArithmeticError("exact survivor input is not passing")
    if first_parent["status"] != (
        "PASS_EXACT_HIDDEN_A1_MW16_FIBRATION_PARAMETER_AND_PUBLIC_SUBGROUP"
    ):
        raise ArithmeticError("first-parent input is not passing")
    if (
        public_rank["claim"] != "rank E(Q) >= 30"
        or public_rank["independence_certificate"]["combined_binary_rank"] != 30
    ):
        raise ArithmeticError("public M30 independence input is not passing")

    sys.path.insert(0, str(ROOT / "elliptic-curves"))
    from latent_lattice.elliptic import EllipticCurve as LatentEllipticCurve
    from latent_lattice.pari import recover_exact_embedding

    table = screen.load_rows(TABLE)
    trace_row = table[PRIORITY_RANK - 1]
    trace_vector = vector(ZZ, screen.parse_vector(trace_row["section_basis_w"]))
    height_gram = matrix(ZZ, model["sections"]["height_gram"])
    if trace_vector * height_gram * trace_vector != 8:
        raise ArithmeticError("selected trace lost norm eight")
    survivor = next(record for record in survivor_document["records"] if record["priority_rank"] == PRIORITY_RANK)
    if survivor["rational_parameter_candidates"] != [
        {"multiplicity": 1, "parameter": shared.qtext(EXPECTED_PARAMETER)}
    ]:
        raise ArithmeticError("second survivor parameter disagrees with the exact screen")

    old_ring = PolynomialRing(QQ, "t")
    old_field = old_ring.fraction_field()
    old_a = old_ring(model["weierstrass_model"]["A_coefficients_low_to_high"])
    old_b = old_ring(model["weierstrass_model"]["B_coefficients_low_to_high"])
    old_curve = EllipticCurve(old_field, [old_a, old_b])
    old_basis = tuple(
        old_curve(
            screen.polynomial_from_record(record["X"], old_ring, QQ),
            screen.polynomial_from_record(record["Y"], old_ring, QQ),
        )
        for record in model["sections"]["records"]
    )
    trace = sum(
        (coefficient * point for coefficient, point in zip(trace_vector, old_basis) if coefficient),
        old_curve(0),
    )
    frame = chord.trace_chord_frame(trace[0], trace[1], old_ring)
    h, nx, ny, m0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    if h.degree() != 2:
        raise ArithmeticError("second survivor is not in the finite-pole chart")

    parameter_ring = PolynomialRing(QQ, "lambda")
    parameter = parameter_ring.gen()
    bivariate_ring = PolynomialRing(parameter_ring, "t")
    hh, nnx, nny, mm0 = map(bivariate_ring, (h, nx, ny, m0))
    slope_numerator = mm0 + parameter * hh**2
    numerator = (
        slope_numerator**4
        - 6 * slope_numerator**2 * nnx
        - 8 * slope_numerator * nny
        - 3 * nnx**2
        - 4 * bivariate_ring(old_a) * hh**4
    )
    quartic, remainder = numerator.quo_rem(hh**6)
    if remainder or quartic.degree() != 4:
        raise ArithmeticError("residual chord did not produce a binary quartic")
    invariant_i, invariant_j = screen.binary_quartic_invariants(quartic, parameter_ring)
    child_a, child_b = -27 * invariant_i, -27 * invariant_j
    child_delta = parameter_ring(-16 * (4 * child_a**3 + 27 * child_b**2))
    if [child_a.degree(), child_b.degree(), child_delta.degree()] != [8, 12, 22]:
        raise ArithmeticError("second A1 child lost its K3 degree profile")
    if child_delta.gcd(child_delta.derivative()).degree() != 0:
        raise ArithmeticError("second child finite discriminant is not squarefree")

    a1, a2, a3, a4, a6 = tuple(QQ(str(value)) for value in target.GENERAL_WEIERSTRASS_COEFFICIENTS)
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
    target_a, target_b = -27 * c4, -54 * c6
    comparison = child_a**3 * target_b**2 - target_a**3 * child_b**2
    factorization = comparison.factor()
    linear = [factor for factor, multiplicity in factorization if factor.degree() == 1 for _ in range(multiplicity)]
    if len(linear) != 1:
        raise ArithmeticError("second trace does not have one rational curve-398 parameter")
    recovered_parameter = -linear[0][0] / linear[0][1]
    if recovered_parameter != EXPECTED_PARAMETER:
        raise ArithmeticError("second curve-398 parameter changed")
    specialized_a = QQ(child_a(recovered_parameter))
    specialized_b = QQ(child_b(recovered_parameter))
    child_curve = EllipticCurve(QQ, [specialized_a, specialized_b])
    target_short = EllipticCurve(QQ, [target_a, target_b])
    if not child_curve.is_isomorphic(target_short):
        raise ArithmeticError("second j-match is a nontrivial quadratic twist")
    child_to_target = child_curve.isomorphism_to(target_short)

    fibre = vector(ZZ, [2, 2] + list(trace_vector))
    old_zero_class = vector(ZZ, [-1, 1] + [0] * 17)
    trace_class = shared.section_class(trace_vector, height_gram)
    ns = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -height_gram)
    if fibre != old_zero_class + trace_class or fibre * ns * fibre or fibre * ns * old_zero_class:
        raise ArithmeticError("second D=O+P trace decomposition failed")
    degree_one = shared.enumerate_degree_one_vectors(height_gram, trace_vector)
    if len(degree_one) != 180:
        raise ArithmeticError(f"second degree-one old-section count changed: {len(degree_one)}")
    selected = tuple(vector(ZZ, row) for row in SECTION_VECTORS)
    if any(value not in degree_one for value in selected):
        raise ArithmeticError("a frozen second-parent section left the complete degree-one shell")
    generic_gram = shared.shioda_gram(selected[1:], selected[0], trace_vector, height_gram)
    if generic_gram.det() != 474 or generic_gram.rank() != 16:
        raise ArithmeticError("second selected MW16 basis is not saturated")

    new_zero_class = shared.section_class(selected[0], height_gram)
    mate = fibre + new_zero_class
    complement = matrix(ZZ, [list(fibre * ns), list(mate * ns)]).right_kernel_matrix()
    transport = matrix(QQ, [list(fibre), list(mate)] + [list(row) for row in complement])
    child_frame = -(complement * ns * complement.transpose())
    root_count = int(pari(matrix(ZZ, child_frame)).qfminim(2)[0])
    if abs(transport.det()) != 1 or child_frame.det() != 948 or root_count != 2:
        raise ArithmeticError("second frame is not primitive with one A1 root")

    fixed_m = m0 + recovered_parameter * h**2
    fixed_quartic = old_ring([QQ(quartic[index](recovered_parameter)) for index in range(5)])
    sum_x = old_ring((fixed_m**2 - nx) // h**2)
    quartic_points = []
    base_maps = []
    for section_vector in selected:
        source_point = sum(
            (coefficient * point for coefficient, point in zip(section_vector, old_basis) if coefficient),
            old_curve(0),
        )
        source_x, source_y = source_point[0], source_point[1]
        base_map = old_field(
            (((source_y + trace[1]) / (source_x - trace[0])) * h - m0) / h**2
        )
        old_parameter = shared.invert_mobius(base_map, recovered_parameter, old_ring)
        x_value = QQ(source_x(old_parameter))
        w_value = (2 * x_value - QQ(sum_x(old_parameter))) / QQ(h(old_parameter))
        if w_value**2 != fixed_quartic(old_parameter):
            raise ArithmeticError("second old section missed the specialized quartic")
        quartic_points.append((old_parameter, w_value))
        base_maps.append(base_map)

    t0, w0 = quartic_points[0]
    shift_ring = PolynomialRing(QQ, "z")
    z = shift_ring.gen()
    shifted = shift_ring(fixed_quartic(t0 + z))
    ee, dd, cc, bb, aa = [QQ(shifted[index]) for index in range(5)]
    if ee != w0**2 or not w0:
        raise ArithmeticError("second pointed quartic origin is invalid")
    a1g = dd / w0
    a2g = cc - dd**2 / (4 * w0**2)
    a3g = 2 * w0 * bb
    a4g = -4 * w0**2 * aa
    a6g = a2g * a4g
    b2g = a1g**2 + 4 * a2g
    b4g = a1g * a3g + 2 * a4g
    b6g = a3g**2 + 4 * a6g
    c4g = b2g**2 - 24 * b4g
    c6g = -b2g**3 + 36 * b2g * b4g - 216 * b6g
    if 81 * (-c4g / 48) != specialized_a or 729 * (-c6g / 864) != specialized_b:
        raise ArithmeticError("second pointed quartic normalization missed the invariant model")

    public_curve = EllipticCurve(QQ, [a1, a2, a3, a4, a6])
    target_to_public = target_short.isomorphism_to(public_curve)
    generic_public_points = []
    for old_parameter, w_value in quartic_points[1:]:
        zz = old_parameter - t0
        if not zz:
            raise ArithmeticError("two second-parent sections meet at the quartic origin")
        x_general = (2 * w0 * (w_value + w0) + dd * zz) / zz**2
        y_general = (
            4 * w0**2 * (w_value + w0)
            + 2 * w0 * dd * zz
            + (2 * w0 * cc - dd**2 / (2 * w0)) * zz**2
        ) / zz**3
        child_point = child_curve(
            9 * (x_general + b2g / 12),
            27 * (y_general + (a1g * x_general + a3g) / 2),
        )
        generic_public_points.append(target_to_public(child_to_target(child_point)))

    latent_curve = LatentEllipticCurve(tuple(Fraction(str(value)) for value in (a1, a2, a3, a4, a6)))
    public_basis = tuple((Fraction(str(x)), Fraction(str(y))) for x, y in target.POINTS)
    generic_public_affine = tuple(
        (Fraction(str(point[0])), Fraction(str(point[1]))) for point in generic_public_points
    )
    embedding = recover_exact_embedding(
        latent_curve, public_basis, generic_public_affine, digits=180, timeout=600.0
    )
    if len(embedding.columns) != 16:
        raise ArithmeticError("second public subgroup embedding has the wrong rank")

    first_rows = [list(map(ZZ, row)) for row in first_parent["public_rank30_embedding"]["matrix_30_by_16_columns"]]
    second_rows = [list(map(ZZ, row)) for row in embedding.columns]
    first_matrix = matrix(ZZ, first_rows)
    second_matrix = matrix(ZZ, second_rows)
    sum_matrix = first_matrix.stack(second_matrix)
    first_rank = int(first_matrix.rank())
    second_rank = int(second_matrix.rank())
    sum_rank = int(sum_matrix.rank())
    intersection_rank = int(first_rank + second_rank - sum_rank)
    smith, smith_left, smith_right = sum_matrix.smith_form(transformation=True)
    if smith_left * sum_matrix * smith_right != smith:
        raise ArithmeticError("two-parent Smith transformation failed exact replay")
    smith_factors = [abs(ZZ(smith[index, index])) for index in range(sum_rank)]
    smith_index = prod(smith_factors) if sum_rank == 30 else None
    if first_rank != 16 or second_rank != 16:
        raise ArithmeticError("a specialized generic parent lost rank")
    second_in_first = first_matrix.transpose().solve_right(second_matrix.transpose()).transpose()
    first_in_second = second_matrix.transpose().solve_right(first_matrix.transpose()).transpose()
    integral_equality = (
        all(value in ZZ for value in second_in_first.list())
        and all(value in ZZ for value in first_in_second.list())
        and second_in_first * first_matrix == second_matrix
        and first_in_second * second_matrix == first_matrix
        and second_in_first * first_in_second == 1
        and first_in_second * second_in_first == 1
        and abs(second_in_first.det()) == 1
        and abs(first_in_second.det()) == 1
    )
    if not integral_equality:
        raise ArithmeticError("the observed equality of the two specialized integral MW16 groups failed replay")

    records = []
    for index, (section_vector, base_map, quartic_point, public_point, row) in enumerate(
        zip(selected[1:], base_maps[1:], quartic_points[1:], generic_public_points, second_rows),
        start=1,
    ):
        source_class = shared.section_class(section_vector, height_gram)
        records.append(
            {
                "basis_index_one_based": index,
                "source_section_basis_coordinates": list(map(int, section_vector)),
                "source_section_height": int(section_vector * height_gram * section_vector),
                "source_section_intersection_with_D": int(shared.ns_intersection(source_class, fibre, height_gram)),
                "new_fibre_component": "nonidentity-O" if source_class * ns * old_zero_class else "identity-P_w",
                "base_map_lambda_of_t": shared.rational_function_record(base_map),
                "specialized_quartic_point": {"t": shared.qtext(quartic_point[0]), "W": shared.qtext(quartic_point[1])},
                "specialized_public_point": shared.point_record(public_point),
                "coordinates_in_public_rank30_points": list(map(int, row)),
                "exact_public_group_law_replay": True,
            }
        )

    result = {
        "schema": "elliptic-curves.icarm-curve398-two-parent-collision.v1",
        "status": "PASS_EXACT_TWO_PARENT_COLLISION",
        "curve_id": 398,
        "ambient_group": {
            "name": "public M30",
            "ordered_basis_source": relative(TARGET),
            "rank": 30,
            "independence_certificate": relative(PUBLIC_RANK),
        },
        "first_parent": {
            "compiled_artifact": relative(FIRST_PARENT),
            "priority_rank": int(first_parent["fibration"]["priority_rank"]),
            "orbit_hex": first_parent["fibration"]["orbit_hex"],
            "rank": first_rank,
            "matrix_16_by_30_rows": [list(map(int, row)) for row in first_rows],
        },
        "second_parent": {
            "priority_rank": PRIORITY_RANK,
            "orbit_mask": int(trace_row["orbit_mask"]),
            "orbit_hex": trace_row["orbit_hex"],
            "trace_section_basis_w": list(map(int, trace_vector)),
            "divisor_class_in_U_plus_M_minus": list(map(int, fibre)),
            "divisor_identity": "D=(2,2,w)=O+P_w",
            "finite_pole_degree": int(h.degree()),
            "equation": "Y^2=X^3+A(lambda)*X+B(lambda)",
            "A_coefficients_low_to_high": shared.poly_record(child_a),
            "B_coefficients_low_to_high": shared.poly_record(child_b),
            "degrees_A_B_Delta": [8, 12, 22],
            "fibre_configuration": "I2 at infinity + 22 I1",
            "finite_discriminant_squarefree": True,
            "child_frame_determinant": int(child_frame.det()),
            "child_frame_norm_two_vector_count_signed": root_count,
            "parameter_recovery": {
                "comparison_polynomial_degree": int(comparison.degree()),
                "factor_degrees_with_multiplicity": [
                    [int(factor.degree()), int(multiplicity)] for factor, multiplicity in factorization
                ],
                "lambda": shared.qtext(recovered_parameter),
                "specialized_child_short_coefficients": [shared.qtext(specialized_a), shared.qtext(specialized_b)],
                "isomorphic_to_curve398_over_Q": True,
                "child_to_curve398_short_isomorphism_u_r_s_t": [shared.qtext(value) for value in child_to_target.tuple()],
            },
            "generic_mw16": {
                "complete_old_degree_one_section_count": len(degree_one),
                "zero_source_section_basis_coordinates": list(map(int, selected[0])),
                "height_gram": [[shared.qtext(value) for value in row] for row in generic_gram.rows()],
                "height_gram_determinant": shared.qtext(generic_gram.det()),
                "rank": 16,
                "saturated": True,
                "records": records,
            },
            "public_rank30_embedding": {
                "orientation": "rows give each generic MW16 point in the ordered public 30-point list",
                "matrix_16_by_30_rows": [list(map(int, row)) for row in second_rows],
                "maximum_absolute_coordinate": embedding.max_abs_coordinate,
                "nonzero_coordinate_count": embedding.nonzero_coordinates,
                "height_dual_numerical_residual_max": embedding.numerical_residual_max,
                "exact_group_law_replay": True,
            },
        },
        "collision": {
            "rank_G1": first_rank,
            "rank_G2": second_rank,
            "rank_intersection": intersection_rank,
            "rank_sum": sum_rank,
            "integral_subgroups_equal": integral_equality,
            "G2_basis_rows_in_G1_basis": [[int(value) for value in row] for row in second_in_first.rows()],
            "G1_basis_rows_in_G2_basis": [[int(value) for value in row] for row in first_in_second.rows()],
            "basis_transition_determinants": [int(second_in_first.det()), int(first_in_second.det())],
            "sum_generator_matrix_32_by_30_rows": [list(map(int, row)) for row in sum_matrix.rows()],
            "smith_diagonal_nonzero": list(map(int, smith_factors)),
            "quotient_free_rank": 30 - sum_rank,
            "quotient_torsion_invariant_factors_nontrivial": [int(value) for value in smith_factors if value > 1],
            "smith_index_in_public_M30": "infinite" if smith_index is None else str(int(smith_index)),
            "smith_transformation_exact_replay": True,
        },
        "inputs": {
            relative(path): digest(path)
            for path in (MODEL, TABLE, SURVIVORS, TARGET, CHORD, SCREEN, FIRST_COMPILER, FIRST_PARENT, PUBLIC_RANK)
        },
        "generation": {
            "script": relative(Path(__file__)),
            "script_sha256": digest(Path(__file__)),
        },
        "software": {"sage_version": SAGE_VERSION, "pari_version": ".".join(map(str, pari.version()))},
        "proof_boundary": (
            "This certifies the second A1/MW16 fibration, its curve-398 specialization, and the exact "
            "relative position of both specialized generic MW16 subgroups inside the displayed public M30. "
            "It does not claim that M30 is saturated in E(Q), nor an unconditional rank upper bound for E(Q)."
        ),
        "reproducing_command": "sage -python elkies-k3/scripts/compile_icarm_curve398_second_a1_mw16_collision.sage --check",
    }
    output_text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.output.read_text() != output_text:
            raise ArithmeticError("stored curve-398 two-parent collision artifact differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(output_text)
    print(
        f"CURVE398COLLISION|G1={first_rank}|G2={second_rank}|intersection={intersection_rank}|"
        f"sum={sum_rank}|smith_index={'infinite' if smith_index is None else smith_index}|"
        f"output={relative(args.output)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
