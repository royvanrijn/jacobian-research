#!/usr/bin/env sage-python
"""Prepare complement-blind MW16 inputs for all nine ICARM A1 presentations.

The complete 11952 norm-eight atlas records nine compiled A1/MW16
presentations through five target curves.  This replay reconstructs each
binary-quartic pencil, selects a deterministic unimodular basis from its
complete shell of old degree-one sections, and specializes those sixteen
sections at the recorded target fibre.

No public point list, public Mordell--Weil basis, displayed complement, or
target-rank outcome is read.  The output is therefore suitable as the sole
input to a point-search detector that is blind to the held-out complements.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
from pathlib import Path
import platform
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
MODEL = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
TABLE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
)
CHORD = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
SCREEN = ROOT / "elkies-k3/scripts/screen_icarm_curve398_norm8_a1_fibrations.sage"
LATTICE = ROOT / "elkies-k3/scripts/compile_icarm_curve398_hidden_a1_mw16.sage"
COMPILED_TEMPLATE = (
    "artifacts/generated-results/"
    "elkies-k3-icarm-curve{curve_id}-11952-norm8-a1-compiled-survivors-v1.json"
)
OUTPUT = (
    ROOT
    / "elliptic-curves/data/icarm_mw16_parent_ladder_blind_inputs_v1.json"
)
CURVE_IDS = (398, 400, 401, 542, 548)


def load_source(name: str, path: Path):
    loader = SourceFileLoader(name, str(path))
    spec = spec_from_loader(name, loader)
    if spec is None:
        raise ImportError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def qtext(value) -> str:
    value = QQ(value)
    return (
        str(value.numerator())
        if value.denominator() == 1
        else f"{value.numerator()}/{value.denominator()}"
    )


def poly_record(poly):
    return [qtext(poly[index]) for index in range(poly.degree() + 1)] if poly else ["0"]


def point_record(point):
    return {"x": qtext(point[0]), "y": qtext(point[1])}


def rational_function_record(value):
    return {
        "numerator_coefficients_low_to_high": poly_record(value.numerator()),
        "denominator_coefficients_low_to_high": poly_record(value.denominator()),
    }


def deterministic_unimodular_sections(
    degree_one, trace_vector, new_zero_vector, height_gram, fibration, lattice
):
    """Select sixteen actual old sections forming the compiled MW quotient."""

    fibre = vector(QQ, fibration["divisor_class_in_U_plus_M_minus"])
    new_zero_class = lattice.section_class(new_zero_vector, height_gram)
    root = vector(
        QQ,
        fibration["generic_mordell_weil"]["a1_root_source_ns_coordinates"],
    )
    quotient_rows = [
        vector(QQ, row)
        for row in fibration["generic_mordell_weil"][
            "quotient_basis_source_ns_coordinates"
        ]
    ]
    marked_basis = matrix(QQ, [fibre, fibre + new_zero_class, root] + quotient_rows)
    if abs(marked_basis.det()) != 1:
        raise ArithmeticError("compiled U+A1+MW marking is not unimodular")
    inverse = marked_basis.inverse()
    quotient_coordinates = []
    for section in degree_one:
        coordinates = lattice.section_class(section, height_gram) * inverse
        if any(value not in ZZ for value in coordinates):
            raise ArithmeticError("an old section is nonintegral in the compiled marking")
        quotient_coordinates.append(vector(ZZ, coordinates[3:]))
    coordinate_matrix = matrix(ZZ, quotient_coordinates)
    if coordinate_matrix.rank() != 16 or coordinate_matrix.elementary_divisors()[:16] != [1] * 16:
        raise ArithmeticError("complete old-section shell does not generate saturated MW16")

    order = sorted(
        range(len(degree_one)),
        key=lambda index: (
            sum(abs(value) for value in quotient_coordinates[index]),
            tuple(quotient_coordinates[index]),
            tuple(degree_one[index]),
        ),
    )
    for shift in range(len(order)):
        shifted = order[shift:] + order[:shift]
        chosen = []
        for index in shifted:
            trial = matrix(
                ZZ, [quotient_coordinates[row] for row in chosen + [index]]
            )
            if trial.rank() > len(chosen):
                chosen.append(index)
            if len(chosen) == 16:
                break
        coordinate_change = matrix(
            ZZ, [quotient_coordinates[index] for index in chosen]
        )
        if abs(coordinate_change.det()) == 1:
            selected = tuple(degree_one[index] for index in chosen)
            generic_gram = lattice.shioda_gram(
                selected,
                new_zero_vector,
                trace_vector,
                height_gram,
            )
            compiled_gram = matrix(
                QQ, fibration["generic_mordell_weil"]["height_gram"]
            )
            if generic_gram != coordinate_change * compiled_gram * coordinate_change.transpose():
                raise ArithmeticError("selected-section height form missed the compiled quotient")
            if generic_gram.det() != 474 or generic_gram.rank() != 16:
                raise ArithmeticError("selected sections are not a saturated MW16 basis")
            return selected, coordinate_change, generic_gram, shift
    raise ArithmeticError("no deterministic unimodular old-section subset was found")


def child_geometry(trace, old_a, old_b, old_ring, chord, screen):
    frame = chord.trace_chord_frame(trace[0], trace[1], old_ring)
    h, nx, ny, m0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    if h.degree() != 2:
        raise ArithmeticError("atlas presentation left the finite-pole chart")
    parameter_ring = PolynomialRing(QQ, "lambda")
    parameter_variable = parameter_ring.gen()
    bivariate_ring = PolynomialRing(parameter_ring, "t")
    hh, nnx, nny, mm0 = map(bivariate_ring, (h, nx, ny, m0))
    slope_numerator = mm0 + parameter_variable * hh**2
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
    invariant_i, invariant_j = screen.binary_quartic_invariants(
        quartic, parameter_ring
    )
    return h, nx, m0, quartic, -27 * invariant_i, -27 * invariant_j


def specialize_sections(
    *,
    parameter,
    selected,
    new_zero_vector,
    trace,
    old_basis,
    old_curve,
    old_a,
    old_ring,
    h,
    nx,
    m0,
    quartic,
    child_a,
    child_b,
    target_short,
    lattice,
):
    fixed_m = m0 + parameter * h**2
    fixed_quartic = old_ring(
        [QQ(quartic[index](parameter)) for index in range(5)]
    )
    sum_x = old_ring((fixed_m**2 - nx) // h**2)
    section_vectors = (new_zero_vector,) + tuple(selected)
    quartic_points = []
    base_maps = []
    for section_vector in section_vectors:
        source_point = sum(
            (
                coefficient * point
                for coefficient, point in zip(section_vector, old_basis)
                if coefficient
            ),
            old_curve(0),
        )
        source_x, source_y = source_point[0], source_point[1]
        base_map = old_curve.base_field()(
            (((source_y + trace[1]) / (source_x - trace[0])) * h - m0) / h**2
        )
        old_parameter = lattice.invert_mobius(base_map, parameter, old_ring)
        x_value = QQ(source_x(old_parameter))
        y_value = QQ(source_y(old_parameter))
        w_value = (2 * x_value - QQ(sum_x(old_parameter))) / QQ(h(old_parameter))
        if w_value**2 != fixed_quartic(old_parameter):
            raise ArithmeticError("old section missed the specialized quartic")
        quartic_points.append((old_parameter, w_value))
        base_maps.append(base_map)

    t0, w0 = quartic_points[0]
    shift_ring = PolynomialRing(QQ, "z")
    z = shift_ring.gen()
    shifted = shift_ring(fixed_quartic(t0 + z))
    ee, dd, cc, bb, aa = [QQ(shifted[index]) for index in range(5)]
    if ee != w0**2 or not w0:
        raise ArithmeticError("pointed quartic origin is invalid")
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
    specialized_a = QQ(child_a(parameter))
    specialized_b = QQ(child_b(parameter))
    if 81 * (-c4g / 48) != specialized_a or 729 * (-c6g / 864) != specialized_b:
        raise ArithmeticError("pointed quartic normalization missed the child model")
    child_curve = EllipticCurve(QQ, [specialized_a, specialized_b])
    if not child_curve.is_isomorphic(target_short):
        raise ArithmeticError("recorded target specialization became a twist")
    child_to_target = child_curve.isomorphism_to(target_short)

    target_points = []
    for old_parameter, w_value in quartic_points[1:]:
        zz = old_parameter - t0
        if not zz:
            raise ArithmeticError("selected section meets the pointed origin")
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
        target_points.append(child_to_target(child_point))
    if len(target_points) != 16 or any(point.curve() != target_short for point in target_points):
        raise ArithmeticError("specialized MW16 point list is incomplete")
    return tuple(target_points), tuple(base_maps), child_to_target


def build():
    screen = load_source("mw16_ladder_screen", SCREEN)
    chord = load_source("mw16_ladder_chord", CHORD)
    lattice = load_source("mw16_ladder_lattice", LATTICE)
    model = json.loads(MODEL.read_text())
    table = screen.load_rows(TABLE)
    height_gram = matrix(ZZ, model["sections"]["height_gram"])
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

    compiled_paths = [
        ROOT / COMPILED_TEMPLATE.format(curve_id=curve_id)
        for curve_id in CURVE_IDS
    ]
    parents = []
    curve_parent_counts = {}
    for curve_id, compiled_path in zip(CURVE_IDS, compiled_paths):
        compiled = json.loads(compiled_path.read_text())
        if compiled.get("status") != "PASS_EXACT_COMPILED_A1_MW16_FIBRATIONS":
            raise ArithmeticError(f"curve {curve_id} compiled atlas status changed")
        if int(compiled["target"]["curve_id"]) != curve_id:
            raise ArithmeticError("compiled atlas target changed")
        target_c4 = QQ(compiled["target"]["c4"])
        target_c6 = QQ(compiled["target"]["c6"])
        target_short = EllipticCurve(QQ, [-27 * target_c4, -54 * target_c6])
        curve_parent_counts[str(curve_id)] = len(compiled["fibrations"])
        for fibration in compiled["fibrations"]:
            priority_rank = int(fibration["priority_rank"])
            trace_vector = vector(
                ZZ, screen.parse_vector(table[priority_rank - 1]["section_basis_w"])
            )
            if list(trace_vector) != fibration["trace_section_basis_w"]:
                raise ArithmeticError("compiled trace and priority table disagree")
            trace = sum(
                (
                    coefficient * point
                    for coefficient, point in zip(trace_vector, old_basis)
                    if coefficient
                ),
                old_curve(0),
            )
            h, nx, m0, quartic, child_a, child_b = child_geometry(
                trace, old_a, old_b, old_ring, chord, screen
            )
            equation = fibration["equation"]
            if poly_record(child_a) != equation["A_coefficients_low_to_high"] or poly_record(child_b) != equation["B_coefficients_low_to_high"]:
                raise ArithmeticError("reconstructed pencil equation differs from atlas")
            specializations = fibration["target_specializations"]
            if len(specializations) != 1:
                raise ArithmeticError("ladder requires one recorded target parameter per parent")
            parameter = QQ(specializations[0]["lambda"])
            new_zero_vector = vector(
                ZZ, fibration["new_zero_source_section_basis_coordinates"]
            )
            degree_one = lattice.enumerate_degree_one_vectors(
                height_gram, trace_vector
            )
            if new_zero_vector not in degree_one:
                raise ArithmeticError("compiled new zero left the degree-one shell")
            selected, coordinate_change, generic_gram, selection_shift = (
                deterministic_unimodular_sections(
                    degree_one,
                    trace_vector,
                    new_zero_vector,
                    height_gram,
                    fibration,
                    lattice,
                )
            )
            target_points, base_maps, child_to_target = specialize_sections(
                parameter=parameter,
                selected=selected,
                new_zero_vector=new_zero_vector,
                trace=trace,
                old_basis=old_basis,
                old_curve=old_curve,
                old_a=old_a,
                old_ring=old_ring,
                h=h,
                nx=nx,
                m0=m0,
                quartic=quartic,
                child_a=child_a,
                child_b=child_b,
                target_short=target_short,
                lattice=lattice,
            )
            parent_id = f"curve{curve_id}-p{priority_rank}"
            parents.append(
                {
                    "parent_id": parent_id,
                    "curve_id": curve_id,
                    "priority_rank": priority_rank,
                    "orbit_hex": fibration["orbit_hex"],
                    "target_parameter": qtext(parameter),
                    "target_short_model": [
                        "0",
                        "0",
                        "0",
                        qtext(target_short.a4()),
                        qtext(target_short.a6()),
                    ],
                    "specialized_generic_points": [
                        point_record(point) for point in target_points
                    ],
                    "generic_height_gram": [
                        [qtext(value) for value in row]
                        for row in generic_gram.rows()
                    ],
                    "generic_rank": 16,
                    "pencil": {
                        "coordinate": "lambda",
                        "A_coefficients_low_to_high": poly_record(child_a),
                        "B_coefficients_low_to_high": poly_record(child_b),
                        "degrees_A_B_Delta": [8, 12, 22],
                        "fibre_configuration": "I2 at infinity + 22 I1",
                    },
                    "source_marking": {
                        "trace_section_basis_w": list(map(int, trace_vector)),
                        "new_zero_source_section_basis_coordinates": list(
                            map(int, new_zero_vector)
                        ),
                        "generic_source_section_basis_coordinates": [
                            list(map(int, section)) for section in selected
                        ],
                        "generic_coordinates_in_compiled_mw_basis": [
                            list(map(int, row)) for row in coordinate_change.rows()
                        ],
                        "complete_old_degree_one_section_count": len(degree_one),
                        "deterministic_cyclic_order_shift": selection_shift,
                        "base_maps_lambda_of_old_t": [
                            rational_function_record(value) for value in base_maps
                        ],
                    },
                    "target_isomorphism": {
                        "child_to_target_short_u_r_s_t": [
                            qtext(value) for value in child_to_target.tuple()
                        ],
                        "exact": True,
                    },
                }
            )
            print(
                f"MW16LADDERINPUT|parent={parent_id}|degree1={len(degree_one)}|"
                f"basis_shift={selection_shift}|det=474|status=PASS",
                flush=True,
            )

    if len(parents) != 9 or curve_parent_counts != {
        "398": 2,
        "400": 2,
        "401": 1,
        "542": 1,
        "548": 3,
    }:
        raise ArithmeticError("the nine-presentation/five-curve atlas partition changed")
    return {
        "schema": "elliptic-curves.icarm-mw16-parent-ladder-blind-inputs.v1",
        "status": "PASS_EXACT_COMPLEMENT_BLIND_NINE_PARENT_INPUTS",
        "observation_unit": {
            "primary": "target_curve",
            "curve_count": 5,
            "curve_ids": list(CURVE_IDS),
            "presentation_count": 9,
            "parent_counts_by_curve": curve_parent_counts,
            "pseudoreplication_rule": (
                "multiple parent presentations of one target curve are retained as "
                "within-curve detector trials and never counted as independent target outcomes"
            ),
        },
        "blindness_boundary": {
            "public_point_lists_loaded": False,
            "public_complement_coordinates_loaded": False,
            "target_rank_lower_bounds_loaded": False,
            "permitted_information": (
                "five target equations, nine exact target parameters, nine pencil equations, "
                "and sixteen specialized generic sections with their exact height forms"
            ),
            "description": (
                "complement-blind replay inputs prepared after the atlas outcomes were known; "
                "this is calibration, not a genuinely prospective outcome-blind sample"
            ),
        },
        "parents": sorted(
            parents, key=lambda row: (row["curve_id"], row["priority_rank"])
        ),
        "inputs": {
            relative(path): digest(path)
            for path in [Path(__file__), MODEL, TABLE, CHORD, SCREEN, LATTICE]
            + compiled_paths
        },
        "software": {
            "python": platform.python_version(),
            "sage": SAGE_VERSION,
        },
        "claim_boundary": [
            "Every stored point is an exact specialization of a section in a saturated generic MW16 basis.",
            "No public exceptional point or displayed target complement is used.",
            "The nine rows are parent presentations nested within five curve-level calibration observations.",
            "The two curve-398 rows are retained as coordinate/presentation trials even though their exact base equivalence is already certified.",
            "No rank beyond sixteen, exact target rank, point-search success, or Selmer statement is asserted.",
        ],
        "reproducing_command": (
            "sage -python elliptic-curves/cas/prepare_icarm_mw16_parent_ladder_inputs.sage --check"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = build()
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not args.output.is_file() or args.output.read_text() != serialized:
            raise ArithmeticError("stored MW16 parent-ladder input differs from replay")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized)
    print(
        f"MW16LADDERINPUT|curves=5|presentations=9|output={relative(args.output)}|"
        "status=PASS_EXACT_COMPLEMENT_BLIND_NINE_PARENT_INPUTS",
        flush=True,
    )


if __name__ == "__main__":
    main()
