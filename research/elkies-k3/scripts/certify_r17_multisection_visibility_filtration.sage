#!/usr/bin/env sage-python
"""Certify the two natural R17 multisection-visibility filtrations.

The literal filtration over all geometrically integral multisections includes
positive-genus curves.  The existing norm-eight genus-one bisection pencil can
be fitted through every displayed exceptional generator at the four published
rank-25--28 controls, so that filtration is already full in degree two.

The finite rational-curve filtration is different.  This script imports the
complete rational-bisection ranks and the explicitly bounded degree-three and
degree-four experiment, preserving their completeness flags.  It deliberately
does not relabel the sampled degree-three/four layers as complete.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
PILOT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json"
)
RELATIONS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_exceptional_specialization_relations_v1.json"
)
RATIONAL_EXPERIMENT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_deep_cover_exceptional_quotients_v1.json"
)
CHORD_SCRIPT = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
PUBLIC_SCRIPTS = {
    "-2/377": ROOT / "elliptic-curves/cas/elkies_rank25.py",
    "-308/251": ROOT / "elliptic-curves/cas/elkies_rank26.py",
    "2456/135": ROOT / "elliptic-curves/cas/elkies_rank27.py",
    "-9529/5471": ROOT / "elliptic-curves/cas/elkies_rank28.py",
}
PARAMETERS = {
    "-2/377": QQ(-2) / QQ(377),
    "-308/251": QQ(-308) / QQ(251),
    "2456/135": QQ(2456) / QQ(135),
    "-9529/5471": QQ(-9529) / QQ(5471),
}
OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-multisection-visibility-filtration-v1.json"
)


def load_script(name: str, path: Path):
    loader = SourceFileLoader(name, str(path))
    spec = spec_from_loader(name, loader)
    if spec is None:
        raise ImportError(f"cannot load {path}")
    module = module_from_spec(spec)
    loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def rational_text(value) -> str:
    value = QQ(value)
    if value.denominator() == 1:
        return str(value.numerator())
    return f"{value.numerator()}/{value.denominator()}"


def polynomial_coefficients(polynomial) -> list[str]:
    return [rational_text(polynomial[index]) for index in range(polynomial.degree() + 1)]


def fit_target(
    *, target, parameter, X, Y, h, Nx, Ny, M0, A, B, discriminant,
    polynomial_ring, function_field, chord_module,
):
    target_point = target["point"]
    target_x, target_y = target_point[0], target_point[1]
    trace_x = QQ(X(parameter))
    trace_y = QQ(Y(parameter))
    h_value = QQ(h(parameter))
    if not h_value or target_x == trace_x:
        raise ArithmeticError(f"degenerate target incidence for {target['label']}")
    target_slope = (target_y + trace_y) / (target_x - trace_x)
    pencil_parameter = (target_slope * h_value - QQ(M0(parameter))) / h_value**2
    M = M0 + pencil_parameter * h**2
    data = chord_module.chord_data_from_slope_numerator(
        h,
        Nx,
        Ny,
        M,
        A,
        B,
        discriminant,
        polynomial_ring,
        function_field,
        expected_q_degree=4,
    )
    q = data["q"]
    if q.degree() != 4 or not q.is_squarefree():
        raise ArithmeticError(f"branch polynomial for {target['label']} is not squarefree quartic")
    if q.gcd(discriminant).degree() or q.gcd(h).degree():
        raise ArithmeticError(f"branch quartic for {target['label']} meets a forbidden divisor")

    sum_x = data["sum_x"]
    product_x = data["product_x"]
    cover_coordinate = (2 * target_x - QQ(sum_x(parameter))) / h_value
    line_value = QQ(M(parameter)) / h_value * target_x - (
        trace_y + QQ(M(parameter)) / h_value * trace_x
    )
    if cover_coordinate**2 != QQ(q(parameter)):
        raise ArithmeticError("target does not give a rational cover witness")
    if line_value != target_y:
        raise ArithmeticError("target y-coordinate misses the fitted line")
    if target_x**2 - QQ(sum_x(parameter)) * target_x + QQ(product_x(parameter)):
        raise ArithmeticError("target x-coordinate misses the residual quadratic")
    if target_y**2 != target_x**3 + A(parameter) * target_x + B(parameter):
        raise ArithmeticError("short target point fails its Weierstrass equation")
    if not (
        data["x0"].degree() <= 4
        and data["x1"].degree() <= 2
        and data["y0"].degree() <= 6
        and data["y1"].degree() <= 4
    ):
        raise ArithmeticError("lifted section exceeds the integral degree bounds")

    factor_degrees = sorted(
        int(factor.degree())
        for factor, exponent in q.factor()
        for _ in range(int(exponent))
    )
    if factor_degrees != [4]:
        raise ArithmeticError(f"branch quartic for {target['label']} is reducible over QQ")
    return {
        "target_label": target["label"],
        "source_public_point_index_one_based": target[
            "source_public_point_index_one_based"
        ],
        "public_point": target["public_point"],
        "short_point": [rational_text(target_x), rational_text(target_y)],
        "pencil_parameter_lambda": rational_text(pencil_parameter),
        "branch_polynomial_q_coefficients_low_to_high": polynomial_coefficients(q),
        "branch_polynomial_factor_degrees_over_Q": factor_degrees,
        "branch_polynomial_squarefree_degree": 4,
        "branch_polynomial_coprime_to_surface_discriminant": True,
        "branch_polynomial_coprime_to_trace_denominator": True,
        "geometrically_integral_genus_one_double_cover": True,
        "cover_rational_witness": {
            "t": rational_text(parameter),
            "s": rational_text(cover_coordinate),
        },
        "exact_target_x_and_y_verified": True,
    }


def build_payload():
    chord_module = load_script("r17_visibility_chords", CHORD_SCRIPT)
    model = json.loads(MODEL.read_text())
    section_data = json.loads(SECTIONS.read_text())
    pilot = json.loads(PILOT.read_text())
    relations = json.loads(RELATIONS.read_text())
    rational_experiment = json.loads(RATIONAL_EXPERIMENT.read_text())

    polynomial_ring = PolynomialRing(QQ, "t")
    function_field = polynomial_ring.fraction_field()
    A = polynomial_ring([QQ(value) for value in model["A_coefficients_low_to_high"]])
    B = polynomial_ring([QQ(value) for value in model["B_coefficients_low_to_high"]])
    discriminant = polynomial_ring(-16 * (4 * A**3 + 27 * B**2))
    if discriminant.degree() != 24:
        raise ArithmeticError("published R17 discriminant degree changed")
    basis_coordinates = chord_module.reconstruct_basis(
        polynomial_ring, A, B, section_data
    )
    generic_curve = EllipticCurve(function_field, [A, B])
    basis = [
        generic_curve(function_field(x_coordinate), function_field(y_coordinate))
        for x_coordinate, y_coordinate in basis_coordinates
    ]

    selected_trace = pilot["traces"][0]
    published_vector = vector(ZZ, selected_trace["published_basis_w"])
    tau = sum(
        (coefficient * point for coefficient, point in zip(published_vector, basis)),
        generic_curve(0),
    )
    X, Y = tau[0], tau[1]
    frame = chord_module.trace_chord_frame(X, Y, polynomial_ring)
    h, Nx, Ny, M0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
    if h.degree() != 2:
        raise ArithmeticError("selected common trace is not in the finite two-pole chart")

    relation_by_parameter = {row["parameter"]: row for row in relations["fibres"]}
    rational_by_parameter = {
        row["parameter"]: row for row in rational_experiment["fibres"]
    }
    literal_fibres = []
    rational_fibres = []
    for parameter_text, parameter in PARAMETERS.items():
        public_module = load_script(
            f"r17_visibility_public_{parameter_text.replace('/', '_').replace('-', 'm')}",
            PUBLIC_SCRIPTS[parameter_text],
        )
        specialized_curve = EllipticCurve(QQ, [A(parameter), B(parameter)])
        public_curve = EllipticCurve(
            QQ, list(public_module.GENERAL_WEIERSTRASS_COEFFICIENTS)
        )
        public_to_short = public_curve.isomorphism_to(specialized_curve)
        relation_fibre = relation_by_parameter[parameter_text]
        targets = []
        for record in relation_fibre["exceptional_basis"]:
            public_x, public_y = (QQ(value) for value in record["public_point"])
            point = public_to_short(public_curve(public_x, public_y))
            targets.append(
                {
                    "label": record["label"],
                    "source_public_point_index_one_based": record[
                        "source_public_point_index_one_based"
                    ],
                    "public_point": [rational_text(public_x), rational_text(public_y)],
                    "point": point,
                }
            )
        fitted = [
            fit_target(
                target=target,
                parameter=parameter,
                X=X,
                Y=Y,
                h=h,
                Nx=Nx,
                Ny=Ny,
                M0=M0,
                A=A,
                B=B,
                discriminant=discriminant,
                polynomial_ring=polynomial_ring,
                function_field=function_field,
                chord_module=chord_module,
            )
            for target in targets
        ]
        dimension = int(relation_fibre["exceptional_rank"])
        if len(fitted) != dimension:
            raise ArithmeticError("genus-one target count does not equal exceptional dimension")
        literal_fibres.append(
            {
                "parameter": parameter_text,
                "displayed_exceptional_dimension": dimension,
                "cumulative_visible_dimensions": [
                    {"maximum_degree": 1, "dimension": 0},
                    {"maximum_degree": 2, "dimension": dimension},
                    {"maximum_degree": 3, "dimension": dimension},
                    {"maximum_degree": 4, "dimension": dimension},
                ],
                "degree_two_fullness_reason": (
                    "Each ordered displayed exceptional basis point lies on its own "
                    "exact member of the common genus-one bisection pencil."
                ),
                "public_to_published_short_weierstrass_isomorphism_u_r_s_t": [
                    rational_text(value) for value in public_to_short.tuple()
                ],
                "targets": fitted,
            }
        )

        rational_row = rational_by_parameter[parameter_text]
        ranks = rational_row["cumulative_captured_exceptional_rank"]
        if [int(row["R_t(D)"]) for row in ranks] not in (
            [5, 5, 5], [3, 3, 3], [2, 2, 2], [1, 1, 1]
        ):
            raise ArithmeticError("imported rational-curve filtration changed")
        rational_fibres.append(
            {
                "parameter": parameter_text,
                "displayed_exceptional_dimension": dimension,
                "tested_cumulative_visible_dimensions": [
                    {
                        "maximum_degree": int(row["maximum_cover_degree"]),
                        "dimension": int(row["R_t(D)"]),
                    }
                    for row in ranks
                ],
                "completeness_by_degree": {
                    "2": "COMPLETE_RATIONAL_BISECTION_TRANSLATION_CLASSES",
                    "3": "INCOMPLETE_FRONTIER_SAMPLE_PLUS_COMPLETE_NORM_26_DEEP_SHELL",
                    "4": "INCOMPLETE_DETERMINISTIC_NORM_34_SAMPLE",
                },
            }
        )

    total_targets = sum(len(row["targets"]) for row in literal_fibres)
    if total_targets != 38:
        raise ArithmeticError(f"expected 38 displayed target fits, got {total_targets}")
    return {
        "schema": "elkies-k3.r17-multisection-visibility-filtration.v1",
        "status": "PASS_EXACT_LITERAL_FILTRATION_AND_FAIL_CLOSED_RATIONAL_EXPERIMENT",
        "definition": {
            "ambient": "(L_t/M_t) tensor_QQ, using the ordered displayed exceptional complement",
            "visible_space": (
                "The QQ-span of classes of rational points P on the fibre E_t that "
                "occur on a geometrically integral multisection C in the declared atlas."
            ),
            "splitting": (
                "C_t has a QQ-rational branch point; for a degree-two cover this is "
                "equivalent to the specialized cover value being a square."
            ),
            "translation_invariance": (
                "Translating C by a generic section changes P by an element of M_t, "
                "so the visible quotient class is unchanged."
            ),
        },
        "literal_all_genus_filtration": {
            "claim": (
                "At all four controls the displayed exceptional quotient is already "
                "visible in degree two through geometrically integral genus-one bisections."
            ),
            "common_trace_published_basis_w": [int(value) for value in published_vector],
            "common_trace_norm": 8,
            "trace_denominator_h_coefficients_low_to_high": polynomial_coefficients(h),
            "fibres": literal_fibres,
        },
        "rational_curve_filtration": {
            "claim": (
                "Only degree two is complete. The tested degree-three and degree-four "
                "layers add no direction, but their frontier universes are incomplete."
            ),
            "fibres": rational_fibres,
            "global_completeness_boundary": {
                "degree_3_norm_20_translation_cosets_total": 18024296,
                "degree_3_norm_20_translation_cosets_with_equations": 138,
                "degree_3_norm_20_translation_cosets_not_constructed": 18024158,
                "degree_3_norm_26_translation_cosets_complete": 320,
                "degree_4_total_translation_cosets": 4**17,
                "degree_4_sampled_cosets": 1025,
                "consequence": (
                    "The conjecture that rational multisections of degree at most four "
                    "span the displayed quotient remains untested, not falsified."
                ),
            },
        },
        "logical_outcome": {
            "unqualified_conjecture": "TRUE_ALREADY_AT_DEGREE_2_ON_DISPLAYED_QUOTIENTS",
            "rational_multisection_conjecture": "OPEN_AFTER_BOUNDED_DEGREE_3_AND_4_TESTS",
            "reason_for_separation": (
                "Allowing moving positive-genus pencils makes target incidence automatic; "
                "restricting to irreducible rational curves restores a finite lattice atlas."
            ),
        },
        "generation": {
            "command": (
                "sage -python elkies-k3/scripts/"
                "certify_r17_multisection_visibility_filtration.sage"
            ),
            "checker_sha256": digest(Path(__file__)),
            "inputs": {
                relative(path): digest(path)
                for path in (
                    MODEL,
                    SECTIONS,
                    PILOT,
                    RELATIONS,
                    RATIONAL_EXPERIMENT,
                    CHORD_SCRIPT,
                    *PUBLIC_SCRIPTS.values(),
                )
            },
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build_payload()
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.is_file() or arguments.output.read_text() != rendered:
            raise SystemExit("stale R17 multisection-visibility filtration certificate")
        terminal = "PASS"
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered)
        terminal = "WROTE"
    print(
        "R17VISIBILITYFILTRATION|targets=38|literal_d2=8,9,10,11|"
        "rational_tested_d4=5,3,2,1|status={}".format(terminal)
    )


if __name__ == "__main__":
    main()
