#!/usr/bin/env sage-python
"""Fit genus-one residual-chord bisections through the rank-28 targets.

For a norm-eight R17 trace section

    tau = (Nx/h^2, Ny/h^3),  deg(h)=2,

regular residual-chord slopes form the pencil

    m = (M0 + lambda*h^2)/h,
    M0*Nx + Ny == 0 (mod h^2),  deg(M0)<4.

The normalized chord discriminant is then a quartic.  This script completely
enumerates and equation-ranks the minimum-norm-eight trace frontier, selects a
short prefix, and solves the single exact incidence equation for ``lambda``
at each of the eleven public rank-28 complement points.  It verifies the
target point, its short-Weierstrass Kummer barcode, the cover, the lifted
section, and every branch/smoothness condition over QQ.

This is a target-fitted positive-control experiment.  It does not search for
a rank-32 specialization or assert that the same fixed pencil members split
at another parameter.

status: ACTIVE_PROOF
claim: exact norm-eight genus-one pencil and 11/11 rank-28 target certificate
inputs: published R17 model/sections, pinned lattice, exact rank-28 points
outputs: artifacts/generated-results/elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json
"""

from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
from pathlib import Path

from sage.all import EllipticCurve, PolynomialRing, QQ, ZZ, QuadraticForm, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
RANK_SCRIPT = ROOT / "elkies-k3/scripts/rank_elkies_2026_bisection_orbits.sage"
CHORD_SCRIPT = ROOT / "elkies-k3/scripts/construct_elkies_2026_bisections.sage"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
SHORT_COORDS = ROOT / "elkies-k3/data/lattice/short_vector_basis_coords.txt"
SHORT_GRAM = ROOT / "elkies-k3/data/lattice/short_vector_basis_gram.txt"
PUBLISHED_IDENTIFICATION = ROOT / "artifacts/generated-results/elkies-2026-published-r17-target.json"
RELATIONS = (
    ROOT
    / "artifacts/generated-results/elliptic-curves"
    / "elkies_2026_exceptional_specialization_relations_v1.json"
)
PUBLIC_POINTS = ROOT / "elliptic-curves/cas/elkies_rank28.py"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json"
)
PARAMETER = QQ(-9529) / QQ(5471)


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


def equation_ranked_norm_eight(rank_module, pari_stack_gb: int):
    pinned = rank_module.load_matrix(PINNED)
    short_coordinates = rank_module.load_matrix(SHORT_COORDS)
    short_gram = rank_module.load_matrix(SHORT_GRAM)
    assert short_gram == short_coordinates * pinned * short_coordinates.transpose()

    identification = json.loads(PUBLISHED_IDENTIFICATION.read_text())
    basis_change = matrix(
        ZZ, identification["pinned_identification"]["basis_change_matrix"]
    )
    assert (
        identification["pinned_identification"]["gram_identity_orientation"]
        == "M^T*Gpub*M=Gpinned"
    )
    short_to_published = short_coordinates * basis_change.transpose()
    section_data = json.loads(SECTIONS.read_text())
    direct_costs, closures = rank_module.section_input_costs(section_data)

    coefficients = [
        short_gram[row, row] // 2 if row == column else short_gram[row, column]
        for row in range(17)
        for column in range(row, 17)
    ]
    pari.allocatemem(pari_stack_gb * 1024**3)
    shells = QuadraticForm(ZZ, 17, coefficients).short_vector_list_up_to_length(5, True)
    lower_masks = {
        rank_module.parity_mask(value)
        for shell_index in (0, 2, 3)
        for value in shells[shell_index]
    }
    best = {}
    multiplicities = Counter()
    for raw_value in shells[4]:
        value = vector(ZZ, raw_value)
        orbit = rank_module.parity_mask(value)
        if orbit in lower_masks:
            continue
        multiplicities[orbit] += 1
        key, score, published, sign = rank_module.canonical_scored_orientation(
            value, short_to_published, direct_costs, closures
        )
        signed_short = sign * value
        candidate = (key, score, published, signed_short)
        if orbit not in best or key < best[orbit][0]:
            best[orbit] = candidate

    if len(best) != 63925:
        raise ArithmeticError(f"expected 63925 norm-eight frontier masks, got {len(best)}")
    ranked = []
    for equation_rank, (orbit, (_key, score, published, short)) in enumerate(
        sorted(best.items(), key=lambda item: (item[1][0], item[0])), start=1
    ):
        pinned_vector = short * short_coordinates
        assert pinned_vector * pinned * pinned_vector == 8
        ranked.append(
            {
                "equation_rank": equation_rank,
                "orbit_mask": orbit,
                "orbit_hex": f"0x{orbit:05x}",
                "published_basis_w": tuple(int(entry) for entry in published),
                "pinned_rank17_w": tuple(int(entry) for entry in pinned_vector),
                "short_basis_w": tuple(int(entry) for entry in short),
                "minimal_representative_count": multiplicities[orbit],
                **score,
            }
        )
    return ranked, {
        "norm_eight_short_representatives": len(shells[4]),
        "minimum_norm_eight_translation_classes": len(best),
        "lower_minimum_masks_excluded": len(lower_masks),
        "group_addition_histogram": {
            str(cost): count
            for cost, count in sorted(
                Counter(row["group_addition_upper_bound"] for row in ranked).items()
            )
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace-limit", type=int, default=1)
    parser.add_argument("--pari-stack-gb", type=int, default=2)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.trace_limit <= 0:
        parser.error("--trace-limit must be positive")

    rank_module = load_script("elkies_r17_rank_bisections", RANK_SCRIPT)
    chord_module = load_script("elkies_r17_construct_bisections", CHORD_SCRIPT)
    public_module = load_script("elkies_rank28_public", PUBLIC_POINTS)
    ranked, enumeration = equation_ranked_norm_eight(
        rank_module, arguments.pari_stack_gb
    )
    selected = ranked[: min(arguments.trace_limit, len(ranked))]

    model = json.loads(MODEL.read_text())
    section_data = json.loads(SECTIONS.read_text())
    relation_payload = json.loads(RELATIONS.read_text())
    target_fibre = next(
        row
        for row in relation_payload["fibres"]
        if row["parameter"] == "-9529/5471"
    )

    polynomial_ring = PolynomialRing(QQ, "t")
    function_field = polynomial_ring.fraction_field()
    A = polynomial_ring([QQ(value) for value in model["A_coefficients_low_to_high"]])
    B = polynomial_ring([QQ(value) for value in model["B_coefficients_low_to_high"]])
    discriminant = polynomial_ring(-16 * (4 * A**3 + 27 * B**2))
    assert discriminant.degree() == 24
    basis_coordinates = chord_module.reconstruct_basis(
        polynomial_ring, A, B, section_data
    )
    generic_curve = EllipticCurve(function_field, [A, B])
    basis = [
        generic_curve(function_field(x_coordinate), function_field(y_coordinate))
        for x_coordinate, y_coordinate in basis_coordinates
    ]

    specialized_curve = EllipticCurve(QQ, [A(PARAMETER), B(PARAMETER)])
    public_curve = EllipticCurve(
        QQ, list(public_module.GENERAL_WEIERSTRASS_COEFFICIENTS)
    )
    public_to_short = public_curve.isomorphism_to(specialized_curve)
    target_records = []
    for record in target_fibre["exceptional_basis"]:
        public_x, public_y = (QQ(value) for value in record["public_point"])
        public_point = public_curve(public_x, public_y)
        short_point = public_to_short(public_point)
        target_records.append(
            {
                "label": record["label"],
                "source_public_point_index_one_based": record[
                    "source_public_point_index_one_based"
                ],
                "public_point": [rational_text(public_x), rational_text(public_y)],
                "short_point": [
                    rational_text(short_point[0]),
                    rational_text(short_point[1]),
                ],
                "point": short_point,
            }
        )

    trace_records = []
    for ranked_trace in selected:
        published_vector = vector(ZZ, ranked_trace["published_basis_w"])
        tau = sum(
            (coefficient * point for coefficient, point in zip(published_vector, basis)),
            generic_curve(0),
        )
        if tau.is_zero():
            raise ArithmeticError("norm-eight trace unexpectedly vanished")
        X, Y = tau[0], tau[1]
        frame = chord_module.trace_chord_frame(X, Y, polynomial_ring)
        h, Nx, Ny, M0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
        if h.degree() != 2:
            raise ArithmeticError(
                "selected pilot trace has a pole at infinity; reciprocal-chart fitting "
                "is outside this finite-pole pilot"
            )

        fitted = []
        for target in target_records:
            target_point = target["point"]
            target_x, target_y = target_point[0], target_point[1]
            trace_x = QQ(X(PARAMETER))
            trace_y = QQ(Y(PARAMETER))
            h_value = QQ(h(PARAMETER))
            if not h_value or target_x == trace_x:
                raise ArithmeticError(
                    f"degenerate target incidence for {target['label']}"
                )
            target_slope = (target_y + trace_y) / (target_x - trace_x)
            pencil_parameter = (
                target_slope * h_value - QQ(M0(PARAMETER))
            ) / h_value**2
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
            if not q.is_irreducible():
                raise ArithmeticError(
                    f"branch quartic for {target['label']} is reducible"
                )
            if q.gcd(discriminant).degree() or q.gcd(h).degree():
                raise ArithmeticError(
                    f"branch quartic for {target['label']} meets a forbidden divisor"
                )

            sum_x = data["sum_x"]
            product_x = data["product_x"]
            cover_coordinate = (2 * target_x - QQ(sum_x(PARAMETER))) / h_value
            line_value = QQ(M(PARAMETER)) / h_value * target_x - (
                trace_y + QQ(M(PARAMETER)) / h_value * trace_x
            )
            if cover_coordinate**2 != QQ(q(PARAMETER)):
                raise ArithmeticError("target does not give a rational cover witness")
            if line_value != target_y:
                raise ArithmeticError("target y-coordinate misses the fitted line")
            if (
                target_x**2
                - QQ(sum_x(PARAMETER)) * target_x
                + QQ(product_x(PARAMETER))
            ):
                raise ArithmeticError("target x-coordinate misses the residual quadratic")
            if target_y**2 != target_x**3 + A(PARAMETER) * target_x + B(PARAMETER):
                raise ArithmeticError("short target point fails its Weierstrass equation")

            # Polynomial degree bounds prove that the lift is disjoint from
            # the zero section on the degree-two base change.
            if not (
                data["x0"].degree() <= 4
                and data["x1"].degree() <= 2
                and data["y0"].degree() <= 6
                and data["y1"].degree() <= 4
            ):
                raise ArithmeticError("lifted section exceeds the integral degree bounds")

            kummer_generator = polynomial_ring.base_ring()["theta"](
                [target_x, -1]
            )
            # The norm of x(Q)-theta in QQ[theta]/(theta^3+A0*theta+B0)
            # is y(Q)^2; literal x equality is the barcode comparison.
            kummer_norm = target_x**3 + A(PARAMETER) * target_x + B(PARAMETER)
            assert kummer_norm == target_y**2
            fitted.append(
                {
                    "target_label": target["label"],
                    "source_public_point_index_one_based": target[
                        "source_public_point_index_one_based"
                    ],
                    "public_point": target["public_point"],
                    "short_point": target["short_point"],
                    "pencil_parameter_lambda": rational_text(pencil_parameter),
                    "slope_numerator_M_coefficients_low_to_high": polynomial_coefficients(M),
                    "branch_polynomial_q_coefficients_low_to_high": polynomial_coefficients(q),
                    "branch_polynomial_degree": int(q.degree()),
                    "branch_polynomial_irreducible_over_Q": True,
                    "branch_polynomial_squarefree": True,
                    "branch_polynomial_coprime_to_surface_discriminant": True,
                    "branch_polynomial_coprime_to_trace_denominator": True,
                    "cover_rational_witness": {
                        "t": rational_text(PARAMETER),
                        "s": rational_text(cover_coordinate),
                    },
                    "residual_quadratic": {
                        "u_coefficients_low_to_high": polynomial_coefficients(-sum_x),
                        "v_coefficients_low_to_high": polynomial_coefficients(product_x),
                    },
                    "lifted_section": {
                        "cover": "s^2=q(t)",
                        "x0_coefficients_low_to_high": polynomial_coefficients(data["x0"]),
                        "x1_coefficients_low_to_high": polynomial_coefficients(data["x1"]),
                        "y0_coefficients_low_to_high": polynomial_coefficients(data["y0"]),
                        "y1_coefficients_low_to_high": polynomial_coefficients(data["y1"]),
                        "trace": "P(t,s)+P(t,-s)=tau(t)",
                    },
                    "kummer_barcode": {
                        "algebra": "QQ[theta]/(theta^3+A(t0)*theta+B(t0))",
                        "target_generator_coefficients": [
                            rational_text(kummer_generator[0]),
                            rational_text(kummer_generator[1]),
                            "0",
                        ],
                        "norm": rational_text(kummer_norm),
                        "norm_square_root": rational_text(target_y),
                        "comparison": "x(P(t0,s0))-theta = x(Q_i)-theta exactly",
                        "equal_modulo_squares": True,
                    },
                    "exact_target_x_and_y_verified": True,
                }
            )

        trace_records.append(
            {
                key: value
                for key, value in ranked_trace.items()
                if key not in {"_sort_key"}
            }
            | {
                "trace_section": {
                    "expression": "sum published_basis_w[j]*P_j",
                    "h_coefficients_low_to_high": polynomial_coefficients(h),
                    "Nx_coefficients_low_to_high": polynomial_coefficients(Nx),
                    "Ny_coefficients_low_to_high": polynomial_coefficients(Ny),
                    "least_slope_M0_coefficients_low_to_high": polynomial_coefficients(M0),
                    "pole_pattern": "two finite simple poles counted by h(t)",
                },
                "coefficient_template": "M(t)=M0(t)+lambda*h(t)^2",
                "targets": fitted,
            }
        )

    success_count = sum(len(record["targets"]) for record in trace_records)
    expected_success_count = len(trace_records) * len(target_records)
    if success_count != expected_success_count:
        raise ArithmeticError("not every selected trace-target pair succeeded")
    payload = {
        "schema": "elkies-k3.r17-rank28-genus-one-bisection-pilot.v1",
        "status": "PASS_EXACT_R17_RANK28_GENUS_ONE_BISECTION_PILOT",
        "parameter": "-9529/5471",
        "inputs": {
            relative(path): digest(path)
            for path in (
                MODEL,
                SECTIONS,
                PINNED,
                SHORT_COORDS,
                SHORT_GRAM,
                PUBLISHED_IDENTIFICATION,
                RELATIONS,
                PUBLIC_POINTS,
                RANK_SCRIPT,
                CHORD_SCRIPT,
                Path(__file__),
            )
        },
        "public_to_published_short_weierstrass_isomorphism_u_r_s_t": [
            rational_text(value) for value in public_to_short.tuple()
        ],
        "complete_norm_eight_frontier": enumeration,
        "selected_trace_count": len(trace_records),
        "target_count": len(target_records),
        "successful_trace_target_pairs": success_count,
        "common_template_group": {
            "trace_shell": "R17 norm 8",
            "arithmetic_genus": 1,
            "degree_over_published_t_line": 2,
            "pole_pattern": "deg(h)=2 with both poles finite",
            "coefficient_template": "M=M0+lambda*h^2",
            "branch_polynomial_degree": 4,
            "targets_succeeded": [record["label"] for record in target_records],
        },
        "generic_cover_consequence": {
            "base_curve_genus": 1,
            "base_curve_has_displayed_QQ_point": True,
            "base_change_degree": 2,
            "base_change_chi": 4,
            "lift_disjoint_from_zero": True,
            "lift_self_intersection": -4,
            "conjugate_lift_intersection": 4,
            "anti_invariant_height": 16,
            "generic_mw_rank_lower_bound_on_each_cover": 18,
            "independent_from_invariant_R17": True,
            "proof": (
                "The four simple branch points lie on smooth source fibres and avoid h=0, "
                "so the two conjugate integral lifts meet transversely exactly four times. "
                "Their difference has height 2*(4-(-4))=16 and is anti-invariant under "
                "the deck involution, hence independent of the invariant R17 subgroup."
            ),
        },
        "interpretation": (
            "One norm-eight divisor class gives a genus-one pencil. Fitting lambda through "
            "each Q_i is therefore the exact selection of the pencil member containing that "
            "point; the 11/11 incidence itself is geometric, while irreducible squarefree "
            "quartic branching over smooth fibres and the height-16 anti-invariant section "
            "are the nontrivial equation-level checks."
        ),
        "proof_boundary": (
            "This exact positive-control certificate uses the equation-cheapest finite-pole "
            "norm-eight trace only. It does not classify all 63,925 frontier pencils, include "
            "the 43 norm-twelve deep classes, prove a rank-32 specialization, or show that "
            "several fixed rational pencil members split together away from the fitted control."
        ),
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "search_elkies_2026_rank28_genus_one_bisections.sage --check"
        ),
        "traces": trace_records,
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.is_file():
            raise FileNotFoundError(f"missing checked output {arguments.output}")
        if arguments.output.read_text() != serialized:
            raise ArithmeticError(f"checked output differs from replay: {arguments.output}")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(serialized)
    print(
        "R17GENUS1BISECTION|traces={}|targets={}|successes={}|status={}|output={}".format(
            len(trace_records),
            len(target_records),
            success_count,
            payload["status"],
            arguments.output,
        )
    )


if __name__ == "__main__":
    main()
