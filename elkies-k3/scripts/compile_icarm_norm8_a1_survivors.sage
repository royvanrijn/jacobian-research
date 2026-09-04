#!/usr/bin/env sage-python
"""Compile exact A1/MW16 pencils from target-screen survivor certificates.

For each QQ-isomorphic survivor this program rebuilds the residual-chord
quartic over QQ, writes its Jacobian Weierstrass equation, proves the singular
fibre configuration, and computes the saturated MW16 height lattice directly
from the rational Neron--Severi marking.  No section basis from another
fibration is assumed to transfer.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader
import json
from pathlib import Path
import sys

from sage.all import (
    EllipticCurve,
    PolynomialRing,
    QQ,
    ZZ,
    block_diagonal_matrix,
    gcd,
    matrix,
    pari,
    vector,
)
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
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
TARGET_SNAPSHOT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-database-sweep-v2.json"
)
CHORD = SCRIPTS / "construct_elkies_2026_bisections.sage"
SCREEN = SCRIPTS / "screen_icarm_norm8_a1_fibrations.sage"
FROZEN_SCREEN = SCRIPTS / "screen_icarm_curve398_norm8_a1_fibrations.sage"
CURVE398_COMPILER = SCRIPTS / "compile_icarm_curve398_hidden_a1_mw16.sage"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def load_source(name: str, path: Path):
    loader = SourceFileLoader(name, str(path))
    spec = spec_from_loader(name, loader)
    if spec is None:
        raise ImportError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


def qtext(value) -> str:
    value = QQ(value)
    return (
        str(value.numerator())
        if value.denominator() == 1
        else f"{value.numerator()}/{value.denominator()}"
    )


def poly_record(poly):
    return [qtext(poly[index]) for index in range(poly.degree() + 1)] if poly else ["0"]


def matrix_record(value):
    return [[qtext(entry) for entry in row] for row in value.rows()]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--curve-id", type=int, required=True)
    parser.add_argument("--exact-survivors", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    exact_path = args.exact_survivors.resolve()
    output_path = args.output.resolve()
    exact = json.loads(exact_path.read_text())
    exact_target = exact.get("target", {})
    exact_curve_id = exact_target.get("curve_id")
    legacy_label_matches = exact_target.get("label") == f"ICARM curve {args.curve_id}"
    if (
        exact_curve_id is not None
        and int(exact_curve_id) != args.curve_id
    ) or (exact_curve_id is None and not legacy_label_matches):
        raise ValueError("exact survivor certificate belongs to another target")
    accepted_exact_statuses = {
        "PASS_EXACT_QQ_ISOMORPHIC_A1_MW16_CANDIDATES",
        "PASS_EXACT_RATIONAL_PARAMETER_CANDIDATES",
    }
    if exact.get("status") not in accepted_exact_statuses:
        raise ValueError("exact survivor certificate has no QQ-isomorphic candidate")

    screen = load_source("icarm_a1_compile_screen", SCREEN)
    helpers = load_source("icarm_a1_compile_helpers", FROZEN_SCREEN)
    chord = load_source("icarm_a1_compile_chord", CHORD)
    lattice_helpers = load_source("icarm_a1_compile_lattice", CURVE398_COMPILER)

    model = json.loads(MODEL.read_text())
    rows = helpers.load_rows(TABLE)
    snapshot = json.loads(TARGET_SNAPSHOT.read_text())
    target_record = screen.target_record(snapshot, args.curve_id)
    target_ainvs, target_c4, target_c6, target_delta = screen.invariants(
        target_record["ainvs"]
    )
    target_a, target_b = -27 * QQ(str(target_c4)), -54 * QQ(str(target_c6))
    target_short = EllipticCurve(QQ, [target_a, target_b])

    old_ring = PolynomialRing(QQ, "t")
    old_field = old_ring.fraction_field()
    old_a = old_ring(model["weierstrass_model"]["A_coefficients_low_to_high"])
    old_b = old_ring(model["weierstrass_model"]["B_coefficients_low_to_high"])
    old_curve = EllipticCurve(old_field, [old_a, old_b])
    old_basis = tuple(
        old_curve(
            helpers.polynomial_from_record(record["X"], old_ring, QQ),
            helpers.polynomial_from_record(record["Y"], old_ring, QQ),
        )
        for record in model["sections"]["records"]
    )
    height_gram = matrix(ZZ, model["sections"]["height_gram"])
    ns_gram = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -height_gram)
    old_zero_class = vector(ZZ, [-1, 1] + [0] * 17)

    compiled = []
    for exact_record in exact["records"]:
        qq_specializations = [
            specialization
            for specialization in exact_record.get("specializations", [])
            if specialization.get("isomorphic_to_target_over_Q")
            or specialization.get("isomorphic_to_curve398_over_Q")
        ]
        if not qq_specializations:
            continue
        priority_rank = int(exact_record["priority_rank"])
        row = rows[priority_rank - 1]
        trace_vector = vector(ZZ, helpers.parse_vector(row["section_basis_w"]))
        if trace_vector * height_gram * trace_vector != 8:
            raise ArithmeticError(f"rank {priority_rank}: trace lost norm eight")
        trace = sum(
            (
                coefficient * point
                for coefficient, point in zip(trace_vector, old_basis)
                if coefficient
            ),
            old_curve(0),
        )
        frame = chord.trace_chord_frame(trace[0], trace[1], old_ring)
        h, nx, ny, m0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
        if h.degree() != 2:
            raise ArithmeticError(f"rank {priority_rank}: exact finite-pole chart unavailable")

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
            raise ArithmeticError(f"rank {priority_rank}: residual chord is not quartic")
        invariant_i, invariant_j = helpers.binary_quartic_invariants(quartic, parameter_ring)
        child_a, child_b = -27 * invariant_i, -27 * invariant_j
        child_delta = parameter_ring(-16 * (4 * child_a**3 + 27 * child_b**2))
        degree_profile = [int(child_a.degree()), int(child_b.degree()), int(child_delta.degree())]
        infinity_orders = [8 - degree_profile[0], 12 - degree_profile[1], 24 - degree_profile[2]]
        finite_squarefree = child_delta.gcd(child_delta.derivative()).degree() == 0
        if degree_profile != [8, 12, 22] or infinity_orders != [0, 0, 2] or not finite_squarefree:
            raise ArithmeticError(
                f"rank {priority_rank}: survivor does not compile to I2(infinity)+22I1"
            )

        fibre = vector(ZZ, [2, 2] + list(trace_vector))
        trace_class = lattice_helpers.section_class(trace_vector, height_gram)
        if (
            fibre != old_zero_class + trace_class
            or fibre * ns_gram * fibre
            or fibre * ns_gram * old_zero_class
        ):
            raise ArithmeticError(f"rank {priority_rank}: D=O+P_w identity failed")
        degree_one = lattice_helpers.enumerate_degree_one_vectors(height_gram, trace_vector)
        if not degree_one:
            raise ArithmeticError(f"rank {priority_rank}: no section for the new pencil")
        new_zero_vector = next(
            (section for section in degree_one if section * height_gram * section == 4),
            degree_one[0],
        )
        new_zero_class = lattice_helpers.section_class(new_zero_vector, height_gram)
        if lattice_helpers.ns_intersection(fibre, new_zero_class, height_gram) != 1:
            raise ArithmeticError(f"rank {priority_rank}: selected new zero is not a section")
        mate = fibre + new_zero_class
        complement = matrix(
            ZZ,
            [list(fibre * ns_gram), list(mate * ns_gram)],
        ).right_kernel_matrix()
        transport = matrix(
            QQ,
            [list(fibre), list(mate)] + [list(row_value) for row_value in complement],
        )
        child_frame = -(complement * ns_gram * complement.transpose())
        minimum = pari(matrix(ZZ, child_frame)).qfminim(2).sage()
        root_count = int(minimum[0])
        root_columns = matrix(ZZ, minimum[2])
        if (
            abs(transport.det()) != 1
            or child_frame.det() != 948
            or root_count != 2
            or root_columns.ncols() != 1
        ):
            raise ArithmeticError(f"rank {priority_rank}: child is not a primitive A1 frame")
        root = vector(ZZ, root_columns.column(0))
        if root * child_frame * root != 2 or gcd(list(root)) != 1:
            raise ArithmeticError(f"rank {priority_rank}: A1 root is not primitive of norm two")
        smith, left, right = matrix(ZZ, [root]).smith_form()
        if smith[0, 0] != 1:
            raise ArithmeticError(f"rank {priority_rank}: primitive root Smith form failed")
        completion = right.inverse()
        if vector(ZZ, completion[0]) != root or abs(completion.det()) != 1:
            raise ArithmeticError(f"rank {priority_rank}: root completion is not unimodular")
        completed_gram = completion * child_frame * completion.transpose()
        cross = matrix(QQ, 16, 1, list(completed_gram[1:, 0]))
        mw_gram = completed_gram[1:, 1:] - cross * cross.transpose() / 2
        if mw_gram.rank() != 16 or mw_gram.det() != 474:
            raise ArithmeticError(f"rank {priority_rank}: saturated MW16 lattice changed")

        comparison = child_a**3 * target_b**2 - target_a**3 * child_b**2
        rational_roots = sorted(
            -factor[0] / factor[1]
            for factor, multiplicity in comparison.factor()
            if factor.degree() == 1
            for unused in range(multiplicity)
        )
        expected_roots = sorted(
            QQ(specialization["parameter"]) for specialization in qq_specializations
        )
        if rational_roots != expected_roots:
            raise ArithmeticError(f"rank {priority_rank}: exact survivor roots changed")
        specializations = []
        for parameter in rational_roots:
            child_curve = EllipticCurve(QQ, [child_a(parameter), child_b(parameter)])
            if not child_curve.is_isomorphic(target_short):
                raise ArithmeticError(f"rank {priority_rank}: target specialization is a twist")
            isomorphism = child_curve.isomorphism_to(target_short)
            specializations.append(
                {
                    "lambda": qtext(parameter),
                    "child_short_coefficients": [
                        qtext(child_curve.a4()),
                        qtext(child_curve.a6()),
                    ],
                    "child_to_target_short_isomorphism_u_r_s_t": [
                        qtext(value) for value in isomorphism.tuple()
                    ],
                    "isomorphic_to_target_over_Q": True,
                }
            )

        mw_marking_in_source_ns = completion[1:, :] * complement
        root_in_source_ns = root * complement
        compiled.append(
            {
                "priority_rank": priority_rank,
                "orbit_mask": int(row["orbit_mask"]),
                "orbit_hex": row["orbit_hex"],
                "trace_section_basis_w": list(map(int, trace_vector)),
                "divisor_class_in_U_plus_M_minus": list(map(int, fibre)),
                "divisor_identity": "D=(2,2,w)=O+P_w",
                "old_fibre_degree": 2,
                "old_zero_degree": 0,
                "complete_old_degree_one_section_count": len(degree_one),
                "new_zero_source_section_basis_coordinates": list(
                    map(int, new_zero_vector)
                ),
                "equation": {
                    "form": "Y^2=X^3+A(lambda)*X+B(lambda)",
                    "A_coefficients_low_to_high": poly_record(child_a),
                    "B_coefficients_low_to_high": poly_record(child_b),
                    "degrees_A_B_Delta": degree_profile,
                    "infinity_orders_c4_c6_Delta": infinity_orders,
                    "finite_discriminant_squarefree": finite_squarefree,
                    "fibre_configuration": "I2 at infinity + 22 I1",
                },
                "generic_mordell_weil": {
                    "rank": 16,
                    "saturated": True,
                    "height_gram": matrix_record(mw_gram),
                    "height_gram_determinant": qtext(mw_gram.det()),
                    "a1_root_child_frame_coordinates": list(map(int, root)),
                    "a1_root_source_ns_coordinates": list(map(int, root_in_source_ns)),
                    "quotient_basis_source_ns_coordinates": [
                        list(map(int, row_value)) for row_value in mw_marking_in_source_ns.rows()
                    ],
                    "child_frame_determinant": int(child_frame.det()),
                    "child_frame_norm_two_vector_count_signed": root_count,
                    "proof": (
                        "The rational rank-19 NS marking splits off a primitive U. Its positive "
                        "frame has determinant 948 and exactly one signed pair of norm-two roots. "
                        "A unimodular completion of that primitive A1 root gives the saturated "
                        "rank-16 quotient height lattice of determinant 474."
                    ),
                },
                "target_specializations": specializations,
            }
        )
        print(
            f"ICARMA1COMPILE|curve={args.curve_id}|rank={priority_rank}"
            f"|degree1={len(degree_one)}|fibres=I2+22I1|mw_rank=16|mw_det=474",
            flush=True,
        )

    if not compiled:
        raise ArithmeticError("no QQ-isomorphic exact survivor was compiled")
    payload = {
        "schema": "elkies-k3.icarm-norm8-a1-compiled-survivors.v1",
        "status": "PASS_EXACT_COMPILED_A1_MW16_FIBRATIONS",
        "target": {
            "curve_id": args.curve_id,
            "label": f"ICARM curve {args.curve_id}",
            "snapshot_rank_lower_bound": int(target_record["snapshot_rank_lower_bound"]),
            "generalized_weierstrass_coefficients": [str(value) for value in target_ainvs],
            "c4": str(target_c4),
            "c6": str(target_c6),
            "discriminant": str(target_delta),
        },
        "compiled_fibration_count": len(compiled),
        "fibrations": compiled,
        "proof_boundary": (
            "This certifies each displayed polynomial A1/MW16 pencil, its saturated generic "
            "MW height lattice from the rational NS marking, and its QQ-isomorphic target "
            "parameter. It does not yet identify explicit generic section equations inside "
            "the target's displayed rational-point subgroup or prove an exact target rank."
        ),
        "inputs": {
            relative(path): digest(path)
            for path in (
                Path(__file__).resolve(),
                MODEL,
                TABLE,
                TARGET_SNAPSHOT,
                exact_path,
                CHORD,
                SCREEN,
                FROZEN_SCREEN,
                CURVE398_COMPILER,
            )
        },
        "software": {
            "sage_version": SAGE_VERSION,
            "pari_version": ".".join(map(str, pari.version())),
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/compile_icarm_norm8_a1_survivors.sage "
            f"--curve-id {args.curve_id} --exact-survivors {relative(exact_path)} "
            f"--output {relative(output_path)} --check"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output_path.is_file() or output_path.read_text() != serialized:
            raise ArithmeticError("stored compiled A1/MW16 artifact differs from replay")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        f"ICARMA1COMPILE|curve={args.curve_id}|fibrations={len(compiled)}"
        f"|status={payload['status']}|output={relative(output_path)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
