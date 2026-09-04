#!/usr/bin/env sage-python
"""Factor exact QQ j-equations for survivors of one ICARM A1/MW16 screen.

The input is either a JSON integer list or the compact/checkpoint ledger from
``screen_icarm_norm8_a1_fibrations.sage``.  Every selected norm-eight trace is
rebuilt over QQ.  Linear factors give all rational j-matching parameters, and
each resulting specialization is tested for QQ-isomorphism to the pinned
target rather than accepted merely up to quadratic twist.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from sage.all import EllipticCurve, PolynomialRing, QQ
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader


DEFAULT_MODEL = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
)
DEFAULT_TABLE = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
)
DEFAULT_TARGET_SNAPSHOT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-r17-norm12-icarm-database-sweep-v2.json"
)
CHORD_SOURCE = SCRIPTS / "construct_elkies_2026_bisections.sage"
SCREEN_SOURCE = SCRIPTS / "screen_icarm_norm8_a1_fibrations.sage"
FROZEN_SCREEN_SOURCE = SCRIPTS / "screen_icarm_curve398_norm8_a1_fibrations.sage"


def load_source(name: str, path: Path):
    loader = SourceFileLoader(name, str(path))
    spec = spec_from_loader(name, loader)
    if spec is None:
        raise ImportError(f"cannot load {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    loader.exec_module(module)
    return module


def rational_roots_from_factors(factorization):
    roots = []
    for factor, multiplicity in factorization:
        if factor.degree() == 1:
            roots.append(
                {
                    "parameter": str(-factor[0] / factor[1]),
                    "multiplicity": int(multiplicity),
                }
            )
    return roots


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    parser.add_argument("--table", type=Path, default=DEFAULT_TABLE)
    parser.add_argument("--target-snapshot", type=Path, default=DEFAULT_TARGET_SNAPSHOT)
    parser.add_argument("--curve-id", type=int, required=True)
    parser.add_argument("--candidate-ranks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    screen = load_source("icarm_norm8_a1_exact_screen", SCREEN_SOURCE)
    helpers = load_source("icarm_norm8_a1_exact_helpers", FROZEN_SCREEN_SOURCE)
    chord = load_source("icarm_norm8_a1_exact_chord", CHORD_SOURCE)

    model_path = arguments.model.resolve()
    table_path = arguments.table.resolve()
    snapshot_path = arguments.target_snapshot.resolve()
    candidates_path = arguments.candidate_ranks.resolve()
    output_path = arguments.output.resolve()
    model = json.loads(model_path.read_text())
    rows = helpers.load_rows(table_path)
    candidate_document = json.loads(candidates_path.read_text())
    candidate_ranks = (
        candidate_document
        if isinstance(candidate_document, list)
        else candidate_document.get("search", {}).get("survivor_priority_ranks")
    )
    if not isinstance(candidate_ranks, list) or any(not isinstance(v, int) for v in candidate_ranks):
        raise ValueError("candidate input has no integer survivor_priority_ranks list")
    if len(set(candidate_ranks)) != len(candidate_ranks):
        raise ValueError("candidate rank list contains duplicates")
    if any(rank < 1 or rank > len(rows) for rank in candidate_ranks):
        raise ValueError("candidate rank lies outside the complete priority table")
    if isinstance(candidate_document, dict):
        candidate_target = candidate_document.get("target")
        prior_curve_id = candidate_target.get("curve_id") if isinstance(candidate_target, dict) else None
        if prior_curve_id is not None and int(prior_curve_id) != arguments.curve_id:
            raise ValueError("candidate ledger belongs to a different ICARM curve")

    snapshot = json.loads(snapshot_path.read_text())
    pinned_target = screen.target_record(snapshot, arguments.curve_id)
    ainvs_fraction, c4_fraction, c6_fraction, delta_fraction = screen.invariants(
        pinned_target["ainvs"]
    )
    a1, a2, a3, a4, a6 = tuple(QQ(str(value)) for value in ainvs_fraction)
    c4, c6 = QQ(str(c4_fraction)), QQ(str(c6_fraction))
    target_a = -27 * c4
    target_b = -54 * c6
    target_short_curve = EllipticCurve(QQ, [target_a, target_b])

    u_ring = PolynomialRing(QQ, "u")
    function_field = u_ring.fraction_field()
    old_a = u_ring(model["weierstrass_model"]["A_coefficients_low_to_high"])
    old_b = u_ring(model["weierstrass_model"]["B_coefficients_low_to_high"])
    old_curve = EllipticCurve(function_field, [old_a, old_b])
    basis = tuple(
        old_curve(
            helpers.polynomial_from_record(section["X"], u_ring, QQ),
            helpers.polynomial_from_record(section["Y"], u_ring, QQ),
        )
        for section in model["sections"]["records"]
    )

    lambda_ring = PolynomialRing(QQ, "lambda")
    lambda_variable = lambda_ring.gen()
    bivariate_ring = PolynomialRing(lambda_ring, "u")
    records = []
    for rank in candidate_ranks:
        row = rows[rank - 1]
        coordinates = helpers.parse_vector(row["section_basis_w"])
        record = {
            "priority_rank": rank,
            "orbit_mask": int(row["orbit_mask"]),
            "orbit_hex": row["orbit_hex"],
            "section_basis_w": list(coordinates),
        }
        trace = sum(
            (coefficient * point for coefficient, point in zip(coordinates, basis) if coefficient),
            old_curve(0),
        )
        frame = chord.trace_chord_frame(trace[0], trace[1], u_ring)
        h, nx, ny, m0 = (frame[key] for key in ("h", "Nx", "Ny", "M0"))
        record["finite_pole_degree"] = int(h.degree())
        if h.degree() != 2:
            record["status"] = "UNKNOWN_POLE_AT_INFINITY_CHART"
            records.append(record)
            continue
        h_symbolic = bivariate_ring(h)
        nx_symbolic = bivariate_ring(nx)
        ny_symbolic = bivariate_ring(ny)
        m_symbolic = bivariate_ring(m0) + lambda_variable * h_symbolic**2
        numerator = (
            m_symbolic**4
            - 6 * m_symbolic**2 * nx_symbolic
            - 8 * m_symbolic * ny_symbolic
            - 3 * nx_symbolic**2
            - 4 * bivariate_ring(old_a) * h_symbolic**4
        )
        quartic, remainder = numerator.quo_rem(h_symbolic**6)
        if remainder or quartic.degree() != 4:
            raise ArithmeticError(f"rank {rank}: residual chord is not quartic")
        invariant_i, invariant_j = helpers.binary_quartic_invariants(quartic, lambda_ring)
        child_a = -27 * invariant_i
        child_b = -27 * invariant_j
        comparison = child_a**3 * target_b**2 - target_a**3 * child_b**2
        factorization = comparison.factor()
        rational_roots = rational_roots_from_factors(factorization)
        specializations = []
        for root_record in rational_roots:
            parameter = QQ(root_record["parameter"])
            specialized_curve = EllipticCurve(
                QQ,
                [QQ(child_a(parameter)), QQ(child_b(parameter))],
            )
            isomorphic = bool(specialized_curve.is_isomorphic(target_short_curve))
            try:
                twist_parameter = specialized_curve.is_quadratic_twist(target_short_curve)
            except (AttributeError, NotImplementedError):
                twist_parameter = None
            specializations.append(
                {
                    "parameter": str(parameter),
                    "child_short_coefficients": [
                        str(specialized_curve.a4()),
                        str(specialized_curve.a6()),
                    ],
                    "isomorphic_to_target_over_Q": isomorphic,
                    "quadratic_twist_parameter": (
                        None if twist_parameter is None else str(twist_parameter)
                    ),
                    "isomorphism_to_target_short": (
                        str(specialized_curve.isomorphism_to(target_short_curve))
                        if isomorphic
                        else None
                    ),
                }
            )
        record.update(
            {
                "status": (
                    "PASS_EXACT_RATIONAL_PARAMETER_CANDIDATE"
                    if rational_roots
                    else "PASS_EXACT_NO_RATIONAL_PARAMETER"
                ),
                "comparison_degree": int(comparison.degree()),
                "factorization": [
                    {
                        "degree": int(factor.degree()),
                        "multiplicity": int(multiplicity),
                        "coefficients_low_to_high": [str(v) for v in factor],
                    }
                    for factor, multiplicity in factorization
                ],
                "rational_parameter_candidates": rational_roots,
                "specializations": specializations,
            }
        )
        records.append(record)
        print(
            f"ICARMA1EXACT|curve={arguments.curve_id}|rank={rank}|degrees="
            f"{[int(f.degree()) for f, _ in factorization]}|rational={len(rational_roots)}"
            f"|qq_isomorphic={sum(item['isomorphic_to_target_over_Q'] for item in specializations)}",
            flush=True,
        )

    hits = [
        record
        for record in records
        if record["status"] == "PASS_EXACT_RATIONAL_PARAMETER_CANDIDATE"
    ]
    qq_hits = [
        (record, specialization)
        for record in hits
        for specialization in record["specializations"]
        if specialization["isomorphic_to_target_over_Q"]
    ]
    unknown = [record for record in records if record["status"].startswith("UNKNOWN_")]
    if qq_hits:
        status = "PASS_EXACT_QQ_ISOMORPHIC_A1_MW16_CANDIDATES"
    elif hits:
        status = "PASS_EXACT_J_MATCHES_BUT_NO_QQ_ISOMORPHIC_CANDIDATE"
    elif unknown:
        status = "UNKNOWN_UNCOMPILED_POLE_AT_INFINITY_SURVIVOR"
    else:
        status = "PASS_EXACT_NO_RATIONAL_PARAMETER_IN_COMPLETE_LAYER"
    payload = {
        "schema": "elkies-k3.icarm-norm8-a1-exact-parameter-factorization.v1",
        "status": status,
        "target": {
            "curve_id": arguments.curve_id,
            "label": f"ICARM curve {arguments.curve_id}",
            "snapshot_rank_lower_bound": int(pinned_target["snapshot_rank_lower_bound"]),
            "generalized_weierstrass_coefficients": [str(value) for value in ainvs_fraction],
            "c4": str(c4_fraction),
            "c6": str(c6_fraction),
            "discriminant": str(delta_fraction),
        },
        "screen_survivor_count": len(candidate_ranks),
        "records": records,
        "rational_j_candidate_count": sum(
            len(record.get("rational_parameter_candidates", [])) for record in hits
        ),
        "qq_isomorphic_candidate_count": len(qq_hits),
        "proof_boundary": (
            "Exact factorization decides every supplied finite-pole survivor. Linear factors "
            "certify rational j-matching parameters; QQ-isomorphism is checked separately. "
            "A complete no-hit applies only to the committed 63,917-class norm-eight layer "
            "on source chart 11952. Any pole-at-infinity survivor remains UNKNOWN."
        ),
        "inputs": {
            screen.relative(path): screen.digest(path)
            for path in (
                Path(__file__).resolve(),
                model_path,
                table_path,
                snapshot_path,
                candidates_path,
                CHORD_SOURCE,
                SCREEN_SOURCE,
                FROZEN_SCREEN_SOURCE,
            )
        },
        "software": {"sage_version": SAGE_VERSION},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/certify_icarm_norm8_a1_survivors.sage "
            f"--curve-id {arguments.curve_id} --model {screen.relative(model_path)} "
            f"--table {screen.relative(table_path)} --target-snapshot {screen.relative(snapshot_path)} "
            f"--candidate-ranks {screen.relative(candidates_path)} "
            f"--output {screen.relative(output_path)}"
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"ICARMA1EXACT|curve={arguments.curve_id}|records={len(records)}"
        f"|rational_j={payload['rational_j_candidate_count']}"
        f"|qq_isomorphic={payload['qq_isomorphic_candidate_count']}"
        f"|status={payload['status']}|output={screen.relative(output_path)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
