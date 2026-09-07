#!/usr/bin/env sage-python
"""Factor the characteristic-zero curve-398 parameter equation for norm-8 traces.

This is the exact second stage of ``screen_icarm_curve398_norm8_a1_fibrations``.
It accepts a JSON list of priority ranks (or a modular-screen ledger), rebuilds
each residual-chord quartic over QQ, and factors the cross-multiplied j-equation.
Only linear factors certify rational candidate parameters; all other factors are
retained so that modular survivors cannot be mistaken for rational parameters.
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

# SourceFileLoader is needed because the reusable implementations have .sage
# filenames rather than importable Python module names.
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
TARGET_SOURCE = ROOT / "elliptic-curves/cas/icarm_curve398.py"
CHORD_SOURCE = SCRIPTS / "construct_elkies_2026_bisections.sage"
SCREEN_SOURCE = SCRIPTS / "screen_icarm_curve398_norm8_a1_fibrations.sage"


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
    parser.add_argument("--candidate-ranks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    screen = load_source("curve398_norm8_screen_exact_helpers", SCREEN_SOURCE)
    chord = load_source("curve398_norm8_chord_exact", CHORD_SOURCE)
    target = load_source("curve398_norm8_target_exact", TARGET_SOURCE)

    model_path = arguments.model.resolve()
    table_path = arguments.table.resolve()
    candidates_path = arguments.candidate_ranks.resolve()
    output_path = arguments.output.resolve()
    model = json.loads(model_path.read_text())
    rows = screen.load_rows(table_path)
    candidate_document = json.loads(candidates_path.read_text())
    candidate_ranks = (
        candidate_document
        if isinstance(candidate_document, list)
        else candidate_document.get("search", {}).get("survivor_priority_ranks")
    )
    if not isinstance(candidate_ranks, list) or any(not isinstance(v, int) for v in candidate_ranks):
        raise ValueError("candidate input must be an integer list or modular-screen ledger")

    a1, a2, a3, a4, a6 = tuple(QQ(str(v)) for v in target.GENERAL_WEIERSTRASS_COEFFICIENTS)
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    c4 = b2**2 - 24 * b4
    c6 = -b2**3 + 36 * b2 * b4 - 216 * b6
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
            screen.polynomial_from_record(section["X"], u_ring, QQ),
            screen.polynomial_from_record(section["Y"], u_ring, QQ),
        )
        for section in model["sections"]["records"]
    )

    lambda_ring = PolynomialRing(QQ, "lambda")
    lambda_variable = lambda_ring.gen()
    bivariate_ring = PolynomialRing(lambda_ring, "u")
    records = []
    for rank in candidate_ranks:
        row = rows[rank - 1]
        coordinates = screen.parse_vector(row["section_basis_w"])
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
        invariant_i, invariant_j = screen.binary_quartic_invariants(quartic, lambda_ring)
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
                    "isomorphic_to_curve398_over_Q": isomorphic,
                    "quadratic_twist_parameter": (
                        None if twist_parameter is None else str(twist_parameter)
                    ),
                    "isomorphism_to_curve398_short": (
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
            f"CURVE398A1EXACT|rank={rank}|degrees="
            f"{[int(f.degree()) for f, _ in factorization]}|rational={len(rational_roots)}",
            flush=True,
        )

    hits = [record for record in records if record["status"] == "PASS_EXACT_RATIONAL_PARAMETER_CANDIDATE"]
    payload = {
        "schema": "elkies-k3.icarm-curve398-norm8-a1-exact-parameter-factorization.v1",
        "status": "PASS_EXACT_RATIONAL_PARAMETER_CANDIDATES" if hits else "PASS_EXACT_NO_RATIONAL_PARAMETER",
        "target": {
            "label": "ICARM curve 398",
            "generalized_weierstrass_coefficients": [str(v) for v in (a1, a2, a3, a4, a6)],
            "c4": str(c4),
            "c6": str(c6),
        },
        "records": records,
        "rational_candidate_count": sum(len(r.get("rational_parameter_candidates", [])) for r in hits),
        "proof_boundary": (
            "Linear factors certify rational parameters with the same j-invariant as curve 398. "
            "They do not yet certify a Q-isomorphism, the A1 fibre configuration, or section transport."
        ),
        "inputs": {
            screen.relative(path): screen.digest(path)
            for path in (model_path, table_path, candidates_path, TARGET_SOURCE, CHORD_SOURCE, SCREEN_SOURCE)
        },
        "software": {"sage_version": SAGE_VERSION},
        "reproducing_command": (
            "sage -python elkies-k3/scripts/certify_icarm_curve398_norm8_a1_survivors.sage "
            f"--model {screen.relative(model_path)} --table {screen.relative(table_path)} "
            f"--candidate-ranks {screen.relative(candidates_path)} --output {screen.relative(output_path)}"
        ),
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"CURVE398A1EXACT|records={len(records)}|rational_candidates={payload['rational_candidate_count']}"
        f"|status={payload['status']}|output={screen.relative(output_path)}",
        flush=True,
    )


if __name__ == "__main__":
    main()
