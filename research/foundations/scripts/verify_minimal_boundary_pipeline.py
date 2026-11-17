#!/usr/bin/env python3
"""Verify the operational minimal-boundary invariant pipeline.

This is an exact regression on finite normalization exports.  It proves that
the extractor is invariant under labels/order and that its chart step is a
coefficient calculation.  It does not compute a canonical normalization from
an arbitrary polynomial map and does not prove the MBP1 gateway conjecture.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.minimal_boundary import (  # noqa: E402
    ExtractionReport,
    Outcome,
    extract_minimal_boundary,
    mutate_record,
)
from jcsearch.minimal_boundary_examples import (  # noqa: E402
    cancellation_record,
    countermodels,
    quadratic_gauge_record,
    weighted_record,
)


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "minimal_boundary_pipeline.json"
)


def serialize_report(report: ExtractionReport) -> dict[str, object]:
    return {
        "record": report.record_name,
        "selected_prime": report.selected_prime,
        "puncture_rank": report.puncture_rank,
        "valuation_rows": [list(row) for row in report.valuation_rows],
        "mbpkg": report.mbpkg.value,
        "chart": (
            None
            if report.chart_certificate is None
            else {
                "mechanism": report.chart_certificate.mechanism,
                "jacobian_unit": str(
                    report.chart_certificate.core_jacobian_unit
                ),
                "extracted_mark": str(
                    report.chart_certificate.extracted_mark
                ),
                "target_shear": str(
                    report.chart_certificate.target_shear
                ),
                "mbp_chart": report.chart_certificate.mbp_chart,
            }
        ),
        "predicates": {
            result.predicate: {
                "outcome": result.outcome.value,
                "reason": result.reason,
                "certificate": dict(result.certificate),
            }
            for result in report.predicates
        },
    }


def outcome_vector(report: ExtractionReport) -> tuple[str, ...]:
    return tuple(result.outcome.value for result in report.predicates)


def renamed_and_reordered(record):
    """Erase suggestive labels and reverse every unordered finite list."""

    primes = tuple(
        replace(prime, label=f"prime-{index}")
        for index, prime in enumerate(reversed(record.primes))
    )
    ledger = tuple(
        replace(row, label=f"row-{index}")
        for index, row in enumerate(reversed(record.ledger))
    )
    return mutate_record(
        record,
        name="opaque-regression-input",
        primes=primes,
        ledger=ledger,
    )


def run_suite() -> dict[str, object]:
    weighted_parameters = tuple(range(3, 9))
    cancellation_parameters = (
        (1, 1),
        (2, 1),
        (3, 1),
        (1, 2),
        (2, 2),
        (1, 3),
    )
    quadratic_parameters = tuple(range(3, 9))

    weighted_reports = tuple(
        extract_minimal_boundary(weighted_record(degree))
        for degree in weighted_parameters
    )
    cancellation_reports = tuple(
        extract_minimal_boundary(cancellation_record(m, r))
        for m, r in cancellation_parameters
    )
    quadratic_reports = tuple(
        extract_minimal_boundary(quadratic_gauge_record(degree))
        for degree in quadratic_parameters
    )

    assert all(report.mbpkg == Outcome.PASS for report in weighted_reports)
    assert all(
        report.chart_certificate is not None
        and report.chart_certificate.mechanism == "positive-section"
        for report in weighted_reports
    )
    assert all(
        report.mbpkg == Outcome.PASS for report in cancellation_reports
    )
    assert all(
        report.chart_certificate is not None
        and report.chart_certificate.mechanism == "reciprocal-integral"
        for report in cancellation_reports
    )

    # Quadratic gauges pass the first seven finite predicates, and the
    # coefficient algorithm extracts S^2.  They fail MBP1's CS definition
    # because that definition lists only weighted/cancellation charts.
    for report in quadratic_reports:
        assert all(
            result.outcome == Outcome.PASS
            for result in report.predicates[:-1]
        )
        assert report.predicates[-1].predicate == "CS"
        assert report.predicates[-1].outcome == Outcome.FAIL
        assert report.chart_certificate is not None
        assert report.chart_certificate.mechanism == "quadratic-incidence"
        assert (
            sp.expand(
                report.chart_certificate.extracted_mark
                - sp.Symbol("S") ** 2
            )
            == 0
        )

    countermodel_reports = tuple(
        extract_minimal_boundary(record) for record in countermodels()
    )
    expected_failures = {
        "countermodel-chart-perturbation": {"CS"},
        "countermodel-imprimitive-conormal": {"PC"},
        "countermodel-contracted-residue": {"NC"},
        "countermodel-nonsaturated-link": {"SAT"},
        "countermodel-spectator-prime": {"SCB", "LC"},
    }
    for report in countermodel_reports:
        actual = {
            result.predicate
            for result in report.predicates
            if result.outcome == Outcome.FAIL
        }
        assert actual == expected_failures[report.record_name]

    # Labels, list order, and the record name are never classification input.
    blindness_inputs = (
        weighted_record(5),
        cancellation_record(2, 2),
        quadratic_gauge_record(6),
    )
    blindness_checks = []
    for record in blindness_inputs:
        original = extract_minimal_boundary(record)
        opaque = extract_minimal_boundary(renamed_and_reordered(record))
        assert outcome_vector(original) == outcome_vector(opaque)
        original_mechanism = (
            original.chart_certificate.mechanism
            if original.chart_certificate is not None
            else None
        )
        opaque_mechanism = (
            opaque.chart_certificate.mechanism
            if opaque.chart_certificate is not None
            else None
        )
        assert original_mechanism == opaque_mechanism
        blindness_checks.append(
            {
                "source": record.name,
                "predicate_vector": list(outcome_vector(original)),
                "mechanism": original_mechanism,
            }
        )

    return {
        "format": "minimal-boundary-invariant-pipeline-v1",
        "arithmetic": "exact characteristic-zero SymPy expressions",
        "input_boundary": (
            "finite exports of canonical normalization data; normalization "
            "is not computed from a bare polynomial map"
        ),
        "weighted": [serialize_report(report) for report in weighted_reports],
        "cancellation": [
            serialize_report(report) for report in cancellation_reports
        ],
        "quadratic_gauge": [
            serialize_report(report) for report in quadratic_reports
        ],
        "countermodels": [
            serialize_report(report) for report in countermodel_reports
        ],
        "formula_blindness": blindness_checks,
        "research_outcome": {
            "operational_predicates": True,
            "quadratic_gauge_mark_extracted": "S**2",
            "missing_implication_proved": False,
            "remaining_gap": (
                "derive the finite export and intrinsic quotient/conormal "
                "marking from an unmarked canonical normalization"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--write-artifact",
        action="store_true",
        help=f"write {ARTIFACT.relative_to(ROOT)}",
    )
    args = parser.parse_args()
    payload = run_suite()
    if args.write_artifact:
        ARTIFACT.parent.mkdir(parents=True, exist_ok=True)
        ARTIFACT.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(f"WROTE {ARTIFACT.relative_to(ROOT)}")
    print("PASS: eight finite MBP predicates run on invariant exports")
    print("PASS: weighted and cancellation grids straighten coefficientwise")
    print("PASS: quadratic gauges extract the primitive S^2 incidence mark")
    print("PASS: perturbed and spectator countermodels fail at their defect layer")
    print("STATUS: the bare-normalization extraction implication remains open")


if __name__ == "__main__":
    main()
