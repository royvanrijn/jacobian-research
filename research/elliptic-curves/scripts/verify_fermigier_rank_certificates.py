#!/usr/bin/env python3
"""Replay Fermigier's generic-section and published E22 rank lower bounds."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import shutil
import subprocess
import sys


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
POINT_DATA = PROGRAM_ROOT / "data" / "fermigier_e22_points.json"
ARTIFACT = (
    REPOSITORY_ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic-curves"
    / "fermigier_rank_certificates_v1.json"
)
sys.path.insert(0, str(PROGRAM_ROOT))
sys.path.insert(0, str(PROGRAM_ROOT / "scripts"))

from ecsearch.fermigier import (  # noqa: E402
    FERMIGIER_E22_RECONSTRUCTION_SHIFT,
    FERMIGIER_REPORTED_PARAMETER,
    fermigier_canonical_coefficients,
    thirteenth_visible_point,
)
from ecsearch.fermigier_rank import (  # noqa: E402
    specialize_fermigier_rank_sections,
)
from ecsearch.rank_certification import (  # noqa: E402
    IndependenceCertificate,
    is_on_weierstrass_curve,
    verify_independence_certificate,
)
from run_fermigier_rank_certificates import build_manifest  # noqa: E402


def _gp_modular_rational(value: Fraction, prime: int) -> str:
    return f"Mod({value.numerator},{prime})/{value.denominator}"


def _append_gp_certificate_replay(
    program: list[str],
    coefficients: tuple[Fraction, ...],
    points: tuple[tuple[Fraction, Fraction], ...],
    certificate: IndependenceCertificate,
    label: str,
) -> None:
    witness_prime = certificate.torsion_witness_prime
    witness_model = ",".join(
        _gp_modular_rational(value, witness_prime) for value in coefficients
    )
    program.append(f"E=ellinit([{witness_model}]);")
    program.append(
        f'if(ellcard(E)!={certificate.torsion_witness_group_order},'
        f'error("{label} torsion witness"));'
    )
    for row_number, row in enumerate(certificate.rows, start=1):
        prime = row.prime
        model = ",".join(
            _gp_modular_rational(value, prime) for value in coefficients
        )
        program.append(f"E=ellinit([{model}]);")
        program.append(f"G=[Mod({row.generator[0]},{prime}),Mod({row.generator[1]},{prime})];")
        program.append(
            f'if(ellcard(E)!={row.group_order}||ellorder(E,G)!={row.group_order},'
            f'error("{label} row {row_number} group"));'
        )
        for point_number, (point, logarithm) in enumerate(
            zip(points, row.logs), start=1
        ):
            x_coordinate = _gp_modular_rational(point[0], prime)
            y_coordinate = _gp_modular_rational(point[1], prime)
            program.append(f"P=[{x_coordinate},{y_coordinate}];")
            program.append(
                f'if(!ellisoncurve(E,P)||ellmul(E,G,{logarithm})!=P,'
                f'error("{label} row {row_number} point {point_number}"));'
            )


def _independent_gp_replay(
    generic_model: tuple[Fraction, ...],
    generic_points: tuple[tuple[Fraction, Fraction], ...],
    generic_certificate: IndependenceCertificate,
    e22_model: tuple[Fraction, ...],
    e22_points: tuple[tuple[Fraction, Fraction], ...],
    e22_certificate: IndependenceCertificate,
) -> str:
    gp = shutil.which("gp")
    if gp is None:
        raise SystemExit("PARI/GP executable 'gp' is required")
    program = ["default(realprecision,80);"]
    _append_gp_certificate_replay(
        program,
        generic_model,
        generic_points,
        generic_certificate,
        "generic",
    )
    _append_gp_certificate_replay(
        program, e22_model, e22_points, e22_certificate, "e22"
    )
    program.append('print("OK ",version());')
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input="\n".join(program) + "\n",
        text=True,
        capture_output=True,
        check=True,
    )
    if "***" in completed.stdout + completed.stderr:
        raise AssertionError(completed.stdout + completed.stderr)
    output = completed.stdout.strip()
    assert output.startswith("OK ")
    return output.removeprefix("OK ")


def main() -> None:
    expected = json.loads(ARTIFACT.read_text())
    maximum_prime = expected["search_bound"]["maximum_reduction_prime"]
    actual = build_manifest(maximum_prime)
    assert actual == expected, "pinned Fermigier rank certificates are stale"

    point_data = json.loads(POINT_DATA.read_text())
    e22_model = tuple(map(Fraction, point_data["weierstrass_coefficients"]))
    e22_points = tuple(
        (Fraction(point[0]), Fraction(point[1])) for point in point_data["points"]
    )
    assert len(e22_points) == 22
    assert all(is_on_weierstrass_curve(e22_model, point) for point in e22_points)

    specialization = specialize_fermigier_rank_sections(
        FERMIGIER_REPORTED_PARAMETER
    )
    assert specialization.quartic_model.shift == FERMIGIER_E22_RECONSTRUCTION_SHIFT
    assert specialization.canonical_model == fermigier_canonical_coefficients(
        FERMIGIER_REPORTED_PARAMETER
    )
    assert len(set(specialization.quartic_points)) == 13
    assert len(set(specialization.canonical_points)) == 13
    assert thirteenth_visible_point(specialization.quartic_model) == (
        Fraction(-46964, 195),
        Fraction(3170976819397626546496, 164775),
    )

    generic_certificate = IndependenceCertificate.from_json_object(
        expected["generic_sections"]["certificate"]
    )
    e22_certificate = IndependenceCertificate.from_json_object(
        expected["published_e22_points"]["certificate"]
    )
    verify_independence_certificate(
        specialization.canonical_model,
        specialization.section_differences,
        generic_certificate,
    )
    verify_independence_certificate(e22_model, e22_points, e22_certificate)
    assert generic_certificate.relation_prime == 5
    assert len(generic_certificate.rows) == 12
    assert e22_certificate.relation_prime == 2
    assert e22_certificate.torsion_witness_prime == 31
    assert e22_certificate.torsion_witness_group_order == 41
    assert len(e22_certificate.rows) == 22

    gp_version = _independent_gp_replay(
        specialization.canonical_model,
        specialization.section_differences,
        generic_certificate,
        e22_model,
        e22_points,
        e22_certificate,
    )
    print(
        "PASS Fermigier rank certificates: 12 generic section differences and "
        f"22 published E22 points are exactly independent (PARI/GP {gp_version} "
        "cross-check); no rank upper bound claimed"
    )


if __name__ == "__main__":
    main()
