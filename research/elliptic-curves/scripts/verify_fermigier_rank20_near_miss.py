#!/usr/bin/env python3
"""Replay the pinned Fermigier rank-at-least-20 near miss."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys


PROGRAM_ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = PROGRAM_ROOT.parent
ARTIFACT = (
    REPOSITORY_ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic-curves"
    / "fermigier_rank20_near_miss_v1.json"
)
sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.fermigier_near_miss import (  # noqa: E402
    build_fermigier_rank20_manifest,
    canonical_ratpoints_output,
)
from ecsearch.fermigier_rank import (  # noqa: E402
    parse_ratpoints_output,
    section_and_point_cloud_differences,
    specialize_fermigier_rank_sections,
)
from ecsearch.rank_certification import (  # noqa: E402
    IndependenceCertificate,
    verify_independence_certificate,
)


def main() -> None:
    expected = json.loads(ARTIFACT.read_text())
    abscissas = tuple(
        Fraction(value) for value in expected["bounded_search"]["abscissas"]
    )
    raw_output = canonical_ratpoints_output(abscissas)
    actual = build_fermigier_rank20_manifest(
        raw_output,
        maximum_reduction_prime=expected["point_cloud"][
            "maximum_reduction_prime"
        ],
    )
    # PARI's version string is provenance, not mathematical output.  Permit an
    # independent replay under a newer PARI release while requiring every
    # model, conductor, discriminant, root number, point, and certificate byte
    # represented in the manifest to agree.
    replay_pari_version = actual["global_curve"]["version"]
    pinned_pari_version = expected["global_curve"]["version"]
    assert replay_pari_version.startswith("[") and pinned_pari_version.startswith("[")
    actual["global_curve"]["version"] = pinned_pari_version
    assert actual == expected, "pinned Fermigier rank-20 near miss is stale"

    specialization = specialize_fermigier_rank_sections(Fraction(28917, 20))
    searched = parse_ratpoints_output(specialization.quartic_model, raw_output)
    cloud = section_and_point_cloud_differences(specialization, searched)
    indices = tuple(expected["point_cloud"]["selected_indices"])
    certificate = IndependenceCertificate.from_json_object(
        expected["point_cloud"]["certificate"]
    )
    verify_independence_certificate(
        specialization.canonical_model,
        tuple(cloud[index] for index in indices),
        certificate,
    )
    assert len(indices) == 20
    assert expected["global_curve"]["below_strict_log_target"] is True
    assert expected["limitations"]["target_status"] == "near miss, not a solution"
    print(
        "PASS Fermigier near miss: u=28917/20 has 20 exactly independent "
        f"points and log(N)={expected['global_curve']['log_conductor'][:20]}... "
        "<182.72; no twenty-first point claimed; "
        f"PARI replay {replay_pari_version} (pinned {pinned_pari_version})"
    )


if __name__ == "__main__":
    main()
