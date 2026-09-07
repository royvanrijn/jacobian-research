#!/usr/bin/env python3
"""Crosswalk the declared Leech coordinate targets to rooted catalogue frames.

The Leech probe emits exact global-catalogue ``surface_key`` objects before
the Co1 embedding quotient.  Compare those keys literally with the current
rooted-Niemeier catalogue and record which rootless MW17 targets already have
rootful presentations.  This is a catalogue crosswalk, not an equation or
Co1-orbit certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
LEECH = (
    ROOT
    / "artifacts/generated-results/elkies-k3-leech-minimal-basis-coordinate-shell-v1.json"
)
CATALOGUE = (
    ROOT / "artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-leech-to-rooted-target-source-crosswalk-v1.json"
)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_key(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def counter_rows(counter: Counter[int]) -> dict[str, int]:
    return {str(key): value for key, value in sorted(counter.items())}


def root_component_count(root_type: str) -> int:
    if root_type == "0":
        return 0
    pieces = root_type.split("+")
    total = 0
    for piece in pieces:
        match = re.fullmatch(r"(?:(\d+))?[ADE]\d+", piece)
        assert match, piece
        total += int(match.group(1) or 1)
    return total


def frame_summary(frame: dict) -> dict:
    backend_ids = sorted(
        {row["backend_id"] for row in frame.get("provenance", [])}
    )
    return {
        "frame_id": frame["frame_id"],
        "mw_rank_for_rho_19": frame["mw_rank_for_rho_19"],
        "root_type": frame["root_type"],
        "root_rank": frame["root_rank"],
        "root_component_count": root_component_count(frame["root_type"]),
        "backend_ids": backend_ids,
        "equation_status": "NOT_ASSERTED_BY_LATTICE_CATALOGUE",
    }


def build(leech: dict, catalogue: dict) -> dict:
    assert leech["schema"] == (
        "elkies-k3.leech-minimal-basis-coordinate-shell.v1"
    )
    assert leech["status"] == (
        "PASS_EXACT_DECLARED_LEECH_COORDINATE_LANGUAGE_PRE_CO1"
    )
    assert catalogue["schema"] == "elkies-k3.rank7-auxiliary-catalogue.v1"
    assert catalogue["status"].startswith("PASS_")
    rooted_by_key = {
        exact_key(surface["surface_key"]): surface
        for surface in catalogue["surfaces"]
    }
    assert len(rooted_by_key) == len(catalogue["surfaces"])

    matches = []
    unmatched = []
    easiest_mw = Counter()
    maximum_mw = Counter()
    match_determinants = Counter()
    unmatched_determinants = Counter()
    for target in leech["preliminary_T_NS_surfaces"]:
        key = exact_key(target["surface_key"])
        rooted = rooted_by_key.get(key)
        if rooted is None:
            unmatched.append(
                {
                    "catalogue_surface_id_if_imported": target[
                        "catalogue_surface_id_if_imported"
                    ],
                    "leech_preliminary_surface_id": target[
                        "preliminary_surface_id"
                    ],
                    "determinant": target["determinant"],
                    "signed_basis_type_ids": target["signed_basis_type_ids"],
                    "coordinate_embedding_multiplicity": target[
                        "coordinate_embedding_multiplicity"
                    ],
                    "leech_rootless_norm_four_pair_distribution": target[
                        "rootless_complement_norm_four_pair_distribution"
                    ],
                    "leech_rootless_norm_four_pairs_minimum": target[
                        "rootless_complement_norm_four_pairs_minimum"
                    ],
                    "leech_rootless_norm_four_pairs_maximum": target[
                        "rootless_complement_norm_four_pairs_maximum"
                    ],
                }
            )
            unmatched_determinants[target["determinant"]] += 1
            continue

        assert rooted["surface_id"] == target["catalogue_surface_id_if_imported"]
        assert rooted["determinant"] == target["determinant"]
        summaries = [frame_summary(frame) for frame in rooted["frames"]]
        selected = min(
            summaries,
            key=lambda row: (
                row["mw_rank_for_rho_19"],
                row["root_component_count"],
                row["frame_id"],
            ),
        )
        minimum_rank = min(row["mw_rank_for_rho_19"] for row in summaries)
        maximum_rank = max(row["mw_rank_for_rho_19"] for row in summaries)
        easiest_mw[minimum_rank] += 1
        maximum_mw[maximum_rank] += 1
        match_determinants[target["determinant"]] += 1
        matches.append(
            {
                "surface_id": rooted["surface_id"],
                "leech_preliminary_surface_id": target[
                    "preliminary_surface_id"
                ],
                "determinant": target["determinant"],
                "leech_mw_rank_for_rho_19": 17,
                "leech_signed_basis_type_ids": target[
                    "signed_basis_type_ids"
                ],
                "leech_coordinate_embedding_multiplicity": target[
                    "coordinate_embedding_multiplicity"
                ],
                "leech_rootless_norm_four_pair_distribution": target[
                    "rootless_complement_norm_four_pair_distribution"
                ],
                "leech_rootless_norm_four_pairs_minimum": target[
                    "rootless_complement_norm_four_pairs_minimum"
                ],
                "leech_rootless_norm_four_pairs_maximum": target[
                    "rootless_complement_norm_four_pairs_maximum"
                ],
                "easiest_catalogued_rooted_frame": selected,
                "easiest_catalogued_rooted_mw_rank": minimum_rank,
                "maximum_catalogued_rooted_mw_rank": maximum_rank,
                "all_catalogued_rooted_backend_ids": sorted(
                    {
                        backend_id
                        for row in summaries
                        for backend_id in row["backend_ids"]
                    }
                ),
                "catalogued_rooted_frame_count": len(summaries),
            }
        )

    matches.sort(
        key=lambda row: (
            row["easiest_catalogued_rooted_mw_rank"],
            row["easiest_catalogued_rooted_frame"]["root_component_count"],
            row["determinant"],
            row["surface_id"],
        )
    )
    unmatched.sort(
        key=lambda row: (
            row["determinant"],
            row["catalogue_surface_id_if_imported"],
        )
    )
    assert len(matches) == 43
    assert len(unmatched) == 107
    assert easiest_mw == Counter({12: 25, 13: 11, 14: 4, 17: 3})
    return {
        "schema": "elkies-k3.leech-to-rooted-target-source-crosswalk.v1",
        "status": "PASS_EXACT_T_NS_KEY_CROSSWALK_PRE_CO1",
        "proof_scope": {
            "proved": (
                "All 150 exact T/NS keys in the declared Leech coordinate "
                "language are compared literally with all current rooted-catalogue "
                "surface keys. Every match records its easiest currently catalogued "
                "rooted frame under the lexicographic (MW rank, root-component "
                "count, frame id) rule."
            ),
            "not_proved": (
                "The Leech embeddings are not yet quotiented by Co1, the coordinate "
                "language is not determinant-band complete, and a rooted lattice "
                "frame is not by itself an explicit equation or a shortest physical "
                "neighbour corridor. The 107 unmatched keys are new relative to the "
                "current catalogue only, not globally new K3 surfaces."
            ),
        },
        "selection_policy": {
            "target": "rootless Leech frame of MW17 at rho=19",
            "source_ranking": [
                "minimum currently catalogued rooted MW rank",
                "minimum root-component support count",
                "frame id tie-break",
            ],
            "unknown_metrics_not_imputed": [
                "equation complexity",
                "field of definition",
                "physical neighbour corridor",
            ],
        },
        "accounting": {
            "leech_preliminary_T_NS_keys": len(
                leech["preliminary_T_NS_surfaces"]
            ),
            "rooted_catalogue_surface_keys": len(catalogue["surfaces"]),
            "exact_key_matches": len(matches),
            "unmatched_leech_keys": len(unmatched),
            "matched_easiest_rooted_mw_rank_distribution": counter_rows(
                easiest_mw
            ),
            "matched_maximum_rooted_mw_rank_distribution": counter_rows(
                maximum_mw
            ),
            "matched_determinant_distribution": counter_rows(
                match_determinants
            ),
            "unmatched_determinant_distribution": counter_rows(
                unmatched_determinants
            ),
        },
        "matches": matches,
        "unmatched_leech_keys": unmatched,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--leech", type=Path, default=LEECH)
    parser.add_argument("--catalogue", type=Path, default=CATALOGUE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    payload = build(
        json.loads(arguments.leech.read_text()),
        json.loads(arguments.catalogue.read_text()),
    )
    payload["inputs"] = {
        str(arguments.leech.resolve().relative_to(ROOT)): digest(arguments.leech),
        str(arguments.catalogue.resolve().relative_to(ROOT)): digest(
            arguments.catalogue
        ),
    }
    payload["reproduce"] = (
        "python3 "
        "elkies-k3/scripts/compare_leech_coordinate_targets_to_rooted_catalogue.py"
    )
    payload["check_command"] = payload["reproduce"] + " --check"
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("Leech/rooted target-source crosswalk is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "LEECHROOTED|targets={}|matches={}|unmatched={}|"
        "easy_mw={}|status=PASS_EXACT_PRE_CO1".format(
            payload["accounting"]["leech_preliminary_T_NS_keys"],
            payload["accounting"]["exact_key_matches"],
            payload["accounting"]["unmatched_leech_keys"],
            payload["accounting"][
                "matched_easiest_rooted_mw_rank_distribution"
            ],
        )
    )


if __name__ == "__main__":
    main()
