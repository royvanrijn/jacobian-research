#!/usr/bin/env sage
"""Build the surface-first rank-seven auxiliary catalogue.

This is the catalogue/merge layer, not an embedding enumerator.  It imports
exact primitive embeddings from independently replayable backends, groups
first by the pair (T, NS), and only then groups partner auxiliaries and frame
lattices.  Backend and determinant-band coverage remain explicit so that a
bounded discovery artifact can never masquerade as the requested all-orbit
census.

status: ACTIVE_SEARCH_INFRASTRUCTURE
claim: exact (T,NS)-first deduplication of the imported one-root and 2C-fixed
  2A7+2D5 shells, the 6A4 double-swap, 4A5+D4 order-four, and
  2A9+D6/3A8/3D8/4A6/4E6/6D4/8A3/12A2 all-residual-class coordinate shells, and the certified
  24A1 Golay-octad design and complete positive seven-octad subfamily.
inputs: artifacts/generated-results/elkies-k3-lattice-foundry-v1.json,
  artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json,
  artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json,
  artifacts/generated-results/elkies-k3-24a1-octad-prefix-orbits-v1.json,
  artifacts/generated-results/elkies-k3-24a1-octad-completion-manifest-v1.json,
  artifacts/generated-results/elkies-k3-24a1-weyl-m24-canonicalization-00000-10547-v2.json,
  artifacts/generated-results/elkies-k3-cross-niemeier-mod2-priority-v1.json,
  artifacts/generated-results/elkies-k3-2a7-2d5-4a-fixed-rank7-v1.json,
  artifacts/generated-results/elkies-k3-2a7-2d5-2c-fixed-high-mw-seed-v1.json,
  artifacts/generated-results/elkies-k3-4d6-swap-fixed-high-mw-seed-v1.json,
  artifacts/generated-results/elkies-k3-6a4-double-swap-fixed-high-mw-seed-v1.json,
  artifacts/generated-results/elkies-k3-4a5-d4-order4-fixed-high-mw-seed-v1.json,
  artifacts/generated-results/elkies-k3-4a6-4e6-fixed-coordinate-shells-v1.json,
  artifacts/generated-results/elkies-k3-8a3-fixed-coordinate-shells-v1.json,
  artifacts/generated-results/elkies-k3-6d4-fixed-coordinate-shells-v1.json,
  artifacts/generated-results/elkies-k3-3d8-fixed-coordinate-shells-v1.json,
  artifacts/generated-results/elkies-k3-2a9-d6-fixed-coordinate-shells-v1.json,
  artifacts/generated-results/elkies-k3-12a2-fixed-coordinate-shells-v1.json,
  artifacts/generated-results/elkies-k3-leech-co0-backend-v1.json,
  artifacts/generated-results/elkies-k3-leech-minimal-basis-coordinate-shell-v1.json
output: artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json
"""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import QQ, ZZ, matrix, pari


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_FOUNDRY = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
DEFAULT_GOLAY = ROOT / "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json"
DEFAULT_NIEMEIER = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
DEFAULT_24A1_PREFIX = ROOT / "artifacts/generated-results/elkies-k3-24a1-octad-prefix-orbits-v1.json"
DEFAULT_24A1_COMPLETION_MANIFEST = (
    ROOT
    / "artifacts/generated-results/elkies-k3-24a1-octad-completion-manifest-v1.json"
)
DEFAULT_MOD2_PRIORITY = (
    ROOT
    / "artifacts/generated-results/elkies-k3-cross-niemeier-mod2-priority-v1.json"
)
DEFAULT_4A_FIXED = (
    ROOT
    / "artifacts/generated-results/elkies-k3-2a7-2d5-4a-fixed-rank7-v1.json"
)
DEFAULT_2C_FIXED_SEED = (
    ROOT
    / "artifacts/generated-results/elkies-k3-2a7-2d5-2c-fixed-high-mw-seed-v1.json"
)
DEFAULT_4D6_SWAP_FIXED_SEED = (
    ROOT
    / "artifacts/generated-results/elkies-k3-4d6-swap-fixed-high-mw-seed-v1.json"
)
DEFAULT_6A4_DOUBLE_SWAP_FIXED_SEED = (
    ROOT
    / "artifacts/generated-results/elkies-k3-6a4-double-swap-fixed-high-mw-seed-v1.json"
)
DEFAULT_4A5_D4_ORDER4_FIXED_SEED = (
    ROOT
    / "artifacts/generated-results/elkies-k3-4a5-d4-order4-fixed-high-mw-seed-v1.json"
)
DEFAULT_4A6_4E6_FIXED_SHELLS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-4a6-4e6-fixed-coordinate-shells-v1.json"
)
DEFAULT_8A3_FIXED_SHELLS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-8a3-fixed-coordinate-shells-v1.json"
)
DEFAULT_6D4_FIXED_SHELLS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-6d4-fixed-coordinate-shells-v1.json"
)
DEFAULT_3D8_FIXED_SHELLS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-3d8-fixed-coordinate-shells-v1.json"
)
DEFAULT_2A9_D6_FIXED_SHELLS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-2a9-d6-fixed-coordinate-shells-v1.json"
)
DEFAULT_3A8_FIXED_SHELLS = (
    ROOT / "artifacts/generated-results/elkies-k3-3a8-fixed-coordinate-shells-v1.json"
)
DEFAULT_12A2_FIXED_SHELLS = (
    ROOT / "artifacts/generated-results/elkies-k3-12a2-fixed-coordinate-shells-v1.json"
)
DEFAULT_LEECH = ROOT / "artifacts/generated-results/elkies-k3-leech-co0-backend-v1.json"
DEFAULT_LEECH_COORDINATE_SHELL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-leech-minimal-basis-coordinate-shell-v1.json"
)
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json"

BANDS = (
    ("D0001-0500", 1, 500),
    ("D0501-1000", 501, 1000),
    ("D1001-2000", 1001, 2000),
    ("D2001-5000", 2001, 5000),
)


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def compact_digest(payload, length=16):
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()[:length]


def gram_digest(gram):
    return hashlib.sha256(
        ("\n".join(" ".join(map(str, row)) for row in gram.rows()) + "\n").encode()
    ).hexdigest()


def negate_discriminant_form_normal_key(normal):
    # The source artifacts store q_M in a canonical normal basis. Since
    # NS=U+M(-1), the same finite-module basis has Gram -q_M. We deliberately
    # retain that basis instead of recomputing an indefinite genus.
    return {
        "invariants": normal["invariants"],
        "quadratic_gram": [
            [str(-QQ(entry)) for entry in row]
            for row in normal["quadratic_gram"]
        ],
        "value_module": normal["value_module"],
        "basis": "negative_of_source_frame_discriminant_normal_basis",
    }


def exact_isometric(left, right):
    return (
        left.nrows() == right.nrows()
        and left.det() == right.det()
        and pari(left).qfisom(pari(right)) != 0
    )


def reduced_gram(gram):
    change = gram.LLL_gram()
    reduced = change.transpose() * gram * change
    assert abs(change.det()) == 1
    assert reduced.det() == gram.det()
    return reduced


def band_for(determinant):
    for band_id, lower, upper in BANDS:
        if lower <= determinant <= upper:
            return band_id
    raise AssertionError(f"determinant outside declared catalogue bands: {determinant}")


def surface_key(ns_discriminant_form, transcendental_gram):
    # For signature (1,18), rank >= length(A)+2 in every imported record, so
    # the discriminant form fixes the NS genus/isometry class.  T is retained
    # literally: different ternary lattices in the same genus are different
    # catalogue surfaces.
    return {
        "ns_discriminant_form_key": ns_discriminant_form,
        "transcendental_gram": rows(transcendental_gram),
    }


def frame_payload(frame, metadata, provenance):
    reduced = reduced_gram(frame)
    payload = {
        "gram": rows(frame),
        "gram_sha256": gram_digest(frame),
        "reduced_gram": rows(reduced),
        "determinant": int(frame.det()),
        "root_type": metadata["root_type"],
        "root_rank": int(metadata["root_rank"]),
        "mw_rank_for_rho_19": int(metadata["mw_rank_for_rho_19"]),
        "provenance": [provenance],
    }
    for optional in (
        "signed_root_count",
        "root_determinant",
        "root_lattice_primitive",
        "root_smith_invariants",
        "torsion_possibilities",
        "determinant_predicted_mw_regulator",
        "rootless_intrinsics",
    ):
        if optional in metadata:
            payload[optional] = metadata[optional]
    return payload


def add_frame(
    surface, frame, metadata, provenance, partner_index, certified_distinct_batch=None
):
    for index, existing in enumerate(surface["frames"]):
        if (
            certified_distinct_batch is not None
            and certified_distinct_batch in existing["_certified_distinct_batches"]
        ):
            # The imported foundry artifact has already performed exact PARI
            # frame-isometry deduplication inside each legacy NS class.
            continue
        if existing["root_type"] != metadata["root_type"]:
            continue
        if existing["root_rank"] != int(metadata["root_rank"]):
            continue
        if exact_isometric(frame, matrix(ZZ, existing["gram"])):
            existing["provenance"].append(provenance)
            if partner_index not in existing["partner_auxiliary_indices"]:
                existing["partner_auxiliary_indices"].append(partner_index)
            if certified_distinct_batch is not None:
                existing["_certified_distinct_batches"].append(
                    certified_distinct_batch
                )
            return index
    payload = frame_payload(frame, metadata, provenance)
    payload["partner_auxiliary_indices"] = [partner_index]
    payload["_certified_distinct_batches"] = (
        [certified_distinct_batch] if certified_distinct_batch is not None else []
    )
    surface["frames"].append(payload)
    return len(surface["frames"]) - 1


def add_partner(surface, auxiliary, provenance):
    for index, existing in enumerate(surface["partner_auxiliaries"]):
        if exact_isometric(auxiliary, matrix(ZZ, existing["gram"])):
            existing["provenance"].append(provenance)
            return index
    reduced = reduced_gram(auxiliary)
    surface["partner_auxiliaries"].append(
        {
            "gram": rows(auxiliary),
            "gram_sha256": gram_digest(auxiliary),
            "reduced_gram": rows(reduced),
            "rank": 7,
            "determinant": int(auxiliary.det()),
            "provenance": [provenance],
        }
    )
    return len(surface["partner_auxiliaries"]) - 1


def find_or_add_surface(surfaces, key):
    digest = compact_digest(key)
    for surface in surfaces:
        if surface["surface_key"] == key:
            return surface
    surface = {
        "surface_id": f"K3-{digest}",
        "surface_key": key,
        "determinant": abs(int(matrix(ZZ, key["transcendental_gram"]).det())),
        "determinant_band": None,
        "partner_auxiliaries": [],
        "frames": [],
        "legacy_ns_ids": [],
    }
    surface["determinant_band"] = band_for(surface["determinant"])
    surfaces.append(surface)
    return surface


def import_foundry(payload, surfaces):
    assert payload["schema"] == "elkies-k3.lattice-foundry-database.v1"
    assert payload["search_specification"]["determinant_bound"] == 5000
    imported_embeddings = 0
    for ns in payload["ns_classes"]:
        assert ns["k3_primitive_embedding_certificate"]["status"] == (
            "PASS_EXACT_TERNARY_DISCRIMINANT_FORM"
        )
        auxiliary = matrix(ZZ, ns["auxiliary_gram"])
        assert auxiliary.nrows() == 7
        assert 0 < auxiliary.det() <= 5000
        for transcendental in ns["ternary_realizations"]:
            ternary = matrix(ZZ, transcendental["gram"])
            ns_form = negate_discriminant_form_normal_key(
                ns["discriminant_form_normal_key"]
            )
            key = surface_key(ns_form, ternary)
            surface = find_or_add_surface(surfaces, key)
            if ns["ns_id"] not in surface["legacy_ns_ids"]:
                surface["legacy_ns_ids"].append(ns["ns_id"])
            partner_provenance = {
                "backend_id": "ROOTED-2A7_2D5",
                "source_artifact": "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json",
                "legacy_ns_id": ns["ns_id"],
                "scope": "complete_only_in_declared_one_root_mutation_shell",
            }
            partner_index = add_partner(surface, auxiliary, partner_provenance)
            for frame_row in ns["frames"]:
                frame = matrix(ZZ, frame_row["gram"])
                embeddings = frame_row["embeddings"]
                imported_embeddings += len(embeddings)
                provenance = {
                    "backend_id": "ROOTED-2A7_2D5",
                    "source_artifact": "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json",
                    "legacy_frame_id": frame_row["frame_id"],
                    "embedding_count": len(embeddings),
                    "embeddings": embeddings,
                    "scope": "complete_only_in_declared_one_root_mutation_shell",
                }
                add_frame(
                    surface,
                    frame,
                    frame_row,
                    provenance,
                    partner_index,
                    certified_distinct_batch=ns["ns_id"],
                )
    return imported_embeddings


def import_golay(payload, surfaces):
    assert payload["schema"] == "elkies-k3.golay-octad-rank17-design.v1"
    assert payload["status"] == "PASS_EXACT_GOLAY_OCTAD_RANK17_DET720_LATTICE_DESIGN"
    auxiliary = matrix(ZZ, payload["auxiliary"]["gram"])
    frame = matrix(ZZ, payload["frame"]["gram"])
    ternary = matrix(ZZ, payload["k3_realizability"]["transcendental_gram"])
    ns_form = negate_discriminant_form_normal_key(
        payload["frame"]["discriminant_form_normal_key"]
    )
    key = surface_key(ns_form, ternary)
    surface = find_or_add_surface(surfaces, key)
    provenance = {
        "backend_id": "ROOTED-24A1",
        "source_artifact": "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json",
        "scope": "one_exact_design_from_a_nonexhaustive_bounded_proposal",
        "auxiliary_basis_in_ambient": payload["auxiliary"]["ambient_basis"],
        "complement_basis_in_ambient": payload["frame"]["ambient_basis"],
    }
    partner_index = add_partner(surface, auxiliary, provenance)
    metadata = {
        "root_type": "0",
        "root_rank": 0,
        "mw_rank_for_rho_19": 17,
        "signed_root_count": 0,
        "root_determinant": 1,
        "root_lattice_primitive": True,
        "torsion_possibilities": [1],
        "determinant_predicted_mw_regulator": str(abs(frame.det())),
        "rootless_intrinsics": {
            "minimum_squared_norm": payload["frame"]["minimum_squared_norm"],
            "norm_four_vectors": payload["frame"]["norm_four_vectors"],
            "norm_four_unoriented_pairs": payload["frame"]["norm_four_unoriented_pairs"],
            "smith_invariants_greater_than_one": payload["frame"]["smith_invariants_greater_than_one"],
        },
    }
    add_frame(surface, frame, metadata, provenance, partner_index)
    return 1


def import_24a1_positive_octad_full_orbits(payload, surfaces):
    assert payload["schema"] == (
        "elkies-k3.24a1-weyl-m24-canonicalization.v2"
    )
    assert payload["status"] == (
        "PASS_EXACT_FULL_WEYL_M24_CANONICALIZATION_OF_DECLARED_INPUT_SHARDS"
    )
    assert payload["parameters"] == {
        "prefix_start_zero_based_inclusive": 0,
        "prefix_stop_zero_based_exclusive": 10547,
        "determinant_bound": 500,
    }
    assert payload["accounting"]["full_weyl_m24_embedding_orbits"] == 24
    assert payload["accounting"][
        "k3_compatible_full_embedding_orbits_by_ternary_genus_gate"
    ] == 18
    imported_records = 0
    for orbit in payload["embedding_orbits"]:
        matches = orbit["ternary_genus_representatives"]
        assert len(matches) == orbit["matching_even_ternary_genera"]
        if not matches:
            continue
        auxiliary = matrix(ZZ, orbit["representative_auxiliary_gram"])
        frame = matrix(ZZ, orbit["representative_frame_gram"])
        assert auxiliary.nrows() == 7 and frame.nrows() == 17
        assert auxiliary.det() == frame.det() == orbit["determinant"]
        root_rank = 17 - orbit["mordell_weil_rank"]
        assert orbit["frame_root_system"] == f"{root_rank}A1"
        ns_form = negate_discriminant_form_normal_key(
            orbit["frame_discriminant_form_normal_key"]
        )
        for ternary_rows in matches:
            ternary = matrix(ZZ, ternary_rows)
            assert abs(ternary.det()) == frame.det()
            key = surface_key(ns_form, ternary)
            surface = find_or_add_surface(surfaces, key)
            provenance = {
                "backend_id": "ROOTED-24A1",
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-24a1-weyl-m24-canonicalization-"
                    "00000-10547-v2.json"
                ),
                "full_weyl_m24_orbit_id": orbit["orbit_id"],
                "intrinsic_auxiliary_class_id": orbit[
                    "intrinsic_auxiliary_class_id"
                ],
                "full_group_stabilizer": orbit["full_group_stabilizer"],
                "representative_input_origin": orbit[
                    "representative_input_origin"
                ],
                "representative_auxiliary_basis_in_ambient": orbit[
                    "representative_auxiliary_basis_in_ambient"
                ],
                "representative_complement_basis_in_ambient": orbit[
                    "representative_complement_basis_in_ambient"
                ],
                "scope": (
                    "complete_positive_seven_octad_generator_subfamily_"
                    "through_determinant_500_under_full_2^24_semidirect_M24;_"
                    "one_exact_representative_per_matching_ternary_genus"
                ),
            }
            partner_index = add_partner(surface, auxiliary, provenance)
            metadata = {
                "root_type": orbit["frame_root_system"],
                "root_rank": root_rank,
                "mw_rank_for_rho_19": orbit["mordell_weil_rank"],
                "signed_root_count": 2 * root_rank,
                "root_determinant": 2**root_rank,
            }
            add_frame(surface, frame, metadata, provenance, partner_index)
            imported_records += 1
    assert imported_records == 18
    return imported_records


def import_2c_fixed_seed(payload, surfaces):
    assert payload["schema"] == (
        "elkies-k3.2a7-2d5-2c-fixed-high-mw-seed.v1"
    )
    assert payload["status"] == (
        "PASS_EXACT_DECLARED_2C_FIXED_HIGH_MW_SEED_SHELL_T_NS_FIRST"
    )
    assert payload["parameters"]["determinant_bound"] == 5000
    assert payload["parameters"]["minimum_mw_rank"] == 12
    assert payload["accounting"]["surface_classes_after_T_NS_first_dedup"] == 73
    imported_orbits = 0
    for source_surface in payload["surfaces_T_NS_first"]:
        key = source_surface["surface_key"]
        surface = find_or_add_surface(surfaces, key)
        partner_index_by_id = {}
        for partner in source_surface["partner_auxiliaries"]:
            auxiliary = matrix(ZZ, partner["gram"])
            provenance = {
                "backend_id": "ROOTED-2A7_2D5",
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-2a7-2d5-2c-fixed-high-mw-seed-v1.json"
                ),
                "source_surface_id": source_surface["surface_id"],
                "source_partner_id": partner["partner_id"],
                "section_orbit_ids": partner["section_orbit_ids"],
                "scope": (
                    "complete_only_in_declared_7_of_16_LLL_coordinate_"
                    "summand_shell_closed_under_Dih_4_section"
                ),
                "representative_auxiliary_basis_in_ambient": partner[
                    "representative_auxiliary_basis_in_ambient"
                ],
            }
            partner_index_by_id[partner["partner_id"]] = add_partner(
                surface, auxiliary, provenance
            )
        for source_frame in source_surface["frames"]:
            frame = matrix(ZZ, source_frame["gram"])
            root_data = source_frame["root_data"]
            metadata = {
                "root_type": root_data["root_type"],
                "root_rank": root_data["root_rank"],
                "mw_rank_for_rho_19": root_data["mw_rank_for_rho_19"],
                "signed_root_count": root_data["signed_root_count"],
                "root_determinant": root_data["root_determinant"],
                "determinant_predicted_mw_regulator": str(
                    QQ(frame.det()) / root_data["root_determinant"]
                ),
            }
            provenance = {
                "backend_id": "ROOTED-2A7_2D5",
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-2a7-2d5-2c-fixed-high-mw-seed-v1.json"
                ),
                "source_surface_id": source_surface["surface_id"],
                "source_frame_id": source_frame["frame_id"],
                "section_orbit_ids": source_frame["section_orbit_ids"],
                "representative_complement_basis_in_ambient": source_frame[
                    "representative_complement_basis_in_ambient"
                ],
                "literal_section_stabilizer": source_frame[
                    "representative_literal_section_stabilizer"
                ],
                "scope": (
                    "complete_only_in_declared_7_of_16_LLL_coordinate_"
                    "summand_shell_closed_under_Dih_4_section"
                ),
            }
            partner_indices = sorted(
                {partner_index_by_id[value] for value in source_frame["partner_ids"]}
            )
            assert partner_indices
            frame_index = add_frame(
                surface,
                frame,
                metadata,
                provenance,
                partner_indices[0],
                certified_distinct_batch=source_surface["surface_id"],
            )
            for partner_index in partner_indices[1:]:
                if partner_index not in surface["frames"][frame_index][
                    "partner_auxiliary_indices"
                ]:
                    surface["frames"][frame_index][
                        "partner_auxiliary_indices"
                    ].append(partner_index)
            imported_orbits += len(source_frame["section_orbit_ids"])
    assert imported_orbits == payload["accounting"]["Dih_4_section_embedding_orbits"]
    return imported_orbits


def import_6a4_double_swap_fixed_seed(payload, surfaces):
    assert payload["schema"] == (
        "elkies-k3.6a4-double-swap-fixed-high-mw-seed.v1"
    )
    assert payload["status"] == (
        "PASS_EXACT_DECLARED_6A4_DOUBLE_SWAP_FIXED_HIGH_MW_SEED_SHELL_T_NS_FIRST"
    )
    assert payload["parameters"]["determinant_bound"] == 5000
    assert payload["parameters"]["minimum_mw_rank"] == 12
    assert payload["residual_group"]["order"] == 240
    assert payload["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 42
    imported_orbits = 0
    for source_surface in payload["surfaces_T_NS_first"]:
        surface = find_or_add_surface(surfaces, source_surface["surface_key"])
        partner_index_by_id = {}
        for partner in source_surface["partner_auxiliaries"]:
            auxiliary = matrix(ZZ, partner["gram"])
            provenance = {
                "backend_id": "ROOTED-6A4",
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-6a4-double-swap-fixed-high-mw-seed-v1.json"
                ),
                "source_surface_id": source_surface["surface_id"],
                "source_partner_id": partner["partner_id"],
                "residual_group_orbit_ids": partner[
                    "residual_group_orbit_ids"
                ],
                "scope": (
                    "complete_only_in_declared_7_of_16_LLL_coordinate_"
                    "summand_shell_closed_under_exact_order_240_residual_group"
                ),
                "representative_auxiliary_basis_in_ambient": partner[
                    "representative_auxiliary_basis_in_ambient"
                ],
            }
            partner_index_by_id[partner["partner_id"]] = add_partner(
                surface, auxiliary, provenance
            )
        for source_frame in source_surface["frames"]:
            frame = matrix(ZZ, source_frame["gram"])
            root_data = source_frame["root_data"]
            metadata = {
                "root_type": root_data["root_type"],
                "root_rank": root_data["root_rank"],
                "mw_rank_for_rho_19": root_data["mw_rank_for_rho_19"],
                "signed_root_count": root_data["signed_root_count"],
                "root_determinant": root_data["root_determinant"],
                "determinant_predicted_mw_regulator": str(
                    QQ(frame.det()) / root_data["root_determinant"]
                ),
            }
            provenance = {
                "backend_id": "ROOTED-6A4",
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-6a4-double-swap-fixed-high-mw-seed-v1.json"
                ),
                "source_surface_id": source_surface["surface_id"],
                "source_frame_id": source_frame["frame_id"],
                "residual_group_orbit_ids": source_frame[
                    "residual_group_orbit_ids"
                ],
                "representative_complement_basis_in_ambient": source_frame[
                    "representative_complement_basis_in_ambient"
                ],
                "literal_residual_stabilizer": source_frame[
                    "representative_literal_residual_stabilizer"
                ],
                "scope": (
                    "complete_only_in_declared_7_of_16_LLL_coordinate_"
                    "summand_shell_closed_under_exact_order_240_residual_group"
                ),
            }
            partner_indices = sorted(
                {
                    partner_index_by_id[value]
                    for value in source_frame["partner_ids"]
                }
            )
            assert partner_indices
            frame_index = add_frame(
                surface,
                frame,
                metadata,
                provenance,
                partner_indices[0],
                certified_distinct_batch=source_surface["surface_id"],
            )
            for partner_index in partner_indices[1:]:
                if partner_index not in surface["frames"][frame_index][
                    "partner_auxiliary_indices"
                ]:
                    surface["frames"][frame_index][
                        "partner_auxiliary_indices"
                    ].append(partner_index)
            imported_orbits += len(source_frame["residual_group_orbit_ids"])
    assert imported_orbits == payload["accounting"][
        "residual_group_embedding_orbits"
    ]
    return imported_orbits


def import_4a5_d4_order4_fixed_seed(payload, surfaces):
    assert payload["schema"] == (
        "elkies-k3.4a5-d4-order4-fixed-high-mw-seed.v1"
    )
    assert payload["status"] == (
        "PASS_EXACT_DECLARED_4A5_D4_ORDER4_FIXED_HIGH_MW_SEED_SHELL_T_NS_FIRST"
    )
    assert payload["parameters"]["determinant_bound"] == 5000
    assert payload["parameters"]["minimum_mw_rank"] == 12
    assert payload["residual_group"]["order"] == 48
    assert payload["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 9
    imported_orbits = 0
    for source_surface in payload["surfaces_T_NS_first"]:
        surface = find_or_add_surface(surfaces, source_surface["surface_key"])
        partner_index_by_id = {}
        for partner in source_surface["partner_auxiliaries"]:
            auxiliary = matrix(ZZ, partner["gram"])
            provenance = {
                "backend_id": "ROOTED-4A5_D4",
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-4a5-d4-order4-fixed-high-mw-seed-v1.json"
                ),
                "source_surface_id": source_surface["surface_id"],
                "source_partner_id": partner["partner_id"],
                "residual_group_orbit_ids": partner[
                    "residual_group_orbit_ids"
                ],
                "scope": (
                    "complete_only_in_declared_7_of_10_LLL_coordinate_"
                    "summand_shell_closed_under_exact_order_48_residual_group"
                ),
                "representative_auxiliary_basis_in_ambient": partner[
                    "representative_auxiliary_basis_in_ambient"
                ],
            }
            partner_index_by_id[partner["partner_id"]] = add_partner(
                surface, auxiliary, provenance
            )
        for source_frame in source_surface["frames"]:
            frame = matrix(ZZ, source_frame["gram"])
            root_data = source_frame["root_data"]
            metadata = {
                "root_type": root_data["root_type"],
                "root_rank": root_data["root_rank"],
                "mw_rank_for_rho_19": root_data["mw_rank_for_rho_19"],
                "signed_root_count": root_data["signed_root_count"],
                "root_determinant": root_data["root_determinant"],
                "determinant_predicted_mw_regulator": str(
                    QQ(frame.det()) / root_data["root_determinant"]
                ),
            }
            if source_frame["rootless_intrinsics"] is not None:
                metadata["rootless_intrinsics"] = source_frame[
                    "rootless_intrinsics"
                ]
            provenance = {
                "backend_id": "ROOTED-4A5_D4",
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-4a5-d4-order4-fixed-high-mw-seed-v1.json"
                ),
                "source_surface_id": source_surface["surface_id"],
                "source_frame_id": source_frame["frame_id"],
                "residual_group_orbit_ids": source_frame[
                    "residual_group_orbit_ids"
                ],
                "representative_complement_basis_in_ambient": source_frame[
                    "representative_complement_basis_in_ambient"
                ],
                "literal_residual_stabilizer": source_frame[
                    "representative_literal_residual_stabilizer"
                ],
                "scope": (
                    "complete_only_in_declared_7_of_10_LLL_coordinate_"
                    "summand_shell_closed_under_exact_order_48_residual_group"
                ),
            }
            partner_indices = sorted(
                {
                    partner_index_by_id[value]
                    for value in source_frame["partner_ids"]
                }
            )
            assert partner_indices
            frame_index = add_frame(
                surface,
                frame,
                metadata,
                provenance,
                partner_indices[0],
                certified_distinct_batch=source_surface["surface_id"],
            )
            for partner_index in partner_indices[1:]:
                if partner_index not in surface["frames"][frame_index][
                    "partner_auxiliary_indices"
                ]:
                    surface["frames"][frame_index][
                        "partner_auxiliary_indices"
                    ].append(partner_index)
            imported_orbits += len(source_frame["residual_group_orbit_ids"])
    assert imported_orbits == payload["accounting"][
        "residual_group_embedding_orbits"
    ]
    return imported_orbits


def import_all_residual_fixed_shells(
    payload,
    surfaces,
    expected_schema,
    expected_status,
    expected_labels,
    source_artifact,
):
    assert payload["schema"] == expected_schema
    assert payload["status"] == expected_status
    assert payload["parameters"]["determinant_bound"] == 5000
    assert payload["parameters"]["minimum_mw_rank"] == 12
    assert {row["ambient_label"] for row in payload["backends"]} == set(
        expected_labels
    )
    imported_orbits = 0
    for backend in payload["backends"]:
        backend_id = backend["backend_id"]
        label = backend["ambient_label"]
        assert backend_id == f"ROOTED-{label}"
        residual_order = backend["residual_group"]["order"]
        scope = (
            "complete_only_in_declared_coordinate_summands_of_pinned_LLL_"
            "fixed_lattice_bases_for_every_nonidentity_residual_matrix_"
            "conjugacy_class_of_fixed_rank_at_least_7_closed_under_exact_"
            f"order_{residual_order}_residual_group"
        )
        for source_surface in backend["surfaces_T_NS_first"]:
            surface = find_or_add_surface(
                surfaces, source_surface["surface_key"]
            )
            partner_index_by_id = {}
            for partner in source_surface["partner_auxiliaries"]:
                auxiliary = matrix(ZZ, partner["gram"])
                provenance = {
                    "backend_id": backend_id,
                    "source_artifact": source_artifact,
                    "source_surface_id": source_surface["surface_id"],
                    "source_partner_id": partner["partner_id"],
                    "residual_group_orbit_ids": partner[
                        "residual_group_orbit_ids"
                    ],
                    "scope": scope,
                    "representative_auxiliary_basis_in_ambient": partner[
                        "representative_auxiliary_basis_in_ambient"
                    ],
                }
                partner_index_by_id[partner["partner_id"]] = add_partner(
                    surface, auxiliary, provenance
                )
            for source_frame in source_surface["frames"]:
                frame = matrix(ZZ, source_frame["gram"])
                root_data = source_frame["root_data"]
                metadata = {
                    "root_type": root_data["root_type"],
                    "root_rank": root_data["root_rank"],
                    "mw_rank_for_rho_19": root_data[
                        "mw_rank_for_rho_19"
                    ],
                    "signed_root_count": root_data["signed_root_count"],
                    "root_determinant": root_data["root_determinant"],
                    "determinant_predicted_mw_regulator": str(
                        QQ(frame.det()) / root_data["root_determinant"]
                    ),
                }
                if source_frame["rootless_intrinsics"] is not None:
                    metadata["rootless_intrinsics"] = source_frame[
                        "rootless_intrinsics"
                    ]
                provenance = {
                    "backend_id": backend_id,
                    "source_artifact": source_artifact,
                    "source_surface_id": source_surface["surface_id"],
                    "source_frame_id": source_frame["frame_id"],
                    "residual_group_orbit_ids": source_frame[
                        "residual_group_orbit_ids"
                    ],
                    "representative_complement_basis_in_ambient": source_frame[
                        "representative_complement_basis_in_ambient"
                    ],
                    "literal_residual_stabilizer": source_frame[
                        "representative_literal_residual_stabilizer"
                    ],
                    "scope": scope,
                }
                partner_indices = sorted(
                    {
                        partner_index_by_id[value]
                        for value in source_frame["partner_ids"]
                    }
                )
                assert partner_indices
                frame_index = add_frame(
                    surface,
                    frame,
                    metadata,
                    provenance,
                    partner_indices[0],
                    certified_distinct_batch=source_surface["surface_id"],
                )
                for partner_index in partner_indices[1:]:
                    if partner_index not in surface["frames"][frame_index][
                        "partner_auxiliary_indices"
                    ]:
                        surface["frames"][frame_index][
                            "partner_auxiliary_indices"
                        ].append(partner_index)
                imported_orbits += len(
                    source_frame["residual_group_orbit_ids"]
                )
    assert imported_orbits == payload["accounting"][
        "k3_compatible_residual_group_embedding_orbits"
    ]
    return imported_orbits


def finalize_surfaces(surfaces):
    surfaces.sort(key=lambda item: (item["determinant"], item["surface_id"]))
    seen_ids = set()
    for surface in surfaces:
        assert surface["surface_id"] not in seen_ids
        seen_ids.add(surface["surface_id"])
        surface["legacy_ns_ids"].sort()
        for partner in surface["partner_auxiliaries"]:
            partner["provenance"].sort(
                key=lambda item: (item["backend_id"], item.get("legacy_ns_id", ""))
            )
        surface["frames"].sort(
            key=lambda item: (
                -item["mw_rank_for_rho_19"],
                item["root_type"],
                item["gram_sha256"],
            )
        )
        for index, frame in enumerate(surface["frames"], start=1):
            frame["frame_id"] = f"{surface['surface_id']}-F{index:03d}"
            frame["partner_auxiliary_indices"].sort()
            del frame["_certified_distinct_batches"]
            frame["provenance"].sort(
                key=lambda item: (item["backend_id"], item.get("legacy_frame_id", ""))
            )


def backend_registry(
    niemeier,
    octad_prefix,
    octad_completions,
    octad_completion_manifest,
    octad_weyl_m24,
    mod2_priority,
    four_a_fixed,
    two_c_fixed_seed,
    four_d6_swap_fixed_seed,
    six_a4_double_swap_fixed_seed,
    four_a5_d4_order4_fixed_seed,
    four_a6_four_e6_fixed_shells,
    eight_a3_fixed_shells,
    six_d4_fixed_shells,
    three_d8_fixed_shells,
    two_a9_d6_fixed_shells,
    three_a8_fixed_shells,
    twelve_a2_fixed_shells,
    leech_foundation,
    leech_coordinate_shell,
    surfaces,
):
    imported = Counter()
    for surface in surfaces:
        for frame in surface["frames"]:
            for provenance in frame["provenance"]:
                imported[provenance["backend_id"]] += 1
    rooted = []
    assert octad_prefix["schema"] == "elkies-k3.24a1-octad-prefix-orbits.v1"
    assert octad_prefix["status"] == (
        "PASS_EXACT_M24_OCTAD_SUBSET_ORBITS_THROUGH_SIZE_5"
    )
    assert octad_completion_manifest["schema"] == (
        "elkies-k3.24a1-octad-completion-manifest.v1"
    )
    assert octad_completion_manifest["status"] == (
        "PASS_EXACT_CONTIGUOUS_24A1_OCTAD_COMPLETION_SHARD_MANIFEST"
    )
    assert len(octad_completions) == len(octad_completion_manifest["shards"])
    completion_ranges = []
    for octad_completion in octad_completions:
        assert octad_completion["schema"] == (
            "elkies-k3.24a1-octad-rank7-completion-shard.v1"
        )
        assert octad_completion["status"] == (
            "PASS_EXACT_DECLARED_24A1_OCTAD_COMPLETION_SHARD"
        )
        assert octad_completion["parameters"]["determinant_bound"] == 500
        completion_ranges.append(
            (
                octad_completion["parameters"][
                    "prefix_start_zero_based_inclusive"
                ],
                octad_completion["parameters"][
                    "prefix_stop_zero_based_exclusive"
                ],
            )
        )
    manifest_ranges = [
        (
            row["prefix_start_zero_based_inclusive"],
            row["prefix_stop_zero_based_exclusive"],
        )
        for row in octad_completion_manifest["shards"]
    ]
    assert sorted(completion_ranges) == manifest_ranges
    frontier_start = octad_completion_manifest["parameters"][
        "prefix_start_zero_based_inclusive"
    ]
    frontier_stop = octad_completion_manifest["parameters"][
        "prefix_stop_zero_based_exclusive"
    ]
    assert frontier_start == 0
    assert octad_weyl_m24["schema"] == (
        "elkies-k3.24a1-weyl-m24-canonicalization.v2"
    )
    assert octad_weyl_m24["status"] == (
        "PASS_EXACT_FULL_WEYL_M24_CANONICALIZATION_OF_DECLARED_INPUT_SHARDS"
    )
    assert octad_weyl_m24["parameters"] == {
        key: octad_completion_manifest["parameters"][key]
        for key in (
            "prefix_start_zero_based_inclusive",
            "prefix_stop_zero_based_exclusive",
            "determinant_bound",
        )
    }
    assert octad_weyl_m24["accounting"][
        "input_shard_local_residual_m24_records"
    ] == (
        sum(
            payload["statistics"]["retained_m24_orbits"]
            for payload in octad_completions
        )
    )
    assert mod2_priority["schema"] == (
        "elkies-k3.cross-niemeier-mod2-priority.v1"
    )
    assert mod2_priority["status"] == (
        "PASS_EXACT_PRIORITY_LEDGER_HEURISTIC_BACKEND_ORDER"
    )
    assert four_a_fixed["schema"] == (
        "elkies-k3.2a7-2d5-4a-fixed-rank7.v1"
    )
    assert four_a_fixed["status"] == (
        "PASS_EXACT_4A_FIXED_CORANK_ONE_RANK7_FAMILY_DET_LE_5000"
    )
    assert four_a_fixed["parameters"]["determinant_bound"] == 5000
    assert four_a_fixed["accounting"][
        "primitive_rank7_embeddings_in_declared_family"
    ] == 336
    assert four_a_fixed["accounting"]["frames_passing_ternary_genus_gate"] == 0
    assert four_a_fixed["accounting"][
        "family_rejected_by_discriminant_length"
    ]
    assert two_c_fixed_seed["schema"] == (
        "elkies-k3.2a7-2d5-2c-fixed-high-mw-seed.v1"
    )
    assert two_c_fixed_seed["status"] == (
        "PASS_EXACT_DECLARED_2C_FIXED_HIGH_MW_SEED_SHELL_T_NS_FIRST"
    )
    assert two_c_fixed_seed["accounting"][
        "Dih_4_section_embedding_orbits"
    ] == 97
    assert two_c_fixed_seed["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 73
    assert four_d6_swap_fixed_seed["schema"] == (
        "elkies-k3.4d6-swap-fixed-high-mw-seed.v1"
    )
    assert four_d6_swap_fixed_seed["status"] == (
        "PASS_EXACT_DECLARED_4D6_SWAP_FIXED_HIGH_MW_SEED_SHELL_T_NS_FIRST"
    )
    assert four_d6_swap_fixed_seed["accounting"][
        "coordinate_subsets_tested"
    ] == 11440
    assert four_d6_swap_fixed_seed["accounting"][
        "high_mw_mod2_accepted_seeds"
    ] == 0
    assert six_a4_double_swap_fixed_seed["schema"] == (
        "elkies-k3.6a4-double-swap-fixed-high-mw-seed.v1"
    )
    assert six_a4_double_swap_fixed_seed["status"] == (
        "PASS_EXACT_DECLARED_6A4_DOUBLE_SWAP_FIXED_HIGH_MW_SEED_SHELL_T_NS_FIRST"
    )
    assert six_a4_double_swap_fixed_seed["residual_group"]["order"] == 240
    assert six_a4_double_swap_fixed_seed["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 42
    assert four_a5_d4_order4_fixed_seed["schema"] == (
        "elkies-k3.4a5-d4-order4-fixed-high-mw-seed.v1"
    )
    assert four_a5_d4_order4_fixed_seed["status"] == (
        "PASS_EXACT_DECLARED_4A5_D4_ORDER4_FIXED_HIGH_MW_SEED_SHELL_T_NS_FIRST"
    )
    assert four_a5_d4_order4_fixed_seed["residual_group"]["order"] == 48
    assert four_a5_d4_order4_fixed_seed["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 9
    assert four_a6_four_e6_fixed_shells["schema"] == (
        "elkies-k3.4a6-4e6-fixed-coordinate-shells.v1"
    )
    assert four_a6_four_e6_fixed_shells["status"] == (
        "PASS_EXACT_DECLARED_FIXED_COORDINATE_SHELLS_T_NS_FIRST"
    )
    a6_e6_by_label = {
        row["ambient_label"]: row
        for row in four_a6_four_e6_fixed_shells["backends"]
    }
    assert a6_e6_by_label["4A6"]["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 9
    assert a6_e6_by_label["4E6"]["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 1
    assert eight_a3_fixed_shells["schema"] == (
        "elkies-k3.8a3-fixed-coordinate-shells.v1"
    )
    assert eight_a3_fixed_shells["status"] == (
        "PASS_EXACT_DECLARED_8A3_FIXED_COORDINATE_SHELLS_T_NS_FIRST"
    )
    eight_a3_backend = eight_a3_fixed_shells["backends"][0]
    assert eight_a3_backend["ambient_label"] == "8A3"
    assert eight_a3_backend["residual_group"]["order"] == 2688
    assert eight_a3_backend["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 435
    assert six_d4_fixed_shells["schema"] == (
        "elkies-k3.6d4-fixed-coordinate-shells.v1"
    )
    assert six_d4_fixed_shells["status"] == (
        "PASS_EXACT_DECLARED_6D4_FIXED_COORDINATE_SHELLS_T_NS_FIRST"
    )
    six_d4_backend = six_d4_fixed_shells["backends"][0]
    assert six_d4_backend["ambient_label"] == "6D4"
    assert six_d4_backend["residual_group"]["order"] == 2160
    assert six_d4_backend["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 218
    assert three_d8_fixed_shells["schema"] == (
        "elkies-k3.3d8-fixed-coordinate-shells.v1"
    )
    assert three_d8_fixed_shells["status"] == (
        "PASS_EXACT_DECLARED_3D8_FIXED_COORDINATE_SHELLS_T_NS_FIRST"
    )
    three_d8_backend = three_d8_fixed_shells["backends"][0]
    assert three_d8_backend["ambient_label"] == "3D8"
    assert three_d8_backend["residual_group"]["order"] == 6
    assert three_d8_backend["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 7
    assert two_a9_d6_fixed_shells["schema"] == (
        "elkies-k3.2a9-d6-fixed-coordinate-shells.v1"
    )
    assert two_a9_d6_fixed_shells["status"] == (
        "PASS_EXACT_DECLARED_2A9_D6_FIXED_COORDINATE_SHELLS_T_NS_FIRST"
    )
    two_a9_d6_backend = two_a9_d6_fixed_shells["backends"][0]
    assert two_a9_d6_backend["ambient_label"] == "2A9_D6"
    assert two_a9_d6_backend["residual_group"]["order"] == 4
    assert two_a9_d6_backend["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 5
    assert three_a8_fixed_shells["schema"] == (
        "elkies-k3.3a8-fixed-coordinate-shells.v1"
    )
    assert three_a8_fixed_shells["status"] == (
        "PASS_EXACT_DECLARED_3A8_FIXED_COORDINATE_SHELLS_T_NS_FIRST"
    )
    three_a8_backend = three_a8_fixed_shells["backends"][0]
    assert three_a8_backend["ambient_label"] == "3A8"
    assert three_a8_backend["residual_group"]["order"] == 12
    assert three_a8_backend["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 25
    assert twelve_a2_fixed_shells["schema"] == (
        "elkies-k3.12a2-fixed-coordinate-shells.v1"
    )
    assert twelve_a2_fixed_shells["status"] == (
        "PASS_EXACT_DECLARED_12A2_FIXED_COORDINATE_SHELLS_T_NS_FIRST"
    )
    twelve_a2_backend = twelve_a2_fixed_shells["backends"][0]
    assert twelve_a2_backend["ambient_label"] == "12A2"
    assert twelve_a2_backend["residual_group"]["order"] == 190080
    assert twelve_a2_backend["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 99
    priority_by_backend = {
        row["backend_id"]: row for row in mod2_priority["backends"]
    }
    for row in niemeier["rooted_niemeier_lattices"]:
        backend_id = f"ROOTED-{row['label']}"
        if row["label"] == "2A7_2D5":
            state = (
                "PARTIAL_IMPORTED_ONE_ROOT_MUTATION_AND_"
                "2C_FIXED_HIGH_MW_SEED_SHELLS"
            )
            group = "W(2A7+2D5) semidirect G^X; residual orbit data exact only in the control shell"
        elif row["label"] == "6A4":
            state = "PARTIAL_IMPORTED_DOUBLE_SWAP_FIXED_HIGH_MW_SEED_SHELL"
            group = (
                "W(6A4) semidirect G^X; the complete 240-element "
                "chamber-preserving residual group is exact in the imported shell"
            )
        elif row["label"] == "4D6":
            state = "PARTIAL_EXACT_SWAP_FIXED_HIGH_MW_NEGATIVE_SHELL"
            group = (
                "W(4D6) semidirect G^X; an exact lifted S4 section is "
                "used in the declared negative shell"
            )
        elif row["label"] == "4A5_D4":
            state = "PARTIAL_IMPORTED_ORDER4_FIXED_HIGH_MW_SEED_SHELL"
            group = (
                "W(4A5+D4) semidirect G^X; the complete 48-element "
                "chamber-preserving residual group is exact in the imported shell"
            )
        elif row["label"] in ("4A6", "4E6"):
            source_backend = a6_e6_by_label[row["label"]]
            state = "PARTIAL_IMPORTED_ALL_RESIDUAL_CLASS_FIXED_COORDINATE_SHELLS"
            group = (
                f"W({row['label']}) semidirect G^X; the complete "
                f"{source_backend['residual_group']['order']}-element "
                "chamber-preserving residual group is exact in the imported shells"
            )
        elif row["label"] == "8A3":
            state = "PARTIAL_IMPORTED_ALL_RESIDUAL_CLASS_FIXED_COORDINATE_SHELLS"
            group = (
                "W(8A3) semidirect G^X; the complete 2688-element "
                "glue-code residual group is exact in the imported shells"
            )
        elif row["label"] == "6D4":
            state = "PARTIAL_IMPORTED_ALL_RESIDUAL_CLASS_FIXED_COORDINATE_SHELLS"
            group = (
                "W(6D4) semidirect G^X; the complete 2160-element "
                "hexacode triality residual group is exact in the imported shells"
            )
        elif row["label"] == "3D8":
            state = "PARTIAL_IMPORTED_ALL_RESIDUAL_CLASS_FIXED_COORDINATE_SHELLS"
            group = (
                "W(3D8) semidirect G^X; the complete 6-element "
                "glue-preserving residual S3 is exact in the imported shells"
            )
        elif row["label"] == "2A9_D6":
            state = "PARTIAL_IMPORTED_ALL_RESIDUAL_CLASS_FIXED_COORDINATE_SHELLS"
            group = (
                "W(2A9+D6) semidirect G^X; the complete cyclic order-four "
                "glue-preserving residual group is exact in the imported shells"
            )
        elif row["label"] == "3A8":
            state = "PARTIAL_IMPORTED_ALL_RESIDUAL_CLASS_FIXED_COORDINATE_SHELLS"
            group = (
                "W(3A8) semidirect G^X; the complete order-twelve "
                "glue-preserving residual group {+/-1} times S3 is exact in "
                "the imported shells"
            )
        elif row["label"] == "12A2":
            state = "PARTIAL_IMPORTED_ALL_RESIDUAL_CLASS_FIXED_COORDINATE_SHELLS"
            group = (
                "W(12A2) semidirect G^X; the complete order-190080 "
                "ternary Golay monomial residual group 2.M12 is exact in "
                "the imported shells"
            )
        elif row["label"] == "3E8":
            state = "PARTIAL_EXACT_ALL_RESIDUAL_CLASS_FIXED_COORDINATE_NEGATIVE_SHELL"
            group = (
                "W(3E8) semidirect G^X; the complete 6-element residual S3 "
                "is exact in the declared negative coordinate shell"
            )
        elif row["label"] == "2D12":
            state = "PARTIAL_EXACT_ALL_RESIDUAL_CLASS_FIXED_COORDINATE_NEGATIVE_SHELL"
            group = (
                "W(2D12) semidirect G^X; the complete 2-element residual S2 "
                "is exact in the declared negative coordinate shell"
            )
        elif row["label"] == "D10_2E7":
            state = "PARTIAL_EXACT_ALL_RESIDUAL_CLASS_FIXED_COORDINATE_NEGATIVE_SHELL"
            group = (
                "W(D10+2E7) semidirect G^X; the complete order-two coupled "
                "D10-diagram/E7-swap residual group is exact in the declared "
                "negative coordinate shell"
            )
        elif row["label"] == "2A12":
            state = "PARTIAL_EXACT_ALL_ELIGIBLE_RESIDUAL_CLASS_FIXED_COORDINATE_NEGATIVE_SHELL"
            group = (
                "W(2A12) semidirect G^X; the complete cyclic order-four "
                "residual group is exact and its sole eligible fixed-coordinate "
                "shell is negative"
            )
        elif row["label"] in (
            "D24",
            "D16_E8",
            "A24",
            "A17_E7",
            "A15_D9",
            "A11_D7_E6",
        ):
            eta_control = mod2_priority["exact_eta_only_controls"][row["label"]]
            if eta_control["residual_group_order"] == 1:
                state = "EXACT_TRIVIAL_RESIDUAL_GROUP_AUXILIARY_ENUMERATION_OPEN"
                group = (
                    f"W({row['label']}); the exhaustive chamber lift test proves "
                    "that G^X is trivial"
                )
            else:
                state = "PARTIAL_EXACT_ETA_FIXED_COORDINATE_NEGATIVE_SHELL"
                group = (
                    f"W({row['label']}) semidirect C2; the complete eta residual "
                    "involution and its declared negative coordinate shell are exact"
                )
        elif row["label"] == "24A1":
            assert frontier_stop == 10547
            state = (
                "PARTIAL_COMPLETE_POSITIVE_SEVEN_OCTAD_"
                "SUBFAMILY_FULL_WEYL_M24"
            )
            group = "2^24 root sign changes semidirect M24, with Golay/M24 support canonicalization required"
        else:
            state = "OPEN_NOT_ENUMERATED"
            group = "W(X) semidirect G^X"
        backend = {
                "backend_id": backend_id,
                "kind": "rooted_niemeier",
                "ambient_label": row["label"],
                "ambient_gram_sha256": gram_digest(matrix(ZZ, row["gram"])),
                "orbit_group": group,
                "state": state,
                "imported_frame_isometry_classes": imported[backend_id],
                "completeness_through_determinant_5000": False,
                "cross_niemeier_mod2_priority": {
                    "tier": priority_by_backend[backend_id]["priority_tier"],
                    "requested_stabilizer_action_types": priority_by_backend[
                        backend_id
                    ]["requested_stabilizer_action_types"],
                    "evidence_status": priority_by_backend[backend_id][
                        "evidence_status"
                    ],
                    "required_gate": "rank_GF2(g_M-I)>0",
                },
            }
        if row["label"] == "24A1":
            backend["octad_prefix_payload_sha256"] = hashlib.sha256(
                json.dumps(octad_prefix, sort_keys=True).encode()
            ).hexdigest()
            backend["complete_octad_subset_orbits_through_size"] = 5
            backend["size_five_octad_orbits"] = octad_prefix["accounting"][-1][
                "orbit_count"
            ]
            backend["rank7_completion_payload_sha256s"] = [
                hashlib.sha256(
                    json.dumps(payload, sort_keys=True).encode()
                ).hexdigest()
                for payload in octad_completions
            ]
            backend["rank7_completion_manifest_sha256"] = hashlib.sha256(
                json.dumps(octad_completion_manifest, sort_keys=True).encode()
            ).hexdigest()
            backend["rank7_completion_prefixes_processed"] = sum(
                payload["statistics"]["prefixes_processed"]
                for payload in octad_completions
            )
            backend["rank7_completion_shard_local_residual_m24_records"] = sum(
                payload["statistics"]["retained_m24_orbits"]
                for payload in octad_completions
            )
            backend["rank7_completion_shard_local_k3_compatible_genus_records"] = sum(
                record["k3_discriminant_gate"]["matching_even_ternary_genera"] > 0
                for payload in octad_completions
                for record in payload["orbits"]
            )
            backend["rank7_completion_weyl_sign_canonicalization"] = (
                f"PASS_EXACT_ON_DECLARED_0_{frontier_stop}_INPUT"
            )
            backend["weyl_m24_payload_sha256"] = hashlib.sha256(
                json.dumps(octad_weyl_m24, sort_keys=True).encode()
            ).hexdigest()
            backend["rank7_full_weyl_m24_embedding_orbits"] = octad_weyl_m24[
                "accounting"
            ]["full_weyl_m24_embedding_orbits"]
            backend["rank7_k3_compatible_full_weyl_m24_orbits"] = (
                octad_weyl_m24["accounting"][
                    "k3_compatible_full_embedding_orbits_by_ternary_genus_gate"
                ]
            )
        if row["label"] == "2A7_2D5":
            backend["exact_4A_fixed_corank_one_family"] = {
                "status": four_a_fixed["status"],
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-2a7-2d5-4a-fixed-rank7-v1.json"
                ),
                "primitive_embeddings": four_a_fixed["accounting"][
                    "primitive_rank7_embeddings_in_declared_family"
                ],
                "auxiliary_isometry_classes": four_a_fixed["accounting"][
                    "auxiliary_isometry_classes"
                ],
                "frame_isometry_classes": four_a_fixed["accounting"][
                    "frame_isometry_classes_within_auxiliary_classes"
                ],
                "mw_rank_distribution": four_a_fixed["accounting"][
                    "mw_rank_distribution"
                ],
                "contains_literal_stabilizer_classes": ["2B", "2C", "4A"],
                "all_requested_classes_nontrivial_mod_2": True,
                "k3_compatible_frames": 0,
                "rejection": (
                    "all discriminant groups have length 7, exceeding the "
                    "maximum length 3 of a rank-three transcendental lattice"
                ),
                "completeness_scope": (
                    "all primitive corank-one sublattices of Fix(4A) through "
                    "determinant 5000; not all rank-seven sublattices of the ambient"
                ),
            }
            backend["exact_2C_fixed_high_mw_seed_shell"] = {
                "status": two_c_fixed_seed["status"],
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-2a7-2d5-2c-fixed-high-mw-seed-v1.json"
                ),
                "coordinate_subsets_tested": two_c_fixed_seed["accounting"][
                    "coordinate_subsets_tested"
                ],
                "Dih_4_section_embedding_orbits": two_c_fixed_seed[
                    "accounting"
                ]["Dih_4_section_embedding_orbits"],
                "surface_classes_after_T_NS_first_dedup": two_c_fixed_seed[
                    "accounting"
                ]["surface_classes_after_T_NS_first_dedup"],
                "partner_auxiliary_isometry_classes": two_c_fixed_seed[
                    "accounting"
                ]["partner_auxiliary_isometry_classes_after_surface_dedup"],
                "frame_isometry_classes": two_c_fixed_seed["accounting"][
                    "frame_isometry_classes_after_surface_dedup"
                ],
                "mw_rank_distribution": two_c_fixed_seed["accounting"][
                    "section_orbit_mw_rank_distribution"
                ],
                "selected_stabilizer_class": "2C",
                "all_selected_actions_nontrivial_mod_2": True,
                "canonicalization": "exact under Dih_4 section; full Weyl quotient open",
                "completeness_scope": (
                    "all 7-of-16 coordinate direct summands of one pinned LLL "
                    "basis of Fix(2C), then Dih_4-section closure; not all "
                    "rank-seven sublattices of Fix(2C)"
                ),
            }
        if row["label"] == "4D6":
            backend["exact_component_swap_fixed_high_mw_seed_shell"] = {
                "status": four_d6_swap_fixed_seed["status"],
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-4d6-swap-fixed-high-mw-seed-v1.json"
                ),
                "coordinate_subsets_tested": four_d6_swap_fixed_seed[
                    "accounting"
                ]["coordinate_subsets_tested"],
                "discriminant_length_rejected": four_d6_swap_fixed_seed[
                    "accounting"
                ]["discriminant_length_rejected"],
                "mw_rank_below_factory_floor_rejected": (
                    four_d6_swap_fixed_seed["accounting"][
                        "mw_rank_below_factory_floor_rejected"
                    ]
                ),
                "high_mw_mod2_accepted_seeds": 0,
                "selected_stabilizer_type": "single D6-component transposition",
                "canonicalization": "exact under a lifted S4 section",
                "completeness_scope": (
                    "all 7-of-16 coordinate direct summands of one pinned LLL "
                    "basis of the selected swap fixed lattice; not all rank-seven "
                    "sublattices of Fix(g) or N(4D6)"
                ),
            }
        if row["label"] == "6A4":
            backend["exact_double_swap_fixed_high_mw_seed_shell"] = {
                "status": six_a4_double_swap_fixed_seed["status"],
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-6a4-double-swap-fixed-high-mw-seed-v1.json"
                ),
                "coordinate_subsets_tested": six_a4_double_swap_fixed_seed[
                    "accounting"
                ]["coordinate_subsets_tested"],
                "residual_group_order": six_a4_double_swap_fixed_seed[
                    "residual_group"
                ]["order"],
                "residual_group_embedding_orbits": (
                    six_a4_double_swap_fixed_seed["accounting"][
                        "residual_group_embedding_orbits"
                    ]
                ),
                "surface_classes_after_T_NS_first_dedup": (
                    six_a4_double_swap_fixed_seed["accounting"][
                        "surface_classes_after_T_NS_first_dedup"
                    ]
                ),
                "partner_auxiliary_isometry_classes": (
                    six_a4_double_swap_fixed_seed["accounting"][
                        "partner_auxiliary_isometry_classes_after_surface_dedup"
                    ]
                ),
                "frame_isometry_classes": six_a4_double_swap_fixed_seed[
                    "accounting"
                ]["frame_isometry_classes_after_surface_dedup"],
                "mw_rank_distribution": six_a4_double_swap_fixed_seed[
                    "accounting"
                ]["orbit_mw_rank_distribution"],
                "selected_stabilizer_class": six_a4_double_swap_fixed_seed[
                    "parameters"
                ]["selected_class"],
                "all_selected_actions_nontrivial_mod_2": True,
                "canonicalization": (
                    "exact under all 240 chamber-preserving residual "
                    "component/diagram automorphisms; full Weyl quotient open"
                ),
                "completeness_scope": (
                    "all 7-of-16 coordinate direct summands of one pinned LLL "
                    "basis of the selected double-swap fixed lattice; not all "
                    "rank-seven sublattices of Fix(g) or N(6A4)"
                ),
            }
        if row["label"] == "4A5_D4":
            backend["exact_order4_fixed_high_mw_seed_shell"] = {
                "status": four_a5_d4_order4_fixed_seed["status"],
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-4a5-d4-order4-fixed-high-mw-seed-v1.json"
                ),
                "coordinate_subsets_tested": four_a5_d4_order4_fixed_seed[
                    "accounting"
                ]["coordinate_subsets_tested"],
                "residual_group_order": four_a5_d4_order4_fixed_seed[
                    "residual_group"
                ]["order"],
                "residual_group_embedding_orbits": (
                    four_a5_d4_order4_fixed_seed["accounting"][
                        "residual_group_embedding_orbits"
                    ]
                ),
                "surface_classes_after_T_NS_first_dedup": (
                    four_a5_d4_order4_fixed_seed["accounting"][
                        "surface_classes_after_T_NS_first_dedup"
                    ]
                ),
                "partner_auxiliary_isometry_classes": (
                    four_a5_d4_order4_fixed_seed["accounting"][
                        "partner_auxiliary_isometry_classes_after_surface_dedup"
                    ]
                ),
                "frame_isometry_classes": four_a5_d4_order4_fixed_seed[
                    "accounting"
                ]["frame_isometry_classes_after_surface_dedup"],
                "mw_rank_distribution": four_a5_d4_order4_fixed_seed[
                    "accounting"
                ]["orbit_mw_rank_distribution"],
                "selected_stabilizer_class": four_a5_d4_order4_fixed_seed[
                    "parameters"
                ]["selected_class"],
                "all_selected_actions_nontrivial_mod_2": True,
                "canonicalization": (
                    "exact under all 48 chamber-preserving residual "
                    "component/diagram automorphisms; full Weyl quotient open"
                ),
                "completeness_scope": (
                    "all 7-of-10 coordinate direct summands of one pinned LLL "
                    "basis of the selected order-four fixed lattice; not all "
                    "rank-seven sublattices of Fix(g) or N(4A5+D4)"
                ),
            }
        if row["label"] in (
            "2A9_D6",
            "3A8",
            "3D8",
            "4A6",
            "4E6",
            "6D4",
            "8A3",
            "12A2",
        ):
            source_backends = {
                "2A9_D6": two_a9_d6_backend,
                "3A8": three_a8_backend,
                "12A2": twelve_a2_backend,
                "3D8": three_d8_backend,
                "4A6": a6_e6_by_label["4A6"],
                "4E6": a6_e6_by_label["4E6"],
                "6D4": six_d4_backend,
                "8A3": eight_a3_backend,
            }
            source_payloads = {
                "2A9_D6": two_a9_d6_fixed_shells,
                "3A8": three_a8_fixed_shells,
                "12A2": twelve_a2_fixed_shells,
                "3D8": three_d8_fixed_shells,
                "4A6": four_a6_four_e6_fixed_shells,
                "4E6": four_a6_four_e6_fixed_shells,
                "6D4": six_d4_fixed_shells,
                "8A3": eight_a3_fixed_shells,
            }
            source_artifacts = {
                "2A9_D6": "elkies-k3-2a9-d6-fixed-coordinate-shells-v1.json",
                "3A8": "elkies-k3-3a8-fixed-coordinate-shells-v1.json",
                "12A2": "elkies-k3-12a2-fixed-coordinate-shells-v1.json",
                "3D8": "elkies-k3-3d8-fixed-coordinate-shells-v1.json",
                "4A6": "elkies-k3-4a6-4e6-fixed-coordinate-shells-v1.json",
                "4E6": "elkies-k3-4a6-4e6-fixed-coordinate-shells-v1.json",
                "6D4": "elkies-k3-6d4-fixed-coordinate-shells-v1.json",
                "8A3": "elkies-k3-8a3-fixed-coordinate-shells-v1.json",
            }
            source_backend = source_backends[row["label"]]
            backend["exact_all_residual_class_fixed_coordinate_shells"] = {
                "status": source_payloads[row["label"]]["status"],
                "source_artifact": (
                    "artifacts/generated-results/"
                    + source_artifacts[row["label"]]
                ),
                "residual_group_order": source_backend["residual_group"][
                    "order"
                ],
                "component_permutation_image_order": source_backend[
                    "residual_group"
                ].get(
                    "component_permutation_image_order",
                    source_backend["residual_group"].get(
                        "a9_component_permutation_image_order"
                    ),
                ),
                "coordinate_subsets_tested": source_backend[
                    "source_probe_accounting"
                ]["coordinate_subsets_tested"],
                "residual_group_embedding_orbits": source_backend[
                    "accounting"
                ]["residual_group_embedding_orbits"],
                "surface_classes_after_T_NS_first_dedup": source_backend[
                    "accounting"
                ]["surface_classes_after_T_NS_first_dedup"],
                "partner_auxiliary_isometry_classes": source_backend[
                    "accounting"
                ]["partner_auxiliary_isometry_classes_after_surface_dedup"],
                "frame_isometry_classes": source_backend["accounting"][
                    "frame_isometry_classes_after_surface_dedup"
                ],
                "mw_rank_distribution": source_backend["accounting"][
                    "orbit_mw_rank_distribution"
                ],
                "nontrivial_mod2_stabilizer_class_coverage": source_backend[
                    "accounting"
                ]["nontrivial_mod2_stabilizer_class_coverage"],
                "canonicalization": (
                    "exact under the complete chamber-preserving residual "
                    "group; full Weyl quotient open"
                ),
                "completeness_scope": (
                    "all coordinate rank-seven summands of one pinned LLL "
                    "basis for every nonidentity residual matrix conjugacy "
                    "class of fixed rank at least seven; not all primitive "
                    "fixed-lattice or ambient auxiliaries"
                ),
            }
        if row["label"] == "3E8":
            backend["exact_all_residual_class_fixed_coordinate_negative_shell"] = {
                "status": (
                    "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_3E8_COORDINATE_SHELL_SCAN"
                ),
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-3e8-fixed-coordinate-shell-probe-v1.json"
                ),
                **mod2_priority["exact_3E8_control"],
                "canonicalization": "vacuous because no seed passes the pre-quotient gates",
                "completeness_scope": (
                    "all 7-of-fixed-rank coordinate summands in both nonidentity "
                    "residual-class pinned LLL bases; not all primitive invariant "
                    "rank-seven sublattices or full Weyl embedding orbits"
                ),
            }
        if row["label"] == "2D12":
            backend["exact_all_residual_class_fixed_coordinate_negative_shell"] = {
                "status": (
                    "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_2D12_COORDINATE_SHELL_SCAN"
                ),
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-2d12-fixed-coordinate-shell-probe-v1.json"
                ),
                **mod2_priority["exact_2D12_control"],
                "canonicalization": "vacuous because no seed passes the pre-quotient gates",
                "completeness_scope": (
                    "all 7-of-12 coordinate summands in the component-swap "
                    "fixed-lattice pinned LLL basis; not all primitive invariant "
                    "rank-seven sublattices or full Weyl embedding orbits"
                ),
            }
        if row["label"] == "D10_2E7":
            backend["exact_all_residual_class_fixed_coordinate_negative_shell"] = {
                "status": (
                    "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_D10_2E7_COORDINATE_SHELL_SCAN"
                ),
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-d10-2e7-fixed-coordinate-shell-probe-v1.json"
                ),
                **mod2_priority["exact_D10_2E7_control"],
                "canonicalization": "vacuous because no seed passes the pre-quotient gates",
                "completeness_scope": (
                    "all 7-of-16 coordinate summands in the coupled-residual-"
                    "involution fixed-lattice pinned LLL basis; not all primitive "
                    "invariant rank-seven sublattices or full Weyl embedding orbits"
                ),
            }
        if row["label"] == "2A12":
            backend["exact_all_eligible_residual_class_fixed_coordinate_negative_shell"] = {
                "status": (
                    "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_2A12_COORDINATE_SHELL_SCAN"
                ),
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-2a12-fixed-coordinate-shell-probe-v1.json"
                ),
                **mod2_priority["exact_2A12_control"],
                "canonicalization": "vacuous because no seed passes the pre-quotient gates",
                "completeness_scope": (
                    "all 7-of-12 coordinate summands for the sole nonidentity "
                    "residual class of fixed rank at least seven; not all "
                    "primitive invariant rank-seven sublattices or full Weyl orbits"
                ),
            }
        if row["label"] in mod2_priority["exact_eta_only_controls"]:
            eta_control = mod2_priority["exact_eta_only_controls"][row["label"]]
            backend["exact_eta_residual_group"] = {
                "status": "PASS_EXACT_ETA_ONLY_SIX_NIEMEIER_RESIDUAL_GROUPS",
                "source_artifact": (
                    "artifacts/generated-results/"
                    "elkies-k3-eta-only-niemeier-residual-groups-v1.json"
                ),
                **eta_control,
                "completeness_scope": (
                    "every product of irreducible-component Dynkin-diagram "
                    "automorphisms is tested for an integral ambient lift"
                ),
            }
            if eta_control["residual_group_order"] > 1:
                backend["exact_eta_fixed_coordinate_negative_shell"] = {
                    "status": (
                        "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_ETA_ONLY_"
                        "COORDINATE_SHELL_SCAN"
                    ),
                    "source_artifact": (
                        "artifacts/generated-results/"
                        "elkies-k3-eta-only-niemeier-fixed-coordinate-shell-probe-v1.json"
                    ),
                    "canonicalization": (
                        "vacuous because no seed passes the pre-quotient gates"
                    ),
                    "completeness_scope": (
                        "all rank-seven coordinate summands of one pinned LLL "
                        "basis for the sole nonidentity eta class; not all "
                        "primitive invariant or ambient rank-seven sublattices"
                    ),
                }
        rooted.append(backend)
    assert leech_foundation["schema"] == "elkies-k3.leech-co0-backend.v1"
    assert leech_foundation["status"] == (
        "PASS_EXACT_LEECH_GRAM_AND_CO0_ATLAS_ACTION_BACKEND_FOUNDATION"
    )
    assert leech_coordinate_shell["schema"] == (
        "elkies-k3.leech-minimal-basis-coordinate-shell.v1"
    )
    assert leech_coordinate_shell["status"] == (
        "PASS_EXACT_DECLARED_LEECH_COORDINATE_LANGUAGE_PRE_CO1"
    )
    assert leech_coordinate_shell["accounting"]["coordinate_subsets_tested"] == 346104
    assert leech_coordinate_shell["accounting"][
        "ternary_compatible_signed_basis_types"
    ] == 194
    assert leech_coordinate_shell["accounting"][
        "preliminary_T_NS_surface_keys"
    ] == 150
    leech = {
        "backend_id": "LEECH-Co0",
        "kind": "leech",
        "ambient_label": "Leech",
        "orbit_group": "Co0",
        "ambient_gram_sha256": leech_foundation["leech_lattice"]["gram_sha256"],
        "co0_generator_sha256": leech_foundation["atlas_representation"][
            "generator_sha256"
        ],
        "state": "EXACT_DECLARED_COORDINATE_LANGUAGE_PRE_CO1_ORBITS_OPEN",
        "imported_frame_isometry_classes": imported["LEECH-Co0"],
        "complements_automatically_rootless": True,
        "completeness_through_determinant_5000": False,
        "exact_minimal_basis_coordinate_shell": {
            "status": leech_coordinate_shell["status"],
            "source_artifact": (
                "artifacts/generated-results/"
                "elkies-k3-leech-minimal-basis-coordinate-shell-v1.json"
            ),
            "coordinate_subsets_tested": leech_coordinate_shell["accounting"][
                "coordinate_subsets_tested"
            ],
            "signed_permutation_basis_types": leech_coordinate_shell[
                "accounting"
            ]["signed_permutation_basis_types"],
            "ternary_compatible_signed_basis_types": leech_coordinate_shell[
                "accounting"
            ]["ternary_compatible_signed_basis_types"],
            "preliminary_T_NS_surface_keys": leech_coordinate_shell["accounting"][
                "preliminary_T_NS_surface_keys"
            ],
            "determinant_range": [
                leech_coordinate_shell["accounting"]["determinant_minimum"],
                leech_coordinate_shell["accounting"]["determinant_maximum"],
            ],
            "catalogue_import_status": (
                "withheld pending exact Co1 embedding quotient"
            ),
            "completeness_scope": (
                "all 7-of-24 coordinate summands of one certified norm-four "
                "unimodular ambient basis; not all primitive Leech auxiliaries"
            ),
        },
        "warning": (
            "The 290 Hoehn-Mason orbits classify fixed-point sublattices, not all "
            "primitive rank-seven sublattices; they are not a completeness shortcut."
        ),
    }
    return rooted, leech


def build_payload(
    foundry,
    golay,
    niemeier,
    octad_prefix,
    octad_completions,
    octad_completion_manifest,
    octad_weyl_m24,
    mod2_priority,
    four_a_fixed,
    two_c_fixed_seed,
    four_d6_swap_fixed_seed,
    six_a4_double_swap_fixed_seed,
    four_a5_d4_order4_fixed_seed,
    four_a6_four_e6_fixed_shells,
    eight_a3_fixed_shells,
    six_d4_fixed_shells,
    three_d8_fixed_shells,
    two_a9_d6_fixed_shells,
    three_a8_fixed_shells,
    twelve_a2_fixed_shells,
    leech_foundation,
    leech_coordinate_shell,
):
    frontier_stop = octad_completion_manifest["parameters"][
        "prefix_stop_zero_based_exclusive"
    ]
    surfaces = []
    imported_foundry_embeddings = import_foundry(foundry, surfaces)
    imported_golay_embeddings = import_golay(golay, surfaces)
    imported_24a1_positive_octad_full_orbits = (
        import_24a1_positive_octad_full_orbits(octad_weyl_m24, surfaces)
    )
    imported_2c_section_orbits = import_2c_fixed_seed(
        two_c_fixed_seed, surfaces
    )
    imported_6a4_residual_orbits = import_6a4_double_swap_fixed_seed(
        six_a4_double_swap_fixed_seed, surfaces
    )
    imported_4a5_d4_residual_orbits = import_4a5_d4_order4_fixed_seed(
        four_a5_d4_order4_fixed_seed, surfaces
    )
    imported_4a6_4e6_residual_orbits = import_all_residual_fixed_shells(
        four_a6_four_e6_fixed_shells,
        surfaces,
        "elkies-k3.4a6-4e6-fixed-coordinate-shells.v1",
        "PASS_EXACT_DECLARED_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        ("4A6", "4E6"),
        "artifacts/generated-results/elkies-k3-4a6-4e6-fixed-coordinate-shells-v1.json",
    )
    imported_8a3_residual_orbits = import_all_residual_fixed_shells(
        eight_a3_fixed_shells,
        surfaces,
        "elkies-k3.8a3-fixed-coordinate-shells.v1",
        "PASS_EXACT_DECLARED_8A3_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        ("8A3",),
        "artifacts/generated-results/elkies-k3-8a3-fixed-coordinate-shells-v1.json",
    )
    imported_6d4_residual_orbits = import_all_residual_fixed_shells(
        six_d4_fixed_shells,
        surfaces,
        "elkies-k3.6d4-fixed-coordinate-shells.v1",
        "PASS_EXACT_DECLARED_6D4_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        ("6D4",),
        "artifacts/generated-results/elkies-k3-6d4-fixed-coordinate-shells-v1.json",
    )
    imported_3d8_residual_orbits = import_all_residual_fixed_shells(
        three_d8_fixed_shells,
        surfaces,
        "elkies-k3.3d8-fixed-coordinate-shells.v1",
        "PASS_EXACT_DECLARED_3D8_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        ("3D8",),
        "artifacts/generated-results/elkies-k3-3d8-fixed-coordinate-shells-v1.json",
    )
    imported_2a9_d6_residual_orbits = import_all_residual_fixed_shells(
        two_a9_d6_fixed_shells,
        surfaces,
        "elkies-k3.2a9-d6-fixed-coordinate-shells.v1",
        "PASS_EXACT_DECLARED_2A9_D6_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        ("2A9_D6",),
        "artifacts/generated-results/elkies-k3-2a9-d6-fixed-coordinate-shells-v1.json",
    )
    imported_3a8_residual_orbits = import_all_residual_fixed_shells(
        three_a8_fixed_shells,
        surfaces,
        "elkies-k3.3a8-fixed-coordinate-shells.v1",
        "PASS_EXACT_DECLARED_3A8_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        ("3A8",),
        "artifacts/generated-results/elkies-k3-3a8-fixed-coordinate-shells-v1.json",
    )
    imported_12a2_residual_orbits = import_all_residual_fixed_shells(
        twelve_a2_fixed_shells,
        surfaces,
        "elkies-k3.12a2-fixed-coordinate-shells.v1",
        "PASS_EXACT_DECLARED_12A2_FIXED_COORDINATE_SHELLS_T_NS_FIRST",
        ("12A2",),
        "artifacts/generated-results/elkies-k3-12a2-fixed-coordinate-shells-v1.json",
    )
    finalize_surfaces(surfaces)
    rooted, leech = backend_registry(
        niemeier,
        octad_prefix,
        octad_completions,
        octad_completion_manifest,
        octad_weyl_m24,
        mod2_priority,
        four_a_fixed,
        two_c_fixed_seed,
        four_d6_swap_fixed_seed,
        six_a4_double_swap_fixed_seed,
        four_a5_d4_order4_fixed_seed,
        four_a6_four_e6_fixed_shells,
        eight_a3_fixed_shells,
        six_d4_fixed_shells,
        three_d8_fixed_shells,
        two_a9_d6_fixed_shells,
        three_a8_fixed_shells,
        twelve_a2_fixed_shells,
        leech_foundation,
        leech_coordinate_shell,
        surfaces,
    )

    frame_rows = [frame for surface in surfaces for frame in surface["frames"]]
    mw_distribution = Counter(frame["mw_rank_for_rho_19"] for frame in frame_rows)
    band_distribution = Counter(surface["determinant_band"] for surface in surfaces)
    assert len(surfaces) == 827
    assert sum(len(surface["partner_auxiliaries"]) for surface in surfaces) == 1074
    assert len(frame_rows) == 1840
    assert mw_distribution == Counter(
        {12: 454, 13: 719, 14: 94, 15: 141, 16: 261, 17: 171}
    )
    assert band_distribution == Counter(
        {
            "D0001-0500": 71,
            "D0501-1000": 280,
            "D1001-2000": 300,
            "D2001-5000": 176,
        }
    )
    cross_backend_6a4_surfaces = [
        surface
        for surface in surfaces
        if "ROOTED-6A4"
        in {
            provenance["backend_id"]
            for frame in surface["frames"]
            for provenance in frame["provenance"]
        }
        and len(
            {
                provenance["backend_id"]
                for frame in surface["frames"]
                for provenance in frame["provenance"]
                if provenance["backend_id"] != "ROOTED-12A2"
            }
        )
        > 1
    ]
    assert {
        "K3-11919eb0b07e0580",
        "K3-ce4de7dfd2f10738",
        "K3-513a666a3ed34e0a",
    } <= {surface["surface_id"] for surface in cross_backend_6a4_surfaces}
    assert len(cross_backend_6a4_surfaces) == 29
    four_a5_d4_surfaces = [
        surface
        for surface in surfaces
        if "ROOTED-4A5_D4"
        in {
            provenance["backend_id"]
            for frame in surface["frames"]
            for provenance in frame["provenance"]
        }
    ]
    assert len(four_a5_d4_surfaces) == 9
    a6_e6_surfaces = [
        surface
        for surface in surfaces
        if {
            provenance["backend_id"]
            for frame in surface["frames"]
            for provenance in frame["provenance"]
        }
        & {"ROOTED-4A6", "ROOTED-4E6"}
    ]
    assert len(a6_e6_surfaces) == 10
    eight_a3_surfaces = [
        surface
        for surface in surfaces
        if "ROOTED-8A3"
        in {
            provenance["backend_id"]
            for frame in surface["frames"]
            for provenance in frame["provenance"]
        }
    ]
    assert len(eight_a3_surfaces) == 435
    eight_a3_backend_sets = Counter(
        tuple(
            sorted(
                {
                    provenance["backend_id"]
                    for frame in surface["frames"]
                    for provenance in frame["provenance"]
                    if provenance["backend_id"] != "ROOTED-12A2"
                }
            )
        )
        for surface in eight_a3_surfaces
    )
    assert eight_a3_backend_sets == Counter(
        {
            ("ROOTED-8A3",): 372,
            ("ROOTED-2A9_D6", "ROOTED-8A3"): 1,
            ("ROOTED-2A7_2D5", "ROOTED-8A3"): 5,
            ("ROOTED-3A8", "ROOTED-8A3"): 5,
            ("ROOTED-6A4", "ROOTED-8A3"): 6,
            ("ROOTED-6D4", "ROOTED-8A3"): 32,
            ("ROOTED-2A7_2D5", "ROOTED-6A4", "ROOTED-8A3"): 1,
            ("ROOTED-2A7_2D5", "ROOTED-6D4", "ROOTED-8A3"): 1,
            ("ROOTED-3A8", "ROOTED-6A4", "ROOTED-8A3"): 5,
            ("ROOTED-3A8", "ROOTED-6D4", "ROOTED-8A3"): 1,
            ("ROOTED-6A4", "ROOTED-6D4", "ROOTED-8A3"): 5,
            (
                "ROOTED-3A8",
                "ROOTED-6A4",
                "ROOTED-6D4",
                "ROOTED-8A3",
            ): 1,
        }
    )
    six_d4_surfaces = [
        surface
        for surface in surfaces
        if "ROOTED-6D4"
        in {
            provenance["backend_id"]
            for frame in surface["frames"]
            for provenance in frame["provenance"]
        }
    ]
    assert len(six_d4_surfaces) == 218
    six_d4_backend_sets = Counter(
        tuple(
            sorted(
                {
                    provenance["backend_id"]
                    for frame in surface["frames"]
                    for provenance in frame["provenance"]
                    if provenance["backend_id"] != "ROOTED-12A2"
                }
            )
        )
        for surface in six_d4_surfaces
    )
    assert six_d4_backend_sets == Counter(
        {
            ("ROOTED-6D4",): 168,
            ("ROOTED-2A7_2D5", "ROOTED-6D4"): 7,
            ("ROOTED-4A5_D4", "ROOTED-6D4"): 1,
            ("ROOTED-4E6", "ROOTED-6D4"): 1,
            ("ROOTED-6A4", "ROOTED-6D4"): 1,
            ("ROOTED-6D4", "ROOTED-8A3"): 32,
            ("ROOTED-2A7_2D5", "ROOTED-6D4", "ROOTED-8A3"): 1,
            ("ROOTED-3A8", "ROOTED-6D4", "ROOTED-8A3"): 1,
            ("ROOTED-6A4", "ROOTED-6D4", "ROOTED-8A3"): 5,
            (
                "ROOTED-3A8",
                "ROOTED-6A4",
                "ROOTED-6D4",
                "ROOTED-8A3",
            ): 1,
        }
    )
    three_d8_surfaces = [
        surface
        for surface in surfaces
        if "ROOTED-3D8"
        in {
            provenance["backend_id"]
            for frame in surface["frames"]
            for provenance in frame["provenance"]
        }
    ]
    assert len(three_d8_surfaces) == 7
    three_d8_backend_sets = Counter(
        tuple(
            sorted(
                {
                    provenance["backend_id"]
                    for frame in surface["frames"]
                    for provenance in frame["provenance"]
                    if provenance["backend_id"] != "ROOTED-12A2"
                }
            )
        )
        for surface in three_d8_surfaces
    )
    assert three_d8_backend_sets == Counter(
        {
            ("ROOTED-3D8",): 6,
            ("ROOTED-2A7_2D5", "ROOTED-3D8", "ROOTED-6A4"): 1,
        }
    )
    two_a9_d6_surfaces = [
        surface
        for surface in surfaces
        if "ROOTED-2A9_D6"
        in {
            provenance["backend_id"]
            for frame in surface["frames"]
            for provenance in frame["provenance"]
        }
    ]
    assert len(two_a9_d6_surfaces) == 5
    two_a9_d6_backend_sets = Counter(
        tuple(
            sorted(
                {
                    provenance["backend_id"]
                    for frame in surface["frames"]
                    for provenance in frame["provenance"]
                    if provenance["backend_id"] != "ROOTED-12A2"
                }
            )
        )
        for surface in two_a9_d6_surfaces
    )
    assert two_a9_d6_backend_sets == Counter(
        {
            ("ROOTED-2A9_D6",): 4,
            ("ROOTED-2A9_D6", "ROOTED-8A3"): 1,
        }
    )
    three_a8_surfaces = [
        surface
        for surface in surfaces
        if "ROOTED-3A8"
        in {
            provenance["backend_id"]
            for frame in surface["frames"]
            for provenance in frame["provenance"]
        }
    ]
    assert len(three_a8_surfaces) == 25
    three_a8_backend_sets = Counter(
        tuple(
            sorted(
                {
                    provenance["backend_id"]
                    for frame in surface["frames"]
                    for provenance in frame["provenance"]
                    if provenance["backend_id"] != "ROOTED-12A2"
                }
            )
        )
        for surface in three_a8_surfaces
    )
    assert three_a8_backend_sets == Counter(
        {
            ("ROOTED-3A8",): 5,
            ("ROOTED-3A8", "ROOTED-6A4"): 8,
            ("ROOTED-3A8", "ROOTED-8A3"): 5,
            ("ROOTED-3A8", "ROOTED-6A4", "ROOTED-8A3"): 5,
            ("ROOTED-3A8", "ROOTED-6D4", "ROOTED-8A3"): 1,
            (
                "ROOTED-3A8",
                "ROOTED-6A4",
                "ROOTED-6D4",
                "ROOTED-8A3",
            ): 1,
        }
    )
    twelve_a2_surfaces = [
        surface
        for surface in surfaces
        if "ROOTED-12A2"
        in {
            provenance["backend_id"]
            for frame in surface["frames"]
            for provenance in frame["provenance"]
        }
    ]
    assert len(twelve_a2_surfaces) == 99
    twelve_a2_backend_sets = Counter(
        tuple(
            sorted(
                {
                    provenance["backend_id"]
                    for frame in surface["frames"]
                    for provenance in frame["provenance"]
                }
            )
        )
        for surface in twelve_a2_surfaces
    )
    assert twelve_a2_backend_sets == Counter(
        {
            ("ROOTED-12A2",): 51,
            ("ROOTED-12A2", "ROOTED-24A1"): 1,
            ("ROOTED-12A2", "ROOTED-2A7_2D5"): 2,
            ("ROOTED-12A2", "ROOTED-3A8", "ROOTED-8A3"): 1,
            ("ROOTED-12A2", "ROOTED-6A4"): 5,
            ("ROOTED-12A2", "ROOTED-6D4"): 5,
            ("ROOTED-12A2", "ROOTED-8A3"): 13,
            ("ROOTED-12A2", "ROOTED-2A7_2D5", "ROOTED-6A4"): 1,
            ("ROOTED-12A2", "ROOTED-2A7_2D5", "ROOTED-6D4"): 1,
            ("ROOTED-12A2", "ROOTED-3A8", "ROOTED-6A4"): 1,
            ("ROOTED-12A2", "ROOTED-6A4", "ROOTED-6D4"): 1,
            ("ROOTED-12A2", "ROOTED-6D4", "ROOTED-8A3"): 6,
            ("ROOTED-12A2", "ROOTED-6A4", "ROOTED-8A3"): 2,
            (
                "ROOTED-12A2",
                "ROOTED-2A7_2D5",
                "ROOTED-6A4",
                "ROOTED-8A3",
            ): 1,
            (
                "ROOTED-12A2",
                "ROOTED-2A7_2D5",
                "ROOTED-6D4",
                "ROOTED-8A3",
            ): 1,
            (
                "ROOTED-12A2",
                "ROOTED-3A8",
                "ROOTED-6A4",
                "ROOTED-8A3",
            ): 3,
            (
                "ROOTED-12A2",
                "ROOTED-3A8",
                "ROOTED-6D4",
                "ROOTED-8A3",
            ): 1,
            (
                "ROOTED-12A2",
                "ROOTED-6A4",
                "ROOTED-6D4",
                "ROOTED-8A3",
            ): 2,
            (
                "ROOTED-12A2",
                "ROOTED-3A8",
                "ROOTED-6A4",
                "ROOTED-6D4",
                "ROOTED-8A3",
            ): 1,
        }
    )
    backend_shards = []
    for backend in rooted + [leech]:
        for band_id, lower, upper in BANDS:
            backend_shards.append(
                {
                    "shard_id": f"{backend['backend_id']}:{band_id}",
                    "backend_id": backend["backend_id"],
                    "determinant_band": [lower, upper],
                    "state": (
                        "PARTIAL_COMPLETE_POSITIVE_SEVEN_OCTAD_"
                        "SUBFAMILY_FULL_WEYL_M24"
                        if backend["backend_id"] == "ROOTED-24A1"
                        and band_id == "D0001-0500"
                        else "OPEN"
                    ),
                    "orbit_completeness_certificate": None,
                }
            )

    return {
        "schema": "elkies-k3.rank7-auxiliary-catalogue.v1",
        "status": "PASS_EXACT_SURFACE_FIRST_IMPORTED_CATALOGUE_FULL_ORBIT_CENSUS_OPEN",
        "proof_scope": {
            "proved": (
                "Every imported record is an exact primitive rank-seven Niemeier "
                "embedding with a saturated rank-17 complement and an exact ternary "
                "discriminant-form mate. Records are deduplicated first by (T,NS), "
                "then by positive-definite partner-auxiliary and frame isometry."
            ),
            "not_proved": (
                "No rooted backend or Leech backend is complete for all primitive "
                "rank-seven embeddings through determinant 5000. The imported "
                "2A7+2D5 mutation shell and 24A1 Golay design retain their original "
                "bounded completeness statements. The exact 4A-fixed corank-one "
                "family is closed and rejected by discriminant length, but does not "
                "cover non-pointwise-fixed auxiliaries. The imported 2C-fixed high-MW "
                "coordinate shell is exact only for its pinned 7-of-16 LLL seed language "
                "and Dih_4-section quotient; it is not a full fixed-lattice or Weyl-orbit "
                "census. The 4D6 component-swap coordinate shell is an exact bounded "
                "negative family, not a backend exclusion. The 6A4 double-swap shell "
                "uses the full 240-element chamber residual group but still covers only "
                "one pinned 7-of-16 fixed-lattice coordinate language, not the full Weyl "
                "quotient or all primitive auxiliaries. The 4A5+D4 order-four shell "
                "uses its full 48-element chamber residual group but covers only one "
                "pinned 7-of-10 fixed-lattice coordinate language. "
                "The 2A9+D6, 3A8, 3D8, 4A6, 4E6, 6D4, 8A3, and 12A2 shells exhaust the coordinate summands of "
                "pinned LLL bases for every eligible residual matrix conjugacy "
                "class, but do not enumerate all primitive fixed-lattice "
                "sublattices or apply the full Weyl quotient. The A24, A17+E7, "
                "A15+D9, and A11+D7+E6 eta-fixed coordinate languages are exact "
                "negative controls, while D24 and D16+E8 have trivial residual "
                "groups; none of these facts closes their all-primitive auxiliary "
                "enumeration. The 24A1 completion shards exhaust all 10,547 "
                "five-octad prefix orbits in the positive seven-octad generator "
                "language; their full Weyl-M24 quotient is exact, but signed/non-octad "
                "generator languages remain open. The Leech minimal-basis coordinate "
                "language exhausts 346,104 primitive summands and gives 150 exact "
                "preliminary T/NS keys, but those records remain unimported until the "
                "full Co1 embedding quotient is certified; it is not an all-primitive "
                "Leech shard. The imported ternary genus "
                "representatives are exact K3 examples, not an enumeration of every "
                "class inside each indefinite ternary genus. Surface IDs catalogue exact "
                "records, not an assertion that the determinant bands are closed."
            ),
        },
        "catalogue_policy": {
            "primary_deduplication": ["transcendental_lattice_T", "Neron_Severi_lattice_NS"],
            "secondary_deduplication": ["partner_auxiliary_K", "fibration_frame_M"],
            "target_mw_ranks": [12, 13, 14, 15, 16, 17],
            "lower_mw_frames_retained": True,
            "determinant_bands": [
                {"band_id": band_id, "lower": lower, "upper": upper}
                for band_id, lower, upper in BANDS
            ],
            "ns_uniqueness_gate": (
                "Imported discriminant length is at most 3; the indefinite even NS "
                "lattice has signature (1,18) and rank >= length(A_NS)+2, so its "
                "genus is a single isometry class."
            ),
            "cross_niemeier_umbral_mod2_priority": {
                "seed_classes_for_A7_2_D5_2": ["2B", "2C", "4A"],
                "negative_control_classes": ["1A", "2A"],
                "required_acceptance_gate": "rank_GF2(g_M-I)>0",
                "policy": (
                    "Prioritize full ambient stabilizers with non-scalar "
                    "component permutations whose induced complement action "
                    "remains nontrivial modulo two; component multiplicity "
                    "alone is only a scheduling heuristic."
                ),
            },
        },
        "accounting": {
            "surface_classes_by_T_NS": len(surfaces),
            "partner_auxiliary_isometry_classes": sum(
                len(surface["partner_auxiliaries"]) for surface in surfaces
            ),
            "frame_isometry_classes": len(frame_rows),
            "frames_with_mw_rank_12_through_17": sum(
                12 <= frame["mw_rank_for_rho_19"] <= 17 for frame in frame_rows
            ),
            "mw_rank_distribution": {
                str(rank): count for rank, count in sorted(mw_distribution.items())
            },
            "surface_determinant_band_distribution": {
                band_id: band_distribution[band_id] for band_id, unused_l, unused_u in BANDS
            },
            "imported_primitive_embedding_records": (
                imported_foundry_embeddings
                + imported_golay_embeddings
                + imported_24a1_positive_octad_full_orbits
                + imported_2c_section_orbits
                + imported_6a4_residual_orbits
                + imported_4a5_d4_residual_orbits
                + imported_4a6_4e6_residual_orbits
                + imported_8a3_residual_orbits
                + imported_6d4_residual_orbits
                + imported_3d8_residual_orbits
                + imported_2a9_d6_residual_orbits
                + imported_3a8_residual_orbits
                + imported_12a2_residual_orbits
            ),
            "6A4_source_surface_classes_before_global_dedup": 42,
            "6A4_new_surface_classes_after_global_T_NS_dedup": 39,
            "6A4_surface_classes_overlapping_prior_backends": 3,
            "6A4_surface_classes_in_any_current_cross_backend_overlap": len(
                cross_backend_6a4_surfaces
            ),
            "4A5_D4_source_surface_classes_before_global_dedup": 9,
            "4A5_D4_new_surface_classes_after_global_T_NS_dedup": 9,
            "4A5_D4_surface_classes_overlapping_prior_backends": 0,
            "4A6_source_surface_classes_before_global_dedup": 9,
            "4A6_new_surface_classes_after_global_T_NS_dedup": 9,
            "4A6_surface_classes_overlapping_prior_backends": 0,
            "4E6_source_surface_classes_before_global_dedup": 1,
            "4E6_new_surface_classes_after_global_T_NS_dedup": 1,
            "4E6_surface_classes_overlapping_prior_backends": 0,
            "8A3_source_surface_classes_before_global_dedup": 435,
            "8A3_new_surface_classes_after_global_T_NS_dedup": 411,
            "8A3_surface_classes_overlapping_prior_backends": 24,
            "8A3_surface_classes_overlapping_2A7_2D5": 7,
            "8A3_surface_classes_overlapping_6A4": 18,
            "6D4_source_surface_classes_before_global_dedup": 218,
            "6D4_new_surface_classes_after_global_T_NS_dedup": 168,
            "6D4_surface_classes_overlapping_prior_backends": 50,
            "6D4_surface_classes_overlapping_8A3": 40,
            "6D4_surface_classes_overlapping_2A7_2D5": 8,
            "6D4_surface_classes_overlapping_6A4": 7,
            "6D4_surface_classes_overlapping_4A5_D4": 1,
            "6D4_surface_classes_overlapping_4E6": 1,
            "6D4_partner_classes_added_after_global_dedup": 251,
            "6D4_frame_classes_added_after_global_dedup": 285,
            "3D8_source_surface_classes_before_global_dedup": 7,
            "3D8_new_surface_classes_after_global_T_NS_dedup": 6,
            "3D8_surface_classes_overlapping_prior_backends": 1,
            "3D8_surface_classes_overlapping_2A7_2D5": 1,
            "3D8_surface_classes_overlapping_6A4": 1,
            "3D8_partner_classes_added_after_global_dedup": 7,
            "3D8_frame_classes_added_after_global_dedup": 7,
            "2A9_D6_source_surface_classes_before_global_dedup": 5,
            "2A9_D6_new_surface_classes_after_global_T_NS_dedup": 4,
            "2A9_D6_surface_classes_overlapping_prior_backends": 1,
            "2A9_D6_surface_classes_overlapping_8A3": 1,
            "2A9_D6_partner_classes_added_after_global_dedup": 5,
            "2A9_D6_frame_classes_added_after_global_dedup": 5,
            "3A8_source_surface_classes_before_global_dedup": 25,
            "3A8_new_surface_classes_after_global_T_NS_dedup": 5,
            "3A8_surface_classes_overlapping_prior_backends": 20,
            "3A8_surface_classes_overlapping_6A4": 14,
            "3A8_surface_classes_overlapping_8A3": 12,
            "3A8_surface_classes_overlapping_6D4": 2,
            "3A8_partner_classes_added_after_global_dedup": 12,
            "3A8_frame_classes_added_after_global_dedup": 64,
            "12A2_source_surface_classes_before_global_dedup": 99,
            "12A2_new_surface_classes_after_global_T_NS_dedup": 52,
            "12A2_surface_classes_overlapping_prior_backends": 47,
            "12A2_surface_classes_overlapping_2A7_2D5": 6,
            "12A2_surface_classes_overlapping_3A8": 7,
            "12A2_surface_classes_overlapping_6A4": 17,
            "12A2_surface_classes_overlapping_6D4": 18,
            "12A2_surface_classes_overlapping_8A3": 31,
            "12A2_partner_classes_added_after_global_dedup": 86,
            "12A2_frame_classes_added_after_global_dedup": 143,
            "rooted_backends": len(rooted),
            "leech_backends": 1,
            "open_backend_band_shards": len(backend_shards),
            "complete_backend_band_shards": 0,
            "leech_pre_co1_coordinate_subsets_tested": leech_coordinate_shell[
                "accounting"
            ]["coordinate_subsets_tested"],
            "leech_pre_co1_ternary_compatible_signed_basis_types": (
                leech_coordinate_shell["accounting"][
                    "ternary_compatible_signed_basis_types"
                ]
            ),
            "leech_pre_co1_preliminary_T_NS_surface_keys": leech_coordinate_shell[
                "accounting"
            ]["preliminary_T_NS_surface_keys"],
        },
        "backends": {"rooted_niemeier": rooted, "leech": leech},
        "enumeration_shards": backend_shards,
        "surfaces": surfaces,
        "literature_boundary": [
            {
                "citation": "Nikulin, Integral symmetric bilinear forms and some of their applications",
                "doi": "10.1070/IM1980v014n01ABEH001060",
                "use": "indefinite even-lattice uniqueness from signature and discriminant form when rank >= length(A)+2",
            },
            {
                "citation": "Nishiyama, Japanese Journal of Mathematics 22 (1996), 293-347",
                "doi": "10.4099/math1924.22.293",
                "use": "primitive auxiliary embeddings into Niemeier lattices classify Jacobian fibrations",
            },
            {
                "citation": "Cheng-Duncan-Harvey, Umbral Moonshine and the Niemeier Lattices",
                "arxiv": "1307.5793",
                "use": "G^X = Aut(N^X)/W(X), including G^(24A1)=M24",
            },
            {
                "citation": "Hoehn-Mason, The 290 fixed-point sublattices of the Leech lattice",
                "doi": "10.1016/j.jalgebra.2015.08.028",
                "use": "Co0 orbit/stabilizer precedent and the fixed-point-only completeness warning",
            },
        ],
        "reproduction": {
            "command": (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elkies-k3/scripts/build_rank7_auxiliary_catalogue.sage"
            ),
            "check_command": (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elkies-k3/scripts/build_rank7_auxiliary_catalogue.sage --check"
            ),
        },
    }


parser = argparse.ArgumentParser()
parser.add_argument("--foundry", type=Path, default=DEFAULT_FOUNDRY)
parser.add_argument("--golay", type=Path, default=DEFAULT_GOLAY)
parser.add_argument("--niemeier", type=Path, default=DEFAULT_NIEMEIER)
parser.add_argument("--octad-prefix", type=Path, default=DEFAULT_24A1_PREFIX)
parser.add_argument(
    "--octad-completion-manifest",
    type=Path,
    default=DEFAULT_24A1_COMPLETION_MANIFEST,
)
parser.add_argument("--octad-completion", type=Path, action="append")
parser.add_argument("--octad-weyl-m24", type=Path)
parser.add_argument("--mod2-priority", type=Path, default=DEFAULT_MOD2_PRIORITY)
parser.add_argument("--four-a-fixed", type=Path, default=DEFAULT_4A_FIXED)
parser.add_argument("--two-c-fixed-seed", type=Path, default=DEFAULT_2C_FIXED_SEED)
parser.add_argument(
    "--four-d6-swap-fixed-seed",
    type=Path,
    default=DEFAULT_4D6_SWAP_FIXED_SEED,
)
parser.add_argument(
    "--six-a4-double-swap-fixed-seed",
    type=Path,
    default=DEFAULT_6A4_DOUBLE_SWAP_FIXED_SEED,
)
parser.add_argument(
    "--four-a5-d4-order4-fixed-seed",
    type=Path,
    default=DEFAULT_4A5_D4_ORDER4_FIXED_SEED,
)
parser.add_argument(
    "--four-a6-four-e6-fixed-shells",
    type=Path,
    default=DEFAULT_4A6_4E6_FIXED_SHELLS,
)
parser.add_argument(
    "--eight-a3-fixed-shells",
    type=Path,
    default=DEFAULT_8A3_FIXED_SHELLS,
)
parser.add_argument(
    "--six-d4-fixed-shells",
    type=Path,
    default=DEFAULT_6D4_FIXED_SHELLS,
)
parser.add_argument(
    "--three-d8-fixed-shells",
    type=Path,
    default=DEFAULT_3D8_FIXED_SHELLS,
)
parser.add_argument(
    "--two-a9-d6-fixed-shells",
    type=Path,
    default=DEFAULT_2A9_D6_FIXED_SHELLS,
)
parser.add_argument(
    "--three-a8-fixed-shells",
    type=Path,
    default=DEFAULT_3A8_FIXED_SHELLS,
)
parser.add_argument(
    "--twelve-a2-fixed-shells",
    type=Path,
    default=DEFAULT_12A2_FIXED_SHELLS,
)
parser.add_argument("--leech", type=Path, default=DEFAULT_LEECH)
parser.add_argument(
    "--leech-coordinate-shell",
    type=Path,
    default=DEFAULT_LEECH_COORDINATE_SHELL,
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

foundry_payload = json.loads(arguments.foundry.read_text())
golay_payload = json.loads(arguments.golay.read_text())
niemeier_payload = json.loads(arguments.niemeier.read_text())
octad_prefix_payload = json.loads(arguments.octad_prefix.read_text())
octad_completion_manifest_bytes = arguments.octad_completion_manifest.read_bytes()
octad_completion_manifest_payload = json.loads(octad_completion_manifest_bytes)
assert octad_completion_manifest_payload["schema"] == (
    "elkies-k3.24a1-octad-completion-manifest.v1"
)
octad_completion_paths = arguments.octad_completion or [
    ROOT / row["artifact"]
    for row in octad_completion_manifest_payload["shards"]
]
if not arguments.octad_completion:
    assert all(
        hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
        for path, row in zip(
            octad_completion_paths,
            octad_completion_manifest_payload["shards"],
        )
    )
octad_completion_payloads = [
    json.loads(path.read_text()) for path in octad_completion_paths
]
octad_weyl_m24_path = arguments.octad_weyl_m24 or (
    ROOT
    / "artifacts/generated-results"
    / (
        "elkies-k3-24a1-weyl-m24-canonicalization-"
        f"{octad_completion_manifest_payload['parameters']['prefix_start_zero_based_inclusive']:05d}-"
        f"{octad_completion_manifest_payload['parameters']['prefix_stop_zero_based_exclusive']:05d}-v2.json"
    )
)
octad_weyl_m24_payload = json.loads(octad_weyl_m24_path.read_text())
mod2_priority_payload = json.loads(arguments.mod2_priority.read_text())
four_a_fixed_payload = json.loads(arguments.four_a_fixed.read_text())
two_c_fixed_seed_payload = json.loads(arguments.two_c_fixed_seed.read_text())
four_d6_swap_fixed_seed_payload = json.loads(
    arguments.four_d6_swap_fixed_seed.read_text()
)
six_a4_double_swap_fixed_seed_payload = json.loads(
    arguments.six_a4_double_swap_fixed_seed.read_text()
)
four_a5_d4_order4_fixed_seed_payload = json.loads(
    arguments.four_a5_d4_order4_fixed_seed.read_text()
)
four_a6_four_e6_fixed_shells_payload = json.loads(
    arguments.four_a6_four_e6_fixed_shells.read_text()
)
eight_a3_fixed_shells_payload = json.loads(
    arguments.eight_a3_fixed_shells.read_text()
)
six_d4_fixed_shells_payload = json.loads(
    arguments.six_d4_fixed_shells.read_text()
)
three_d8_fixed_shells_payload = json.loads(
    arguments.three_d8_fixed_shells.read_text()
)
two_a9_d6_fixed_shells_payload = json.loads(
    arguments.two_a9_d6_fixed_shells.read_text()
)
three_a8_fixed_shells_payload = json.loads(
    arguments.three_a8_fixed_shells.read_text()
)
twelve_a2_fixed_shells_payload = json.loads(
    arguments.twelve_a2_fixed_shells.read_text()
)
leech_payload = json.loads(arguments.leech.read_text())
leech_coordinate_shell_payload = json.loads(
    arguments.leech_coordinate_shell.read_text()
)
payload = build_payload(
    foundry_payload,
    golay_payload,
    niemeier_payload,
    octad_prefix_payload,
    octad_completion_payloads,
    octad_completion_manifest_payload,
    octad_weyl_m24_payload,
    mod2_priority_payload,
    four_a_fixed_payload,
    two_c_fixed_seed_payload,
    four_d6_swap_fixed_seed_payload,
    six_a4_double_swap_fixed_seed_payload,
    four_a5_d4_order4_fixed_seed_payload,
    four_a6_four_e6_fixed_shells_payload,
    eight_a3_fixed_shells_payload,
    six_d4_fixed_shells_payload,
    three_d8_fixed_shells_payload,
    two_a9_d6_fixed_shells_payload,
    three_a8_fixed_shells_payload,
    twelve_a2_fixed_shells_payload,
    leech_payload,
    leech_coordinate_shell_payload,
)
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"

if arguments.check:
    if not arguments.output.exists() or arguments.output.read_text() != encoded:
        raise SystemExit("rank-seven auxiliary catalogue artifact is stale")
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded)

print(
    "RANK7CATALOGUE|surfaces={}|partners={}|frames={}|mw12_17={}|"
    "backends=24|open_shards=96|status=PASS".format(
        payload["accounting"]["surface_classes_by_T_NS"],
        payload["accounting"]["partner_auxiliary_isometry_classes"],
        payload["accounting"]["frame_isometry_classes"],
        payload["accounting"]["frames_with_mw_rank_12_through_17"],
    )
)
