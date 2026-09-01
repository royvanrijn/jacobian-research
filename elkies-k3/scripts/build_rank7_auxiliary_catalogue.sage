#!/usr/bin/env sage
"""Build the surface-first rank-seven auxiliary catalogue.

This is the catalogue/merge layer, not an embedding enumerator.  It imports
exact primitive embeddings from independently replayable backends, groups
first by the pair (T, NS), and only then groups partner auxiliaries and frame
lattices.  Backend and determinant-band coverage remain explicit so that a
bounded discovery artifact can never masquerade as the requested all-orbit
census.

status: ACTIVE_SEARCH_INFRASTRUCTURE
claim: exact (T,NS)-first deduplication of the imported one-root 2A7+2D5
  shell and the certified 24A1 Golay-octad design.
inputs: artifacts/generated-results/elkies-k3-lattice-foundry-v1.json,
  artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json,
  artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json,
  artifacts/generated-results/elkies-k3-24a1-octad-prefix-orbits-v1.json,
  artifacts/generated-results/elkies-k3-24a1-octad-rank7-completion-00000-00250-v1.json,
  artifacts/generated-results/elkies-k3-24a1-octad-rank7-completion-00250-00500-v1.json,
  artifacts/generated-results/elkies-k3-24a1-weyl-m24-canonicalization-00000-00500-v2.json,
  artifacts/generated-results/elkies-k3-cross-niemeier-mod2-priority-v1.json,
  artifacts/generated-results/elkies-k3-leech-co0-backend-v1.json
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
DEFAULT_24A1_COMPLETIONS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-24a1-octad-rank7-completion-00000-00250-v1.json",
    ROOT
    / "artifacts/generated-results/elkies-k3-24a1-octad-rank7-completion-00250-00500-v1.json",
)
DEFAULT_24A1_WEYL_M24 = (
    ROOT
    / "artifacts/generated-results/elkies-k3-24a1-weyl-m24-canonicalization-00000-00500-v2.json"
)
DEFAULT_MOD2_PRIORITY = (
    ROOT
    / "artifacts/generated-results/elkies-k3-cross-niemeier-mod2-priority-v1.json"
)
DEFAULT_LEECH = ROOT / "artifacts/generated-results/elkies-k3-leech-co0-backend-v1.json"
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
    octad_weyl_m24,
    mod2_priority,
    leech_foundation,
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
    assert len(octad_completions) == 2
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
    assert sorted(completion_ranges) == [(0, 250), (250, 500)]
    assert octad_weyl_m24["schema"] == (
        "elkies-k3.24a1-weyl-m24-canonicalization.v2"
    )
    assert octad_weyl_m24["status"] == (
        "PASS_EXACT_FULL_WEYL_M24_CANONICALIZATION_OF_DECLARED_INPUT_SHARDS"
    )
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
    priority_by_backend = {
        row["backend_id"]: row for row in mod2_priority["backends"]
    }
    for row in niemeier["rooted_niemeier_lattices"]:
        backend_id = f"ROOTED-{row['label']}"
        if row["label"] == "2A7_2D5":
            state = "PARTIAL_IMPORTED_ONE_ROOT_MUTATION_SHELL"
            group = "W(2A7+2D5) semidirect G^X; residual orbit data exact only in the control shell"
        elif row["label"] == "24A1":
            state = "PARTIAL_POSITIVE_OCTAD_FULL_WEYL_M24_500_OF_10547"
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
                "PASS_EXACT_ON_DECLARED_0_250_INPUT"
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
        rooted.append(backend)
    assert leech_foundation["schema"] == "elkies-k3.leech-co0-backend.v1"
    assert leech_foundation["status"] == (
        "PASS_EXACT_LEECH_GRAM_AND_CO0_ATLAS_ACTION_BACKEND_FOUNDATION"
    )
    leech = {
        "backend_id": "LEECH-Co0",
        "kind": "leech",
        "ambient_label": "Leech",
        "orbit_group": "Co0",
        "ambient_gram_sha256": leech_foundation["leech_lattice"]["gram_sha256"],
        "co0_generator_sha256": leech_foundation["atlas_representation"][
            "generator_sha256"
        ],
        "state": "AMBIENT_READY_EMBEDDING_ORBITS_OPEN",
        "imported_frame_isometry_classes": imported["LEECH-Co0"],
        "complements_automatically_rootless": True,
        "completeness_through_determinant_5000": False,
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
    octad_weyl_m24,
    mod2_priority,
    leech_foundation,
):
    surfaces = []
    imported_foundry_embeddings = import_foundry(foundry, surfaces)
    imported_golay_embeddings = import_golay(golay, surfaces)
    finalize_surfaces(surfaces)
    rooted, leech = backend_registry(
        niemeier,
        octad_prefix,
        octad_completions,
        octad_weyl_m24,
        mod2_priority,
        leech_foundation,
        surfaces,
    )

    frame_rows = [frame for surface in surfaces for frame in surface["frames"]]
    mw_distribution = Counter(frame["mw_rank_for_rho_19"] for frame in frame_rows)
    band_distribution = Counter(surface["determinant_band"] for surface in surfaces)
    backend_shards = []
    for backend in rooted + [leech]:
        for band_id, lower, upper in BANDS:
            backend_shards.append(
                {
                    "shard_id": f"{backend['backend_id']}:{band_id}",
                    "backend_id": backend["backend_id"],
                    "determinant_band": [lower, upper],
                    "state": (
                        "PARTIAL_POSITIVE_OCTAD_PREFIXES_0_500_FULL_WEYL_M24"
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
                "bounded completeness statements. The 24A1 completion shard covers "
                "only positive seven-octad generators for prefixes 0:500; their "
                "full Weyl-M24 quotient is exact, but the remaining prefixes and "
                "signed/non-octad generator languages remain open. Surface IDs catalogue exact "
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
                imported_foundry_embeddings + imported_golay_embeddings
            ),
            "rooted_backends": len(rooted),
            "leech_backends": 1,
            "open_backend_band_shards": len(backend_shards),
            "complete_backend_band_shards": 0,
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
parser.add_argument("--octad-completion", type=Path, action="append")
parser.add_argument("--octad-weyl-m24", type=Path, default=DEFAULT_24A1_WEYL_M24)
parser.add_argument("--mod2-priority", type=Path, default=DEFAULT_MOD2_PRIORITY)
parser.add_argument("--leech", type=Path, default=DEFAULT_LEECH)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

foundry_payload = json.loads(arguments.foundry.read_text())
golay_payload = json.loads(arguments.golay.read_text())
niemeier_payload = json.loads(arguments.niemeier.read_text())
octad_prefix_payload = json.loads(arguments.octad_prefix.read_text())
octad_completion_paths = arguments.octad_completion or list(DEFAULT_24A1_COMPLETIONS)
octad_completion_payloads = [
    json.loads(path.read_text()) for path in octad_completion_paths
]
octad_weyl_m24_payload = json.loads(arguments.octad_weyl_m24.read_text())
mod2_priority_payload = json.loads(arguments.mod2_priority.read_text())
leech_payload = json.loads(arguments.leech.read_text())
payload = build_payload(
    foundry_payload,
    golay_payload,
    niemeier_payload,
    octad_prefix_payload,
    octad_completion_payloads,
    octad_weyl_m24_payload,
    mod2_priority_payload,
    leech_payload,
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
