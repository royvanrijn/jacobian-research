#!/usr/bin/env python3
"""Audit the committed BCR3--BCR12 ledgers without recomputation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCES = {
    "scripts/audit_hvc38_cross_construction_frontier.py":
        "18e8744f1375d50f45cfb6a2834760044cd974f87faa24e7b5334a02b281403b",
    "scripts/audit_k12_coordinate_pair_frontier.py":
        "47ed579e5103001a56840e17fda4c53a0eeac7b9cfe07e2142a7969ccb57170d",
    "scripts/audit_k12_parameterized_completion.py":
        "f46464fec86512084d6284f96c1b89c1964b533ecd6240fdbd0a2a385909e221",
    "scripts/audit_hvc38_gap_closure.py":
        "1f06861322fc5e05f38cb2af5794aee71f151cee1f41c9f0c0ef151b08dbe0d0",
    "scripts/audit_k12_z8_cubic_completion.py":
        "743bb5e34c5c309cdc97a00c9e2f28217cd2488ee91be0ebb3661faf47a53618",
    "scripts/audit_k12_single_defect_quartic_completion.py":
        "8503d3778aee82df901a9f9cf77dd13ece8392ff9ff9f8ce488d972cb98ec79b",
    "scripts/audit_hvc38_maximal_block_closure.py":
        "bfc4f6e10b9c3860421e18c0554988d3399993cf49d0ea005e2a7b7e75f2ee35",
    "scripts/verify_k12_tensor_module_frontier.py":
        "1f576bdf1ba2204659f43f0d9e94cad6ae150a7534ad86bb99ef66ee5545b726",
    "scripts/verify_k12_cubic_graph_bilinear_obstruction.py":
        "274427fc21c12c35554d46bd7da57ce057e25454272c7af2375798b5539ae416",
    "scripts/verify_k12_graph_cubic_completion_obstruction.py":
        "ce76adfa02e87cda72c81eeaeca5a9b5a71275e6149c6282cdb0b12619a8e677",
    "scripts/audit_macfarlane_g20_dimension_reduction.py":
        "d29ad30d8f4d14890dc21f24c70ced9a4cee2735cb05d8902abe53370518a404",
    "scripts/verify_macfarlane_f12_reduction.py":
        "60e1cfdf3132763d4285a296873bf82f6deacc807316600c4a60597c12a1c6ff",
    "scripts/search_k12_cubic_graph_bilinear_completions.py":
        "f68687bacd5d6d038a7ad331c0c0e8ecbacc3f2371a2a5606c4ac6a07a34dc30",
    "scripts/verify_k12_constant_graph_cubic_obstruction.py":
        "418c4bd6546510adb33eff474e379cb5a538dfc910486ad3b46beaa17ff8ae44",
}
ARTIFACTS = {
    "hvc38_cross_construction_frontier.json":
        "e7efcfc81808c3caee11e18e18cf8e650e9dd58fc6ed66995a4e39ca1016cd2d",
    "k12_coordinate_pair_frontier.json":
        "40d4c1633a5e1b7b6f33ae1ad8d45e334e9dcdc4422feb4babfe810309aa8c3f",
    "k12_parameterized_completion_frontier.json":
        "8d57ff7d90d92941179f537fdd634dd9341a7b313005dc56e0799bb0e215199c",
    "hvc38_gap_closure.json":
        "97770dd794d921498fede1951003236742d9f05c41aec30c41c2ef0e0d8b767c",
    "k12_z8_cubic_completion_frontier.json":
        "17c4899437ffd41623836c06095e60a055cae8456f750515ffcd9f2ecca6e0bc",
    "k12_single_defect_quartic_completion_frontier.json":
        "65f32980a3055619dbfe45e58fde3bec29c06493d5c0bfcee84b58ea23643336",
    "hvc38_maximal_block_closure.json":
        "2ff1efc977ff2609ab36707ef79223a898a2c651bda639d9a8ec377b3edf1414",
    "k12_tensor_module_frontier.json":
        "7eb546812451dd5770b63e9e2e6a022841a7f7be679c3bf26edaaf6086f224cc",
    "k12_cubic_graph_bilinear_obstruction.json":
        "086a81c68ba2ffbbae5450636a8f999cb32ab5ea66fdad95412a9ae79fc5b124",
    "k12_graph_cubic_completion_obstruction.json":
        "3525645bff60df3be172d91bc5c72e273ce9742a68e8ef3d1e410d9eaa21d953",
}


def load(name: str) -> dict[str, object]:
    return json.loads(
        (ROOT / "artifacts" / "generated-results" / name).read_text(
            encoding="utf-8"
        )
    )


def audit_frontiers() -> None:
    cross = load("hvc38_cross_construction_frontier.json")
    assert cross["format"] == "hvc38-cross-construction-frontier-v1"
    assert cross["public_construction"]["degree_three_dimension"] == 11
    assert cross["public_construction"]["cubic_output_rank"] == 7
    assert cross["local_independent_construction"]["degree_three_dimension"] == 12
    assert cross["local_independent_construction"]["cubic_output_rank"] == 6
    assert len(cross["public_G11_quadratic_pivot_completion_obstructions"]) == 7
    assert cross["local_F12_coordinated_source_shear"][
        "degree_preserving_kernel_dimension"
    ] == 14
    assert cross["local_F12_coordinated_source_shear"][
        "rank_at_most_five_linear_rank"
    ] == 5
    assert cross["local_F12_coordinated_source_shear"][
        "rank_at_most_five_augmented_rank"
    ] == 6
    assert "Exact bounded obstructions only" in cross["scope"]

    coordinate = load("k12_coordinate_pair_frontier.json")
    assert coordinate["format"] == "k12-coordinate-pair-frontier-v1"
    assert coordinate["status"] == "exact bounded obstruction over Q"
    assert coordinate["literal_triangular_components"] == list(range(4, 13))
    assert [row["omitted_component"] for row in coordinate["deletions"]] == list(
        range(4, 13)
    )
    linear_rows = coordinate["linear_target_coordinate_graph_audit"]
    assert [row["source_pivot_variable"] for row in linear_rows] == list(
        range(1, 13)
    )
    assert all(not row["nonzero_pivot_possible"] for row in linear_rows[:3])
    assert all(row["raw_degree_three_ideal_is_unit"] for row in linear_rows[3:])
    assert [
        row["omitted_component"]
        for row in coordinate["quartic_target_screens_for_closest_cases"]
    ] == [11, 12]
    assert "not a lower bound" in coordinate["scope"]

    parameterized = load("k12_parameterized_completion_frontier.json")
    assert parameterized["format"] == "k12-parameterized-completion-frontier-v1"
    quadratic = parameterized[
        "quadratic_graph_families_quadratic_target_completion"
    ]
    cubic = parameterized["single_defect_families_cubic_target_completion"]
    assert [row["source_pivot"] for row in quadratic] == [7, 8, 9, 10, 11, 12]
    assert [row["source_pivot"] for row in cubic] == [7, 9, 10, 11, 12]
    assert "z8 family" in parameterized["scope"]
    assert "not a dimension-eleven lower bound" in parameterized["scope"]

    gap = load("hvc38_gap_closure.json")
    assert gap["format"] == "hvc38-gap-closure-v1"
    public = gap["nonlinear_pivots"]["public_d_pivot"]
    local = gap["nonlinear_pivots"]["local_z8_pivot"]
    assert [row["target_degree"] for row in public["degree_records"]] == list(
        range(1, 9)
    )
    assert [row["target_degree"] for row in local["degree_records"]] == list(
        range(1, 9)
    )
    assert public["degree_records"][-1]["nullity"] == 3
    assert local["degree_records"][-1]["nullity"] == 5
    combined = gap["coordinated_source_target"]
    assert combined["high_degree_kernel_dimension"] == 36
    assert combined["finite_triangular_family"]["combined_family_parameters"] == 17
    assert combined["finite_triangular_family"][
        "degree_three_plus_rank_drop_groebner_basis"
    ] == ["1"]
    assert "bounded by target degree eight" in gap["scope"]


def audit_completion_closures() -> None:
    z8 = load("k12_z8_cubic_completion_frontier.json")
    assert z8["format"] == "k12-z8-cubic-completion-frontier-v1"
    assert (z8["source_pivot"], z8["maximum_target_degree"]) == (8, 3)
    assert (z8["parameter_count"], z8["full_target_basis_count"]) == (5, 285)
    assert (z8["high_degree_row_count"], z8["nonzero_high_degree_columns"]) == (
        54977,
        277,
    )
    assert len(z8["certificates"]) == 3
    assert z8["cover_ideal_is_unit"] is True
    assert z8["cover_groebner_basis"] == ["1"]
    assert "does not cover target degree at least four" in z8["scope"]

    quartic = load("k12_single_defect_quartic_completion_frontier.json")
    assert quartic["format"] == (
        "k12-single-defect-quartic-completion-frontier-v1"
    )
    families = quartic["families"]
    assert [row["source_pivot"] for row in families] == [7, 9, 10, 11, 12]
    assert all(row["maximum_target_degree"] == 4 for row in families)
    assert all(row["full_target_basis_count"] == 1000 for row in families)
    assert all(row["nonzero_high_degree_columns"] == 990 for row in families)
    assert "multi-defect z8 family at target degree four is not covered" in (
        quartic["scope"]
    )

    maximal = load("hvc38_maximal_block_closure.json")
    assert maximal["format"] == "hvc38-maximal-block-closure-v1"
    assert maximal["status"] == "exact bounded obstruction theorem"
    records = maximal["records"]
    assert len(records) == 6
    assert len({tuple(row["source_coordinates"]) for row in records}) == 6
    assert all(row["kernel_dimension_source_degree_at_least_3"] == 0 for row in records)
    for row in records:
        family = row["combined_full_kernel_family"]
        if "degree_three_plus_rank_drop_unit_ideal" in family:
            assert family["degree_three_plus_rank_drop_unit_ideal"] is True
            assert family["degree_three_plus_rank_drop_groebner_basis"] == ["1"]
        else:
            assert family["high_degree_equations"] == 0
            assert family["selected_cubic_minor"] == "6"
    assert "not a dimension-38 minimality theorem" in maximal["scope"]


def audit_tensor_and_full_graph_rows() -> None:
    tensor = load("k12_tensor_module_frontier.json")
    assert tensor["format"] == "k12-tensor-module-frontier-v1"
    assert tensor["K12"]["cubic_tensor"]["output_flattening_rank"] == 6
    assert tensor["K12"]["cubic_tensor"]["input_directional_flattening_rank"] == 12
    assert tensor["G19"]["cubic_tensor"]["output_flattening_rank"] == 18
    assert tensor["G19"]["cubic_tensor"]["input_directional_flattening_rank"] == 19
    assert tensor["consequences"]["pure_cube_summand_lower_bounds"] == {
        "G19_cubic_tensor": 19,
        "K12_cubic_tensor": 12,
    }
    assert "exclude only constant linear source quotients" in tensor["scope"]

    expected = {
        "k12_cubic_graph_bilinear_obstruction.json": (
            "k12-cubic-graph-bilinear-obstruction-v1",
            2,
        ),
        "k12_graph_cubic_completion_obstruction.json": (
            "k12-graph-cubic-completion-obstruction-v1",
            3,
        ),
    }
    for name, (format_name, maximum_degree) in expected.items():
        artifact = load(name)
        assert artifact["format"] == format_name
        assert artifact["status"] == "exact bounded obstruction over Q"
        assert artifact["desired_restricted_degree"] == 3
        assert artifact["maximum_target_completion_degree"] == maximum_degree
        assert [row["source_pivot"] for row in artifact["constant_minor_families"]] == [
            4,
            9,
            10,
            11,
            12,
        ]
        assert [row["source_pivot"] for row in artifact["stratified_families"]] == [
            5,
            6,
            7,
            8,
        ]
        assert "linear target graph coordinates" in artifact["scope"]
        assert "not a dimension-eleven lower bound" in artifact["scope"]


def audit_fail_closed_writes() -> None:
    producer_sources = tuple(SOURCES)[:10]
    for relative_path in producer_sources:
        source = (ROOT / relative_path).read_text(encoding="utf-8")
        assert '"--write"' in source
        assert "if args.write:" in source
        assert "assert OUTPUT.read_text() == serialized" in source
        assert "is stale; regenerate with --write" in source


def main() -> None:
    expected_paths = {
        **SOURCES,
        **{
            f"artifacts/generated-results/{name}": digest
            for name, digest in ARTIFACTS.items()
        },
    }
    for relative_path, expected_hash in expected_paths.items():
        actual_hash = hashlib.sha256((ROOT / relative_path).read_bytes()).hexdigest()
        assert actual_hash == expected_hash, (relative_path, actual_hash, expected_hash)
    audit_frontiers()
    audit_completion_closures()
    audit_tensor_and_full_graph_rows()
    audit_fail_closed_writes()
    print(
        "PASS: BCR3--BCR12 sources and ten committed ledgers match pinned "
        "hashes, exact route partitions, and fail-closed write behavior"
    )
    print(
        "BOUNDARY: these are coordinate-, degree-, stage-, source-block-, or "
        "tensor-architecture obstructions; none is a dimension-11/36/38 "
        "minimality theorem"
    )
    print(
        "WRITE SAFETY: ordinary checker runs compare committed bytes; artifact "
        "replacement now requires explicit --write"
    )
    print("NO RECOMPUTATION: no SymPy, Singular, or search calculation was run")


if __name__ == "__main__":
    main()
