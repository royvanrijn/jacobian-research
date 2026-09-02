#!/usr/bin/env sage
"""Build the cross-Niemeier component-permutation/mod-2 priority ledger.

This is a deterministic experiment scheduler, not an embedding census.  It
uses repeated root-component types to identify backends whose component
permutation envelope can contain non-scalar involutions or order-four
elements.  For A7^2 D5^2 it additionally consumes the exact Dih_4 umbral
section and the existing six-frame negative control: the section contains
2B, 2C, and 4A classes, while all observed stabilizer images contain only
1A/2A and act trivially on M/2M.

Future embeddings are accepted into the high-priority experiment only after
their *full ambient* stabilizer is computed, its induced complement action is
reduced modulo two, and rank(g-I)>0 is certified.  Repeated components alone
are therefore a scheduling signal, never a mathematical claim about the
actual stabilizer.

status: EXACT_PRIORITY_LEDGER_HEURISTIC_BACKEND_ORDER
"""

import argparse
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_NIEMEIER = (
    ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
)
DEFAULT_UMBRAL = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-umbral-orbits-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-cross-niemeier-mod2-priority-v1.json"
)
DEFAULT_4D6_SWAP = (
    ROOT
    / "artifacts/generated-results/elkies-k3-4d6-swap-fixed-high-mw-seed-v1.json"
)
DEFAULT_6A4_DOUBLE_SWAP = (
    ROOT
    / "artifacts/generated-results/elkies-k3-6a4-double-swap-fixed-high-mw-seed-v1.json"
)
DEFAULT_4A5_D4_ORDER4 = (
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
DEFAULT_3E8_FIXED_SHELL_PROBE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-3e8-fixed-coordinate-shell-probe-v1.json"
)
DEFAULT_3D8_FIXED_SHELLS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-3d8-fixed-coordinate-shells-v1.json"
)
DEFAULT_2D12_FIXED_SHELL_PROBE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-2d12-fixed-coordinate-shell-probe-v1.json"
)
DEFAULT_D10_2E7_FIXED_SHELL_PROBE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-d10-2e7-fixed-coordinate-shell-probe-v1.json"
)
DEFAULT_2A12_FIXED_SHELL_PROBE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-2a12-fixed-coordinate-shell-probe-v1.json"
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
DEFAULT_ETA_ONLY_RESIDUAL_GROUPS = (
    ROOT
    / "artifacts/generated-results/elkies-k3-eta-only-niemeier-residual-groups-v1.json"
)
DEFAULT_ETA_ONLY_FIXED_SHELL_PROBE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-eta-only-niemeier-fixed-coordinate-shell-probe-v1.json"
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def component_type(component):
    return f"{component['family']}{component['rank']}"


def build(
    niemeier,
    umbral,
    four_d6_swap,
    six_a4_double_swap,
    four_a5_d4_order4,
    four_a6_four_e6_fixed_shells,
    eight_a3_fixed_shells,
    six_d4_fixed_shells,
    three_e8_fixed_shell_probe,
    three_d8_fixed_shells,
    two_d12_fixed_shell_probe,
    d10_two_e7_fixed_shell_probe,
    two_a12_fixed_shell_probe,
    two_a9_d6_fixed_shells,
    three_a8_fixed_shells,
    twelve_a2_fixed_shells,
    eta_only_residual_groups,
    eta_only_fixed_shell_probe,
):
    assert niemeier["schema"] == "elkies-k3.rooted-niemeier-catalog.v1"
    assert umbral["schema"] == "elkies-k3.lattice-foundry-umbral-orbits.v1"
    assert umbral["status"] == (
        "PASS_EXACT_AMBIENT_STABILIZERS_D2_ORBITS_AND_SAMPLED_D3_ORBITS"
    )
    assert four_d6_swap["schema"] == (
        "elkies-k3.4d6-swap-fixed-high-mw-seed.v1"
    )
    assert four_d6_swap["accounting"]["coordinate_subsets_tested"] == 11440
    assert four_d6_swap["accounting"]["high_mw_mod2_accepted_seeds"] == 0
    assert six_a4_double_swap["schema"] == (
        "elkies-k3.6a4-double-swap-fixed-high-mw-seed.v1"
    )
    assert six_a4_double_swap["residual_group"]["order"] == 240
    assert six_a4_double_swap["accounting"][
        "residual_group_embedding_orbits"
    ] == 161
    assert four_a5_d4_order4["schema"] == (
        "elkies-k3.4a5-d4-order4-fixed-high-mw-seed.v1"
    )
    assert four_a5_d4_order4["residual_group"]["order"] == 48
    assert four_a5_d4_order4["accounting"][
        "residual_group_embedding_orbits"
    ] == 39
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
    assert a6_e6_by_label["4A6"]["residual_group"]["order"] == 24
    assert a6_e6_by_label["4A6"]["accounting"][
        "residual_group_embedding_orbits"
    ] == 86
    assert a6_e6_by_label["4E6"]["residual_group"]["order"] == 48
    assert a6_e6_by_label["4E6"]["accounting"][
        "residual_group_embedding_orbits"
    ] == 45
    assert eight_a3_fixed_shells["schema"] == (
        "elkies-k3.8a3-fixed-coordinate-shells.v1"
    )
    assert eight_a3_fixed_shells["status"] == (
        "PASS_EXACT_DECLARED_8A3_FIXED_COORDINATE_SHELLS_T_NS_FIRST"
    )
    eight_a3_backend = eight_a3_fixed_shells["backends"][0]
    assert eight_a3_backend["residual_group"]["order"] == 2688
    assert eight_a3_backend["accounting"][
        "residual_group_embedding_orbits"
    ] == 1162
    assert six_d4_fixed_shells["schema"] == (
        "elkies-k3.6d4-fixed-coordinate-shells.v1"
    )
    assert six_d4_fixed_shells["status"] == (
        "PASS_EXACT_DECLARED_6D4_FIXED_COORDINATE_SHELLS_T_NS_FIRST"
    )
    six_d4_backend = six_d4_fixed_shells["backends"][0]
    assert six_d4_backend["residual_group"]["order"] == 2160
    assert six_d4_backend["accounting"][
        "residual_group_embedding_orbits"
    ] == 466
    assert three_e8_fixed_shell_probe["schema"] == (
        "elkies-k3.3e8-fixed-coordinate-shell-probe.v1"
    )
    assert three_e8_fixed_shell_probe["status"] == (
        "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_3E8_COORDINATE_SHELL_SCAN"
    )
    assert three_e8_fixed_shell_probe["residual_group_order"] == 6
    assert three_e8_fixed_shell_probe["accounting"][
        "coordinate_subsets_tested"
    ] == 11448
    assert three_e8_fixed_shell_probe["accounting"][
        "high_mw_mod2_accepted_seeds_before_residual_dedup"
    ] == 0
    assert three_d8_fixed_shells["schema"] == (
        "elkies-k3.3d8-fixed-coordinate-shells.v1"
    )
    assert three_d8_fixed_shells["status"] == (
        "PASS_EXACT_DECLARED_3D8_FIXED_COORDINATE_SHELLS_T_NS_FIRST"
    )
    three_d8_backend = three_d8_fixed_shells["backends"][0]
    assert three_d8_backend["residual_group"]["order"] == 6
    assert three_d8_backend["accounting"][
        "residual_group_embedding_orbits"
    ] == 40
    assert three_d8_backend["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 7
    assert two_d12_fixed_shell_probe["schema"] == (
        "elkies-k3.2d12-fixed-coordinate-shell-probe.v1"
    )
    assert two_d12_fixed_shell_probe["status"] == (
        "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_2D12_COORDINATE_SHELL_SCAN"
    )
    assert two_d12_fixed_shell_probe["residual_group_order"] == 2
    assert two_d12_fixed_shell_probe["accounting"][
        "coordinate_subsets_tested"
    ] == 792
    assert two_d12_fixed_shell_probe["accounting"][
        "high_mw_mod2_accepted_seeds_before_residual_dedup"
    ] == 0
    assert d10_two_e7_fixed_shell_probe["schema"] == (
        "elkies-k3.d10-2e7-fixed-coordinate-shell-probe.v1"
    )
    assert d10_two_e7_fixed_shell_probe["status"] == (
        "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_D10_2E7_COORDINATE_SHELL_SCAN"
    )
    assert d10_two_e7_fixed_shell_probe["residual_group_order"] == 2
    assert d10_two_e7_fixed_shell_probe["accounting"][
        "coordinate_subsets_tested"
    ] == 11440
    assert d10_two_e7_fixed_shell_probe["accounting"][
        "high_mw_mod2_accepted_seeds_before_residual_dedup"
    ] == 0
    assert two_a12_fixed_shell_probe["schema"] == (
        "elkies-k3.2a12-fixed-coordinate-shell-probe.v1"
    )
    assert two_a12_fixed_shell_probe["status"] == (
        "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_2A12_COORDINATE_SHELL_SCAN"
    )
    assert two_a12_fixed_shell_probe["residual_group_order"] == 4
    assert two_a12_fixed_shell_probe["accounting"][
        "coordinate_subsets_tested"
    ] == 792
    assert two_a12_fixed_shell_probe["accounting"][
        "high_mw_mod2_accepted_seeds_before_residual_dedup"
    ] == 0
    assert two_a9_d6_fixed_shells["schema"] == (
        "elkies-k3.2a9-d6-fixed-coordinate-shells.v1"
    )
    assert two_a9_d6_fixed_shells["status"] == (
        "PASS_EXACT_DECLARED_2A9_D6_FIXED_COORDINATE_SHELLS_T_NS_FIRST"
    )
    two_a9_d6_backend = two_a9_d6_fixed_shells["backends"][0]
    assert two_a9_d6_backend["residual_group"]["order"] == 4
    assert two_a9_d6_backend["accounting"][
        "residual_group_embedding_orbits"
    ] == 32
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
    assert three_a8_backend["residual_group"]["order"] == 12
    assert three_a8_backend["accounting"][
        "residual_group_embedding_orbits"
    ] == 189
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
    assert twelve_a2_backend["residual_group"]["order"] == 190080
    assert twelve_a2_backend["accounting"][
        "residual_group_embedding_orbits"
    ] == 214
    assert twelve_a2_backend["accounting"][
        "surface_classes_after_T_NS_first_dedup"
    ] == 99
    assert eta_only_residual_groups["schema"] == (
        "elkies-k3.eta-only-niemeier-residual-groups.v1"
    )
    assert eta_only_residual_groups["status"] == (
        "PASS_EXACT_ETA_ONLY_SIX_NIEMEIER_RESIDUAL_GROUPS"
    )
    assert eta_only_fixed_shell_probe["schema"] == (
        "elkies-k3.eta-only-niemeier-fixed-coordinate-shell-probe.v1"
    )
    assert eta_only_fixed_shell_probe["status"] == (
        "PASS_EXACT_PRE_RESIDUAL_QUOTIENT_ETA_ONLY_COORDINATE_SHELL_SCAN"
    )
    eta_group_by_label = {
        row["ambient_label"]: row for row in eta_only_residual_groups["backends"]
    }
    eta_probe_by_label = {
        row["ambient_label"]: row for row in eta_only_fixed_shell_probe["backends"]
    }
    assert set(eta_group_by_label) == set(eta_probe_by_label) == {
        "D24",
        "D16_E8",
        "A24",
        "A17_E7",
        "A15_D9",
        "A11_D7_E6",
    }
    assert eta_only_fixed_shell_probe["accounting"][
        "coordinate_subsets_tested"
    ] == 35112
    assert eta_only_fixed_shell_probe["accounting"][
        "high_mw_mod2_accepted_seeds_before_residual_dedup"
    ] == 0
    exact_section_classes = Counter(row["class"] for row in umbral["group_section"])
    assert exact_section_classes == Counter(
        {"1A": 1, "2A": 1, "2B": 2, "2C": 2, "4A": 2}
    )
    observed_classes = sorted(
        {
            label
            for target in umbral["targets"]
            for label in target["full_ambient_stabilizer"][
                "umbral_image_classes"
            ]
        }
    )
    assert observed_classes == ["1A", "2A"]
    for target in umbral["targets"]:
        for action in target["full_ambient_stabilizer"]["induced_actions"]:
            if "2A" in action["compatible_umbral_classes"]:
                assert action["fixed_degree_two_rational_cosets"] == target[
                    "degree_two_rational_cosets"
                ]["count"]

    rows = []
    for ambient in niemeier["rooted_niemeier_lattices"]:
        multiplicities = Counter(
            component_type(component) for component in ambient["root_components"]
        )
        repeated = {
            key: value for key, value in sorted(multiplicities.items()) if value >= 2
        }
        swap_types = [key for key, value in repeated.items() if value >= 2]
        four_cycle_types = [key for key, value in repeated.items() if value >= 4]
        envelope_order = math.prod(math.factorial(value) for value in multiplicities.values())
        if ambient["label"] == "2A7_2D5":
            tier = 0
            reason = (
                "Exact Dih_4 section contains 2B, 2C, and 4A, but the six "
                "existing full stabilizers realize only the mod-2-trivial 1A/2A control."
            )
            requested = ["2B", "2C", "4A"]
            evidence = "EXACT_GROUP_SECTION_AND_NEGATIVE_CONTROL"
        elif ambient["label"] == "6A4":
            tier = 1
            reason = (
                "The exact order-240 chamber residual group contains a literal "
                "fixed-rank-16 double component swap; its declared coordinate "
                "shell gives 161 mod-2-nontrivial high-MW orbits."
            )
            requested = [
                six_a4_double_swap["parameters"]["selected_class"],
                "other non-scalar order-2/order-4 residual classes",
            ]
            evidence = "EXACT_RESIDUAL_GROUP_AND_POSITIVE_FIXED_SHELL"
        elif ambient["label"] == "4D6":
            tier = 1
            reason = (
                "An exact lifted S4 section contains component transpositions. "
                "The first 11,440-coordinate fixed shell has no MW12--17 seed, "
                "so other invariant languages remain the target."
            )
            requested = [
                "non-pointwise or non-coordinate component-swap invariant auxiliary",
                "order-4 component permutation",
            ]
            evidence = "EXACT_S4_SECTION_AND_NEGATIVE_FIXED_COORDINATE_CONTROL"
        elif ambient["label"] == "4A5_D4":
            tier = 1
            reason = (
                "The exact order-48 chamber residual group contains an order-four "
                "double A5-component swap. Its declared fixed coordinate shell "
                "gives 39 mod-2-nontrivial MW13/MW17 orbits."
            )
            requested = [
                four_a5_d4_order4["parameters"]["selected_class"],
                "other non-pointwise order-four invariant auxiliaries",
            ]
            evidence = "EXACT_RESIDUAL_GROUP_AND_POSITIVE_ORDER4_FIXED_SHELL"
        elif ambient["label"] == "4A6":
            tier = 1
            reason = (
                "The exact order-24 chamber residual group has A4 component "
                "image. Exhausting every eligible matrix-conjugacy-class "
                "coordinate shell gives 86 mod-2-nontrivial order-three "
                "orbits and nine local K3 surface classes."
            )
            requested = [
                "4A6-C03 and 4A6-C04 order-three stabilizers",
                "non-coordinate invariant auxiliaries for the negative order-two/order-four classes",
            ]
            evidence = "EXACT_RESIDUAL_GROUP_AND_POSITIVE_ORDER3_FIXED_SHELL"
        elif ambient["label"] == "4E6":
            tier = 1
            reason = (
                "The exact order-48 chamber residual group has S4 component "
                "image. Exhausting every eligible matrix-conjugacy-class "
                "coordinate shell gives 45 mod-2-nontrivial residual orbits "
                "and one local K3 surface class; the order-four shell is negative."
            )
            requested = [
                "4E6-C04 order-three stabilizer",
                "non-coordinate invariant auxiliaries for the negative order-four class",
            ]
            evidence = "EXACT_RESIDUAL_GROUP_AND_POSITIVE_ORDER3_FIXED_SHELL"
        elif ambient["label"] == "8A3":
            tier = 1
            reason = (
                "The exact order-2688 glue-code residual group has component "
                "image of order 1344. Exhausting every eligible fixed-coordinate "
                "shell gives 1162 residual orbits and 435 local K3 surface classes."
            )
            requested = [
                "8A3-C02/C03 involution stabilizers",
                "8A3-C05 order-three and C06 order-four stabilizers",
                "8A3-C09 order-six stabilizers",
                "non-coordinate invariant auxiliaries",
            ]
            evidence = "EXACT_GLUE_RESIDUAL_GROUP_AND_POSITIVE_ALL_CLASS_FIXED_SHELLS"
        elif ambient["label"] == "6D4":
            tier = 1
            reason = (
                "The exact order-2160 hexacode residual group has triality "
                "kernel 3 and component image S6. Every eligible fixed-coordinate "
                "shell gives 466 residual orbits and 218 local K3 surface classes."
            )
            requested = [
                "6D4-C02/C03 involution stabilizers",
                "6D4-C06 order-three and C08 order-four stabilizers",
                "6D4-C10 order-five and C12 order-six stabilizers",
                "non-coordinate invariant auxiliaries",
            ]
            evidence = "EXACT_HEXACODE_RESIDUAL_GROUP_AND_POSITIVE_ALL_CLASS_FIXED_SHELLS"
        elif ambient["label"] == "12A2":
            tier = 1
            reason = (
                "The intrinsic ternary Golay glue has full monomial residual "
                "group 2.M12 of order 190080 and component image M12. Every "
                "eligible fixed-coordinate shell gives 214 residual orbits "
                "and 99 local K3 surface classes."
            )
            requested = [
                "12A2-C02 involution stabilizers",
                "12A2-C05 order-three and C07 order-four stabilizers",
                "12A2-C10 order-five and C11 order-six stabilizers",
                "non-coordinate invariant auxiliaries",
            ]
            evidence = "EXACT_TERNARY_GOLAY_RESIDUAL_GROUP_AND_POSITIVE_ALL_CLASS_FIXED_SHELLS"
        elif ambient["label"] == "3E8":
            tier = 2
            reason = (
                "The exact S3 chamber residual group permutes the three E8 "
                "components. Its complete nonidentity fixed-coordinate shell "
                "has no determinant/length/MW12--17/mod-2-qualified seed, so "
                "non-coordinate invariant languages remain the target."
            )
            requested = [
                "non-coordinate transposition-invariant auxiliary",
                "non-coordinate 3-cycle-invariant auxiliary",
            ]
            evidence = "EXACT_S3_RESIDUAL_GROUP_AND_NEGATIVE_ALL_CLASS_FIXED_COORDINATE_CONTROL"
        elif ambient["label"] == "3D8":
            tier = 2
            reason = (
                "The exact order-six glue stabilizer is the natural S3 "
                "component-permutation group. Its all-class fixed-coordinate "
                "shell gives 40 mod-2-nontrivial residual orbits and seven "
                "local K3 surface classes."
            )
            requested = [
                "3D8-C02 transposition stabilizers",
                "non-coordinate transposition- and 3-cycle-invariant auxiliaries",
            ]
            evidence = "EXACT_GLUE_RESIDUAL_GROUP_AND_POSITIVE_ALL_CLASS_FIXED_SHELLS"
        elif ambient["label"] == "3A8":
            tier = 2
            reason = (
                "The exact order-twelve glue stabilizer is {+/-1} times S3. "
                "Its all-class fixed-coordinate shell gives 189 mod-2-"
                "nontrivial residual orbits and 25 local K3 surface classes."
            )
            requested = [
                "3A8-C02 transposition stabilizers",
                "non-coordinate invariant auxiliaries for the negative central, signed-transposition, and 3-cycle classes",
            ]
            evidence = "EXACT_GLUE_RESIDUAL_GROUP_AND_POSITIVE_ALL_CLASS_FIXED_SHELLS"
        elif ambient["label"] == "2D12":
            tier = 2
            reason = (
                "The exact order-two glue stabilizer is the natural component "
                "swap. Its complete fixed-coordinate shell has no "
                "determinant/length/MW12--17/mod-2-qualified seed, so "
                "non-coordinate invariant languages remain the target."
            )
            requested = ["non-coordinate component-swap-invariant auxiliary"]
            evidence = "EXACT_S2_RESIDUAL_GROUP_AND_NEGATIVE_ALL_CLASS_FIXED_COORDINATE_CONTROL"
        elif ambient["label"] == "D10_2E7":
            tier = 2
            reason = (
                "Glue retains exactly the simultaneous D10 diagram involution "
                "and E7-component swap. Its complete fixed-coordinate shell "
                "has no determinant/length/MW12--17/mod-2-qualified seed, so "
                "non-coordinate invariant languages remain the target."
            )
            requested = [
                "non-coordinate auxiliary invariant under the coupled D10/E7 involution"
            ]
            evidence = "EXACT_COUPLED_RESIDUAL_INVOLUTION_AND_NEGATIVE_ALL_CLASS_FIXED_COORDINATE_CONTROL"
        elif ambient["label"] == "2A12":
            tier = 2
            reason = (
                "The exact glue-preserving residual group is cyclic of order "
                "four. Its sole eligible order-two fixed-coordinate shell has "
                "no determinant/length/MW12--17/mod-2-qualified seed, so "
                "non-coordinate invariant languages remain the target."
            )
            requested = [
                "non-coordinate auxiliary invariant under the order-two diagram action"
            ]
            evidence = "EXACT_C4_RESIDUAL_GROUP_AND_NEGATIVE_ALL_ELIGIBLE_CLASS_FIXED_COORDINATE_CONTROL"
        elif ambient["label"] == "2A9_D6":
            tier = 2
            reason = (
                "The exact glue-preserving residual group is cyclic of order "
                "four. Its order-four fixed-coordinate shells give 32 residual "
                "orbits and five local K3 surface classes; the order-two shell "
                "is negative."
            )
            requested = [
                "2A9_D6-C03/C04 order-four stabilizers",
                "non-coordinate invariant auxiliaries",
            ]
            evidence = "EXACT_C4_RESIDUAL_GROUP_AND_POSITIVE_ALL_CLASS_FIXED_SHELLS"
        elif ambient["label"] in eta_group_by_label:
            eta_group = eta_group_by_label[ambient["label"]]
            eta_probe = eta_probe_by_label[ambient["label"]]
            tier = 3
            if eta_group["residual_group"]["order"] == 1:
                reason = (
                    "The exhaustive component Dynkin-diagram lift test proves "
                    "that the chamber residual group is trivial. There is no "
                    "nonidentity residual fixed-coordinate language; systematic "
                    "non-symmetry auxiliary enumeration remains open."
                )
                requested = ["systematic primitive rank-seven auxiliary enumeration"]
                evidence = "EXACT_TRIVIAL_ETA_RESIDUAL_GROUP"
            else:
                assert eta_probe["accounting"]["conjugacy_classes_scanned"] == 1
                assert eta_probe["accounting"][
                    "high_mw_mod2_accepted_seeds_before_residual_dedup"
                ] == 0
                reason = (
                    "The complete residual group is generated by the eta "
                    "diagram involution. Its sole nonidentity fixed-coordinate "
                    "shell has no determinant/length/MW12--17/mod-2-qualified "
                    "seed, so non-coordinate invariant languages remain open."
                )
                requested = ["non-coordinate eta-invariant auxiliary"]
                evidence = "EXACT_ETA_RESIDUAL_GROUP_AND_NEGATIVE_FIXED_COORDINATE_CONTROL"
        elif four_cycle_types:
            tier = 1
            reason = (
                "At least four identical root components give the permutation "
                "envelope order-four and non-scalar involution candidates."
            )
            requested = [
                "order-4 component permutation",
                "non-scalar order-2 component permutation",
            ]
            evidence = "COMPONENT_ENVELOPE_HEURISTIC_ACTUAL_GX_AND_STABILIZER_OPEN"
        elif swap_types:
            tier = 2
            reason = (
                "Repeated identical root components give non-scalar swap candidates."
            )
            requested = ["non-scalar order-2 component permutation"]
            evidence = "COMPONENT_ENVELOPE_HEURISTIC_ACTUAL_GX_AND_STABILIZER_OPEN"
        else:
            tier = 3
            reason = (
                "No repeated root-component type; component permutations are not "
                "the first mod-2 symmetry route."
            )
            requested = []
            evidence = "LOW_COMPONENT_PERMUTATION_PRIORITY"
        rows.append(
            {
                "backend_id": f"ROOTED-{ambient['label']}",
                "ambient_label": ambient["label"],
                "priority_tier": tier,
                "component_type_multiplicities": dict(sorted(multiplicities.items())),
                "repeated_component_types": repeated,
                "full_symmetric_component_permutation_envelope_order": envelope_order,
                "swap_component_types": swap_types,
                "four_cycle_component_types": four_cycle_types,
                "requested_stabilizer_action_types": requested,
                "reason": reason,
                "evidence_status": evidence,
                "acceptance_gate": {
                    "full_ambient_stabilizer_required": True,
                    "chamber_section_intersection_alone_is_insufficient": True,
                    "induced_complement_action_required": True,
                    "mod2_moved_dimension_test": "rank_GF2(g_M - I) > 0",
                    "then_compute": (
                        "fixed-point and orbit distribution on the rational subset of M/2M"
                    ),
                },
            }
        )
    rows.sort(key=lambda row: (row["priority_tier"], row["ambient_label"]))
    tier_distribution = Counter(row["priority_tier"] for row in rows)
    assert len(rows) == 23
    return {
        "schema": "elkies-k3.cross-niemeier-mod2-priority.v1",
        "status": "PASS_EXACT_PRIORITY_LEDGER_HEURISTIC_BACKEND_ORDER",
        "proof_scope": {
            "proved": (
                "The root-component multiplicities of all 23 rooted Niemeier "
                "backends are classified into a deterministic experiment order. "
                "For A7^2 D5^2, the exact Dih_4 section and six-frame negative "
                "control are imported and checked. For 4D6 and 6A4, exact "
                "component sections and their first bounded fixed shells are checked. "
                "For 4A5+D4, the exact residual group and first order-four fixed "
                "shell are checked. For 4A6 and 4E6, every eligible residual "
                "matrix-conjugacy-class coordinate shell is checked. The same "
                "all-class check is complete for the declared 6D4 and 8A3 "
                "coordinate languages. The complete nonidentity 3E8 "
                "fixed-coordinate language is an exact negative control. The "
                "same all-class coordinate check is positive for 3D8, while "
                "the complete 2D12 and D10+2E7 involution-coordinate languages "
                "and the sole eligible 2A12 class language are negative; the "
                "all-class 2A9+D6, 3A8, and 12A2 checks are positive. The six "
                "eta-only residual groups are exact; D24 and D16+E8 are trivial, "
                "and all four nontrivial fixed-coordinate controls are negative."
            ),
            "not_proved": (
                "A repeated-component envelope does not prove that its permutations "
                "preserve Niemeier glue, stabilize a future auxiliary, or act "
                "nontrivially on M/2M. Every retained experiment must pass the "
                "declared full-stabilizer and induced-action gates."
            ),
        },
        "selection_policy": {
            "primary_signal": (
                "full stabilizer contains a non-scalar component permutation"
            ),
            "required_exact_signal": "rank_GF2(g_M-I)>0 on the complement",
            "seed_classes_for_A7_2_D5_2": ["2B", "2C", "4A"],
            "negative_control_classes": ["1A", "2A"],
            "use_after_gate": (
                "prioritize orbit-resolved rational bisection cosets and source searches"
            ),
        },
        "exact_A7_2_D5_2_control": {
            "section_class_distribution": dict(sorted(exact_section_classes.items())),
            "observed_full_stabilizer_image_classes": observed_classes,
            "frames_checked": len(umbral["targets"]),
            "mod2_result": (
                "all observed complement images are generated by +/-I and act trivially on M/2M"
            ),
        },
        "exact_4D6_control": {
            "section": "lifted S4 component section",
            "coordinate_subsets_tested": four_d6_swap["accounting"][
                "coordinate_subsets_tested"
            ],
            "high_mw_mod2_accepted_seeds": 0,
            "scope": "selected component-transposition fixed-lattice coordinate shell only",
        },
        "exact_6A4_control": {
            "residual_group_order": six_a4_double_swap["residual_group"][
                "order"
            ],
            "component_permutation_image_order": six_a4_double_swap[
                "residual_group"
            ]["component_permutation_image_order"],
            "selected_class": six_a4_double_swap["parameters"][
                "selected_class"
            ],
            "high_mw_mod2_accepted_seeds": six_a4_double_swap["accounting"][
                "high_mw_mod2_accepted_seeds"
            ],
            "residual_group_embedding_orbits": six_a4_double_swap["accounting"][
                "residual_group_embedding_orbits"
            ],
            "surface_classes_after_T_NS_first_dedup": six_a4_double_swap[
                "accounting"
            ]["surface_classes_after_T_NS_first_dedup"],
        },
        "exact_4A5_D4_control": {
            "residual_group_order": four_a5_d4_order4["residual_group"][
                "order"
            ],
            "component_permutation_image_order": four_a5_d4_order4[
                "residual_group"
            ]["component_permutation_image_order"],
            "selected_class": four_a5_d4_order4["parameters"][
                "selected_class"
            ],
            "high_mw_mod2_accepted_seeds": four_a5_d4_order4["accounting"][
                "high_mw_mod2_accepted_seeds"
            ],
            "residual_group_embedding_orbits": four_a5_d4_order4["accounting"][
                "residual_group_embedding_orbits"
            ],
            "surface_classes_after_T_NS_first_dedup": four_a5_d4_order4[
                "accounting"
            ]["surface_classes_after_T_NS_first_dedup"],
        },
        "exact_4A6_4E6_control": {
            label: {
                "residual_group_order": backend["residual_group"]["order"],
                "component_permutation_image_order": backend[
                    "residual_group"
                ]["component_permutation_image_order"],
                "coordinate_subsets_tested": backend[
                    "source_probe_accounting"
                ]["coordinate_subsets_tested"],
                "high_mw_mod2_accepted_seeds": backend["accounting"][
                    "high_mw_mod2_accepted_seeds_before_residual_dedup"
                ],
                "residual_group_embedding_orbits": backend["accounting"][
                    "residual_group_embedding_orbits"
                ],
                "surface_classes_after_T_NS_first_dedup": backend[
                    "accounting"
                ]["surface_classes_after_T_NS_first_dedup"],
                "positive_stabilizer_classes": sorted(
                    backend["accounting"][
                        "nontrivial_mod2_stabilizer_class_coverage"
                    ]
                ),
            }
            for label, backend in sorted(a6_e6_by_label.items())
        },
        "exact_8A3_control": {
            "residual_group_order": eight_a3_backend["residual_group"]["order"],
            "component_permutation_image_order": eight_a3_backend[
                "residual_group"
            ]["component_permutation_image_order"],
            "coordinate_subsets_tested": eight_a3_backend[
                "source_probe_accounting"
            ]["coordinate_subsets_tested"],
            "high_mw_mod2_accepted_seeds": eight_a3_backend["accounting"][
                "high_mw_mod2_accepted_seeds_before_residual_dedup"
            ],
            "residual_group_embedding_orbits": eight_a3_backend["accounting"][
                "residual_group_embedding_orbits"
            ],
            "surface_classes_after_T_NS_first_dedup": eight_a3_backend[
                "accounting"
            ]["surface_classes_after_T_NS_first_dedup"],
            "positive_stabilizer_classes": sorted(
                eight_a3_backend["accounting"][
                    "nontrivial_mod2_stabilizer_class_coverage"
                ]
            ),
        },
        "exact_6D4_control": {
            "residual_group_order": six_d4_backend["residual_group"]["order"],
            "component_permutation_image_order": six_d4_backend[
                "residual_group"
            ]["component_permutation_image_order"],
            "component_kernel_order": six_d4_backend["residual_group"][
                "component_kernel_order"
            ],
            "coordinate_subsets_tested": six_d4_backend[
                "source_probe_accounting"
            ]["coordinate_subsets_tested"],
            "high_mw_mod2_accepted_seeds": six_d4_backend["accounting"][
                "high_mw_mod2_accepted_seeds_before_residual_dedup"
            ],
            "residual_group_embedding_orbits": six_d4_backend["accounting"][
                "residual_group_embedding_orbits"
            ],
            "surface_classes_after_T_NS_first_dedup": six_d4_backend[
                "accounting"
            ]["surface_classes_after_T_NS_first_dedup"],
            "positive_stabilizer_classes": sorted(
                six_d4_backend["accounting"][
                    "nontrivial_mod2_stabilizer_class_coverage"
                ]
            ),
        },
        "exact_3E8_control": {
            "residual_group_order": three_e8_fixed_shell_probe[
                "residual_group_order"
            ],
            "component_permutation_image_order": three_e8_fixed_shell_probe[
                "component_permutation_image_order"
            ],
            "conjugacy_classes_scanned": three_e8_fixed_shell_probe[
                "accounting"
            ]["conjugacy_classes_scanned"],
            "coordinate_subsets_tested": three_e8_fixed_shell_probe[
                "accounting"
            ]["coordinate_subsets_tested"],
            "high_mw_mod2_accepted_seeds": three_e8_fixed_shell_probe[
                "accounting"
            ]["high_mw_mod2_accepted_seeds_before_residual_dedup"],
            "scope": (
                "every 7-of-fixed-rank coordinate summand for both "
                "nonidentity residual matrix conjugacy classes"
            ),
        },
        "exact_3D8_control": {
            "residual_group_order": three_d8_backend["residual_group"]["order"],
            "component_permutation_image_order": three_d8_backend[
                "residual_group"
            ]["component_permutation_image_order"],
            "coordinate_subsets_tested": three_d8_backend[
                "source_probe_accounting"
            ]["coordinate_subsets_tested"],
            "high_mw_mod2_accepted_seeds": three_d8_backend["accounting"][
                "high_mw_mod2_accepted_seeds_before_residual_dedup"
            ],
            "residual_group_embedding_orbits": three_d8_backend["accounting"][
                "residual_group_embedding_orbits"
            ],
            "surface_classes_after_T_NS_first_dedup": three_d8_backend[
                "accounting"
            ]["surface_classes_after_T_NS_first_dedup"],
            "positive_stabilizer_classes": sorted(
                three_d8_backend["accounting"][
                    "nontrivial_mod2_stabilizer_class_coverage"
                ]
            ),
        },
        "exact_2D12_control": {
            "residual_group_order": two_d12_fixed_shell_probe[
                "residual_group_order"
            ],
            "component_permutation_image_order": two_d12_fixed_shell_probe[
                "component_permutation_image_order"
            ],
            "conjugacy_classes_scanned": two_d12_fixed_shell_probe[
                "accounting"
            ]["conjugacy_classes_scanned"],
            "coordinate_subsets_tested": two_d12_fixed_shell_probe[
                "accounting"
            ]["coordinate_subsets_tested"],
            "high_mw_mod2_accepted_seeds": two_d12_fixed_shell_probe[
                "accounting"
            ]["high_mw_mod2_accepted_seeds_before_residual_dedup"],
            "scope": (
                "all 7-of-12 coordinate summands of the component-swap "
                "fixed-lattice pinned LLL basis"
            ),
        },
        "exact_D10_2E7_control": {
            "residual_group_order": d10_two_e7_fixed_shell_probe[
                "residual_group_order"
            ],
            "d10_diagram_image_order": d10_two_e7_fixed_shell_probe[
                "d10_diagram_image_order"
            ],
            "e7_component_permutation_image_order": (
                d10_two_e7_fixed_shell_probe[
                    "e7_component_permutation_image_order"
                ]
            ),
            "conjugacy_classes_scanned": d10_two_e7_fixed_shell_probe[
                "accounting"
            ]["conjugacy_classes_scanned"],
            "coordinate_subsets_tested": d10_two_e7_fixed_shell_probe[
                "accounting"
            ]["coordinate_subsets_tested"],
            "high_mw_mod2_accepted_seeds": d10_two_e7_fixed_shell_probe[
                "accounting"
            ]["high_mw_mod2_accepted_seeds_before_residual_dedup"],
            "scope": (
                "all 7-of-16 coordinate summands of the coupled-residual-"
                "involution fixed-lattice pinned LLL basis"
            ),
        },
        "exact_2A12_control": {
            "residual_group_order": two_a12_fixed_shell_probe[
                "residual_group_order"
            ],
            "component_permutation_image_order": two_a12_fixed_shell_probe[
                "component_permutation_image_order"
            ],
            "component_diagram_kernel_order": two_a12_fixed_shell_probe[
                "component_diagram_kernel_order"
            ],
            "conjugacy_classes_scanned": two_a12_fixed_shell_probe[
                "accounting"
            ]["conjugacy_classes_scanned"],
            "coordinate_subsets_tested": two_a12_fixed_shell_probe[
                "accounting"
            ]["coordinate_subsets_tested"],
            "high_mw_mod2_accepted_seeds": two_a12_fixed_shell_probe[
                "accounting"
            ]["high_mw_mod2_accepted_seeds_before_residual_dedup"],
            "scope": (
                "all 7-of-12 coordinate summands for the sole nonidentity "
                "residual class of fixed rank at least seven"
            ),
        },
        "exact_2A9_D6_control": {
            "residual_group_order": two_a9_d6_backend["residual_group"]["order"],
            "a9_component_permutation_image_order": two_a9_d6_backend[
                "residual_group"
            ]["a9_component_permutation_image_order"],
            "coordinate_subsets_tested": two_a9_d6_backend[
                "source_probe_accounting"
            ]["coordinate_subsets_tested"],
            "high_mw_mod2_accepted_seeds": two_a9_d6_backend["accounting"][
                "high_mw_mod2_accepted_seeds_before_residual_dedup"
            ],
            "residual_group_embedding_orbits": two_a9_d6_backend["accounting"][
                "residual_group_embedding_orbits"
            ],
            "surface_classes_after_T_NS_first_dedup": two_a9_d6_backend[
                "accounting"
            ]["surface_classes_after_T_NS_first_dedup"],
            "positive_stabilizer_classes": sorted(
                two_a9_d6_backend["accounting"][
                    "nontrivial_mod2_stabilizer_class_coverage"
                ]
            ),
        },
        "exact_3A8_control": {
            "residual_group_order": three_a8_backend["residual_group"]["order"],
            "component_permutation_image_order": three_a8_backend[
                "residual_group"
            ]["component_permutation_image_order"],
            "coordinate_subsets_tested": three_a8_backend[
                "source_probe_accounting"
            ]["coordinate_subsets_tested"],
            "high_mw_mod2_accepted_seeds": three_a8_backend["accounting"][
                "high_mw_mod2_accepted_seeds_before_residual_dedup"
            ],
            "residual_group_embedding_orbits": three_a8_backend["accounting"][
                "residual_group_embedding_orbits"
            ],
            "surface_classes_after_T_NS_first_dedup": three_a8_backend[
                "accounting"
            ]["surface_classes_after_T_NS_first_dedup"],
            "positive_stabilizer_classes": sorted(
                three_a8_backend["accounting"][
                    "nontrivial_mod2_stabilizer_class_coverage"
                ]
            ),
        },
        "exact_12A2_control": {
            "residual_group_order": twelve_a2_backend["residual_group"]["order"],
            "component_permutation_image_order": twelve_a2_backend[
                "residual_group"
            ]["component_permutation_image_order"],
            "central_diagram_kernel_order": twelve_a2_backend[
                "residual_group"
            ]["central_diagram_kernel_order"],
            "coordinate_subsets_tested": twelve_a2_backend[
                "source_probe_accounting"
            ]["coordinate_subsets_tested"],
            "high_mw_mod2_accepted_seeds": twelve_a2_backend["accounting"][
                "high_mw_mod2_accepted_seeds_before_residual_dedup"
            ],
            "residual_group_embedding_orbits": twelve_a2_backend["accounting"][
                "residual_group_embedding_orbits"
            ],
            "surface_classes_after_T_NS_first_dedup": twelve_a2_backend[
                "accounting"
            ]["surface_classes_after_T_NS_first_dedup"],
            "positive_stabilizer_classes": sorted(
                twelve_a2_backend["accounting"][
                    "nontrivial_mod2_stabilizer_class_coverage"
                ]
            ),
        },
        "exact_eta_only_controls": {
            label: {
                "residual_group_order": eta_group_by_label[label][
                    "residual_group"
                ]["order"],
                "component_diagram_kernel_order": eta_group_by_label[label][
                    "residual_group"
                ]["component_diagram_kernel_order"],
                "root_lattice_index": eta_group_by_label[label]["root_lattice"][
                    "index_in_niemeier"
                ],
                "candidate_chamber_maps_tested": eta_group_by_label[label][
                    "residual_group"
                ]["candidate_maps_tested"],
                "conjugacy_classes_scanned": eta_probe_by_label[label][
                    "accounting"
                ]["conjugacy_classes_scanned"],
                "coordinate_subsets_tested": eta_probe_by_label[label][
                    "accounting"
                ]["coordinate_subsets_tested"],
                "high_mw_mod2_accepted_seeds": eta_probe_by_label[label][
                    "accounting"
                ]["high_mw_mod2_accepted_seeds_before_residual_dedup"],
            }
            for label in sorted(eta_group_by_label)
        },
        "accounting": {
            "rooted_backends": len(rows),
            "priority_tier_distribution": {
                str(key): value for key, value in sorted(tier_distribution.items())
            },
        },
        "backends": rows,
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--niemeier", type=Path, default=DEFAULT_NIEMEIER)
parser.add_argument("--umbral", type=Path, default=DEFAULT_UMBRAL)
parser.add_argument("--four-d6-swap", type=Path, default=DEFAULT_4D6_SWAP)
parser.add_argument(
    "--six-a4-double-swap", type=Path, default=DEFAULT_6A4_DOUBLE_SWAP
)
parser.add_argument(
    "--four-a5-d4-order4", type=Path, default=DEFAULT_4A5_D4_ORDER4
)
parser.add_argument(
    "--four-a6-four-e6-fixed-shells",
    type=Path,
    default=DEFAULT_4A6_4E6_FIXED_SHELLS,
)
parser.add_argument(
    "--eight-a3-fixed-shells", type=Path, default=DEFAULT_8A3_FIXED_SHELLS
)
parser.add_argument(
    "--six-d4-fixed-shells", type=Path, default=DEFAULT_6D4_FIXED_SHELLS
)
parser.add_argument(
    "--three-e8-fixed-shell-probe",
    type=Path,
    default=DEFAULT_3E8_FIXED_SHELL_PROBE,
)
parser.add_argument(
    "--three-d8-fixed-shells",
    type=Path,
    default=DEFAULT_3D8_FIXED_SHELLS,
)
parser.add_argument(
    "--two-d12-fixed-shell-probe",
    type=Path,
    default=DEFAULT_2D12_FIXED_SHELL_PROBE,
)
parser.add_argument(
    "--d10-two-e7-fixed-shell-probe",
    type=Path,
    default=DEFAULT_D10_2E7_FIXED_SHELL_PROBE,
)
parser.add_argument(
    "--two-a12-fixed-shell-probe",
    type=Path,
    default=DEFAULT_2A12_FIXED_SHELL_PROBE,
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
parser.add_argument(
    "--eta-only-residual-groups",
    type=Path,
    default=DEFAULT_ETA_ONLY_RESIDUAL_GROUPS,
)
parser.add_argument(
    "--eta-only-fixed-shell-probe",
    type=Path,
    default=DEFAULT_ETA_ONLY_FIXED_SHELL_PROBE,
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
result = build(
    json.loads(arguments.niemeier.read_text()),
    json.loads(arguments.umbral.read_text()),
    json.loads(arguments.four_d6_swap.read_text()),
    json.loads(arguments.six_a4_double_swap.read_text()),
    json.loads(arguments.four_a5_d4_order4.read_text()),
    json.loads(arguments.four_a6_four_e6_fixed_shells.read_text()),
    json.loads(arguments.eight_a3_fixed_shells.read_text()),
    json.loads(arguments.six_d4_fixed_shells.read_text()),
    json.loads(arguments.three_e8_fixed_shell_probe.read_text()),
    json.loads(arguments.three_d8_fixed_shells.read_text()),
    json.loads(arguments.two_d12_fixed_shell_probe.read_text()),
    json.loads(arguments.d10_two_e7_fixed_shell_probe.read_text()),
    json.loads(arguments.two_a12_fixed_shell_probe.read_text()),
    json.loads(arguments.two_a9_d6_fixed_shells.read_text()),
    json.loads(arguments.three_a8_fixed_shells.read_text()),
    json.loads(arguments.twelve_a2_fixed_shells.read_text()),
    json.loads(arguments.eta_only_residual_groups.read_text()),
    json.loads(arguments.eta_only_fixed_shell_probe.read_text()),
)
result["inputs"] = {
    str(arguments.niemeier.relative_to(ROOT)): digest(arguments.niemeier),
    str(arguments.umbral.relative_to(ROOT)): digest(arguments.umbral),
    str(arguments.four_d6_swap.relative_to(ROOT)): digest(
        arguments.four_d6_swap
    ),
    str(arguments.six_a4_double_swap.relative_to(ROOT)): digest(
        arguments.six_a4_double_swap
    ),
    str(arguments.four_a5_d4_order4.relative_to(ROOT)): digest(
        arguments.four_a5_d4_order4
    ),
    str(arguments.four_a6_four_e6_fixed_shells.relative_to(ROOT)): digest(
        arguments.four_a6_four_e6_fixed_shells
    ),
    str(arguments.eight_a3_fixed_shells.relative_to(ROOT)): digest(
        arguments.eight_a3_fixed_shells
    ),
    str(arguments.six_d4_fixed_shells.relative_to(ROOT)): digest(
        arguments.six_d4_fixed_shells
    ),
    str(arguments.three_e8_fixed_shell_probe.relative_to(ROOT)): digest(
        arguments.three_e8_fixed_shell_probe
    ),
    str(arguments.three_d8_fixed_shells.relative_to(ROOT)): digest(
        arguments.three_d8_fixed_shells
    ),
    str(arguments.two_d12_fixed_shell_probe.relative_to(ROOT)): digest(
        arguments.two_d12_fixed_shell_probe
    ),
    str(arguments.d10_two_e7_fixed_shell_probe.relative_to(ROOT)): digest(
        arguments.d10_two_e7_fixed_shell_probe
    ),
    str(arguments.two_a12_fixed_shell_probe.relative_to(ROOT)): digest(
        arguments.two_a12_fixed_shell_probe
    ),
    str(arguments.two_a9_d6_fixed_shells.relative_to(ROOT)): digest(
        arguments.two_a9_d6_fixed_shells
    ),
    str(arguments.three_a8_fixed_shells.relative_to(ROOT)): digest(
        arguments.three_a8_fixed_shells
    ),
    str(arguments.twelve_a2_fixed_shells.relative_to(ROOT)): digest(
        arguments.twelve_a2_fixed_shells
    ),
    str(arguments.eta_only_residual_groups.relative_to(ROOT)): digest(
        arguments.eta_only_residual_groups
    ),
    str(arguments.eta_only_fixed_shell_probe.relative_to(ROOT)): digest(
        arguments.eta_only_fixed_shell_probe
    ),
}
encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not arguments.output.exists() or arguments.output.read_text() != encoded:
        raise SystemExit("cross-Niemeier mod-2 priority artifact is stale")
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded)
print(
    "MOD2PRIORITY|backends={}|tiers={}|seed=2B,2C,4A|status=PASS".format(
        result["accounting"]["rooted_backends"],
        result["accounting"]["priority_tier_distribution"],
    )
)
