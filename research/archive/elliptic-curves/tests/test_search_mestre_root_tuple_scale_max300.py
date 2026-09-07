#!/usr/bin/env python3
"""Focused checks for the standalone exact max-root-300 Mestre continuation."""

from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
GENERATED = ROOT / "artifacts" / "generated-results"
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))

from search_mestre_root_tuple_scale import tuple_digest  # noqa: E402


SOURCE = CAS / "enumerate_mestre_root_tuples_scale_max300.cpp"
SCRIPT = CAS / "search_mestre_root_tuple_scale_max300.py"
CENSUS = GENERATED / "elliptic_mestre_root_tuple_scale_max300_census.json"
ARTIFACT = GENERATED / "elliptic_mestre_root_tuple_scale_max300.json"
CERTIFICATE_SCRIPT = CAS / "certify_mestre_02136217261290_t2_rank15.py"
CERTIFICATE = (
    GENERATED / "elliptic_mestre_02136217261290_t2_rank15_certificate.json"
)
EXPECTED_SOURCE_SHA256 = (
    "9b5c939638ca3b1193088434f3c5c386915d049068c9d6c0d13febd059a02d7d"
)
EXPECTED_SCRIPT_SHA256 = (
    "922cd33621e882dbb5483b041c547f568d3f4fdfea2bffc8cab0d3741a3445b4"
)
EXPECTED_CENSUS_SHA256 = (
    "c5a68905977f059182efc1233e7301c039b2164b45bdfef1f8fd106b13d263ea"
)
EXPECTED_ARTIFACT_SHA256 = (
    "e4a7be774ac0cae3c636c70bde7490e7ced7e313971dfba3c9017e48d730fca7"
)
EXPECTED_CERTIFICATE_SCRIPT_SHA256 = (
    "bacc7f735c86d81adaac1bb0d5c1916f8120f79b2f12f33c2b720418358bc06e"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "35abefefab42b19f49fad074f0c2cd65b039e8f36c398fbe7b46f68a0c2f09ea"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class MestreRootTupleScaleMax300Test(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.census = json.loads(CENSUS.read_text())
        cls.data = json.loads(ARTIFACT.read_text())
        cls.certificate = json.loads(CERTIFICATE.read_text())

    def test_new_files_and_artifacts_are_pinned(self) -> None:
        self.assertEqual(sha256(SOURCE), EXPECTED_SOURCE_SHA256)
        self.assertEqual(sha256(SCRIPT), EXPECTED_SCRIPT_SHA256)
        self.assertEqual(sha256(CENSUS), EXPECTED_CENSUS_SHA256)
        self.assertEqual(sha256(ARTIFACT), EXPECTED_ARTIFACT_SHA256)
        self.assertEqual(
            sha256(CERTIFICATE_SCRIPT), EXPECTED_CERTIFICATE_SCRIPT_SHA256
        )
        self.assertEqual(sha256(CERTIFICATE), EXPECTED_CERTIFICATE_SHA256)
        self.assertEqual(
            self.census["result_sha256"],
            "a15f9c65701676ecfdeb4ef6473393254f98aa94b30c9ee70b491f6df11cc5f7",
        )
        self.assertEqual(
            self.data["result_sha256"],
            "a567485ccc11cf56a11ddf83a32c809b5fd2670a820bdcf3f13a754d75b669ce",
        )

    def test_exact_census_and_record_for_record_max200_prefix(self) -> None:
        census = self.census["census"]
        self.assertEqual(
            census["affine_normalized_primitive_reflection_quotient_count"],
            9_456_367_899,
        )
        self.assertEqual(census["degree_five_obstruction_zero_count"], 931_824)
        self.assertEqual(census["reflection_obstruction_zero_count"], 922_407)
        self.assertEqual(census["nonreflection_obstruction_zero_count"], 9_417)
        self.assertEqual(census["nonreflection_generically_nonsingular_count"], 2_329)
        self.assertEqual(census["nonreflection_generically_singular_count"], 7_088)
        self.assertEqual(
            census["genuinely_new_diameter_201_to_300_family_count"], 1_291
        )
        self.assertEqual(
            census["genuinely_new_diameter_201_to_300_family_sha256"],
            "5bf617b68c8aea2a78202ffc75418cb918b498f53ea121caa4ebb70a99f7008e",
        )
        self.assertEqual(
            census["independent_burnside_mobius_count"],
            {
                "primitive_unquotiented_count": 18_911_813_391,
                "primitive_reflection_fixed_count": 922_407,
                "primitive_reflection_orbit_count": 9_456_367_899,
            },
        )
        prefix = self.census["exact_max200_prefix_recovery"]
        self.assertEqual(prefix["counts"], [1_225_592_478, 275_387, 271_600, 3_787])
        self.assertTrue(prefix["record_for_record_equal_to_frozen_enumerator"])
        self.assertTrue(prefix["nonreflection_records_equal_frozen_census"])
        new_roots = tuple(
            tuple(roots)
            for roots in self.census["tuple_populations"][
                "genuinely_new_nonsingular_roots"
            ]
        )
        self.assertEqual(len(new_roots), 1_291)
        self.assertTrue(all(201 <= roots[-1] <= 300 for roots in new_roots))
        self.assertEqual(
            tuple_digest(new_roots),
            "5bf617b68c8aea2a78202ffc75418cb918b498f53ea121caa4ebb70a99f7008e",
        )

    def test_complete_panel_and_fixed_leader_population(self) -> None:
        population = self.data["complete_panel_screen"]["population"]
        self.assertEqual(population["new_family_count"], 1_291)
        self.assertEqual(population["proposed_panel_fiber_count"], 10_328)
        self.assertEqual(population["admissible_panel_fiber_count"], 9_563)
        self.assertEqual(population["inadmissible_panel_fiber_count"], 765)
        self.assertEqual(population["exact_mod3_certificate_count"], 9_563)
        self.assertEqual(population["maximum_visible_certified_rank_lower_bound"], 11)
        self.assertEqual(
            population["visible_rank_lower_bound_histogram"],
            {"5": 6, "6": 158, "7": 882, "8": 45, "9": 845, "10": 5_619, "11": 2_008},
        )
        selection = self.data["rank_aware_diversity_selection"]
        self.assertTrue(selection["population_closed_before_conductor_calls"])
        self.assertEqual(selection["selected_family_count"], 64)
        self.assertEqual(
            selection["selected_identifier_sha256"],
            "a3ad2e369f61e2a0e6a47c04266cd7c1885a2b82e06f3b0676ae4f765bdf6ae0",
        )
        self.assertEqual(
            Counter(selection["diversity_decile_counts"].values()), Counter({3: 10})
        )
        leaders = self.data["leader_followup"]["population"]
        self.assertEqual(leaders["selected_leaders"], 64)
        self.assertEqual(leaders["conductor_completed"], 64)
        self.assertEqual(leaders["point_search_completed"], 64)
        self.assertEqual(leaders["subtarget_conductors"], 64)
        self.assertEqual(leaders["maximum_stable_numerical_rank"], 15)
        self.assertEqual(
            leaders["stable_numerical_rank_histogram"],
            {"11": 20, "12": 24, "13": 17, "14": 2, "15": 1},
        )

    def test_standalone_exact_rank15_certificate(self) -> None:
        certificate = self.certificate
        self.assertEqual(
            certificate["status"], "certified exact algebraic rank lower bound 15"
        )
        self.assertEqual(certificate["curve"]["family_roots"], [0, 2, 136, 217, 261, 290])
        self.assertEqual(certificate["curve"]["parameter_T"], "2")
        self.assertEqual(
            certificate["curve"]["conductor"],
            "27535464408096664363840114552671696329422686810",
        )
        self.assertEqual(certificate["curve"]["root_number"], -1)
        self.assertEqual(
            certificate["input_selection"]["point_sha256"],
            "43063eff53a2764cda3b950ee257e87cac9f498be3025be582b742ca7303583c",
        )
        exact = certificate["exact_finite_reduction_certificate"]
        self.assertEqual(exact["certified_algebraic_rank_lower_bound"], 15)
        self.assertEqual(exact["combined_exact_rank_over_F3"], 15)
        self.assertEqual(
            exact["certificate_primes"],
            [17, 23, 29, 37, 41, 43, 83, 101, 103, 107, 131, 149, 173, 193, 199],
        )
        self.assertEqual(
            certificate["result_sha256"],
            "c1de5071cf9ac8bb993345804bb0ab6f96656c72912c294f8e5fe097d002a77b",
        )


if __name__ == "__main__":
    unittest.main()
