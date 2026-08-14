from __future__ import annotations

import importlib.util
import json
from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves" / "cas"
SCRIPT = CAS / "search_mestre_rank15_490_9_outside_neighborhood.py"
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "elliptic_mestre_rank15_490_9_outside_neighborhood.json"
)
sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location(
    "search_mestre_rank15_490_9_outside_neighborhood", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class MestreRank15OutsideNeighborhoodTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.anchor = MODULE.validate_pinned_anchor(ROOT)
        cls.population = MODULE.generate_population()

    def test_anchor_and_complete_prior_box_are_exactly_separated(self) -> None:
        self.assertEqual(tuple(self.anchor["curve"]["roots"]), MODULE.ROOTS)
        self.assertEqual(Q(self.anchor["curve"]["parameter"]), MODULE.T0)
        self.assertEqual(
            self.anchor["claim"]["certified_algebraic_rank_lower_bound"], 15
        )
        self.assertEqual(len(self.population), MODULE.EXPECTED_GENERATOR_COUNT)
        self.assertEqual(
            MODULE.generator_digest(self.population),
            MODULE.EXPECTED_GENERATOR_SHA256,
        )
        self.assertNotIn(MODULE.T0, self.population)
        self.assertTrue(
            all(not MODULE.belongs_to_prior_box(value) for value in self.population)
        )
        self.assertTrue(all(value > 0 for value in self.population))

    def test_farey_and_crt_mutations_replay_exactly(self) -> None:
        self.assertEqual(490 * 7 - 9 * 381, 1)
        self.assertEqual(490 * 2 - 9 * 109, -1)
        source_counts = {
            family: sum(
                any(source.startswith(family) for source in sources)
                for sources in self.population.values()
            )
            for family in (
                "gauss-near",
                "farey-left",
                "farey-right",
                "discriminant-power",
                "local-trace",
            )
        }
        self.assertEqual(
            source_counts,
            {
                "gauss-near": 19_222,
                "farey-left": 1_988,
                "farey-right": 1_987,
                "discriminant-power": 6_495,
                "local-trace": 6_602,
            },
        )
        checked_power = checked_trace = 0
        for parameter, sources in self.population.items():
            for source in sources:
                if source.startswith("discriminant-power:M="):
                    modulus = int(source.split(":", 2)[1].split("=", 1)[1])
                    self.assertEqual(
                        (9 * parameter.numerator - 490 * parameter.denominator)
                        % modulus,
                        0,
                    )
                    checked_power += 1
                elif source.startswith("local-trace:M="):
                    modulus = int(source.split(":", 2)[1].split("=", 1)[1])
                    self.assertEqual(
                        (9 * parameter.numerator - 490 * parameter.denominator)
                        % modulus,
                        0,
                    )
                    checked_trace += 1
        self.assertGreater(checked_power, 6_400)
        self.assertGreater(checked_trace, 6_500)

    def test_generated_artifact_closes_conductor_and_point_stages(self) -> None:
        data = json.loads(ARTIFACT.read_text())
        self.assertEqual(
            data["reproduction"]["script_sha256"], MODULE.sha256_file(SCRIPT)
        )
        self.assertEqual(data["scope"]["generator_count"], 35_507)
        self.assertEqual(
            data["scope"]["generator_sha256"], MODULE.EXPECTED_GENERATOR_SHA256
        )
        self.assertEqual(
            data["population"]["feature_population_sha256"],
            "f2ed41e124aee6c91bc62e25742b6c04904d0ead6fd3b9ea168c90b38fb27e64",
        )
        selection = data["conductor_selection"]
        self.assertTrue(selection["all_proxy_below_190_retained"])
        self.assertEqual(selection["proxy_below_190_count"], 15)
        self.assertEqual(selection["selected_count"], 34)
        self.assertEqual(
            selection["selected_sha256"],
            "7df96e5d8f0fead6c22efe4372183bb7e2cfa24356b8fb1c233531579fdad437",
        )
        conductors = data["conductors"]
        self.assertEqual(len(conductors["records"]), 34)
        self.assertEqual(
            sum(record["status"] == "timeout" for record in conductors["records"].values()),
            1,
        )
        self.assertEqual(conductors["strict_subtarget_count"], 9)
        self.assertEqual(
            set(conductors["strict_subtarget_parameters"]),
            {
                "33044/607",
                "33880/989",
                "34214/625",
                "37114/535",
                "37510/689",
                "39472/725",
                "39970/739",
                "41104/755",
                "56420/1403",
            },
        )
        protocol = data["point_protocol"]
        self.assertEqual(len(protocol["H50000_panel"]), 14)
        self.assertEqual(len(protocol["H250000_panel"]), 6)
        self.assertEqual(protocol["H1000000_panel"], ["37510/689", "39970/739"])
        self.assertTrue(
            all(
                record["stable_numerical_rank"] == 11
                for record in protocol["stage_records"]["H50000"].values()
            )
        )
        self.assertTrue(
            all(
                record["stable_numerical_rank"] == 11
                for record in protocol["stage_records"]["H250000"].values()
            )
        )
        self.assertEqual(
            {
                parameter: record["stable_numerical_rank"]
                for parameter, record in protocol["stage_records"]["H1000000"].items()
            },
            {"37510/689": 12, "39970/739": 11},
        )
        self.assertEqual(protocol["maximum_stable_numerical_rank"], 12)
        self.assertEqual(protocol["certified_signals"], [])
        self.assertFalse(data["outcome"]["rank_signal_at_least_16"])
        self.assertFalse(data["outcome"]["breakthrough_found"])


if __name__ == "__main__":
    unittest.main()
