from __future__ import annotations

from decimal import Decimal
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
SCRIPT = CAS / "search_nagao_section7_anomalous_root_crt.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic_nagao_section7_anomalous_root_crt.json"
)
if str(CAS) not in sys.path:
    sys.path.insert(0, str(CAS))
SPEC = importlib.util.spec_from_file_location("section7_anomalous_root_crt", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SEARCH = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = SEARCH
SPEC.loader.exec_module(SEARCH)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class Section7AnomalousRootCrtTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.population, cls.population_audit = SEARCH.generate_population()
        cls.artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))

    def test_anomalous_roots_and_crt_classes_are_exact(self) -> None:
        self.assertEqual(SEARCH.anomalous_root_audit(), {
            "11": [
                {"residue": 5, "forced_valuation": 3},
                {"residue": 6, "forced_valuation": 3},
            ],
            "19": [
                {"residue": 1, "forced_valuation": 3},
                {"residue": 18, "forced_valuation": 3},
            ],
            "43": [
                {"residue": 1, "forced_valuation": 3},
                {"residue": 42, "forced_valuation": 3},
            ],
            "47": [
                {"residue": 23, "forced_valuation": 3},
                {"residue": 24, "forced_valuation": 3},
            ],
        })
        for primes in SEARCH.profile_primes():
            classes = SEARCH.signless_crt_classes(primes)
            self.assertEqual(len(classes), 2)
            self.assertEqual({modulus for _, modulus in classes}, {primes[0] * primes[1]})

    def test_population_is_new_and_pinned(self) -> None:
        audit = self.population_audit
        self.assertEqual(len(self.population), 80_883)
        self.assertEqual(audit["multi_profile_parameter_count"], 30_321)
        self.assertEqual(
            audit["population_sha256"],
            "b961f0727d2bc2e55939215e2c5144956f72873736e57ca18c71e1a179299e94",
        )
        base_counts = {
            key: value["shell_population_count"]
            for key, value in audit["profiles"].items()
            if not key.endswith("low-height-extension")
        }
        self.assertEqual(base_counts, {
            "root-11-19": 3_674,
            "root-11-43": 10_758,
            "root-11-47": 11_870,
            "root-19-43": 8_447,
            "root-19-47": 14_493,
            "root-43-47": 16_772,
        })
        self.assertTrue(
            all(
                candidate.parameter.denominator > 1_000
                and candidate.height <= 15_000
                for candidate in self.population
            )
        )
        by_parameter = {candidate.parameter: candidate for candidate in self.population}
        expected = {
            Fraction(137, 1022),
            Fraction(875, 1046),
            Fraction(241, 1058),
            Fraction(823, 1210),
            Fraction(2941, 1402),
        }
        self.assertTrue(expected <= set(by_parameter))
        self.assertIn(
            "root-19-43-low-height-extension",
            by_parameter[Fraction(2941, 1402)].profiles,
        )

    def test_complete_prior_population_exclusion_is_exact(self) -> None:
        exclusions, audit = SEARCH.exact_exclusion_population()
        self.assertEqual(len(exclusions), 6_488)
        self.assertEqual(
            audit["finite_outside_box_exclusion_sha256"],
            "4ca379527b3564156f12de3eee45297ed1fd5840b94746267b35d997f413dda2",
        )
        self.assertEqual(
            audit["exact_auxiliary_population"]["exact_union_count"], 6_431
        )
        self.assertEqual(
            audit["exact_auxiliary_population"]["exact_union_sha256"],
            "017d765c220846d4dc943d910f7bb0e401a1c574809531098c43459f61298e17",
        )
        replay = audit["full_neighborhood_deterministic_replay_checkpoint"]
        self.assertEqual(replay["full_count"], 1_574)
        self.assertEqual(replay["outside_global_rectangle_count"], 84)
        self.assertEqual(
            replay["outside_global_rectangle_sha256"],
            "5864d87b437ff465cdbd0f9ea98371f09d99e5fa2b6e48df83d5f4fb132861f7",
        )
        self.assertFalse(
            any(candidate.parameter in exclusions for candidate in self.population)
        )

    def test_artifact_pins_scoring_conductors_and_subtargets(self) -> None:
        data = self.artifact
        self.assertEqual(data["reproduction"]["script_sha256"], sha256_file(SCRIPT))
        self.assertEqual(data["population"]["novel_population_count"], 80_883)
        self.assertEqual(data["population"]["exact_prior_parameter_overlaps_removed"], 0)
        self.assertEqual(
            data["population"]["novel_population_sha256"],
            "b961f0727d2bc2e55939215e2c5144956f72873736e57ca18c71e1a179299e94",
        )
        scoring = data["leakage_free_scoring"]
        self.assertEqual(scoring["exact_discovery_trace_primes"], list(SEARCH.DISCOVERY_PRIMES))
        self.assertEqual(scoring["exact_held_forward_trace_primes"], list(SEARCH.HELD_FORWARD_PRIMES))
        self.assertEqual(scoring["proxy_discovery"]["retained_for_exact_discovery_trace_scoring"], 4_096)
        self.assertEqual(scoring["discovery_frontier_count"], 576)
        self.assertEqual(
            scoring["discovery_frontier_sha256"],
            "a3981d414db8ab2b4f63e37c273caeac56abd359c3a4acf339900a4f0c889567",
        )

        conductor = data["conductor_and_root_parity"]
        self.assertEqual(
            {key: conductor[key] for key in ("selected", "completed", "timeouts", "errors")},
            {"selected": 64, "completed": 61, "timeouts": 3, "errors": 0},
        )
        subtarget = {
            record["constructor_parameter_T"]: (
                Decimal(record["conductor"]["log_conductor"]),
                record["conductor"]["root_number"],
            )
            for record in conductor["records"]
            if record.get("conductor")
            and record["conductor"]["below_strict_log_conductor_target"]
        }
        self.assertEqual(set(subtarget), {
            "137/1022", "2941/1402", "875/1046", "241/1058", "823/1210"
        })
        self.assertEqual(
            {parameter for parameter, (_, root) in subtarget.items() if root == -1},
            {"137/1022", "241/1058"},
        )
        self.assertTrue(all(log_n < Decimal("182.72") for log_n, _ in subtarget.values()))

    def test_staged_point_frontier_is_honestly_negative(self) -> None:
        data = self.artifact
        stages = data["point_stages"]
        self.assertEqual(
            [stage["quartic_naive_height_bound"] for stage in stages],
            [50_000, 250_000, 1_000_000],
        )
        self.assertEqual([stage["population_searched"] for stage in stages], [24, 10, 3])
        self.assertTrue(all(stage["completed"] == stage["population_searched"] for stage in stages))
        self.assertTrue(all(stage["timeouts"] == 0 and stage["errors"] == 0 for stage in stages))
        h50 = {record["constructor_parameter_T"]: record for record in stages[0]["ranked_population"]}
        for parameter in ("137/1022", "2941/1402", "875/1046", "241/1058", "823/1210"):
            self.assertIn(parameter, h50)
            self.assertEqual(h50[parameter]["height_rank"]["stable_numerical_rank"], 12)
            self.assertEqual(
                h50[parameter]["point_search"]["new_distinct_jacobian_sign_pairs_beyond_21_predeclared"],
                0,
            )
        final = stages[-1]["ranked_population"]
        self.assertEqual(final[0]["constructor_parameter_T"], "4409/1052")
        self.assertEqual(final[0]["height_rank"]["stable_numerical_rank"], 13)
        self.assertEqual(data["maximum_stable_numerical_rank"], 13)
        self.assertEqual(data["exact_checkpoints_numerical_rank_at_least_21"], [])
        self.assertEqual(data["target"]["certified_hits"], [])


if __name__ == "__main__":
    unittest.main()
