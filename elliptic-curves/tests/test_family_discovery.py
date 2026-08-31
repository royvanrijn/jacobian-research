from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
PROGRAM_ROOT = ROOT / "elliptic-curves"
CAS_ROOT = PROGRAM_ROOT / "cas"
SCRIPTS_ROOT = PROGRAM_ROOT / "scripts"
for path in (PROGRAM_ROOT, CAS_ROOT, SCRIPTS_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from discover_record_families import (  # noqa: E402
    load_mestre_families,
    load_polynomial_families,
    load_targets,
)
from ecsearch.family_discovery import (  # noqa: E402
    PolynomialWeierstrassFamily,
    discover_target_families,
    fraction_valuation,
)


SPECIFICATION = (
    ROOT / "elliptic-curves/data/family-discovery/icarm_273_282_302.json"
)
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic-curves/"
    "icarm_273_282_302_family_discovery_v1.json"
)


class GeneratedConstructionSpaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.specification = json.loads(SPECIFICATION.read_text())

    def test_fermigier_family_is_emitted_by_parameter_generator(self) -> None:
        families, summary = load_mestre_families(
            self.specification["six_root_mestre"]
        )
        family = next(
            candidate
            for candidate in families
            if candidate.roots == (0, 29, 658, 722, 981, 1036)
        )
        parameter_pairs = {
            (
                origin["generator_parameters"]["u"],
                origin["generator_parameters"]["v"],
            )
            for origin in family.metadata["origins"]
            if origin.get("generator_kind") == "fermigier-six-root"
        }
        self.assertEqual(parameter_pairs, {("3", "5"), ("5", "3")})
        self.assertEqual(summary["generator_parameter_pairs_tested"], 25)
        self.assertEqual(summary["generator_degenerate_parameter_pairs"], 15)

    def test_modular_exclusion_retains_a_possible_root_at_infinity(self) -> None:
        family = PolynomialWeierstrassFamily(
            identifier="short-test",
            coefficient_polynomials={
                "a1": (0,),
                "a2": (0,),
                "a3": (0,),
                "a4": (0, 1),
                "a6": (1,),
            },
        )
        # For target j=3 the cubic leading coefficient vanishes modulo 5,
        # while the affine reduction is a nonzero constant.  This cannot be
        # used as a no-root witness because [1:0] remains possible.
        self.assertIsNone(family.modular_roots(Fraction(3), 5))

    def test_curve282_is_rediscovered_in_both_family_coordinates(self) -> None:
        target = load_targets(self.specification)[1]
        polynomial_family = load_polynomial_families(
            self.specification["polynomial_families"]
        )[0]
        mestre_families, _ = load_mestre_families(
            self.specification["six_root_mestre"]
        )
        mestre_family = next(
            family
            for family in mestre_families
            if family.roots == (0, 29, 658, 722, 981, 1036)
        )
        result = discover_target_families(
            target,
            (polynomial_family, mestre_family),
            modular_primes=(101, 103, 107),
        )
        matches = {
            (match["family_id"], match["parameter_name"]): match
            for match in result["q_isomorphism_matches"]
        }
        self.assertEqual(
            matches[("fermigier-mestre-rank12", "u")]["parameter"],
            "11671/42",
        )
        mestre_key = ("six-root-mestre:0,29,658,722,981,1036", "T")
        self.assertEqual(matches[mestre_key]["parameter"], "11671/21")
        self.assertEqual(
            matches[("fermigier-mestre-rank12", "u")][
                "q_isomorphism_invariant_scale"
            ],
            "882",
        )
        self.assertEqual(matches[mestre_key]["q_isomorphism_invariant_scale"], "147")


class PinnedDiscoveryArtifactTests(unittest.TestCase):
    def test_pinned_scope_and_results(self) -> None:
        artifact = json.loads(ARTIFACT.read_text())
        self.assertEqual(
            artifact["schema"], "elliptic-curves.generated-family-discovery.v1"
        )
        self.assertEqual(artifact["construction_space"]["total_family_count"], 2334)
        self.assertEqual(
            artifact["construction_space"]["six_root_generation_summary"],
            {
                "duplicate_family_emission_count": 6,
                "emitted_family_count_before_deduplication": 2339,
                "generator_degenerate_parameter_pairs": 15,
                "generator_parameter_pairs_tested": 25,
            },
        )
        by_target = {record["target"]: record for record in artifact["targets"]}
        self.assertEqual(
            {
                target: record["exact_factorization_survivor_count"]
                for target, record in by_target.items()
            },
            {
                "ICARM curve 273": 113,
                "ICARM curve 282": 114,
                "ICARM curve 302": 146,
            },
        )
        self.assertFalse(by_target["ICARM curve 273"]["q_isomorphism_matches"])
        self.assertFalse(by_target["ICARM curve 302"]["q_isomorphism_matches"])
        matches = by_target["ICARM curve 282"]["q_isomorphism_matches"]
        self.assertEqual(len(matches), 2)
        for match in matches:
            scale = Fraction(match["q_isomorphism_invariant_scale"])
            for prime_text, target_value in match[
                "target_repeated_prime_valuations"
            ].items():
                prime = int(prime_text)
                source_value = match[
                    "source_discriminant_valuations_at_target_repeated_primes"
                ][prime_text]
                self.assertEqual(
                    target_value,
                    source_value + 12 * fraction_valuation(scale, prime),
                )


if __name__ == "__main__":
    unittest.main()
