"""Check real retained proof routes against the catalogue's promised dependencies."""

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from dependency_graph import DependencyGraph
from dependency_sources import missing_source_closure


class RetainedDependencyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT / "MATH_STATUS.json").read_text())
        cls.graph = DependencyGraph(cls.data, ROOT)

    def test_every_principal_scope_has_a_route_and_local_proof_closure(self):
        for topic in self.data["topics"]:
            self.assertTrue(topic["scope_notes"], topic["id"])
            for result in topic["results"]:
                self.assertTrue(self.graph.incoming[result["id"]], result["id"])
                self.assertIn(result["dependency_coverage"], {"partial", "complete"})
        self.assertEqual(missing_source_closure(self.graph), [])

    def test_included_arithmetic_proof_changes_reach_the_stated_consequence(self):
        path = "../unpublished/common-arithmetic-fibers/sections/01-keller-fibers.tex"
        affected = self.graph.source_impact(path)
        self.assertIn("ARITHMETIC-COMMON-ARITHMETIC-FIBERS", affected)
        self.assertIn("ARITHMETIC-JC-NEGATION", affected)
        changed = deepcopy(self.data)
        for node in changed["graph"]["nodes"]:
            for source in node.get("sources", []):
                if source["path"] == path:
                    source["sha256"] = "0" * 64
        graph = DependencyGraph(changed, ROOT)
        self.assertEqual(graph.evaluate()["ARITHMETIC-COMMON-ARITHMETIC-FIBERS"]["state"], "stale")
        # An independent cubic collision may continue to support the same negated
        # conjecture, so its conclusion must not be unconditionally downgraded here.
        self.assertTrue(any(route["state"] == "stale"
                            for route in graph.evaluate()["ARITHMETIC-JC-NEGATION"]["routes"]))

    def test_known_tang_and_trace_inputs_cannot_be_bypassed(self):
        cases = [
            ("HESS-TANG", "HESSIAN-RICCI-FLAT"),
            ("HESS-LW", "HESSIAN-MOMENT-TRACE"),
            ("HESS-LY", "HESSIAN-MOMENT-TRACE"),
            ("HESS-KN", "HESSIAN-MOMENT-TRACE"),
        ]
        for premise, conclusion in cases:
            with self.subTest(premise=premise):
                scenario = self.graph.simulate({premise: "stale"})
                self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_wirtinger_verification_reaches_named_equalities_without_blocking_lift_bounds(self):
        equality_claims = (
            "MRC-ALL-CHARACTERISTIC-LIFTS-THM-PANEL",
            "MRC-FINAL-FIVE-DOUBLE-COVERS-SHARP-FOUR-THM-FIVE",
            "MRC-FINAL-FIVE-DOUBLE-COVERS-SHARP-FOUR-COROLLARY-2",
            "MRC-FOUR-PIECE-DOUBLE-COVERS-SHARP-FOUR-THM-KNOTS",
            "MRC-THREE-PIECE-DOUBLE-COVERS-SHARP-FOUR-THEOREM-1",
            "MRC-MUTATION-COMMON-DOUBLE-COVERS-THM-BLOCKS",
            "MRC-ORIGINAL-SIXTEEN-KNOWN-GRAPH-COVER-FAMILY-THM-FAMILY",
            "MRC-GENERAL-COVER-THM-TWENTYFOUR",
        )
        # The published upper-bound theorem is now cached and its coloring
        # definition checked. Revoking that input must still block equalities,
        # while leaving purely algebraic lower bounds available.
        actual = self.graph.evaluate()
        self.assertEqual(actual["EXT-WIRTINGER-BRIDGE"]["state"], "proved")
        scenario = self.graph.simulate({"EXT-WIRTINGER-BRIDGE": "stale"})
        for conclusion in equality_claims:
            with self.subTest(conclusion=conclusion):
                self.assertEqual(actual[conclusion]["state"], "proved")
                self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for lower_bound in ("MRC-ALL-CHARACTERISTIC-LIFTS-THM-FIELD",
                            "MRC-GENERAL-COVER-THM-GENERAL",
                            "MRC-RESTRICTED-CAPACITY-THM-SEVEN"):
            with self.subTest(lower_bound=lower_bound):
                self.assertEqual(actual[lower_bound]["state"], "proved")
                self.assertEqual(scenario["nodes"][lower_bound]["state"], "proved")

    def test_seven_sheet_data_cannot_be_bypassed_or_invalidate_uniform_bound(self):
        changed = deepcopy(self.data)
        source = next(n for n in changed["graph"]["nodes"]
                      if n["id"] == "SRC-MRC-SEVEN-SHEET-CERTIFICATES")
        for pin in source["sources"]:
            if pin["path"].endswith("seven_sheet_cover_inputs.json"):
                pin["sha256"] = "0" * 64
        graph = DependencyGraph(changed, ROOT)
        states = graph.evaluate()
        self.assertEqual(states["MRC-SEVEN-SHEET-EXAMPLES"]["state"], "stale")
        self.assertEqual(states["MRC-GENERAL-COVER-THM-GENERAL"]["state"], "proved")
        self.assertEqual(states["MRC-RESTRICTED-CAPACITY-THM-SEVEN"]["state"], "proved")

    def test_twisted_family_import_does_not_underwrite_the_new_equalities(self):
        scenario = self.graph.simulate({"EXT-BBK-TWISTED-BRIDGE": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-TWISTED-FAMILY-SCOPE"]["state"], "stale")
        self.assertEqual(scenario["nodes"]["MRC-SEVEN-SHEET-EXAMPLES"]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-FISHNET-AND-PLAT-SCOPE"]["state"], "proved")

    def test_generalized_sum_rank_import_is_separate_from_finite_exclusion(self):
        scenario = self.graph.simulate({"EXT-GRUSHKO": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-NORMALIZED-COVER-ADDITIVITY"]["state"], "stale")
        self.assertEqual(scenario["nodes"]["MRC-NO-SHARP-FINITE-COXETER"]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-SEVEN-SHEET-EXAMPLES"]["state"], "proved")

    def test_scoped_family_exclusions_do_not_close_global_priority_or_mrc(self):
        actual = self.graph.evaluate()
        for result in ("MRC-FISHNET-AND-PLAT-SCOPE", "MRC-TWISTED-FAMILY-SCOPE",
                       "MRC-NO-SHARP-FINITE-COXETER", "MRC-NORMALIZED-COVER-ADDITIVITY"):
            self.assertEqual(actual[result]["state"], "proved")
        self.assertEqual(actual["GOAL-MRC"]["state"], "open")
        self.assertEqual(actual["GAP-MERIDIONAL-RANK-3"]["state"], "open")

    def test_canonical_coxeter_import_does_not_underwrite_field_obstruction_or_equality(self):
        actual = self.graph.evaluate()
        self.assertEqual(actual["MRC-NO-SHARP-COXETER"]["state"], "proved")
        self.assertEqual(actual["GAP-MRC-RATIONAL-REPLACEMENTS"]["state"], "proved")
        scenario = self.graph.simulate({"EXT-COXETER-CANONICAL": "stale"})
        for conclusion in ("MRC-NO-SHARP-COXETER", "MRC-UNIVERSAL-REFLECTION-COR-BBK"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-UNIVERSAL-REFLECTION-THM-REFLECTION",
                           "MRC-SEVEN-SHEET-EXAMPLES", "MRC-GENERAL-COVER-THM-GENERAL"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        self.assertEqual(actual["GOAL-MRC"]["state"], "open")
        self.assertEqual(actual["GAP-MERIDIONAL-RANK-3"]["state"], "open")
        # The independent four-face topology argument still closes the narrow
        # replacement comparison when the Coxeter import is revoked.
        self.assertEqual(scenario["nodes"]["GAP-MRC-RATIONAL-REPLACEMENTS"]["state"], "proved")

    def test_optional_reflection_discovery_certificates_are_not_a_field_scope_premise(self):
        scenario = self.graph.simulate({"SRC-MRC-REFLECTION-ALGEBRA": "stale"})
        # The retained identities15 are useful independent controls; the direct
        # coefficient proof and exact diagram identity establish the theorem.
        self.assertEqual(scenario["nodes"]["MRC-NO-SHARP-COXETER"]["state"], "proved")
        scenario = self.graph.simulate({"SRC-MRC-REFLECTION-DIAGRAM": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-NO-SHARP-COXETER"]["state"], "stale")
        self.assertEqual(scenario["nodes"]["MRC-UNIVERSAL-REFLECTION-LEM-SPAN"]["state"], "proved")

    def test_classical_filling_import_blocks_topology_but_not_coxeter_or_equalities(self):
        scenario = self.graph.simulate({"EXT-MARTELLI-M6-FILLINGS": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-RATIONAL-REPLACEMENT-TOPOLOGY"]["state"], "stale")
        self.assertEqual(scenario["nodes"]["MRC-NO-SHARP-COXETER"]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-SEVEN-SHEET-EXAMPLES"]["state"], "proved")
        self.assertEqual(scenario["nodes"]["GAP-MRC-RATIONAL-REPLACEMENTS"]["state"], "proved")
        both = self.graph.simulate({"EXT-MARTELLI-M6-FILLINGS": "stale",
                                    "EXT-COXETER-CANONICAL": "stale"})
        self.assertEqual(both["nodes"]["GAP-MRC-RATIONAL-REPLACEMENTS"]["state"], "stale")

    def test_exterior_hyperbolicity_is_separate_from_cover_topology_and_reflections(self):
        scenario = self.graph.simulate({"EXT-MENASCO-ALTERNATING-HYPERBOLIC": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-FOUR-ALTERNATING-HYPERBOLIC"]["state"], "stale")
        for conclusion in ("MRC-SEED-COVER-TOPOLOGY", "MRC-NO-SHARP-COXETER",
                           "MRC-RATIONAL-REPLACEMENT-TOPOLOGY", "MRC-SEVEN-SHEET-EXAMPLES"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_finite_affine_route_survives_cover_group_revocation(self):
        scenario = self.graph.simulate({"SRC-MRC-SEVEN-SHEET-CERTIFICATES": "stale"})
        for conclusion in ("MRC-TERNARY-AFFINE-QUOTIENTS",
                           "MRC-GENERAL-COVER-THM-TWENTYFOUR",
                           "MRC-RANK3-THM-RANK",
                           "MRC-REFLECTIVE-AUGMENTATION-THM-UPPER"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        # The older principal also asserts58 specific cover groups, so that
        # broader data claim cannot inherit just the alternative equalities.
        self.assertEqual(scenario["nodes"]["MRC-SEVEN-SHEET-EXAMPLES"]["state"], "stale")

    def test_affine_controls_are_separate_from_uniform_realization_and_old_covers(self):
        scenario = self.graph.simulate({"SRC-MRC-TERNARY-AFFINE-CONTROLS": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-TERNARY-AFFINE-QUOTIENTS"]["state"], "stale")
        self.assertEqual(scenario["nodes"]["MRC-FINITE-AFFINE-COVER-REALIZATION"]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-SEVEN-SHEET-EXAMPLES"]["state"], "proved")
        scenario = self.graph.simulate({"SRC-MRC-TERNARY-AFFINE-DIAGRAMS": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-TERNARY-AFFINE-THM-KNOTS"]["state"], "stale")
        self.assertEqual(scenario["nodes"]["MRC-TERNARY-AFFINE-THM-GROUP"]["state"], "proved")

    def test_ordinary_rank_import_is_not_a_meridional_or_augmentation_premise(self):
        scenario = self.graph.simulate({"EXT-SNAPPY332-CELL-AND-PD": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-RANK3-LEM-UPPER"]["state"], "stale")
        for conclusion in ("MRC-TERNARY-AFFINE-QUOTIENTS",
                           "MRC-REFLECTIVE-AUGMENTATION-GAP",
                           "MRC-NO-SHARP-COXETER"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_deck_character_argument_does_not_import_cord_identification(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-REFLECTIVE-CORD-PRESENTATION": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-REFLECTIVE-AUGMENTATION-THM-UPPER"]["state"], "stale")
        self.assertEqual(scenario["nodes"]["MRC-REFLECTIVE-AUGMENTATION-THM-DECK"]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-TERNARY-AFFINE-QUOTIENTS"]["state"], "proved")

    def test_real_multistep_progress_stays_hypothetical_and_preserves_assurance(self):
        before = deepcopy(self.data)
        # Fine decomposition exposes a retained-input verification gap. A
        # hypothetical proof must carry that assumption too, rather than bypass it.
        overrides = {"GAP-AUTOMATIC-T3-PERFECTION": "proved"}
        for blocker in self.graph.blockers("HESS-LZ"):
            if blocker["type"] == "node":
                overrides[blocker["id"]] = "proved"
        result = self.graph.simulate(overrides)
        for conclusion in ("GOAL-JC2", "GOAL-HC4"):
            self.assertEqual(result["nodes"][conclusion]["state"], "conditional")
            self.assertIn("GAP-AUTOMATIC-T3-PERFECTION",
                          result["nodes"][conclusion]["assumptions"])
        self.assertEqual(self.data, before)
        self.assertFalse(self.graph.evaluate()["GOAL-HC4"]["hypothetical"])

    def test_reverse_hessian_equivalence_propagates_without_extra_tang_context(self):
        result = self.graph.simulate({"GOAL-HC4": "proved"})
        self.assertEqual(result["nodes"]["GOAL-JC2"]["state"], "conditional")
        self.assertEqual(result["nodes"]["GOAL-JC2"]["assumptions"], ["GOAL-HC4"])


if __name__ == "__main__":
    unittest.main()
