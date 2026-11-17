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

    def test_a6_equalities_do_not_import_the_cover_family_exclusion(self):
        scenario = self.graph.simulate({"SRC-MRC-A6-TOPOLOGY-CERTIFICATES": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-A6-COVER-TOPOLOGY"]["state"], "stale")
        for conclusion in ("MRC-A6-AFFINE-EXAMPLES", "MRC-SIXTY-EIGHT-NAMED-EQUALITIES",
                           "MRC-A6-AFFINE-THM-GROUPS", "MRC-UNIVERSAL-REFLECTION-THM-A6EIGHT"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_a6_class_rank_is_separate_from_geometric_bridge_import(self):
        scenario = self.graph.simulate({"EXT-WIRTINGER-BRIDGE": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-A6-AFFINE-THM-SHARP"]["state"], "stale")
        for conclusion in ("MRC-A6-AFFINE-LEM-QUOTIENT", "MRC-A6-AFFINE-THM-GROUPS",
                           "MRC-TERNARY-AFFINE-THM-ALLPRIMES"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_uniform_integral_control_is_not_a_named_knot_premise(self):
        scenario = self.graph.simulate({"SRC-MRC-UNIFORM-AFFINE-CERTIFICATES": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-TERNARY-AFFINE-THM-ALLPRIMES"]["state"], "stale")
        for conclusion in ("MRC-TERNARY-AFFINE-THM-KNOTS", "MRC-A6-AFFINE-EXAMPLES",
                           "MRC-SIXTY-EIGHT-NAMED-EQUALITIES"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_kch_correspondence_only_underwrites_the_augmentation_consequence(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-KCH-AUGMENTATION-CORRESPONDENCE": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-INVOLUTORY-KCH-COR-AUGMENTATION"]["state"], "stale")
        for conclusion in ("MRC-INVOLUTORY-KCH-THM-SIXTEEN", "MRC-NO-SHARP-COXETER",
                           "MRC-A7-BINARY-EXAMPLES"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_generic_kch_correspondence_only_supports_augmentation_outputs(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-KCH-AUGMENTATION-CORRESPONDENCE": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-GENERIC-KCH-THM-GENERIC"]["state"], "proved")
        for conclusion in ("MRC-GENERIC-KCH-COR-AUGMENTATION", "MRC-GENERIC-KCH-COR-GAP",
                           "MRC-GENERIC-KCH-COR-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        self.assertEqual(scenario["nodes"]["MRC-INVOLUTION-SUMS-SHARP-FAMILIES"]["state"], "proved")

    def test_generic_augmentation_upper_does_not_require_positive_attainment(self):
        scenario = self.graph.simulate({"MRC-REFLECTIVE-AUGMENTATION-THM-POSITIVE": "stale"})
        for conclusion in ("MRC-GENERIC-KCH-THM-GENERIC", "MRC-GENERIC-KCH-COR-AUGMENTATION"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        for conclusion in ("MRC-GENERIC-KCH-COR-GAP", "MRC-GENERIC-KCH-COR-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_augmentation_sum_import_is_separate_from_seed_and_cover_sums(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-AUGMENTATION-SUM": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-GENERIC-KCH-COR-SUMS"]["state"], "stale")
        for conclusion in ("MRC-GENERIC-KCH-COR-GAP", "MRC-GENERIC-KCH-COR-AUGMENTATION",
                           "MRC-INVOLUTION-SUMS-SHARP-FAMILIES"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_generic_sixteen_upper_is_independent_of_new_a6_attainment(self):
        scenario = self.graph.simulate({"SRC-MRC-A6-RANK-THREE-CERTIFICATES": "stale"})
        for conclusion in ("MRC-GENERIC-KCH-EXT-THM-ALL16", "MRC-GENERIC-KCH-EXT-COR-SEMISIMPLE",
                           "MRC-FINITE-COMPLEX-COR-KNOTS", "MRC-GENERIC-KCH-EXT-COR-AUGMENTATION"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        for conclusion in ("MRC-A6-RANK-THREE-THM-ATTAINMENT", "MRC-GENERIC-KCH-EXT-COR-GAP",
                           "MRC-GENERIC-KCH-EXT-COR-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_generic_sixteen_group_bounds_do_not_need_augmentation_correspondence(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-KCH-AUGMENTATION-CORRESPONDENCE": "stale"})
        for conclusion in ("MRC-GENERIC-KCH-EXT-THM-ALL16", "MRC-GENERIC-KCH-EXT-COR-SEMISIMPLE",
                           "MRC-FINITE-COMPLEX-COR-KNOTS", "MRC-A6-RANK-THREE-COR-REFLECTIVE"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        for conclusion in ("MRC-GENERIC-KCH-EXT-COR-AUGMENTATION", "MRC-GENERIC-KCH-EXT-COR-GAP",
                           "MRC-A6-RANK-THREE-COR-AUGMENTATION"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_finite_reflection_class_reduction_does_not_import_classification(self):
        scenario = self.graph.simulate({"EXT-SHEPHARD-TODD-LOW-RANK-GENERATION": "stale"})
        for conclusion in ("MRC-FINITE-COMPLEX-LEM-CLASS", "MRC-GENERIC-KCH-EXT-THM-ALL16",
                           "MRC-GENERIC-KCH-EXT-COR-SEMISIMPLE"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-FINITE-COMPLEX-THM-LOWRANK"]["state"], "stale")
        # The new global capacity proof is a second route to the knot bound.
        self.assertEqual(scenario["nodes"]["MRC-FINITE-COMPLEX-COR-KNOTS"]["state"], "proved")

    def test_reduced_core_is_a_separate_premise_of_both_generic_proofs(self):
        scenario = self.graph.simulate({"MRC-GENERIC-KCH-CORE-373": "stale"})
        for conclusion in ("MRC-GENERIC-KCH-THM-GENERIC", "MRC-GENERIC-KCH-EXT-LEM-CORES",
                           "MRC-GENERIC-KCH-EXT-THM-ALL16"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-A6-RANK-THREE-THM-ATTAINMENT", "MRC-FINITE-COMPLEX-THM-LOWRANK"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_same_parameter_trefoil_input_only_blocks_augmentation_sum_families(self):
        scenario = self.graph.simulate({"MRC-TREFOIL-MINUS-ONE-AUGMENTATION": "stale"})
        for conclusion in ("MRC-GENERIC-KCH-COR-SUMS", "MRC-GENERIC-KCH-EXT-COR-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-GENERIC-KCH-EXT-COR-GAP", "MRC-INVOLUTION-SUMS-SHARP-FAMILIES"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_finite_complex_knot_application_is_separate_from_class_theorem(self):
        scenario = self.graph.simulate({"MRC-GENERIC-KCH-EXT-COR-SEMISIMPLE": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-FINITE-COMPLEX-COR-KNOTS"]["state"], "proved")
        both = self.graph.simulate({"MRC-GENERIC-KCH-EXT-COR-SEMISIMPLE": "stale",
                                    "MRC-COMPLEX-DOMINATION-LEM-CAPACITY": "stale"})
        self.assertEqual(both["nodes"]["MRC-FINITE-COMPLEX-COR-KNOTS"]["state"], "stale")
        for conclusion in ("MRC-FINITE-COMPLEX-THM-LOWRANK", "MRC-GENERIC-KCH-EXT-THM-ALL16",
                           "MRC-A6-RANK-THREE-COR-AUGMENTATION"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_paired_span_and_form_bounds_do_not_require_irreducibility_consequence(self):
        scenario = self.graph.simulate({"MRC-GENERIC-KCH-EXT-THM-ALL16": "stale"})
        for conclusion in ("MRC-GENERIC-KCH-EXT-THM-PAIREDSPANS", "MRC-GENERIC-KCH-EXT-COR-FORMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        for conclusion in ("MRC-GENERIC-KCH-EXT-COR-SEMISIMPLE", "MRC-GENERIC-KCH-EXT-COR-AUGMENTATION"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_paired_span_revocation_does_not_revoke_positive_augmentation_witnesses(self):
        scenario = self.graph.simulate({"MRC-GENERIC-KCH-EXT-THM-PAIREDSPANS": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-GENERIC-KCH-EXT-COR-FORMS"]["state"], "stale")
        self.assertEqual(scenario["nodes"]["MRC-A6-RANK-THREE-THM-ATTAINMENT"]["state"], "proved")

    def test_simple_a8_route_preserves_numeric_twentyfour_without_old_affine_or_cover_data(self):
        scenario = self.graph.simulate({"SRC-MRC-TERNARY-AFFINE-CONTROLS": "stale",
                                       "SRC-MRC-SEVEN-SHEET-CERTIFICATES": "stale"})
        for conclusion in ("MRC-A8-THM-SHARP", "MRC-GENERAL-COVER-THM-TWENTYFOUR"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-SEVEN-SHEET-EXAMPLES"]["state"], "stale")

    def test_simple_a8_group_bound_is_separate_from_bridge_import(self):
        scenario = self.graph.simulate({"EXT-WIRTINGER-BRIDGE": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-A8-THM-SHARP"]["state"], "stale")
        self.assertEqual(scenario["nodes"]["MRC-A8-THM-GROUP"]["state"], "proved")
        scenario = self.graph.simulate({"SRC-MRC-A8-CERTIFICATES": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-A8-LEM-TRANSITIVE"]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-GENERAL-COVER-THM-TWENTYFOUR"]["state"], "proved")

    def test_g31_explicit_group_proof_is_separate_from_its_published_name(self):
        scenario = self.graph.simulate({"EXT-KONISHI-MINABE-G31-MODEL": "stale"})
        for conclusion in ("MRC-G31-THM-QUOTIENT", "MRC-G31-COR-CLASSRANK", "MRC-G31-COR-DOMINANCE"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-G31-ROUTE-DOMINATED"]["state"], "stale")

    def test_g31_augmentation_import_is_separate_from_group_class_bound(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-KCH-AUGMENTATION-CORRESPONDENCE": "stale"})
        for conclusion in ("MRC-G31-THM-QUOTIENT", "MRC-G31-COR-CLASSRANK", "MRC-A8-THM-GROUP"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-G31-COR-DOMINANCE"]["state"], "stale")

    def test_new17_analytic_import_is_not_a_finite_or_crossing_premise(self):
        scenario = self.graph.simulate({"EXT-NEW17-KRAWCZYK-INTERVAL": "stale"})
        for conclusion in ("MRC-NEW17-THM-HYPERBOLIC", "MRC-NEW17-CLOSED-REPRESENTATION-VOLUME",
                           "MRC-NEW17-THM-NONGRAPH"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-NEW17-THM-SHARP", "MRC-NEW17-ALEXANDER", "MRC-NEW17-THM-CROSSING",
                           "MRC-NEW17-COR-TREFOILS", "MRC-NEW17-REFLECTION-COR-COXETER"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_new17_actual_cover_and_cusped_hyperbolicity_are_separate(self):
        scenario = self.graph.simulate({"EXT-NEW17-COVER-PERIPHERAL": "stale"})
        for conclusion in ("MRC-NEW17-ACTUAL-DOUBLE-COVER", "MRC-NEW17-THM-NONGRAPH"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-NEW17-THM-HYPERBOLIC", "MRC-NEW17-THM-SHARP",
                           "MRC-NEW17-COR-TREFOILS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_new17_graph_volume_vanishing_has_an_independent_all_parameter_route(self):
        scenario = self.graph.simulate({"EXT-LTY-GRAPH-REPRESENTATION-VOLUME": "stale"})
        for conclusion in ("MRC-NEW17-THM-NONGRAPH", "MRC-ALL-COVERS-THM-ALLCOVERS",
                           "MRC-NEW17-CLOSED-REPRESENTATION-VOLUME", "MRC-NEW17-THM-HYPERBOLIC",
                           "MRC-NEW17-THM-SHARP", "MRC-NEW17-REFLECTION-COR-COXETER"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        both = self.graph.simulate({"EXT-LTY-GRAPH-REPRESENTATION-VOLUME": "stale",
                                    "MRC-ALL-COVERS-ALL-NONGRAPH": "stale"})
        self.assertEqual(both["nodes"]["MRC-NEW17-THM-NONGRAPH"]["state"], "stale")

    def test_new17_yamada_only_supplies_the_crossing_minimality_route(self):
        scenario = self.graph.simulate({"EXT-YAMADA-SEIFERT-BRAID": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-NEW17-THM-CROSSING"]["state"], "stale")
        for conclusion in ("MRC-NEW17-ALEXANDER", "MRC-NEW17-THM-SHARP", "MRC-NEW17-COR-TREFOILS",
                           "MRC-NEW17-THM-NONGRAPH"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_new17_universal_certificate_does_not_underwrite_numeric_or_topology_proofs(self):
        scenario = self.graph.simulate({"SRC-MRC-NEW17-REFLECTION-CERTIFICATES": "stale"})
        for conclusion in ("MRC-NEW17-REFLECTION-COR-COXETER", "MRC-NEW17-REFLECTION-COR-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-NEW17-THM-SHARP", "MRC-NEW17-THM-CROSSING", "MRC-NEW17-THM-HYPERBOLIC",
                           "MRC-NEW17-THM-NONGRAPH", "MRC-NEW17-COR-TREFOILS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_new17_sum_cover_bound_is_separate_from_coxeter_sum_obstruction(self):
        scenario = self.graph.simulate({"MRC-GENERAL-COVER-THM-GENERAL": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-NEW17-COR-TREFOILS"]["state"], "stale")
        for conclusion in ("MRC-NEW17-REFLECTION-COR-SUMS", "MRC-NEW17-THM-SHARP",
                           "MRC-NEW17-THM-NONGRAPH"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_new17_symmetric_augmentation_comparison_does_not_underwrite_coxeter_or_mr(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-MAXIMAL-BRAID": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-NEW17-FAMILY-COR-SYMMETRIC"]["state"], "stale")
        for conclusion in ("MRC-NEW17-REFLECTION-COR-COXETER", "MRC-NEW17-THM-SHARP",
                           "MRC-NEW17-FAMILY-PROP-RASMUSSEN"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_new17_abe_formula_is_separate_from_meridional_and_coxeter_bounds(self):
        scenario = self.graph.simulate({"EXT-ABE-HOMOGENEOUS-RASMUSSEN": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-NEW17-FAMILY-PROP-RASMUSSEN"]["state"], "stale")
        for conclusion in ("MRC-NEW17-THM-SHARP", "MRC-NEW17-REFLECTION-COR-COXETER",
                           "MRC-NEW17-FAMILY-PROP-COXETERFAMILIES"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_new17_comparison_imports_both_seed_and_sum_coxeter_bounds(self):
        scenario = self.graph.simulate({"MRC-NEW17-REFLECTION-COR-SUMS": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-NEW17-COR-COMPARISON"]["state"], "stale")
        for conclusion in ("MRC-NEW17-REFLECTION-COR-COXETER", "MRC-NEW17-COR-TREFOILS",
                           "MRC-NEW17-THM-SHARP"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_new17_generic_upper_uses_only_the_independent_392_core(self):
        scenario = self.graph.simulate({"MRC-GENERIC-KCH-CORE-373": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-GENERIC-KCH-EXT-LEM-CORES"]["state"], "stale")
        for conclusion in ("MRC-GENERIC-KCH-CORE-392", "MRC-NEW17-GENERIC-THM-PAIREDSPANS",
                           "MRC-NEW17-GENERIC-COR-AUGMENTATION", "MRC-NEW17-GENERIC-EXACT-AUGMENTATION"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_new17_392_revocation_preserves_independent_positive_and_symmetric_proofs(self):
        scenario = self.graph.simulate({"MRC-GENERIC-KCH-CORE-392": "stale"})
        for conclusion in ("MRC-NEW17-GENERIC-THM-PAIREDSPANS", "MRC-NEW17-GENERIC-COR-AUGMENTATION"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-NEW17-GENERIC-THM-ATTAINMENT", "MRC-NEW17-GENERIC-POSITIVE-CLASS",
                           "MRC-NEW17-REFLECTION-COR-COXETER", "MRC-NEW17-THM-SHARP"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_new17_augmentation_import_is_not_a_representation_or_finite_class_premise(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-KCH-AUGMENTATION-CORRESPONDENCE": "stale"})
        for conclusion in ("MRC-NEW17-GENERIC-THM-PAIREDSPANS", "MRC-NEW17-GENERIC-COR-FORMS",
                           "MRC-NEW17-GENERIC-COR-FINITECLASS", "MRC-NEW17-GENERIC-THM-ATTAINMENT",
                           "MRC-NEW17-GENERIC-COR-FINITEATTAINMENT", "MRC-NEW17-COR-TREFOILS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        for conclusion in ("MRC-NEW17-GENERIC-COR-AUGMENTATION", "MRC-NEW17-GENERIC-EXACT-AUGMENTATION",
                           "MRC-NEW17-GENERIC-COR-GAP", "MRC-NEW17-GENERIC-AUGMENTATION-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_new17_generic_upper_does_not_require_positive_attainment(self):
        scenario = self.graph.simulate({"MRC-NEW17-GENERIC-THM-ATTAINMENT": "stale"})
        for conclusion in ("MRC-NEW17-GENERIC-THM-PAIREDSPANS", "MRC-NEW17-GENERIC-COR-AUGMENTATION",
                           "MRC-NEW17-GENERIC-COR-FINITECLASS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        for conclusion in ("MRC-NEW17-GENERIC-EXACT-AUGMENTATION", "MRC-NEW17-GENERIC-COR-GAP",
                           "MRC-NEW17-GENERIC-AUGMENTATION-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_new17_augmentation_sum_is_separate_from_meridional_cover_sums(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-AUGMENTATION-SUM": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-NEW17-GENERIC-AUGMENTATION-SUMS"]["state"], "stale")
        for conclusion in ("MRC-NEW17-GENERIC-EXACT-AUGMENTATION", "MRC-NEW17-COR-TREFOILS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        scenario = self.graph.simulate({"MRC-NEW17-THM-SHARP": "stale"})
        for conclusion in ("MRC-NEW17-GENERIC-EXACT-AUGMENTATION", "MRC-NEW17-GENERIC-AUGMENTATION-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-NEW17-GENERIC-COR-SUMS"]["state"], "stale")

    def test_new17_finite_class_attainment_is_separate_from_the_classification_upper(self):
        scenario = self.graph.simulate({"MRC-FINITE-COMPLEX-THM-LOWRANK": "stale"})
        for conclusion in ("MRC-NEW17-GENERIC-COR-FINITECLASS", "MRC-NEW17-GENERIC-COR-FINITEATTAINMENT"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        scenario = self.graph.simulate({"MRC-FINITE-COMPLEX-THM-LOWRANK": "stale",
                                       "MRC-COMPLEX-DOMINATION-LEM-CAPACITY": "stale"})
        for conclusion in ("MRC-NEW17-GENERIC-COR-FINITECLASS", "MRC-NEW17-GENERIC-COR-FINITEATTAINMENT"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-NEW17-GENERIC-POSITIVE-CLASS", "MRC-NEW17-GENERIC-THM-ATTAINMENT",
                           "MRC-NEW17-GENERIC-EXACT-AUGMENTATION"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_global_complex_capacity_is_separate_from_augmentation_conversion(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-KCH-AUGMENTATION-CORRESPONDENCE": "stale"})
        for conclusion in ("MRC-COMPLEX-DOMINATION-LEM-CAPACITY", "MRC-COMPLEX-DOMINATION-COR-FINITECAPACITY",
                           "MRC-COMPLEX-DOMINATION-LEM-ADJOINT"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        for conclusion in ("MRC-COMPLEX-DOMINATION-THM-DOMINATION", "MRC-COMPLEX-DOMINATION-COR-KNOTS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_global_complex_theorem_has_no_specific_new17_premise(self):
        scenario = self.graph.simulate({"MRC-NEW17-GENERIC-EXACT-AUGMENTATION": "stale"})
        for conclusion in ("MRC-COMPLEX-DOMINATION-LEM-CAPACITY", "MRC-COMPLEX-DOMINATION-THM-DOMINATION",
                           "MRC-COMPLEX-DOMINATION-COR-FINITECAPACITY"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-COMPLEX-DOMINATION-COR-KNOTS"]["state"], "stale")

    def test_global_complex_classification_import_does_not_underwrite_adjoint_algebra(self):
        scenario = self.graph.simulate({"EXT-ST54-GLOBAL-GENERATION": "stale"})
        for conclusion in ("MRC-COMPLEX-DOMINATION-LEM-ACTIVE", "MRC-COMPLEX-DOMINATION-LEM-DETERMINANT",
                           "MRC-COMPLEX-DOMINATION-LEM-ADJOINT", "MRC-NEW17-GENERIC-POSITIVE-CLASS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-COMPLEX-DOMINATION-LEM-CAPACITY"]["state"], "stale")

    def test_global_complex_capacity_needs_the_exceptional_group_quotient(self):
        scenario = self.graph.simulate({"MRC-G31-THM-QUOTIENT": "stale"})
        for conclusion in ("MRC-COMPLEX-DOMINATION-LEM-CAPACITY", "MRC-COMPLEX-DOMINATION-THM-DOMINATION"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-COMPLEX-DOMINATION-LEM-ADJOINT", "MRC-NEW17-GENERIC-COR-FINITECLASS",
                           "MRC-NEW17-THM-SHARP"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_global_capacity_and_lowrank_routes_fail_independently(self):
        scenario = self.graph.simulate({"MRC-FINITE-COMPLEX-THM-LOWRANK": "stale"})
        for conclusion in ("MRC-FINITE-COMPLEX-COR-KNOTS", "MRC-NEW17-GENERIC-COR-FINITECLASS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        scenario = self.graph.simulate({"MRC-COMPLEX-DOMINATION-LEM-CAPACITY": "stale"})
        for conclusion in ("MRC-FINITE-COMPLEX-COR-KNOTS", "MRC-NEW17-GENERIC-COR-FINITECLASS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        scenario = self.graph.simulate({"MRC-COMPLEX-DOMINATION-LEM-CAPACITY": "stale",
                                       "MRC-FINITE-COMPLEX-THM-LOWRANK": "stale"})
        for conclusion in ("MRC-FINITE-COMPLEX-COR-KNOTS", "MRC-NEW17-GENERIC-COR-FINITECLASS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_rank_one_quotient_reduction_is_independent_of_the_finite_panel(self):
        scenario = self.graph.simulate({"SRC-MRC-RANKONE-BINARY-CERTIFICATES": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-RANKONE-BINARY-LEM-REDUCTION"]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-RANKONE-BINARY-COR-BOUNDARY"]["state"], "stale")
        self.assertEqual(scenario["nodes"]["MRC-A7-BINARY-EXAMPLES"]["state"], "proved")

    def test_projective_numeric_bound_is_separate_from_known_family_membership(self):
        scenario = self.graph.simulate({"MRC-GENERAL-COVER-THM-GENERAL": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-SL3-PROJECTIVE-THM-EQUALITIES"]["state"], "stale")
        for conclusion in ("MRC-SL3-PROJECTIVE-THM-GRAPH", "MRC-SL3-PROJECTIVE-COR-KNOWN",
                           "MRC-A7-BINARY-EXAMPLES", "MRC-INVOLUTORY-KCH-THM-SIXTEEN"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_finite_lifts_are_separate_from_augmentation_conversion(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-KCH-AUGMENTATION-CORRESPONDENCE": "stale"})
        for conclusion in ("MRC-FINITE-LIFTS-LEM-UNIQUE", "MRC-FINITE-LIFTS-LEM-CLASS",
                           "MRC-FINITE-LIFTS-THM-LIFTS", "MRC-FINITE-LIFTS-FORTY-EIGHT-POSITIVE"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-FINITE-LIFTS-COR-SIXTEEN"]["state"], "stale")

    def test_finite_lift_lower_bounds_do_not_require_generic_upper(self):
        scenario = self.graph.simulate({"MRC-GENERIC-KCH-EXT-THM-ALL16": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-FINITE-LIFTS-FORTY-EIGHT-POSITIVE"]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-FINITE-LIFTS-COR-SIXTEEN"]["state"], "stale")

    def test_prime_family_visual_primeness_is_separate_from_sharp_group_bound(self):
        scenario = self.graph.simulate({"EXT-FLO-HOMOGENEOUS-VISUAL-PRIMENESS": "stale"})
        for conclusion in ("MRC-PRIME-GAP-PRIMENESS", "MRC-ALT-PRIME-PRIMENESS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-PRIME-GAP-SHARP", "MRC-PRIME-GAP-THM-PAIRED",
                           "MRC-PRIME-GAP-AUGMENTATION", "MRC-ALT-PRIME-SHARP",
                           "MRC-ALT-PRIME-COR-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_prime_family_crossing_minimality_is_separate_from_mr_and_primeness(self):
        scenario = self.graph.simulate({"EXT-YAMADA-SEIFERT-BRAID": "stale"})
        for conclusion in ("MRC-PRIME-GAP-CROSSING", "MRC-ALT-PRIME-CROSSING"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-PRIME-GAP-SHARP", "MRC-PRIME-GAP-PRIMENESS",
                           "MRC-PRIME-GAP-LEM-ALEXANDER", "MRC-ALT-PRIME-SHARP"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_uniform_prime_gap_uses_only_the_independent_392_core(self):
        scenario = self.graph.simulate({"MRC-GENERIC-KCH-CORE-373": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-PRIME-GAP-THM-PAIRED"]["state"], "proved")
        scenario = self.graph.simulate({"MRC-GENERIC-KCH-CORE-392": "stale"})
        for conclusion in ("MRC-PRIME-GAP-THM-PAIRED", "MRC-PRIME-GAP-AUGMENTATION",
                           "MRC-PRIME-GAP-COR-COXETER", "MRC-PRIME-GAP-FINITECLASS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-PRIME-GAP-SHARP", "MRC-PRIME-GAP-POSITIVE",
                           "MRC-PRIME-GAP-LEM-LIFT", "MRC-ALT-PRIME-THM-FAMILY"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_singular_gram_extension_is_an_essential_separate_coxeter_premise(self):
        scenario = self.graph.simulate({"MRC-PRIME-GAP-LEM-LIFT": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-PRIME-GAP-COR-COXETER"]["state"], "stale")
        for conclusion in ("MRC-PRIME-GAP-AUGMENTATION", "MRC-PRIME-GAP-FINITECLASS",
                           "MRC-PRIME-GAP-SHARP", "MRC-NEW17-REFLECTION-COR-COXETER"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_prime_family_finite_class_bound_has_a_pure_group_capacity_route(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-KCH-AUGMENTATION-CORRESPONDENCE": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-PRIME-GAP-AUGMENTATION"]["state"], "stale")
        for conclusion in ("MRC-PRIME-GAP-FINITECLASS", "MRC-PRIME-GAP-COR-COXETER",
                           "MRC-PRIME-GAP-POSITIVE", "MRC-ALT-PRIME-THM-FAMILY"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_uniform_prime_sum_numeric_equality_has_an_independent_simple_target_route(self):
        scenario = self.graph.simulate({"MRC-GENERAL-COVER-THM-GENERAL": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-PRIME-GAP-SUM-MERIDIONAL"]["state"], "stale")
        for conclusion in ("MRC-ALT-PRIME-UNIFORM-TWIST-SUMS", "MRC-PRIME-GAP-COR-SUMS",
                           "MRC-ALT-PRIME-COR-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_general_a8_seed_sum_method_does_not_require_prime_or_genus_imports(self):
        scenario = self.graph.simulate({"EXT-FLO-HOMOGENEOUS-VISUAL-PRIMENESS": "stale",
                                       "EXT-CROMWELL-HOMOGENEOUS-GENUS": "stale",
                                       "EXT-YAMADA-SEIFERT-BRAID": "stale"})
        for conclusion in ("MRC-ALT-PRIME-COR-SUMS", "MRC-ALT-PRIME-THM-TARGETS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        self.assertEqual(scenario["nodes"]["MRC-ALT-PRIME-THM-FAMILY"]["state"], "stale")

    def test_alternating_target_written_ordinary_argument_is_separate_from_literal_controls(self):
        scenario = self.graph.simulate({"SRC-MRC-ALT-PRIME-CERTIFICATES": "stale"})
        for conclusion in ("MRC-ALT-PRIME-LEM-ORDINARY", "MRC-ALT-PRIME-LEM-SEESAWS",
                           "MRC-PRIME-GAP-SHARP"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        for conclusion in ("MRC-ALT-PRIME-THM-TARGETS", "MRC-ALT-PRIME-LEM-CLOSURE"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_cofinite_topology_has_an_independent_all_parameter_route(self):
        scenario = self.graph.simulate({"EXT-BC-ONE-BOUNDARY-NONGRAPH-FILLING": "stale"})
        for conclusion in ("MRC-COFINITE-LEM-PERSISTENCE", "MRC-COFINITE-THM-ONEFILLING"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-COFINITE-THM-COFINITE", "MRC-COFINITE-LEM-TANGLE",
                           "MRC-COFINITE-SPHERE-VOLUME-ADDITIVITY", "MRC-COFINITE-LEM-SUMVOLUME",
                           "MRC-PRIME-GAP-THM-MAIN", "MRC-ALL-COVERS-THM-ALLCOVERS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        both = self.graph.simulate({"EXT-BC-ONE-BOUNDARY-NONGRAPH-FILLING": "stale",
                                    "MRC-ALL-COVERS-ALL-NONGRAPH": "stale"})
        self.assertEqual(both["nodes"]["MRC-COFINITE-THM-COFINITE"]["state"], "stale")

    def test_cofinite_initial_volume_is_not_a_general_or_all_parameter_topology_premise(self):
        scenario = self.graph.simulate({"MRC-NEW17-CLOSED-REPRESENTATION-VOLUME": "stale"})
        for conclusion in ("MRC-NEW17-CLOSED-REPRESENTATION-VOLUME", "MRC-COFINITE-LEM-SUMVOLUME"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-NEW17-THM-NONGRAPH", "MRC-COFINITE-THM-COFINITE",
                           "MRC-COFINITE-LEM-PERSISTENCE", "MRC-COFINITE-THM-ONEFILLING",
                           "MRC-COFINITE-GRAPH-TYPE-VOLUME-ZERO", "MRC-PRIME-GAP-THM-MAIN",
                           "MRC-ALL-COVERS-THM-ALLCOVERS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        both = self.graph.simulate({"MRC-NEW17-CLOSED-REPRESENTATION-VOLUME": "stale",
                                    "MRC-ALL-COVERS-ALL-NONGRAPH": "stale"})
        for conclusion in ("MRC-NEW17-THM-NONGRAPH", "MRC-COFINITE-THM-COFINITE"):
            self.assertEqual(both["nodes"][conclusion]["state"], "stale")

    def test_cofinite_topology_is_separate_from_augmentation_gap_consequences(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-KCH-AUGMENTATION-CORRESPONDENCE": "stale"})
        for conclusion in ("MRC-COFINITE-THM-COFINITE", "MRC-COFINITE-TREFOIL-TOPOLOGY"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        for conclusion in ("MRC-COFINITE-PRIME-GAP-KNOTS", "MRC-COFINITE-GAP-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_original_cgls_context_is_not_an_indispensable_cofinite_premise(self):
        scenario = self.graph.simulate({"SRC-MRC-CGLS-ORIGINAL-CONTEXT": "stale"})
        for conclusion in ("MRC-COFINITE-LEM-PERSISTENCE", "MRC-COFINITE-THM-ONEFILLING",
                           "MRC-COFINITE-THM-COFINITE", "MRC-COFINITE-PRIME-GAP-KNOTS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_fixed_rational_tangle_binding_has_an_independent_all_parameter_route(self):
        scenario = self.graph.simulate({"EXT-MZ-RATIONAL-TANGLE-FILLINGS": "stale"})
        self.assertEqual(scenario["nodes"]["MRC-COFINITE-LEM-TANGLE"]["state"], "stale")
        for conclusion in ("MRC-COFINITE-THM-COFINITE", "MRC-COFINITE-LEM-SUMVOLUME",
                           "MRC-COFINITE-THM-ONEFILLING", "MRC-PRIME-GAP-THM-MAIN",
                           "MRC-ALL-COVERS-THM-BINDING", "MRC-ALL-COVERS-THM-ALLCOVERS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        both = self.graph.simulate({"EXT-MZ-RATIONAL-TANGLE-FILLINGS": "stale",
                                    "MRC-ALL-COVERS-ALL-NONGRAPH": "stale"})
        self.assertEqual(both["nodes"]["MRC-COFINITE-THM-COFINITE"]["state"], "stale")

    def test_all_parameter_binding_is_separate_from_fixed_parent_geometry(self):
        for premise in ("MRC-ALL-COVERS-THM-BINDING", "SRC-MRC-ALL-COVERS-TOPOLOGY-CERTIFICATES",
                        "EXT-SNAPPY332-ORIENTED-CURVE-BASIS", "EXT-SNAPPY332-MARKED-PERMANENT-FILL"):
            with self.subTest(premise=premise):
                scenario = self.graph.simulate({premise: "stale"})
                for conclusion in ("MRC-ALL-COVERS-THM-BINDING", "MRC-ALL-COVERS-THM-ALLCOVERS",
                                   "MRC-ALL-COVERS-ALL-NONGRAPH", "MRC-ALL-COVERS-COR-FAMILY",
                                   "MRC-ALL-COVERS-GAP-SUMS"):
                    self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
                for conclusion in ("MRC-ALL-COVERS-COMPLETE-PARENT", "MRC-ALL-COVERS-EMBEDDED-CUSP-LENGTH",
                                   "MRC-ALL-COVERS-LEM-GEOMETRY", "MRC-NEW17-THM-NONGRAPH",
                                   "MRC-COFINITE-THM-COFINITE", "MRC-PRIME-GAP-SHARP"):
                    self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_all_parameter_six_theorem_is_not_a_rank_or_fixed_geometry_premise(self):
        scenario = self.graph.simulate({"EXT-LACKENBY-SIX-THEOREM": "stale"})
        for conclusion in ("MRC-ALL-COVERS-THM-ALLCOVERS", "MRC-ALL-COVERS-ALL-NONGRAPH",
                           "MRC-ALL-COVERS-COR-FAMILY", "MRC-ALL-COVERS-GAP-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-ALL-COVERS-THM-BINDING", "MRC-ALL-COVERS-COMPLETE-PARENT",
                           "MRC-ALL-COVERS-EMBEDDED-CUSP-LENGTH", "MRC-ALL-COVERS-LEM-GEOMETRY",
                           "MRC-PRIME-GAP-SHARP", "MRC-PRIME-GAP-AUGMENTATION",
                           "MRC-NEW17-CLOSED-REPRESENTATION-VOLUME", "MRC-NEW17-THM-NONGRAPH",
                           "MRC-COFINITE-THM-COFINITE"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_all_parameter_cusp_import_is_separate_from_complete_shapes_and_ranks(self):
        scenario = self.graph.simulate({"EXT-SNAPPY332-EMBEDDED-CUSP": "stale"})
        for conclusion in ("MRC-ALL-COVERS-EMBEDDED-CUSP-LENGTH", "MRC-ALL-COVERS-LEM-GEOMETRY",
                           "MRC-ALL-COVERS-THM-ALLCOVERS", "MRC-ALL-COVERS-COR-FAMILY",
                           "MRC-ALL-COVERS-GAP-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-ALL-COVERS-THM-BINDING", "MRC-ALL-COVERS-COMPLETE-PARENT",
                           "MRC-PRIME-GAP-SHARP", "MRC-PRIME-GAP-AUGMENTATION",
                           "MRC-NEW17-THM-NONGRAPH", "MRC-COFINITE-THM-COFINITE"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_all_parameter_new_analytic_imports_preserve_the_old_qualitative_route(self):
        for premise in ("EXT-SNAPPY332-VERIFIED-COMPLETE", "SRC-MRC-ALL-COVERS-GEOMETRY-CERTIFICATES"):
            with self.subTest(premise=premise):
                scenario = self.graph.simulate({premise: "stale"})
                for conclusion in ("MRC-ALL-COVERS-COMPLETE-PARENT", "MRC-ALL-COVERS-EMBEDDED-CUSP-LENGTH",
                                   "MRC-ALL-COVERS-LEM-GEOMETRY", "MRC-ALL-COVERS-THM-ALLCOVERS",
                                   "MRC-ALL-COVERS-ALL-NONGRAPH"):
                    self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
                for conclusion in ("MRC-ALL-COVERS-THM-BINDING", "MRC-NEW17-CLOSED-REPRESENTATION-VOLUME",
                                   "MRC-NEW17-THM-NONGRAPH", "MRC-COFINITE-LEM-PERSISTENCE",
                                   "MRC-COFINITE-THM-ONEFILLING", "MRC-COFINITE-THM-COFINITE"):
                    self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_all_parameter_theorem_does_not_use_old_seed_volume_or_qualitative_fillings(self):
        scenario = self.graph.simulate({"MRC-NEW17-CLOSED-REPRESENTATION-VOLUME": "stale",
                                       "EXT-BC-ONE-BOUNDARY-NONGRAPH-FILLING": "stale",
                                       "EXT-MZ-RATIONAL-TANGLE-FILLINGS": "stale"})
        for conclusion in ("MRC-ALL-COVERS-THM-BINDING", "MRC-ALL-COVERS-LEM-GEOMETRY",
                           "MRC-ALL-COVERS-THM-ALLCOVERS", "MRC-ALL-COVERS-ALL-NONGRAPH",
                           "MRC-ALL-COVERS-COR-FAMILY", "MRC-ALL-COVERS-GAP-SUMS",
                           "MRC-NEW17-THM-NONGRAPH", "MRC-COFINITE-THM-COFINITE"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
        for conclusion in ("MRC-NEW17-CLOSED-REPRESENTATION-VOLUME", "MRC-COFINITE-LEM-SUMVOLUME",
                           "MRC-COFINITE-LEM-PERSISTENCE", "MRC-COFINITE-THM-ONEFILLING",
                           "MRC-COFINITE-LEM-TANGLE"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")

    def test_all_parameter_augmentation_import_affects_gap_consequences_not_cover_topology(self):
        scenario = self.graph.simulate({"EXT-CORNWELL-KCH-AUGMENTATION-CORRESPONDENCE": "stale"})
        for conclusion in ("MRC-ALL-COVERS-COR-FAMILY", "MRC-ALL-COVERS-GAP-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-ALL-COVERS-THM-BINDING", "MRC-ALL-COVERS-LEM-GEOMETRY",
                           "MRC-ALL-COVERS-THM-ALLCOVERS", "MRC-ALL-COVERS-ALL-NONGRAPH",
                           "MRC-ALL-COVERS-TREFOIL-TOPOLOGY", "MRC-PRIME-GAP-SHARP",
                           "MRC-NEW17-THM-NONGRAPH", "MRC-COFINITE-THM-COFINITE"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

    def test_all_parameter_shared_krawczyk_and_cover_imports_are_real_premises(self):
        for premise in ("EXT-NEW17-KRAWCZYK-INTERVAL", "EXT-NEW17-COVER-PERIPHERAL"):
            with self.subTest(premise=premise):
                scenario = self.graph.simulate({premise: "stale"})
                for conclusion in ("MRC-ALL-COVERS-THM-ALLCOVERS", "MRC-ALL-COVERS-ALL-NONGRAPH",
                                   "MRC-ALL-COVERS-COR-FAMILY", "MRC-ALL-COVERS-GAP-SUMS",
                                   "MRC-NEW17-THM-NONGRAPH"):
                    self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
                for conclusion in ("MRC-PRIME-GAP-SHARP", "MRC-PRIME-GAP-AUGMENTATION"):
                    self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")
                if premise == "EXT-NEW17-KRAWCZYK-INTERVAL":
                    self.assertEqual(scenario["nodes"]["MRC-ALL-COVERS-THM-BINDING"]["state"], "proved")
                    self.assertEqual(scenario["nodes"]["MRC-ALL-COVERS-COMPLETE-PARENT"]["state"], "stale")
                else:
                    self.assertEqual(scenario["nodes"]["MRC-ALL-COVERS-THM-BINDING"]["state"], "stale")
                    self.assertEqual(scenario["nodes"]["MRC-ALL-COVERS-LEM-GEOMETRY"]["state"], "proved")

    def test_all_parameter_revocation_preserves_old_volume_and_general_filling_proofs(self):
        scenario = self.graph.simulate({"MRC-ALL-COVERS-THM-ALLCOVERS": "stale"})
        for conclusion in ("MRC-ALL-COVERS-ALL-NONGRAPH", "MRC-ALL-COVERS-COR-FAMILY",
                           "MRC-ALL-COVERS-GAP-SUMS"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "stale")
        for conclusion in ("MRC-NEW17-CLOSED-REPRESENTATION-VOLUME", "MRC-NEW17-THM-NONGRAPH",
                           "MRC-COFINITE-LEM-SUMVOLUME", "MRC-COFINITE-LEM-PERSISTENCE",
                           "MRC-COFINITE-THM-ONEFILLING", "MRC-COFINITE-THM-COFINITE"):
            self.assertEqual(scenario["nodes"][conclusion]["state"], "proved")

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
