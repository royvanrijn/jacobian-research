"""Pruning is exact-scope, theorem-bound and invalidated by changed evidence."""

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"elliptic-curves/cas"))
from research_runtime.pruning import PruningRegistry,SearchRequest,known_constraints,theorem_fingerprint


class PruningTests(unittest.TestCase):
    def test_actual_product_and_ct_exclusions_without_chart_confusion(self):
        registry=PruningRegistry(ROOT)
        rules=known_constraints(ROOT)["constraints"]
        for rule in rules:
            scope=tuple(rule["scope"].items())
            if rule["kind"]=="rank_upper":
                request=SearchRequest("nontorsion_section",scope,height="1000000000000")
                self.assertFalse(registry.decision(request)["search_allowed"])
            else:
                self.assertFalse(registry.decision(SearchRequest("kummer_class",scope,class_mask=1,class_dimension=20))["search_allowed"])
                for mask in rule["basis_masks"]:
                    self.assertTrue(registry.decision(SearchRequest("kummer_class",scope,class_mask=mask,class_dimension=20))["search_allowed"])
                # Pointed charts are birational models of the whole elliptic
                # curve, so their integer labels are NOT Kummer classes.
                self.assertTrue(registry.decision(SearchRequest("birational_chart",scope,class_mask=1))["search_allowed"])
        self.assertTrue(registry.decision(SearchRequest("nontorsion_section",(("surface","different constant twist"),)))["search_allowed"])

    def test_new_theorem_and_retraction_are_seen_immediately(self):
        with TemporaryDirectory() as directory:
            root=Path(directory);proof=root/"proof.json";proof.write_text('{"rank":0}')
            rule={"scope":{"surface":"s"},"kind":"rank_upper","upper":0,
                  "evidence":{"proof.json":sha256(proof.read_bytes()).hexdigest()}}
            entry={"id":"R0","state":"unknown","scope":"rank zero on s","canonical_source":"proof.json",
                   "search_constraints":[rule]}
            status=root/"MATH_STATUS.json"
            status.write_text(json.dumps({"entries":[entry]}))
            registry=PruningRegistry(root)
            request=SearchRequest("nontorsion_section",(("surface","s"),))
            self.assertTrue(registry.decision(request)["search_allowed"])
            entry["state"]="proved";status.write_text(json.dumps({"entries":[entry]}))
            self.assertFalse(registry.decision(request)["search_allowed"])
            proof.write_text('{"rank":1}')
            self.assertTrue(registry.decision(request)["search_allowed"])
            self.assertIn("changed",registry.inactive[0]["reason"])

    def test_stale_mathematical_statement_never_excludes(self):
        with TemporaryDirectory() as directory:
            root=Path(directory);proof=root/"proof";proof.write_text("witness")
            entry={"id":"T","state":"proved","scope":"old statement","canonical_source":"proof"}
            record={"theorem_id":"T","theorem_fingerprint":theorem_fingerprint(entry),
                "evidence":{"proof":sha256(proof.read_bytes()).hexdigest()},"kind":"rank_upper","upper":0,"scope":{"surface":"s"}}
            registry_path=root/"constraints.json";registry_path.write_text(json.dumps({"constraints":[record]}))
            (root/"MATH_STATUS.json").write_text(json.dumps({"entries":[{**entry,"scope":"corrected statement"}]}))
            registry=PruningRegistry(root,registry_path)
            self.assertTrue(registry.decision(SearchRequest("nontorsion_section",(("surface","s"),)))["search_allowed"])


if __name__=="__main__":unittest.main()
