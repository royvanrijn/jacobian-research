#!/usr/bin/env python3
"""Focused exact specialization-relation replay for the rational component."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
SCRIPT = CAS / "audit_mestre_transverse_component_relations.py"


def load() -> object:
    sys.path.insert(0, str(CAS))
    spec = importlib.util.spec_from_file_location("mestre_component_relations", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class MestreTransverseComponentRelationsTest(unittest.TestCase):
    def test_exact_visible_subgroup_relations(self) -> None:
        audit = load()
        result = audit.replay(include_pari=False)
        self.assertEqual(result["visible_mod3_rank"], 9)
        self.assertEqual(result["augmented_mod3_rank"], 9)
        self.assertEqual(
            result["exact_affine_relations"],
            [
                "P1=V(0,-)-V(0,+)+V(1,+)",
                "P2=-V(0,-)-V(1,+)-V((42-z)/6,-)",
            ],
        )


if __name__ == "__main__":
    unittest.main()
