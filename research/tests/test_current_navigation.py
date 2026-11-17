"""Detect stale human-facing navigation without replaying mathematics."""

from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import research


class CurrentNavigationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name) / "research"
        self.root.mkdir()
        change = patch.object(research, "ROOT", self.root)
        change.start()
        self.addCleanup(change.stop)
        self.data = {"topics": [{
            "id": "moments", "title": "Moments", "summary": "Current mathematics.",
            "status": "The larger question remains open.", "readme": "moments/README.md",
            "results": [{"id": "two-pair", "title": "Two-pair counterexample",
                         "state": "proved", "scope": "A scoped theorem.",
                         "canonical_source": "proof.md"}],
            "open_questions": [{"question": "Does the larger statement hold?",
                                "source": "proof.md"}],
        }]}

    def write_navigation(self, data):
        for path, text in research.navigation_outputs(data).items():
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)

    def test_stale_research_index_or_topic_guide_is_rejected(self):
        for relative in ("README.md", "moments/README.md"):
            with self.subTest(page=relative):
                self.write_navigation(self.data)
                (self.root / relative).write_text("Outdated status.\n")
                with self.assertRaisesRegex(AssertionError, "is stale"):
                    research.check_navigation(self.data)

    def test_changed_result_state_cannot_leave_the_topic_table_behind(self):
        self.write_navigation(self.data)
        changed = deepcopy(self.data)
        changed["topics"][0]["results"][0]["state"] = "partial"
        for path, text in research.navigation_outputs(changed).items():
            if path != self.root / "moments/README.md":
                path.write_text(text)
        with self.assertRaisesRegex(AssertionError, "moments/README.md is stale"):
            research.check_navigation(changed)
        self.write_navigation(changed)
        research.check_navigation(changed)

    def test_rendering_cannot_target_a_published_file(self):
        self.data["topics"][0]["readme"] = "../results/paper/main.tex"
        with self.assertRaisesRegex(AssertionError, "topic guide must be its own README"):
            research.navigation_outputs(self.data)


if __name__ == "__main__":
    unittest.main()
