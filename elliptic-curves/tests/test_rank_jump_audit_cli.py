from pathlib import Path
import json
import subprocess
import sys
import tempfile
import unittest

CAS = Path(__file__).resolve().parents[1]/"cas"
SCRIPT = CAS/"audit_rank_jump.py"


class AuditCommandTests(unittest.TestCase):
    def run_command(self, *args):
        return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                              capture_output=True, text=True, timeout=10)

    def test_family_from_retained_anchor_format(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)/"anchor.json"
            source.write_text(json.dumps({"anchor": {"id": "test", "base_polynomial_ascending": ["10", "-7", "0", "1"]}}))
            result = self.run_command("family", "--input", source)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(json.loads(result.stdout)["arithmetic_generic_rank"], 0)

    def test_outputs_are_immutable(self):
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp)/"null.json"
            args = ("ct-null", "--dimension", "13", "--output", target)
            self.assertEqual(self.run_command(*args).returncode, 0)
            original = target.read_bytes()
            self.assertNotEqual(self.run_command(*args).returncode, 0)
            self.assertEqual(target.read_bytes(), original)

    def test_mask_outputs_do_not_mix_oracle_into_request(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source, request, oracle = tmp/"input.json", tmp/"request.json", tmp/"oracle.json"
            source.write_text(json.dumps({"curve": [-7, 10], "points": [[1, 2], [2, 2]], "metric_gram": [[4, 1], [1, 6]]}))
            result = self.run_command("mask", "--input", source, "--withhold", "1", "--search-input", request, "--oracle", oracle)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertNotIn("withheld_points", json.loads(request.read_text()))
            self.assertEqual(json.loads(oracle.read_text())["withheld_points"], [{"x": "2", "y": "2"}])

    def test_bad_mask_must_not_create_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            source = tmp/"input.json"
            source.write_text(json.dumps({"curve": [-7, 10], "points": [[1, 2], [2, 2]], "metric_gram": [[4, 1], [1, 6]]}))
            request, oracle = tmp/"request.json", tmp/"oracle.json"
            result = self.run_command("mask", "--input", source, "--withhold", "0,1", "--search-input", request, "--oracle", oracle)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(request.exists())
            self.assertFalse(oracle.exists())

    def test_ct_null_label(self):
        result = self.run_command("ct-null", "--dimension", "17")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("not an elliptic-curve sampling law", json.loads(result.stdout)["model"])


if __name__ == "__main__":
    unittest.main()
