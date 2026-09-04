from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elliptic-curves/cas/certify_record_f1_from_class_quotient.py"


def load_module():
    spec = importlib.util.spec_from_file_location("record_f1_class_quotient", SOURCE)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RecordF1ClassQuotientTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.module = load_module()

    def test_exact_change_of_generator(self):
        self.assertEqual(
            self.module.TARGETS[356]["pressure_cubic"],
            "z^3 + 4*z^2 - "
            f"{16 * self.module.A356}*z + "
            f"{64 * self.module.B356}",
        )
        self.assertEqual(
            self.module.TARGETS[385]["pressure_cubic"],
            "z^3 - 3*z^2 + "
            f"{3 - self.module.A385}*z + "
            f"{self.module.A385 + self.module.B385 - 1}",
        )

    def test_binary_rank(self):
        self.assertEqual(
            self.module.f2_rank([[1, 0, 1], [0, 1, 1], [1, 1, 0]], 3),
            2,
        )

    def test_one_sided_bad_span_is_subtracted_in_quotient(self):
        self.assertEqual(self.module.quotient_dimension_upper(22, 10), 12)
        with self.assertRaises(ArithmeticError):
            self.module.quotient_dimension_upper(21, 22)
        source = SOURCE.read_text()
        self.assertIn("C'/B' -> C/B is still surjective", source)
        self.assertNotIn("no subtraction is allowed", source)

    def test_strict_marker_parser(self):
        protocol = self.module.PROTOCOL
        log = "\n".join(
            [
                f"{protocol}|stage=class_quotient|status=PASS|cyc=[2, 4, 3]|signature=[1, 1]",
                f"{protocol}|stage=rational_prime|status=PASS|p=2|count=1",
                f"{protocol}|stage=prime_ideal|status=PASS|p=2|index=1|e=1|f=3|coordinates=[1, 3, 2]~",
                f"{protocol}|stage=complete|status=PASS",
            ]
        )
        parsed = self.module.parse_gp_log(log, [2])
        self.assertEqual(parsed["class_group_invariants"], [2, 4, 3])
        self.assertEqual(parsed["prime_ideals"][0]["class_coordinates"], [1, 3, 2])

    def test_parser_rejects_missing_terminal(self):
        with self.assertRaises(ArithmeticError):
            self.module.parse_gp_log("", [2])


if __name__ == "__main__":
    unittest.main()
