from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from explicit_formula_rank_diagnostic import (  # noqa: E402
    DELTA,
    E20_MODEL,
    PRIME_LIMIT,
    gp_program,
)
from explicit_formula_section7_rank17 import (  # noqa: E402
    CANDIDATES,
    calibrated_bounds,
    load_exact_inputs,
)
from nagao_1994 import short_jacobian_coefficients  # noqa: E402
from nagao_1994_section7 import SECTION7_CONSTRUCTION  # noqa: E402


SCRIPT = CAS / "explicit_formula_section7_rank17.py"
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/elliptic_nagao_section7_rank17_explicit_formula.json"
)


class ExplicitFormulaSection7Rank17Tests(unittest.TestCase):
    def test_exact_inputs_are_the_two_certified_rank17_fibers(self) -> None:
        inputs = load_exact_inputs()
        self.assertEqual(set(inputs), {"T599", "T426"})
        self.assertEqual(inputs["T599"]["constructor_parameter_T"], "599/2")
        self.assertEqual(inputs["T426"]["constructor_parameter_T"], "426")
        for candidate in inputs.values():
            self.assertEqual(candidate["certified_algebraic_rank_lower_bound"], 17)
            self.assertEqual(candidate["root_number"], -1)

    def test_gp_program_has_all_three_exact_models(self) -> None:
        curves = (
            ("E20", E20_MODEL),
            *(
                (
                    label,
                    short_jacobian_coefficients(SECTION7_CONSTRUCTION, parameter),
                )
                for label, parameter in CANDIDATES
            ),
        )
        program = gp_program(curves)
        self.assertIn(f"D={DELTA};", program)
        self.assertIn(f"LIM={PRIME_LIMIT};", program)
        for label in ("E20", "T599", "T426"):
            self.assertIn(f"ROW|{label}|", program)

    def test_pinned_calibrated_values(self) -> None:
        records = {
            "E20": {
                "log_conductor": "170.0876648422488811354314776110703548643",
                "prime_sum": "-53.5271982465901051823721048180945488965",
            },
            "T599": {
                "log_conductor": "124.0610122569483449003220429890761973897",
                "prime_sum": "-53.0005817587059842090900824239944506477",
            },
            "T426": {
                "log_conductor": "138.8258223152913849024419185729231226275",
                "prime_sum": "-50.5088175321621536892883507397483669490",
            },
        }
        result = calibrated_bounds(records)
        self.assertGreater(float(result["T599"]["calibrated_upper_value"]), 17.95)
        self.assertLess(float(result["T599"]["calibrated_upper_value"]), 17.96)
        self.assertGreater(float(result["T426"]["calibrated_upper_value"]), 18.73)
        self.assertLess(float(result["T426"]["calibrated_upper_value"]), 18.74)

    def test_generated_artifact_matches_script(self) -> None:
        artifact = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        self.assertEqual(
            artifact["script_sha256"], hashlib.sha256(SCRIPT.read_bytes()).hexdigest()
        )
        self.assertEqual(artifact["method"]["reference_published_upper_value"], "21.70")
        for label in ("T599", "T426"):
            self.assertLess(
                float(artifact["comparisons"][label]["calibrated_upper_value"]),
                19,
            )
            self.assertEqual(artifact["exact_inputs"][label]["root_number"], -1)
        self.assertIn(
            "no algebraic-rank upper bound",
            artifact["interpretation"]["unconditional"],
        )


if __name__ == "__main__":
    unittest.main()
