from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest


TOOLS_DIRECTORY = Path(__file__).resolve().parents[1] / "tools"
REPOSITORY = Path(__file__).resolve().parents[2]
GENERATOR = TOOLS_DIRECTORY / "generate_u42_magma.py"
PARI_RUNNER = TOOLS_DIRECTORY / "run_u42_pari_rank.py"
PROBE_ARTIFACT = (
    REPOSITORY
    / "artifacts/generated-results/elliptic_nagao_u42_magma_probe.json"
)
TOOLCHAIN_ARTIFACT = (
    REPOSITORY
    / "artifacts/generated-results/elliptic_nagao_u42_descent_toolchain.json"
)
if str(TOOLS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(TOOLS_DIRECTORY))

from generate_u42_magma import (  # noqa: E402
    EXPECTED_MINIMAL,
    SOURCE_SHA256,
    load_points,
    probe_code,
    setup_code,
)
from run_u42_pari_rank import (  # noqa: E402
    gp_program as pari_gp_program,
    load_points as pari_load_points,
)


class NagaoU42MagmaGeneratorTests(unittest.TestCase):
    def test_pinned_basis_is_complete(self) -> None:
        points = load_points()
        self.assertEqual(len(points), 17)
        self.assertEqual(
            points[0],
            ("5609374266179619/9604", "-31644788272829862624/2401"),
        )
        self.assertEqual(
            SOURCE_SHA256,
            "4fea0207fd637988bcc1147143657cbec5c2404cb81b4c4a487e2dde20cc43b8",
        )

    def test_setup_verifies_before_mapping_to_the_minimal_model(self) -> None:
        code = setup_code(load_points())
        equation_check = code.index("P[2]^2 eq P[1]^3")
        minimal_model = code.index("MinimalModel(Eshort)")
        independence = code.index("IsLinearlyIndependent(Pmin)")
        self.assertLess(equation_check, minimal_model)
        self.assertLess(minimal_model, independence)
        self.assertIn("Pmin := [ short_to_min(P) : P in Pshort ]", code)
        self.assertIn(
            ", ".join(str(value) for value in EXPECTED_MINIMAL),
            code,
        )

    def test_unconditional_selmer_probe_is_explicit_and_guarded(self) -> None:
        code = setup_code(load_points()) + probe_code("twoselmer")
        self.assertIn("TwoSelmerGroup(Emin : Bound := -1)", code)
        self.assertIn("procedure RunTwoSelmer()", code)
        self.assertLess(len(code.encode()), 50_000)

    def test_quotient_descent_factors_out_known_points(self) -> None:
        code = setup_code(load_points()) + probe_code("twodescent")
        self.assertIn("RemoveTorsion := true", code)
        self.assertIn("RemoveGens := known", code)
        self.assertIn("WithMaps := false", code)
        self.assertIn("procedure RunTwoDescent()", code)
        self.assertLess(len(code.encode()), 50_000)

    def test_rankbounds_probe_has_no_unknown_verbose_flag(self) -> None:
        code = probe_code("rankbounds")
        self.assertIn("RankBounds(Emin : Effort := 1)", code)
        self.assertNotIn('SetVerbose("EC"', code)

    def test_probe_artifact_pins_generator_and_does_not_claim_an_upper_bound(self) -> None:
        artifact = json.loads(PROBE_ARTIFACT.read_text())
        generator_sha256 = hashlib.sha256(GENERATOR.read_bytes()).hexdigest()
        self.assertEqual(artifact["generator"]["sha256"], generator_sha256)
        self.assertEqual(
            artifact["verification"]["rigorous_rank_lower_bound"], 17
        )
        self.assertEqual(
            artifact["upper_bound_probes"]["twoselmer"]["result"],
            "no Selmer group or upper bound returned",
        )

    def test_pari_runner_supplies_the_same_exact_basis(self) -> None:
        self.assertEqual(pari_load_points(), load_points())
        program = pari_gp_program(8_000_000_000)
        self.assertIn("EXACT_POINTS_ON_CURVE", program)
        self.assertIn("R=ellrank(E,0,P);", program)
        self.assertIn("default(parisizemax,8000000000);", program)

    def test_toolchain_artifact_records_bounded_failures_only(self) -> None:
        artifact = json.loads(TOOLCHAIN_ARTIFACT.read_text())
        runner_sha256 = hashlib.sha256(PARI_RUNNER.read_bytes()).hexdigest()
        self.assertEqual(artifact["pari_ellrank"]["runner_sha256"], runner_sha256)
        self.assertEqual(artifact["pari_ellrank"]["status"], "strict timeout")
        self.assertIsNone(artifact["summary"]["certified_rank_upper_bound"])
        self.assertIsNone(artifact["summary"]["new_rank_claim"])


if __name__ == "__main__":
    unittest.main()
