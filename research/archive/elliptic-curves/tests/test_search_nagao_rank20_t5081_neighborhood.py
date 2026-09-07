from __future__ import annotations

from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
GENERATED = ROOT / "artifacts/generated-results"
sys.path.insert(0, str(CAS))

from search_nagao_rank20_t5081_neighborhood import (  # noqa: E402
    CALIBRATION_PARAMETER,
    DISCRIMINANT_POLYNOMIAL,
    EXPECTED_TRACE_DESIGN,
    INVARIANT_I,
    INVARIANT_J,
    PROXY_LIMIT,
    QUADRATIC_COMPANION_RELATIONS,
    SAVING_PRIMES,
    build_beam_strata,
    build_residue_tables,
    build_trace_beams,
    companion_quartic_points,
    exact_predeclared_seeds,
    gauss_shell,
    generate_candidates,
    homogenized_discriminant,
    integer_valuation,
    learn_discriminant_root_balls,
    learn_local_trace_fingerprint,
    outside_prior_box,
    polynomial_value,
    quadratic_companion_quartic_points,
    smooth_even_denominators,
)
from certify_nagao_rank20_t5081 import CONSTRUCTION  # noqa: E402
from nagao_1994 import primitive_quartic_coefficients, quartic_value  # noqa: E402


SCRIPT = CAS / "search_nagao_rank20_t5081_neighborhood.py"
ARTIFACT = GENERATED / "elliptic_nagao_rank20_t5081_neighborhood.json"


class NagaoRank20T5081NeighborhoodTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tables = build_residue_tables(200)
        cls.fingerprint = learn_local_trace_fingerprint(cls.tables)
        cls.root_balls = learn_discriminant_root_balls()

    def test_learned_invariants_and_discriminant_are_pinned(self) -> None:
        self.assertEqual((len(INVARIANT_I), len(INVARIANT_J)), (9, 13))
        self.assertTrue(all(INVARIANT_I[index] == 0 for index in range(1, 9, 2)))
        self.assertTrue(all(INVARIANT_J[index] == 0 for index in range(1, 13, 2)))
        self.assertEqual(len(DISCRIMINANT_POLYNOMIAL), 21)
        self.assertEqual(
            hashlib.sha256(",".join(map(str, DISCRIMINANT_POLYNOMIAL)).encode()).hexdigest(),
            "347a8f04bf8184a2e969b89035ba02c2a7bce2cb497a5bd8f5800e78edf4d55c",
        )
        for parameter in (Q(1), Q(17, 3), CALIBRATION_PARAMETER):
            invariant_i, invariant_j = CONSTRUCTION.primitive_binary_invariants(parameter)
            self.assertEqual(polynomial_value(INVARIANT_I, parameter), invariant_i)
            self.assertEqual(polynomial_value(INVARIANT_J, parameter), invariant_j)

    def test_trace_fingerprint_is_automatically_pinned(self) -> None:
        self.assertEqual(self.fingerprint, EXPECTED_TRACE_DESIGN)
        self.assertEqual([row[0] for row in self.fingerprint], [53, 109, 151, 163, 197])
        self.assertTrue(all(row[2] < 0 for row in self.fingerprint))

    def test_discriminant_root_ball_unions_are_exact(self) -> None:
        expected_counts = {2: 2, 3: 2, 5: 12, 7: 14, 13: 6, 17: 5, 23: 8}
        expected_targets = {2: 6, 3: 14, 5: 7, 7: 6, 13: 2, 17: 2, 23: 2}
        self.assertEqual({prime: len(value) for prime, value in self.root_balls.items()}, expected_counts)
        lead_discriminant = homogenized_discriminant(CALIBRATION_PARAMETER)
        for prime in SAVING_PRIMES:
            self.assertEqual(integer_valuation(lead_discriminant, prime), expected_targets[prime])
            self.assertTrue(
                any(
                    (
                        CALIBRATION_PARAMETER.numerator
                        - ball.residue * CALIBRATION_PARAMETER.denominator
                    )
                    % ball.modulus
                    == 0
                    for ball in self.root_balls[prime]
                )
            )
            self.assertTrue(
                all(ball.forced_valuation >= expected_targets[prime] for ball in self.root_balls[prime])
            )

    def test_all_nine_companion_sections_are_exact_and_predeclared(self) -> None:
        for parameter in (Q(1), CALIBRATION_PARAMETER):
            quartic = primitive_quartic_coefficients(CONSTRUCTION, parameter)
            companions = companion_quartic_points(parameter) + quadratic_companion_quartic_points(parameter)
            self.assertEqual(len(companions), 9)
            self.assertTrue(
                all(y_value**2 == quartic_value(quartic, x_value) for x_value, y_value in companions)
            )
            seed_quartic, seed_images, _ = exact_predeclared_seeds(parameter)
            self.assertEqual((len(seed_quartic), len(seed_images)), (21, 21))
        self.assertEqual(len(QUADRATIC_COMPANION_RELATIONS), 3)
        self.assertTrue(all(len(relation) == 12 for relation in QUADRATIC_COMPANION_RELATIONS))

    def test_gauss_shell_replays_class_up_to_even_sign(self) -> None:
        shell = gauss_shell(37, 101, limit=6)
        self.assertGreaterEqual(len(shell), 2)
        for parameter, _, _ in shell:
            plus = (parameter.numerator - 37 * parameter.denominator) % 101
            minus = (parameter.numerator + 37 * parameter.denominator) % 101
            self.assertTrue(plus == 0 or minus == 0)

    def test_population_exclusion_and_smooth_even_stratum(self) -> None:
        self.assertFalse(outside_prior_box(Q(9999, 100)))
        self.assertTrue(outside_prior_box(Q(9998, 101)))
        self.assertTrue(outside_prior_box(Q(10001, 100)))
        denominators = smooth_even_denominators()
        self.assertEqual(len(denominators), 35)
        self.assertTrue(all(101 <= denominator <= 2000 and denominator % 16 == 0 for denominator in denominators))

    def test_small_generation_is_bounded_new_and_leak_free(self) -> None:
        full, prefix, _ = build_trace_beams(self.fingerprint, self.tables, width=30)
        strata = build_beam_strata(full, prefix, self.root_balls, root_width=8)
        candidates, audit = generate_candidates(
            strata, proxy_limit=PROXY_LIMIT, max_survivors=40
        )
        self.assertGreater(len(candidates), 0)
        self.assertLessEqual(len(candidates), 40)
        self.assertEqual(audit["retained_after_stratified_cap"], len(candidates))
        self.assertTrue(all(candidate.parameter > 0 for candidate in candidates))
        self.assertTrue(all(candidate.parameter != CALIBRATION_PARAMETER for candidate in candidates))
        self.assertTrue(all(outside_prior_box(candidate.parameter) for candidate in candidates))

    def test_generated_artifact_is_pinned_when_present(self) -> None:
        if not ARTIFACT.exists():
            self.skipTest("the bounded section-7 neighborhood artifact is absent")
        data = json.loads(ARTIFACT.read_text())
        self.assertEqual(
            data["reproducibility"]["script_sha256"],
            hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
        )
        self.assertEqual(data["calibration_only"]["constructor_parameter_T"], "5081/47")
        self.assertTrue(data["calibration_only"]["excluded_before_every_population_selection"])
        self.assertEqual(
            [row["prime"] for row in data["local_design"]["learned_trace_fingerprint"]],
            [53, 109, 151, 163, 197],
        )
        self.assertEqual(data["point_triage"]["predeclared_sections_decontaminated"], 21)


if __name__ == "__main__":
    unittest.main()
