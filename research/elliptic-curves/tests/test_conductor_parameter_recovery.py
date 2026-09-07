from __future__ import annotations

from pathlib import Path
import shutil
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
PROGRAM_ROOT = ROOT / "elliptic-curves"
if str(PROGRAM_ROOT) not in sys.path:
    sys.path.insert(0, str(PROGRAM_ROOT))

from ecsearch.conductor_engineering import (  # noqa: E402
    ProjectiveRootBall,
    combine_projective_rows,
    discover_projective_root_balls,
    enumerate_projective_lattice,
    exact_j_match,
    factor_over_known_primes,
    select_repeated_prime_constraints,
    weierstrass_invariant_polynomials,
    weierstrass_invariants,
)
from ecsearch.crt_lattice import (  # noqa: E402
    gauss_reduce_linear_congruence_lattice,
)


class ProjectiveFingerprintTests(unittest.TestCase):
    def test_factorization_and_automatic_repeated_power_selection(self) -> None:
        value = -(2**5) * 3 * (5**2) * 101
        factors, cofactor = factor_over_known_primes(value, (2, 3, 5, 101))
        self.assertEqual(factors, ((2, 5), (3, 1), (5, 2), (101, 1)))
        self.assertEqual(cofactor, 1)
        self.assertEqual(
            select_repeated_prime_constraints(
                value,
                maximum_prime=11,
                minimum_valuation=2,
                excluded_primes=(2,),
            ),
            ((5, 2),),
        )

    def test_affine_and_infinity_charts_are_disjoint(self) -> None:
        # H(a,b)=b+5a has no affine root mod 5, but has the infinity
        # class b/a=0 mod 5.
        balls, profile = discover_projective_root_balls((1, 5), 5, 1)
        self.assertEqual(profile["affine_target_root_count"], 0)
        self.assertEqual(profile["infinity_target_root_count"], 1)
        self.assertEqual(
            balls,
            (ProjectiveRootBall(5, 1, "infinity", 1, 0),),
        )

    def test_multiple_roots_compress_to_their_actual_cost(self) -> None:
        balls, profile = discover_projective_root_balls((0, 0, 1), 3, 3)
        self.assertEqual(profile["affine_target_root_count"], 3)
        self.assertIn(ProjectiveRootBall(3, 3, "affine", 2, 0), balls)

    def test_mixed_chart_kernel_contains_the_declared_parameter(self) -> None:
        choices = (
            ProjectiveRootBall(5, 1, "affine", 1, 2),
            ProjectiveRootBall(7, 1, "infinity", 1, 0),
        )
        coefficient_a, coefficient_b, modulus = combine_projective_rows(choices)
        basis = gauss_reduce_linear_congruence_lattice(
            coefficient_a, coefficient_b, modulus
        )
        pairs = set(
            enumerate_projective_lattice(
                basis, modulus, coefficient_radius=4, height_cap=100
            )
        )
        self.assertIn((9, 7), pairs)
        self.assertTrue(all(choice.matches(9, 7) for choice in choices))


class ExactFiberVerificationTests(unittest.TestCase):
    def test_polynomial_weierstrass_invariants_and_j_match(self) -> None:
        model = {
            "a1": (0,),
            "a2": (0,),
            "a3": (0,),
            "a4": (0, 0, -1),
            "a6": (0, 0, 1),
        }
        family = weierstrass_invariant_polynomials(model)
        target = weierstrass_invariants((0, 0, 0, -25, 25))
        self.assertTrue(
            exact_j_match(
                5,
                1,
                family_c4=family["c4"],
                family_discriminant=family["discriminant"],
                target_c4=target["c4"],
                target_discriminant=target["discriminant"],
            )
        )
        self.assertFalse(
            exact_j_match(
                4,
                1,
                family_c4=family["c4"],
                family_discriminant=family["discriminant"],
                target_c4=target["c4"],
                target_discriminant=target["discriminant"],
            )
        )


@unittest.skipUnless(shutil.which("gp"), "PARI/GP is required for the pinned replay")
class Icarm282ReplayTests(unittest.TestCase):
    def test_pinned_recovery_is_current(self) -> None:
        specification = (
            ROOT
            / "elliptic-curves/data/conductor-engineering/icarm_curve282_fermigier.json"
        )
        artifact = (
            ROOT
            / "artifacts/generated-results/elliptic-curves/"
            "icarm_curve282_conductor_parameter_recovery_v1.json"
        )
        completed = subprocess.run(
            [
                sys.executable,
                str(ROOT / "elliptic-curves/scripts/recover_conductor_parameter.py"),
                str(specification),
                "--output",
                str(artifact),
                "--check",
            ],
            text=True,
            capture_output=True,
            check=True,
            timeout=60,
        )
        self.assertIn("exact_j_matches=['11671/42']", completed.stdout)


if __name__ == "__main__":
    unittest.main()
