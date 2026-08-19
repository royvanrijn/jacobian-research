from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from elliptic_candidate_record import point_sequence_sha256  # noqa: E402
from search_fermigier_rank20_direction import (  # noqa: E402
    AnchorData,
    EXPECTED_STRONG_BASIS_SHA256,
    canonical_sign,
    load_anchor,
    rational_square_root,
    scan_frontier,
    score_base_point,
)


Q = Fraction


class FermigierRank20DirectionTests(unittest.TestCase):
    def test_pinned_anchor_loads_the_stronger_basis(self) -> None:
        anchor = load_anchor()
        self.assertEqual(len(anchor.basis), 20)
        self.assertEqual(anchor.basis_sha256, EXPECTED_STRONG_BASIS_SHA256)
        self.assertGreaterEqual(len(anchor.known_signless_points), 115)
        self.assertEqual(anchor.model[:3], (Q(0), Q(0), Q(0)))

    def test_exact_helpers(self) -> None:
        self.assertEqual(rational_square_root(Q(81, 100)), Q(9, 10))
        with self.assertRaises(ValueError):
            rational_square_root(Q(2, 3))
        self.assertEqual(canonical_sign((Q(7, 3), Q(-11, 5))), (Q(7, 3), Q(11, 5)))

    def test_interrupted_frontier_resumes_from_consistent_gray_state(self) -> None:
        model = (Q(0), Q(0), Q(0), Q(-1), Q(1))
        basis = ((Q(0), Q(1)), (Q(1), Q(1)))
        anchor = AnchorData(
            model=model,
            basis=basis,
            basis_sha256=point_sequence_sha256(basis),
            known_signless_points=frozenset(canonical_sign(point) for point in basis),
            known_quartic_abscissas=(),
            candidate_record={},
        )
        original = score_base_point
        calls = 0

        def interrupt_second(coefficient_a, current_basis, base_point):
            nonlocal calls
            calls += 1
            if calls == 2:
                raise KeyboardInterrupt
            return original(coefficient_a, current_basis, base_point)

        with tempfile.TemporaryDirectory() as directory:
            checkpoint = Path(directory) / "frontier.json"
            with patch(
                "search_fermigier_rank20_direction.score_base_point",
                side_effect=interrupt_second,
            ):
                with self.assertRaises(KeyboardInterrupt):
                    scan_frontier(
                        anchor,
                        checkpoint,
                        retain_count=3,
                        global_count=2,
                        maximum_count=3,
                        progress_interval=1,
                    )
            partial = json.loads(checkpoint.read_text())
            self.assertEqual(partial["status"], "partial")
            self.assertEqual(partial["state"]["processed_integer"], 1)
            self.assertEqual(partial["state"]["previous_gray"], 1)

            selected = scan_frontier(
                anchor,
                checkpoint,
                retain_count=3,
                global_count=2,
                maximum_count=3,
                progress_interval=1,
            )
            completed = json.loads(checkpoint.read_text())
            self.assertEqual(completed["status"], "complete")
            self.assertEqual(completed["state"]["processed_integer"], 3)
            self.assertGreaterEqual(len(selected), 2)


if __name__ == "__main__":
    unittest.main()
