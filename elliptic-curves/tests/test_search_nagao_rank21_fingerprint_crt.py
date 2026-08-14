from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
CAS = ROOT / "elliptic-curves/cas"
sys.path.insert(0, str(CAS))

from search_nagao_rank21_fingerprint_crt import (  # noqa: E402
    CALIBRATION_PARAMETER,
    DESIGN_PRIMES,
    POWER_PRIMES,
    TRACE_THRESHOLDS,
    build_power_strata,
    build_trace_beam,
    calibration_verification,
    cheap_power_symbol_union,
    gauss_representatives,
    in_old_positive_rectangle,
    leakage_free_tables,
    power_symbol_union,
    select_point_survivors,
    smooth_stratum_closure,
    trace_symbol_union,
)


class NagaoRank21FingerprintCRTTests(unittest.TestCase):
    def test_trace_unions_and_calibration_are_pinned(self) -> None:
        expected_sizes = {11: 6, 19: 10, 31: 20, 41: 20, 47: 12, 59: 39}
        for prime, threshold in TRACE_THRESHOLDS:
            symbols = trace_symbol_union(prime, threshold)
            self.assertEqual(len(symbols), expected_sizes[prime])
            self.assertTrue(all(symbol.ellap <= threshold for symbol in symbols))
        records = calibration_verification()
        self.assertEqual(CALIBRATION_PARAMETER, Q(6793, 64))
        self.assertEqual([record["prime"] for record in records], [11, 19, 31, 41, 47, 59])
        self.assertTrue(all(not record["used_for_generation"] for record in records))

    def test_gauss_representatives_replay_the_crt_class(self) -> None:
        representatives = gauss_representatives(37, 101, limit=4)
        self.assertGreaterEqual(len(representatives), 2)
        for parameter, _, _ in representatives:
            self.assertEqual(
                (parameter.numerator - 37 * parameter.denominator) % 101,
                0,
            )
        self.assertTrue(in_old_positive_rectangle(Q(-85, 2)))
        self.assertFalse(in_old_positive_rectangle(Q(-10001, 2)))

    def test_power_symbols_are_exact_multiple_root_balls(self) -> None:
        unions = {prime: power_symbol_union(prime) for prime in POWER_PRIMES}
        self.assertTrue(all(symbol.forced_valuation >= 2 for values in unions.values() for symbol in values))
        self.assertTrue(any(symbol.modulus == 13**2 and symbol.forced_valuation == 5 for symbol in unions[13]))
        self.assertTrue(any(symbol.modulus == 83**2 and symbol.forced_valuation == 3 for symbol in unions[83]))

    def test_small_beam_and_power_strata_are_exactly_bounded(self) -> None:
        states, audit = build_trace_beam(width=40)
        self.assertEqual(len(states), 40)
        self.assertEqual(len(audit), len(TRACE_THRESHOLDS))
        strata = build_power_strata(states, single_width=12, double_width=6, triple_width=3)
        self.assertEqual(len(strata["trace-only"]), 40)
        self.assertTrue(all(len(value) > 0 for value in strata.values()))
        self.assertNotIn("power-13-37-83", strata)
        self.assertIn("even-denominator-v2-22", strata)

    def test_cheap_power_unions_are_not_single_residue_clones(self) -> None:
        self.assertEqual(len(cheap_power_symbol_union(5)), 4)
        self.assertEqual(len(cheap_power_symbol_union(23)), 12)
        self.assertEqual(len(cheap_power_symbol_union(7)), 16)

    def test_leakage_tables_and_closed_smooth_stratum(self) -> None:
        tables = leakage_free_tables(200)
        self.assertFalse(set(tables) & set(DESIGN_PRIMES))
        records = smooth_stratum_closure()
        self.assertEqual(len(records), 20)
        self.assertTrue(all(not record["point_search_in_this_lane"] for record in records))
        self.assertGreater(
            min(record["radical_proxy"]["log_radical_upper_proxy"] for record in records),
            190,
        )

    def test_point_selection_preserves_declared_b200_tail(self) -> None:
        # A structural smoke test: empty input is valid and deterministic.
        self.assertEqual(select_point_survivors((), keep_count=4), ())


if __name__ == "__main__":
    unittest.main()
