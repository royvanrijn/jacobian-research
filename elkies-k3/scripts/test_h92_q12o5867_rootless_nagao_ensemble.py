#!/usr/bin/env python3
"""Tests for the fresh-prime q12/orbit5867 ensemble reranker."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parents[1]
sys.path.insert(0, str(SCRIPTS))
SCRIPT = SCRIPTS / "rerank_h92_q12o5867_rootless_nagao_ensemble.py"
SPEC = importlib.util.spec_from_file_location("q12o5867_ensemble", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
ENSEMBLE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ENSEMBLE
SPEC.loader.exec_module(ENSEMBLE)


class FreshPrimeEnsembleTests(unittest.TestCase):
    def test_round_robin_blocks_are_complete_and_balanced(self) -> None:
        blocks = ENSEMBLE.round_robin_blocks((11, 13, 17, 19, 23, 29, 31), 3)
        self.assertEqual(blocks, ((11, 19, 31), (13, 23), (17, 29)))
        self.assertEqual(sorted(prime for block in blocks for prime in block), [11, 13, 17, 19, 23, 29, 31])

    def test_prime_standardization_ignores_bad_fibres(self) -> None:
        symbol = ENSEMBLE.LocalSymbol
        table = (
            symbol(11, 0, 10, 2, 1, True, False),
            symbol(11, 1, 11, 1, 3, True, False),
            symbol(11, 2, None, None, 0, False, True),
        )
        mean, standard_deviation = ENSEMBLE.prime_standardization(table)
        self.assertEqual(mean, 2.0)
        self.assertEqual(standard_deviation, 1.0)

    def test_pinned_crt_populations_reconstruct_exactly(self) -> None:
        _, box_rows = ENSEMBLE.load_box_population(ENSEMBLE.DEFAULT_BOX)
        old_pairs = {row[:2] for row in box_rows}
        observed = []
        for path in ENSEMBLE.DEFAULT_CRT:
            document, rows = ENSEMBLE.reconstruct_crt_population(path, old_pairs)
            observed.append(
                (
                    len(rows),
                    ENSEMBLE.population_sha256(rows),
                    document["construction"]["novel_parameter_population_sha256"],
                )
            )
        self.assertEqual([row[0] for row in observed], [49197, 45753])
        self.assertTrue(all(row[1] == row[2] for row in observed))


if __name__ == "__main__":
    unittest.main()
