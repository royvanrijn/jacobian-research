#!/usr/bin/env python3
"""Unit and smoke tests for the q12/orbit5867 projective Nagao sieve."""

from __future__ import annotations

import importlib.util
import json
from math import gcd
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from time import perf_counter
import unittest
import statistics
from math import sqrt


SCRIPT = Path(__file__).with_name("search_h92_q12o5867_rootless_nagao.py")
CPP = Path(__file__).with_name("scan_h92_q12o5867_rootless_nagao.cpp")
SPEC = importlib.util.spec_from_file_location("q12o5867_nagao", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SIEVE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = SIEVE
SPEC.loader.exec_module(SIEVE)


class ProjectiveNagaoSieveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = SIEVE.load_family_model()

    def test_complete_projective_table_and_trace_identity(self) -> None:
        self.assertEqual(self.model.coordinate, "elkies_2026_published_t")
        table = SIEVE.residue_table(self.model, 19)
        self.assertEqual(len(table), 20)
        self.assertEqual(table[-1].label, "infinity")
        self.assertTrue(any(symbol.good_reduction for symbol in table))
        for symbol in table:
            self.assertEqual(symbol.good_reduction, not symbol.singular_mod_prime)
            if symbol.good_reduction:
                self.assertEqual(symbol.point_count, 20 - symbol.trace)

    def test_raw_q12_chart_remains_available_as_regression_input(self) -> None:
        raw_model = SIEVE.load_family_model(
            SIEVE.ROOT / "artifacts/local/elkies-k3/q12o5867-smooth-rr-qq.json"
        )
        self.assertEqual(raw_model.coordinate, "q12o5867_raw_u")
        self.assertEqual((raw_model.a_degree, raw_model.b_degree), (8, 12))

    def test_projective_reduction_including_infinity(self) -> None:
        self.assertEqual(SIEVE.projective_index(3, 2, 7), 5)
        self.assertEqual(SIEVE.projective_index(1, 7, 7), 7)
        with self.assertRaises(ValueError):
            SIEVE.projective_index(7, 7, 7)

    def test_primitive_enumeration_is_canonical(self) -> None:
        pairs = [(row.numerator, row.denominator) for row in SIEVE.primitive_parameters(2, 2)]
        self.assertEqual(pairs[0], (1, 0))
        self.assertIn((0, 1), pairs)
        self.assertIn((-1, 2), pairs)
        self.assertIn((1, 2), pairs)
        self.assertNotIn((2, 2), pairs)
        self.assertEqual(len(pairs), len(set(pairs)))

    def test_actual_family_staged_smoke_and_benchmark(self) -> None:
        blocks, rejected = SIEVE.build_residue_tables(self.model, ((19, 41), (43, 61)))
        self.assertEqual(rejected, ())
        started = perf_counter()
        survivors, stages = SIEVE.run_staged_sieve(
            numerator_bound=40,
            denominator_bound=30,
            table_blocks=blocks,
            keep_per_bucket=(8, 4),
            bucket_width=10,
        )
        elapsed = perf_counter() - started
        self.assertTrue(survivors)
        self.assertEqual(len(stages), 2)
        self.assertGreater(stages[0]["population_scored"], stages[0]["retained_count"])
        self.assertLessEqual(stages[1]["retained_count"], 4 * stages[1]["height_bucket_count"])
        self.assertTrue(all(len(candidate.block_score_units) == 2 for candidate in survivors))
        print(
            "SMOKE_BENCHMARK "
            f"population={stages[0]['population_scored']} "
            f"stage1_rate={stages[0]['parameters_per_second']:.0f}/s "
            f"elapsed={elapsed:.4f}s"
        )

    def test_cpp_hot_loop_matches_python_survivor_order_and_scores(self) -> None:
        compiler = shutil.which("g++")
        if compiler is None:
            self.skipTest("g++ is unavailable")
        blocks, rejected = SIEVE.build_residue_tables(
            self.model, SIEVE.DEFAULT_PRIME_BLOCKS
        )
        self.assertEqual(rejected, ())
        python_survivors, _ = SIEVE.run_staged_sieve(
            numerator_bound=40,
            denominator_bound=30,
            table_blocks=blocks,
            keep_per_bucket=(8, 4, 2),
            bucket_width=10,
        )
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            table_path = temporary / "tables.txt"
            binary_path = temporary / "scan"
            output_path = temporary / "result.json"
            SIEVE.export_cpp_tables(table_path, self.model, blocks)
            subprocess.run(
                [compiler, "-O3", "-std=c++17", str(CPP), "-o", str(binary_path)],
                check=True,
            )
            subprocess.run(
                [
                    str(binary_path),
                    str(table_path),
                    "40",
                    "30",
                    "10",
                    "8,4,2",
                    "100",
                    str(output_path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            cpp = json.loads(output_path.read_text())
        expected = [
            (
                [candidate.numerator, candidate.denominator],
                list(candidate.block_score_units),
                candidate.total_score_units,
            )
            for candidate in python_survivors
        ]
        observed = [
            (
                record["projective_pair"],
                record["block_score_units_1e12"],
                record["total_score_units_1e12"],
            )
            for record in cpp["finalists"]
        ]
        self.assertEqual(observed, expected)

    def test_cpp_skew_chart_scores_the_actual_projective_parameter(self) -> None:
        compiler = shutil.which("g++")
        if compiler is None:
            self.skipTest("g++ is unavailable")
        blocks, rejected = SIEVE.build_residue_tables(
            self.model, SIEVE.DEFAULT_PRIME_BLOCKS
        )
        self.assertEqual(rejected, ())
        scale = 17
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            table_path = temporary / "tables.txt"
            binary_path = temporary / "scan"
            output_path = temporary / "result.json"
            SIEVE.export_cpp_tables(table_path, self.model, blocks)
            subprocess.run(
                [compiler, "-O3", "-std=c++17", str(CPP), "-o", str(binary_path)],
                check=True,
            )
            subprocess.run(
                [
                    str(binary_path),
                    str(table_path),
                    "20",
                    "15",
                    "5",
                    "8,4,2",
                    "100",
                    str(output_path),
                    str(scale),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            result = json.loads(output_path.read_text())
        self.assertEqual(result["search"]["parameter_scale"], scale)
        self.assertEqual(
            result["schema"], "h92-q12o5867-rootless-projective-nagao-cpp-skew-v2"
        )
        inverse_cache = {}
        for record in result["finalists"]:
            chart_numerator, chart_denominator = record["chart_projective_pair"]
            common = gcd(abs(scale * chart_numerator), chart_denominator)
            expected_pair = [
                scale * chart_numerator // common,
                chart_denominator // common,
            ]
            self.assertEqual(record["projective_pair"], expected_pair)
            candidate = SIEVE.Candidate(
                numerator=expected_pair[0],
                denominator=expected_pair[1],
                height=max(abs(expected_pair[0]), expected_pair[1]),
            )
            for block in blocks:
                candidate = SIEVE.score_block(candidate, block, inverse_cache)
            self.assertEqual(
                record["block_score_units_1e12"],
                list(candidate.block_score_units),
            )

    def test_cpp_complete_worst_block_control_ranks_match_python(self) -> None:
        compiler = shutil.which("g++")
        if compiler is None:
            self.skipTest("g++ is unavailable")
        requested = ((19, 31, 43), (23, 37, 47), (29, 41, 53))
        blocks, rejected = SIEVE.build_residue_tables(self.model, requested)
        self.assertEqual(rejected, ())
        controls = ((0, 1), (1, 1), (-1, 1))

        standards = {}
        for block in blocks:
            for prime, table in block.items():
                values = [
                    symbol.contribution_units
                    for symbol in table
                    if symbol.good_reduction
                ]
                standards[prime] = (
                    statistics.fmean(values),
                    statistics.pstdev(values),
                )

        def python_key(pair):
            numerator, denominator = pair
            signals = []
            good = bad = 0
            for block in blocks:
                total = 0.0
                for prime, table in block.items():
                    index = SIEVE.projective_index(numerator, denominator, prime)
                    symbol = table[index]
                    if symbol.good_reduction:
                        mean, deviation = standards[prime]
                        total += (symbol.contribution_units - mean) / deviation
                        good += 1
                    else:
                        bad += 1
                signals.append(total / sqrt(len(block)))
            height = 1 if denominator == 0 else max(abs(numerator), denominator)
            return (
                min(signals),
                statistics.fmean(signals),
                good,
                -bad,
                -height,
                -denominator,
                -numerator,
            )

        population = [
            (candidate.numerator, candidate.denominator)
            for candidate in SIEVE.primitive_parameters(20, 20)
        ]
        expected_ranks = [
            1 + sum(python_key(pair) > python_key(control) for pair in population)
            for control in controls
        ]

        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            table_path = temporary / "tables.txt"
            binary_path = temporary / "scan"
            output_path = temporary / "result.json"
            SIEVE.export_cpp_tables(table_path, self.model, blocks)
            subprocess.run(
                [compiler, "-O3", "-std=c++17", str(CPP), "-o", str(binary_path)],
                check=True,
            )
            completed = subprocess.run(
                [
                    str(binary_path),
                    str(table_path),
                    "20",
                    "20",
                    "1",
                    "1,1,1",
                    "10",
                    str(output_path),
                    "1",
                    ",".join(f"{a}/{b}" for a, b in controls),
                ],
                capture_output=True,
                text=True,
            )
            self.assertIn(completed.returncode, (0, 3))
            result = json.loads(output_path.read_text())
        self.assertEqual(result["population_count"], len(population))
        self.assertEqual(
            [record["population_rank"] for record in result["positive_controls"]],
            expected_ranks,
        )
        self.assertEqual(result["scoring"]["primary_ranking_key"], "minimum block signal")

    def test_cpp_rank_region_scores_complete_disjoint_shell(self) -> None:
        compiler = shutil.which("g++")
        if compiler is None:
            self.skipTest("g++ is unavailable")
        requested = ((19, 31, 43), (23, 37, 47), (29, 41, 53))
        blocks, rejected = SIEVE.build_residue_tables(self.model, requested)
        self.assertEqual(rejected, ())
        standards = {}
        for block in blocks:
            for prime, table in block.items():
                values = [
                    symbol.contribution_units
                    for symbol in table
                    if symbol.good_reduction
                ]
                standards[prime] = (
                    statistics.fmean(values),
                    statistics.pstdev(values),
                )

        def score(pair):
            numerator, denominator = pair
            signals = []
            good = bad = 0
            for block in blocks:
                total = 0.0
                for prime, table in block.items():
                    symbol = table[SIEVE.projective_index(numerator, denominator, prime)]
                    if symbol.good_reduction:
                        mean, deviation = standards[prime]
                        total += (symbol.contribution_units - mean) / deviation
                        good += 1
                    else:
                        bad += 1
                signals.append(total / sqrt(len(block)))
            height = max(abs(numerator), denominator)
            frozen_key = (
                min(signals),
                statistics.fmean(signals),
                good,
                -bad,
                -height,
                -denominator,
                -numerator,
            )
            ordinary_key = (statistics.fmean(signals), *frozen_key)
            return frozen_key, ordinary_key

        population = [
            (candidate.numerator, candidate.denominator)
            for candidate in SIEVE.primitive_parameters(30, 30)
            if candidate.denominator and candidate.height >= 11
        ]
        expected_frozen = sorted(population, key=lambda pair: score(pair)[0], reverse=True)[:20]
        expected_ordinary = sorted(population, key=lambda pair: score(pair)[1], reverse=True)[:5]
        with tempfile.TemporaryDirectory() as temporary_directory:
            temporary = Path(temporary_directory)
            table_path = temporary / "tables.txt"
            binary_path = temporary / "scan"
            output_path = temporary / "result.json"
            SIEVE.export_cpp_tables(table_path, self.model, blocks)
            subprocess.run(
                [compiler, "-O3", "-std=c++17", str(CPP), "-o", str(binary_path)],
                check=True,
            )
            subprocess.run(
                [
                    str(binary_path), str(table_path), "30", "30", "1", "1,1,1",
                    "20", str(output_path), "1", "--rank-region", "11", "5",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            result = json.loads(output_path.read_text())
        self.assertEqual(result["population_count"], len(population))
        self.assertEqual(
            [tuple(row["score"]["projective_pair"]) for row in result["ranked_prefix"]],
            expected_frozen,
        )
        self.assertEqual(
            [
                tuple(row["score"]["projective_pair"])
                for row in result["ordinary_nagao_control_prefix"]
            ],
            expected_ordinary,
        )
        self.assertEqual(len(result["random_control_lane"]), 5)


if __name__ == "__main__":
    unittest.main()
