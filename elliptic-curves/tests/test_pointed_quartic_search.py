from fractions import Fraction as Q
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import ast
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"elliptic-curves/cas"))
import pointed_quartic_search as pq
import run_pointed_quartic_search as jobs
import mw16_sensitivity_backend as regression
import pointed_quartic_migration as migration


class PointedQuarticSearchTests(unittest.TestCase):
    def test_same_chart_after_a_known_point_observation_has_distinct_checkpoint(self):
        # This is the state change seen when a timed-out chart is retried:
        # an additional known-point observation changes no search coordinates.
        from research_runtime.finite_reduction import ReductionCache
        from research_runtime.store import FactStore
        from research_runtime.search_state import raw_state
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);cache=ReductionCache(FactStore(root/'facts'))
            state=raw_state((0,0,0,-7,10),[(1,2),(2,2)],cache=cache,prime_bound=31)
            observed=state.adjoin(state.basis[0],cache=cache)
            self.assertEqual(state.basis,observed.basis)
            self.assertNotEqual(state.key,observed.key)
            first=pq.PointedQuarticSearch(state=state,centre={'coefficients':[1,0]})
            second=pq.PointedQuarticSearch(state=observed,centre={'coefficients':[1,0]})
            self.assertEqual(first.input_record(),second.input_record())
            a=first.search(19,2,checkpoint_dir=root/'charts')
            b=second.search(19,2,checkpoint_dir=root/'charts')
            self.assertEqual(a.curve_points,b.curve_points)
            self.assertEqual(len(list((root/'charts').glob('*.json'))),2)
            with patch.object(pq,'search_box',side_effect=AssertionError('unexpected re-enumeration')):
                self.assertEqual(second.search(19,2,checkpoint_dir=root/'charts'),b)

    model = (0, 0, 0, -1, 1)
    p = (Q(0), Q(1))

    def search(self, policy=None):
        return pq.PointedQuarticSearch(self.model, [self.p], {"coefficients": [1]}, policy)

    def test_arbitrary_subgroup_sizes_and_explicit_centres(self):
        expected = self.search().search(19, 2).curve_points
        # These deliberately dependent generators test the API's dimension
        # independence; no fake rank claim is attached to their cardinality.
        for size in (1, 16, 17, 18, 32):
            s = pq.PointedQuarticSearch(self.model, [self.p]*size,
                {"coefficients": [1]+[0]*(size-1), "point": self.p})
            self.assertEqual(s.search(19, 2).curve_points, expected)
        s = pq.PointedQuarticSearch(self.model, [], {"point": self.p})
        self.assertEqual(s.search(19, 2).curve_points, expected)

    def test_general_weierstrass_transport(self):
        # E: y^2+xy+y=x^3-x has (0,0). Complete square and cube over Q,
        # then compare every recovered point with the same short-model box.
        curve = (1, 0, 1, -1, 0)
        s = pq.PointedQuarticSearch(curve, [(0, 0)], {"coefficients": [1]}, "raw")
        short = pq.PointedQuarticSearch(s.short_model, [s.to_short((0, 0))], {"coefficients": [1]}, "raw")
        found = s.search(31, 2)
        self.assertTrue(found.curve_points)
        self.assertTrue(all(s.on_curve(p) for p in found.curve_points))
        self.assertEqual({s.to_short(p) for p in found.curve_points}, set(short.search(31, 2).curve_points))
        self.assertEqual(s.verify_record(found.record), found)

    def test_coordinates_match_preserved_sensitivity_control(self):
        p = pq.base.linear_combination(self.model, [self.p], [7])
        for spec in ("gauss", "raw", "metric:1/16", "metric:16", "gauss:1,2,1", "gauss:2,1,-1"):
            with self.subTest(spec=spec):
                new = pq.PointedQuarticSearch(self.model, [p], {"coefficients": [1]}, spec).search(35, 2)
                old = regression.run(model=self.model, points=[p], representative=[1], mask=1,
                    specification=spec, height=35, seconds=2)
                self.assertEqual(new.record["finite_curve_points"], old["finite_curve_points"])

    def test_general_pgl2_and_projective_poles_against_brute_force(self):
        for matrix in ((0, 1, 1, 0), (Q(2, 3), -1, Q(5, 7), 2), (2, 2, 0, 2)):
            s = self.search({"kind": "metric", "weight": "16", "matrix": matrix})
            r = s.search(17, 2)
            expected = set()
            for d in range(1, 18):
                for n in range(-17, 18):
                    v = sum(c*n**i*d**(4-i) for i, c in enumerate(s.coefficients))
                    if gcd(n, d) == 1 and v >= 0 and isqrt(v)**2 == v:
                        expected.add((n, d, isqrt(v)))
            got = {tuple(map(int, p)) for p in r.record["primitive_square_hits"] if int(p[1])}
            self.assertEqual(got, expected)
            self.assertTrue(r.record["infinity_checked"])
            self.assertEqual(s.verify_record(r.record), r)
        self.assertEqual(pq.CoordinatePolicy(matrix=(2, 2, 0, 2)), pq.CoordinatePolicy(matrix=(1, 1, 0, 1)))

    def test_checkpoint_is_bound_to_sources_inputs_policy_and_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            s = self.search()
            original = s.search(19, 2, checkpoint_dir=directory)
            with patch.object(pq, "search_box", side_effect=RuntimeError("fresh search")):
                self.assertEqual(s.search(19, 2, checkpoint_dir=directory), original)
                for changed, height in ((self.search("gauss"), 19), (s, 20)):
                    with self.assertRaises(RuntimeError):
                        changed.search(height, 2, checkpoint_dir=directory)
            path = next(Path(directory).glob("*.json"))
            saved = json.loads(path.read_text())
            # Recompute the outer digest too: exact chart/map validation must
            # catch tampering that a checksum alone would miss.
            saved["record"]["finite_curve_points"] = []
            saved["record_sha256"] = sha256(pq.canonical(saved["record"]).encode()).hexdigest()
            path.write_text(json.dumps(saved))
            with self.assertRaises(ArithmeticError):
                s.search(19, 2, checkpoint_dir=directory)

    def test_timeout_cannot_be_resumed_as_complete_and_shards_match(self):
        s = self.search()
        with tempfile.TemporaryDirectory() as directory:
            timeout = s.search(10000, 1e-9, checkpoint_dir=directory)
            self.assertEqual(timeout.record["status"], "bounded_search_timeout")
            self.assertEqual(timeout.record["integer_pairs_covered"], 0)
            with patch.object(pq, "search_box", side_effect=RuntimeError("retry")):
                with self.assertRaises(RuntimeError):
                    s.search(10000, 1e-9, checkpoint_dir=directory)
        one = s.search(29, 2, denominator_end=13)
        two = s.search(29, 2, denominator_start=14)
        self.assertEqual(set(one.curve_points)|set(two.curve_points), set(s.search(29, 2).curve_points))
        self.assertEqual(one.record["integer_pairs_covered"]+two.record["integer_pairs_covered"], 29*59)

    def test_invalid_inputs_fail_before_search(self):
        bad = [((0, 0), [self.p], {"coefficients": [1]}, None),
               (self.model, [(0, 0)], {"coefficients": [1]}, None),
               (self.model, [self.p], {"coefficients": [0]}, None),
               (self.model, [self.p], {"coefficients": [Q(1, 2)]}, None),
               (self.model, [self.p], {"coefficients": [1], "point": (1, 1)}, None),
               (self.model, [self.p], {"coefficients": [1]}, "red"),
               (self.model, [self.p], {"coefficients": [1]}, {"matrix": (1, 1, 1, 1)})]
        for args in bad:
            with self.assertRaises(ValueError):
                pq.PointedQuarticSearch(*args)
        for seconds in (0, float("inf"), float("nan")):
            with self.assertRaises(ValueError):
                self.search().search(19, seconds)

    def test_actual_mw18_specialization_uses_same_service(self):
        document = json.loads(jobs.MW18.read_text())
        row = document["candidates"][0]
        manifest = jobs.mw18_jobs(document, row["candidate_id"],
            [{"coefficients": [1]+[0]*17}], "metric:16", 5, 2)
        with tempfile.TemporaryDirectory() as directory:
            result = jobs.execute(manifest, Path(directory)/"result.json", Path(directory)/"charts")
        self.assertEqual(result["status"], "COMPLETE")
        self.assertEqual(result["results"][0]["search"]["backend"], pq.BACKEND_NAME)
        self.assertEqual(len(result["results"][0]["search"]["input"]["subgroup"]), 18)

    def test_large_prospective_quartic_canary(self):
        source = ROOT/"artifacts/generated-results/elliptic-curves/icarm_mw16_nagao_finalist_specializations_h300_v1.json"
        row = json.loads(source.read_text())["candidates"][0]
        search = pq.PointedQuarticSearch(row["raw_short_model"], row["raw_generic_points"],
                                        {"coefficients": [1]+[0]*15})
        self.assertGreater(search.chart_record()["maximum_coefficient_bits"], 3000)
        result = search.search(17, 2)
        self.assertEqual(result.record["status"], "bounded_search_complete")
        self.assertEqual(result.record["integer_pairs_covered"], 17*35)

    def test_active_raw_wrappers_have_no_pari_search(self):
        for filename, function in (("run_mw17_jump_v2.sage", "run_quartic_search_raw"),
                                   ("run_a1_mw16_target_free_parameter_search.sage", "direct_integral_quartic_search")):
            source = (ROOT/"elliptic-curves/cas"/filename).read_text()
            node = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == function)
            segment = ast.get_source_segment(source, node)
            self.assertIn("from pointed_quartic_search import run_quartic_search", segment)
            self.assertNotIn("hyperellratpoints", segment)

    def test_frozen_sources_and_runtime_checkpoints_are_separate(self):
        campaign = json.loads((ROOT/"artifacts/generated-results/elkies-k3-mw17-jump-v2-campaign-v1.json").read_text())
        migration.validate_frozen_sources(campaign["implementation_hashes"])
        migration.require_runtime({"runtime_search": migration.runtime_search()})
        with self.assertRaises(ArithmeticError):
            migration.require_runtime({"campaign_sha256": "historical"})
        with self.assertRaises(ArithmeticError):
            migration.validate_frozen_sources({"elliptic-curves/cas/half_lattice_pointed_sieve.py": "wrong"})

    def test_sensitivity_checker_accepts_universal_rational_coordinates(self):
        import verify_mw16_sensitivity as exact
        with tempfile.TemporaryDirectory() as directory:
            record = pq.checkpoint(directory, model=self.model, points=[self.p],
                representative=[1], mask=1, specification="gauss:2,1,-1", height=19, seconds=2)
            self.assertEqual(exact.check_chart(self.model, [self.p], record), set(map(pq.point, record["finite_curve_points"])))


if __name__ == "__main__":
    unittest.main()
