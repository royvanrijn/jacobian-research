import importlib.util
import os
import itertools
import math
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

CAS = Path(__file__).resolve().parents[1] / "cas"
sys.path.insert(0, str(CAS))

from curve302_short_vector_core import (
    enumerate_primitive_directions, exact_ldl, in_rational_span,
    intersection_saturated, lattice_contains, one_step_core_hit, primitive,
    qnorm, rank_interval, rational_rank, saturation_basis, saturation_index,
)


class ShortVectorCoreTests(unittest.TestCase):
    def test_primitive_sign_and_gcd(self):
        self.assertEqual(primitive((-4, 2, 0)), (2, -1, 0))
        self.assertEqual(primitive((0, -3, 6)), (0, 1, -2))
        self.assertEqual(primitive((0, 0)), (0, 0))

    def test_ldl_reconstructs(self):
        q = ((F(4), F(1), F(1)), (F(1), F(3), F(0)), (F(1), F(0), F(2)))
        L, D = exact_ldl(q)
        for i in range(3):
            for j in range(3):
                value = sum(L[i][k] * D[k] * L[j][k] for k in range(3))
                self.assertEqual(value, q[i][j])

    def test_exact_enumeration_matches_bruteforce(self):
        q = ((F(3), F(1)), (F(1), F(2)))
        bound = F(12)
        rows, _ = enumerate_primitive_directions(q, bound, max_directions=1000, max_nodes=10000)
        got = {(norm, vector) for norm, vector in rows}
        brute = set()
        for x, y in itertools.product(range(-5, 6), repeat=2):
            if not (x or y) or math.gcd(abs(x), abs(y)) != 1:
                continue
            first = x if x else y
            if first < 0:
                continue
            v = (x, y); norm = qnorm(q, v)
            if norm <= bound:
                brute.add((norm, v))
        self.assertEqual(got, brute)

    def test_enumeration_boundary_vector(self):
        q = ((F(2), F(1)), (F(1), F(2)))
        rows, _ = enumerate_primitive_directions(q, F(6), max_directions=100, max_nodes=1000)
        self.assertIn((F(6), (1, 1)), rows)

    def test_resource_cap_is_fail_closed(self):
        q = ((F(1), F(0)), (F(0), F(1)))
        with self.assertRaisesRegex(ValueError, "direction cap"):
            enumerate_primitive_directions(q, F(20), max_directions=1, max_nodes=10000)

    def test_saturation_index_and_basis(self):
        rows = ((2, 0, 0), (0, 2, 0))
        self.assertEqual(saturation_index(rows, 3), 4)
        sat = saturation_basis(rows, 3)
        self.assertEqual(sat, ((1, 0, 0), (0, 1, 0)))
        self.assertEqual(saturation_index(sat, 3), 1)

    def test_saturation_skew_line(self):
        self.assertEqual(saturation_basis(((6, 9, 3),), 3), ((2, 3, 1),))

    def test_intersection_saturated(self):
        a = ((1, 0, 0), (0, 1, 0))
        b = ((1, 0, 0), (0, 0, 1))
        self.assertEqual(intersection_saturated((a, b), 3), ((1, 0, 0),))

    def test_lattice_containment(self):
        plane = ((1, 0, 0), (0, 1, 0))
        line = ((1, 1, 0),)
        other = ((0, 0, 1),)
        self.assertTrue(lattice_contains(plane, line, 3))
        self.assertFalse(lattice_contains(plane, other, 3))

    def test_rank_interval_exact_tie(self):
        rows = [(F(1), (1, 0)), (F(1), (0, 1)), (F(2), (1, -1))]
        rows.sort(key=lambda row: (row[0], row[1]))
        lo, hi, norm = rank_interval(rows, (0, 1))
        self.assertEqual((lo, hi, norm), (1, 2, F(1)))

    def test_one_step_core_hit(self):
        prefix = ((1, 0, 0),)
        core = ((1, 0, 0), (0, 1, 0))
        self.assertTrue(one_step_core_hit(prefix, (0, 1, 1), core, 3) is False)
        self.assertTrue(one_step_core_hit(prefix, (0, 1, 0), core, 3))

    def test_rational_span(self):
        self.assertTrue(in_rational_span((2, 2), ((1, 1),), 2))
        self.assertFalse(in_rational_span((1, 0), ((1, 1),), 2))


class ControllerIntegrationTests(unittest.TestCase):
    def test_full_synthetic_controller_and_check(self):
        import hashlib, json, subprocess, tempfile
        from pathlib import Path
        runner = CAS / "run_curve302_short_vector_core.py"
        names = tuple([f"recovered-local-{i:02d}" for i in range(1, 5)] +
                      [f"recovered-strict-{i:02d}" for i in range(1, 4)] +
                      [f"residual-strict-{i:02d}" for i in range(1, 8)])
        with tempfile.TemporaryDirectory(prefix="short-core-test-") as tmp:
            tmp = Path(tmp); source = tmp / "source"; source.mkdir(); output = tmp / "out"
            rel = {"status":"PASS_QUOTIENT_RELATION_ANALYSIS", "direction_ids":list(names),
                   "schur_quotient":[["1" if i == j else "0" for j in range(14)] for i in range(14)]}
            runs = []
            for seed in range(14):
                order = [j for j in range(14) if j != seed]
                if seed == 0: order = order[:11]
                stages = []
                for epoch, j in enumerate(order):
                    word = [0]*14; word[j] = 1
                    stages.append({"epoch":epoch, "new":[{"integral":True, "denominator":1,
                                  "quotient_word":word, "primitive_quotient_word":word}]})
                runs.append({"seed":names[seed], "final_rank":18+len(order), "stages":stages})
            traj = {"status":"PASS_ALL_14_SEEDED_TRAJECTORIES_RECONCILED", "direction_ids":list(names), "runs":runs}
            (source/"quotient-relations.json").write_text(json.dumps(rel, sort_keys=True)+"\n")
            (source/"trajectories.json").write_text(json.dumps(traj, sort_keys=True)+"\n")
            def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()
            report = {"status":"PASS_THREE_CLOSURE_EXPERIMENTS", "outputs":{
                "quotient-relations":digest(source/"quotient-relations.json"),
                "trajectories":digest(source/"trajectories.json")}}
            (source/"REPORT.json").write_text(json.dumps(report, sort_keys=True)+"\n")
            env = {**os.environ, "PYTHONPATH":str(CAS)}
            subprocess.run([sys.executable, str(runner), "run", "--source", str(source), "--folder", str(output),
                            "--max-directions", "1000", "--max-nodes", "100000", "--stage-seconds", "60"],
                           check=True, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            subprocess.run([sys.executable, str(runner), "check", "--folder", str(output)],
                           check=True, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            result = json.loads((output/"REPORT.json").read_text())
            self.assertEqual(result["status"], "PASS_THREE_SHORT_VECTOR_CORE_EXPERIMENTS")
            enum = json.loads((output/"enumeration.json").read_text())
            self.assertEqual(enum["direction_count"], 14)
            ranks = json.loads((output/"ranks-basins.json").read_text())
            self.assertEqual(ranks["rank_summary"]["coverage"], 180)

if __name__ == "__main__":
    unittest.main()
