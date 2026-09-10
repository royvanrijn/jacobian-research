"""Run with python -m unittest discover -s research/elliptic-curves/cas -p test_euclidean_seed_sieve.py."""
from fractions import Fraction as F
from pathlib import Path
import json
import random
import tempfile
import unittest
from unittest.mock import patch

import euclidean_seed_sieve as s
import run_euclidean_seed_foundry as runner


def normal_form(q=(2, 1, 1), salt=0):
    # Universal polynomial normal form (8), independent of any record fibre.
    h = s.poly((1, 1, 0, 1)); m = s.poly((1+salt, 3, 0, 0, 0, 1))
    b = s.poly((2, 0, 0, 0, 1)); k = s.poly((4, 1, 0, 1)); q = s.poly(q)
    h2 = s.mul(h, h)
    N = s.sub(s.mul(m, m), s.mul(h2, b))
    V = s.add(s.sub(s.scale(s.mul(s.mul(m, m), m), -1), s.scale(s.mul(s.mul(h2, h2), k), F(1, 2))),
              s.scale(s.mul(s.mul(h2, m), b), F(3, 2)))
    A = s.sub(s.sub(s.mul(m, k), s.scale(s.mul(b, b), F(3, 4))), s.scale(s.mul(h2, q), F(1, 4)))
    B = s.scale(s.add(s.sub(s.add(s.mul(s.mul(m, m), q), s.mul(s.mul(b, b), b)),
                           s.scale(s.mul(s.mul(m, b), k), 2)),
                     s.mul(h2, s.sub(s.mul(k, k), s.mul(b, q)))), F(1, 4))
    return (A, B, h, N, V), {'m': m, 'b': b, 'k': k, 'q': q}


class ArithmeticTests(unittest.TestCase):
    def test_trim_and_zero(self):
        self.assertEqual(s.poly([]), (0,))
        self.assertEqual(s.poly([1, 0, 0]), (1,))
        self.assertEqual(s.mul((0,), (1, 2)), (0,))

    def test_division_identities(self):
        rng = random.Random(101)
        for _ in range(100):
            a = s.poly(F(rng.randint(-20, 20), rng.randint(1, 7)) for _ in range(12))
            b = s.poly([rng.randint(-5, 5) for _ in range(4)]+[1])
            q, r = s.divrem(a, b)
            self.assertEqual(s.add(s.mul(q, b), r), a)
            self.assertTrue(r == (0,) or len(r) < len(b))

    def test_exact_div_failure(self):
        with self.assertRaises(ArithmeticError): s.exact_div((1, 0, 1), (0, 1))
        with self.assertRaises(ZeroDivisionError): s.divrem((1,), (0,))

    def test_inverse(self):
        a, mod = (1, 1), (1, 0, 1)
        self.assertEqual(s.divrem(s.mul(a, s.inverse_mod(a, mod)), mod)[1], (1,))
        with self.assertRaises(ArithmeticError): s.inverse_mod((0, 1), (0, 0, 1))

    def test_large_square_receipt(self):
        a, b = 10**600+3, 10**500+7
        value = F(a*a, b*b)
        c = s.square_certificate(value)
        self.assertEqual(F(c['root'])**2, value)
        self.assertEqual(s.square_certificate(F(a*a+1, b*b))['status'], 'NONSQUARE')
        self.assertEqual(s.square_certificate(-1)['status'], 'NONSQUARE')
        self.assertEqual(s.square_certificate(0)['root'], '0')
        self.assertEqual(s.square_certificate(F(1, 2))['status'], 'NONSQUARE')

    def test_normal_form_recovery(self):
        for salt in range(15):
            args, wanted = normal_form(salt=salt)
            result = s.conic_from_trace(*args)
            for name, coefficients in wanted.items():
                self.assertEqual(result[name], list(map(str, coefficients)))

    def test_wrong_trace_rejected(self):
        args, _ = normal_form()
        with self.assertRaises(ArithmeticError): s.conic_from_trace(*args[:-1], s.add(args[-1], (1,)))

    def test_wrong_chart_rejected(self):
        args, _ = normal_form()
        with self.assertRaises(ArithmeticError): s.conic_from_trace(args[0], args[1], (1, 1), *args[3:])

    def test_degenerate_cover_rejected(self):
        args, _ = normal_form(q=(1, 2, 1))
        with self.assertRaises(ArithmeticError): s.conic_from_trace(*args)

    def test_linear_cover(self):
        args, _ = normal_form(q=(1, 1))
        self.assertEqual(s.conic_from_trace(*args)['q'], ['1', '1'])

    def test_ramification_is_not_new_direction(self):
        args, _ = normal_form(q=(0, 1, 1))
        conic = s.conic_from_trace(*args)
        self.assertEqual(s.incidence(conic, 0)['status'], 'RAMIFIED_NO_NEW_DIRECTION')

    def test_specialized_maps_and_denominator_transport(self):
        args, _ = normal_form(q=(1, 0, 3))
        conic = s.conic_from_trace(*args)
        # 1+3*(1/2)^2 is nonsquare; at t=1 both signs map to the curve.
        self.assertEqual(s.incidence(conic, F(1, 2))['status'], 'EXACT_NONSPLIT')
        hit = s.incidence(conic, 1)
        self.assertEqual(hit['status'], 'SPLIT_REQUIRES_ADMISSION')
        x, y = map(F, hit['point'])
        self.assertEqual(y*y, x**3+s.evaluate(args[0], 1)*x+s.evaluate(args[1], 1))
        for t in (F(1, 3), F(2, 3), F(-5, 4)):
            # Build a cover whose q(t) is a square at a nonintegral address.
            args, _ = normal_form(q=(1-3*t*t, 0, 3))
            hit = s.incidence(s.conic_from_trace(*args), t)
            x, y = map(F, hit['point']); d = t.denominator
            X, Y = x*d**4, y*d**6
            self.assertEqual(Y*Y, X**3+s.evaluate(args[0], t)*d**8*X+s.evaluate(args[1], t)*d**12)

    def test_independent_sympy_oracle(self):
        try: import sympy as S
        except ImportError: self.skipTest('optional independent symbolic implementation unavailable')
        t = S.symbols('t')
        def expr(p): return sum(S.Rational(x.numerator, x.denominator)*t**i for i, x in enumerate(p))
        args, _ = normal_form(); A, B, h, N, V = map(expr, args)
        self.assertEqual(S.expand(V*V-N**3-A*N*h**4-B*h**6), 0)
        result = s.conic_from_trace(*args)
        m = S.rem(-V*S.invert(N, h*h, t), h*h, t)
        self.assertEqual(S.expand(m-expr(s.poly(result['m']))), 0)
        q = S.cancel((m**4-6*N*m*m-8*V*m-3*N*N-4*A*h**4)/h**6)
        self.assertEqual(S.expand(q-expr(s.poly(result['q']))), 0)


class ProtocolTests(unittest.TestCase):
    def test_addresses_unique_and_stable(self):
        self.assertEqual([s.address(i) for i in range(6)], ['1', '-1', '1/2', '-1/2', '2', '-2'])
        self.assertEqual(len({s.address(i) for i in range(10000)}), 10000)
        self.assertEqual(s.address(10000), s.address(10000))
        for v in (-1, 1.5, True):
            with self.assertRaises(ValueError): s.address(v)

    def test_orbit_windows_fixed_and_disjoint(self):
        rows = [{'category': 'rational', 'minimum_norm': '10', 'orbit_mask': str(i), 'parent_MW17_w': '1 -1'} for i in range(80)]
        rows += [{'category': 'elliptic', 'minimum_norm': '8', 'orbit_mask': '999', 'parent_MW17_w': '0 1'}]
        a = s.ordered_orbits(rows, 0, 30); b = s.ordered_orbits(rows[::-1], 0, 30)
        self.assertEqual(a, b)
        self.assertFalse({r['mask'] for r in a} & {r['mask'] for r in s.ordered_orbits(rows, 30, 30)})
        with self.assertRaises(ValueError): s.ordered_orbits(rows, 79, 2)
        with self.assertRaises(ValueError): s.ordered_orbits(rows+rows[:1], 0, 1)

    def test_immutable_receipts(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)/'receipt.json'; runner.save(p, {'x': [1, 2]}); runner.save(p, {'x': (1, 2)})
            with self.assertRaises(ArithmeticError): runner.save(p, {'x': [3]})

    def test_failed_worker_never_admitted(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d); runner.save(p/'result.json', {'status': 'PASS_EXACT_CONIC_MAP'})
            runner.save(p/'supervision-receipt.json', {'outcome': 'timeout', 'returncode': -9})
            self.assertFalse(runner.successful(p, 'result.json'))

    def test_success_output_tampering_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d); runner.save(p/'result.json', {'x': 1})
            runner.save(p/'supervision-receipt.json', {'outcome': 'completed', 'returncode': 0,
                'output_sha256': {'result.json': runner.sha(p/'result.json')}})
            self.assertTrue(runner.successful(p, 'result.json'))
            (p/'result.json').write_text('{"x": 2}')
            with self.assertRaises(ArithmeticError): runner.successful(p, 'result.json')

    def test_frozen_source_tampering_fails_closed(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d); root = p/'runtime'; root.mkdir(); f = root/'code.py'; f.write_text('x=1')
            runner.save(p/'plan.json', {'root': str(root)})
            runner.save(p/'manifest.json', {'plan_sha256': runner.sha(p/'plan.json'),
                'files': {'code.py': runner.sha(f)}, 'executables': {}})
            self.assertEqual(runner.guard(p)[1], root)
            f.write_text('x=2')
            with self.assertRaises(ArithmeticError): runner.guard(p)

    def test_generic_projection_drops_record_fields(self):
        zero = {'numerator': ['0'], 'denominator': ['1']}
        reduced = {'a_invariants': [zero, zero, zero,
            {'numerator': ['2', '4'], 'denominator': ['2']}, {'numerator': ['3'], 'denominator': ['1']}],
            'basis_weierstrass_coordinates': [[zero, zero]]*17, 'exceptional_points': ['NOT AN INPUT']}
        parent = runner.normalized_parent(reduced, [[int(i == j) for i in range(17)] for j in range(17)])
        self.assertEqual(parent['A_coefficients_low_to_high'], ['1', '2'])
        self.assertNotIn('exceptional_points', parent)
        self.assertEqual(parent['sections'][0]['X']['denominator_coefficients_low_to_high'], ['1'])




class ControllerIntegrationTests(unittest.TestCase):
    """Fake only the unavailable CAS; exercise real plans, files and dispatch."""
    def scenario(self, folder, *, fail_trace=False):
        args, _ = normal_form(q=(1, 0, 3))
        conic = s.conic_from_trace(*args)
        root = folder/'runtime/research'
        parent = {'A_coefficients_low_to_high': list(map(str, args[0])),
                  'B_coefficients_low_to_high': list(map(str, args[1]))}
        runner.save(root/'seed-inputs/parent.json', parent)
        plan = {'root': str(root), 'sage': '/test/sage', 'workers': 2, 'trace_seconds': 25,
                'admission_seconds': 900, 'cascade_seconds': 7200, 'rss_bytes': 100000,
                'point_calls': 64, 'height': 125000, 'orbits': [{'mask': 9, 'word': [1]*17}],
                'parameters': [{'index': 0, 'parameter': '1', 'control': True},
                               {'index': 2, 'parameter': '1/2', 'control': False}],
                'claim_boundary': 'Test only'}
        runner.save(folder/'plan.json', plan)
        runner.save(folder/'manifest.json', {'plan_sha256': runner.sha(folder/'plan.json'),
            'files': {'seed-inputs/parent.json': runner.sha(root/'seed-inputs/parent.json')}, 'executables': {}})
        calls = []
        def fake_bounded(plan, job, command, seconds):
            if (job/'supervision-receipt.json').exists(): return runner.read(job/'supervision-receipt.json')
            calls.append(command)
            r = runner.read(job/'request.json')
            if '_build' in command:
                if not fail_trace:
                    runner.save(job/'result.json', {'status': 'PASS_EXACT_CONIC_MAP', 'mask': 9,
                        'parent_sha256': r['parent_sha256'], 'request_sha256': runner.sha(job/'request.json'), 'conic': conic})
            elif '_admit' in command:
                self.assertEqual(r['incidence_sha256'], runner.sha(root/r['incidence']))
                self.assertEqual(len(r['candidates']), 1)  # Not two directions for +/- W.
                runner.save(job/'seed.json', {'rank_lower_bound': 18})
                runner.save(job/'admission.json', {'status': 'PASS_CERTIFIED_SEED'})
            else:
                self.assertEqual(r['allowance'], 64)
                self.assertEqual(r['bank_index'], 0)
                if job.name == 'seeded': self.assertEqual(r['packet_sha256'], runner.sha(root/r['packet']))
                else: self.assertNotIn('packet', r)
                runner.save(job/'result.json', {'status': 'TEST_ONLY', 'rank_lower_bound': 19})
            receipt = {'outcome': 'timeout' if fail_trace and '_build' in command else 'completed',
                       'returncode': -9 if fail_trace and '_build' in command else 0,
                       'output_sha256': {p.name: runner.sha(p) for p in job.glob('*.json') if p.name != 'request.json'}}
            runner.save(job/'supervision-receipt.json', receipt)
            return receipt
        return calls, fake_bounded

    def test_seeded_and_control_requests_and_resume(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d); calls, fake = self.scenario(folder)
            with patch.object(runner, 'bounded', fake), patch('builtins.print'):
                runner.run(folder); n = len(calls); runner.run(folder)
            self.assertEqual(len(calls), n)
            self.assertEqual(n, 4)  # Construct, admit, seeded cascade, raw control.
            rows = runner.read(folder/'STATUS.json')['slots']
            self.assertIn('seeded', rows[0]); self.assertIn('control', rows[0])
            self.assertEqual(rows[1]['square_splits'], 0)
            self.assertNotIn('seeded', rows[1])

    def test_trace_timeout_preserved_and_control_still_runs(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d); calls, fake = self.scenario(folder, fail_trace=True)
            with patch.object(runner, 'bounded', fake), patch('builtins.print'): runner.run(folder)
            result = runner.read(folder/'STATUS.json')
            self.assertEqual(result['trace_outcomes'], {'UNRESOLVED_timeout': 1})
            self.assertIn('control', result['slots'][0])
            self.assertNotIn('seeded', result['slots'][0])
            self.assertEqual(len(calls), 2)

    def test_stop_leaves_unstarted_exposures_explicit(self):
        with tempfile.TemporaryDirectory() as d:
            folder = Path(d); calls, fake = self.scenario(folder); (folder/'STOP').touch()
            with patch.object(runner, 'bounded', fake), patch('builtins.print'): runner.run(folder)
            self.assertEqual(calls, [])
            self.assertEqual(runner.read(folder/'STATUS.json')['trace_outcomes'], {'NOT_COMPLETED': 1})


if __name__ == '__main__': unittest.main()
