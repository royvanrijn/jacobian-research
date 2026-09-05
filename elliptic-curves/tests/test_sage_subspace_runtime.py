"""Tiny exact controls for the concrete shared subspace backend."""
import copy
import gzip
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"cas"))
from research_runtime.arithmetic import TwoTorsionContext
from research_runtime.finite_reduction import ReductionCache
from research_runtime.mw_state import MWState
from research_runtime.sage_arithmetic import SageArithmetic
from research_runtime.sage_subspace import SageSubspaceBackend, PointClassWitness
from research_runtime.store import FactStore
from research_runtime.subspace import SubspaceDescent


class SageSubspaceTests(unittest.TestCase):
    def setUp(self):
        self.directory = TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.arithmetic = SageArithmetic(FactStore(self.directory.name))
        algebra = TwoTorsionContext((10, -7, 0, 1))
        self.context = self.arithmetic.prepare_congruent((0, 0, 0, -7, 10), algebra,
            (0, 1, 0), factor_primes=[2, 83], discover=True)
        self.arithmetic.field(algebra, factor_primes=[2, 83], discover=True)
        self.cache = ReductionCache(self.arithmetic.store)
        self.state = MWState.empty(self.context, cache=self.cache,
            primes=[5, 7, 11, 13, 17, 19, 23, 29], no_two_torsion_prime=3)
        for point in [(1, 2), (2, 2)]:
            self.state = self.state.adjoin(point, cache=self.cache)
        self.global_witness = PointClassWitness(self.state, self.cache)
        self.backend = SageSubspaceBackend(self.arithmetic, self.context, self.global_witness)
        self.descent = SubspaceDescent(self.context, self.global_witness.classes, self.backend)

    def test_local_features_reduce_rational_coefficients_without_global_setup(self):
        algebra=TwoTorsionContext(('1/25','-1/5','0','1'))
        with patch.object(self.arithmetic,'field',side_effect=AssertionError('order requested')),patch.object(self.arithmetic,'bnf',side_effect=AssertionError('BNF requested')):
            record=self.arithmetic.fast_features(algebra,primes=[2,3,5])
        self.assertEqual(record['local_features'][-1]['status'],'NONINTEGRAL_PRESENTATION')
        self.assertTrue(all(r['completion_factor_degrees'] is not None for r in record['local_features'][:2]))

    def test_raw_rational_curve_upgrades_to_an_integral_labelled_field(self):
        from fractions import Fraction
        from research_runtime.search_state import raw_state
        state=raw_state([0,0,0,'-7/625','10/15625'],[['1/25','2/125']],cache=self.cache,prime_bound=43)
        context=self.arithmetic.prepare_context(state.arithmetic,factor_primes=[2,5,83],discover=True)
        self.assertTrue(all(Fraction(c).denominator==1 for c in context.two_torsion.polynomial))
        self.arithmetic.field(context.two_torsion,factor_primes=[2,5,83],discover=True)
        self.assertEqual(state.with_arithmetic(context).reductions,state.reductions)
        self.assertNotEqual(context.two_torsion.key,state.arithmetic.two_torsion.key)

    def test_state_roundtrip_and_wrong_signatures(self):
        self.assertEqual(MWState.from_record(self.state.record(), cache=self.cache), self.state)
        record = self.state.record()
        record["state"]["reductions"]["columns"][0] ^= 1
        with self.assertRaises((ArithmeticError, ValueError)):
            MWState.from_record(record, cache=self.cache)

    def test_two_point_classes_covers_ct_and_readonly_replay(self):
        with patch.object(self.arithmetic, "bnf", side_effect=AssertionError("BNF requested")):
            witness = self.descent.run()
        self.assertEqual(witness["admissible_masks"], [1, 2])
        self.assertEqual(witness["radical"]["radical_dimension"], 2)
        self.assertEqual(witness["ct"]["matrix"], [[0, 0], [0, 0]])
        with patch.object(self.backend, "local_map", side_effect=AssertionError("local search")), \
             patch.object(self.backend, "cover", side_effect=AssertionError("conic solving")), \
             patch.object(self.backend, "ct_pairing", side_effect=AssertionError("pairing discovery")), \
             patch.object(self.arithmetic, "factor_integer", side_effect=AssertionError("refactoring")), \
             patch.object(self.arithmetic, "square_root", side_effect=AssertionError("root search")):
            self.assertEqual(self.descent.run(retained=witness), witness)
        for location in ("local", "cover", "ct"):
            bad = copy.deepcopy(witness)
            if location == "local":
                bad["local_maps"][0]["quotient_rows"][0].append(1)
            elif location == "cover":
                bad["covers"][0]["d_over_quartic_y"] = "0"
            else:
                bad["ct"]["pairs"][0]["local_terms"][0]["hilbert_symbol"] *= -1
            with self.assertRaises((ArithmeticError, ValueError)):
                self.descent.run(retained=bad)

    def test_published_nonzero_pairing_with_generic_replayer(self):
        from sage.all import QQ
        root = Path(__file__).resolve().parents[2]
        evidence = root/"artifacts/generated-results/elliptic-curves/fixed_cubic_u_minus1_cassels_tate_evidence_v1.json.gz"
        with gzip.open(evidence, "rt") as handle:
            control = json.load(handle)["published_control"]
        self.backend.I, self.backend.J = QQ(control["I"]), QQ(control["J"])
        quartics = [self.backend.R(list(map(QQ, q))) for q in control["quartics"]]
        places, factors = self.backend._pair_support(quartics)
        self.assertEqual(places, [2, 3, 5, 7, 571, "infinity"])
        record = {"masks": [1, 2, 3], "gamma": control["gamma"],
                  "square_root_phi_coefficients": control["square_root_phi_coefficients"],
                  "support_factors": factors, "local_terms": control["local_terms"], "value": 1}
        result = self.backend._pair(None, (1, 2, 3), [{"quartic": q} for q in control["quartics"]], retained=record)
        self.assertEqual(result["value"], 1)

    def test_portable_command_replay_without_discovery(self):
        from run_arithmetic_pipeline import worker
        from research_runtime.supervisor import Limits
        request = json.loads((Path(__file__).resolve().parents[1]/"data/runtime_subspace_control.json").read_text())
        limits = Limits(30, 1_073_741_824)
        witness = worker(request, self.arithmetic.store, limits)
        with TemporaryDirectory() as directory, \
             patch("research_runtime.sage_subspace.local_candidates", side_effect=AssertionError("local search")), \
             patch("research_runtime.sage_subspace.quartic_local_witness", side_effect=AssertionError("CT search")), \
             patch.object(SageSubspaceBackend, "cover", side_effect=AssertionError("conic solving")), \
             patch.object(SageSubspaceBackend, "search_cover", side_effect=AssertionError("point search")):
            replay = worker(request, FactStore(directory), limits, retained=witness)
        self.assertEqual(replay["descent"], witness["descent"])
        self.assertEqual(replay["final_search_state"], witness["final_search_state"])

    def test_general_weierstrass_state(self):
        context = self.arithmetic.prepare((0, 0, 1, -1, 0), factor_primes=[37], discover=True)
        state = MWState.empty(context, cache=self.cache, primes=[5, 7, 11, 13, 17, 19], no_two_torsion_prime=3)
        state = state.adjoin((0, 0), cache=self.cache)
        self.assertEqual(state.rank, 1)
        self.assertTrue(state.verify(self.cache))

    def test_complete_selmer_cold_command_and_retained_integrity(self):
        from run_arithmetic_pipeline import worker
        from research_runtime.supervisor import Limits
        request = {"mode":"complete-selmer", "requirement":"unconditional-upper-bound",
                   "curve":{"model":[0,0,1,-1,0], "factor_primes":[37]}}
        with TemporaryDirectory() as directory:
            witness = worker(request,FactStore(directory),Limits(10,600_000_000))
        self.assertEqual(witness["selmer"]["full_selmer_dimension"],1)
        with TemporaryDirectory() as directory, patch.object(SageArithmetic,"bnf",side_effect=AssertionError("BNF replay")):
            replay = worker(request,FactStore(directory),Limits(10,600_000_000),retained=witness)
        self.assertEqual(replay["status"],"RETAINED_RESULT_INTEGRITY_ONLY")

    def test_scheduling_does_not_prepare_or_factor_a_number_field(self):
        from run_arithmetic_pipeline import worker
        from research_runtime.supervisor import Limits
        request = {"mode": "features", "curve": {"model": [0, 0, 0, -7, 10]}}
        with TemporaryDirectory() as directory, \
             patch.object(SageArithmetic, "prepare", side_effect=AssertionError("global setup")), \
             patch.object(SageArithmetic, "field", side_effect=AssertionError("maximal order")), \
             patch.object(SageArithmetic, "bnf", side_effect=AssertionError("BNF")):
            result = worker(request, FactStore(directory), Limits(30, 1_073_741_824))
        self.assertEqual(result["status"], "LOCAL_FEATURES_ONLY")
        self.assertFalse(result["features"]["full_bnf_requested"])


if __name__ == "__main__":
    unittest.main()
