import importlib.util
import py_compile
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
CAS=ROOT/'elliptic-curves/cas'
sys.path.insert(0,str(CAS))

import det1092_v4_bootstrap_contract as c
import det1092_v4_bootstrap_worker as w


class V4ContractTests(unittest.TestCase):
    def test_fixed_budget_and_exported_survivor_contract(self):
        self.assertEqual(c.FRESH_CHARTS,512)
        self.assertEqual(c.SURVIVOR_SHELLS,(8,10,12))
        self.assertEqual(c.EXPORTED_SURVIVOR_COUNTS,{8:63922,10:40917,12:139})
        self.assertEqual(sum(c.EXPORTED_SURVIVOR_COUNTS.values()),104978)
        self.assertLess(sum(c.EXPORTED_SURVIVOR_COUNTS.values()),1<<17)

    def test_full_quotient_is_generated_from_masks_not_tsv_rows(self):
        prior={1,7,131071}
        fresh=c.fresh_masks(prior)
        self.assertEqual(len(fresh),(1<<17)-1-len(prior))
        self.assertNotIn(0,fresh)
        self.assertTrue(prior.isdisjoint(fresh))
        self.assertEqual(fresh[0],2)
        self.assertEqual(fresh[-1],131070)

    def test_parity_mask_is_sign_and_even_translation_invariant(self):
        word=[0,1,-2,3]+[0]*13
        shifted=[x+2*(i-3) for i,x in enumerate(word)]
        negated=[-x for x in word]
        self.assertEqual(c.parity_mask(word),c.parity_mask(shifted))
        self.assertEqual(c.parity_mask(word),c.parity_mask(negated))
        self.assertNotEqual(c.parity_mask(word),0)

    def test_hash_domain_is_case_and_parity_bound(self):
        a=c.hash_key('case-a',17)
        self.assertEqual(a,c.hash_key('case-a',17))
        self.assertNotEqual(a,c.hash_key('case-b',17))
        self.assertNotEqual(a,c.hash_key('case-a',18))

    def test_even_positions_deterministic_and_in_range(self):
        order=list(range(1000))
        out=w.even_positions(order,128)
        self.assertEqual(len(out),128)
        self.assertEqual(out,w.even_positions(order,128))
        self.assertEqual(len(set(out)),128)
        self.assertTrue(all(0<=x<1000 for x in out))

    def test_even_positions_keeps_small_population(self):
        self.assertEqual(w.even_positions([7,2,9],10),[7,2,9])

    def test_old_v4_evidence_namespace_is_not_reused(self):
        self.assertEqual(c.D.name,'det1092-v4-wide-bootstrap-v2')
        self.assertEqual(c.AUTO.name,'det1092-v4-wide-controller-v2')
        self.assertEqual(c.DOMAIN,'det1092-v4-wide-bootstrap-v2')

    def test_claim_does_not_promote_null_to_upper_bound(self):
        text=c.CLAIM.lower()
        self.assertIn('not a rank upper bound',text)
        self.assertIn('tsv is used only',text)
        self.assertIn('does not retune or execute',text)

    def test_worker_does_not_require_survivor_tsv_to_cover_full_quotient(self):
        source=(CAS/'det1092_v4_bootstrap_worker.py').read_text()
        self.assertNotIn('orbit table is not the complete parity quotient',source)
        self.assertIn('c.fresh_masks(prior_masks)',source)
        self.assertIn('exported_survivors(ctx)',source)
        self.assertIn("'complete_parity_classes':1<<17",source)

    def test_exact_cvp_censor_has_deterministic_reserve(self):
        source=(CAS/'det1092_v4_bootstrap_worker.py').read_text()
        self.assertIn('cvp_censored.append(mask)',source)
        self.assertIn('candidate_order=initial+[m for m in hashed if m not in chosen]',source)
        self.assertIn("'cvp-reserve'",source)

    def test_all_v4_python_sources_compile(self):
        for name in ('det1092_v4_bootstrap_contract.py','det1092_v4_bootstrap_worker.py',
                     'det1092_v4_bootstrap_replay.py','run_det1092_v4_bootstrap.py'):
            py_compile.compile(str(CAS/name),doraise=True)

    def test_controller_imports_without_sage(self):
        name='v4_controller_import_test'
        spec=importlib.util.spec_from_file_location(name,CAS/'run_det1092_v4_bootstrap.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        self.assertTrue(callable(module.status))
        self.assertTrue(callable(module.launch))

    def test_source_manifest_is_complete(self):
        sources=c.own_sources()
        for suffix in ('det1092_v4_bootstrap_contract.py','det1092_v4_bootstrap_worker.py',
                       'det1092_v4_bootstrap_replay.py','run_det1092_v4_bootstrap.py',
                       'visibility_lattice_v2.py'):
            self.assertIn('elliptic-curves/cas/'+suffix,sources)
        self.assertTrue(all((ROOT/path).is_file() for path in sources))


if __name__=='__main__':unittest.main()
