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
    def test_fixed_budget_and_shells(self):
        self.assertEqual(c.FRESH_CHARTS,512)
        self.assertEqual(c.SHELLS,(4,6,8,10,12))
        self.assertNotIn(0,c.SHELLS)

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

    def test_claim_does_not_promote_null_to_upper_bound(self):
        text=c.CLAIM.lower()
        self.assertIn('not a rank upper bound',text)
        self.assertIn('does not retune or execute',text)

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
                       'det1092_v4_bootstrap_replay.py','run_det1092_v4_bootstrap.py'):
            self.assertIn('elliptic-curves/cas/'+suffix,sources)
        self.assertTrue(all((ROOT/path).is_file() for path in sources))


if __name__=='__main__':unittest.main()
