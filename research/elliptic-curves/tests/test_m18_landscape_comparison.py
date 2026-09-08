"""M18 boundary, immutable labels, and missing-order-statistic regressions."""
import copy
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

CAS=Path(__file__).resolve().parents[1]/'cas'
sys.path.insert(0,str(CAS))
import m18_landscape_comparison as m
import m18_landscape_outcomes as labels

class LandscapeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        protocol=m.read(m.OUT/'protocol.json')
        cls.case=protocol['cases'][0]['case'];cls.folder=m.SNAP/cls.case
        cls.seed=m.read(cls.folder/'seed-input.json')
        cls.selection=m.read(cls.folder/'selection.json')

    def test_rejects_later_basis(self):
        seed=copy.deepcopy(self.seed);seed['points'].append(seed['points'][0])
        with self.assertRaisesRegex(ValueError,'M18'):m.features(seed,self.selection,self.folder,replay=False)
        selection=copy.deepcopy(self.selection);selection['basis'][17]=selection['basis'][0]
        with self.assertRaisesRegex(ValueError,'boundary'):m.features(self.seed,selection,self.folder,replay=False)

    def test_only_npz_reads_and_absent_128(self):
        # The feature function is allowed only its M18 values and 32 frozen NPZs.
        # A read of a transcript, label or higher-rank packet fails immediately.
        allowed={str(p.resolve()) for p in self.folder.glob('anchor-*-full.npz')}
        import builtins
        real_open=builtins.open;real_path_open=Path.open
        def guard(path):
            if isinstance(path,(str,Path)):
                self.assertIn(str(Path(path).resolve()),allowed)
        def guarded_open(path,*args,**kwargs):guard(path);return real_open(path,*args,**kwargs)
        def guarded_path_open(path,*args,**kwargs):guard(path);return real_path_open(path,*args,**kwargs)
        with patch('builtins.open',guarded_open),patch.object(Path,'open',guarded_path_open):
            result=m.features(self.seed,self.selection,self.folder,replay=False)
        expected=m.read(m.OUT/'features'/f'{self.case}.json')
        self.assertEqual(result['features'],expected['features'])
        self.assertIsNone(result['features']['gap_128'])
        self.assertEqual(result['features']['cosets'],64)

    def test_gram_tampering_fails(self):
        selection=copy.deepcopy(self.selection);selection['rounded_gram'][0][1]+=1
        with self.assertRaisesRegex(ValueError,'asymmetric'):m.features(self.seed,selection,self.folder,replay=False)

    def test_labels_require_seal(self):
        with tempfile.TemporaryDirectory() as d,patch.object(labels,'OUT',Path(d)):
            with self.assertRaises(FileNotFoundError):labels.joined()

    def test_labels_reject_feature_tampering(self):
        with tempfile.TemporaryDirectory() as d,patch.object(labels,'OUT',Path(d)):
            out=Path(d);(out/'protocol.json').write_text('{}')
            (out/'feature.json').write_text('modified')
            m.write(out/'feature-seal.json',{'protocol_sha256':m.sha(out/'protocol.json'),'features':{'feature.json':'incorrect'}})
            with self.assertRaisesRegex(ValueError,'features changed'):labels.joined()

    def test_nearest_rank_and_entropy(self):
        self.assertEqual(m.quantile([4,1,3,2],1,2),2)
        self.assertEqual(m.entropy([2,2]),1)
        self.assertEqual(m.entropy([4]),0)
        self.assertIsNone(m.entropy([]))

if __name__=='__main__':unittest.main()
