"""Rational input normalization and independent conductor replay regressions."""
import sys
from pathlib import Path
from fractions import Fraction
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from run_inventory_conductor_resume_v2 import integral_packet
from v3_warm_support import atomic, read
try:
    from inventory_conductor_worker_v2 import run, verify_transport
    SAGE=True
except ModuleNotFoundError as error:
    if not error.name.startswith('sage'):raise
    SAGE=False


class Normalization(unittest.TestCase):
    def test_integral_unchanged(self):
        source=dict(curve=['0','-1','1','-10','-20'])
        result=integral_packet(source)
        self.assertEqual(result['curve'],source['curve'])
        self.assertEqual(result['integral_model_scale'],'1')

    @unittest.skipUnless(SAGE,'Sage required')
    def test_rational_replay_and_corruption(self):
        with tempfile.TemporaryDirectory() as name:
            d=Path(name)
            for scale in (2,6):
                curve=[str(Fraction(a,scale**w)) for a,w in zip([0,-1,1,-10,-20],(1,2,3,4,6))]
                source=dict(id=str(scale),curve=curve,rank_lower_bound=0,known_primes=[],factor_hints=[])
                packet=integral_packet(source)
                inp=d/f'input-{scale}.json';out=d/f'out-{scale}.json'
                atomic(inp,packet);run(inp,out);run(inp,out,True)
                self.assertEqual(read(out)['exact_conductor'],'11')
                packet['integral_model_scale']='1'
                with self.assertRaisesRegex(ValueError,'transport coefficients differ'):
                    verify_transport(packet)


if __name__=='__main__':unittest.main()
