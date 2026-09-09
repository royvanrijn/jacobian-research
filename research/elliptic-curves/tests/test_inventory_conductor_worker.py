"""Run with Sage -python; exact and nonminimal arithmetic regressions."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
try:
    from inventory_conductor_worker import run
    AVAILABLE=True
except ModuleNotFoundError as error:
    if not error.name.startswith('sage'):raise
    AVAILABLE=False
from v3_warm_support import atomic, read


@unittest.skipUnless(AVAILABLE,'Sage required')
class Conductors(unittest.TestCase):
    def test_minimal_and_nonminimal(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d)
            for scale in (1,2):
                packet=d/f'input-{scale}.json';out=d/f'out-{scale}.json'
                curve=[str(a*scale**w) for a,w in zip([0,-1,1,-10,-20],[1,2,3,4,6])]
                atomic(packet,dict(id=str(scale),curve=curve,rank_lower_bound=0,known_primes=[],factor_hints=[]))
                run(packet,out);run(packet,out,True)
                self.assertEqual(read(out)['exact_conductor'],'11')
            saved=read(out)
            # A genuine partial certificate must replay as UNKNOWN.
            saved['local_data']=saved['local_data'][:1]
            q=saved['local_data'][0]
            remaining=int(saved['discriminant_abs'])//int(q['prime'])**q['discriminant_valuation']
            divisor=int(q['prime'])**q['conductor_exponent']
            saved.update(status='UNKNOWN',remaining_cofactor=str(remaining),conductor_divisor=str(divisor),
                         conductor_upper_bound=str(divisor*remaining),exact_conductor=None)
            atomic(out,saved);run(packet,out,True)
            saved['status']='EXACT';saved['exact_conductor']='11';atomic(out,saved)
            with self.assertRaisesRegex(Exception,'independent local replay differs'):run(packet,out,True)


if __name__=='__main__':unittest.main()
