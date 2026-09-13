import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('tangent_replay_test',ROOT/'elkies-k3/scripts/verify_r17_shared_ordinate_tangents.py')
V=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(V)


class TangentGate(unittest.TestCase):
    def copied(self):
        tmp=tempfile.TemporaryDirectory();self.addCleanup(tmp.cleanup);path=Path(tmp.name)
        for name in ('input.json','result.json'):shutil.copyfile(V.DEFAULT/name,path/name)
        return path

    def test_missing_pair_is_rejected(self):
        path=self.copied();result=json.loads((path/'result.json').read_text());result['rows'].pop()
        (path/'result.json').write_text(json.dumps(result))
        with self.assertRaisesRegex(ValueError,'complete unique basis-pair coverage'):V.verify(path)

    def test_altered_branch_witness_is_rejected(self):
        path=self.copied();result=json.loads((path/'result.json').read_text())
        result['rows'][0]['witness']['odd_part'][0]+=1
        (path/'result.json').write_text(json.dumps(result))
        with self.assertRaisesRegex(ValueError,'altered finite branch witness'):V.verify(path)

    def test_square_factors_and_colliding_branch_points(self):
        # x*(x-1)^2*(x-2)^3 has odd part x*(x-2), not its radical.
        p=101;mul=lambda a,b:V.V.ff_mul(a,b,p)
        f=mul(mul([0,1],mul([-1%p,1],[-1%p,1])),mul(mul([-2%p,1],[-2%p,1]),[-2%p,1]))
        self.assertEqual(V.odd_part(f,p),[0,-2%p,1])
        # The two characteristic-zero branch points of (x-1)*(x+4)
        # collide at5; the bound must allow this loss under reduction.
        self.assertEqual(V.odd_part([1,3,1],5),[1])

    def test_nonintegral_input_is_rejected(self):
        packet=json.loads((V.DEFAULT/'input.json').read_text());packet['A'][0]='1/1009'
        with self.assertRaisesRegex(ValueError,'nonintegral generic coefficient'):
            V.finite_witness(packet,[0,1],1009)


if __name__=='__main__':unittest.main()
