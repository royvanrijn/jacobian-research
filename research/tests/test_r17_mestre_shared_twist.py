"""Failure controls for the section and root-field proof obligations."""
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('mestre_replay',ROOT/'elkies-k3/scripts/verify_r17_mestre_shared_twist.py')
V=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(V)


class MestreReplay(unittest.TestCase):
    def test_reducible_and_inseparable_modular_polynomials_are_rejected(self):
        self.assertFalse(V.irreducible([1,0,1],5))
        self.assertFalse(V.irreducible([1,0,1],2))
        self.assertTrue(V.irreducible([1,0,1],3))

    def test_composite_modulus_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'nonprime'):
            V.reduction([V.Q(1),V.Q(1)],9)

    def test_degree_drop_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'degree drop'):
            V.reduction([V.Q(1),V.Q(29)],29)

    def replay_mutation(self,filename,mutate):
        with tempfile.TemporaryDirectory() as name:
            path=Path(name)
            for f in ('input.json','result.json','branch-gate.json'):shutil.copyfile(V.DEFAULT/f,path/f)
            data=json.loads((path/filename).read_text());mutate(data);(path/filename).write_text(json.dumps(data))
            return V.verify(path)

    def test_altered_new_section_is_rejected(self):
        def mutate(d):d['points'][1]['y_radical_numerator'][0]='2'
        with self.assertRaisesRegex(ValueError,'section formula'):
            self.replay_mutation('result.json',mutate)

    def test_wrong_local_root_is_rejected(self):
        def mutate(d):d['root_field_certificates'][0]['degree_one_prime']['root']=52
        with self.assertRaisesRegex(ValueError,'degree-one local place'):
            self.replay_mutation('branch-gate.json',mutate)

    def test_producer_height_hint_is_not_used(self):
        def mutate(d):d['height_certificate_inputs']['height_matrix']=[[0,0],[0,0]]
        result=self.replay_mutation('result.json',mutate)
        self.assertEqual(result['height_matrix'],[[24,0],[0,24]])
        self.assertFalse(result['goal_complete'])


if __name__=='__main__':unittest.main()
