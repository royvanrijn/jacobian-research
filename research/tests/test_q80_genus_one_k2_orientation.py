"""Targeted failures for the finite nodal hypotheses and exact local controls."""
import copy
import importlib.util
import os
from pathlib import Path
import resource
import unittest

ROOT=Path(__file__).resolve().parents[1]
PATH=ROOT/'elkies-k3/scripts/verify_q80_genus_one_k2_orientation.py'
spec=importlib.util.spec_from_file_location('q80_orientation',PATH)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
OUT=Path(os.environ.get('Q80_ORIENTATION_TEST_DIR',str(mod.DEFAULT)))

class OrientationControls(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        resource.setrlimit(resource.RLIMIT_CPU,(20,25))
        resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
        cls.source=mod.read(ROOT/mod.SOURCE)
        cls.nodal=mod.read(OUT/'nodal-hypotheses.json')
        cls.symbolic=mod.read(OUT/'symbolic-identities.json')
        cls.controls=mod.read(OUT/'local-controls.json')
        cls.packet=mod.read(OUT/'input.json')

    def test_literal_nodal_and_control_records(self):
        mod.check_nodal(self.nodal,self.source)
        mod.check_symbolic(self.symbolic)
        mod.check_controls(self.controls,self.packet)

    def test_wrong_simple_root_derivative_is_rejected(self):
        row=copy.deepcopy(self.nodal);row['rho_prime']=[6,0]
        with self.assertRaisesRegex(ValueError,'first-derivative'):
            mod.check_nodal(row,self.source)

    def test_nonunit_nodal_slope_is_rejected(self):
        row=copy.deepcopy(self.nodal);row['V_prime']=[0,0]
        with self.assertRaisesRegex(ValueError,'first-derivative'):
            mod.check_nodal(row,self.source)

    def test_unused_root_character_cannot_be_changed(self):
        row=copy.deepcopy(self.nodal);row['norm_codes_node_simple'][1]=1
        with self.assertRaisesRegex(ValueError,'character record'):
            mod.check_nodal(row,self.source)

    def test_symbolic_coefficient_corruption_is_rejected(self):
        row=copy.deepcopy(self.symbolic);row['toy_Q_coefficients'][0][0]='0'
        with self.assertRaisesRegex(ValueError,'coefficient record'):
            mod.check_symbolic(row)

    def test_even_contact_does_not_force_split_branch_field(self):
        row=copy.deepcopy(self.controls);row['rows'][1]['branch_field_over_Qp']='split'
        with self.assertRaisesRegex(ValueError,'control and Hensel'):
            mod.check_controls(row,self.packet)

    def test_insufficient_hensel_precision_is_rejected(self):
        packet=copy.deepcopy(self.packet);packet['hensel_unit_digits']=3
        with self.assertRaisesRegex(ValueError,'retained branch precision'):
            mod.check_controls(self.controls,packet)

    def test_same_orientation_cannot_replace_opposite_offsets(self):
        row=copy.deepcopy(self.controls);row['rows'][2]['valuations']['node_offsets_sum']=1
        with self.assertRaisesRegex(ValueError,'control and Hensel'):
            mod.check_controls(row,self.packet)

    def test_invalid_square_root_seed_is_rejected(self):
        unit=self.controls['rows'][2]['disc_unit_mod_p']
        bad=next(i for i in range(1,131) if i*i%131!=unit)
        with self.assertRaisesRegex(ValueError,'square-root seed'):
            mod.sqrt_hensel(unit,16,bad)

if __name__=='__main__':unittest.main()
