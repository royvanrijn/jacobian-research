"""Reject changes to the finite NS0031 witnesses without importing a CAS."""

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'elkies-k3/scripts'))
import verify_ns0031_marking_arithmetic as replay


class NS0031ReplayTests(unittest.TestCase):
    def setUp(self):
        self.packet=json.loads(replay.PACKET.read_text())
        self.certificate=json.loads(replay.CERTIFICATE.read_text())

    def test_optimized_python_cannot_report_an_unchecked_pass(self):
        result=subprocess.run([sys.executable,'-O',replay.__file__],capture_output=True,text=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('requires assertions',result.stderr)
        self.assertNotIn('PASS',result.stdout)

    def test_portable_arithmetic_needs_no_original_catalogues(self):
        self.assertIn('PASS NS0031 finite arithmetic',replay.verify(self.packet,self.certificate))
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(AssertionError,'missing original input'):
                replay.check_source_projection(self.packet,Path(directory))

    def test_wrong_frobenius_trace_or_curve_is_rejected(self):
        for field,value in [('frobenius_trace',6),('minimal_ainvariants',[1,1,1,-8,7])]:
            with self.subTest(field=field):
                certificate=deepcopy(self.certificate)
                certificate['x0_37_rational_point_gate']['excluded_lifts'][0][field]=value
                with self.assertRaises(AssertionError):
                    replay.verify(self.packet,certificate)

    def test_bad_cartan_signature_and_clifford_model_are_rejected(self):
        certificate=deepcopy(self.certificate)
        certificate['norm_one_modular_curve']['cusp_widths']=[4,4,147,149]
        with self.assertRaises(AssertionError):replay.verify(self.packet,certificate)
        certificate=deepcopy(self.certificate)
        certificate['even_clifford_order']['split_basis_before_conjugation'][3][1][0]='147'
        with self.assertRaises(AssertionError):replay.verify(self.packet,certificate)

    def test_changed_projection_is_not_accepted_as_the_original_input(self):
        packet=deepcopy(self.packet)
        packet['sources'][0]['fields']['surface_key/transcendental_gram'][1][1]=76
        with self.assertRaises(AssertionError):replay.verify(packet,self.certificate)
        packet=deepcopy(self.packet)
        packet['sources'][0]['sha256']='0'*64
        with self.assertRaises(AssertionError):replay.verify(packet,self.certificate)


if __name__=='__main__':unittest.main()
