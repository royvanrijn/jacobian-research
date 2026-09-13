"""Projective-boundary failure controls for the family comparison."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('projective_node_replay',ROOT/'elkies-k3/scripts/verify_r17_one_node_complete_pairs.py')
V=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(V)


class ProjectiveNodeReplay(unittest.TestCase):
    def test_source_includes_parameter_line_at_infinity(self):
        C=[[0],[0],[0],[0],[1,0,1]];p=5
        points=list(V.source_image(C,p))
        self.assertEqual(len(points),p*p+p+1)
        self.assertEqual(points[-1],(1,0,1,0,0,0,0))

    def test_zero_source_value_remains_a_base_point(self):
        C=[[0],[0],[0],[0],[1,0,1]]
        points=list(V.source_image(C,5))
        self.assertIsNone(points[0])
        self.assertEqual(sum(q is None for q in points),1)

    def test_target_includes_both_infinite_parameters(self):
        B=[[int(i==j) for j in range(5)] for i in range(5)]
        points=list(V.target_image(B,5))
        self.assertEqual(len(points),36)
        self.assertEqual(points[-1],(0,0,0,0,1,0,0))

    def test_nonintegral_source_coefficients_are_deferred(self):
        self.assertIsNone(V.reduce_channels([['1/101'],[],[],[],['1']],101))

    def test_projective_scalars_cannot_hide_a_match(self):
        self.assertEqual(V.normalize([1,2,0,3,0,1,4],101),V.normalize([7,14,0,21,0,7,28],101))

    def test_missing_completion_certificate_cannot_close_the_last_pair(self):
        with tempfile.TemporaryDirectory() as name:
            with self.assertRaisesRegex(ValueError,'final pair remains UNKNOWN'):
                V.verify_continuation(Path(name),{}, {(27,3)})


if __name__=='__main__':unittest.main()
