"""Exact model, point-chart and resource-failure regressions (sage -python)."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from sage.all import EllipticCurve, QQ, matrix, vector
import run_fixed_field_radical_covers as m


class RadicalModels(unittest.TestCase):
    def setUp(self):
        self.row={'reduced_quartic_model':[['1','-1','0','0','1'],['0']*3]}
        self.E=EllipticCurve([0,0,0,-4,1])
        self.var=matrix(QQ,[[1,1,0,0],[0,1,0,0],[0,0,1,0],[0,0,1,1]])
        self.eq=matrix(QQ,[[1,2],[0,1]])
        old=m.initial_quadrics(self.row['reduced_quartic_model'][0])
        new=[sum((self.eq[i,j]*self.var*old[j]*self.var.transpose()
                  for j in range(2)),matrix(QQ,4)) for i in range(2)]
        self.model={'quadric_matrices':[[str(z) for z in a.list()] for a in new],
                    'equation_transform':list(map(str,self.eq.list())),
                    'variable_transform':list(map(str,self.var.list()))}

    def test_exact_transformation_and_minimality(self):
        qs,var=m.verify_quadric_model(self.row,self.model,self.E)
        for original in ([0,0,1,1],[1,0,0,1]):
            new=vector(QQ,original)*var.inverse()
            self.assertTrue(all(new*q*new==0 for q in qs))
            s,t,y=m.quadric_to_quartic(new*var)
            self.assertEqual(y*y,m.base.homogeneous(self.row['reduced_quartic_model'][0],s,t,4))

    def test_reject_changed_map(self):
        bad=copy.deepcopy(self.model)
        bad['variable_transform'][0]='2'
        with self.assertRaises(AssertionError):m.verify_quadric_model(self.row,bad,self.E)

    def test_reject_wrong_jacobian(self):
        with self.assertRaises(AssertionError):
            m.verify_quadric_model(self.row,self.model,EllipticCurve([0,0,0,-4,-1]))

    def test_cross_terms_cannot_be_silently_dropped(self):
        self.row['reduced_quartic_model'][1][0]='1'
        with self.assertRaises(AssertionError):m.verify_quadric_model(self.row,self.model,self.E)
        with self.assertRaises(AssertionError):m.magma_input(self.row)

    def test_infinity_and_projective_scaling(self):
        self.assertEqual(m.quadric_to_quartic([7,0,0,-14]),(1,0,-2))
        self.assertEqual(m.quadric_to_quartic([8,12,18,24]),(QQ(2)/3,1,QQ(4)/3))
        with self.assertRaises(AssertionError):m.quadric_to_quartic([0,0,0,1])

    def test_done_does_not_hide_resource_error(self):
        raw='<calculator><results><line>System Error: User memory limit has been reached</line><line>DONE</line></results></calculator>'
        self.assertFalse(m.parse_xml(raw)[0])
        self.assertEqual(m.descent_status(raw),'MEMORY_LIMIT_NO_DESCENT_RESULT')
        with self.assertRaises(AssertionError):m.parse_model(raw)

    def test_timeout_is_not_empty_selmer_set(self):
        raw='<calculator><headers><warning>The computation exceeded the time limit</warning></headers><results><line>FOUR_DESCENT_COUNT 0</line><line>DONE</line></results></calculator>'
        self.assertEqual(m.descent_status(raw),'TIME_LIMIT_NO_DESCENT_RESULT')

    def test_output_fractions_and_forbidden_syntax(self):
        self.assertEqual(m.literal('POINTS [ [ -3/7, 1, 0, 5 ] ]','POINTS'),[[QQ(-3)/7,1,0,5]])
        with self.assertRaises(AssertionError):m.literal('POINTS [danger()]','POINTS')


if __name__=='__main__':unittest.main()
