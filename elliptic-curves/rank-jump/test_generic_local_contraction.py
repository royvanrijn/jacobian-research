import itertools
import unittest
import generic_local_contraction as local
import generic_contraction_consequences as consequences
import retrospective as r


def elements(basis):
    return {local.xor_selected(mask,basis) for mask in range(1<<len(basis))}


class GenericContractionTests(unittest.TestCase):
    def test_right_inverse_and_failed_gate(self):
        # The two local products have a nonzero common direction.
        result=local.linear_certificate([1,2,4],[1,2,4],[4,8,16])
        self.assertTrue(result['generic_surjectivity_certified'])
        self.assertEqual(result['local_change_dimension'],2)
        failed=local.linear_certificate([1,4],[1,2,4],[4,8,16])
        self.assertFalse(failed['generic_surjectivity_certified'])
        self.assertEqual(failed['original_basis_generic_correction_masks'],[])

    def test_all_binary_lagrangians_dimension_four(self):
        pair=lambda x,y:(((x&3)&(y>>2)).bit_count()+((y&3)&(x>>2)).bit_count())%2
        lagrangians={frozenset(elements([x,y])) for x,y in itertools.combinations(range(1,16),2)
                     if r.rank([x,y])==2 and pair(x,y)==0}
        self.assertEqual(len(lagrangians),15)
        original={0,1,2,3};twist={0,4,8,12}
        for boundary in lagrangians:
            a=boundary&original;b=boundary&twist
            for generic in [frozenset({0})]+[frozenset({0,x}) for x in a if x]+[a]:
                e=(len(generic)).bit_length()-1
                drop=len(a).bit_length()-len(b).bit_length()
                self.assertGreaterEqual(drop,2*e-2)
                self.assertLessEqual(drop,2)
                if e==2:
                    self.assertTrue(b<=a);self.assertEqual(drop,2)

    def test_counting_baseline_exhaustively(self):
        spaces={frozenset(elements([x,y])) for x,y in itertools.combinations(range(1,16),2) if r.rank([x,y])==2}
        self.assertEqual(len(spaces),consequences.gaussian(4,2))
        # Projection onto the first two coordinates.
        onto=sum(len({x&3 for x in s})==4 for s in spaces)
        self.assertEqual(onto,consequences.surjective_count(4,2,2))

    def test_coordinates_allow_dependent_inputs_and_reject_outside(self):
        columns=[1,2,3,0]
        self.assertEqual(local.xor_selected(local.coordinates(3,columns),columns),3)
        with self.assertRaises(ValueError):local.coordinates(4,columns)


if __name__=='__main__':unittest.main()
