"""Tests distinguishing local and simultaneous generic images."""
import unittest
import retrospective as r
import bad_prime_support as b


def example(second_generic=False):
    return {'id':'synthetic','generic_dimension':2 if second_generic else 1,'witness_dimension':2,
      'local':[{'prime':3,'point_signature_rows':[[0,1],[0,1]],'point_kummer_dimension':1,'local_reduction':{}},
               {'prime':5,'point_signature_rows':[[0,1],[0,0]],'point_kummer_dimension':1,'local_reduction':{}}]}


class BadPrimeSupportTests(unittest.TestCase):
    def test_individual_surjectivity_does_not_imply_product_surjectivity(self):
        row=b.characterize(example())
        self.assertTrue(all(p['generic_image_is_full_point_image'] for p in row['local']))
        self.assertFalse(row['generic_image_is_full_product_point_image'])
        self.assertEqual(row['joint_quotient_support_dimension'],1)
        self.assertIsNone(row['simultaneous_generic_corrections'][0]['generic_correction_mask'])

    def test_extra_generic_direction_completes_product(self):
        row=b.characterize(example(True))
        self.assertTrue(row['generic_image_is_full_product_point_image'])
        self.assertEqual(row['joint_generic_image_dimension'],2)

    def test_reported_corrections_cancel_every_local_matrix(self):
        for record in r.read(b.INPUT)['cases']:
            row=b.characterize(record);m=record['generic_dimension']
            for correction in row['simultaneous_generic_corrections']:
                mask=correction['generic_correction_mask'];self.assertIsNotNone(mask)
                for local in record['local']:
                    sig=list(map(r.pack,local['point_signature_rows']))
                    v=sig[m+correction['quotient_index']]
                    for i in range(m):
                        if mask>>i&1:v^=sig[i]
                    self.assertEqual(v,0)


if __name__=='__main__':unittest.main()
