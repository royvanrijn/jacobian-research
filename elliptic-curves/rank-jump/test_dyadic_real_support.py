"""Exact real-component regression and seven-place corrections."""
import unittest
import retrospective as r
import bad_prime_support as bad
import dyadic_real_support as dr


class DyadicRealSupportTests(unittest.TestCase):
    def test_exact_real_components(self):
        # y^2=x^3-2x has a bounded oval and an unbounded component.
        self.assertEqual(dr.real_bit(-2,0,-1),1)
        self.assertEqual(dr.real_bit(-2,0,0),1)
        self.assertEqual(dr.real_bit(-2,0,2),0)
        # Negative discriminant: the real elliptic curve is connected.
        self.assertEqual(dr.real_bit(-1,1,0),0)
        with self.assertRaises(ValueError):dr.real_bit(0,0,1)

    def test_all_simultaneous_corrections_cancel_the_seven_images(self):
        for record in r.read(dr.INPUT)['cases']:
            combined=dr.combine(record);out=bad.characterize(combined);m=record['generic_dimension']
            for correction in out['simultaneous_generic_corrections']:
                mask=correction['generic_correction_mask'];self.assertIsNotNone(mask)
                for local in combined['local']:
                    sig=list(map(r.pack,local['point_signature_rows']))
                    bit=sig[m+correction['quotient_index']]
                    for i in range(m):
                        if mask>>i&1:bit^=sig[i]
                    self.assertEqual(bit,0)

    def test_sage_rational_interoperability_when_available(self):
        try:from sage.all import QQ
        except ImportError:self.skipTest('optional Sage rational conversion regression')
        self.assertEqual(dr.real_bit(QQ(-2),QQ(0),QQ(-1)),1)
        self.assertEqual(dr.real_bit(QQ(-2),QQ(0),QQ(2)),0)


if __name__=='__main__':unittest.main()
