"""Synthetic status tests: no BNF or expensive research calculation."""
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cas'))
from validate_11952_class_quotient import classify


def provisional(g):
    return ('CERTIFY_ORDER|[]\nPROVISIONAL_CYCLIC_FACTORS|' + json.dumps([2]*g) +
            '\nPROVISIONAL_2RANK|' + str(g) + '\nDONE|provisional\n')


class QuotientGates(unittest.TestCase):
    def test_conditional_not_unconditional(self):
        r = classify('provisional', provisional(16), 'completed', True)
        self.assertEqual(r['grh_conditional_exact_rank'], 25)
        self.assertIsNone(r['unconditional_exact_rank'])

    def test_not_sufficient(self):
        r = classify('provisional', provisional(17), 'completed', True)
        self.assertEqual(r['status'], 'NOT_SUFFICIENT')
        self.assertEqual(r['grh_conditional_rank_upper'], 26)
        self.assertIsNone(r['grh_conditional_exact_rank'])

    def test_inconsistent_low_dimension(self):
        r = classify('provisional', provisional(15), 'completed', True)
        self.assertEqual(r['status'], 'INCONSISTENT_WITH_CERTIFIED_LOWER_BOUND')
        self.assertIsNone(r['grh_conditional_rank_upper'])

    def test_error_after_done_is_not_success(self):
        for error in ['FAIL|failure', '  *** bnfinit: stack overflow', '  *** Warning: unknown warning']:
            r = classify('provisional', provisional(16)+error, 'completed', True)
            self.assertIsNone(r['grh_conditional_rank_upper'])

    def test_only_known_diagnostics_are_permitted(self):
        log = ('  ***   Warning: new maximum stack size = 8589934592 (8192.000 Mbytes).\n'
               '*** Bach constant: 1.5852739882908610802\n') + provisional(16)
        self.assertEqual(classify('provisional', log, 'completed', True)['status'], 'PROVISIONAL_GRH_SUFFICIENT')

    def test_no_binary_or_no_done(self):
        self.assertIsNone(classify('provisional', provisional(16), 'completed', False)['grh_conditional_rank_upper'])
        self.assertIsNone(classify('provisional', provisional(16).replace('DONE|provisional', ''), 'completed', True)['grh_conditional_rank_upper'])

    def test_resource_failure(self):
        for outcome in ['strict_wall_timeout', 'strict_rss_limit']:
            r = classify('provisional', provisional(16), outcome, True)
            self.assertEqual(r['status'], 'RESOURCE_LIMIT')
            self.assertIsNone(r['grh_conditional_exact_rank'])

    def test_duplicate_markers_rejected(self):
        with self.assertRaises(ValueError):
            classify('provisional', provisional(16)+'PROVISIONAL_2RANK|16', 'completed', True)

    def test_successful_quotient(self):
        p = classify('provisional', provisional(16), 'completed', True)
        log = ('PROVISIONAL_CYCLIC_FACTORS|'+json.dumps([2]*16)+
               '\nCLASSGROUP_QUOTIENT_CERTIFIED|1\nG_UPPER|16\nDONE|quotient')
        r = classify('quotient', log, 'completed', True, p)
        self.assertEqual(r['unconditional_exact_rank'], 25)
        self.assertEqual(r['unconditional_class_2rank_upper'], 16)
        for bad in [log.replace('CERTIFIED|1','CERTIFIED|0'), log+'\nFAIL|bad']:
            r = classify('quotient', bad, 'completed', True, p)
            self.assertIsNone(r['unconditional_rank_upper'])
            self.assertEqual(r['grh_conditional_exact_rank'], 25)

    def test_quotient_timeout_preserves_only_conditional_result(self):
        p = classify('provisional', provisional(16), 'completed', True)
        r = classify('quotient', '', 'strict_wall_timeout', True, p)
        self.assertEqual(r['grh_conditional_exact_rank'], 25)
        self.assertIsNone(r['unconditional_exact_rank'])


if __name__ == '__main__':
    unittest.main()
