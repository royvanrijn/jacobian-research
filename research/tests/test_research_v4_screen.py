"""Guard frozen rank-screen coverage and UNKNOWN without running descent."""

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import audit_v4_base_rank_screen as audit


class V4RecordTests(unittest.TestCase):
    def setUp(self):
        self.screen = json.loads((audit.ROOT / audit.SCREEN).read_text())
        self.shortlist = json.loads((audit.ROOT / audit.SHORTLIST).read_text())

    def test_original_records_keep_62_intervals_and_two_unknowns(self):
        summary = audit.validate(self.screen, self.shortlist)
        self.assertEqual((summary['completed'], summary['exact_rank_count'], summary['timeouts']), (62, 17, 2))

    def test_missing_duplicate_and_misattributed_cases_fail(self):
        for mode in ('omitted', 'duplicate', 'misattributed'):
            with self.subTest(mode=mode):
                screen = deepcopy(self.screen)
                if mode == 'omitted': screen['results'].pop()
                elif mode == 'duplicate': screen['results'][2] = deepcopy(screen['results'][1])
                else: screen['results'][2]['pair_key'] = screen['results'][1]['pair_key']
                with self.assertRaises(ValueError): audit.validate(screen, self.shortlist)

    def test_timeout_cannot_gain_bounds_or_become_completed(self):
        for updates in ({'rank_lower_bound': 0}, {'rank_upper_bound': None},
                        {'status': 'completed', 'rank_status': 'EXACT_INTERVAL',
                         'rank_lower_bound': 0, 'rank_upper_bound': 0, 'pari_effort': 0}):
            with self.subTest(updates=updates):
                screen = deepcopy(self.screen)
                screen['results'][2].update(updates)
                with self.assertRaises(ValueError): audit.validate(screen, self.shortlist)

    def test_invalid_intervals_and_invented_summaries_fail(self):
        for mode in ('interval', 'summary', 'ranked-list'):
            with self.subTest(mode=mode):
                screen = deepcopy(self.screen)
                if mode == 'interval': screen['results'][0]['rank_lower_bound'] = 4
                elif mode == 'summary': screen['summary']['completed'] = 64
                else: screen['completed_pairs_ranked_by_upper_then_lower_bound'].reverse()
                with self.assertRaises(ValueError): audit.validate(screen, self.shortlist)

    def test_stale_run_identity_or_changed_effort_fails(self):
        for key, value in [('script_sha256', '0' * 64), ('pari_effort', 1), ('timeout_seconds', 100)]:
            with self.subTest(key=key):
                screen = deepcopy(self.screen)
                screen['results'][2]['run_key'][key] = value
                with self.assertRaises(ValueError): audit.validate(screen, self.shortlist)

    def test_all_selected_products_have_retained_closures(self):
        docs = [json.loads((audit.ROOT / path).read_text()) for path in (audit.CAMPAIGN, audit.SINGLE, audit.SWEEP)]
        audit.validate_product_closures(self.screen, *docs)

    def test_wrong_product_identity_reopened_queue_and_geometric_promotion_fail(self):
        original = [json.loads((audit.ROOT / path).read_text()) for path in (audit.CAMPAIGN, audit.SINGLE, audit.SWEEP)]
        for mode in ('identity', 'queue', 'geometric'):
            with self.subTest(mode=mode):
                docs = deepcopy(original)
                sweep = docs[2]
                if mode == 'identity': sweep['targets'][0]['pair_key'] = docs[1]['pair_key']
                elif mode == 'queue': sweep['explicit_section_solving_queue'] = [sweep['targets'][0]['pair_key']]
                else: sweep['targets'][0]['rank_over_QQbar_u'] = dict(lower=0, upper=0, status='PROVED')
                with self.assertRaises(ValueError): audit.validate_product_closures(self.screen, *docs)


if __name__ == '__main__':
    unittest.main()
