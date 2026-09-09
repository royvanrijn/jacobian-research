"""README curation must not confuse UNKNOWN conductors with large ones."""
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'cas'))
from render_main_readme_curves import retention_reason, conductor_benchmarks, select_section, short_introduction


class SelectionTests(unittest.TestCase):
    def test_concise_introduction(self):
        text=short_introduction(dict(count=321,conductor_status_counts=dict(EXACT=195,UNKNOWN=126)))
        self.assertIn('195 exact conductors',text)
        self.assertIn('certified lower bounds',text)
        self.assertIn('Methods and selection',text)
        self.assertNotIn('structural examples',text)
        self.assertLess(len(text.split()),75)

    def row(self, rank=21, status='UNKNOWN', conductor=None, identifier='ordinary'):
        return dict(id=identifier, local_search_rank_lower_bound=rank,
                    conductor_status=status, conductor=conductor)

    def test_threshold_and_exception(self):
        self.assertTrue(retention_reason(self.row(22), {}))
        self.assertIsNone(retention_reason(self.row(), {}))
        self.assertTrue(retention_reason(self.row(identifier='det1092-small-conic'), {}))

    def test_exact_comparison_only(self):
        self.assertTrue(retention_reason(self.row(status='EXACT', conductor='100'), {21:100}))
        self.assertIsNone(retention_reason(self.row(status='EXACT', conductor='101'), {21:100}))
        self.assertIsNone(retention_reason(self.row(conductor='1'), {21:100}))
        self.assertIsNone(retention_reason(self.row(status='EXACT', conductor='1'), {}))

    def test_benchmark_is_at_least_rank(self):
        b = conductor_benchmarks([dict(rank_lower_bound=22, conductor='90'),
                                  dict(rank_lower_bound=20, conductor='10')])
        self.assertEqual(b[21], 90)
        self.assertEqual(b[20], 10)

    def test_filter_does_not_mutate_archive(self):
        rows = [self.row(), self.row(22, identifier='high')]
        section = '| [ordinary](data/research_curves/ordinary.md) | x |\n| [high](data/research_curves/high.md) | y |'
        result, kept = select_section(section, rows, {})
        self.assertEqual(kept, {'high'})
        self.assertNotIn('ordinary', result)
        self.assertEqual(len(rows), 2)
        with self.assertRaises(ArithmeticError):
            select_section('', rows, {})


if __name__ == '__main__':
    unittest.main()
