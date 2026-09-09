"""Displayed minima are per rank/column, include ties and never promote dashes."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from curve_table_highlights import highlight_rank_minima


class Highlights(unittest.TestCase):
    def test_groups_ties_and_missing(self):
        text='\n'.join([
            '| Curve | a-invariants | Rank | log N | Naive height | Faltings height | log abs(Δ) |',
            '| A | x | ≥ 28 | 318.98 | 420.19 | 32.92 | 409.64 |',
            '| B | x | ≥ 27 | 275.84 | 350.00 | 26.00 | — |',
            '| C | x | ≥ 27 | 280.00 | 340.00 | 26.00 | 331.00 |'])
        result=highlight_rank_minima(text)
        self.assertIn('| A | x | ≥ 28 | **318.98** | **420.19** | **32.92** | **409.64** |',result)
        self.assertIn('| B | x | ≥ 27 | **275.84** | 350.00 | **26.00** | — |',result)
        self.assertIn('| C | x | ≥ 27 | 280.00 | **340.00** | **26.00** | **331.00** |',result)
        self.assertEqual(result,highlight_rank_minima(result))

    def test_filtered_subset_recomputes(self):
        row='| B | x | ≥ 18 | 100.00 | — | -1.00 | — |'
        self.assertEqual(highlight_rank_minima(row),
                         '| B | x | ≥ 18 | **100.00** | — | **-1.00** | — |')

    def test_stale_bold_removed(self):
        rows='| A | x | ≥ 27 | **5.00** | — | — | — |\n| B | x | ≥ 27 | 4.00 | — | — | — |'
        result=highlight_rank_minima(rows)
        self.assertNotIn('**5.00**',result)
        self.assertIn('**4.00**',result)


if __name__=='__main__':unittest.main()
