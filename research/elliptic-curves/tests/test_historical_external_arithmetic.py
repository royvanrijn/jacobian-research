"""Cheap cohort-isolation, completeness and missing-data regressions."""
from pathlib import Path
import tempfile
import unittest

import historical_external_arithmetic as lane
import report_historical_external_arithmetic as report
import wide_arithmetic_profile_core as core


def roster():
    return [{"id": str(i), "family": "11952", "parameter": t, "rank_lower_bound": rank,
             "rank_provenance": "EXISTING_CERTIFICATE"} for i, (t, rank) in enumerate(lane.REQUIRED.items())]


def profile(key, rank, value, kind=lane.PROSPECTIVE):
    result = {"curve_key": key, "id": key, "family": "11952", "t": key,
              "final_rank_lower_bound": rank, "rank_bucket": core.rank_bucket(rank),
              "initial_rank_lower_bound": None, "improved_since_initial": None,
              "base_status": "PASS", "local_status": "PASS" if value is not None else "UNKNOWN_TIMEOUT",
              **lane.labels(kind)}
    if value is not None:
        result["bk_local_term"] = value
    return result


class SelectionTests(unittest.TestCase):
    def test_discovers_all_qualifying_not_only_required(self):
        rows = roster() + [{**roster()[0], "id": "extra", "parameter": "3", "rank_lower_bound": 29},
                           {**roster()[0], "id": "low", "parameter": "4", "rank_lower_bound": 26},
                           {**roster()[0], "id": "other", "family": "1092", "parameter": "5", "rank_lower_bound": 31}]
        self.assertEqual({r["id"] for r in lane.select(rows)}, {str(i) for i in range(5)} | {"extra"})

    def test_missing_required_fails(self):
        with self.assertRaisesRegex(RuntimeError, "mandatory"):
            lane.select(roster()[:-1])

    def test_duplicate_rational_parameter_fails(self):
        rows = roster()
        rows.append({**rows[0], "id": "duplicate", "parameter": "-4896/22"})
        with self.assertRaisesRegex(RuntimeError, "duplicate historical parameters"):
            lane.select(rows)

    def test_claimed_28_cannot_drop_to_27(self):
        rows = roster()
        rows[-1]["rank_lower_bound"] = 27
        with self.assertRaisesRegex(RuntimeError, "mandatory"):
            lane.select(rows)


class AnalysisTests(unittest.TestCase):
    def test_historical_never_enters_prospective_diagnostics(self):
        p = [profile("p1", 17, 4), profile("p2", 25, 9)]
        h = [profile("h1", 28, 1000000, lane.HISTORICAL)]
        expected = core.summarize_profiles(p)
        self.assertEqual(report.prospective_summary(p, h), expected)
        h[0]["bk_local_term"] = -1000000
        self.assertEqual(report.prospective_summary(p, h), expected)

    def test_cohort_contamination_fails(self):
        with self.assertRaisesRegex(RuntimeError, "cohort contamination"):
            report.prospective_summary([profile("p", 17, 4, lane.HISTORICAL)], [])

    def test_selection_label_contamination_fails(self):
        p = profile("p", 17, 4)
        p["selection_mode"] = lane.HISTORICAL[1]
        with self.assertRaisesRegex(RuntimeError, "cohort contamination"):
            report.validate_cohorts([p], [])

    def test_duplicate_and_overlap_fail(self):
        p = profile("p", 17, 4)
        with self.assertRaisesRegex(RuntimeError, "duplicate"):
            report.validate_cohorts([p, p], [])
        with self.assertRaisesRegex(RuntimeError, "overlapping"):
            report.validate_cohorts([p], [profile("p", 28, 8, lane.HISTORICAL)])

    def test_unknown_not_imputed_as_zero(self):
        group = report.summarize_group("test", [profile("p1", 17, 8), profile("p2", 19, None)], lane.PROSPECTIVE)
        self.assertEqual(group["median_bk_local_term"], 8)
        self.assertEqual(group["n_bk_local_term"], 1)
        self.assertEqual(group["count"], 2)
        self.assertEqual(group["local_status_counts"]["UNKNOWN_TIMEOUT"], 1)

    def test_all_unknown_stays_unknown(self):
        group = report.summarize_group("test", [profile("p1", 17, None)], lane.PROSPECTIVE)
        self.assertIsNone(group["median_bk_local_term"])
        self.assertEqual(group["n_bk_local_term"], 0)
        self.assertEqual(report.fmt(group["median_bk_local_term"]), "UNKNOWN")

    def test_markdown_escapes_discriminant_absolute_value_pipes(self):
        self.assertEqual(report.table(["log2 |D_K|"], [[3]])[0], "| log2 \\|D_K\\| |")

    def test_group_cannot_mix_cohorts(self):
        with self.assertRaisesRegex(RuntimeError, "mixed"):
            report.summarize_group("test", [profile("p", 17, 5), profile("h", 28, 8, lane.HISTORICAL)], lane.PROSPECTIVE)

    def test_comparison_has_separate_requested_groups(self):
        p = [profile("p17", 17, 5), profile("p23", 23, 7), profile("p24", 24, 9), profile("p25", 25, 9)]
        p[-1]["t"] = "921/653"
        groups = report.comparison(p, [profile("h", 28, 10, lane.HISTORICAL)])
        self.assertEqual(groups[0]["cohort"], "historical_external")
        self.assertEqual(groups[1]["count"], 1)
        self.assertEqual(groups[2]["count"], 3)
        self.assertEqual(groups[6]["count"], 4)
        self.assertTrue(all(g["cohort"] == "prospective_broad_2080" for g in groups[1:]))

    def test_labeled_unknown_retains_status_not_arithmetic(self):
        row = {**roster()[0], "curve_key": "test", "t": "-2448/11", "final_rank_lower_bound": 27,
               "local_search_rank_lower_bound": 27}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            for stage in ("base", "local"):
                (path / stage).mkdir()
                core.write_immutable(path / stage / "test.json", {"status": "UNKNOWN_TIMEOUT", "bk_local_term": 999})
            got = report.labeled_profile(row, path, lane.HISTORICAL)
        self.assertEqual(got["local_status"], "UNKNOWN_TIMEOUT")
        self.assertEqual(got["selection_mode"], "retrospective_known_high_rank")
        self.assertNotIn("bk_local_term", got)

    def test_tree_snapshot_catches_added_changed_removed_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            core.write_immutable(path / "one.json", {"v": 1})
            initial = lane.tree_hashes(path)
            core.write_immutable(path / "two.json", {"v": 2})
            self.assertNotEqual(initial, lane.tree_hashes(path))


if __name__ == "__main__":
    unittest.main()
