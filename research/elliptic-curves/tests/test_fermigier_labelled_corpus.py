from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import tempfile
import unittest


ECSEARCH = Path(__file__).resolve().parents[1] / "ecsearch"
sys.path.insert(0, str(ECSEARCH))

from fermigier_labelled_corpus import (  # noqa: E402
    canonical_parameter,
    iter_path,
    parameter_record,
    parse_legacy_log,
    split_bucket,
)


class FermigierLabelledCorpusTests(unittest.TestCase):
    def test_coordinate_normalization_and_sign_deduplication(self) -> None:
        self.assertEqual(canonical_parameter("19754/39", coordinate="u"), (39508, 39))
        self.assertEqual(canonical_parameter("-39508/39", coordinate="T"), (39508, 39))
        self.assertEqual(canonical_parameter(Fraction(28917, 20), coordinate="u"), (28917, 10))
        record = parameter_record((39508, 39))
        self.assertEqual(record["adapter_u"], "19754/39")
        self.assertEqual(record["normalized_T_squared"], "1560882064/1521")

    def test_json_path_wildcards_are_deterministic(self) -> None:
        document = {"stages": [{"rows": [{"t": "1/2"}]}, {"rows": [{"t": "2/3"}]}]}
        found = list(iter_path(document, "stages.*.rows.*"))
        self.assertEqual([path for path, _ in found], ["stages.0.rows.0", "stages.1.rows.0"])

    def test_legacy_rank_is_explicitly_uncertified(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "legacy.log"
            path.write_text(
                "=== h20k === H=20000 input=1 promote_rank>=16\n"
                "#0001 u=-19754/39 rank>=18 qpts=60 12.5s\n"
                "#0002 u=1/3 score=1 rank>=? qpts=? status=error 0.0s\n"
            )
            found = list(parse_legacy_log(path, 5_000))
        self.assertEqual(found[0][0], (39508, 39))
        self.assertEqual(found[0][1]["height_bound"], 20_000)
        self.assertEqual(found[0][1]["reported_rank_floor_uncertified"], 18)
        self.assertEqual(found[1][1]["status"], "error")
        self.assertIsNone(found[1][1]["reported_rank_floor_uncertified"])

    def test_split_is_stable(self) -> None:
        value = split_bucket((17, 29), "corpus-test")
        self.assertEqual(value, split_bucket((17, 29), "corpus-test"))
        self.assertIn(value, {"train", "validation", "internal_test"})


if __name__ == "__main__":
    unittest.main()
