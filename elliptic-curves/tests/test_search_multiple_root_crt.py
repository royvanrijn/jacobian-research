from __future__ import annotations

from pathlib import Path
import sys
import unittest


CAS_DIRECTORY = Path(__file__).resolve().parents[1] / "cas"
if str(CAS_DIRECTORY) not in sys.path:
    sys.path.insert(0, str(CAS_DIRECTORY))

from search_multiple_root_crt import (  # noqa: E402
    DEFAULT_GROUP_ORDER,
    choose_groups,
    crt_classes,
    enumerate_candidates,
    exact_group_certificate,
)


class MultipleRootCRTSearchTests(unittest.TestCase):
    def test_default_groups_have_144_unique_crt_classes(self) -> None:
        groups = choose_groups(DEFAULT_GROUP_ORDER, ())
        classes = crt_classes(groups)
        self.assertEqual(len(classes), 144)
        self.assertEqual({item["crt_modulus"] for item in classes}, {6_441_589})
        self.assertEqual(len({item["crt_residue"] for item in classes}), 144)

    def test_group_selection_and_certificates(self) -> None:
        groups = choose_groups(DEFAULT_GROUP_ORDER, (17, 19))
        self.assertEqual(tuple(group.prime for group in groups), (7, 11, 37))
        for group in groups:
            certificate = exact_group_certificate(group)
            self.assertEqual(certificate["prime"], group.prime)
            self.assertTrue(
                all(
                    item["minimal_c4_valuation"] == 0
                    and item["minimal_trace"] == 1
                    for item in certificate["residue_certificates"]
                )
            )
        with self.assertRaises(ValueError):
            choose_groups((7, 7), ())
        with self.assertRaises(ValueError):
            choose_groups((7,), (7,))

    def test_default_bounded_enumeration_pins_shortest_candidate(self) -> None:
        groups = choose_groups(DEFAULT_GROUP_ORDER, ())
        records, counts = enumerate_candidates(
            crt_classes(groups),
            groups,
            coefficient_radius=12,
            representatives_per_class=12,
        )
        self.assertEqual(counts["crt_classes_visited"], 144)
        self.assertEqual(counts["unique_nonsingular_representatives"], 864)
        self.assertEqual(records[0]["t"], "70/223")
        self.assertEqual(records[0]["height"], 223)
        self.assertEqual(
            records[0]["h_valuations"],
            {"7": 18, "11": 5, "17": 5, "19": 4, "37": 3},
        )


if __name__ == "__main__":
    unittest.main()
