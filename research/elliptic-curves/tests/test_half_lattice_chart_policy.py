from __future__ import annotations

from copy import deepcopy
import sys
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))

from half_lattice_chart_policy import (  # noqa: E402
    MISS_INFERENCES,
    bind_ordering,
    policy_document,
    validate_ordering,
)


def fixture() -> dict:
    return {
        "basis_records": [
            {"x": "1", "y": "2"},
            {"x": "3", "y": "4"},
        ],
        "height_gram_rows": [["2", "1/2"], ["1/2", "3"]],
        "generic_coordinate_rows": [[1, 0]],
        "quotient_coordinate_rows": [[0, 1]],
        "chart_universe_id": "example:M2/2M2:all",
        "ordered_chart_ids": ["mask:3", "mask:1", "mask:2"],
        "heuristics": [
            "legacy_half_lattice_depth",
            "old_deep_43",
            "quotient_hamming_weight",
        ],
    }


class HalfLatticeChartPolicyTests(unittest.TestCase):
    def test_semantics_forbid_arithmetic_or_selmer_meaning(self) -> None:
        policy = policy_document()
        self.assertFalse(policy["quartic_is_nontrivial_2_covering"])
        self.assertFalse(policy["quartic_represents_a_selmer_class"])
        self.assertTrue(all(value is False for value in MISS_INFERENCES.values()))
        meanings = str(policy["heuristic_meanings"])
        self.assertIn("search charts", meanings)
        self.assertIn("Selmer", meanings)
        self.assertIn("basis-invariant", meanings)

    def test_exact_state_and_order_validate(self) -> None:
        data = fixture()
        certificate = bind_ordering(**data)
        validate_ordering(
            certificate,
            **{key: value for key, value in data.items() if key != "heuristics"},
        )
        self.assertFalse(certificate["calibration_transfers_to_a_changed_state"])

    def test_lattice_enlargement_invalidates_order(self) -> None:
        data = fixture()
        certificate = bind_ordering(**data)
        changed = deepcopy(data)
        changed["basis_records"].append({"x": "5", "y": "6"})
        changed["height_gram_rows"] = [
            ["2", "1/2", "0"],
            ["1/2", "3", "0"],
            ["0", "0", "5"],
        ]
        with self.assertRaisesRegex(ValueError, "stale chart ordering"):
            validate_ordering(
                certificate,
                **{
                    key: value
                    for key, value in changed.items()
                    if key != "heuristics"
                },
            )

    def test_basis_or_quotient_coordinate_change_invalidates_order(self) -> None:
        data = fixture()
        certificate = bind_ordering(**data)
        for field, replacement in (
            ("basis_records", list(reversed(data["basis_records"]))),
            ("quotient_coordinate_rows", [[1, 1]]),
            ("chart_universe_id", "example:changed-basis"),
        ):
            changed = deepcopy(data)
            changed[field] = replacement
            with self.subTest(field=field), self.assertRaisesRegex(
                ValueError, "stale chart ordering"
            ):
                validate_ordering(
                    certificate,
                    **{
                        key: value
                        for key, value in changed.items()
                        if key != "heuristics"
                    },
                )

    def test_changed_chart_order_requires_recomputation(self) -> None:
        data = fixture()
        certificate = bind_ordering(**data)
        changed = deepcopy(data)
        changed["ordered_chart_ids"] = list(reversed(data["ordered_chart_ids"]))
        with self.assertRaisesRegex(ValueError, "ordered chart identities changed"):
            validate_ordering(
                certificate,
                **{
                    key: value
                    for key, value in changed.items()
                    if key != "heuristics"
                },
            )


if __name__ == "__main__":
    unittest.main()
