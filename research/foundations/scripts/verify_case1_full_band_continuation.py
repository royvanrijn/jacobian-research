#!/usr/bin/env python3
"""Verify the pinned Case-1 full-band continuation ledger.

The complete exact replay is intentionally a specialized long calculation.
This verifier checks its immutable source hashes, all deterministic ledger
invariants, and replays the first omitted layer exactly against the pinned
hashes.  Use ``continue_case1_full_bands.py --stop-layer -11`` for the full
eight-layer characteristic-zero replay.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from continue_case1_full_bands import (  # noqa: E402
    continue_descent,
    deterministic_ledger,
    support_p,
    support_q,
)


LEDGER = (
    ROOT
    / "artifacts/generated-results/case1_full_band_continuation.json"
)
LEDGER_SHA256 = (
    "75e24fc937c162b877336e4bac3ac2c2"
    "ee5f3c432c55d6ccf48d3cc3f7827f60"
)


assert hashlib.sha256(LEDGER.read_bytes()).hexdigest() == LEDGER_SHA256
ledger = json.loads(LEDGER.read_text())
assert ledger["schema_version"] == 1
assert ledger["initial_layer"] == -4
assert ledger["terminal_layer"] == -11
assert ledger["initial_parameter_count"] == 6
assert ledger["final_parameter_count"] == 6
assert ledger["derived_P_bands"] == [-8, -7, -6]
assert ledger["derived_Q_bands"] == list(range(-12, -4))

layers = ledger["layers"]
assert [item["layer"] for item in layers] == list(range(-4, -12, -1))
assert [item["column_count"] for item in layers] == [11, 9, 7, 5, 4, 3, 2, 1]
assert [item["compatibility_count"] for item in layers] == [6, 7, 8, 9, 9, 9, 9, 9]
assert sum(item["compatibility_count"] for item in layers) == 66

for item in layers:
    layer = item["layer"]
    expected_p_support = support_p(layer - 2)
    expected_q_support = support_q(layer - 1)
    assert item["P_support"] == expected_p_support
    assert item["Q_support"] == expected_q_support
    assert item["column_count"] == len(expected_p_support) + len(expected_q_support)
    assert item["rank"] == item["column_count"]
    assert item["nullity"] == 0
    assert item["parameters_before"] == item["parameters_after"] == 6
    assert item["compatibility_shape"]["coefficient_slots"] == item[
        "compatibility_count"
    ]
    if item["P_shape"] is not None:
        assert item["P_shape"]["coefficient_slots"] == len(expected_p_support)
        assert item["P_shape"]["nonzero_coefficients"] == len(expected_p_support)
    if item["Q_shape"] is not None:
        assert item["Q_shape"]["coefficient_slots"] == len(expected_q_support)
        assert item["Q_shape"]["nonzero_coefficients"] == len(expected_q_support)

first_report, _, _ = continue_descent(-4)
first_ledger = deterministic_ledger(first_report)
for field in (
    "checkpoint_sha256",
    "exact_core_sha256",
    "initial_layer",
    "initial_parameter_count",
    "final_parameter_count",
):
    assert first_ledger[field] == ledger[field]
assert first_ledger["layers"] == layers[:1]

print("PASS: pinned source and full Case-1 continuation ledger are consistent")
print("PASS: all eight lower-band systems have full column rank and no new moduli")
print("PASS: first omitted layer replays exactly against its band and compatibility hashes")
