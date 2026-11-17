#!/usr/bin/env python3
"""Replay the complete four-cubic finite-field obstruction."""

from __future__ import annotations

from pathlib import Path
import runpy


SCRIPT_DIRECTORY = Path(__file__).parent

runpy.run_path(
    str(SCRIPT_DIRECTORY / "verify_hc4_meng_four_cubic_rank_gate.py")
)
runpy.run_path(
    str(SCRIPT_DIRECTORY / "verify_hc4_meng_four_cubic_rank_zero.py")
)

print(
    "PASS: cubic support four is excluded over the certificate field "
    "F_1000003"
)
