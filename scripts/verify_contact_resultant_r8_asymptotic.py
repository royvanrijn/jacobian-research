#!/usr/bin/env python3
"""Run the exact endpoint and branch-at-infinity certificate at r=8."""

from __future__ import annotations

import os
import runpy
from pathlib import Path


os.environ["CONTACT_R"] = "8"
runpy.run_path(
    Path(__file__).with_name("verify_contact_resultant_r7_asymptotic.py"),
    run_name="__main__",
)
