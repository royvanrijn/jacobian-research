#!/usr/bin/env python3
"""Run the effective cancellation contact-resultant certificate at r=7."""

from __future__ import annotations

import os
import runpy
from pathlib import Path


os.environ["CONTACT_R"] = "7"
runpy.run_path(
    Path(__file__).with_name("verify_contact_resultant_r6_effective.py"),
    run_name="__main__",
)
