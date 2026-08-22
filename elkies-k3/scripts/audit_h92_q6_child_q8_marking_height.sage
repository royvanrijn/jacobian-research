#!/usr/bin/env sage -python
"""Regression entry point for the repaired H92 q6-child q8 marking.

The historical audit exposed the withdrawn height-96 / degree-46 point.  The
canonical regression now runs the corrected 2-cover marking certificate,
which simultaneously verifies:

  * withdrawn point = 2 * corrected point;
  * corrected MW coordinate (-2,-2,0) has height 24;
  * corrected O-intersection/collision degree is 10;
  * II* and IV* specializations are on the identity component.
"""

from pathlib import Path

TARGET = Path(__file__).with_name("derive_h92_q6_child_q8_marking_2cover.sage")
exec(compile(TARGET.read_text(), str(TARGET), "exec"))
