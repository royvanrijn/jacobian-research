#!/usr/bin/env sage -python
"""Compatibility entry point for the corrected H92 q6-child q8 marking.

The previous implementation double-counted the binary-quartic 2-covering map
and produced a section of height 96 / collision degree 46.  The canonical
marking is now derived by ``derive_h92_q6_child_q8_marking_2cover.sage`` and
has height 24 / collision degree 10.
"""

from pathlib import Path

TARGET = Path(__file__).with_name("derive_h92_q6_child_q8_marking_2cover.sage")
exec(compile(TARGET.read_text(), str(TARGET), "exec"))
