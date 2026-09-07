"""Canonical loader for the eleven exact hidden generic newfamily sections.

The source of truth is ``hidden_sections_data.py`` committed next to this file.
The old Sage ``.sobj`` reconstruction is intentionally not required.
"""

from __future__ import annotations

from sage.all import FractionField, PolynomialRing, QQ

from hidden_sections_data import SECTIONS


def _q(value):
    if isinstance(value, tuple):
        return QQ(value[0]) / QQ(value[1])
    return QQ(value)


def _poly(ring, coefficients):
    return ring([_q(value) for value in coefficients])


def load_hidden_sections():
    """Return ``[(x_0(T),y_0(T)),...,(x_10(T),y_10(T))]`` over ``QQ(T)``."""
    if len(SECTIONS) != 11:
        raise RuntimeError(f"expected 11 hidden section records, found {len(SECTIONS)}")

    ring = PolynomialRing(QQ, "T")
    field = FractionField(ring)
    result = []
    for expected, record in enumerate(SECTIONS):
        if record["index"] != expected:
            raise RuntimeError(
                f"hidden sections out of order: expected U{expected}, got U{record['index']}"
            )
        xn = _poly(ring, record["x_num"])
        xd = _poly(ring, record["x_den"])
        yn = _poly(ring, record["y_num"])
        yd = _poly(ring, record["y_den"])
        result.append((field(xn) / field(xd), field(yn) / field(yd)))
    return result
