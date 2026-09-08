#!/usr/bin/env sage -python
"""Git-only entry point for robust newfamily extra-point discovery.

This reuses the calibrated v2 search logic but replaces the old ``.sobj``
section source with the committed ``hidden_sections_data.py`` representation.
It therefore requires no ``/tmp/newfamily_hidden_sections_complete.sobj``.

The v2 parent launches isolated child processes via ``Path(__file__)``.  Since
that lookup happens in the imported v2 module's globals, this shim also points
``v2.__file__`` at itself so child processes re-enter v3 and retain the git-only
section loader.
"""

from __future__ import annotations

from pathlib import Path

from sage.all import EllipticCurve, QQ, ZZ

import search_unseeded_extra_points_v2 as v2
from screen_seeded_rational_candidates_fast import ROOTS, extra_integral_scale, load_builder
from hidden_sections import load_hidden_sections


def build_integral_specialization_git(numerator: int, denominator: int, _legacy_source=None):
    a = ZZ(numerator)
    b = ZZ(denominator)
    t = QQ(a) / QQ(b)
    sections = load_hidden_sections()
    family = load_builder()(ROOTS)

    A = QQ(family["Amin"](t))
    B = QQ(family["Bmin"](t))
    Ah = A * b**8
    Bh = B * b**12
    c = extra_integral_scale(Ah, Bh)
    Aint_q = Ah * c**4
    Bint_q = Bh * c**6
    if Aint_q.denominator() != 1 or Bint_q.denominator() != 1:
        raise RuntimeError("failed to integralize homogeneous short model")

    E = EllipticCurve(QQ, [0, 0, 0, ZZ(Aint_q), ZZ(Bint_q)])
    if E.discriminant() == 0:
        raise RuntimeError("singular specialization")

    xscale = b**4 * c**2
    yscale = b**6 * c**3
    known = [
        E([QQ(xf(t)) * xscale, QQ(yf(t)) * yscale])
        for xf, yf in sections
    ]
    if len(set(known)) != 11:
        raise RuntimeError("hidden sections collide at specialization")
    return E, known, c


# v2 resolves this symbol dynamically from its module globals inside run_single.
v2.build_integral_specialization = build_integral_specialization_git

# v2.run_parent launches isolated children using Path(__file__).  Point that
# module-global name at this shim so the child also uses the git-only loader.
v2.__file__ = str(Path(__file__).resolve())


if __name__ == "__main__":
    raise SystemExit(v2.main())
