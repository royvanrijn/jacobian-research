#!/usr/bin/env sage -python
"""Git-only exact verifier for v2/v3 Schur-complement rank-gain hits.

This reuses the baseline-first verifier from ``batch_verify_v2_rank_gain_hits``
but replaces its legacy ``.sobj`` specialization builder with the committed
``hidden_sections_data.py`` representation.  Child verifier processes recurse
through this shim, so no ``/tmp/newfamily_hidden_sections_complete.sobj`` is
required.
"""

from __future__ import annotations

from pathlib import Path

import batch_verify_v2_rank_gain_hits as v2
from search_unseeded_extra_points_v3 import build_integral_specialization_git


# The imported verifier resolves this symbol dynamically in verify_specialization.
v2.build_integral_specialization = build_integral_specialization_git

# Its parent launches isolated children via Path(__file__). Point that module
# global at this shim so children retain the git-only specialization builder.
v2.__file__ = str(Path(__file__).resolve())


if __name__ == "__main__":
    raise SystemExit(v2.main())
