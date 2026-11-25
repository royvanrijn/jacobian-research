#!/usr/bin/env python3
"""Exact retained complements for the terminal-thirteen lower bridge."""

from fractions import Fraction as Q

import verify_multivariate_terminal_thirteen_retained_fallbacks as replay


# Outer endpoints, lower and upper caps, and expected retained lists.
BANDS = {
    "63-64": (Q(63, 100), Q(16, 25), Q(1, 40), Q(1, 30), 13),
    "64-66": (Q(16, 25), Q(33, 50), Q(1, 40), Q(1, 30), 13),
    "63-66-tight": (Q(63, 100), Q(33, 50), Q(1, 50), Q(1, 40), 38),
    "68-72": (Q(17, 25), Q(18, 25), Q(1, 50), Q(1, 36), 32),
    "72-76": (Q(18, 25), Q(19, 25), Q(1, 50), Q(1, 34), 16),
    "68-72-tight": (Q(17, 25), Q(18, 25), Q(1, 60), Q(1, 50), 7),
    "72-76-tight": (Q(18, 25), Q(19, 25), Q(1, 60), Q(1, 50), 3),
}


if __name__ == "__main__":
    replay.BANDS = BANDS
    replay.main()
