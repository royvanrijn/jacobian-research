#!/usr/bin/env python3
"""Compatibility entry point for the disproved cubic-survivor candidate.

The exact calculation now lives in
``verify_two_pair_sic_bidegree33_Q_residual_slice.py`` and proves that the
candidate is supported on the excluded ``u=0`` boundary.
"""

from verify_two_pair_sic_bidegree33_Q_residual_slice import main


if __name__ == "__main__":
    main()
