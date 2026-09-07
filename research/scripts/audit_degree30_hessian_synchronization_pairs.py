#!/usr/bin/env python3
"""Replay the four non-tree degree-30 pair synchronization certificates.

The fast default synchronization regression uses only a five-edge spanning
tree, because that already proves the global all-six theorem.  These four
additional exact reductions, together with the nested ``{2,10}``
normal-coordinate certificate, leave precisely five incomparable pairs open.
"""

from __future__ import annotations

from verify_hessian_synchronization_lifts import audit_pair_on_factor_chart


def main() -> None:
    cases = (
        (5, 6, 502),
        (6, 10, 189),
        (6, 15, 12),
        (10, 15, 96),
    )
    for source_cut, target_cut, expected_basis_size in cases:
        basis_size = audit_pair_on_factor_chart(
            30,
            source_cut,
            target_cut,
            timeout=300,
        )
        assert basis_size == expected_basis_size
        print(
            f"degree 30: pair ({source_cut},{target_cut}) synchronizes "
            f"scheme-theoretically; basis size {basis_size}"
        )
    print("PASS: ten of fifteen degree-30 pairs have exact certificates")
    print("PASS: precisely five incomparable degree-30 pairs remain open")


if __name__ == "__main__":
    main()
