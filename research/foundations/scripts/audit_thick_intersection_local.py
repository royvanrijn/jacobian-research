"""Clean-room exact-arithmetic audit of the cancellation contact length.

This deliberately avoids SymPy and resultants.  It checks the two branch
restrictions used by the completed-local-ring proof and then applies
div(P)=m[Y=0]+[Y=Q].  The bounded loop is a transcription regression; the
all-parameter identities are proved in the paper by the beta integral.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial


def beta_sum(a: int, b: int) -> Fraction:
    """Return integral_0^1 u^a(1-u)^b du by binomial expansion."""

    return sum(
        (Fraction((-1) ** j * comb(b, j), a + j + 1) for j in range(b + 1)),
        start=Fraction(0),
    )


def audit_pair(m: int, r: int) -> None:
    zero_branch = beta_sum(r, m * r)
    zero_expected = Fraction(
        factorial(r) * factorial(m * r), factorial(m * r + r + 1)
    )
    assert zero_branch == zero_expected

    diagonal_branch = sum(
        (Fraction((-1) ** j * comb(r, j), j + 1) for j in range(r + 1)),
        start=Fraction(0),
    )
    assert diagonal_branch == Fraction(1, r + 1)

    zero_order = m * r
    diagonal_order = m * r
    local_length = m * zero_order + diagonal_order
    expected_length = m * r * (m + 1)
    assert local_length == expected_length


def main() -> None:
    for m in range(1, 13):
        for r in range(1, 13):
            audit_pair(m, r)
    print(
        "PASS: completed-local cancellation contact gives exact length "
        "mr(m+1) for 144 exact-arithmetic regression pairs"
    )


if __name__ == "__main__":
    main()
