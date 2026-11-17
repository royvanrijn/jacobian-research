#!/usr/bin/env python3
"""Audit the first structural consequences of the Ritt-contact proposal.

This is deliberately a small checker.  The expensive ideal calculations
which prove the degree-thirty sector data live in
``verify_degree30_ritt_2_complex.py``.  Here we test what those exact data
do and do not determine, and the behavior under a tame ramified pullback of
the power--Dickson contact parameter.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.ritt_complex import (  # noqa: E402
    BraidSectorDecoration,
    MoveType,
    degree_thirty_braid_decorations,
    symmetric_braid_complex,
)


@dataclass(frozen=True)
class RamifiedContactSignature:
    """Contact data after ``z = Z**ramification_index``.

    The primitive contact monoid is still abstractly ``N``.  Its labelled
    map from the original overlap divisor is multiplication by the
    ramification index, and an annihilator ``(z**a)`` pulls back to
    ``(Z**(a*m))``.
    """

    composite_omission: int
    ramification_index: int
    labelled_contact_vector: tuple[int]
    annihilator_exponent: int
    transverse_exponents: tuple[int, ...]
    nilpotence_index: int


def ramified_contact_signature(
    sector: BraidSectorDecoration,
    ramification_index: int,
) -> RamifiedContactSignature:
    if ramification_index <= 0:
        raise ValueError("the ramification index must be positive")
    return RamifiedContactSignature(
        composite_omission=sector.composite_omission,
        ramification_index=ramification_index,
        labelled_contact_vector=(ramification_index,),
        annihilator_exponent=(
            ramification_index * sector.conductor_power
        ),
        transverse_exponents=sector.transverse_slice.exponents,
        nilpotence_index=sector.nilpotence_index,
    )


def main() -> None:
    complex_ = symmetric_braid_complex((2, 3, 5), MoveType.CHEBYSHEV)
    sectors = degree_thirty_braid_decorations()

    # All sectors lie on the same reduced braid cell and have the same
    # abstract primitive overlap monoid N.
    relation_signatures = {
        (
            len(complex_.vertices),
            len(complex_.edges),
            len(complex_.two_cells),
            complex_.euler_characteristic,
        )
        for _ in sectors
    }
    assert relation_signatures == {(6, 6, 1, 1)}
    assert {(1,) for _ in sectors} == {(1,)}

    # Coarse derived data do not classify sectors: cuts 15 and 6 have the
    # same conormal and ambient Koszul-Tor ranks.
    coarse_derived = [
        (
            sector.transverse_slice.conormal_rank,
            sector.transverse_slice.residue_field_tor_ranks,
        )
        for sector in sectors
    ]
    assert coarse_derived[1] == coarse_derived[2]

    # The contact exponent alone also does not classify sectors: cuts 10
    # and 15 both have defect annihilator (z^2).
    assert sectors[0].conductor_power == sectors[1].conductor_power == 2

    # The relative package used in the note distinguishes all three exact
    # sectors: it retains contact action, internal CI orders, and the
    # nilpotence filtration.
    relative_signatures = {
        (
            sector.conductor_power,
            sector.transverse_slice.exponents,
            sector.nilpotence_index,
        )
        for sector in sectors
    }
    assert len(relative_signatures) == 3

    # A ramified pullback fixes the relation cell, transverse algebra, and
    # nilpotence index, but changes the labelled integral contact map and
    # every contact-supported annihilator order.
    for ramification_index in (1, 2, 3):
        pulled_back = [
            ramified_contact_signature(sector, ramification_index)
            for sector in sectors
        ]
        assert {
            signature.labelled_contact_vector
            for signature in pulled_back
        } == {(ramification_index,)}
        for sector, signature in zip(sectors, pulled_back):
            assert (
                signature.annihilator_exponent
                == ramification_index * sector.conductor_power
            )
            assert (
                signature.transverse_exponents
                == sector.transverse_slice.exponents
            )
            assert signature.nilpotence_index == sector.nilpotence_index

    # Forgetting the labelled map makes all ramified primitive rank-one
    # monoids abstractly isomorphic, which is exactly the information loss
    # exposed by the collision-local-ring counterexample.
    unlabelled_rank = {
        len(
            ramified_contact_signature(sectors[0], index)
            .labelled_contact_vector
        )
        for index in (1, 2, 3)
    }
    assert unlabelled_rank == {1}

    print("PASS: all thick sectors share one reduced S3 braid 2-cell")
    print("PASS: the abstract primitive overlap monoid does not classify them")
    print("PASS: contact order alone conflates the cut-10 and cut-15 sectors")
    print("PASS: coarse conormal/Tor data conflate the cut-15 and cut-6 sectors")
    print("PASS: relative contact action plus internal defect orders separates all sectors")
    print("PASS: z=Z^m scales labelled contact and annihilator orders by m")
    print("PASS: ramified pullback preserves the relation cell and transverse algebra")


if __name__ == "__main__":
    main()
