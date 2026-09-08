#!/usr/bin/env sage
"""Certify the class-group images of curve 273's residual endpoint ideals.

The full-ideal chain terminates at five degree-one ideals in the cubic
2-division field ``Q(theta)``.  This script moves to the reduced defining
polynomial returned by Sage's ``optimized_representation`` and computes the
class group with ``proof=True``.  It reports the exact class coordinates and
their images modulo squares; these are the finite obstruction that the
large-prime relation collection cannot see.

No Selmer or Mordell--Weil rank assertion is made here.  In particular, a
zero endpoint class modulo squares only says that the global ideal-class
obstruction has been removed; the 2-cover local conditions remain separate.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
import time

from sage.all import NumberField, PolynomialRing, QQ, ZZ, proof


sys.path.insert(0, str(Path(__file__).resolve().parent))

from curve273_full_ideal_chain import T9
from icarm_curve273 import short_coefficients
from run_curve273_mod2_relations import BAD_RATIONAL_PRIMES


PROTOCOL = "R30ENDCLASS"


def sage_q(value):
    numerator = value.numerator
    denominator = value.denominator
    if callable(numerator):
        numerator = numerator()
    if callable(denominator):
        denominator = denominator()
    return QQ(ZZ(numerator)) / QQ(ZZ(denominator))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--s-unit-basis",
        action="store_true",
        help=(
            "construct the global K(S,2) envelope for all primes above the "
            "2-division discriminant support"
        ),
    )
    args = parser.parse_args()
    proof.all(True)
    coefficients = short_coefficients()
    A = ZZ(sage_q(coefficients[3]))
    B = ZZ(sage_q(coefficients[4]))
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    original = NumberField(x**3 + A * x + B, "theta")
    reduced, reduced_to_original, original_to_reduced = (
        original.optimized_representation()
    )
    alpha = reduced.gen()

    # For this field Sage gives theta |--> -9*alpha + 3.  Keep the map
    # symbolic rather than relying on a positional interpretation of the
    # optimized-representation tuple.
    theta_image = original_to_reduced(original.gen())
    print(
        f"{PROTOCOL}|stage=input|original_disc_bits="
        f"{abs(ZZ(original.discriminant())).nbits()}"
        f"|reduced_polynomial={reduced.defining_polynomial()}"
        f"|theta_image={theta_image}",
        flush=True,
    )

    started = time.monotonic()
    classes = reduced.class_group(proof=True)
    orders = tuple(ZZ(order) for order in classes.gens_orders())
    print(
        f"{PROTOCOL}|stage=class_group|order={classes.order()}"
        f"|invariants={','.join(map(str, orders))}"
        f"|seconds={time.monotonic()-started:.3f}",
        flush=True,
    )

    # A degree-one ideal (q, theta-r) transports to
    # (q, theta_image-r).  Its residue is checked against the reduced
    # defining polynomial before querying its class.
    for index, (q, residue) in enumerate(T9, 1):
        q = ZZ(q)
        residue = ZZ(residue)
        transported = reduced.ideal(q, theta_image - residue)
        assert transported.is_prime() and transported.norm() == q
        coordinates = tuple(ZZ(value) for value in classes(transported).exponents())
        parity = tuple(
            coordinate % 2
            for coordinate, order in zip(coordinates, orders)
            if order % 2 == 0
        )
        print(
            f"{PROTOCOL}|endpoint={index}|q={q}|q_bits={q.nbits()}"
            f"|theta_residue={residue}|class_coordinates="
            f"{','.join(map(str, coordinates)) or 'empty'}"
            f"|class_mod_2={','.join(map(str, parity)) or 'empty'}",
            flush=True,
        )

    if args.s_unit_basis:
        # For an irreducible cubic 2-division polynomial, a global 2-descent
        # squareclass is unramified away from the primes above 2*disc(f).
        # BAD_RATIONAL_PRIMES is the exact rational support reconstructed from
        # that discriminant in the independent mod-2 relation audit.
        S = tuple(
            prime
            for rational_prime in BAD_RATIONAL_PRIMES
            for prime in reduced.primes_above(ZZ(rational_prime))
        )
        assert len(S) == len(set(S))
        space, generators, _, _ = reduced.selmer_space(S, 2, proof=True)
        print(
            f"{PROTOCOL}|stage=s_unit_envelope|S_size={len(S)}"
            f"|dimension={space.dimension()}",
            flush=True,
        )
        for index, generator in enumerate(generators, 1):
            print(
                f"{PROTOCOL}|s_unit_generator={index}|value={generator}",
                flush=True,
            )

    print(f"{PROTOCOL}|status=PASS", flush=True)


if __name__ == "__main__":
    main()
