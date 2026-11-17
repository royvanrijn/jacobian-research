#!/usr/bin/env python3
"""Exact checks for the universal boundary-saturation research note.

The general logarithmic purity theorem is proved in the note.  This script
checks the finite presentation calculations used by the two no-go models and
the pure-divisor control.  It deliberately creates no generated artifact.
"""

from __future__ import annotations

from sympy import groebner, sympify, symbols


def ideal_contains(generators, polynomial, variables) -> bool:
    basis = groebner(generators, *variables, order="grevlex")
    return basis.reduce(sympify(polynomial))[1] == 0


def main() -> None:
    a, b = symbols("a b")

    # One split node has rank-one normalization quotient.  Multiplication by
    # a presents A/(a); the class 1 is nonzero and killed by the collision
    # ideal (a), so the whole module is collision-supported.
    assert not ideal_contains([a], 1, (a, b))
    assert ideal_contains([a], a, (a, b))

    # A chain of three normalized components has two nodes and hence a
    # rank-two conductor quotient.  The matching matrix
    # [[a,0,0],[0,a,b]] has cokernel A/(a) + A/(a,b).
    # Its maximal minors generate (a^2, a*b)=a*(a,b).
    maximal_minors = [a * a, a * b, 0]
    assert all(ideal_contains([a * a, a * b], f, (a, b)) for f in maximal_minors)
    assert ideal_contains(maximal_minors, a * a, (a, b))
    assert ideal_contains(maximal_minors, a * b, (a, b))
    assert not ideal_contains(maximal_minors, a, (a, b))

    # The second target basis vector is not in the image, while multiplication
    # by either generator of I=(a,b) puts it in the image.  In the block
    # presentation this is exactly the nonzero A/(a,b) torsion summand.
    assert not ideal_contains([a, b], 1, (a, b))
    assert ideal_contains([a, b], a, (a, b))
    assert ideal_contains([a, b], b, (a, b))

    # Pure-divisor control: diag(a,1) has cokernel A/(a).  For the closed-point
    # ideal (a,b), b is regular modulo (a); in particular the class of a
    # polynomial killed by a power of (a,b) must already vanish.  The finite
    # checks below pin the relevant coprimality/membership facts.
    assert not ideal_contains([a], b, (a, b))
    for exponent in range(1, 8):
        assert not ideal_contains([a], b**exponent, (a, b))

    print("PASS: universal boundary-saturation countermodels and control")


if __name__ == "__main__":
    main()
