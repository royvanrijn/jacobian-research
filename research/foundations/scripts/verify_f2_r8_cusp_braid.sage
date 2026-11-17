#!/usr/bin/env sage-python
"""Certify the four algebraic ``r=8`` cusp-pair braid obstructions.

The exact cusp locus is the quartic field derived by
``verify_f2_geometric_degree_six_stein_reduction.py``.  Every complex
embedding is checked because braid type is not automatically preserved by an
arbitrary Galois embedding.  SageMath and its optional SIROCCO package are
required.
"""

import argparse
import importlib.machinery
import importlib.util
from pathlib import Path

from sage.all import NumberField, PolynomialRing, QQ, QQbar
from sage.schemes.curves.zariski_vankampen import braid_monodromy


def load_helpers():
    path = Path(__file__).with_name("verify_f2_r6_cusp_braid.sage")
    loader = importlib.machinery.SourceFileLoader("f2_braid_helpers", str(path))
    spec = importlib.util.spec_from_loader(loader.name, loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def coefficient_row(field):
    a = field.gen()
    c = -(
        51940000 * a**3
        + 65374350 * a**2
        + 15840099 * a
        - 898128
    ) / 4361202
    d = -(
        53410000 * a**3
        + 72989725 * a**2
        + 21185232 * a
        - 301824
    ) / 7268670
    x = -(
        7840000 * a**3
        + 22668000 * a**2
        + 25088409 * a
        + 4472496
    ) / 5814936
    p0 = (
        42140000 * a**3
        + 55884050 * a**2
        + 19349013 * a
        + 295344
    ) / 2422890
    A = a + x
    E = 5 * x / 3
    C = (4 * x + 3 * c + 5 * p0) / 3
    D = (5 * A**2 - 15 * A * a + 10 * a**2 + 12 * p0 + 9 * d) / 9
    q0 = (
        2 * A**2
        - 8 * A * a
        + 6 * A * c
        + 10 * A * p0
        + 6 * a**2
        - 6 * a * c
        - 15 * a * p0
    ) / 9
    return a, 1, c, d, A, E, C, D, p0, q0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--embedding",
        choices=("all", "0", "1", "2", "3"),
        default="all",
    )
    arguments = parser.parse_args()

    polynomial_ring = PolynomialRing(QQ, names=("a",))
    variable = polynomial_ring.gen()
    quartic = (
        196000000 * variable**4
        + 260940000 * variable**3
        + 82362825 * variable**2
        - 2390688 * variable
        + 20736
    )
    assert len(list(quartic.factor())) == 1
    roots = sorted(
        quartic.roots(QQbar, multiplicities=False),
        key=lambda root: (root.real(), root.imag()),
    )
    assert len(roots) == 4

    helpers = load_helpers()
    for index, root in enumerate(roots):
        if arguments.embedding != "all" and int(arguments.embedding) != index:
            continue
        field = NumberField(quartic, "a", embedding=root)
        first, second = helpers.implicit_pair(coefficient_row(field), field)
        braids, component_by_strand, vertical, degree = braid_monodromy(
            first * second,
            arrangement=(first, second),
        )
        assert degree == 6
        assert len(braids) == 12
        assert vertical == {}
        assert sorted(component_by_strand.values()) == [0, 0, 0, 1, 1, 1]
        words = [tuple(braid.Tietze()) for braid in braids]
        helpers.transposition_audit(
            words,
            component_by_strand,
            f"r=8 embedding {index}",
            expect_equal_disjoint=True,
        )


if __name__ == "__main__":
    main()
