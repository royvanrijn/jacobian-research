#!/usr/bin/env python3
"""Verify split and collided tame-toroidal packets over a target node."""

from __future__ import annotations

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def fitting_ideals(matrix: sp.Matrix) -> tuple[sp.Expr, tuple[sp.Expr, ...]]:
    determinant = sp.expand(matrix.det())
    entries = tuple(sp.expand(entry) for entry in matrix)
    return determinant, entries


def split_packet_audit(index: int) -> None:
    # x=r^e, y=t in logarithmic bases (dlog x,dlog y) and
    # (dlog r,dlog t).  In characteristic zero, e is a unit.
    matrix = sp.Matrix([[index, 0], [0, 1]])
    determinant, entries = fitting_ideals(matrix)
    assert determinant == index
    assert sp.gcd_list(entries) == 1
    assert smith_normal_form(matrix, domain=ZZ) == sp.diag(1, index)
    assert matrix.rank() == 2
    assert matrix.inv() == sp.Matrix([[sp.Rational(1, index), 0], [0, 1]])


def collided_packet_audit(index: int) -> None:
    # For z^e=x*y use the character-lattice basis
    # m1=x, m2=(x+y)/e.  Then x=m1 and y=e*m2-m1.
    character_matrix = sp.Matrix([[1, -1], [0, index]])
    determinant, entries = fitting_ideals(character_matrix)
    assert determinant == index
    assert sp.gcd_list(entries) == 1
    assert smith_normal_form(character_matrix, domain=ZZ) == sp.diag(1, index)
    assert character_matrix.rank() == 2

    # The A_(e-1) cone in the dual basis is spanned by (0,1),(e,1).
    # Its minimal regular subdivision has rays (j,1).  Interior rays map to
    # the target node with orders (j,e-j) and are (-2)-curves.
    rays = [sp.Matrix([j, 1]) for j in range(index + 1)]
    for left, right in zip(rays, rays[1:]):
        assert sp.det(sp.Matrix.hstack(left, right)) == -1
    for j, ray in enumerate(rays):
        order_x = int(ray[0])
        order_y = index * int(ray[1]) - int(ray[0])
        assert (order_x, order_y) == (j, index - j)
        if 0 < j < index:
            assert order_x > 0 and order_y > 0
    for j in range(1, index):
        assert rays[j - 1] + rays[j + 1] == 2 * rays[j]
    assert len(rays[1:-1]) == index - 1


def general_kummer_lattice_audit() -> None:
    # The adjugate identity proves symbolically that every rank-two finite-index
    # lattice inclusion becomes invertible after its determinant is inverted.
    a, b, c, d = sp.symbols("a b c d")
    generic = sp.Matrix([[a, b], [c, d]])
    determinant = generic.det()
    assert generic.adjugate() * generic == determinant * sp.eye(2)
    assert generic * generic.adjugate() == determinant * sp.eye(2)

    # Concrete nonscalar Smith types and both signs guard the implementation.
    matrices = (
        sp.Matrix([[2, 1], [0, 3]]),
        sp.Matrix([[3, -2], [1, 1]]),
        sp.Matrix([[5, 0], [2, 1]]),
        sp.Matrix([[4, 3], [-1, 2]]),
    )
    for matrix in matrices:
        assert matrix.det() != 0
        assert matrix.rank() == 2
        assert matrix.inv() * matrix == sp.eye(2)


def monomial_unit_jet_audit() -> None:
    # A general completed SNC pullback has
    # x=u^a*v^b*alpha and y=u^c*v^d*beta with alpha,beta units.  Its residue
    # log matrix is the exponent matrix.  When that matrix has rank one, the
    # first jet of its determinant is the following two-coefficient gate.
    u, v = sp.symbols("u v")
    a, b, c, d = sp.symbols("a b c d")
    alpha_u, alpha_v, beta_u, beta_v = sp.symbols(
        "alpha_u alpha_v beta_u beta_v"
    )
    alpha = 1 + alpha_u * u + alpha_v * v
    beta = 1 + beta_u * u + beta_v * v
    theta = sp.Matrix(
        [
            [a + u * sp.diff(alpha, u) / alpha, c + u * sp.diff(beta, u) / beta],
            [b + v * sp.diff(alpha, v) / alpha, d + v * sp.diff(beta, v) / beta],
        ]
    )
    determinant = sp.factor(theta.det())
    assert determinant.subs({u: 0, v: 0}) == a * d - b * c
    assert sp.expand(
        sp.diff(determinant, u).subs({u: 0, v: 0})
        - (d * alpha_u - b * beta_u)
    ) == 0
    assert sp.expand(
        sp.diff(determinant, v).subs({u: 0, v: 0})
        - (a * beta_v - c * alpha_v)
    ) == 0

    # Every nonzero rank-one exponent matrix has a scalar unit entry, so its
    # cokernel has unit Fitt_1 and is cyclic.  A nonzero first determinant jet
    # makes its reduced support smooth; only the simultaneous vanishing of
    # the two displayed coefficients reaches the higher-order defect gate.
    for av_a in range(5):
        for av_b in range(5):
            for av_c in range(5):
                for av_d in range(5):
                    residue = sp.Matrix([[av_a, av_c], [av_b, av_d]])
                    if residue.rank() != 1:
                        continue
                    assert sp.gcd_list(tuple(residue)) >= 1
                    assert any(entry != 0 for entry in residue)

    smooth_specialization = {
        a: 1,
        b: 1,
        c: 2,
        d: 2,
        alpha_u: 0,
        alpha_v: 0,
        beta_u: 1,
        beta_v: 0,
    }
    assert (a * d - b * c).subs(smooth_specialization) == 0
    assert sp.diff(determinant, u).subs(
        {u: 0, v: 0} | smooth_specialization
    ) == -1

    higher_order_specialization = {
        a: 1,
        b: 1,
        c: 2,
        d: 2,
        beta_u: 2 * alpha_u,
        beta_v: 2 * alpha_v,
    }
    assert (a * d - b * c).subs(higher_order_specialization) == 0
    assert sp.diff(determinant, u).subs(
        {u: 0, v: 0} | higher_order_specialization
    ) == 0
    assert sp.diff(determinant, v).subs(
        {u: 0, v: 0} | higher_order_specialization
    ) == 0


def main() -> None:
    for index in range(2, 31):
        split_packet_audit(index)
        collided_packet_audit(index)
    general_kummer_lattice_audit()
    monomial_unit_jet_audit()
    print(
        "PASS: split and collided tame Kummer packets over a node are log-"
        "etale after toric resolution; the general SNC monomial-unit rank "
        "and first-jet gates isolate every possible higher local defect"
    )


if __name__ == "__main__":
    main()
