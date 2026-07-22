#!/usr/bin/env python3
"""Regression tests for the superelliptic de Rham kernel."""

import importlib.util
from pathlib import Path

import sympy as sp


MODULE = Path(__file__).with_name("superelliptic_derham.py")
SPEC = importlib.util.spec_from_file_location("superelliptic_derham", MODULE)
assert SPEC and SPEC.loader
ENGINE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ENGINE)

SuperellipticDeRham = ENGINE.SuperellipticDeRham
weighted_wronskian_compatibility = ENGINE.weighted_wronskian_compatibility

t = sp.symbols("t")

# The (72,108) curve has genus three and a six-dimensional compact character.
a2, a3, a4, a5, a6, a7 = sp.symbols("a2:8")
A72 = t + a2 * t**2 + a3 * t**3 + a4 * t**4 + a5 * t**5 + a6 * t**6 + a7 * t**7 + t**8
# Generic squarefreeness is automatic on the Wronskian solution locus and is
# deliberately not recomputed over QQ(a2,...,a7) in this fast regression.
curve72 = SuperellipticDeRham(t, A72, 2, check_squarefree=False)
assert curve72.genus == 3
assert curve72.affine_h1_dimension == 7
assert curve72.compact_h1_dimension == 6
assert curve72.character_dimension(3, compact=False) == 7
assert curve72.character_dimension(3, compact=True) == 6

# Reproduce the infinity relation and basis in WEIGHTED_WRONSKIAN_FIRST_BLOCK.
residues = curve72.residue_coefficients(3)
expected = (
    0,
    0,
    0,
    1,
    -a7 / 2,
    3 * a7**2 / 8 - a6 / 2,
    -a5 / 2 + 3 * a7 * a6 / 4 - 5 * a7**3 / 16,
)
assert all(sp.expand(x - y) == 0 for x, y in zip(residues, expected))
expected_basis = (
    1,
    t,
    t**2,
    t**4 + a7 * t**3 / 2,
    t**5 + (a6 / 2 - 3 * a7**2 / 8) * t**3,
    t**6 + (a5 / 2 - 3 * a7 * a6 / 4 + 5 * a7**3 / 16) * t**3,
)
assert all(
    sp.expand(x - y) == 0 for x, y in zip(curve72.compact_basis(3), expected_basis)
)

# The support-constrained block solves 11 D coefficients and leaves 6 equations.
dvars = sp.symbols("d2:13")
solved, compatibility = weighted_wronskian_compatibility(
    t, A72, t**2, 2, 3, dvars, range(2, 13)
)
assert len(solved) == 11
assert len(compatibility) == 6

# Hermite reduction recognizes an exact weighted Wronskian over QQ.
Aq = t**8 + 2 * t**7 - t**3 + t + 1
curveq = SuperellipticDeRham(t, Aq, 2)
Dq = t**12 - 3 * t**7 + 2 * t**2
Rq = sp.expand(2 * Aq * sp.diff(Dq, t) - 3 * sp.diff(Aq, t) * Dq)
reduced = curveq.reduce_weighted_wronskian(Rq, 3)
assert reduced.is_exact
assert reduced.remainder == 0

# A genuine superelliptic example: y^3=A of degree four has genus three.
A3 = t**4 + t + 1
curve3 = SuperellipticDeRham(t, A3, 3)
assert curve3.genus == 3
assert curve3.compact_h1_dimension == 6
assert curve3.character_dimension(1) == 3
assert curve3.character_dimension(2) == 3
D3 = t**7 + 2 * t + 1
R3 = sp.expand(3 * A3 * sp.diff(D3, t) - 2 * sp.diff(A3, t) * D3)
assert curve3.reduce_weighted_wronskian(R3, 2).is_exact

# The trivial character descends to an ordinary rational differential on the
# quotient instead of entering one of the nontrivial eigenspace reducers.
trivial = curve3.reduce_weighted_wronskian(1 + t, 3)
assert trivial.character == 0
assert sp.cancel(trivial.differential - (1 + t) / (3 * A3**2)) == 0

# When gcd(a,deg A)>1, residue directions occur in selected characters only.
A6 = t**6 + t + 1
curve6 = SuperellipticDeRham(t, A6, 4)
assert curve6.delta == 2
assert curve6.character_dimension(1) == 5
assert curve6.character_dimension(2) == 4
assert curve6.character_dimension(3) == 5
assert sum(curve6.character_dimension(r) for r in range(1, 4)) == 14
assert curve6.compact_h1_dimension == 14

# The Gauss--Manin layer differentiates in the parameter directions and then
# reuses Hermite reduction.  On the two-parameter elliptic family, the
# connection is flat and its only denominator is the discriminant divisor.
u, v = sp.symbols("u v")
elliptic_A = t**3 + u * t + v
elliptic = SuperellipticDeRham(t, elliptic_A, 2, check_squarefree=False)
connection = elliptic.gauss_manin_connection((u, v), 1)
assert connection.basis == (1, t)
Mu, Mv = connection.matrices
curvature = (
    Mv.diff(u) - Mu.diff(v) + Mu * Mv - Mv * Mu
).applyfunc(sp.factor)
assert curvature == sp.zeros(2)
discriminant = sp.factor(sp.discriminant(elliptic_A, t))
for matrix in connection.matrices:
    for value in matrix:
        denominator = sp.factor(sp.denom(sp.cancel(value)))
        assert sp.rem(sp.Poly(discriminant, u), sp.Poly(denominator, u)) == 0 or sp.rem(
            sp.Poly(denominator, u), sp.Poly(-discriminant, u)
        ) == 0

# The actual plane-block rank also works without specializing to genus one.
# This one-parameter degree-eight slice gives a dense 6x6 connection whose
# denominators are scalar multiples of the slice discriminant.
plane_u = sp.symbols("plane_u")
plane_A = t**8 + plane_u * t**7 + t
plane_curve = SuperellipticDeRham(t, plane_A, 2, check_squarefree=False)
plane_connection = plane_curve.gauss_manin_connection((plane_u,), 1)
plane_matrix = plane_connection.matrices[0]
assert plane_curve.genus == 3
assert plane_matrix.shape == (6, 6)
assert all(value != 0 for value in plane_matrix)
plane_discriminant = sp.factor(sp.discriminant(plane_A, t))
for value in plane_matrix:
    denominator = sp.factor(sp.denom(sp.cancel(value)))
    assert not sp.cancel(denominator / plane_discriminant).has(plane_u)

print("PASS: superelliptic Hermite/de Rham reduction")
print("PASS: (72,108) gives 11 solved coefficients and 6 compact obstructions")
print("PASS: character dimensions sum to compact H^1_deRham")
print("PASS: trivial character descends to a rational quotient differential")
print("PASS: elliptic Gauss--Manin matrices are flat with discriminant poles")
print("PASS: genus-three plane slice has a dense 6x6 Gauss--Manin matrix")
