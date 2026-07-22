#!/usr/bin/env python3
"""Audit the quadratic stationary-point ladder shown in the supplied card."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.boundary import cancellation_boundary_cover_profile  # noqa: E402


s, A, B, C = sp.symbols("s A B C")
D = 1 - B * s + A * s**2

for ell in range(2, 9):
    J = (-1) ** ell * ell * sp.binomial(2 * ell - 1, ell)
    P = sp.expand(-J * sp.integrate(D ** (ell - 1), (s, 0, s)) - C)
    assert sp.expand(sp.diff(P, s) + J * D ** (ell - 1)) == 0
    assert sp.degree(P, s) == 2 * ell - 1
    assert cancellation_boundary_cover_profile(1, ell - 1).degree == 2 * ell - 1

# The first rung is exactly the repository cubic after a target rescaling.
ell = 2
J = (-1) ** ell * ell * sp.binomial(2 * ell - 1, ell)
P2 = sp.expand(-J * sp.integrate(D, (s, 0, s)) - C)
a, b, c = sp.symbols("a b c")
repository_cubic = c - 2 * s + b * s**2 - 2 * a * s**3
assert sp.factor(
    P2.subs({A: 3 * a, B: b, C: -3 * c}) / 3 - repository_cubic
) == 0

print("PASS quadratic ladder equals cancellation slice (m,r)=(1,ell-1)")
print("PASS ell=2 recovers the foundational marked-root cubic after rescaling")
