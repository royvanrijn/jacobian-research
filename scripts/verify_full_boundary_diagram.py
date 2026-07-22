#!/usr/bin/env python3
"""Exact regression for the full cancellation boundary-prime diagram."""

from __future__ import annotations

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.boundary import cancellation_full_boundary_diagram  # noqa: E402


w, v = sp.symbols("w v")

for m in range(1, 6):
    for r in range(1, 6):
        diagram = cancellation_full_boundary_diagram(m, r)
        K = diagram.branch_polynomial
        L = diagram.critical_coefficient

        # Reconstruct both defining integrals independently of the library.
        assert sp.factor(
            sp.diff(w ** (r + 1) * K, w)
            - w**r * (1 - w) ** (m * r)
        ) == 0
        p = lambda value: value * (1 - value) ** m
        independent_L = sp.cancel(
            sp.integrate((p(w) - p(v)) ** r, (v, 0, w)) / w ** (r + 1)
        )
        assert sp.factor(L - independent_L) == 0
        assert sp.degree(K, w) == m * r
        assert sp.degree(L, w) == m * r
        assert sp.Poly(K, w).coeff_monomial(1) == sp.Rational(1, r + 1)
        assert sp.Poly(K, w).coeff_monomial(w) == -sp.Rational(m * r, r + 2)
        assert sp.Poly(L, w).coeff_monomial(1) == sp.Rational(1, r + 1)
        assert sp.Poly(L, w).coeff_monomial(w) == -sp.Rational(
            m * r * (r + 3), (r + 1) * (r + 2)
        )

        # K is a fractional-linear transform of the cancellation parameter
        # polynomial M.  Hence every proved irreducibility range for M also
        # proves the nonzero contact resultant without computing it.
        q = sp.Symbol("q")
        n = m * r
        M = sum(
            (-1) ** j * sp.binomial(n + r + 1, j) * q ** (n - j)
            for j in range(n + 1)
        )
        transformed_K = sp.factor(
            sp.cancel((1 - q) ** n * K.subs(w, -q / (1 - q)))
        )
        assert sp.rem(sp.Poly(M, q), sp.Poly(transformed_K, q)) == 0
        assert diagram.branch_discriminant != 0
        assert diagram.contact_resultant != 0
        assert diagram.contact_certified

        boundary_count = m * r - 1
        assert len(diagram.boundary_primes) == 1 + boundary_count
        assert len(diagram.cluster_primes) == m * r
        assert diagram.critical_cluster_total_length == m * m * r
        assert diagram.critical_boundary_total_length == m * boundary_count

        critical_edges = [
            edge for edge in diagram.intersections if edge.left == "E_Delta"
        ]
        cluster_edges = [
            edge for edge in diagram.intersections if edge.left != "E_Delta"
        ]
        assert len(critical_edges) == boundary_count
        assert all(edge.contact_length == m for edge in critical_edges)
        assert len(cluster_edges) == boundary_count * (boundary_count - 1) // 2
        assert all(edge.contact_length == 1 for edge in cluster_edges)

# Concrete values independently reproduced by the normalization experiments.
assert cancellation_full_boundary_diagram(2, 1).critical_cluster_total_length == 4
assert cancellation_full_boundary_diagram(1, 2).critical_cluster_total_length == 2
assert cancellation_full_boundary_diagram(3, 1).critical_cluster_total_length == 9
assert cancellation_full_boundary_diagram(2, 2).critical_cluster_total_length == 8

print("PASS full cancellation boundary diagram for 1 <= m,r <= 5")
print("PASS every critical/K-branch contact has length m")
print("PASS every distinct K-branch intersection is the reduced central stratum")
