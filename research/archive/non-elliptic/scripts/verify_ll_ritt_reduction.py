#!/usr/bin/env python3
"""Verify the LL profile and degree-12/30 Ritt coverage reduction."""

from __future__ import annotations

import pathlib
import sys

import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.ll_ritt import (
    ll_signature,
    prime_words,
    proper_outer_cuts,
    selected_relation_graph_types,
)
from jcsearch.ritt_complex import dickson


def labels(edge: frozenset[object]) -> frozenset[str]:
    return frozenset(vertex.label for vertex in edge)


# A two-factor (a,b) composition has the generic profile
# b^(a-1) 1^(b-1), and the count recovers both Riemann--Hurwitz and the
# known N-a-b marked codimension in the N-3 dimensional seed space.
for degree in (12, 30):
    for outer in proper_outer_cuts(degree):
        inner = degree // outer
        signature = ll_signature((outer, inner))
        assert signature.expanded_transposition_multiplicities == (
            (inner,) * (outer - 1) + (1,) * (inner - 1)
        )
        assert signature.ramification_point_count == degree - 1
        assert (degree - 3) - signature.marked_dimension == degree - outer - inner


# Every complete factor word has sum(d_i-1) generic branch values, hence the
# same marked dimension as its normalized factor chart.
for degree in (12, 30):
    for word in prime_words(degree):
        degree_word = tuple(occurrence.degree for occurrence in word)
        signature = ll_signature(degree_word)
        assert signature.branch_value_count == sum(
            factor_degree - 1 for factor_degree in degree_word
        )
        assert signature.marked_dimension == sum(
            factor_degree - 1 for factor_degree in degree_word
        ) - 1


# In degree 12, every choice of prime refinements carrying all four cuts
# forces the same relation graph: two ordered quadratic occurrences joined
# bidirectionally through the cubic.  The entire graph is one SCC.
cuts12 = proper_outer_cuts(12)
assert cuts12 == (2, 3, 4, 6)
types12 = selected_relation_graph_types(12, cuts12)
assert sum(types12.values()) == 4
assert len(types12) == 1
graph12 = next(iter(types12))
assert len(graph12.strongly_connected_components()) == 1
assert {labels(edge) for edge in graph12.two_way_edges} == {
    frozenset(("2", "3")),
    frozenset(("2_2", "3")),
}
assert {
    (left.label, right.label) for left, right in graph12.one_way_edges
} == {("2", "2_2")}


# In degree 30 there are 2^6 choices of prime refinements, but all 64 force
# the same fully bidirected triangle on degrees 2,3,5.  Hence there is one
# strongly connected collision block before any coefficient calculation.
cuts30 = proper_outer_cuts(30)
assert cuts30 == (2, 3, 5, 6, 10, 15)
types30 = selected_relation_graph_types(30, cuts30)
assert sum(types30.values()) == 64
assert len(types30) == 1
graph30 = next(iter(types30))
assert len(graph30.strongly_connected_components()) == 1
assert {labels(edge) for edge in graph30.two_way_edges} == {
    frozenset(("2", "3")),
    frozenset(("2", "5")),
    frozenset(("3", "5")),
}
assert not graph30.one_way_edges


# The tame multi-collision normal form for these forced SCCs is the Dickson
# family.  Its generic finite branch locus has two values.  The power
# degeneration is z=0 and is therefore contained in the Dickson closure.
w, target = sp.symbols("w target")
for degree in (12, 30):
    polynomial = dickson(degree, w, sp.Integer(1))
    critical_value_resultant = sp.factor(
        sp.resultant(sp.diff(polynomial, w), polynomial - target, w)
    )
    assert sp.degree(critical_value_resultant, target) == degree - 1
    assert sp.degree(
        sp.Poly(critical_value_resultant, target).sqf_part(),
        target,
    ) == 2
    translation, parameter = sp.symbols(
        f"translation_{degree} parameter_{degree}"
    )
    family = sp.expand(
        dickson(degree, w + translation, parameter)
        - dickson(degree, translation, parameter)
    )
    assert sp.expand(
        family.subs(parameter, 0)
        - ((w + translation) ** degree - translation**degree)
    ) == 0


print("PASS: generic (a,b) composition has LL profile b^(a-1) 1^(b-1)")
print("PASS: LL branch counts recover the exact codimension N-a-b")
print("PASS: all degree-12 prime-refinement choices force one relation graph")
print("PASS: all 64 degree-30 choices force the fully bidirected 2,3,5 block")
print("PASS: the forced degree-12/30 Dickson families have two branch values")
print("PASS: the power family is the z=0 boundary of the Dickson family")
