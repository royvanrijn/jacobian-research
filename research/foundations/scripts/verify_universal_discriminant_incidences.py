#!/usr/bin/env python3
"""Exact audit of the universal saturated discriminant incidences."""

import sys
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import (  # noqa: E402
    incidence_elimination_generators,
    universal_discriminant_incidences,
    universal_primitive,
)


W = sp.symbols("W")
expected_contacts = {
    "ordinary_bitangent": (2, 2),
    "higher_cusp": (4,),
    "cusp_branch": (3, 2),
    "tritangent": (2, 2, 2),
}

# The admissible coefficient chart and its three weighted open conditions are
# explicit in every degree.  Each incidence carries one exact Rabinowitsch
# equation representing saturation by its regularity and distinctness factor.
for degree in range(3, 11):
    model = universal_primitive(degree, W, prefix=f"h{degree}_")
    H = model.polynomial
    assert len(model.coefficients) == degree - 2
    assert H.subs(W, 0) == 0
    assert sp.diff(H, W).subs(W, 0) == 0
    assert H.subs(W, 1) == 0
    assert sp.Poly(model.open_factor, *model.coefficients).is_zero is False

    incidences = universal_discriminant_incidences(model, prefix=f"i{degree}")
    assert set(incidences) == set(expected_contacts)
    for name, incidence in incidences.items():
        assert incidence.contact_partition == expected_contacts[name]
        assert incidence.parameters == model.coefficients
        assert incidence.saturated_equations[-1] == (
            1 - incidence.gate * incidence.saturation_factor
        )

    # The gate turns every forbidden diagonal into the unit equation.
    bitangent = incidences["ordinary_bitangent"]
    r, u = bitangent.marked_points
    assert bitangent.saturated_equations[-1].subs(r, u) == 1
    tritangent = incidences["tritangent"]
    r, u, v = tritangent.marked_points
    assert tritangent.saturated_equations[-1].subs(r, u) == 1
    assert tritangent.saturated_equations[-1].subs(r, v) == 1
    assert tritangent.saturated_equations[-1].subs(u, v) == 1

# Each universal incidence vanishes identically on its contact-factor
# parameterization, while its saturation factor remains generically nonzero.
minimal_models = {
    "ordinary_bitangent": universal_primitive(4, W, "b_", admissible=False),
    "higher_cusp": universal_primitive(4, W, "c_", admissible=False),
    "cusp_branch": universal_primitive(5, W, "cb_", admissible=False),
    "tritangent": universal_primitive(6, W, "tt_", admissible=False),
}
for name, model in minimal_models.items():
    incidence = universal_discriminant_incidences(model, prefix=f"contact_{name}")[name]
    points = incidence.marked_points
    contact_polynomial = sp.prod(
        (W - point) ** multiplicity
        for point, multiplicity in zip(points, incidence.contact_partition)
    )
    coefficient_substitution = {
        coefficient: sp.Poly(contact_polynomial, W).coeff_monomial(W**power)
        for power, coefficient in zip(range(2, model.degree + 1), model.coefficients)
    }
    assert all(
        sp.factor(equation.subs(coefficient_substitution)) == 0
        for equation in incidence.equations
    )
    assert sp.factor(incidence.saturation_factor.subs(coefficient_substitution)) != 0

# A real elimination check: the saturated degree-four higher-cusp incidence
# recovers exactly the discriminant hypersurface of H'' in coefficient space.
quartic = universal_primitive(4, W, "q4_")
higher_cusp = universal_discriminant_incidences(quartic, "q4_inc")["higher_cusp"]
eliminated = incidence_elimination_generators(higher_cusp)
assert len(eliminated) == 1
cusp_discriminant = sp.factor(sp.discriminant(sp.diff(quartic.polynomial, W, 2), W))
ratio = sp.cancel(eliminated[0] / cusp_discriminant)
assert not ratio.free_symbols and ratio != 0

print("PASS: universal admissible coefficient charts and open factors")
print("PASS: Rabinowitsch gates saturate every diagonal and cusp factor")
print("PASS: contact factorizations satisfy all four universal incidences")
print("PASS: elimination recovers the quartic higher-cusp coefficient locus")
