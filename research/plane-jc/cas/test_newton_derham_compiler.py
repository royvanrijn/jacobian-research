#!/usr/bin/env python3
"""Exact regressions for the Newton-chain/de Rham compiler IR."""

from dataclasses import replace
from functools import reduce
from operator import mul

import sympy as sp

from newton_derham_compiler import (
    certify_local_system_reuse,
    compile_weighted_wronskian,
    frontier_75_125_record,
    normalized_72_108_block,
    repeated_tail_96_144_record,
)


block = normalized_72_108_block()
compiled = compile_weighted_wronskian(block)
certificate = compiled.certificate

assert compiled.genus == 3
assert compiled.character == 1
assert compiled.affine_dimension == 7
assert compiled.compact_dimension == 6
assert len(certificate.solved_coefficients) == 11
assert len(certificate.compatibility) == 6
assert certificate.residual_degrees == tuple(range(13, 19))
assert certificate.tail_infinity_orders == (5, 4, 3, 2, 1, 0)
assert certificate.tail_is_second_kind
assert certificate.identity_holds
assert certificate.tail_basis_certified
assert all(
    sp.expand(coordinate + value / 2) == 0
    for coordinate, value in zip(
        certificate.de_rham_coordinates, certificate.compatibility
    )
)

# The low-coefficient matrix is triangular with diagonal 2*k-3.  This is the
# exact independence certificate for the six tail classes: an exact linear
# combination would have a primitive of degree <=12, whose thirteen low
# coefficients force that primitive to vanish.
expected_determinant = reduce(mul, (2 * k - 3 for k in range(13)), 1)
assert certificate.low_operator_determinant == expected_determinant

# Verify the representative identity independently of the compiler flag.
t = block.t
D = certificate.primitive
residual = sp.expand(
    block.covering_exponent * block.A * sp.diff(D, t)
    - block.primitive_weight * sp.diff(block.A, t) * D
    - block.R
)
assert sp.expand(residual - certificate.residual) == 0

# The old block has an executable local-system and supported-problem
# fingerprint.  The proposed repeated-tail row honestly remains front-end
# incomplete and records the lower-side source reconciliation.
assert compiled.local_system_fingerprint[:5] == (2, 8, 1, 2, 6)

# A fingerprint is only a screening invariant.  Reuse requires an exact
# cyclic-curve identity after an explicit base/coordinate change.  The
# forcing term may change without changing the Gauss--Manin local system.
different_forcing = replace(block, R=t**3)
reuse = certify_local_system_reuse(block, different_forcing)
assert reuse.curve_identity == 0
assert reuse.covering_exponent == 2
assert reuse.character == 1

# Exercise a nontrivial coordinate rescaling: s=2*t and y_candidate=3*y.
s = sp.symbols("s")
scaled_curve = replace(
    block,
    t=s,
    A=9 * block.A.subs(t, s / 2),
)
scaled_reuse = certify_local_system_reuse(
    block,
    scaled_curve,
    t_scale=2,
    y_scale=3,
)
assert scaled_reuse.curve_identity == 0

try:
    certify_local_system_reuse(block, replace(block, A=block.A + 1))
except ValueError as error:
    assert "does not identify the curve families" in str(error)
else:
    raise AssertionError("a fingerprint-compatible but distinct curve must not certify")

try:
    certify_local_system_reuse(
        block,
        block,
        base_substitution=((sp.Symbol("a2"), t),),
    )
except ValueError as error:
    assert "may not depend on the fiber coordinate" in str(error)
else:
    raise AssertionError("a fiber-dependent substitution is not a base map")

repeated = repeated_tail_96_144_record()
assert not repeated.frontend_complete
assert len(repeated.missing_frontend_data) == 1
assert "impossible last lower corner (8,4)" in repeated.source_reconciliation
assert "q1=d0=4" in repeated.source_reconciliation
assert "remaining triple-root partition" in repeated.source_reconciliation
assert "open counts 1,6,3,0" in repeated.source_reconciliation
assert "source-excluded" in repeated.missing_frontend_data[0]
frontier = frontier_75_125_record()
assert not frontier.frontend_complete
assert frontier.multiplicities == (3, 5)
assert len(frontier.missing_frontend_data) == 3
assert "support control" in frontier.missing_frontend_data[0]
assert "exhaustive gamma branches" in frontier.missing_frontend_data[1]
assert "bracket layer 4" in frontier.missing_frontend_data[2]
assert "j=0 degree-(50,75)" in frontier.source_reconciliation

try:
    compile_weighted_wronskian(block.__class__(
        chain=repeated,
        t=block.t,
        A=block.A,
        R=block.R,
        covering_exponent=block.covering_exponent,
        primitive_weight=block.primitive_weight,
        primitive_exponents=block.primitive_exponents,
        full_primitive_bounds=block.full_primitive_bounds,
    ))
except ValueError as error:
    assert "lacks front-end data" in str(error)
else:
    raise AssertionError("an incomplete chain must not compile")

print("PASS: (72,108) support equations are certified de Rham tail coordinates")
print("PASS: the six tail classes have an exact triangular independence certificate")
print("PASS: local-system reuse requires an exact curve-family identity")
print("PASS: compiler IR refuses to invent missing frontier Laurent bands")
