#!/usr/bin/env python3
"""Exact regressions for the Newton-chain/de Rham compiler IR."""

from functools import reduce
from operator import mul

import sympy as sp

from newton_derham_compiler import (
    compile_weighted_wronskian,
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
repeated = repeated_tail_96_144_record()
assert not repeated.frontend_complete
assert len(repeated.missing_frontend_data) == 5
assert "impossible last lower corner (8,4)" in repeated.source_reconciliation

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
print("PASS: compiler IR refuses to invent missing (96,144) Laurent bands")
