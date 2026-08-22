#!/usr/bin/env sage
"""Regression tests for elliptic_neighbor_compiler_field_generic.sage."""

from pathlib import Path

from sage.all import PolynomialRing, QQ, QuadraticField, matrix, vector

HERE = Path(__file__).resolve().parent
load(str(HERE / "elliptic_neighbor_compiler_field_generic.sage"))

K = QuadraticField(-3, "j")
j = K.gen()
ambient = ("a", "b", "c", "d")

images = {
    "a": (K(1), K(0)),
    "b": (K(0), K(1)),
    "c": (K(0), K(0)),
    "d": (K(0), K(0)),
}
block = quotient_condition(
    "K-block",
    ambient,
    lambda name: images[name],
    ("r0", "r1"),
    "number-field regression",
    coefficient_field=K,
)
assert block["matrix"].base_ring() is K
compiled = compile_resolved_conditions(ambient, (block,), complete=True)
assert compiled["condition_matrix"].base_ring() is K
assert compiled["rank"] == 2
assert compiled["kernel_dimension"] == 2
assert compiled["h0_certified"]

# Mixed QQ/K blocks must promote to K without changing rank.
qq_block = {
    "name": "QQ-zero",
    "matrix": matrix(QQ, 1, 4, [0, 0, 0, 0]),
    "quotient_basis": ("zero",),
    "provenance": "QQ compatibility regression",
}
mixed = compile_resolved_conditions(ambient, (qq_block, block), complete=True)
assert mixed["condition_matrix"].base_ring() is K
assert mixed["rank"] == 2
assert mixed["kernel_dimension"] == 2

# Resolved-chart quotient must preserve the local coefficient field.
R = PolynomialRing(K, "t")
t = R.gen()
resolved = resolved_chart_quotient_condition(
    "resolved-K",
    (R(1), t, t**2),
    R,
    lambda value: value,
    R.ideal(t**2),
    (R(1), t),
    "univariate K[t]/(t^2) regression",
)
assert resolved["matrix"].base_ring() is K
assert resolved["matrix"].rank() == 2
resolved_compiled = compile_resolved_conditions(
    (R(1), t, t**2), (resolved,), complete=False
)
assert resolved_compiled["condition_matrix"].base_ring() is K
assert resolved_compiled["rank"] == 2
assert resolved_compiled["kernel_dimension"] == 1

# Historical default remains QQ.
qq_default = quotient_condition(
    "QQ-default",
    ("x", "y"),
    lambda name: (1,) if name == "x" else (0,),
    ("r",),
    "default-field regression",
)
assert qq_default["matrix"].base_ring() is QQ
qq_compiled = compile_resolved_conditions(("x", "y"), (qq_default,))
assert qq_compiled["condition_matrix"].base_ring() is QQ

print(
    "ELLIPTICNEIGHBORFIELD|K_rank=2|K_nullity=2|mixed_field=K|"
    "resolved_rank=2|QQ_default=1|status=PASS_FIELD_GENERIC_COMPILER",
    flush=True,
)
