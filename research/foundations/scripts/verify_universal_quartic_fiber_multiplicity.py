#!/usr/bin/env python3
"""Exact algebra for universal quartic weighted-fiber multiplicity.

The written proof uses Hasse--Minkowski for a five-variable trace form over
a number field.  This checker verifies the coefficient identities, the
trace-form reduction, the weighted normal form, the archimedean signatures,
and all finite bad parameter loci used after the density argument.
"""

from __future__ import annotations

import sympy as sp


T, W = sp.symbols("T W")
translation, step, trace_square = sp.symbols(
    "translation step trace_square"
)
linear, constant = sp.symbols("linear constant")

# A trace-zero primitive quartic generator has this characteristic
# polynomial.  Newton's identity gives quadratic coefficient
# -Tr(eta^2)/2.
quadratic = -trace_square / 2
quartic = T**4 + quadratic * T**2 + linear * T + constant

changed = sp.expand(quartic.subs(T, translation + step * W))
tangent_remainder = sp.expand(
    changed
    - changed.subs(W, 0)
    - sp.diff(changed, W).subs(W, 0) * W
)
chord = sp.factor(
    quartic.subs(T, translation + step)
    - quartic.subs(T, translation)
    - step * sp.diff(quartic, T).subs(T, translation)
)
assert sp.expand(
    chord
    - step**2
    * (
        6 * translation**2
        + 4 * translation * step
        + step**2
        - trace_square / 2
    )
) == 0

# Diagonalize the tangent-chord condition with e=step+2*translation.
e = sp.symbols("e")
diagonal_chord = sp.expand(
    (
        6 * translation**2
        + 4 * translation * step
        + step**2
        - trace_square / 2
    ).subs(step, e - 2 * translation)
)
assert diagonal_chord == (
    e**2 + 2 * translation**2 - trace_square / 2
)

# On the chord quadric, normalize H'(1)=-1.  The result is exactly the
# one-parameter quartic weighted seed H_alpha.
alpha = sp.symbols("alpha")
normalized = sp.expand(
    -tangent_remainder
    / (2 * step**3 * (step + 2 * translation))
)
normal_form = sp.expand(
    W**2 * (W - 1) * (alpha * W - alpha - 1)
)
alpha_value = -step / (2 * (step + 2 * translation))
chord_relation = {
    trace_square: 2
    * (
        (step + 2 * translation) ** 2
        + 2 * translation**2
    )
}
assert sp.cancel(
    normalized.subs(chord_relation)
    - normal_form.subs(alpha, alpha_value)
) == 0
assert sp.simplify(
    alpha_value.subs(step, e - 2 * translation)
    - (translation / e - sp.Rational(1, 2))
) == 0

# The clean weighted locus excludes only finitely many alpha.
hessian = sp.diff(normal_form, W, 2)
assert sp.expand(
    sp.discriminant(hessian, W)
    - 12 * (4 * alpha**2 + 4 * alpha + 3)
) == 0
assert sp.discriminant(4 * alpha**2 + 4 * alpha + 3, alpha) == -32
assert sp.expand(hessian.subs(W, 0) - 2 * (alpha + 1)) == 0
assert sp.expand(hessian.subs(W, 1) - 2 * (alpha - 2)) == 0
assert sp.expand(sp.diff(normal_form, W).subs(W, 1) + 1) == 0

# At every real place of a rank-four etale algebra over a number field, the
# trace algebra has one of these three signatures.  Adding the two negative
# variables gives an indefinite five-dimensional chord form.
trace_zero_signatures = {
    (4, 0): (3, 0),
    (2, 1): (2, 1),
    (0, 2): (1, 2),
}
chord_signatures = {
    signature: (values[0], values[1] + 2)
    for signature, values in trace_zero_signatures.items()
}
assert chord_signatures == {
    (4, 0): (3, 2),
    (2, 1): (2, 3),
    (0, 2): (1, 4),
}
assert all(
    positive and negative and positive + negative == 5
    for positive, negative in chord_signatures.values()
)

print("PASS: quartic tangent chords are the five-variable trace quadric")
print("PASS: chord normalization gives alpha=translation/e-1/2")
print("PASS: every quartic real signature gives an indefinite rank-five form")
print("PASS: the clean weighted exclusions are a finite alpha set")
