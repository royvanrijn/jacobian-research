#!/usr/bin/env sage -python
"""Map the exact degree-10 residual origin to the H92 hyperplane section.

For the 3-neighbor ternary cubic, the marked degree-29 divisor has a unique
degree-10 residual point R.  Taking the tangent at R gives 2R+Q in the plane
hyperplane class, so Q is the class H-3R on the Jacobian.  This script maps Q
to the final H92 model without using the H21/P1 q=6 construction.
"""

from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[2]
RECOVERY = ROOT / "elkies-k3/scripts/recover_h92_p2_coordinates.sage"
globals()["__file__"] = str(RECOVERY)
exec(compile(RECOVERY.read_text(), str(RECOVERY), "exec"))

residual = json.loads((
    ROOT / "artifacts/generated-results/elkies-k3-h92-p2-residual-crt.json"
).read_text())
if not residual["complete"]:
    raise RuntimeError("the degree-10 residual point is not fully CRT-lifted")


def polynomial(coefficients, variable):
    return sum(Kw(QQ(value)) * variable**index for index, value in enumerate(coefficients))


target_base = Kw.gen() / (2 * Kw(s))
residual_z = polynomial(residual["residual_z"]["numerator"], target_base) / polynomial(residual["residual_z"]["denominator"], target_base)
residual_x = polynomial(residual["residual_x"]["numerator"], target_base) / polynomial(residual["residual_x"]["denominator"], target_base)

affine = PolynomialRing(Kw, names=("z", "X"))
z, X = affine.gens()
cubic = affine(sum(
    Kw(coefficient) * Kw.gen()**powers[2] * z**powers[0] * X**powers[1]
    for powers, coefficient in three_neighbor_cubic.dict().items()
))
assert cubic(residual_z, residual_x) == 0

partial_z = cubic.derivative(z)(residual_z, residual_x)
partial_x = cubic.derivative(X)(residual_z, residual_x)
assert partial_x
tangent_x = residual_x - partial_z / partial_x * (z - residual_z)
intersection = affine(cubic.subs({X: tangent_x}))
univariate = PolynomialRing(Kw, "z")(intersection)
quotient, remainder = univariate.quo_rem((univariate.parent().gen() - residual_z)**2)
assert not remainder and quotient.degree() == 1
q_z = -quotient[0] / quotient[1]
q_x = tangent_x.subs({z: q_z})
assert cubic(q_z, q_x) == 0

projective = PolynomialRing(Kw, names=("z", "X", "Z"))
zp, Xp, Zp = projective.gens()
form = projective(sum(
    coefficient * zp**powers[0] * Xp**powers[1] * Zp**(3 - powers[0] - powers[1])
    for powers, coefficient in cubic.dict().items()
))
cubic_to_jacobian = EllipticCurve_from_cubic(
    form, (residual_z, residual_x, Kw(1)), morphism=True
)

_, h92_formulas = anchor.parse_h92(ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt")
A1, A, B1, B, B2 = tuple(QQ(value(r, s)) for value in h92_formulas)
h92 = EllipticCurve([
    0, 0, 0,
    A1 * target_base**3 + A * target_base**4,
    B1 * target_base**5 + B * target_base**6 + B2 * target_base**7,
])
isomorphism = cubic_to_jacobian.codomain().isomorphism_to(h92)
candidate = isomorphism(cubic_to_jacobian((q_z, q_x, Kw(1))))
assert candidate

candidate_x, candidate_y = candidate.xy()
assert candidate_y**2 == candidate_x**3 + (
    A1 * target_base**3 + A * target_base**4
) * candidate_x + (
    B1 * target_base**5 + B * target_base**6 + B2 * target_base**7
)
print(
    "H92P2MAP|status=PASS|"
    f"x_degrees={candidate_x.numerator().degree()},{candidate_x.denominator().degree()}|"
    f"y_degrees={candidate_y.numerator().degree()},{candidate_y.denominator().degree()}",
    flush=True,
)
print("H92P2MAP|x=" + str(candidate_x), flush=True)
print("H92P2MAP|y=" + str(candidate_y), flush=True)
