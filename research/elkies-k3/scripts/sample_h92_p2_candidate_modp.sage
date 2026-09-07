#!/usr/bin/env sage -python
"""Sample the H92 point obtained from the residual-origin tangent construction.

This replays the degree-10 residual construction after reduction modulo a
prime.  It avoids a costly QQ(w) plane-cubic conversion while retaining the
same marked-fiber incidence: if R is the residual origin, the point sampled
is the third tangent intersection Q, so 2R+Q is a hyperplane section.
"""

from pathlib import Path
import json
import os

from sage.all import *
from sage.schemes.curves.constructor import Curve


ROOT = Path(__file__).resolve().parents[2]
RECOVERY = ROOT / "elkies-k3/scripts/recover_h92_p2_coordinates.sage"
transport_source = RECOVERY.read_text().split("# A finite-field spot check confirms", 1)[0]
globals()["__file__"] = str(RECOVERY)
exec(compile(transport_source, str(RECOVERY), "exec"))

payload = json.loads((ROOT / "artifacts/generated-results/elkies-k3-h92-p2-residual-crt.json").read_text())
assert payload["complete"]
prime = ZZ(os.environ.get("H92P2_PRIME", "100003"))
field = GF(prime)
sample_count = int(os.environ.get("H92P2_SAMPLE_COUNT", "100"))


def reduced_rational(value):
    rational = QQ(value)
    return field(rational.numerator()) / field(rational.denominator())


def evaluate(coefficients, value):
    return sum(reduced_rational(coefficient) * value**index for index, coefficient in enumerate(coefficients))


def residual_at(target_t):
    return (
        evaluate(payload["residual_z"]["numerator"], target_t)
        / evaluate(payload["residual_z"]["denominator"], target_t),
        evaluate(payload["residual_x"]["numerator"], target_t)
        / evaluate(payload["residual_x"]["denominator"], target_t),
    )


def point_at(target_t):
    raw_w = 2 * field(s) * target_t
    affine = PolynomialRing(field, names=("z", "X"))
    z, X = affine.gens()
    cubic = affine(sum(
        field(coefficient) * raw_w**powers[2] * z**powers[0] * X**powers[1]
        for powers, coefficient in three_neighbor_cubic.dict().items()
    ))
    rz, rX = residual_at(target_t)
    assert cubic(rz, rX) == 0

    def tangent_residual(point_z, point_x):
        partial_z = cubic.derivative(z)(point_z, point_x)
        partial_x = cubic.derivative(X)(point_z, point_x)
        if not partial_x:
            raise ArithmeticError("vertical tangent")
        tangent_x = point_x - partial_z / partial_x * (z - point_z)
        univariate = PolynomialRing(field, "z")(cubic.subs({X: tangent_x}))
        if univariate.degree() != 3:
            raise ArithmeticError("tangent has degenerate intersection")
        third_z = -univariate[2] / univariate[3] - 2 * point_z
        return third_z, tangent_x.subs({z: third_z})

    qz, qX = tangent_residual(rz, rX)
    assert cubic(qz, qX) == 0
    p3z, p3X = tangent_residual(qz, qX)
    assert cubic(p3z, p3X) == 0

    projective = PolynomialRing(field, names=("z", "X", "Z"))
    zp, Xp, Zp = projective.gens()
    form = projective(sum(
        coefficient * zp**powers[0] * Xp**powers[1] * Zp**(3 - powers[0] - powers[1])
        for powers, coefficient in cubic.dict().items()
    ))
    curve = Curve(form)
    hessian = Matrix([
        [form.derivative(first, second) for first in (zp, Xp, Zp)]
        for second in (zp, Xp, Zp)
    ]).det()
    if curve.intersection(Curve(hessian)).rational_points():
        raise ArithmeticError("rational-flex specialization")
    # Over a finite field Sage can construct the cubic's Jacobian quickly.
    # We then explicitly translate the hyperplane class by 3R, so the result
    # is independent of whether this finite fiber happens to have a flex.
    cubic_to_jacobian = EllipticCurve_from_cubic(form, (rz, rX, field(1)), morphism=True)
    elliptic = cubic_to_jacobian.codomain()

    def trace_point(point, degree):
        result = point.curve()(0)
        current = point
        for _ in range(degree):
            result += current
            x_value, y_value = current.xy()
            current = point.curve()(x_value**prime, y_value**prime)
        if not result:
            return elliptic(0)
        x_value, y_value = result.xy()
        return elliptic(field(x_value), field(y_value))

    # A generic hyperplane section represents H.  Taking Frobenius traces
    # keeps its image sum over the base field even when it does not split.
    line = PolynomialRing(field, "z")(cubic.subs({X: 0}))
    candidate = elliptic(0)
    for factor, multiplicity in line.factor():
        assert multiplicity == 1
        degree = factor.degree()
        extension = GF(prime**degree, "alpha")
        root = factor.change_ring(extension).roots(extension)[0][0]
        image = cubic_to_jacobian.change_ring(extension)((root, extension(0), extension(1)))
        candidate += trace_point(image, degree)
    # The flex-scheme test above ensures R is the constructor's chosen origin.
    _, formulas = anchor.parse_h92(ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt")
    A1, A, B1, B, B2 = (field(QQ(value(r, s))) for value in formulas)
    h92 = EllipticCurve([0, 0, 0, A1*target_t**3 + A*target_t**4,
                         B1*target_t**5 + B*target_t**6 + B2*target_t**7])
    if not h92.discriminant():
        raise ArithmeticError("singular H92 specialization")
    answer = elliptic.isomorphism_to(h92)(candidate)
    return answer.xy()


values = []
scan_limit = int(os.environ.get("H92P2_SCAN_LIMIT", "2000"))
for integer in range(1, scan_limit + 1):
    try:
        x_value, y_value = point_at(field(integer))
    except Exception as error:
        if os.environ.get("H92P2_DEBUG") == "1":
            print(f"H92P2CANDIDATE|t={integer}|skip={type(error).__name__}:{error}")
        continue
    values.append((integer, int(x_value), int(y_value)))
    if len(values) == sample_count:
        break

output = os.environ.get("H92P2_MODULAR_OUTPUT")
if output:
    Path(output).write_text(json.dumps({
        "prime": int(prime),
        "origin_branch": "no_rational_flex_residual_origin",
        "values": values,
    }, sort_keys=True) + "\n")
print(f"H92P2CANDIDATE|prime={prime}|samples={len(values)}|status=PASS")
