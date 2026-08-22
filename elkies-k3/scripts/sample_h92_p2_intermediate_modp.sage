#!/usr/bin/env sage -python
"""Produce one fast modular record for the marked intermediate H92 section.

The canonical transport construction is in
``recover_h92_p2_coordinates.sage``.  This worker loads only its exact
two-neighbor prefix, evaluates the marked degree-three divisor on enough split
finite fibers, and interpolates the intermediate Jacobian coordinates.  It
deliberately avoids the source file's 2,000-fiber diagnostic sweep.
"""

from pathlib import Path
import json
import os


ROOT = Path(__file__).resolve().parents[2]
RECOVERY = ROOT / "elkies-k3/scripts/recover_h92_p2_coordinates.sage"
transport_source = RECOVERY.read_text().split("axis_candidates =", 1)[0]
globals()["__file__"] = str(RECOVERY)
exec(compile(transport_source, str(RECOVERY), "exec"))

prime = int(os.environ["H92P2_PRIME"])
field = GF(prime)
needed = int(os.environ.get("H92P2_NEEDED_SAMPLES", "70"))
scan_limit = int(os.environ.get("H92P2_SCAN_LIMIT", "1200"))
output = Path(os.environ["H92P2_MODULAR_OUTPUT"])


def reduce_function(value, base_value):
    numerator = value.numerator()
    denominator = value.denominator()
    top = sum(field(coefficient) * base_value**index
              for index, coefficient in enumerate(numerator.list()))
    bottom = sum(field(coefficient) * base_value**index
                 for index, coefficient in enumerate(denominator.list()))
    if not bottom:
        raise ZeroDivisionError
    return top / bottom


def specialize_univariate(polynomial, base_value):
    ring = PolynomialRing(field, "z")
    variable = ring.gen()
    return ring(sum(
        reduce_function(coefficient, base_value) * variable**index
        for index, coefficient in enumerate(polynomial.list())
    ))


def marked_divisor_class(base_value):
    f0 = specialize_univariate(quartic, base_value)
    g0 = specialize_univariate(source_divisor, base_value)
    roots = g0.roots(field)
    if len(roots) != 3 or any(multiplicity != 1 for _, multiplicity in roots):
        raise ValueError("non-split marked degree-three fiber")
    numerator = specialize_univariate(source_v_numerator, base_value)
    denominator = specialize_univariate(source_v_denominator, base_value)
    h0 = (numerator * denominator.inverse_mod(g0)).mod(g0)
    complement = (f0 - h0**2) // g0
    if complement.degree() != 1:
        raise ValueError("no rational complement")
    qx = -complement[0] / complement[1]
    qy = h0(qx)
    z0 = f0.parent().gen()
    shifted = f0(z0 + qx)
    qa, qb, qc, qd = shifted[4], shifted[3], shifted[2], shifted[1]
    long_curve = EllipticCurve([
        qd / qy, qc - qd**2 / (4 * qy**2), 2 * qy * qb,
        -4 * qy**2 * qa,
        -(qc - qd**2 / (4 * qy**2)) * 4 * qy**2 * qa,
    ])
    answer = long_curve(0)
    for root, _ in roots:
        delta = root - qx
        y_value = h0(root)
        x_value = (2 * qy * (y_value + qy) + qd * delta) / delta**2
        y_image = (
            4 * qy**2 * (y_value + qy)
            + 2 * qy * (qd * delta + qc * delta**2)
            - qd**2 * delta**2 / (2 * qy)
        ) / delta**3
        answer += long_curve(x_value, y_image)
    jacobian = EllipticCurve([
        0, 0, 0,
        reduce_function(weierstrass_2neighbor[1], base_value),
        reduce_function(weierstrass_2neighbor[0], base_value),
    ])
    return long_curve.isomorphism_to(jacobian)(answer)


values = []
for integer in range(1, scan_limit + 1):
    try:
        point = marked_divisor_class(field(integer))
    except (ArithmeticError, ValueError, ZeroDivisionError):
        continue
    if point:
        values.append((field(integer), point.xy()[0], point.xy()[1]))
    if len(values) >= needed:
        break
if len(values) < needed:
    raise RuntimeError(f"only {len(values)} split marked fibers found")


def interpolate(index, numerator_degree, denominator_degree):
    rows = [
        [-entry[0]**power for power in range(numerator_degree + 1)]
        + [entry[index] * entry[0]**power
           for power in range(denominator_degree + 1)]
        for entry in values
    ]
    kernel = matrix(field, rows).right_kernel()
    if kernel.dimension() != 1:
        raise RuntimeError(f"unexpected interpolation kernel {kernel.dimension()}")
    relation = kernel.basis()[0]
    denominator = tuple(relation[numerator_degree + 1:])
    scale = denominator[-1]
    if not scale:
        raise RuntimeError("bad denominator normalization")
    return (
        tuple(int(value / scale) for value in relation[:numerator_degree + 1]),
        tuple(int(value / scale) for value in denominator),
    )


x = interpolate(1, 22, 18)
y = interpolate(2, 33, 27)
payload = {"prime": prime, "samples": len(values), "x": x, "y": y}
output.write_text(json.dumps(payload, sort_keys=True) + "\n")
print(
    f"H92P2FAST|prime={prime}|samples={len(values)}|"
    "x_degrees=22,18|y_degrees=33,27|status=PASS",
    flush=True,
)
