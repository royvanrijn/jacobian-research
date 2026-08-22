#!/usr/bin/env sage -python
"""Evaluate the marked H92 degree-21 divisor as a final Jacobian point.

This is a one-fiber diagnostic for the explicit H92 P2 transport.  It is
deliberately independent of the H21/P1 q=6 construction.
"""

from pathlib import Path
import json
import os


ROOT = Path(__file__).resolve().parents[2]
RECOVERY = ROOT / "elkies-k3/scripts/recover_h92_p2_coordinates.sage"
globals()["__file__"] = str(RECOVERY)
exec(compile(RECOVERY.read_text(), str(RECOVERY), "exec"))

payload = json.loads((
    ROOT / "artifacts/generated-results/elkies-k3-h92-p2-intermediate-crt.json"
).read_text())
if any(value is None for part in payload["x"].values() for value in part):
    raise RuntimeError("the intermediate x coordinate has not fully CRT-lifted")


def lift_polynomial(coefficients):
    return sum(Ku(QQ(value)) * Ku.gen()**index for index, value in enumerate(coefficients))


intermediate_x = lift_polynomial(payload["x"]["numerator"]) / lift_polynomial(payload["x"]["denominator"])
intermediate_y = (
    intermediate_x**3
    + weierstrass_2neighbor[1] * intermediate_x
    + weierstrass_2neighbor[0]
).sqrt()

# Fix the square-root branch by the original marked divisor, not by a choice
# made during rational reconstruction.
orientation_field = spot_base
for old_base in orientation_field:
    try:
        marked = marked_divisor_class(orientation_field, old_base)
        if reduce_function(intermediate_y, orientation_field, old_base) != marked.xy()[1]:
            intermediate_y = -intermediate_y
        break
    except (ArithmeticError, ValueError, ZeroDivisionError):
        continue
else:
    raise RuntimeError("no good orientation specialization")

field = GF(int(os.environ.get("H92P2_PRIME", "100003")))

first_x = (
    4 * intermediate_x
    - shift0 - shift1 * Ku.gen() - shift2 * Ku.gen()**2 - shift3 * Ku.gen()**3
) / normalization_scale
first_y = 8 * intermediate_y / (s + r + r**2)**3
assert first_y**2 == Ku(p_first(first_x))
marked_w_function = new_base_scale * (
    first_y + Ku(transport_A) * first_x + Ku(transport_B)
) / Ku.gen()**6
plane_check = sum(
    Ku(coefficient) * Ku.gen()**powers[0] * (first_x / Ku.gen()**3)**powers[1]
    * marked_w_function**powers[2]
    for powers, coefficient in three_neighbor_cubic.dict().items()
)
assert plane_check == 0


def reduce_polynomial_qq(polynomial, field):
    ring = PolynomialRing(field, "u")
    variable = ring.gen()
    return ring(sum(field(coefficient) * variable**index for index, coefficient in enumerate(polynomial.list())))


def reduce_rational_mod(value, modulus, field):
    numerator = reduce_polynomial_qq(value.numerator(), field)
    denominator = reduce_polynomial_qq(value.denominator(), field)
    return (numerator * denominator.inverse_mod(modulus)).mod(modulus)


def image_on_extended_codomain(morphism, point, extension):
    # ``morphism`` has already been base-changed by the caller.  A second
    # base extension introduces an unrelated copy of the extension generator
    # and invalidates otherwise correct projective coordinates.
    codomain = morphism.codomain()
    domain_point = morphism.domain()(point)
    return morphism(domain_point)


def trace_point(point, prime, degree):
    curve = point.curve()
    result = curve(0)
    current = point
    for _ in range(degree):
        result += current
        x_value, y_value = current.xy()
        current = curve(x_value**prime, y_value**prime)
    return result


def image_trace_from_factor(morphism, factor, point_from_root, prime):
    degree = factor.degree()
    extension = GF(prime**degree, "alpha")
    root = factor.change_ring(extension).roots(extension)[0][0]
    point = point_from_root(extension, root)
    image = image_on_extended_codomain(morphism.change_ring(extension), point, extension)
    traced = trace_point(image, prime, degree)
    if not traced:
        return morphism.codomain()(0)
    x_value, y_value = traced.xy()
    return morphism.codomain()(morphism.codomain().base_ring()(x_value), morphism.codomain().base_ring()(y_value))


def final_point_at(field, target_t):
    """Return D-7H on the final H92 fiber at target_t over field."""
    prime = field.characteristic()
    raw_w = 2 * field(s) * target_t
    affine = PolynomialRing(field, names=("z", "X"))
    z_value, X_value = affine.gens()
    specialized = affine(sum(
        field(coefficient) * raw_w**powers[2] * z_value**powers[0] * X_value**powers[1]
        for powers, coefficient in three_neighbor_cubic.dict().items()
    ))
    projective = PolynomialRing(field, names=("z", "X", "Z"))
    z_projective, X_projective, Z_projective = projective.gens()
    form = projective(sum(
        coefficient * z_projective**powers[0] * X_projective**powers[1] * Z_projective**(3 - powers[0] - powers[1])
        for powers, coefficient in specialized.dict().items()
    ))
    u_ring = PolynomialRing(field, "u")
    u_value = u_ring.gen()
    numerator = reduce_polynomial_qq(marked_w_function.numerator(), field)
    denominator = reduce_polynomial_qq(marked_w_function.denominator(), field)
    divisor_polynomial = numerator - raw_w * denominator
    assert divisor_polynomial.degree() == 29

    # A degree-10 plane form has 30 independent restrictions to a smooth
    # plane cubic.  Its unique restriction vanishing on the rational
    # degree-29 marked divisor has one residual point R.  Since
    # 10*3 - 29 = 1, R is a rational point and trivializes the torsor.
    reduced_x = reduce_rational_mod(first_x, divisor_polynomial, field)
    reduced_x *= u_value.inverse_mod(divisor_polynomial)**3
    basis = [
        (power_z, power_x)
        for power_x in range(3)
        for power_z in range(11 - power_x)
    ]
    evaluations = [
        (u_value**power_z * reduced_x**power_x).mod(divisor_polynomial)
        for power_z, power_x in basis
    ]
    coefficient_rows = [
        [evaluation[power] for evaluation in evaluations]
        for power in range(29)
    ]
    kernel = matrix(field, coefficient_rows).right_kernel()
    if kernel.dimension() != 1:
        raise RuntimeError(f"unexpected degree-10 restriction kernel {kernel.dimension()}")
    relation = kernel.basis()[0]
    residual_form = affine(sum(
        relation[index] * z_value**power_z * X_value**power_x
        for index, (power_z, power_x) in enumerate(basis)
    ))
    residual_polynomial = specialized.resultant(residual_form, X_value)
    z_ring = PolynomialRing(field, "z")
    residual_polynomial = z_ring(residual_polynomial)
    divisor_in_z = z_ring(divisor_polynomial)
    quotient, remainder = residual_polynomial.quo_rem(divisor_in_z)
    if remainder or quotient.degree() != 1:
        raise RuntimeError("degree-10 residual did not isolate one affine point")
    residual_z = -quotient[0] / quotient[1]
    residual_x_candidates = specialized.subs({z_value: residual_z}).univariate_polynomial().roots(field)
    residual_x = next(
        candidate
        for candidate, multiplicity in residual_x_candidates
        if multiplicity == 1 and residual_form.subs({z_value: residual_z, X_value: candidate}) == 0
    )
    base_point = (residual_z, residual_x, field(1))
    if os.environ.get("H92P2_RESIDUAL_ONLY") == "1":
        return residual_z, residual_x
    print("H92P2FINAL|stage=degree10_residual", flush=True)
    cubic_to_jacobian = EllipticCurve_from_cubic(form, base_point, morphism=True)
    jacobian = cubic_to_jacobian.codomain()

    # Any non-tangent line gives the same degree-three hyperplane class.
    line_polynomial = specialized.subs({X_value: 0}).univariate_polynomial()
    print(f"H92P2FINAL|stage=hyperplane_factor|degree={line_polynomial.degree()}", flush=True)
    hyperplane_sum = jacobian(0)
    for factor, multiplicity in line_polynomial.factor():
        assert multiplicity == 1
        hyperplane_sum += image_trace_from_factor(
            cubic_to_jacobian,
            factor,
            lambda extension, root: (root, extension(0), extension(1)),
            prime,
        )

    h92 = EllipticCurve([
        0, 0, 0,
        field(A1_target) * target_t**3 + field(A_target) * target_t**4,
        field(B1_target) * target_t**5 + field(B_target) * target_t**6 + field(B2_target) * target_t**7,
    ])
    print("H92P2FINAL|stage=target_isomorphism", flush=True)
    isomorphism = jacobian.isomorphism_to(h92)
    return isomorphism(hyperplane_sum)


if os.environ.get("H92P2_PROBE", "1") != "0":
    modular_values = []
    for target_value in range(1, int(os.environ.get("H92P2_SAMPLE_COUNT", "6")) + 1):
        probe = final_point_at(field, field(target_value))
        if os.environ.get("H92P2_RESIDUAL_ONLY") == "1":
            modular_values.append((field(target_value), probe[0], probe[1]))
        print(
            f"H92P2FINAL|prime={field.characteristic()}|t={target_value}|"
            f"point={probe}|status=PASS"
        )
    modular_output = os.environ.get("H92P2_RESIDUAL_OUTPUT")
    if modular_output:
        def interpolate(values):
            for denominator_degree in range(31):
                for numerator_degree in range(31):
                    if len(values) < numerator_degree + denominator_degree + 1:
                        continue
                    relation_space = matrix(field, [
                        [-input_value**index for index in range(numerator_degree + 1)]
                        + [output_value * input_value**index for index in range(denominator_degree + 1)]
                        for input_value, output_value in values
                    ]).right_kernel()
                    if relation_space.dimension() == 1:
                        relation = relation_space.basis()[0]
                        denominator = tuple(relation[numerator_degree + 1:])
                        scale = denominator[-1]
                        return (
                            tuple(int(value / scale) for value in relation[:numerator_degree + 1]),
                            tuple(int(value / scale) for value in denominator),
                        )
            raise RuntimeError("residual interpolation degree cap exhausted")
        residual_z = interpolate([(value[0], value[1]) for value in modular_values])
        residual_x = interpolate([(value[0], value[2]) for value in modular_values])
        Path(modular_output).write_text(json.dumps({
            "prime": int(field.characteristic()),
            "residual_z": residual_z,
            "residual_x": residual_x,
        }, sort_keys=True) + "\n")
else:
    numerator = marked_w_function.numerator()
    denominator = marked_w_function.denominator()
    print(
        "H92P2FINAL|marked_map="
        f"numerator_degree={numerator.degree()}|denominator_degree={denominator.degree()}|"
        f"numerator_factor_degrees={[factor.degree() for factor, _ in numerator.factor()]}|"
        f"denominator_factor_degrees={[factor.degree() for factor, _ in denominator.factor()]}",
        flush=True,
    )

if os.environ.get("H92P2_EXACT_RESIDUAL") == "1":
    # Lift the same degree-10 restriction calculation from a finite fiber to
    # QQ(w).  This is the torsor-trivialization gate needed before asking for
    # a global plane-cubic-to-H92 map.
    exact_u_ring = PolynomialRing(Kw, "u")
    exact_u = exact_u_ring.gen()
    exact_divisor = exact_u_ring(marked_w_function.numerator()) - Kw.gen() * exact_u_ring(marked_w_function.denominator())
    exact_x = exact_u_ring(first_x.numerator()) * exact_u_ring(first_x.denominator()).inverse_mod(exact_divisor)
    exact_x *= exact_u.inverse_mod(exact_divisor)**3
    exact_basis = [
        (power_z, power_x)
        for power_x in range(3)
        for power_z in range(11 - power_x)
    ]
    exact_evaluations = [
        (exact_u**power_z * exact_x**power_x).mod(exact_divisor)
        for power_z, power_x in exact_basis
    ]
    exact_matrix = matrix(Kw, [
        [evaluation[power] for evaluation in exact_evaluations]
        for power in range(29)
    ])
    exact_kernel = exact_matrix.right_kernel()
    print(
        "H92P2FINAL|exact_degree10_kernel="
        f"{exact_kernel.dimension()}|divisor_degree={exact_divisor.degree()}",
        flush=True,
    )
