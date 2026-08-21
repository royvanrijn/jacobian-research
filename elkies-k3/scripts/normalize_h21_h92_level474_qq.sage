#!/usr/bin/env sage -python
"""Normalize the exact H21/H92 level-474 component over QQ.

Starting from the certified degree-21 component ``C(r,s)``, this script takes
the quotient by its exact involution, lowers the resulting rational plane curve
by an explicit Cremona transformation, and parametrizes the reduced model.

The expensive step is the adjoint-ideal calculation.  It has a strict timeout
and may be cached in an untracked JSON file.  Before the rational-normal-curve
inversion, the adjoint lattice is saturated and LLL-reduced; this preserves its
QQ-span while avoiding the enormous coefficients returned by Singular.
"""

from sage.all import PowerSeriesRing, PolynomialRing, QQ, ZZ, gcd, matrix, prod

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time


FACTOR_SCHEMA = "elkies-k3.h21-h92-level474-qq-factor.v1"
ADJOINT_SCHEMA = "elkies-k3.h21-h92-level474-adjoint-cache.v1"
OUTPUT_SCHEMA = "elkies-k3.h21-h92-level474-normalization.v1"


def stage(name, **values):
    payload = "|".join(f"{key}={value}" for key, value in values.items())
    print(
        f"H21H92NORM|stage={name}" + (f"|{payload}" if payload else ""),
        flush=True,
    )


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def primitive(polynomial):
    coefficients = [ZZ(value) for value in polynomial.coefficients()]
    content = gcd(coefficients)
    result = polynomial // content
    leading = result.monomial_coefficient(max(result.monomials()))
    return -result if leading < 0 else result


def singular_polynomial(polynomial):
    return str(polynomial).replace("**", "^")


def parse_singular_polynomial(value, ring):
    names = "".join(re.escape(name) for name in ring.variable_names())
    value = re.sub(rf"([{names}])([0-9]+)", r"\1^\2", value)
    value = re.sub(rf"([0-9])([{names}])", r"\1*\2", value)
    value = re.sub(rf"([{names}])(?=[{names}])", r"\1*", value)
    return ring(value)


def build_models(record):
    plane = PolynomialRing(QQ, names=("r", "s"))
    r, s = plane.gens()
    component = sum(
        QQ(item["coefficient"]) * r ** int(item["r"]) * s ** int(item["s"])
        for item in record["factor"]["coefficients"]
    )

    inverted = component(-1 / s, -1 / r)
    if (r**13 * s**13 * inverted + component).numerator():
        raise ArithmeticError("the recorded component does not satisfy the involution")

    invariant_ring = PolynomialRing(QQ, names=("t", "m"))
    t, m = invariant_ring.gens()
    invariant_field = invariant_ring.fraction_field()
    extension_ring = PolynomialRing(invariant_field, "S")
    S = extension_ring.gen()
    remainder = extension_ring(component(invariant_field(t) * S, S)).mod(
        S**2 - invariant_field(m) * S - 1 / invariant_field(t)
    )
    coefficients = remainder.list() + [invariant_field.zero()] * (2 - len(remainder))
    quotient = gcd(
        invariant_field(coefficients[0]).numerator(),
        invariant_field(coefficients[1]).numerator(),
    )
    quotient = primitive(quotient)

    # a=t*m=r-1/s is a better invariant than m.  Then b=t-a collapses
    # the 73-term degree-20 affine equation to 37 terms of bidegree (5,8).
    a_ring = PolynomialRing(QQ, names=("t", "a"))
    ta, a = a_ring.gens()
    raw_a = a_ring(quotient(ta, a / ta) * ta**8)
    if raw_a % ta**7:
        raise ArithmeticError("unexpected invariant-coordinate denominator")
    a_model = primitive(raw_a // ta**7)
    b_ring = PolynomialRing(QQ, names=("t", "b"))
    tb, b = b_ring.gens()
    b_model = primitive(b_ring(a_model(tb, tb - b)))
    if b_model(-1, -2) or b_model.derivative(tb)(-1, -2) == 0:
        raise ArithmeticError("the CM-24 quotient point is missing or singular")

    projective = b_model.homogenize("z").parent()
    tp, bp, z = projective.gens()
    homogeneous = projective(b_model.homogenize("z"))

    # Send [1:0:0], [0:1:0], and the rational node [-1:46:1] to the
    # coordinate points, then apply [X:Y:Z] -> [YZ:XZ:XY].
    linear = homogeneous(tp - z, bp + 46 * z, z)
    cremona_pullback = linear(bp * z, tp * z, tp * bp)
    exceptional = tp**8 * bp**5 * z**2
    if cremona_pullback % exceptional:
        raise ArithmeticError("Cremona exceptional multiplicities changed")
    reduced = primitive(cremona_pullback // exceptional)
    if reduced.total_degree() != 11:
        raise ArithmeticError("the reduced quotient is not degree 11")

    return component, quotient, a_model, b_model, reduced


def compute_adjoint(reduced, timeout):
    singular = Path(os.environ["SAGE_LOCAL"]) / "bin" / "Singular"
    program = f'''\
LIB "paraplanecurves.lib";
ring R=0,(t,b,z),dp;
poly f={singular_polynomial(reduced)};
timer=1;
ideal AI=adjointIdeal(f,list(4,"rattestyes"));
"ADJOINT_DONE|"+string(timer)+"|"+string(size(AI));
int k;
for(k=1;k<=size(AI);k++)
{{
  "ADJOINT|"+string(k)+"|"+string(AI[k]);
}}
'''
    started = time.time()
    try:
        completed = subprocess.run(
            [str(singular), "-q"],
            input=program,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        raise TimeoutError(
            f"Singular adjoint computation exceeded its {timeout}s hard limit"
        ) from error
    if completed.returncode:
        print(completed.stdout, end="", flush=True)
        raise RuntimeError("Singular adjoint computation failed")
    values = {}
    for line in completed.stdout.splitlines():
        if line.startswith("ADJOINT|"):
            _, index, polynomial = line.split("|", 2)
            values[int(index)] = parse_singular_polynomial(
                polynomial, reduced.parent()
            )
    if tuple(sorted(values)) != tuple(range(1, 11)):
        print(completed.stdout, end="", flush=True)
        raise RuntimeError("Singular returned no complete degree-9 adjoint basis")
    stage("adjoint", cache="miss", seconds=f"{time.time()-started:.2f}")
    return tuple(values[index] for index in range(1, 11))


def short_adjoint_basis(adjoint):
    ring = adjoint[0].parent()
    variables = ring.gens()
    monomials = tuple(
        variables[0] ** i * variables[1] ** j * variables[2] ** (9 - i - j)
        for i in range(10)
        for j in range(10 - i)
    )
    rows = matrix(
        ZZ,
        [[ZZ(polynomial.monomial_coefficient(term)) for term in monomials] for polynomial in adjoint],
    )
    saturated = rows.row_module(ZZ).saturation().basis_matrix()
    short = saturated.LLL()
    if short.rank() != 10 or short.row_space(QQ) != rows.row_space(QQ):
        raise ArithmeticError("LLL changed the adjoint QQ-span")
    polynomials = tuple(
        sum((row[index] * term for index, term in enumerate(monomials)), ring.zero())
        for row in short.rows()
    )
    original_bits = max(abs(value).nbits() for value in rows.list() if value)
    short_bits = max(abs(value).nbits() for value in short.list() if value)
    stage(
        "adjoint_lll",
        original_max_bits=original_bits,
        short_max_bits=short_bits,
    )
    return polynomials, original_bits, short_bits


def local_branch_series(reduced, precision):
    t, _, _ = reduced.parent().gens()
    series_ring = PowerSeriesRing(QQ, "w", default_prec=precision)
    w = series_ring.gen()
    z_series = series_ring(-54) + w
    t_series = series_ring(-54)
    derivative = reduced.derivative(t)
    for _ in range(8):
        correction = reduced(t_series, 1, z_series) / derivative(
            t_series, 1, z_series
        )
        t_series = (t_series - correction).add_bigoh(precision)
        if correction.valuation() >= precision - 1:
            break
    if reduced(t_series, 1, z_series).valuation() < precision - 1:
        raise ArithmeticError("Newton series did not reach the requested precision")
    return t_series, z_series


def osculating_projection(reduced, adjoint, precision=32):
    """Recover a degree-one parameter from the complete adjoint series."""

    t_series, z_series = local_branch_series(reduced, precision)
    section_series = tuple(
        value(t_series, 1, z_series).add_bigoh(precision) for value in adjoint
    )
    jets = matrix(
        QQ,
        [[section[order] for section in section_series] for order in range(10)],
    )
    kernel_8 = jets[:8].right_kernel()
    kernel_9 = jets[:9].right_kernel()
    kernel_10 = jets.right_kernel()
    if (kernel_8.dimension(), kernel_9.dimension(), kernel_10.dimension()) != (
        2,
        1,
        0,
    ):
        raise ArithmeticError(
            "the adjoint jets do not have the O(9) osculating dimensions"
        )
    order_9_vector = kernel_9.basis()[0]
    order_8_vector = next(
        value
        for value in kernel_8.basis()
        if value not in kernel_9
    )
    ring = reduced.parent()
    order_9 = sum(
        (coefficient * value for coefficient, value in zip(order_9_vector, adjoint)),
        ring.zero(),
    )
    order_8 = sum(
        (coefficient * value for coefficient, value in zip(order_8_vector, adjoint)),
        ring.zero(),
    )
    numerator_series = order_9(t_series, 1, z_series)
    denominator_series = order_8(t_series, 1, z_series)
    if numerator_series.valuation() != 9 or denominator_series.valuation() != 8:
        raise ArithmeticError("the osculating sections have the wrong orders")
    stage(
        "osculating_projection",
        jet_rank=jets.rank(),
        numerator_order=numerator_series.valuation(),
        denominator_order=denominator_series.valuation(),
    )
    return order_9, order_8


def invert_projection_by_series(reduced, projection, precision=32):
    """Invert a birational projection using one smooth rational branch.

    Generic elimination retains every adjoint base point and is needlessly
    expensive here.  The reduced curve has the smooth rational point
    ``[-54:1:-54]``, the Cremona image of ``(t,b)=(0,-8)``.  A local inverse
    series there determines the global degree-at-most-11 rational inverse by
    Pade reconstruction, after which both defining identities are checked
    exactly.
    """

    t, b, z = reduced.parent().gens()
    point = (QQ(-54), QQ(1), QQ(-54))
    if reduced(*point) or all(value(*point) == 0 for value in reduced.gradient()):
        raise ArithmeticError("the selected inversion point is not smooth")
    numerator, denominator = projection
    t_series, z_series = local_branch_series(reduced, precision)

    parameter_series = (
        numerator(t_series, 1, z_series)
        / denominator(t_series, 1, z_series)
    ).add_bigoh(precision)
    parameter_at_point = parameter_series[0]
    centered = parameter_series - parameter_at_point
    if centered.valuation() != 1:
        raise ArithmeticError("the selected projection is ramified at the inversion point")
    inverse_series = centered.reverse(precision=precision)
    t_in_parameter = t_series(inverse_series).add_bigoh(precision)
    z_in_parameter = z_series(inverse_series).add_bigoh(precision)

    parameter_ring = PolynomialRing(QQ, "x")
    x = parameter_ring.gen()
    parameter_field = parameter_ring.fraction_field()
    if parameter_at_point:
        raise ArithmeticError("the osculating parameter should vanish at the base point")
    t_candidate = parameter_field(t_in_parameter.pade(11, 11))
    z_candidate = parameter_field(z_in_parameter.pade(11, 11))

    common_denominator = (
        t_candidate.denominator() * z_candidate.denominator()
    )
    projective = (
        t_candidate.numerator() * z_candidate.denominator(),
        common_denominator,
        z_candidate.numerator() * t_candidate.denominator(),
    )
    common = gcd(gcd(projective[0], projective[1]), projective[2])
    projective = tuple(value // common for value in projective)

    def evaluate_homogeneous(polynomial):
        return sum(
            (
                coefficient
                * prod(value**exponent for value, exponent in zip(projective, exponents))
                for exponents, coefficient in polynomial.dict().items()
            ),
            parameter_ring.zero(),
        )

    if evaluate_homogeneous(reduced):
        raise ArithmeticError("the Pade inverse does not lie on the reduced curve")
    if evaluate_homogeneous(numerator) - x * evaluate_homogeneous(denominator):
        raise ArithmeticError("the Pade inverse does not invert the osculating parameter")
    stage(
        "parametrization",
        local_parameter=parameter_at_point,
        projective_degree=max(value.degree() for value in projective),
    )
    return projective


def identify_published_model(discriminant):
    parameter_ring = discriminant.parent().ring()
    x = parameter_ring.gen()
    parameter_field = parameter_ring.fraction_field()
    published = -27 * x**6 + 198 * x**4 - 171 * x**2 + 576
    # Numerically discovered once from the two branch sextics, then promoted
    # only after the following exact square identity.  The second choice is
    # obtained by x_published -> -x_published.
    A = ZZ(164590478411323)
    C = ZZ(18542755847723)
    D = ZZ(52351461697088)
    published_x = parameter_field(A * x / (C * x + D))
    ratio = parameter_field(published(published_x) / discriminant)
    if not ratio.is_square():
        raise ArithmeticError("the recovered double cover is not the published twist")
    published_y_multiplier = ratio.sqrt()
    if published_y_multiplier**2 * discriminant != published(published_x):
        raise ArithmeticError("the published-model square identity failed")
    stage(
        "published_model",
        x_map=f"{A}*x/({C}*x+{D})",
        y_multiplier_degrees=(
            f"{published_y_multiplier.numerator().degree()},"
            f"{published_y_multiplier.denominator().degree()}"
        ),
    )
    return published, published_x, published_y_multiplier


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--factor", required=True, type=Path)
    parser.add_argument("--adjoint-cache", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--adjoint-timeout", type=int, default=60)
    arguments = parser.parse_args()

    record = json.loads(arguments.factor.read_text())
    if record.get("schema") != FACTOR_SCHEMA or record.get("status") != "PASS_CHARACTERISTIC_ZERO_FACTOR":
        raise ValueError("input is not the passing characteristic-zero factor")
    component, quotient, a_model, b_model, reduced = build_models(record)
    stage(
        "models",
        component_degree=component.total_degree(),
        quotient_degree=quotient.total_degree(),
        bidegree=f"{b_model.degree(b_model.parent().gen(0))},{b_model.degree(b_model.parent().gen(1))}",
        reduced_degree=reduced.total_degree(),
    )

    adjoint = None
    if arguments.adjoint_cache and arguments.adjoint_cache.is_file():
        cache = json.loads(arguments.adjoint_cache.read_text())
        if (
            cache.get("schema") == ADJOINT_SCHEMA
            and cache.get("reduced_equation") == str(reduced)
        ):
            adjoint = tuple(reduced.parent()(value) for value in cache["basis"])
            stage("adjoint", cache="hit", path=arguments.adjoint_cache)
    if adjoint is None:
        adjoint = compute_adjoint(reduced, arguments.adjoint_timeout)
        if arguments.adjoint_cache:
            arguments.adjoint_cache.parent.mkdir(parents=True, exist_ok=True)
            arguments.adjoint_cache.write_text(
                json.dumps(
                    {
                        "schema": ADJOINT_SCHEMA,
                        "reduced_equation": str(reduced),
                        "basis": [str(value) for value in adjoint],
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n"
            )
    short, original_bits, short_bits = short_adjoint_basis(adjoint)
    projection = osculating_projection(reduced, short)
    parameter = invert_projection_by_series(reduced, projection)

    reduced_t, reduced_b, reduced_z = parameter
    original_t_numerator = reduced_b * (reduced_z - reduced_t)
    original_b_numerator = reduced_t * (reduced_z + 46 * reduced_b)
    common_denominator = reduced_t * reduced_b
    parameter_ring = reduced_t.parent()
    original_t = parameter_ring.fraction_field()(original_t_numerator / common_denominator)
    original_b = parameter_ring.fraction_field()(original_b_numerator / common_denominator)
    original_a = original_t - original_b
    discriminant = original_a**2 + 4 * original_t

    if b_model(original_t, original_b):
        raise ArithmeticError("the recovered parameter does not satisfy the quotient")
    if reduced(reduced_t, reduced_b, reduced_z):
        raise ArithmeticError("the recovered parameter does not satisfy the reduced model")

    numerator = discriminant.numerator()
    denominator = discriminant.denominator()
    numerator_factorization = numerator.factor()
    denominator_factorization = denominator.factor()
    squarefree_numerator = prod(
        (factor for factor, exponent in numerator_factorization if exponent % 2),
        parameter_ring.one(),
    )
    squarefree_denominator = prod(
        (factor for factor, exponent in denominator_factorization if exponent % 2),
        parameter_ring.one(),
    )
    stage(
        "double_cover",
        numerator_degree=numerator.degree(),
        denominator_degree=denominator.degree(),
        squarefree_degrees=f"{squarefree_numerator.degree()},{squarefree_denominator.degree()}",
    )
    published, published_x, published_y_multiplier = identify_published_model(
        discriminant
    )

    output = {
        "schema": OUTPUT_SCHEMA,
        "status": "PASS_LEVEL474_NORMALIZATION",
        "input": {
            "factor": str(arguments.factor),
            "sha256": sha256(arguments.factor),
        },
        "involution": {
            "map": ["-1/s", "-1/r"],
            "identity": "r^13*s^13*C(-1/s,-1/r)=-C(r,s)",
            "invariants": {"t": "r/s", "m": "s-1/r", "a": "t*m", "b": "t-a"},
        },
        "models": {
            "quotient_t_m": str(quotient),
            "quotient_t_a": str(a_model),
            "quotient_t_b": str(b_model),
            "reduced_degree_11": str(reduced),
            "cremona_inverse": {
                "t": "B*(Z-T)/(T*B)",
                "b": "T*(Z+46*B)/(T*B)",
            },
            "cremona_forward": {
                "T": "(b-46*z)*z",
                "B": "(t+z)*z",
                "Z": "(t+z)*(b-46*z)",
            },
        },
        "adjoint": {
            "dimension": len(adjoint),
            "degree": 9,
            "original_max_coefficient_bits": original_bits,
            "short_max_coefficient_bits": short_bits,
            "short_basis": [str(value) for value in short],
        },
        "parameter": {
            "birational_projection": [str(value) for value in projection],
            "reduced_projective": [str(value) for value in parameter],
            "t": str(original_t),
            "b": str(original_b),
            "a": str(original_a),
        },
        "double_cover": {
            "equation": "Y^2=a^2+4*t, with Y=r+1/s on C(r,s)=0",
            "discriminant": str(discriminant),
            "numerator_factorization": str(numerator_factorization),
            "denominator_factorization": str(denominator_factorization),
            "squarefree_numerator": str(squarefree_numerator),
            "squarefree_denominator": str(squarefree_denominator),
        },
        "published_level_474": {
            "equation": f"y^2={published}",
            "x": str(published_x),
            "y": f"({published_y_multiplier})*Y",
            "exact_identity": "published_y_multiplier^2*(a^2+4*t)=published_sextic(published_x)",
            "sign_variant": "negating published_x gives the second map because the sextic is even",
        },
        "proof_boundary": (
            "This gives an exact birational normalization of the H21/H92 "
            "component to the published level-474 genus-two curve.  Transporting "
            "the resulting labeled H3 Kumar family through the K3 neighbor chain "
            "is the next construction-level step."
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    stage("artifact", path=arguments.output, sha256=sha256(arguments.output))
    stage("complete", status=output["status"])


if __name__ == "__main__":
    main()
