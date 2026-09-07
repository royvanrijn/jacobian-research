#!/usr/bin/env python3
"""Adapt the exact mapped-fibre worker to a common-producer quadratic prime."""

import argparse
import contextlib
import hashlib
import io
import json
import re
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CORE = ROOT / "elkies-k3/scripts/sample_q80_third_q12_weierstrass_mod19_quadratic.sage"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, required=True)
parser.add_argument("--new-base", type=int, required=True)
parser.add_argument("--new-base-r", type=int, default=0)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
args.input = args.input.resolve()
args.output = args.output.resolve()


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


pencil = json.loads(args.input.read_text())
if pencil.get("status") != "PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER":
    raise ValueError("common-producer resolved pencil is not certified")
specialization = pencil["specialization"]
prime = int(specialization["prime"])
match = re.fullmatch(r"r\^2 \+ (\d+)\*r \+ (\d+)", specialization["extension_modulus"])
if match is None:
    raise ValueError("cannot parse quadratic extension modulus")
linear, constant = map(int, match.groups())

local_root = ROOT / "artifacts/local/elkies-k3"
local_root.mkdir(parents=True, exist_ok=True)
with tempfile.TemporaryDirectory(dir=local_root) as temporary_directory:
    core_output = Path(temporary_directory) / "sample.json"
    source = CORE.read_text()
    input_pattern = re.compile(
        r"INPUT = ROOT / \"artifacts/generated-results/q80-third-q12-um2-p19-resolved-pencil.json\"\n"
        r"DEFAULT_OUTPUT = ROOT / \"artifacts/generated-results/q80-third-q12-um2-p19-weierstrass-sample.json\""
    )
    replacement = f"INPUT = Path({str(args.input)!r})\nDEFAULT_OUTPUT = Path({str(core_output)!r})"
    source, replacements = input_pattern.subn(replacement, source, count=1)
    if replacements != 1:
        raise ArithmeticError("immutable sample core input contract changed")
    source = source.replace(
        '"PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_MOD19_QUADRATIC"',
        '"PASS_EXACT_RESOLVED_THIRD_Q12_PENCIL_COMMON_PRODUCER"',
        1,
    )
    if source.count("base_finite = GF(19)") != 1:
        raise ArithmeticError("immutable sample core prime contract changed")
    source = source.replace("base_finite = GF(19)", f"base_finite = GF({prime})", 1)
    modulus_old = 'finite = GF(19**2, "r", modulus=m**2 + 12 * m + 3)'
    modulus_new = f'finite = GF({prime}**2, "r", modulus=m**2 + {linear} * m + {constant})'
    if source.count(modulus_old) != 1:
        raise ArithmeticError("immutable sample core modulus contract changed")
    source = source.replace(modulus_old, modulus_new, 1)
    source = source.replace("(xi + 6) * (xi + 16) ** 2", "(xi + 6) * (xi - 3) ** 2", 1)

    # Sage's default normalization over a non-prime finite constant field
    # uses a p-th-power algorithm.  Its intermediate degrees therefore grow
    # with p (and become prohibitive already at p=61).  Singular's
    # Grauert--Remmert ``normal`` procedure returns module generators for the
    # same integral closure without that characteristic-dependent blow-up.
    # The post-processing below is the same reversed-Hermite reduction used
    # by Sage's FunctionField_integral._maximal_order_basis.
    maximal_order_block = (
        "finite_order = curve_function.maximal_order()\n"
        "integral_basis = finite_order.basis()"
    )
    fast_normalization_block = r'''
def fast_integral_basis(function_field):
    from sage.all import lcm, matrix
    from sage.libs.singular.function import lib, singular_function
    from sage.rings.function_field.hermite_form_polynomial import reversed_hermite_form

    model, from_model, unused_to_model = function_field.monic_integral_model("z_fast")
    constants = model.constant_base_field()
    rational = model.base_field()
    bivariate, (y_variable, x_variable) = PolynomialRing(
        constants, names="y_fast,x_fast", order="degrevlex"
    ).objgens()
    coefficients = model.polynomial().list()
    equation_flat_over_constants = sum(
        coefficients[index].numerator().subs(x_variable) * y_variable**index
        for index in range(len(coefficients))
    )
    flattened = PolynomialRing(
        constants.prime_subfield(), names="yy_fast,xx_fast,zz_fast"
    )
    yy_fast, xx_fast, zz_fast = flattened.gens()
    equation_flat = flattened.zero()
    for monomial in equation_flat_over_constants.monomials():
        coefficient = equation_flat_over_constants.monomial_coefficient(monomial)
        equation_flat += flattened(coefficient.polynomial("zz_fast")) * flattened(monomial)
    constant_modulus = flattened(constants.polynomial("zz_fast"))

    lib("normal.lib")
    normal = singular_function("normal")
    normalization = normal(flattened.ideal([constant_modulus, equation_flat]), "isPrim")
    if len(normalization) != 2 or len(normalization[1]) != 1:
        raise ArithmeticError("Grauert--Remmert normalization returned unexpected components")
    module_generators = list(normalization[1][0])
    if not module_generators:
        raise ArithmeticError("Grauert--Remmert normalization returned no module generators")

    unflatten = flattened.hom(
        [y_variable, x_variable, constants.gen()], bivariate
    )
    polynomials = [unflatten(value) for value in module_generators]
    x_base = rational.gen()
    y_model = model.gen()
    field_generators = []
    for value in polynomials:
        polynomial_in_y = value.polynomial(bivariate.gen(0))
        reconstructed = model.zero()
        for index in range(polynomial_in_y.degree() + 1):
            reconstructed += polynomial_in_y[index].subs(x_base) * y_model**index
        field_generators.append(reconstructed)

    inverse_denominator = ~field_generators[-1]
    spanning = []
    for value in field_generators:
        integral = inverse_denominator * value
        for unused in range(model.degree()):
            spanning.append(integral)
            integral *= y_model
    vector_space, from_vector, to_vector = model.free_module()
    spanning_vectors = [to_vector(value) for value in spanning]
    common_denominator = lcm(value.denominator() for value in spanning_vectors)
    basis_matrix = matrix(
        [
            [coefficient.numerator() for coefficient in common_denominator * value]
            for value in spanning_vectors
        ]
    )
    reversed_hermite_form(basis_matrix)
    model_basis = tuple(
        from_vector(row) / common_denominator
        for row in basis_matrix
        if not row.is_zero()
    )
    if len(model_basis) != model.degree() or model_basis[0] != 1:
        raise ArithmeticError("fast normalization did not return a normalized integral basis")
    return tuple(from_model(value) for value in model_basis)


integral_basis = fast_integral_basis(curve_function)
'''.strip()
    if source.count(maximal_order_block) != 1:
        raise ArithmeticError("immutable sample core maximal-order contract changed")
    source = source.replace(maximal_order_block, fast_normalization_block, 1)

    saved_argv = sys.argv
    sys.argv = [
        str(CORE),
        "--new-base",
        str(args.new_base),
        "--new-base-r",
        str(args.new_base_r),
        "--output",
        str(core_output),
    ]
    namespace = {"__file__": str(CORE), "__name__": "__main__"}
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(source, str(CORE), "exec"), namespace)
    finally:
        sys.argv = saved_argv
    sample = json.loads(core_output.read_text())

if sample.get("status") != "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_MOD19_QUADRATIC":
    raise ArithmeticError("adapted mapped-fibre worker did not pass")
sample["schema"] = "elkies-k3.q80-third-q12-weierstrass-sample-modp2.v2"
sample["status"] = "PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_COMMON_PRODUCER"
sample["specialization"] = {
    "u": specialization["u"],
    "prime": prime,
    "extension_modulus": specialization["extension_modulus"],
    "new_base": args.new_base % prime,
    "new_base_coefficients_1_r": [args.new_base % prime, args.new_base_r % prime],
}
sample["infinity"]["simple_branch"] = "xi=-6"
sample["infinity"]["double_branch"] = "xi=3"
sample["input"] = {"path": str(args.input.relative_to(ROOT)), "sha256": sha256(args.input)}
sample["worker"] = {
    "core": {"path": str(CORE.relative_to(ROOT)), "sha256": sha256(CORE)},
    "adapter": {
        "path": str(Path(__file__).resolve().relative_to(ROOT)),
        "sha256": sha256(Path(__file__).resolve()),
    },
}
sample["normalization"] = {
    "algorithm": "Singular normal.lib normal (Grauert--Remmert module normalization)",
    "reason": "avoids the characteristic-p power-map degree growth in Sage's non-prime-constant-field default",
    "integral_basis_rank": 3,
}
sample["claim_boundary"] = (
    "This exactly converts one smooth common-producer pencil member to its "
    "elliptic Jacobian with forward and inverse maps. It does not interpolate "
    "the generic child or assert a characteristic-zero lift."
)
sample["reproduce"] = (
    "sage -python elkies-k3/scripts/sample_q80_third_q12_weierstrass_modp2.py "
    f"--input {args.input} --new-base {args.new_base} --new-base-r {args.new_base_r} "
    f"--output {args.output}"
)
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(sample, indent=2, sort_keys=True) + "\n")
print(
    f"Q80THIRDQ12COMMONSAMPLE|prime={prime}|T={args.new_base % prime}+{args.new_base_r % prime}r|"
    "RR=1,1,2,3|maps=both|status=PASS_EXACT_THIRD_Q12_WEIERSTRASS_SAMPLE_COMMON_PRODUCER"
)
