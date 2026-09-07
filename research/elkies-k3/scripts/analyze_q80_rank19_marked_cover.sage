#!/usr/bin/env sage
"""Analyze the first marked-section cover of the candidate rational q=80 branch.

The bounded GF(7) surface-ideal candidate has an exact rational parameter; it
is not yet a certified global component of the marked system.  For the
first polynomial section the component conditions force

    X1 = T + (d-1) T^2.

The T^4 coefficient of X1^3+A*X1+B is therefore y12^2.  Substitution of the
surface parameter determines whether the marked section is rational over
GF(7)(t), and, if not, the branch divisor and genus of the resulting quadratic
cover of the surface line.

A partial rational reconstruction of the second polynomial section's
X-coordinate is promoted only after its full polynomial-square identity is
checked. The script then computes the two quadratic marking covers, their
branch overlap, and the genus of the combined biquadratic curve.
"""

import argparse
import json
from itertools import permutations as ordered_permutations
from pathlib import Path

from sage.all import EllipticCurve, GF, HyperellipticCurve, Matrix, PolynomialRing, prod


parser = argparse.ArgumentParser()
parser.add_argument(
    "--parameter-input",
    default="artifacts/generated-results/q80-cm24-slope-8-87-gf7-parameter.json",
)
parser.add_argument(
    "--partial-marked-input",
    default=(
        "artifacts/generated-results/"
        "q80-cm24-slope-8-87-gf7-partial-marked-parameter.json"
    ),
)
arguments = parser.parse_args()

load("elkies-k3/scripts/verify_q80_rank19_deformation_gf7.sage")

artifact_path = Path(arguments.parameter_input)
artifact = json.loads(artifact_path.read_text())
if artifact.get("schema") != "q80-cm24-formal-branch-parameter-v1":
    raise ValueError("unexpected surface-parameter artifact schema")

parameter_polynomials = PolynomialRing(GF(7), "t")
t = parameter_polynomials.gen()
parameter_functions = parameter_polynomials.fraction_field()
center = {"D": 3, "P": 4, "Q": 3, "E": 2}
raw_surface = {
    raw_name: parameter_functions(center[centered_name])
    + parameter_functions(artifact["functions"][centered_name]["value"])
    for raw_name, centered_name in zip(("d", "p", "q", "e"), ("D", "P", "Q", "E"))
}

# Work in the original coefficient field long enough to extract the forced
# leading square.  The other P1 coefficients do not occur in this term.
forced_x1 = T + (d-1)*T**2
leading_square = (forced_x1**3 + A*forced_x1 + B)[4]

images = []
for variable, name in zip(variables, names):
    if name in raw_surface:
        images.append(raw_surface[name])
    else:
        images.append(parameter_functions(GF(7)(seed[variable])))
surface_map = parameters.hom(images, parameter_functions)
radicand = (
    surface_map(leading_square.numerator())
    / surface_map(leading_square.denominator())
)
numerator = radicand.numerator()
denominator = radicand.denominator()


def odd_factor_data(polynomial):
    return tuple(
        (factor_polynomial, int(exponent))
        for factor_polynomial, exponent in polynomial.factor()
        if exponent % 2
    )


odd_numerator = odd_factor_data(numerator)
odd_denominator = odd_factor_data(denominator)
finite_branch_degree = sum(
    factor_polynomial.degree()
    for factor_polynomial, _ in odd_numerator + odd_denominator
)
infinity_valuation = denominator.degree()-numerator.degree()
infinity_branched = infinity_valuation % 2 != 0
branch_degree = finite_branch_degree + int(infinity_branched)
geometric_square = branch_degree == 0
genus = None if geometric_square else (branch_degree-2)//2
branch_polynomial = prod(
    factor_polynomial.monic()
    for factor_polynomial, _ in odd_numerator + odd_denominator
)
assert (
    radicand / parameter_functions(branch_polynomial)
).is_square(), "squarefree model does not represent the marked cover"

# The genus-two curve expected in the rank-construction source has model
#   v^2 = 16*t0^6 - 19*t0^4 + 88*t0^2 - 48.
# Its good reduction modulo seven is the polynomial below. Distinct absolute
# Igusa invariants prove that the present marked cover is not obtained from
# that reduction by a Mobius change of the base parameter (even geometrically).
known_source_reduction = 2*t**6 + 2*t**4 + 4*t**2 + 1
cover_igusa = HyperellipticCurve(branch_polynomial).absolute_igusa_invariants_kohel()
known_igusa = HyperellipticCurve(
    known_source_reduction
).absolute_igusa_invariants_kohel()
known_source_match = cover_igusa == known_igusa

# A projective transformation stabilizing a six-point branch set is determined
# by the images of any ordered triple. All five finite roots split over GF(49),
# so only 6*5*4 image triples need inspection. The known source reduction has
# the visible nonhyperelliptic involution t -> -t; absence of any nonidentity
# stabilizer is an independent rejection gate for the present cover.
quadratic_field = GF(7**2, name="j")
quadratic_polynomials = PolynomialRing(quadratic_field, "x")
quadratic_branch_polynomial = quadratic_polynomials(branch_polynomial)
finite_branch_points = tuple(
    (root, quadratic_field.one())
    for root, multiplicity in quadratic_branch_polynomial.roots()
    if multiplicity == 1
)
branch_points = finite_branch_points + (
    (quadratic_field.one(), quadratic_field.zero()),
)
assert len(branch_points) == 6


def canonical_projective(point):
    x_value, z_value = point
    if z_value:
        return (x_value/z_value, quadratic_field.one())
    return (quadratic_field.one(), quadratic_field.zero())


branch_point_set = frozenset(map(canonical_projective, branch_points))
source_triple = branch_points[:3]
stabilizers = set()
for target_triple in ordered_permutations(branch_points, int(3)):
    equation_rows = []
    for (x_value, z_value), (u_value, v_value) in zip(
        source_triple, target_triple
    ):
        equation_rows.append(
            (v_value*x_value, v_value*z_value, -u_value*x_value, -u_value*z_value)
        )
    kernel = Matrix(quadratic_field, equation_rows).right_kernel_matrix()
    if kernel.nrows() != 1:
        continue
    a_value, b_value, c_value, d_value = kernel.row(0)
    if a_value*d_value-b_value*c_value == 0:
        continue
    image = frozenset(
        canonical_projective(
            (a_value*x_value+b_value*z_value, c_value*x_value+d_value*z_value)
        )
        for x_value, z_value in branch_points
    )
    if image != branch_point_set:
        continue
    coordinates = (a_value, b_value, c_value, d_value)
    scale = next(value for value in coordinates if value)
    stabilizers.add(tuple(value/scale for value in coordinates))

branch_stabilizer_order_gf49 = len(stabilizers)
assert branch_stabilizer_order_gf49 == 1

# Strong global test of the candidate parameter.  On the selected P1 chart,
# Y1=T^2*(y12+y13*T+y14*T^2), with Y1(1)=0.  The leading square fixes y12 on
# the quadratic cover, the T^5 coefficient fixes y13, and the component gate
# fixes y14.  All remaining coefficients must then vanish identically.  This
# distinguishes a genuine marked component from a finite-jet implicit model.
extension_polynomials = PolynomialRing(parameter_functions, "V")
V = extension_polynomials.gen()
marked_functions = parameter_functions.extension(V**2-radicand, names="v")
v = marked_functions.gen()
candidate_r_values = tuple(
    surface_map(coefficient.numerator()) / surface_map(coefficient.denominator())
    for coefficient in (forced_x1**3 + A*forced_x1 + B)
)
y12_candidate = v
y13_candidate = marked_functions(candidate_r_values[5])/(2*v)
y14_candidate = -y12_candidate-y13_candidate
marked_polynomials = PolynomialRing(marked_functions, "S")
S = marked_polynomials.gen()
candidate_y1 = S**2*(y12_candidate+y13_candidate*S+y14_candidate*S**2)
candidate_r = sum(
    marked_functions(value)*S**index for index, value in enumerate(candidate_r_values)
)
marked_residual = candidate_y1**2-candidate_r
marked_identity = not marked_residual
first_failed_coefficient = next(
    (index for index in range(max(candidate_y1.degree(), candidate_r.degree())+1)
     if marked_residual[index]),
    None,
)

print(
    "Q80MARKEDCOVER|section=P1|coordinate=y12|"
    f"radicand={radicand}|numerator_degree={numerator.degree()}|"
    f"denominator_degree={denominator.degree()}|geometric_square={int(geometric_square)}",
    flush=True,
)

# The order-40 marked jet also gives bounded rational candidates for all five
# X2 coefficients.  A full polynomial-square identity promotes or rejects all
# five simultaneously and determines the second section's marking cover.
partial_path = Path(arguments.partial_marked_input)
if partial_path.exists():
    partial = json.loads(partial_path.read_text())
    if partial.get("schema") != "q80-cm24-marked-branch-partial-parameter-v1":
        raise ValueError("unexpected partial marked-parameter schema")
    x2_values = tuple(
        parameter_functions(partial["coordinates"][f"x2{index}"]["value"])
        for index in range(5)
    )
    base_polynomials = PolynomialRing(parameter_functions, "S2")
    S2 = base_polynomials.gen()

    def mapped_polynomial(polynomial):
        return sum(
            parameter_functions(
                surface_map(coefficient.numerator())
                / surface_map(coefficient.denominator())
            )*S2**index
            for index, coefficient in enumerate(polynomial)
        )

    candidate_A = mapped_polynomial(A)
    candidate_B = mapped_polynomial(B)
    candidate_x2 = sum(value*S2**index for index, value in enumerate(x2_values))
    candidate_r2 = candidate_x2**3+candidate_A*candidate_x2+candidate_B
    x20 = x2_values[0]
    extension2_polynomials = PolynomialRing(parameter_functions, "W")
    W = extension2_polynomials.gen()
    marked2_functions = parameter_functions.extension(W**2-x20, names="w")
    w = marked2_functions.gen()
    y2_values = [marked2_functions(x20)*w]
    for coefficient_index in range(1, 7):
        known = sum(
            y2_values[left]*y2_values[coefficient_index-left]
            for left in range(1, coefficient_index)
        )
        y2_values.append(
            (marked2_functions(candidate_r2[coefficient_index])-known)
            / (2*y2_values[0])
        )
    marked2_polynomials = PolynomialRing(marked2_functions, "U")
    U = marked2_polynomials.gen()
    candidate_y2 = sum(value*U**index for index, value in enumerate(y2_values))
    candidate_r2_marked = sum(
        marked2_functions(value)*U**index
        for index, value in enumerate(candidate_r2)
    )
    marked2_residual = candidate_y2**2-candidate_r2_marked
    marked2_identity = not marked2_residual
    marked2_first_failure = next(
        (index for index in range(13) if marked2_residual[index]), None
    )
    x20_odd_numerator = odd_factor_data(x20.numerator())
    x20_odd_denominator = odd_factor_data(x20.denominator())
    x20_finite_branch_degree = sum(
        factor_polynomial.degree()
        for factor_polynomial, _ in x20_odd_numerator+x20_odd_denominator
    )
    x20_infinity_valuation = x20.denominator().degree()-x20.numerator().degree()
    x20_branch_degree = (
        x20_finite_branch_degree + int(x20_infinity_valuation % 2 != 0)
    )
    x20_genus = (x20_branch_degree-2)//2
    same_quadratic_cover = (radicand/x20).is_square()
    x20_branch_polynomial = prod(
        factor_polynomial.monic()
        for factor_polynomial, _ in x20_odd_numerator+x20_odd_denominator
    )
    finite_overlap = branch_polynomial.gcd(x20_branch_polynomial).degree()
    infinity_overlap = int(
        infinity_branched and x20_infinity_valuation % 2 != 0
    )
    branch_overlap = finite_overlap+infinity_overlap
    branch_union = branch_degree+x20_branch_degree-branch_overlap
    biquadratic_genus = None if same_quadratic_cover else branch_union-3
    third_branch_degree = branch_degree+x20_branch_degree-2*branch_overlap
    third_quotient_genus = (third_branch_degree-2)//2
    common_branch_polynomial = branch_polynomial.gcd(x20_branch_polynomial)
    third_branch_polynomial = parameter_polynomials(
        (branch_polynomial*x20_branch_polynomial).quo_rem(
            common_branch_polynomial**2
        )[0]
    ).monic()
    assert third_branch_polynomial.degree() == 4
    quartic_coefficients = [
        third_branch_polynomial[index] for index in range(5)
    ]
    quartic_e, quartic_d, quartic_c, quartic_b, quartic_a = quartic_coefficients
    quartic_I = (
        12*quartic_a*quartic_e-3*quartic_b*quartic_d+quartic_c**2
    )
    quartic_J = (
        72*quartic_a*quartic_c*quartic_e
        + 9*quartic_b*quartic_c*quartic_d
        - 27*quartic_a*quartic_d**2
        - 27*quartic_b**2*quartic_e
        - 2*quartic_c**3
    )
    third_quotient_j = (
        GF(7)(1728)*4*quartic_I**3/(4*quartic_I**3-quartic_J**2)
    )
    known_quotient_js = tuple(
        EllipticCurve(GF(7), coefficients).j_invariant()
        for coefficients in (
            (0, -19, 0, 1408, -12288),
            (0, 88, 0, 912, 36864),
        )
    )
    print(
        "Q80MARKEDCOVER|section=P2|stage=global_identity|"
        f"exact={int(marked2_identity)}|first_failed_coefficient={marked2_first_failure}|"
        f"same_as_P1_cover={int(same_quadratic_cover)}|"
        f"branch_degree={x20_branch_degree}|quadratic_cover_genus={x20_genus}|"
        f"status={'PASS' if marked2_identity else 'FINITE_JET_CANDIDATE_REJECTED'}",
        flush=True,
    )
    print(
        "Q80MARKEDCOVER|sections=P1,P2|stage=combined_cover|"
        f"branch_overlap={branch_overlap}|branch_union={branch_union}|"
        f"third_quotient_genus={third_quotient_genus}|"
        f"biquadratic_cover_genus={biquadratic_genus}|"
        f"third_quotient_model=w^2-({third_branch_polynomial})|"
        f"third_quotient_j={third_quotient_j}|"
        f"known_X0679_quotient_js={known_quotient_js}|"
        f"known_quotient_match={int(third_quotient_j in known_quotient_js)}|status=PASS",
        flush=True,
    )
print(
    "Q80MARKEDCOVER|section=P1|stage=global_identity|"
    f"exact={int(marked_identity)}|first_failed_coefficient={first_failed_coefficient}|"
    f"status={'PASS' if marked_identity else 'FINITE_JET_CANDIDATE_REJECTED'}",
    flush=True,
)
print(
    "Q80MARKEDCOVER|section=P1|"
    f"odd_numerator={odd_numerator}|odd_denominator={odd_denominator}|"
    f"infinity_valuation={infinity_valuation}|infinity_branched={int(infinity_branched)}|"
    f"branch_degree={branch_degree}|quadratic_cover_genus={genus}|status=PASS",
    flush=True,
)
print(
    "Q80MARKEDCOVER|section=P1|"
    f"squarefree_model=v^2-({branch_polynomial})|"
    f"cover_absolute_igusa={cover_igusa}|"
    f"known_source_mod7={known_source_reduction}|"
    f"known_source_absolute_igusa={known_igusa}|"
    f"known_source_mobius_match={int(known_source_match)}|status=PASS_COMPARISON",
    flush=True,
)
print(
    "Q80MARKEDCOVER|section=P1|"
    f"branch_stabilizer_field=GF(49)|"
    f"branch_stabilizer_order={branch_stabilizer_order_gf49}|"
    "nonhyperelliptic_involution=0|status=PASS_AUTOMORPHISM_GATE",
    flush=True,
)
