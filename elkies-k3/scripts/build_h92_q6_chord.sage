#!/usr/bin/env sage -python
"""Build and verify the exact H92 q=6 chord pencil.

The pinned point is the ancillary marked ``P1``.  Therefore the chord

    q = (y-y(P1))/(x-x(P1))

has its section pole at ``-P1`` and spans the neighbor pencil associated to
``D=O+(-P1)-F``.  Eliminating ``x`` gives an exact double-cover equation over
``QQ(q)``; its squarefree old-base degree is the equation-level genus gate.
"""

from sage.all import *

import argparse
import hashlib
import itertools
import json
from importlib.machinery import SourceFileLoader
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"
SECTION = ROOT / "artifacts/generated-results/elkies-k3-h92-p1-lift.json"
SECTION_SHA256 = "c323bf6346bb239934a5a2d8b1a3f4067e70e993d2e4eb32aaa30f469fca6397"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h92-q6-chord.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rational_polynomial(ring, coefficients):
    return ring([QQ(value) for value in coefficients])


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--section", type=Path, default=SECTION)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--scan-infinity-scales", action="store_true")
parser.add_argument("--search-modp", type=int, default=0)
parser.add_argument("--scale-min", type=int, default=0)
parser.add_argument("--scale-max", type=int, default=6)
parser.add_argument("--polynomial-scale", type=int, default=0)
parser.add_argument("--scan-simple-hensel", action="store_true")
parser.add_argument("--triple-diagnostic", action="store_true")
parser.add_argument("--scan-triple-jets", action="store_true")
parser.add_argument("--test-e7-depressed", action="store_true")
parser.add_argument("--test-section-jet", action="store_true")
parser.add_argument("--e8-chain", action="store_true")
parser.add_argument("--full-cusp", action="store_true")
parser.add_argument("--test-current-section-jet", action="store_true")
parser.add_argument("--test-e7-cusp-plane", action="store_true")
parser.add_argument("--smith-module", action="store_true")
parser.add_argument("--staged-popov", action="store_true")
parser.add_argument("--staged-overmodule", action="store_true")
parser.add_argument("--test-local-infinity-overmodule", action="store_true")
parser.add_argument("--exact-staged-smith", action="store_true")
parser.add_argument("--solve-e7-kernel", action="store_true")
args = parser.parse_args()

assert digest(args.section) == SECTION_SHA256
section = json.loads(args.section.read_text())
assert section["target_model"] == "h92"
assert section["exact_weierstrass_square"]

anchor = SourceFileLoader("h92_q6_anchor", str(ANCHOR)).load_module()
h92_ring, h92_formulas = anchor.parse_h92(H92)
r92, s92 = anchor.EXPECTED_H92
A1, A, B1, B, B2 = tuple(QQ(value(r92, s92)) for value in h92_formulas)

base_ring = PolynomialRing(QQ, "u")
u0 = base_ring.gen()
base_field = base_ring.fraction_field()
u = base_field.gen()
x_numerator = rational_polynomial(
    base_ring, section["x_entrance_base"]["numerator_coefficients"]
)
x_denominator = rational_polynomial(
    base_ring, section["x_entrance_base"]["denominator_coefficients"]
)
y_numerator = rational_polynomial(
    base_ring, section["y_entrance_base"]["numerator_coefficients"]
)
y_denominator = rational_polynomial(
    base_ring, section["y_entrance_base"]["denominator_coefficients"]
)
xP = base_field(x_numerator) / base_field(x_denominator)
yP = base_field(y_numerator) / base_field(y_denominator)
T = 1 / u
a = A1 * T**3 + A * T**4
b = B1 * T**5 + B * T**6 + B2 * T**7
assert yP**2 == xP**3 + a * xP + b

# Saturate only at the four smooth P.O intersection roots.  The full square
# root of the affine x-denominator is scalar*u^2*Z4; its u^2 factor is the
# E8 minimal-model scaling, not two further intersections with O.
z_polynomial = x_denominator.sqrt()
assert z_polynomial**2 == x_denominator
h_polynomial = rational_polynomial(
    base_ring, section["structured_denominator"]["Z4_coefficients"]
)
assert h_polynomial.degree() == 4
assert (z_polynomial / (u0**2 * h_polynomial)) in QQ
principal_function = base_field(h_polynomial) * yP / xP
assert gcd(principal_function.denominator(), h_polynomial) == 1
principal_part = (
    principal_function.numerator()
    * principal_function.denominator().inverse_mod(h_polynomial)
).mod(h_polynomial)

# The compensated chord is recorded symbolically as a basis (1,q).  Since
# the raw numerator cancels at P1 and not at -P1, its section pole is -P1.
q_ring = PolynomialRing(QQ, "q")
q = q_ring.gen()
parameter_field = q_ring.fraction_field()
old_base_ring = PolynomialRing(parameter_field, "u")
u_parameter = old_base_ring.gen()
old_base_field = old_base_ring.fraction_field()


def transport(value):
    numerator = old_base_ring(
        [parameter_field(coefficient) for coefficient in value.numerator().list()]
    )
    denominator = old_base_ring(
        [parameter_field(coefficient) for coefficient in value.denominator().list()]
    )
    return old_base_field(numerator) / old_base_field(denominator)


xP_parameter = transport(xP)
yP_parameter = transport(yP)
a_parameter = transport(a)
h_parameter = old_base_ring(
    [parameter_field(coefficient) for coefficient in h_polynomial.list()]
)
principal_parameter = old_base_ring(
    [parameter_field(coefficient) for coefficient in principal_part.list()]
)
e8_cusp_value = u**2 * (
    base_field(h_polynomial) * yP / xP - base_field(principal_part)
)
e8_series_ring = PowerSeriesRing(QQ, "e", default_prec=3)
e = e8_series_ring.gen()
e8_cusp_series = e8_series_ring(e8_cusp_value.numerator()(e)) / e8_series_ring(
    e8_cusp_value.denominator()(e)
)
e8_chain_jet = base_field(e8_cusp_series[0] + e8_cusp_series[1] * u)
if (
    args.staged_popov
    or args.staged_overmodule
    or args.test_local_infinity_overmodule
    or args.exact_staged_smith
    or args.solve_e7_kernel
):
    # Keep precisely the already justified finite saturation: the degree-four
    # smooth P.O factor h, followed by the depth-two repeated E8 cusp jet.
    # Unlike the full Smith kernel below, this module must not be saturated by
    # the remaining degree-18 vertical factor.
    y_projective = base_field(yP) * base_field(z_polynomial**3)
    assert y_projective.denominator() in QQ
    y_projective = base_ring(y_projective)
    one_numerator = vector(
        base_field,
        [-base_field(x_numerator * z_polynomial), base_field(z_polynomial**3), 0],
    )
    chord_numerator = vector(
        base_field,
        [-base_field(y_projective), 0, base_field(z_polynomial**3)],
    )
    staged_numerator = (
        base_field(h_polynomial) * chord_numerator
        - (base_field(principal_part) + e8_chain_jet / u**2) * one_numerator
    )
    assert all(entry.denominator() in QQ for entry in staged_numerator)
    staged_numerator = vector(base_ring, staged_numerator)
    assert all(entry % u0**2 == 0 for entry in staged_numerator)
    e8_saturated_numerator = vector(
        base_ring, [entry // u0**2 for entry in staged_numerator]
    )
    staged_module = matrix(
        base_ring,
        [vector(base_ring, one_numerator), e8_saturated_numerator],
    )

    def shifted_degree(row, shifts):
        return max(
            (-Infinity if not value else value.degree()) + shift
            for value, shift in zip(row, shifts)
        )

    print(
        "H92Q6STAGED|column_degrees="
        + ",".join(
            str(max(entry.degree() for entry in row if entry))
            for row in staged_module.rows()
        ),
        flush=True,
    )
    staged_minors = [
        staged_module[:, columns].det()
        for columns in ((0, 1), (0, 2), (1, 2))
    ]
    staged_minor_gcd = gcd(staged_minors)
    print(
        f"H92Q6STAGED|minors_gcd_degree={staged_minor_gcd.degree()}|"
        f"minors_gcd_factors={staged_minor_gcd.factor()}",
        flush=True,
    )
    if args.exact_staged_smith or args.solve_e7_kernel:
        entry_gcd = gcd(list(staged_module.list()))
        assert entry_gcd.degree() == 4
        assert (entry_gcd / h_polynomial) in QQ
        normalized_module = matrix(
            base_ring,
            [[entry // entry_gcd for entry in row] for row in staged_module.rows()],
        )
        smith_diagonal, smith_left, smith_right = normalized_module.smith_form()
        assert smith_left * normalized_module * smith_right == smith_diagonal
        expected_second = u0**4 * h_polynomial**2
        assert smith_diagonal[0, 0] in QQ and smith_diagonal[0, 0]
        assert (smith_diagonal[1, 1] / expected_second) in QQ
        partial_diagonal = matrix(base_ring, 2, 3)
        partial_diagonal[0, 0] = smith_diagonal[0, 0]
        partial_diagonal[1, 1] = smith_diagonal[1, 1] / u0**2
        e8_module = partial_diagonal * smith_right.inverse()
        assert all(entry in base_ring for entry in e8_module.list())
        e8_module = matrix(base_ring, e8_module)
        e8_minors = [
            e8_module[:, columns].det()
            for columns in ((0, 1), (0, 2), (1, 2))
        ]
        e8_minor_gcd = gcd(e8_minors)
        assert (e8_minor_gcd / (u0**2 * h_polynomial**2)) in QQ
        shifts = (0, -2, -3)
        e8_popov, e8_transform = e8_module.popov_form(
            shifts=list(shifts), transformation=True
        )
        assert e8_transform * e8_module == e8_popov
        e8_degrees = [shifted_degree(row, shifts) for row in e8_popov.rows()]
        primitive_diagonal = matrix(base_ring, 2, 3)
        primitive_diagonal[0, 0] = smith_diagonal[0, 0]
        primitive_diagonal[1, 1] = smith_diagonal[0, 0]
        primitive_module = primitive_diagonal * smith_right.inverse()
        primitive_module = matrix(base_ring, primitive_module)
        primitive_popov, primitive_transform = primitive_module.popov_form(
            shifts=list(shifts), transformation=True
        )
        primitive_degrees = [
            shifted_degree(row, shifts) for row in primitive_popov.rows()
        ]
        assert sorted(primitive_degrees) == [6, 7]
        # Track the Popov functions in the normalized staged basis (1,R/u^2).
        invertible_columns = None
        for columns in ((0, 1), (0, 2), (1, 2)):
            minor = normalized_module[:, columns]
            if minor.det():
                invertible_columns = columns
                break
        assert invertible_columns is not None
        basis_change = (
            e8_module[:, invertible_columns]
            * normalized_module[:, invertible_columns].inverse()
        )
        assert basis_change * normalized_module == e8_module
        popov_functions = e8_transform * basis_change
        primitive_basis_change = (
            primitive_module[:, invertible_columns]
            * normalized_module[:, invertible_columns].inverse()
        )
        assert primitive_basis_change * normalized_module == primitive_module
        primitive_functions = primitive_transform * primitive_basis_change
        e7_valuations = (
            (3, 2, 2),
            (6, 4, 3),
            (9, 6, 4),
            (5, 3, 2),
            (7, 5, 3),
            (5, 4, 2),
            (3, 3, 1),
        )

        def toric_row_valuation(row, valuation):
            y_weight, x_weight, t_weight = valuation
            coordinate_weights = (0, x_weight, y_weight)
            return min(
                coordinate_weight - coefficient.degree() * t_weight
                for coefficient, coordinate_weight in zip(row, coordinate_weights)
                if coefficient
            )

        denominator_valuations = tuple(
            toric_row_valuation(normalized_module.row(0), valuation)
            for valuation in e7_valuations
        )
        popov_valuations = tuple(
            tuple(
                toric_row_valuation(row, valuation) - denominator_value
                for valuation, denominator_value in zip(
                    e7_valuations, denominator_valuations
                )
            )
            for row in e8_popov.rows()
        )
        print(
            f"H92Q6EXACTSMITH|entry_gcd_degree={entry_gcd.degree()}|"
            f"smith={smith_diagonal[0,0]},{smith_diagonal[1,1]}|"
            f"residual_minor_gcd={e8_minor_gcd.factor()}|"
            f"popov_degrees={e8_degrees}|"
            f"function_coefficient_degrees="
            f"{[[(-Infinity if not value else value.numerator().degree()-value.denominator().degree()) for value in row] for row in popov_functions.rows()]}|"
            f"E7_function_valuations={popov_valuations}",
            flush=True,
        )
        if e8_degrees != [10, 13]:
            raise ArithmeticError("partial E8 Smith saturation has wrong degrees")
        # Apply the exact resolved III* length-six quotient.  The coordinate
        # specific reversals are governed by the same shifts (0,-2,-3): at
        # cutoff C the coefficient rows reverse in degrees C,C+2,C+3.
        local_ring = PolynomialRing(QQ, names=("Z", "U", "Y"), order="degrevlex")
        Z_local, U_local, Y_local = local_ring.gens()
        local_ideal = local_ring.ideal(
            (
                U_local**2,
                Z_local * U_local,
                Z_local**4,
                Z_local * Y_local,
                U_local * Y_local,
                Y_local**2,
            )
        )
        local_groebner = local_ideal.groebner_basis()
        local_basis = (
            local_ring(1),
            Z_local,
            Z_local**2,
            Z_local**3,
            U_local,
            Y_local,
        )
        c2_local = -B1 / A1
        c3_local = -(c2_local**3 + A * c2_local + B) / A1
        x_local = c2_local * Z_local**2 + c3_local * Z_local**3 + U_local

        def reversed_coefficient(polynomial, reversal_degree):
            if not polynomial:
                return local_ring.zero()
            if polynomial.degree() > reversal_degree:
                raise ArithmeticError("row exceeds the proposed reversed cutoff")
            return sum(
                QQ(value) * Z_local ** (reversal_degree - degree)
                for degree, value in enumerate(polynomial.list())
            )

        def local_residue(row, cutoff):
            value = (
                reversed_coefficient(row[0], cutoff)
                + reversed_coefficient(row[1], cutoff + 2) * x_local
                + reversed_coefficient(row[2], cutoff + 3) * Y_local
            )
            remainder = local_ring(value).reduce(local_groebner)
            return vector(
                QQ,
                [remainder.monomial_coefficient(monomial) for monomial in local_basis],
            )

        def formal_cusp_x(max_degree):
            answer = base_ring.zero()
            for degree in range(2, max_degree + 1):
                residual = (
                    answer**3
                    + (A1 * u0**3 + A * u0**4) * answer
                    + B1 * u0**5
                    + B * u0**6
                    + B2 * u0**7
                )
                coefficient = residual[degree + 3]
                answer += (-coefficient / A1) * u0**degree
            residual = (
                answer**3
                + (A1 * u0**3 + A * u0**4) * answer
                + B1 * u0**5
                + B * u0**6
                + B2 * u0**7
            )
            assert residual.valuation() >= max_degree + 4
            return answer

        denominator_row_normalized = normalized_module.row(0)
        denominator_cutoff = shifted_degree(denominator_row_normalized, (0, -2, -3))

        def common_local_lift(row, cutoff, x_chart):
            return (
                Z_local**3 * reversed_coefficient(row[0], cutoff)
                + Z_local
                * reversed_coefficient(row[1], cutoff + 2)
                * x_chart
                + reversed_coefficient(row[2], cutoff + 3) * Y_local
            )

        denominator_cusp = formal_cusp_x(denominator_cutoff + 6)
        denominator_x_chart = local_ring(
            sum(
                QQ(value) * Z_local**degree
                for degree, value in enumerate(denominator_cusp.list())
            )
        ) + U_local
        denominator_lift = common_local_lift(
            denominator_row_normalized, denominator_cutoff, denominator_x_chart
        )
        valuation_vectors = (
            vector(ZZ, (2, 2, 4, 3, 1, 2, 3)),
            vector(ZZ, (2, 4, 6, 4, 3, 3, 5)),
            vector(ZZ, (3, 5, 9, 6, 3, 5, 7)),
        )
        vZ_local, vU_local, vY_local = valuation_vectors
        denominator_lift_valuation = vector(
            ZZ,
            [
                min(
                    monomial[0] * vZ_local[index]
                    + monomial[1] * vU_local[index]
                    + monomial[2] * vY_local[index]
                    for monomial in denominator_lift.dict()
                )
                for index in range(7)
            ],
        )
        denominator_laurent_valuation = (
            denominator_lift_valuation
            - (denominator_cutoff + 3) * vZ_local
        )
        # Here the common denominator after cancelling the smooth h factor is
        # z^3*(x-xP)/h.  At the E7 place z has base order -6, h has base
        # order -4, and x-xP=U+O(Z^3) has exceptional valuation vU.  Thus its
        # exact divisorial valuation is vU-14*vZ.  This independent local
        # identity guards against silently changing the common
        # trivialization in the reversed-row calculation.
        assert denominator_laurent_valuation == vU_local - 14 * vZ_local
        print(
            f"H92Q6EXACTSMITH|denominator_cutoff={denominator_cutoff}|"
            f"denominator_laurent_valuation={tuple(denominator_laurent_valuation)}",
            flush=True,
        )

        def correct_local_remainders(
            rows,
            cutoff,
            extra_cycle=False,
            cycle_multiple=None,
            fiber_gate=False,
        ):
            """Return the exact shifted-complete-ideal obstruction.

            ``common_local_lift`` represents ``Z^(C+3) N`` for a numerator
            over the fixed denominator ``z^3*(x-xP)/h``.  Regularity of the
            rational function therefore requires membership in the shifted
            denominator cycle; the q=6 E7 correction adds ``vZ+vU``.  Since
            these lifts are linear in U and Y, complete-ideal membership is
            exactly the componentwise monomial-valuation test below; no
            reduction by the surface equation is needed.
            """
            shift = cutoff + 3
            cusp = formal_cusp_x(shift + 3)
            cusp_local = local_ring(
                sum(QQ(value) * Z_local**degree for degree, value in enumerate(cusp.list()))
            )
            x_chart = cusp_local + U_local
            target_cycle = shift * vZ_local + denominator_laurent_valuation
            if cycle_multiple is None:
                cycle_multiple = ZZ(bool(extra_cycle))
            target_cycle += ZZ(cycle_multiple) * (vZ_local + vU_local)
            if fiber_gate:
                target_cycle += vZ_local
            remainders = []
            for row in rows:
                common_lift = common_local_lift(row, cutoff, x_chart)
                remainder = local_ring.zero()
                for monomial, coefficient in local_ring(common_lift).dict().items():
                    z_degree, u_degree, y_degree = monomial
                    assert y_degree <= 1 and u_degree <= 1
                    valuation = (
                        z_degree * vZ_local
                        + u_degree * vU_local
                        + y_degree * vY_local
                    )
                    if min(valuation - target_cycle) < 0:
                        remainder += coefficient * (
                            Z_local**z_degree
                            * U_local**u_degree
                            * Y_local**y_degree
                        )
                remainders.append(remainder)
            monomials = sorted(
                set().union(*(set(remainder.monomials()) for remainder in remainders)),
                key=str,
            )
            residue_matrix = (
                matrix(
                    QQ,
                    [
                        [
                            remainder.monomial_coefficient(monomial)
                            for remainder in remainders
                        ]
                        for monomial in monomials
                    ],
                )
                if monomials
                else matrix(QQ, 0, len(rows))
            )
            return residue_matrix

        corrected_partial_data = []
        for cutoff in range(10, 19):
            ambient_rows = []
            for row, degree in zip(e8_popov.rows(), e8_degrees):
                for power in range(max(0, cutoff - degree + 1)):
                    ambient_rows.append(u0**power * row)
            residue_matrix = correct_local_remainders(ambient_rows, cutoff)
            rank = matrix(GF(19), residue_matrix).rank()
            corrected_partial_data.append(
                (cutoff, len(ambient_rows), rank, len(ambient_rows) - rank)
            )
        corrected_primitive_data = []
        for cutoff in range(7, 15):
            ambient_rows = []
            for row, degree in zip(primitive_popov.rows(), primitive_degrees):
                for power in range(max(0, cutoff - degree + 1)):
                    ambient_rows.append(u0**power * row)
            residue_matrix = correct_local_remainders(ambient_rows, cutoff)
            rank = matrix(GF(19), residue_matrix).rank()
            corrected_primitive_data.append(
                (cutoff, len(ambient_rows), rank, len(ambient_rows) - rank)
            )
        print(
            f"H92Q6EXACTSMITH|correct_common_lift_partial={corrected_partial_data}|"
            f"correct_common_lift_primitive={corrected_primitive_data}",
            flush=True,
        )

        # The corrected common-lift condition cuts the primitive closure to a
        # four-dimensional local space at cutoff 10.  Intersect that space
        # back with the deliberately retained partial E8 lattice.  This is the
        # small local-global gate: it must be done after the E7 condition, not
        # by taking the globally over-saturated primitive module as H^0(D).
        e8_invertible_columns = next(
            columns
            for columns in ((0, 1), (0, 2), (1, 2))
            if e8_module[:, columns].det()
        )
        intersection_prime = args.search_modp or 19
        intersection_field = GF(intersection_prime)
        intersection_ring = PolynomialRing(intersection_field, "u19")

        def reduce_polynomial_mod_intersection_prime(polynomial):
            return intersection_ring(
                [intersection_field(QQ(value)) for value in polynomial.list()]
            )

        intersection_rational_field = intersection_ring.fraction_field()

        def reduce_rational_mod_intersection_prime(value):
            return intersection_rational_field(
                reduce_polynomial_mod_intersection_prime(value.numerator())
            ) / intersection_rational_field(
                reduce_polynomial_mod_intersection_prime(value.denominator())
            )

        e8_module_mod_prime = matrix(
            intersection_ring,
            [
                [reduce_polynomial_mod_intersection_prime(value) for value in row]
                for row in e8_module.rows()
            ],
        )
        intersection_data = []
        first_intersection_cutoff = None
        factor_membership = []
        cutoff10_saturation_grid = []
        cutoff10_u2_functions = None
        extra_cycle_chord_functions = None
        extra_cycle_space_data = []
        extra_cycle_spaces = {}
        extra_factor_candidate_spaces = []
        first_intersection_functions = None
        cutoff12_intersection_functions = None
        cutoff13_intersection_functions = None
        for intersection_cutoff in range(7, 19):
            primitive_ambient_rows = []
            for numerator_row, degree in zip(
                primitive_popov.rows(), primitive_degrees
            ):
                for power in range(max(0, intersection_cutoff - degree + 1)):
                    primitive_ambient_rows.append(u0**power * numerator_row)
            primitive_residues = correct_local_remainders(
                primitive_ambient_rows, intersection_cutoff
            )
            primitive_residues_mod_prime = matrix(
                intersection_field, primitive_residues
            )
            primitive_rank_mod_prime = primitive_residues_mod_prime.rank()
            primitive_kernel_mod_prime = (
                identity_matrix(intersection_field, len(primitive_ambient_rows))
                if primitive_residues.nrows() == 0
                else primitive_residues_mod_prime.right_kernel().basis_matrix()
            )
            primitive_ambient_rows_mod_prime = matrix(
                intersection_ring,
                [
                    [
                        reduce_polynomial_mod_intersection_prime(value)
                        for value in row
                    ]
                    for row in primitive_ambient_rows
                ],
            )
            primitive_local_rows_mod_prime = (
                primitive_kernel_mod_prime * primitive_ambient_rows_mod_prime
            )
            primitive_e8_coefficients_mod_prime = (
                primitive_local_rows_mod_prime[:, e8_invertible_columns]
                * e8_module_mod_prime[:, e8_invertible_columns].inverse()
            )
            assert (
                primitive_e8_coefficients_mod_prime * e8_module_mod_prime
                == primitive_local_rows_mod_prime
            )
            membership_rows = []
            for coordinate in range(2):
                common_denominator = lcm(
                    [
                        value.denominator()
                        for value in primitive_e8_coefficients_mod_prime.column(
                            coordinate
                        )
                    ]
                )
                remainder_columns = [
                    intersection_ring(common_denominator * value).mod(
                        common_denominator
                    )
                    for value in primitive_e8_coefficients_mod_prime.column(
                        coordinate
                    )
                ]
                for degree in range(common_denominator.degree()):
                    membership_rows.append(
                        [remainder[degree] for remainder in remainder_columns]
                    )
            membership_matrix_mod_prime = (
                matrix(intersection_field, membership_rows)
                if membership_rows
                else matrix(
                    intersection_field, 0, primitive_kernel_mod_prime.nrows()
                )
            )
            membership_rank_mod_prime = membership_matrix_mod_prime.rank()
            finite_dimension_mod_prime = (
                membership_matrix_mod_prime.ncols() - membership_rank_mod_prime
            )
            if intersection_cutoff == 10:
                for coordinate in range(2):
                    coordinate_denominator = lcm(
                        [
                            value.denominator()
                            for value in primitive_e8_coefficients_mod_prime.column(
                                coordinate
                            )
                        ]
                    )
                    for factor, exponent in coordinate_denominator.factor():
                        modulus = factor**exponent
                        factor_rows = []
                        factor_remainders = [
                            intersection_ring(coordinate_denominator * value).mod(
                                modulus
                            )
                            for value in primitive_e8_coefficients_mod_prime.column(
                                coordinate
                            )
                        ]
                        for degree in range(modulus.degree()):
                            factor_rows.append(
                                [
                                    remainder[degree]
                                    for remainder in factor_remainders
                                ]
                            )
                        factor_matrix = matrix(intersection_field, factor_rows)
                        factor_rank = factor_matrix.rank()
                        factor_membership.append(
                            (
                                coordinate,
                                factor,
                                exponent,
                                factor_rank,
                                factor_matrix.ncols() - factor_rank,
                            )
                        )
                        if (
                            coordinate == 1
                            and factor == intersection_ring.gen()
                            and exponent == 2
                            and factor_matrix.ncols() - factor_rank == 2
                        ):
                            u2_kernel = factor_matrix.right_kernel().basis_matrix()
                            primitive_ambient_functions_mod_prime = []
                            for function_row, degree in zip(
                                primitive_functions.rows(), primitive_degrees
                            ):
                                for power in range(
                                    max(0, intersection_cutoff - degree + 1)
                                ):
                                    primitive_ambient_functions_mod_prime.append(
                                        vector(
                                            intersection_rational_field,
                                            [
                                                intersection_ring.gen()**power
                                                * reduce_rational_mod_intersection_prime(
                                                    value
                                                )
                                                for value in function_row
                                            ],
                                        )
                                    )
                            cutoff10_u2_functions = (
                                u2_kernel
                                * primitive_kernel_mod_prime
                                * matrix(
                                    intersection_rational_field,
                                    primitive_ambient_functions_mod_prime,
                                )
                            )
                residual_denominator = lcm(
                    [
                        value.denominator()
                        for value in primitive_e8_coefficients_mod_prime.column(1)
                    ]
                )
                residual_numerators = [
                    intersection_ring(residual_denominator * value)
                    for value in primitive_e8_coefficients_mod_prime.column(1)
                ]
                h_mod_intersection = reduce_polynomial_mod_intersection_prime(
                    h_polynomial
                ).monic()
                for u_exponent in range(3):
                    for h_exponent in range(3):
                        modulus = (
                            intersection_ring.gen()**u_exponent
                            * h_mod_intersection**h_exponent
                        )
                        if modulus.degree() == 0:
                            condition_matrix = matrix(
                                intersection_field,
                                0,
                                primitive_kernel_mod_prime.nrows(),
                            )
                        else:
                            condition_remainders = [
                                numerator.mod(modulus)
                                for numerator in residual_numerators
                            ]
                            condition_matrix = matrix(
                                intersection_field,
                                [
                                    [
                                        remainder[degree]
                                        for remainder in condition_remainders
                                    ]
                                    for degree in range(modulus.degree())
                                ],
                            )
                        condition_rank = condition_matrix.rank()
                        cutoff10_saturation_grid.append(
                            (
                                u_exponent,
                                h_exponent,
                                condition_rank,
                                condition_matrix.ncols() - condition_rank,
                            )
                        )
                extra_cycle_residues = correct_local_remainders(
                    primitive_ambient_rows,
                    intersection_cutoff,
                    extra_cycle=True,
                )
                extra_cycle_kernel = matrix(
                    intersection_field, extra_cycle_residues
                ).right_kernel().basis_matrix()
                assert extra_cycle_kernel.nrows() == 4
                extra_cycle_local_rows = (
                    extra_cycle_kernel * primitive_ambient_rows_mod_prime
                )
                extra_cycle_e8_coefficients = (
                    extra_cycle_local_rows[:, e8_invertible_columns]
                    * e8_module_mod_prime[:, e8_invertible_columns].inverse()
                )
                extra_cycle_denominator = lcm(
                    [
                        value.denominator()
                        for value in extra_cycle_e8_coefficients.column(1)
                    ]
                )
                extra_cycle_numerators = [
                    intersection_ring(extra_cycle_denominator * value)
                    for value in extra_cycle_e8_coefficients.column(1)
                ]
                u2_modulus = intersection_ring.gen()**2
                extra_cycle_u2_matrix = matrix(
                    intersection_field,
                    [
                        [
                            numerator.mod(u2_modulus)[degree]
                            for numerator in extra_cycle_numerators
                        ]
                        for degree in range(2)
                    ],
                )
                extra_cycle_u2_kernel = (
                    extra_cycle_u2_matrix.right_kernel().basis_matrix()
                )
                assert extra_cycle_u2_kernel.nrows() == 2
                extra_cycle_ambient_functions = []
                for function_row, degree in zip(
                    primitive_functions.rows(), primitive_degrees
                ):
                    for power in range(
                        max(0, intersection_cutoff - degree + 1)
                    ):
                        extra_cycle_ambient_functions.append(
                            vector(
                                intersection_rational_field,
                                [
                                    intersection_ring.gen()**power
                                    * reduce_rational_mod_intersection_prime(value)
                                    for value in function_row
                                ],
                            )
                        )
                extra_cycle_function_matrix = matrix(
                    intersection_rational_field,
                    extra_cycle_ambient_functions,
                )
                finite_factor_gates = []
                for coordinate in range(2):
                    coordinate_denominator = lcm(
                        [
                            value.denominator()
                            for value in extra_cycle_e8_coefficients.column(
                                coordinate
                            )
                        ]
                    )
                    coordinate_numerators = [
                        intersection_ring(coordinate_denominator * value)
                        for value in extra_cycle_e8_coefficients.column(coordinate)
                    ]
                    for factor, exponent in coordinate_denominator.factor():
                        modulus = factor**exponent
                        gate = matrix(
                            intersection_field,
                            [
                                [
                                    numerator.mod(modulus)[degree]
                                    for numerator in coordinate_numerators
                                ]
                                for degree in range(modulus.degree())
                            ],
                        )
                        finite_factor_gates.append(
                            (f"c{coordinate}:{factor}^{exponent}", gate)
                        )
                for mask in range(1 << len(finite_factor_gates)):
                    selected_rows = []
                    selected_labels = []
                    for gate_index, (label, gate) in enumerate(
                        finite_factor_gates
                    ):
                        if (mask >> gate_index) & 1:
                            selected_labels.append(label)
                            selected_rows.extend(gate.rows())
                    combined_gate = (
                        matrix(intersection_field, selected_rows)
                        if selected_rows
                        else matrix(
                            intersection_field,
                            0,
                            extra_cycle_kernel.nrows(),
                        )
                    )
                    combined_kernel = combined_gate.right_kernel().basis_matrix()
                    if combined_kernel.nrows() != 2:
                        continue
                    candidate_space = (
                        combined_kernel
                        * extra_cycle_kernel
                        * extra_cycle_function_matrix
                    )
                    if candidate_space.rank() == 2:
                        extra_factor_candidate_spaces.append(
                            (tuple(selected_labels), candidate_space)
                        )
                extra_cycle_chord_functions = (
                    extra_cycle_u2_kernel
                    * extra_cycle_kernel
                    * extra_cycle_function_matrix
                )
                extra_cycle_spaces[intersection_cutoff] = (
                    extra_cycle_chord_functions
                )
                extra_cycle_space_data.append(
                    (
                        intersection_cutoff,
                        extra_cycle_kernel.nrows(),
                        extra_cycle_u2_matrix.rank(),
                        extra_cycle_u2_kernel.nrows(),
                        extra_cycle_chord_functions.rank(),
                    )
                )
            if 8 <= intersection_cutoff <= 14 and intersection_cutoff != 10:
                extra_residues = correct_local_remainders(
                    primitive_ambient_rows,
                    intersection_cutoff,
                    extra_cycle=True,
                )
                extra_residues_mod_prime = matrix(
                    intersection_field, extra_residues
                )
                extra_kernel = extra_residues_mod_prime.right_kernel().basis_matrix()
                extra_local_rows = extra_kernel * primitive_ambient_rows_mod_prime
                extra_coefficients = (
                    extra_local_rows[:, e8_invertible_columns]
                    * e8_module_mod_prime[:, e8_invertible_columns].inverse()
                )
                extra_denominator = lcm(
                    [value.denominator() for value in extra_coefficients.column(1)]
                )
                extra_numerators = [
                    intersection_ring(extra_denominator * value)
                    for value in extra_coefficients.column(1)
                ]
                extra_u2_modulus = intersection_ring.gen()**2
                extra_u2_matrix = matrix(
                    intersection_field,
                    [
                        [
                            numerator.mod(extra_u2_modulus)[degree]
                            for numerator in extra_numerators
                        ]
                        for degree in range(2)
                    ],
                )
                extra_u2_kernel = extra_u2_matrix.right_kernel().basis_matrix()
                extra_ambient_functions = []
                for function_row, degree in zip(
                    primitive_functions.rows(), primitive_degrees
                ):
                    for power in range(
                        max(0, intersection_cutoff - degree + 1)
                    ):
                        extra_ambient_functions.append(
                            vector(
                                intersection_rational_field,
                                [
                                    intersection_ring.gen()**power
                                    * reduce_rational_mod_intersection_prime(value)
                                    for value in function_row
                                ],
                            )
                        )
                extra_function_matrix = matrix(
                    intersection_rational_field, extra_ambient_functions
                )
                cutoff_factor_gates = []
                for coordinate in range(2):
                    coordinate_denominator = lcm(
                        [
                            value.denominator()
                            for value in extra_coefficients.column(coordinate)
                        ]
                    )
                    coordinate_numerators = [
                        intersection_ring(coordinate_denominator * value)
                        for value in extra_coefficients.column(coordinate)
                    ]
                    for factor, exponent in coordinate_denominator.factor():
                        modulus = factor**exponent
                        gate = matrix(
                            intersection_field,
                            [
                                [
                                    numerator.mod(modulus)[degree]
                                    for numerator in coordinate_numerators
                                ]
                                for degree in range(modulus.degree())
                            ],
                        )
                        cutoff_factor_gates.append(
                            (
                                f"C{intersection_cutoff}:c{coordinate}:"
                                f"{factor}^{exponent}",
                                gate,
                            )
                        )
                for mask in range(1 << len(cutoff_factor_gates)):
                    selected_rows = []
                    selected_labels = []
                    for gate_index, (label, gate) in enumerate(
                        cutoff_factor_gates
                    ):
                        if (mask >> gate_index) & 1:
                            selected_labels.append(label)
                            selected_rows.extend(gate.rows())
                    combined_gate = (
                        matrix(intersection_field, selected_rows)
                        if selected_rows
                        else matrix(
                            intersection_field, 0, extra_kernel.nrows()
                        )
                    )
                    combined_kernel = combined_gate.right_kernel().basis_matrix()
                    if combined_kernel.nrows() != 2:
                        continue
                    candidate_space = (
                        combined_kernel
                        * extra_kernel
                        * extra_function_matrix
                    )
                    if candidate_space.rank() == 2:
                        extra_factor_candidate_spaces.append(
                            (tuple(selected_labels), candidate_space)
                        )
                extra_functions = (
                    extra_u2_kernel
                    * extra_kernel
                    * extra_function_matrix
                )
                extra_cycle_spaces[intersection_cutoff] = extra_functions
                extra_cycle_space_data.append(
                    (
                        intersection_cutoff,
                        extra_kernel.nrows(),
                        extra_u2_matrix.rank(),
                        extra_u2_kernel.nrows(),
                        extra_functions.rank(),
                    )
                )
            intersection_data.append(
                (
                    intersection_cutoff,
                    len(primitive_ambient_rows),
                    primitive_rank_mod_prime,
                    primitive_kernel_mod_prime.nrows(),
                    membership_matrix_mod_prime.nrows(),
                    membership_rank_mod_prime,
                    finite_dimension_mod_prime,
                )
            )
            if finite_dimension_mod_prime == 2 and first_intersection_cutoff is None:
                first_intersection_cutoff = intersection_cutoff
                finite_kernel_mod_prime = (
                    membership_matrix_mod_prime.right_kernel().basis_matrix()
                )
                primitive_ambient_functions_mod_prime = []
                for function_row, degree in zip(
                    primitive_functions.rows(), primitive_degrees
                ):
                    for power in range(
                        max(0, intersection_cutoff - degree + 1)
                    ):
                        primitive_ambient_functions_mod_prime.append(
                            vector(
                                intersection_rational_field,
                                [
                                    intersection_ring.gen()**power
                                    * reduce_rational_mod_intersection_prime(value)
                                    for value in function_row
                                ],
                            )
                        )
                first_intersection_functions = (
                    finite_kernel_mod_prime
                    * primitive_kernel_mod_prime
                    * matrix(
                        intersection_rational_field,
                        primitive_ambient_functions_mod_prime,
                    )
                )
            if intersection_cutoff in (12, 13) and finite_dimension_mod_prime == 3:
                stable_kernel_mod_prime = (
                    membership_matrix_mod_prime.right_kernel().basis_matrix()
                )
                stable_ambient_functions_mod_prime = []
                for function_row, degree in zip(
                    primitive_functions.rows(), primitive_degrees
                ):
                    for power in range(
                        max(0, intersection_cutoff - degree + 1)
                    ):
                        stable_ambient_functions_mod_prime.append(
                            vector(
                                intersection_rational_field,
                                [
                                    intersection_ring.gen()**power
                                    * reduce_rational_mod_intersection_prime(value)
                                    for value in function_row
                                ],
                            )
                        )
                stable_intersection_functions = (
                    stable_kernel_mod_prime
                    * primitive_kernel_mod_prime
                    * matrix(
                        intersection_rational_field,
                        stable_ambient_functions_mod_prime,
                    )
                )
                if intersection_cutoff == 12:
                    cutoff12_intersection_functions = stable_intersection_functions
                else:
                    cutoff13_intersection_functions = stable_intersection_functions
        print(
            f"H92Q6INTERSECTION|prime={intersection_prime}|"
            f"data={intersection_data}|"
            f"first_dimension_two_cutoff={first_intersection_cutoff}|"
            f"cutoff10_factor_membership={factor_membership}|"
            f"cutoff10_saturation_grid={cutoff10_saturation_grid}|"
            f"cutoff11_function_rank="
            f"{None if first_intersection_functions is None else first_intersection_functions.rank()}|"
            f"cutoff12_function_rank="
            f"{None if cutoff12_intersection_functions is None else cutoff12_intersection_functions.rank()}|"
            f"cutoff13_function_rank="
            f"{None if cutoff13_intersection_functions is None else cutoff13_intersection_functions.rank()}|"
            f"extra_cycle_chord_rank="
            f"{None if extra_cycle_chord_functions is None else extra_cycle_chord_functions.rank()}|"
            f"extra_cycle_spaces={sorted(extra_cycle_space_data)}",
            flush=True,
        )
        goodprime_functions = None
        goodprime_mode = None
        extra_cycle_genus_scan = []
        if extra_cycle_chord_functions is not None:
            sample_xP = reduce_rational_mod_intersection_prime(xP)
            sample_yP = reduce_rational_mod_intersection_prime(yP)
            sample_a = reduce_rational_mod_intersection_prime(a)
            sample_h = intersection_rational_field(
                reduce_polynomial_mod_intersection_prime(h_polynomial)
            )
            sample_principal = intersection_rational_field(
                reduce_polynomial_mod_intersection_prime(principal_part)
            )
            sample_e8_jet = reduce_rational_mod_intersection_prime(
                e8_chain_jet
            )

            def specialized_odd_degree(candidate_functions, parameter_value):
                (a0_value, b0_value), (a1_value, b1_value) = (
                    candidate_functions.rows()
                )
                denominator_value = (
                    intersection_field(parameter_value) * b0_value - b1_value
                )
                if not denominator_value:
                    return None
                s_value = (
                    a1_value
                    - intersection_field(parameter_value) * a0_value
                ) / denominator_value
                r_value = intersection_rational_field(
                    intersection_ring.gen()**2
                ) * s_value
                raw_value = (
                    r_value
                    + sample_e8_jet / intersection_ring.gen()**2
                    + sample_principal
                ) / sample_h
                cover_value = (
                    raw_value**4
                    - 6 * sample_xP * raw_value**2
                    + 8 * sample_yP * raw_value
                    - 3 * sample_xP**2
                    - 4 * sample_a
                )
                parity = {}
                for polynomial_value in (
                    cover_value.numerator(),
                    cover_value.denominator(),
                ):
                    for factor, exponent in polynomial_value.factor():
                        monic_factor = factor.monic()
                        parity[monic_factor] = (
                            parity.get(monic_factor, 0) + exponent
                        ) % 2
                return sum(
                    factor.degree() for factor, odd in parity.items() if odd
                )

            u4_diagonal = matrix(base_ring, 2, 3)
            u4_diagonal[0, 0] = smith_diagonal[0, 0]
            u4_diagonal[1, 1] = smith_diagonal[1, 1] / u0**4
            u4_module = matrix(base_ring, u4_diagonal * smith_right.inverse())
            u4_popov, u4_transform = u4_module.popov_form(
                shifts=list(shifts), transformation=True
            )
            u4_degrees = [shifted_degree(row, shifts) for row in u4_popov.rows()]
            u4_basis_change = (
                u4_module[:, invertible_columns]
                * normalized_module[:, invertible_columns].inverse()
            )
            assert u4_basis_change * normalized_module == u4_module
            u4_functions = u4_transform * u4_basis_change
            u4_fiber_gate_data = []
            u4_fiber_gate_spaces = {}
            for u4_cutoff in range(min(u4_degrees), max(u4_degrees) + 9):
                u4_rows = []
                u4_function_rows = []
                for numerator_row, function_row, degree in zip(
                    u4_popov.rows(), u4_functions.rows(), u4_degrees
                ):
                    for power in range(max(0, u4_cutoff - degree + 1)):
                        u4_rows.append(u0**power * numerator_row)
                        u4_function_rows.append(
                            vector(
                                intersection_rational_field,
                                [
                                    intersection_ring.gen()**power
                                    * reduce_rational_mod_intersection_prime(value)
                                    for value in function_row
                                ],
                            )
                        )
                u4_residues = matrix(
                    intersection_field,
                    correct_local_remainders(
                        u4_rows, u4_cutoff, fiber_gate=True
                    ),
                )
                u4_kernel = u4_residues.right_kernel().basis_matrix()
                u4_candidate_functions = (
                    u4_kernel
                    * matrix(intersection_rational_field, u4_function_rows)
                )
                pencil_degrees = None
                if u4_candidate_functions.rank() == 2:
                    u4_candidate_functions = (
                        u4_candidate_functions.row_space().basis_matrix()
                    )
                    pencil_degrees = tuple(
                        specialized_odd_degree(
                            u4_candidate_functions, parameter_value
                        )
                        for parameter_value in (1, 2, 3)
                    )
                    u4_fiber_gate_spaces[u4_cutoff] = u4_candidate_functions
                u4_fiber_gate_data.append(
                    (
                        u4_cutoff,
                        len(u4_rows),
                        u4_residues.rank(),
                        u4_kernel.nrows(),
                        u4_candidate_functions.rank(),
                        pencil_degrees,
                    )
                )

            primitive_fiber_data = []
            primitive_fiber_factor_spaces = []
            for primitive_cutoff in range(7, 15):
                primitive_rows = []
                primitive_function_rows = []
                for numerator_row, function_row, degree in zip(
                    primitive_popov.rows(), primitive_functions.rows(), primitive_degrees
                ):
                    for power in range(max(0, primitive_cutoff - degree + 1)):
                        primitive_rows.append(u0**power * numerator_row)
                        primitive_function_rows.append(
                            vector(
                                intersection_rational_field,
                                [
                                    intersection_ring.gen()**power
                                    * reduce_rational_mod_intersection_prime(value)
                                    for value in function_row
                                ],
                            )
                        )
                primitive_fiber_residues = matrix(
                    intersection_field,
                    correct_local_remainders(
                        primitive_rows, primitive_cutoff, fiber_gate=True
                    ),
                )
                primitive_fiber_kernel = (
                    primitive_fiber_residues.right_kernel().basis_matrix()
                )
                primitive_rows_mod = matrix(
                    intersection_ring,
                    [
                        [
                            reduce_polynomial_mod_intersection_prime(value)
                            for value in row
                        ]
                        for row in primitive_rows
                    ],
                )
                local_rows = primitive_fiber_kernel * primitive_rows_mod
                e8_coefficients = (
                    local_rows[:, e8_invertible_columns]
                    * e8_module_mod_prime[:, e8_invertible_columns].inverse()
                )
                assert e8_coefficients * e8_module_mod_prime == local_rows
                local_function_matrix = matrix(
                    intersection_rational_field, primitive_function_rows
                )
                raw_functions = (
                    primitive_fiber_kernel * local_function_matrix
                )
                factor_gates = []
                for coordinate in range(2):
                    coordinate_denominator = lcm(
                        [
                            value.denominator()
                            for value in e8_coefficients.column(coordinate)
                        ]
                    )
                    coordinate_numerators = [
                        intersection_ring(coordinate_denominator * value)
                        for value in e8_coefficients.column(coordinate)
                    ]
                    for factor, exponent in coordinate_denominator.factor():
                        modulus = factor**exponent
                        gate = matrix(
                            intersection_field,
                            [
                                [
                                    numerator.mod(modulus)[degree]
                                    for numerator in coordinate_numerators
                                ]
                                for degree in range(modulus.degree())
                            ],
                        )
                        factor_gates.append(
                            (
                                f"PF{primitive_cutoff}:c{coordinate}:"
                                f"{factor}^{exponent}",
                                gate,
                            )
                        )
                two_space_count = 0
                for mask in range(1 << len(factor_gates)):
                    selected_rows = []
                    selected_labels = []
                    for gate_index, (label, gate) in enumerate(factor_gates):
                        if (mask >> gate_index) & 1:
                            selected_labels.append(label)
                            selected_rows.extend(gate.rows())
                    combined_gate = (
                        matrix(intersection_field, selected_rows)
                        if selected_rows
                        else matrix(
                            intersection_field,
                            0,
                            primitive_fiber_kernel.nrows(),
                        )
                    )
                    combined_kernel = combined_gate.right_kernel().basis_matrix()
                    candidate_space = combined_kernel * raw_functions
                    if candidate_space.rank() != 2:
                        continue
                    candidate_space = candidate_space.row_space().basis_matrix()
                    if candidate_space.nrows() != 2:
                        continue
                    two_space_count += 1
                    primitive_fiber_factor_spaces.append(
                        (tuple(selected_labels), candidate_space)
                    )
                primitive_fiber_data.append(
                    (
                        primitive_cutoff,
                        len(primitive_rows),
                        primitive_fiber_residues.rank(),
                        primitive_fiber_kernel.nrows(),
                        raw_functions.rank(),
                        len(factor_gates),
                        two_space_count,
                    )
                )

            partial_extra_space_data = []
            partial_extra_spaces = {}
            for partial_cutoff in range(10, 19):
                partial_rows = []
                partial_function_rows = []
                for numerator_row, function_row, degree in zip(
                    e8_popov.rows(), popov_functions.rows(), e8_degrees
                ):
                    for power in range(max(0, partial_cutoff - degree + 1)):
                        partial_rows.append(u0**power * numerator_row)
                        partial_function_rows.append(
                            vector(
                                intersection_rational_field,
                                [
                                    intersection_ring.gen()**power
                                    * reduce_rational_mod_intersection_prime(value)
                                    for value in function_row
                                ],
                            )
                        )
                partial_extra_residues = matrix(
                    intersection_field,
                    correct_local_remainders(
                        partial_rows, partial_cutoff, extra_cycle=True
                    ),
                )
                partial_extra_kernel = (
                    partial_extra_residues.right_kernel().basis_matrix()
                )
                partial_extra_functions = (
                    partial_extra_kernel
                    * matrix(
                        intersection_rational_field, partial_function_rows
                    )
                )
                pencil_degrees = None
                if (
                    partial_extra_functions.nrows() == 2
                    and partial_extra_functions.rank() == 2
                ):
                    pencil_degrees = tuple(
                        specialized_odd_degree(
                            partial_extra_functions, parameter_value
                        )
                        for parameter_value in (1, 2, 3)
                    )
                    partial_extra_spaces[partial_cutoff] = (
                        partial_extra_functions
                    )
                partial_extra_space_data.append(
                    (
                        partial_cutoff,
                        len(partial_rows),
                        partial_extra_residues.rank(),
                        partial_extra_kernel.nrows(),
                        partial_extra_functions.rank(),
                        pencil_degrees,
                    )
                )

            partial_fiber_gate_data = []
            partial_fiber_gate_spaces = {}
            for partial_cutoff in range(10, 19):
                partial_rows = []
                partial_function_rows = []
                for numerator_row, function_row, degree in zip(
                    e8_popov.rows(), popov_functions.rows(), e8_degrees
                ):
                    for power in range(max(0, partial_cutoff - degree + 1)):
                        partial_rows.append(u0**power * numerator_row)
                        partial_function_rows.append(
                            vector(
                                intersection_rational_field,
                                [
                                    intersection_ring.gen()**power
                                    * reduce_rational_mod_intersection_prime(value)
                                    for value in function_row
                                ],
                            )
                        )
                fiber_residues = matrix(
                    intersection_field,
                    correct_local_remainders(
                        partial_rows, partial_cutoff, fiber_gate=True
                    ),
                )
                fiber_kernel = fiber_residues.right_kernel().basis_matrix()
                fiber_functions = (
                    fiber_kernel
                    * matrix(
                        intersection_rational_field, partial_function_rows
                    )
                )
                pencil_degrees = None
                if fiber_functions.rank() == 2:
                    fiber_functions = fiber_functions.row_space().basis_matrix()
                    pencil_degrees = tuple(
                        specialized_odd_degree(fiber_functions, parameter_value)
                        for parameter_value in (1, 2, 3)
                    )
                    if fiber_functions.nrows() == 2:
                        partial_fiber_gate_spaces[partial_cutoff] = fiber_functions
                partial_fiber_gate_data.append(
                    (
                        partial_cutoff,
                        len(partial_rows),
                        fiber_residues.rank(),
                        fiber_kernel.nrows(),
                        fiber_functions.rank(),
                        pencil_degrees,
                    )
                )

            partial_negative_cycle_data = []
            partial_negative_cycle_spaces = {}
            for partial_cutoff in range(10, 19):
                partial_rows = []
                partial_function_rows = []
                for numerator_row, function_row, degree in zip(
                    e8_popov.rows(), popov_functions.rows(), e8_degrees
                ):
                    for power in range(max(0, partial_cutoff - degree + 1)):
                        partial_rows.append(u0**power * numerator_row)
                        partial_function_rows.append(
                            vector(
                                intersection_rational_field,
                                [
                                    intersection_ring.gen()**power
                                    * reduce_rational_mod_intersection_prime(value)
                                    for value in function_row
                                ],
                            )
                        )
                negative_residues = matrix(
                    intersection_field,
                    correct_local_remainders(
                        partial_rows, partial_cutoff, cycle_multiple=-1
                    ),
                )
                negative_kernel = negative_residues.right_kernel().basis_matrix()
                negative_functions = (
                    negative_kernel
                    * matrix(
                        intersection_rational_field, partial_function_rows
                    )
                )
                pencil_degrees = None
                if negative_functions.rank() == 2:
                    negative_functions = negative_functions.row_space().basis_matrix()
                    pencil_degrees = tuple(
                        specialized_odd_degree(negative_functions, parameter_value)
                        for parameter_value in (1, 2, 3)
                    )
                    partial_negative_cycle_spaces[partial_cutoff] = negative_functions
                partial_negative_cycle_data.append(
                    (
                        partial_cutoff,
                        len(partial_rows),
                        negative_residues.rank(),
                        negative_kernel.nrows(),
                        negative_functions.rank(),
                        pencil_degrees,
                    )
                )

            regular_cutoff = 13
            regular_rows = []
            regular_function_rows = []
            for numerator_row, function_row, degree in zip(
                e8_popov.rows(), popov_functions.rows(), e8_degrees
            ):
                for power in range(max(0, regular_cutoff - degree + 1)):
                    regular_rows.append(u0**power * numerator_row)
                    regular_function_rows.append(
                        vector(
                            intersection_rational_field,
                            [
                                intersection_ring.gen()**power
                                * reduce_rational_mod_intersection_prime(value)
                                for value in function_row
                            ],
                        )
                    )
            regular_residues = matrix(
                intersection_field,
                correct_local_remainders(regular_rows, regular_cutoff),
            )
            regular_kernel = regular_residues.right_kernel().basis_matrix()
            regular_functions = (
                regular_kernel
                * matrix(intersection_rational_field, regular_function_rows)
            )
            assert regular_functions.nrows() == 3
            regular_plane_histogram = {}
            regular_plane_candidates = []
            projective_normals = [
                vector(intersection_field, (1, aa, bb))
                for aa in intersection_field
                for bb in intersection_field
            ] + [
                vector(intersection_field, (0, 1, aa))
                for aa in intersection_field
            ] + [vector(intersection_field, (0, 0, 1))]
            for normal in projective_normals:
                plane = matrix(intersection_field, [normal]).right_kernel().basis_matrix()
                candidate_functions = plane * regular_functions
                if candidate_functions.rank() != 2:
                    continue
                degrees = []
                for parameter_value in (1, 2, 3):
                    degree = specialized_odd_degree(
                        candidate_functions, parameter_value
                    )
                    degrees.append(degree)
                    if degree not in (3, 4):
                        break
                prefix = tuple(degrees)
                regular_plane_histogram[prefix] = (
                    regular_plane_histogram.get(prefix, 0) + 1
                )
                if len(degrees) == 3:
                    regular_plane_candidates.append((tuple(normal), candidate_functions))

            old_basis_0, old_basis_1, old_basis_2 = (
                cutoff13_intersection_functions.rows()
            )
            direct_space_degrees = []
            finite_factor_space_data = []
            finite_factor_genus_candidates = []
            for labels, factor_space in extra_factor_candidate_spaces:
                degrees = tuple(
                    specialized_odd_degree(factor_space, parameter_value)
                    for parameter_value in (1, 2, 3)
                )
                finite_factor_space_data.append((labels, degrees))
                if all(degree in (3, 4) for degree in degrees):
                    finite_factor_genus_candidates.append((labels, factor_space))
            primitive_fiber_space_data = []
            primitive_fiber_genus_candidates = []
            for labels, factor_space in primitive_fiber_factor_spaces:
                degrees = tuple(
                    specialized_odd_degree(factor_space, parameter_value)
                    for parameter_value in (1, 2, 3)
                )
                primitive_fiber_space_data.append((labels, degrees))
                if all(degree in (3, 4) for degree in degrees):
                    primitive_fiber_genus_candidates.append((labels, factor_space))
            if regular_plane_candidates:
                normal, regular_plane = regular_plane_candidates[0]
                goodprime_functions = regular_plane
                goodprime_mode = f"regular_cutoff13_plane_{normal}"
            if primitive_fiber_genus_candidates and goodprime_functions is None:
                labels, primitive_fiber_space = primitive_fiber_genus_candidates[0]
                goodprime_functions = primitive_fiber_space
                goodprime_mode = "primitive_fiber_gates_" + ",".join(labels)
            for u4_cutoff, u4_space in sorted(u4_fiber_gate_spaces.items()):
                if goodprime_functions is not None:
                    break
                degrees = tuple(
                    specialized_odd_degree(u4_space, parameter_value)
                    for parameter_value in (1, 2, 3)
                )
                if all(degree in (3, 4) for degree in degrees):
                    goodprime_functions = u4_space
                    goodprime_mode = f"u4_fiber_gate_cutoff_{u4_cutoff}"
                    break
            for partial_cutoff, partial_space in sorted(
                partial_fiber_gate_spaces.items()
            ):
                if goodprime_functions is not None:
                    break
                degrees = tuple(
                    specialized_odd_degree(partial_space, parameter_value)
                    for parameter_value in (1, 2, 3)
                )
                if all(degree in (3, 4) for degree in degrees):
                    goodprime_functions = partial_space
                    goodprime_mode = f"direct_fiber_gate_cutoff_{partial_cutoff}"
                    break
            for partial_cutoff, partial_space in sorted(
                partial_negative_cycle_spaces.items()
            ):
                if goodprime_functions is not None:
                    break
                degrees = tuple(
                    specialized_odd_degree(partial_space, parameter_value)
                    for parameter_value in (1, 2, 3)
                )
                if all(degree in (3, 4) for degree in degrees):
                    goodprime_functions = partial_space
                    goodprime_mode = f"direct_partial_negative_cycle_cutoff_{partial_cutoff}"
                    break
            for partial_cutoff, partial_space in sorted(
                partial_extra_spaces.items()
            ):
                degrees = tuple(
                    specialized_odd_degree(partial_space, parameter_value)
                    for parameter_value in (1, 2, 3)
                )
                if all(degree in (3, 4) for degree in degrees):
                    goodprime_functions = partial_space
                    goodprime_mode = f"direct_partial_extra_cutoff_{partial_cutoff}"
                    break
            for chord_cutoff, chord_space in sorted(extra_cycle_spaces.items()):
                if goodprime_functions is not None:
                    break
                if chord_space.nrows() != 2 or chord_space.rank() != 2:
                    continue
                degrees = tuple(
                    specialized_odd_degree(chord_space, parameter_value)
                    for parameter_value in (1, 2, 3)
                )
                direct_space_degrees.append((chord_cutoff, degrees))
                if all(degree in (3, 4) for degree in degrees):
                    goodprime_functions = chord_space
                    goodprime_mode = f"direct_extra_cycle_cutoff_{chord_cutoff}"
                    break
            if finite_factor_genus_candidates and goodprime_functions is None:
                labels, factor_space = finite_factor_genus_candidates[0]
                goodprime_functions = factor_space
                goodprime_mode = "finite_factor_gates_" + ",".join(labels)
            projective_old_rows = [
                ("infinity_infinity", old_basis_2)
            ] + [
                (
                    f"infinity_{int(coefficient)}",
                    old_basis_1 + coefficient * old_basis_2,
                )
                for coefficient in intersection_field
            ] + [
                (
                    f"{int(coefficient_1)}_{int(coefficient_2)}",
                    old_basis_0
                    + coefficient_1 * old_basis_1
                    + coefficient_2 * old_basis_2,
                )
                for coefficient_1 in intersection_field
                for coefficient_2 in intersection_field
            ]
            degree_prefix_counts = {}
            tested_pairs = 0
            for chord_cutoff, chord_space in ((10, extra_cycle_spaces[10]),):
                if goodprime_functions is not None:
                    break
                chord_basis_0, chord_basis_1 = chord_space.rows()
                projective_chord_rows = [("infinity", chord_basis_1)] + [
                    (
                        int(coefficient),
                        chord_basis_0 + coefficient * chord_basis_1,
                    )
                    for coefficient in intersection_field
                ]
                for old_label, old_row in projective_old_rows:
                    for chord_label, chord_row in projective_chord_rows:
                        candidate_functions = matrix(
                            intersection_rational_field, [old_row, chord_row]
                        )
                        if not candidate_functions.det():
                            continue
                        tested_pairs += 1
                        sample_degrees = []
                        for parameter_value in (1, 2, 3):
                            degree = specialized_odd_degree(
                                candidate_functions, parameter_value
                            )
                            sample_degrees.append(degree)
                            if degree not in (3, 4):
                                break
                        prefix = tuple(sample_degrees)
                        degree_prefix_counts[prefix] = (
                            degree_prefix_counts.get(prefix, 0) + 1
                        )
                        if len(sample_degrees) == 3:
                            extra_cycle_genus_scan.append(
                                ((chord_cutoff, old_label, chord_label), prefix)
                            )
                            goodprime_functions = candidate_functions
                            goodprime_mode = (
                                f"cutoff_{chord_cutoff}_old_{old_label}_"
                                f"plus_chord_{chord_label}"
                            )
                            break
                    if goodprime_functions is not None:
                        break
                if goodprime_functions is not None:
                    break
            print(
                f"H92Q6GENUSSCAN|prime={intersection_prime}|"
                f"tested_pairs={tested_pairs}|"
                f"degree_prefix_counts={degree_prefix_counts}|"
                f"partial_extra_spaces={partial_extra_space_data}|"
                f"partial_fiber_gate_spaces={partial_fiber_gate_data}|"
                f"u4_degrees={u4_degrees}|u4_fiber_gate_spaces={u4_fiber_gate_data}|"
                f"primitive_fiber_data={primitive_fiber_data}|"
                f"primitive_fiber_spaces={primitive_fiber_space_data}|"
                f"primitive_fiber_genus_labels="
                f"{[labels for labels, unused_space in primitive_fiber_genus_candidates]}|"
                f"partial_negative_cycle_spaces={partial_negative_cycle_data}|"
                f"regular_plane_histogram={regular_plane_histogram}|"
                f"regular_plane_candidates="
                f"{[normal for normal, unused_space in regular_plane_candidates]}|"
                f"finite_factor_spaces={finite_factor_space_data}|"
                f"finite_factor_genus_labels="
                f"{[labels for labels, unused_space in finite_factor_genus_candidates]}|"
                f"direct_space_degrees={direct_space_degrees}|"
                f"survivor_data={extra_cycle_genus_scan}|"
                f"survivor={goodprime_mode}",
                flush=True,
            )
            if goodprime_functions is None:
                raise SystemExit(0)
        if goodprime_functions is not None:
            function_determinant_mod_prime = goodprime_functions.det()
            assert function_determinant_mod_prime
            new_parameter_ring = PolynomialRing(
                intersection_field, "r19"
            )
            r19 = new_parameter_ring.gen()
            new_parameter_field = new_parameter_ring.fraction_field()
            new_old_base_ring = PolynomialRing(new_parameter_field, "u19")
            u19_parameter = new_old_base_ring.gen()
            new_old_base_field = new_old_base_ring.fraction_field()

            def transport_intersection_rational(value):
                numerator = new_old_base_ring(
                    [new_parameter_field(coefficient) for coefficient in value.numerator()]
                )
                denominator = new_old_base_ring(
                    [
                        new_parameter_field(coefficient)
                        for coefficient in value.denominator()
                    ]
                )
                return new_old_base_field(numerator) / new_old_base_field(
                    denominator
                )

            (a0_mod, b0_mod), (a1_mod, b1_mod) = (
                goodprime_functions.rows()
            )
            s_mod = (
                transport_intersection_rational(a1_mod)
                - r19 * transport_intersection_rational(a0_mod)
            ) / (
                r19 * transport_intersection_rational(b0_mod)
                - transport_intersection_rational(b1_mod)
            )
            r_compensated_mod = u19_parameter**2 * s_mod
            xP_mod = transport_intersection_rational(
                reduce_rational_mod_intersection_prime(xP)
            )
            yP_mod = transport_intersection_rational(
                reduce_rational_mod_intersection_prime(yP)
            )
            a_mod = transport_intersection_rational(
                reduce_rational_mod_intersection_prime(a)
            )
            h_mod = new_old_base_ring(
                reduce_polynomial_mod_intersection_prime(h_polynomial)
            )
            principal_mod = new_old_base_ring(
                reduce_polynomial_mod_intersection_prime(principal_part)
            )
            e8_jet_mod = transport_intersection_rational(
                reduce_rational_mod_intersection_prime(e8_chain_jet)
            )
            raw_chord_mod = (
                r_compensated_mod
                + e8_jet_mod / u19_parameter**2
                + principal_mod
            ) / h_mod
            cover_mod = (
                raw_chord_mod**4
                - 6 * xP_mod * raw_chord_mod**2
                + 8 * yP_mod * raw_chord_mod
                - 3 * xP_mod**2
                - 4 * a_mod
            )
            numerator_factorization_mod = new_old_base_ring(
                cover_mod.numerator()
            ).factor()
            denominator_factorization_mod = new_old_base_ring(
                cover_mod.denominator()
            ).factor()
            parity_mod = {}
            for factorization_mod in (
                numerator_factorization_mod,
                denominator_factorization_mod,
            ):
                for factor, exponent in factorization_mod:
                    monic_factor = factor.monic()
                    parity_mod[monic_factor] = (
                        parity_mod.get(monic_factor, 0) + exponent
                    ) % 2
            odd_degree_mod = sum(
                factor.degree() for factor, odd in parity_mod.items() if odd
            )
            squareclass_unit_mod = (
                numerator_factorization_mod.unit()
                / denominator_factorization_mod.unit()
            )
            quartic_mod = new_old_base_ring.one()
            for factor, odd in parity_mod.items():
                if odd:
                    quartic_mod *= factor
            quartic_mod = quartic_mod.monic()
            assert quartic_mod.degree() == odd_degree_mod
            quartic_coefficients_mod = list(quartic_mod.list()) + [
                new_parameter_field.zero()
            ] * 5
            e_mod, d_mod, c_mod, b_mod, aa_mod = quartic_coefficients_mod[:5]
            invariant_i_mod = (
                12 * aa_mod * e_mod - 3 * b_mod * d_mod + c_mod**2
            )
            invariant_j_mod = (
                72 * aa_mod * c_mod * e_mod
                + 9 * b_mod * c_mod * d_mod
                - 27 * aa_mod * d_mod**2
                - 27 * b_mod**2 * e_mod
                - 2 * c_mod**3
            )
            jacobian_a_mod = squareclass_unit_mod**2 * (-27 * invariant_i_mod)
            jacobian_b_mod = squareclass_unit_mod**3 * (-27 * invariant_j_mod)
            jacobian_delta_mod = squareclass_unit_mod**6 * (
                4 * invariant_i_mod**3 - invariant_j_mod**2
            )

            def kodaira_data_mod(ord_a, ord_b, ord_delta):
                if ord_delta == 1:
                    return 0, 0, 1, 1, "I1"
                if ord_a == 0 or ord_b == 0:
                    n_value = int(ord_delta)
                    return n_value - 1, n_value * (n_value - 1), n_value, n_value, f"I{n_value}"
                if ord_delta == 2:
                    return 0, 0, 1, 2, "II"
                if ord_delta == 3:
                    return 1, 2, 2, 3, "III"
                if ord_delta == 4:
                    return 2, 6, 3, 4, "IV"
                if ord_delta == 6 and ord_a >= 2 and ord_b >= 3:
                    return 4, 24, 4, 6, "I0*"
                if ord_delta >= 7 and ord_a == 2 and ord_b == 3:
                    n_value = int(ord_delta - 6)
                    rank = n_value + 4
                    return (
                        rank,
                        2 * rank * (rank - 1),
                        4,
                        n_value + 6,
                        f"I{n_value}*",
                    )
                if ord_delta == 8:
                    return 6, 72, 3, 8, "IV*"
                if ord_delta == 9:
                    return 7, 126, 2, 9, "III*"
                if ord_delta == 10:
                    return 8, 240, 1, 10, "II*"
                raise ArithmeticError((ord_a, ord_b, ord_delta))

            jacobian_factors_mod = set()
            for value in (
                jacobian_a_mod,
                jacobian_b_mod,
                jacobian_delta_mod,
            ):
                for polynomial_value in (value.numerator(), value.denominator()):
                    if not polynomial_value:
                        continue
                    jacobian_factors_mod.update(
                        factor for factor, unused_exponent in polynomial_value.factor()
                    )
            finite_signature_mod = []
            root_rank_mod = 0
            root_count_mod = 0
            root_determinant_mod = 1
            euler_mod = 0
            for factor in sorted(jacobian_factors_mod, key=str):
                raw_orders = (
                    int(jacobian_a_mod.valuation(factor)),
                    int(jacobian_b_mod.valuation(factor)),
                    int(jacobian_delta_mod.valuation(factor)),
                )
                scaling = min(
                    raw_orders[0] // 4,
                    raw_orders[1] // 6,
                    raw_orders[2] // 12,
                )
                orders = (
                    raw_orders[0] - 4 * scaling,
                    raw_orders[1] - 6 * scaling,
                    raw_orders[2] - 12 * scaling,
                )
                if orders[2] == 0:
                    continue
                rank, count, determinant, local_euler, kind = kodaira_data_mod(
                    *orders
                )
                factor_degree = int(factor.degree())
                root_rank_mod += factor_degree * rank
                root_count_mod += factor_degree * count
                root_determinant_mod *= determinant**factor_degree
                euler_mod += factor_degree * local_euler
                finite_signature_mod.append(
                    (
                        str(factor),
                        factor_degree,
                        raw_orders,
                        scaling,
                        orders,
                        kind,
                    )
                )
            infinity_raw_orders = tuple(
                int(value.denominator().degree() - value.numerator().degree())
                for value in (
                    jacobian_a_mod,
                    jacobian_b_mod,
                    jacobian_delta_mod,
                )
            )
            infinity_scaling = min(
                infinity_raw_orders[0] // 4,
                infinity_raw_orders[1] // 6,
                infinity_raw_orders[2] // 12,
            )
            infinity_orders = (
                infinity_raw_orders[0] - 4 * infinity_scaling,
                infinity_raw_orders[1] - 6 * infinity_scaling,
                infinity_raw_orders[2] - 12 * infinity_scaling,
            )
            infinity_kind = "smooth"
            if infinity_orders[2] > 0:
                rank, count, determinant, local_euler, infinity_kind = (
                    kodaira_data_mod(*infinity_orders)
                )
                root_rank_mod += rank
                root_count_mod += count
                root_determinant_mod *= determinant
                euler_mod += local_euler
            print(
                f"H92Q6GOODPRIME|prime={intersection_prime}|"
                f"mode={goodprime_mode}|"
                f"kernel_dimension=2|"
                f"function_determinant={function_determinant_mod_prime}|"
                f"functions={goodprime_functions}|"
                f"numerator_factors="
                + ",".join(
                    f"{factor.degree()}^{exponent}"
                    for factor, exponent in numerator_factorization_mod
                )
                + "|denominator_factors="
                + ",".join(
                    f"{factor.degree()}^{exponent}"
                    for factor, exponent in denominator_factorization_mod
                )
                + f"|odd_degree={odd_degree_mod}|"
                f"squareclass_unit={squareclass_unit_mod}|"
                f"finite_signature={finite_signature_mod}|"
                f"infinity_orders={infinity_orders}|"
                f"infinity_fiber={infinity_kind}|"
                f"root_data={(root_rank_mod, root_count_mod, root_determinant_mod)}|"
                f"euler={euler_mod}",
                flush=True,
            )
        raise SystemExit(0)

        e7_cutoff_data = []
        solved_kernel_functions = None
        solved_kernel_rows = None
        for cutoff in range(10, 19):
            ambient_rows = []
            for row, degree in zip(e8_popov.rows(), e8_degrees):
                for power in range(max(0, cutoff - degree + 1)):
                    ambient_rows.append(u0**power * row)
            residue_matrix = matrix(
                QQ,
                6,
                len(ambient_rows),
                [
                    value
                    for row in ambient_rows
                    for value in local_residue(row, cutoff)
                ],
            ).transpose() if False else matrix(
                QQ,
                [list(local_residue(row, cutoff)) for row in ambient_rows],
            ).transpose()
            e7_cutoff_data.append(
                (cutoff, len(ambient_rows), residue_matrix.rank(), residue_matrix.right_kernel().dimension())
            )
            if cutoff == 14:
                kernel_basis = residue_matrix.right_kernel().basis_matrix()
                assert kernel_basis.nrows() == 2
                ambient_row_matrix = matrix(base_ring, ambient_rows)
                solved_kernel_rows = kernel_basis * ambient_row_matrix
                ambient_function_rows = []
                for function_row, degree in zip(popov_functions.rows(), e8_degrees):
                    for power in range(max(0, cutoff - degree + 1)):
                        ambient_function_rows.append(
                            vector(base_field, [u**power * value for value in function_row])
                        )
                solved_kernel_functions = kernel_basis * matrix(
                    base_field, ambient_function_rows
                )
        print(f"H92Q6EXACTSMITH|E7_cutoffs={e7_cutoff_data}", flush=True)
        primitive_cutoff_data = []
        primitive_kernel_functions = None
        local_global_kernel_functions = None
        local_global_data = []
        for cutoff in range(7, 15):
            ambient_rows = []
            for row, degree in zip(primitive_popov.rows(), primitive_degrees):
                for power in range(max(0, cutoff - degree + 1)):
                    ambient_rows.append(u0**power * row)
            residue_matrix = matrix(
                QQ,
                [list(local_residue(row, cutoff)) for row in ambient_rows],
            ).transpose()
            primitive_cutoff_data.append(
                (
                    cutoff,
                    len(ambient_rows),
                    residue_matrix.rank(),
                    residue_matrix.right_kernel().dimension(),
                )
            )
            if cutoff == 9:
                kernel_basis = residue_matrix.right_kernel().basis_matrix()
                assert kernel_basis.nrows() == 2
                ambient_function_rows = []
                for function_row, degree in zip(
                    primitive_functions.rows(), primitive_degrees
                ):
                    for power in range(max(0, cutoff - degree + 1)):
                        ambient_function_rows.append(
                            vector(
                                base_field,
                                [u**power * value for value in function_row],
                            )
                        )
                primitive_kernel_functions = kernel_basis * matrix(
                    base_field, ambient_function_rows
                )
            if cutoff >= 10:
                kernel_basis = residue_matrix.right_kernel().basis_matrix()
                candidate_rows = kernel_basis * matrix(base_ring, ambient_rows)
                candidate_coefficients = (
                    candidate_rows[:, invertible_columns]
                    * e8_module[:, invertible_columns].inverse()
                )
                assert candidate_coefficients * e8_module == candidate_rows
                membership_rows = []
                for coordinate in range(2):
                    common_denominator = lcm(
                        [
                            value.denominator()
                            for value in candidate_coefficients.column(coordinate)
                        ]
                    )
                    remainder_columns = []
                    for value in candidate_coefficients.column(coordinate):
                        polynomial = base_ring(common_denominator * value)
                        remainder_columns.append(polynomial.mod(common_denominator))
                    for degree in range(common_denominator.degree()):
                        membership_rows.append(
                            [QQ(polynomial[degree]) for polynomial in remainder_columns]
                        )
                membership_matrix = matrix(QQ, membership_rows)
                membership_rank_modular = matrix(GF(19), membership_matrix).rank()
                finite_dimension = membership_matrix.ncols() - membership_rank_modular
                local_global_data.append(
                    (
                        cutoff,
                        kernel_basis.nrows(),
                        membership_rank_modular,
                        finite_dimension,
                    )
                )
                if finite_dimension == 2 and local_global_kernel_functions is None:
                    print(
                        f"H92Q6EXACTSMITH|local_global_candidate_cutoff={cutoff}|"
                        f"membership_shape={membership_matrix.dimensions()}",
                        flush=True,
                    )
                    finite_kernel = membership_matrix.right_kernel().basis_matrix()
                    assert finite_kernel.nrows() == 2
                    ambient_function_rows = []
                    for function_row, degree in zip(
                        primitive_functions.rows(), primitive_degrees
                    ):
                        for power in range(max(0, cutoff - degree + 1)):
                            ambient_function_rows.append(
                                vector(
                                    base_field,
                                    [u**power * value for value in function_row],
                                )
                            )
                    local_global_kernel_functions = (
                        finite_kernel
                        * kernel_basis
                        * matrix(base_field, ambient_function_rows)
                    )
        print(
            f"H92Q6EXACTSMITH|primitive_degrees={primitive_degrees}|"
            f"primitive_E7_cutoffs={primitive_cutoff_data}",
            flush=True,
        )
        print(f"H92Q6EXACTSMITH|local_global_cutoffs={local_global_data}", flush=True)
        assert primitive_kernel_functions is not None
        assert solved_kernel_functions is not None
        kernel_determinant = solved_kernel_functions.det()
        assert kernel_determinant
        print(
            f"H92Q6E7KERNEL|dimension=2|function_determinant={kernel_determinant}|"
            f"function_degree_profiles="
            f"{[[value.numerator().degree()-value.denominator().degree() if value else -Infinity for value in row] for row in solved_kernel_functions.rows()]}",
            flush=True,
        )
        exact_smith_popov_functions = popov_functions
        e7_kernel_functions = (
            local_global_kernel_functions
            if local_global_kernel_functions is not None
            else primitive_kernel_functions
        )
        if not args.search_modp and not args.solve_e7_kernel:
            raise SystemExit(0)
    else:
        exact_smith_popov_functions = None
    if args.test_local_infinity_overmodule:
        shifts = (0, -2, -3)
        staged_popov, staged_transform = staged_module.popov_form(
            shifts=list(shifts), transformation=True
        )
        assert staged_transform * staged_module == staged_popov
        assert staged_transform[1, 1] in QQ and staged_transform[1, 1]
        correction = base_ring(staged_transform[1, 0] / staged_transform[1, 1])
        assert correction.degree() <= 1
        # The one-step local transform at t=1/u replaces the second function
        # by (R/u^2+correction)/u.  Thus R=u^3*r-correction*u^2.
        print(
            f"H92Q6LOCALINF|correction={correction}|"
            f"transform_degrees="
            f"{[[(-1 if not value else value.degree()) for value in row] for row in staged_transform.rows()]}",
            flush=True,
        )
        # Defer the squareclass computation until bivariate_numerator exists.
        local_infinity_correction = correction
    else:
        local_infinity_correction = None
    if args.staged_overmodule:
        # Diagnostic only: test whether the remaining one-step transform is
        # represented by literal division by u in this cleared chart.  A
        # failure means the resolved III* local lattice (at t=1/u), rather
        # than a finite-place polynomial division, is required.
        equation_rows = []
        equation_rhs = []
        for entry0, entry1 in zip(staged_module.row(0), staged_module.row(1)):
            for degree in range(1):
                equation_rows.append(
                    [
                        QQ(entry0[degree - power]) if degree >= power else QQ(0)
                        for power in range(2)
                    ]
                )
                equation_rhs.append(-QQ(entry1[degree]))
        equation_matrix = matrix(QQ, equation_rows)
        equation_rhs = vector(QQ, equation_rhs)
        consistent = equation_matrix.rank() == equation_matrix.augment(
            matrix(QQ, len(equation_rhs), 1, equation_rhs)
        ).rank()
        solution = equation_matrix.solve_right(equation_rhs) if consistent else None
        print(
            f"H92Q6OVERMODULE|divisibility_equations={len(equation_rows)}|"
            f"rank={equation_matrix.rank()}|solution={solution}",
            flush=True,
        )
        if solution is None or equation_matrix.rank() != 2:
            raise ArithmeticError("literal one-step division is unavailable")
        c_value = base_ring(
            sum(solution[index] * u0**index for index in range(2))
        )
        transformed_numerator = staged_module.row(1) + c_value * staged_module.row(0)
        assert all(entry % u0 == 0 for entry in transformed_numerator)
        overmodule_generator = vector(base_ring, [entry // u0 for entry in transformed_numerator])
        overmodule = matrix(base_ring, [staged_module.row(0), overmodule_generator])
        shifts = (0, -2, -3)
        final_popov, final_transform = overmodule.popov_form(
            shifts=list(shifts), transformation=True
        )
        assert final_transform * overmodule == final_popov
        final_degrees = [shifted_degree(row, shifts) for row in final_popov.rows()]
        print(
            f"H92Q6OVERMODULE|c={c_value}|degrees={final_degrees}|"
            f"transform_degrees="
            f"{[[(-1 if not value else value.degree()) for value in row] for row in final_transform.rows()]}",
            flush=True,
        )
        raise SystemExit(0)
    for shifts in ((0, 4, 6), (0, -2, -3), (0, 0, 0)):
        popov_basis, popov_transform = staged_module.popov_form(
            shifts=list(shifts), transformation=True
        )
        assert popov_transform * staged_module == popov_basis
        denominator_degree = shifted_degree(staged_module.row(0), shifts)
        degrees = [shifted_degree(row, shifts) for row in popov_basis.rows()]
        dimensions = []
        for cutoff in range(min(degrees) - 2, denominator_degree + 3):
            dimension = sum(max(0, cutoff - degree + 1) for degree in degrees)
            if dimension <= 6:
                dimensions.append((cutoff, dimension))
        print(
            f"H92Q6STAGED|shifts={shifts}|denominator_degree={denominator_degree}|"
            f"popov_degrees={degrees}|small_dimensions={dimensions}|"
            f"transform_degrees="
            f"{[[(-1 if not value else value.degree()) for value in row] for row in popov_transform.rows()]}",
            flush=True,
        )
    if args.staged_popov:
        raise SystemExit(0)
if args.smith_module:
    y_projective = base_field(yP) * base_field(z_polynomial**3)
    assert y_projective.denominator() in QQ
    y_projective = base_ring(y_projective)
    numerator_module = matrix(
        base_ring,
        [
            [-x_numerator * z_polynomial, -y_projective],
            [z_polynomial**3, 0],
            [0, z_polynomial**3],
        ],
    )
    primitive_normal = matrix(
        base_ring, [[z_polynomial**3, x_numerator * z_polynomial, y_projective]]
    )
    smith_diagonal, smith_left, smith_right = primitive_normal.smith_form()
    assert smith_left * primitive_normal * smith_right == smith_diagonal
    assert smith_diagonal[0, 0] in QQ and smith_diagonal[0, 0]
    saturated_basis = smith_right[:, 1:]
    assert primitive_normal * saturated_basis == 0
    saturated_minors = [
        saturated_basis.matrix_from_rows(rows).det()
        for rows in ((0, 1), (0, 2), (1, 2))
    ]
    assert gcd(saturated_minors) in QQ
    popov_basis, popov_transform = saturated_basis.transpose().popov_form(
        shifts=[0, 4, 6], transformation=True
    )
    assert popov_transform * saturated_basis.transpose() == popov_basis

    def shifted_degree(row):
        return max(
            (-Infinity if not value else value.degree()) + shift
            for value, shift in zip(row, (0, 4, 6))
        )

    denominator_row = vector(
        base_ring, [-x_numerator * z_polynomial, z_polynomial**3, 0]
    )
    denominator_degree = shifted_degree(denominator_row)
    popov_degrees = [shifted_degree(row) for row in popov_basis.rows()]
    cutoff = denominator_degree - 1
    global_rows = []
    for row, degree in zip(popov_basis.rows(), popov_degrees):
        for power in range(max(0, cutoff - degree + 1)):
            global_rows.append(u0**power * row)
    print(
        f"H92Q6SMITH|smith={smith_diagonal}|popov_degrees={popov_degrees}|"
        f"denominator_degree={denominator_degree}|cutoff={cutoff}|"
        f"h0={len(global_rows)}",
        flush=True,
    )
    print(f"H92Q6SMITH|popov_basis={popov_basis}", flush=True)
    if len(global_rows) != 2:
        raise ArithmeticError("the -F weighted cutoff did not have dimension two")

    coefficients = []
    for row in global_rows:
        a_coefficient = base_field(row[1]) / base_field(z_polynomial**3)
        b_coefficient = base_field(row[2]) / base_field(z_polynomial**3)
        assert (
            base_field(row[0])
            == -a_coefficient * base_field(x_numerator * z_polynomial)
            - b_coefficient * base_field(y_projective)
        )
        coefficients.append((a_coefficient, b_coefficient))
    (a0, b0), (a1, b1) = coefficients
    determinant = a0 * b1 - a1 * b0
    assert determinant
    transformed_raw_chord = (
        transport(a1) - q * transport(a0)
    ) / (q * transport(b0) - transport(b1))
    transformed_cover = (
        transformed_raw_chord**4
        - 6 * xP_parameter * transformed_raw_chord**2
        + 8 * yP_parameter * transformed_raw_chord
        - 3 * xP_parameter**2
        - 4 * a_parameter
    )
    numerator_factorization = old_base_ring(
        transformed_cover.numerator()
    ).factor()
    denominator_factorization = old_base_ring(
        transformed_cover.denominator()
    ).factor()
    parity = {}
    representatives = {}
    for factorization_value in (numerator_factorization, denominator_factorization):
        for factor, exponent in factorization_value:
            key = factor.monic()
            parity[key] = (parity.get(key, 0) + exponent) % 2
            representatives[key] = factor
    odd_factors = [representatives[key] for key, value in parity.items() if value]
    odd_degree = sum(factor.degree() for factor in odd_factors)
    print(
        f"H92Q6SMITH|coefficients={coefficients}|determinant={determinant}|"
        f"numerator_factors="
        + ",".join(f"{f.degree()}^{e}" for f, e in numerator_factorization)
        + "|denominator_factors="
        + ",".join(f"{f.degree()}^{e}" for f, e in denominator_factorization)
        + f"|odd_degree={odd_degree}",
        flush=True,
    )
    if odd_degree not in (3, 4):
        raise ArithmeticError("Smith/Popov pencil failed the genus-one gate")
    raise SystemExit(0)
if args.full_cusp:
    raw_chord = q / h_parameter + transport(yP / xP)
elif args.e8_chain:
    # R=(u^2*q-e8_chain_jet)/u^2, hence q=R+jet/u^2.
    raw_chord = (
        q + transport(e8_chain_jet / u**2) + principal_parameter
    ) / h_parameter
else:
    raw_chord = (q + principal_parameter) / h_parameter

# Put y=yP+q*(x-xP), divide the cubic identity by x-xP, and take the
# discriminant of the resulting quadratic in x-xP.
double_cover = (
    raw_chord**4
    - 6 * xP_parameter * raw_chord**2
    + 8 * yP_parameter * raw_chord
    - 3 * xP_parameter**2
    - 4 * a_parameter
)
numerator = old_base_ring(double_cover.numerator())
denominator = old_base_ring(double_cover.denominator())
denominator_factorization = denominator.factor()
print(
    "H92Q6CHORD|denominator_factors="
    + ",".join(
        f"{factor.degree()}^{exponent}"
        for factor, exponent in denominator_factorization
    )
    + f"|denominator_unit={denominator_factorization.unit()}",
    flush=True,
)
bivariate_ring = PolynomialRing(QQ, names=("uu", "qq"))
uu, qq = bivariate_ring.gens()
bivariate_numerator = bivariate_ring.zero()
for u_degree, q_coefficient in enumerate(numerator.list()):
    q_numerator = q_coefficient.numerator()
    q_denominator = q_coefficient.denominator()
    if q_denominator not in QQ:
        raise ArithmeticError("unexpected q-dependent denominator")
    bivariate_numerator += sum(
        QQ(coefficient) / QQ(q_denominator) * uu**u_degree * qq**q_degree
        for q_degree, coefficient in enumerate(q_numerator.list())
    )
if args.test_local_infinity_overmodule:
    correction_bivariate = sum(
        QQ(value) * uu**degree
        for degree, value in enumerate(local_infinity_correction.list())
    )
    transformed_bivariate = bivariate_numerator(
        uu, uu**3 * qq - correction_bivariate * uu**2
    )
    if args.search_modp:
        finite_field = GF(args.search_modp)
        finite_ring = PolynomialRing(finite_field, names=("v", "r"))
        v_finite, r_finite = finite_ring.gens()
        finite_cover = finite_ring(bivariate_numerator)
        hits = []
        for constant in finite_field:
            candidate = finite_cover(
                v_finite, v_finite**3 * r_finite - constant * v_finite**2
            )
            parameter_ring = PolynomialRing(finite_field["r"].fraction_field(), "v")
            candidate_univariate = parameter_ring(candidate)
            odd_degree = sum(
                factor.degree()
                for factor, exponent in candidate_univariate.factor()
                if exponent % 2
            )
            if odd_degree <= 6:
                hits.append((int(constant), odd_degree))
        linear_denominator_hits = []
        for center in finite_field:
            for constant in finite_field:
                candidate = finite_cover(
                    v_finite,
                    v_finite**2
                    * ((v_finite - center) * r_finite - constant),
                )
                candidate_univariate = parameter_ring(candidate)
                odd_degree = sum(
                    factor.degree()
                    for factor, exponent in candidate_univariate.factor()
                    if exponent % 2
                )
                if odd_degree <= 6:
                    linear_denominator_hits.append(
                        (int(center), int(constant), odd_degree)
                    )
        # Cutoff 17 has the three functions (1,u,S), S=R/u^2.  A single
        # resolved III* row selects a hyperplane.  Enumerate every such
        # hyperplane and retain the full rational squareclass in u.
        staged_hyperplane_hits = []
        parameter_field_finite = finite_field["r"].fraction_field()
        r_parameter = parameter_field_finite.gen()
        u_ring_finite = PolynomialRing(parameter_field_finite, "v")
        u_parameter_finite = u_ring_finite.gen()
        projective_rows = []
        for pivot in range(3):
            for tail in itertools.product(finite_field, repeat=2 - pivot):
                row = [finite_field(0)] * 3
                row[pivot] = finite_field(1)
                for index, value in enumerate(tail, start=pivot + 1):
                    row[index] = value
                projective_rows.append(vector(finite_field, row))
        assert len(projective_rows) == args.search_modp**2 + args.search_modp + 1
        for functional in projective_rows:
            kernel_basis = matrix(finite_field, [functional]).right_kernel().basis()
            f0 = kernel_basis[0]
            f1 = kernel_basis[1]
            s_denominator = r_parameter * f0[2] - f1[2]
            if not s_denominator:
                continue
            s_value = (
                f1[0]
                + f1[1] * u_parameter_finite
                - r_parameter * (f0[0] + f0[1] * u_parameter_finite)
            ) / s_denominator
            r_value = u_parameter_finite**2 * s_value
            candidate_value = u_ring_finite.zero()
            # Evaluate the original bivariate polynomial F(u,R).
            candidate_fraction = u_ring_finite.fraction_field().zero()
            for monomial, value in finite_cover.dict().items():
                candidate_fraction += (
                    parameter_field_finite(value)
                    * u_parameter_finite**monomial[0]
                    * r_value**monomial[1]
                )
            odd_degree = 0
            factor_signature = []
            for polynomial in (
                u_ring_finite(candidate_fraction.numerator()),
                u_ring_finite(candidate_fraction.denominator()),
            ):
                for factor, exponent in polynomial.factor():
                    factor_signature.append((factor.degree(), int(exponent)))
                    if exponent % 2:
                        odd_degree += factor.degree()
            if odd_degree <= 6:
                staged_hyperplane_hits.append(
                    {
                        "functional": tuple(int(value) for value in functional),
                        "basis": tuple(
                            tuple(int(value) for value in row)
                            for row in kernel_basis
                        ),
                        "odd_degree": int(odd_degree),
                        "factors": tuple(factor_signature),
                    }
                )
        print(
            f"H92Q6LOCALINF|prime={args.search_modp}|constant_hits={hits}|"
            f"linear_denominator_hits={linear_denominator_hits}|"
            f"staged_hyperplane_hits={staged_hyperplane_hits}",
            flush=True,
        )
        raise SystemExit(0)
    transformed_cover = old_base_ring.zero()
    for monomial, value in transformed_bivariate.dict().items():
        transformed_cover += (
            parameter_field(value * q**monomial[1]) * u_parameter**monomial[0]
        )
    transformed_factorization = transformed_cover.factor()
    transformed_odd_degree = sum(
        factor.degree()
        for factor, exponent in transformed_factorization
        if exponent % 2
    )
    print(
        "H92Q6LOCALINF|factors="
        + ",".join(
            f"{factor.degree()}^{exponent}"
            for factor, exponent in transformed_factorization
        )
        + f"|unit={transformed_factorization.unit()}|"
        f"odd_degree={transformed_odd_degree}",
        flush=True,
    )
    raise SystemExit(0)
if args.solve_e7_kernel:
    if not args.e8_chain:
        raise SystemExit("--solve-e7-kernel requires --e8-chain")
    (a0_function, b0_function), (a1_function, b1_function) = tuple(
        tuple(value for value in row) for row in e7_kernel_functions.rows()
    )
    transformed_S = (
        transport(a1_function) - q * transport(a0_function)
    ) / (q * transport(b0_function) - transport(b1_function))
    transformed_R = u_parameter**2 * transformed_S
    transformed_raw_chord = (
        transformed_R
        + transport(e8_chain_jet / u**2)
        + principal_parameter
    ) / h_parameter
    transformed_cover = (
        transformed_raw_chord**4
        - 6 * xP_parameter * transformed_raw_chord**2
        + 8 * yP_parameter * transformed_raw_chord
        - 3 * xP_parameter**2
        - 4 * a_parameter
    )
    numerator_factorization = old_base_ring(transformed_cover.numerator()).factor()
    denominator_factorization = old_base_ring(transformed_cover.denominator()).factor()
    parity = {}
    for factorization_value in (numerator_factorization, denominator_factorization):
        for factor, exponent in factorization_value:
            key = factor.monic()
            parity[key] = (parity.get(key, 0) + int(exponent)) % 2
    odd_factors = [factor for factor, odd in parity.items() if odd]
    odd_degree = sum(factor.degree() for factor in odd_factors)
    squareclass_unit = (
        numerator_factorization.unit() / denominator_factorization.unit()
    )
    print(
        "H92Q6E7KERNEL|numerator_factors="
        + ",".join(
            f"{factor.degree()}^{exponent}"
            for factor, exponent in numerator_factorization
        )
        + "|denominator_factors="
        + ",".join(
            f"{factor.degree()}^{exponent}"
            for factor, exponent in denominator_factorization
        )
        + f"|odd_degree={odd_degree}|squareclass_unit={squareclass_unit}",
        flush=True,
    )
    if odd_degree not in (3, 4):
        raise ArithmeticError("resolved E7 kernel failed the genus-one gate")
    raise SystemExit(0)
if args.exact_staged_smith and args.search_modp:
    prime = ZZ(args.search_modp)
    finite_field = GF(prime)
    finite_base_ring = PolynomialRing(finite_field, "v")
    v_finite = finite_base_ring.gen()
    finite_base_field = finite_base_ring.fraction_field()

    def reduce_polynomial(polynomial):
        return finite_base_ring([finite_field(value) for value in polynomial.list()])

    def reduce_rational_function(value):
        return finite_base_field(reduce_polynomial(value.numerator())) / finite_base_field(
            reduce_polynomial(value.denominator())
        )

    function_matrix = matrix(
        finite_base_field,
        [
            [reduce_rational_function(value) for value in row]
            for row in exact_smith_popov_functions.rows()
        ],
    )
    print(f"H92Q6E7SEARCH|prime={prime}|popov_functions={function_matrix}", flush=True)
    A0, B0 = function_matrix.row(0)
    A1_function, B1_function = function_matrix.row(1)
    assert not B0 and A0 and B1_function
    finite_cover_ring = PolynomialRing(finite_field, names=("v", "R"))
    v_cover, R_cover = finite_cover_ring.gens()
    finite_cover = finite_cover_ring(bivariate_numerator)
    finite_original_denominator = finite_base_ring.zero()
    for degree, coefficient in enumerate(denominator.list()):
        assert coefficient.denominator() in QQ
        assert coefficient.numerator().degree() <= 0
        finite_original_denominator += (
            finite_field(QQ(coefficient)) * v_finite**degree
        )

    def rational_squareclass_degree(candidate):
        parity = {}
        for polynomial in (
            candidate.numerator(),
            candidate.denominator(),
            finite_original_denominator,
        ):
            for factor, exponent in finite_base_ring(polynomial).factor():
                key = factor.monic()
                parity[key] = (parity.get(key, 0) + int(exponent)) % 2
        return sum(factor.degree() for factor, odd in parity.items() if odd)

    def squareclass_degree_at(c_values, parameter_value):
        correction = sum(c_values[index] * v_finite**index for index in range(3))
        S_value = (
            parameter_value * v_finite**3 * A0
            - A1_function
            - correction * A0
        ) / B1_function
        R_value = finite_base_field(v_finite**2) * S_value
        candidate = finite_base_field.zero()
        for monomial, value in finite_cover.dict().items():
            candidate += (
                finite_field(value)
                * v_finite**monomial[0]
                * R_value**monomial[1]
            )
        return rational_squareclass_degree(candidate)

    sample_values = [finite_field(value) for value in (1, 2, 3, 5)]
    survivors = []
    centered_survivors = []
    tested = 0
    for c_values in itertools.product(finite_field, repeat=3):
        tested += 1
        degrees = []
        for parameter_value in sample_values:
            degree = squareclass_degree_at(c_values, parameter_value)
            degrees.append(degree)
            if degree not in (3, 4):
                break
        if len(degrees) == len(sample_values):
            survivors.append(
                (tuple(int(value) for value in c_values), tuple(int(value) for value in degrees))
            )
        centered_correction = sum(
            c_values[index] * v_finite ** (index + 1) for index in range(3)
        )
        centered_degrees = []
        for parameter_value in sample_values:
            S_value = (
                parameter_value * A0
                - A1_function
                - centered_correction * A0
            ) / B1_function
            R_value = finite_base_field(v_finite**2) * S_value
            candidate = finite_base_field.zero()
            for monomial, value in finite_cover.dict().items():
                candidate += (
                    finite_field(value)
                    * v_finite**monomial[0]
                    * R_value**monomial[1]
                )
            degree = rational_squareclass_degree(candidate)
            centered_degrees.append(degree)
            if degree not in (3, 4):
                break
        if len(centered_degrees) == len(sample_values):
            centered_survivors.append(
                (
                    tuple(int(value) for value in c_values),
                    tuple(int(value) for value in centered_degrees),
                )
            )
        if tested % 1000 == 0:
            print(
                f"H92Q6E7SEARCH|prime={prime}|tested={tested}|"
                f"scaled_survivors={len(survivors)}|"
                f"centered_survivors={len(centered_survivors)}",
                flush=True,
            )
    print(
        f"H92Q6E7SEARCH|prime={prime}|tested={tested}|"
        f"scaled_survivors={survivors}|centered_survivors={centered_survivors}",
        flush=True,
    )
    raise SystemExit(0)
q_coefficient_degrees = []
q_coefficient_valuations = []
for q_degree in range(5):
    coefficient = bivariate_numerator.coefficient({qq: q_degree})
    q_coefficient_degrees.append(-1 if not coefficient else coefficient.degree(uu))
    q_coefficient_valuations.append(
        -1
        if not coefficient
        else min(monomial[0] for monomial in coefficient.dict())
    )
print(
    "H92Q6CHORD|q_coefficient_u_degrees="
    + ",".join(map(str, q_coefficient_degrees)),
    flush=True,
)
if args.e8_chain:
    print(
        f"H92Q6CHORD|e8_cusp_series={e8_cusp_series}|"
        f"e8_chain_jet={e8_chain_jet}",
        flush=True,
    )
print(
    "H92Q6CHORD|q_coefficient_u_valuations="
    + ",".join(map(str, q_coefficient_valuations)),
    flush=True,
)
edge_ring_qq = PolynomialRing(QQ, "Q")
Q_edge_qq = edge_ring_qq.gen()
edge_polynomial_qq = sum(
    QQ(bivariate_numerator.monomial_coefficient(uu ** (2 * q_degree) * qq**q_degree))
    * Q_edge_qq**q_degree
    for q_degree in range(5)
)
print(
    f"H92Q6CHORD|u0_edge_qq={edge_polynomial_qq}|"
    f"factorization={edge_polynomial_qq.factor()}",
    flush=True,
)
if args.test_e7_cusp_plane:
    if not args.e8_chain:
        raise SystemExit("--test-e7-cusp-plane requires --e8-chain")
    current_cusp_value = (
        base_field(h_polynomial) * yP / xP
        - base_field(principal_part)
        - e8_chain_jet / u**2
    )
    e7_cusp_jet = (
        current_cusp_value.numerator() // current_cusp_value.denominator()
    )
    assert e7_cusp_jet.degree() <= 4
    e7_cusp_bivariate = sum(
        QQ(value) * uu**degree
        for degree, value in enumerate(e7_cusp_jet.list())
    )
    print(f"H92Q6CHORD|e7_cusp_jet={e7_cusp_jet}", flush=True)
    for division_order in range(1, 9):
        transformed_bivariate = bivariate_numerator(
            uu, uu**division_order * qq + e7_cusp_bivariate
        )
        transformed_cover = old_base_ring.zero()
        for monomial, value in transformed_bivariate.dict().items():
            transformed_cover += (
                parameter_field(value * q**monomial[1])
                * u_parameter**monomial[0]
            )
        transformed_factorization = transformed_cover.factor()
        transformed_odd_degree = sum(
            factor.degree()
            for factor, exponent in transformed_factorization
            if exponent % 2
        )
        print(
            f"H92Q6CHORD|e7_division_order={division_order}|"
            f"factors="
            + ",".join(
                f"{factor.degree()}^{exponent}"
                for factor, exponent in transformed_factorization
            )
            + f"|odd_degree={transformed_odd_degree}",
            flush=True,
        )
    raise SystemExit(0)
if args.test_current_section_jet:
    tangent_slope = (3 * xP**2 + a) / (2 * yP)
    current_section_value = (
        base_field(h_polynomial) * tangent_slope - base_field(principal_part)
    )
    if args.e8_chain:
        current_section_value -= e8_chain_jet / u**2
    section_jet = (
        current_section_value.numerator() // current_section_value.denominator()
    )
    section_jet_bivariate = sum(
        QQ(value) * uu**degree
        for degree, value in enumerate(section_jet.list())
    )
    transformed_bivariate = bivariate_numerator(
        uu, qq + section_jet_bivariate
    )
    transformed_cover = old_base_ring.zero()
    for monomial, value in transformed_bivariate.dict().items():
        transformed_cover += (
            parameter_field(value * q**monomial[1]) * u_parameter**monomial[0]
        )
    transformed_factorization = transformed_cover.factor()
    transformed_odd_degree = sum(
        factor.degree()
        for factor, exponent in transformed_factorization
        if exponent % 2
    )
    print(
        f"H92Q6CHORD|current_section_jet={section_jet}|factors="
        + ",".join(
            f"{factor.degree()}^{exponent}"
            for factor, exponent in transformed_factorization
        )
        + f"|odd_degree={transformed_odd_degree}",
        flush=True,
    )
    raise SystemExit(0)
if args.test_section_jet:
    coefficient_ring = PolynomialRing(QQ, "u")
    u_coefficient = coefficient_ring.gen()
    q_coefficients_qq = []
    for q_degree in range(5):
        coefficient = bivariate_numerator.coefficient({qq: q_degree})
        univariate = coefficient_ring(
            sum(
                QQ(value) * u_coefficient**monomial[0]
                for monomial, value in coefficient.dict().items()
            )
        )
        assert univariate % u_coefficient ** (2 * q_degree) == 0
        q_coefficients_qq.append(
            univariate // u_coefficient ** (2 * q_degree)
        )
    tangent_slope = (3 * xP**2 + a) / (2 * yP)
    q_on_P = base_field(h_polynomial) * tangent_slope - base_field(principal_part)
    Q_on_P = u**2 * q_on_P
    section_jet = Q_on_P.numerator() // Q_on_P.denominator()
    assert (Q_on_P - section_jet).numerator().degree() < Q_on_P.denominator().degree()
    transformed_Q = q + old_base_ring(section_jet)
    transformed_cover = sum(
        old_base_ring(coefficient) * transformed_Q**degree
        for degree, coefficient in enumerate(q_coefficients_qq)
    )
    transformed_numerator = old_base_ring(transformed_cover.numerator())
    transformed_factorization = transformed_numerator.factor()
    transformed_odd_degree = sum(
        factor.degree()
        for factor, exponent in transformed_factorization
        if exponent % 2
    )
    print(
        f"H92Q6CHORD|section_jet={section_jet}|factors="
        + ",".join(
            f"{factor.degree()}^{exponent}"
            for factor, exponent in transformed_factorization
        )
        + f"|odd_degree={transformed_odd_degree}",
        flush=True,
    )
    raise SystemExit(0)
if args.test_e7_depressed:
    coefficient_ring = PolynomialRing(QQ, "u")
    u_coefficient = coefficient_ring.gen()
    q_coefficients_qq = []
    for q_degree in range(5):
        coefficient = bivariate_numerator.coefficient({qq: q_degree})
        univariate = coefficient_ring(
            sum(
                QQ(value) * u_coefficient**monomial[0]
                for monomial, value in coefficient.dict().items()
            )
        )
        assert univariate % u_coefficient ** (2 * q_degree) == 0
        q_coefficients_qq.append(
            univariate // u_coefficient ** (2 * q_degree)
        )
    assert q_coefficients_qq[4].degree() == 0
    e7_shift = q_coefficients_qq[3] / (4 * q_coefficients_qq[4])
    assert e7_shift.denominator() in QQ and e7_shift.numerator().degree() == 5
    transformed_Q = q - old_base_ring(e7_shift)
    transformed_cover = sum(
        old_base_ring(coefficient) * transformed_Q**degree
        for degree, coefficient in enumerate(q_coefficients_qq)
    )
    transformed_numerator = old_base_ring(transformed_cover.numerator())
    transformed_factorization = transformed_numerator.factor()
    transformed_odd_degree = sum(
        factor.degree()
        for factor, exponent in transformed_factorization
        if exponent % 2
    )
    print(
        f"H92Q6CHORD|e7_shift={e7_shift}|factors="
        + ",".join(
            f"{factor.degree()}^{exponent}"
            for factor, exponent in transformed_factorization
        )
        + f"|odd_degree={transformed_odd_degree}",
        flush=True,
    )
    raise SystemExit(0)
if args.scan_simple_hensel:
    simple_factors = [
        factor
        for factor, exponent in edge_polynomial_qq.factor()
        if factor.degree() == 1 and exponent == 1
    ]
    assert len(simple_factors) == 1
    simple_root = -simple_factors[0][0] / simple_factors[0][1]
    coefficient_ring = PolynomialRing(QQ, "u")
    u_coefficient = coefficient_ring.gen()
    q_coefficients_qq = []
    for q_degree in range(5):
        coefficient = bivariate_numerator.coefficient({qq: q_degree})
        univariate = coefficient_ring(
            sum(
                QQ(value) * u_coefficient**monomial[0]
                for monomial, value in coefficient.dict().items()
            )
        )
        assert univariate % u_coefficient ** (2 * q_degree) == 0
        q_coefficients_qq.append(
            univariate // u_coefficient ** (2 * q_degree)
        )
    series_ring = PowerSeriesRing(QQ, "v", default_prec=12)
    v_series = series_ring.gen()
    Q_series = series_ring(simple_root)
    for unused_iteration in range(5):
        value = sum(
            series_ring(coefficient(v_series)) * Q_series**degree
            for degree, coefficient in enumerate(q_coefficients_qq)
        )
        derivative = sum(
            degree * series_ring(coefficient(v_series)) * Q_series ** (degree - 1)
            for degree, coefficient in enumerate(q_coefficients_qq)
            if degree
        )
        Q_series -= value / derivative
    print(
        f"H92Q6CHORD|simple_root={simple_root}|series={Q_series}",
        flush=True,
    )
    for jet_order in range(1, 11):
        jet = old_base_ring(
            [parameter_field(Q_series[index]) for index in range(jet_order)]
        )
        transformed_Q = jet + u_parameter**jet_order * q
        transformed_cover = sum(
            old_base_ring(coefficient) * transformed_Q**degree
            for degree, coefficient in enumerate(q_coefficients_qq)
        )
        transformed_numerator = old_base_ring(transformed_cover.numerator())
        transformed_factorization = transformed_numerator.factor()
        transformed_odd_degree = sum(
            factor.degree()
            for factor, exponent in transformed_factorization
            if exponent % 2
        )
        print(
            f"H92Q6CHORD|simple_jet={jet_order}|factors="
            + ",".join(
                f"{factor.degree()}^{exponent}"
                for factor, exponent in transformed_factorization
            )
            + f"|odd_degree={transformed_odd_degree}",
            flush=True,
        )
    raise SystemExit(0)
if args.scan_infinity_scales:
    for scale_exponent in range(1, 9):
        scaled_raw_chord = (
            u_parameter**scale_exponent * q + principal_parameter
        ) / h_parameter
        scaled_cover = (
            scaled_raw_chord**4
            - 6 * xP_parameter * scaled_raw_chord**2
            + 8 * yP_parameter * scaled_raw_chord
            - 3 * xP_parameter**2
            - 4 * a_parameter
        )
        scaled_numerator = old_base_ring(scaled_cover.numerator())
        scaled_factorization = scaled_numerator.factor()
        scaled_odd_degree = sum(
            factor.degree()
            for factor, exponent in scaled_factorization
            if exponent % 2
        )
        print(
            f"H92Q6CHORD|infinity_scale={scale_exponent}|"
            f"factors="
            + ",".join(
                f"{factor.degree()}^{exponent}"
                for factor, exponent in scaled_factorization
            )
            + f"|odd_degree={scaled_odd_degree}",
            flush=True,
        )
    raise SystemExit(0)
if args.search_modp:
    prime = ZZ(args.search_modp)
    finite_field = GF(prime)
    finite_ring = PolynomialRing(finite_field, names=("u", "q"))
    u_finite, q_finite = finite_ring.gens()
    finite_cover = finite_ring(bivariate_numerator)
    finite_univariate = PolynomialRing(finite_field, "v")
    v_finite = finite_univariate.gen()
    cover_coefficients = [
        finite_univariate(
            sum(
                coefficient * v_finite**monomial[0]
                for monomial, coefficient in finite_cover.dict().items()
                if monomial[1] == degree
            )
        )
        for degree in range(5)
    ]
    edge_polynomial_ring = PolynomialRing(finite_field, "Q")
    Q_edge = edge_polynomial_ring.gen()
    edge_polynomial = sum(
        coefficient[2 * degree] * Q_edge**degree
        for degree, coefficient in enumerate(cover_coefficients)
        if coefficient.valuation() == 2 * degree
    )
    print(
        f"H92Q6RR|prime={prime}|u0_edge={edge_polynomial}|"
        f"roots={edge_polynomial.roots()}",
        flush=True,
    )
    if args.scan_triple_jets:
        def jet_odd_squareclass_degree(polynomial):
            if not polynomial:
                return -1
            return sum(
                factor.degree()
                for factor, exponent in polynomial.squarefree_decomposition()
                if exponent % 2
            )

        repeated_roots = [
            root for root, multiplicity in edge_polynomial.roots() if multiplicity == 3
        ]
        assert len(repeated_roots) == 1
        jet = finite_univariate(repeated_roots[0])
        jet_ring = PolynomialRing(finite_field, names=("v", "S"))
        v_jet, S_jet = jet_ring.gens()
        Q_coefficients = [
            coefficient // v_finite ** (2 * degree)
            for degree, coefficient in enumerate(cover_coefficients)
        ]
        for jet_order in range(1, 11):
            substituted_Q = jet_ring(jet) + v_jet**jet_order * S_jet
            transformed = sum(
                jet_ring(coefficient) * substituted_Q**degree
                for degree, coefficient in enumerate(Q_coefficients)
            )
            coefficient_valuations = []
            for degree in range(5):
                coefficient = transformed.coefficient({S_jet: degree})
                if coefficient:
                    coefficient_valuations.append(
                        min(monomial[0] for monomial in coefficient.dict())
                    )
            common_valuation = min(coefficient_valuations)
            transformed //= v_jet**common_valuation
            edge = edge_polynomial_ring(
                sum(
                    value * Q_edge**monomial[1]
                    for monomial, value in transformed.dict().items()
                    if monomial[0] == 0
                )
            )
            sample_degrees = []
            for sample in (1, 2, 3):
                sample_Q = jet + v_finite**jet_order * finite_field(sample)
                sample_polynomial = sum(
                    coefficient * sample_Q**degree
                    for degree, coefficient in enumerate(Q_coefficients)
                )
                sample_degrees.append(jet_odd_squareclass_degree(sample_polynomial))
            print(
                f"H92Q6RR|triple_jet={jet_order}|center={jet}|"
                f"common_valuation={common_valuation}|edge={edge}|"
                f"edge_factorization={edge.factor()}|"
                f"sample_odd_degrees={sample_degrees}",
                flush=True,
            )
            repeated = sorted(
                edge.roots(), key=lambda item: (item[1], int(item[0])), reverse=True
            )
            if not repeated or repeated[0][1] < 2:
                break
            jet += repeated[0][0] * v_finite**jet_order
        raise SystemExit(0)
    if args.triple_diagnostic:
        repeated_roots = [
            root for root, multiplicity in edge_polynomial.roots() if multiplicity == 3
        ]
        assert len(repeated_roots) == 1
        repeated_root = repeated_roots[0]
        translated_ring = PolynomialRing(finite_field, names=("v", "R"))
        v_translated, R_translated = translated_ring.gens()
        translated = translated_ring.zero()
        for degree, coefficient in enumerate(cover_coefficients):
            reduced_coefficient = coefficient // v_finite ** (2 * degree)
            translated += translated_ring(reduced_coefficient) * (
                R_translated + repeated_root
            ) ** degree
        translated_valuations = []
        for degree in range(5):
            coefficient = translated.coefficient({R_translated: degree})
            translated_valuations.append(
                -1
                if not coefficient
                else min(monomial[0] for monomial in coefficient.dict())
            )
        print(
            f"H92Q6RR|triple_root={repeated_root}|"
            f"translated_valuations={translated_valuations}|"
            f"translated={translated}",
            flush=True,
        )
        raise SystemExit(0)

    def rref_planes():
        for pivot0 in range(3):
            for pivot1 in range(pivot0 + 1, 4):
                free_positions = [
                    (row, column)
                    for row, pivot in enumerate((pivot0, pivot1))
                    for column in range(pivot + 1, 4)
                    if column not in (pivot0, pivot1)
                ]
                for values in itertools.product(finite_field, repeat=len(free_positions)):
                    rows = [[finite_field.zero()] * 4 for _ in range(2)]
                    rows[0][pivot0] = rows[1][pivot1] = finite_field.one()
                    for (row, column), value in zip(free_positions, values):
                        rows[row][column] = value
                    yield rows

    def odd_squareclass_degree(polynomial):
        if not polynomial:
            return -1
        return sum(
            factor.degree()
            for factor, exponent in polynomial.squarefree_decomposition()
            if exponent % 2
        )

    specialization_values = [finite_field(value) for value in (1, 2, 3)]
    tested = 0
    hits = []
    for infinity_scale in range(args.scale_min, args.scale_max + 1):
        for rows in rref_planes():
            a = v_finite**args.polynomial_scale * (
                rows[0][0] + rows[0][1] * v_finite
            )
            b_row = v_finite**infinity_scale * (
                rows[0][2] + rows[0][3] * v_finite
            )
            c = v_finite**args.polynomial_scale * (
                rows[1][0] + rows[1][1] * v_finite
            )
            d_row = v_finite**infinity_scale * (
                rows[1][2] + rows[1][3] * v_finite
            )
            determinant = a * d_row - b_row * c
            if not determinant:
                continue
            tested += 1
            degrees = []
            for parameter_value in specialization_values:
                q_numerator = c - parameter_value * a
                q_denominator = parameter_value * b_row - d_row
                if not q_denominator:
                    degrees.append(-1)
                    break
                transformed = sum(
                    cover_coefficients[degree]
                    * q_numerator**degree
                    * q_denominator ** (4 - degree)
                    for degree in range(5)
                )
                degrees.append(odd_squareclass_degree(transformed))
                if degrees[-1] not in (3, 4):
                    break
            if len(degrees) == len(specialization_values) and all(
                degree in (3, 4) for degree in degrees
            ):
                hit = {
                    "infinity_scale": infinity_scale,
                    "rows": [[int(value) for value in row] for row in rows],
                    "sample_odd_degrees": degrees,
                }
                hits.append(hit)
                print(f"H92Q6RR|hit={hit}", flush=True)
        print(
            f"H92Q6RR|prime={prime}|poly_scale={args.polynomial_scale}|"
            f"q_scale={infinity_scale}|"
            f"tested={tested}|hits={len(hits)}",
            flush=True,
        )
    print(
        f"H92Q6RR|prime={prime}|tested={tested}|hits={len(hits)}|"
        "status=PASS_BOUNDED_MODP_RR_SEARCH",
        flush=True,
    )
    raise SystemExit(0)
factorization = numerator.factor()
factor_data = [
    {
        "degree_in_u": int(factor.degree()),
        "exponent": int(exponent),
        "factor": str(factor),
    }
    for factor, exponent in factorization
]
odd_part = prod(
    factor for factor, exponent in factorization if exponent % 2
)
odd_degree = ZZ(odd_part.degree())
assert denominator.is_square()
print(
    "H92Q6CHORD|raw_factors="
    + ",".join(f"{factor.degree()}^{exponent}" for factor, exponent in factorization),
    flush=True,
)
assert odd_degree in (3, 4)

payload = {
    "schema": "elkies-k3.h92-q6-chord.v1",
    "status": "PASS_EXACT_H92_Q6_CHORD",
    "input": {
        "section": str(args.section.relative_to(ROOT)),
        "sha256": SECTION_SHA256,
    },
    "section": {
        "label": "P1",
        "orientation": "ancillary split-component marking",
        "x": str(xP),
        "y": str(yP),
    },
    "divisor": {
        "class": "D=O+(-P1)-F",
        "chord_section_pole": "-P1",
        "raw_chord": "m=(y-y(P1))/(x-x(P1))",
        "full_affine_denominator_root": str(z_polynomial),
        "smooth_PO_saturation_polynomial_h": str(h_polynomial),
        "principal_part_A": str(principal_part),
        "basis": ["1", "q=z*m-A"],
        "dimension": 2,
    },
    "double_cover": {
        "identity": (
            "w^2=m^4-6*xP*m^2+8*yP*m-3*xP^2-4*a, "
            "m=(q+A)/z"
        ),
        "numerator_degree_in_u": int(numerator.degree()),
        "denominator_degree_in_u": int(denominator.degree()),
        "denominator_is_square": True,
        "factorization": factor_data,
        "odd_part": str(odd_part),
        "odd_part_degree_in_u": int(odd_degree),
        "generic_fiber_genus": 1,
    },
    "proof_boundary": (
        "This certifies the exact marked section, chord map, and genus-one "
        "pencil over QQ.  A minimal Weierstrass model and Kodaira-fiber "
        "classification of the q-pencil are separate downstream steps."
    ),
}

args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    f"H92Q6CHORD|numerator_degree={numerator.degree()}|"
    f"odd_degree={odd_degree}|genus=1|basis_dimension=2",
    flush=True,
)
print("H92Q6CHORD|status=PASS_EXACT_H92_Q6_CHORD", flush=True)
