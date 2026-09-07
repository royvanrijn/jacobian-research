#!/usr/bin/env sage -python
"""Lift odd-I4 q4/orbit164 polynomial-section candidates to QQ.

The modular input fixes one reducible-fibre node.  Add its exact x/y
incidence equations to the thirteen Weierstrass coefficient equations, choose
twelve independent rows, and Newton lift each regular branch.  Rational
reconstruction followed by literal substitution retains only QQ sections.
No Groebner basis is used.
"""

import argparse
import hashlib
import json
import time
from pathlib import Path

from sage.all import GF, NumberField, PolynomialRing, QQ, Qp, ZZ, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
MODEL = LOCAL / "q4o164-compact-weierstrass-qq.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--modular", type=Path,
    default=LOCAL / "q4o164-odd-finite-i4-sections-mod41.json",
)
parser.add_argument(
    "--output", type=Path,
    default=LOCAL / "q4o164-odd-finite-i4-sections-qq.json",
)
parser.add_argument("--precision", type=int, default=180)
parser.add_argument(
    "--field-orbit",
    default="",
    help=(
        "comma-separated modular indices for conjugate representatives; "
        "recover their common number field by symmetric functions"
    ),
)
parser.add_argument(
    "--quadratic-pair-scan", action="store_true",
    help="score all non-QQ branch pairs by rational trace and norm recognition",
)
parser.add_argument(
    "--rational-sum-scan", action="store_true",
    help=(
        "scan signed sums of the non-QQ inverse-pair representatives and "
        "retain sums that reconstruct and verify over QQ(t)"
    ),
)
args = parser.parse_args()
MODULAR = args.modular if args.modular.is_absolute() else ROOT / args.modular
OUTPUT = args.output if args.output.is_absolute() else ROOT / args.output
PRECISION = int(args.precision)
FIELD_ORBIT = [int(value) for value in args.field_orbit.split(",") if value]
if FIELD_ORBIT and len(FIELD_ORBIT) < 2:
    parser.error("--field-orbit requires at least two modular indices")

started = time.monotonic()
model = json.loads(MODEL.read_text())
modular = json.loads(MODULAR.read_text())
assert model["status"] == "PASS_EXACT_QQ_Q4O164_COMPACT_WEIERSTRASS_NORMALIZATION"
assert modular["status"] == "PASS_MODP_Q4O164_ODD_I4_SECTION_SCAN"
PRIME = ZZ(modular["prime"])
selected = int(modular["selected_fibre"]["index"])

RQ = PolynomialRing(QQ, "t")
tq = RQ.gen()
A_QQ = RQ([QQ(value) for value in model["compact_model"]["A_coefficients_low_to_high"]])
B_QQ = RQ([QQ(value) for value in model["compact_model"]["B_coefficients_low_to_high"]])


def reduce_qq(value, field):
    value = QQ(value)
    return field(value.numerator()) / field(value.denominator())


supports_QQ = [
    None if row["support"] == "infinity" else QQ(row["support"])
    for row in model["compact_model"]["reducible_fibres"]
]
support_QQ = supports_QQ[selected]
assert support_QQ is not None
RX = PolynomialRing(QQ, "x")
xvar = RX.gen()
cubic = xvar**3 + A_QQ(support_QQ) * xvar + B_QQ(support_QQ)
repeated = cubic.gcd(cubic.derivative())
assert repeated.degree() == 1
node_QQ = QQ(-repeated[0] / repeated[1])

F = GF(PRIME)
RF = PolynomialRing(F, "t")
tf = RF.gen()
A_F = RF([reduce_qq(value, F) for value in A_QQ])
B_F = RF([reduce_qq(value, F) for value in B_QQ])
support_F = reduce_qq(support_QQ, F)
node_F = reduce_qq(node_QQ, F)
assert int(node_F) == int(modular["selected_fibre"]["node"])

K = Qp(PRIME, prec=PRECISION, type="capped-rel")
RT = PolynomialRing(K, "t")
t = RT.gen()
A = RT([K(value) for value in A_QQ])
B = RT([K(value) for value in B_QQ])
support = K(support_QQ)
node = K(node_QQ)


def polynomials(values, ring):
    return ring(list(values[:5])), ring(list(values[5:]))


def residual(values):
    X, Y = polynomials(values, RT)
    equation = Y**2 - X**3 - A * X - B
    return vector(K, [equation[index] for index in range(13)] + [
        X(support) - node, Y(support),
    ])


def jacobian(values, ring, surface_A, fibre_support):
    X, Y = polynomials(values, ring)
    dx = -3 * X**2 - surface_A
    dy = 2 * Y
    zero = ring.base_ring().zero()
    rows = [[
        dx[degree-shift] if 0 <= degree-shift <= dx.degree() else zero
        for shift in range(5)
    ] + [
        dy[degree-shift] if 0 <= degree-shift <= dy.degree() else zero
        for shift in range(7)
    ] for degree in range(13)]
    rows.append([fibre_support**degree for degree in range(5)] + [zero] * 7)
    rows.append([zero] * 5 + [fibre_support**degree for degree in range(7)])
    return matrix(ring.base_ring(), rows)


def minimum_valuation(values):
    nonzero = [value.valuation() for value in values if value]
    return min(nonzero) if nonzero else PRECISION


def rational_reconstruct(value):
    if not value:
        return QQ.zero()
    valuation_shift = max(0, -int(value.valuation())) if value else 0
    integral_value = value * PRIME**valuation_shift
    usable_precision = min(
        PRECISION - 12,
        int(integral_value.precision_absolute()) - 4,
    )
    if usable_precision < 16:
        raise ArithmeticError("insufficient p-adic precision for rational reconstruction")
    modulus = PRIME ** usable_precision
    residue = ZZ(integral_value.lift()) % modulus
    return QQ(residue.rational_reconstruction(modulus)) / PRIME**valuation_shift


def coefficient_bits(value):
    value = QQ(value)
    return max(abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())


lifted = []
nonrational_dependencies = []
padic_branches = {}
for modular_index, row in enumerate(modular["sections"]):
    seed = vector(F, row["x_coefficients_low_to_high"] + row["y_coefficients_low_to_high"])
    seed = vector(F, list(seed[:5]) + list(seed[5:]) + [F.zero()] * (12-len(seed)))
    XF, YF = polynomials(seed, RF)
    assert YF**2 == XF**3 + A_F * XF + B_F
    assert XF(support_F) == node_F and YF(support_F) == 0
    JF = jacobian(seed, RF, A_F, support_F)
    assert JF.rank() == 12
    pivot_rows = list(map(int, JF.transpose().pivots()))
    determinant = int(matrix(F, [JF.row(index) for index in pivot_rows]).det())
    assert determinant

    values = vector(K, [K(value).add_bigoh(1) for value in seed])
    known_precision = 1
    iterations = []
    while known_precision < PRECISION:
        working_precision = min(2 * known_precision, PRECISION)
        values = vector(K, [K(value.lift()).add_bigoh(working_precision) for value in values])
        full = residual(values)
        chosen = vector(K, [full[index] for index in pivot_rows])
        square = matrix(K, [jacobian(values, RT, A, support).row(index) for index in pivot_rows])
        correction = square.solve_right(-chosen)
        values += correction
        iterations.append({
            "working_precision_p_adic_digits": working_precision,
            "minimum_full_residual_valuation_after": int(minimum_valuation(residual(values))),
            "minimum_correction_valuation": int(minimum_valuation(correction)),
        })
        known_precision = working_precision

    try:
        reconstructed = [rational_reconstruct(value) for value in values]
    except ArithmeticError:
        padic_branches[modular_index] = values
        nonrational_dependencies.append({
            "modular_candidate_index": modular_index,
            "dependencies_degree_at_most_four": [
                str(value.algebraic_dependency(4)) for value in values
            ],
        })
        continue
    X, Y = polynomials(reconstructed, RQ)
    if Y**2 != X**3 + A_QQ * X + B_QQ:
        continue
    if X(support_QQ) != node_QQ or Y(support_QQ) != 0:
        continue
    lifted.append({
        "modular_candidate_index": modular_index,
        "selected_independent_equation_rows": pivot_rows,
        "selected_jacobian_determinant_modp": determinant,
        "iterations": iterations,
        "x_coefficients_low_to_high": [str(value) for value in X.list()],
        "y_coefficients_low_to_high": [str(value) for value in Y.list()],
        "maximum_x_rational_bits": max(map(coefficient_bits, X)),
        "maximum_y_rational_bits": max(map(coefficient_bits, Y)),
        "exact_compact_weierstrass_identity": True,
        "exact_selected_node_incidence": True,
    })
    print(
        f"Q4O164ODDLIFT|modular_index={modular_index}|bits="
        f"{lifted[-1]['maximum_x_rational_bits']}/{lifted[-1]['maximum_y_rational_bits']}|"
        f"elapsed={time.monotonic()-started:.3f}",
        flush=True,
    )

number_field_lifts = []
quadratic_pair_scores = []
rational_section_sums = []
if args.quadratic_pair_scan:
    branch_indices = sorted(padic_branches)
    for left_position, left_index in enumerate(branch_indices):
        for right_index in branch_indices[left_position+1:]:
            recognized_coordinates = 0
            for left, right in zip(padic_branches[left_index], padic_branches[right_index]):
                try:
                    rational_reconstruct(left + right)
                    rational_reconstruct(left * right)
                    recognized_coordinates += 1
                except ArithmeticError:
                    pass
            quadratic_pair_scores.append({
                "modular_indices": [left_index, right_index],
                "rational_trace_and_norm_coordinates": recognized_coordinates,
                "total_coordinates": 12,
            })
    quadratic_pair_scores.sort(
        key=lambda row: (-row["rational_trace_and_norm_coordinates"], row["modular_indices"])
    )
    print(
        "Q4O164ODDLIFT|quadratic_pair_scores="
        + ";".join(
            f"{row['modular_indices'][0]},{row['modular_indices'][1]}:"
            f"{row['rational_trace_and_norm_coordinates']}"
            for row in quadratic_pair_scores[:12]
        ),
        flush=True,
    )

if args.rational_sum_scan:
    fraction_field = RT.fraction_field()

    def add_points(left, right, curve_A):
        if left is None:
            return right
        if right is None:
            return left
        x_left, y_left = left
        x_right, y_right = right
        if x_left == x_right:
            if y_left == -y_right:
                return None
            slope = (3*x_left**2 + curve_A) / (2*y_left)
        else:
            slope = (y_right-y_left) / (x_right-x_left)
        x_sum = slope**2 - x_left - x_right
        y_sum = -y_left + slope*(x_left-x_sum)
        return x_sum, y_sum

    def reconstruct_fraction(value, exact_ring):
        numerator = value.numerator()
        denominator = value.denominator()
        scale = denominator.leading_coefficient()
        numerator /= scale
        denominator /= scale
        numerator_QQ = exact_ring([rational_reconstruct(entry) for entry in numerator])
        denominator_QQ = exact_ring([rational_reconstruct(entry) for entry in denominator])
        return exact_ring.fraction_field()(numerator_QQ / denominator_QQ)

    # One representative from each inverse pair.  A rational odd signed sum
    # can expose the arithmetic ninth direction even when no individual
    # polynomial representative is defined over QQ.
    representatives = []
    for index in sorted(padic_branches):
        row = modular["sections"][index]
        inverse_seen = any(
            row["x_coefficients_low_to_high"] == modular["sections"][old]["x_coefficients_low_to_high"]
            for old in representatives
        )
        if not inverse_seen:
            representatives.append(index)
    assert len(representatives) == 6
    points = {}
    branch_polynomials = {}
    for index in representatives:
        X, Y = polynomials(padic_branches[index], RT)
        branch_polynomials[index] = X, Y
        points[index] = fraction_field(X), fraction_field(Y)

    screen_supports = [K(value) for value in (0, 1, 2)]
    screen_points = {
        support_value: {
            index: (X(support_value), Y(support_value))
            for index, (X, Y) in points.items()
        }
        for support_value in screen_supports
    }

    def signed_sum(coefficients, point_table, curve_A):
        point = None
        for coefficient, index in zip(coefficients, representatives):
            if not coefficient:
                continue
            X, Y = point_table[index]
            point = add_points(point, (X, Y if coefficient > 0 else -Y), curve_A)
        return point

    from itertools import product
    survivor_words = []
    for coefficients in product((-1, 0, 1), repeat=len(representatives)):
        if not any(coefficients):
            continue
        first_nonzero = next(value for value in coefficients if value)
        if first_nonzero < 0:
            continue
        screen_passes = True
        for support_value in screen_supports:
            screen_point = signed_sum(
                coefficients, screen_points[support_value], A(support_value)
            )
            if screen_point is None:
                screen_passes = False
                break
            try:
                rational_reconstruct(screen_point[0])
                rational_reconstruct(screen_point[1])
            except ArithmeticError:
                screen_passes = False
                break
        if not screen_passes:
            continue
        survivor_words.append(coefficients)

    print(
        f"Q4O164ODDLIFT|signed_words=364|"
        f"scalar_screen_survivors={len(survivor_words)}|stage=SCALAR_SCREEN",
        flush=True,
    )

    def interpolate_rational(samples, numerator_degree, denominator_degree):
        unknowns = numerator_degree + 1 + denominator_degree
        if len(samples) < unknowns:
            raise ArithmeticError("not enough regular specializations for interpolation")
        rows = []
        right = []
        for support_value, value in samples[:unknowns]:
            rows.append(
                [support_value**degree for degree in range(numerator_degree+1)]
                + [-value*support_value**degree for degree in range(denominator_degree)]
            )
            right.append(value*support_value**denominator_degree)
        solution = matrix(K, rows).solve_right(vector(K, right))
        solution_entries = list(solution)
        numerator = RT(solution_entries[:numerator_degree+1])
        denominator = RT(solution_entries[numerator_degree+1:] + [K.one()])
        for support_value, value in samples[unknowns:]:
            if minimum_valuation([numerator(support_value)-value*denominator(support_value)]) < PRECISION-20:
                raise ArithmeticError("rational interpolation fails held-out specialization")
        numerator_QQ = RQ([rational_reconstruct(value) for value in numerator])
        denominator_QQ = RQ([rational_reconstruct(value) for value in denominator])
        return RQ.fraction_field()(numerator_QQ/denominator_QQ)

    for coefficients in survivor_words:
        samples_x = []
        samples_y = []
        for residue in range(int(PRIME)):
            support_value = K(residue)
            table = {
                index: (X(support_value), Y(support_value))
                for index, (X, Y) in branch_polynomials.items()
            }
            try:
                point = signed_sum(coefficients, table, A(support_value))
            except (ArithmeticError, ZeroDivisionError):
                continue
            if point is None:
                continue
            samples_x.append((support_value, point[0]))
            samples_y.append((support_value, point[1]))
        try:
            # The marked q8 horizontal has normalized degrees X/Y/Z=12/18/4;
            # Y degree 17 is allowed when the leading coefficient cancels.
            X_QQ = interpolate_rational(samples_x, 12, 8)
            Y_QQ = interpolate_rational(samples_y, 18, 12)
        except (ArithmeticError, ZeroDivisionError, ValueError, TypeError):
            continue
        if Y_QQ**2 != X_QQ**3 + RQ.fraction_field()(A_QQ)*X_QQ + B_QQ:
            continue
        record = {
            "representative_indices": representatives,
            "signed_coefficients": list(coefficients),
            "x_numerator_coefficients_low_to_high": [str(value) for value in X_QQ.numerator()],
            "x_denominator_coefficients_low_to_high": [str(value) for value in X_QQ.denominator()],
            "y_numerator_coefficients_low_to_high": [str(value) for value in Y_QQ.numerator()],
            "y_denominator_coefficients_low_to_high": [str(value) for value in Y_QQ.denominator()],
            "degrees_x_numerator_y_numerator_x_denominator_y_denominator": [
                X_QQ.numerator().degree(), Y_QQ.numerator().degree(),
                X_QQ.denominator().degree(), Y_QQ.denominator().degree(),
            ],
            "exact_compact_weierstrass_identity": True,
        }
        rational_section_sums.append(record)
        print(
            f"Q4O164ODDLIFT|rational_sum={coefficients}|"
            f"degrees={record['degrees_x_numerator_y_numerator_x_denominator_y_denominator']}|"
            "status=PASS_EXACT_QQ_SIGNED_SUM",
            flush=True,
        )

    print(
        f"Q4O164ODDLIFT|signed_words=364|"
        f"scalar_screen_survivors={len(survivor_words)}|"
        f"exact_rational_sums={len(rational_section_sums)}|"
        "status=PASS_BOUNDED_RATIONAL_SUM_SCAN",
        flush=True,
    )

if FIELD_ORBIT:
    assert all(index in padic_branches for index in FIELD_ORBIT)
    orbit_values = [padic_branches[index] for index in FIELD_ORBIT]
    field_degree = len(FIELD_ORBIT)

    # Choose a section coefficient whose three residues are distinct.  It is
    # then a primitive element for this split p-adic cubic orbit.
    primitive_coordinate = next(
        coordinate for coordinate in range(12)
        if len({int(modular["sections"][index][
            "x_coefficients_low_to_high" if coordinate < 5
            else "y_coefficients_low_to_high"
        ][coordinate if coordinate < 5 else coordinate-5]) for index in FIELD_ORBIT}) == field_degree
    )
    alphas = [values[primitive_coordinate] for values in orbit_values]
    SP = PolynomialRing(K, "z")
    z = SP.gen()
    minimal_padic = SP.one()
    for alpha in alphas:
        minimal_padic *= z - alpha
    minimal_coefficients = [rational_reconstruct(value) for value in minimal_padic]
    SQ = PolynomialRing(QQ, "z")
    minimal = SQ(minimal_coefficients)
    assert minimal.degree() == field_degree and minimal.is_monic() and minimal.is_irreducible()
    assert all(minimum_valuation(minimal(alpha)) >= PRECISION-12 for alpha in alphas)

    vandermonde = matrix(K, [
        [alpha**power for power in range(field_degree)] for alpha in alphas
    ])
    coordinate_expressions = []
    for coordinate in range(12):
        interpolated = vandermonde.solve_right(
            vector(K, [values[coordinate] for values in orbit_values])
        )
        exact = [rational_reconstruct(value) for value in interpolated]
        assert all(
            minimum_valuation([sum(K(exact[power]) * alpha**power for power in range(field_degree))
                               - values[coordinate]]) >= PRECISION-12
            for alpha, values in zip(alphas, orbit_values)
        )
        coordinate_expressions.append(exact)

    NF = NumberField(minimal, "a")
    a = NF.gen()
    exact_values = [
        sum(NF(expression[power]) * a**power for power in range(field_degree))
        for expression in coordinate_expressions
    ]
    RNF = PolynomialRing(NF, "t")
    tn = RNF.gen()
    XNF, YNF = polynomials(exact_values, RNF)
    ANF = RNF(A_QQ)
    BNF = RNF(B_QQ)
    assert YNF**2 == XNF**3 + ANF * XNF + BNF
    assert XNF(NF(support_QQ)) == NF(node_QQ) and YNF(NF(support_QQ)) == 0

    # The chosen prime splits completely.  Check that a -> alpha mod p for
    # each root reproduces every modular branch used in the interpolation.
    reduction_checks = []
    minimal_F = PolynomialRing(F, "z")(minimal)
    for orbit_position, modular_index in enumerate(FIELD_ORBIT):
        alpha_F = F(int(alphas[orbit_position].lift()))
        assert minimal_F(alpha_F) == 0
        reduced = [
            sum(reduce_qq(expression[power], F) * alpha_F**power for power in range(field_degree))
            for expression in coordinate_expressions
        ]
        expected = modular["sections"][modular_index]
        assert [int(value) for value in reduced[:5]] == expected["x_coefficients_low_to_high"]
        assert [int(value) for value in reduced[5:]] == expected["y_coefficients_low_to_high"]
        reduction_checks.append(modular_index)

    number_field_lifts.append({
        "modular_representative_indices": FIELD_ORBIT,
        "primitive_coordinate_index": primitive_coordinate,
        "minimal_polynomial_low_to_high": [str(value) for value in minimal],
        "coordinate_expressions_in_power_basis_low_to_high": [
            [str(value) for value in expression] for expression in coordinate_expressions
        ],
        "x_coefficients_low_to_high": [str(value) for value in XNF.list()],
        "y_coefficients_low_to_high": [str(value) for value in YNF.list()],
        "exact_compact_weierstrass_identity": True,
        "exact_selected_node_incidence": True,
        "modular_reduction_checks": reduction_checks,
    })
    print(
        f"Q4O164ODDLIFT|field_orbit={','.join(map(str, FIELD_ORBIT))}|"
        f"primitive_coordinate={primitive_coordinate}|minimal={minimal}|"
        "status=PASS_EXACT_NUMBER_FIELD_ORBIT_LIFT",
        flush=True,
    )

payload = {
    "schema": "elkies-k3.q4o164-odd-i4-sections-qq.v3",
    "status": "PASS_EXACT_QQ_Q4O164_ODD_I4_SECTION_LIFTS",
    "prime": int(PRIME),
    "precision_p_adic_digits": PRECISION,
    "selected_fibre": {
        "index": selected,
        "label": modular["selected_fibre"]["label"],
        "support": str(support_QQ),
        "node": str(node_QQ),
    },
    "exact_rational_lifts": lifted,
    "exact_number_field_orbit_lifts": number_field_lifts,
    "quadratic_pair_scores": quadratic_pair_scores,
    "exact_rational_signed_sums": rational_section_sums,
    "nonrational_padic_branches": nonrational_dependencies,
    "method": {
        "large_Groebner_required": False,
        "rational_sum_scan_requested": bool(args.rational_sum_scan),
        "rational_sum_scan_word_count_up_to_sign": 364 if args.rational_sum_scan else 0,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "Each retained branch is an exact QQ polynomial section with the selected node "
        "incidence. Each retained number-field orbit is an exact section over the displayed "
        "number field, reconstructed simultaneously from the specified p-adic conjugates. "
        "Each retained signed sum is verified exactly over QQ(t). Failure to retain a sum "
        "is only a bounded degree/precision/height construction result, not non-existence. "
        "Independence and marked-lattice identification remain separate gates."
    ),
    "inputs": {
        "paths": [str(path.relative_to(ROOT)) for path in (MODEL, MODULAR)],
        "sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (MODEL, MODULAR)
        },
    },
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "Q4O164ODDLIFT|modular_candidates={}|exact_lifts={}|status={}|output={}".format(
        len(modular["sections"]), len(lifted), payload["status"], OUTPUT,
    ),
    flush=True,
)
