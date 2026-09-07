#!/usr/bin/env sage
"""Compile the fixed-corridor 4A1 fibration on the exact 3A1 equation.

For the exact q114 section use

    D = O + P + 3*C1 + 2*C2 - 3F,   P.O=21.

The chord ambient has dimension 57.  Smooth collision gives dimension 15.
The two ordinary resolved-node blocks and the marked-node block have measured
individual ranks 5, 3, 5 and stacked rank 13, leaving h0=2.  Only exact
truncated local algebra and linear algebra are used.
"""

import hashlib
import json
import sys
import time
from pathlib import Path

from sage.all import GF, PolynomialRing, PowerSeriesRing, QQ, ZZ, matrix


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
TARGET = LOCAL / "fixed-reverse-4a1-horizontal-from-3a1-qq.json"
SURFACE = LOCAL / "fixed-reverse-3a1-rr-qq.json"
POINTING = LOCAL / "fixed-reverse-3a1-pointing-qq.json"
OUTPUT = LOCAL / "fixed-reverse-4a1-rr-qq.json"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficients(poly):
    return [str(value) for value in poly.list()]


def rational_bits(values):
    answer = 0
    for value in values:
        value = QQ(value)
        answer = max(answer, abs(ZZ(value.numerator())).nbits(), ZZ(value.denominator()).nbits())
    return int(answer)


started = time.monotonic()
target_artifact = json.loads(TARGET.read_text())
surface_artifact = json.loads(SURFACE.read_text())
pointing_artifact = json.loads(POINTING.read_text())
assert target_artifact["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_HORIZONTAL_ON_3A1"
assert surface_artifact["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_3A1_RR_JACOBIAN"
assert pointing_artifact["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_3A1_POINTING"
assert target_artifact["fixed_edge"]["P_dot_O"] == 21

R = PolynomialRing(QQ, "u")
u = R.gen()
K = R.fraction_field()
A = R(surface_artifact["child"]["minimal_A_coefficients_low_to_high"])
B = R(surface_artifact["child"]["minimal_B_coefficients_low_to_high"])
section = target_artifact["section"]
X = R(section["x_numerator_coefficients_low_to_high"])
Y = R(section["y_numerator_coefficients_low_to_high"])
Z = R(section["Z_coefficients_low_to_high"])
assert (X.degree(), Y.degree(), Z.degree()) == (46, 69, 21)
assert Y ** 2 == X ** 3 + A * X * Z ** 4 + B * Z ** 6
# Remove the irrelevant common projective scale before polynomial Euclid.
# This keeps the same affine point X/Z^2,Y/Z^3 and substantially reduces the
# coefficient content seen by inverse_mod.
projective_scale = Z.leading_coefficient()
X = R(X / projective_scale ** 2)
Y = R(Y / projective_scale ** 3)
Z = R(Z / projective_scale)
assert Z.is_monic()
assert Y ** 2 == X ** 3 + A * X * Z ** 4 + B * Z ** 6

# For D=O+P+kF, the compact chord bounds are
# deg(AA)<=2(P.O)+k and deg(BB)<=P.O-2+k.  Here k=-3.
aa_degree = 39
bb_degree = 16
ambient = [(u ** degree, R.zero()) for degree in range(aa_degree + 1)] + [
    (R.zero(), u ** degree) for degree in range(bb_degree + 1)
]
collision_modulus = Z ** 2
assert len(ambient) == 57 and collision_modulus.degree() == 42
SMOOTH_CHECKPOINT = LOCAL / "fixed-reverse-4a1-smooth-kernel-qq.json"
checkpoint_data = json.loads(SMOOTH_CHECKPOINT.read_text()) if SMOOTH_CHECKPOINT.exists() else None
checkpoint_valid = bool(
    checkpoint_data
    and checkpoint_data.get("status") == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_SMOOTH_KERNEL"
    and checkpoint_data.get("target_sha256") == sha256(TARGET)
    and checkpoint_data.get("surface_sha256") == sha256(SURFACE)
)
if checkpoint_valid:
    smooth_kernel = matrix(QQ, checkpoint_data["rows"])
    pivot_pair = tuple(map(int, checkpoint_data["endpoint_pivot_pair"]))
    assert smooth_kernel.dimensions() == (15, 57)
else:
    # Solve the congruence without forming a dense 42-by-57 matrix.  X is a unit
    # modulo Z^2, so BB determines AA uniquely.  Only degrees 40 and 41 of the
    # remainder must vanish; all lower coefficients are the resulting AA.
    X_inverse = X.inverse_mod(collision_modulus)
    bb_monomials = [u ** degree for degree in range(bb_degree + 1)]
    aa_remainders = [R((BB * Y * X_inverse) % collision_modulus) for BB in bb_monomials]
    endpoint_conditions = matrix(QQ, [
        [remainder[degree] for remainder in aa_remainders]
        for degree in (40, 41)
    ])
    assert endpoint_conditions.dimensions() == (2, 17)
    # A generic dense rational-kernel call clears enormous common denominators and
    # invokes IML.  Solve the two equations directly by choosing a nonzero 2-by-2
    # pivot minor and assigning each of the other 15 coordinates in turn.
    pivot_pair = next(
        (left, right)
        for left in range(17) for right in range(left + 1, 17)
        if endpoint_conditions[0, left] * endpoint_conditions[1, right]
           - endpoint_conditions[0, right] * endpoint_conditions[1, left]
    )
    pivot_left, pivot_right = pivot_pair
    pivot_det = (
        endpoint_conditions[0, pivot_left] * endpoint_conditions[1, pivot_right]
        - endpoint_conditions[0, pivot_right] * endpoint_conditions[1, pivot_left]
    )
    bb_kernel_rows = []
    for free in range(17):
        if free in pivot_pair:
            continue
        row = [QQ.zero()] * 17
        row[free] = QQ.one()
        row[pivot_left] = (
            -endpoint_conditions[0, free] * endpoint_conditions[1, pivot_right]
            + endpoint_conditions[0, pivot_right] * endpoint_conditions[1, free]
        ) / pivot_det
        row[pivot_right] = (
            -endpoint_conditions[0, pivot_left] * endpoint_conditions[1, free]
            + endpoint_conditions[0, free] * endpoint_conditions[1, pivot_left]
        ) / pivot_det
        assert endpoint_conditions * matrix(QQ, 17, 1, row) == 0
        bb_kernel_rows.append(row)
    bb_kernel = matrix(QQ, bb_kernel_rows)
    assert bb_kernel.dimensions() == (15, 17)
    smooth_rows = []
    for row in bb_kernel.rows():
        BB = sum((row[index] * bb_monomials[index] for index in range(17)), R.zero())
        AA = R((BB * Y * X_inverse) % collision_modulus)
        assert AA.degree() <= aa_degree
        assert (AA * X - BB * Y) % collision_modulus == 0
        smooth_rows.append(
            [AA[degree] for degree in range(40)]
            + [BB[degree] for degree in range(17)]
        )
    smooth_kernel = matrix(QQ, smooth_rows)
    assert smooth_kernel.dimensions() == (15, 57)
    SMOOTH_CHECKPOINT.write_text(json.dumps({
        "schema": "elkies-k3.fixed-reverse-4a1-smooth-kernel-qq.v1",
        "status": "PASS_EXACT_QQ_FIXED_REVERSE_4A1_SMOOTH_KERNEL",
        "target_sha256": sha256(TARGET),
        "surface_sha256": sha256(SURFACE),
        "projective_scale": str(projective_scale),
        "endpoint_pivot_pair": list(map(int, pivot_pair)),
        "rows": [[str(value) for value in row] for row in smooth_kernel.rows()],
    }, indent=2, sort_keys=True) + "\n")

# -------------------------------------------------------------------------
# Complete resolved A1-node quotients.  Modulo m^k (k<=3), the completed
# split-node relation is ab=kappa*tau^2; all omitted terms have total degree
# at least k.  The quotient bases therefore have intrinsic lengths k^2.
# -------------------------------------------------------------------------
Hx = K(X) / K(Z ** 2)
Hy = K(Y) / K(Z ** 3)


def series_at(value, support, precision):
    PS = PowerSeriesRing(QQ, "tau", default_prec=precision)
    tau = PS.gen()
    value = K(value)

    def evaluate(poly):
        answer = PS.zero()
        for coefficient in reversed(R(poly).list()):
            answer = answer * (support + tau) + coefficient
        return answer

    return evaluate(value.numerator()) / evaluate(value.denominator())


def reduce_monomial(key, coefficient, node_lambda, kappa, order, cubic_corrections=()):
    """Reduce xi^i*y^j*tau^l by y^2=lambda*xi^2+kappa*tau^2."""
    xi_power, y_power, tau_power = key
    terms = [(xi_power, y_power, tau_power, coefficient)]
    while any(term[1] >= 2 for term in terms):
        reduced = []
        for i, j, ell, value in terms:
            if j < 2:
                reduced.append((i, j, ell, value))
                continue
            reduced.append((i + 2, j - 2, ell, value * node_lambda))
            reduced.append((i, j - 2, ell + 2, value * kappa))
            for xi_increment, tau_increment, correction in cubic_corrections:
                reduced.append((
                    i + xi_increment, j - 2, ell + tau_increment,
                    value * correction,
                ))
        terms = reduced
    return [
        ((i, j, ell), value)
        for i, j, ell, value in terms
        if i + j + ell < order and value
    ]


def alg_add(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = answer.get(key, QQ.zero()) + value
        if not answer[key]:
            del answer[key]
    return answer


def alg_scale(value, scalar):
    return {key: scalar * coefficient for key, coefficient in value.items() if scalar * coefficient}


def alg_mul(left, right, node_lambda, kappa, order, cubic_corrections=()):
    answer = {}
    for left_key, left_value in left.items():
        for right_key, right_value in right.items():
            reduced_terms = reduce_monomial(
                tuple(left_key[index] + right_key[index] for index in range(3)),
                left_value * right_value, node_lambda, kappa, order,
                cubic_corrections,
            )
            for key, value in reduced_terms:
                answer[key] = answer.get(key, QQ.zero()) + value
                if not answer[key]:
                    del answer[key]
    return answer


def alg_from_series(value, order):
    return {
        (0, 0, degree): QQ(value[degree])
        for degree in range(order) if value[degree]
    }


def alg_inverse(value, node_lambda, kappa, order, cubic_corrections=()):
    constant = value.get((0, 0, 0), QQ.zero())
    assert constant
    normalized_tail = alg_add(alg_scale(value, 1 / constant), {(0, 0, 0): -1})
    answer = {(0, 0, 0): QQ.one()}
    power = {(0, 0, 0): QQ.one()}
    for exponent in range(1, order):
        power = alg_mul(
            power, normalized_tail, node_lambda, kappa, order,
            cubic_corrections,
        )
        answer = alg_add(answer, alg_scale(power, (-1) ** exponent))
    return alg_scale(answer, 1 / constant)


def node_condition_matrix(support, order):
    precision = order + 2
    A_series = series_at(A, support, precision)
    B_series = series_at(B, support, precision)
    c0 = QQ(-3 * B(support) / (2 * A(support)))
    assert c0 ** 3 + A(support) * c0 + B(support) == 0
    assert 3 * c0 ** 2 + A(support) == 0
    c_series = (-A_series / 3).sqrt()
    if c_series[0] != c0:
        c_series = -c_series
    node_lambda = QQ(3 * c0)
    node_error = c_series ** 3 + A_series * c_series + B_series
    assert node_error.valuation() == 2
    kappa = QQ(node_error[2])

    one = {(0, 0, 0): QQ.one()}
    xi_variable = {(1, 0, 0): QQ.one()}
    y_variable = {(0, 1, 0): QQ.one()}
    x_local = alg_add(alg_from_series(c_series, order), xi_variable)
    y_local = y_variable
    hx_local = alg_from_series(series_at(Hx, support, precision), order)
    hy_local = alg_from_series(series_at(Hy, support, precision), order)
    denominator = alg_add(x_local, alg_scale(hx_local, -1))
    assert denominator.get((0, 0, 0), QQ.zero())
    chord = alg_mul(
        alg_add(y_local, hy_local),
        alg_inverse(denominator, node_lambda, kappa, order),
        node_lambda, kappa, order,
    )
    z_local = alg_from_series(series_at(Z, support, precision), order)
    z_chord = alg_mul(z_local, chord, node_lambda, kappa, order)

    if order == 3:
        quotient_basis = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
                          (2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1), (0, 0, 2)]
    elif order == 2:
        quotient_basis = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    else:
        raise ValueError("only the q114 orders 3 and 2 are supported")
    residues = []
    for AA, BB in ambient:
        aa_local = alg_from_series(series_at(AA, support, precision), order)
        bb_local = alg_from_series(series_at(BB, support, precision), order)
        numerator = alg_add(
            aa_local,
            alg_scale(
                alg_mul(bb_local, z_chord, node_lambda, kappa, order), -1,
            ),
        )
        residues.append([numerator.get(key, QQ.zero()) for key in quotient_basis])
    return matrix(QQ, list(zip(*residues))), {
        "support": str(support), "order": order, "quotient_length": order ** 2,
        "node_x": str(c0), "tangent_square": str(node_lambda), "kappa": str(kappa),
        "quotient_basis": [list(key) for key in quotient_basis],
    }


def marked_node_condition_matrix(support, order=3):
    """Invariant marked-node block from the cleared chord numerator.

    Put G=(x-Hx)AA + BB*Z*(y+Hy).  The old nonidentity component at this
    node is the prescribed new zero, so restriction to the exceptional conic
    must retain the full two-dimensional H0(E,O_E(1)).  Modulo m^4 its two
    exact leading generators are

        tau^2*(xi-alpha0*tau),  tau^2*(y+eta0*tau).

    Keeping only the first generator contracts the proposed zero and produces
    the rejected bidegree-(4,6) genus-two residual.
    """
    local_order = 5
    precision = local_order + 2
    A_series = series_at(A, support, precision)
    B_series = series_at(B, support, precision)
    c0 = QQ(-3 * B(support) / (2 * A(support)))
    c_series = (-A_series / 3).sqrt()
    if c_series[0] != c0:
        c_series = -c_series
    hx_series = series_at(Hx, support, precision)
    hy_series = series_at(Hy, support, precision)
    z_series = series_at(Z, support, precision)
    assert (hx_series - c_series).valuation() == 1
    assert hy_series.valuation() == 1
    assert z_series[0]
    alpha0 = QQ((hx_series - c_series)[1])
    eta0 = QQ(hy_series[1])
    node_error = c_series ** 3 + A_series * c_series + B_series
    assert node_error.valuation() == 2
    kappa = QQ(node_error[2])
    node_lambda = QQ(3 * c0)
    # Modulo m^4 the completed node relation is not just its quadratic
    # tangent cone.  In x=c(tau)+xi coordinates one has exactly
    #
    #   y^2 = 3*c(tau)*xi^2 + xi^3 + node_error(tau).
    #
    # The three total-degree-three corrections below are invisible in the
    # order-two/three unmarked blocks but essential in this marked order-four
    # quotient.
    cubic_corrections = (
        (3, 0, QQ.one()),
        (2, 1, QQ(3 * c_series[1])),
        (0, 3, QQ(node_error[3])),
        (2, 2, QQ(3 * c_series[2])),
        (0, 4, QQ(node_error[4])),
    )
    assert eta0 ** 2 == node_lambda * alpha0 ** 2 + kappa

    xi_variable = {(1, 0, 0): QQ.one()}
    y_variable = {(0, 1, 0): QQ.one()}
    x_local = alg_add(alg_from_series(c_series, local_order), xi_variable)
    hx_local = alg_from_series(hx_series, local_order)
    hy_local = alg_from_series(hy_series, local_order)
    z_local = alg_from_series(z_series, local_order)
    chord_denominator = alg_add(x_local, alg_scale(hx_local, -1))
    chord_numerator = alg_add(y_variable, hy_local)
    quotient_basis = [
        (xi_degree, y_degree, tau_degree)
        for total_degree in range(local_order)
        for y_degree in (0, 1)
        for xi_degree in range(total_degree - y_degree + 1)
        for tau_degree in [total_degree - y_degree - xi_degree]
    ]
    # The admissible residue is a submodule, not a list of coordinate
    # directions.  The two degree-three generators cut out the marked point
    # on the exceptional conic.  Close them under multiplication by the local
    # maximal ideal through degree four.  Conditions are the annihilator of
    # this exact submodule in the 25-dimensional order-five local algebra.
    tau_squared = {(0, 0, 2): QQ.one()}
    zero_generators = [
        alg_mul(
            tau_squared, chord_denominator, node_lambda, kappa, local_order,
            cubic_corrections,
        ),
        alg_mul(
            tau_squared, chord_numerator, node_lambda, kappa, local_order,
            cubic_corrections,
        ),
    ]
    local_variables = [
        {(1, 0, 0): QQ.one()},
        {(0, 1, 0): QQ.one()},
        {(0, 0, 1): QQ.one()},
    ]
    allowed = list(zero_generators)
    for generator in zero_generators:
        for variable in local_variables:
            allowed.append(alg_mul(
                generator, variable, node_lambda, kappa, local_order,
                cubic_corrections,
            ))
    allowed_matrix = matrix(QQ, [
        [value.get(key, QQ.zero()) for key in quotient_basis]
        for value in allowed
    ])
    quotient_annihilator = allowed_matrix.right_kernel_matrix(
        algorithm="generic", basis="computed"
    )

    residues = []
    for AA, BB in ambient:
        aa_local = alg_from_series(series_at(AA, support, precision), local_order)
        bb_local = alg_from_series(series_at(BB, support, precision), local_order)
        cleared = alg_add(
            alg_mul(
                chord_denominator, aa_local, node_lambda, kappa, local_order,
                cubic_corrections,
            ),
            alg_scale(
                alg_mul(
                    alg_mul(
                        bb_local, z_local, node_lambda, kappa, local_order,
                        cubic_corrections,
                    ),
                    chord_numerator, node_lambda, kappa, local_order,
                    cubic_corrections,
                ),
                -1,
            ),
        )
        residues.append([cleared.get(key, QQ.zero()) for key in quotient_basis])
    raw_residue_matrix = matrix(QQ, list(zip(*residues)))
    conditions = quotient_annihilator * raw_residue_matrix
    return conditions, {
        "support": str(support), "exceptional_jet_order": order,
        "raw_jet_rows": len(quotient_basis),
        "node_x": str(c0), "tangent_square": str(node_lambda),
        "kappa": str(kappa), "marked_tangent": [str(alpha0), str(eta0), "1"],
        "quotient_basis": [list(key) for key in quotient_basis],
        "allowed_submodule_dimension": int(
            len(quotient_basis) - quotient_annihilator.nrows()
        ),
        "condition_row_count": int(quotient_annihilator.nrows()),
        "condition_rows": "annihilator of the m-closed order-five reflexive module tau^2*(x-Hx,y+Hy)",
        "allowed_top_residual": "exact two-generator marked-section module, including its degree-four infinitely-near jet",
    }


support_1 = QQ(pointing_artifact["effective_horizontal_components"][0]["child_I2_support"])
support_3 = QQ(pointing_artifact["effective_horizontal_components"][1]["child_I2_support"])
remaining_factor = R(pointing_artifact["remaining_vertical_components"]["child_I2_support_factor_coefficients_low_to_high"])
assert remaining_factor.degree() == 1
support_2 = QQ(-remaining_factor[0] / remaining_factor[1])
node_1, node_1_record = node_condition_matrix(support_1, 3)
node_2, node_2_record = node_condition_matrix(support_2, 2)
node_3, node_3_record = marked_node_condition_matrix(support_3, 3)
restricted_1 = node_1 * smooth_kernel.transpose()
restricted_2 = node_2 * smooth_kernel.transpose()
restricted_3 = node_3 * smooth_kernel.transpose()
restricted_12 = restricted_1.stack(restricted_2)
restricted = restricted_12.stack(restricted_3)
rank_prime = 131
rank_1_modp = restricted_1.change_ring(GF(rank_prime)).rank()
rank_2_modp = restricted_2.change_ring(GF(rank_prime)).rank()
rank_3_modp = restricted_3.change_ring(GF(rank_prime)).rank()
rank_12_modp = restricted_12.change_ring(GF(rank_prime)).rank()
rank_stacked_modp = restricted.change_ring(GF(rank_prime)).rank()
search_marked_plane_modp = "--search-marked-plane-modp" in sys.argv
print(
    f"FIXEDREVERSE4A1RR_NODE_STAGE|p={rank_prime}|"
    f"individual={rank_1_modp},{rank_2_modp},{rank_3_modp}|"
    f"stack12={rank_12_modp}|stack123={rank_stacked_modp}",
    flush=True,
)
restricted_12_modp = restricted_12.change_ring(GF(rank_prime))
marked_modp = restricted_3.change_ring(GF(rank_prime))
marked_increment_rows = []
running_modp = restricted_12_modp
running_rank = rank_12_modp
for row_index in range(marked_modp.nrows()):
    trial = running_modp.stack(matrix(GF(rank_prime), [marked_modp.row(row_index)]))
    trial_rank = trial.rank()
    if trial_rank > running_rank:
        marked_increment_rows.append(row_index)
        running_modp = trial
        running_rank = trial_rank
print(
    f"FIXEDREVERSE4A1RR_MARKED_PIVOTS|rows={marked_increment_rows}|"
    f"basis={node_3_record.get('quotient_basis')}",
    flush=True,
)
assert (rank_1_modp, rank_2_modp) == (5, 3)
assert rank_12_modp == 8
assert rank_stacked_modp == (12 if search_marked_plane_modp else 13)
# This is a tiny 17-by-15 field kernel with very large rational entries.
# Clearing all denominators routes Sage through IML and is much slower than
# direct pivot elimination over QQ.
resolved_coordinates = restricted.right_kernel_matrix(
    algorithm="generic", basis="computed"
)
kernel = resolved_coordinates * smooth_kernel
assert kernel.nrows() == (3 if search_marked_plane_modp else 2)

pairs = []
for row in kernel.rows():
    AA = sum((row[index] * ambient[index][0] for index in range(57)), R.zero())
    BB = sum((row[index] * ambient[index][1] for index in range(57)), R.zero())
    assert (AA * X - BB * Y) % collision_modulus == 0
    # Nullspace bases are only defined up to a scalar.  Clear the common
    # rational content now; carrying it into a fourth power makes otherwise
    # modest bivariate multiplication prohibitively expensive.
    pair_values = AA.list() + BB.list()
    common_denominator = ZZ.one()
    for value in pair_values:
        common_denominator = common_denominator.lcm(QQ(value).denominator())
    integral_values = [ZZ(QQ(value) * common_denominator) for value in pair_values]
    common_numerator = ZZ.zero()
    for value in integral_values:
        common_numerator = common_numerator.gcd(value)
    primitive_scale = QQ(common_denominator) / abs(common_numerator)
    AA = R(primitive_scale * AA)
    BB = R(primitive_scale * BB)
    assert all(QQ(value).denominator() == 1 for value in AA.list() + BB.list())
    assert (AA * X - BB * Y) % collision_modulus == 0
    pairs.append((AA, BB))

# Search the planes in the exact three-space modulo one good prime.  For every
# plane, form the ternary chord quartic once, strip its univariate square
# content, and retain only binary-quartic residuals.  This is a construction
# selector; a retained plane must still be lifted and replayed exactly.
if search_marked_plane_modp:
    p = rank_prime
    Fp = GF(p)
    Rp = PolynomialRing(Fp, "up")
    up = Rp.gen()

    def mod_rational(value):
        value = QQ(value)
        return Fp(value.numerator()) / Fp(value.denominator())

    def mod_poly(poly):
        return Rp([mod_rational(value) for value in R(poly).list()])

    Ternary = PolynomialRing(Rp, names=("z0", "z1", "z2"))
    z0, z1, z2 = Ternary.gens()
    coordinates = (z0, z1, z2)
    aa_ternary = sum(
        (Ternary(mod_poly(pairs[index][0])) * coordinates[index] for index in range(3)),
        Ternary.zero(),
    )
    bb_ternary = sum(
        (Ternary(mod_poly(pairs[index][1])) * coordinates[index] for index in range(3)),
        Ternary.zero(),
    )
    Xp, Yp, Zp, Ap = map(mod_poly, (X, Y, Z, A))
    raw_ternary = (
        aa_ternary ** 4
        - 6 * Ternary(Xp) * aa_ternary ** 2 * bb_ternary ** 2
        + 8 * Ternary(Yp) * aa_ternary * bb_ternary ** 3
        - 3 * Ternary(Xp ** 2) * bb_ternary ** 4
        - 4 * Ternary(Ap * Zp ** 4) * bb_ternary ** 4
    )
    collision_power = Zp ** 4
    ternary_coefficients = {}
    for exponent, value in raw_ternary.dict().items():
        quotient, remainder = Rp(value).quo_rem(collision_power)
        assert not remainder
        ternary_coefficients[tuple(map(int, exponent))] = quotient
    assert set(sum(exponent) for exponent in ternary_coefficients) == {4}

    def maximal_square_divisor(poly):
        layers = [Rp(poly)]
        while True:
            repeated = layers[-1].gcd(layers[-1].derivative())
            if repeated.is_constant():
                break
            layers.append(repeated)
        answer = Rp.one()
        for index in range(1, len(layers), 2):
            denominator = layers[index + 1] if index + 1 < len(layers) else Rp.one()
            quotient, remainder = layers[index].quo_rem(denominator)
            assert not remainder
            answer *= quotient
        return answer.monic()

    def score_plane(first, second):
        coefficients_by_degree = [Rp.zero() for unused in range(5)]
        for exponent, value in ternary_coefficients.items():
            i, j, k = exponent
            # Substitute z_r = first_r*S + second_r*T and collect S-degree.
            partial = {(0, 0): Fp.one()}
            for power, left, right in zip((i, j, k), first, second):
                updated = {}
                for (s_degree, t_degree), scalar in partial.items():
                    for left_count in range(power + 1):
                        right_count = power - left_count
                        contribution = (
                            scalar * Fp(ZZ(power).binomial(left_count))
                            * left ** left_count * right ** right_count
                        )
                        key = (s_degree + left_count, t_degree + right_count)
                        updated[key] = updated.get(key, Fp.zero()) + contribution
                partial = updated
            for (s_degree, t_degree), scalar in partial.items():
                assert s_degree + t_degree == 4
                coefficients_by_degree[s_degree] += scalar * value
        content = Rp.zero()
        for value in coefficients_by_degree:
            content = content.gcd(value)
        content = content.monic()
        square = maximal_square_divisor(content)
        reduced = []
        for value in coefficients_by_degree:
            quotient, remainder = value.quo_rem(square ** 2)
            assert not remainder
            reduced.append(quotient)
        odd_content = Rp.zero()
        for value in reduced:
            odd_content = odd_content.gcd(value)
        odd_content = odd_content.monic()
        moving = []
        for value in reduced:
            quotient, remainder = value.quo_rem(odd_content)
            assert not remainder
            moving.append(quotient)
        return {
            "square_degree": int(square.degree()),
            "odd_content_degree": int(odd_content.degree()),
            "moving_degree": max(int(value.degree()) for value in moving),
        }

    winners = []
    scanned = 0
    # Main affine chart: z2=alpha*z0+beta*z1.
    for alpha in Fp:
        for beta in Fp:
            score = score_plane((Fp.one(), Fp.zero(), alpha), (Fp.zero(), Fp.one(), beta))
            scanned += 1
            if score["moving_degree"] == 4 and score["odd_content_degree"] in (0, 2):
                winners.append({"chart": "z2", "alpha": int(alpha), "beta": int(beta), **score})
        if int(alpha) % 16 == 0:
            print(
                f"FIXEDREVERSE4A1RR_PLANE_SEARCH_PROGRESS|p={p}|alpha={int(alpha)}|"
                f"scanned={scanned}|winners={len(winners)}",
                flush=True,
            )
    # Boundary chart z1=alpha*z0, plus z0=0.
    for alpha in Fp:
        score = score_plane((Fp.one(), alpha, Fp.zero()), (Fp.zero(), Fp.zero(), Fp.one()))
        scanned += 1
        if score["moving_degree"] == 4 and score["odd_content_degree"] in (0, 2):
            winners.append({"chart": "z1", "alpha": int(alpha), **score})
    score = score_plane((Fp.zero(), Fp.one(), Fp.zero()), (Fp.zero(), Fp.zero(), Fp.one()))
    scanned += 1
    if score["moving_degree"] == 4 and score["odd_content_degree"] in (0, 2):
        winners.append({"chart": "z0", **score})
    print(
        "FIXEDREVERSE4A1RR_PLANE_SEARCH|"
        f"p={p}|scanned={scanned}|winners={json.dumps(winners, sort_keys=True)}",
        flush=True,
    )
    raise SystemExit(0)

# Reduce the two-row integral pencil lattice before any fourth powers are
# formed.  This is only a GL_2(ZZ) change of the new base coordinates and can
# sharply lower coefficient growth without changing the H0 plane.
pre_lll_pair_bits = [rational_bits(AA.list() + BB.list()) for AA, BB in pairs]
pair_lattice = matrix(ZZ, [
    [ZZ(AA[degree]) for degree in range(aa_degree + 1)]
    + [ZZ(BB[degree]) for degree in range(bb_degree + 1)]
    for AA, BB in pairs
])
reduced_pair_lattice = pair_lattice.LLL()
pivot_columns = next(
    (left, right)
    for left in range(pair_lattice.ncols())
    for right in range(left + 1, pair_lattice.ncols())
    if matrix(QQ, pair_lattice.matrix_from_columns([left, right])).det()
)
old_pivot = matrix(QQ, pair_lattice.matrix_from_columns(pivot_columns))
new_pivot = matrix(QQ, reduced_pair_lattice.matrix_from_columns(pivot_columns))
pair_change = new_pivot * old_pivot.inverse()
assert pair_change * pair_lattice == reduced_pair_lattice
assert abs(pair_change.det()) == 1
pairs = [
    (
        R(list(row[:aa_degree + 1])),
        R(list(row[aa_degree + 1:])),
    )
    for row in reduced_pair_lattice.rows()
]
for AA, BB in pairs:
    assert (AA * X - BB * Y) % collision_modulus == 0
AA0, BB0 = pairs[0]
AA1, BB1 = pairs[1]
print(
    "FIXEDREVERSE4A1RR_KERNEL_STAGE|"
    f"pre_lll_pair_bits={pre_lll_pair_bits}|"
    f"primitive_pair_bits={[rational_bits(AA.list() + BB.list()) for AA, BB in pairs]}|"
    f"degrees={[(AA.degree(), BB.degree()) for AA, BB in pairs]}",
    flush=True,
)
if "--stop-after-kernel" in sys.argv:
    raise SystemExit(0)

# Compile the chord radicand with s as the new base parameter.  Its s-degree is
# only four, whereas the u-coefficients are enormous.  A nested ring QQ[u][s]
# therefore delegates every costly coefficient product to FLINT's optimized
# univariate arithmetic and avoids both QQ(s) gcd normalization and generic
# bivariate schoolbook multiplication.
S = PolynomialRing(QQ, "s")
s = S.gen()
Pencil = PolynomialRing(R, "v")
v = Pencil.gen()
SU = PolynomialRing(QQ, names=("s", "u"))
ss, uu = SU.gens()


def lift_pencil(poly):
    return Pencil([R(poly)])


aa = Pencil([AA0, AA1])
bb = Pencil([BB0, BB1])
Xu, Yu, Zu, Au = map(lift_pencil, (X, Y, Z, A))
pair_fingerprint = hashlib.sha256(json.dumps([
    [coefficients(AA), coefficients(BB)] for AA, BB in pairs
], separators=(",", ":")).encode()).hexdigest()
RADICAND_CHECKPOINT = LOCAL / "fixed-reverse-4a1-after-collision-qq.json"
radicand_checkpoint = json.loads(RADICAND_CHECKPOINT.read_text()) if RADICAND_CHECKPOINT.exists() else None
if (
    radicand_checkpoint
    and radicand_checkpoint.get("status") == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_AFTER_COLLISION"
    and radicand_checkpoint.get("pair_fingerprint") == pair_fingerprint
):
    after_collision = Pencil([
        R(values) for values in radicand_checkpoint["coefficients_in_s_then_u_low_to_high"]
    ])
else:
    raw = aa ** 4 - 6 * Xu * aa ** 2 * bb ** 2 + 8 * Yu * aa * bb ** 3 - 3 * Xu ** 2 * bb ** 4 - 4 * Au * bb ** 4 * Zu ** 4
    after_collision, remainder = raw.quo_rem(Zu ** 4)
    assert not remainder
    RADICAND_CHECKPOINT.write_text(json.dumps({
        "schema": "elkies-k3.fixed-reverse-4a1-after-collision-qq.v1",
        "status": "PASS_EXACT_QQ_FIXED_REVERSE_4A1_AFTER_COLLISION",
        "pair_fingerprint": pair_fingerprint,
        "coefficients_in_s_then_u_low_to_high": [
            coefficients(R(value)) for value in after_collision.list()
        ],
    }, indent=2, sort_keys=True) + "\n")

# Strip the forced square by multivariate gcd rather than factoring the large
# bivariate radicand.  This is exact in QQ[s,u].
after = sum(
    (
        QQ(coefficient) * ss ** s_degree * uu ** u_degree
        for s_degree, u_polynomial in enumerate(after_collision.list())
        for u_degree, coefficient in enumerate(R(u_polynomial).list())
    ),
    SU.zero(),
)
# Extract the maximal square divisor with multiplicities.  The one-gcd shortcut
# gcd(f,f_s,f_u)^2 is valid only when every repeated factor has multiplicity
# exactly two; the vertical coefficients 3 and 2 on this edge create higher
# even powers.  Successive derivative-gcd layers d_i carry exponent m-i for a
# factor of multiplicity m, so (d_1/d_2)(d_3/d_4)... has exponent floor(m/2).
derivative_layers = [after]
while True:
    current = derivative_layers[-1]
    repeated_layer = current.gcd(current.derivative(ss)).gcd(current.derivative(uu))
    if repeated_layer.is_constant():
        break
    derivative_layers.append(repeated_layer)
    print(
        "FIXEDREVERSE4A1RR_SQUARE_LAYER|"
        f"index={len(derivative_layers) - 1}|degrees=({repeated_layer.degree(ss)},{repeated_layer.degree(uu)})",
        flush=True,
    )
square_factor = SU.one()
for index in range(1, len(derivative_layers), 2):
    denominator = derivative_layers[index + 1] if index + 1 < len(derivative_layers) else SU.one()
    layer_quotient, remainder = derivative_layers[index].quo_rem(denominator)
    assert not remainder
    square_factor *= layer_quotient
squarefree_residual, remainder = after.quo_rem(square_factor ** 2)
assert not remainder
# The resolved divisor also has a fixed vertical branch factor with odd
# multiplicity one.  It is the content in QQ[u] common to all five pencil
# coefficients.  Removing this fixed component turns the residual sextic into
# the moving binary quartic; compute it by gcd, not factorization.
residual_by_s_degree = [R.zero() for unused in range(squarefree_residual.degree(ss) + 1)]
for (s_degree, u_degree), coefficient in squarefree_residual.dict().items():
    residual_by_s_degree[s_degree] += QQ(coefficient) * u ** u_degree
fixed_branch_factor = R.zero()
for value in residual_by_s_degree:
    fixed_branch_factor = fixed_branch_factor.gcd(value)
fixed_branch_factor = fixed_branch_factor.monic()
assert fixed_branch_factor.degree() in (0, 2)
fixed_branch_factor_SU = sum(
    (QQ(value) * uu ** degree for degree, value in enumerate(fixed_branch_factor.list())),
    SU.zero(),
)
quartic, remainder = squarefree_residual.quo_rem(fixed_branch_factor_SU)
assert not remainder
print(
    "FIXEDREVERSE4A1RR_RESIDUAL_DIAGNOSTIC|"
    f"squarefree_degrees=({squarefree_residual.degree(ss)},{squarefree_residual.degree(uu)})|"
    f"fixed_branch_degree={fixed_branch_factor.degree()}|"
    f"quotient_degrees=({quartic.degree(ss)},{quartic.degree(uu)})",
    flush=True,
)
for diagnostic_prime in (131, 137):
    try:
        DiagnosticRing = PolynomialRing(GF(diagnostic_prime), names=("s", "u"))
        diagnostic_s, diagnostic_u = DiagnosticRing.gens()
        diagnostic_residual = DiagnosticRing(squarefree_residual)
        diagnostic_factor_degrees = [
            (int(factor.degree(diagnostic_s)), int(factor.degree(diagnostic_u)), int(multiplicity))
            for factor, multiplicity in diagnostic_residual.factor()
        ]
        print(
            "FIXEDREVERSE4A1RR_RESIDUAL_MODP|"
            f"p={diagnostic_prime}|factor_degrees={diagnostic_factor_degrees}",
            flush=True,
        )
    except (ZeroDivisionError, TypeError, ValueError):
        print(
            f"FIXEDREVERSE4A1RR_RESIDUAL_MODP|p={diagnostic_prime}|bad_denominator=1",
            flush=True,
        )
assert quartic.degree(uu) == 4
assert quartic.gcd(quartic.derivative(ss)).gcd(quartic.derivative(uu)).is_constant()

U_over_S = PolynomialRing(S, "u")
quartic_univariate = U_over_S(quartic)
square_factor_univariate = U_over_S(square_factor)
quartic_coefficients = [S(quartic_univariate[degree]) for degree in range(5)]
e, d, c, b, a = quartic_coefficients
I = S(12 * a * e - 3 * b * d + c ** 2)
J = S(72 * a * c * e + 9 * b * c * d - 27 * a * d ** 2 - 27 * b ** 2 * e - 2 * c ** 3)
A_child = S(-27 * I)
B_child = S(-27 * J)

# Remove any finite fourth/sixth common scaling.  Factoring only gcd(A,B) is
# tiny compared with factoring the chord radicand or the surface ideal.
removed_scalings = []
for factor, unused in A_child.gcd(B_child).factor():
    order = min(A_child.valuation(factor) // 4, B_child.valuation(factor) // 6)
    if order:
        A_child //= factor ** (4 * order)
        B_child //= factor ** (6 * order)
        removed_scalings.append((factor, int(order)))
assert A_child.degree() <= 8 and B_child.degree() <= 12
print(
    "FIXEDREVERSE4A1RR_QUARTIC_STAGE|"
    f"square_factor_total_degree={square_factor.total_degree()}|"
    f"quartic_degrees_s_u=({quartic.degree(ss)},{quartic.degree(uu)})|"
    f"preminimal_degrees=({A_child.degree()},{B_child.degree()})|"
    f"removed_scalings={[(factor.degree(), order) for factor, order in removed_scalings]}",
    flush=True,
)
assert not (A_child.degree() <= 4 and B_child.degree() <= 6)

Delta = S(-16 * (4 * A_child ** 3 + 27 * B_child ** 2))
assert Delta.degree() <= 24

# Semistable classification without factoring the degree-24 discriminant.
# The repeated part gives the expected I2 supports; the residual factor is
# squarefree and disjoint from A and B, hence gives geometric I1 fibres.
repeated = Delta.gcd(Delta.derivative()).monic()
expected_i2_count = 4
assert repeated.degree() == expected_i2_count
residual, remainder = Delta.quo_rem(repeated ** 2)
assert not remainder and residual.degree() == 24 - 2 * expected_i2_count
assert repeated.is_squarefree()
assert residual.gcd(residual.derivative()).is_constant()
assert repeated.gcd(residual).is_constant()
assert repeated.gcd(A_child).is_constant() and repeated.gcd(B_child).is_constant()
assert residual.gcd(A_child).is_constant() and residual.gcd(B_child).is_constant()
infinity_orders = [int(8 - A_child.degree()), int(12 - B_child.degree()), int(24 - Delta.degree())]
assert infinity_orders == [0, 0, 0]

quartic_values = []
for coefficient in quartic_coefficients:
    quartic_values.extend(coefficient.list())
jacobian_values = A_child.list() + B_child.list() + Delta.list()
payload = {
    "schema": "elkies-k3.fixed-reverse-4a1-resolved-rr-qq.v1",
    "status": "PASS_EXACT_QQ_FIXED_REVERSE_4A1_RR_JACOBIAN",
    "reproducing_command": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/compile_h92_fixed_reverse_4a1_rr_qq.sage"
    ),
    "inputs": {
        "target": {"path": str(TARGET.relative_to(ROOT)), "sha256": sha256(TARGET)},
        "surface": {"path": str(SURFACE.relative_to(ROOT)), "sha256": sha256(SURFACE)},
        "pointing": {"path": str(POINTING.relative_to(ROOT)), "sha256": sha256(POINTING)},
    },
    "divisor": {
        "class": "O+P+3*C1+2*C2-3F", "P_dot_O": 21,
        "fibre_twist": -3, "vertical_support_coefficients": [3, 2, 0],
    },
    "smooth_RR": {
        "ambient_dimension": 57,
        "collision_modulus_degree": 42,
        "collision_condition_rank": 42,
        "collision_solver": "X inverse modulo Z^2 plus a 2-by-17 endpoint kernel",
        "post_collision_dimension": int(smooth_kernel.nrows()),
        "resolved_node_condition_ranks": [5, 3, 5],
        "resolved_node_stacked_rank": 13,
        "rank_witness_prime": rank_prime,
        "total_condition_rank": 55,
        "h0": int(kernel.nrows()),
        "AA_degree_bound": aa_degree,
        "BB_degree_bound": bb_degree,
        "basis_pairs": [
            {"AA_coefficients_low_to_high": coefficients(AA), "BB_coefficients_low_to_high": coefficients(BB)}
            for AA, BB in pairs
        ],
        "resolved_node_quotients": [node_1_record, node_2_record, node_3_record],
        "complete_resolved_vertical_rows": True,
    },
    "binary_quartic": {
        "coefficients_in_old_u_low_to_high": [coefficients(value) for value in quartic_coefficients],
        "square_factor_coefficients_in_old_u_low_to_high": [
            coefficients(S(square_factor_univariate[degree]))
            for degree in range(square_factor_univariate.degree() + 1)
        ],
        "square_factor_total_degree": int(square_factor.total_degree()),
        "removed_fixed_branch_factor_coefficients_in_old_u_low_to_high": coefficients(fixed_branch_factor),
        "maximum_rational_bits": rational_bits(quartic_values),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": coefficients(A_child),
        "minimal_B_coefficients_low_to_high": coefficients(B_child),
        "discriminant_coefficients_low_to_high": coefficients(Delta),
        "degrees_A_B_Delta": [int(A_child.degree()), int(B_child.degree()), int(Delta.degree())],
        "finite_fibres": [
            {"factor_coefficients_low_to_high": coefficients(repeated), "factor_degree": int(repeated.degree()), "orders_A_B_Delta": [0, 0, 2], "kodaira": "I2", "root_rank_contribution": int(repeated.degree())},
            {"factor_coefficients_low_to_high": coefficients(residual), "factor_degree": int(residual.degree()), "orders_A_B_Delta": [0, 0, 1], "kodaira": "I1", "root_rank_contribution": 0},
        ],
        "infinity": {"orders_A_B_Delta": infinity_orders, "kodaira": "smooth"},
        "euler_number": 24,
        "root_rank": expected_i2_count,
        "ADE": "4A1",
        "MW_rank_if_rho19": 13,
        "removed_nonminimal_finite_scalings": [
            {"factor_coefficients_low_to_high": coefficients(factor), "order": order}
            for factor, order in removed_scalings
        ],
        "maximum_A_B_Delta_rational_bits": rational_bits(jacobian_values),
    },
    "method": {
        "exact_resolved_A1_node_ideal_quotients": True,
        "resolved_rank_sequence": [57, 15, 10, 7, 2],
        "multivariate_gcd_square_stripping": True,
        "full_discriminant_factorization": False,
        "groebner_or_surface_elimination": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": "Exact QQ resolved h0=2, binary quartic, minimal semistable Jacobian and fibre classification. Pointing and component marking are separate gates.",
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE4A1RR|edge=3A1/MW14 reverse to 4A1/MW13|ambient=57|collision_rank={}|resolved_ranks=5+3+5|h0={}|quartic=4|degrees={}|I2={}|I1={}|"
    "root_rank={}|bits={}|seconds={:.3f}|status={}|output={}".format(
        42, kernel.nrows(), payload["child"]["degrees_A_B_Delta"],
        expected_i2_count, residual.degree(), expected_i2_count,
        payload["child"]["maximum_A_B_Delta_rational_bits"], payload["method"]["runtime_seconds"],
        payload["status"], OUTPUT,
    ), flush=True,
)
