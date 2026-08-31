#!/usr/bin/env sage
"""Compile the physical-nef fixed reverse q4/orbit114 pencil over QQ.

After the exact four-reflection physical reduction, the divisor is

    D = O + P - C1 - 9F,   P.O=21,

where C1 is the first old nonidentity I2 component.  Thus the compact chord
ambient has dimension 45, smooth collision leaves dimension three, and the
single exceptional-component vanishing condition leaves h0=2.  No Groebner
basis or generic surface elimination is used.
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
NEF_AUDIT = LOCAL / "fixed-reverse-4a1-physical-nef-audit.json"
SURFACE = LOCAL / "fixed-reverse-3a1-rr-qq.json"
POINTING = LOCAL / "fixed-reverse-3a1-pointing-qq.json"
OUTPUT = LOCAL / "fixed-reverse-4a1-rr-qq.json"


def read_json(path):
    return json.loads(path.read_text())


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coefficients(poly):
    return [str(value) for value in poly.list()]


def rational_bits(values):
    answer = 0
    for value in values:
        value = QQ(value)
        answer = max(
            answer,
            abs(ZZ(value.numerator())).nbits(),
            ZZ(value.denominator()).nbits(),
        )
    return int(answer)


started = time.monotonic()
target = read_json(TARGET)
nef_audit = read_json(NEF_AUDIT)
surface = read_json(SURFACE)
pointing = read_json(POINTING)
assert target["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_HORIZONTAL_ON_3A1"
assert nef_audit["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_PHYSICAL_NEF"
assert surface["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_3A1_RR_JACOBIAN"
assert pointing["status"] == "PASS_EXACT_QQ_FIXED_REVERSE_3A1_POINTING"
assert nef_audit["physical_component_reduction"]["q"] == 24

R = PolynomialRing(QQ, "u")
u = R.gen()
A = R(surface["child"]["minimal_A_coefficients_low_to_high"])
B = R(surface["child"]["minimal_B_coefficients_low_to_high"])
section = target["section"]
X = R(section["x_numerator_coefficients_low_to_high"])
Y = R(section["y_numerator_coefficients_low_to_high"])
Z = R(section["Z_coefficients_low_to_high"])
assert (X.degree(), Y.degree(), Z.degree()) == (46, 69, 21)
assert Y ** 2 == X ** 3 + A * X * Z ** 4 + B * Z ** 6
projective_scale = Z.leading_coefficient()
X = R(X / projective_scale ** 2)
Y = R(Y / projective_scale ** 3)
Z = R(Z / projective_scale)
assert Z.is_monic()
assert Y ** 2 == X ** 3 + A * X * Z ** 4 + B * Z ** 6

# D=O+P-C1-9F.  First compile O+P-9F; subtracting C1 is the final
# one-dimensional exceptional-component vanishing quotient.
aa_degree = 33
bb_degree = 10
ambient_dimension = aa_degree + bb_degree + 2
collision_modulus = Z ** 2
assert ambient_dimension == 45 and collision_modulus.degree() == 42

SMOOTH_CHECKPOINT = LOCAL / "fixed-reverse-4a1-physical-smooth-kernel-qq.json"
smooth_checkpoint = read_json(SMOOTH_CHECKPOINT) if SMOOTH_CHECKPOINT.exists() else None
if (
    smooth_checkpoint
    and smooth_checkpoint.get("status") == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_PHYSICAL_SMOOTH_KERNEL"
    and smooth_checkpoint.get("target_sha256") == sha256(TARGET)
    and smooth_checkpoint.get("surface_sha256") == sha256(SURFACE)
):
    smooth_pairs = [
        (R(record["AA"]), R(record["BB"]))
        for record in smooth_checkpoint["pairs"]
    ]
    endpoint_rank = 8
else:
    X_inverse = X.inverse_mod(collision_modulus)
    bb_monomials = [u ** degree for degree in range(bb_degree + 1)]
    aa_remainders = [R((BB * Y * X_inverse) % collision_modulus) for BB in bb_monomials]
    endpoint_conditions = matrix(QQ, [
        [remainder[degree] for remainder in aa_remainders]
        for degree in range(aa_degree + 1, collision_modulus.degree())
    ])
    assert endpoint_conditions.dimensions() == (8, 11)
    endpoint_rank = endpoint_conditions.rank()
    bb_kernel = endpoint_conditions.right_kernel_matrix(
        algorithm="generic", basis="computed"
    )
    assert bb_kernel.dimensions() == (3, 11)
    smooth_pairs = []
    for row in bb_kernel.rows():
        BB = sum((row[index] * bb_monomials[index] for index in range(11)), R.zero())
        AA = R((BB * Y * X_inverse) % collision_modulus)
        assert AA.degree() <= aa_degree
        assert (AA * X - BB * Y) % collision_modulus == 0
        smooth_pairs.append((AA, BB))
    SMOOTH_CHECKPOINT.write_text(json.dumps({
        "schema": "elkies-k3.fixed-reverse-4a1-physical-smooth-kernel-qq.v1",
        "status": "PASS_EXACT_QQ_FIXED_REVERSE_4A1_PHYSICAL_SMOOTH_KERNEL",
        "target_sha256": sha256(TARGET),
        "surface_sha256": sha256(SURFACE),
        "pairs": [
            {"AA": coefficients(AA), "BB": coefficients(BB)}
            for AA, BB in smooth_pairs
        ],
    }, indent=2, sort_keys=True) + "\n")
assert endpoint_rank == 8 and len(smooth_pairs) == 3

# Use a short integral basis of the entire three-space before selecting its
# physical hyperplane.  This makes the modular plane normals suitable for
# rational reconstruction and lowers the eventual quartic coefficient cost.
primitive_smooth_rows = []
for AA, BB in smooth_pairs:
    entries = AA.list() + [QQ.zero()] * (aa_degree + 1 - len(AA.list()))
    entries += BB.list() + [QQ.zero()] * (bb_degree + 1 - len(BB.list()))
    denominator = ZZ.one()
    for value in entries:
        denominator = denominator.lcm(QQ(value).denominator())
    integral = [ZZ(QQ(value) * denominator) for value in entries]
    content = ZZ.zero()
    for value in integral:
        content = content.gcd(value)
    primitive_smooth_rows.append([value // abs(content) for value in integral])
smooth_lattice = matrix(ZZ, primitive_smooth_rows).LLL()
smooth_pairs = [
    (R(list(row[:aa_degree + 1])), R(list(row[aa_degree + 1:])))
    for row in smooth_lattice.rows()
]
assert all((AA * X - BB * Y) % collision_modulus == 0 for AA, BB in smooth_pairs)


if "--probe-local-functionals-modp" in sys.argv:
    p = 131
    Fp = GF(p)

    def normalize_row(row):
        reduced = [Fp(QQ(value).numerator()) / Fp(QQ(value).denominator()) for value in row]
        pivot = next((value for value in reversed(reduced) if value), None)
        return None if pivot is None else [int(value / pivot) for value in reduced]

    records = []
    supports = [
        QQ(record["child_I2_support"])
        for record in pointing["effective_horizontal_components"]
    ]
    remaining_factor = R(pointing["remaining_vertical_components"][
        "child_I2_support_factor_coefficients_low_to_high"
    ])
    supports.insert(1, QQ(-remaining_factor[0] / remaining_factor[1]))
    for support_index, local_support in enumerate(supports, 1):
        PS = PowerSeriesRing(QQ, "tau", default_prec=7)
        tau = PS.gen()

        def at_series(poly):
            answer = PS.zero()
            for coefficient in reversed(R(poly).list()):
                answer = answer * (local_support + tau) + coefficient
            return answer

        A_series = at_series(A)
        c0 = QQ(-3 * B(local_support) / (2 * A(local_support)))
        c_series = (-A_series / 3).sqrt()
        if c_series[0] != c0:
            c_series = -c_series
        hx_series = at_series(X) / at_series(Z) ** 2
        hy_series = at_series(Y) / at_series(Z) ** 3
        z_series = at_series(Z)
        den_series = c_series - hx_series
        cleared_series = [
            at_series(AA) * den_series - at_series(BB) * z_series * hy_series
            for AA, BB in smooth_pairs
        ]
        cleared_plus_series = [
            at_series(AA) * den_series + at_series(BB) * z_series * hy_series
            for AA, BB in smooth_pairs
        ]
        regular_chord_series = [value / den_series for value in cleared_series]
        bbz_series = [at_series(BB) * z_series for unused, BB in smooth_pairs]
        for order in range(5):
            records.append({
                "support": support_index,
                "functional": "cleared_constant_tau_{}".format(order),
                "normal": normalize_row([value[order] for value in cleared_series]),
            })
            records.append({
                "support": support_index,
                "functional": "cleared_plus_tau_{}".format(order),
                "normal": normalize_row([value[order] for value in cleared_plus_series]),
            })
            records.append({
                "support": support_index,
                "functional": "AA_tau_{}".format(order),
                "normal": normalize_row([at_series(AA)[order] for AA, unused in smooth_pairs]),
            })
            records.append({
                "support": support_index,
                "functional": "BB_tau_{}".format(order),
                "normal": normalize_row([at_series(BB)[order] for unused, BB in smooth_pairs]),
            })
            records.append({
                "support": support_index,
                "functional": "BBZ_tau_{}".format(order),
                "normal": normalize_row([value[order] for value in bbz_series]),
            })
            records.append({
                "support": support_index,
                "functional": "regular_chord_tau_{}".format(order),
                "normal": normalize_row([value[order] for value in regular_chord_series]),
            })
        alpha = QQ((hx_series - c_series)[1])
        eta = QQ(hy_series[1])
        # Leading cleared numerator AA*(xi-alpha*tau)-BB*Z*(y+eta*tau).
        for tangent_name, tangent in (
            ("plus_section_tangent", (alpha, eta, QQ.one())),
            ("minus_section_tangent", (alpha, -eta, QQ.one())),
        ):
            xi_value, y_value, tau_value = tangent
            row = []
            for AA, BB in smooth_pairs:
                row.append(
                    AA(local_support) * (xi_value - alpha * tau_value)
                    - BB(local_support) * Z(local_support)
                    * (y_value + eta * tau_value)
                )
            records.append({
                "support": support_index,
                "functional": tangent_name,
                "normal": normalize_row(row),
            })
        tangent_square_modp = Fp(QQ(3 * c0).numerator()) / Fp(QQ(3 * c0).denominator())
        if tangent_square_modp.is_square():
            tangent_root_modp = tangent_square_modp.sqrt()
            for sign in (-1, 1):
                branch_slope = Fp(sign) * tangent_root_modp
                row = []
                for AA, BB in smooth_pairs:
                    aa0 = Fp(QQ(AA(local_support)).numerator()) / Fp(QQ(AA(local_support)).denominator())
                    bb0 = Fp(QQ(BB(local_support)).numerator()) / Fp(QQ(BB(local_support)).denominator())
                    z0 = Fp(QQ(Z(local_support)).numerator()) / Fp(QQ(Z(local_support)).denominator())
                    row.append(aa0 - bb0 * z0 * branch_slope)
                pivot = next(value for value in reversed(row) if value)
                records.append({
                    "support": support_index,
                    "functional": "old_fibre_branch_{}".format("plus" if sign == 1 else "minus"),
                    "normal": [int(value / pivot) for value in row],
                })
    print(
        "FIXEDREVERSE4A1PHYSRR_LOCAL_NORMALS|p={}|records={}".format(
            p, json.dumps(records, sort_keys=True)
        ),
        flush=True,
    )
    raise SystemExit(0)


search_plane_arguments = [
    argument for argument in sys.argv if argument.startswith("--search-plane-modp")
]
if search_plane_arguments:
    p = (
        int(search_plane_arguments[0].split("=", 1)[1])
        if "=" in search_plane_arguments[0]
        else 131
    )
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
    aa_ternary = sum((
        Ternary(mod_poly(smooth_pairs[index][0])) * coordinates[index]
        for index in range(3)
    ), Ternary.zero())
    bb_ternary = sum((
        Ternary(mod_poly(smooth_pairs[index][1])) * coordinates[index]
        for index in range(3)
    ), Ternary.zero())
    Xp, Yp, Zp, Ap = map(mod_poly, (X, Y, Z, A))
    raw_ternary = (
        aa_ternary ** 4
        - 6 * Ternary(Xp) * aa_ternary ** 2 * bb_ternary ** 2
        + 8 * Ternary(Yp) * aa_ternary * bb_ternary ** 3
        - 3 * Ternary(Xp ** 2) * bb_ternary ** 4
        - 4 * Ternary(Ap * Zp ** 4) * bb_ternary ** 4
    )
    ternary_coefficients = {}
    for exponent, value in raw_ternary.dict().items():
        quotient, remainder = Rp(value).quo_rem(Zp ** 4)
        assert not remainder
        ternary_coefficients[tuple(map(int, exponent))] = quotient

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
            partial = {(0, 0): Fp.one()}
            for power, left, right in zip(exponent, first, second):
                updated = {}
                for (left_degree, right_degree), scalar in partial.items():
                    for left_count in range(power + 1):
                        right_count = power - left_count
                        key = (left_degree + left_count, right_degree + right_count)
                        updated[key] = updated.get(key, Fp.zero()) + (
                            scalar * ZZ(power).binomial(left_count)
                            * left ** left_count * right ** right_count
                        )
                partial = updated
            for (left_degree, unused), scalar in partial.items():
                coefficients_by_degree[left_degree] += scalar * value
        content = Rp.zero()
        for value in coefficients_by_degree:
            content = content.gcd(value)
        square = maximal_square_divisor(content.monic())
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
        return int(square.degree()), int(odd_content.degree()), max(
            int(value.degree()) for value in moving
        )

    winners = []
    scanned = 0
    for alpha in Fp:
        for beta in Fp:
            score = score_plane(
                (Fp.one(), Fp.zero(), alpha),
                (Fp.zero(), Fp.one(), beta),
            )
            scanned += 1
            if score[2] == 4 and score[1] in (0, 2):
                winners.append({
                    "chart": "z2", "alpha": int(alpha), "beta": int(beta),
                    "square_degree": score[0], "odd_content_degree": score[1],
                })
    for alpha in Fp:
        score = score_plane(
            (Fp.one(), alpha, Fp.zero()),
            (Fp.zero(), Fp.zero(), Fp.one()),
        )
        scanned += 1
        if score[2] == 4 and score[1] in (0, 2):
            winners.append({
                "chart": "z1", "alpha": int(alpha),
                "square_degree": score[0], "odd_content_degree": score[1],
            })
    score = score_plane(
        (Fp.zero(), Fp.one(), Fp.zero()),
        (Fp.zero(), Fp.zero(), Fp.one()),
    )
    scanned += 1
    if score[2] == 4 and score[1] in (0, 2):
        winners.append({
            "chart": "z0", "square_degree": score[0],
            "odd_content_degree": score[1],
        })
    print(
        "FIXEDREVERSE4A1PHYSRR_PLANE_SEARCH|p={}|scanned={}|winners={}".format(
            p, scanned, json.dumps(winners, sort_keys=True)
        ),
        flush=True,
    )
    raise SystemExit(0)

# At the first old I2 support, P misses the node, so the chord denominator is
# a unit.  The quartic uses the opposite-point chord trivialization: its
# resolved restriction is AA+BB*Z*(y+Hy)/(x-Hx).  The plus sign is fixed by
# the exact three-plane selector (and distinguishes C1 from C2).
selected_component = 2 if "--component=2" in sys.argv else 1
if selected_component == 1:
    support = QQ(pointing["effective_horizontal_components"][0]["child_I2_support"])
else:
    remaining_factor_for_component = R(pointing["remaining_vertical_components"][
        "child_I2_support_factor_coefficients_low_to_high"
    ])
    support = QQ(-remaining_factor_for_component[0] / remaining_factor_for_component[1])
node_x = QQ(-3 * B(support) / (2 * A(support)))
assert node_x ** 3 + A(support) * node_x + B(support) == 0
assert 3 * node_x ** 2 + A(support) == 0
Hx = X / Z ** 2
Hy = Y / Z ** 3
assert node_x - Hx(support)
exceptional_row = matrix(QQ, [[
    AA(support) + BB(support) * Z(support) * Hy(support) / (node_x - Hx(support))
    for AA, BB in smooth_pairs
]])
rank_prime = 131
assert exceptional_row.rank() == 1
assert exceptional_row.change_ring(GF(rank_prime)).rank() == 1
resolved_coordinates = exceptional_row.right_kernel_matrix(
    algorithm="generic", basis="computed"
)
assert resolved_coordinates.dimensions() == (2, 3)
pairs = []
for row in resolved_coordinates.rows():
    AA = sum((row[index] * smooth_pairs[index][0] for index in range(3)), R.zero())
    BB = sum((row[index] * smooth_pairs[index][1] for index in range(3)), R.zero())
    assert (AA * X - BB * Y) % collision_modulus == 0
    assert AA(support) + BB(support) * Z(support) * Hy(support) / (node_x - Hx(support)) == 0
    pair_values = AA.list() + BB.list()
    denominator = ZZ.one()
    for value in pair_values:
        denominator = denominator.lcm(QQ(value).denominator())
    integral = [ZZ(QQ(value) * denominator) for value in pair_values]
    content = ZZ.zero()
    for value in integral:
        content = content.gcd(value)
    scale = QQ(denominator) / abs(content)
    pairs.append((R(scale * AA), R(scale * BB)))

pre_lll_bits = [rational_bits(AA.list() + BB.list()) for AA, BB in pairs]
pair_lattice = matrix(ZZ, [
    [ZZ(AA[degree]) for degree in range(aa_degree + 1)]
    + [ZZ(BB[degree]) for degree in range(bb_degree + 1)]
    for AA, BB in pairs
])
reduced_lattice = pair_lattice.LLL()
pairs = [
    (R(list(row[:aa_degree + 1])), R(list(row[aa_degree + 1:])))
    for row in reduced_lattice.rows()
]
for AA, BB in pairs:
    assert (AA * X - BB * Y) % collision_modulus == 0
    assert AA(support) + BB(support) * Z(support) * Hy(support) / (node_x - Hx(support)) == 0
print(
    "FIXEDREVERSE4A1PHYSRR_KERNEL|ambient=45|collision_rank=42|"
    "component_rank=1|h0=2|pre_lll_bits={}|bits={}|degrees={}".format(
        pre_lll_bits,
        [rational_bits(AA.list() + BB.list()) for AA, BB in pairs],
        [(AA.degree(), BB.degree()) for AA, BB in pairs],
    ),
    flush=True,
)

# Compile the binary quartic in the nested ring QQ[u][v].  This keeps the
# expensive coefficient arithmetic in FLINT univariate polynomials.
Pencil = PolynomialRing(R, "v")
v = Pencil.gen()
AA0, BB0 = pairs[0]
AA1, BB1 = pairs[1]
aa = Pencil([AA0, AA1])
bb = Pencil([BB0, BB1])
Xu, Yu, Zu, Au = [Pencil([poly]) for poly in (X, Y, Z, A)]
pair_fingerprint = hashlib.sha256(json.dumps([
    [coefficients(AA), coefficients(BB)] for AA, BB in pairs
], separators=(",", ":")).encode()).hexdigest()
RADICAND_CHECKPOINT = LOCAL / (
    "fixed-reverse-4a1-physical-after-collision-c{}-qq.json".format(selected_component)
)
checkpoint = read_json(RADICAND_CHECKPOINT) if RADICAND_CHECKPOINT.exists() else None
if (
    checkpoint
    and checkpoint.get("status") == "PASS_EXACT_QQ_FIXED_REVERSE_4A1_PHYSICAL_AFTER_COLLISION"
    and checkpoint.get("pair_fingerprint") == pair_fingerprint
):
    after_collision = Pencil([
        R(values) for values in checkpoint["coefficients_in_v_then_u_low_to_high"]
    ])
else:
    raw = (
        aa ** 4
        - 6 * Xu * aa ** 2 * bb ** 2
        + 8 * Yu * aa * bb ** 3
        - 3 * Xu ** 2 * bb ** 4
        - 4 * Au * bb ** 4 * Zu ** 4
    )
    after_collision, remainder = raw.quo_rem(Zu ** 4)
    assert not remainder
    RADICAND_CHECKPOINT.write_text(json.dumps({
        "schema": "elkies-k3.fixed-reverse-4a1-physical-after-collision-qq.v1",
        "status": "PASS_EXACT_QQ_FIXED_REVERSE_4A1_PHYSICAL_AFTER_COLLISION",
        "pair_fingerprint": pair_fingerprint,
        "coefficients_in_v_then_u_low_to_high": [
            coefficients(R(value)) for value in after_collision.list()
        ],
    }, indent=2, sort_keys=True) + "\n")

VU = PolynomialRing(QQ, names=("v", "u"))
vv, uu = VU.gens()
after = sum((
    QQ(coefficient) * vv ** v_degree * uu ** u_degree
    for v_degree, u_polynomial in enumerate(after_collision.list())
    for u_degree, coefficient in enumerate(R(u_polynomial).list())
), VU.zero())

# Extract the maximal square divisor by derivative-gcd layers.  This is exact
# and avoids factoring the large bivariate radicand.
derivative_layers = [after]
while True:
    current = derivative_layers[-1]
    repeated_layer = current.gcd(current.derivative(vv)).gcd(current.derivative(uu))
    if repeated_layer.is_constant():
        break
    derivative_layers.append(repeated_layer)
    print(
        "FIXEDREVERSE4A1PHYSRR_SQUARE_LAYER|index={}|degrees=({}, {})".format(
            len(derivative_layers) - 1,
            repeated_layer.degree(vv),
            repeated_layer.degree(uu),
        ),
        flush=True,
    )
square_factor = VU.one()
for index in range(1, len(derivative_layers), 2):
    denominator = (
        derivative_layers[index + 1]
        if index + 1 < len(derivative_layers)
        else VU.one()
    )
    quotient, remainder = derivative_layers[index].quo_rem(denominator)
    assert not remainder
    square_factor *= quotient
squarefree, remainder = after.quo_rem(square_factor ** 2)
assert not remainder

# A fixed odd vertical branch, if present, is the univariate content of the
# five pencil coefficients.  Strip it by gcd and retain the moving quartic.
by_v_degree = [R.zero() for unused in range(squarefree.degree(vv) + 1)]
for (v_degree, u_degree), coefficient in squarefree.dict().items():
    by_v_degree[v_degree] += QQ(coefficient) * u ** u_degree
fixed_branch = R.zero()
for value in by_v_degree:
    fixed_branch = fixed_branch.gcd(value)
fixed_branch = fixed_branch.monic()
fixed_branch_vu = sum((
    QQ(value) * uu ** degree
    for degree, value in enumerate(fixed_branch.list())
), VU.zero())
quartic, remainder = squarefree.quo_rem(fixed_branch_vu)
assert not remainder
print(
    "FIXEDREVERSE4A1PHYSRR_RESIDUAL|squarefree=({}, {})|fixed_branch={}|quartic=({}, {})".format(
        squarefree.degree(vv), squarefree.degree(uu), fixed_branch.degree(),
        quartic.degree(vv), quartic.degree(uu),
    ),
    flush=True,
)
assert quartic.degree(vv) == 4 and quartic.degree(uu) == 4
assert quartic.gcd(quartic.derivative(vv)).gcd(quartic.derivative(uu)).is_constant()

S = PolynomialRing(QQ, "v")
U_over_S = PolynomialRing(S, "u")
quartic_univariate = U_over_S(quartic)
square_factor_univariate = U_over_S(square_factor)
quartic_coefficients = [S(quartic_univariate[degree]) for degree in range(5)]
e, d, c, b, a = quartic_coefficients
I = S(12 * a * e - 3 * b * d + c ** 2)
J = S(72 * a * c * e + 9 * b * c * d - 27 * a * d ** 2 - 27 * b ** 2 * e - 2 * c ** 3)
A_child = S(-27 * I)
B_child = S(-27 * J)
removed_scalings = []
for factor, unused in A_child.gcd(B_child).factor():
    order = min(A_child.valuation(factor) // 4, B_child.valuation(factor) // 6)
    if order:
        A_child //= factor ** (4 * order)
        B_child //= factor ** (6 * order)
        removed_scalings.append((factor, int(order)))
assert A_child.degree() <= 8 and B_child.degree() <= 12
Delta = S(-16 * (4 * A_child ** 3 + 27 * B_child ** 2))
assert Delta.degree() <= 24

# Exact semistable fibre classification by gcd layers only.
repeated = Delta.gcd(Delta.derivative()).monic()
expected_i2_count = 4
print(
    "FIXEDREVERSE4A1PHYSRR_JACOBIAN_DIAGNOSTIC|degrees={}|repeated_degree={}|"
    "removed_scalings={}".format(
        (A_child.degree(), B_child.degree(), Delta.degree()),
        repeated.degree(),
        [(factor.degree(), order) for factor, order in removed_scalings],
    ),
    flush=True,
)
assert repeated.degree() == expected_i2_count
residual, remainder = Delta.quo_rem(repeated ** 2)
assert not remainder and residual.degree() == 16
assert repeated.is_squarefree()
assert residual.gcd(residual.derivative()).is_constant()
assert repeated.gcd(residual).is_constant()
assert repeated.gcd(A_child).is_constant() and repeated.gcd(B_child).is_constant()
assert residual.gcd(A_child).is_constant() and residual.gcd(B_child).is_constant()
infinity_orders = [
    int(8 - A_child.degree()),
    int(12 - B_child.degree()),
    int(24 - Delta.degree()),
]
assert infinity_orders == [0, 0, 0]

quartic_values = []
for coefficient in quartic_coefficients:
    quartic_values.extend(coefficient.list())
jacobian_values = A_child.list() + B_child.list() + Delta.list()
payload = {
    "schema": "elkies-k3.fixed-reverse-4a1-physical-resolved-rr-qq.v1",
    "status": "PASS_EXACT_QQ_FIXED_REVERSE_4A1_RR_JACOBIAN",
    "reproducing_command": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/compile_h92_fixed_reverse_4a1_physical_rr_qq.sage"
    ),
    "inputs": {
        name: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
        for name, path in (
            ("target", TARGET),
            ("physical_nef_audit", NEF_AUDIT),
            ("surface", SURFACE),
            ("pointing", POINTING),
        )
    },
    "divisor": {
        "class": "O+P-C1-9F",
        "P_dot_O": 21,
        "fibre_twist": -9,
        "subtracted_old_nonidentity_components": [1],
        "physical_fibre_in_3A1_coordinates": nef_audit[
            "physical_component_reduction"
        ]["reduced_fibre_in_3A1_coordinates"],
    },
    "smooth_RR": {
        "ambient_dimension": ambient_dimension,
        "collision_modulus_degree": int(collision_modulus.degree()),
        "collision_condition_rank": int(endpoint_rank),
        "post_collision_dimension": len(smooth_pairs),
        "exceptional_component_vanishing_rank": int(exceptional_row.rank()),
        "rank_witness_prime": rank_prime,
        "total_condition_rank": 43,
        "h0": len(pairs),
        "AA_degree_bound": aa_degree,
        "BB_degree_bound": bb_degree,
        "basis_pairs": [
            {
                "AA_coefficients_low_to_high": coefficients(AA),
                "BB_coefficients_low_to_high": coefficients(BB),
            }
            for AA, BB in pairs
        ],
        "maximum_basis_rational_bits": max(
            rational_bits(AA.list() + BB.list()) for AA, BB in pairs
        ),
    },
    "binary_quartic": {
        "coefficients_in_old_u_low_to_high": [
            coefficients(value) for value in quartic_coefficients
        ],
        "square_factor_coefficients_in_old_u_low_to_high": [
            coefficients(S(square_factor_univariate[degree]))
            for degree in range(square_factor_univariate.degree() + 1)
        ],
        "square_factor_total_degree": int(square_factor.total_degree()),
        "removed_fixed_branch_factor_coefficients_low_to_high": coefficients(fixed_branch),
        "maximum_rational_bits": rational_bits(quartic_values),
    },
    "child": {
        "minimal_A_coefficients_low_to_high": coefficients(A_child),
        "minimal_B_coefficients_low_to_high": coefficients(B_child),
        "discriminant_coefficients_low_to_high": coefficients(Delta),
        "degrees_A_B_Delta": [
            int(A_child.degree()), int(B_child.degree()), int(Delta.degree())
        ],
        "finite_fibres": [
            {
                "factor_coefficients_low_to_high": coefficients(repeated),
                "factor_degree": int(repeated.degree()),
                "orders_A_B_Delta": [0, 0, 2],
                "kodaira": "I2",
                "root_rank_contribution": int(repeated.degree()),
            },
            {
                "factor_coefficients_low_to_high": coefficients(residual),
                "factor_degree": int(residual.degree()),
                "orders_A_B_Delta": [0, 0, 1],
                "kodaira": "I1",
                "root_rank_contribution": 0,
            },
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
        "physical_affine_weyl_reduction": True,
        "smooth_collision_kernel": True,
        "single_exceptional_component_vanishing": True,
        "resolved_rank_sequence": [45, 3, 2],
        "multivariate_gcd_square_stripping": True,
        "full_discriminant_factorization": False,
        "groebner_or_surface_elimination": False,
        "runtime_seconds": time.monotonic() - started,
    },
    "proof_boundary": (
        "Exact QQ physical-nef h0=2, binary quartic, minimal semistable "
        "Jacobian and fibre classification. Exact child zero/component "
        "pointing is a separate gate."
    ),
}
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "FIXEDREVERSE4A1PHYSRR|ambient=45|collision_rank=42|component_rank=1|"
    "h0=2|quartic=4|degrees={}|I2=4|I1=16|bits={}|seconds={:.3f}|status={}|output={}".format(
        payload["child"]["degrees_A_B_Delta"],
        payload["child"]["maximum_A_B_Delta_rational_bits"],
        payload["method"]["runtime_seconds"],
        payload["status"],
        OUTPUT,
    ),
    flush=True,
)
