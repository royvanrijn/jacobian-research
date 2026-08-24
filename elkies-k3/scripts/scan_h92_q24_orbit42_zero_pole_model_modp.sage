#!/usr/bin/env sage -python
"""Scan degree-(4,6) sections of the cached exact q24 D12 model mod p.

The exact model cache is produced by
``recover_h92_q24_orbit42_zero_pole_smallprime.sage``.  This lightweight
stage permits good-prime selection without repeating the expensive exact
normalization.  It branches over the nonzero leading parameter and records
the full coefficient-Jacobian rank of every section.
"""

import argparse
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, Integers, PolynomialRing, matrix, vector


ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "artifacts/local/elkies-k3"
DEFAULT_MODEL = LOCAL / "q24-orbit42-zero-pole-seeds-mod-43.json"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--prime", type=int, required=True)
parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
parser.add_argument("--output", type=Path)
args = parser.parse_args()

p = ZZ(args.prime)
if not p.is_prime() or p in (2, 3):
    raise ValueError("prime must be an odd prime other than 3")
F = GF(p)

model_path = args.model.resolve()
model = json.loads(model_path.read_text())
if model.get("status") != "PASS_Q42_ZERO_POLE_SMALLPRIME_SEEDS":
    raise ValueError("model cache is not a passing small-prime seed artifact")
exact = model["exact_model"]


def red_q(value):
    value = QQ(value)
    denominator = ZZ(value.denominator())
    if denominator % p == 0:
        raise ZeroDivisionError(f"bad reduction denominator at p={p}")
    return F(ZZ(value.numerator())) / F(denominator)


R = PolynomialRing(F, "u")
u = R.gen()
A = R([red_q(v) for v in exact["A_coefficients_low_to_high"]])
B = R([red_q(v) for v in exact["B_coefficients_low_to_high"]])
if (A.degree(), B.degree()) != (6, 9):
    raise ArithmeticError("cached exact D12 model has wrong reduced degrees")

print(
    "Q42ZPSCAN|stage=MODEL|"
    f"prime={p}|Adeg=6|Bdeg=9|status=PASS",
    flush=True,
)


def solve_dx3():
    lead_ring = PolynomialRing(F, "c")
    c = lead_ring.gen()
    lead_equation = c**3 + F(A[6]) * c + F(B[9])
    leads = sorted({
        F(-factor[0] / factor[1])
        for factor, exponent in lead_equation.factor()
        if factor.degree() == 1
    }, key=int)
    answer = []
    for lead in leads:
        SR = PolynomialRing(
            F,
            names=("ell", "x2", "x1", "x0", "inv"),
            order="degrevlex",
        )
        ell, x2, x1, x0, inv = SR.gens()
        K = SR.fraction_field()
        U = PolynomialRing(K, "z")
        z = U.gen()
        AA = U([K(v) for v in A.list()])
        BB = U([K(v) for v in B.list()])
        x = K(lead) * z**3 + K(x2) * z**2 + K(x1) * z + K(x0)
        rhs = x**3 + AA * x + BB
        if rhs[9] or rhs.degree() > 8:
            raise ArithmeticError("dx3 leading equation failed")
        ys = {4: K(ell)}
        equations = [SR((K(ell)**2 - K(rhs[8])).numerator())]
        for degree in range(7, 3, -1):
            j = degree - 4
            known = sum(
                ys[i] * ys[degree - i]
                for i in ys
                if (degree - i) in ys and i != 4 and (degree - i) != 4
            )
            ys[j] = (K(rhs[degree]) - known) / (K(2) * ys[4])
        y = sum(ys[i] * z**i for i in range(5))
        residual = y**2 - rhs
        equations.extend(SR(K(residual[k]).numerator()) for k in range(4))
        equations.append(inv * ell - 1)
        solutions = SR.ideal(equations).variety()
        for solution in solutions:
            values = {g: F(solution[g]) for g in SR.gens()}
            xx = R([values[x0], values[x1], values[x2], lead])
            rhs_value = xx**3 + A * xx + B
            yy = [F.zero()] * 5
            yy[4] = values[ell]
            for degree in range(7, 3, -1):
                j = degree - 4
                known = sum(
                    yy[i] * yy[degree - i]
                    for i in range(5)
                    if 0 <= degree - i < 5 and i != 4 and (degree - i) != 4
                )
                yy[j] = (rhs_value[degree] - known) / (F(2) * yy[4])
            yy = R(yy)
            if yy**2 != rhs_value:
                raise ArithmeticError("modular dx3 section identity failed")
            answer.append((xx, yy))
    print(
        "Q42ZPSCAN|stage=DX3|"
        f"leads={len(leads)}|solutions={len(answer)}|status=PASS",
        flush=True,
    )
    return answer


def solve_dx4():
    answer = []
    for s_integer in range(1, int(p)):
        s_value = F(s_integer)
        SR = PolynomialRing(F, names=("x3", "x2", "x1", "x0"), order="degrevlex")
        x3, x2, x1, x0 = SR.gens()
        K = SR.fraction_field()
        U = PolynomialRing(K, "z")
        z = U.gen()
        AA = U([K(v) for v in A.list()])
        BB = U([K(v) for v in B.list()])
        x = (
            K(s_value**2) * z**4 + K(x3) * z**3 + K(x2) * z**2
            + K(x1) * z + K(x0)
        )
        rhs = x**3 + AA * x + BB
        ys = {6: K(s_value**3)}
        for degree in range(11, 5, -1):
            j = degree - 6
            known = sum(
                ys[i] * ys[degree - i]
                for i in ys
                if (degree - i) in ys and i != 6 and (degree - i) != 6
            )
            ys[j] = (K(rhs[degree]) - known) / (K(2) * ys[6])
        y = sum(ys[i] * z**i for i in range(7))
        residual = y**2 - rhs
        equations = [SR(K(residual[k]).numerator()) for k in range(6)]
        solutions = SR.ideal(equations).variety()
        for solution in solutions:
            values = {g: F(solution[g]) for g in SR.gens()}
            xx = R([
                values[x0], values[x1], values[x2], values[x3], s_value**2
            ])
            rhs_value = xx**3 + A * xx + B
            yy = [F.zero()] * 7
            yy[6] = s_value**3
            for degree in range(11, 5, -1):
                j = degree - 6
                known = sum(
                    yy[i] * yy[degree - i]
                    for i in range(7)
                    if 0 <= degree - i < 7 and i != 6 and (degree - i) != 6
                )
                yy[j] = (rhs_value[degree] - known) / (F(2) * yy[6])
            yy = R(yy)
            if yy**2 != rhs_value:
                raise ArithmeticError("modular dx4 section identity failed")
            answer.append((xx, yy))
    return answer


def padded(poly, length):
    return [poly[i] if i <= poly.degree() else F.zero() for i in range(length)]


def coefficient_jacobian(x_value, y_value):
    derivative_x = -3 * x_value**2 - A
    derivative_y = 2 * y_value
    columns = []
    for power in range(5):
        columns.append(padded(u**power * derivative_x, 13))
    for power in range(7):
        columns.append(padded(u**power * derivative_y, 13))
    return matrix(F, columns).transpose()


def degree_adapted_hensel_test(x_value, y_value, target_depth=4):
    dx = int(x_value.degree())
    dy = int(y_value.degree())
    equation_degree = max(2 * dy, 3 * dx, A.degree() + dx, B.degree())
    derivative_x = -3 * x_value**2 - A
    derivative_y = 2 * y_value
    columns = []
    for power in range(dx + 1):
        columns.append(padded(u**power * derivative_x, equation_degree + 1))
    for power in range(dy + 1):
        columns.append(padded(u**power * derivative_y, equation_degree + 1))
    jacobian = matrix(F, columns).transpose()
    unknowns = dx + dy + 2
    rank = int(jacobian.rank())
    if rank != unknowns:
        return rank, 1

    rows = tuple(map(int, jacobian.transpose().pivots()))
    square = jacobian.matrix_from_rows(rows)
    if not square.is_invertible():
        raise ArithmeticError("degree-adapted modular minor is singular")

    values = (
        [ZZ(int(v)) for v in x_value.list()]
        + [ZZ(int(v)) for v in y_value.list()]
    )
    achieved = 1
    for level in range(1, target_depth):
        modulus = p**(level + 1)
        coefficient_ring = Integers(modulus)
        polynomial_ring = PolynomialRing(coefficient_ring, "u")

        def red_mod(value):
            value = QQ(value)
            denominator = coefficient_ring(ZZ(value.denominator()))
            if not denominator.is_unit():
                raise ZeroDivisionError("bad prime-power model denominator")
            return coefficient_ring(ZZ(value.numerator())) / denominator

        coefficient_a = polynomial_ring([
            red_mod(v) for v in exact["A_coefficients_low_to_high"]
        ])
        coefficient_b = polynomial_ring([
            red_mod(v) for v in exact["B_coefficients_low_to_high"]
        ])
        x_lift = polynomial_ring(values[:dx + 1])
        y_lift = polynomial_ring(values[dx + 1:])
        residual = (
            y_lift**2 - x_lift**3 - coefficient_a * x_lift - coefficient_b
        )
        rhs = []
        for row in rows:
            coefficient = (
                ZZ(residual[row]) if row <= residual.degree() else ZZ.zero()
            )
            if coefficient % p**level:
                return rank, achieved
            rhs.append(F(-(coefficient // p**level)))
        correction = square.solve_right(vector(F, rhs))
        values = [
            value + p**level * ZZ(delta)
            for value, delta in zip(values, correction)
        ]
        x_check = polynomial_ring(values[:dx + 1])
        y_check = polynomial_ring(values[dx + 1:])
        check = (
            y_check**2 - x_check**3
            - coefficient_a * x_check - coefficient_b
        )
        if check:
            return rank, achieved
        achieved = level + 1
    return rank, achieved


sections = solve_dx3() + solve_dx4()
unique = {}
for x_value, y_value in sections:
    unique[(tuple(x_value.list()), tuple(y_value.list()))] = (x_value, y_value)
sections = list(unique.values())

records = []
for index, (x_value, y_value) in enumerate(sections):
    rank = int(coefficient_jacobian(x_value, y_value).rank())
    adapted_rank, hensel_depth = degree_adapted_hensel_test(x_value, y_value)
    records.append({
        "index": index,
        "x_coefficients_low_to_high": [int(v) for v in x_value.list()],
        "y_coefficients_low_to_high": [int(v) for v in y_value.list()],
        "coefficient_jacobian_rank": rank,
        "degree_adapted_jacobian_rank": adapted_rank,
        "hensel_depth": hensel_depth,
        "lifts_mod_p2": hensel_depth >= 2,
    })

payload = {
    "schema": "elkies-k3.h3-q24-orbit42-zero-pole-model-modp.v1",
    "status": "PASS_Q42_ZERO_POLE_MODEL_MODP_SCAN",
    "prime": int(p),
    "model_cache": str(model_path.relative_to(ROOT)),
    "section_count": len(records),
    "isolated_count": sum(r["coefficient_jacobian_rank"] == 12 for r in records),
    "degree_adapted_isolated_count": sum(
        r["degree_adapted_jacobian_rank"]
        == len(r["x_coefficients_low_to_high"])
        + len(r["y_coefficients_low_to_high"])
        for r in records
    ),
    "p2_liftable_count": sum(r["lifts_mod_p2"] for r in records),
    "depth4_liftable_count": sum(r["hensel_depth"] >= 4 for r in records),
    "rank_histogram": {
        str(rank): sum(r["coefficient_jacobian_rank"] == rank for r in records)
        for rank in sorted({r["coefficient_jacobian_rank"] for r in records})
    },
    "sections": records,
    "proof_boundary": (
        "Modular section enumeration and Jacobian ranks only. No lift, MW "
        "identification, orbit42 pencil, or A11 child is claimed."
    ),
}

OUTPUT = (
    args.output.resolve()
    if args.output
    else LOCAL / f"q24-orbit42-zero-pole-model-mod-{p}.json"
)
OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(f"OUTPUT|{OUTPUT}", flush=True)
print(
    "Q42ZPSCAN_RESULT|"
    f"prime={p}|sections={len(records)}|isolated={payload['isolated_count']}|"
    f"adapted={payload['degree_adapted_isolated_count']}|"
    f"p2={payload['p2_liftable_count']}|"
    f"depth4={payload['depth4_liftable_count']}|"
    "ranks=" + ",".join(
        f"{rank}:{count}" for rank, count in payload["rank_histogram"].items()
    ) + "|status=PASS_Q42_ZERO_POLE_MODEL_MODP_SCAN",
    flush=True,
)
