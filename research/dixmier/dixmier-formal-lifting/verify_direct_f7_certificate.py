"""Verify the F7 generators against direct bracket expansion and replay the exact certificate."""
from fractions import Fraction
from hashlib import sha256
from math import comb
from pathlib import Path
from zipfile import ZipFile

import flint


def read_pinned_helpers():
    archive = (
        Path(__file__).resolve().parents[2]
        / "jacobian"
        / "planar-degree-frontier"
        / "sources"
        / "niedbala-giraudin-2611.00896v1.zip"
    )
    with ZipFile(archive) as zipped:
        ancillary = zipped.read("anc/verify_F7F8.py")
    expected = "3e76cdb1062474cbd7101e964d8404e5ec48b12101d68e7cc8852dd7841a7a83"
    assert sha256(ancillary).hexdigest() == expected
    source = ancillary.decode()
    helpers = {"__name__": "pinned_f7_helpers"}
    exec(source.split("OK = []", 1)[0], helpers)
    definitions = source.index("def machine(")
    end = source.index("\nfor name, kmax, N, DEG, MOD, wt in ", definitions)
    exec(source[definitions:end], helpers)
    return helpers


def verify_fixed_data(system):
    machine = system["m"]
    field = machine.F

    def edge_coefficients(power):
        coefficients = [field.zero] * (2 * power + 1)
        H_power = field.one
        for j in range(power + 1):
            for i in range(power + 1):
                value = field.smul(comb(power, i) * comb(power, j), H_power)
                coefficients[i + j] = field.add(coefficients[i + j], value)
            H_power = field.mul(H_power, field.e)
        return coefficients

    assert system["Psupp"] == {
        (b, a) for a in range(7) for b in range(2*a + 1) if 7*a-5*b <= 2
    }
    assert system["Qsupp"] == {
        (b, a) for a in range(22) for b in range(2*a + 1) if 7*a-5*b <= 7
    }
    P_fixed = {(1, 1): field.one}
    P_fixed.update({(8+i, 6): c for i, c in enumerate(edge_coefficients(2))})
    Q_fixed = {(k, k): field.zero for k in (1, 2, 3)}
    Q_fixed.update({
        (7*j, 1+5*j): field.const(c)
        for j, c in enumerate((Fraction(35,128), Fraction(35,16), Fraction(35,8), Fraction(7,2)))
    })
    Q_fixed.update({(28+i, 21): c for i, c in enumerate(edge_coefficients(7))})
    assert system["Pfix"] == P_fixed
    assert system["Qfix"] == Q_fixed
    assert system["J"] == field.const(Fraction(-35, 128))


def direct_conditions(system):
    machine = system["m"]
    field = machine.F
    weight = system["W"]
    P_fixed, Q_fixed = system["Pfix"], system["Qfix"]
    P = {key: machine.c(coef) for key, coef in P_fixed.items()}
    Q = {key: machine.c(coef) for key, coef in Q_fixed.items()}

    context = flint.fmpq_mpoly_ctx.get(machine.names + ["s", "v"], "degrevlex")
    s_index, v_index = machine.nv, machine.nv + 1

    def expand(values):
        terms = {}
        for (b, a), value in values.items():
            for exponent, coefficient in value.to_dict().items():
                for r in range(b // 2 + 1):
                    shifted = list(exponent)
                    shifted[1] += r
                    key = tuple(shifted) + (b - r, a - r)
                    terms[key] = terms.get(key, 0) + coefficient * comb(b // 2, r)
        return context.from_dict({key: value for key, value in terms.items() if value})

    def bracket_at(left, right, target):
        bracket = left.derivative(s_index) * right.derivative(v_index)
        bracket -= left.derivative(v_index) * right.derivative(s_index)
        coefficients = {}
        for monomial, coefficient in bracket.to_dict().items():
            key = monomial[-2:]
            if weight(*key) == target:
                coefficients.setdefault(key, {})[monomial[:-2]] = coefficient
        reduced = {
            key: machine.red(machine.ctx.from_dict(poly))
            for key, poly in coefficients.items()
        }
        return {key: value for key, value in reduced.items() if value}

    def as_field(polynomial):
        coefficients = [Fraction(0)] * field.d
        for exponents, number in polynomial.to_dict().items():
            assert not any(exponents[1:])
            coefficients[exponents[0]] = Fraction(int(number.p), int(number.q))
        return tuple(coefficients)

    P_top = expand({key: c for key, c in P.items() if weight(*key) == 2})
    Q_top = expand({key: c for key, c in Q.items() if weight(*key) == 7})
    assert bracket_at(P_top, Q_top, 7) == {(0, 1): machine.c(Fraction(35, 128))}

    conditions = []
    for defect in range(1, 31):
        target = 7 - defect
        known = bracket_at(expand(P), expand(Q), target)
        columns = []
        for name, poly, fixed, support, top_weight in (
            ("P", P, P_fixed, system["Psupp"], 2),
            ("Q", Q, Q_fixed, system["Qsupp"], 7),
        ):
            for key in sorted(support):
                if weight(*key) != top_weight - defect or key in fixed or key == (0, 0):
                    continue
                basis_element = expand({key: machine.c(1)})
                if name == "P":
                    column = bracket_at(basis_element, Q_top, target)
                else:
                    column = bracket_at(P_top, basis_element, target)
                columns.append((poly, name, key, column))

        rows = sorted(set(known).union(*(column.keys() for _, _, _, column in columns)))
        matrix = [
            [as_field(column.get(key, machine.zero)) for _, _, _, column in columns]
            for key in rows
        ]
        rhs = [-known.get(key, machine.zero) for key in rows]
        pivots = []
        for j in range(len(columns)):
            i = len(pivots)
            pivot = next((k for k in range(i, len(rows)) if not field.iszero(matrix[k][j])), None)
            if pivot is None:
                continue
            matrix[i], matrix[pivot] = matrix[pivot], matrix[i]
            rhs[i], rhs[pivot] = rhs[pivot], rhs[i]
            factor = field.inv(matrix[i][j])
            matrix[i] = [field.mul(factor, v) for v in matrix[i]]
            rhs[i] = machine.smul(factor, rhs[i])
            for k in range(len(rows)):
                factor = matrix[k][j]
                if k == i or field.iszero(factor):
                    continue
                matrix[k] = [
                    field.sub(x, field.mul(factor, y))
                    for x, y in zip(matrix[k], matrix[i])
                ]
                rhs[k] -= machine.smul(factor, rhs[i])
            pivots.append(j)

        free = [j for j in range(len(columns)) if j not in pivots]
        values = {j: machine.newpar(columns[j][1:3]) for j in free}
        for i, j in enumerate(pivots):
            value = rhs[i]
            for k in free:
                value -= machine.smul(matrix[i][k], values[k])
            values[j] = value
        for j, (poly, _, key, _) in enumerate(columns):
            poly[key] = values[j]
        conditions.extend(value for value in rhs[len(pivots):] if value)
    return conditions


def main():
    helpers = read_pinned_helpers()
    pinned = helpers["run_machine"](helpers["machine"]("F7", 30), 30)
    system = helpers["machine"]("F7", 30)
    verify_fixed_data(system)
    conditions = direct_conditions(system)
    assert system["m"].npar == pinned.npar == 5
    assert len(conditions) == len(pinned.conds) == 17
    remaining = list(pinned.conds)
    for condition in conditions:
        remaining.remove(condition)
    assert not remaining
    print("All 17 F7 compatibility generators agree with direct bracket expansion.", flush=True)
    ok, details = helpers["certificate"](pinned, 20, 14, 5, [2, 2, 4, 1, 2, 3])
    assert ok
    print(f"Exact F7 unit-ideal identity verified ({details}).", flush=True)


if __name__ == "__main__":
    main()
