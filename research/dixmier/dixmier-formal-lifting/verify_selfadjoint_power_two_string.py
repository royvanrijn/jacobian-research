#!/usr/bin/env python3
"""Exact bounded controls for the written selfadjoint power-of-two theorem.

The universal theorem uses adjoint parity, the retained coprime theorem,
finite flow subtraction and valuation arguments. This checker does not
search all mates or prove those universal statements by finite testing.
"""

from __future__ import annotations

if not __debug__:
    raise RuntimeError("This verifier requires assertions; do not use -O or -OO.")

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform

import sympy as sp


x, h = sp.symbols("x h")


def clean(op):
    out = {}
    for degree, coefficient in op.items():
        coefficient = sp.cancel(sp.expand(coefficient))
        if coefficient != 0:
            out[degree] = coefficient
    return out


def add(*operators):
    out = {}
    for op in operators:
        for degree, coefficient in op.items():
            out[degree] = out.get(degree, 0) + coefficient
    return clean(out)


def scale(op, coefficient):
    return clean({degree: coefficient * value for degree, value in op.items()})


def dx(value):
    return sp.diff(value, x)


def multiply(left, right, deriv=dx):
    out = {}
    for i, a in left.items():
        for j, b in right.items():
            db = b
            for k in range(i + 1):
                degree = i + j - k
                out[degree] = out.get(degree, 0) + sp.binomial(i, k) * a * db
                db = deriv(db)
    return clean(out)


def power(op, n, deriv=dx):
    assert n >= 0
    out = {0: sp.Integer(1)}
    while n:
        if n & 1:
            out = multiply(out, op, deriv)
        op = multiply(op, op, deriv)
        n //= 2
    return out


def commutator(left, right, deriv=dx):
    return add(multiply(left, right, deriv), scale(multiply(right, left, deriv), -1))


def adjoint(op, deriv=dx):
    out = {}
    for i, coefficient in op.items():
        derivative = coefficient
        for k in range(i + 1):
            degree = i - k
            out[degree] = out.get(degree, 0) + (
                (-1) ** i * sp.binomial(i, k) * derivative
            )
            derivative = deriv(derivative)
    return clean(out)


def order(op):
    return max(op, default=-1)


def reduce_h(op):
    return clean({degree: value.subs(h, 0) for degree, value in op.items()})


def verify():
    counts = {
        "generic_quartic_identities": 0,
        "adjoint_parity_controls": 0,
        "power_two_positive_pairs": 0,
        "normalization_controls": 0,
        "specialization_controls": 0,
        "negative_boundary_controls": 0,
        "global_polynomial_landing_controls": 0,
    }
    one, ee, xx = {0: sp.Integer(1)}, {1: sp.Integer(1)}, {0: x}

    # Generic quartic adjoint criterion and first-order leading comparison.
    u2, u3, u4, aa = (sp.Function(name)(x) for name in ("u2", "u3", "u4", "a"))
    pp = {4: sp.Integer(1), 2: u2, 1: u3, 0: u4}
    defect = add(adjoint(pp), scale(pp, -1))
    assert defect == clean({
        1: 2 * (dx(u2) - u3),
        0: dx(dx(u2)) - dx(u3),
    })
    assert add(adjoint({4: 1, 2: u2, 1: dx(u2), 0: u4}),
               scale({4: 1, 2: u2, 1: dx(u2), 0: u4}, -1)) == {}
    tt = {1: aa, 0: dx(aa) / 2}
    assert adjoint(tt) == scale(tt, -1)
    assert commutator(pp, tt).get(4) == 4 * dx(aa)
    counts["generic_quartic_identities"] += 2

    # Exact antiinvolution and skewadjoint parity, with nonconstant coefficients.
    test = {4: x + 1, 3: x**2, 1: 1 / x, 0: x**3}
    assert adjoint(adjoint(test)) == clean(test)
    assert adjoint(multiply(test, pp)) == multiply(adjoint(pp), adjoint(test))
    for n in range(1, 11):
        qq = {n: x**2 + 1, 0: x**3}
        skew = scale(add(qq, scale(adjoint(qq), -1)), sp.Rational(1, 2))
        assert skew and order(skew) % 2 == 1
        assert adjoint(skew) == scale(skew, -1)
        counts["adjoint_parity_controls"] += 1
    assert add({0: x}, scale(adjoint({0: x}), -1)) == {}

    # Actual mates in every sampled power-of-two order; the theorem is not
    # inferred from this finite list.
    for n in (2, 4, 8, 16):
        pp = add(power(ee, n), scale(power(ee, 2), 3), xx)
        qq = add(power(pp, 2), scale(pp, 5), scale(ee, -1))
        assert order(pp) == n and order(qq) == 2 * n
        assert adjoint(pp) == pp
        assert commutator(pp, qq) == one
        skew = scale(add(qq, scale(adjoint(qq), -1)), sp.Rational(1, 2))
        assert skew == scale(ee, -1)
        assert add(qq, ee) == add(power(pp, 2), scale(pp, 5))
        assert commutator(pp, add(power(pp, 3), pp, one)) == {}
        counts["power_two_positive_pairs"] += 1

    # A nonmonic rational selfadjoint pair. E=rho D+delta(rho)/2
    # with rho=x^2; normalized derivative of w=-1/x is one.
    rho = x**2
    normalized = {1: rho, 0: dx(rho) / 2}
    assert adjoint(normalized) == scale(normalized, -1)
    ww = -1 / x
    assert sp.cancel(rho * dx(ww)) == 1
    assert dx(ww) != 1
    pp = add(power(normalized, 4), {0: ww})
    qq = scale(normalized, -1)
    assert adjoint(pp) == pp and pp[4] == x**8
    assert commutator(pp, qq) == one
    counts["normalization_controls"] += 1

    # The original derivation need not have the same time: delta=x d/dx,
    # E=x^{-1}D, and delta_E(x)=1 whereas delta(x)=x.
    delta = lambda value: x * sp.diff(value, x)
    normalized = {1: 1 / x}
    pp = add(power(normalized, 4, delta), xx)
    assert commutator(pp, scale(normalized, -1), delta) == one
    assert delta(x) == x and sp.cancel(delta(x) / x) == 1
    counts["normalization_controls"] += 1

    # Nonunit leading coefficient and a commutator that vanishes on reduction.
    pp = add(scale(power(ee, 4), h), power(ee, 2), xx)
    qq = add(power(pp, 2), scale(ee, -h))
    assert commutator(pp, qq) == {0: h}
    assert order(pp) == 4 and order(qq) == 8
    assert order(reduce_h(pp)) == 2 and order(reduce_h(qq)) == 4
    assert reduce_h(qq) == power(reduce_h(pp), 2)
    counts["specialization_controls"] += 1

    # A stable maximal ideal suffices even if its uniformizer is not constant.
    euler_h = lambda value: h * sp.diff(value, h)
    pp = {2: sp.Integer(1), 4: h}
    qq = add(scale(power(pp, 2, euler_h), 1 + h), {1: h})
    assert euler_h(h) == h
    assert order(reduce_h(pp)) == 2 and order(reduce_h(qq)) == 4
    counts["specialization_controls"] += 1

    # If the maximal ideal is not differential, reduction is not multiplicative.
    dh = lambda value: sp.diff(value, h)
    qq = power(pp, 2, dh)
    assert reduce_h(qq) == {4: sp.Integer(1), 5: sp.Integer(2)}
    assert order(reduce_h(qq)) == 5
    assert reduce_h(qq) != power(reduce_h(pp), 2, dh)
    counts["negative_boundary_controls"] += 1

    # The low-order-remainder hypothesis is also necessary for coefficient
    # integrality: a polar scalar can cancel a remainder of higher order.
    pp, qq = power(ee, 2), power(ee, 3)
    remainder = add(qq, scale(pp, -1 / h))
    assert add(scale(pp, 1 / h), remainder) == qq
    assert order(remainder) >= order(pp)
    counts["negative_boundary_controls"] += 1

    # A common-power leading symbol is not an actual operator power.
    ll = add(power(ee, 2), xx)
    pp = add(power(ll, 2), ee)
    assert commutator(pp, ll) == one
    assert pp.get(2) == 2 * x and pp.get(1) == 3
    assert pp.get(1) != dx(pp.get(2))
    assert add(pp, scale(adjoint(pp), -1)) == scale(ee, 2)
    assert add(pp, scale(power(ll, 2), -1)) == ee
    qq = add(ll, power(pp, 2))
    assert order(pp) == 4 and order(qq) == 8
    assert commutator(pp, qq) == one
    counts["negative_boundary_controls"] += 1

    # Selfadjointness alone does not imply the polynomial-centralizer claim.
    pp = power(ee, 4)
    assert adjoint(pp) == pp and commutator(pp, ee) == {}
    assert order(ee) % order(pp) != 0
    counts["negative_boundary_controls"] += 1

    # Squares outside the older infinity hypotheses still lie in the new
    # selfadjoint locus, and their nonconstant U violates the exact criterion.
    for potential, correction in ((x**2, x**10), (1 / x, x**3)):
        ll = add(power(ee, 2), {0: potential})
        pp = add(power(ll, 2), {0: correction})
        assert adjoint(pp) == pp
        assert pp.get(2) == 2 * potential and dx(pp.get(2)) != 0
        counts["negative_boundary_controls"] += 1

    # Actual polynomial landing after a non-polynomial-preservation claim
    # would be insufficient; here the subtraction gives explicit generators.
    normalized = add(ee, {0: x**2})
    pp = add(power(normalized, 4), xx)
    qq = add(power(pp, 2), scale(normalized, -1))
    q1 = add(qq, scale(power(pp, 2), -1))
    assert commutator(pp, qq) == one and order(q1) == 1
    recovered_x = add(pp, scale(power(q1, 4), -1))
    recovered_e = add(scale(q1, -1), scale(power(recovered_x, 2), -1))
    assert recovered_x == xx and recovered_e == ee
    counts["global_polynomial_landing_controls"] += 1
    return counts


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path,
                        default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo_root.resolve()
    directory = Path("research/dixmier/dixmier-formal-lifting")
    ledger_path = directory / "selfadjoint_power_two_string_source_ledger.json"
    ledger = json.loads((repo / ledger_path).read_text())
    for source in ledger["sources"]:
        assert digest(repo / source["path"]) == source["sha256"], source["path"]
    counts = verify()
    proof = directory / "SELFADJOINT_POWER_OF_TWO_STRING_CLASSIFICATION.md"
    checker = directory / Path(__file__).name
    receipt = {
        "schema": "selfadjoint_power_two_string_verification_v1",
        "status": "PASS",
        "arithmetic": "exact SymPy rational and symbolic differential expressions",
        "python": platform.python_version(),
        "sympy": sp.__version__,
        "counts": counts,
        "source_pins_checked": len(ledger["sources"]),
        "packet_pins": {str(path): digest(repo / path)
                        for path in (proof, checker, ledger_path)},
        "scope": "Bounded exact controls. The written proof, imported all-coprime "
                 "theorem, algebraic normalization/descent and differential-ideal "
                 "specialization argument establish the unrestricted claims.",
        "not_claimed": [
            "All even or all composite differential orders are classified",
            "Every quartic can be normalized to a selfadjoint quartic",
            "An original derivation has an algebraic primitive of one",
            "Selfadjointness alone forces centralizer C[P]",
            "A finite-control proof of the universal theorem",
            "External human review or formal verification",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "counts": counts,
                      "receipt": str(args.output)}, sort_keys=True))


if __name__ == "__main__":
    main()
