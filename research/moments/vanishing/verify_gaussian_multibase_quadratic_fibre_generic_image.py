"""Bounded exact controls for the multibase quadratic-fibre generic image.

The all-degree proof is in the companion note. These finite controls compare
original Wick contractions with normalized selected-field marker formulas.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import hashlib
import json
from math import factorial
from pathlib import Path
import platform
import sympy as sp

z1, z2, w1, w2, y, t = sp.symbols("z1 z2 w1 w2 y t")
zs = (z1, z2)
variables = (z1, z2, w1, w2, y)


def zero(expr):
    assert sp.cancel(sp.expand(expr)) == 0, expr


def trunc(expr, order):
    return sp.expand(sp.series(expr, t, 0, order+1).removeO())


def wick(expr):
    answer = 0
    for (a, b, c, d, e), coefficient in sp.Poly(sp.expand(expr), *variables).terms():
        if a == c and b == d and e % 2 == 0:
            # Ordinary covariance 2.
            answer += coefficient * factorial(a) * factorial(b) * sp.Rational(
                factorial(e), factorial(e//2))
    return sp.expand(answer)


def raw_periods(phase, markers, order):
    values = [sp.Integer(0) for _ in markers]
    power = sp.Integer(1)
    for n in range(order+1):
        if n:
            power = sp.expand(power*phase)
        for j, marker in enumerate(markers):
            values[j] += wick(marker*power)*t**n/factorial(n)
    return [trunc(value, order) for value in values]


h = sp.Matrix([1+z1**2, 1+z1])
J = sp.eye(2)-t*h.jacobian(zs)
Jinv = J.inv()
detJ = J.det()
D = 1-4*t*(1+z2)
S = 2/D
u = 1+z1
kk = z2**2
Rhat = t*kk+t**2*u**2/D
phase = w1*h[0]+w2*h[1]+(1+z2)*y**2+u*y+kk


def X(j, value):
    return sum(Jinv[i, j]*sp.diff(value, zs[i]) for i in range(2))


def B(j, value):
    correction = X(j, Rhat)-X(j, detJ)/detJ-X(j, D)/(2*D)
    return sp.cancel(X(j, value)+correction*value)


@lru_cache(None)
def normalized_monomial(a, b, c, d, e):
    mean = t*S*u
    moment = sum(
        sp.Rational(factorial(e), 2**j*factorial(j)*factorial(e-2*j))
        * S**j * mean**(e-2*j) for j in range(e//2+1))
    value = sp.cancel(z1**a*z2**b*moment)
    for _ in range(d):
        value = B(1, value)
    for _ in range(c):
        value = B(0, value)
    return value


def normalized(marker):
    return sp.cancel(sum(
        coefficient*normalized_monomial(*powers)
        for powers, coefficient in sp.Poly(marker, *variables).terms()))


def field_reduce(expr):
    # zeta_2=t(1+zeta_1), zeta_1=t(1+zeta_1^2).
    expression = sp.cancel(expr.subs(z2, t*(1+z1)))
    numerator, denominator = sp.fraction(expression)
    q = sp.Poly(z1-t*(1+z1**2), z1, domain=sp.QQ.frac_field(t))
    n = sp.Poly(numerator, z1, domain=q.domain).rem(q)
    den = sp.Poly(denominator, z1, domain=q.domain).rem(q)
    assert not den.is_zero
    return sp.cancel((n*sp.invert(den, q)).rem(q).as_expr())


def controls():
    order = 3
    root = sp.Matrix([0, 0])
    for _ in range(order+1):
        root = sp.Matrix([trunc(t*value.subs({z1:root[0], z2:root[1]}, simultaneous=True), order)
                          for value in h])
    selected = {z1:root[0], z2:root[1]}
    zero(trunc(root[0]-t*(1+root[0]**2), order))
    zero(trunc(root[1]-t*(1+root[0]), order))
    f = trunc(sp.exp(trunc(Rhat.subs(selected, simultaneous=True), order))
              /(detJ.subs(selected, simultaneous=True)
                *sp.sqrt(D.subs(selected, simultaneous=True))), order)
    markers = [sp.Integer(1), z1, z2, w1, w2, y, y**2, w1*y, z1*w2]
    rows = raw_periods(phase, markers, order)
    for marker, row in zip(markers, rows):
        ratio = trunc(normalized(marker).subs(selected, simultaneous=True), order)
        zero(row-trunc(f*ratio, order))

    # Exact selected-field relation and generic rank-two monomial frame.
    zero(field_reduce(normalized(z2)-t*(1+normalized(z1))))
    q = sp.Poly(z1-t*(1+z1**2), z1, domain=sp.QQ.frac_field(t))
    assert q.degree() == 2 and q.is_irreducible
    zero(field_reduce(normalized(z1))-z1)
    zero(field_reduce(normalized(sp.Integer(1)))-1)

    # Original Gaussian Stein identities transported through the
    # normalized source: these equalities are exact in the selected field.
    amplitude = sp.Integer(1)
    gaussian_duals = [w1, w2, z1, z2, y/2]
    for coordinate, dual in zip(variables, gaussian_duals):
        stein = sp.diff(amplitude, coordinate)-dual*amplitude+t*sp.diff(phase, coordinate)*amplitude
        zero(field_reduce(normalized(stein)))

    velocity = Jinv*h
    delta = lambda value: sp.diff(value, t)+sum(velocity[i]*sp.diff(value, zs[i]) for i in range(2))
    tau = delta(Rhat)-delta(detJ)/detJ-delta(D)/(2*D)
    zero(field_reduce(normalized(phase)-tau))
    zero(field_reduce(B(0, B(1, sp.Integer(1)))-B(1, B(0, sp.Integer(1)))))

    # Full root scheme can have a positive-dimensional component.
    hx = sp.Matrix([z1*z2, z2**2])
    Jx = sp.eye(2)-t*hx.jacobian(zs)
    assert Jx.subs({z1:0,z2:0}) == sp.eye(2)
    for relation in sp.Matrix(zs)-t*hx:
        zero(relation.subs(z2, 1/t))
    other_det = sp.cancel(Jx.det().subs(z2, 1/t))
    assert other_det == 0

    # Nonscalar C: order D^-1*C is essential, S is symmetric.
    Cmat = sp.Matrix([[2, 1], [1, 3]])
    Amat = sp.Matrix([[z1, 1+z2], [1+z2, 2*z1]])
    Dmat = sp.eye(2)-2*t*Cmat*Amat
    Smat = Dmat.inv()*Cmat
    for entry in Smat-Smat.T:
        zero(entry)
    for entry in (Cmat.inv()-2*t*Amat)*Smat-sp.eye(2):
        zero(entry)
    wrong = Cmat*Dmat.inv()
    assert any(sp.cancel(entry) != 0 for entry in wrong-Smat)

    return {
        "status":"passed",
        "scope":"Bounded exact controls only; unrestricted selected-field theorem is established by the companion proof.",
        "two_base_phase":str(phase),
        "ordinary_covariance":2,
        "selected_root_through_t3":[str(value) for value in root],
        "selected_generic_rank":2,
        "original_frame":["Phi[1]","Phi[Z1]"],
        "original_Wick_markers":[str(marker) for marker in markers],
        "original_Wick_orders":[0,order],
        "exact_Stein_coordinates":[str(v) for v in variables],
        "exact_connection_identity":True,
        "variable_A_B_marking_commutator":True,
        "extraneous_positive_dimensional_component":"h=(z1*z2,z2^2): selected0; other line z2=1/t",
        "nonscalar_covariance_order_control":"S=(I-2tCA)^-1 C is symmetric; C(I-2tCA)^-1 differs",
        "independent_replay":False,"formal_verification":False,"external_review":False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = controls()
    directory = Path(__file__).resolve().parent
    names = [
        "GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_GENERIC_IMAGE.md",
        Path(__file__).name,
        "GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md",
        "GAUSSIAN_W_AFFINE_C_ZERO_GENERIC_IMAGE.md",
    ]
    report.update(python=platform.python_version(), sympy=sp.__version__,
                  sha256={name:hashlib.sha256((directory/name).read_bytes()).hexdigest()
                          for name in names})
    text=json.dumps(report,indent=2)+"\n"
    if args.output is not None:
        args.output.write_text(text,encoding="utf-8",newline="\n")
    print(text,end="")


if __name__ == "__main__":
    main()
