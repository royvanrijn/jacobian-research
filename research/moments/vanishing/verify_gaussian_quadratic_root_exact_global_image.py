"""Bounded exact controls for P=W(1+Z^2), its global image and literal kernel.

The all-degree proof is in the companion note. Finite controls compare the
original Wick map, algebraic source coordinates, regularizations, Bezout
identities and tangent normal forms. No unrestricted inference uses a prefix.
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

z, w, t = sp.symbols("z w t")
Delta = 1-4*t**2
phase = w*(1+z**2)


def zero(expr):
    assert sp.cancel(sp.expand(expr)) == 0, expr


def falling(x, n):
    return sp.prod(x-j for j in range(n))


@lru_cache(None)
def coordinates(a, b):
    prefactor = sp.Rational(1, 2**a)*(-4)**b*t**(b-a)
    even = sum(sp.binomial(a, 2*k)*falling(sp.Rational(2*k-1, 2), b)
               *Delta**(k-b) for k in range(a//2+1))
    odd = sum(sp.binomial(a, 2*k+1)*falling(k, b)*Delta**(k-b)
              for k in range(b, (a-1)//2+1))
    return sp.cancel(prefactor*even), sp.expand(-prefactor*odd)


def marker_coordinates(marker):
    pair = [sp.Integer(0), sp.Integer(0)]
    for (a, b), coefficient in sp.Poly(marker, z, w).terms():
        aa, bb = coordinates(a, b)
        pair[0] += coefficient*aa
        pair[1] += coefficient*bb
    return tuple(sp.cancel(value) for value in pair)


def source_coefficient(a, b, n):
    degree = n+a-b
    if degree < 0 or degree % 2:
        return sp.Integer(0)
    q = degree//2
    return sp.cancel(sp.Rational((-4)**b, 2**a)*(-4)**q/factorial(q)*sum(
        (-1)**j*sp.binomial(a, j)*falling(sp.Rational(j-1, 2), b)
        *falling(sp.Rational(j-1, 2)-b, q) for j in range(a+1)))


def wick(expr):
    total = 0
    for (a, b), coefficient in sp.Poly(sp.expand(expr), z, w).terms():
        if a == b:
            total += coefficient*factorial(a)
    return sp.expand(total)


def delta(marker):
    return sp.expand(w*(1-z**2)*marker-(1+z**2)*sp.diff(marker, z)
                     +2*z*w*sp.diff(marker, w))


def normal_form(marker):
    value = sp.expand(marker)
    certificate = sp.Integer(0)
    while True:
        divisible = [(a+b, a, b, c) for (a, b), c in sp.Poly(value, z, w).terms()
                     if a >= 2 and b >= 1]
        if not divisible:
            break
        _, a, b, coefficient = max(divisible, key=lambda row: row[:3])
        quotient = coefficient*z**(a-2)*w**(b-1)
        value = sp.expand(value+delta(quotient))
        certificate += quotient
    zero(value-marker-delta(certificate))
    return value, sp.expand(certificate)


def f_coefficient(n):
    return sp.binomial(n, n//2) if n >= 0 and n % 2 == 0 else sp.Integer(0)


def q_coordinates(j):
    if j == 0:
        return sp.Integer(1), sp.Integer(0)
    return t**(-j), -sum(f_coefficient(n)*t**(n-j) for n in range(j))


def is_laurent_at_origin(expr):
    _, denominator = sp.fraction(sp.cancel(expr))
    return len(sp.Poly(denominator, t).terms()) == 1


def checks():
    # Original Wick moments and exact monomial source coefficients.
    order, max_a, max_b = 12, 6, 4
    powers = [sp.expand(phase**n) for n in range(order+1)]
    monomials = 0
    for a in range(max_a+1):
        for b in range(max_b+1):
            aa, bb = coordinates(a, b)
            assert is_laurent_at_origin(bb)
            for term in sp.Add.make_args(sp.expand(bb)):
                if term:
                    assert term.as_powers_dict().get(t, 0) <= -b-1
            for n in range(min(0, b-a), 0):
                zero(source_coefficient(a, b, n))
            for n, power in enumerate(powers):
                raw = wick(z**a*w**b*power)/factorial(n)
                zero(raw-source_coefficient(a, b, n))
                lower = sp.Rational(n+b-a, 2)
                direct = (sp.Rational(factorial(n+b), factorial(n))
                          *sp.binomial(n, int(lower))
                          if lower.is_integer and 0 <= lower <= n else sp.Integer(0))
                zero(raw-direct)
            monomials += 1

    # Exact triangular regularizations; every q_j is in literal V.
    regularized = []
    for a in range(1, 11):
        H = sp.expand(sum(sp.binomial(a, 2*k)*Delta**k
                          for k in range(a//2+1))/2**a)
        assert H.subs(t, 0) == sp.Rational(1, 2)
        target = [sp.Integer(0), sp.Integer(0)]
        for k in range(a//2+1):
            coefficient = H.coeff(t, 2*k)
            qa, qb = q_coordinates(a-2*k)
            target[0] += coefficient*qa
            target[1] += coefficient*qb
        actual = coordinates(a, 0)
        zero(actual[0]-target[0])
        zero(actual[1]-target[1])
        # Rational projection has leading pole exactly -1/(2t^a).
        zero(sp.expand(actual[1]).coeff(t, -a)+sp.Rational(1, 2))
        regularized.append(a)

    # Removing t^b from the W^b period uses a polynomial Bezout identity.
    bezout_orders = []
    for b in range(1, 7):
        u, v, gcd = sp.gcdex(t**b, Delta**b, t)
        zero(gcd-1)
        zero(u*t**b+v*Delta**b-1)
        constant = 2**b*sp.prod(2*j-1 for j in range(1, b+1))
        aa, bb = coordinates(0, b)
        zero(aa-constant*t**b/Delta**b)
        zero(bb)
        zero(u*aa/constant+v-Delta**(-b))
        bezout_orders.append(b)
    za, zb = coordinates(1, 0)
    zero(1-2*t*za)
    zero(-2*t*zb-1)

    # Representative reverse lattice decomposition with both finite and
    # origin denominators; polynomial plus q_j plus Delta^-2*f.
    a0 = (1+3*t)/Delta**2
    polynomial = 1+t**3
    lambdas = {1:sp.Integer(2), 4:sp.Integer(-3)}
    aa = a0+sum(c*t**(-j) for j, c in lambdas.items())
    bb = polynomial+sum(c*q_coordinates(j)[1] for j, c in lambdas.items())
    target_a = a0+sum(c*q_coordinates(j)[0] for j, c in lambdas.items())
    target_b = polynomial+sum(c*q_coordinates(j)[1] for j, c in lambdas.items())
    zero(aa-target_a)
    zero(bb-target_b)
    reverse_sample = str(sp.factor(aa))

    # Degree-preserving exact tangent reduction and fixed-kernel controls.
    reduction_count = 0
    for degree in range(8):
        for a in range(degree+1):
            b = degree-a
            marker = z**a*w**b
            nf, certificate = normal_form(marker)
            assert sp.Poly(nf, z, w).total_degree() <= degree
            assert all(i < 2 or j < 1 for (i, j), _ in sp.Poly(nf, z, w).terms())
            left, right = marker_coordinates(marker), marker_coordinates(nf)
            zero(left[0]-right[0])
            zero(left[1]-right[1])
            reduction_count += 1
    for degree in range(5):
        for a in range(degree+1):
            kernel_pair = marker_coordinates(delta(z**a*w**(degree-a)))
            zero(kernel_pair[0])
            zero(kernel_pair[1])
    filtered = []
    for D in range(1, 8):
        basis = ([z**a for a in range(D+1)]
                 +[w**b for b in range(1, D+1)]
                 +[z*w**b for b in range(1, D)])
        assert len(basis) == 3*D
        filtered.append(dict(degree=D, normal_form_basis_count=len(basis)))
    for b in range(1, 8):
        aa, bb = coordinates(1, b)
        wa, wb = coordinates(0, b)
        zero(aa-wa/(2*t))
        zero(bb)
        zero(wb)
        evaluation = sp.Matrix([[sp.Rational(1, 2)**b, sp.Rational(1, 2)**(b-1)],
                                [sp.Rational(-1, 2)**b, sp.Rational(-1, 2)**(b-1)]])
        zero(evaluation.det()-sp.Rational(-1, 4)**(b-1))

    # The naive ring span contains f*f=Delta^-1, which is not in M because
    # the exact rational intersection is K[t], not K[t,Delta^-1].
    assert not is_laurent_at_origin(1/Delta)
    assert is_laurent_at_origin(1+t+t**3)
    factors = sp.factor_list(Delta, t)[1]
    assert len(factors) == 2 and all(exponent == 1 for _, exponent in factors)

    return dict(
        status="passed",
        original_Wick_control=dict(a_range=[0,max_a], b_range=[0,max_b],
                                   time_orders=[0,order], marker_count=monomials),
        exact_negative_Laurent_projection_bounds=True,
        exact_triangular_q_indices=regularized,
        exact_Bezout_orders=bezout_orders,
        exact_time_amplitude="Phi[1-2tZ]=1",
        reverse_lattice_sample=dict(f_coordinate=reverse_sample,
                                   description="(1+3t)Delta^-2*f plus 2q1-3q4 and 1+t^3"),
        degree_preserving_monomial_reductions=reduction_count,
        exact_tangent_kernel_controls="delta applied to every monomial of degree<=4",
        filtered_basis_controls=filtered,
        exact_two_pole_independence_orders=[1,7],
        nonrational_f_control="Delta has simple zeros at ±1/2",
        rejected_naive_ring_element="f^2=Delta^-1; rational nonpolynomial, outside M",
        scope=("Bounded exact controls only. The all-degree source, kernel independence "
               "and global lattice equality are proved in the companion note. "
               "Two active Witt variables only; no inert ordinary-marker kernel is asserted."),
        independent_replay=False, formal_verification=False, external_review=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = checks()
    directory = Path(__file__).resolve().parent
    names = ["GAUSSIAN_QUADRATIC_ROOT_EXACT_GLOBAL_IMAGE.md", Path(__file__).name,
             "GAUSSIAN_W_AFFINE_C_ZERO_GENERIC_IMAGE.md",
             "GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_GENERIC_IMAGE.md"]
    report.update(python=platform.python_version(), sympy=sp.__version__,
                  sha256={name:hashlib.sha256((directory/name).read_bytes()).hexdigest()
                          for name in names})
    text=json.dumps(report,indent=2)+"\n"
    if args.output is not None:
        args.output.write_text(text,encoding="utf-8",newline="\n")
    print(text,end="")


if __name__ == "__main__":
    main()
