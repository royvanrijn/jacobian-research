#!/usr/bin/env python3
"""Audit the rank-one fixed-quintic common-resolvent elliptic slice.

On the slice

    kappa = (5/4) A,  R = 4,

the cube condition is

    y^2 = x^3 + 30464/15,  x = 4*Pi,  y = 12*A.

Every rational point fails the Hasse local condition over Q_5.  The exact
proof is a valuation calculation followed by the Newton polygon of the
cubic after translating T=2A+Y.  The bounded Mordell--Weil enumeration is
only a regression for that all-points proof.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess

import sympy as sp


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--multiple-bound",
        type=int,
        default=12,
        help="enumerate nonzero multiples n*P with |n| at most this bound",
    )
    parser.add_argument("--padic-precision", type=int, default=30)
    parser.add_argument("--show-failures", action="store_true")
    return parser.parse_args()


def verify_symbolic_reconstruction() -> None:
    A, Pi, T, Y = sp.symbols("A Pi T Y")
    V = (A**2 + 48) / 4
    W = Pi**3
    d = A**2 - V - 5
    e = (4 * W - V * d) / A
    q = T**2 + A * T + V
    h = T**3 - A * T**2 + d * T + e
    product = sp.Poly(sp.expand(q * h), T)

    assert product.coeff_monomial(T**4) == 0
    assert product.coeff_monomial(T**3) == -5
    assert sp.factor(product.coeff_monomial(T) - 4 * Pi**3) == 0
    assert sp.factor(sp.discriminant(q, T) + 48) == 0

    relation = {Pi**3: (135 * A**2 - 1904) / 60}
    discriminant_quotient = sp.factor(
        sp.discriminant(h, T) / sp.discriminant(q, T)
    )
    expected_square = (
        (5 * A**2 - 204) * (15 * A**2 - 272) / (960 * A)
    ) ** 2
    assert sp.factor(
        discriminant_quotient.subs(relation) - expected_square
    ) == 0

    translated = sp.Poly(
        sp.factor(sp.expand(h.subs(T, 2 * A + Y)).subs(relation)),
        Y,
    )
    assert translated.coeff_monomial(Y**3) == 1
    assert translated.coeff_monomial(Y**2) == 5 * A
    assert translated.coeff_monomial(Y) == 35 * A**2 / 4 - 17
    assert sp.factor(
        translated.coeff_monomial(1)
        - (1275 * A**4 - 7140 * A**2 + 18496) / (240 * A)
    ) == 0
    print("PASS: exact factor and translated-cubic reconstruction")


def verify_five_adic_obstruction() -> None:
    # From 135*A^2=60*Pi^3+1904, v_5(Pi)>=0 would give
    # 1+2*v_5(A)=0, which is impossible.  For negative valuation,
    # comparison of the unique lowest terms gives
    # v_5(Pi)=-2m and v_5(A)=-3m with m>=1.
    assert sp.factorint(135)[5] == 1
    assert sp.factorint(60)[5] == 1
    assert 1904 % 5 == 4

    m = sp.symbols("m", integer=True, positive=True)
    coefficient_valuations = (
        1 - 9 * m,
        1 - 6 * m,
        1 - 3 * m,
        0,
    )
    endpoint_slope = 3 * m - sp.Rational(1, 3)
    endpoint_line = tuple(
        coefficient_valuations[0] + index * endpoint_slope
        for index in range(4)
    )
    assert sp.factor(
        coefficient_valuations[1] - endpoint_line[1]
    ) == sp.Rational(1, 3)
    assert sp.factor(
        coefficient_valuations[2] - endpoint_line[2]
    ) == sp.Rational(2, 3)
    assert sp.fraction(sp.together(endpoint_slope))[1] == 3

    # The quadratic discriminant -48 is the nonsquare unit 2 modulo 5.
    assert (-48) % 5 == 2
    assert {value**2 % 5 for value in range(5)} == {0, 1, 4}
    print(
        "PASS: every rational point on the rank-one slice fails over Q_5 "
        "(quadratic nonsplit; translated cubic slope denominator 3)"
    )


def main() -> None:
    args = parse_args()
    if args.multiple_bound <= 0:
        raise SystemExit("multiple-bound must be positive")
    if args.padic_precision <= 0:
        raise SystemExit("padic-precision must be positive")

    verify_symbolic_reconstruction()
    verify_five_adic_obstruction()

    gp = shutil.which("gp")
    if gp is None:
        raise SystemExit("PARI/GP executable 'gp' is required")

    program = rf"""
E=ellinit([0,0,0,0,30464/15]);
P=[-687701486576/167997515625,3050290856818742032/68857981716796875];
if(ellisoncurve(E,P)==0,error("bad rank-one generator"));
ntested=0;
nirred=0;
nlocalfail=0;
check(n,Q)={{
  my(Pi,A,V,d,e,q,h,qq,hh,discquot);
  if(#Q==1,return());
  Pi=Q[1]/4; A=Q[2]/12;
  if(135*A^2-60*Pi^3-1904,error("elliptic reconstruction failed"));
  V=(A^2+48)/4;
  d=A^2-V-5;
  e=(4*Pi^3-V*d)/A;
  q=x^2+A*x+V;
  h=x^3-A*x^2+d*x+e;
  ntested++;
  if(!polisirreducible(q)||!polisirreducible(h),return());
  nirred++;
  discquot=poldisc(h)/poldisc(q);
  if(!issquare(discquot),error("common resolvent identity failed"));
  qq=q/content(q);
  hh=h/content(h);
  if(#polrootspadic(qq,5,{args.padic_precision})
     ||#polrootspadic(hh,5,{args.padic_precision}),
    error("the exact five-adic obstruction failed"));
  nlocalfail++;
  if({1 if args.show_failures else 0},
    print("LOCAL_FAIL n,p,A,Pi=",[n,5,A,Pi]))
}};
for(n=-{args.multiple_bound},{args.multiple_bound},if(n,Q=ellmul(E,P,n);check(n,Q)));
print("SEARCH_COMPLETE tested=",ntested," irreducible=",nirred," hasse=0 local_fail=",nlocalfail);
"""
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if "***" in completed.stdout or "***" in completed.stderr:
        raise SystemExit(completed.stdout + completed.stderr)
    print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="")


if __name__ == "__main__":
    main()
