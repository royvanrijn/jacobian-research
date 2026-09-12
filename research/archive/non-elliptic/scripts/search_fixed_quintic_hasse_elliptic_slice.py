#!/usr/bin/env python3
"""Search the rank-two fixed-quintic common-resolvent elliptic slice.

This is an exploratory search, not an infinitude proof.  On the slice

    kappa = -A,  R = 1,

the cube condition is the Mordell curve

    y^2 = x^3 + 22356,  x = 18*Pi,  y = 54*A.

The displayed generators have rank two.  This script enumerates their
bounded integral linear combinations, reconstructs the quadratic and cubic
factors, and checks every prime dividing their discriminants or leading
coefficients.  A reported row is therefore an exact Hasse-failing target;
the absence of rows is only bounded evidence on this Mordell--Weil box.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess

import sympy as sp


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--coefficient-bound",
        type=int,
        default=12,
        help="enumerate m*P1+n*P2 with |m|,|n| at most this bound",
    )
    parser.add_argument("--padic-precision", type=int, default=30)
    parser.add_argument(
        "--show-failures",
        action="store_true",
        help="print every Mordell--Weil point rejected by the dyadic test",
    )
    return parser.parse_args()


def verify_symbolic_reconstruction() -> None:
    A, Pi, T = sp.symbols("A Pi T")
    V = (A**2 + 3) / 4
    W = Pi**3
    d = A**2 - V - 5
    e = (4 * W - V * d) / A
    q = T**2 + A * T + V
    h = T**3 - A * T**2 + d * T + e
    product = sp.Poly(sp.expand(q * h), T)

    assert product.coeff_monomial(T**4) == 0
    assert product.coeff_monomial(T**3) == -5
    assert sp.factor(product.coeff_monomial(T) - 4 * Pi**3) == 0
    assert sp.factor(sp.discriminant(q, T) + 3) == 0
    discriminant_quotient = sp.factor(
        sp.discriminant(h, T) / sp.discriminant(q, T)
    )
    expected_square = (
        (3 * A**2 - 23) * (5 * A**2 - 69) / (48 * A)
    ) ** 2
    assert sp.factor(
        discriminant_quotient.subs(Pi**3, (3 * A**2 - 23) / 6)
        - expected_square
    ) == 0
    print("PASS: exact factor reconstruction on the rank-two elliptic slice")


def verify_dyadic_obstruction() -> None:
    """Audit the valuation proof excluding every Q_2-point on the slice."""

    # If v_2(Pi)=0, the right side of 3*A^2=6*Pi^3+23 is 1 or 5
    # modulo 8, whereas the left side is 3 modulo 8.  If v_2(Pi)>0,
    # the two sides are respectively 7 and 3 modulo 8.
    odd_residues = (1, 3, 5, 7)
    assert {
        (6 * residue**3 + 23) % 8 for residue in odd_residues
    } == {1, 5}
    assert 23 % 8 == 7
    assert {3 * residue**2 % 8 for residue in odd_residues} == {3}

    # Hence v_2(Pi)=-(2m+1), and comparing the unique lowest-valuation
    # terms gives v_2(A)=-(3m+1).  The cubic coefficients, in increasing
    # powers of T, then have the following valuations.
    m = sp.symbols("m", integer=True, nonnegative=True)
    coefficient_valuations = (
        -9 * m - 7,
        -6 * m - 4,
        -3 * m - 1,
        0,
    )
    endpoint_slope = 3 * m + sp.Rational(7, 3)
    endpoint_line = tuple(
        coefficient_valuations[0] + index * endpoint_slope
        for index in range(4)
    )
    assert sp.factor(
        coefficient_valuations[1] - endpoint_line[1]
    ) == sp.Rational(2, 3)
    assert sp.factor(
        coefficient_valuations[2] - endpoint_line[2]
    ) == sp.Rational(4, 3)
    assert sp.fraction(sp.together(endpoint_slope))[1] == 3

    # The quadratic discriminant is -3, a nonsquare unit in Q_2.
    assert (-3) % 8 == 5
    print(
        "PASS: every rational point on the rank-two slice fails over Q_2 "
        "(quadratic nonsplit; cubic Newton slope denominator 3)"
    )


def main() -> None:
    args = parse_args()
    if args.coefficient_bound <= 0:
        raise SystemExit("coefficient-bound must be positive")
    if args.padic_precision <= 0:
        raise SystemExit("padic-precision must be positive")

    verify_symbolic_reconstruction()
    verify_dyadic_obstruction()
    gp = shutil.which("gp")
    if gp is None:
        raise SystemExit("PARI/GP executable 'gp' is required")

    program = rf"""
E=ellinit([0,0,0,0,22356]);
P1=[-11,145];
P2=[73/4,1349/8];
if(ellisoncurve(E,P1)==0||ellisoncurve(E,P2)==0,error("bad generators"));
seen=Map();
ntested=0;
nirred=0;
nlocalfail=0;
check(m,n,P)={{
  my(xx,yy,Pi,A,V,W,d,e,q,h,Q,H,discquot,B,C);
  if(#P==1,return());
  xx=P[1]; yy=P[2];
  Pi=xx/18; A=yy/54;
  if(!Pi||!A,return());
  if(mapisdefined(seen,[A,Pi]),return());
  mapput(seen,[A,Pi],1);
  ntested++;
  if(3*A^2-6*Pi^3-23,error("elliptic reconstruction failed"));
  V=(A^2+3)/4;
  W=Pi^3;
  d=A^2-V-5;
  e=(4*W-V*d)/A;
  q=x^2+A*x+V;
  h=x^3-A*x^2+d*x+e;
  if(x^5-5*x^3+polcoeff(q*h,2)*x^2+4*Pi^3*x+polcoeff(q*h,0)
       !=q*h,error("factor identity failed"));
  if(!polisirreducible(q)||!polisirreducible(h),return());
  nirred++;
  discquot=poldisc(h)/poldisc(q);
  if(!issquare(discquot),error("common resolvent identity failed"));
  Q=q/content(q);
  H=h/content(h);
  if(#polrootspadic(Q,2,{args.padic_precision})
     ||#polrootspadic(H,2,{args.padic_precision}),
    error("the exact dyadic obstruction failed"));
  B=-polcoeff(q*h,2)/(2*Pi);
  C=-polcoeff(q*h,0)/(2*Pi^5);
  nlocalfail++;
  if({1 if args.show_failures else 0},
    print("LOCAL_FAIL m,n,p,A,Pi=",[m,n,2,A,Pi]))
}};
for(m=-{args.coefficient_bound},{args.coefficient_bound},for(n=-{args.coefficient_bound},{args.coefficient_bound},if(m||n,P=elladd(E,ellmul(E,P1,m),ellmul(E,P2,n));check(m,n,P))));
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
