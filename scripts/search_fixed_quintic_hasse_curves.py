#!/usr/bin/env python3
"""Search rational curve slices of the fixed-quintic Hasse incidence.

This is an experiment, not an infinitude proof.  It first verifies the exact
normalized conic parametrization from FIXED_QUINTIC_MODULI_DOMINANCE.md and
then asks PARI/GP to search proportional conic parameters K=cA.  Candidates
must have irreducible quadratic and cubic factors and a cubic root over
Q_2, Q_3, and Q_5.  The remaining ramified primes still require an exact
candidate-by-candidate audit.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess

import sympy as sp


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r-denominator", type=int, default=4)
    parser.add_argument("--r-numerator", type=int, default=30)
    parser.add_argument("--a-denominator", type=int, default=4)
    parser.add_argument("--a-numerator", type=int, default=30)
    parser.add_argument("--pi-bound", type=int, default=30)
    return parser.parse_args()


def verify_parametrization() -> None:
    A, V, W, R, K = sp.symbols("A V W R K")
    cubic_discriminant_numerator = (
        3 * A**8
        - 24 * A**6 * V
        - 50 * A**6
        + 70 * A**4 * V**2
        + 270 * A**4 * V
        + 56 * A**4 * W
        + 275 * A**4
        - 76 * A**2 * V**3
        - 510 * A**2 * V**2
        - 288 * A**2 * V * W
        - 750 * A**2 * V
        - 360 * A**2 * W
        - 500 * A**2
        + 27 * V**4
        + 270 * V**3
        + 216 * V**2 * W
        + 675 * V**2
        + 1080 * V * W
        + 432 * W**2
    )
    H = 2 * A**2 - 3 * V - 15
    assert sp.factor(
        sp.discriminant(cubic_discriminant_numerator, W)
        + 256 * A**2 * H**3
    ) == 0

    specialized_v = (A**2 + 3 * R**2) / 4
    specialized_h = sp.factor(H.subs(V, specialized_v))
    linear_coefficient = sp.factor(
        sp.diff(cubic_discriminant_numerator, W)
        .subs({V: specialized_v, W: 0})
    )
    parametrized_w = sp.factor(
        (
            8
            * A
            * specialized_h
            * (specialized_h / K - K)
            - linear_coefficient
        )
        / 864
    )
    parametrized_y = sp.factor(
        R
        * A
        * specialized_h
        * (K + specialized_h / K)
        / 3
    )
    square_equation = -(
        A**2 - 4 * V
    ) * cubic_discriminant_numerator
    assert sp.factor(
        square_equation.subs({V: specialized_v, W: parametrized_w})
        - parametrized_y**2
    ) == 0

    c = sp.symbols("c")
    proportional_w = sp.factor(parametrized_w.subs(K, c * A))
    a4_coefficient = sp.Poly(
        sp.together(proportional_w).as_numer_denom()[0],
        A,
    ).coeff_monomial(A**4)
    assert sp.expand(
        a4_coefficient + 5 * (c + 1) * (4 * c - 5)
    ) == 0

    assert sp.factor(
        proportional_w.subs(c, -1)
        - (3 * R**2 + 5) * (3 * A**2 - 3 * R**2 - 20) / 48
    ) == 0
    assert sp.factor(
        proportional_w.subs(c, sp.Rational(5, 4))
        - (
            135 * A**2 * R**2
            - 99 * R**4
            - 420 * R**2
            + 1600
        )
        / 960
    ) == 0

    print("PASS exact normalized discriminant and conic parametrization")
    print("PASS elliptic degenerations occur at K/A=-1 and K/A=5/4")


def main() -> None:
    args = parse_args()
    for name, value in vars(args).items():
        if value <= 0:
            raise SystemExit(f"{name} must be positive")

    verify_parametrization()
    gp = shutil.which("gp")
    if gp is None:
        raise SystemExit("PARI/GP executable 'gp' is required")

    program = rf"""
Eminus=ellinit([0,0,0,0,22356]);
rminus=ellrank(Eminus);
if(rminus[1]!=2||rminus[2]!=2,error("unexpected rank interval for K/A=-1, R=1: ",rminus));
Eplus=ellinit([0,0,0,0,30464/15]);
rplus=ellrank(Eplus);
if(rplus[1]!=1||rplus[2]!=1,error("unexpected rank interval for K/A=5/4, R=4: ",rplus));
print("PASS PARI rank intervals: K/A=-1,R=1 has rank 2; K/A=5/4,R=4 has rank 1");
candidate_count=0;
small_prime_count=0;
check(R,A,P)={{
  my(qa,qb,qc,DD,ss,c,V,d,e,q,h,r2,r3,r5);
  qa=20*A^4-36*A^2*R^2-240*A^2;
  qb=-5*A^4-270*A^2*R^2-180*A^2
     +1728*P^3+243*R^4+1620*R^2;
  qc=-25*A^4+90*A^2*R^2+600*A^2
     -81*R^4-1080*R^2-3600;
  if(!qa,return());
  DD=qb^2-4*qa*qc;
  ss=0;
  if(DD<0||!issquare(DD,&ss),return());
  forstep(sign=-1,1,2,
    c=(-qb+sign*ss)/(2*qa);
    if(!c,next());
    V=(A^2+3*R^2)/4;
    d=A^2-V-5;
    e=(4*P^3-V*d)/A;
    q=x^2+A*x+V;
    h=x^3-A*x^2+d*x+e;
    if(polisirreducible(q)&&polisirreducible(h),
      candidate_count++;
      r2=#polrootspadic(h,2,15);
      r3=#polrootspadic(h,3,15);
      r5=#polrootspadic(h,5,15);
      if(r2&&r3&&r5,
        print("SMALL_PRIME_CANDIDATE c,R,A,Pi=",[c,R,A,P]);
        small_prime_count++
      )
    )
  )
}};
for(rd=1,{args.r_denominator},for(rn=1,{args.r_numerator},if(gcd(rn,rd)==1,R=rn/rd;for(ad=1,{args.a_denominator},for(an=-{args.a_numerator},{args.a_numerator},if(an&&gcd(abs(an),ad)==1,A=an/ad;for(P=-{args.pi_bound},{args.pi_bound},if(P,check(R,A,P)))))))));
print("SEARCH_COMPLETE candidates=",candidate_count," small_prime_candidates=",small_prime_count);
"""
    completed = subprocess.run(
        [gp, "-q", "-f"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    if "***" in completed.stderr:
        raise SystemExit(completed.stderr)
    print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="")


if __name__ == "__main__":
    main()
