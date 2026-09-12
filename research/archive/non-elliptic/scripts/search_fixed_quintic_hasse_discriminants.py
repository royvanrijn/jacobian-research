#!/usr/bin/env python3
"""Search Hasse fibers while varying the shared quadratic resolvent.

This is a bounded experiment, not an infinitude proof.  In rational
normalized quadratic-cubic coordinates it imposes

    disc(q) = D R^2

for squarefree nonsquare integers D.  It then requires the cubic
discriminant to have the same square class, both factors to be irreducible,
and at least one factor to have a root over Q_p at every prime visible in
the two polynomial discriminants and coefficient denominators.

PARI's p-adic root test is a search filter.  A reported target still needs
an independent exact Hensel audit before it becomes a certificate.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--d-bound", type=int, default=30)
    parser.add_argument("--d-value", type=int)
    parser.add_argument("--r-bound", type=int, default=16)
    parser.add_argument("--r-denominator", type=int, default=1)
    parser.add_argument("--a-bound", type=int, default=40)
    parser.add_argument("--a-denominator", type=int, default=1)
    parser.add_argument("--pi-bound", type=int, default=30)
    parser.add_argument("--pi-denominator", type=int, default=1)
    parser.add_argument("--padic-precision", type=int, default=20)
    parser.add_argument("--show-failures", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    for name, value in vars(args).items():
        if (
            value is not None
            and name not in {"d_value", "show_failures"}
            and value <= 0
        ):
            raise SystemExit(f"{name} must be positive")
    if args.d_value in (0, 1):
        raise SystemExit("d_value must be a squarefree nonsquare integer")

    gp = shutil.which("gp")
    if gp is None:
        raise SystemExit("PARI/GP executable 'gp' is required")

    parameter_loop = (
        f"for(rd=1,{args.r_denominator},"
        f"for(rn=1,{args.r_bound},"
        "if(gcd(rn,rd)==1,R=rn/rd;"
        f"for(ad=1,{args.a_denominator},"
        f"for(an=-{args.a_bound},{args.a_bound},"
        "if(an&&gcd(abs(an),ad)==1,A=an/ad;"
        f"for(pd=1,{args.pi_denominator},"
        f"for(pn=-{args.pi_bound},{args.pi_bound},"
        "if(pn&&gcd(abs(pn),pd)==1,P=pn/pd;"
        "check(D,R,A,P))))))))))"
    )
    if args.d_value is None:
        search_loop = (
            f"for(D=-{args.d_bound},{args.d_bound},"
            "if(D&&D!=1&&issquarefree(D),"
            f"{parameter_loop}));"
        )
    else:
        search_loop = (
            f"D={args.d_value};"
            'if(!issquarefree(D),error("d_value must be squarefree"));'
            f"{parameter_loop};"
        )

    program = rf"""
candidate_count=0;
local_count=0;

Mvalue(A,V,W)={{
  3*A^8-24*A^6*V-50*A^6
  +70*A^4*V^2+270*A^4*V+56*A^4*W+275*A^4
  -76*A^2*V^3-510*A^2*V^2-288*A^2*V*W
  -750*A^2*V-360*A^2*W-500*A^2
  +27*V^4+270*V^3+216*V^2*W+675*V^2
  +1080*V*W+432*W^2
}};

badprimelist(q,h)={{
  my(dq,dh,n);
  dq=poldisc(q);
  dh=poldisc(h);
  n=abs(numerator(dq)*denominator(dq)
        *numerator(dh)*denominator(dh)
        *denominator(content(q))*denominator(content(h)));
  if(n<=1,return([]));
  factor(n)[,1]
}};

firstuncovered(q,h)={{
  my(ps=badprimelist(q,h));
  for(i=1,#ps,
    if(!#polrootspadic(q,ps[i],{args.padic_precision})
       &&!#polrootspadic(h,ps[i],{args.padic_precision}),
       return(ps[i]))
  );
  0
}};

check(D,R,A,P)={{
  my(V,W,M,s,d,e,q,h,badp,beta,B,C);
  V=(A^2-D*R^2)/4;
  W=P^3;
  M=Mvalue(A,V,W);
  if(!issquare(-D*M,&s),return());
  d=A^2-V-5;
  e=(4*W-V*d)/A;
  q=x^2+A*x+V;
  h=x^3-A*x^2+d*x+e;
  if(!polisirreducible(q)||!polisirreducible(h),return());
  candidate_count++;
  badp=firstuncovered(q,h);
  if(badp,
    if({1 if args.show_failures else 0},
      print("LOCAL_FAILURE p,D,R,A,Pi=",[badp,D,R,A,P]));
    return()
  );
  beta=e+A*(d-V);
  B=-beta/(2*P);
  C=-V*e/(2*P^5);
  print("LOCAL_CANDIDATE D,R,A,Pi=",[D,R,A,P]);
  print("  target Pi,B,C=",[P,B,C]);
  print("  q=",q);
  print("  h=",h);
  print("  bad_primes=",badprimelist(q,h));
  local_count++
}};

{search_loop}
print("SEARCH_COMPLETE common_resolvent_irreducible=",candidate_count," all_bad_prime_candidates=",local_count);
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
