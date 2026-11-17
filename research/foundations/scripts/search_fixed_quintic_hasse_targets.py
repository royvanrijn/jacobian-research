#!/usr/bin/env python3
"""Bounded PARI/GP search for low-height Hasse fibers in the fixed quintic.

The search uses the quadratic-times-cubic incidence from
``FIXED_QUINTIC_MODULI_DOMINANCE.md``.  For fixed ``(Pi,a)``, equality of
the quadratic and cubic discriminant squareclasses is a genus-two
hyperelliptic equation in ``b``.  PARI's ``hyperellratpoints`` enumerates
its rational points in the requested height box.  Every candidate is then
tested for irreducibility and for roots over all possibly ramified ``Q_p``.

This is a bounded search, not a proof of global height minimality.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pi-denominator", type=int, default=8)
    parser.add_argument("--pi-numerator", type=int, default=16)
    parser.add_argument("--pi-absolute", type=int, default=4)
    parser.add_argument("--a-denominator", type=int, default=8)
    parser.add_argument("--a-numerator", type=int, default=120)
    parser.add_argument("--b-height", type=int, default=10_000)
    parser.add_argument(
        "--target-height",
        type=int,
        default=458_080,
        help="report only targets of strictly smaller projective height",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    gp = shutil.which("gp")
    if gp is None:
        raise SystemExit("PARI/GP executable 'gp' is required")

    for name, value in vars(args).items():
        if value <= 0:
            raise SystemExit(f"{name} must be positive")

    program = rf"""
Nfun(p,a,b)=3*a^8*p^12-24*a^6*b*p^12-50*a^6*p^8+70*a^4*b^2*p^12+270*a^4*b*p^8+56*a^4*p^7+275*a^4*p^4-76*a^2*b^3*p^12-510*a^2*b^2*p^8-288*a^2*b*p^7-750*a^2*b*p^4-360*a^2*p^3-500*a^2+27*b^4*p^12+270*b^3*p^8+216*b^2*p^7+675*b^2*p^4+1080*b*p^3+432*p^2;
height3(p,b,c)={{my(d=lcm([denominator(p),denominator(b),denominator(c)]),v=[d,d*p,d*b,d*c],g=0);for(i=1,4,g=gcd(g,abs(v[i])));v=v/g;[vecmax(abs(v)),v]}};
localok(q,h)={{my(Q=q/content(q),H=h/content(h),n=abs(poldisc(Q)*poldisc(H)*pollead(Q)*pollead(H)),fa=factor(n),pp);for(i=1,matsize(fa)[1],pp=fa[i,1];if(#polrootspadic(Q,pp,20)==0&&#polrootspadic(H,pp,20)==0,return(0)));1}};
evalpoint(p,a,b,cap)={{my(B,C,q,h,ht);B=-(p^5*a^4-3*p^5*a^2*b+p^5*b^2-5*p*a^2+5*p*b+4)/(2*a);C=b*(p^5*a^2*b-p^5*b^2-5*p*b-4)/(2*a);ht=height3(p,B,C);if(ht[1]>=cap,return());q=x^2+a*x+b;h=p^5*x^3-p^5*a*x^2+(p^5*(a^2-b)-5*p)*x+(-p^5*a^2*b+p^5*b^2+5*p*b+4)/a;if(polisirreducible(q)&&polisirreducible(h)&&localok(q,h),print("VALID height=",ht[1]," projective=",ht[2]," Pi,a,b,B,C=",[p,a,b,B,C]))}};
searchone(p,a,H,cap)={{my(F,v,bb);F=-Nfun(p,a,x)*(a^2-4*x);v=hyperellratpoints(F,H);for(i=1,#v,if(v[i][2]>0,bb=v[i][1];if(a^2-4*bb!=0,evalpoint(p,a,bb,cap))))}};
for(q=1,{args.pi_denominator},for(P=-{args.pi_numerator},{args.pi_numerator},if(P!=0&&gcd(abs(P),q)==1&&abs(P/q)<={args.pi_absolute},p=P/q;for(d=1,{args.a_denominator},for(A=-{args.a_numerator},{args.a_numerator},if(A!=0&&gcd(abs(A),d)==1,a=A/d;searchone(p,a,{args.b_height},{args.target_height})))))));
print("SEARCH_COMPLETE");
"""

    print(
        "Searching "
        f"|num(Pi)|<={args.pi_numerator}, den(Pi)<={args.pi_denominator}, "
        f"|Pi|<={args.pi_absolute}; "
        f"|num(a)|<={args.a_numerator}, den(a)<={args.a_denominator}; "
        f"H(b)<={args.b_height}; H(target)<{args.target_height}"
    )
    completed = subprocess.run(
        [gp, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, end="")


if __name__ == "__main__":
    main()
