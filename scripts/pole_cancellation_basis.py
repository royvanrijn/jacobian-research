#!/usr/bin/env python3
"""Compute the exact sparse kernel of pole-cancellation constraints."""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.laurent import forbidden_pole_equations, linear_certificate

x,y,t,r=sp.symbols("x y t r")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--degree",type=int,default=6)
    parser.add_argument("--pole-window",type=int,default=3)
    args=parser.parse_args()
    support=[]
    for i in range(args.degree+1):
        for j in range(-args.degree, args.pole_window+1):
            support.append((i,j))
    coeffs=sp.symbols(f"c0:{len(support)}")
    expression=sum(c*t**i*r**j for c,(i,j) in zip(coeffs,support))
    pullback=sp.expand(expression.subs({t:y+1/x,r:2/x}, simultaneous=True))
    equations=forbidden_pole_equations(pullback,x,(y,))
    cert,matrix,rhs=linear_certificate(equations,coeffs)
    nullity=len(coeffs)-cert.rank
    print("chart: t=y+1/x, r=2/x")
    print("support size:",len(support),"pole equations:",len(equations))
    print("exact rank:",cert.rank,"cancellation-kernel dimension:",nullity)
    # Verify that every ordinary monomial x^a y^b in the requested degree has
    # a Laurent representative 2^a r^-a (t-r/2)^b inside the window when possible.
    checked=0
    for total in range(args.degree+1):
        for b in range(total+1):
            a=total-b
            representative=sp.expand(2**a*r**(-a)*(t-r/2)**b)
            exponents={(int(term.as_powers_dict().get(t,0)),int(term.as_powers_dict().get(r,0)))
                       for term in sp.Add.make_args(representative)}
            if exponents.issubset(set(support)):
                assert sp.cancel(representative.subs({t:y+1/x,r:2/x}, simultaneous=True)-x**a*y**b)==0
                checked+=1
    print("verified polynomial monomial representatives:",checked)

if __name__=="__main__": main()

