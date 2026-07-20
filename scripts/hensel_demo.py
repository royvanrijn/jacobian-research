#!/usr/bin/env python3
"""Exercise Hensel lifting plus rational reconstruction on an exact system."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.arithmetic import hensel_newton, rational_reconstruct

a,b=sp.symbols("a b")
polys=[3*a-1, 5*b-2]
p=1009
root=[pow(3,-1,p),2*pow(5,-1,p)%p]
lifted,modulus=hensel_newton(polys,(a,b),root,p,levels=3)
reconstructed=[rational_reconstruct(v,modulus) for v in lifted]
assert reconstructed==[sp.Rational(1,3),sp.Rational(2,5)]
assert all(f.subs({a:reconstructed[0],b:reconstructed[1]})==0 for f in polys)
print("PASS: lifted modulus",modulus,"and reconstructed",reconstructed)
