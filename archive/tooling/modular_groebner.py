#!/usr/bin/env python3
"""Run the nonlinear 3D rediscovery residue with msolve F4."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.msolve import run,available

a20,a11,a30,b10,b01,b21,t,r,c=sp.symbols("a20 a11 a30 b10 b01 b21 t r c")
A=a20*t**2+a11*r*t+a30*c*t**3
B=b10*t+b01*r+b21*c*t**2
res=sp.Poly(sp.expand(sp.Matrix([A,B,c]).jacobian((t,r,c)).det()-r/2),t,r,c)
equations=[]
for equation in res.coeffs()+[b01-1,b10-4,b21+3]:
    _denominator, integral = sp.Poly(equation, a20,a11,a30,b10,b01,b21, domain=sp.QQ).clear_denoms()
    equations.append(integral.as_expr())
variables=(a20,a11,a30,b10,b01,b21)
assert available()
print("F4 backend: msolve")
for prime in (1000003,1000033,1000037):
    result=run(equations,variables,prime=prime)
    assert result.returncode==0 and not result.empty
    print(prime,"F4 completed; reduced basis present:","#Reduced Groebner basis" in result.output)
