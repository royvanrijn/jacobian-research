#!/usr/bin/env python3
"""Validation ladder for charts, cancellations, fake modular roots, and filters."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.triangular import TriangularChart,u,v,x,y
from jcsearch.filters import generic_resultant_degree,line_restriction_metrics,candidate_metrics
from jcsearch.msolve import available,run
from jcsearch.weighted import construct,w
from jcsearch.toric import catalog,u as tu,v as tv

for p in (0,v,v**2+v,v**4):
    chart=TriangularChart(sp.sympify(p));assert chart.jacobian==x**-2
    # Planted cancellation construction: A=1/v, B=u-p(v).
    A,B=1/v,u-chart.p
    assert sp.factor(sp.Matrix([A,B]).jacobian((u,v)).det()-v**-2)==0
    assert chart.pullback(A)==x and chart.pullback(B)==y

# Random-looking triangular automorphism and generic degree recognition.
P=x+(y+x**2)**3;Q=y+x**2
assert sp.factor(sp.Matrix([P,Q]).jacobian((x,y)).det())==1
assert generic_resultant_degree(P,Q)==1
assert candidate_metrics(P,Q)["degree_flags"]["birational"]
assert line_restriction_metrics(x,y)["y=0"]["linear_coordinate"]

# Characteristic-p fake: x-x^p has derivative 1 mod p, but not over Q.
p=5;fake=x-x**p
assert sp.diff(fake,x)!=1 and sp.Poly(sp.diff(fake,x)-1,x,modulus=p).is_zero
assert available()
sat=run([x*(y-1),x*(x-2)],(x,y),prime=1000003,saturate=x)
assert sat.returncode==0 and not sat.empty and "x" in sat.output and "y" in sat.output
for seed in (2*w-3*w**2,w-2*w**3):
    assert construct(seed)
for chart in catalog():
    det=chart.a*chart.d-chart.b*chart.c
    X=tu**(chart.d//det)*tv**(-chart.b//det);Y=tu**(-chart.c//det)*tv**(chart.a//det)
    assert sp.factor(sp.Matrix([X,Y]).jacobian((tu,tv)).det()-chart.reciprocal_outer_jacobian)==0
    assert chart.pullback(X)==x and chart.pullback(Y)==y
print("PASS: planted chart cancellations")
print("PASS: polynomial automorphism and generic-degree filter")
print("PASS: characteristic-p fake is rejected over Q")
print("PASS: msolve F4 backend available")
print("PASS: msolve F4SAT component removal")
print("PASS: two weighted-lift 3D seed deformations")
print("PASS: planted inverses for three two-divisor toric charts")
