#!/usr/bin/env python3
"""Verify K_n = polynomials of total degree <= n for translated chart boxes."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.translated import TranslatedTwoDivisorChart,cancellation_system,kernel_expressions,x,y

chart=TranslatedTwoDivisorChart()
for n in range(1,5):
    support=tuple((i,j) for i in range(-n,n+1) for j in range(-n,n+1))
    system=cancellation_system(support);basis=kernel_expressions(system)
    expected=(n+1)*(n+2)//2
    assert len(basis)==expected
    pullbacks=[sp.Poly(chart.pullback(f),x,y) for f in basis]
    assert all(p.total_degree()<=n for p in pullbacks)
    coefficient_matrix=[]
    monomials=[x**a*y**b for total in range(n+1) for a in range(total+1) for b in [total-a]]
    for p in pullbacks:
      coefficient_matrix.append([p.coeff_monomial(m) for m in monomials])
    assert sp.Matrix(coefficient_matrix).rank()==expected
    print(f"PASS n={n}: box size {(2*n+1)**2}, kernel dimension {expected}, full degree<={n} polynomial space")
