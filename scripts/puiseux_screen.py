#!/usr/bin/env python3
"""Puiseux-first screen for explicit candidate polynomials."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import sympy as sp
from jcsearch.filters import puiseux_leading_scan
x,y=sp.symbols("x y")
# Controls: identity has no bad branch; a deliberately non-Keller pair may.
for name,P,Q in (("identity",x,y),("triangular",x+y**3,y)):
 print(name,puiseux_leading_scan(P,Q))

