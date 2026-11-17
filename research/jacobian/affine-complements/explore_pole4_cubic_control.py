"""Bounded exploratory lift of one base-point-free cubic symbol family.

This is a control family, not a complete classification at pole pair (3,4).
"""
from pathlib import Path
import sympy as S,json,hashlib,time
HERE=Path(__file__).resolve().parent
start=time.monotonic()
z,w=S.symbols("z w")
t0,a0,b0,c0,k,ell=S.symbols("t0 a0 b0 c0 k ell")
g1=w/z
g2=w*w/z
g3=(z+w*w)*w/z
g4=(z+w*w)**2/z
T=g4+t0
A=3*g3+a0
B=3*g2+b0
C=S.Rational(8,3)*g1+S.Rational(2,3)*k*A+ell*T+c0
P=T**3+T*A+B
Q=T**4+S.Rational(4,3)*T*T*A+S.Rational(4,3)*T*B+S.Rational(2,9)*A*A+k*T*T+C
J=S.cancel(S.diff(P,z)*S.diff(Q,w)-S.diff(P,w)*S.diff(Q,z))
pole4=S.cancel(z**4*J)
num,den=S.fraction(pole4)
assert S.cancel(den).subs(z,0)!=0
coefficient=S.cancel(pole4.subs(z,0))
equations=S.Poly(coefficient,w).all_coeffs()
equations=[S.factor(v) for v in equations if v!=0]
report=dict(as_of="2026-11-14",scope="One cubic-symbol control, not the full pole-four family",
            pole_four_coefficient=str(S.factor(coefficient)),
            equations=list(map(str,equations)),
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            elapsed_seconds=round(time.monotonic()-start,6))
(HERE/"pole4_cubic_control.json").write_text(json.dumps(report,indent=2)+"\n")
print("pole-four coefficient",S.factor(coefficient))
print("necessary equations",equations)

