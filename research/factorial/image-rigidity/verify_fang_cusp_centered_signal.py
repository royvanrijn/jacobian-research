"""Exact local coefficients of the d=3 cusp centered insertion.

Run: .venv/bin/python research/factorial/image-rigidity/verify_fang_cusp_centered_signal.py
"""

import sympy as s

t, S = s.symbols('t S')
phi = s.log(S) + s.log(16 + S**3) - S
f = 1 / (16 + S**3)
assert s.factor(s.diff(phi, S)) == -(S - 2)**3 * (S + 2) / (S * (S**3 + 16))
phicoeff = s.series(phi.subs(S, 2 + t), t, 0, 7).removeO().expand()
fcoeff = s.series(f.subs(S, 2 + t), t, 0, 4).removeO().expand()
assert phicoeff.coeff(t, 4) == -s.Rational(1, 48)
assert phicoeff.coeff(t, 5) == s.Rational(1, 80)
assert fcoeff.coeff(t, 0) == s.Rational(1, 24)
assert fcoeff.coeff(t, 1) == -s.Rational(1, 48)
assert fcoeff.coeff(t, 2) == 0
moment_coefficient = fcoeff.coeff(t, 1) * phicoeff.coeff(t, 5) * s.Rational(3, 4) * 48
assert moment_coefficient == -s.Rational(3, 320)
print('PASS d=3 cusp derivative, quartic/quintic phase, insertion and nonzero centered signal')
