"""Complete highest-cubic chart: gamma3=1, gamma2=0, then eliminate mu2.

All cubic phases are covered up to oriented isometry and phase amplitude.
Positive weights are3 minus the circular charge. No purity classification
or global safe-cover membership is asserted by this adapter alone.
"""
import sympy as S

from conformal_euler_cubic import harmonic_basis, z, w, t

LW, LT, LZ = S.symbols('ell_w ell_t ell_z')
BM2, BM1, B0, B1, B2 = S.symbols('bminus2 bminus1 b0 b1 b2')
CM3, CM2, CM1, C0, C1 = S.symbols('cminus3 cminus2 cminus1 c0 c1')
VARIABLES13 = (LW, LT, LZ, BM2, BM1, B0, B1, B2, CM3, CM2, CM1, C0, C1)
WEIGHTS13 = (4, 3, 2, 5, 4, 3, 2, 1, 6, 5, 4, 3, 2)
VARIABLES = tuple(v for v in VARIABLES13 if v != CM3)
WEIGHTS = tuple(k for v, k in zip(VARIABLES13, WEIGHTS13) if v != CM3)
BASE = (LZ, B1, B2, C0, C1)
DEPENDENT = (LW, LT, BM2, BM1, B0, CM2, CM1)
DEPENDENT_WEIGHTS = (4, 3, 5, 4, 3, 5, 4)


def phase13():
    basis = harmonic_basis()
    return S.expand(LW*w+LT*t+LZ*z
                    +sum(v*basis[2][j] for v, j in zip((BM2, BM1, B0, B1, B2), range(-2, 3)))
                    +z**3+sum(v*basis[3][j] for v, j in zip((CM3, CM2, CM1, C0, C1), range(-3, 2))))


def eliminated_cminus3():
    return S.expand(-S.Rational(35, 4)*(S.Rational(2, 3)*LZ*LW+LT*LT/3
                      +S.Rational(4, 15)*B2*BM2+S.Rational(2, 15)*B1*BM1+B0*B0/5
                      +C1*CM1/21+C0*C0/7))


def phase():
    return S.expand(phase13().subs(CM3, eliminated_cminus3()))


def central_phase():
    return S.expand(phase().subs({v: 0 for v in BASE}))
