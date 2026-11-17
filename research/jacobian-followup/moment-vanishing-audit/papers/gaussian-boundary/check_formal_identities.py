"""Finite universal algebra only; does not certify the written formal frame.

The first check was run inline during revision 2. This file preserves that
transcription for review; it was not rerun after the successful inline check.
It deliberately does not overwrite the retained check receipt.
"""
import sympy as s


def check():
    a, B, C, F, dA, dB, dC = s.symbols('a B C F dA dB dC', nonzero=True)
    A = a**3 / 3
    m2 = -(2*B*F + C)/(3*A)
    m3 = -(2*B*m2 + C*F + 1)/(3*A)
    m4 = -(2*B*m3 + C*m2 + 2*F)/(3*A)
    ip = s.factor(dA*m3 + dB*m2 + dC*F)
    fp = s.factor(dA*m4 + dB*m3 + dC*m2 - F*ip)
    X = (B**2 - 3*A*C)/a**4
    airy = -a*F - B/a**2

    def derivative(value):
        return (s.diff(value, a)*dA/a**2 + s.diff(value, B)*dB
                + s.diff(value, C)*dC + s.diff(value, F)*fp)

    xp = s.factor(derivative(X))
    assert s.factor(derivative(airy) - xp*(X-airy**2)) == 0
    assert s.factor(s.diff(ip, F) + a*xp) == 0

    z, t, aa, T = s.symbols('z t a T')
    cs = s.symbols('g0:10')
    G = (cs[0] + cs[1]*z + cs[2]*T + cs[3]*z**2 + cs[4]*z*T
         + cs[5]*T**2 + cs[6]*z**3 + cs[7]*z**2*T
         + cs[8]*z*T**2 + cs[9]*T**3)
    V = G.subs(T, -aa-z**2)
    L = s.diff(G, T).subs(T, -aa-z**2)
    b0 = z*L - (aa+z**2)**2/2
    b1 = z**2*s.diff(G, T, 2).subs(T, -aa-z**2)/2 + z*(aa+z**2)
    b2 = cs[9]*z**3 - z**2/2
    phi = t*G.subs(T, z/t-aa-z**2) - (z/t-aa-z**2)**2/2
    assert s.expand(phi-(t*V+b0+b1/t+b2/t**2)) == 0
    assert s.expand(V).coeff(z, 6) == -cs[9]
    assert s.expand(V).coeff(z, 5) == cs[8]
    assert s.expand(V).coeff(z, 4) == cs[5]-cs[7]-3*aa*cs[9]
    assert s.expand(b0.subs({cs[9]: 0, cs[8]: 0,
                            cs[7]: cs[5]})).coeff(z, 4) == -s.Rational(1, 2)

    zz, v = s.symbols('zeta v')
    for n in (1, 2, 3):
        H = v*zz**n - zz**4/2
        relation = zz**(4-n) - s.Rational(n, 2)*v
        action = (1-s.Rational(n, 4))*v*zz**n
        assert s.rem(s.expand(H-action), relation, zz) == 0
        assert s.rem(s.expand(s.diff(H, zz, 2)-2*(n-4)*zz**2),
                     relation, zz) == 0


if __name__ == '__main__':
    check()
    print('Finite universal identities pass; completeness and analytic estimates '
          'require the written proofs. No receipt or theorem assurance changed.')
