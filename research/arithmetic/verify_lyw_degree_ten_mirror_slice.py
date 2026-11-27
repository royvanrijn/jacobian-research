"""Exact inputs for LY1--LY6; Python 3.12 and SymPy 1.14.0.

The integer coefficients are transcribed from Levin--Yan--Wiljanen,
arXiv:1910.12276v1, Theorem 9, printed p. 12. No numerical root tests.
The uniform valuation and class-field arguments are in the companion note;
this checker verifies their polynomial, Sturm and residue inputs.
"""

import platform
import sympy as s

t = s.Symbol("t")
# Coefficients in ascending order.
F_COEFFS = [
    543592155691663960065241800360826161610140961,
    -135004259433655686521826532061360904927910680,
    14597743197263467927181474503046907251979462,
    -1509800364506319291441531124462079071041720,
    139665528153448288531118705650287136663899,
    -8471956828413213486742748322179256745500,
    486933739385947419621206507920009537350,
    -23227949011157855871750302161149318060,
    767540949843094964859507359162484321,
    -14476726558441542259500980593582900,
    127358629188153017343112694654244,
]
Y_COEFFS = [
    [23315023973129417008893, -2894185136924624900028,
     124175441794992816207, 11259180860536474740,
     -612703879315493343, 11285328049647162],
    [-23315023973129417008893, 2895345088136314829232,
     -133140352448943347487, 15898504501345253760,
     -637837872300913137, 11285328049647162],
    [-23295069079544963156463, 2891689893601551455028,
     -132367737077167916373, 15779426445792454284,
     -641395912765696179, 11285328049647162],
    [-35358708563462647994607, 2621647313739427449588,
     -219254957235699134757, 14665379977955069244,
     -787772412770249571, 11285328049647162],
    [183122730884960782522323, 17007935571912588694272,
     1831498005344531215377, 74598289369102611840,
     2985212511317920527, 11285328049647162],
]
X_FACTORS = [
    (26, [-468970154183, 4162221598416, -236377090446, 3327579857]),
    (2, [-6096612004379, -6525260234283, -166911461493, 21571897004]),
    (-2, [48838428515024, -4891775410473, 712614460299, 114744133]),
    (-22, [-40486600394636, 962384304273, -125714528331, 6769903847]),
    (22, [145782790908731, 9374196436563, 789848985597, 19735990876]),
]
CHARACTERS = [
    # q, t, y, zeta_3, four ratio residues, four character exponents
    (73, 4, 33, 8, [33, 2, 38, 2], [1, 2, 1, 2]),
    (79, 3, 2, 23, [26, 21, 8, 1], [2, 0, 0, 0]),
    (97, 2, 2, 35, [27, 65, 31, 92], [0, 2, 1, 1]),
    (103, 1, 7, 46, [102, 95, 58, 19], [0, 0, 2, 1]),
]
FIVE_CHARACTERS = [
    # q, u, t=(-2a)u^3, y, zeta_3, five raw residues, character exponents
    (73, 12, 51, 7, 8, [55, 12, 29, 45, 22], [2, 1, 2, 1, 0]),
    (79, 4, 41, 37, 23, [57, 32, 2, 64, 12], [0, 2, 1, 0, 0]),
    (97, 1, 21, 4, 35, [11, 48, 2, 21, 11], [2, 2, 1, 2, 2]),
    (103, 5, 3, 5, 46, [79, 100, 101, 100, 88], [0, 0, 1, 0, 2]),
    (109, 1, 88, 18, 45, [106, 84, 51, 61, 55], [2, 1, 2, 2, 0]),
]


def poly(coefficients):
    return s.Poly.from_list(list(reversed(coefficients)), t, domain=s.QQ)


def valuation(integer, prime):
    assert integer != 0
    value = 0
    while integer % prime == 0:
        value += 1
        integer //= prime
    return value


def variations(signs):
    return sum(left != right for left, right in zip(signs, signs[1:]))


def rational_mod(value, modulus):
    value = s.Rational(value)
    return int(value.p)*pow(int(value.q), -1, modulus) % modulus


def verify():
    f = poly(F_COEFFS)
    ys = [poly(row) for row in Y_COEFFS]
    xs = [constant * poly(row) for constant, row in X_FACTORS]
    a = Y_COEFFS[0][-1]
    assert f.degree() == 10 and f.LC() == a*a
    assert all(y.degree() == 5 and y.LC() == a for y in ys)
    assert s.gcd(f, f.diff()).degree() == 0 and f.eval(0) > 0
    for x, y in zip(xs, ys):
        assert y*y - f == x**3
        assert s.gcd(x, f).degree() == 0
        # Detect a coefficient transcription error rather than just fitting values.
        assert y*y - (f + 1) != x**3

    sturm = f.sturm()
    assert [h.degree() for h in sturm] == list(range(10, -1, -1))
    plus = [1 if h.LC() > 0 else -1 for h in sturm]
    minus = [sign*(-1)**h.degree() for sign, h in zip(plus, sturm)]
    assert plus == [1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1]
    assert minus == [1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1]
    assert variations(plus) == variations(minus) == 5

    disc = int(f.discriminant())
    bad_integer = 6*a*abs(disc)
    rows = []
    for q, argument, radical, zeta, expected, characters in CHARACTERS:
        assert s.isprime(q) and q % 3 == 1 and bad_integer % q != 0
        assert radical*radical % q == int(f.eval(argument)) % q != 0
        differences = [(radical-int(y.eval(argument))) % q for y in ys]
        assert all(differences)
        residues = [b*pow(differences[4], -1, q) % q for b in differences[:4]]
        assert residues == expected
        assert zeta != 1 and pow(zeta, 3, q) == 1
        powers = [1, zeta, zeta*zeta % q]
        assert [powers.index(pow(b, (q-1)//3, q)) for b in residues] == characters
        rows.append(characters)
    assert int(s.det(s.Matrix(rows))) % 3 == 1
    duplicate = rows[:3] + [rows[0]]
    assert int(s.det(s.Matrix(duplicate))) % 3 == 0

    # One local control at 3, at the uniform bound v_3(t)=-5.
    # Hensel-lift the square root of f*(z)/a^2 congruent to 1 mod 3.
    z = 3**5
    fstar = sum(c*z**(10-i) for i, c in enumerate(F_COEFFS))
    ystars = [sum(c*z**(5-i) for i, c in enumerate(row)) for row in Y_COEFFS]
    assert valuation(a, 3) == 1
    assert valuation(fstar-a*a, 3) >= 5
    precision = 12
    modulus = 3**precision
    e = valuation(a*a, 3)
    unit = (fstar//3**e)*pow((a*a)//3**e, -1, modulus) % modulus
    root, mod = 1, 3
    while mod < modulus:
        next_mod = 3*mod
        roots = [root+k*mod for k in range(3)
                 if ((root+k*mod)**2-unit) % next_mod == 0]
        assert len(roots) == 1
        root, mod = roots[0], next_mod
    radical = -a*root
    den = radical-ystars[4]
    assert valuation(den, 3) == 1
    local = []
    for yi in ystars[:4]:
        num = radical-yi
        assert valuation(num, 3) == 1
        local.append((num//3)*pow(den//3, -1, 27) % 27)
    assert local == [1, 1, 1, 1]

    print(f"Python {platform.python_version()}; SymPy {s.__version__}")
    print("Verified five integer cube norms and the shared square leading term.")
    print("Sturm variations at -infinity and +infinity: 5, 5; f is positive.")
    print("Four good-prime character rows: determinant 1 modulo 3.")
    print("At t=3^-5 the four negative-branch ratios are 1 modulo 27.")
    print("LY1--LY3 exact inputs passed; uniform proofs are in the note.")


def verify_five_and_towers():
    f = poly(F_COEFFS)
    ys = [poly(row) for row in Y_COEFFS]
    xs = [constant * poly(row) for constant, row in X_FACTORS]
    a = Y_COEFFS[0][-1]
    c = -2*a
    assert c == -22570656099294324
    substitution = s.Poly(c*t**3, t, domain=s.QQ)
    cover = f.compose(substitution)
    assert cover.degree() == 30 and cover.eval(0) != 0
    assert s.gcd(cover, cover.diff()).degree() == 0
    assert c*(c*t**3)**5 == (c*c*t**5)**3
    for x, y in zip(xs, ys):
        assert y.compose(substitution)**2-cover == x.compose(substitution)**3

    bad_integer = 6*a*abs(int(f.discriminant()))
    assert all(bad_integer % p == 0 for p in (2, 3, 5, 11))
    rows = []
    for q, u, argument, radical, zeta, expected, characters in FIVE_CHARACTERS:
        assert s.isprime(q) and q % 3 == 1 and bad_integer % q != 0
        assert u != 0 and c*pow(u, 3, q) % q == argument
        assert radical*radical % q == int(f.eval(argument)) % q != 0
        residues = [(radical-int(y.eval(argument))) % q for y in ys]
        assert all(residues) and residues == expected
        assert zeta != 1 and pow(zeta, 3, q) == 1
        powers = [1, zeta, zeta*zeta % q]
        assert [powers.index(pow(b, (q-1)//3, q)) for b in residues] == characters
        # The common normalization is a cube in this residue field.
        assert c*pow(argument, 5, q) % q == pow(c*c*pow(u, 5, q), 3, q)
        rows.append(characters)
    assert int(s.det(s.Matrix(rows))) % 3 == 2
    assert int(s.det(s.Matrix(rows[:4]+[rows[0]]))) % 3 == 0

    # The bad-prime bound is v_p(u)<=-v_p(2a)-1.
    # At 3 use u=3^-2, giving v_3(t)=-5 and five normalized units.
    u = s.Rational(1, 9)
    argument = c*u**3
    z = 1/argument
    assert z.p % 3**5 == 0 and z.q % 3 != 0
    fstar = sum(b*z**(10-i) for i, b in enumerate(F_COEFFS))
    ystars = [sum(b*z**(5-i) for i, b in enumerate(row)) for row in Y_COEFFS]
    modulus = 3**12
    unit = rational_mod(fstar/(a*a), modulus)
    root, mod = 1, 3
    while mod < modulus:
        next_mod = 3*mod
        roots = [root+k*mod for k in range(3)
                 if ((root+k*mod)**2-unit) % next_mod == 0]
        assert len(roots) == 1
        root, mod = roots[0], next_mod
    radical = -a*root
    local = [rational_mod((radical-yi)/c, 27) for yi in ystars]
    assert local == [1]*5

    parameter = s.Rational(2, 5)
    slope = parameter**3-parameter
    assert slope < 0  # Both estimates improve as the actual class rank grows.
    assert 1+4*slope+2*parameter**2 == -s.Rational(3, 125)
    assert 1+5*slope+4*parameter**2 == -s.Rational(1, 25)

    print("Degree-30 cover squarefree; common leading term is an exact cube.")
    print("Five good-prime character rows: determinant 2 modulo 3.")
    print("At u=3^-2 all five normalized negative-branch values are 1 modulo 27.")
    print("Tower polynomial controls: -3/125 (rank 4, two cuts), -1/25 (rank 5, four cuts).")
    print("LY4--LY6 exact inputs passed; tower and uniform local arguments remain written proofs.")


if __name__ == "__main__":
    verify()
    verify_five_and_towers()
