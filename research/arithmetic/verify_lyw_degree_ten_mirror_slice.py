"""Exact inputs for LY1--LY7; Python 3.12 and SymPy 1.14.0.

The integer coefficients are transcribed from Levin--Yan--Wiljanen,
arXiv:1910.12276v1, Theorem 9, printed p. 12. No numerical root tests.
The uniform valuation and class-field arguments are in the companion note;
this checker verifies their polynomial, Sturm and residue inputs. The optional
unit probe searches bounded Pell prefixes on three full-slice LY4 fibres.
"""

import argparse
import hashlib
import json
import math
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
    assert all(pow(-3, (q-1)//2, q) == q-1 for q in (5, 11))
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


def unit_residues(radicand, numerator, coefficient, denominator=1):
    """Residues at both split places, for an already certified integral unit.

    The caller must certify integrality and norm +/-1. A noncube residue is
    sufficient even for a unit of a nonmaximal order; cube residues for such a
    unit do not determine the fundamental unit's localization.
    """
    if denominator % 3 == 0:
        raise ValueError("the unit's displayed denominator must be prime to three")
    roots = [r for r in range(9) if r*r % 9 == radicand % 9]
    if len(roots) != 2 or radicand % 3 != 1:
        raise ValueError("three must split in the displayed quadratic algebra")
    values = [(numerator+coefficient*r)*pow(denominator, -1, 9) % 9
              for r in roots]
    if any(value % 3 == 0 for value in values):
        raise ValueError("the displayed element is not a local unit")
    return values


def pell_prefix(radicand, steps):
    """First `steps` convergents to sqrt(radicand), using exact complete quotients.

    Q_(j+1)=1 certifies a Pell unit with norm (-1)^(j+1). Only the
    convergent residues modulo 9 are needed to test its local cube class.
    Exhausting the budget says nothing about maximal-order units.
    """
    if steps < 1:
        raise ValueError("steps must be positive")
    a0 = math.isqrt(radicand)
    if a0*a0 == radicand:
        raise ValueError("a quadratic nonsquare radicand is required")
    p, q, a = 0, 1, a0
    h0, h1, k0, k1 = 0, 1, 1, 0
    for j in range(steps):
        h0, h1 = h1, (a*h1+h0) % 9
        k0, k1 = k1, (a*k1+k0) % 9
        p_next = a*q-p
        quotient, remainder = divmod(radicand-p_next*p_next, q)
        assert remainder == 0 and quotient > 0
        p, q = p_next, quotient
        if q == 1:
            residues = unit_residues(radicand, h1, k1) if radicand % 3 == 1 else None
            return {"steps": j+1, "pell_unit_found": True,
                    "norm": (-1)**(j+1), "unit_residues_mod9": residues,
                    "unit_noncube_at_three": None if residues is None else
                        any(value not in (1, 8) for value in residues)}
        a = (a0+p)//q
    return {"steps": steps, "pell_unit_found": False,
            "unit_noncube_at_three": None}


def verify_unit_criterion():
    # Modulo 9 the cubes of 3-adic units are exactly +/-1.
    assert {pow(r, 3, 9) for r in range(9) if r % 3} == {1, 8}
    # 3+sqrt(10) is an integral norm -1 unit and a noncube at both places.
    assert 3*3-10 == -1
    assert unit_residues(10, 3, 1) == [4, 2]
    assert pell_prefix(10, 1)["unit_noncube_at_three"] is True
    # Its cube 117+37sqrt(10) is locally cubic. A negative test on an
    # arbitrary unit therefore cannot be used to exclude separation.
    assert 117*117-10*37*37 == -1
    assert unit_residues(10, 117, 37) == [1, 8]
    # Half-integral units must be supported too: (3+sqrt(13))/2 has norm -1.
    assert 3*3-13 == -4
    assert unit_residues(13, 3, 1, 2) == [7, 5]
    assert pell_prefix(10, 1)["norm"] == -1
    assert pell_prefix(3, 2)["norm"] == 1
    assert pell_prefix(13, 1)["pell_unit_found"] is False
    for arguments in ((3, 2, 1, 1), (10, 3, 1, 3), (10, 1, 1, 1)):
        try:
            unit_residues(*arguments)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid local-unit data were accepted")
    print("Unit controls: a noncube, its locally cubic cube, and a half-integral unit.")


def probe_ly4_units(steps=2000, fibres=3):
    """A bounded search on the full LY4 slice, not just the local u=1/9 control."""
    if fibres < 1 or steps < 1:
        raise ValueError("fibres and steps must be positive")
    a = Y_COEFFS[0][-1]
    c = -2*a
    bad_integer = 6*a*abs(int(poly(F_COEFFS).discriminant()))
    factors = s.factorint(bad_integer, limit=10000)
    # Do not silently treat a remaining composite cofactor as one bad prime.
    assert all(s.isprime(p) for p in factors)
    assert math.prod(int(p)**int(e) for p, e in factors.items()) == bad_integer
    bad_primes = sorted(int(p) for p in factors)
    denominator = math.prod(p**(valuation(2*a, p)+1) for p in bad_primes)
    from sympy.ntheory.modular import crt
    primes = [row[0] for row in FIVE_CHARACTERS]
    numerator0, modulus = map(int, crt(primes, [row[1]*denominator % row[0]
                                               for row in FIVE_CHARACTERS]))
    records, offset = [], 0
    while len(records) < fibres:
        numerator = numerator0+offset*modulus
        offset += 1
        if math.gcd(numerator, denominator) != 1:
            continue
        assert all(numerator % p != 0 and
                   valuation(denominator, p) >= valuation(2*a, p)+1
                   for p in bad_primes)
        assert all(numerator*pow(denominator, -1, row[0]) % row[0] == row[1]
                   for row in FIVE_CHARACTERS)
        top, bottom = c*numerator**3, denominator**3
        radicand = sum(b*top**i*bottom**(10-i) for i, b in enumerate(F_COEFFS))
        assert radicand > 0
        # Remove certified square factors only. No squarefree/fundamental
        # discriminant claim is made for the remaining radicand.
        for p in bad_primes:
            e = valuation(radicand, p)
            radicand //= p**(2*(e//2))
        assert radicand % 3 == 1
        witness = next((int(p) for p in s.primerange(5, 1000)
                        if pow(radicand % int(p), (int(p)-1)//2, int(p)) == int(p)-1), None)
        assert witness is not None  # Exact certificate that this fibre is quadratic.
        digits = str(radicand)
        record = {"u_numerator": numerator, "u_denominator": denominator,
                  "radicand_digits": len(digits),
                  "radicand_sha256": hashlib.sha256(digits.encode("ascii")).hexdigest(),
                  "nonsquare_witness_prime": witness,
                  "full_LY4_slice": True,
                  **pell_prefix(radicand, steps)}
        records.append(record)
    print(json.dumps({"python": platform.python_version(), "sympy": s.__version__,
                      "bad_primes": bad_primes, "fibres": records}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-probe", action="store_true")
    parser.add_argument("--steps", type=int, default=2000)
    parser.add_argument("--fibres", type=int, default=3)
    args = parser.parse_args()
    if args.unit_probe:
        probe_ly4_units(args.steps, args.fibres)
    else:
        verify()
        verify_five_and_towers()
        verify_unit_criterion()
