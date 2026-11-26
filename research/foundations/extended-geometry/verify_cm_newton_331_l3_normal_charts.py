"""Conductor-18 weighted boundary: exclude v*w != 0 on L3."""
from math import comb

from sage.all import GF, QQ, PolynomialRing, matrix, vector

from verify_cm_newton_full_width_reconstruction_certificates import reduced_residual_rows


def homogeneous(polynomial, degree):
    return polynomial.parent()({e: c for e, c in polynomial.dict().items() if sum(e) == degree})


def l3_jets(rows):
    parameter = PolynomialRing(QQ, "w")
    w = parameter.gen()
    field = parameter.fraction_field()
    ring = PolynomialRing(field, names=("L", "H", "tau"), order="degrevlex")
    l, h, tau = ring.gens()
    order = 7
    expansions = []
    for k, row in enumerate(rows, 7):
        blocks = [ring.zero() for _ in range(order + 1)]
        for (r, ell, h_power, v, w_power), coefficient in row.dict().items():
            tau_power, remainder = divmod(k - 1 - (2*r + ell + 2*h_power + v + 2*w_power), 2)
            assert not remainder and tau_power >= 0
            if r + tau_power > order:
                continue
            for a in range(min(ell, order - r - tau_power) + 1):
                for b in range(min(h_power, order - r - tau_power - a) + 1):
                    blocks[r] += (
                        coefficient * comb(ell, a) * comb(h_power, b) * 8**(ell - a)
                        * 17**(h_power - b) * w**(w_power + h_power - b)
                        * l**a * h**b * tau**tau_power
                    )
        expansions.append(blocks)

    def truncate(f, degree):
        return ring({e: c for e, c in f.dict().items() if sum(e) <= degree})

    def evaluate(blocks, r, degree):
        result = ring.zero()
        power = ring.one()
        for block in blocks:
            result += truncate(block * power, degree)
            power = truncate(power * r, degree)
            if not power:
                break
        return result

    pivot = expansions[0][1].constant_coefficient()
    assert pivot == -QQ(12829696)/194828780125125*w**2
    solution = ring.zero()
    for degree in range(1, order + 1):
        correction = -evaluate(expansions[0], solution, degree)/pivot
        assert correction == homogeneous(correction, degree)
        solution += correction
    residuals = [evaluate(blocks, solution, order) for blocks in expansions]
    assert not residuals[0] and all(not truncate(f, 2) for f in residuals)
    return field, residuals[1:]


def verify_l3_nonendpoint_chart(rows):
    exact, jets = l3_jets(rows)
    w = exact.gen()
    reduction_ring = PolynomialRing(GF(331), "w")
    reduced = reduction_ring.fraction_field()
    w_mod = reduced.gen()
    finite_chart = PolynomialRing(GF(331), names=("W", "A", "U", "V"), order="degrevlex")
    W, A, U, V = finite_chart.gens()

    def modulo(coefficient):
        coefficient = exact(coefficient)
        numerator = reduction_ring(coefficient.numerator())
        denominator = reduction_ring(coefficient.denominator())
        assert denominator and denominator == denominator.leading_coefficient() * reduction_ring.gen()**denominator.degree()
        return reduced(numerator) / reduced(denominator)

    def modulo_polynomial(f, target):
        return target({e: modulo(c) for e, c in f.dict().items()})

    def chart_substitution(f, values):
        answer = finite_chart.zero()
        for exponent, coefficient in f.dict().items():
            c = modulo(coefficient)
            denominator = c.denominator()
            term = sum(value*W**i for i, value in c.numerator().dict().items())
            term *= A**denominator.degree()/denominator.leading_coefficient()
            for value, power in zip(values, exponent, strict=True):
                term *= value**power
            answer += term
        return answer

    def membership(rows, values, identities):
        generators = [chart_substitution(f, values) for f in rows] + [W*A-1]
        ideal = finite_chart.ideal(generators)
        for identity in identities:
            certificate = identity.lift(ideal)
            assert sum(p*c for p, c in zip(generators, certificate, strict=True)) == identity

    def exactly_lift(leading, expected_modulus, expected_first, substitute):
        # Solve the characteristic-zero initial equations, then check all
        # residual identities in the resulting finite algebra.
        ring = PolynomialRing(exact, names=("Y", "Z"), order="lex")
        y, z = ring.gens()
        mapped = [f.parent().hom(substitute(f, y, z), ring)(f) for f in leading]
        linear, univariate = list(ring.ideal(mapped).groebner_basis())
        assert linear.degree(y) == 1 and linear.monomial_coefficient(y) == 1
        assert univariate.degree(y) == 0
        z_ring = PolynomialRing(exact, "z")
        z_map = ring.hom([0, z_ring.gen()], z_ring)
        modulus = z_map(univariate)
        assert modulus.is_monic()
        assert modulo_polynomial(modulus, expected_modulus.parent()) == expected_modulus
        algebra = z_ring.quotient(modulus, "z")
        root = algebra.gen()
        first = -ring.hom([0, root], algebra)(linear)
        evaluate = ring.hom([first, root], algebra)
        assert all(not evaluate(f) for f in mapped)
        polynomial = modulo_polynomial(first.lift(), expected_modulus.parent())
        assert polynomial == expected_first

    def obstruction(leading, first_jets, modulus, first, values, pivot_norm, next_norm):
        ring = PolynomialRing(reduced, names=("x", "z"))
        x, z = ring.gens()
        algebra = modulus.parent().quotient(modulus, "z")
        root = algebra.gen()
        evaluate = ring.hom(values(first(root), root), algebra)
        def dehomogenize(f):
            return sum(modulo(c)*x**a*z**b for (a,b,_),c in f.dict().items())
        maps = [dehomogenize(f) for f in leading]
        assert all(not evaluate(g) for g in maps)
        pivot = matrix(algebra, [
            [evaluate(g.derivative(var)) for var in ring.gens()] for g in maps[:2]
        ])
        norm = lambda value: modulus.resultant(value.lift())
        assert norm(pivot.det()) == pivot_norm
        mapped_jet = [dehomogenize(f) for f in first_jets]
        correction = pivot.inverse() * vector(algebra, [-evaluate(f) for f in mapped_jet[:2]])
        extra = evaluate(mapped_jet[2]) + sum(
            evaluate(maps[2].derivative(var))*delta
            for var, delta in zip(ring.gens(), correction, strict=True)
        )
        assert norm(extra) == next_norm
        # A monomial determinant makes the inverse of the pivot regular at
        # every specialization with w a unit. Confirm regularity directly.
        for f in [extra, *correction]:
            for c in f.lift():
                den = c.denominator()
                assert den == den.leading_coefficient()*reduction_ring.gen()**den.degree()

    # First normal chart. The displayed constraints also identify the
    # only possible direction at its boundary tau=0.
    cubics = [homogeneous(f, 3) for f in jets]
    membership(cubics, (finite_chart.one(), finite_chart.zero(), finite_chart.zero()),
               [finite_chart.one()])
    membership(cubics, (U, finite_chart.one(), finite_chart.zero()),
               [(U-GF(331)(8)**-1)**2])
    ring = PolynomialRing(reduced, "z")
    z = ring.gen()
    cubic = z**3-136*z**2-27*z+75
    cubic_y = 98*z**2 + (8/w_mod-101)*z + 62
    assert cubic.is_squarefree()
    membership(cubics, (A*U, V, finite_chart.one()), [
        U**3-136*U**2-27*U+75,
        V-98*U**2-(8*A-101)*U-62,
    ])
    exactly_lift(cubics, cubic, cubic_y, lambda f, y, z: [z/w, y, 1])
    obstruction(cubics[:3], [homogeneous(f, 4) for f in jets[:3]],
                cubic, cubic_y, lambda y, z: [z/w_mod, y], 191*w_mod**6, reduced(68))

    # Second normal chart. Set H=s, L=s(1/8+X), tau=sZ. Dividing by
    # s^3 and cancelling the common quadratic leaves cubics.
    chart = PolynomialRing(exact, names=("X", "Z", "s"), order="degrevlex")
    x, z_old, s = chart.gens()
    second = []
    for f in jets:
        g = sum(c*(QQ(1)/8+x)**a*z_old**d*s**(a+b+d-3) for (a,b,d), c in f.dict().items())
        assert not g.constant_coefficient() and not homogeneous(g, 1)
        second.append(g)
    pivot = second[0]
    quadratic = homogeneous(pivot, 2)
    head = quadratic.monomial_coefficient(x**2)
    assert head == -QQ(41156608)/139765022353299915*w
    cancelled = []
    for f in second[1:]:
        coefficient = homogeneous(f, 2).monomial_coefficient(x**2) / head
        g = f-coefficient*pivot
        assert not homogeneous(g, 2)
        cancelled.append(g)

    leading = [quadratic, *(homogeneous(f, 3) for f in cancelled)]
    for values in [(finite_chart.one(), U, finite_chart.zero()),
                   (U, finite_chart.one(), finite_chart.zero())]:
        membership(leading, values, [finite_chart.one()])
    sextic = z**6-15*z**5-45*z**4-110*z**3+74*z**2+124*z-76
    sextic_x = -41*z**5+74*z**4+18*z**3+140*z**2-81*z+32
    assert sextic.is_squarefree()
    membership(leading, (U, V, finite_chart.one()), [
        V**6-15*V**5-45*V**4-110*V**3+74*V**2+124*V-76,
        U+41*V**5-74*V**4-18*V**3-140*V**2+81*V-32,
    ])
    exactly_lift(leading, sextic, sextic_x, lambda f, y, z: [y, z, 1])
    first_jets = [homogeneous(pivot, 3), *(homogeneous(f, 4) for f in cancelled)]
    obstruction(leading[:3], first_jets[:3], sextic, sextic_x, lambda x, z: [x, z],
                245*w_mod**12, 158*w_mod**6)
    print("At 331: both normal charts exclude the L3 locus with v and w nonzero.")


if __name__ == "__main__":
    verify_l3_nonendpoint_chart(reduced_residual_rows(5, 16, even_linear=True, primitive=False))

