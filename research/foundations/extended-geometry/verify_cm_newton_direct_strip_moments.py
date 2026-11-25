"""Replay the fixed-strip Frobenius and moment calculations in SageMath."""

from pathlib import Path
import sys

from sage.all import (
    CyclotomicField,
    GF,
    PolynomialRing,
    QQ,
    ZZ,
    euler_phi,
    lcm,
    matrix,
)


sys.path.insert(0, str(Path(__file__).resolve().parent))

from verify_cm_newton_mu_cyclic_path import verify_distinct_powers


def relative_frobenius(prime, zeta):
    field = GF(prime**2, "a")
    generator = field.multiplicative_generator()
    logs = {}
    value = field.one()
    for index in range(prime**2 - 1):
        logs[value] = index
        value *= generator

    traces = []
    for elements, stride in ((map(field, range(prime)), prime + 1), (iter(field), 1)):
        counts = [0] * 10
        for u in elements:
            value = u**5 * (u - 1)**3 * (u - 3)**3
            if value:
                index = logs[value]
                assert index % stride == 0
                counts[(index // stride) % 10] += 1
        traces.append(-sum(count * zeta**i for i, count in enumerate(counts)))
    trace, second_power_trace = traces
    return trace, (trace**2 - second_power_trace) / 2


def non_cm_replay():
    field = CyclotomicField(10)
    zeta = field.gen()
    relative = PolynomialRing(field, "X")
    x = relative.gen()
    rationals = PolynomialRing(QQ, "X")
    exponent = lcm([j for j in range(1, 8193) if euler_phi(j) <= 64])
    records = (
        (11, 137, zeta**3 + zeta**2 - zeta, -9*zeta**3 + 3*zeta**2 - 7*zeta - 3,
         11**2 * 149**2),
        (31, 103, -zeta**3 + 2*zeta - 3, -17*zeta**3 + 20*zeta**2 - 26*zeta + 38,
         31**2 * 389**2),
    )
    for p, ell, wanted_trace, wanted_det, wanted_norm in records:
        trace, determinant = relative_frobenius(p, zeta)
        assert (trace, determinant) == (wanted_trace, wanted_det)
        assert (trace**2 - 4*determinant).norm() == wanted_norm
        polynomial = relative.one()
        for character in (1, 3, 7, 9):
            conjugate = field.hom([zeta**character], field)
            polynomial *= x**2 - conjugate(trace) * x + conjugate(determinant)
        polynomial = rationals(polynomial)
        assert polynomial.is_irreducible()
        powered = verify_distinct_powers(polynomial, ell, exponent)
        print(f"p={p}: irreducible degree eight, discriminant norm {wanted_norm}, "
              f"separable power polynomial modulo {ell}: {powered}")


def angular_polynomial(base):
    ring = PolynomialRing(base, "z")
    z = ring.gen()
    return 3**8 * (z - 1)**3 * (z - 3)**3


def exact_moment_replay():
    ring = PolynomialRing(ZZ, ("x", "t"))
    x, t = ring.gens()
    shifted_phase = ring.zero()
    for degree, coefficient in enumerate(angular_polynomial(ZZ).list()):
        angle = degree - 1
        shifted_phase += coefficient * x**degree * t**(10 - 2 * angle)
    # x times the lower Laurent term from lifting x^(-1)t^12 by -2yt^10.
    shifted_phase -= angular_polynomial(ZZ)[0] * t**10
    prime = 7
    phases = (
        (1, shifted_phase + 2*x*t**8 + 3*x**2*t**6 + x**5*t),
        (2, shifted_phase**2 + x**2*t**19),
    )
    for scale, phase in phases:
        assert all(j + 2*(i-scale) <= 10*scale and 9*j + 20*(i-scale) >= 70*scale
                   for i, j in phase.dict())
        assert (9*prime**2 - 1) % (10*scale) == 0
        moment = (9*prime**2 - 1) // (10*scale)
        contribution = QQ.zero()
        for (i, j), coefficient in (phase**moment).dict().items():
            if i == scale*moment and j % 2 == 0:
                assert 9*j >= 70*scale*moment and j <= 10*scale*moment
                contribution += QQ(coefficient) / (j + 1)
        normalized = 9*prime**2 * contribution
        assert normalized.denominator() % prime != 0
        assert GF(prime)(normalized) == 4
        print(f"p={prime}, k={scale}: full rational strip moment has normalized residue 4.")


def selected_entries():
    for p, expected in ((7, 4), (17, 4), (37, 9)):
        m = (9*p**2 - 1) // 10
        value = (angular_polynomial(GF(p))**m)[m]
        assert value == expected
        print(f"p={p}, m={m}: selected entry {value}.")


def order_five_span():
    def residues(a, b, c):
        scale = (b-a)**6
        value = QQ(a-c)/(a-b)**2
        derivative = QQ(1)/(a-b)**2 - 2*QQ(a-c)/(a-b)**3
        return scale*derivative, scale*(a*derivative+value)

    rows = [residues(0, 1, 3), residues(0, 3, 1)]
    assert rows == [(-5, -3), (27, -81)]
    assert matrix(QQ, rows).det() == 486
    print("Order-five source classes have pairing determinant 486.")


def main():
    if not __debug__:
        raise SystemExit("Run with assertions enabled.")
    non_cm_replay()
    selected_entries()
    exact_moment_replay()
    order_five_span()


if __name__ == "__main__":
    main()
