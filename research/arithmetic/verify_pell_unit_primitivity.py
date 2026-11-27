"""Exact polynomial controls for PF1–PF4."""

from sage.all import QQ, PolynomialRing

from verify_quartic_pell_two_class_directions import P, t


def dickson(variable, sign, degree):
    previous, current = variable.parent()(2), variable
    for _ in range(1, degree):
        previous, current = current, variable*current-sign*previous
    return current


def verify_dickson_trace():
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    fraction = ring.fraction_field()
    u = fraction(x)
    for sign in (-1, 1):
        for degree in (2, 3, 5, 7, 11):
            polynomial = dickson(x, sign, degree)
            substitution = u+sign/u
            assert polynomial(substitution) == u**degree+(sign/u)**degree
            derivative = polynomial.derivative()(substitution)
            assert derivative*(1-sign/u**2) == degree*(u**degree-(sign/u)**degree)/u
            critical_value_power = 4*sign**degree
            assert (polynomial**2-critical_value_power) % polynomial.derivative() == 0


def verify_critical_trace_shifts():
    ring = PolynomialRing(QQ, "x")
    x = ring.gen()
    for degree in (3, 5, 7, 11):
        for sign in (-1, 1):
            critical_fiber = dickson(x, 1, degree)-2*sign
            factors = dict(critical_fiber.factor())
            assert factors.pop(x-2*sign) == 1
            assert sum(f.degree() for f in factors) == (degree-1)//2
            assert set(factors.values()) == {2}
    for degree in (2, 4):
        lower_trace = dickson(x, 1, degree)
        assert lower_trace.is_squarefree()
        assert dickson(x, 1, 2*degree)+2 == lower_trace**2


def verify_trace_representations():
    a, b, c = P[4], P[2], P[0]
    assert 2*P == 2*a*(t*t+b/(2*a))**2+2*c-b*b/(2*a)

    ring = PolynomialRing(QQ, names=("B", "c", "power"))
    B, c, power = ring.gens()
    fraction = ring.fraction_field()
    denominator = fraction(B*c-1)
    n = (power-c*c)/(2*denominator)
    offset = 2-B*B*c*c/denominator
    assert 2*(1+B*B*n) == B*B*power/denominator+offset
    assert offset-2 == -(B*c)**2/denominator
    assert offset+2 == -(B*c-2)**2/denominator


if __name__ == "__main__":
    verify_dickson_trace()
    verify_critical_trace_shifts()
    verify_trace_representations()
    print("PF1–PF4 exact trace controls passed")
