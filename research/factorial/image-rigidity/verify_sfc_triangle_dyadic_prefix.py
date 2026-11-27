"""Exact residue certificates for the dyadic prefix windows."""

from math import isqrt, lcm

from sage.all import GF, PolynomialRing


def check_all_frobenius_exponents():
    """Prove G_{3(25*2^u+1)}=T^2+T+1 for every u>=0."""
    field = GF(4, "w")
    omega = field.gen()
    assert omega**2 + omega + 1 == 0
    ring = PolynomialRing(field, "z")
    z = ring.gen()
    cube_roots = (ring(1), ring(omega), ring(omega**2))
    period = 480
    for xi in cube_roots:
        for eta in cube_roots:
            polynomial = z**25 * (1 + z)**25 + xi * (1 + z)**25 + eta * z**25
            assert polynomial.degree() == 50
            assert polynomial.gcd(z * (1 + z)) == 1

            image = z
            for _ in range(period):
                image = image**2 % polynomial
            reduced = polynomial.gcd(image - z)
            # Every root of the degree-50 polynomial now lies in F_(2^480).
            assert reduced**2 % polynomial == 0

            initial = (xi * (z**25).inverse_mod(reduced)) % reduced
            root = initial
            for _ in range(period):
                candidates = reduced.gcd(root - z)
                assert (z**2 + z + 1) % candidates == 0
                root = root**2 % reduced
            assert root == initial


def check_gcd_examples():
    binary_ring = PolynomialRing(GF(2), "T")
    x = binary_ring.gen()
    left = sum(x**i for i in range(16) if i != 8)
    right = x**15 + x**8 + 1
    assert left * (x**23 + 1) + right * ((1 + x)**23 + 1) == 1
    assert pow(2, 11, 23) == 1 and 2 % 23 != 1
    assert pow(2, 10, 25) == 24
    assert (x**827 + 1).gcd((1 + x)**827 + 1) == 1
    assert (x**621 + 1).gcd((1 + x)**621 + 1) == x**2 + x + 1


def check_nonprime_power_progression(multiplier, start, period, primes, valuations):
    for shift, (prime, exponent) in enumerate(zip(primes, valuations, strict=True)):
        assert all(prime % d != 0 for d in range(2, isqrt(prime) + 1))
        argument = 2 * multiplier * 2**start - 1 - 2 * shift
        assert argument > prime**exponent
        assert argument % prime**exponent == 0
        assert argument % prime**(exponent + 1) != 0
        assert pow(2, period, prime**(exponent + 1)) == 1


if __name__ == "__main__":
    check_all_frobenius_exponents()
    check_gcd_examples()
    check_nonprime_power_progression(23, 8, 26311740, (3, 61, 79, 3, 7, 5), (1,) * 6)
    primes = (191, 461, 3, 43, 23, 3)
    valuations = (1, 1, 1, 1, 1, 2)
    period = lcm(*((p - 1) * p**e for p, e in zip(primes, valuations, strict=True)))
    assert period == 22932289486260
    check_nonprime_power_progression(25, 14, period, primes, valuations)
    print("Verified both prefix families, the period-480 residue certificate, and the n=826 gcds.")
