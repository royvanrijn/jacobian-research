#!/usr/bin/env python3
"""The explicit Mestre--Fermigier family that produced the rank-22 curve.

Starting with roots ``(0,55,314,378,1007,1036)``, write

    q(X-T) q(X+T) = g(X,T)^2 - r(X,T).

After removing the square factor ``(50616*T)^2``, the genus-one quartic
``Y^2 = R_T(X)`` has the coefficients implemented below.  The parameter in
this normalized tuple model is twice the parameter printed by Fermigier;
thus his ``19754/39`` is ``T=39508/39`` here.  PARI independently reduces
that specialization to Fermigier's published E_22 model.

The associated Jacobian model is computed from the classical binary-quartic
invariants I and J:

    y^2 = x^3 - 27 I(T) x - 27 J(T).

This model is intended for exact local scoring and conductor computations.
The binary-quartic covariant map below also transports the thirteen visible
quartic points to exact rational points on this model.  Their existence is
checked exactly; their Mordell--Weil independence still requires a separate
certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd

from ek_k3 import (
    Q,
    fraction_mod,
    hensel_lift_simple_root,
    legendre_symbol,
    polynomial_derivative,
    polynomial_eval,
    rational_square_root,
)


ROOTS = (0, 55, 314, 378, 1007, 1036)
PUBLISHED_PARAMETER = Q(19754, 39)
NORMALIZED_RECORD_PARAMETER = 2 * PUBLISHED_PARAMETER

# H(T) = disc_X(R_T) / 16, in ascending coefficient order.  Odd
# coefficients vanish.  Keeping this explicit makes Hensel replay
# dependency-free and guards against silently changing the family.
DISCRIMINANT_FACTOR_COEFFICIENTS = (
    178751927593521952355531210961215773383808056363395226527454720000,
    0,
    -18656650379213797389666658813347865821574479769276654052771500,
    0,
    529759946663021294253196844469955118570259147664813128103,
    0,
    -2794585724783814536180462100993070634116448921747400,
    0,
    -74195486391707382079823737680731330188240761295,
    0,
    1197167655825977922608410702757393878691200,
    0,
    -7102859996723716046015855835033267744,
    0,
    19079756640074226760165837126400,
    0,
    -17584994369128502534584064,
    0,
    -11924229808341504000,
    0,
    20185251840000,
)


@dataclass(frozen=True)
class FermigierLocalData:
    prime: int
    residue: int
    good_reduction: bool
    point_count: int
    trace: int
    split_multiplicative: bool | None


@dataclass(frozen=True)
class FermigierPowerRoot:
    prime: int
    exponent: int
    modulus: int
    residue: int
    split_multiplicative: bool


class FermigierMestreFamily:
    """Exact normalized tuple model of Fermigier's one-parameter family."""

    excluded_primes = (2, 3)

    @staticmethod
    def quartic_coefficients(t: Fraction) -> tuple[Fraction, ...]:
        """Return ``(a,b,c,d,e)`` for ``aX^4+bX^3+cX^2+dX+e``."""

        return (
            t**2 + 1149050,
            -30 * (62 * t**2 + 68377393),
            -(2 * t**4 - 1718550 * t**2 - 1195214262641),
            30 * (62 * t**4 - 21690305 * t**2 - 8594794400346),
            t**6
            - 879500 * t**4
            + 102302344648 * t**2
            + 18103855887324900,
        )

    @classmethod
    def binary_invariants(cls, t: Fraction) -> tuple[Fraction, Fraction]:
        a, b, c, d, e = cls.quartic_coefficients(t)
        invariant_i = 12 * a * e - 3 * b * d + c**2
        invariant_j = (
            72 * a * c * e
            + 9 * b * c * d
            - 27 * a * d**2
            - 27 * b**2 * e
            - 2 * c**3
        )
        return invariant_i, invariant_j

    @classmethod
    def coefficients(cls, t: Fraction) -> tuple[Fraction, ...]:
        invariant_i, invariant_j = cls.binary_invariants(t)
        return (Q(0), Q(0), Q(0), -27 * invariant_i, -27 * invariant_j)

    @staticmethod
    def discriminant_factor(t: Fraction) -> Fraction:
        answer = Q(0)
        for coefficient in reversed(DISCRIMINANT_FACTOR_COEFFICIENTS):
            answer = answer * t + coefficient
        return answer

    @classmethod
    def invariants(cls, t: Fraction) -> dict[str, Fraction]:
        invariant_i, invariant_j = cls.binary_invariants(t)
        # For the short Jacobian model.  Its discriminant is a fixed
        # 2,3-unit times H(T), so valuations agree at every p >= 5.
        c4 = 1296 * invariant_i
        c6 = 23328 * invariant_j
        discriminant = 16 * 3**9 * (4 * invariant_i**3 - invariant_j**2)
        return {
            "binary_i": invariant_i,
            "binary_j": invariant_j,
            "c4": c4,
            "c6": c6,
            "quartic_discriminant": 16 * cls.discriminant_factor(t),
            "weierstrass_discriminant": discriminant,
        }

    @classmethod
    def quartic_value(cls, t: Fraction, x: Fraction) -> Fraction:
        answer = Q(0)
        for coefficient in cls.quartic_coefficients(t):
            answer = answer * x + coefficient
        return answer

    @classmethod
    def quartic_covariants_at(
        cls, t: Fraction, x: Fraction
    ) -> tuple[Fraction, Fraction]:
        """Return the binary-quartic covariants ``(g(x,1), h(x,1))``.

        For the homogenization

        ``U=a*X^4+b*X^3*Y+c*X^2*Y^2+d*X*Y^3+e*Y^4``, this uses

        ``g=(U_XY^2-U_XX*U_YY)/144`` and
        ``h=(U_X*g_Y-U_Y*g_X)/8``.

        All divisions are performed in :class:`fractions.Fraction`; no
        floating approximation enters the covariant calculation.
        """

        t = Q(t)
        x = Q(x)
        a, b, c, d, e = cls.quartic_coefficients(t)

        # Coefficients of g in descending binary-quartic order.  Writing
        # these out once lets us evaluate g and its two first derivatives
        # without a symbolic-algebra dependency.
        g0 = b**2 / Q(16) - a * c / Q(6)
        g1 = b * c / Q(12) - a * d / Q(2)
        g2 = c**2 / Q(12) - b * d / Q(8) - a * e
        g3 = c * d / Q(12) - b * e / Q(2)
        g4 = d**2 / Q(16) - c * e / Q(6)

        u_x = 4 * a * x**3 + 3 * b * x**2 + 2 * c * x + d
        u_y = b * x**3 + 2 * c * x**2 + 3 * d * x + 4 * e
        g_value = g0 * x**4 + g1 * x**3 + g2 * x**2 + g3 * x + g4
        g_x = 4 * g0 * x**3 + 3 * g1 * x**2 + 2 * g2 * x + g3
        g_y = g1 * x**3 + 2 * g2 * x**2 + 3 * g3 * x + 4 * g4
        h_value = (u_x * g_y - u_y * g_x) / Q(8)
        return g_value, h_value

    @classmethod
    def quartic_point_to_jacobian(
        cls,
        t: Fraction,
        point: tuple[Fraction, Fraction],
    ) -> tuple[Fraction, Fraction]:
        """Map an affine quartic point to the short Jacobian model.

        For ``point=(x,z)`` on ``z^2=R_t(x)`` with ``z != 0``, the map is

        ``(x,z) -> (36*g(x,1)/z^2, 108*h(x,1)/z^3)``.

        The source and target equations are both checked exactly.  A point
        with ``z=0`` lies outside this affine formula and is rejected rather
        than being represented by a sentinel for the point at infinity.
        """

        t = Q(t)
        x, z = (Q(coordinate) for coordinate in point)
        if z == 0:
            raise ValueError("the affine quartic-to-Jacobian map requires z != 0")
        if z**2 != cls.quartic_value(t, x):
            raise ValueError("the supplied point is not on the quartic")

        g_value, h_value = cls.quartic_covariants_at(t, x)
        jacobian_x = 36 * g_value / z**2
        jacobian_y = 108 * h_value / z**3
        _, _, _, coefficient_a, coefficient_b = cls.coefficients(t)
        expected_square = (
            jacobian_x**3 + coefficient_a * jacobian_x + coefficient_b
        )
        if jacobian_y**2 != expected_square:
            raise AssertionError("the binary-quartic covariant identity failed")
        return jacobian_x, jacobian_y

    @staticmethod
    def square_approximant(t: Fraction, x: Fraction) -> Fraction:
        """Return g(X,T) for the unscaled identity q=g^2-r."""

        return (
            x**6
            - 2790 * x**5
            + (2726125 - 3 * t**2) * x**4
            + (5580 * t**2 - 1106081640) * x**3
            + (3 * t**4 - 3892050 * t**2 + 176868664084) * x**2
            + (-2790 * t**4 + 1034713080 * t**2 - 6810411651120) * x
            - t**6
            + 1165925 * t**4
            - 128370083212 * t**2
        )

    @classmethod
    def visible_quartic_points(
        cls, t: Fraction
    ) -> tuple[tuple[Fraction, Fraction], ...]:
        """Return the twelve points coming from the roots of q.

        They are points on the normalized quartic.  Their existence does not
        by itself certify their specialized Mordell--Weil independence.
        """

        if t == 0:
            raise ValueError("the normalized visible-point formula excludes T=0")
        points: list[tuple[Fraction, Fraction]] = []
        for root in ROOTS:
            for sign in (-1, 1):
                x = Q(root) + sign * t
                y = cls.square_approximant(t, x) / (50616 * t)
                if y**2 != cls.quartic_value(t, x):
                    raise AssertionError("Mestre identity failed at a visible point")
                points.append((x, y))
        return tuple(points)

    @classmethod
    def extra_quartic_point(cls, t: Fraction) -> tuple[Fraction, Fraction]:
        """Return Mestre's thirteenth point for the ``u=3,v=5`` member.

        The affine change carrying Fermigier's ``u=3,v=5`` roots to
        ``ROOTS`` sends the published abscissa ``A+B*t`` to
        ``1256/5 - (17/35)T``.
        """

        x = Q(1256, 5) - Q(17, 35) * t
        y = rational_square_root(cls.quartic_value(t, x))
        if y is None:
            raise AssertionError("Mestre's extra section failed the curve equation")
        return x, y

    @classmethod
    def known_quartic_points(
        cls, t: Fraction
    ) -> tuple[tuple[Fraction, Fraction], ...]:
        """Return all thirteen visible points used by the generic-rank construction."""

        return cls.visible_quartic_points(t) + (cls.extra_quartic_point(t),)

    @classmethod
    def known_jacobian_points(
        cls, t: Fraction
    ) -> tuple[tuple[Fraction, Fraction], ...]:
        """Return exact Jacobian images of all thirteen known quartic points.

        Each image is checked against ``y^2=x^3-27*I(t)*x-27*J(t)`` by
        :meth:`quartic_point_to_jacobian`.  This is a point certificate, not
        an independence or specialized-rank certificate.
        """

        return tuple(
            cls.quartic_point_to_jacobian(t, point)
            for point in cls.known_quartic_points(t)
        )

    @staticmethod
    def _check_prime(prime: int) -> None:
        if prime < 5:
            raise ValueError("the normalized local model is restricted to p >= 5")

    @classmethod
    def _short_coefficients_mod(cls, residue: int, prime: int) -> tuple[int, int]:
        cls._check_prime(prime)
        invariant_i, invariant_j = cls.binary_invariants(Q(residue))
        return (
            fraction_mod(-27 * invariant_i, prime),
            fraction_mod(-27 * invariant_j, prime),
        )

    @classmethod
    def local_data(cls, residue: int, prime: int) -> FermigierLocalData:
        """Count points on the short Jacobian model over F_p, including bad fibers."""

        cls._check_prime(prime)
        residue %= prime
        coefficient_a, coefficient_b = cls._short_coefficients_mod(residue, prime)
        character_sum = 0
        for x in range(prime):
            rhs = (x**3 + coefficient_a * x + coefficient_b) % prime
            character_sum += legendre_symbol(rhs, prime)
        point_count = prime + 1 + character_sum
        trace = -character_sum
        discriminant_zero = cls.discriminant_factor(Q(residue)) % prime == 0
        split = None
        if discriminant_zero:
            # A nodal cubic with c4 != 0 has a_p = +1 (split) or -1
            # (non-split).  The point count distinguishes the two cases.
            invariant_i, _ = cls.binary_invariants(Q(residue))
            if fraction_mod(invariant_i, prime) != 0 and trace in (-1, 1):
                split = trace == 1
        return FermigierLocalData(
            prime,
            residue,
            not discriminant_zero,
            point_count,
            trace,
            split,
        )

    @classmethod
    def power_roots(
        cls, prime: int, exponent: int, *, split_only: bool = False
    ) -> tuple[FermigierPowerRoot, ...]:
        """Find and Hensel-lift clean roots of the degree-20 discriminant factor."""

        cls._check_prime(prime)
        if exponent < 1:
            raise ValueError("the lifting exponent must be positive")
        derivative = polynomial_derivative(DISCRIMINANT_FACTOR_COEFFICIENTS)
        roots: list[FermigierPowerRoot] = []
        for residue in range(prime):
            if polynomial_eval(DISCRIMINANT_FACTOR_COEFFICIENTS, residue, prime) != 0:
                continue
            if polynomial_eval(derivative, residue, prime) == 0:
                continue
            local = cls.local_data(residue, prime)
            if local.split_multiplicative is None:
                continue
            if split_only and not local.split_multiplicative:
                continue
            lifted = hensel_lift_simple_root(
                DISCRIMINANT_FACTOR_COEFFICIENTS,
                residue,
                prime,
                exponent,
            )
            roots.append(
                FermigierPowerRoot(
                    prime,
                    exponent,
                    prime**exponent,
                    lifted,
                    bool(local.split_multiplicative),
                )
            )
        return tuple(roots)

    @classmethod
    def verify_power_constraint(
        cls, numerator: int, denominator: int, prime: int, exponent: int
    ) -> int:
        if gcd(denominator, prime) != 1:
            raise ValueError("the rational denominator is not a p-adic unit")
        value = cls.discriminant_factor(Q(numerator, denominator))
        if value == 0:
            raise ValueError("the specialization is singular")

        def integer_valuation(integer: int) -> int:
            integer = abs(integer)
            result = 0
            while integer % prime == 0:
                integer //= prime
                result += 1
            return result

        actual = integer_valuation(value.numerator) - integer_valuation(value.denominator)
        if actual < exponent:
            raise AssertionError(
                f"expected v_{prime}(H(T)) >= {exponent}, obtained {actual}"
            )
        return actual


def is_rational_square(value: Fraction) -> bool:
    return rational_square_root(value) is not None
