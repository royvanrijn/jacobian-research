"""Complete H2=z*t, H1=ell_z*z family; cubic coefficients unrestricted.

This is the intrinsic ker(B) linear-vector boundary, including H1=0.
The exact input is selected from the full linear-coefficient moment
generator by its uniquely recoverable power of the linear w coefficient.
No purity classification is asserted by this generator.
"""
from functools import lru_cache
import math

import sympy as S

from nilpotent_cyclic_slice import CHARGES, CUBIC, LINEAR
from nilpotent_cyclic_slice import exact_moment as cyclic_moment

VARIABLES = (LINEAR[0], *CUBIC)


@lru_cache(maxsize=None)
def exact_moment(order):
    result = {}
    for powers, coefficient in cyclic_moment(order).items():
        if powers[1]:
            continue
        linear_z = powers[0]
        cubic_count = sum(powers[2:])
        charge = sum(j*k for j, k in zip(CHARGES, powers[2:]))
        numerator = order-cubic_count-2*linear_z-charge
        assert numerator % 2 == 0
        quadratic_count = numerator//2
        linear_w = charge+quadratic_count+linear_z
        assert quadratic_count >= 0 and linear_w >= 0
        if linear_w:
            continue
        target = (linear_z, *powers[2:])
        assert target not in result
        result[target] = coefficient
    return result


def raw_expression(order):
    return S.expand(sum(S.Rational(c.numerator, c.denominator)
                        *math.prod(v**k for v, k in zip(VARIABLES, powers))
                        for powers, c in exact_moment(order).items()))


def primitive_expression(order):
    return S.Poly(raw_expression(order), *VARIABLES).clear_denoms()[1].primitive()[1].as_expr()
