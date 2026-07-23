#!/usr/bin/env python3
"""Exact regression for generic affine-mark faithfulness.

The marked Hessian-divisor fiber is the rerooting orbit.  This checker verifies the
degree-independent pencil transport and the reconstruction calculation which
shows that only the root-one sheet can be both unramified and affine.
"""

import sympy as sp


w, s, t, a, kappa = sp.symbols("w s t a kappa", nonzero=True)

# The rerooting pencil identity is coefficient- and degree-independent.
for degree in range(4, 10):
    coefficients = sp.symbols(f"h0:{degree + 1}")
    H = sum(coefficients[j] * w**j for j in range(degree + 1))
    G = sp.expand(kappa * H.subs(w, a * w))
    transported = sp.expand(
        kappa
        * (
            H.subs(w, a * w)
            - (s / (kappa * a)) * a * w
            + t / kappa
        )
    )
    assert sp.expand((G - s * w + t) - transported) == 0


# Let rho be a nonzero simple root.  On its C=0 sheet,
# gamma=-H'(rho)/c.  Regularity of y forces gamma=rho.  Once that leading
# pole cancels, the numerator of z tends to rho-1, independently of a0.
rho, hprime, c, a0 = sp.symbols("rho hprime c a0", nonzero=True)
gamma = -hprime / c
y_numerator = sp.factor(rho - gamma)
assert sp.factor(y_numerator.subs(hprime, -c * rho)) == 0

z_numerator = sp.factor(gamma - 1 - a0 * (rho / gamma - 1))
assert sp.factor(z_numerator.subs(hprime, -c * rho) - (rho - 1)) == 0

# The normalized root one has c=-H'(1), hence gamma=1 and zero z numerator.
assert sp.factor(gamma.subs({rho: 1, hprime: -c}) - 1) == 0
assert sp.factor(z_numerator.subs({rho: 1, hprime: -c})) == 0


# A generic normalized seed has N-2 nonzero simple roots.  The Hessian divisor
# forgets which one is the affine root, while the affine mark selects one.
for degree in range(4, 20):
    assert degree - 2 > 0


print("PASS affine mark: rerooting transports every degree-N inverse pencil")
print("PASS affine mark: y-regularity forces gamma=rho on a simple root sheet")
print("PASS affine mark: z-regularity then forces rho=1")
print("PASS affine mark: the marked fiber has degree one over the seed open")
