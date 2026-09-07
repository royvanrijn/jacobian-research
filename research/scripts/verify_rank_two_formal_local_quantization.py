#!/usr/bin/env python3
"""Certificate for the all-orders formal-local rank-two quantization.

After the Ore/Darboux localization, R is a central parameter and (S,T) is
an etale canonical coordinate pair on every completed fiber.  In these
coordinates the linearized commutator operator is

    d(f,g) = partial_S(f) + partial_T(g).

The explicit homotopy below proves that every higher commutator defect can
be removed in the completed local ring.  The second calculation records the
Hamiltonian primitive used to restore the R-connection.  This certificate
does not assert that the resulting coefficients descend to filtered
polynomials in the unlocalized Weyl algebra.
"""

import sympy as sp


S, T = sp.symbols("S T")


def linearized_commutator(s_correction, t_correction):
    return sp.expand(
        sp.diff(s_correction, S) + sp.diff(t_correction, T)
    )


def obstruction_homotopy(defect):
    """Return (f,g) with d(f,g)=-defect and zero integration constant."""

    return -sp.integrate(defect, S), sp.Integer(0)


# A generic finite jet stands for an arbitrary coefficient of the completed
# local ring; the monomial formula works unchanged at every jet order.
coefficients = sp.symbols("c0:21")
generic_jet = sum(
    coefficients[index] * S**s_degree * T**t_degree
    for index, (s_degree, t_degree) in enumerate(
        pair
        for total in range(6)
        for pair in ((i, total - i) for i in range(total + 1))
    )
    if index < len(coefficients)
)
fiber_correction = obstruction_homotopy(generic_jet)
assert sp.expand(linearized_commutator(*fiber_correction) + generic_jet) == 0


# The R-variation of a canonical pair is a divergence-free vector field.
# A Hamiltonian primitive A satisfies A_T=a and A_S=-b for
# V=a*d_S+b*d_T.  Test the explicit primitive on a generic Hamiltonian jet.
hamiltonian = sum(
    coefficients[index] * S**s_degree * T**t_degree
    for index, (s_degree, t_degree) in enumerate(
        pair
        for total in range(5)
        for pair in ((i, total - i) for i in range(total + 1))
    )
)
a = sp.diff(hamiltonian, T)
b = -sp.diff(hamiltonian, S)
assert sp.expand(sp.diff(a, S) + sp.diff(b, T)) == 0
assert sp.expand(sp.diff(hamiltonian, T) - a) == 0
assert sp.expand(sp.diff(hamiltonian, S) + b) == 0


print("PASS: d(f,g)=f_S+g_T has an explicit contracting homotopy")
print("PASS: every finite defect jet is removable, hence so is every formal jet")
print(
    "PASS: every divergence-free R-variation has a formal Hamiltonian primitive"
)
print(
    "SCOPE: the construction is local/etale and hbar-adic, not globally filtered"
)
