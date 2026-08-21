#!/usr/bin/env sage -python
"""Certify the resolved E7 local flag for the first H92 q=6 pencil.

Near the III* fiber ``t=0`` the short H92 equation is

    y^2 = x^3 + (A1*t^3+A*t^4)*x + B1*t^5+B*t^6+B2*t^7.

Putting ``x=x0(t)+U`` along the formal branch with ``y=U=0`` gives the
semiquasihomogeneous E7 coordinates ``(Z,U,Y)=(t,x-x0,y)``.  Through the
length-six quotient needed below,

    x0 = c2*Z^2+c3*Z^3,
    c2=-B1/A1,
    c3=-(c2^3+A*c2+B)/A1.

For the canonical E7 singularity ``Y^2=U^3+U*Z^3``, the relevant complete
ideal is

    I=(U^2,Z*U,Z^4,Z*Y,U*Y).

Its quotient basis is ``1,Z,Z^2,Z^3,U,Y``.  Extracting these six residue
coefficients gives the resolved III* functionals for a numerator already
written in the common local trivialization.  For a global row with shifted
coefficient bounds ``(C,C+2,C+3)``, the common polynomial lift is

    Z^3*rev(N0,C) + Z*rev(Nx,C+2)*x + rev(Ny,C+3)*y.

It must be tested against the shifted cycle ``c+(C+3)*vZ``.  Applying the
three reversals directly to ``1,x,y`` would mix different powers of ``Z``
and is not a section of the same line bundle.  This script certifies both the
length-six flag and a representative shifted complete ideal; it does not
claim that the global two-dimensional q=6 pencil has been recovered.
"""

from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import Matrix, PolynomialRing, QQ, ZZ, vector


ROOT = Path(__file__).resolve().parents[2]
ANCHOR = ROOT / "elkies-k3/scripts/verify_h3_noncm_q6_source_anchor.sage"
H92 = ROOT / "artifacts/local/humbert-inputs/92/igusa92.txt"

anchor = SourceFileLoader("h92_e7_flag_anchor", str(ANCHOR)).load_module()
h92_ring, h92_formulas = anchor.parse_h92(H92)
r92, s92 = anchor.EXPECTED_H92
A1, A, B1, B, B2 = tuple(QQ(value(r92, s92)) for value in h92_formulas)
assert A1

# The E7 graph order is the one returned by the minimal resolution of
# Y^2+U^3+U*Z^3.  The positive Cartan matrix has edges
# 1-4, 2-5, 2-7, 3-4, 3-6, 3-7.
e7_cartan = Matrix(
    QQ,
    [
        [2, 0, 0, -1, 0, 0, 0],
        [0, 2, 0, 0, -1, 0, -1],
        [0, 0, 2, -1, 0, -1, -1],
        [-1, 0, -1, 2, 0, 0, 0],
        [0, -1, 0, 0, 2, 0, 0],
        [0, 0, -1, 0, 0, 2, 0],
        [0, -1, -1, 0, 0, 0, 2],
    ],
)
assert e7_cartan.det() == 2

vZ = vector(QQ, (2, 2, 4, 3, 1, 2, 3))
vU = vector(QQ, (2, 4, 6, 4, 3, 3, 5))
vY = vector(QQ, (3, 5, 9, 6, 3, 5, 7))
e1 = vector(QQ, (1, 0, 0, 0, 0, 0, 0))
e5 = vector(QQ, (0, 0, 0, 0, 1, 0, 0))
e6 = vector(QQ, (0, 0, 0, 0, 0, 1, 0))
assert e7_cartan * vZ == e1
assert e7_cartan * vU == 2 * e5
assert e7_cartan * vY == e5 + e6

# The cycle has boundary equal to the affine branch plus twice the branch
# meeting the minuscule E7 component.  These five monomials generate its
# integrally closed local ideal.
cycle = vZ + vU
generator_valuations = {
    "U2": 2 * vU,
    "ZU": vZ + vU,
    "Z4": 4 * vZ,
    "ZY": vZ + vY,
    "UY": vU + vY,
}
for valuation in generator_valuations.values():
    assert min(valuation - cycle) >= 0
assert e7_cartan * cycle == e1 + 2 * e5
assert cycle * e7_cartan * cycle == 12

# Work directly in the six-dimensional quotient.  Y^2 is redundant modulo
# the surface equation and the displayed generators, but including it makes
# the finite quotient computation transparent.
local_ring = PolynomialRing(QQ, names=("Z", "U", "Y"), order="degrevlex")
Z, U, Y = local_ring.gens()
complete_ideal = local_ring.ideal(
    (U**2, Z * U, Z**4, Z * Y, U * Y, Y**2)
)
assert complete_ideal.vector_space_dimension() == 6
quotient_basis = (local_ring(1), Z, Z**2, Z**3, U, Y)


def shifted_complete_ideal(shift):
    """Return I(c+shift*vZ) using its minimal monomial generators.

    The normal-form monomials have Y exponent zero or one.  Componentwise
    exceptional valuations select the complete ideal.  Equality between the
    quotient colength and the rational-double-point cycle formula certifies
    that no additional cancellation generator is missing.
    """
    shifted_cycle = cycle + shift * vZ
    candidates = []
    for y_exponent in range(2):
        for z_exponent in range(shift + 5):
            for u_exponent in range(shift + 3):
                valuation = (
                    z_exponent * vZ
                    + u_exponent * vU
                    + y_exponent * vY
                )
                if min(valuation - shifted_cycle) >= 0:
                    candidates.append(
                        Z**z_exponent * U**u_exponent * Y**y_exponent
                    )
    generators = []
    for monomial in sorted(set(candidates), key=lambda value: (value.degree(), str(value))):
        if not any(monomial.quo_rem(generator)[1] == 0 for generator in generators):
            generators.append(monomial)
    ideal = local_ring.ideal((Y**2 - U**3 - U * Z**3, *generators))
    expected_length = ZZ(shifted_cycle * e7_cartan * shifted_cycle) // 2
    assert expected_length == shift**2 + 4 * shift + 6
    assert ideal.vector_space_dimension() == expected_length
    return ideal, tuple(generators), expected_length


# C=14 is the first cutoff at which the earlier mixed-reversal diagnostic
# produced a spurious two-dimensional kernel.  Its correct common lift has
# shift C+3=17 and must be reduced in this length-363 quotient, not the
# unshifted length-six quotient above.
shifted_ideal_17, shifted_generators_17, shifted_length_17 = shifted_complete_ideal(17)
assert len(shifted_generators_17) == 39


def residue_vector(polynomial):
    """Return coefficients on 1,Z,Z^2,Z^3,U,Y modulo the E7 ideal."""
    remainder = local_ring(polynomial).reduce(complete_ideal.groebner_basis())
    answer = vector(
        QQ,
        (remainder.monomial_coefficient(monomial) for monomial in quotient_basis),
    )
    assert remainder == sum(value * monomial for value, monomial in zip(answer, quotient_basis))
    return answer


c2 = -B1 / A1
c3 = -(c2**3 + A * c2 + B) / A1
x0_jet = c2 * Z**2 + c3 * Z**3

# The translated H92 equation is zero in the quotient.  Equivalently, all
# omitted analytic-coordinate corrections lie in the complete ideal and do
# not alter the six residue rows.
translated_equation = (
    Y**2
    - (x0_jet + U)**3
    - (A1 * Z**3 + A * Z**4) * (x0_jet + U)
    - (B1 * Z**5 + B * Z**6 + B2 * Z**7)
)
assert residue_vector(translated_equation) == vector(QQ, 6)

# This is the concrete local chart map to use on any reversed/global
# numerator: t=Z, x=c2*Z^2+c3*Z^3+U, y=Y.  The columns below display the six
# lambda functionals on the smallest ambient monomials.
ambient = (local_ring(1), Z, Z**2, Z**3, x0_jet + U, Y)
residue_matrix = Matrix(QQ, [residue_vector(value) for value in ambient]).transpose()
assert residue_matrix.rank() == 6

print(
    "H92Q6E7FLAG|normal_form=Y2-U3-UZ3|"
    "graph_edges=1-4,2-5,2-7,3-4,3-6,3-7|"
    f"vZ={tuple(vZ)}|vU={tuple(vU)}|vY={tuple(vY)}",
    flush=True,
)
print(
    "H92Q6E7FLAG|cycle=vZ+vU|boundary=e1+2e5|"
    "ideal=(U2,ZU,Z4,ZY,UY)|basis=1,Z,Z2,Z3,U,Y|length=6",
    flush=True,
)
print(
    f"H92Q6E7FLAG|c2={c2}|c3={c3}|"
    f"residue_matrix={tuple(tuple(row) for row in residue_matrix.rows())}",
    flush=True,
)
print(
    "H92Q6E7FLAG|common_lift="
    "Z3*rev(N0,C)+Z*rev(Nx,C+2)*x+rev(Ny,C+3)*y|"
    "shifted_cycle=c+(C+3)vZ|"
    f"C14_length={shifted_length_17}|"
    f"C14_generators={len(shifted_generators_17)}",
    flush=True,
)
print("H92Q6E7FLAG|status=PASS_RESOLVED_E7_LOCAL_FLAG", flush=True)
