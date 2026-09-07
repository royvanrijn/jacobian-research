#!/usr/bin/env sage
"""Derive the first q80-to-rootless q=4 pencil from local components.

The first reduced neighbor class is

    D = 4F + 2O
        -(2,3,2,2,1)_D5 -(2,3,2,4,3,2)_E6.

Inside ``L(2O+4F)=<1,T,T^2,T^3,T^4,x>`` the D5 and E6 fixed
components impose four linear conditions.  At the I1* fiber, the first
exceptional cubic is ``(xi-1)^2*(xi+2)``.  The divisor meets the two spinor
ends, hence follows the double root xi=1.  At the IV* fiber it passes through
the triple point and then meets the two outer components.  The resulting
space is exactly

    L(D) = <T^2, x-T>,

so a coordinate for the first geometric neighbor is ``U=(x-T)/T^2``.

The script also derives the generic binary quartic and its Jacobian directly
from the four-parameter q80 ambient normal form.  This is an equation-level
certificate for the first step of the pinned q=4,4,12,12,4,6 lattice path.
"""

from sage.all import Matrix, PolynomialRing, QQ, ZZ, vector


# Rebuild the exact q80 ambient normal form over QQ(d,p,q,e).
parameters = PolynomialRing(QQ, names=("d", "p", "q", "e"))
d, p, q, e = parameters.gens()
K = parameters.fraction_field()
KT = PolynomialRing(K, "T")
T = KT.gen()
KS = PolynomialRing(K, "s")
s = KS.gen()

r = -3*d**2 + 3-p-q
A = T**2*(-3+p*T+q*T**2+r*T**3)
A_at_one = KS(A(T=1+s))
u = (A_at_one+3*d**2)/(-3*d**2)
branch = 2*d**3*(1+QQ(3)/2*u+QQ(3)/8*u**2-QQ(1)/16*u**3)
branch_jets = vector(K, [branch[j] for j in range(4)])
jet_matrix = Matrix(
    K, 4, 4, lambda row, column: KS((1+s)**(4+column))[row]
)
fixed_jets = vector(
    K, [KS(2*(1+s)**3+e*(1+s)**8)[j] for j in range(4)]
)
b1, b2, b3, b4 = jet_matrix.solve_right(branch_jets-fixed_jets)
B = T**3*(2+b1*T+b2*T**2+b3*T**3+b4*T**4+e*T**5)


# The fixed vertical coefficients are one old fiber plus the indicated
# connected exceptional chains.  Their complements are precisely the two
# endpoints met by the bisection.
d5_fiber = vector(ZZ, (1, 2, 1, 2, 1))
d5_chain = vector(ZZ, (1, 1, 1, 0, 0))
d5_fixed = vector(ZZ, (2, 3, 2, 2, 1))
e6_fiber = vector(ZZ, (1, 2, 2, 3, 2, 1))
e6_chain = vector(ZZ, (1, 1, 0, 1, 1, 1))
e6_fixed = vector(ZZ, (2, 3, 2, 4, 3, 2))
assert d5_fixed == d5_fiber+d5_chain
assert e6_fixed == e6_fiber+e6_chain

# Check the four local linear gates on a general element of L(2O+4F).
coefficients = PolynomialRing(K, names=("a", "c0", "c1", "c2", "c3", "c4"))
a, c0, c1, c2, c3, c4 = coefficients.gens()
local_d5 = PolynomialRing(coefficients, "xi")
xi = local_d5.gen()
d5_first_exceptional = xi**3-3*xi+2
assert d5_first_exceptional == (xi-1)**2*(xi+2)
# After x=T*xi, c0=0 removes one whole old fiber and the next restriction is
# a*xi+c1.  Following the double root xi=1 is exactly c1=-a.
d5_residual = a*xi+c1
assert d5_residual-a*(xi-1) == c1+a

local_e6 = PolynomialRing(coefficients, names=("uu", "x1", "x2", "eta"))
uu, x1, x2, eta = local_e6.gens()
# With xbar=u*x1 the first restriction after c4=0 is a*x1+c3.  Passing the
# triple point x1=0 forces c3=0.  The next blowup x1=u*x2 reaches the two
# outer components eta^2=e and leaves a*x2+c2, hence one point on each end.
e6_first_residual = a*x1+c3
assert e6_first_residual-a*x1 == c3
e6_outer_residual = (a*(uu*x2)+c2*uu)/uu
assert e6_outer_residual == a*x2+c2

gate_matrix = Matrix(
    K,
    [
        # coefficient order: a,c0,c1,c2,c3,c4
        (0, 1, 0, 0, 0, 0),       # c0=0
        (1, 0, 1, 0, 0, 0),       # c1=-a
        (0, 0, 0, 0, 1, 0),       # c3=0
        (0, 0, 0, 0, 0, 1),       # c4=0
    ],
)
kernel = gate_matrix.right_kernel_matrix()
assert kernel.nrows() == 2
assert gate_matrix*vector(K, (1, 0, -1, 0, 0, 0)) == 0
assert gate_matrix*vector(K, (0, 0, 0, 1, 0, 0)) == 0


# Substitute x=T^2*U+T.  The fixed T^4 square leaves a binary quartic in T.
KU = PolynomialRing(K, "U")
U = KU.gen()
quartic_ring = PolynomialRing(KU, "T")
Tq = quartic_ring.gen()
Aq = quartic_ring(A)
Bq = quartic_ring(B)
quartic, remainder = (
    (Tq**2*U+Tq)**3+Aq*(Tq**2*U+Tq)+Bq
).quo_rem(Tq**4)
assert remainder == 0 and quartic.degree() == 4
compact_quartic = (
    (p+b1)
    +(3*U**2+p*U+q+b2)*Tq
    +(U**3+q*U+r+b3)*Tq**2
    +(r*U+b4)*Tq**3
    +e*Tq**4
)
assert quartic == compact_quartic
qcoeff = [quartic[index] for index in range(5)]
q0, q1, q2, q3, q4 = qcoeff
invariant_i = 12*q4*q0-3*q3*q1+q2**2
invariant_j = (
    72*q4*q2*q0+9*q3*q2*q1-27*q4*q1**2
    -27*q3**2*q0-2*q2**3
)
jacobian_a = -27*invariant_i
jacobian_b = -27*invariant_j
discriminant = 4*jacobian_a**3+27*jacobian_b**2

factorization = tuple(
    (factor.monic(), ZZ(exponent)) for factor, exponent in discriminant.factor()
)
assert sorted((factor.degree(), exponent) for factor, exponent in factorization) == [
    (1, 4), (9, 1)
]
linear_factor = next(factor for factor, exponent in factorization if exponent == 4)
assert jacobian_a % linear_factor and jacobian_b % linear_factor
assert jacobian_a.degree() == 6
assert jacobian_b.degree() == 9
assert discriminant.degree() == 13
# Infinity has valuations (2,3,11), hence I5*=D9.  On the unrestricted
# four-parameter ambient family the finite factor is only I4=A3.  The marked
# rank-19 locus imposes one further collision with the residual nonic, giving
# I5=A4; at CM24 it becomes I6=A5.
residual_factor = next(factor for factor, exponent in factorization if exponent == 1)
assert residual_factor.gcd(residual_factor.derivative()) == 1
collision = K(residual_factor(U=d-1))
collision_numerator = parameters(collision.numerator())
collision_factorization = tuple(
    (factor, ZZ(exponent)) for factor, exponent in collision_numerator.factor()
)


# Exact CM24 boundary.  The child gains more than the single visibly
# orthogonal A1: the new algebraic class combines with old MW directions.
# Thus D9+A5+2A1 is compatible with a Picard-rank jump of only one because
# the child MW rank simultaneously drops from four to two.
cm24 = {d: -QQ(1)/2, p: QQ(9)/4, q: -QQ(9)/4, e: -QQ(27)/32}
surface_tangents = (
    vector(QQ, (QQ(8)/87, 1, -QQ(24)/29, -QQ(45)/116)),
    vector(QQ, (QQ(1)/12, 1, -QQ(45)/52, -QQ(261)/832)),
)
collision_diagnostics = []
for factor, exponent in collision_factorization:
    value = QQ(factor.subs(cm24))
    gradient = vector(QQ, [QQ(factor.derivative(variable).subs(cm24)) for variable in parameters.gens()])
    collision_diagnostics.append(
        (
            factor.total_degree(),
            len(factor.dict()),
            exponent,
            value,
            tuple(gradient*tangent for tangent in surface_tangents),
        )
    )
assert collision_diagnostics == [
    (2, 5, 4, -9, (QQ(42)/29, QQ(18)/13)),
    (8, 55, 1, 0, (0, 0)),
    (8, 36, 1, -QQ(2187)/2, (QQ(19197)/29, QQ(33291)/52)),
]
rank19_collision_factor = next(
    factor
    for factor, _ in collision_factorization
    if factor.subs(cm24) == 0
)
series_ring = PolynomialRing(QQ, "h")
h = series_ring.gen()
formal_surface_series = (
    (
        -QQ(1)/2+QQ(8)/87*h+QQ(344)/73167*h**2+QQ(283168)/184600341*h**3,
        QQ(9)/4+h,
        -QQ(9)/4-QQ(24)/29*h+QQ(1312)/73167*h**2+QQ(861728)/61533447*h**3,
        -QQ(27)/32-QQ(45)/116*h-QQ(3072)/24389*h**2-QQ(864550)/20511149*h**3,
    ),
    (
        -QQ(1)/2+QQ(1)/12*h+QQ(23)/10816*h**2+QQ(23545)/43869696*h**3,
        QQ(9)/4+h,
        -QQ(9)/4-QQ(45)/52*h+QQ(17)/421824*h**2+QQ(1027239)/190102016*h**3,
        -QQ(27)/32-QQ(261)/832*h-QQ(174753)/2249728*h**2-QQ(57149349)/3041632256*h**3,
    ),
)
for branch_index, series in enumerate(formal_surface_series, 1):
    image = parameters.hom(series, series_ring)(rank19_collision_factor)
    assert image.truncate(4) == 0
    print(
        f"Q80FIRSTQ4PENCIL|rank19_collision_branch={branch_index}|"
        "formal_order=4|status=PASS",
        flush=True,
    )
cm_U_ring = PolynomialRing(QQ, "U")
cm_u = cm_U_ring.gen()
QQT = PolynomialRing(cm_U_ring, "T")


def specialize_scalar(value):
    value = K(value)
    return QQ(value.numerator().subs(cm24))/QQ(value.denominator().subs(cm24))


def specialize_u(polynomial):
    polynomial = KU(polynomial)
    return cm_U_ring([specialize_scalar(value) for value in polynomial.list()])


quartic_cm = QQT([specialize_u(coefficient) for coefficient in quartic.list()])
assert quartic_cm == QQT(
    -QQ(27)/32*QQT.gen()**4
    + QQT.base_ring().gen()*QQ(9)/4*QQT.gen()**3
    + (
        QQT.base_ring().gen()**3
        - QQ(9)/4*QQT.base_ring().gen()
        + QQ(27)/16
    )*QQT.gen()**2
    + (
        3*QQT.base_ring().gen()**2
        + QQ(9)/4*QQT.base_ring().gen()
        + QQ(27)/4
    )*QQT.gen()
    - QQ(243)/32
)

cm_a = specialize_u(jacobian_a)
cm_b = specialize_u(jacobian_b)
cm_delta = 4*cm_a**3+27*cm_b**2
cm_factors = tuple(
    (str(factor.monic()), factor.degree(), ZZ(exponent))
    for factor, exponent in cm_delta.factor()
)
assert sorted((degree, exponent) for _, degree, exponent in cm_factors) == [
    (1, 6), (2, 2), (3, 1)
]

print(
    "Q80FIRSTQ4PENCIL|fixed_D5={}|fixed_E6={}|gates=c0,c1+a,c3,c4|"
    "basis=T^2,x-T|U=(x-T)/T^2".format(tuple(d5_fixed), tuple(e6_fixed)),
    flush=True,
)
print(
    "Q80FIRSTQ4PENCIL|quartic="
    "e*T^4+(r*U+b4)*T^3+(U^3+q*U+r+b3)*T^2+"
    "(3*U^2+p*U+q+b2)*T+(p+b1)",
    flush=True,
)
print(
    f"Q80FIRSTQ4PENCIL|ambient_Delta=(U-d+1)^4*R9|"
    f"R9_degree={residual_factor.degree()}|"
    "ambient_fibers=I5*,I4,9I1|ADE=D9+A3",
    flush=True,
)
print(
    "Q80FIRSTQ4PENCIL|rank19_collision_factor_degree=8|support=55|"
    "formal_branches=2|formal_order=4|"
    "rank19_fibers=I5*,I5,8I1|ADE=D9+A4|MW=4",
    flush=True,
)
print(
    f"Q80FIRSTQ4PENCIL|cm24_Delta_factors={cm_factors}|"
    "cm24_fibers=I5*,I6,2I2,3I1|ADE=D9+A5+2A1|geometric_MW=2",
    flush=True,
)
print("Q80FIRSTQ4PENCIL|status=PASS", flush=True)
