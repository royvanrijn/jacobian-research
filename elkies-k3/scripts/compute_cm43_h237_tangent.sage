#!/usr/bin/env sage
"""Audit the attempted H237 tangent at the exact CM(-43) Humbert-8 point.

At the CM point the generic level-79 section is

    Q79 = 4*P1 - 5*P2 + P3,

where P1,P2 have height 5/2 and P3 is the generic height-four Humbert
section.  Write

    x(Q79)=Nx/h^2,  y(Q79)=Ny/h^3,

with degrees 120,180,58.  Linearizing the cleared section identity in the two
Humbert parameters ``(r,s)`` gives a 361 by 362 finite-field system.  At the
CM(-43) point this first-order system is deliberately singular: the section
block has a two-dimensional kernel and both parameter directions lift to
first order.  Projecting the second-order residual to the three-dimensional
cokernel leaves three quadrics in four first-order kernel coordinates.  A
third-order calculation on their reduced radical is included together with
the essential audit of the freedom to add kernel vectors to the second
correction.  That audit shows that the apparent canonical-slice cubic is
fully gauge-absorbed and is therefore *not* a Humbert-237 tangent equation.
"""

from sage.all import *
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--p", type=int, default=101)
parser.add_argument("--print-cone", action="store_true")
parser.add_argument("--third-order", action="store_true")
args = parser.parse_args()

EXPECTED_CANONICAL_SLICE_CUBIC = (
    QQ(1),
    -QQ(32777647185971477137326047735483125)
    / QQ(5537010211609548283434042242558544),
    -QQ(62256599453976317430685478929919140625)
    / QQ(4252423842516133081677344442284961792),
    QQ(16726681288628536079919155997729849853515625)
    / QQ(446402445291973586382160910173306149076992),
)

p = ZZ(args.p)
assert p.is_prime() and p not in (2, 3, 5, 7, 19, 29, 37, 43, 241, 349)
field = GF(p)
RT = PolynomialRing(field, "T")
T = RT.gen()
KT = RT.fraction_field()


def red(value):
    value = QQ(value)
    assert value.denominator() % p
    return field(value.numerator()) / field(value.denominator())


r = red(-QQ(1225)/722)
s = red(-QQ(93312)/442225)

A1 = 2*r*s**2
A = -(9*r*s + 4*r**2 + 4*r + 1)/3
B1 = r*s**2*(3*s + 8*r - 2)/3
B0 = -(
    54*r**2*s + 81*r*s - 16*r**3 - 24*r**2 - 12*r - 2
)/27
B2 = r**2
a4 = RT(A1*T**3 + A*T**4)
a6 = RT(B1*T**5 + B0*T**6 + B2*T**7)
curve = EllipticCurve(KT, [0, 0, 0, a4, a6])

# The two rational height-5/2 sections.
c0 = red(QQ(1194481)/442225)
W = red(QQ(663613890625)/34828517376)
section_data = (
    (
        red(QQ(684775)/93312),
        tuple(map(red, (
            QQ(914233879)/294079625,
            QQ(26371835)/1181952,
            QQ(557834834375)/11609505792,
            QQ(540596465650390625)/6499837226778624,
        ))),
    ),
    (
        red(-QQ(1765225)/93312),
        tuple(map(red, (
            -QQ(1085766121)/294079625,
            QQ(57241835)/1181952,
            -QQ(1437996415625)/11609505792,
            QQ(540596465650390625)/6499837226778624,
        ))),
    ),
)
points = []
for c, d in section_data:
    x = RT(T**2*(c0+c*T+W*T**2))
    y = RT(T**3*sum(d[index]*T**index for index in range(4)))
    points.append(curve((x, y)))

# The generic height-four section from the second D9+E7 quartic point.
quartic_b = -3*s*T**2
quartic_c = 3*s**2*T + (2*r+1)*T**2
quartic_d = -s**3 - 2*r*s*T + r*T**2 - 2*s*T
x_long = quartic_d**2/(4*s**2) - quartic_c
x3 = RT(x_long + quartic_c/3)
y3 = RT(-quartic_b*s - quartic_d*x_long/(2*s))
point3 = curve((x3, y3))

Q79 = 4*points[0] - 5*points[1] + point3
Nx = RT(Q79[0].numerator())
Ny = RT(Q79[1].numerator())
h = RT(Q79[0].denominator().sqrt())
assert Q79[0].denominator() == h**2
assert Q79[1].denominator() == h**3
assert (Nx.degree(), Ny.degree(), h.degree()) == (120, 180, 58)
assert Ny**2 == Nx**3 + a4*Nx*h**4 + a6*h**6

# Parameter derivatives of the Kumar coefficients, evaluated at (r,s).
a4_r = RT(2*s**2*T**3 - (9*s+8*r+4)/3*T**4)
a4_s = RT(4*r*s*T**3 - 3*r*T**4)
a6_r = RT(
    (s**2*(3*s+16*r-2)/3)*T**5
    - (108*r*s+81*s-48*r**2-48*r-12)/27*T**6
    + 2*r*T**7
)
a6_s = RT(
    (r*(9*s**2+16*r*s-4*s)/3)*T**5
    - (54*r**2+81*r)/27*T**6
)

degree_bound = 360


def coefficients(poly):
    poly = RT(poly)
    return [poly[index] for index in range(degree_bound+1)]


# Linearization of
#
#   Ny^2-Nx^3-a4*Nx*h^4-a6*h^6 = 0.
columns = []
for index in range(181):
    columns.append(coefficients(2*Ny*T**index))
for index in range(121):
    columns.append(coefficients(-(3*Nx**2+a4*h**4)*T**index))
for index in range(58):
    columns.append(coefficients(
        -(4*a4*Nx*h**3+6*a6*h**5)*T**index
    ))
columns.append(coefficients(-(a4_r*Nx*h**4+a6_r*h**6)))
columns.append(coefficients(-(a4_s*Nx*h**4+a6_s*h**6)))

matrix_full = matrix(field, columns).transpose()
matrix_section = matrix_full[:, :360]
assert matrix_full.dimensions() == (361, 362)
section_rank = matrix_section.rank()
full_rank = matrix_full.rank()
print(
    f"CM43H237TAN|p={p}|diagnostic_section_rank={section_rank}"
    f"|diagnostic_full_rank={full_rank}",
    flush=True,
)
assert section_rank == 358
assert full_rank == 358
kernel = matrix_full.right_kernel_matrix()
cokernel = matrix_full.left_kernel_matrix()
assert kernel.nrows() == 4
assert cokernel.nrows() == 3
assert matrix(field, [row[-2:] for row in kernel]).rank() == 2


def second_parameter_coefficients(dr, ds):
    """Coefficient of epsilon^2 in a_i(r+epsilon*dr,s+epsilon*ds)."""
    qa4 = RT(
        2*(r*ds**2 + 2*s*dr*ds)*T**3
        - (9*dr*ds + 4*dr**2)/3*T**4
    )
    Rrs = PolynomialRing(field, names=("rr", "ss"))
    rr, ss = Rrs.gens()
    b1 = rr*ss**2*(3*ss + 8*rr - 2)/3
    b0 = -(
        54*rr**2*ss + 81*rr*ss - 16*rr**3
        - 24*rr**2 - 12*rr - 2
    )/27

    def quadratic_coefficient(poly):
        values = {rr: r, ss: s}
        return field(
            poly.derivative(rr, 2).subs(values)*dr**2/2
            + poly.derivative(rr).derivative(ss).subs(values)*dr*ds
            + poly.derivative(ss, 2).subs(values)*ds**2/2
        )

    qa6 = RT(
        quadratic_coefficient(b1)*T**5
        + quadratic_coefficient(b0)*T**6
        + dr**2*T**7
    )
    return qa4, qa6


def split_kernel_vector(vector_data):
    dy = RT(sum(vector_data[index]*T**index for index in range(181)))
    dx = RT(sum(vector_data[181+index]*T**index for index in range(121)))
    dh = RT(sum(vector_data[302+index]*T**index for index in range(58)))
    return dy, dx, dh, vector_data[-2], vector_data[-1]


def second_order_residual_vector(vector_data):
    dy, dx, dh, dr, ds = split_kernel_vector(vector_data)
    da4 = a4_r*dr + a4_s*ds
    da6 = a6_r*dr + a6_s*ds
    qa4, qa6 = second_parameter_coefficients(dr, ds)
    residual = RT(
        dy**2 - 3*Nx*dx**2
        - a4*(6*Nx*h**2*dh**2 + 4*dx*h**3*dh)
        - da4*(dx*h**4 + 4*Nx*h**3*dh)
        - qa4*Nx*h**4
        - 15*a6*h**4*dh**2
        - 6*da6*h**5*dh
        - qa6*h**6
    )
    return vector(field, coefficients(residual))


def second_order_obstruction(vector_data):
    residual_vector = second_order_residual_vector(vector_data)
    return vector(field, [row*residual_vector for row in cokernel])


# Recover the three obstruction quadrics by polarization on a kernel basis.
Rc = PolynomialRing(field, names=("c0", "c1", "c2", "c3"))
cs = Rc.gens()
diagonal = [second_order_obstruction(kernel[index]) for index in range(4)]
quadrics = [sum(diagonal[i][j]*cs[i]**2 for i in range(4)) for j in range(3)]
for i in range(4):
    for k in range(i+1, 4):
        mixed = (
            second_order_obstruction(kernel[i]+kernel[k])
            - diagonal[i] - diagonal[k]
        )
        for j in range(3):
            quadrics[j] += mixed[j]*cs[i]*cs[k]

dr_linear = sum(kernel[index][-2]*cs[index] for index in range(4))
ds_linear = sum(kernel[index][-1]*cs[index] for index in range(4))
assert all(poly.degree() == 2 for poly in quadrics)
if args.print_cone:
    print(
        f"CM43H237TAN|p={p}|quadrics=" + ";".join(map(str, quadrics)),
        flush=True,
    )
    print(
        f"CM43H237TAN|p={p}|parameter_projection="
        f"dr:{dr_linear};ds:{ds_linear}",
        flush=True,
    )


def affine_solutions(normalizer):
    ideal = Rc.ideal(quadrics + [normalizer-1])
    if args.print_cone:
        print(
            f"CM43H237TAN|p={p}|chart={normalizer}=1|groebner="
            + ";".join(map(str, ideal.groebner_basis())),
            flush=True,
        )
    if ideal.dimension() != 0:
        return None, ideal.dimension()
    return ideal.variety(), 0


solutions_ds, dimension_ds = affine_solutions(ds_linear)
solutions_dr, dimension_dr = affine_solutions(dr_linear)
print(
    f"CM43H237TAN|p={p}|diagnostic_cone_dimensions="
    f"ds:{dimension_ds},dr:{dimension_dr}",
    flush=True,
)
nonisolated = dimension_ds != 0 or dimension_dr != 0
if nonisolated:
    print(
        f"CM43H237TAN|p={p}|section_matrix=361x360|section_rank=358"
        f"|full_rank=358|kernel=4|cokernel=3",
        flush=True,
    )
    print(
        f"CM43H237TAN|p={p}|second_order_projective_dimension="
        f"{max(dimension_ds, dimension_dr)}"
        "|status=NONISOLATED_SECOND_ORDER",
        flush=True,
    )
    if args.third_order:
        # The reduced quadratic cone is c0=c1=0 at every tested good prime.
        # Work only on its two-dimensional radical, with coordinates (u,v).
        Ru = PolynomialRing(field, names=("u", "v"))
        u, v = Ru.gens()
        basis_u, basis_v = kernel[2], kernel[3]
        q_uu = second_order_residual_vector(basis_u)
        q_vv = second_order_residual_vector(basis_v)
        q_uv = (
            second_order_residual_vector(basis_u+basis_v) - q_uu - q_vv
        )
        q2 = vector(
            Ru,
            [q_uu[i]*u**2 + q_uv[i]*u*v + q_vv[i]*v**2 for i in range(361)],
        )
        assert not any(cokernel.change_ring(Ru)*q2)

        # Choose a fixed right inverse of the rank-358 linearized matrix and
        # solve for a canonical second correction.  Changing the right inverse
        # changes the lift, not the projected third-order obstruction locus.
        pivot_columns = tuple(matrix_full.pivots())
        column_basis = matrix_full.matrix_from_columns(pivot_columns)
        pivot_rows = tuple(column_basis.transpose().pivots())
        square = matrix_full.matrix_from_rows_and_columns(
            pivot_rows, pivot_columns
        )
        assert square.dimensions() == (358, 358) and square.is_invertible()
        selected_rhs = vector(Ru, [-q2[row] for row in pivot_rows])
        # Invert over the finite field first.  Asking Sage to echelonize the
        # same 358-square matrix over Ru is needlessly expensive.
        selected_solution = square.inverse().change_ring(Ru)*selected_rhs
        correction2 = vector(Ru, [0]*362)
        for column, value in zip(pivot_columns, selected_solution):
            correction2[column] = value
        assert matrix_full.change_ring(Ru)*correction2 == -q2

        correction1 = vector(
            Ru, [u*basis_u[index] + v*basis_v[index] for index in range(362)]
        )
        Re = PolynomialRing(Ru, "epsilon")
        epsilon = Re.gen()
        RTe = PolynomialRing(Re, "TT")
        TT = RTe.gen()

        def expanded_polynomial(base, start, length):
            return RTe(sum(
                (
                    Re(base[index])
                    + epsilon*Re(correction1[start+index])
                    + epsilon**2*Re(correction2[start+index])
                )*TT**index
                for index in range(length)
            ))

        Nye = expanded_polynomial(Ny, 0, 181)
        Nxe = expanded_polynomial(Nx, 181, 121)
        he = expanded_polynomial(h, 302, 58) + TT**58
        # h is monic and correction2 has no leading coefficient variable;
        # expanded_polynomial included the original lower coefficients only.
        assert h[58] == 1
        re = Re(r) + epsilon*Re(correction1[-2]) + epsilon**2*Re(correction2[-2])
        se = Re(s) + epsilon*Re(correction1[-1]) + epsilon**2*Re(correction2[-1])
        a4e = RTe(
            2*re*se**2*TT**3
            - (9*re*se + 4*re**2 + 4*re + 1)/3*TT**4
        )
        a6e = RTe(
            re*se**2*(3*se + 8*re - 2)/3*TT**5
            - (
                54*re**2*se + 81*re*se - 16*re**3
                - 24*re**2 - 12*re - 2
            )/27*TT**6
            + re**2*TT**7
        )
        expanded_identity = Nye**2 - Nxe**3 - a4e*Nxe*he**4 - a6e*he**6
        for order in range(3):
            assert all(expanded_identity[index][order] == 0 for index in range(361))
        residual3 = vector(
            Ru, [expanded_identity[index][3] for index in range(361)]
        )
        third_obstruction_vector = cokernel.change_ring(Ru)*residual3
        third_obstructions = tuple(
            polynomial for polynomial in third_obstruction_vector if polynomial
        )

        # A second correction is defined only modulo ker(matrix_full).  Its
        # effect on the third obstruction is the polarized quadratic residual
        # B(v1,k).  Audit this gauge before interpreting the canonical-slice
        # cubic as a genuine parameter obstruction.
        lifting_columns = []
        for kernel_index in range(4):
            kernel_vector = kernel[kernel_index]
            effect_u = cokernel*vector(
                field,
                second_order_residual_vector(basis_u+kernel_vector)
                - second_order_residual_vector(basis_u)
                - second_order_residual_vector(kernel_vector),
            )
            effect_v = cokernel*vector(
                field,
                second_order_residual_vector(basis_v+kernel_vector)
                - second_order_residual_vector(basis_v)
                - second_order_residual_vector(kernel_vector),
            )
            lifting_columns.append(
                vector(Ru, [effect_u[row]*u + effect_v[row]*v for row in range(3)])
            )
        lifting_matrix = matrix(Ru, lifting_columns).transpose()
        lifting_field = Ru.fraction_field()
        lifting_rank = lifting_matrix.change_ring(lifting_field).rank()
        augmented_rank = lifting_matrix.augment(
            matrix(Ru, 3, 1, list(third_obstruction_vector))
        ).change_ring(lifting_field).rank()
        print(
            f"CM43H237TAN3|p={p}|second_correction_gauge_rank={lifting_rank}"
            f"|augmented_rank={augmented_rank}",
            flush=True,
        )
        third_ideal = Ru.ideal(third_obstructions)
        print(
            f"CM43H237TAN3|p={p}|equations={len(third_obstructions)}"
            f"|dimension={third_ideal.dimension()}|degrees="
            + ",".join(map(str, sorted(set(poly.total_degree() for poly in third_obstructions)))),
            flush=True,
        )
        if args.print_cone:
            print(
                f"CM43H237TAN3|p={p}|obstructions="
                + ";".join(map(str, third_obstructions)),
                flush=True,
            )

        # Express the binary cubic in the invariant parameter tangent
        # coordinates (D,S)=(dr,ds), rather than in the prime-dependent
        # echelon kernel coordinates (u,v).  Normalize projectively by the
        # first nonzero coefficient in D^3,D^2*S,D*S^2,S^3 order.
        Rd = PolynomialRing(field, names=("D", "S"))
        Dvar, Svar = Rd.gens()
        parameter_matrix = matrix(field, (
            (basis_u[-2], basis_v[-2]),
            (basis_u[-1], basis_v[-1]),
        ))
        assert parameter_matrix.is_invertible()
        parameter_inverse = parameter_matrix.inverse()
        u_expression = parameter_inverse[0, 0]*Dvar + parameter_inverse[0, 1]*Svar
        v_expression = parameter_inverse[1, 0]*Dvar + parameter_inverse[1, 1]*Svar
        projected_cubic = Rd(sum(
            coefficient*u_expression**exponents[0]*v_expression**exponents[1]
            for exponents, coefficient in third_obstructions[0].dict().items()
        ))
        projected_coefficients = [
            projected_cubic.monomial_coefficient(monomial)
            for monomial in (Dvar**3, Dvar**2*Svar, Dvar*Svar**2, Svar**3)
        ]
        leading = next(value for value in projected_coefficients if value)
        projected_coefficients = [value/leading for value in projected_coefficients]
        expected_match = "normalization_ramified"
        if all(
            value.denominator() % p
            for value in EXPECTED_CANONICAL_SLICE_CUBIC
        ):
            expected_reduction = [
                red(value) for value in EXPECTED_CANONICAL_SLICE_CUBIC
            ]
            assert projected_coefficients == expected_reduction
            expected_match = "1"
        print(
            f"CM43H237TAN3|p={p}|canonical_slice_cubic_coefficients="
            + ",".join(str(int(value)) for value in projected_coefficients)
            + f"|expected_match={expected_match}",
            flush=True,
        )

        Rz = PolynomialRing(field, "z")
        zeta = Rz.gen()

        def affine_binary(poly):
            return Rz(sum(
                coefficient*zeta**exponents[0]
                for exponents, coefficient in poly.dict().items()
            ))

        affine_gcd = Rz(0)
        for polynomial in third_obstructions:
            affine_gcd = gcd(affine_gcd, affine_binary(polynomial))
        tangent_uv = [(root, field(1)) for root, _ in affine_gcd.roots()]
        if third_obstructions and all(poly(u=1, v=0) == 0 for poly in third_obstructions):
            tangent_uv.append((field(1), field(0)))
        tangent_directions = set()
        for u_value, v_value in tangent_uv:
            dr_value = basis_u[-2]*u_value + basis_v[-2]*v_value
            ds_value = basis_u[-1]*u_value + basis_v[-1]*v_value
            if ds_value:
                tangent_directions.add((dr_value/ds_value, field(1)))
            elif dr_value:
                tangent_directions.add((field(1), field(0)))
        assert lifting_rank == augmented_rank == 1
        print(
            f"CM43H237TAN3|p={p}|canonical_slice_gcd={affine_gcd}"
            f"|slice_directions={len(tangent_directions)}|dr:ds="
            + ";".join(
                f"{int(direction[0])}:{int(direction[1])}"
                for direction in sorted(
                    tangent_directions,
                    key=lambda item: (int(item[1]), int(item[0])),
                )
            )
            + "|status=GAUGE_ABSORBED_NOT_A_TANGENT",
            flush=True,
        )
else:
    directions = set()
    for solution in solutions_ds + solutions_dr:
        dr_value = dr_linear(solution)
        ds_value = ds_linear(solution)
        if (dr_value, ds_value) == (0, 0):
            continue
        if ds_value:
            directions.add((dr_value/ds_value, field(1)))
        else:
            directions.add((field(1), field(0)))

    direction_text = ";".join(
        f"{int(direction[0])}:{int(direction[1])}"
        for direction in sorted(
            directions, key=lambda item: (int(item[1]), int(item[0]))
        )
    )
    print(
        f"CM43H237TAN|p={p}|section_matrix=361x360|section_rank=358"
        f"|full_rank=358|kernel=4|cokernel=3",
        flush=True,
    )
    print(
        f"CM43H237TAN|p={p}|projected_directions={len(directions)}"
        f"|dr:ds={direction_text}|status=PASS",
        flush=True,
    )
