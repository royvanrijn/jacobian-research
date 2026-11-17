#!/usr/bin/env python3
"""Exact audit of the normalized-boundary A4 Keller assembly."""

import itertools
import shutil
import subprocess

import sympy as sp


lam, r, c, S, T = sp.symbols("lam r c S T")


def target_cubic(first, second, third):
    return sp.expand(
        first**3
        - 3 * first * second**2
        + 2 * second**3
        - 9 * first * second * third
        + 9 * second**2 * third
        - 27 * first * third**2
        + 27 * second * third**2
        + 27 * third**3
    )


f1 = 2 * r**3 - 3 * r**2 + 3 * r - 1
f2 = -r**3 + 3 * r - 1
f3 = r * (r - 1)
a = 2 * r - 1
b = 2 * r**2 + r - 5


# ---------------------------------------------------------------------------
# 1. Polynomial normalized boundary and unimodular ambient completion
# ---------------------------------------------------------------------------

assert sp.factor(target_cubic(lam * f1, lam * f2, lam * f3)) == 0
assert sp.expand(a * f2 + b * f3) == 1

completion_matrix = sp.Matrix([
    [f1, 1, 0],
    [f2, 0, b],
    [f3, 0, -a],
])
assert sp.factor(completion_matrix.det()) == 1

p = sp.expand(f1 * lam + S)
q = sp.expand(f2 * lam + b * T)
rho = sp.expand(f3 * lam - a * T)
new_target_cubic = target_cubic(p, q, rho)
assert new_target_cubic.subs({S: 0, T: 0}) == 0


# ---------------------------------------------------------------------------
# 2. The two inverse-mask numerators are not divisible
# ---------------------------------------------------------------------------

difference = sp.expand(p - q)
incidence = sp.expand(
    27 * rho**2
    - difference**2
    - 3 * (difference - 3 * rho) * q
)
numerator_1 = sp.expand(c - difference * r)
numerator_2 = sp.expand(-incidence * c + 27 * rho**3 * r)

variables = (lam, r, c, S, T)
_, remainder_1 = sp.div(
    numerator_1,
    new_target_cubic,
    *variables,
)
_, remainder_2 = sp.div(
    numerator_2,
    new_target_cubic,
    *variables,
)
assert remainder_1 != 0
assert remainder_2 != 0


# ---------------------------------------------------------------------------
# 3. The rational Jacobian quotient is nonconstant
# ---------------------------------------------------------------------------

old_target_cubic = target_cubic(lam, r, c)
assert sp.factor(old_target_cubic) != 0
assert sp.factor(
    new_target_cubic - old_target_cubic
) != 0

# At the old mask-zero section the prospective denominator vanishes, while
# the numerator is a nonzero polynomial.  Hence their ratio is not constant.
assert new_target_cubic.subs({S: 0, T: 0}) == 0
assert old_target_cubic.subs({S: 0, T: 0}) == old_target_cubic


# ---------------------------------------------------------------------------
# 4. The first nonautomorphic log-crepant incidence still misses divisibility
# ---------------------------------------------------------------------------

# Homogeneity gives a particularly simple solution of the corrected
# log-Jacobian equation: scale the three cubic coordinates by the first mask.
# This map is birational over S != 0 and has determinant S^3.
scaled_alpha = sp.Matrix([
    S * lam,
    S * r,
    S * c,
    S,
    T,
])
scaled_jacobian = sp.factor(
    scaled_alpha.jacobian((lam, r, c, S, T)).det()
)
scaled_target_cubic = target_cubic(S * lam, S * r, S * c)

assert scaled_jacobian == S**3
assert sp.factor(scaled_target_cubic - S**3 * old_target_cubic) == 0
assert sp.factor(
    scaled_target_cubic - old_target_cubic * scaled_jacobian
) == 0

# Inverting beta_2 after this incidence has first numerator
# S-S(lam-r)T.  Cancelling the single visible S leaves a factor congruent to
# one modulo S, whereas the denominator has S^3.  Thus this exact solution of
# the log equation cannot factor polynomially, even before the A4 lift is
# substituted.
old_difference = sp.expand(lam - r)
old_incidence = sp.expand(
    27 * c**2
    - old_difference**2
    - 3 * (old_difference - 3 * c) * r
)
scaled_difference = sp.expand(S * old_difference)
scaled_incidence = sp.expand(S**2 * old_incidence)
scaled_numerator_1 = sp.expand(S - scaled_difference * T)
scaled_numerator_2 = sp.expand(
    -scaled_incidence * S + 27 * (S * c)**3 * T
)

assert sp.expand(
    scaled_numerator_1 - S * (1 - old_difference * T)
) == 0
assert sp.expand(
    scaled_numerator_2
    - S**3 * (-old_incidence + 27 * c**3 * T)
) == 0
assert sp.expand(
    sp.cancel(scaled_numerator_1 / S).subs(S, 0)
    - (1 - old_difference * T)
) == 0
assert sp.expand(
    sp.cancel(scaled_numerator_1 / S).subs({S: 0, T: 0})
) == 1


# The same obstruction is universal for this radial base contraction, even
# when both mask outputs depend on all five variables.  The proof in the
# note uses the following exact determinant identity and the m-adic order at
# m=(lam,r,c,S).  SymPy Functions keep A and C genuinely arbitrary here.
all_variables = (lam, r, c, S, T)
arbitrary_A = sp.Function("arbitrary_A")(*all_variables)
arbitrary_C = sp.Function("arbitrary_C")(*all_variables)
radial_with_masks = sp.Matrix([
    S * lam,
    S * r,
    S * c,
    arbitrary_A,
    arbitrary_C,
])
radial_mask_jacobian = sp.factor(
    radial_with_masks.jacobian(all_variables).det()
)
euler_A = sum(
    variable * sp.diff(arbitrary_A, variable)
    for variable in (lam, r, c)
)
euler_C = sum(
    variable * sp.diff(arbitrary_C, variable)
    for variable in (lam, r, c)
)
radial_mask_expected = S**2 * (
    (S * sp.diff(arbitrary_A, S) - euler_A)
    * sp.diff(arbitrary_C, T)
    - sp.diff(arbitrary_A, T)
    * (S * sp.diff(arbitrary_C, S) - euler_C)
)
assert sp.simplify(radial_mask_jacobian - radial_mask_expected) == 0


# ---------------------------------------------------------------------------
# 5. The quadratic polar shortcut is not the missing nonradial base map
# ---------------------------------------------------------------------------

# The inverse-mask matrix is a two-generator affine-modification ledger: its
# determinant is exactly the cubic boundary.  A successful coupling can
# therefore be searched as a mixed Rees module for this matrix after the pure
# A4 lift, rather than as two unrelated scalar divisibility conditions.
inverse_mask_matrix = sp.Matrix([
    [27 * c**3, old_difference],
    [old_incidence, 1],
])
assert sp.expand(inverse_mask_matrix.det() - old_target_cubic) == 0
numerator_matrix = sp.Matrix([
    [1, -old_difference],
    [-old_incidence, 27 * c**3],
])
assert sp.expand(
    numerator_matrix * inverse_mask_matrix
    - old_target_cubic * sp.eye(2)
) == sp.zeros(2)

# Feeding the two source masks directly through the adjugate makes both
# inverse numerators divisible, but it cannot be the missing coupling.  The
# first three coordinates remain fixed, while the mask block has determinant
# B; hence the full Jacobian is B and the log-crepant identity would read
# B=B^2 already on the mask-zero section.
adjugate_mask_incidence = sp.Matrix([
    lam,
    r,
    c,
    27 * c**3 * S + old_difference * T,
    old_incidence * S + T,
])
adjugate_mask_jacobian = sp.factor(
    adjugate_mask_incidence.jacobian(all_variables).det()
)
assert adjugate_mask_jacobian == old_target_cubic
assert sp.expand(
    old_target_cubic
    - old_target_cubic * adjugate_mask_jacobian
) != 0

# The contraction module proposed in the note can be closed completely for
# the normalized nonradial triple.  A component of a pulled-back irreducible
# boundary which meets the etale locus maps dominantly to that boundary; it
# therefore detects every target function modulo the boundary equation.  In
# particular, a nonzero contraction class can occur only when every factor
# of the pulled-back boundary is supported on the critical divisor.
root_u, root_v, root_w, root_z1, root_z2 = sp.symbols(
    "root_u root_v root_w root_z1 root_z2"
)
root_H = (
    8 * root_u**3
    - 6 * root_u * root_v**2
    - 18 * root_u * root_v
    - 54 * root_u
    - 2 * root_v**3
    - 9 * root_v**2
    - 27 * root_v
    - 27
)
root_K = (
    4 * root_u**2
    + 4 * root_u * root_v
    + 6 * root_u
    + root_v**2
    + 3 * root_v
    + 9
)
root_M = root_u**2 + 2 * root_v**2 + 6 * root_v + 18
root_L = (
    root_u**3
    - 3 * root_u * root_v**2
    - 9 * root_u * root_v
    - 27 * root_u
    + 2 * root_v**3
    + 9 * root_v**2
    + 27 * root_v
    + 27
)
root_N1 = sp.expand(root_M * root_K)
root_N2 = (
    8 * root_u**3 * root_v
    + 12 * root_u**2 * root_v**2
    + 36 * root_u**2 * root_v
    + 108 * root_u**2
    + 6 * root_u * root_v**3
    + 36 * root_u * root_v**2
    + 108 * root_u * root_v
    + 162 * root_u
    + root_v**4
    + 9 * root_v**3
    + 27 * root_v**2
    + 54 * root_v
)
pure_lift = sp.Matrix([
    root_w * root_N1,
    root_w * root_N2,
    root_w * root_H,
    root_w * root_L * root_z1 / 4,
    root_z2,
])
pure_lift_jacobian = sp.factor(
    pure_lift.jacobian(
        (root_u, root_v, root_w, root_z1, root_z2)
    ).det()
)
assert pure_lift_jacobian == root_w**3 * root_K**3 * root_L**2

pure_lift_substitution = dict(zip((lam, r, c, S, T), pure_lift))
normalized_boundary_pullback = new_target_cubic.subs(
    pure_lift_substitution,
    simultaneous=True,
)
assert sp.expand(
    normalized_boundary_pullback.subs(root_w, 0)
    + 133 * root_z2**3
) == 0

# Since the specialization is a nonzero constant in root_u and root_v,
# none of root_w, root_K, or root_L divides the pulled-back boundary.  It is
# therefore coprime to the complete critical divisor of the pure lift.
assert sp.Poly(old_target_cubic, lam, r, c).is_irreducible

# The old boundary is an irreducible ordinary nodal cubic.  In the affine
# chart P=1, put Q=1+x and R=y.  Its tangent cone at [1:1:0] is the product
# of the two distinct branch directions x and x+3*y.  These exact local
# data feed the written basepoint/ramification proof excluding a quadratic
# realization B(p,q,rho)=v*B^2.  The remaining intersection calculation is
# recorded in the total-transform basis (H,E1,E2): a dominant quadratic
# rational map with two simple, infinitely-near basepoints has pullback
# line D=2H-E1-E2 and ramification class 3H-2E1-2E2, whereas the strict
# nodal cubic has class 3H-2E1-E2.  Their difference -E2 has negative
# intersection with the nef class D.
node_x, node_y = sp.symbols("node_x node_y")
node_point = {lam: 1, r: 1, c: 0}
assert old_target_cubic.subs(node_point) == 0
assert all(
    sp.diff(old_target_cubic, variable).subs(node_point) == 0
    for variable in (lam, r, c)
)
node_chart = sp.Poly(
    sp.expand(
        old_target_cubic.subs(
            {lam: 1, r: 1 + node_x, c: node_y},
            simultaneous=True,
        )
    ),
    node_x,
    node_y,
)
node_tangent_cone = sp.expand(sum(
    coefficient * node_x**monomial[0] * node_y**monomial[1]
    for monomial, coefficient in node_chart.terms()
    if sum(monomial) == 2
))
assert sp.factor(node_tangent_cone) == 3 * node_x * (node_x + 3 * node_y)

# Diagonal intersection form in the total-transform basis.
pullback_line_class = sp.Matrix([2, -1, -1])
ramification_class = sp.Matrix([3, -2, -2])
strict_boundary_class = sp.Matrix([3, -2, -1])
residual_ramification_class = ramification_class - strict_boundary_class
intersection_form = sp.diag(1, -1, -1)
assert residual_ramification_class == sp.Matrix([0, 0, -1])
assert (
    pullback_line_class.T
    * intersection_form
    * residual_ramification_class
)[0] == -1

# A second use of the nodal normalization produces a genuinely different,
# dominant base rechart.  Moving from the normalization surface in the node
# direction monomializes the boundary, but every target factor still has an
# etale source component.  The written valuation argument therefore gives a
# zero descended coupling module for this candidate as well.
node_chord_base = sp.Matrix([
    c + lam * f1,
    c + lam * f2,
    lam * f3,
])
node_chord_boundary = sp.factor(target_cubic(*node_chord_base))
node_chord_jacobian = sp.factor(
    node_chord_base.jacobian((c, lam, r)).det()
)
assert node_chord_boundary == 27 * lam**2 * c * r**3 * (r - 1)**3
assert node_chord_jacobian == 3 * lam * r**2 * (r - 1)**2

node_chord_pullback = sp.factor(
    node_chord_boundary.subs(pure_lift_substitution, simultaneous=True)
)
assert node_chord_pullback == sp.factor(
    27
    * (root_w * root_N1)**2
    * (root_w * root_H)
    * (root_w * root_N2)**3
    * (root_w * root_N2 - 1)**3
)

# Etale witnesses for lambda=0, c=0, r=0, and r=1, respectively.
# For lambda, use the noncritical component M of N1=M*K.
node_chord_etale_factors = (
    root_M,
    root_H,
    root_N2,
    root_w * root_N2 - 1,
)
assert all(
    sp.Poly(factor, root_u, root_v, root_w).is_irreducible
    for factor in node_chord_etale_factors
)
assert all(
    sp.gcd(factor, root_w * root_K * root_L) == 1
    for factor in node_chord_etale_factors
)

# The same rejection holds for all 5P3 coordinate placements of
# (mu,lambda,t).  Every coordinate hyperplane has the displayed etale
# witness, and every shifted coordinate iota(t)-1 is itself irreducible and
# coprime to the critical divisor.
coordinate_etale_witnesses = (
    root_M,
    root_N2,
    root_H,
    root_z1,
    root_z2,
)
assert all(
    sp.Poly(factor, root_u, root_v, root_w, root_z1, root_z2).is_irreducible
    and sp.gcd(factor, root_w * root_K * root_L) == 1
    for factor in coordinate_etale_witnesses
)
shifted_coordinate_witnesses = tuple(
    sp.together(coordinate - 1).as_numer_denom()[0]
    for coordinate in pure_lift
)
assert all(
    sp.Poly(factor, root_u, root_v, root_w, root_z1, root_z2).is_irreducible
    and sp.gcd(factor, root_w * root_K * root_L) == 1
    for factor in shifted_coordinate_witnesses
)
node_chord_coordinate_assignments = tuple(itertools.permutations(range(5), 3))
assert len(node_chord_coordinate_assignments) == 60

# The old cubic boundary passes the exceptional-support test, and its
# contraction module has a concrete minimal class.  The conic below is the
# image of K=0 in the cubic base.  Its pullback has exact K-order three.
coupling_conic = r**2 + 3 * r * c + 9 * c**2
root_rho = root_v**2 + 3 * root_v + 9
assert sp.expand(
    old_target_cubic
    - (
        lam**3
        + (2 * r + 3 * c - 3 * lam) * coupling_conic
    )
) == 0
assert sp.expand(
    coupling_conic.subs(pure_lift_substitution, simultaneous=True)
    - root_w**2 * root_rho * root_K**3
) == 0

old_boundary_pullback = old_target_cubic.subs(
    pure_lift_substitution,
    simultaneous=True,
)
assert sp.expand(
    old_boundary_pullback - root_w**3 * root_K**3 * root_L**2
) == 0

minimal_coupling_class = sp.expand(S**2 * coupling_conic)
minimal_coupling_pullback = minimal_coupling_class.subs(
    pure_lift_substitution,
    simultaneous=True,
)
minimal_coupling_quotient = (
    root_w * root_rho * root_z1**2 / 16
)
assert sp.expand(
    minimal_coupling_pullback
    - old_boundary_pullback * minimal_coupling_quotient
) == 0
assert sp.rem(minimal_coupling_class, old_target_cubic, lam) != 0


def homogeneous_base_forms(degree):
    """Return the monomial basis of homogeneous base forms of one degree."""

    return tuple(
        lam**lam_degree * r**r_degree * c ** (degree - lam_degree - r_degree)
        for lam_degree in range(degree + 1)
        for r_degree in range(degree - lam_degree + 1)
    )


def assert_base_factor_kernel(degree, l_exponent, expected_generators):
    """Certify a complete homogeneous K^3*L^e divisibility kernel."""

    forms = homogeneous_base_forms(degree)
    factor = sp.expand(root_K**3 * root_L**l_exponent)
    remainders = tuple(
        sp.Poly(
            sp.rem(
                form.subs(
                    {
                        lam: root_N1,
                        r: root_N2,
                        c: root_H,
                    },
                    simultaneous=True,
                ),
                factor,
                root_u,
            ),
            root_u,
            root_v,
        )
        for form in forms
    )
    monomials = tuple(
        sorted(set().union(*(set(remainder.monoms()) for remainder in remainders)))
    )
    matrix = sp.Matrix([
        [
            remainder.coeff_monomial(monomial)
            for remainder in remainders
        ]
        for monomial in monomials
    ])
    if expected_generators:
        expected = sp.Matrix.hstack(*(
            sp.Matrix([
                sp.Poly(generator, lam, r, c).coeff_monomial(form)
                for form in forms
            ])
            for generator in expected_generators
        ))
        assert matrix * expected == sp.zeros(matrix.rows, expected.cols)
        assert expected.rank() == len(expected_generators)
        assert matrix.rank() + expected.rank() == matrix.cols
    else:
        assert matrix.rank() == matrix.cols


# Decomposing by z1, z2, and W reduces the complete target-degree-four
# contraction calculation to these small homogeneous kernels.  The five
# old-boundary multiples are the trivial classes; S^2*coupling_conic is the
# unique additional class.
for base_degree in range(3):
    assert_base_factor_kernel(base_degree, 2, ())
assert_base_factor_kernel(3, 2, (old_target_cubic,))
assert_base_factor_kernel(
    4,
    2,
    (
        lam * old_target_cubic,
        r * old_target_cubic,
        c * old_target_cubic,
    ),
)
for base_degree in range(3):
    assert_base_factor_kernel(base_degree, 1, ())
assert_base_factor_kernel(3, 1, (old_target_cubic,))
for base_degree in range(2):
    assert_base_factor_kernel(base_degree, 0, ())
assert_base_factor_kernel(2, 0, (coupling_conic,))

# Globally, contraction of (K^3) from the root chart to the cubic base is
# exactly (coupling_conic,P^3).  Singular supplies an independent elimination
# certificate for the nonreduced ideal.  The L=0 curve maps nonconstantly
# into the irreducible projective cubic B=0, hence dominantly; every
# homogeneous base form whose pullback has an L factor is therefore a
# multiple of B.  Together with the independent W,z1,z2 gradings this proves
# J_B=(B,S^2*coupling_conic,S^2*P^3) in the full five-variable target ring.
root_projective_tangent = sp.expand(
    (
        sp.diff(root_N1, root_u) * root_H
        - root_N1 * sp.diff(root_H, root_u)
    )
    * sp.diff(root_L, root_v)
    - (
        sp.diff(root_N1, root_v) * root_H
        - root_N1 * sp.diff(root_H, root_v)
    )
    * sp.diff(root_L, root_u)
)
assert sp.rem(root_projective_tangent, root_L, root_u) != 0
assert sp.Poly(root_L, root_u, root_v).is_irreducible


def local_singular_expression(expression):
    return str(sp.expand(expression)).replace("**", "^")


contraction_program = f"""
option(redSB);
ring e=0,(root_u,root_v,lam,r,c),dp;
poly H={local_singular_expression(root_H)};
poly K={local_singular_expression(root_K)};
poly M={local_singular_expression(root_M)};
poly N1=M*K;
poly N2={local_singular_expression(root_N2)};
poly CK=r^2+3*r*c+9*c^2;
ideal graph=lam-N1,r-N2,c-H;
ideal contracted=std(eliminate(graph+K^3,root_u*root_v));
ideal expected=CK,lam^3;
expected=std(expected);
ideal forward=reduce(contracted,expected);
ideal reverse=reduce(expected,contracted);
int equal_ideals=1;
for (int i=1; i<=size(forward); i++)
{{
    if (forward[i]!=0) {{ equal_ideals=0; }}
}}
for (int j=1; j<=size(reverse); j++)
{{
    if (reverse[j]!=0) {{ equal_ideals=0; }}
}}
if ((size(contracted)==2) && (equal_ideals==1))
{{
    print("PASS: contraction of (K^3) is (C_K,P^3)");
}}
else
{{
    print("FAIL: K^3 contraction");
}}
quit;
"""
if shutil.which("Singular") is None:
    raise RuntimeError("Singular is required for the contraction audit")
contraction_completed = subprocess.run(
    ["Singular", "-q"],
    input=contraction_program,
    text=True,
    capture_output=True,
    check=True,
)
print(contraction_completed.stdout)
assert "FAIL:" not in contraction_completed.stdout
assert contraction_completed.stdout.count("PASS:") == 1

# A genuinely different coupling boundary would have to pull back to a pure
# W,K,L monomial.  The W grading first forces it to be homogeneous in the
# three base variables, so the bounded search is the linear identity
#
#     F(N1,N2,H) = scalar*K^a*L^b.
#
# Exhaustive linear algebra over F_101 is a rigorous characteristic-zero
# exclusion: the homogeneous substitution columns remain independent after
# reduction, and any primitive rational identity would reduce to one of the
# enumerated modular identities.  Exact characteristic-zero identities
# supply the two expected powers of B.
exceptional_prime = 101


def modular_polynomial_dictionary(expression):
    return {
        monomial: int(coefficient) % exceptional_prime
        for monomial, coefficient in sp.Poly(
            expression, root_u, root_v
        ).terms()
        if int(coefficient) % exceptional_prime
    }


def modular_leading_monomial(polynomial):
    return max(
        polynomial,
        key=lambda monomial: (
            monomial[0] + monomial[1],
            monomial[0],
            monomial[1],
        ),
    )


def modular_independent_basis(columns):
    basis = {}
    for column in columns:
        working = dict(column)
        while working:
            pivot = modular_leading_monomial(working)
            pivot_coefficient = working[pivot]
            if pivot not in basis:
                inverse = pow(
                    pivot_coefficient, -1, exceptional_prime
                )
                basis[pivot] = {
                    monomial: coefficient * inverse % exceptional_prime
                    for monomial, coefficient in working.items()
                }
                break
            for monomial, coefficient in basis[pivot].items():
                new_coefficient = (
                    working.get(monomial, 0)
                    - pivot_coefficient * coefficient
                ) % exceptional_prime
                if new_coefficient:
                    working[monomial] = new_coefficient
                else:
                    working.pop(monomial, None)
        else:
            raise AssertionError(
                "homogeneous substitution columns lost independence mod 101"
            )
    return basis


def modular_basis_contains(polynomial, basis):
    working = dict(polynomial)
    while working:
        pivot = modular_leading_monomial(working)
        pivot_coefficient = working[pivot]
        if pivot not in basis:
            return False
        for monomial, coefficient in basis[pivot].items():
            new_coefficient = (
                working.get(monomial, 0)
                - pivot_coefficient * coefficient
            ) % exceptional_prime
            if new_coefficient:
                working[monomial] = new_coefficient
            else:
                working.pop(monomial, None)
    return True


maximum_exceptional_degree = 6
root_n1_powers = [sp.Integer(1)]
root_n2_powers = [sp.Integer(1)]
root_h_powers = [sp.Integer(1)]
for _ in range(maximum_exceptional_degree):
    root_n1_powers.append(sp.expand(root_n1_powers[-1] * root_N1))
    root_n2_powers.append(sp.expand(root_n2_powers[-1] * root_N2))
    root_h_powers.append(sp.expand(root_h_powers[-1] * root_H))
root_k_powers = [sp.Integer(1)]
root_l_powers = [sp.Integer(1)]
for _ in range(2 * maximum_exceptional_degree):
    root_k_powers.append(sp.expand(root_k_powers[-1] * root_K))
for _ in range(4 * maximum_exceptional_degree // 3):
    root_l_powers.append(sp.expand(root_l_powers[-1] * root_L))

exceptional_hits = {}
for degree in range(1, maximum_exceptional_degree + 1):
    substitution_columns = []
    for lam_degree in range(degree + 1):
        for r_degree in range(degree - lam_degree + 1):
            c_degree = degree - lam_degree - r_degree
            substitution_columns.append(
                modular_polynomial_dictionary(
                    root_n1_powers[lam_degree]
                    * root_n2_powers[r_degree]
                    * root_h_powers[c_degree]
                )
            )
    modular_basis = modular_independent_basis(substitution_columns)
    degree_hits = []
    for k_exponent in range(2 * degree + 1):
        maximum_l_exponent = (4 * degree - 2 * k_exponent) // 3
        for l_exponent in range(maximum_l_exponent + 1):
            if k_exponent == 0 and l_exponent == 0:
                continue
            target = modular_polynomial_dictionary(
                root_k_powers[k_exponent]
                * root_l_powers[l_exponent]
            )
            if modular_basis_contains(target, modular_basis):
                degree_hits.append((k_exponent, l_exponent))
    exceptional_hits[degree] = tuple(degree_hits)

assert exceptional_hits == {
    1: (),
    2: (),
    3: ((3, 2),),
    4: (),
    5: (),
    6: ((6, 4),),
}
assert sp.expand(
    old_target_cubic.subs(
        {lam: root_N1, r: root_N2, c: root_H},
        simultaneous=True,
    )
    - root_K**3 * root_L**2
) == 0
print(
    "PASS: through degree six the only pure K/L boundary identities"
    " are B and B^2"
)

# A representative pair (A,C)=(X*T+g,T) makes both inverse-mask quotients
# polynomial.  It is an explicit coupled polynomial factorization, but its
# mask Jacobian is 2*S*coupling_conic, so it is neither Keller nor a recovery
# of the original A4 cover without additional base feedback.
coupled_mask_A = sp.expand(old_difference * T + minimal_coupling_class)
coupled_mask_C = T
coupled_incidence = sp.Matrix([
    lam,
    r,
    c,
    coupled_mask_A,
    coupled_mask_C,
])
coupled_incidence_jacobian = sp.factor(
    coupled_incidence.jacobian(all_variables).det()
)
assert coupled_incidence_jacobian == 2 * S * coupling_conic

coupled_numerator_1 = sp.expand(
    coupled_mask_A - old_difference * coupled_mask_C
)
coupled_numerator_2 = sp.expand(
    -old_incidence * coupled_mask_A + 27 * c**3 * coupled_mask_C
)
coupled_quotient_1 = sp.cancel(
    coupled_numerator_1.subs(
        pure_lift_substitution,
        simultaneous=True,
    )
    / old_boundary_pullback
)
coupled_quotient_2 = sp.cancel(
    coupled_numerator_2.subs(
        pure_lift_substitution,
        simultaneous=True,
    )
    / old_boundary_pullback
)
pulled_old_incidence = old_incidence.subs(
    pure_lift_substitution,
    simultaneous=True,
)
assert sp.expand(
    coupled_quotient_1 - minimal_coupling_quotient
) == 0
assert sp.expand(
    coupled_quotient_2
    - (root_z2 - pulled_old_incidence * minimal_coupling_quotient)
) == 0
coupled_composite_jacobian = sp.factor(
    coupled_incidence_jacobian.subs(
        pure_lift_substitution,
        simultaneous=True,
    )
)
assert coupled_composite_jacobian == (
    root_w**3 * root_K**3 * root_L * root_rho * root_z1 / 2
)

# Twisting the normal parameter of the node-chord map by the old boundary
# produces the first genuinely different boundary with a nonzero descended
# coupling class.  For the natural coordinate choice lambda=P, t=Q, the
# old class multiplied by the remaining chord factor gives both polynomial
# inverse masks.  Its simple C=T realization is still not log crepant.
twisted_chord_h = sp.expand(lam**2 * r**3 * (r - 1)**3)
twisted_chord_base = sp.Matrix([
    old_target_cubic + lam * f1,
    old_target_cubic + lam * f2,
    lam * f3,
])
twisted_chord_p, twisted_chord_q, twisted_chord_rho = twisted_chord_base
twisted_chord_boundary = sp.factor(target_cubic(*twisted_chord_base))
assert sp.factor(
    twisted_chord_boundary - 27 * old_target_cubic * twisted_chord_h
) == 0

twisted_chord_difference = sp.expand(twisted_chord_p - twisted_chord_q)
twisted_chord_U = sp.expand(
    27 * twisted_chord_rho**2
    - twisted_chord_difference**2
    - 3
    * (twisted_chord_difference - 3 * twisted_chord_rho)
    * twisted_chord_q
)
twisted_chord_g = sp.expand(
    27 * twisted_chord_h * minimal_coupling_class
)
twisted_chord_C = T
twisted_chord_A = sp.expand(
    twisted_chord_difference * twisted_chord_C + twisted_chord_g
)

# Chain the three small Jacobians instead of expanding one large determinant:
# chi contributes 3*lambda*t^2*(t-1)^2, (B,P,Q) contributes B_R, and the
# mask block contributes d_S(27*h*S^2*C_K).
twisted_chord_base_jacobian = sp.factor(
    3
    * lam
    * r**2
    * (r - 1)**2
    * sp.Matrix([old_target_cubic, lam, r])
    .jacobian((lam, r, c))
    .det()
)
assert twisted_chord_base_jacobian == sp.factor(
    3 * lam * r**2 * (r - 1)**2 * sp.diff(old_target_cubic, c)
)
twisted_chord_mask_jacobian = sp.factor(sp.diff(twisted_chord_g, S))
assert sp.factor(
    twisted_chord_mask_jacobian
    - 54 * twisted_chord_h * S * coupling_conic
) == 0
twisted_chord_incidence_jacobian = sp.factor(
    twisted_chord_base_jacobian * twisted_chord_mask_jacobian
)
twisted_chord_log_ratio = sp.factor(
    twisted_chord_incidence_jacobian / (27 * twisted_chord_h)
)
assert twisted_chord_log_ratio == sp.factor(
    6
    * lam
    * r**2
    * (r - 1)**2
    * S
    * coupling_conic
    * sp.diff(old_target_cubic, c)
)

twisted_chord_h_pullback = twisted_chord_h.subs(
    pure_lift_substitution,
    simultaneous=True,
)
twisted_chord_boundary_pullback = sp.expand(
    27 * old_boundary_pullback * twisted_chord_h_pullback
)
twisted_chord_numerator_1 = sp.expand(
    twisted_chord_A - twisted_chord_difference * twisted_chord_C
)
twisted_chord_numerator_2 = sp.expand(
    -twisted_chord_U * twisted_chord_A
    + 27 * twisted_chord_rho**3 * twisted_chord_C
)
assert sp.factor(
    twisted_chord_numerator_1 - twisted_chord_g
) == 0
assert sp.factor(
    twisted_chord_numerator_2
    - (
        twisted_chord_boundary * twisted_chord_C
        - twisted_chord_U * twisted_chord_g
    )
) == 0
twisted_chord_g_pullback = sp.expand(
    27 * twisted_chord_h_pullback * minimal_coupling_pullback
)
pulled_twisted_chord_U = twisted_chord_U.subs(
    pure_lift_substitution,
    simultaneous=True,
)
assert sp.expand(
    twisted_chord_g_pullback
    - twisted_chord_boundary_pullback * minimal_coupling_quotient
) == 0
twisted_chord_quotient_1 = minimal_coupling_quotient
twisted_chord_quotient_2 = sp.expand(
    root_z2 - pulled_twisted_chord_U * minimal_coupling_quotient
)
twisted_chord_composite_jacobian = sp.factor(
    twisted_chord_log_ratio.subs(
        pure_lift_substitution,
        simultaneous=True,
    )
)
assert not twisted_chord_composite_jacobian.is_constant()

# The same simple mask choice C=T is non-Keller for every ordered coordinate
# placement of (lambda,t) among P,Q,R,S.  This is a bounded twelve-case
# audit, not an exclusion of non-coordinate choices or general mask feedback.
chord_simple_coordinate_pool = (lam, r, c, S)
twisted_chord_coordinate_ratios = []
for chord_lambda, chord_t in itertools.permutations(
    chord_simple_coordinate_pool, 2
):
    chord_transverse_jacobian = sp.factor(
        sp.Matrix([
            old_target_cubic,
            chord_lambda,
            chord_t,
            minimal_coupling_class,
            T,
        ]).jacobian(all_variables).det()
    )
    chord_ratio = sp.factor(
        3
        * chord_lambda
        * chord_t**2
        * (chord_t - 1)**2
        * chord_transverse_jacobian
    )
    assert chord_ratio != 0
    assert set(chord_ratio.free_symbols) & set(all_variables)
    twisted_chord_coordinate_ratios.append(sp.factor(chord_ratio))
assert len(twisted_chord_coordinate_ratios) == 12

# For arbitrary masks the exterior-form log equation is stronger.  Etale
# valuations first force g=h*g0.  Its left side is therefore divisible by
# h=lambda^2*t^3*(t-1)^3, while the required right side has only the factor
# 9*lambda*t*(t-1).  This rejects all 5P2 coordinate choices, including
# those involving T, without choosing the two mask outputs.
twisted_chord_wedge_obstructions = []
for chord_lambda, chord_t in itertools.permutations(all_variables, 2):
    chord_h = chord_lambda**2 * chord_t**3 * (chord_t - 1)**3
    chord_required_wedge = 9 * chord_lambda * chord_t * (chord_t - 1)
    chord_wedge_quotient = sp.cancel(chord_required_wedge / chord_h)
    assert sp.denom(chord_wedge_quotient) != 1
    twisted_chord_wedge_obstructions.append(chord_wedge_quotient)
assert len(twisted_chord_wedge_obstructions) == 20

# A tempting nonradial solution of the log equation would be a quadratic
# polar/adjugate map of the cubic.  It would require both the Hessian
# determinant and B(grad B) to be supported on B.  The rational smooth point
# below lies on B, but both quantities are nonzero there, so this shortcut is
# excluded without expanding the large sextic B(grad B).
base_variables = (lam, r, c)
polar_map = sp.Matrix([
    sp.diff(old_target_cubic, variable)
    for variable in base_variables
])
polar_boundary = target_cubic(*polar_map)
polar_jacobian = sp.hessian(old_target_cubic, base_variables).det()
smooth_boundary_point = {lam: 9, r: -3, c: 2}
assert old_target_cubic.subs(smooth_boundary_point) == 0
assert polar_map.subs(smooth_boundary_point) == sp.Matrix([162, 54, -648])
assert polar_boundary.subs(smooth_boundary_point) == -8533918944
assert polar_jacobian.subs(smooth_boundary_point) == 1259712

# The first adjoint-divisor completion suggested by free-divisor theory is
# B times an arbitrary hyperplane.  If this homogeneous quartic were free,
# its two Jacobian syzygy degrees would sum to three, so at least one
# non-Euler syzygy would have linear coefficients.  The complete coefficient
# matrix for such a syzygy has full rank for every nonzero hyperplane: its
# maximal-minor ideal in the hyperplane parameters is the ninth power of the
# irrelevant ideal.
ell_0, ell_1, ell_2 = sp.symbols("ell_0 ell_1 ell_2")
adjoint_linear_form = ell_0 * lam + ell_1 * r + ell_2 * c
adjoint_quartic = sp.expand(old_target_cubic * adjoint_linear_form)
degree_four_monomials = tuple(
    lam**lam_power * r**r_power * c ** (4 - lam_power - r_power)
    for lam_power in range(5)
    for r_power in range(5 - lam_power)
)
linear_syzygy_columns = tuple(
    sp.Poly(
        sp.expand(
            coefficient * sp.diff(adjoint_quartic, derivative_variable)
        ),
        lam,
        r,
        c,
    )
    for derivative_variable in base_variables
    for coefficient in base_variables
)
linear_syzygy_matrix = sp.Matrix([
    [
        column.coeff_monomial(monomial)
        for column in linear_syzygy_columns
    ]
    for monomial in degree_four_monomials
])
assert linear_syzygy_matrix.shape == (15, 9)


def singular_expression(expression):
    return str(sp.expand(expression)).replace("**", "^")


linear_syzygy_entries = ",\n".join(
    singular_expression(linear_syzygy_matrix[row, column])
    for row in range(15)
    for column in range(9)
)
singular_program = f"""
option(redSB);
ring q=0,(ell_0,ell_1,ell_2),dp;
matrix M[15][9]={linear_syzygy_entries};
ideal I=minor(M,9);
ideal G=std(I);
ideal E=maxideal(9);
ideal GE=reduce(G,E);
ideal EG=reduce(E,G);
int equal_ideals=1;
for (int i=1; i<=size(GE); i++)
{{
    if (GE[i]!=0) {{ equal_ideals=0; }}
}}
for (int j=1; j<=size(EG); j++)
{{
    if (EG[j]!=0) {{ equal_ideals=0; }}
}}
if ((size(G)==55) && (equal_ideals==1))
{{
    print("PASS: the adjoint-hyperplane minor ideal is (ell_0,ell_1,ell_2)^9");
}}
else
{{
    print("FAIL: adjoint-hyperplane minor ideal");
}}
quit;
"""
if shutil.which("Singular") is None:
    raise RuntimeError("Singular is required for the adjoint-divisor audit")
singular_completed = subprocess.run(
    ["Singular", "-q"],
    input=singular_program,
    text=True,
    capture_output=True,
    check=True,
)
print(singular_completed.stdout)
assert "FAIL:" not in singular_completed.stdout
assert singular_completed.stdout.count("PASS:") == 1


print("PASS: the normalized target boundary has a unimodular A5 completion")
print("PASS: both two-mask inverse numerators have nonzero remainders")
print("PASS: the assembled rational candidate has nonconstant Jacobian")
print("PASS: homogeneous scaling is a birational nonautomorphic log rechart")
print("OBSTRUCTION: its first inverse mask has unavoidable residue 1 mod S")
print("OBSTRUCTION: arbitrary polynomial masks cannot repair the radial base")
print("PASS: the inverse masks form a determinant-B mixed Rees matrix")
print("OBSTRUCTION: its direct adjugate mask block has Jacobian B, not one")
print("PASS: the normalized pulled boundary meets the etale locus")
print("OBSTRUCTION: its descended contraction module J_pi/(B_pi) is zero")
print("PASS: B is an irreducible ordinary nodal cubic with two tangent branches")
print("OBSTRUCTION: the quadratic B^2 basepoint ledger leaves residual class -E2")
print("PASS: the node-chord base rechart monomializes B and is dominant")
print("OBSTRUCTION: all four node-chord boundary factors retain etale witnesses")
print("OBSTRUCTION: all 60 coordinate placements of the node-chord parameters fail")
print("PASS: S^2*(Q^2+3*Q*R+9*R^2) is the unique new degree-four coupling class")
print("PASS: J_B=(B,S^2*C_K,S^2*P^3) is the complete old-boundary contraction")
print("PASS: its displayed mask pair gives two exact polynomial quotients")
print("OBSTRUCTION: that coupled polynomial map has nonconstant Jacobian")
print("PASS: twisting the node chord by B gives a different-boundary coupling")
print("OBSTRUCTION: all 20 coordinate B-twists fail even with arbitrary masks")
print("OBSTRUCTION: the complete module closes every boundary-preserving incidence")
print("OBSTRUCTION: the quadratic polar/adjugate shortcut misses B at a smooth point")
print("OBSTRUCTION: no single adjoint hyperplane makes B a Saito-free divisor")
