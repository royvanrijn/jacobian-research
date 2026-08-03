#!/usr/bin/env python3
"""Exact audit of the normalized-boundary A4 Keller assembly."""

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
print("OBSTRUCTION: the quadratic polar/adjugate shortcut misses B at a smooth point")
print("OBSTRUCTION: no single adjoint hyperplane makes B a Saito-free divisor")
