#!/usr/bin/env python3
"""Exact Ferrand-norm comparison for the first higher transfer blocks.

Over a characteristic-zero base, divided powers agree with symmetric
tensors.  For a finite free algebra B/A of rank h and a B-algebra C, Ferrand's
norm is

    Gamma_A^h(C) tensor_{Gamma_A^h(B)} A.

The tensor product uses the determinant norm Gamma_A^h(B) -> A.  This script
constructs that quotient directly from the universal normic-law relations

    gamma^h(bm) = Nm_{B/A}(b) gamma^h(m)

and uses the componentwise multiplication on symmetric tensors.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import combinations_with_replacement, permutations, product
import shutil
import subprocess
import sys

import sympy as sp
from sympy.polys.matrices import DomainMatrix


def companion_matrix(coefficients: tuple[sp.Expr, ...]) -> sp.Matrix:
    """Multiplication by z in A[z]/(z^h + ... + coefficients[0])."""

    degree = len(coefficients)
    matrix = sp.zeros(degree)
    for column in range(degree - 1):
        matrix[column + 1, column] = 1
    for row, coefficient in enumerate(coefficients):
        matrix[row, degree - 1] = -coefficient
    return matrix


class FerrandNorm:
    """Symmetric-tensor model of a monogenic Ferrand norm in characteristic 0."""

    def __init__(
        self,
        coefficients: tuple[sp.Expr, ...],
        parameters: tuple[sp.Symbol, ...],
        preferred_free_labels: tuple[tuple[int, ...], ...],
    ) -> None:
        self.degree = len(coefficients)
        self.parameters = parameters
        self.c_dimension = 2 * self.degree
        self.orbits = tuple(
            combinations_with_replacement(
                range(self.c_dimension), self.degree
            )
        )
        self.orbit_index = {
            orbit: index for index, orbit in enumerate(self.orbits)
        }
        self.z_matrix = companion_matrix(coefficients)
        self.preferred_free = tuple(
            self.orbit_index[label] for label in preferred_free_labels
        )
        self._construct_normic_quotient()

    def _construct_normic_quotient(self) -> None:
        h = self.degree
        b = sp.symbols(f"b0:{h}")
        m = sp.symbols(f"m0:{2 * h}")
        multiplication_by_b = sum(
            (b[power] * self.z_matrix**power for power in range(h)),
            start=sp.zeros(h),
        )
        determinant_norm = sp.factor(multiplication_by_b.det())
        bm: list[sp.Expr] = []
        for epsilon_degree in range(2):
            block = sp.Matrix(
                m[epsilon_degree * h : (epsilon_degree + 1) * h]
            )
            bm.extend(multiplication_by_b * block)

        polynomial_variables = (*b, *m)
        coefficient_vectors: dict[tuple[int, ...], dict[int, sp.Expr]] = {}
        for row, orbit in enumerate(self.orbits):
            left = sp.prod(bm[index] for index in orbit)
            right = determinant_norm * sp.prod(m[index] for index in orbit)
            relation = sp.Poly(
                sp.expand(left - right), *polynomial_variables
            )
            for monomial, coefficient in relation.terms():
                coefficient_vectors.setdefault(monomial, {})[row] = coefficient

        relation_rows = [
            [entries.get(index, 0) for index in range(len(self.orbits))]
            for entries in coefficient_vectors.values()
        ]
        pivot_indices = tuple(
            index
            for index in range(len(self.orbits))
            if index not in self.preferred_free
        )
        column_order = (*pivot_indices, *self.preferred_free)
        coefficient_field = (
            sp.QQ.frac_field(*self.parameters) if self.parameters else sp.QQ
        )
        reordered_rows = [
            [row[index] for index in column_order] for row in relation_rows
        ]
        domain_matrix = DomainMatrix.from_list_sympy(
            len(reordered_rows), len(self.orbits), reordered_rows
        ).convert_to(coefficient_field)
        reduced, pivots = domain_matrix.rref()

        expected_rank = len(self.orbits) - len(self.preferred_free)
        assert pivots == tuple(range(expected_rank))
        relation_matrix = reduced.to_Matrix()[:expected_rank, :]
        self.pivot_indices = pivot_indices
        self.reduction_block = relation_matrix[:, expected_rank:]

        # The preferred bases work over the polynomial base, not merely over
        # its fraction field: every reduction coefficient has denominator 1.
        for entry in self.reduction_block:
            assert sp.cancel(entry).as_numer_denom()[1] == 1
        if self.parameters:
            self._certify_integral_basis(relation_rows)

    def _certify_integral_basis(
        self, relation_rows: list[list[sp.Expr]]
    ) -> None:
        """Use a polynomial-module standard basis to exclude hidden torsion."""

        singular = shutil.which("Singular")
        assert singular is not None, (
            "Singular is required for the integral certificate"
        )
        column_order = (*self.preferred_free, *self.pivot_indices)
        new_position = {
            old_position: new_position + 1
            for new_position, old_position in enumerate(column_order)
        }

        def singular_expression(expression: sp.Expr) -> str:
            return str(sp.expand(expression)).replace("**", "^")

        module_generators: list[str] = []
        for row in relation_rows:
            terms = [
                f"({singular_expression(coefficient)})"
                f"*gen({new_position[index]})"
                for index, coefficient in enumerate(row)
                if coefficient
            ]
            module_generators.append("+".join(terms))

        parameter_names = ",".join(str(parameter) for parameter in self.parameters)
        free_rank = len(self.preferred_free)
        relation_rank = len(self.pivot_indices)
        singular_program = "\n".join(
            (
                f"ring R=0,({parameter_names}),(C,dp);",
                "module M=" + ",\n".join(module_generators) + ";",
                "option(redSB);",
                "module G=std(M);",
                "int ok=1;",
                f"if (size(G)!={relation_rank}) {{ ok=0; }}",
                "int i;",
                f"for (i=1;i<={relation_rank};i++)",
                "{",
                f"  if (leadmonom(G[i])!=gen(i+{free_rank})) {{ ok=0; }}",
                "}",
                'if (ok==1) { "INTEGRAL_BASIS_PASS"; }',
            )
        )
        completed = subprocess.run(
            (singular, "-q"),
            input=singular_program,
            text=True,
            capture_output=True,
            check=True,
        )
        assert "INTEGRAL_BASIS_PASS" in completed.stdout
        assert "?" not in completed.stderr

    @lru_cache(maxsize=None)
    def _reduce_z_power(self, power: int) -> tuple[sp.Expr, ...]:
        if power < self.degree:
            answer = [sp.Integer(0)] * self.degree
            answer[power] = sp.Integer(1)
            return tuple(answer)
        previous = sp.Matrix(self._reduce_z_power(power - 1))
        return tuple(sp.expand(entry) for entry in self.z_matrix * previous)

    @lru_cache(maxsize=None)
    def _multiply_c_basis(
        self, left: int, right: int
    ) -> tuple[sp.Expr, ...]:
        h = self.degree
        left_epsilon, left_z = divmod(left, h)
        right_epsilon, right_z = divmod(right, h)
        if left_epsilon + right_epsilon >= 2:
            return (sp.Integer(0),) * self.c_dimension
        answer = [sp.Integer(0)] * self.c_dimension
        offset = (left_epsilon + right_epsilon) * h
        for z_power, coefficient in enumerate(
            self._reduce_z_power(left_z + right_z)
        ):
            answer[offset + z_power] = coefficient
        return tuple(answer)

    @staticmethod
    @lru_cache(maxsize=None)
    def _distinct_permutations(
        orbit: tuple[int, ...]
    ) -> tuple[tuple[int, ...], ...]:
        return tuple(set(permutations(orbit)))

    @lru_cache(maxsize=None)
    def _multiply_orbit_basis(
        self, left: int, right: int
    ) -> tuple[sp.Expr, ...]:
        ordered_coefficients: dict[tuple[int, ...], sp.Expr] = {}
        for left_ordering in self._distinct_permutations(self.orbits[left]):
            for right_ordering in self._distinct_permutations(self.orbits[right]):
                factor_expansions = [
                    self._multiply_c_basis(a, b)
                    for a, b in zip(left_ordering, right_ordering)
                ]
                nonzero_terms = [
                    tuple(
                        (index, coefficient)
                        for index, coefficient in enumerate(expansion)
                        if coefficient
                    )
                    for expansion in factor_expansions
                ]
                for selected_terms in product(*nonzero_terms):
                    output_indices = tuple(index for index, _ in selected_terms)
                    coefficient = sp.prod(
                        coefficient for _, coefficient in selected_terms
                    )
                    if coefficient:
                        ordered_coefficients[output_indices] = (
                            ordered_coefficients.get(output_indices, 0)
                            + coefficient
                        )

        answer = [sp.Integer(0)] * len(self.orbits)
        for orbit, index in self.orbit_index.items():
            answer[index] = sp.expand(ordered_coefficients.get(orbit, 0))
        return tuple(answer)

    def reduce(self, vector: sp.Matrix) -> sp.Matrix:
        pivot_part = sp.Matrix([vector[index] for index in self.pivot_indices])
        free_part = sp.Matrix([vector[index] for index in self.preferred_free])
        return sp.Matrix(
            [
                sp.factor(entry)
                for entry in free_part - self.reduction_block.T * pivot_part
            ]
        )

    def multiply(self, left: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
        answer = sp.zeros(len(self.preferred_free), 1)
        for left_position, left_coefficient in enumerate(left):
            if not left_coefficient:
                continue
            left_index = self.preferred_free[left_position]
            for right_position, right_coefficient in enumerate(right):
                if not right_coefficient:
                    continue
                right_index = self.preferred_free[right_position]
                orbit_product = sp.Matrix(
                    self._multiply_orbit_basis(left_index, right_index)
                )
                answer += (
                    left_coefficient
                    * right_coefficient
                    * self.reduce(orbit_product)
                )
        return sp.Matrix([sp.factor(entry) for entry in answer])

    def basis_vector(self, position: int) -> sp.Matrix:
        return sp.eye(len(self.preferred_free))[:, position]

    def power(self, element: sp.Matrix, exponent: int) -> sp.Matrix:
        answer = self.basis_vector(0)
        for _ in range(exponent):
            answer = self.multiply(answer, element)
        return answer


def assert_zero(vector: sp.Matrix) -> None:
    assert all(sp.factor(entry) == 0 for entry in vector)


def verify_rees_cones_with_singular() -> None:
    """Certify cubic and quartic Rees saturations."""

    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for the Rees certificate"
    program = r'''
LIB "elim.lib";

ring R3=0,(T,Y,X,h,p,q,r),dp;
ideal J3=
T^3*h-6*T^2*r^2-6*X^2*p*r^3+12*X*Y*r^3,
T^2*X*h-2*T^2*p*r-2*X^2*p^2*r^2+4*X*Y*p*r^2,
-2*T^2*q+T*X^2*h-2*X^2*p*q*r+4*X*Y*q*r,
-4*T^2*p+X^2*Y*h-4*X^2*p^2*r+8*X*Y*p*r,
-6*T^2+X^3*h-6*X^2*p*r+12*X*Y*r,
2*T*Y+X^2*p*q-X^2*r-2*X*Y*q,
2*T*X+X^2*p^2-X^2*q-2*X*Y*p+Y^2;
ideal K3=sat(J3,ideal(h));
ideal C3=subst(K3,h,0);
ideal F3=
X^4,
2*p*X^3-3*X^2*Y,
3*X^2*T-q*X^3,
Y^2+2*T*X-2*p*X*Y+(p^2-q)*X^2,
2*T*Y-2*q*X*Y+(p*q-r)*X^2,
T^2-2*r*X*Y+p*r*X^2;
ideal GC3=std(C3);
ideal GF3=std(F3);
int ok3=1;
ideal R3a=reduce(C3,GF3);
ideal R3b=reduce(F3,GC3);
int i;
for (i=1;i<=size(R3a);i++) { if (R3a[i]!=0) { ok3=0; } }
for (i=1;i<=size(R3b);i++) { if (R3b[i]!=0) { ok3=0; } }

ring R4=0,(D,C,B,A,h),dp;
ideal I4=
A^2*C*D,
D^3,
C*D^2,
B*D^2,
D*(A*C-2*B^2),
B^3-3*A^2*D,
A*D^2,
A*B*D,
A*(A*C+B^2),
A^2*B-2*D^2,
A^3-12*C*D,
C^2+2*B*D,
A*D+B*C;
ideal J4=subst(I4,D,h*D);
J4=subst(J4,C,h*C);
J4=subst(J4,B,h*B);
J4=subst(J4,A,h*A);
J4[1]=J4[1]/h^4;
for (i=2;i<=9;i++) { J4[i]=J4[i]/h^3; }
for (i=10;i<=13;i++) { J4[i]=J4[i]/h^2; }
ideal K4=sat(J4,ideal(h));
ideal C4=subst(K4,h,0);
ideal F4=
C*B+D*A,
C^2+2*D*B,
D*C,
D^2,
B^2*A+C*A^2,
D*B*A,
B^3-3*D*A^2,
D*B^2,
B*A^3,
C*A^3,
D*A^3,
A^5;
ideal GC4=std(C4);
ideal GF4=std(F4);
int ok4=1;
ideal R4a=reduce(C4,GF4);
ideal R4b=reduce(F4,GC4);
for (i=1;i<=size(R4a);i++) { if (R4a[i]!=0) { ok4=0; } }
for (i=1;i<=size(R4b);i++) { if (R4b[i]!=0) { ok4=0; } }

if ((ok3==1) and (ok4==1)) { "REES_CONES_PASS"; }
'''
    completed = subprocess.run(
        (singular, "-q"),
        input=program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "REES_CONES_PASS" in completed.stdout, (
        completed.stdout + completed.stderr
    )
    assert "?" not in completed.stderr

    # Repeat the quartic elimination over the full monic-quartic base.  It
    # is much smaller to eliminate the coefficients of U triangularly in
    # SymPy and ask Singular only for the Rees saturation.
    z = sp.symbols("z")
    p4, q4, r4, s4 = sp.symbols("p q r s")
    a4, b4, c4, d4, rees = sp.symbols("A B C D h")
    u4 = sp.symbols("u0:12")
    collision = z**4 + p4 * z**3 + q4 * z**2 + r4 * z + s4
    v4 = sp.expand(
        collision**2 + a4 * z**3 + b4 * z**2 + c4 * z + d4
    )
    u_polynomial = z**12 + sum(u4[index] * z**index for index in range(12))
    difference = sp.Poly(sp.expand(u_polynomial**2 - v4**3), z)
    solved: dict[sp.Symbol, sp.Expr] = {}
    for degree, variable in zip(range(23, 11, -1), reversed(u4)):
        equation = sp.expand(
            difference.coeff_monomial(z**degree).subs(solved)
        )
        roots = sp.solve(equation, variable, dict=False)
        assert len(roots) == 1
        solved[variable] = sp.factor(roots[0])
    remaining = [
        sp.factor(difference.coeff_monomial(z**degree).subs(solved))
        for degree in range(11, -1, -1)
    ]
    transfer = sp.groebner(
        remaining,
        d4,
        c4,
        b4,
        a4,
        order="grevlex",
        domain=sp.QQ.poly_ring(p4, q4, r4, s4),
    )
    assert len(transfer.polys) == 13

    transverse = (a4, b4, c4, d4)

    def transverse_degree(term: sp.Expr) -> int:
        powers = term.as_powers_dict()
        return sum(int(powers.get(variable, 0)) for variable in transverse)

    def rees_homogenization(expression: sp.Expr) -> sp.Expr:
        terms = sp.Add.make_args(sp.expand(expression))
        order = min(transverse_degree(term) for term in terms)
        homogenized = sp.expand(
            sum(
                term * rees ** (transverse_degree(term) - order)
                for term in terms
            )
        )
        # Singular parses integral polynomials more robustly than expressions
        # containing trailing rational divisions.
        return sp.Poly(
            homogenized,
            a4,
            b4,
            c4,
            d4,
            rees,
            p4,
            q4,
            r4,
            s4,
            domain=sp.QQ,
        ).clear_denoms()[1].as_expr()

    def singular_expression(expression: sp.Expr) -> str:
        return str(sp.expand(expression)).replace("**", "^")

    rees_generators = ",\n".join(
        singular_expression(rees_homogenization(polynomial.as_expr()))
        for polynomial in transfer.polys
    )
    degree_five_monomials = ",\n".join(
        singular_expression(
            a4**a_exponent
            * b4**b_exponent
            * c4**c_exponent
            * d4**d_exponent
        )
        for a_exponent in range(6)
        for b_exponent in range(6 - a_exponent)
        for c_exponent in range(6 - a_exponent - b_exponent)
        for d_exponent in (5 - a_exponent - b_exponent - c_exponent,)
    )
    relative_program = r'''
LIB "elim.lib";
ring R=0,(D,C,B,A,h,p,q,r,s),dp;
ideal J=REES_GENERATORS;
ideal K=sat(J,ideal(h));
ideal Cone=subst(K,h,0);
ideal Universal=std(Cone);
ideal DegreeFive=DEGREE_FIVE_MONOMIALS;
ideal DegreeFiveRemainders=reduce(DegreeFive,Universal);
int finite_ok=1;
int i;
for (i=1;i<=size(DegreeFiveRemainders);i++)
{
  if (DegreeFiveRemainders[i]!=0) { finite_ok=0; }
}

ring Q=(0,p,q,r,s),(D,C,B,A),dp;
ideal Generic=imap(R,Cone);
Generic=std(Generic);
ideal ExpectedBasis=
A^4,A^3,B*A^2,C*A^2,D*A^2,A^2,B*A,C*A,D*A,A,
B^2,D*B,B,C,D,1;
ideal ReducedBasis=reduce(ExpectedBasis,Generic);
int generic_ok=(vdim(Generic)==16);
for (i=1;i<=size(ExpectedBasis);i++)
{
  if (ReducedBasis[i]!=ExpectedBasis[i]) { generic_ok=0; }
}

setring R;
ideal Special=Cone;
Special=subst(Special,p,0);
Special=subst(Special,q,0);
Special=subst(Special,r,0);
Special=subst(Special,s,0);
ideal ExpectedSpecial=
C*B+D*A,
C^2+2*D*B,
D*C,
D^2,
B^2*A+C*A^2,
D*B*A,
B^3-3*D*A^2,
D*B^2,
B*A^3,
C*A^3,
D*A^3,
A^5;
ideal GS=std(Special);
ideal GE=std(ExpectedSpecial);
ideal Left=reduce(Special,GE);
ideal Right=reduce(ExpectedSpecial,GS);
int special_ok=1;
for (i=1;i<=size(Left);i++) { if (Left[i]!=0) { special_ok=0; } }
for (i=1;i<=size(Right);i++) { if (Right[i]!=0) { special_ok=0; } }

ring L=0,(D,C,B,A,l),(dp(4),dp(1));
map DoubleRootFamily=R,D,C,B,A,0,-1-l,l,0,0;
ideal Family=DoubleRootFamily(Cone);
ideal GF=std(Family);
ideal ExpectedLeading=
C*B,C^2,D*C,D^2,B^2*A,D*B*A,B^3,D*B^2,B*A^3,C*A^3,D*A^3,A^5;
int family_ok=(size(GF)==12);
for (i=1;i<=12;i++)
{
  if (leadmonom(GF[i])!=ExpectedLeading[i]) { family_ok=0; }
}

setring R;
if ((finite_ok==1) and (generic_ok==1) and (special_ok==1))
{
  "RELATIVE_QUARTIC_CONE_PASS";
}
setring L;
if (family_ok==1) { "QUARTIC_DOUBLE_ROOT_FAMILY_PASS"; }
'''
    relative_program = relative_program.replace(
        "REES_GENERATORS", rees_generators
    ).replace("DEGREE_FIVE_MONOMIALS", degree_five_monomials)
    relative_completed = subprocess.run(
        (singular, "-q"),
        input=relative_program,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "RELATIVE_QUARTIC_CONE_PASS" in relative_completed.stdout, (
        relative_completed.stdout + relative_completed.stderr
    )
    assert "QUARTIC_DOUBLE_ROOT_FAMILY_PASS" in relative_completed.stdout, (
        relative_completed.stdout + relative_completed.stderr
    )
    assert "?" not in relative_completed.stderr


# ---------------------------------------------------------------------------
# The universal quadratic collision algebra.
# ---------------------------------------------------------------------------

p, q = sp.symbols("p q")
quadratic = FerrandNorm(
    coefficients=(q, p),
    parameters=(p, q),
    preferred_free_labels=((0, 0), (0, 2), (0, 3), (2, 2)),
)
one2, alpha2, beta2, pair2 = (
    quadratic.basis_vector(index) for index in range(4)
)

assert quadratic.multiply(alpha2, alpha2) == 2 * pair2
assert quadratic.multiply(alpha2, beta2) == -p * pair2
assert quadratic.multiply(beta2, beta2) == 2 * q * pair2
for element in (alpha2, beta2, pair2):
    assert_zero(quadratic.multiply(pair2, element))

# The transfer coordinates are X=alpha and Y=-beta.
X2 = alpha2
Y2 = -beta2
assert_zero(quadratic.power(X2, 3))
assert_zero(2 * quadratic.multiply(X2, Y2) - p * quadratic.power(X2, 2))
assert_zero(quadratic.power(Y2, 2) - q * quadratic.power(X2, 2))


# ---------------------------------------------------------------------------
# The universal cubic Ferrand norm.
# ---------------------------------------------------------------------------

p, q, r = sp.symbols("p q r")
cubic = FerrandNorm(
    coefficients=(r, q, p),
    parameters=(p, q, r),
    preferred_free_labels=(
        (0, 0, 0),
        (0, 0, 3),
        (0, 1, 3),
        (0, 3, 3),
        (1, 1, 3),
        (1, 3, 3),
        (2, 3, 3),
        (3, 3, 3),
    ),
)
one3 = cubic.basis_vector(0)
alpha = cubic.basis_vector(1)
beta = cubic.basis_vector(2)
tau = cubic.basis_vector(4)

alpha2 = cubic.power(alpha, 2)
alpha3 = cubic.power(alpha, 3)
alpha4 = cubic.power(alpha, 4)
alpha_beta = cubic.multiply(alpha, beta)
alpha_tau = cubic.multiply(alpha, tau)
beta2 = cubic.power(beta, 2)
beta_tau = cubic.multiply(beta, tau)
tau2 = cubic.power(tau, 2)

ferrand_relations = (
    alpha4,
    3 * cubic.multiply(alpha2, beta) + 2 * p * alpha3,
    3 * cubic.multiply(alpha2, tau) - q * alpha3,
    beta2
    + 2 * alpha_tau
    + 2 * p * alpha_beta
    + (p**2 - q) * alpha2,
    2 * beta_tau - 2 * q * alpha_beta - (p * q - r) * alpha2,
    tau2 + 2 * r * alpha_beta + p * r * alpha2,
)
for relation in ferrand_relations:
    assert_zero(relation)

# These eight monomials are a basis of the norm quotient.
norm_monomials = (
    one3,
    alpha,
    alpha2,
    alpha3,
    beta,
    alpha_beta,
    tau,
    alpha_tau,
)
norm_basis_matrix = sp.Matrix.hstack(*norm_monomials)
assert sp.factor(norm_basis_matrix.det()) in (12, -12)

a, b, t = sp.symbols("a b t")
ferrand_presentation = (
    a**4,
    3 * a**2 * b + 2 * p * a**3,
    3 * a**2 * t - q * a**3,
    b**2 + 2 * a * t + 2 * p * a * b + (p**2 - q) * a**2,
    2 * b * t - 2 * q * a * b - (p * q - r) * a**2,
    t**2 + 2 * r * a * b + p * r * a**2,
)
coefficient_domain = sp.QQ.poly_ring(p, q, r)
ferrand_groebner = sp.groebner(
    ferrand_presentation, t, b, a, order="grevlex", domain=coefficient_domain
)
assert len(ferrand_groebner.polys) == 6
for monomial in (a**3, a**2, a * b, a * t, a, b, t, sp.Integer(1)):
    assert ferrand_groebner.reduce(monomial)[1] == monomial


# ---------------------------------------------------------------------------
# Independent recomputation of the degree-three square/cube transfer block.
# ---------------------------------------------------------------------------

Z = sp.symbols("Z")
U_coefficients = sp.symbols("u0:9")
S = Z**3 + p * Z**2 + q * Z + r
X, Y, T = sp.symbols("X Y T")
V = sp.expand(S**2 + X * Z**2 + Y * Z + T)
U = Z**9 + sum(
    U_coefficients[index] * Z**index for index in range(9)
)
difference = sp.Poly(sp.expand(U**2 - V**3), Z)

solution: dict[sp.Symbol, sp.Expr] = {}
for degree, variable in zip(
    range(17, 8, -1), reversed(U_coefficients)
):
    equation = sp.expand(
        difference.coeff_monomial(Z**degree).subs(solution)
    )
    roots = sp.solve(equation, variable, dict=False)
    assert len(roots) == 1
    solution[variable] = sp.factor(roots[0])

remaining = [
    sp.factor(difference.coeff_monomial(Z**degree).subs(solution))
    for degree in range(8, -1, -1)
]
transfer_relations = (
    T**3 - 6 * T**2 * r**2 - 6 * X**2 * p * r**3 + 12 * X * Y * r**3,
    T**2 * X
    - 2 * T**2 * p * r
    - 2 * X**2 * p**2 * r**2
    + 4 * X * Y * p * r**2,
    T * X**2
    - 2 * T**2 * q
    - 2 * X**2 * p * q * r
    + 4 * X * Y * q * r,
    X**2 * Y
    - 4 * T**2 * p
    - 4 * X**2 * p**2 * r
    + 8 * X * Y * p * r,
    X**3 - 6 * T**2 - 6 * X**2 * p * r + 12 * X * Y * r,
    2 * T * Y + X**2 * p * q - X**2 * r - 2 * X * Y * q,
    Y**2 + 2 * T * X + X**2 * p**2 - X**2 * q - 2 * X * Y * p,
)
transfer_groebner = sp.groebner(
    transfer_relations, T, Y, X, order="grevlex", domain=coefficient_domain
)
affine_groebner = sp.groebner(
    remaining[:7], T, Y, X, order="grevlex", domain=coefficient_domain
)
assert all(
    transfer_groebner.reduce(equation)[1] == 0
    for equation in remaining
)
assert all(
    affine_groebner.reduce(equation)[1] == 0
    for equation in transfer_relations
)
transfer_basis = (
    sp.Integer(1), T, T**2, X, T * X, X**2, Y, X * Y
)
assert all(
    transfer_groebner.reduce(monomial)[1] == monomial
    for monomial in transfer_basis
)

# Rees saturation is checked below in Singular.  Under the cofactor-sign map
# alpha=X, beta=-Y, tau=T, the six cubic Ferrand relations are exactly the
# defining equations of the relative I=(X,Y,T) normal cone of Z_3.
relative_cubic_cone_relations = tuple(
    sp.expand(relation.subs({a: X, b: -Y, t: T}))
    for relation in ferrand_presentation
)


# ---------------------------------------------------------------------------
# The cubic norm and transfer fibers are not isomorphic.
# ---------------------------------------------------------------------------

special_norm_relations = tuple(
    relation.subs({p: 0, q: 0, r: 0})
    for relation in ferrand_presentation
)
special_norm = sp.groebner(
    special_norm_relations, t, b, a, order="grevlex", domain=sp.QQ
)
assert special_norm.reduce(t**2)[1] == 0
assert special_norm.reduce(t)[1] == t

special_transfer_relations = tuple(
    relation.subs({p: 0, q: 0, r: 0})
    for relation in transfer_relations
)
special_transfer = sp.groebner(
    special_transfer_relations, T, Y, X, order="grevlex", domain=sp.QQ
)
assert special_transfer.reduce(T**2)[1] == T**2
assert special_transfer.reduce(X**3 - 6 * T**2)[1] == 0

# The Ferrand fiber is the specialization of the relative normal-cone
# isomorphism.  These are the initial forms for m=(X,Y,T); X^4 is the extra
# initial relation obtained from X(X^3-6T^2) and T^2X.
tangent_cone_relations = (
    T**2,
    T * Y,
    Y**2 + 2 * T * X,
    T * X**2,
    X**2 * Y,
    X**4,
)
tangent_cone = sp.groebner(
    tangent_cone_relations, T, Y, X, order="grevlex", domain=sp.QQ
)
renamed_special_norm = tuple(
    relation.subs({a: X, b: Y, t: T})
    for relation in special_norm_relations
)
assert all(
    tangent_cone.reduce(relation)[1] == 0
    for relation in renamed_special_norm
)
assert all(
    special_norm.reduce(
        relation.subs({X: a, Y: b, T: t})
    )[1]
    == 0
    for relation in tangent_cone_relations
)

# Modulo m^3, a square of aX+bY+cT has coefficients
# a^2, 2ab, 2(ac-b^2) in the basis X^2,XY,TX.  Hence a=b=0.
# Any possible square-zero tangent lift is therefore cT+h with h in m^2.
# But T annihilates m^2, m^4=0, and T^2 is nonzero.
for monomial in (X**2, X * Y, T * X, T**2):
    assert special_transfer.reduce(T * monomial)[1] == 0
for left in (X**2, X * Y, T * X, T**2):
    for right in (X**2, X * Y, T * X, T**2):
        assert special_transfer.reduce(left * right)[1] == 0


# ---------------------------------------------------------------------------
# The maximally collided quartic norm and transfer tangent cone.
# ---------------------------------------------------------------------------

A, B, C, D = sp.symbols("A B C D")
quartic_u = sp.symbols("w0:12")
quartic_V = Z**8 + A * Z**3 + B * Z**2 + C * Z + D
quartic_U = Z**12 + sum(
    quartic_u[index] * Z**index for index in range(12)
)
quartic_difference = sp.Poly(sp.expand(quartic_U**2 - quartic_V**3), Z)
quartic_solution: dict[sp.Symbol, sp.Expr] = {}
for degree, variable in zip(range(23, 11, -1), reversed(quartic_u)):
    equation = sp.expand(
        quartic_difference.coeff_monomial(Z**degree).subs(quartic_solution)
    )
    roots = sp.solve(equation, variable, dict=False)
    assert len(roots) == 1
    quartic_solution[variable] = sp.factor(roots[0])
quartic_remaining = [
    sp.factor(
        quartic_difference.coeff_monomial(Z**degree).subs(quartic_solution)
    )
    for degree in range(11, -1, -1)
]
quartic_transfer_relations = (
    A**2 * C * D,
    D**3,
    C * D**2,
    B * D**2,
    D * (A * C - 2 * B**2),
    B**3 - 3 * A**2 * D,
    A * D**2,
    A * B * D,
    A * (A * C + B**2),
    A**2 * B - 2 * D**2,
    A**3 - 12 * C * D,
    C**2 + 2 * B * D,
    A * D + B * C,
)
quartic_transfer = sp.groebner(
    quartic_transfer_relations, D, C, B, A, order="grevlex", domain=sp.QQ
)
quartic_affine = sp.groebner(
    quartic_remaining[:10], D, C, B, A, order="grevlex", domain=sp.QQ
)
assert all(
    quartic_transfer.reduce(equation)[1] == 0
    for equation in quartic_remaining
)
assert all(
    quartic_affine.reduce(equation)[1] == 0
    for equation in quartic_transfer_relations
)

quartic_preferred_labels = (
    (0, 0, 0, 0),
    (0, 0, 0, 4),
    (0, 0, 1, 4),
    (0, 0, 4, 4),
    (0, 1, 1, 4),
    (0, 1, 4, 4),
    (0, 2, 4, 4),
    (0, 4, 4, 4),
    (1, 1, 1, 4),
    (1, 1, 4, 4),
    (1, 2, 4, 4),
    (1, 4, 4, 4),
    (2, 2, 4, 4),
    (2, 4, 4, 4),
    (3, 4, 4, 4),
    (4, 4, 4, 4),
)
quartic = FerrandNorm(
    coefficients=(sp.Integer(0),) * 4,
    parameters=(),
    preferred_free_labels=quartic_preferred_labels,
)
one4 = quartic.basis_vector(0)
quartic_A = quartic.basis_vector(1)
quartic_B = -quartic.basis_vector(2)
quartic_C = quartic.basis_vector(4)
quartic_D = -quartic.basis_vector(8)


def quartic_power(element: sp.Matrix, exponent: int) -> sp.Matrix:
    return quartic.power(element, exponent)


quartic_cone_vectors = (
    quartic.multiply(quartic_C, quartic_B)
    + quartic.multiply(quartic_D, quartic_A),
    quartic_power(quartic_C, 2)
    + 2 * quartic.multiply(quartic_D, quartic_B),
    quartic.multiply(quartic_D, quartic_C),
    quartic_power(quartic_D, 2),
    quartic.multiply(quartic_power(quartic_B, 2), quartic_A)
    + quartic.multiply(quartic_C, quartic_power(quartic_A, 2)),
    quartic.multiply(
        quartic.multiply(quartic_D, quartic_B), quartic_A
    ),
    quartic_power(quartic_B, 3)
    - 3 * quartic.multiply(quartic_D, quartic_power(quartic_A, 2)),
    quartic.multiply(quartic_D, quartic_power(quartic_B, 2)),
    quartic.multiply(quartic_B, quartic_power(quartic_A, 3)),
    quartic.multiply(quartic_C, quartic_power(quartic_A, 3)),
    quartic.multiply(quartic_D, quartic_power(quartic_A, 3)),
    quartic_power(quartic_A, 5),
)
for relation in quartic_cone_vectors:
    assert_zero(relation)

quartic_norm_monomials = (
    one4,
    quartic_power(quartic_A, 4),
    quartic_power(quartic_A, 3),
    quartic.multiply(quartic_B, quartic_power(quartic_A, 2)),
    quartic.multiply(quartic_C, quartic_power(quartic_A, 2)),
    quartic.multiply(quartic_D, quartic_power(quartic_A, 2)),
    quartic_power(quartic_A, 2),
    quartic.multiply(quartic_B, quartic_A),
    quartic.multiply(quartic_C, quartic_A),
    quartic.multiply(quartic_D, quartic_A),
    quartic_A,
    quartic_power(quartic_B, 2),
    quartic.multiply(quartic_D, quartic_B),
    quartic_B,
    quartic_C,
    quartic_D,
)
quartic_norm_basis_matrix = sp.Matrix.hstack(*quartic_norm_monomials)
assert abs(quartic_norm_basis_matrix.det()) == 4608

quartic_cone_relations = (
    C * B + D * A,
    C**2 + 2 * D * B,
    D * C,
    D**2,
    B**2 * A + C * A**2,
    D * B * A,
    B**3 - 3 * D * A**2,
    D * B**2,
    B * A**3,
    C * A**3,
    D * A**3,
    A**5,
)
quartic_cone = sp.groebner(
    quartic_cone_relations, D, C, B, A, order="grevlex", domain=sp.QQ
)
quartic_cone_basis = (
    sp.Integer(1),
    A,
    B,
    C,
    D,
    A**2,
    A * B,
    B**2,
    A * C,
    A * D,
    B * D,
    A**3,
    A**2 * B,
    A**2 * C,
    A**2 * D,
    A**4,
)
assert all(
    quartic_cone.reduce(monomial)[1] == monomial
    for monomial in quartic_cone_basis
)


def verify_quartic_collision_strata() -> None:
    """Exact optional certificate on the moduli-bearing double-root stratum.

    The normal form S=z^2(z-1)(z-lambda) covers every quartic root type
    2+1+1 up to affine change of z.  Its lambda=0 and lambda=1 fibers cover
    3+1 and 2+2.  The 4 fiber is certified above, while the split locus is
    the formal tensor-product case.
    """

    lam = sp.symbols("lambda")
    family = FerrandNorm(
        coefficients=(0, 0, lam, -(1 + lam)),
        parameters=(lam,),
        preferred_free_labels=quartic_preferred_labels,
    )
    family_one = family.basis_vector(0)
    family_a = family.basis_vector(1)
    family_b = -family.basis_vector(2)
    family_c = family.basis_vector(4)
    family_d = -family.basis_vector(8)

    def family_power(element: sp.Matrix, exponent: int) -> sp.Matrix:
        return family.power(element, exponent)

    family_monomials = (
        family_power(family_a, 4),
        family_power(family_a, 3),
        family.multiply(family_b, family_power(family_a, 2)),
        family.multiply(family_c, family_power(family_a, 2)),
        family.multiply(family_d, family_power(family_a, 2)),
        family_power(family_a, 2),
        family.multiply(family_b, family_a),
        family.multiply(family_c, family_a),
        family.multiply(family_d, family_a),
        family_a,
        family_power(family_b, 2),
        family.multiply(family_d, family_b),
        family_b,
        family_c,
        family_d,
        family_one,
    )
    determinant = sp.factor(sp.Matrix.hstack(*family_monomials).det())
    assert determinant in (sp.Integer(4608), sp.Integer(-4608))


verify_rees_cones_with_singular()

print("PASS: the quadratic Ferrand norm is exactly the universal Z_2 block")
print("PASS: the universal cubic Ferrand norm is free of rank 8 with six relations")
print("PASS: the independently eliminated transfer Z_3 is free of rank 8 with seven relations")
print("PASS: their triple-root fibers are not isomorphic (square-zero tangent obstruction)")
print("PASS: the full relative cubic Ferrand norm is gr_(X,Y,T)(Z_3)")
print("PASS: the quadruple-root Ferrand norm is the rank-16 tangent cone gr_m(Z_4)")
print("PASS: the relative quartic transfer cone is finite with the correct collision-stratum bases")
if "--quartic-strata" in sys.argv:
    verify_quartic_collision_strata()
    print("PASS: the quartic double-root family has an integral cofactor basis")
