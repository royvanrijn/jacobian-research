#!/usr/bin/env python3
"""Verify generic Schur rigidity on a two-parameter sextic surface.

The sextic is

  (x^6+y^6+z^6)/30
  + mu*x^2*y^2*z^2
  + nu*sum_{i != j} x_i^4*x_j^2.

On nu != 0, six boundary coefficients solve the quadratic Schur quotient.
After clearing the resulting nu^2 denominator, 114 intrinsic equations
remain in the fifteen quartic coefficients.  Over Q(mu,nu), their exact
Groebner basis makes every quartic coefficient nilpotent of exponent three.

This proves rigidity at the generic point of the parameter surface.  It
does not claim that the specialization locus inside nu != 0 is empty.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess

import sympy as sp


x, y, z, mu, nu = sp.symbols("x y z mu nu")
quartic_coefficients = sp.symbols("s0:15")

quartic_monomials = [
    x**i * y**j * z ** (4 - i - j)
    for i in range(5)
    for j in range(5 - i)
]
quartic = sum(
    coefficient * monomial
    for coefficient, monomial in zip(
        quartic_coefficients, quartic_monomials
    )
)

mixed_42 = sum(
    left**4 * right**2
    for left in (x, y, z)
    for right in (x, y, z)
    if left != right
)
h6 = (
    (x**6 + y**6 + z**6) / 30
    + mu * x**2 * y**2 * z**2
    + nu * mixed_42
)
hessian = sp.hessian(h6, (x, y, z))
hessian_determinant = sp.expand(hessian.det())
gradient = sp.Matrix(
    [sp.diff(quartic, variable) for variable in (x, y, z)]
)
schur_numerator = sp.expand(
    (gradient.T * hessian.adjugate() * gradient)[0]
)

quadratic_monomials = (
    x**2,
    y**2,
    z**2,
    x * y,
    x * z,
    y * z,
)
degree_14_monomials = [
    x**i * y**j * z ** (14 - i - j)
    for i in range(15)
    for j in range(15 - i)
]
quotient_matrix = sp.Matrix(
    [
        [
            sp.Poly(
                hessian_determinant * quadratic_monomial,
                x,
                y,
                z,
            ).coeff_monomial(degree_14_monomial)
            for quadratic_monomial in quadratic_monomials
        ]
        for degree_14_monomial in degree_14_monomials
    ]
)
numerator_vector = sp.Matrix(
    [
        sp.Poly(schur_numerator, x, y, z).coeff_monomial(
            monomial
        )
        for monomial in degree_14_monomials
    ]
)

# These are z^14, y*z^13, y^2*z^12, x*z^13, x*y*z^12,
# and x^2*z^12.  Their quotient matrix has determinant 4096*nu^12.
pivot_rows = (0, 1, 2, 15, 16, 29)
pivot_matrix = quotient_matrix[list(pivot_rows), :]
assert sp.factor(pivot_matrix.det()) == 4096 * nu**12

s = quartic_coefficients

# To avoid repeated rational-matrix simplification, store 2*nu^2 times
# the six exact quotient solutions.
quotient_numerators = (
    -mu * s[1] ** 2
    + 640 * nu**3 * s[0] ** 2
    - 96 * nu**2 * s[0] * s[9]
    + 2 * nu**2 * s[5] ** 2
    + 2 * nu * s[1] * s[10]
    + 6 * nu * s[12] * s[5]
    - 6 * nu * s[5] ** 2
    + nu * s[6] ** 2
    + 4 * nu * s[9] ** 2,
    -mu * s[5] ** 2
    + 640 * nu**3 * s[0] ** 2
    - 96 * nu**2 * s[0] * s[2]
    + 2 * nu**2 * s[1] ** 2
    - 6 * nu * s[1] ** 2
    + 6 * nu * s[1] * s[3]
    + 4 * nu * s[2] ** 2
    + 2 * nu * s[5] * s[7]
    + nu * s[6] ** 2,
    nu * (32 * nu * s[0] ** 2 + s[1] ** 2 + s[5] ** 2),
    -4
    * (
        mu * s[1] * s[5]
        + 24 * nu**2 * s[0] * s[6]
        - nu**2 * s[1] * s[5]
        - nu * s[1] * s[7]
        - nu * s[10] * s[5]
        - nu * s[2] * s[6]
        - nu * s[6] * s[9]
    ),
    -2
    * nu
    * (8 * nu * s[0] * s[5] - s[1] * s[6] - 2 * s[5] * s[9]),
    -2
    * nu
    * (8 * nu * s[0] * s[1] - 2 * s[1] * s[2] - s[5] * s[6]),
)

for pivot_position, row in enumerate(pivot_rows):
    pivot_remainder = sp.expand(
        2 * nu**2 * numerator_vector[row]
        - sum(
            quotient_matrix[row, column]
            * quotient_numerators[column]
            for column in range(6)
        )
    )
    assert pivot_remainder == 0, pivot_position

intrinsic_equations = []
for row in range(len(degree_14_monomials)):
    equation = sp.expand(
        2 * nu**2 * numerator_vector[row]
        - sum(
            quotient_matrix[row, column]
            * quotient_numerators[column]
            for column in range(6)
        )
    )
    if equation != 0:
        intrinsic_equations.append(equation)
assert len(intrinsic_equations) == 114


def singular_expression(expression: sp.Expr) -> str:
    return str(sp.expand(expression)).replace("**", "^")


singular = shutil.which("Singular")
if singular is None:
    raise RuntimeError("Singular is required for the exact generic check")

parser = argparse.ArgumentParser()
parser.add_argument(
    "--extract-denominators",
    action="store_true",
    help="print the factored coefficient denominator of the generic basis",
)
parser.add_argument(
    "--basis-profile",
    action="store_true",
    help="print the generic quotient dimension and degree profile",
)
parser.add_argument(
    "--fitting-profile",
    action="store_true",
    help="compute the degree-two fraction-free Fitting pivot",
)
arguments = parser.parse_args()

if arguments.fitting_profile:
    from sympy.polys.matrices import DomainMatrix

    quadratic_coefficient_monomials = [
        s[left] * s[right]
        for left in range(len(s))
        for right in range(left, len(s))
    ]
    quadratic_coefficient_matrix = sp.Matrix(
        [
            [
                sp.Poly(
                    equation,
                    *s,
                ).coeff_monomial(monomial)
                for monomial in quadratic_coefficient_monomials
            ]
            for equation in intrinsic_equations
        ]
    )
    row_contents = []
    for row in range(quadratic_coefficient_matrix.rows):
        content = sp.Integer(0)
        for entry in quadratic_coefficient_matrix.row(row):
            content = sp.gcd(content, entry)
        content = sp.factor(content)
        row_contents.append(content)
        if content != 1:
            quadratic_coefficient_matrix.row_op(
                row, lambda entry, _: sp.cancel(entry / content)
            )
    content_counts = {
        content: row_contents.count(content)
        for content in set(row_contents)
    }
    print(
        "ROW_CONTENTS:",
        ", ".join(
            f"{sp.sstr(content)}:{multiplicity}"
            for content, multiplicity in sorted(
                content_counts.items(), key=lambda item: sp.sstr(item[0])
            )
        ),
    )
    polynomial_domain = sp.QQ.poly_ring(mu, nu)
    fitting_matrix = DomainMatrix.from_Matrix(
        quadratic_coefficient_matrix
    ).convert_to(polynomial_domain)
    fitting_columns = ()
    fitting_rows = ()
    fitting_point = None
    for candidate_mu, candidate_nu in (
        (1, 1),
        (1, 2),
        (2, 1),
        (2, 3),
        (0, 1),
    ):
        numerical_matrix = DomainMatrix.from_Matrix(
            quadratic_coefficient_matrix.subs(
                {mu: candidate_mu, nu: candidate_nu}
            )
        ).convert_to(sp.QQ)
        _, candidate_columns = numerical_matrix.rref()
        print(
            "SPECIALIZED_RANK:",
            candidate_mu,
            candidate_nu,
            len(candidate_columns),
        )
        if len(candidate_columns) > len(fitting_columns):
            fitting_columns = candidate_columns
            fitting_point = (candidate_mu, candidate_nu)
            selected_columns = numerical_matrix[:, list(candidate_columns)]
            _, fitting_rows = selected_columns.transpose().rref()
    assert len(fitting_columns) == 99
    assert len(fitting_rows) == 99
    fitting_minor = fitting_matrix[
        list(fitting_rows), list(fitting_columns)
    ]
    fitting_determinant = fitting_minor.det().as_expr()
    print(
        "DEGREE_TWO_FITTING:",
        len(fitting_columns),
        sp.factor(fitting_determinant),
    )
    raise SystemExit

denominator_program = ""
standard_basis_program = "ideal G=slimgb(I);"
if arguments.extract_denominators:
    standard_basis_program = """
matrix transformation;
ideal G=liftstd(I,transformation,"slimgb");
"""
    all_quartic_variables = "*".join(map(str, quartic_coefficients))
    denominator_program = f"""
matrix coefficient_matrix;
int coefficient_index;
number coefficient_denominator;
for (i=1;i<=size(G);i++)
{{
  coefficient_matrix=coef(G[i],{all_quartic_variables});
  for (
    coefficient_index=1;
    coefficient_index<=ncols(coefficient_matrix);
    coefficient_index++
  )
  {{
    coefficient_denominator=denominator(
      leadcoef(coefficient_matrix[2,coefficient_index])
    );
    print("DENOMINATOR "+string(coefficient_denominator));
  }}
}}
int transformation_row;
int transformation_column;
poly transformation_entry;
for (
  transformation_row=1;
  transformation_row<=nrows(transformation);
  transformation_row++
)
{{
  for (
    transformation_column=1;
    transformation_column<=ncols(transformation);
    transformation_column++
  )
  {{
    transformation_entry=transformation[
      transformation_row,transformation_column
    ];
    if (transformation_entry!=0)
    {{
      coefficient_matrix=coef(
        transformation_entry,{all_quartic_variables}
      );
      for (
        coefficient_index=1;
        coefficient_index<=ncols(coefficient_matrix);
        coefficient_index++
      )
      {{
        coefficient_denominator=denominator(
          leadcoef(coefficient_matrix[2,coefficient_index])
        );
        print("DENOMINATOR "+string(coefficient_denominator));
      }}
    }}
  }}
}}
"""

profile_program = ""
if arguments.basis_profile:
    profile_program = """
ideal standard_monomials=kbase(G);
int maximum_standard_degree=0;
int standard_index;
for (standard_index=1;standard_index<=size(standard_monomials);standard_index++)
{
  if (deg(standard_monomials[standard_index])>maximum_standard_degree)
  {
    maximum_standard_degree=deg(standard_monomials[standard_index]);
  }
}
print(
  "PROFILE "
  +string(vdim(G))+" "
  +string(size(standard_monomials))+" "
  +string(maximum_standard_degree)
);
intvec hilbert_function=hilb(G,1);
print("HILBERT "+string(hilbert_function));
"""

program = f"""
ring rr=(0,mu,nu),({",".join(map(str, quartic_coefficients))}),dp;
option(redSB);
ideal I={",".join(map(singular_expression, intrinsic_equations))};
{standard_basis_program}
ideal M={",".join(map(str, quartic_coefficients))};
ideal GM=std(M);
print(
  "GENERIC "
  +string(size(I))+" "
  +string(size(G))+" "
  +string(size(reduce(I,GM)))
);
int i;
poly cube;
for (i=1;i<=size(M);i++)
{{
  cube=M[i]^3;
  print("CUBE "+string(i)+" "+string(reduce(cube,G)==0));
}}
{denominator_program}
{profile_program}
"""
completed = subprocess.run(
    [singular, "-q"],
    input=program,
    text=True,
    capture_output=True,
    check=True,
    timeout=900 if arguments.extract_denominators else 180,
)
if completed.stderr.strip():
    raise RuntimeError(completed.stderr)

generic_marker = re.search(
    r"(?m)^GENERIC (\d+) (\d+) (\d+)$", completed.stdout
)
assert generic_marker is not None
assert tuple(map(int, generic_marker.groups())) == (114, 117, 0)
cube_markers = re.findall(
    r"(?m)^CUBE (\d+) ([01])$", completed.stdout
)
assert len(cube_markers) == 15
assert all(success == "1" for _, success in cube_markers)

if arguments.extract_denominators:
    denominator_strings = set(
        re.findall(r"(?m)^DENOMINATOR (.+)$", completed.stdout)
    )
    coefficient_denominator = sp.Integer(1)
    for denominator_string in denominator_strings:
        denominator = sp.sympify(
            denominator_string.replace("^", "**"),
            locals={"mu": mu, "nu": nu},
        )
        coefficient_denominator = sp.lcm(
            coefficient_denominator,
            sp.Poly(denominator, mu, nu),
        ).as_expr()
    print(
        "GENERIC_BASIS_DENOMINATOR:",
        sp.factor(coefficient_denominator),
    )
if arguments.basis_profile:
    profile_marker = re.search(
        r"(?m)^PROFILE (\d+) (\d+) (\d+)$", completed.stdout
    )
    assert profile_marker is not None
    print("GENERIC_QUOTIENT_PROFILE:", " ".join(profile_marker.groups()))
    hilbert_marker = re.search(r"(?m)^HILBERT (.+)$", completed.stdout)
    assert hilbert_marker is not None
    print("GENERIC_HILBERT_DATA:", hilbert_marker.group(1))

print("PASS: six quotient pivots have determinant 4096*nu^12")
print("PASS: quotient elimination leaves 114 intrinsic equations")
print("PASS: the generic Groebner basis has 117 elements")
print("PASS: all fifteen quartic coefficient cubes reduce to zero")
print("SCOPE: generic rigidity only; exceptional nu!=0 fibers remain open")
