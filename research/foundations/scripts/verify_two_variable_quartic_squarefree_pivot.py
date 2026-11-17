#!/usr/bin/env python3
"""Exact pivot decomposition of the squarefree binary-quartic rank-one fiber.

For A=u*v*(u-v)*(u-lambda*v), eliminate the first moment and work on the
projective chart e=1.  The second moment is linear in a.  This checker:

* derives its pivot exactly;
* proves over Q(lambda) that the pivot-zero boundary has no common zero of
  moments two through six;
* identifies the expected three-point ideal after pivot elimination as

      (8*c-3*d^2, d*(d-4)*(d-4*lambda));

* verifies that, after imposing 8*c=3*d^2, moments three and four contain
  the displayed cubic factor, while moments five and six contain its square.
* proves the generic radical statement by a resultant decomposition.  The
  pairwise resultant gcd is the expected cubic times a sixth power of one
  quadratic.  The quadratic component is supported on p=0 because p^3
  reduces to zero there, while the three linear d-branches force the
  expected c-values.

Thus the generic squarefree six-moment fiber consists exactly of the
annihilator sections.  Exceptional lambda-specializations are not
classified here; collided-root strata are covered by the separate
two-root and three-root theorems.

The checker also closes two distinguished algebraic specializations:
lambda^2+4*lambda+1=0, where the pivot meets an annihilator section, and
lambda^2-lambda+1=0, the equianharmonic orbit.  The harmonic orbit is
covered by the separate lambda=2 anchor checker.
"""

from __future__ import annotations

import json
from math import factorial
from pathlib import Path
import shutil
import subprocess

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_variable_quartic_squarefree_pivot.json"
)

x, y, u, v = sp.symbols("x y u v")
a, b, c, d, e, lam = sp.symbols("a b c d e lam")
P = a * x**4 + b * x**3 * y + c * x**2 * y**2 + d * x * y**3 + e * y**4
A = u * v * (u - v) * (u - lam * v)


def apolar_moment(order: int) -> sp.Expr:
    symbol_power = sp.Poly(sp.expand(A**order), u, v)
    polynomial_power = sp.Poly(sp.expand(P**order), x, y)
    return sp.expand(
        sum(
            coefficient
            * polynomial_power.coeff_monomial(x**x_order * y**y_order)
            * factorial(x_order)
            * factorial(y_order)
            for (x_order, y_order), coefficient in symbol_power.terms()
        )
    )


def singular_expression(expression: sp.Expr) -> str:
    return sp.sstr(expression).replace("**", "^")


def cleared_polynomial(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> sp.Expr:
    numerator = sp.fraction(sp.cancel(expression))[0]
    polynomial = sp.Poly(numerator, *variables, domain=sp.QQ)
    return sp.expand(polynomial.clear_denoms()[1].primitive()[1].as_expr())


def main() -> None:
    moments = [apolar_moment(order) for order in range(1, 7)]
    assert moments[0] == 6 * b - 4 * c * lam - 4 * c + 6 * d * lam
    b_solution = sp.Rational(2, 3) * c * (lam + 1) - d * lam
    chart_moments = [
        sp.cancel(moment.subs(b, b_solution).subs(e, 1))
        for moment in moments[1:]
    ]

    second = sp.Poly(chart_moments[0], a)
    pivot = sp.factor(second.coeff_monomial(a) / 576)
    constant = second.coeff_monomial(1)
    assert pivot == 5 * c - 5 * d * lam - 5 * d + 2 * lam**2 + 8 * lam + 2

    # Boundary pivot=0: solve for c, retain the constant part of mu_2,
    # and adjoin mu_3,...,mu_6.  The resulting ideal in Q(lam)[a,d] is
    # the unit ideal.
    c_boundary = sp.solve(pivot, c)[0]
    boundary_polynomials = [
        cleared_polynomial(
            constant.subs(c, c_boundary),
            (a, d, lam),
        )
    ]
    boundary_polynomials.extend(
        cleared_polynomial(
            moment.subs(c, c_boundary),
            (a, d, lam),
        )
        for moment in chart_moments[1:]
    )
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    boundary_program = f"""
ring r=(0,lam),(a,d),dp;
option(redSB);
ideal I={",".join(singular_expression(value) for value in boundary_polynomials)};
ideal G=std(I);
if (size(G)!=1 || G[1]!=1)
{{
  print("BOUNDARY_NONUNIT");
  exit(1);
}}
print("BOUNDARY_UNIT");
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=boundary_program,
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    assert completed.stdout.strip() == "BOUNDARY_UNIT", completed.stdout

    # The unit certificate can specialize only where the following
    # necessary resultant gcd vanishes.  The squarefree part consists of
    # one S_3-orbit; lambda=0,1 are collided-root values.
    pivot_orbit_factors = (
        lam**2 + 4 * lam + 1,
        lam**2 - 6 * lam + 6,
        6 * lam**2 - 6 * lam + 1,
    )
    pivot_exception_polynomial = (
        lam**4
        * (lam - 1) ** 4
        * sp.prod(pivot_orbit_factors)
    )
    boundary_exception_program = """
ring r=0,(a,d,lam),lp;
"""
    for order, polynomial in enumerate(boundary_polynomials, start=2):
        boundary_exception_program += (
            f"poly f{order}={singular_expression(polynomial)};\n"
        )
    boundary_exception_program += f"""
poly r34=resultant(f3,f4,a);
poly r35=resultant(f3,f5,a);
poly r36=resultant(f3,f6,a);
poly s34=resultant(f2,r34,d);
poly s35=resultant(f2,r35,d);
poly s36=resultant(f2,r36,d);
poly exceptional=gcd(gcd(s34,s35),s36);
ideal expected={singular_expression(pivot_exception_polynomial)};
if (deg(exceptional)!=14 || reduce(exceptional,expected)!=0)
{{
  print("BAD_BOUNDARY_EXCEPTIONAL_GCD");
  exit(1);
}}
print("BOUNDARY_EXCEPTIONAL_GCD");
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=boundary_exception_program,
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    assert completed.stdout.strip() == "BOUNDARY_EXCEPTIONAL_GCD", (
        completed.stdout
    )

    # Principal-open pivot: eliminate a with mu_2.  The three expected
    # points in (c,d) are (0,0), (6,4), and (6*lam^2,4*lam).
    a_solution = sp.cancel(-constant / second.coeff_monomial(a))
    expected_relation = 8 * c - 3 * d**2
    expected_cubic = d * (d - 4) * (d - 4 * lam)
    expected_points = (
        {c: 0, d: 0},
        {c: 6, d: 4},
        {c: 6 * lam**2, d: 4 * lam},
    )
    for point in expected_points:
        assert sp.expand(expected_relation.subs(point)) == 0
        assert sp.expand(expected_cubic.subs(point)) == 0

    multiplicities: list[int] = []
    for order, moment in enumerate(chart_moments[1:], start=3):
        restricted = sp.factor(
            sp.cancel(
                moment.subs(a, a_solution).subs(
                    c,
                    sp.Rational(3, 8) * d**2,
                )
            )
        )
        numerator = sp.Poly(sp.fraction(restricted)[0], d, domain=sp.QQ[lam])
        multiplicity = 0
        divisor = sp.Poly(expected_cubic, d, domain=sp.QQ[lam])
        while numerator.rem(divisor).is_zero:
            numerator = numerator.exquo(divisor)
            multiplicity += 1
        multiplicities.append(multiplicity)
        assert multiplicity >= (1 if order <= 4 else 2)
    assert multiplicities == [1, 1, 2, 2]

    # Generic radical decomposition.  Any common (c,d)-zero makes every
    # pairwise c-resultant vanish.  Their gcd over Q(lambda)[d] is the
    # expected cubic times q(d)^6.  The q-component is supported on p=0,
    # while the three expected d-branches force c=0,6,6*lambda^2.
    extra_quadratic = (
        (25 * lam**4 - 50 * lam**3 + 75 * lam**2 - 50 * lam + 25)
        * d**2
        + (
            100 * lam**5
            - 350 * lam**4
            + 100 * lam**3
            + 100 * lam**2
            - 350 * lam
            + 100
        )
        * d
        - 44 * lam**6
        - 68 * lam**5
        + 386 * lam**4
        + 208 * lam**3
        + 386 * lam**2
        - 68 * lam
        - 44
    )
    eliminated_polynomials = [
        sp.Poly(
            sp.fraction(sp.cancel(moment.subs(a, a_solution)))[0],
            c,
            d,
            domain=sp.QQ[lam],
        ).primitive()[1].as_expr()
        for moment in chart_moments[1:]
    ]
    decomposition_program = """
ring r=(0,lam),(c,d),dp;
option(redSB);
"""
    for order, polynomial in enumerate(eliminated_polynomials, start=3):
        decomposition_program += (
            f"poly f{order}={singular_expression(polynomial)};\n"
        )
    decomposition_program += f"""
poly p={singular_expression(pivot)};
poly cubic={singular_expression(expected_cubic)};
poly q={singular_expression(extra_quadratic)};
poly common=resultant(f3,f4,c);
common=gcd(common,resultant(f3,f5,c));
common=gcd(common,resultant(f3,f6,c));
common=gcd(common,resultant(f4,f5,c));
common=gcd(common,resultant(f4,f6,c));
common=gcd(common,resultant(f5,f6,c));
ideal E=cubic*q^6;
if (deg(common)!=15 || reduce(common,E)!=0)
{{
  print("BAD_RESULTANT_GCD");
  exit(1);
}}
ideal J0=q,f3,f4,f5,f6;
ideal J=std(J0);
if (reduce(p^3,J)!=0)
{{
  print("BAD_PIVOT_COMPONENT");
  exit(1);
}}
ideal I00=d,f3,f4,f5,f6;
ideal I0=std(I00);
if (reduce(c,I0)!=0)
{{
  print("BAD_D_ZERO_BRANCH");
  exit(1);
}}
ideal I40=d-4,f3,f4,f5,f6;
ideal I4=std(I40);
if (reduce(c-6,I4)!=0)
{{
  print("BAD_D_FOUR_BRANCH");
  exit(1);
}}
ideal IL0=d-4*lam,f3,f4,f5,f6;
ideal IL=std(IL0);
if (reduce(c-6*lam^2,IL)!=0)
{{
  print("BAD_D_LAMBDA_BRANCH");
  exit(1);
}}
print("GENERIC_DECOMPOSITION");
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=decomposition_program,
        text=True,
        capture_output=True,
        check=True,
        timeout=180,
    )
    assert completed.stdout.strip() == "GENERIC_DECOMPOSITION", (
        completed.stdout
    )

    # The h=0, g!=0 specialization chart admits a much smaller exact
    # certificate than the four-variable saturation.  Substitute
    # c=3*d^2/8 first.  The four moment numerators contain the invertible
    # factors g,g,g^2,g^2, respectively; after exact division, the
    # residual ideal together with the Rabinowitsch equation is the unit
    # ideal over Q.
    restricted_residuals: list[sp.Expr] = []
    for index, polynomial in enumerate(eliminated_polynomials):
        restricted = sp.Poly(
            sp.expand(
                polynomial.subs(c, sp.Rational(3, 8) * d**2)
            ),
            d,
            lam,
            domain=sp.QQ,
        )
        restricted = sp.Poly(
            restricted.clear_denoms()[1].primitive()[1],
            d,
            lam,
            domain=sp.QQ,
        )
        divisor_power = 1 if index < 2 else 2
        residual = restricted.exquo(
            sp.Poly(expected_cubic**divisor_power, d, lam)
        )
        restricted_residuals.append(residual.primitive()[1].as_expr())
    restricted_pivot = cleared_polynomial(
        pivot.subs(c, sp.Rational(3, 8) * d**2),
        (d, lam),
    )
    h_zero_program = "ring r=0,(z,d,lam),lp;\noption(redSB);\n"
    for order, polynomial in enumerate(restricted_residuals, start=3):
        h_zero_program += (
            f"poly r{order}={singular_expression(polynomial)};\n"
        )
    h_zero_program += f"""
poly p={singular_expression(restricted_pivot)};
poly g={singular_expression(expected_cubic)};
ideal I=r3,r4,r5,r6,z*p*g-1;
ideal G=std(I);
if (size(G)!=1 || G[1]!=1)
{{
  print("BAD_H_ZERO_G_OPEN_CHART");
  exit(1);
}}
print("H_ZERO_G_OPEN_UNIT");
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=h_zero_program,
        text=True,
        capture_output=True,
        check=True,
        timeout=60,
    )
    assert completed.stdout.strip() == "H_ZERO_G_OPEN_UNIT", (
        completed.stdout
    )

    # The remaining h!=0, p!=0 chart can be written without the
    # Rabinowitsch variable.  If q=lambda^4*(lambda-1)^4, then
    #
    #   q in (f3,f4,f5,f6,z*p*h-1)
    #
    # follows from q*(p*h)^N in (f3,f4,f5,f6).  Modular reconnaissance at
    # three good primes finds the same least exponent N=5 and the same
    # degree-order basis size 87.  At the first prime, an exact finite-field
    # lift has the displayed four multiplier profiles.  These computations
    # narrow the rational target but do not certify its Q-membership.
    squarefree_degeneracy = lam**4 * (lam - 1) ** 4
    saturation_primes = (101, 103, 107)
    saturation_exponents: dict[str, int] = {}
    saturation_basis_sizes: dict[str, int] = {}
    saturation_lift_profile: list[list[int]] = []
    for prime in saturation_primes:
        saturation_program = f"""
ring r={prime},(c,d,lam),dp;
option(redSB);
"""
        for order, polynomial in enumerate(eliminated_polynomials, start=3):
            saturation_program += (
                f"poly f{order}={singular_expression(polynomial)};\n"
            )
        saturation_program += f"""
poly p={singular_expression(pivot)};
poly h={singular_expression(expected_relation)};
poly q={singular_expression(squarefree_degeneracy)};
poly multiplier=p*h;
ideal F=f3,f4,f5,f6;
ideal G=std(F);
int exponent=-1;
poly powered=q;
int candidate;
for (candidate=0;candidate<=6;candidate++)
{{
  if (reduce(powered,G)==0)
  {{
    exponent=candidate;
    break;
  }}
  powered=powered*multiplier;
}}
print("SATURATION_EXPONENT");
print(exponent);
print("BASIS_SIZE");
print(size(G));
"""
        if prime == saturation_primes[0]:
            saturation_program += """
ideal T=powered;
matrix L=lift(F,T);
poly reconstructed=0;
int row;
for (row=1;row<=size(F);row++)
{
  reconstructed=reconstructed+F[row]*L[row,1];
}
if (reconstructed-powered!=0)
{
  print("BAD_FINITE_FIELD_SATURATION_LIFT");
  exit(1);
}
print("LIFT_PROFILE");
for (row=1;row<=size(F);row++)
{
  print(deg(L[row,1]));
  print(size(L[row,1]));
}
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=saturation_program,
            text=True,
            capture_output=True,
            check=True,
            timeout=300,
        )
        output_lines = completed.stdout.strip().splitlines()
        assert output_lines[:4] == [
            "SATURATION_EXPONENT",
            "5",
            "BASIS_SIZE",
            "87",
        ], (prime, completed.stdout)
        saturation_exponents[str(prime)] = int(output_lines[1])
        saturation_basis_sizes[str(prime)] = int(output_lines[3])
        if prime == saturation_primes[0]:
            assert output_lines[4:] == [
                "LIFT_PROFILE",
                "34",
                "5356",
                "29",
                "3679",
                "27",
                "3037",
                "22",
                "1853",
            ], completed.stdout
            saturation_lift_profile = [
                [int(output_lines[index]), int(output_lines[index + 1])]
                for index in range(5, len(output_lines), 2)
            ]

    # On each expected d-branch, shift c by its expected value and remove
    # the common powers t^(1,1,2,2).  The gcd of all pairwise t-resultants
    # gives a complete necessary specialization polynomial for an
    # additional c-value on that branch.
    branch_quartics = (
        22 * lam**4 - 54 * lam**3 + lam**2 - 54 * lam + 22,
        22 * lam**4 - 34 * lam**3 - 29 * lam**2 + 126 * lam - 63,
        63 * lam**4 - 126 * lam**3 + 29 * lam**2 + 34 * lam - 22,
    )
    branch_specs = (
        ("d_zero", 0, 0, pivot_orbit_factors[0], branch_quartics[0]),
        ("d_four", 4, 6, pivot_orbit_factors[1], branch_quartics[1]),
        (
            "d_four_lambda",
            4 * lam,
            6 * lam**2,
            pivot_orbit_factors[2],
            branch_quartics[2],
        ),
    )
    branch_exception_data: dict[str, str] = {}
    t = sp.symbols("t")
    for label, d_value, c_value, pivot_factor, quartic_factor in branch_specs:
        shifted_polynomials: list[sp.Expr] = []
        removed_powers: list[int] = []
        for polynomial in eliminated_polynomials:
            shifted = sp.Poly(
                sp.expand(
                    polynomial.subs(
                        {
                            d: d_value,
                            c: c_value + t,
                        }
                    )
                ),
                t,
                domain=sp.QQ[lam],
            )
            removed_power = min(
                monomial[0] for monomial, _ in shifted.terms()
            )
            removed_powers.append(removed_power)
            shifted_polynomials.append(
                sp.cancel(shifted.as_expr() / t**removed_power)
            )
        assert removed_powers == [1, 1, 2, 2]
        expected_branch_exception = (
            lam**4
            * (lam - 1) ** 4
            * pivot_factor**2
            * quartic_factor**6
        )
        branch_program = "ring r=0,(t,lam),lp;\n"
        for order, polynomial in enumerate(
            shifted_polynomials,
            start=3,
        ):
            branch_program += (
                f"poly f{order}={singular_expression(polynomial)};\n"
            )
        pair_index = 0
        for left in range(3, 7):
            for right in range(left + 1, 7):
                branch_program += (
                    f"poly r{pair_index}="
                    f"resultant(f{left},f{right},t);\n"
                )
                pair_index += 1
        branch_program += "poly exceptional=r0;\n"
        for index in range(1, pair_index):
            branch_program += (
                f"exceptional=gcd(exceptional,r{index});\n"
            )
        branch_program += f"""
ideal expected={singular_expression(expected_branch_exception)};
if (deg(exceptional)!=36 || reduce(exceptional,expected)!=0)
{{
  print("BAD_BRANCH_EXCEPTIONAL_GCD");
  exit(1);
}}
print("BRANCH_EXCEPTIONAL_GCD");
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=branch_program,
            text=True,
            capture_output=True,
            check=True,
            timeout=60,
        )
        assert completed.stdout.strip() == "BRANCH_EXCEPTIONAL_GCD", (
            label,
            completed.stdout,
        )
        branch_exception_data[label] = str(expected_branch_exception)

    # Two distinguished exceptional S_3-orbits.  In each quadratic
    # coefficient field the homogeneous projective fiber has degree four,
    # equals the four expected annihilator points set-theoretically, and
    # admits uniform eighth-power radical certificates.
    homogeneous_restricted = [
        sp.Poly(
            sp.fraction(sp.cancel(moment.subs(b, b_solution)))[0],
            a,
            c,
            d,
            e,
            domain=sp.QQ[lam],
        ).primitive()[1].as_expr()
        for moment in moments[1:]
    ]
    exceptional_fibers = {
        "pivot_annihilator_orbit": "lam^2+4*lam+1",
        "equianharmonic_orbit": "lam^2-lam+1",
        "branch_quartic_orbit": (
            "22*lam^4-54*lam^3+lam^2-54*lam+22"
        ),
    }
    for label, minimal_polynomial in exceptional_fibers.items():
        fiber_program = f"""
ring r=(0,lam),(a,c,d,e),dp;
minpoly={minimal_polynomial};
option(redSB);
ideal I={",".join(
    singular_expression(value) for value in homogeneous_restricted
)};
ideal GI=std(I);
ideal J1=c,d,e;
ideal J2=a,c,d;
ideal J3=a-e,c-6*e,d-4*e;
ideal J4=a-lam^4*e,c-6*lam^2*e,d-4*lam*e;
ideal J=intersect(intersect(J1,J2),intersect(J3,J4));
ideal GJ=std(J);
if (dim(GI)!=1 || mult(GI)!=4 || dim(GJ)!=1 || mult(GJ)!=4)
{{
  print("BAD_EXCEPTIONAL_FIBER");
  exit(1);
}}
if (size(reduce(I,GJ))!=0)
{{
  print("BAD_EXCEPTIONAL_CONTAINMENT");
  exit(1);
}}
int generator;
for (generator=1;generator<=size(GJ);generator++)
{{
  if (reduce(GJ[generator]^8,GI)!=0)
  {{
    print("BAD_EXCEPTIONAL_POWER");
    exit(1);
  }}
}}
print("EXCEPTIONAL_FIBER");
"""
        completed = subprocess.run(
            [singular, "-q"],
            input=fiber_program,
            text=True,
            capture_output=True,
            check=True,
            timeout=60,
        )
        assert completed.stdout.strip() == "EXCEPTIONAL_FIBER", (
            label,
            completed.stdout,
        )

    artifact = {
        "format": "two-variable-quartic-squarefree-pivot-v2",
        "field": "characteristic zero",
        "symbol": "u*v*(u-v)*(u-lambda*v)",
        "chart": "e=1 after eliminating b with mu_1",
        "mu_2_pivot": str(pivot),
        "pivot_boundary": {
            "status": "excluded over Q(lambda)",
            "certificate": "unit ideal from mu_2,...,mu_6",
        },
        "expected_affine_ideal": [
            "8*c-3*d^2",
            "d*(d-4)*(d-4*lambda)",
        ],
        "expected_points": [
            ["0", "0"],
            ["6", "4"],
            ["6*lambda^2", "4*lambda"],
        ],
        "cubic_factor_multiplicities_in_mu_3_through_mu_6": multiplicities,
        "generic_radical_decomposition": {
            "pairwise_c_resultant_gcd": (
                "d*(d-4)*(d-4*lambda)*q_lambda(d)^6"
            ),
            "extra_quadratic": str(extra_quadratic),
            "extra_component": "p^3=0, hence absent on p nonzero",
            "remaining_branches": [
                ["d=0", "c=0"],
                ["d=4", "c=6"],
                ["d=4*lambda", "c=6*lambda^2"],
            ],
            "status": (
                "generic squarefree fiber equals the annihilator sections"
            ),
        },
        "closed_special_squarefree_orbits": {
            label: {
                "minimal_polynomial": minimal_polynomial,
                "projective_fiber_degree": 4,
                "radical": "four annihilator sections",
                "radical_power_certificate": 8,
            }
            for label, minimal_polynomial in exceptional_fibers.items()
        },
        "harmonic_orbit": (
            "closed separately by the exact lambda=2 anchor certificate"
        ),
        "pivot_boundary_exceptional_gcd": str(
            pivot_exception_polynomial
        ),
        "expected_branch_exceptional_gcds": branch_exception_data,
        "h_zero_g_nonzero_chart": {
            "status": "exact unit ideal over Q",
            "reduction": (
                "substitute c=3*d^2/8 and divide the four moment "
                "numerators by the invertible powers g,g,g^2,g^2"
            ),
        },
        "h_nonzero_p_nonzero_modular_reconnaissance": {
            "target": (
                "lambda^4*(lambda-1)^4"
                "*(p*(8*c-3*d^2))^5 in (f3,f4,f5,f6)"
            ),
            "primes": list(saturation_primes),
            "least_saturation_exponents": saturation_exponents,
            "degree_order_basis_sizes": saturation_basis_sizes,
            "prime_101_lift_profile_degree_and_terms": (
                saturation_lift_profile
            ),
            "status": (
                "finite-field evidence only; exact Q lift remains open"
            ),
        },
        "remaining_exact_gate": (
            "lift lambda^4*(lambda-1)^4"
            "*(p*(8*c-3*d^2))^5 into (f3,f4,f5,f6) over Q"
        ),
        "written_source": (
            "extended-geometry/TWO_PAIR_SIC_BIDEGREE44_RANK_FRONTIER.md"
        ),
    }
    OUTPUT.write_text(json.dumps(artifact, indent=2) + "\n")
    print("PASS squarefree quartic: the mu_2 pivot boundary is empty")
    print("PASS squarefree quartic: expected affine ideal has three sections")
    print("PASS squarefree quartic: cubic factor multiplicities are 1,1,2,2")
    print(
        "PASS squarefree quartic: generic fiber is exactly the "
        "annihilator sections"
    )
    print(
        "PASS squarefree quartic: h=0 and g nonzero is an exact empty "
        "chart"
    )
    print(
        "PASS squarefree quartic: pivot-annihilator and equianharmonic "
        "and branch-quartic fibers have only the four annihilator points"
    )
    print(
        "PASS squarefree quartic: good-prime saturation exponent is "
        "consistently five (modular evidence only)"
    )
    print(f"PASS wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
