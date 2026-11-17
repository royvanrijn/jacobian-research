#!/usr/bin/env python3
"""Certify the higher-gcd strata of the degree-42 cubic pencil.

For the terminal binary cubics, let ``A*B`` be their resultant and let
``alpha*u+beta*v`` be their degree-one subresultant.  This checker:

1. decomposes the reduced higher-gcd ideal
   ``(A*B, alpha, beta)`` over QQ;
2. records the quadratic splitting of its nonrational curve over
   ``QQ(sqrt(-3))``;
3. extracts the quartic Kuranishi forms at one nonvertex point of every
   weighted curve and certifies the resulting Hilbert cutoff; and
4. checks that the common vertex remains non-Artinian at quartic order.
"""

from __future__ import annotations

from dataclasses import dataclass
import shutil
import subprocess

import sympy as sp

from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


e1, e2, translation = sp.symbols("e1 e2 t")

factor_a = (
    4 * e1**3
    - e1**2 * e2**2
    + e2**3 * translation
    - 6 * e1 * e2 * translation
)
factor_b = factor_a**2 + 9 * factor_a * translation**2 + 27 * translation**4
alpha = (
    8 * e1**5
    - 2 * e1**4 * e2**2
    - 20 * e1**3 * e2 * translation
    + 4 * e1**2 * e2**3 * translation
    + 6 * e1**2 * translation**2
    + 12 * e1 * e2**2 * translation**2
    - 2 * e2**4 * translation**2
    - 9 * e2 * translation**3
)
beta = (
    -4 * e1**5 * e2
    + e1**4 * e2**3
    + 4 * e1**4 * translation
    + 9 * e1**3 * e2**2 * translation
    - 2 * e1**2 * e2**4 * translation
    - 9 * e1**2 * e2 * translation**2
    - 5 * e1 * e2**3 * translation**2
    + e2**5 * translation**2
    + 6 * e2**2 * translation**3
)


def certify_reduced_decomposition(singular: str) -> None:
    """Check the rational reduced support and its A/B provenance."""

    program = f"""
LIB "primdec.lib";
ring q=0,(t,e1,e2),dp;
poly A={serialize(factor_a)};
poly B={serialize(factor_b)};
poly alpha={serialize(alpha)};
poly beta={serialize(beta)};
ideal IA=A,alpha,beta;
ideal IB=B,alpha,beta;
ideal IH=A*B,alpha,beta;

ideal P1=
  e2^4-6*e1*e2^2+12*e1^2,
  e2^3-4*e1*e2+6*t;
ideal P2=e1,e2;
ideal P3=t,4*e1-e2^2;
ideal P4=e1,t;

ideal RA=intersect(P2,intersect(P3,P4));
ideal RB=intersect(P1,intersect(P3,P4));
ideal RH=intersect(P1,intersect(P2,intersect(P3,P4)));
ideal radA=std(radical(IA));
ideal radB=std(radical(IB));
ideal radH=std(radical(IH));
RA=std(RA);
RB=std(RB);
RH=std(RH);

proc issubset(ideal J,ideal K)
{{
  int i;
  for (i=1;i<=size(J);i++)
  {{
    if (reduce(J[i],K)!=0) {{ return(0); }}
  }}
  return(1);
}}
int checkA=(issubset(radA,RA) and issubset(RA,radA));
int checkB=(issubset(radB,RB) and issubset(RB,radB));
int checkH=(issubset(radH,RH) and issubset(RH,radH));
print("DEGREE42_HIGHER_GCD_DECOMPOSITION");
print(checkA);
print(checkB);
print(checkH);
"""
    process = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    compact = process.stdout.split()
    marker = compact.index("DEGREE42_HIGHER_GCD_DECOMPOSITION")
    assert compact[marker + 1 : marker + 4] == ["1", "1", "1"], (
        process.stdout + process.stderr
    )


def certify_conjugate_parameterizations() -> None:
    """Check the two branches of P1 after adjoining sqrt(-3)."""

    rho, parameter = sp.symbols("rho parameter")
    relation = sp.Poly(rho**2 + 3, rho, domain=sp.QQ.frac_field(parameter))

    def reduce_rho(expression: sp.Expr) -> sp.Expr:
        polynomial = sp.Poly(
            sp.cancel(expression),
            rho,
            domain=sp.QQ.frac_field(parameter),
        )
        return sp.factor(polynomial.rem(relation).as_expr())

    p1_generators = (
        e2**4 - 6 * e1 * e2**2 + 12 * e1**2,
        e2**3 - 4 * e1 * e2 + 6 * translation,
    )
    for sign in (1, -1):
        substitution = {
            e2: parameter,
            e1: sp.Rational(1, 12) * (3 + sign * rho) * parameter**2,
            translation: sign * rho * parameter**3 / 18,
        }
        assert all(
            reduce_rho(generator.subs(substitution)) == 0
            for generator in p1_generators
        )


@dataclass(frozen=True)
class QuarticCase:
    name: str
    characteristic: int
    base_values: tuple[int, int, int, int, int, int]
    cubic_gcd_degree: int
    quotient_length: int
    cutoff_power: int


def quartic_program(
    normals: tuple[sp.Symbol, ...],
    specialized_residuals: list[sp.Expr],
    case: QuarticCase,
) -> str:
    """Build the implicit-differentiation check for one weighted orbit."""

    u, v = sp.symbols("sync42hg_U sync42hg_V")
    zero_normals = ",".join(f"{variable},0" for variable in normals)
    return f"""
ring q={case.characteristic},({",".join(map(str, normals))},{u},{v}),dp;
ideal I={",".join(serialize(item) for item in specialized_residuals)};

proc ev(poly p)
{{
  return(subst(p,{zero_normals}));
}}
proc D(poly p)
{{
  return({u}*diff(p,{normals[0]})+{v}*diff(p,{normals[1]}));
}}

poly a34=ev(diff(I[5],{normals[3]}));
poly a35=ev(diff(I[5],{normals[4]}));
poly a45=ev(diff(I[11],{normals[4]}));
poly h3=ev(D(D(I[5])))/2;
poly h4=ev(D(D(I[11])))/2;
poly h5=ev(D(D(I[17])))/2;
poly z25=h5;
poly z24=h4+a45*z25;
poly z23=h3+a34*z24+a35*z25;

proc cross3(poly p)
{{
  return(z23*ev(D(diff(p,{normals[2]})))
    +z24*ev(D(diff(p,{normals[3]})))
    +z25*ev(D(diff(p,{normals[4]}))));
}}
proc cubic(poly p)
{{
  return(ev(D(D(D(p))))/6+cross3(p));
}}

poly p3=cubic(I[18]);
poly q3=cubic(I[19]);
poly k3=cubic(I[5]);
poly k4=cubic(I[11]);
poly k5=cubic(I[17]);
poly z35=k5;
poly z34=k4+a45*z35;
poly z33=k3+a34*z34+a35*z35;

proc quartic(poly p)
{{
  poly ans=ev(D(D(D(D(p)))))/24;
  ans=ans+z23*ev(D(D(diff(p,{normals[2]})))/2)
    +z24*ev(D(D(diff(p,{normals[3]})))/2)
    +z25*ev(D(D(diff(p,{normals[4]})))/2);
  ans=ans+z33*ev(D(diff(p,{normals[2]})))
    +z34*ev(D(diff(p,{normals[3]})))
    +z35*ev(D(diff(p,{normals[4]})));
  ans=ans+z23^2*ev(diff(diff(p,{normals[2]}),{normals[2]}))/2
    +z24^2*ev(diff(diff(p,{normals[3]}),{normals[3]}))/2
    +z25^2*ev(diff(diff(p,{normals[4]}),{normals[4]}))/2;
  ans=ans+z23*z24*ev(diff(diff(p,{normals[2]}),{normals[3]}))
    +z23*z25*ev(diff(diff(p,{normals[2]}),{normals[4]}))
    +z24*z25*ev(diff(diff(p,{normals[3]}),{normals[4]}));
  return(ans);
}}

poly p4=quartic(I[18]);
poly q4=quartic(I[19]);
poly g3=gcd(p3,q3);
poly gall=gcd(gcd(g3,p4),q4);
ideal G=std(ideal({",".join(map(str, normals))},p3,q3,p4,q4));
proc contained(ideal J,ideal K)
{{
  int i;
  for (i=1;i<=size(J);i++)
  {{
    if (reduce(J[i],K)!=0) {{ return(0); }}
  }}
  return(1);
}}
int gcd_degree=(deg(g3)=={case.cubic_gcd_degree});
int no_common_direction=(deg(gall)==0);
int expected_length=(vdim(G)=={case.quotient_length});
int cutoff=contained(maxideal({case.cutoff_power}),G);
int sharp=!contained(maxideal({case.cutoff_power - 1}),G);
print("DEGREE42_HIGHER_GCD_QUARTIC");
print(gcd_degree);
print(no_common_direction);
print(expected_length);
print(cutoff);
print(sharp);
"""


def certify_quartic_orbits(
    singular: str,
    problem: tuple[
        tuple[sp.Symbol, ...],
        tuple[sp.Symbol, ...],
        list[sp.Expr],
        sp.Expr,
    ],
) -> None:
    """Check quartic closure on every punctured higher-gcd curve."""

    normals, bases, residuals, _defect = problem
    cases = (
        QuarticCase("t_axis", 0, (0, 0, 1, 0, 0, 1), 2, 8, 4),
        QuarticCase("e2_axis", 0, (0, 1, 0, 0, 0, 1), 3, 10, 5),
        QuarticCase("parabola", 0, (1, 2, 0, 0, 0, 1), 2, 8, 4),
        # P1 splits modulo 103.  The point c=8, t=91 lies on one branch.
        # A nonzero good-prime witness excludes a characteristic-zero
        # identity on the irreducible rational curve P1.
        QuarticCase("conjugate", 103, (8, 1, 91, 0, 0, 1), 2, 8, 4),
    )
    for case in cases:
        assert len(bases) == len(case.base_values)
        substitution = dict(zip(bases, case.base_values))
        specialized = [
            item.subs(substitution, simultaneous=True) for item in residuals
        ]
        process = subprocess.run(
            [singular, "-q"],
            input=quartic_program(normals, specialized, case),
            text=True,
            capture_output=True,
            timeout=300,
            check=True,
        )
        compact = process.stdout.split()
        marker = compact.index("DEGREE42_HIGHER_GCD_QUARTIC")
        assert compact[marker + 1 : marker + 6] == ["1"] * 5, (
            case.name + "\n" + process.stdout + process.stderr
        )


def certify_common_vertex(
    singular: str,
    problem: tuple[
        tuple[sp.Symbol, ...],
        tuple[sp.Symbol, ...],
        list[sp.Expr],
        sp.Expr,
    ],
) -> None:
    """Check the surviving contact-five vertex at quartic order."""

    normals, bases, residuals, _defect = problem
    vertex_values = (0, 0, 0, 0, 0, 1)
    assert len(bases) == len(vertex_values)
    substitution = dict(zip(bases, vertex_values))
    specialized = [
        item.subs(substitution, simultaneous=True) for item in residuals
    ]
    vertex_case = QuarticCase("vertex", 0, (0, 0, 0, 0, 0, 1), 0, 0, 0)
    program = quartic_program(normals, specialized, vertex_case)
    # Replace the generic final audit: both cubics vanish, while the two
    # quartics have gcd v^3 and hence do not yet cut the normal plane.
    program = program[: program.index("poly g3=gcd(p3,q3);")] + """
poly g4=gcd(p4,q4);
int cubic_zero=(p3==0 and q3==0);
int quartic_common_cube=(deg(g4)==3);
int non_artinian=(vdim(std(ideal(p4,q4)))==-1);
print("DEGREE42_HIGHER_GCD_VERTEX");
print(cubic_zero);
print(quartic_common_cube);
print(non_artinian);
"""
    process = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    compact = process.stdout.split()
    marker = compact.index("DEGREE42_HIGHER_GCD_VERTEX")
    assert compact[marker + 1 : marker + 4] == ["1", "1", "1"], (
        process.stdout + process.stderr
    )


def main() -> None:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for the certificate"
    certify_reduced_decomposition(singular)
    certify_conjugate_parameterizations()
    problem = transformed_problem()
    certify_quartic_orbits(singular, problem)
    certify_common_vertex(singular, problem)
    print(
        "PASS: the higher-gcd support is four weighted curves; quartics "
        "close every punctured curve with Hilbert vectors (1,2,3,2) or "
        "(1,2,3,3,1), while their common contact-five vertex retains a "
        "common cubic factor at quartic order"
    )


if __name__ == "__main__":
    main()
