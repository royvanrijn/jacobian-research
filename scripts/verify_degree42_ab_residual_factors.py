#!/usr/bin/env python3
"""Certify symbolic common-line residual factors on the degree-42 A/B strata."""

from __future__ import annotations

import argparse
import shutil
import subprocess

import sympy as sp

from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


def parse_singular_rational(value: str, symbols: dict[str, sp.Symbol]) -> sp.Expr:
    value = value.strip().replace("^", "**")
    return sp.factor(sp.sympify(value, locals=symbols))


def symbolic_residual(case: str) -> sp.Expr:
    normals, bases, residuals, _defect = transformed_problem()
    e1, e2, translation, w0, w1, w2 = bases
    factor_b = (
        e1**4 * e2**4
        - 2 * e1**2 * e2**5 * translation
        + e2**6 * translation**2
        - 8 * e1**5 * e2**2
        + 20 * e1**3 * e2**3 * translation
        - 12 * e1 * e2**4 * translation**2
        + 16 * e1**6
        - 48 * e1**4 * e2 * translation
        + 27 * e1**2 * e2**2 * translation**2
        + 9 * e2**3 * translation**3
        + 36 * e1**3 * translation**2
        - 54 * e1 * e2 * translation**3
        + 27 * translation**4
    )
    alpha_b = (
        8 * e1**5
        - 2 * e1**4 * e2**2
        - 20 * e1**3 * e2 * translation
        + 4 * e1**2 * e2**3 * translation
        + 6 * e1**2 * translation**2
        + 12 * e1 * e2**2 * translation**2
        - 2 * e2**4 * translation**2
        - 9 * e2 * translation**3
    )
    beta_b = (
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
    if case == "A":
        t_on_a = sp.cancel(
            e1**2 * (e2**2 - 4 * e1) / (e2 * (e2**2 - 6 * e1))
        )
        substitutions = {translation: t_on_a, w0: 0, w1: 0}
        coefficients = (e1, e2, w2)
        du = 2 * e1 - e2**2
        dv = -e2
        linear_u_coefficient = e2
        minpoly = ""
        relation = None
    elif case == "B":
        substitutions = {e2: 1, w0: 0, w1: 0}
        coefficients = (e1, w2, translation)
        du = beta_b.subs(e2, 1)
        dv = -alpha_b.subs(e2, 1)
        linear_u_coefficient = alpha_b.subs(e2, 1)
        minpoly = ""
        relation = factor_b.subs(e2, 1)
    elif case == "BN":
        r, q = sp.symbols("sync42ab_r sync42ab_q")
        x = (1 - r + q * (r - 1) * (2 * r - 1)) / 2
        t_normalized = -(1 + q) * (r - 1) ** 2 * (2 * r - 1) / 2
        substitutions = {
            e1: x,
            e2: 1,
            translation: t_normalized,
            w0: 0,
            w1: 0,
        }
        coefficients = (q, r, w2)
        du = r
        dv = 1
        linear_u_coefficient = 1
        minpoly = ""
        relation = q**2 + 3
        assert sp.rem(
            sp.together(factor_b.subs(substitutions)).as_numer_denom()[0],
            relation,
            q,
        ) == 0
        assert sp.rem(
            sp.together((alpha_b * r + beta_b).subs(substitutions))
            .as_numer_denom()[0],
            relation,
            q,
        ) == 0
    else:
        raise ValueError(case)
    specialized = [
        item.subs(substitutions, simultaneous=True)
        for item in residuals
    ]
    epsilon = sp.Symbol("sync42ab_eps")
    variables = normals + (epsilon,)
    zero_normals = ",".join(f"{variable},0" for variable in normals)
    marker = f"DEGREE42_{case}_SYMBOLIC_RHO"

    program = f"""
ring q=(0,{",".join(map(str, coefficients))}),({",".join(map(str, variables))}),dp;
{minpoly}
ideal I={",".join(serialize(item) for item in specialized)};

proc ev(poly p)
{{
  return(subst(p,{zero_normals}));
}}
proc D(poly p)
{{
  return(({serialize(du)})*diff(p,{normals[0]})
    +({serialize(dv)})*diff(p,{normals[1]}));
}}
proc DE(poly p)
{{
  return((({serialize(du)})+{epsilon})*diff(p,{normals[0]})
    +({serialize(dv)})*diff(p,{normals[1]}));
}}

poly a34=ev(diff(I[5],{normals[3]}));
poly a35=ev(diff(I[5],{normals[4]}));
poly a45=ev(diff(I[11],{normals[4]}));

poly eh3=ev(DE(DE(I[5])))/2;
poly eh4=ev(DE(DE(I[11])))/2;
poly eh5=ev(DE(DE(I[17])))/2;
poly ez25=eh5;
poly ez24=eh4+a45*ez25;
poly ez23=eh3+a34*ez24+a35*ez25;
proc ecubic(poly p)
{{
  return(ev(DE(DE(DE(p))))/6
    +ez23*ev(DE(diff(p,{normals[2]})))
    +ez24*ev(DE(diff(p,{normals[3]})))
    +ez25*ev(DE(diff(p,{normals[4]}))));
}}
poly ep3=ecubic(I[18]);
poly eq3=ecubic(I[19]);
poly ca2=subst(diff(ep3,{epsilon}),{epsilon},0)/({serialize(linear_u_coefficient)});
poly cb2=subst(diff(eq3,{epsilon}),{epsilon},0)/({serialize(linear_u_coefficient)});

poly h3=ev(D(D(I[5])))/2;
poly h4=ev(D(D(I[11])))/2;
poly h5=ev(D(D(I[17])))/2;
poly z25=h5;
poly z24=h4+a45*z25;
poly z23=h3+a34*z24+a35*z25;
proc cubic(poly p)
{{
  return(ev(D(D(D(p))))/6
    +z23*ev(D(diff(p,{normals[2]})))
    +z24*ev(D(diff(p,{normals[3]})))
    +z25*ev(D(diff(p,{normals[4]}))));
}}
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
poly rho=ca2*q4-cb2*p4;
print("{marker}");
print(rho);
"""
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    process = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=900,
        check=True,
    )
    assert marker in process.stdout, process.stdout + process.stderr
    value = "".join(process.stdout.split(marker, 1)[1].split())
    result = parse_singular_rational(
        value,
        {str(item): item for item in coefficients},
    )
    if case == "B":
        numerator, denominator = sp.fraction(sp.cancel(result))
        assert relation is not None
        numerator = sp.rem(numerator, relation, translation)
        result = sp.factor(sp.cancel(numerator / denominator))
    elif case == "BN":
        assert relation is not None
        numerator, denominator = sp.fraction(sp.cancel(result))
        numerator = sp.rem(numerator, relation, coefficients[0])
        denominator = sp.rem(denominator, relation, coefficients[0])
        q = coefficients[0]
        conjugate = denominator.subs(q, -q)
        numerator = sp.rem(sp.expand(numerator * conjugate), relation, q)
        denominator = sp.rem(sp.expand(denominator * conjugate), relation, q)
        result = sp.factor(sp.cancel(numerator / denominator))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        choices=("A", "B", "BN", "both"),
        default="both",
    )
    args = parser.parse_args()
    if args.case in ("A", "both"):
        rho_a = symbolic_residual("A")
        _, bases, _, _ = transformed_problem()
        e1, e2, _translation, _w0, _w1, w2 = bases
        p_a = (
            5376 * e1**8
            - 7168 * e1**7 * e2**2
            + 3808 * e1**6 * e2**4
            - 1008 * e1**5 * e2**6
            + 133 * e1**4 * e2**8
            - 7 * e1**3 * e2**10
            - 360 * e1**3 * e2**4 * w2
            + 192 * e1**2 * e2**6 * w2
            - 34 * e1 * e2**8 * w2
            + 2 * e2**10 * w2
        )
        expected_a = (
            -75
            * e1**2
            * w2
            * (4 * e1 - e2**2)
            * p_a
            / (512 * e2 * (6 * e1 - e2**2) ** 3)
        )
        assert sp.cancel(rho_a - expected_a) == 0
        print("A_RESIDUAL_FACTORIZATION")
        print(rho_a)
    if args.case == "B":
        rho_b = symbolic_residual("B")
        print("B_RESIDUAL_FACTORIZATION_MOD_B")
        print(rho_b)
    if args.case in ("BN", "both"):
        rho_bn = symbolic_residual("BN")
        q, r = sp.symbols("sync42ab_q sync42ab_r")
        _, bases, _, _ = transformed_problem()
        w2 = bases[-1]
        p_b = (
            112 * r**7
            - 560 * r**6
            + 1176 * r**5
            - 1344 * r**4
            + 903 * r**3
            - 357 * r**2
            + 77 * r
            - 7
            + (3 - 5 * r) * w2
        )
        expected_bn = (
            sp.Rational(75, 1024)
            * w2
            * (q - 1)
            * (r - 1) ** 2
            * (2 * r - 1)
            * p_b
        )
        assert sp.cancel(rho_bn - expected_bn) == 0
        print("B_NORMALIZED_RESIDUAL_FACTORIZATION")
        print(rho_bn)


if __name__ == "__main__":
    main()
