#!/usr/bin/env python3
"""Certify nonzero fifth common-line residuals on the degree-42 A/B graphs.

On either residual graph the terminal cubic equations share one reduced
linear factor and their invariant quartic residual vanishes.  This checker
solves the three pivot equations through fourth order, follows one terminal
equation through its quartic correction, and evaluates the invariant fifth
coefficient of the other terminal equation.

The A witness is exact over QQ.  The B witness is over the good prime 103;
a nonzero value there excludes an identically zero characteristic-zero
residual on the corresponding geometric normalization.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp

from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


def quintic_residual_value(
    characteristic: int,
    base_values: tuple[sp.Expr | int, ...],
    direction: tuple[int, int],
) -> tuple[str, str, str, str]:
    """Return ``(c1, p5, q5, rho5)`` at one residual-graph point."""

    normals, bases, residuals, _defect = transformed_problem()
    specialized = [
        item.subs(dict(zip(bases, base_values)), simultaneous=True)
        for item in residuals
    ]
    u, v, radial = sp.symbols(
        "sync42ab5_U sync42ab5_V sync42ab5_R"
    )
    variables = normals + (u, v, radial)
    zero_normals = ",".join(f"{variable},0" for variable in normals)
    du, dv = direction
    marker = "DEGREE42_AB_QUINTIC"

    program = f"""
ring q={characteristic},({",".join(map(str, variables))}),dp;
ideal I={",".join(serialize(item) for item in specialized)};

proc ev(poly p)
{{
  return(subst(p,{zero_normals}));
}}
proc D(poly p)
{{
  return({u}*diff(p,{normals[0]})+{v}*diff(p,{normals[1]}));
}}
proc path(poly p, poly X1, poly X2, poly X3, poly X4, poly X5)
{{
  p=subst(p,{normals[0]},X1);
  p=subst(p,{normals[1]},X2);
  p=subst(p,{normals[2]},X3);
  p=subst(p,{normals[3]},X4);
  p=subst(p,{normals[4]},X5);
  return(p);
}}
proc hc(poly p, int k)
{{
  int j;
  int fac=1;
  for (j=1; j<=k; j++)
  {{
    p=diff(p,{radial});
    fac=fac*j;
  }}
  return(subst(p,{radial},0)/fac);
}}

// The constant pivot Jacobian is upper-unitriangular.
poly a34=ev(diff(I[5],{normals[3]}));
poly a35=ev(diff(I[5],{normals[4]}));
poly a45=ev(diff(I[11],{normals[4]}));

// Cubic terminal forms after solving the quadratic pivot coefficients.
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
poly p3=cubic(I[18]);
poly q3=cubic(I[19]);
poly p3u=subst(diff(p3,{u}),{u},{du},{v},{dv});
poly q3u=subst(diff(q3,{u}),{u},{du},{v},{dv});

// Solve the pivot equations along the straight common-line arc.
poly X1={du}*{radial};
poly X2={dv}*{radial};
poly h23=hc(path(I[5],X1,X2,0,0,0),2);
poly h24=hc(path(I[11],X1,X2,0,0,0),2);
poly h25=hc(path(I[17],X1,X2,0,0,0),2);
poly zz25=h25;
poly zz24=h24+a45*zz25;
poly zz23=h23+a34*zz24+a35*zz25;
poly X3=zz23*{radial}^2;
poly X4=zz24*{radial}^2;
poly X5=zz25*{radial}^2;

poly h33=hc(path(I[5],X1,X2,X3,X4,X5),3);
poly h34=hc(path(I[11],X1,X2,X3,X4,X5),3);
poly h35=hc(path(I[17],X1,X2,X3,X4,X5),3);
poly zz35=h35;
poly zz34=h34+a45*zz35;
poly zz33=h33+a34*zz34+a35*zz35;
poly Y3=X3+zz33*{radial}^3;
poly Y4=X4+zz34*{radial}^3;
poly Y5=X5+zz35*{radial}^3;
poly p4=hc(path(I[18],X1,X2,Y3,Y4,Y5),4);
poly q4=hc(path(I[19],X1,X2,Y3,Y4,Y5),4);

// Follow p=0 through quartic order by bending transversely in X1.
poly c1=-p4/p3u;
X1={du}*{radial}+c1*{radial}^2;
X2={dv}*{radial};
X3=zz23*{radial}^2;
X4=zz24*{radial}^2;
X5=zz25*{radial}^2;

// Re-solve cubic pivot coefficients after the transverse bend.
h33=hc(path(I[5],X1,X2,X3,X4,X5),3);
h34=hc(path(I[11],X1,X2,X3,X4,X5),3);
h35=hc(path(I[17],X1,X2,X3,X4,X5),3);
zz35=h35;
zz34=h34+a45*zz35;
zz33=h33+a34*zz34+a35*zz35;
Y3=X3+zz33*{radial}^3;
Y4=X4+zz34*{radial}^3;
Y5=X5+zz35*{radial}^3;
poly p4zero=hc(path(I[18],X1,X2,Y3,Y4,Y5),4);
poly q4zero=hc(path(I[19],X1,X2,Y3,Y4,Y5),4);

// Solve the fourth pivot coefficients; fifth pivot terms cannot enter the
// terminal order because both terminal equations have zero differential.
poly h43=hc(path(I[5],X1,X2,Y3,Y4,Y5),4);
poly h44=hc(path(I[11],X1,X2,Y3,Y4,Y5),4);
poly h45=hc(path(I[17],X1,X2,Y3,Y4,Y5),4);
poly zz45=h45;
poly zz44=h44+a45*zz45;
poly zz43=h43+a34*zz44+a35*zz45;
poly Z3=Y3+zz43*{radial}^4;
poly Z4=Y4+zz44*{radial}^4;
poly Z5=Y5+zz45*{radial}^4;
poly p5=hc(path(I[18],X1,X2,Z3,Z4,Z5),5);
poly q5=hc(path(I[19],X1,X2,Z3,Z4,Z5),5);
poly rho5=p3u*q5-q3u*p5;

print("{marker}");
print(p3u);
print(q3u);
print(c1);
print(p4zero);
print(q4zero);
print(p5);
print(q5);
print(rho5);
"""
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for the certificate"
    process = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=600,
        check=True,
    )
    compact = process.stdout.split()
    marker_index = compact.index(marker)
    values = compact[marker_index + 1 : marker_index + 9]
    assert len(values) == 8, process.stdout + process.stderr
    p3u, _q3u, c1, p4zero, q4zero, p5, q5, rho5 = values
    assert p3u != "0", process.stdout + process.stderr
    assert p4zero == "0", process.stdout + process.stderr
    assert q4zero == "0", process.stdout + process.stderr
    assert rho5 != "0", process.stdout + process.stderr
    return c1, p5, q5, rho5


def main() -> None:
    e1, e2, translation = sp.symbols("e1 e2 translation")
    factor_a = (
        4 * e1**3
        - e1**2 * e2**2
        + e2**3 * translation
        - 6 * e1 * e2 * translation
    )
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

    # A=0 at (e1,e2,t)=(1,1,3/5), and P_A=0 fixes w2=567/100.
    a_point = {e1: 1, e2: 1, translation: sp.Rational(3, 5)}
    assert factor_a.subs(a_point) == 0
    assert factor_b.subs(a_point) == sp.Rational(2187, 625)
    a_result = quintic_residual_value(
        0,
        (1, 1, sp.Rational(3, 5), 0, 0, sp.Rational(567, 100)),
        (1, -1),
    )

    # On the B-normalization modulo 103 take q=10 (q^2=-3), r=0.
    # Then e1=t=(1+q)/2=57, w2=Q_B(0)/(-3)=7/3=71.
    prime = 103
    b_point = {e1: 57, e2: 1, translation: 57}
    assert int(factor_b.subs(b_point)) % prime == 0
    assert int(factor_a.subs(b_point)) % prime == 73
    assert int(alpha_b.subs(b_point)) % prime == 36
    assert int(beta_b.subs(b_point)) % prime == 0
    b_result = quintic_residual_value(
        prime,
        (57, 1, 57, 0, 0, sp.Rational(7, 3)),
        (0, 1),
    )

    print(
        "PASS: the invariant fifth common-line residual is generically "
        "nonzero on both degree-42 residual graphs; "
        f"A has (c1,p5,q5,rho5)={a_result}, "
        f"B has {b_result} modulo {prime}"
    )


if __name__ == "__main__":
    main()
