#!/usr/bin/env python3
"""Test the first residual quartic on the degree-42 A/B divisors.

The terminal binary cubics share one reduced linear factor on either
resultant divisor.  After blowing up the two-dimensional normal slice, the
common-line residual theorem identifies the first remaining coefficient as

    rho = a_2(p) q_4(p) - b_2(p) p_4(p),

where p_3=L*a_2 and q_3=L*b_2.  A nonzero value at one point of an
irreducible divisor proves that its quartic restriction is not identically
zero.  The A witness is characteristic zero; the B witness is at the good
prime 103, which is enough to exclude a characteristic-zero polynomial
identity.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp

from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


def residual_value(
    characteristic: int,
    base_values: tuple[sp.Expr | int, ...],
    direction: tuple[int, int],
    linear_u_coefficient: int,
) -> str:
    normals, bases, residuals, _defect = transformed_problem()
    specialized = [
        item.subs(dict(zip(bases, base_values)), simultaneous=True)
        for item in residuals
    ]
    u, v = sp.symbols("sync42ab_U sync42ab_V")
    variables = normals + (u, v)
    zero_normals = ",".join(f"{variable},0" for variable in normals)
    du, dv = direction
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
poly p3d=subst(p3,{u},{du},{v},{dv});
poly q3d=subst(q3,{u},{du},{v},{dv});
poly p4d=subst(p4,{u},{du},{v},{dv});
poly q4d=subst(q4,{u},{du},{v},{dv});
poly a2=subst(diff(p3,{u}),{u},{du},{v},{dv})/{linear_u_coefficient};
poly b2=subst(diff(q3,{u}),{u},{du},{v},{dv})/{linear_u_coefficient};
poly rho=a2*q4d-b2*p4d;
print("DEGREE42_AB_QUARTIC");
print(p3d);
print(q3d);
print(p4d);
print(q4d);
print(a2);
print(b2);
print(rho);
"""
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for the certificate"
    process = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    compact = process.stdout.split()
    marker = compact.index("DEGREE42_AB_QUARTIC")
    values = compact[marker + 1 : marker + 8]
    assert len(values) == 7, process.stdout + process.stderr
    assert values[0:2] == ["0", "0"], process.stdout + process.stderr
    assert values[-1] != "0", process.stdout + process.stderr
    return values[-1]


def main() -> None:
    # A=0, B=2187/625 != 0, and e2(6e1-e2^2)=5.
    e1, e2, t = sp.symbols("e1 e2 t")
    factor_a = 4 * e1**3 - e1**2 * e2**2 + e2**3 * t - 6 * e1 * e2 * t
    factor_b = (
        e1**4 * e2**4
        - 2 * e1**2 * e2**5 * t
        + e2**6 * t**2
        - 8 * e1**5 * e2**2
        + 20 * e1**3 * e2**3 * t
        - 12 * e1 * e2**4 * t**2
        + 16 * e1**6
        - 48 * e1**4 * e2 * t
        + 27 * e1**2 * e2**2 * t**2
        + 9 * e2**3 * t**3
        + 36 * e1**3 * t**2
        - 54 * e1 * e2 * t**3
        + 27 * t**4
    )
    a_point = {e1: 1, e2: 1, t: sp.Rational(3, 5)}
    assert factor_a.subs(a_point) == 0
    assert factor_b.subs(a_point) == sp.Rational(2187, 625)
    rho_a = residual_value(
        0,
        (1, 1, sp.Rational(3, 5), 0, 0, 1),
        (1, -1),
        1,
    )

    # Modulo 103: B(1,1,21)=0, A=1, alpha_B=9, beta_B=3.
    b_point = {e1: 1, e2: 1, t: 21}
    assert int(factor_b.subs(b_point)) % 103 == 0
    assert int(factor_a.subs(b_point)) % 103 == 1
    rho_b = residual_value(
        103,
        (1, 1, 21, 0, 0, 1),
        (3, -9),
        9,
    )
    print(
        "PASS: the common-line residual quartic is generically nonzero "
        f"on A=0 (rho={rho_a}) and B=0 (rho={rho_b} mod 103)"
    )


if __name__ == "__main__":
    main()
