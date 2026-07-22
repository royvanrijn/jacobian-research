#!/usr/bin/env python3
"""Search the three-real-Gaussian cubic-Hermite correction exactly.

For

    S = a He_1(T) + b He_2(T) + c He_3(T),
    (a,b,c) = z (v_1,v_2,v_3) / h,

put I_k=E(T^k exp(S)).  Gaussian integration by parts gives

    3 c I_(k+2) - (1-2b) I_(k+1)
      + (a-3c) I_k + k I_(k-1) = 0.

Thus the connection has rank two.  Prescribing I_0=R=(h-z*h')/h makes
the logarithmic-derivative equation linear in q=I_1/I_0.  Eliminating q
from that equation and the remaining Riccati equation gives one rational
identity.  This script clears its known powers of h before expansion and
asks Singular whether its coefficient ideal is empty on a chosen chart of
the genuine-cubic locus.

The chart index r means

    v_3[0] = ... = v_3[r-1] = 0,  v_3[r] != 0.

The inequation is implemented with a Rabinowitsch variable s.  Clearing the
connection denominators creates spurious components unless this chart
saturation is included.

By default all coefficients of the compatibility polynomial are used.  The
``--low-coefficients`` option retains only that many lowest-z coefficients;
this is a search prefilter, not an all-orders certificate unless the full
coefficient ideal is subsequently checked.

The current seed is the first nonlinear canonical bridge seed
h=1+z^2-z^3.  Computations are modular reconnaissance unless ``--prime 0``
is requested.  A unit ideal modulo a prime is strong obstruction evidence,
but the characteristic-zero certificate should still be reconstructed.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys


def polynomial(prefix: str, degree: int, start: int = 0) -> str:
    terms = []
    for exponent in range(start, degree + 1):
        coefficient = f"{prefix}{exponent}"
        terms.append(coefficient if exponent == 0 else f"{coefficient}*z^{exponent}")
    return "+".join(terms) or "0"


def singular_program(
    degree: int,
    chart: int,
    prime: int,
    low_coefficients: int | None,
    algorithm: str,
) -> str:
    if not 0 <= chart <= degree:
        raise ValueError("chart must lie between zero and degree")

    parameters = [
        *(f"a{i}" for i in range(degree + 1)),
        *(f"b{i}" for i in range(degree + 1)),
        *(f"c{i}" for i in range(chart, degree + 1)),
    ]
    chart_constraints = [f"s*c{chart}-1"]
    a1 = "a1" if degree >= 1 else "0"
    b1 = "b1" if degree >= 1 else "0"
    c0 = "c0" if chart == 0 else "0"
    c1 = "c1" if degree >= 1 and chart <= 1 else "0"
    # The z^2 and z^3 coefficients of Z=D/h select the Gaussian formal
    # solution of the singular connection.  Besides being mathematically
    # necessary initial data, these low-degree equations are very effective
    # Groebner preconditioners on the difficult c0 != 0 chart.
    branch_constraints = [
        f"a0^2+2*b0^2+6*({c0})^2+4",
        (
            f"3*a0*({a1})+6*b0*({b1})+18*({c0})*({c1})"
            f"+3*a0^2*b0+18*a0*b0*({c0})+4*b0^3"
            f"+54*b0*({c0})^2-9"
        ),
    ]
    coefficient_selection = (
        "int start=1;"
        if low_coefficients is None
        else f"int start=nc-{low_coefficients - 1}; if (start<1) {{ start=1; }}"
    )

    # The capitalized polynomial names below are numerators after clearing
    # the displayed common powers of h in the module docstring derivation.
    return f"""ring r={prime},(z,{','.join(parameters)},s),dp;
option(redSB);
poly h=1+z^2-z^3;
poly hp=diff(h,z);
poly D=h-z*hp;
poly v1={polynomial('a', degree)};
poly v2={polynomial('b', degree)};
poly v3={polynomial('c', degree, chart)};

poly A0=z*v1;
poly B0=z*v2;
poly C0=z*v3;
poly Ap=h-2*B0;
poly Cp=3*C0;
poly Ep=A0-3*C0;
poly Bp=Ap^2-Cp*Ep;
poly Kp=Ap*Ep+Cp*h;

poly ad=diff(A0,z)*h-A0*hp;
poly bd=diff(B0,z)*h-B0*hp;
poly cd=diff(C0,z)*h-C0*hp;
poly Lp=diff(D,z)*h-D*hp;

poly Mp=ad*Cp^2+bd*Ap*Cp+cd*(Bp-3*Cp^2);
poly Np=bd*(Ep*Cp+Cp^2)+cd*Kp;
poly QN=D*Np+Lp*Cp^2*h;
poly QM=D*Mp;

poly Fp=Ap*Bp-Ep*Cp*Ap-2*Cp^2*h;
poly Hp=Fp-3*Ap*Cp^2;
poly Gp=-Ap*Kp+Ep^2*Cp;
poly Jp=Gp+3*Ep*Cp^2;
poly Tb=ad*Ap*Cp^2+bd*Cp*(Bp-Cp^2)+cd*Hp;
poly TQN=D*Tb-Lp*Cp^3*h;
poly T0p=-ad*Ep*Cp^2-bd*Cp*Kp+cd*Jp;

poly Compat=D*h^2*Cp^3*(diff(QN,z)*QM-QN*diff(QM,z))
  -QM*(TQN*QN+D*T0p*QM);
matrix MC=coef(Compat,z);
int nc=ncols(MC);
{coefficient_selection}
ideal I;
poly f;
for (int j=start;j<=nc;j++) {{
  f=MC[2,j];
  // c_chart is a unit on this chart.  Remove its pure common powers before
  // the Groebner calculation; this is equivalent to the chart saturation
  // and dramatically lowers the degrees of the coefficient generators.
  while ((f!=0) && (subst(f,c{chart},0)==0)) {{ f=f/c{chart}; }}
  I[j-start+1]=f;
}}
I=I,{','.join(branch_constraints + chart_constraints)};

print("compatibility_terms");
print(size(Compat));
print("available_coefficient_equations");
print(nc);
print("used_coefficient_equations");
print(nc-start+1);
print("starting_groebner_basis");
ideal G={algorithm}(I);
print("basis_size");
print(size(G));
print("normal_form_one");
print(reduce(1,G));
quit;
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--degree", type=int, default=1)
    parser.add_argument("--chart", type=int, default=0)
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--low-coefficients", type=int)
    parser.add_argument("--algorithm", choices=("std", "slimgb"), default="slimgb")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--emit", action="store_true")
    args = parser.parse_args()

    program = singular_program(
        args.degree,
        args.chart,
        args.prime,
        args.low_coefficients,
        args.algorithm,
    )
    if args.emit:
        print(program)
        return

    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required; use --emit to print the program")

    try:
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            timeout=args.timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as error:
        partial = error.stdout or ""
        if isinstance(partial, bytes):
            partial = partial.decode(errors="replace")
        if partial:
            print(partial, end="")
        print(f"TIMEOUT after {args.timeout} seconds", file=sys.stderr)
        raise SystemExit(124)

    print(completed.stdout, end="")
    if completed.stderr:
        print(completed.stderr, file=sys.stderr, end="")
    raise SystemExit(completed.returncode)


if __name__ == "__main__":
    main()
