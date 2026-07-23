#!/usr/bin/env python3
"""Factor the degree-42 common-line quartic residuals on A and B.

The cubic resultant divisors have one common linear factor away from the
higher-gcd locus.  This checker derives the generic cubic and quartic
Kuranishi forms from the original transported residual system, forms the
cleared common-line residual, and factors it in the A and geometric B
coordinate rings.
"""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp

from verify_degree42_higher_gcd_strata import QuarticCase, quartic_program
from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


def generic_quartic_prefix(
    algebraic: bool,
) -> tuple[str, tuple[sp.Symbol, ...], tuple[sp.Symbol, ...]]:
    """Return a Singular program through the generic quartic pair."""

    normals, bases, residuals, _defect = transformed_problem()
    e1, e2, translation, w0, w1, _w2 = bases
    specialized = [
        item.subs({w0: 0, w1: 0}, simultaneous=True) for item in residuals
    ]
    dummy = QuarticCase("generic", 0, (0, 0, 0, 0, 0, 0), 1, 1, 1)
    program = quartic_program(normals, specialized, dummy)
    old_ring = (
        f'ring q=0,({",".join(map(str, normals))},'
        "sync42hg_U,sync42hg_V),dp;"
    )
    coefficient_field = "(0,rho)" if algebraic else "0"
    new_ring = (
        f'ring q={coefficient_field},'
        f'({",".join(map(str, normals + bases))},'
        "sync42hg_U,sync42hg_V),(dp(5),lp(3),dp(3),dp(2));"
    )
    if algebraic:
        new_ring += "\nminpoly=rho^2+3;"
    assert old_ring in program
    program = program.replace(old_ring, new_ring, 1)
    program = program[: program.index("poly g3=gcd(p3,q3);")]
    # The block order makes e1^3 the leading term of A-k*t^2 while the
    # five eliminated normal variables remain in the first block.
    assert (e1, e2, translation) == bases[:3]
    return program, normals, bases


def run_singular(program: str, marker: str, count: int) -> list[str]:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for the certificate"
    process = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=900,
        check=True,
    )
    compact = process.stdout.split()
    location = compact.index(marker)
    values = compact[location + 1 : location + 1 + count]
    assert len(values) == count, process.stdout + process.stderr
    return values


def certify_rational_factorizations() -> None:
    """Check the A factorization and the rational norm on B."""

    program, _normals, bases = generic_quartic_prefix(algebraic=False)
    e1, e2, translation, _w0, _w1, w2 = bases
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
    residual_a = (
        28 * e1**2 * e2**3 * translation
        + 336 * e1**2 * translation**2
        - 7 * e1 * e2**5 * translation
        - 280 * e1 * e2**2 * translation**2
        - 10 * e1 * e2**2 * w2
        + 35 * e2**4 * translation**2
        + 2 * e2**4 * w2
    )
    program += f"""
poly A={serialize(factor_a)};
poly B={serialize(factor_b)};
poly alpha={serialize(alpha)};
poly beta={serialize(beta)};
poly RA={serialize(residual_a)};

poly Au={e2}^2-2*{e1};
poly Av={e2};
poly Adp=subst(diff(p3,sync42hg_U),
  sync42hg_U,Au,sync42hg_V,Av);
poly Adq=subst(diff(q3,sync42hg_U),
  sync42hg_U,Au,sync42hg_V,Av);
poly Ap4=subst(p4,sync42hg_U,Au,sync42hg_V,Av);
poly Aq4=subst(q4,sync42hg_U,Au,sync42hg_V,Av);
poly rhoA=Adp*Aq4-Adq*Ap4;
ideal GA=std(ideal(A));
poly expectedA=-75/512*{translation}*{e2}^3*{w2}*RA;
int factorA=(reduce(rhoA-expectedA,GA)==0);

list facRA=factorize(RA);
ideal factorsRA=facRA[1];
int irreducibleRA=(size(factorsRA)==2);
int linearRA=(
  diff(RA,{w2})!=0 and diff(diff(RA,{w2}),{w2})==0
);

poly witnessA=subst(A,{e1},1,{e2},1,{translation},3/5);
poly witnessRA=subst(RA,{e1},1,{e2},1,{translation},3/5,{w2},567/100);
poly witnessAlpha=subst(alpha,{e1},1,{e2},1,{translation},3/5);
poly witnessBeta=subst(beta,{e1},1,{e2},1,{translation},3/5);
int newA=(
  witnessA==0 and witnessRA==0 and
  witnessAlpha!=0 and witnessBeta!=0
);

poly Bu=beta;
poly Bv=-alpha;
poly Bdp=subst(diff(p3,sync42hg_U),
  sync42hg_U,Bu,sync42hg_V,Bv);
poly Bdq=subst(diff(q3,sync42hg_U),
  sync42hg_U,Bu,sync42hg_V,Bv);
poly Bp4=subst(p4,sync42hg_U,Bu,sync42hg_V,Bv);
poly Bq4=subst(q4,sync42hg_U,Bu,sync42hg_V,Bv);
poly rhoBQ=Bdp*Bq4-Bdq*Bp4;
poly remBQ=reduce(rhoBQ,std(ideal(B)));
poly normResidual=remBQ/({translation}^11*{w2});
int exactNorm=(
  remBQ-{translation}^11*{w2}*normResidual==0
);
list facNorm=factorize(normResidual);
ideal factorsNorm=facNorm[1];
int irreducibleNorm=(size(factorsNorm)==2);

print("DEGREE42_AB_RATIONAL_FACTORIZATION");
print(factorA);
print(irreducibleRA);
print(linearRA);
print(newA);
print(exactNorm);
print(irreducibleNorm);
"""
    values = run_singular(
        program,
        "DEGREE42_AB_RATIONAL_FACTORIZATION",
        6,
    )
    assert values == ["1"] * 6, values


def certify_geometric_b_factorization() -> None:
    """Factor the residual on one B branch over QQ(sqrt(-3))."""

    program, _normals, bases = generic_quartic_prefix(algebraic=True)
    e1, e2, translation, _w0, _w1, w2 = bases
    factor_a = (
        4 * e1**3
        - e1**2 * e2**2
        + e2**3 * translation
        - 6 * e1 * e2 * translation
    )
    program += f"""
poly A={serialize(factor_a)};
number k=(-9+3*rho)/2;
poly C=A-k*{translation}^2;
poly alphaReduced=
  (2*k+6)*{e1}^2-(2*k+9)*{e2}*{translation};
poly betaReduced=
  -(k+3)*{e1}^2*{e2}
  +k*{e1}*{translation}
  +(k+6)*{e2}^2*{translation};

poly Bdp=subst(diff(p3,sync42hg_U),
  sync42hg_U,betaReduced,sync42hg_V,-alphaReduced);
poly Bdq=subst(diff(q3,sync42hg_U),
  sync42hg_U,betaReduced,sync42hg_V,-alphaReduced);
poly Bp4=subst(p4,
  sync42hg_U,betaReduced,sync42hg_V,-alphaReduced);
poly Bq4=subst(q4,
  sync42hg_U,betaReduced,sync42hg_V,-alphaReduced);
ideal GC=std(ideal(C));
poly remB=reduce(Bdp*Bq4-Bdq*Bp4,GC);
poly RB=remB/({translation}*{w2});
int exactB=(remB-{translation}*{w2}*RB==0);
int linearB=(
  diff(RB,{w2})!=0 and diff(diff(RB,{w2}),{w2})==0
);
int nonzeroConstant=(
  reduce(subst(RB,{w2},0),GC)!=0
);
int nonzeroSlope=(
  reduce(diff(RB,{w2}),GC)!=0
);
list facRB=factorize(RB);
ideal factorsRB=facRB[1];
int irreducibleRepresentative=(size(factorsRB)==2);

print("DEGREE42_B_GEOMETRIC_FACTORIZATION");
print(exactB);
print(linearB);
print(nonzeroConstant);
print(nonzeroSlope);
print(irreducibleRepresentative);
"""
    values = run_singular(
        program,
        "DEGREE42_B_GEOMETRIC_FACTORIZATION",
        5,
    )
    assert values == ["1"] * 5, values


def certify_b_off_subresultant_witness() -> None:
    """Check an off-higher-gcd zero of the B residual modulo 103."""

    normals, bases, residuals, _defect = transformed_problem()
    e1, e2, translation, w0, w1, w2 = bases
    specialized = [
        item.subs(
            {e1: 1, e2: 1, translation: 21, w0: 0, w1: 0},
            simultaneous=True,
        )
        for item in residuals
    ]
    dummy = QuarticCase("b_witness", 103, (0, 0, 0, 0, 0, 0), 1, 1, 1)
    program = quartic_program(normals, specialized, dummy)
    old_ring = (
        f'ring q=103,({",".join(map(str, normals))},'
        "sync42hg_U,sync42hg_V),dp;"
    )
    new_ring = (
        f'ring q=103,({",".join(map(str, normals))},{w2},'
        "sync42hg_U,sync42hg_V),(dp(5),dp(1),dp(2));"
    )
    assert old_ring in program
    program = program.replace(old_ring, new_ring, 1)
    program = program[: program.index("poly g3=gcd(p3,q3);")]
    program += f"""
poly dp3=subst(diff(p3,sync42hg_U),
  sync42hg_U,3,sync42hg_V,-9);
poly dq3=subst(diff(q3,sync42hg_U),
  sync42hg_U,3,sync42hg_V,-9);
poly ep4=subst(p4,sync42hg_U,3,sync42hg_V,-9);
poly eq4=subst(q4,sync42hg_U,3,sync42hg_V,-9);
poly rho=dp3*eq4-dq3*ep4;
int residualShape=(
  rho==9*{w2}*(64+86*{w2})
);
int newZero=(subst(rho,{w2},28)==0);
print("DEGREE42_B_OFF_SUBRESULTANT_WITNESS");
print(residualShape);
print(newZero);
"""
    values = run_singular(
        program,
        "DEGREE42_B_OFF_SUBRESULTANT_WITNESS",
        2,
    )
    assert values == ["1", "1"], values


def main() -> None:
    certify_rational_factorizations()
    certify_geometric_b_factorization()
    certify_b_off_subresultant_witness()
    print(
        "PASS: rho_A=-75/512*t*e2^3*w2*R_A modulo A; the rational "
        "B norm has t^11*w2 times one irreducible representative, and "
        "each geometric B branch has t*w2*R_B with R_B linear in w2; "
        "both R_A and R_B have zeros away from the higher-gcd locus"
    )


if __name__ == "__main__":
    main()
