#!/usr/bin/env python3
"""Branchwise corrected-moment calculation on the s0 common boundary.

This checker continues after the exact mu_2, A, and B pivots.  It replaces
mu_4 and mu_5 by their generic finite-fiber Groebner basis on a chosen
L/Q stratum while preserving the base content of every relation.  On the
deepest L=Q=0 stratum it can also pseudo-reduce corrected later moments in
the exact rank-fifteen algebra without dividing by a base polynomial.
The exported system is equivalent on the explicitly recorded principal
open.

The script is intentionally branchwise.  It does not prove a corrected
zero-fiber theorem: the J-divisor still uses the separate exact rank-five
checker, and exact radical reconstruction of the later-moment ideals remains
open.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

import sympy as sp

from explore_two_pair_sic_bidegree33_full_anchor import (
    PARAMETERS,
    chart_expression,
    moment_terms,
    prepare_s0_branch_for_msolve,
)
from verify_two_pair_sic_bidegree33_boundary_generic_quotient import (
    exact_chart_expression,
    exact_moment_terms,
    substitute,
)


ROOT = Path(__file__).resolve().parents[1]
CORRECTED_ORDERS = (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14)
L_ADAPTED = "L"
Q_ADAPTED = "Q"
J_ADAPTED = "((99*Q+155*t0^2)^2+30420*L^2)"


def boundary_polynomials(
    singular: str,
    orders: tuple[int, ...],
    prime: int,
    timeout: int,
) -> list[str]:
    if prime:
        expressions = [
            chart_expression(moment_terms(order, prime), 0, prime)
            for order in orders
        ]
    else:
        expressions = [
            exact_chart_expression(exact_moment_terms(order))
            for order in orders
        ]
    variables, polynomials = prepare_s0_branch_for_msolve(
        singular,
        expressions,
        prime,
        "s0-boundary",
        timeout,
    )
    assert variables == (
        "s1",
        "s2",
        "s3",
        "s5",
        "t0",
        "t1",
        "t2",
        "t4",
    )
    assert len(polynomials) == len(orders) - 1
    return polynomials


def branch_data(
    branch: str,
) -> tuple[tuple[str, ...], tuple[tuple[str, str], ...], str, list[str]]:
    adapted = (
        ("t1", "(s1*t0-L)"),
        ("s2", "(s1^2-(13/3)*t0^2-Q)"),
    )
    if branch == "generic":
        return (
            ("s1", "s3", "t0", "L", "Q", "t2"),
            adapted,
            f"({L_ADAPTED})*({Q_ADAPTED})*({J_ADAPTED})",
            [],
        )
    if branch == "L":
        return (
            ("s1", "s3", "t0", "Q", "t2"),
            adapted + (("L", "(0)"),),
            "Q*(99*Q+155*t0^2)",
            [],
        )
    if branch == "Q":
        return (
            ("s1", "s3", "t0", "L", "t2"),
            adapted + (("Q", "(0)"),),
            f"L*({J_ADAPTED.replace('Q', '(0)')})",
            [],
        )
    if branch == "J":
        raise NotImplementedError(
            "Use verify_two_pair_sic_bidegree33_boundary_generic_quotient.py "
            "for the exact rank-five J-divisor algebra. Exporting corrected "
            "later moments requires its fraction-free quadratic-pair "
            "reduction; treating alpha as transcendental gives the wrong "
            "generic rank."
        )
    if branch == "LQ":
        return (
            ("s1", "s3", "t0", "t2"),
            adapted + (("L", "(0)"), ("Q", "(0)")),
            "1",
            [],
        )
    raise AssertionError(branch)


def finite_algebra_export(
    singular: str,
    polynomials: list[str],
    orders: tuple[int, ...],
    branch: str,
    prime: int,
    timeout: int,
) -> dict[str, object]:
    base_variables, replacements, prescribed_open, extra_equations = (
        branch_data(branch)
    )
    restricted = [
        substitute(polynomial, replacements) for polynomial in polynomials
    ]
    coefficient_specification = (
        str(prime) if prime else "0"
    ) + "," + ",".join(base_variables)
    ordinary_variables = base_variables + ("s5", "t4", "uinv")
    declarations = [
        f"poly p{order}={polynomial};"
        for order, polynomial in zip(orders[1:], restricted)
    ]
    later_orders = orders[4:]
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring fiber=({coefficient_specification}),(s5,t4),dp;
{chr(10).join(declarations)}
ideal G=std(p4,p5);
proc polynomialLcm(poly left, poly right)
{{
  return(left*right/gcd(left,right));
}}
proc clearPreservingContent(poly input)
{{
  poly commonDenominator=1;
  poly cursor=input;
  number coefficient;
  poly currentTerm;
  while(cursor!=0)
  {{
    coefficient=leadcoef(cursor);
    currentTerm=leadmonom(cursor);
    commonDenominator=polynomialLcm(
      commonDenominator,
      denominator(coefficient)
    );
    cursor=cursor-coefficient*currentTerm;
  }}
  return(list(commonDenominator*input,commonDenominator));
}}
poly denominatorPolynomial=1;
int polynomialIndex;
list cleared=clearPreservingContent(p3);
ideal exportedFiber=cleared[1];
denominatorPolynomial=polynomialLcm(denominatorPolynomial,cleared[2]);
for (polynomialIndex=1;polynomialIndex<=size(G);polynomialIndex++)
{{
  cleared=clearPreservingContent(G[polynomialIndex]);
  exportedFiber[size(exportedFiber)+1]=cleared[1];
  denominatorPolynomial=polynomialLcm(denominatorPolynomial,cleared[2]);
}}
{chr(10).join(
    f'cleared=clearPreservingContent(p{order}); '
    f'exportedFiber[size(exportedFiber)+1]=cleared[1]; '
    f'denominatorPolynomial=polynomialLcm('
    f'denominatorPolynomial,cleared[2]);'
    for order in later_orders
)}
print("META_GSIZE "+string(size(G)));
print("META_VDIM "+string(vdim(G)));
for (polynomialIndex=1;polynomialIndex<=size(G);polynomialIndex++)
{{
print("META_LEADEXP "+string(leadexp(G[polynomialIndex])));
}}
ring ordinary={prime},({",".join(ordinary_variables)}),dp;
poly openFactor=imap(fiber,denominatorPolynomial)
  *({prescribed_open});
print("EXPORT_OPEN "+string(openFactor));
ideal exported=imap(fiber,exportedFiber);
for (polynomialIndex=1;polynomialIndex<=size(exported);polynomialIndex++)
{{
  print("EXPORT_POLY "+string(exported[polynomialIndex]));
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    gsize_marker = re.search(r"(?m)^META_GSIZE (\d+)$", completed.stdout)
    vdim_marker = re.search(r"(?m)^META_VDIM (\d+)$", completed.stdout)
    assert gsize_marker is not None and vdim_marker is not None, (
        completed.stdout[-2000:]
    )
    gsize = int(gsize_marker.group(1))
    open_marker = re.search(r"(?m)^EXPORT_OPEN (.*)$", completed.stdout)
    assert open_marker is not None
    exported = re.findall(r"(?m)^EXPORT_POLY (.*)$", completed.stdout)
    assert len(exported) == 1 + gsize + len(later_orders), (
        len(exported),
        1 + gsize + len(later_orders),
        completed.stdout[-4000:],
        completed.stderr[-2000:],
    )
    exported.extend(extra_equations)
    open_factor = open_marker.group(1)
    if open_factor != "1":
        exported.append(f"uinv*({open_factor})-1")
    else:
        ordinary_variables = tuple(
            variable for variable in ordinary_variables if variable != "uinv"
        )
    return {
        "branch": branch,
        "base_variables": list(base_variables),
        "ordinary_variables": list(ordinary_variables),
        "groebner_basis_size": gsize,
        "quotient_length": int(vdim_marker.group(1)),
        "leading_exponents": re.findall(
            r"(?m)^META_LEADEXP ([0-9,]+)$",
            completed.stdout,
        ),
        "later_orders": list(later_orders),
        "later_corrected_moments_reduced": False,
        "open_factor": open_factor,
        "polynomials": exported,
    }


def run_msolve(
    msolve: str,
    export: dict[str, object],
    prime: int,
    timeout: int,
    linear_algebra: int,
    eliminate: int,
) -> dict[str, object]:
    variables = export["ordinary_variables"]
    polynomials = export["polynomials"]
    with tempfile.TemporaryDirectory(prefix="sic33-corrected-boundary-") as path:
        input_path = Path(path) / "system.ms"
        output_path = Path(path) / "result.ms"
        input_path.write_text(
            ",".join(variables)
            + "\n"
            + str(prime)
            + "\n"
            + ",\n".join(polynomials)
            + "\n"
        )
        started = time.monotonic()
        try:
            command = [
                    msolve,
                    "-f",
                    str(input_path),
                    "-o",
                    str(output_path),
                    "-t",
                    "4",
                    "-v",
                    "1",
                    "-l",
                    str(linear_algebra),
                    "-g",
                    "1",
                ]
            if eliminate:
                command.extend(["-e", str(eliminate)])
            completed = subprocess.run(
                command,
                text=True,
                capture_output=True,
                check=False,
                timeout=timeout,
            )
        except subprocess.TimeoutExpired:
            return {
                "status": "timeout",
                "seconds": float(timeout),
                "returncode": None,
                "stdout_tail": "",
                "stderr_tail": "",
                "result_tail": "",
            }
        elapsed = time.monotonic() - started
        result = output_path.read_text() if output_path.exists() else ""
    if completed.returncode != 0:
        status = f"solver-error-{completed.returncode}"
    elif result.rstrip().endswith("[1]:"):
        status = "unit"
    elif result.rstrip().endswith("[-1]:"):
        status = "unit"
    else:
        status = "nonunit"
    return {
        "status": status,
        "seconds": round(elapsed, 3),
        "returncode": completed.returncode,
        "stdout_tail": completed.stdout[-1000:],
        "stderr_tail": completed.stderr[-1000:],
        "result_tail": result[-1000:],
    }


def run_modstd(
    singular: str,
    export: dict[str, object],
    timeout: int,
) -> dict[str, object]:
    variables = export["ordinary_variables"]
    polynomials = export["polynomials"]
    started = time.monotonic()
    try:
        completed = subprocess.run(
            [singular, "-q"],
            input=f"""
LIB "modstd.lib";
ring exact=0,({",".join(variables)}),dp;
option(redSB);
ideal I={",".join(polynomials)};
ideal G=modStd(I,1);
print(
  "EXACT "+string(dim(G))+" "+string(size(G))+" "
  +string(vdim(G))+" "+string(G[1]==1)
);
""",
            text=True,
            capture_output=True,
            check=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return {
            "status": "timeout",
            "seconds": float(timeout),
            "algorithm": "Singular modStd with exactness=1",
            "stdout_tail": "",
            "stderr_tail": "",
        }
    elapsed = time.monotonic() - started
    marker = re.search(
        r"(?m)^EXACT (-?\d+) (\d+) (-?\d+) ([01])$",
        completed.stdout,
    )
    if marker is None:
        raise AssertionError(completed.stdout[-8000:])
    dimension, size, quotient_length, unit = marker.groups()
    return {
        "status": "unit" if unit == "1" else "nonunit",
        "seconds": round(elapsed, 3),
        "dimension": int(dimension),
        "basis_size": int(size),
        "quotient_length": int(quotient_length),
        "algorithm": "Singular modStd with exactness=1",
        "stdout_tail": completed.stdout[-1000:],
        "stderr_tail": completed.stderr[-1000:],
    }


def deepest_three_pivot_certificate(
    singular: str,
    polynomials: list[str],
    timeout: int,
) -> dict[str, object]:
    replacements = (
        ("t1", "(s1*t0-L)"),
        ("s2", "(s1^2-(13/3)*t0^2-Q)"),
        ("L", "(0)"),
        ("Q", "(0)"),
    )
    p3, p4, p5 = [
        substitute(polynomial, replacements)
        for polynomial in polynomials[:3]
    ]
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring deepest=(0,s1,s3,t0),(t2,s5,t4),dp;
option(redSB);
poly p3={p3};
poly p4={p4};
poly p5={p5};
ideal I=p3,p4,p5;
ideal G=std(I);
print("META "+string(size(G))+" "+string(vdim(G)));
print("P3LEAD "+string(leadexp(p3))+" "+string(leadcoef(p3)));
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEADEXP "+string(leadexp(G[basisIndex])));
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    marker = re.search(r"(?m)^META (\d+) (\d+)$", completed.stdout)
    p3_marker = re.search(
        r"(?m)^P3LEAD ([0-9,]+) (.*)$",
        completed.stdout,
    )
    assert marker is not None and marker.groups() == ("11", "15"), (
        completed.stdout[-4000:],
        completed.stderr[-1000:],
    )
    assert p3_marker is not None and p3_marker.group(1) == "3,0,0"
    leading_exponents = re.findall(
        r"(?m)^LEADEXP ([0-9,]+)$",
        completed.stdout,
    )
    assert leading_exponents == [
        "1,1,1",
        "2,0,1",
        "1,2,0",
        "2,1,0",
        "3,0,0",
        "0,0,4",
        "0,1,3",
        "1,0,3",
        "0,2,2",
        "0,3,1",
        "0,4,0",
    ]
    return {
        "coefficient_field": "QQ(s1,s3,t0)",
        "fiber_variables": ["t2", "s5", "t4"],
        "moments": [3, 4, 5],
        "groebner_basis_size": 11,
        "quotient_length": 15,
        "mu3_leading_exponent": [3, 0, 0],
        "mu3_leading_coefficient": p3_marker.group(2),
        "leading_exponents": leading_exponents,
        "scope": (
            "exact characteristic-zero finite algebra after the mu3 pivot; "
            "later corrected moments are not yet reduced fraction-free"
        ),
    }


def deepest_fraction_free_normal_forms(
    singular: str,
    polynomials: list[str],
    orders: tuple[int, ...],
    prime: int,
    timeout: int,
) -> dict[str, object]:
    """Pseudo-reduce later moments in the L=Q=0 rank-fifteen algebra.

    The Groebner basis is first cleared coefficient by coefficient.  The
    reduction then cancels a leading term by cross multiplication, rather
    than division by the leading coefficient.  Thus no base divisor is
    discarded.  The returned normal forms are deliberately summarized by
    support and step count: their coefficients are far too large to be a
    useful JSON certificate.
    """

    replacements = (
        ("t1", "(s1*t0-L)"),
        ("s2", "(s1^2-(13/3)*t0^2-Q)"),
        ("L", "(0)"),
        ("Q", "(0)"),
    )
    restricted = [
        substitute(polynomial, replacements) for polynomial in polynomials
    ]
    declarations = "\n".join(
        f"poly p{order}={polynomial};"
        for order, polynomial in zip(orders[1:], restricted)
    )
    later_orders = orders[4:]
    reductions = "\n".join(
        f"""
reductionData=ffReduce(p{order},clearedBasis);
remainderPolynomial=reductionData[1];
commonContent=0;
cursor=remainderPolynomial;
while(cursor!=0)
{{
  coefficient=leadcoef(cursor);
  commonContent=gcd(
    commonContent,
    numerator(coefficient)
  );
  cursor=cursor-lead(cursor);
}}
print(
  "FFNF {order} "+string(reductionData[2])+" "
  +string(size(remainderPolynomial))+" "
  +string(deg(remainderPolynomial))+" "
  +string(commonContent==1)
);
"""
        for order in later_orders
    )
    coefficient_specification = (
        str(prime) if prime else "0"
    ) + ",s1,s3,t0"
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring deepest=({coefficient_specification}),(t2,s5,t4),dp;
option(redSB);
{declarations}
ideal threePivot=p3,p4,p5;
ideal groebnerBasis=std(threePivot);
proc polynomialLcm(poly left, poly right)
{{
  return(left*right/gcd(left,right));
}}
proc clearPreservingContent(poly input)
{{
  poly commonDenominator=1;
  poly cursor=input;
  number coefficient;
  poly currentTerm;
  while(cursor!=0)
  {{
    coefficient=leadcoef(cursor);
    currentTerm=leadmonom(cursor);
    commonDenominator=polynomialLcm(
      commonDenominator,
      denominator(coefficient)
    );
    cursor=cursor-coefficient*currentTerm;
  }}
  return(list(commonDenominator*input,commonDenominator));
}}
ideal clearedBasis;
poly denominatorSupport=1;
list cleared;
int basisIndex;
for (
  basisIndex=1;
  basisIndex<=size(groebnerBasis);
  basisIndex++
)
{{
  cleared=clearPreservingContent(groebnerBasis[basisIndex]);
  clearedBasis[basisIndex]=cleared[1];
  denominatorSupport=polynomialLcm(
    denominatorSupport,
    cleared[2]
  );
}}
proc ffReduce(poly input, ideal divisors)
{{
  poly active=input;
  poly remainder=0;
  poly quotientMonomial;
  number activeCoefficient;
  number divisorCoefficient;
  int divisorIndex;
  int found;
  int steps=0;
  while(active!=0)
  {{
    found=0;
    for (
      divisorIndex=1;
      divisorIndex<=size(divisors);
      divisorIndex++
    )
    {{
      if (
        leadmonom(active)
        /leadmonom(divisors[divisorIndex])!=0
      )
      {{
        activeCoefficient=leadcoef(active);
        divisorCoefficient=leadcoef(divisors[divisorIndex]);
        quotientMonomial=
          leadmonom(active)
          /leadmonom(divisors[divisorIndex]);
        active=
          divisorCoefficient*active
          -activeCoefficient
           *quotientMonomial
           *divisors[divisorIndex];
        remainder=divisorCoefficient*remainder;
        found=1;
        steps++;
        break;
      }}
    }}
    if(found==0)
    {{
      remainder=remainder+lead(active);
      active=active-lead(active);
    }}
  }}
  return(list(remainder,steps));
}}
print(
  "FFALG "+string(size(groebnerBasis))+" "
  +string(vdim(groebnerBasis))
);
list reductionData;
poly remainderPolynomial;
poly cursor;
poly commonContent;
number coefficient;
{reductions}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    algebra_marker = re.search(r"(?m)^FFALG (\d+) (\d+)$", completed.stdout)
    assert algebra_marker is not None
    assert algebra_marker.groups() == ("11", "15")
    normal_forms: dict[str, object] = {}
    for order, steps, support, fiber_degree, primitive in re.findall(
        r"(?m)^FFNF (\d+) (\d+) (\d+) (-?\d+) ([01])$",
        completed.stdout,
    ):
        normal_forms[order] = {
            "pseudo_reduction_steps": int(steps),
            "support_size": int(support),
            "maximum_fiber_degree": int(fiber_degree),
            "primitive_base_content": primitive == "1",
        }
    assert set(normal_forms) == {str(order) for order in later_orders}, (
        completed.stdout[-4000:],
        completed.stderr[-2000:],
    )
    return {
        "algorithm": "content-preserving fraction-free pseudo-normal form",
        "coefficient_field": (
            f"F_{prime}(s1,s3,t0)" if prime else "QQ(s1,s3,t0)"
        ),
        "groebner_basis_size": 11,
        "quotient_length": 15,
        "normal_forms": normal_forms,
        "scope": (
            "exact reductions in the finite algebra; exceptional leading-"
            "coefficient support and the later radical are not yet resolved"
        ),
    }


def t0_zero_branch_certificate(
    msolve: str,
    polynomials: list[str],
    orders: tuple[int, ...],
    branch: str,
    timeout: int,
) -> dict[str, object]:
    """Exclude one t0=0 branch divisor over QQ."""

    assert 10 in orders
    available = dict(zip(orders[1:], polynomials))
    adapted = (
        ("t1", "(s1*t0-L)"),
        ("s2", "(s1^2-(13/3)*t0^2-Q)"),
        ("t0", "(0)"),
    )
    extra_equations: list[str] = []
    if branch == "generic":
        replacements = adapted
        open_factor = "L*Q*(9801*Q^2+30420*L^2)"
        variables = ["t2", "s5", "t4", "s1", "s3", "L", "Q"]
    elif branch == "L":
        replacements = adapted + (("L", "(0)"),)
        open_factor = "Q"
        variables = ["t2", "s5", "t4", "s1", "s3", "Q"]
    elif branch == "Q":
        replacements = adapted + (("Q", "(0)"),)
        open_factor = "L"
        variables = ["t2", "s5", "t4", "s1", "s3", "L"]
    elif branch == "J":
        replacements = adapted
        open_factor = "L*Q"
        extra_equations.append("9801*Q^2+30420*L^2")
        variables = ["t2", "s5", "t4", "s1", "s3", "L", "Q"]
    elif branch == "LQ":
        replacements = adapted + (("L", "(0)"), ("Q", "(0)"))
        open_factor = "1"
        variables = ["t2", "s5", "t4", "s1", "s3"]
    else:
        raise AssertionError(branch)
    equations = [
        substitute(available[order], replacements)
        for order in range(3, 11)
    ]
    equations.extend(extra_equations)
    if open_factor != "1":
        variables.append("uinv")
        equations.append(f"uinv*({open_factor})-1")
    with tempfile.TemporaryDirectory(
        prefix=f"sic33-corrected-{branch}-t0-zero-"
    ) as path:
        input_path = Path(path) / "system.ms"
        output_path = Path(path) / "result.ms"
        input_path.write_text(
            ",".join(variables)
            + "\n"
            "0\n"
            + ",\n".join(equations)
            + "\n"
        )
        started = time.monotonic()
        completed = subprocess.run(
            [
                msolve,
                "-f",
                str(input_path),
                "-o",
                str(output_path),
                "-t",
                "4",
                "-v",
                "1",
                "-l",
                "2",
                "-g",
                "1",
            ],
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
        )
        elapsed = time.monotonic() - started
        result = output_path.read_text() if output_path.exists() else ""
    unit = (
        completed.returncode == 0
        and result.rstrip().endswith("[1]:")
        and "#length of basis:      1 element" in result
    )
    assert unit, (
        completed.returncode,
        completed.stdout[-2000:],
        completed.stderr[-2000:],
        result[-2000:],
    )
    return {
        "characteristic": 0,
        "branch": branch,
        "divisor": "t0=0",
        "open_factor": open_factor,
        "variables": variables,
        "moments": list(range(3, 11)),
        "unit_ideal": True,
        "basis_size": 1,
        "backend": "msolve exact sparse linear algebra 2",
        "seconds": round(elapsed, 3),
        "scope": (
            "exact exclusion of this t0=0 branch divisor; the t0-open "
            "and the full corrected radical remain unresolved"
        ),
    }


def t0_open_rank_six_certificate(
    singular: str,
    timeout: int,
) -> dict[str, object]:
    """Construct the generic rank-six fibre after normalizing t0 to one.

    The residual torus sends every point of the s0-chart with t0 nonzero to
    this chart.  We retain u=s0^{-1}, use mu2 to remove t4, and use the
    homogenized A=B=0 equations to remove t3 and s4.  The remaining mu4,
    mu5 fibre is taken in (s6,s5); this reverse order has substantially
    smaller leading-coefficient factors than (s5,s6).
    """

    symbols = sp.symbols(" ".join(PARAMETERS))
    (
        s0,
        s1,
        s2,
        s3,
        s4,
        s5,
        s6,
        t0,
        t1,
        t2,
        t3,
        t4,
    ) = symbols
    u = sp.symbols("u")

    def full_moment(order: int) -> sp.Expr:
        return sum(
            coefficient
            * sp.prod(
                variable**exponent
                for variable, exponent in zip(
                    symbols,
                    exponent_tuple,
                )
            )
            for exponent_tuple, coefficient in exact_moment_terms(order).items()
        ).subs(t0, 1)

    moments = {
        order: full_moment(order)
        for order in (2, 3, 4, 5, 6)
    }
    mu2_coefficient = sp.diff(moments[2], t4)
    assert mu2_coefficient == 336
    t4_value = sp.cancel(
        -(moments[2] - mu2_coefficient * t4) / mu2_coefficient
    )
    homogenized_a = (
        6 * s1**2 * t1
        - 3 * s1 * s2
        - 3 * s0 * s1 * t2
        - 3 * s0 * s2 * t1
        + 2 * s0 * s3
        - 3 * s0
        + s0**2 * t3
    )
    homogenized_b = (
        12 * s0 * s1 * s3
        + 28 * s1 * t1
        - 18 * s0 * s1
        - 9 * s0 * s2**2
        - 14 * s2
        - 3 * s0**2 * s4
        - 2 * s0 * t2
        - 12 * s0 * t1**2
    )
    t3_value = sp.cancel(
        -(homogenized_a - s0**2 * t3) / s0**2
    )
    s4_value = sp.cancel(
        (homogenized_b + 3 * s0**2 * s4) / (3 * s0**2)
    )
    transformed: dict[int, sp.Poly] = {}
    cleared_denominators: dict[str, str] = {}
    polynomial_variables = (s6, s5, s1, s2, s3, t1, t2, u)
    for order in (3, 4, 5, 6):
        expression = sp.cancel(
            moments[order]
            .subs(t4, t4_value)
            .subs(t3, t3_value)
            .subs(s4, s4_value)
            .subs(s0, 1 / u)
        )
        numerator, denominator = sp.fraction(expression)
        transformed[order] = sp.Poly(
            numerator,
            *polynomial_variables,
            domain=sp.QQ,
        )
        cleared_denominators[str(order)] = str(denominator)
    assert transformed[3].degree(s6) == 0
    assert transformed[3].degree(s5) == 0
    assert len(transformed[3].terms()) == 73

    def singular_expression(polynomial: sp.Poly) -> str:
        return sp.sstr(polynomial.as_expr()).replace("**", "^")

    declarations = "\n".join(
        f"poly p{order}={singular_expression(transformed[order])};"
        for order in (4, 5, 6)
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring t0open=(0,s1,s2,s3,t1,t2,u),(s6,s5),dp;
option(redSB);
{declarations}
ideal G=std(p4,p5);
poly r6=reduce(p6,G);
poly x=s1^2*u-s2;
poly ell=s1*u-t1;
poly Qfactor=3*x-13*u;
poly Afactor=99*x-274*u;
poly Jfactor=(99*x-274*u)^2+30420*ell^2;
poly Kfactor=351*x-901*u;
poly Hfactor=(99*x-274*u)*Kfactor+121680*ell^2;
print(
  "META "+string(size(G))+" "+string(vdim(G))+" "
  +string(size(r6))+" "+string(deg(r6))
);
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEAD "+string(leadexp(G[basisIndex])));
}}
poly cursor=r6;
while(cursor!=0)
{{
  print("NFEXP "+string(leadexp(cursor)));
  cursor=cursor-lead(cursor);
}}
print(
  "LC "+string(leadcoef(G[1])/Kfactor)+" "
  +string(leadcoef(G[2])/Hfactor)+" "
  +string(
    leadcoef(G[3])/(Qfactor*Kfactor*Jfactor*Hfactor)
  )
);
print(
  "IDENT "+string(Kfactor-(4*Afactor-15*Qfactor)==0)+" "
  +string(Hfactor-(4*Jfactor-15*Afactor*Qfactor)==0)
);
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    marker = re.search(r"(?m)^META (\d+) (\d+) (\d+) (\d+)$", completed.stdout)
    assert marker is not None and marker.groups() == ("3", "6", "6", "3"), (
        completed.stdout[-4000:],
        completed.stderr[-1000:],
    )
    leading_exponents = re.findall(
        r"(?m)^LEAD ([0-9,]+)$",
        completed.stdout,
    )
    assert leading_exponents == ["2,0", "1,2", "0,4"]
    normal_form_exponents = re.findall(
        r"(?m)^NFEXP ([0-9,]+)$",
        completed.stdout,
    )
    assert set(normal_form_exponents) == {
        "0,0",
        "1,0",
        "0,1",
        "1,1",
        "0,2",
        "0,3",
    }
    leading_factor_constants = re.search(
        r"(?m)^LC (\S+) (\S+) (\S+)$",
        completed.stdout,
    )
    assert (
        leading_factor_constants is not None
        and leading_factor_constants.groups() == ("311040", "324", "122472")
    ), completed.stdout[-2000:]
    identity_marker = re.search(r"(?m)^IDENT ([01]) ([01])$", completed.stdout)
    assert identity_marker is not None and identity_marker.groups() == ("1", "1")

    k_substitution = {
        s2: s1**2 * u - sp.Rational(901, 351) * u,
    }

    def specialized_polynomial(
        order: int,
        replacements: dict[sp.Symbol, sp.Expr],
        variables: tuple[sp.Symbol, ...],
    ) -> sp.Poly:
        polynomial = sp.Poly(
            sp.expand(
                transformed[order].as_expr().subs(replacements)
            ),
            *variables,
            domain=sp.QQ,
        )
        return polynomial.clear_denoms(convert=True)[1]

    k_variables = (s6, s5, s1, s3, t1, t2, u)
    k_polynomials = {
        order: specialized_polynomial(
            order,
            k_substitution,
            k_variables,
        )
        for order in (4, 5, 6)
    }
    k_declarations = "\n".join(
        f"poly p{order}={singular_expression(k_polynomials[order])};"
        for order in (4, 5, 6)
    )
    k_completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring kbranch=(0,s1,s3,t1,t2,u),(s6,s5),dp;
option(redSB);
{k_declarations}
ideal G=std(p4,p5);
poly r6=reduce(p6,G);
poly ell=s1*u-t1;
poly Jfactor=(-775*u/39)^2+30420*ell^2;
print(
  "META "+string(size(G))+" "+string(vdim(G))+" "
  +string(size(r6))+" "+string(deg(r6))
);
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEAD "+string(leadexp(G[basisIndex])));
}}
print(
  "KLC "+string(leadcoef(G[1])/(u*ell))+" "
  +string(leadcoef(G[2])/ell^3)+" "
  +string(leadcoef(G[3])/(u^2*ell^2*Jfactor))
);
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    k_marker = re.search(
        r"(?m)^META (\d+) (\d+) (\d+) (\d+)$",
        k_completed.stdout,
    )
    assert k_marker is not None and k_marker.groups() == ("3", "6", "6", "3")
    k_leading_exponents = re.findall(
        r"(?m)^LEAD ([0-9,]+)$",
        k_completed.stdout,
    )
    assert k_leading_exponents == ["1,1", "3,0", "0,4"]
    k_factor_constants = re.search(
        r"(?m)^KLC (-?\d+) (-?\d+) (-?\d+)$",
        k_completed.stdout,
    )
    assert k_factor_constants is not None, k_completed.stdout[-3000:]

    kh_substitution = {
        **k_substitution,
        t1: s1 * u,
    }
    kh_variables = (s6, s5, s1, s3, t2, u)
    kh_polynomials = {
        order: specialized_polynomial(
            order,
            kh_substitution,
            kh_variables,
        )
        for order in (4, 5, 6)
    }
    kh_declarations = "\n".join(
        f"poly p{order}={singular_expression(kh_polynomials[order])};"
        for order in (4, 5, 6)
    )
    kh_completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring khbranch=(0,s1,s3,t2,u),(s6,s5),dp;
option(redSB);
{kh_declarations}
ideal G=std(p4,p5);
poly r6=reduce(p6,G);
print(
  "META "+string(size(G))+" "+string(vdim(G))+" "
  +string(size(r6))+" "+string(deg(r6))
);
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print(
    "LEAD "+string(leadexp(G[basisIndex]))+" "
    +string(leadcoef(G[basisIndex]))
  );
}}
poly cursor=r6;
while(cursor!=0)
{{
  print("NFEXP "+string(leadexp(cursor)));
  cursor=cursor-lead(cursor);
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    kh_marker = re.search(
        r"(?m)^META (\d+) (\d+) (\d+) (\d+)$",
        kh_completed.stdout,
    )
    assert kh_marker is not None and kh_marker.groups() == ("2", "6", "6", "3")
    kh_leading_data = re.findall(
        r"(?m)^LEAD ([0-9,]+) (.*)$",
        kh_completed.stdout,
    )
    assert [exponent for exponent, _ in kh_leading_data] == ["0,2", "3,0"]
    assert re.fullmatch(
        r"\d+\*u\^3",
        kh_leading_data[0][1].strip("()"),
    )
    assert re.fullmatch(r"\d+", kh_leading_data[1][1].strip("()"))
    kh_normal_form_exponents = re.findall(
        r"(?m)^NFEXP ([0-9,]+)$",
        kh_completed.stdout,
    )
    assert kh_normal_form_exponents == [
        "2,1",
        "2,0",
        "1,1",
        "1,0",
        "0,1",
        "0,0",
    ]

    h_completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring hsource=(0,s1,s3,t2,u,r),(s6,s5,s2,t1),dp;
poly p4={singular_expression(transformed[4])};
poly p5={singular_expression(transformed[5])};
poly D=1287*r^2+40560;
poly ell=-775*r*u/D;
poly q=-155*u/33+r*ell;
poly x=(q+13*u)/3;
p4=subst(subst(p4,s2,s1^2*u-x),t1,s1*u-ell);
p5=subst(subst(p5,s2,s1^2*u-x),t1,s1*u-ell);
ring hbranch=(0,s1,s3,t2,u,r),(s6,s5),dp;
option(redSB);
poly q4=imap(hsource,p4);
poly q5=imap(hsource,p5);
ideal G=std(q4,q5);
print(
  "META "+string(size(G))+" "+string(vdim(G))+" "+string(dim(G))
);
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEAD "+string(leadexp(G[basisIndex])));
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    h_marker = re.search(
        r"(?m)^META (\d+) (\d+) (-?\d+)$",
        h_completed.stdout,
    )
    assert h_marker is not None and h_marker.groups() == ("2", "6", "0")
    h_leading_exponents = re.findall(
        r"(?m)^LEAD ([0-9,]+)$",
        h_completed.stdout,
    )
    assert h_leading_exponents == ["2,0", "0,3"]
    return {
        "chart": "t0=1 with u=s0^-1",
        "characteristic": 0,
        "mu2_t4_pivot_coefficient": 336,
        "eliminated_variables": ["t4", "t3", "s4"],
        "base_variables": ["s1", "s2", "s3", "t1", "t2", "u"],
        "fiber_variables": ["s6", "s5"],
        "mu3_base_term_count": 73,
        "mu3_independent_of_fiber": True,
        "cleared_denominators": cleared_denominators,
        "groebner_basis_size": 3,
        "quotient_length": 6,
        "leading_exponents": leading_exponents,
        "standard_basis": ["1", "s6", "s5", "s6*s5", "s5^2", "s5^3"],
        "mu6_normal_form": {
            "support_size": 6,
            "maximum_fiber_degree": 3,
            "support_exponents": normal_form_exponents,
        },
        "adapted_base_coordinates": {
            "x": "s1^2*u-s2",
            "ell": "s1*u-t1",
        },
        "leading_coefficient_factors": {
            "A": "99*x-274*u",
            "Q": "3*x-13*u",
            "J": "(99*x-274*u)^2+30420*ell^2",
            "K": "351*x-901*u",
            "H": "(99*x-274*u)*K+121680*ell^2",
            "identities": [
                "K=4*A-15*Q",
                "H=4*J-15*A*Q",
            ],
            "constants": [311040, 324, 122472],
            "factorization": ["K", "H", "Q*K*J*H"],
        },
        "K_divisor": {
            "substitution": "s2=s1^2*u-(901/351)*u",
            "H": "121680*ell^2",
            "groebner_basis_size": 3,
            "quotient_length": 6,
            "leading_exponents": k_leading_exponents,
            "leading_coefficient_factors": [
                "u*ell",
                "ell^3",
                "u^2*ell^2*J",
            ],
            "mu6_normal_form_support_size": 6,
            "scope": "exact on the generic K=0 divisor away from ell*J=0",
        },
        "K_H_intersection": {
            "radical_substitutions": [
                "s2=s1^2*u-(901/351)*u",
                "t1=s1*u",
            ],
            "derivation": "K=0 and H=121680*ell^2",
            "groebner_basis_size": 2,
            "quotient_length": 6,
            "leading_exponents": ["0,2", "3,0"],
            "standard_basis": [
                "1",
                "s6",
                "s6^2",
                "s5",
                "s6*s5",
                "s6^2*s5",
            ],
            "mu6_normal_form": {
                "support_size": 6,
                "support_exponents": kh_normal_form_exponents,
            },
            "scope": (
                "exact finite algebra on the reduced K=H=0 intersection; "
                "later common-root equations remain unresolved"
            ),
        },
        "H_divisor": {
            "parameter": "r",
            "denominator": "D=1287*r^2+40560",
            "parametrization": {
                "ell": "-775*r*u/D",
                "Q": "-155*u/33+r*ell",
                "x": "(Q+13*u)/3",
                "s2": "s1^2*u-x",
                "t1": "s1*u-ell",
            },
            "omitted_base_point": "A=ell=0, contained in J=0",
            "groebner_basis_size": 2,
            "quotient_length": 6,
            "leading_exponents": h_leading_exponents,
            "standard_basis": [
                "1",
                "s6",
                "s5",
                "s6*s5",
                "s5^2",
                "s6*s5^2",
            ],
            "scope": (
                "exact generic finite algebra on H=0; later-moment "
                "normal forms and lower-dimensional specializations "
                "remain unresolved"
            ),
        },
        "scope": (
            "exact generic finite algebra and exact mu6 normal form on the "
            "t0-open chart; coefficient-denominator strata and the later "
            "radical remain unresolved"
        ),
    }


def fixed_chart_expression_exact(
    terms: dict[tuple[int, ...], int],
    fixed_index: int,
) -> str:
    """Serialize an exact moment after setting one parameter equal to one."""

    combined: dict[tuple[int, ...], int] = {}
    for exponents, coefficient in terms.items():
        reduced = exponents[:fixed_index] + exponents[fixed_index + 1 :]
        combined[reduced] = combined.get(reduced, 0) + coefficient
    variable_names = PARAMETERS[:fixed_index] + PARAMETERS[fixed_index + 1 :]
    serialized: list[str] = []
    for exponents, coefficient in sorted(combined.items()):
        if coefficient == 0:
            continue
        factors: list[str] = []
        for variable, exponent in zip(variable_names, exponents):
            if exponent == 1:
                factors.append(variable)
            elif exponent > 1:
                factors.append(f"{variable}^{exponent}")
        monomial = "*".join(factors)
        if not monomial:
            serialized.append(str(coefficient))
        elif coefficient == 1:
            serialized.append(monomial)
        elif coefficient == -1:
            serialized.append(f"-{monomial}")
        else:
            serialized.append(f"{coefficient}*{monomial}")
    return "+".join(serialized).replace("+-", "-") or "0"


def t0_open_direct_export(
    orders: tuple[int, ...],
    prime: int,
) -> dict[str, object]:
    """Keep the common-root incidence before any rank-six localization.

    Setting ``t0=1`` is valid on the residual-torus quotient of the
    ``s0*t0 != 0`` chart.  The equation ``s0*uinv=1`` retains precisely
    the original ``s0`` principal open.  Keeping A, B, mu2, and every later
    moment as equations avoids the invalid stronger condition that every
    coordinate of a later moment vanish in the rank-six algebra.
    """

    assert orders and orders[0] == 2
    fixed_index = PARAMETERS.index("t0")
    if prime:
        moments = [
            chart_expression(
                moment_terms(order, prime),
                fixed_index,
                prime,
            )
            for order in orders
        ]
    else:
        moments = [
            fixed_chart_expression_exact(
                exact_moment_terms(order),
                fixed_index,
            )
            for order in orders
        ]
    homogenized_a = (
        "6*s1^2*t1-3*s1*s2-3*s0*s1*t2-3*s0*s2*t1"
        "+2*s0*s3-3*s0+s0^2*t3"
    )
    homogenized_b = (
        "12*s0*s1*s3+28*s1*t1-18*s0*s1-9*s0*s2^2"
        "-14*s2-3*s0^2*s4-2*s0*t2-12*s0*t1^2"
    )
    variables = [
        "s6",
        "s5",
        "t4",
        "t3",
        "s4",
        "s0",
        "s1",
        "s2",
        "s3",
        "t1",
        "t2",
        "uinv",
    ]
    polynomials = [
        homogenized_a,
        homogenized_b,
        *moments,
        "s0*uinv-1",
    ]
    return {
        "branch": "t0-open common-root incidence",
        "characteristic": prime,
        "ordinary_variables": variables,
        "polynomials": polynomials,
        "orders": list(orders),
        "equation_count": len(polynomials),
        "open_factor": "s0",
        "scope": (
            "specialization-safe direct incidence on s0*t0 != 0; "
            "no rank-six leading coefficient has been inverted"
        ),
    }


def t0_open_reduced_export(
    singular: str,
    orders: tuple[int, ...],
    prime: int,
    timeout: int,
) -> dict[str, object]:
    """Apply the three linear pivots while retaining ``s0*uinv-1``.

    Unlike the rational-function-field rank-six presentation, this
    polynomial presentation does not invert Q, J, K, H, or any leading
    coefficient.  It is therefore valid simultaneously on all their
    specializations.
    """

    assert orders and orders[0] == 2
    fixed_index = PARAMETERS.index("t0")
    if prime:
        moment_expressions = {
            order: chart_expression(
                moment_terms(order, prime),
                fixed_index,
                prime,
            )
            for order in orders
            if order >= 3
        }
    else:
        moment_expressions = {
            order: fixed_chart_expression_exact(
                exact_moment_terms(order),
                fixed_index,
            )
            for order in orders
            if order >= 3
        }
    declarations = "\n".join(
        f"poly p{order}={expression};"
        for order, expression in moment_expressions.items()
    )
    transformations = "\n".join(
        f"""
poly q{order}=subst(p{order},t4,t4Value);
q{order}=subst(q{order},t3,t3Value);
q{order}=subst(q{order},s4,s4Value);
q{order}=reduce(q{order},inverseBasis);
print("EXPORT {order} "+string(q{order}));
"""
        for order in moment_expressions
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring source={prime},(
  s6,s5,t4,t3,s4,s0,s1,s2,s3,t1,t2,uinv
),dp;
{declarations}
poly aWithoutT3=
  6*s1^2*t1-3*s1*s2-3*s0*s1*t2-3*s0*s2*t1
  +2*s0*s3-3*s0;
poly bWithoutS4=
  12*s0*s1*s3+28*s1*t1-18*s0*s1-9*s0*s2^2
  -14*s2-2*s0*t2-12*s0*t1^2;
poly t3Value=-aWithoutT3*uinv^2;
poly s4Value=(bWithoutS4*uinv^2)/3;
poly t4Value=(
  3*s0*s6-18*s1*s5+45*s2*s4-30*s3^2
  +56*t1*t3-42*t2^2-70
)/14;
ideal inverseBasis=std(s0*uinv-1);
{transformations}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    assert "\n   ? " not in completed.stdout, (
        completed.stdout[-8000:],
        completed.stderr[-2000:],
    )
    exported = {
        int(order): polynomial
        for order, polynomial in re.findall(
            r"(?m)^EXPORT (\d+) (.*)$",
            completed.stdout,
        )
    }
    assert set(exported) == set(moment_expressions), (
        completed.stdout[-4000:],
        completed.stderr[-2000:],
    )
    for order, polynomial in exported.items():
        assert not re.search(r"\b(?:t3|s4|t4)\b", polynomial), (
            order,
            polynomial[-2000:],
        )
    variables = [
        "s6",
        "s5",
        "s0",
        "s1",
        "s2",
        "s3",
        "t1",
        "t2",
        "uinv",
    ]
    polynomials = [
        exported[order] for order in moment_expressions
    ] + ["s0*uinv-1"]
    return {
        "branch": "t0-open specialization-safe reduced incidence",
        "characteristic": prime,
        "ordinary_variables": variables,
        "polynomials": polynomials,
        "orders": list(orders),
        "eliminated_variables": ["t3", "s4", "t4"],
        "equation_count": len(polynomials),
        "open_factor": "s0",
        "scope": (
            "the A, B, and mu2 pivots are applied polynomially modulo "
            "s0*uinv-1; Q, J, K, H and all fibre leading coefficients "
            "remain uninverted"
        ),
    }


def t0_open_fitting_export(
    singular: str,
    orders: tuple[int, ...],
    prime: int,
    timeout: int,
) -> dict[str, object]:
    """Export a common-root Fitting system from the generic rank-six algebra."""

    assert prime == 0 or prime > 7
    multiplication_orders = tuple(order for order in orders if order >= 6)
    assert multiplication_orders and multiplication_orders[0] == 6
    reduced = t0_open_reduced_export(
        singular,
        orders,
        prime,
        timeout,
    )
    available = dict(zip(orders[1:], reduced["polynomials"][:-1]))
    localized = {
        order: substitute(
            polynomial,
            (("s0", "(1/u)"), ("uinv", "(u)")),
        )
        for order, polynomial in available.items()
    }
    declarations = "\n".join(
        f"poly p{order}={polynomial};"
        for order, polynomial in localized.items()
    )
    matrix_declarations = "\n".join(
        f"matrix M{order}[6][6];"
        for order in multiplication_orders
    )
    normal_form_declarations = "\n".join(
        f"poly r{order}=reduce(p{order},G);"
        for order in multiplication_orders
    )
    multiplication_matrix_fills = "\n".join(
        f"""
z=r{order};
while(z!=0)
{{
  basisRow=coordinateIndex(leadmonom(z));
  if(basisRow==0)
  {{
    print("COORDINATE_ERROR {order} "+string(leadmonom(z)));
    exit;
  }}
  if(basisRow==1) {{ M{order}=M{order}+leadcoef(z)*B1; }}
  if(basisRow==2) {{ M{order}=M{order}+leadcoef(z)*B2; }}
  if(basisRow==3) {{ M{order}=M{order}+leadcoef(z)*B3; }}
  if(basisRow==4) {{ M{order}=M{order}+leadcoef(z)*B4; }}
  if(basisRow==5) {{ M{order}=M{order}+leadcoef(z)*B5; }}
  if(basisRow==6) {{ M{order}=M{order}+leadcoef(z)*B6; }}
  z=z-lead(z);
}}
"""
        for order in multiplication_orders
    )
    parameter_orders = multiplication_orders[1:]
    parameter_variables = [f"z{order}" for order in parameter_orders]
    fitting_ring_variables = ",".join(parameter_variables) or "zaux"
    mapped_matrices = "\n".join(
        f"matrix N{order}=imap(fiber,M{order});"
        for order in multiplication_orders
    )
    matrix_expression = "N6" + "".join(
        f"+z{order}*N{order}" for order in parameter_orders
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring fiber=({prime},s1,s2,s3,t1,t2,u),(s6,s5),dp;
option(redSB);
{declarations}
ideal G=std(p4,p5);
print("ALGEBRA "+string(size(G))+" "+string(vdim(G)));
{normal_form_declarations}
{matrix_declarations}
proc coordinateIndex(poly termValue)
{{
  if(termValue==1) {{ return(1); }}
  if(termValue==s6) {{ return(2); }}
  if(termValue==s5) {{ return(3); }}
  if(termValue==s6*s5) {{ return(4); }}
  if(termValue==s5^2) {{ return(5); }}
  if(termValue==s5^3) {{ return(6); }}
  return(0);
}}
poly basisMonomial;
poly z;
int basisColumn;
int basisRow;
matrix B1[6][6];
matrix B2[6][6];
matrix B3[6][6];
matrix B4[6][6];
matrix B5[6][6];
matrix B6[6][6];
for(basisColumn=1;basisColumn<=6;basisColumn++)
{{
  B1[basisColumn,basisColumn]=1;
  if(basisColumn==1) {{ basisMonomial=1; }}
  if(basisColumn==2) {{ basisMonomial=s6; }}
  if(basisColumn==3) {{ basisMonomial=s5; }}
  if(basisColumn==4) {{ basisMonomial=s6*s5; }}
  if(basisColumn==5) {{ basisMonomial=s5^2; }}
  if(basisColumn==6) {{ basisMonomial=s5^3; }}
  z=reduce(s6*basisMonomial,G);
  while(z!=0)
  {{
    basisRow=coordinateIndex(leadmonom(z));
    if(basisRow==0)
    {{
      print("COORDINATE_ERROR s6 "+string(leadmonom(z)));
      exit;
    }}
    B2[basisRow,basisColumn]=leadcoef(z);
    z=z-lead(z);
  }}
  z=reduce(s5*basisMonomial,G);
  while(z!=0)
  {{
    basisRow=coordinateIndex(leadmonom(z));
    if(basisRow==0)
    {{
      print("COORDINATE_ERROR s5 "+string(leadmonom(z)));
      exit;
    }}
    B3[basisRow,basisColumn]=leadcoef(z);
    z=z-lead(z);
  }}
}}
B4=B2*B3;
B5=B3*B3;
B6=B5*B3;
{multiplication_matrix_fills}
ring fitting=({prime},s1,s2,s3,t1,t2,u),(
  {fitting_ring_variables}
),dp;
{mapped_matrices}
poly baseEquation=imap(fiber,p3);
matrix combination[6][6];
combination={matrix_expression};
poly determinantPolynomial=det(combination);
poly determinantCursor=determinantPolynomial;
while(determinantCursor!=0)
{{
  print(
    "FITTING "+string(leadexp(determinantCursor))+" "
    +string(numerator(leadcoef(determinantCursor)))
  );
  determinantCursor=determinantCursor-lead(determinantCursor);
}}
print("BASE "+string(numerator(leadcoef(baseEquation))));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    assert "\n   ? " not in completed.stdout, (
        completed.stdout[-8000:],
        completed.stderr[-2000:],
    )
    algebra = re.search(r"(?m)^ALGEBRA (\d+) (\d+)$", completed.stdout)
    assert algebra is not None and algebra.groups() == ("3", "6"), (
        completed.stdout[-4000:],
        completed.stderr[-2000:],
    )
    fitting_values = [
        {
            "parameter_exponents": [
                int(exponent) for exponent in exponents.split(",")
            ],
            "polynomial": polynomial,
        }
        for exponents, polynomial in re.findall(
            r"(?m)^FITTING ([0-9,]+) (.*)$",
            completed.stdout,
        )
    ]
    assert fitting_values, completed.stdout[-4000:]
    base_marker = re.search(r"(?m)^BASE (.*)$", completed.stdout)
    assert base_marker is not None
    inverse_variable = "finv"
    open_factor = (
        "u"
        "*(3*(s1^2*u-s2)-13*u)"
        "*((99*(s1^2*u-s2)-274*u)^2+30420*(s1*u-t1)^2)"
        "*(351*(s1^2*u-s2)-901*u)"
        "*((99*(s1^2*u-s2)-274*u)"
        "*(351*(s1^2*u-s2)-901*u)+121680*(s1*u-t1)^2)"
    )
    polynomials = [base_marker.group(1)] + [
        datum["polynomial"] for datum in fitting_values
    ] + [f"{inverse_variable}*({open_factor})-1"]
    return {
        "branch": "generic t0-open rank-six common-root Fitting locus",
        "characteristic": prime,
        "ordinary_variables": [
            "s1",
            "s2",
            "s3",
            "t1",
            "t2",
            "u",
            inverse_variable,
        ],
        "polynomials": polynomials,
        "orders": list(orders),
        "quotient_length": 6,
        "fitting_determinant": {
            "parameter_orders": list(parameter_orders),
            "parameter_variables": parameter_variables,
            "determinant_degree": 6,
            "nonzero_coefficient_count": len(fitting_values),
            "maximum_possible_coefficient_count": math.comb(
                len(parameter_variables) + 6,
                6,
            ),
        },
        "open_factor": open_factor,
        "scope": (
            "exact common-root Fitting equations on the generic "
            "Q*J*K*H*u open; exceptional divisors are not covered"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--branch",
        choices=("generic", "L", "Q", "J", "LQ"),
        required=True,
    )
    parser.add_argument("--prime", type=int, default=47)
    parser.add_argument("--timeout", type=int, default=600)
    parser.add_argument(
        "--through",
        choices=(6, 7, 8, 9, 10, 11, 12, 14),
        type=int,
        default=8,
    )
    parser.add_argument("--export-only", action="store_true")
    parser.add_argument(
        "--backend",
        choices=("msolve", "modstd"),
        default="msolve",
    )
    parser.add_argument(
        "--linear-algebra",
        choices=(1, 2, 42, 44),
        type=int,
        default=44,
    )
    parser.add_argument("--fiber-first", action="store_true")
    parser.add_argument("--eliminate", type=int, default=0)
    parser.add_argument("--replay-prime", type=int, default=0)
    parser.add_argument("--include-branch-table", action="store_true")
    parser.add_argument("--deepest-ffnf", action="store_true")
    parser.add_argument(
        "--deepest-ffnf-through",
        choices=(6, 7, 8, 9, 10, 11, 12, 14),
        type=int,
        default=0,
    )
    parser.add_argument("--deepest-t0-zero", action="store_true")
    parser.add_argument("--t0-zero-branch-table", action="store_true")
    parser.add_argument("--t0-open-rank-six", action="store_true")
    parser.add_argument("--t0-open-direct", action="store_true")
    parser.add_argument("--t0-open-reduced", action="store_true")
    parser.add_argument("--t0-open-fitting", action="store_true")
    parser.add_argument(
        "--t0-open-certificate-only",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--output", default="")
    arguments = parser.parse_args()
    orders = tuple(
        order
        for order in CORRECTED_ORDERS
        if order <= arguments.through or (
            arguments.through == 14 and order == 14
        )
    )
    if arguments.prime:
        assert 3 * max(orders) < arguments.prime
    singular = shutil.which("Singular")
    msolve = shutil.which("msolve")
    assert singular is not None and msolve is not None
    if (
        arguments.t0_open_direct
        or arguments.t0_open_reduced
        or arguments.t0_open_fitting
    ):
        export = (
            t0_open_fitting_export(
                singular,
                orders,
                arguments.prime,
                arguments.timeout,
            )
            if arguments.t0_open_fitting
            else (
                t0_open_reduced_export(
                    singular,
                    orders,
                    arguments.prime,
                    arguments.timeout,
                )
                if arguments.t0_open_reduced
                else t0_open_direct_export(orders, arguments.prime)
            )
        )
        result = None
        if not arguments.export_only:
            if arguments.backend == "modstd":
                assert arguments.prime == 0
                result = run_modstd(singular, export, arguments.timeout)
            else:
                result = run_msolve(
                    msolve,
                    export,
                    arguments.prime,
                    arguments.timeout,
                    arguments.linear_algebra,
                    arguments.eliminate,
                )
        summary = {
            key: value for key, value in export.items() if key != "polynomials"
        }
        summary["exported_polynomial_terms"] = [
            polynomial.count("+") + polynomial.count("-") + 1
            for polynomial in export["polynomials"]
        ]
        summary["solve"] = result
        summary["reproduction_command"] = " ".join(sys.argv)
        print(json.dumps(summary, indent=2, sort_keys=True))
        return
    if arguments.t0_open_certificate_only:
        assert arguments.prime == 0
        print(
            json.dumps(
                t0_open_rank_six_certificate(
                    singular,
                    arguments.timeout,
                ),
                sort_keys=True,
            )
        )
        return
    t0_open_result = None
    if arguments.t0_open_rank_six:
        completed = subprocess.run(
            [
                sys.executable,
                str(Path(__file__).resolve()),
                "--branch",
                "LQ",
                "--prime",
                "0",
                "--through",
                "6",
                "--export-only",
                "--t0-open-certificate-only",
                "--timeout",
                str(arguments.timeout),
            ],
            text=True,
            capture_output=True,
            check=True,
            timeout=arguments.timeout,
        )
        t0_open_result = json.loads(completed.stdout)
    polynomials = boundary_polynomials(
        singular,
        orders,
        arguments.prime,
        arguments.timeout,
    )
    export = finite_algebra_export(
        singular,
        polynomials,
        orders,
        arguments.branch,
        arguments.prime,
        arguments.timeout,
    )
    if arguments.fiber_first:
        variables = export["ordinary_variables"]
        fiber_variables = (
            ("t2", "s5", "t4")
            if arguments.branch == "LQ"
            else ("s5", "t4")
        )
        export["ordinary_variables"] = [
            variable for variable in fiber_variables if variable in variables
        ] + [
            variable
            for variable in variables
            if variable not in fiber_variables
        ]
    result = None
    if not arguments.export_only:
        if arguments.backend == "modstd":
            assert arguments.prime == 0
            result = run_modstd(singular, export, arguments.timeout)
        else:
            result = run_msolve(
                msolve,
                export,
                arguments.prime,
                arguments.timeout,
                arguments.linear_algebra,
                arguments.eliminate,
            )
    summary = {key: value for key, value in export.items() if key != "polynomials"}
    summary["corrected_moment_set"] = list(range(1, 13)) + [14]
    summary["constructed_orders"] = list(orders)
    summary["exported_polynomial_terms"] = [
        polynomial.count("+") + polynomial.count("-") + 1
        for polynomial in export["polynomials"]
    ]
    if arguments.branch == "LQ" and arguments.prime == 0:
        summary["deepest_three_pivot"] = deepest_three_pivot_certificate(
            singular,
            polynomials,
            arguments.timeout,
        )
    if arguments.deepest_ffnf:
        assert arguments.branch == "LQ"
        fraction_free_through = (
            arguments.deepest_ffnf_through or arguments.through
        )
        assert fraction_free_through <= arguments.through
        fraction_free_orders = tuple(
            order
            for order in orders
            if order <= fraction_free_through
            or (fraction_free_through == 14 and order == 14)
        )
        fraction_free_polynomials = polynomials[
            : len(fraction_free_orders) - 1
        ]
        try:
            summary["deepest_fraction_free_normal_forms"] = (
                deepest_fraction_free_normal_forms(
                    singular,
                    fraction_free_polynomials,
                    fraction_free_orders,
                    arguments.prime,
                    arguments.timeout,
                )
            )
        except subprocess.TimeoutExpired:
            summary["deepest_fraction_free_normal_forms"] = {
                "status": "timeout",
                "seconds": float(arguments.timeout),
                "requested_later_orders": list(fraction_free_orders[4:]),
                "scope": (
                    "bounded fraction-free reduction attempt; no "
                    "mathematical conclusion"
                ),
            }
    if arguments.deepest_t0_zero:
        assert arguments.branch == "LQ"
        assert arguments.prime == 0
        summary["deepest_t0_zero"] = t0_zero_branch_certificate(
            msolve,
            polynomials,
            orders,
            "LQ",
            arguments.timeout,
        )
    if arguments.t0_zero_branch_table:
        assert arguments.prime == 0
        summary["t0_zero_branch_table"] = {
            branch: t0_zero_branch_certificate(
                msolve,
                polynomials,
                orders,
                branch,
                arguments.timeout,
            )
            for branch in ("generic", "L", "Q", "J", "LQ")
        }
    if arguments.t0_open_rank_six:
        assert arguments.prime == 0
        assert t0_open_result is not None
        summary["t0_open_rank_six"] = t0_open_result
    summary["solve"] = result
    if arguments.include_branch_table:
        assert arguments.prime == 0
        branch_table: dict[str, object] = {}
        for branch in ("generic", "L", "Q", "LQ"):
            branch_export = (
                export
                if branch == arguments.branch
                else finite_algebra_export(
                    singular,
                    polynomials,
                    orders,
                    branch,
                    0,
                    arguments.timeout,
                )
            )
            branch_table[branch] = {
                key: value
                for key, value in branch_export.items()
                if key not in ("polynomials", "ordinary_variables")
            }
            branch_table[branch]["exported_polynomial_terms"] = [
                polynomial.count("+") + polynomial.count("-") + 1
                for polynomial in branch_export["polynomials"]
            ]
        branch_table["J"] = {
            "quotient_length": 5,
            "standard_basis": ["1", "t4", "t4^2", "s5", "s5*t4"],
            "leading_exponents": ["2,0", "1,2", "0,3"],
            "exact_fraction_free_normal_forms": {
                "6": {
                    "support_size": 5,
                    "pseudo_reduction_steps": 5,
                },
                "7": {
                    "support_size": 5,
                    "pseudo_reduction_steps": 10,
                },
            },
            "source": (
                "verify_two_pair_sic_bidegree33_boundary_generic_quotient.py"
            ),
            "later_corrected_moments_reduced": [6, 7],
        }
        summary["branch_table"] = branch_table
    if arguments.replay_prime:
        assert arguments.prime == 0
        assert 3 * max(orders) < arguments.replay_prime
        replay_polynomials = boundary_polynomials(
            singular,
            orders,
            arguments.replay_prime,
            arguments.timeout,
        )
        replay_export = finite_algebra_export(
            singular,
            replay_polynomials,
            orders,
            arguments.branch,
            arguments.replay_prime,
            arguments.timeout,
        )
        summary["modular_replay"] = {
            "prime": arguments.replay_prime,
            "solve": run_msolve(
                msolve,
                replay_export,
                arguments.replay_prime,
                arguments.timeout,
                arguments.linear_algebra,
                arguments.eliminate,
            ),
            "scope": "finite-field reconnaissance only",
        }
    summary["scope"] = (
        "exact finite-algebra branch exports through the constructed orders; "
        "no semistable point and no characteristic-zero radical equality"
    )
    summary["reproduction_command"] = " ".join(sys.argv)
    if arguments.output:
        output = Path(arguments.output)
        if not output.is_absolute():
            output = ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
