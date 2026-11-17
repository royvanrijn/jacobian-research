#!/usr/bin/env python3
"""Modular coefficient-ideal test on the reduced SIC(2) (3,3) boundary.

On the s0=1, A=B=0 chart, the exact mu_2, A, and B pivots leave six
base variables and the two fiber variables (s5,t4).  The certified generic
quotient by (mu_4,mu_5) has rank six.  This script reduces the later moments
in that quotient, extracts their fiber coefficients, and computes the
resulting base ideal on the L*Q*J principal open.

The output is modular evidence only.  In particular, a unit ideal must be
reconstructed over characteristic zero before it can be promoted to a
nullcone certificate.

The coefficient-ideal modes test whether a later moment vanishes as an
element of the whole finite quotient; that is stronger than vanishing at
one quotient point and is not a common-root test.  The ``--trace-norm``
mode keeps mu_3, mu_4, and mu_5 in a rank-twelve finite algebra after
L=1 and is the appropriate starting point for trace, norm, characteristic
polynomial, and Fitting computations for a common mu_6,mu_7 zero.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from explore_two_pair_sic_bidegree33_full_anchor import (  # noqa: E402
    chart_expression,
    moment_terms,
    prepare_s0_branch_for_msolve,
)
from verify_two_pair_sic_bidegree33_boundary_generic_quotient import (  # noqa: E402
    T2_NUMERATOR_ADAPTED,
    exact_chart_expression,
    exact_moment_terms,
    substitute,
)


BASE_VARIABLES = ("s1", "s2", "s3", "t0", "t1", "t2")
LINEAR = "s1*t0-t1"
QUADRATIC = "s1^2-s2-(13/3)*t0^2"
J_DIVISOR = (
    "9801*s1^4-19602*s1^2*s2-23832*s1^2*t0^2"
    "-60840*s1*t0*t1+9801*s2^2+54252*s2*t0^2"
    "+75076*t0^4+30420*t1^2"
)


def quotient_coefficients(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    polynomials: list[str],
    timeout: int,
) -> dict[int, list[str]]:
    names = [f"p{order}" for order in orders]
    declarations = "\n".join(
        f"poly {name}={polynomial};"
        for name, polynomial in zip(names, polynomials)
    )
    reductions = []
    for order in orders:
        if order in (3, 4, 5):
            continue
        reductions.append(
            f"""
poly r{order}=reduce(p{order},G);
poly z{order}=r{order};
number c{order};
while (z{order}!=0)
{{
  c{order}=leadcoef(z{order});
  print(
    "COEFFICIENT {order} "
    +string(leadexp(z{order}))
    +" "+string(numerator(c{order}))
  );
  z{order}=z{order}-c{order}*leadmonom(z{order});
}}
"""
        )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring fiber=({prime},{",".join(BASE_VARIABLES)}),(s5,t4),dp;
{declarations}
ideal G=std(p4,p5);
print("QUOTIENT "+string(vdim(G))+" "+string(size(G)));
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("LEADING "+string(leadexp(G[basisIndex])));
}}
{"".join(reductions)}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    quotient = re.search(r"(?m)^QUOTIENT (\d+) (\d+)$", completed.stdout)
    if quotient is None:
        raise AssertionError(completed.stdout[-4000:])
    print(
        f"QUOTIENT length={quotient.group(1)} basis={quotient.group(2)}",
        flush=True,
    )
    leading = re.findall(r"(?m)^LEADING ([0-9,]+)$", completed.stdout)
    print(f"FIBER_INITIAL {','.join(leading)}", flush=True)
    coefficients: dict[int, list[str]] = {order: [] for order in orders}
    for order, _exponent, coefficient in re.findall(
        r"(?m)^COEFFICIENT (\d+) ([0-9,]+) (.*)$",
        completed.stdout,
    ):
        coefficients[int(order)].append(coefficient)
    for order in orders:
        if order not in (3, 4, 5):
            print(
                f"REMAINDER order={order} "
                f"coefficients={len(coefficients[order])}",
                flush=True,
            )
    return coefficients


def base_ideal(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    base_polynomial: str,
    coefficients: dict[int, list[str]],
    timeout: int,
) -> None:
    generators = [base_polynomial]
    for order in orders:
        if order in (3, 4, 5):
            continue
        generators.extend(coefficients[order])
    generator_text = ",".join(generators)
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring base={prime},({",".join(BASE_VARIABLES)},u),dp;
poly L={LINEAR};
poly Q={QUADRATIC};
poly J={J_DIVISOR};
ideal I={generator_text},u*L*Q*J-1;
ideal G=std(I);
print(
  "BASE "+string(dim(G))+" "+string(size(G))+" "+string(G[1]==1)
);
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("BASE_LEADING "+string(leadexp(G[basisIndex])));
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    marker = re.search(r"(?m)^BASE (-?\d+) (\d+) ([01])$", completed.stdout)
    if marker is None:
        raise AssertionError(completed.stdout[-4000:])
    dimension, size, unit = marker.groups()
    print(
        f"BASE dimension={dimension} basis={size} unit={unit}",
        flush=True,
    )
    leading = re.findall(r"(?m)^BASE_LEADING ([0-9,]+)$", completed.stdout)
    print(f"BASE_INITIAL {','.join(leading)}", flush=True)


def adapted_coefficients(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    polynomials: list[str],
    timeout: int,
) -> tuple[str, dict[int, list[str]]]:
    """Use the known t2 pivot and return coefficients in adapted variables."""

    replacements = (
        ("t1", "(s1*t0-L)"),
        ("s2", "(s1^2-(13/3)*t0^2-Q)"),
        ("t2", "tt"),
    )
    adapted = [
        substitute(polynomial, replacements)
        for polynomial in polynomials
    ]
    names = [f"p{order}" for order in orders]
    declarations = "\n".join(
        f"poly {name}={polynomial};"
        for name, polynomial in zip(names, adapted)
    )
    reductions = []
    for order in orders:
        if order in (3, 4, 5):
            continue
        reductions.append(
            f"""
poly r{order}=reduce(p{order},G);
poly z{order}=r{order};
number c{order};
while (z{order}!=0)
{{
  c{order}=leadcoef(z{order});
  print(
    "ADAPTED_COEFFICIENT {order} "
    +string(leadexp(z{order}))
    +" "+string(numerator(c{order}))
  );
  z{order}=z{order}-c{order}*leadmonom(z{order});
}}
"""
        )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring fiber=({prime},s1,s3,t0,L,Q),(s5,t4),dp;
number tt=({T2_NUMERATOR_ADAPTED})/(93366*L*Q);
{declarations}
ideal G=std(p4,p5);
print("ADAPTED_P3 "+string(numerator(leadcoef(p3))));
print("ADAPTED_QUOTIENT "+string(vdim(G))+" "+string(size(G)));
{"".join(reductions)}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    p3_marker = re.search(r"(?m)^ADAPTED_P3 \((.*)\)$", completed.stdout)
    quotient = re.search(
        r"(?m)^ADAPTED_QUOTIENT (\d+) (\d+)$",
        completed.stdout,
    )
    if p3_marker is None or quotient is None:
        raise AssertionError(completed.stdout[-4000:])
    print(
        f"ADAPTED_QUOTIENT length={quotient.group(1)} "
        f"basis={quotient.group(2)}",
        flush=True,
    )
    coefficients: dict[int, list[str]] = {order: [] for order in orders}
    for order, _exponent, coefficient in re.findall(
        r"(?m)^ADAPTED_COEFFICIENT (\d+) ([0-9,]+) (.*)$",
        completed.stdout,
    ):
        coefficients[int(order)].append(coefficient)
    for order in orders:
        if order not in (3, 4, 5):
            print(
                f"ADAPTED_REMAINDER order={order} "
                f"coefficients={len(coefficients[order])}",
                flush=True,
            )
    metadata_declarations = [f"poly c3={p3_marker.group(1)};"]
    metadata_prints = [
        'print("ADAPTED_COEFFICIENT_META 3 0 "+string(deg(c3))+" "+string(size(c3)));'
    ]
    for order in orders:
        for index, coefficient in enumerate(coefficients[order]):
            name = f"c{order}_{index}"
            metadata_declarations.append(f"poly {name}={coefficient};")
            metadata_prints.append(
                f'print("ADAPTED_COEFFICIENT_META {order} {index} "'
                f'+string(deg({name}))+" "+string(size({name})));'
            )
    metadata = subprocess.run(
        [singular, "-q"],
        input=f"""
ring metadata={prime},(s1,s3,t0,L,Q),dp;
{"".join(metadata_declarations)}
{"".join(metadata_prints)}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    for order, index, degree, terms in re.findall(
        r"(?m)^ADAPTED_COEFFICIENT_META (\d+) (\d+) (\d+) (\d+)$",
        metadata.stdout,
    ):
        print(
            f"ADAPTED_COEFFICIENT_META order={order} index={index} "
            f"degree={degree} terms={terms}",
            flush=True,
        )
    return p3_marker.group(1), coefficients


def adapted_base_ideal(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    base_polynomial: str,
    coefficients: dict[int, list[str]],
    timeout: int,
) -> None:
    generators = [base_polynomial]
    for order in orders:
        if order in (3, 4, 5):
            continue
        generators.extend(coefficients[order])
    adapted_j = "(99*Q+155*t0^2)^2+30420*L^2"
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring base={prime},(s1,s3,t0,L,Q,u),dp;
ideal I={",".join(generators)},u*L*Q*({adapted_j})-1;
ideal G=std(I);
print(
  "ADAPTED_BASE "+string(dim(G))+" "+string(size(G))+" "
  +string(G[1]==1)
);
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("ADAPTED_BASE_LEADING "+string(leadexp(G[basisIndex])));
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    marker = re.search(
        r"(?m)^ADAPTED_BASE (-?\d+) (\d+) ([01])$",
        completed.stdout,
    )
    if marker is None:
        raise AssertionError(completed.stdout[-4000:])
    dimension, size, unit = marker.groups()
    print(
        f"ADAPTED_BASE dimension={dimension} basis={size} unit={unit}",
        flush=True,
    )
    leading = re.findall(
        r"(?m)^ADAPTED_BASE_LEADING ([0-9,]+)$",
        completed.stdout,
    )
    print(f"ADAPTED_BASE_INITIAL {','.join(leading)}", flush=True)


def adapted_base_msolve(
    msolve: str,
    prime: int,
    orders: tuple[int, ...],
    base_polynomial: str,
    coefficients: dict[int, list[str]],
    timeout: int,
    threads: int,
    linear_algebra: int,
    coefficient_counts: dict[int, int],
    normalize_linear: bool,
) -> None:
    generators = [base_polynomial]
    for order in orders:
        if order in (3, 4, 5):
            continue
        generators.extend(
            coefficients[order][: coefficient_counts.get(order, len(coefficients[order]))]
        )
    saturation = "L*Q*((99*Q+155*t0^2)^2+30420*L^2)"
    variables = "s1,s3,t0,L,Q"
    if normalize_linear:
        replacements = (("L", "(1)"),)
        generators = [
            substitute(generator, replacements)
            for generator in generators
        ]
        saturation = substitute(saturation, replacements)
        variables = "s1,s3,t0,Q"
    msolve_options = ["-S"]
    if prime == 0:
        variables += ",u"
        generators.append(f"u*({saturation})-1")
        msolve_options = []
    else:
        generators.append(saturation)
    with tempfile.TemporaryDirectory(prefix="sic33-boundary-base-") as directory:
        input_path = Path(directory) / "base.ms"
        output_path = Path(directory) / "base.out"
        input_path.write_text(
            variables
            + "\n"
            f"{prime}\n"
            + ",\n".join(generators)
            + "\n"
        )
        command = [
            msolve,
            "-f",
            str(input_path),
            "-o",
            str(output_path),
            "-t",
            str(threads),
            "-v",
            "1",
            "-l",
            str(linear_algebra),
            "-g",
            "1",
        ] + msolve_options
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
            timeout=timeout,
        )
        if completed.returncode != 0:
            raise AssertionError(
                f"msolve exit {completed.returncode}\n{completed.stdout[-2000:]}\n"
                f"{completed.stderr[-2000:]}"
            )
        result = output_path.read_text().strip()
    if result.rstrip().endswith("[1]:") and "#length of basis:      1 element" in result:
        status = "unit"
    elif result in ("[-1]", "[-1]:"):
        status = "unit"
    elif result.startswith("[1,") and ",-1,[]" in result.replace(" ", ""):
        status = "positive-dimensional"
    else:
        status = "finite-nonempty"
    leading_lines = [
        line
        for line in completed.stdout.splitlines()
        if "leading" in line.lower() or line.startswith("[")
    ]
    print(f"ADAPTED_MSOLVE status={status}", flush=True)
    if status != "unit":
        print(
            "ADAPTED_MSOLVE_RESULT "
            + result.replace("\n", " ")[:1000],
            flush=True,
        )
    if leading_lines:
        print(
            "ADAPTED_MSOLVE_INITIAL "
            + " | ".join(leading_lines[-10:]),
            flush=True,
        )


def generic_common_root(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    polynomials: list[str],
    timeout: int,
    lift: bool,
    normalize_linear: bool,
) -> None:
    """Keep one fiber point and extract the generic resultant denominator."""

    available = dict(zip(orders, polynomials))
    required = (3, 4, 5, 6)
    assert all(order in available for order in required)
    ring_specification = (
        f"ring common=({prime},s1,s2,s3,t0,t1,t2),(s5,t4),dp;"
    )
    minpoly_program = ""
    if normalize_linear:
        raise ValueError(
            "Singular does not support this multivariate algebraic "
            "coefficient-field shortcut; use --trace-norm instead"
        )
    declarations = "\n".join(
        f"poly p{order}={available[order]};"
        for order in (4, 5, 6)
    )
    basis_program = """
matrix transformation;
ideal I456=p4,p5,p6;
ideal G456=liftstd(I456,transformation,"slimgb");
"""
    denominator_program = """
int rowIndex;
int columnIndex;
poly z;
number c;
for (rowIndex=1;rowIndex<=nrows(transformation);rowIndex++)
{
  for (columnIndex=1;columnIndex<=ncols(transformation);columnIndex++)
  {
    z=transformation[rowIndex,columnIndex];
    while (z!=0)
    {
      c=leadcoef(z);
      print(
        "COMMON_DENOMINATOR "+string(rowIndex)+" "
        +string(columnIndex)+" "+string(denominator(c))
      );
      z=z-c*leadmonom(z);
    }
  }
}
"""
    if not lift:
        basis_program = """
ideal I456=p4,p5,p6;
ideal G456=slimgb(I456);
"""
        denominator_program = ""
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
{ring_specification}
{minpoly_program}
{declarations}
ideal I45=p4,p5;
ideal G45=std(I45);
{basis_program}
print(
  "COMMON_45 "+string(dim(G45))+" "+string(vdim(G45))+" "
  +string(size(G45))
);
print(
  "COMMON_456 "+string(dim(G456))+" "+string(size(G456))+" "
  +string(G456[1]==1)
);
{denominator_program}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    first = re.search(
        r"(?m)^COMMON_45 (-?\d+) (-?\d+) (\d+)$",
        completed.stdout,
    )
    second = re.search(
        r"(?m)^COMMON_456 (-?\d+) (\d+) ([01])$",
        completed.stdout,
    )
    if first is None or second is None:
        raise AssertionError(completed.stdout[-4000:])
    print(
        f"COMMON_45 dimension={first.group(1)} "
        f"length={first.group(2)} basis={first.group(3)}",
        flush=True,
    )
    print(
        f"COMMON_456 dimension={second.group(1)} "
        f"basis={second.group(2)} unit={second.group(3)}",
        flush=True,
    )
    denominators = re.findall(
        r"(?m)^COMMON_DENOMINATOR \d+ \d+ (.*)$",
        completed.stdout,
    )
    unique_denominators = list(dict.fromkeys(denominators))
    print(
        f"COMMON_DENOMINATORS total={len(denominators)} "
        f"unique={len(unique_denominators)}",
        flush=True,
    )
    for index, denominator in enumerate(unique_denominators):
        print(
            f"COMMON_DENOMINATOR_UNIQUE index={index} value={denominator}",
            flush=True,
        )


def multiplication_fitting(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    polynomials: list[str],
    timeout: int,
    normalize_linear: bool,
    norm_only: bool,
    matrix_only: bool,
) -> None:
    """Construct selected Fitting minors in the rank-six fiber quotient."""

    available = dict(zip(orders, polynomials))
    assert all(order in available for order in (4, 5, 6))
    p7_expression = available.get(7, "0") if not norm_only else "0"
    maximum_replacement = -1 if matrix_only else (0 if norm_only else 6)
    matrix_entry_program = ""
    if matrix_only:
        matrix_entry_program = """
for (rowIndex=1;rowIndex<=6;rowIndex++)
{
  for (columnIndex=1;columnIndex<=6;columnIndex++)
  {
    print(
      "FITTING_MATRIX_ENTRY "+string(rowIndex)+" "+string(columnIndex)+" "
      +string(numerator(leadcoef(M6[rowIndex,columnIndex])))+" @ "
      +string(denominator(leadcoef(M6[rowIndex,columnIndex])))
    );
  }
}
"""
    ring_specification = (
        f"ring fit=({prime},s1,s2,s3,t0,t1,t2),(s5,t4),dp;"
    )
    minpoly_program = ""
    if normalize_linear:
        raise ValueError(
            "use --trace-norm for the normalized common-root quotient"
        )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
{ring_specification}
{minpoly_program}
poly p4={available[4]};
poly p5={available[5]};
poly p6={available[6]};
poly p7={p7_expression};
ideal G=std(p4,p5);
poly r6=reduce(p6,G);
poly r7=reduce(p7,G);

proc coordinateIndex(poly m)
{{
  if (m==1) {{ return(1); }}
  if (m==t4) {{ return(2); }}
  if (m==t4^2) {{ return(3); }}
  if (m==t4^3) {{ return(4); }}
  if (m==s5) {{ return(5); }}
  if (m==s5*t4) {{ return(6); }}
  return(0);
}}

poly basisMonomial;
poly z;
poly m;
number c;
int rowIndex;
int columnIndex;
matrix M6[6][6];
matrix M7[6][6];
for (columnIndex=1;columnIndex<=6;columnIndex++)
{{
  if (columnIndex==1) {{ basisMonomial=1; }}
  if (columnIndex==2) {{ basisMonomial=t4; }}
  if (columnIndex==3) {{ basisMonomial=t4^2; }}
  if (columnIndex==4) {{ basisMonomial=t4^3; }}
  if (columnIndex==5) {{ basisMonomial=s5; }}
  if (columnIndex==6) {{ basisMonomial=s5*t4; }}
  z=reduce(r6*basisMonomial,G);
  while (z!=0)
  {{
    m=leadmonom(z); c=leadcoef(z); rowIndex=coordinateIndex(m);
    if (rowIndex==0) {{ print("COORDINATE_ERROR_6 "+string(m)); exit; }}
    M6[rowIndex,columnIndex]=c;
    z=z-c*m;
  }}
  z=reduce(r7*basisMonomial,G);
  while (z!=0)
  {{
    m=leadmonom(z); c=leadcoef(z); rowIndex=coordinateIndex(m);
    if (rowIndex==0) {{ print("COORDINATE_ERROR_7 "+string(m)); exit; }}
    M7[rowIndex,columnIndex]=c;
    z=z-c*m;
  }}
}}

matrix S[6][6];
number minorValue;
int replacement;
{matrix_entry_program}
for (replacement=0;replacement<={maximum_replacement};replacement++)
{{
  S=M6;
  if (replacement>0)
  {{
    for (rowIndex=1;rowIndex<=6;rowIndex++)
    {{
      S[rowIndex,replacement]=M7[rowIndex,replacement];
    }}
  }}
  minorValue=det(S);
  print(
    "FITTING_MINOR "+string(replacement)+" "
    +string(numerator(minorValue))
  );
}}
print("FITTING_QUOTIENT "+string(vdim(G))+" "+string(size(G)));
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    quotient = re.search(
        r"(?m)^FITTING_QUOTIENT (\d+) (\d+)$",
        completed.stdout,
    )
    if quotient is None:
        raise AssertionError(completed.stdout[-4000:])
    print(
        f"FITTING_QUOTIENT length={quotient.group(1)} "
        f"basis={quotient.group(2)}",
        flush=True,
    )
    minors = re.findall(r"(?m)^FITTING_MINOR (\d+) (.*)$", completed.stdout)
    assert len(minors) == maximum_replacement + 1
    for index, polynomial in minors:
        print(
            f"FITTING_MINOR_META index={index} "
            f"sha256={hashlib.sha256(polynomial.encode()).hexdigest()} "
            f"characters={len(polynomial)}",
            flush=True,
        )
    if matrix_only:
        entries = re.findall(
            r"(?m)^FITTING_MATRIX_ENTRY (\d+) (\d+) (.*) @ (.*)$",
            completed.stdout,
        )
        if len(entries) != 36:
            raise AssertionError(completed.stdout[-8000:])
        output = (
            ROOT
            / "artifacts"
            / "generated-results"
            / f"two_pair_sic_bidegree33_boundary_m6_matrix_mod{prime}.json"
        )
        output.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "prime": prime,
            "normalization": "s0=1, A=B=mu2=0, L=1, mu3(t2)=0",
            "basis": ["1", "t4", "t4^2", "t4^3", "s5", "s5*t4"],
            "entries": [
                {
                    "row": int(row),
                    "column": int(column),
                    "numerator": numerator,
                    "denominator": denominator,
                }
                for row, column, numerator, denominator in entries
            ],
            "scope": (
                "exact modular multiplication matrix for mu6 in the "
                "rank-six boundary quotient; determinant not yet computed"
            ),
        }
        output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"FITTING_MATRIX_WROTE {output.relative_to(ROOT)}", flush=True)


def trace_norm_export(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    polynomials: list[str],
    timeout: int,
) -> None:
    """Export rank-twelve multiplication matrices after mu3,mu4,mu5."""

    available = dict(zip(orders, polynomials))
    assert all(order in available for order in (3, 4, 5, 6))
    p7_expression = available.get(7, "0")
    maximum_matrix_order = 7 if 7 in available else 6
    replacements = (
        ("t1", "(s1*t0-L)"),
        ("s2", "(s1^2-(13/3)*t0^2-Q)"),
        ("L", "(1)"),
    )
    available = {
        order: substitute(polynomial, replacements)
        for order, polynomial in available.items()
    }
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring tn=({prime},s1,t0,Q,t2),(s5,t4,s3),dp;
poly p3={available[3]};
poly p4={available[4]};
poly p5={available[5]};
poly p6={available[6]};
poly p7={p7_expression};
ideal I345=p3,p4,p5;
ideal G=slimgb(I345);
poly r6=reduce(p6,G);
poly r7=reduce(p7,G);

proc coordinateIndex(poly m)
{{
  if (m==1) {{ return(1); }}
  if (m==t4) {{ return(2); }}
  if (m==t4^2) {{ return(3); }}
  if (m==t4^3) {{ return(4); }}
  if (m==s5) {{ return(5); }}
  if (m==s5*t4) {{ return(6); }}
  if (m==s3) {{ return(7); }}
  if (m==s3*t4) {{ return(8); }}
  if (m==s3*t4^2) {{ return(9); }}
  if (m==s3*t4^3) {{ return(10); }}
  if (m==s3*s5) {{ return(11); }}
  if (m==s3*s5*t4) {{ return(12); }}
  return(0);
}}

poly basisMonomial;
poly z;
poly m;
number c;
int rowIndex;
int columnIndex;
matrix MT4[12][12];
matrix MS5[12][12];
matrix MS3[12][12];
matrix M6[12][12];
matrix M7[12][12];
for (columnIndex=1;columnIndex<=12;columnIndex++)
{{
  if (columnIndex==1) {{ basisMonomial=1; }}
  if (columnIndex==2) {{ basisMonomial=t4; }}
  if (columnIndex==3) {{ basisMonomial=t4^2; }}
  if (columnIndex==4) {{ basisMonomial=t4^3; }}
  if (columnIndex==5) {{ basisMonomial=s5; }}
  if (columnIndex==6) {{ basisMonomial=s5*t4; }}
  if (columnIndex==7) {{ basisMonomial=s3; }}
  if (columnIndex==8) {{ basisMonomial=s3*t4; }}
  if (columnIndex==9) {{ basisMonomial=s3*t4^2; }}
  if (columnIndex==10) {{ basisMonomial=s3*t4^3; }}
  if (columnIndex==11) {{ basisMonomial=s3*s5; }}
  if (columnIndex==12) {{ basisMonomial=s3*s5*t4; }}
  z=reduce(t4*basisMonomial,G);
  while (z!=0)
  {{
    m=leadmonom(z); c=leadcoef(z); rowIndex=coordinateIndex(m);
    if (rowIndex==0) {{ print("TRACE_NORM_COORDINATE_ERROR_T4 "+string(m)); exit; }}
    MT4[rowIndex,columnIndex]=c;
    z=z-c*m;
  }}
  z=reduce(s5*basisMonomial,G);
  while (z!=0)
  {{
    m=leadmonom(z); c=leadcoef(z); rowIndex=coordinateIndex(m);
    if (rowIndex==0) {{ print("TRACE_NORM_COORDINATE_ERROR_S5 "+string(m)); exit; }}
    MS5[rowIndex,columnIndex]=c;
    z=z-c*m;
  }}
  z=reduce(s3*basisMonomial,G);
  while (z!=0)
  {{
    m=leadmonom(z); c=leadcoef(z); rowIndex=coordinateIndex(m);
    if (rowIndex==0) {{ print("TRACE_NORM_COORDINATE_ERROR_S3 "+string(m)); exit; }}
    MS3[rowIndex,columnIndex]=c;
    z=z-c*m;
  }}
}}

matrix coefficient6[1][12];
matrix coefficient7[1][12];
z=r6;
while (z!=0)
{{
  m=leadmonom(z); c=leadcoef(z); rowIndex=coordinateIndex(m);
  coefficient6[1,rowIndex]=c;
  z=z-c*m;
}}
z=r7;
while (z!=0)
{{
  m=leadmonom(z); c=leadcoef(z); rowIndex=coordinateIndex(m);
  coefficient7[1,rowIndex]=c;
  z=z-c*m;
}}

print("TRACE_NORM_QUOTIENT "+string(vdim(G))+" "+string(size(G)));
for (rowIndex=1;rowIndex<=size(G);rowIndex++)
{{
  print("TRACE_NORM_LEADING "+string(leadexp(G[rowIndex])));
}}
for (columnIndex=1;columnIndex<=12;columnIndex++)
{{
  print(
    "TRACE_NORM_COEFFICIENT 6 "+string(columnIndex)+" "
    +string(numerator(leadcoef(coefficient6[1,columnIndex])))+" @ "
    +string(denominator(leadcoef(coefficient6[1,columnIndex])))
  );
  print(
    "TRACE_NORM_COEFFICIENT 7 "+string(columnIndex)+" "
    +string(numerator(leadcoef(coefficient7[1,columnIndex])))+" @ "
    +string(denominator(leadcoef(coefficient7[1,columnIndex])))
  );
}}
for (rowIndex=1;rowIndex<=12;rowIndex++)
{{
  for (columnIndex=1;columnIndex<=12;columnIndex++)
  {{
    print(
      "TRACE_NORM_GENERATOR t4 "+string(rowIndex)+" "+string(columnIndex)+" "
      +string(numerator(leadcoef(MT4[rowIndex,columnIndex])))+" @ "
      +string(denominator(leadcoef(MT4[rowIndex,columnIndex])))
    );
    print(
      "TRACE_NORM_GENERATOR s5 "+string(rowIndex)+" "+string(columnIndex)+" "
      +string(numerator(leadcoef(MS5[rowIndex,columnIndex])))+" @ "
      +string(denominator(leadcoef(MS5[rowIndex,columnIndex])))
    );
    print(
      "TRACE_NORM_GENERATOR s3 "+string(rowIndex)+" "+string(columnIndex)+" "
      +string(numerator(leadcoef(MS3[rowIndex,columnIndex])))+" @ "
      +string(denominator(leadcoef(MS3[rowIndex,columnIndex])))
    );
  }}
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    quotient = re.search(
        r"(?m)^TRACE_NORM_QUOTIENT (\d+) (\d+)$",
        completed.stdout,
    )
    if quotient is None:
        raise AssertionError(completed.stdout[-8000:])
    leading = re.findall(r"(?m)^TRACE_NORM_LEADING ([0-9,]+)$", completed.stdout)
    coefficients = re.findall(
        r"(?m)^TRACE_NORM_COEFFICIENT ([67]) (\d+) (.*) @ (.*)$",
        completed.stdout,
    )
    generators = re.findall(
        r"(?m)^TRACE_NORM_GENERATOR (t4|s5|s3) "
        r"(\d+) (\d+) (.*) @ (.*)$",
        completed.stdout,
    )
    if len(coefficients) != 24 or len(generators) != 432:
        raise AssertionError(completed.stdout[-8000:])
    output = (
        ROOT
        / "artifacts"
        / "generated-results"
        / (
            "two_pair_sic_bidegree33_boundary_finite_algebra_"
            f"through_mu{maximum_matrix_order}_mod{prime}.json"
        )
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "prime": prime,
        "normalization": "s0=1, A=B=mu2=0, L=1",
        "quotient_generators": ["mu3", "mu4", "mu5"],
        "quotient_length": int(quotient.group(1)),
        "groebner_basis_size": int(quotient.group(2)),
        "leading_exponents_s5_t4_s3": [
            [int(value) for value in exponent.split(",")]
            for exponent in leading
        ],
        "basis": [
            "1", "t4", "t4^2", "t4^3", "s5", "s5*t4",
            "s3", "s3*t4", "s3*t4^2", "s3*t4^3", "s3*s5", "s3*s5*t4",
        ],
        "normal_form_coefficients": {
            order: [
                {
                    "basis_index": int(index),
                    "numerator": numerator,
                    "denominator": denominator,
                }
                for coefficient_order, index, numerator, denominator in coefficients
                if coefficient_order == order
            ]
            for order in ("6", "7")
        },
        "generator_matrices": {
            generator: [
                {
                    "row": int(row),
                    "column": int(column),
                    "numerator": numerator,
                    "denominator": denominator,
                }
                for entry_generator, row, column, numerator, denominator in generators
                if entry_generator == generator
            ]
            for generator in ("t4", "s5", "s3")
        },
        "scope": (
            "exact modular generator matrices and moment normal forms in the "
            "rank-twelve "
            f"common-root quotient through mu{maximum_matrix_order}; "
            "assemble multiplication matrices and compute trace, norm, "
            "characteristic polynomials, and Fitting minors externally"
        ),
    }
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"TRACE_NORM_QUOTIENT length={quotient.group(1)} "
        f"basis={quotient.group(2)}",
        flush=True,
    )
    print(f"TRACE_NORM_WROTE {output.relative_to(ROOT)}", flush=True)


def trace_norm_samples(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    polynomials: list[str],
    timeout: int,
    sample_count: int,
) -> None:
    """Compute exact multiplication invariants at deterministic base points."""

    available = dict(zip(orders, polynomials))
    assert all(order in available for order in (3, 4, 5, 6, 7))
    adapted = {
        order: substitute(
            polynomial,
            (
                ("t1", "(s1*t0-L)"),
                ("s2", "(s1^2-(13/3)*t0^2-Q)"),
                ("L", "(1)"),
            ),
        )
        for order, polynomial in available.items()
    }
    samples: list[dict[str, object]] = []
    candidate_index = 0
    while len(samples) < sample_count:
        block, offset = divmod(candidate_index, prime - 1)
        point = {
            "s1": 1 + (2 + 3 * offset + 5 * block + 7 * block**2) % (prime - 1),
            "t0": 1 + (3 + 5 * offset + 7 * block + 11 * block**2) % (prime - 1),
            "Q": 1 + (5 + 7 * offset + 11 * block + 13 * block**2) % (prime - 1),
            "t2": 1 + (
                7 + 11 * offset + 13 * block + 17 * block**2
            ) % (prime - 1),
        }
        candidate_index += 1
        base_declarations = "\n".join(
            (
                f"poly a{order}={polynomial};\n"
                f"poly q{order}=subst(subst(subst(subst("
                f"a{order},s1,{point['s1']}),t0,{point['t0']}),"
                f"Q,{point['Q']}),t2,{point['t2']});"
            )
            for order, polynomial in adapted.items()
        )
        fiber_declarations = "\n".join(
            f"poly p{order}=imap(base,q{order});"
            for order in adapted
        )
        completed = subprocess.run(
            [singular, "-q"],
            input=f"""
ring base={prime},(s1,t0,Q,t2,s5,t4,s3),dp;
{base_declarations}
ring tn={prime},(s5,t4,s3),dp;
{fiber_declarations}
ideal I345=p3,p4,p5;
ideal G=slimgb(I345);
poly remainder6=reduce(p6,G);
poly remainder7=reduce(p7,G);
print("SAMPLE_QUOTIENT "+string(vdim(G))+" "+string(size(G)));
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("SAMPLE_LEADING "+string(leadexp(G[basisIndex])));
}}

proc coordinateIndex(poly m)
{{
  if (m==1) {{ return(1); }}
  if (m==t4) {{ return(2); }}
  if (m==t4^2) {{ return(3); }}
  if (m==t4^3) {{ return(4); }}
  if (m==s5) {{ return(5); }}
  if (m==s5*t4) {{ return(6); }}
  if (m==s3) {{ return(7); }}
  if (m==s3*t4) {{ return(8); }}
  if (m==s3*t4^2) {{ return(9); }}
  if (m==s3*t4^3) {{ return(10); }}
  if (m==s3*s5) {{ return(11); }}
  if (m==s3*s5*t4) {{ return(12); }}
  return(0);
}}

poly basisMonomial;
poly z;
poly m;
number c;
poly trace6=0;
poly trace7=0;
int rowIndex;
int columnIndex;
matrix M6[12][12];
matrix M7[12][12];
for (columnIndex=1;columnIndex<=12;columnIndex++)
{{
  if (columnIndex==1) {{ basisMonomial=1; }}
  if (columnIndex==2) {{ basisMonomial=t4; }}
  if (columnIndex==3) {{ basisMonomial=t4^2; }}
  if (columnIndex==4) {{ basisMonomial=t4^3; }}
  if (columnIndex==5) {{ basisMonomial=s5; }}
  if (columnIndex==6) {{ basisMonomial=s5*t4; }}
  if (columnIndex==7) {{ basisMonomial=s3; }}
  if (columnIndex==8) {{ basisMonomial=s3*t4; }}
  if (columnIndex==9) {{ basisMonomial=s3*t4^2; }}
  if (columnIndex==10) {{ basisMonomial=s3*t4^3; }}
  if (columnIndex==11) {{ basisMonomial=s3*s5; }}
  if (columnIndex==12) {{ basisMonomial=s3*s5*t4; }}
  z=reduce(remainder6*basisMonomial,G);
  while (z!=0)
  {{
    m=leadmonom(z); c=leadcoef(z); rowIndex=coordinateIndex(m);
    if (rowIndex==0) {{ print("SAMPLE_COORDINATE_ERROR_6 "+string(m)); exit; }}
    M6[rowIndex,columnIndex]=c;
    z=z-c*m;
  }}
  z=reduce(remainder7*basisMonomial,G);
  while (z!=0)
  {{
    m=leadmonom(z); c=leadcoef(z); rowIndex=coordinateIndex(m);
    if (rowIndex==0) {{ print("SAMPLE_COORDINATE_ERROR_7 "+string(m)); exit; }}
    M7[rowIndex,columnIndex]=c;
    z=z-c*m;
  }}
  trace6=trace6+M6[columnIndex,columnIndex];
  trace7=trace7+M7[columnIndex,columnIndex];
}}
matrix joint[12][24];
for (rowIndex=1;rowIndex<=12;rowIndex++)
{{
  for (columnIndex=1;columnIndex<=12;columnIndex++)
  {{
    joint[rowIndex,columnIndex]=M6[rowIndex,columnIndex];
    joint[rowIndex,columnIndex+12]=M7[rowIndex,columnIndex];
  }}
}}
ideal Ireduced=I345,remainder6,remainder7;
ideal Greduced=std(Ireduced);
ideal Iall=I345,p6,p7;
ideal Gall=std(Iall);
print(
  "SAMPLE_INVARIANTS "+string(trace6)+" "+string(det(M6))+" "
  +string(trace7)+" "+string(det(M7))+" "+string(rank(M6))+" "
  +string(rank(M7))+" "+string(rank(joint))+" "
  +string(vdim(Greduced))+" "+string(Greduced[1]==1)+" "
  +string(vdim(Gall))+" "+string(Gall[1]==1)
);

ring cp={prime},(x),dp;
matrix N6=imap(tn,M6);
matrix N7=imap(tn,M7);
matrix X[12][12];
for (rowIndex=1;rowIndex<=12;rowIndex++) {{ X[rowIndex,rowIndex]=x; }}
poly characteristic6=det(X-N6);
poly characteristic7=det(X-N7);
matrix C6=coef(characteristic6,x);
matrix C7=coef(characteristic7,x);
for (columnIndex=1;columnIndex<=ncols(C6);columnIndex++)
{{
  print(
    "SAMPLE_CHAR6 "+string(leadexp(C6[1,columnIndex]))+" "
    +string(C6[2,columnIndex])
  );
}}
for (columnIndex=1;columnIndex<=ncols(C7);columnIndex++)
{{
  print(
    "SAMPLE_CHAR7 "+string(leadexp(C7[1,columnIndex]))+" "
    +string(C7[2,columnIndex])
  );
}}
""",
            text=True,
            capture_output=True,
            check=True,
            timeout=timeout,
        )
        quotient = re.search(
            r"(?m)^SAMPLE_QUOTIENT (\d+) (\d+)$",
            completed.stdout,
        )
        if quotient is None:
            raise AssertionError(completed.stdout[-8000:])
        if int(quotient.group(1)) != 12:
            print(
                f"TRACE_SAMPLE_SKIP point={point} length={quotient.group(1)}",
                flush=True,
            )
            continue
        invariant = re.search(
            r"(?m)^SAMPLE_INVARIANTS "
            r"(\S+) (\S+) (\S+) (\S+) (\d+) (\d+) (\d+) "
            r"(\d+) ([01]) (\d+) ([01])$",
            completed.stdout,
        )
        if invariant is None:
            raise AssertionError(completed.stdout[-8000:])

        def characteristic(order: int) -> list[int]:
            terms = re.findall(
                rf"(?m)^SAMPLE_CHAR{order} (\d+) (\S+)$",
                completed.stdout,
            )
            coefficients = [0] * 13
            for degree, coefficient in terms:
                coefficients[int(degree)] = int(coefficient) % prime
            return coefficients

        (
            trace6,
            norm6,
            trace7,
            norm7,
            rank6,
            rank7,
            joint_rank,
            reduced_common_length,
            reduced_unit,
            common_length,
            unit,
        ) = invariant.groups()
        if bool(int(reduced_unit)) != bool(int(unit)):
            raise AssertionError(
                "reduced and direct common-quotient unit tests disagree"
            )
        if bool(int(unit)) != (int(joint_rank) == 12):
            coordinate_errors = re.findall(
                r"(?m)^SAMPLE_COORDINATE_ERROR_[67] (.*)$",
                completed.stdout,
            )
            raise AssertionError(
                "joint multiplication rank and common-quotient test disagree: "
                f"point={point}, ranks=({rank6},{rank7},{joint_rank}), "
                f"lengths=({reduced_common_length},{common_length}), "
                f"coordinate_errors={coordinate_errors}"
            )
        sample = {
            "point": point,
            "quotient_length": int(quotient.group(1)),
            "groebner_basis_size": int(quotient.group(2)),
            "leading_exponents_s5_t4_s3": [
                [int(value) for value in exponent.split(",")]
                for exponent in re.findall(
                    r"(?m)^SAMPLE_LEADING ([0-9,]+)$",
                    completed.stdout,
                )
            ],
            "mu6": {
                "trace": int(trace6) % prime,
                "norm": int(norm6) % prime,
                "multiplication_rank": int(rank6),
                "characteristic_coefficients_ascending": characteristic(6),
            },
            "mu7": {
                "trace": int(trace7) % prime,
                "norm": int(norm7) % prime,
                "multiplication_rank": int(rank7),
                "characteristic_coefficients_ascending": characteristic(7),
            },
            "joint_multiplication_rank": int(joint_rank),
            "reduced_common_quotient_length": int(reduced_common_length),
            "common_quotient_length": int(common_length),
            "common_quotient_is_unit": bool(int(unit)),
        }
        samples.append(sample)
        if (
            sample["mu6"]["norm"] == 0
            or sample["mu7"]["norm"] == 0
            or int(joint_rank) < 12
        ):
            print(
                f"TRACE_SAMPLE_EXCEPTIONAL point={point} "
                f"norm6={sample['mu6']['norm']} "
                f"norm7={sample['mu7']['norm']} "
                f"joint_rank={joint_rank} common_length={common_length}",
                flush=True,
            )

    output = (
        ROOT
        / "artifacts"
        / "generated-results"
        / f"two_pair_sic_bidegree33_boundary_trace_norm_samples_mod{prime}.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "prime": prime,
        "normalization": "s0=1, A=B=mu2=0, L=1",
        "quotient_generators": ["mu3", "mu4", "mu5"],
        "basis": [
            "1", "t4", "t4^2", "t4^3", "s5", "s5*t4",
            "s3", "s3*t4", "s3*t4^2", "s3*t4^3", "s3*s5", "s3*s5*t4",
        ],
        "samples": samples,
        "scope": (
            "exact finite-field evaluations of multiplication traces, norms, "
            "characteristic polynomials, and the joint Fitting rank for "
            "mu6 and mu7 in the rank-twelve quotient"
        ),
    }
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"TRACE_SAMPLE_WROTE {output.relative_to(ROOT)}", flush=True)


def trace_norm_basis_probe(
    singular: str,
    prime: int,
    orders: tuple[int, ...],
    polynomials: list[str],
    timeout: int,
    algorithm: str,
    variable_order: str,
) -> None:
    """Benchmark the rank-twelve mu3,mu4,mu5 basis without matrices."""

    available = dict(zip(orders, polynomials))
    assert all(order in available for order in (3, 4, 5))
    replacements = (
        ("t1", "(s1*t0-L)"),
        ("s2", "(s1^2-(13/3)*t0^2-Q)"),
        ("L", "(1)"),
    )
    available = {
        order: substitute(polynomial, replacements)
        for order, polynomial in available.items()
    }
    variables = {
        "s3-s5-t4": "s3,s5,t4",
        "s5-t4-s3": "s5,t4,s3",
        "t4-s5-s3": "t4,s5,s3",
    }[variable_order]
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring tn=({prime},s1,t0,Q,t2),({variables}),dp;
option(redSB);
poly p3={available[3]};
poly p4={available[4]};
poly p5={available[5]};
ideal I345=p3,p4,p5;
ideal G={algorithm}(I345);
print(
  "TRACE_NORM_BASIS "+string(dim(G))+" "+string(vdim(G))+" "
  +string(size(G))
);
int basisIndex;
for (basisIndex=1;basisIndex<=size(G);basisIndex++)
{{
  print("TRACE_NORM_BASIS_LEADING "+string(leadexp(G[basisIndex])));
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    marker = re.search(
        r"(?m)^TRACE_NORM_BASIS (-?\d+) (-?\d+) (\d+)$",
        completed.stdout,
    )
    if marker is None:
        raise AssertionError(completed.stdout[-8000:])
    leading = re.findall(
        r"(?m)^TRACE_NORM_BASIS_LEADING ([0-9,]+)$",
        completed.stdout,
    )
    print(
        f"TRACE_NORM_BASIS algorithm={algorithm} order={variable_order} "
        f"dimension={marker.group(1)} length={marker.group(2)} "
        f"basis={marker.group(3)}",
        flush=True,
    )
    print(f"TRACE_NORM_BASIS_INITIAL {','.join(leading)}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=47)
    parser.add_argument("--orders", default="2,3,4,5,6")
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument("--adapted", action="store_true")
    parser.add_argument(
        "--backend",
        choices=("singular", "msolve"),
        default="singular",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--linear-algebra",
        choices=(1, 2, 42, 44),
        type=int,
        default=2,
    )
    parser.add_argument(
        "--coefficient-counts",
        default="",
        help=(
            "optional comma-separated order:count limits for the adapted "
            "msolve base ideal, for example 6:3,7:1"
        ),
    )
    parser.add_argument(
        "--normalize-l",
        action="store_true",
        help="use weighted homogeneity and L!=0 to impose L=1",
    )
    parser.add_argument("--common-root", action="store_true")
    parser.add_argument("--no-lift", action="store_true")
    parser.add_argument("--multiplication", action="store_true")
    parser.add_argument("--norm-only", action="store_true")
    parser.add_argument("--matrix-only", action="store_true")
    parser.add_argument("--trace-norm", action="store_true")
    parser.add_argument(
        "--trace-samples",
        type=int,
        default=0,
        help=(
            "with --trace-norm, specialize the adapted base at this many "
            "deterministic points and compute exact mu6/mu7 invariants"
        ),
    )
    parser.add_argument("--basis-only", action="store_true")
    parser.add_argument(
        "--trace-algorithm",
        choices=("std", "slimgb"),
        default="slimgb",
    )
    parser.add_argument(
        "--trace-variable-order",
        choices=("s3-s5-t4", "s5-t4-s3", "t4-s5-s3"),
        default="s3-s5-t4",
    )
    arguments = parser.parse_args()
    orders = tuple(int(value) for value in arguments.orders.split(",") if value)
    coefficient_counts = {
        int(order): int(count)
        for item in arguments.coefficient_counts.split(",")
        if item
        for order, count in (item.split(":", 1),)
    }
    assert orders[:4] == (2, 3, 4, 5)
    assert arguments.prime == 0 or 3 * max(orders) < arguments.prime
    singular = shutil.which("Singular")
    assert singular is not None
    msolve = shutil.which("msolve")
    if arguments.backend == "msolve":
        assert msolve is not None

    if arguments.prime == 0:
        expressions = [
            exact_chart_expression(exact_moment_terms(order))
            for order in orders
        ]
    else:
        expressions = [
            chart_expression(
                moment_terms(order, arguments.prime),
                0,
                arguments.prime,
            )
            for order in orders
        ]
    variables, polynomials = prepare_s0_branch_for_msolve(
        singular,
        expressions,
        arguments.prime,
        "s0-boundary",
        arguments.timeout,
    )
    assert variables == ("s1", "s2", "s3", "s5", "t0", "t1", "t2", "t4")
    reduced_orders = orders[1:]
    if arguments.trace_norm and arguments.trace_samples:
        trace_norm_samples(
            singular,
            arguments.prime,
            reduced_orders,
            polynomials,
            arguments.timeout,
            arguments.trace_samples,
        )
    elif arguments.trace_norm and arguments.basis_only:
        trace_norm_basis_probe(
            singular,
            arguments.prime,
            reduced_orders,
            polynomials,
            arguments.timeout,
            arguments.trace_algorithm,
            arguments.trace_variable_order,
        )
    elif arguments.trace_norm:
        trace_norm_export(
            singular,
            arguments.prime,
            reduced_orders,
            polynomials,
            arguments.timeout,
        )
    elif arguments.multiplication:
        multiplication_fitting(
            singular,
            arguments.prime,
            reduced_orders,
            polynomials,
            arguments.timeout,
            arguments.normalize_l,
            arguments.norm_only,
            arguments.matrix_only,
        )
    elif arguments.common_root:
        generic_common_root(
            singular,
            arguments.prime,
            reduced_orders,
            polynomials,
            arguments.timeout,
            not arguments.no_lift,
            arguments.normalize_l,
        )
    elif arguments.adapted:
        base_polynomial, coefficients = adapted_coefficients(
            singular,
            arguments.prime,
            reduced_orders,
            polynomials,
            arguments.timeout,
        )
        if arguments.backend == "msolve":
            adapted_base_msolve(
                msolve,
                arguments.prime,
                reduced_orders,
                base_polynomial,
                coefficients,
                arguments.timeout,
                arguments.threads,
                arguments.linear_algebra,
                coefficient_counts,
                arguments.normalize_l,
            )
        else:
            adapted_base_ideal(
                singular,
                arguments.prime,
                reduced_orders,
                base_polynomial,
                coefficients,
                arguments.timeout,
            )
    else:
        coefficients = quotient_coefficients(
            singular,
            arguments.prime,
            reduced_orders,
            polynomials,
            arguments.timeout,
        )
        base_ideal(
            singular,
            arguments.prime,
            reduced_orders,
            polynomials[0],
            coefficients,
            arguments.timeout,
        )
    if arguments.trace_norm:
        print(
            "EVIDENCE ONLY: modular finite-algebra invariants on L=1; "
            "no characteristic-zero theorem is promoted"
        )
    else:
        print(
            "EVIDENCE ONLY: modular L*Q*J-open coefficient ideal; "
            "no characteristic-zero theorem is promoted"
        )


if __name__ == "__main__":
    main()
