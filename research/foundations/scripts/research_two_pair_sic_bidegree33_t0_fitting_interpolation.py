#!/usr/bin/env python3
"""Bounded modular samples of the t0-open determinant pencil.

The transformed base equation mu_3 is quadratic in s3.  At a five-base
specialization where it has two distinct roots, evaluating the rank-six
fiber at both roots recovers the two coordinates A+B*s3 of any rational
function modulo mu_3.  This script computes those paired evaluations for all
seven coefficients of

    det(M_mu6 + z*M_mu7) = c0 + c1*z + ... + c6*z^6.

Directional mode reconstructs c0 and c1 as before.  Random mode instead scans
bounded deterministic shards for specializations where the whole determinant
pencil vanishes.  This is a candidate search, not a global zero-fiber
certificate.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import combinations
import json
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys

from verify_two_pair_sic_bidegree33_corrected_boundary import (
    CORRECTED_ORDERS,
    t0_open_localized_export,
)
from explore_two_pair_sic_bidegree33_full_anchor import (
    PARAMETERS,
    moment_terms,
)


ROOT = Path(__file__).resolve().parents[1]
BASE_VARIABLES = ("s1", "s2", "t1", "t2", "u")
FIBER_VARIABLES = ("s6", "s5")
ALL_VARIABLES = (*FIBER_VARIABLES, "s1", "s2", "s3", "t1", "t2", "u")
DEFAULT_BASE = {
    "s1": 5,
    "s2": 7,
    "t1": 11,
    "t2": 13,
    "u": 17,
}
STRATA = (
    "generic",
    "Q",
    "J",
    "K",
    "H",
    "KH",
    "QJH",
    "JH",
    "JK",
    "a2",
    "discriminant",
)
BORDER_NAMES = (
    "pivot_A",
    "pivot_B",
    "A_M6_c6",
    "A_M7_c1",
    "A_M7_c2",
    "A_M7_c3",
    "A_M7_c4",
    "A_M7_c5",
    "A_M7_c6",
    "A_M8_c1",
)

Monomial = tuple[int, ...]
SparsePolynomial = dict[Monomial, int]


def parse_sparse_mod(polynomial: str, prime: int) -> SparsePolynomial:
    """Parse the expanded Singular serialization used by the exporter."""

    terms = re.findall(r"[+-]?[^+-]+", polynomial)
    assert terms and "".join(terms) == polynomial
    result: dict[Monomial, int] = defaultdict(int)
    variable_index = {
        variable: index for index, variable in enumerate(ALL_VARIABLES)
    }
    for term in terms:
        sign = -1 if term.startswith("-") else 1
        body = term[1:] if term[:1] in "+-" else term
        coefficient = sign
        exponents = [0] * len(ALL_VARIABLES)
        for factor in body.split("*"):
            if re.fullmatch(r"\d+", factor):
                coefficient = coefficient * int(factor) % prime
                continue
            match = re.fullmatch(r"([a-z]\w*)(?:\^(\d+))?", factor)
            assert match is not None, factor
            variable = match.group(1)
            assert variable in variable_index, variable
            exponents[variable_index[variable]] += int(match.group(2) or 1)
        key = tuple(exponents)
        result[key] = (result[key] + coefficient) % prime
    return {monomial: value for monomial, value in result.items() if value}


def evaluate_base_to_bivariate(
    polynomial: SparsePolynomial,
    base: dict[str, int],
    s3_value: int,
    prime: int,
) -> dict[tuple[int, int], int]:
    """Evaluate all coordinates except the two fiber variables."""

    values = {**base, "s3": s3_value}
    result: dict[tuple[int, int], int] = defaultdict(int)
    for exponents, coefficient in polynomial.items():
        scalar = coefficient
        for index, variable in enumerate(ALL_VARIABLES[2:], start=2):
            scalar = (
                scalar
                * pow(values[variable] % prime, exponents[index], prime)
            ) % prime
        fiber = (exponents[0], exponents[1])
        result[fiber] = (result[fiber] + scalar) % prime
    return {monomial: value for monomial, value in result.items() if value}


def evaluate_mu3_coefficients(
    polynomial: SparsePolynomial,
    base: dict[str, int],
    prime: int,
) -> list[int]:
    """Return the coefficients of mu_3 as a polynomial in s3."""

    coefficients = [0, 0, 0]
    for exponents, coefficient in polynomial.items():
        assert exponents[0] == exponents[1] == 0
        s3_exponent = exponents[4]
        assert s3_exponent <= 2
        scalar = coefficient
        for index, variable in enumerate(ALL_VARIABLES[2:], start=2):
            if variable == "s3":
                continue
            scalar = (
                scalar
                * pow(base[variable] % prime, exponents[index], prime)
            ) % prime
        coefficients[s3_exponent] = (
            coefficients[s3_exponent] + scalar
        ) % prime
    return coefficients


def square_root_mod(value: int, prime: int) -> int | None:
    """Return one square root for primes congruent to three modulo four."""

    assert prime % 4 == 3
    value %= prime
    root = pow(value, (prime + 1) // 4, prime)
    return root if root * root % prime == value else None


def quadratic_roots(
    coefficients: list[int],
    prime: int,
) -> tuple[int, int] | None:
    constant, linear, quadratic = coefficients
    if quadratic == 0:
        return None
    discriminant = (linear * linear - 4 * quadratic * constant) % prime
    root = square_root_mod(discriminant, prime)
    if root is None or root == 0:
        return None
    denominator_inverse = pow(2 * quadratic, -1, prime)
    left = (-linear + root) * denominator_inverse % prime
    right = (-linear - root) * denominator_inverse % prime
    assert left != right
    return tuple(sorted((left, right)))


def generic_factors(base: dict[str, int], prime: int) -> dict[str, int]:
    s1 = base["s1"] % prime
    s2 = base["s2"] % prime
    t1 = base["t1"] % prime
    u = base["u"] % prime
    x = (s1 * s1 * u - s2) % prime
    ell = (s1 * u - t1) % prime
    a_factor = (99 * x - 274 * u) % prime
    q_factor = (3 * x - 13 * u) % prime
    j_factor = (a_factor * a_factor + 30420 * ell * ell) % prime
    k_factor = (351 * x - 901 * u) % prime
    h_factor = (
        a_factor * k_factor + 121680 * ell * ell
    ) % prime
    return {
        "u": u,
        "Q": q_factor,
        "J": j_factor,
        "K": k_factor,
        "H": h_factor,
    }


def random_stratum_base(
    generator: random.Random,
    prime: int,
    stratum: str,
    mu3: SparsePolynomial | None = None,
) -> dict[str, int] | None:
    """Generate one base point on a selected coefficient divisor."""

    assert stratum in STRATA
    if stratum == "generic":
        return {
            variable: generator.randrange(prime)
            for variable in BASE_VARIABLES
        }
    s1 = generator.randrange(prime)
    t2 = generator.randrange(prime)
    u = generator.randrange(1, prime)
    if stratum == "a2":
        s2 = generator.randrange(prime)
        t1 = generator.randrange(prime)
        t2 = (
            2 * s1 * t1 * u
            - s2 * u
            - 8 * u * u * pow(9, -1, prime)
        ) % prime
        return {
            "s1": s1,
            "s2": s2,
            "t1": t1,
            "t2": t2,
            "u": u,
        }
    if stratum == "discriminant":
        assert mu3 is not None
        s2 = generator.randrange(prime)
        t1 = generator.randrange(prime)
        candidates = []
        for candidate_t2 in range(prime):
            base = {
                "s1": s1,
                "s2": s2,
                "t1": t1,
                "t2": candidate_t2,
                "u": u,
            }
            constant, linear, quadratic = evaluate_mu3_coefficients(
                mu3,
                base,
                prime,
            )
            if quadratic and (
                linear * linear - 4 * quadratic * constant
            ) % prime == 0:
                candidates.append(base)
        return generator.choice(candidates) if candidates else None
    if stratum in {"Q", "K"}:
        t1 = generator.randrange(prime)
        x = (
            13 * u * pow(3, -1, prime)
            if stratum == "Q"
            else 901 * u * pow(351, -1, prime)
        ) % prime
        return {
            "s1": s1,
            "s2": (s1 * s1 * u - x) % prime,
            "t1": t1,
            "t2": t2,
            "u": u,
        }
    if stratum == "KH":
        x = 901 * u * pow(351, -1, prime) % prime
        return {
            "s1": s1,
            "s2": (s1 * s1 * u - x) % prime,
            "t1": s1 * u % prime,
            "t2": t2,
            "u": u,
        }
    if stratum == "QJH":
        alpha = square_root_mod(-30420, prime)
        if alpha is None:
            raise ValueError(
                f"-30420 is not a square modulo {prime}; choose another prime"
            )
        x = 13 * u * pow(3, -1, prime) % prime
        a_factor = 155 * u % prime
        ell = a_factor * pow(alpha, -1, prime) % prime
        return {
            "s1": s1,
            "s2": (s1 * s1 * u - x) % prime,
            "t1": (s1 * u - ell) % prime,
            "t2": t2,
            "u": u,
        }
    if stratum == "JH":
        x = 274 * u * pow(99, -1, prime) % prime
        return {
            "s1": s1,
            "s2": (s1 * s1 * u - x) % prime,
            "t1": s1 * u % prime,
            "t2": t2,
            "u": u,
        }
    if stratum == "JK":
        alpha = square_root_mod(-30420, prime)
        if alpha is None:
            raise ValueError(
                f"-30420 is not a square modulo {prime}; choose another prime"
            )
        x = 901 * u * pow(351, -1, prime) % prime
        a_factor = (99 * x - 274 * u) % prime
        ell = a_factor * pow(alpha, -1, prime) % prime
        return {
            "s1": s1,
            "s2": (s1 * s1 * u - x) % prime,
            "t1": (s1 * u - ell) % prime,
            "t2": t2,
            "u": u,
        }
    if stratum == "J":
        alpha = square_root_mod(-30420, prime)
        if alpha is None:
            raise ValueError(
                f"-30420 is not a square modulo {prime}; choose another prime"
            )
        ell = generator.randrange(prime)
        a_factor = alpha * ell % prime
        x = (a_factor + 274 * u) * pow(99, -1, prime) % prime
        return {
            "s1": s1,
            "s2": (s1 * s1 * u - x) % prime,
            "t1": (s1 * u - ell) % prime,
            "t2": t2,
            "u": u,
        }
    assert stratum == "H"
    parameter = generator.randrange(prime)
    denominator = (1287 * parameter * parameter + 40560) % prime
    if denominator == 0:
        return None
    ell = (
        -775 * parameter * u * pow(denominator, -1, prime)
    ) % prime
    q_factor = (
        -155 * u * pow(33, -1, prime) + parameter * ell
    ) % prime
    x = (q_factor + 13 * u) * pow(3, -1, prime) % prime
    return {
        "s1": s1,
        "s2": (s1 * s1 * u - x) % prime,
        "t1": (s1 * u - ell) % prime,
        "t2": t2,
        "u": u,
    }


def serialize_bivariate(
    polynomial: dict[tuple[int, int], int],
    prime: int,
) -> str:
    terms = []
    for (s6_exponent, s5_exponent), coefficient in sorted(
        polynomial.items(),
        reverse=True,
    ):
        coefficient %= prime
        if not coefficient:
            continue
        factors = [str(coefficient)]
        if s6_exponent:
            factors.append("s6" if s6_exponent == 1 else f"s6^{s6_exponent}")
        if s5_exponent:
            factors.append("s5" if s5_exponent == 1 else f"s5^{s5_exponent}")
        terms.append("*".join(factors))
    return "+".join(terms) or "0"


def singular_procedure(
    weights: list[list[int]],
    emit_all_borders: bool = False,
) -> str:
    coefficient_expressions = [
        "+".join(
            f"{weight}*leadcoef(d{value})"
            for value, weight in enumerate(coefficient_weights)
            if weight
        )
        or "0"
        for coefficient_weights in weights
    ]
    coefficient_declarations = "\n  ".join(
        f"number c{index}={expression};"
        for index, expression in enumerate(coefficient_expressions)
    )
    coefficient_print = '+ " " +'.join(
        f"string(c{index})" for index in range(7)
    )
    candidate_condition = (
        "(c0==0) && (c1==0) && (c2==0) && (c3==0)"
        " && (c4==0) && (c5==0) && (c6==0)"
    )
    border_output = ""
    m7_borders = ""
    if emit_all_borders:
        m7_borders = "\n".join(
            f"""
    for(basisRow=1;basisRow<=6;basisRow++)
    {{
      borderMatrix[basisRow,6]=M7[basisRow,{column}];
    }}
    borderValue=det(borderMatrix);
    print(
      "BORDER "+string(sampleId)+" A_M7_c{column} "
      +string(borderValue)
    );
"""
            for column in range(1, 7)
        )
    border_output = f"""
    matrix pivotA[5][5];
    matrix pivotB[5][5];
    matrix borderMatrix[6][6];
    poly borderValue;
    for(basisRow=1;basisRow<=5;basisRow++)
    {{
      for(basisColumn=1;basisColumn<=5;basisColumn++)
      {{
        pivotA[basisRow,basisColumn]=M6[basisRow,basisColumn];
      }}
      for(basisColumn=1;basisColumn<=4;basisColumn++)
      {{
        pivotB[basisRow,basisColumn]=M6[basisRow,basisColumn];
      }}
      pivotB[basisRow,5]=M7[basisRow,1];
    }}
    print(
      "BORDER "+string(sampleId)+" pivot_A "+string(det(pivotA))
    );
    print(
      "BORDER "+string(sampleId)+" pivot_B "+string(det(pivotB))
    );
    for(basisRow=1;basisRow<=6;basisRow++)
    {{
      for(basisColumn=1;basisColumn<=5;basisColumn++)
      {{
        borderMatrix[basisRow,basisColumn]=M6[basisRow,basisColumn];
      }}
      borderMatrix[basisRow,6]=M6[basisRow,6];
    }}
    borderValue=det(borderMatrix);
    print(
      "BORDER "+string(sampleId)+" A_M6_c6 "+string(borderValue)
    );
{m7_borders}
    for(basisRow=1;basisRow<=6;basisRow++)
    {{
      borderMatrix[basisRow,6]=M8[basisRow,1];
    }}
    borderValue=det(borderMatrix);
    print(
      "BORDER "+string(sampleId)+" A_M8_c1 "+string(borderValue)
    );
    matrix rankMatrix67[6][12];
    matrix rankMatrix678[6][18];
    for(basisRow=1;basisRow<=6;basisRow++)
    {{
      for(basisColumn=1;basisColumn<=6;basisColumn++)
      {{
        rankMatrix67[basisRow,basisColumn]
          =M6[basisRow,basisColumn];
        rankMatrix67[basisRow,6+basisColumn]
          =M7[basisRow,basisColumn];
        rankMatrix678[basisRow,basisColumn]
          =M6[basisRow,basisColumn];
        rankMatrix678[basisRow,6+basisColumn]
          =M7[basisRow,basisColumn];
        rankMatrix678[basisRow,12+basisColumn]
          =M8[basisRow,basisColumn];
      }}
    }}
    print(
      "RANK "+string(sampleId)+" "+string(rank(M6))+" "
      +string(rank(rankMatrix67))+" "+string(rank(rankMatrix678))
    );
"""
    if not emit_all_borders:
        border_output = ""
    matrix_output = f"""
    if({candidate_condition})
    {{
      for(basisRow=1;basisRow<=6;basisRow++)
      {{
        for(basisColumn=1;basisColumn<=6;basisColumn++)
        {{
          print(
            "MATRIX "+string(sampleId)+" 6 "
            +string(basisRow)+" "+string(basisColumn)+" "
            +string(M6[basisRow,basisColumn])
          );
          print(
            "MATRIX "+string(sampleId)+" 7 "
            +string(basisRow)+" "+string(basisColumn)+" "
            +string(M7[basisRow,basisColumn])
          );
          print(
            "MATRIX "+string(sampleId)+" 8 "
            +string(basisRow)+" "+string(basisColumn)+" "
            +string(M8[basisRow,basisColumn])
          );
        }}
      }}
    }}
"""
    m8_condition = "1" if emit_all_borders else candidate_condition
    m8_and_outputs = f"""
  if({m8_condition})
  {{
    poly r8=reduce(p8,G);
    z=r8;
    while(z!=0)
    {{
      basisRow=coordinateIndex(leadmonom(z));
      if(basisRow==1) {{ M8=M8+leadcoef(z)*B1; }}
      if(basisRow==2) {{ M8=M8+leadcoef(z)*B2; }}
      if(basisRow==3) {{ M8=M8+leadcoef(z)*B3; }}
      if(basisRow==4) {{ M8=M8+leadcoef(z)*B4; }}
      if(basisRow==5) {{ M8=M8+leadcoef(z)*B5; }}
      if(basisRow==6) {{ M8=M8+leadcoef(z)*B6; }}
      z=z-lead(z);
    }}
{border_output}
{matrix_output}
  }}
"""
    return f"""
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
proc fittingSample(poly p4, poly p5, poly p6, poly p7, poly p8, int sampleId)
{{
  ideal G=std(p4,p5);
  if((dim(G)!=0) || (vdim(G)!=6))
  {{
    print("BAD "+string(sampleId)+" "+string(dim(G))+" "+string(vdim(G)));
    return();
  }}
  poly r6=reduce(p6,G);
  poly r7=reduce(p7,G);
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
  matrix M6[6][6];
  matrix M7[6][6];
  matrix M8[6][6];
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
      if(basisRow==0) {{ print("COORDINATE_ERROR"); return(); }}
      B2[basisRow,basisColumn]=leadcoef(z);
      z=z-lead(z);
    }}
    z=reduce(s5*basisMonomial,G);
    while(z!=0)
    {{
      basisRow=coordinateIndex(leadmonom(z));
      if(basisRow==0) {{ print("COORDINATE_ERROR"); return(); }}
      B3[basisRow,basisColumn]=leadcoef(z);
      z=z-lead(z);
    }}
  }}
  B4=B2*B3;
  B5=B3*B3;
  B6=B5*B3;
  z=r6;
  while(z!=0)
  {{
    basisRow=coordinateIndex(leadmonom(z));
    if(basisRow==1) {{ M6=M6+leadcoef(z)*B1; }}
    if(basisRow==2) {{ M6=M6+leadcoef(z)*B2; }}
    if(basisRow==3) {{ M6=M6+leadcoef(z)*B3; }}
    if(basisRow==4) {{ M6=M6+leadcoef(z)*B4; }}
    if(basisRow==5) {{ M6=M6+leadcoef(z)*B5; }}
    if(basisRow==6) {{ M6=M6+leadcoef(z)*B6; }}
    z=z-lead(z);
  }}
  z=r7;
  while(z!=0)
  {{
    basisRow=coordinateIndex(leadmonom(z));
    if(basisRow==1) {{ M7=M7+leadcoef(z)*B1; }}
    if(basisRow==2) {{ M7=M7+leadcoef(z)*B2; }}
    if(basisRow==3) {{ M7=M7+leadcoef(z)*B3; }}
    if(basisRow==4) {{ M7=M7+leadcoef(z)*B4; }}
    if(basisRow==5) {{ M7=M7+leadcoef(z)*B5; }}
    if(basisRow==6) {{ M7=M7+leadcoef(z)*B6; }}
    z=z-lead(z);
  }}
  poly d0=det(M6);
  poly d1=det(M6+M7);
  poly d2=det(M6+2*M7);
  poly d3=det(M6+3*M7);
  poly d4=det(M6+4*M7);
  poly d5=det(M6+5*M7);
  poly d6=det(M6+6*M7);
  {coefficient_declarations}
  print(
    "SAMPLE "+string(sampleId)+" "+{coefficient_print}
  );
{m8_and_outputs}
}}
"""


def vandermonde_weights(prime: int, coefficient: int) -> list[int]:
    """Return weights extracting one coefficient from values at 0,...,6."""

    augmented = [
        [pow(value, degree, prime) for value in range(7)]
        + [1 if degree == coefficient else 0]
        for degree in range(7)
    ]
    for column in range(7):
        pivot = next(
            row for row in range(column, 7)
            if augmented[row][column] % prime
        )
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        inverse = pow(augmented[column][column], -1, prime)
        augmented[column] = [
            value * inverse % prime for value in augmented[column]
        ]
        for row in range(7):
            if row == column:
                continue
            scalar = augmented[row][column]
            augmented[row] = [
                (left - scalar * right) % prime
                for left, right in zip(
                    augmented[row], augmented[column], strict=True
                )
            ]
    weights = [augmented[row][-1] for row in range(7)]
    assert all(
        sum(weights[value] * pow(value, degree, prime) for value in range(7))
        % prime
        == (1 if degree == coefficient else 0)
        for degree in range(7)
    )
    return weights


def run_samples(
    singular: str,
    prime: int,
    samples: list[dict[str, object]],
    timeout: int,
    pivot_scout: bool = False,
) -> tuple[
    dict[int, tuple[int, ...]],
    dict[int, dict[int, list[list[int]]]],
    dict[int, dict[str, int]],
    dict[int, dict[str, int]],
]:
    weights = [
        vandermonde_weights(prime, coefficient)
        for coefficient in range(7)
    ]
    blocks = [
        f"ring sample0={prime},(s6,s5),dp;",
        singular_procedure(weights, pivot_scout),
    ]
    for sample_id, sample in enumerate(samples):
        moments = sample["moments"]
        assert isinstance(moments, dict)
        blocks.append(f"ring sample{sample_id + 1}={prime},(s6,s5),dp;")
        blocks.extend(
            f"poly p{order}={moments[str(order)]};"
            for order in (4, 5, 6, 7)
        )
        blocks.append(f"poly p8={moments.get('8', '0')};")
        blocks.append(f"fittingSample(p4,p5,p6,p7,p8,{sample_id});")
    completed = subprocess.run(
        [singular, "-q"],
        input="\n".join(blocks),
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    assert "\n   ? " not in completed.stdout, (
        completed.stdout[-4000:],
        completed.stderr[-2000:],
    )
    bad = re.findall(r"(?m)^BAD .*$", completed.stdout)
    assert not bad, bad[:10]
    values = {}
    for match in re.findall(
        r"(?m)^SAMPLE (\d+)((?: -?\d+){7})$",
        completed.stdout,
    ):
        sample, serialized_values = match
        values[int(sample)] = tuple(
            int(value) % prime for value in serialized_values.split()
        )
    assert len(values) == len(samples), (
        len(values),
        len(samples),
        completed.stdout[-4000:],
    )
    matrices: dict[int, dict[int, list[list[int]]]] = {}
    for sample, order, row, column, value in re.findall(
        r"(?m)^MATRIX (\d+) ([678]) ([1-6]) ([1-6]) (-?\d+)$",
        completed.stdout,
    ):
        sample_index = int(sample)
        moment_order = int(order)
        sample_matrices = matrices.setdefault(
            sample_index,
            {
                current_order: [[0] * 6 for _ in range(6)]
                for current_order in (6, 7, 8)
            },
        )
        sample_matrices[moment_order][int(row) - 1][int(column) - 1] = (
            int(value) % prime
        )
    assert all(
        set(sample_matrices) == {6, 7, 8}
        for sample_matrices in matrices.values()
    )
    borders: dict[int, dict[str, int]] = defaultdict(dict)
    for sample, name, value in re.findall(
        r"(?m)^BORDER (\d+) ([A-Za-z0-9_]+) (-?\d+)$",
        completed.stdout,
    ):
        assert name in BORDER_NAMES
        borders[int(sample)][name] = int(value) % prime
    if pivot_scout:
        assert set(borders) == set(range(len(samples)))
        assert all(set(values) == set(BORDER_NAMES) for values in borders.values())
    else:
        assert not borders
    ranks = {
        int(sample): {
            "M6": int(rank_m6),
            "M6_M7": int(rank_m67),
            "M6_M7_M8": int(rank_m678),
        }
        for sample, rank_m6, rank_m67, rank_m678 in re.findall(
            r"(?m)^RANK (\d+) ([0-6]) ([0-6]) ([0-6])$",
            completed.stdout,
        )
    }
    if pivot_scout:
        assert set(ranks) == set(range(len(samples)))
    else:
        assert not ranks
    return values, matrices, dict(borders), ranks


def matrix_rank_mod(matrix: list[list[int]], prime: int) -> int:
    reduced = [
        [value % prime for value in row]
        for row in matrix
    ]
    row_count = len(reduced)
    column_count = len(reduced[0]) if reduced else 0
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (
                row
                for row in range(pivot_row, row_count)
                if reduced[row][column]
            ),
            None,
        )
        if pivot is None:
            continue
        reduced[pivot_row], reduced[pivot] = (
            reduced[pivot],
            reduced[pivot_row],
        )
        inverse = pow(reduced[pivot_row][column], -1, prime)
        reduced[pivot_row] = [
            value * inverse % prime for value in reduced[pivot_row]
        ]
        for row in range(row_count):
            if row == pivot_row:
                continue
            scalar = reduced[row][column]
            if scalar:
                reduced[row] = [
                    (left - scalar * right) % prime
                    for left, right in zip(
                        reduced[row],
                        reduced[pivot_row],
                        strict=True,
                    )
                ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def determinant_mod(matrix: list[list[int]], prime: int) -> int:
    assert matrix and len(matrix) == len(matrix[0])
    reduced = [
        [value % prime for value in row]
        for row in matrix
    ]
    determinant = 1
    for column in range(len(reduced)):
        pivot = next(
            (
                row
                for row in range(column, len(reduced))
                if reduced[row][column]
            ),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            reduced[column], reduced[pivot] = (
                reduced[pivot],
                reduced[column],
            )
            determinant = -determinant
        pivot_value = reduced[column][column]
        determinant = determinant * pivot_value % prime
        inverse = pow(pivot_value, -1, prime)
        for row in range(column + 1, len(reduced)):
            scalar = reduced[row][column] * inverse % prime
            for later_column in range(column, len(reduced)):
                reduced[row][later_column] = (
                    reduced[row][later_column]
                    - scalar * reduced[column][later_column]
                ) % prime
    return determinant % prime


def block_matrix(
    matrices: dict[int, list[list[int]]],
    orders: tuple[int, ...],
) -> list[list[int]]:
    return [
        [
            value
            for order in orders
            for value in matrices[order][row]
        ]
        for row in range(6)
    ]


def pivot_signature(
    matrices: dict[int, list[list[int]]],
    prime: int,
) -> dict[str, object]:
    b7 = block_matrix(matrices, (6, 7))
    b8 = block_matrix(matrices, (6, 7, 8))
    rank_b7 = matrix_rank_mod(b7, prime)
    rank_b8 = matrix_rank_mod(b8, prime)
    assert rank_b7 == 5 and rank_b8 == 6, (rank_b7, rank_b8)
    pivot_rows = None
    pivot_columns = None
    pivot_determinant = 0
    for rows in combinations(range(6), 5):
        for columns in combinations(range(12), 5):
            determinant = determinant_mod(
                [
                    [b7[row][column] for column in columns]
                    for row in rows
                ],
                prime,
            )
            if determinant:
                pivot_rows = rows
                pivot_columns = columns
                pivot_determinant = determinant
                break
        if pivot_rows is not None:
            break
    assert pivot_rows is not None and pivot_columns is not None
    completion_column = None
    completion_determinant = 0
    for column in range(12, 18):
        determinant = determinant_mod(
            [
                [
                    b8[row][current_column]
                    for current_column in (*pivot_columns, column)
                ]
                for row in range(6)
            ],
            prime,
        )
        if determinant:
            completion_column = column - 12
            completion_determinant = determinant
            break
    assert completion_column is not None
    return {
        "rank_B7": rank_b7,
        "rank_B8": rank_b8,
        "pivot_rows_zero_based": list(pivot_rows),
        "pivot_columns_B7_zero_based": list(pivot_columns),
        "pivot_determinant": pivot_determinant,
        "mu8_completion_column_zero_based": completion_column,
        "completion_determinant": completion_determinant,
        "completion_to_pivot_ratio": (
            completion_determinant * pow(pivot_determinant, -1, prime)
        ) % prime,
    }


def replay_common_pencil_zeros(
    singular: str,
    prime: int,
    samples: list[dict[str, object]],
    candidate_sample_indices: list[int],
    timeout: int,
) -> list[dict[str, object]]:
    """Directly replay candidate fibers with the four-generator ideal."""

    if not candidate_sample_indices:
        return []
    blocks = []
    for replay_id, sample_index in enumerate(candidate_sample_indices):
        sample = samples[sample_index]
        moments = sample["moments"]
        assert isinstance(moments, dict)
        blocks.append(f"ring replay{replay_id}={prime},(s6,s5),dp;")
        blocks.extend(
            f"poly p{order}={moments[str(order)]};"
            for order in (4, 5, 6, 7)
        )
        blocks.extend(
            [
                "option(redSB);",
                "ideal I=p4,p5,p6,p7;",
                "ideal C=std(I);",
                (
                    'print("REPLAY '
                    + str(replay_id)
                    + ' "+string(dim(C))+" "+string(vdim(C))'
                    '+ " "+string(size(C)));'
                ),
                "int generatorIndex;",
                "for(generatorIndex=1;generatorIndex<=size(C);generatorIndex++)",
                "{",
                (
                    '  print("GENERATOR '
                    + str(replay_id)
                    + ' "+string(generatorIndex)+" "+string(C[generatorIndex]));'
                ),
                "}",
            ]
        )
    completed = subprocess.run(
        [singular, "-q"],
        input="\n".join(blocks),
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    assert "\n   ? " not in completed.stdout, (
        completed.stdout[-4000:],
        completed.stderr[-2000:],
    )
    headers = {
        int(replay_id): {
            "dimension": int(dimension),
            "vector_space_dimension": int(vector_space_dimension),
            "standard_basis_size": int(size),
            "standard_basis": [],
        }
        for replay_id, dimension, vector_space_dimension, size in re.findall(
            r"(?m)^REPLAY (\d+) (-?\d+) (-?\d+) (\d+)$",
            completed.stdout,
        )
    }
    for replay_id, generator_index, polynomial in re.findall(
        r"(?m)^GENERATOR (\d+) (\d+) (.+)$",
        completed.stdout,
    ):
        assert int(generator_index) == len(
            headers[int(replay_id)]["standard_basis"]
        ) + 1
        headers[int(replay_id)]["standard_basis"].append(polynomial)
    assert len(headers) == len(candidate_sample_indices), (
        headers,
        completed.stdout[-4000:],
    )
    result = []
    for replay_id, sample_index in enumerate(candidate_sample_indices):
        datum = headers[replay_id]
        assert datum["dimension"] == 0
        assert datum["vector_space_dimension"] > 0
        assert len(datum["standard_basis"]) == datum["standard_basis_size"]
        result.append(
            {
                "sample_index": sample_index,
                **datum,
            }
        )
    return result


def run_direct_samples(
    singular: str,
    prime: int,
    samples: list[dict[str, object]],
    timeout: int,
) -> list[dict[str, object]]:
    """Compute direct common ideals without assuming a rank-six basis."""

    blocks = []
    for sample_id, sample in enumerate(samples):
        moments = sample["moments"]
        assert isinstance(moments, dict)
        blocks.append(f"ring direct{sample_id}={prime},(s6,s5),dp;")
        blocks.extend(
            f"poly p{order}={moments[str(order)]};"
            for order in (4, 5, 6, 7)
        )
        blocks.extend(
            [
                "option(redSB);",
                "ideal G45=std(p4,p5);",
                "ideal I47=p4,p5,p6,p7;",
                "ideal G47=std(I47);",
                (
                    'print("DIRECT '
                    + str(sample_id)
                    + ' "+string(dim(G45))+" "+string(vdim(G45))'
                    '+ " "+string(dim(G47))+" "+string(vdim(G47))'
                    '+ " "+string(size(G47)));'
                ),
                "int directGenerator;",
                "for(directGenerator=1;directGenerator<=size(G47);directGenerator++)",
                "{",
                (
                    '  print("DIRECT_GENERATOR '
                    + str(sample_id)
                    + ' "+string(directGenerator)+" "'
                    '+string(G47[directGenerator]));'
                ),
                "}",
            ]
        )
    completed = subprocess.run(
        [singular, "-q"],
        input="\n".join(blocks),
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    assert "\n   ? " not in completed.stdout, (
        completed.stdout[-4000:],
        completed.stderr[-2000:],
    )
    records = {
        int(sample): {
            "sample_index": int(sample),
            "mu4_mu5_dimension": int(dimension45),
            "mu4_mu5_length": int(length45),
            "through_mu7_dimension": int(dimension47),
            "through_mu7_length": int(length47),
            "through_mu7_standard_basis_size": int(size47),
            "through_mu7_standard_basis": [],
        }
        for (
            sample,
            dimension45,
            length45,
            dimension47,
            length47,
            size47,
        ) in re.findall(
            r"(?m)^DIRECT (\d+) (-?\d+) (-?\d+) (-?\d+) (-?\d+) (\d+)$",
            completed.stdout,
        )
    }
    for sample, generator, polynomial in re.findall(
        r"(?m)^DIRECT_GENERATOR (\d+) (\d+) (.+)$",
        completed.stdout,
    ):
        record = records[int(sample)]
        assert int(generator) == len(
            record["through_mu7_standard_basis"]
        ) + 1
        record["through_mu7_standard_basis"].append(polynomial)
    assert set(records) == set(range(len(samples))), (
        len(records),
        len(samples),
        completed.stdout[-4000:],
    )
    result = []
    for sample in range(len(samples)):
        record = records[sample]
        assert len(record["through_mu7_standard_basis"]) == record[
            "through_mu7_standard_basis_size"
        ]
        result.append(record)
    return result


def unique_fiber_point(
    standard_basis: list[str],
    prime: int,
) -> dict[str, int] | None:
    """Read a reduced length-one basis ``s5-a,s6-b``."""

    point = {}
    for polynomial in standard_basis:
        match = re.fullmatch(r"(s[56])(?:([+-])(\d+))?", polynomial)
        if match is None:
            continue
        variable, sign, magnitude = match.groups()
        constant = 0
        if magnitude is not None:
            constant = int(magnitude) * (-1 if sign == "-" else 1)
        point[variable] = -constant % prime
    if set(point) != set(FIBER_VARIABLES):
        return None
    return point


def full_parameter_point(
    base: dict[str, int],
    s3_value: int,
    fiber: dict[str, int],
    prime: int,
) -> dict[str, int]:
    """Undo the three linear pivots at one localized t0-open point."""

    s0 = pow(base["u"], -1, prime)
    s1 = base["s1"] % prime
    s2 = base["s2"] % prime
    s3 = s3_value % prime
    s5 = fiber["s5"] % prime
    s6 = fiber["s6"] % prime
    t1 = base["t1"] % prime
    t2 = base["t2"] % prime
    u = base["u"] % prime
    a_without_t3 = (
        6 * s1 * s1 * t1
        - 3 * s1 * s2
        - 3 * s0 * s1 * t2
        - 3 * s0 * s2 * t1
        + 2 * s0 * s3
        - 3 * s0
    ) % prime
    t3 = -a_without_t3 * u * u % prime
    b_without_s4 = (
        12 * s0 * s1 * s3
        + 28 * s1 * t1
        - 18 * s0 * s1
        - 9 * s0 * s2 * s2
        - 14 * s2
        - 2 * s0 * t2
        - 12 * s0 * t1 * t1
    ) % prime
    s4 = b_without_s4 * u * u * pow(3, -1, prime) % prime
    t4 = (
        3 * s0 * s6
        - 18 * s1 * s5
        + 45 * s2 * s4
        - 30 * s3 * s3
        + 56 * t1 * t3
        - 42 * t2 * t2
        - 70
    ) * pow(14, -1, prime) % prime
    return {
        "s0": s0,
        "s1": s1,
        "s2": s2,
        "s3": s3,
        "s4": s4,
        "s5": s5,
        "s6": s6,
        "t0": 1,
        "t1": t1,
        "t2": t2,
        "t3": t3,
        "t4": t4,
    }


def evaluate_moment_at_point(
    order: int,
    point: dict[str, int],
    prime: int,
) -> int:
    """Evaluate a corrected scalar moment without expanding a chart."""

    result = 0
    for exponents, coefficient in moment_terms(order, prime).items():
        term = coefficient
        for variable, exponent in zip(PARAMETERS, exponents, strict=True):
            term = term * pow(point[variable], exponent, prime) % prime
        result = (result + term) % prime
    return result


def paired_coordinates(
    left_root: int,
    right_root: int,
    left_value: int,
    right_value: int,
    prime: int,
) -> tuple[int, int]:
    slope = (
        (left_value - right_value)
        * pow(left_root - right_root, -1, prime)
    ) % prime
    intercept = (left_value - slope * left_root) % prime
    return intercept, slope


def polynomial_trim(polynomial: list[int], prime: int) -> list[int]:
    while len(polynomial) > 1 and polynomial[-1] % prime == 0:
        polynomial.pop()
    return [coefficient % prime for coefficient in polynomial]


def polynomial_add(
    left: list[int],
    right: list[int],
    prime: int,
    right_sign: int = 1,
) -> list[int]:
    result = [0] * max(len(left), len(right))
    for index, coefficient in enumerate(left):
        result[index] += coefficient
    for index, coefficient in enumerate(right):
        result[index] += right_sign * coefficient
    return polynomial_trim(result, prime)


def polynomial_scale(
    polynomial: list[int],
    scalar: int,
    prime: int,
) -> list[int]:
    return polynomial_trim(
        [coefficient * scalar % prime for coefficient in polynomial],
        prime,
    )


def polynomial_multiply(
    left: list[int],
    right: list[int],
    prime: int,
) -> list[int]:
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_coefficient in enumerate(left):
        for right_index, right_coefficient in enumerate(right):
            result[left_index + right_index] = (
                result[left_index + right_index]
                + left_coefficient * right_coefficient
            ) % prime
    return polynomial_trim(result, prime)


def polynomial_divmod(
    dividend: list[int],
    divisor: list[int],
    prime: int,
) -> tuple[list[int], list[int]]:
    remainder = polynomial_trim(dividend[:], prime)
    divisor = polynomial_trim(divisor[:], prime)
    quotient = [0] * max(1, len(remainder) - len(divisor) + 1)
    leading_inverse = pow(divisor[-1], -1, prime)
    while len(remainder) >= len(divisor) and remainder != [0]:
        degree = len(remainder) - len(divisor)
        scalar = remainder[-1] * leading_inverse % prime
        quotient[degree] = scalar
        for index, coefficient in enumerate(divisor):
            remainder[degree + index] = (
                remainder[degree + index] - scalar * coefficient
            ) % prime
        polynomial_trim(remainder, prime)
    return (
        polynomial_trim(quotient, prime),
        polynomial_trim(remainder, prime),
    )


def polynomial_evaluate(
    polynomial: list[int],
    value: int,
    prime: int,
) -> int:
    result = 0
    for coefficient in reversed(polynomial):
        result = (result * value + coefficient) % prime
    return result


def interpolate_polynomial(
    arguments: list[int],
    values: list[int],
    prime: int,
) -> tuple[list[int], list[int]]:
    """Return the interpolant and product of the linear sample factors."""

    divided_differences = values[:]
    for width in range(1, len(arguments)):
        for index in range(len(arguments) - 1, width - 1, -1):
            divided_differences[index] = (
                (divided_differences[index] - divided_differences[index - 1])
                * pow(
                    (arguments[index] - arguments[index - width]) % prime,
                    -1,
                    prime,
                )
            ) % prime
    interpolant = [0]
    basis = [1]
    for argument, coefficient in zip(
        arguments,
        divided_differences,
        strict=True,
    ):
        interpolant = polynomial_add(
            interpolant,
            polynomial_scale(basis, coefficient, prime),
            prime,
        )
        basis = polynomial_multiply(
            basis,
            [(-argument) % prime, 1],
            prime,
        )
    return interpolant, basis


def rational_reconstruct(
    arguments: list[int],
    values: list[int],
    prime: int,
) -> tuple[list[int], list[int]]:
    """Find the minimum-total-degree rational interpolant."""

    interpolant, modulus = interpolate_polynomial(
        arguments,
        values,
        prime,
    )
    left_remainder, right_remainder = modulus, interpolant
    left_denominator, right_denominator = [0], [1]
    candidates = []
    while right_remainder != [0]:
        if all(
            polynomial_evaluate(right_denominator, argument, prime)
            and (
                polynomial_evaluate(right_remainder, argument, prime)
                - value
                * polynomial_evaluate(
                    right_denominator,
                    argument,
                    prime,
                )
            )
            % prime
            == 0
            for argument, value in zip(arguments, values, strict=True)
        ):
            inverse = pow(right_denominator[-1], -1, prime)
            numerator = polynomial_scale(
                right_remainder,
                inverse,
                prime,
            )
            denominator = polynomial_scale(
                right_denominator,
                inverse,
                prime,
            )
            candidates.append(
                (
                    len(numerator) + len(denominator) - 2,
                    numerator,
                    denominator,
                )
            )
        quotient, remainder = polynomial_divmod(
            left_remainder,
            right_remainder,
            prime,
        )
        next_denominator = polynomial_add(
            left_denominator,
            polynomial_multiply(quotient, right_denominator, prime),
            prime,
            right_sign=-1,
        )
        left_remainder, right_remainder = right_remainder, remainder
        left_denominator, right_denominator = (
            right_denominator,
            next_denominator,
        )
    assert candidates
    _, numerator, denominator = min(
        candidates,
        key=lambda candidate: candidate[0],
    )
    return numerator, denominator


def reconstruct_coordinates(
    pairs: list[dict[str, object]],
    training_count: int,
    prime: int,
) -> dict[str, object]:
    training = pairs[:training_count]
    held_out = pairs[training_count:]
    arguments = [int(pair["variable_value"]) for pair in training]
    result = {}
    for coefficient in ("c0", "c1"):
        for coordinate in ("constant", "s3"):
            values = [
                int(pair[f"{coefficient}_mod_mu3"][coordinate])
                for pair in training
            ]
            numerator, denominator = rational_reconstruct(
                arguments,
                values,
                prime,
            )
            held_out_verified = all(
                polynomial_evaluate(
                    denominator,
                    int(pair["variable_value"]),
                    prime,
                )
                and (
                    polynomial_evaluate(
                        numerator,
                        int(pair["variable_value"]),
                        prime,
                    )
                    - int(pair[f"{coefficient}_mod_mu3"][coordinate])
                    * polynomial_evaluate(
                        denominator,
                        int(pair["variable_value"]),
                        prime,
                    )
                )
                % prime
                == 0
                for pair in held_out
            )
            assert held_out_verified
            result[f"{coefficient}_{coordinate}"] = {
                "numerator_degree": len(numerator) - 1,
                "denominator_degree": len(denominator) - 1,
                "numerator_coefficients_ascending": numerator,
                "denominator_coefficients_ascending": denominator,
                "held_out_count": len(held_out),
                "held_out_verified": True,
            }
    return result


def reconstruct_bordered_coordinates(
    pairs: list[dict[str, object]],
    training_count: int,
    prime: int,
) -> dict[str, object]:
    """Reconstruct directional residues of the pivot and bordered minors."""

    training = pairs[:training_count]
    held_out = pairs[training_count:]
    arguments = [int(pair["variable_value"]) for pair in training]
    result = {}
    for name in BORDER_NAMES:
        for coordinate in ("constant", "s3"):
            values = [
                int(pair["bordered_mod_mu3"][name][coordinate])
                for pair in training
            ]
            numerator, denominator = rational_reconstruct(
                arguments,
                values,
                prime,
            )
            held_out_verified = all(
                polynomial_evaluate(
                    denominator,
                    int(pair["variable_value"]),
                    prime,
                )
                and (
                    polynomial_evaluate(
                        numerator,
                        int(pair["variable_value"]),
                        prime,
                    )
                    - int(pair["bordered_mod_mu3"][name][coordinate])
                    * polynomial_evaluate(
                        denominator,
                        int(pair["variable_value"]),
                        prime,
                    )
                )
                % prime
                == 0
                for pair in held_out
            )
            assert held_out_verified
            result[f"{name}_{coordinate}"] = {
                "numerator_degree": len(numerator) - 1,
                "denominator_degree": len(denominator) - 1,
                "numerator_coefficients_ascending": numerator,
                "denominator_coefficients_ascending": denominator,
                "held_out_count": len(held_out),
                "held_out_verified": True,
            }
    return result


def monic_polynomial(
    polynomial: list[int],
    prime: int,
) -> list[int]:
    polynomial = polynomial_trim(polynomial[:], prime)
    return polynomial_scale(
        polynomial,
        pow(polynomial[-1], -1, prime),
        prime,
    )


def known_denominator_model(
    mu3: SparsePolynomial,
    varied_variable: str,
    prime: int,
    reconstruction: dict[str, object],
    base_template: dict[str, int],
) -> dict[str, object]:
    """Verify the observed a2,Q,J denominator formula on this line."""

    arguments = list(range(12))
    a2_values = []
    q_values = []
    j_values = []
    for value in arguments:
        base = {key: entry % prime for key, entry in base_template.items()}
        base[varied_variable] = value
        a2_values.append(
            evaluate_mu3_coefficients(mu3, base, prime)[2]
        )
        factors = generic_factors(base, prime)
        q_values.append(factors["Q"])
        j_values.append(factors["J"])
    a2, _ = interpolate_polynomial(arguments, a2_values, prime)
    q_factor, _ = interpolate_polynomial(arguments, q_values, prime)
    j_factor, _ = interpolate_polynomial(arguments, j_values, prime)
    models = {
        "c0": (41, 3, 3),
        "c1": (42, 4, 4),
    }
    result = {
        "factor_degrees": {
            "mu3_s3_quadratic_coefficient": len(a2) - 1,
            "Q": len(q_factor) - 1,
            "J": len(j_factor) - 1,
        },
        "models": {},
    }
    for coefficient, exponents in models.items():
        model = [1]
        for factor, exponent in zip(
            (a2, q_factor, j_factor),
            exponents,
            strict=True,
        ):
            for _ in range(exponent):
                model = polynomial_multiply(model, factor, prime)
        model = monic_polynomial(model, prime)
        coordinate_matches = {}
        for coordinate in ("constant", "s3"):
            datum = reconstruction[f"{coefficient}_{coordinate}"]
            observed = monic_polynomial(
                list(datum["denominator_coefficients_ascending"]),
                prime,
            )
            coordinate_matches[coordinate] = observed == model
        assert all(coordinate_matches.values())
        result["models"][coefficient] = {
            "formula": (
                "a2^41*Q^3*J^3"
                if coefficient == "c0"
                else "a2^42*Q^4*J^4"
            ),
            "degree": len(model) - 1,
            "coordinate_matches": coordinate_matches,
        }
    result["verified"] = True
    return result


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=1019)
    parser.add_argument("--variable", choices=BASE_VARIABLES, default="s1")
    parser.add_argument("--sample-count", type=int, default=16)
    parser.add_argument("--training-count", type=int, default=0)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument(
        "--base-values",
        default=",".join(
            f"{variable}={DEFAULT_BASE[variable]}"
            for variable in BASE_VARIABLES
        ),
    )
    parser.add_argument("--timeout", type=int, default=20)
    parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="scan deterministic random five-base specializations",
    )
    parser.add_argument(
        "--max-attempts",
        type=int,
        default=0,
        help="candidate bases tried before stopping (default: 20*sample-count)",
    )
    parser.add_argument(
        "--retain-pairs",
        action="store_true",
        help="retain all paired residues in random-mode output",
    )
    parser.add_argument("--stratum", choices=STRATA, default="generic")
    parser.add_argument(
        "--pivot-scout",
        action="store_true",
        help="evaluate the two pivots, B7 borders, and first mu8 border",
    )
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    prime = arguments.prime
    assert prime > 42 and prime % 4 == 3
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    base_template = {}
    for assignment in arguments.base_values.split(","):
        variable, value = assignment.split("=", 1)
        assert variable in BASE_VARIABLES and variable not in base_template
        base_template[variable] = int(value) % prime
    assert set(base_template) == set(BASE_VARIABLES)
    export = t0_open_localized_export(
        singular,
        (2, 3, 4, 5, 6, 7, 8),
        prime,
        arguments.timeout,
    )
    parsed = {
        order: parse_sparse_mod(polynomial, prime)
        for order, polynomial in zip(
            (3, 4, 5, 6, 7, 8),
            export["polynomials"][:-1],
            strict=True,
        )
    }
    samples = []
    pair_metadata = []
    candidate = arguments.start
    attempts = 0
    maximum_attempts = arguments.max_attempts or 20 * arguments.sample_count
    generator = (
        random.Random(arguments.random_seed)
        if arguments.random_seed is not None
        else None
    )
    if arguments.stratum != "generic":
        assert generator is not None
        assert not arguments.training_count
        assert not arguments.pivot_scout
    while (
        len(pair_metadata) < arguments.sample_count
        and attempts < maximum_attempts
        and (generator is not None or candidate < prime)
    ):
        if generator is None:
            base = dict(base_template)
            base[arguments.variable] = candidate
            candidate += 1
        else:
            base = random_stratum_base(
                generator,
                prime,
                arguments.stratum,
                parsed[3],
            )
        attempts += 1
        if base is None:
            continue
        factors = generic_factors(base, prime)
        if arguments.stratum == "generic":
            if any(value == 0 for value in factors.values()):
                continue
        else:
            stratum_zero_factors = {
                "Q": {"Q"},
                "J": {"J"},
                "K": {"K"},
                "H": {"H"},
                "KH": {"K", "H"},
                "QJH": {"Q", "J", "H"},
                "JH": {"J", "H"},
                "JK": {"J", "K"},
                "a2": set(),
                "discriminant": set(),
            }[arguments.stratum]
            if any(
                factors[factor] != 0 for factor in stratum_zero_factors
            ) or any(
                value == 0
                for factor, value in factors.items()
                if factor not in stratum_zero_factors
            ):
                continue
        mu3_coefficients = evaluate_mu3_coefficients(
            parsed[3],
            base,
            prime,
        )
        if arguments.stratum == "a2":
            constant, linear, quadratic = mu3_coefficients
            assert quadratic == 0
            roots = (
                (-constant * pow(linear, -1, prime)) % prime,
            ) if linear else None
        elif arguments.stratum == "discriminant":
            constant, linear, quadratic = mu3_coefficients
            assert quadratic
            assert (
                linear * linear - 4 * quadratic * constant
            ) % prime == 0
            roots = (
                (-linear * pow(2 * quadratic, -1, prime)) % prime,
            )
        else:
            roots = quadratic_roots(mu3_coefficients, prime)
        if roots is None:
            continue
        sample_indices = []
        for root in roots:
            sampled_orders = (
                (4, 5, 6, 7, 8)
                if arguments.pivot_scout
                else (4, 5, 6, 7)
            )
            moments = {
                str(order): serialize_bivariate(
                    evaluate_base_to_bivariate(
                        parsed[order],
                        base,
                        root,
                        prime,
                    ),
                    prime,
                )
                for order in sampled_orders
            }
            sample_indices.append(len(samples))
            samples.append(
                {
                    "base": base,
                    "s3": root,
                    "moments": moments,
                }
            )
        pair_metadata.append(
            {
                **(
                    {"variable_value": base[arguments.variable]}
                    if generator is None
                    else {"random_pair_index": len(pair_metadata)}
                ),
                "base": base,
                "s3_roots": list(roots),
                "sample_indices": sample_indices,
                "generic_factors": factors,
            }
        )
    assert len(pair_metadata) == arguments.sample_count, (
        len(pair_metadata),
        arguments.sample_count,
        attempts,
        maximum_attempts,
    )
    if arguments.stratum != "generic":
        direct_records = run_direct_samples(
            singular,
            prime,
            samples,
            arguments.timeout,
        )
        direct_candidates = []
        for record in direct_records:
            if record["through_mu7_length"] <= 0:
                continue
            sample = samples[record["sample_index"]]
            candidate_record = {
                **record,
                "base": sample["base"],
                "s3": sample["s3"],
            }
            if record["through_mu7_length"] == 1:
                fiber = unique_fiber_point(
                    record["through_mu7_standard_basis"],
                    prime,
                )
                assert fiber is not None
                point = full_parameter_point(
                    sample["base"],
                    sample["s3"],
                    fiber,
                    prime,
                )
                moment_values = {}
                first_nonzero = None
                for order in CORRECTED_ORDERS:
                    moment_values[str(order)] = evaluate_moment_at_point(
                        order,
                        point,
                        prime,
                    )
                    if order >= 8 and moment_values[str(order)]:
                        first_nonzero = order
                        break
                assert all(
                    moment_values[str(order)] == 0
                    for order in (2, 3, 4, 5, 6, 7)
                )
                candidate_record["higher_moment_evaluation"] = {
                    "full_parameter_point": point,
                    "moment_values": moment_values,
                    "first_nonzero_order": first_nonzero,
                    "all_corrected_moments_zero": first_nonzero is None,
                }
            direct_candidates.append(candidate_record)
        length_distribution: dict[str, int] = defaultdict(int)
        for record in direct_records:
            length_distribution[str(record["mu4_mu5_length"])] += 1
        generic_length = max(
            record["mu4_mu5_length"] for record in direct_records
        )
        length_drop_samples = []
        for record in direct_records:
            if record["mu4_mu5_length"] >= generic_length:
                continue
            sample = samples[record["sample_index"]]
            coefficients = evaluate_mu3_coefficients(
                parsed[3],
                sample["base"],
                prime,
            )
            constant, linear, quadratic = coefficients
            length_drop_samples.append(
                {
                    "sample_index": record["sample_index"],
                    "base": sample["base"],
                    "s3": sample["s3"],
                    "mu4_mu5_length": record["mu4_mu5_length"],
                    "mu3_coefficients_ascending": coefficients,
                    "mu3_discriminant": (
                        linear * linear - 4 * quadratic * constant
                    ) % prime,
                }
            )
        payload = {
            "format": (
                "two-pair-sic-bidegree33-t0-specialized-direct-scout-v1"
            ),
            "status": "bounded modular specialized direct ideals; not a proof",
            "prime": prime,
            "stratum": arguments.stratum,
            "random_seed": arguments.random_seed,
            "attempted_base_count": attempts,
            "accepted_paired_base_count": len(pair_metadata),
            "evaluated_mu3_root_count": len(samples),
            "mu4_mu5_length_distribution": dict(length_distribution),
            "generic_mu4_mu5_length": generic_length,
            "length_drop_sample_count": len(length_drop_samples),
            "length_drop_samples": length_drop_samples,
            "common_through_mu7_point_count": len(direct_candidates),
            "common_through_mu7_points": direct_candidates,
            "reproduction_command": " ".join(sys.argv),
        }
        if arguments.retain_pairs:
            payload["direct_samples"] = [
                {
                    "sample_index": record["sample_index"],
                    "base": samples[record["sample_index"]]["base"],
                    "s3": samples[record["sample_index"]]["s3"],
                    "mu4_mu5_length": record["mu4_mu5_length"],
                    "through_mu7_length": record["through_mu7_length"],
                }
                for record in direct_records
            ]
        if arguments.output is not None:
            output = arguments.output
            if not output.is_absolute():
                output = ROOT / output
            output.parent.mkdir(parents=True, exist_ok=True)
            output.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        print(json.dumps(payload, indent=2, sort_keys=True))
        return
    values, _, border_values, rank_values = run_samples(
        singular,
        prime,
        samples,
        arguments.timeout,
        arguments.pivot_scout,
    )
    pairs = []
    common_pencil_zero_samples = []
    selected_pivot_complement_samples = []
    rank_at_most_four_samples = []
    candidate_sample_indices = []
    for metadata in pair_metadata:
        left_index, right_index = metadata["sample_indices"]
        left_root, right_root = metadata["s3_roots"]
        pencil_coordinates = []
        for coefficient in range(7):
            constant, slope = paired_coordinates(
                left_root,
                right_root,
                values[left_index][coefficient],
                values[right_index][coefficient],
                prime,
            )
            pencil_coordinates.append(
                {"constant": constant, "s3": slope}
            )
        bordered_coordinates = {}
        if arguments.pivot_scout:
            for name in BORDER_NAMES:
                constant, slope = paired_coordinates(
                    left_root,
                    right_root,
                    border_values[left_index][name],
                    border_values[right_index][name],
                    prime,
                )
                bordered_coordinates[name] = {
                    "constant": constant,
                    "s3": slope,
                }
        for sample_index, root in zip(
            (left_index, right_index),
            (left_root, right_root),
            strict=True,
        ):
            if arguments.pivot_scout:
                rank_record = rank_values[sample_index]
                if (
                    border_values[sample_index]["pivot_A"] == 0
                    and border_values[sample_index]["pivot_B"] == 0
                ):
                    selected_pivot_complement_samples.append(
                        {
                            "base": metadata["base"],
                            "s3": root,
                            "sample_index": sample_index,
                            "ranks": rank_record,
                            "nonzero_restoring_borders": [
                                name
                                for name in BORDER_NAMES[2:]
                                if border_values[sample_index][name] != 0
                            ],
                        }
                    )
                if rank_record["M6_M7"] <= 4:
                    rank_at_most_four_samples.append(
                        {
                            "base": metadata["base"],
                            "s3": root,
                            "sample_index": sample_index,
                            "ranks": rank_record,
                        }
                    )
            if all(value == 0 for value in values[sample_index]):
                candidate_sample_indices.append(sample_index)
                common_pencil_zero_samples.append(
                    {
                        "base": metadata["base"],
                        "s3": root,
                        "sample_index": sample_index,
                        "pencil_coefficients": list(values[sample_index]),
                    }
                )
        pairs.append(
            {
                **metadata,
                "c0_mod_mu3": pencil_coordinates[0],
                "c1_mod_mu3": pencil_coordinates[1],
                "pencil_mod_mu3": {
                    f"c{coefficient}": coordinates
                    for coefficient, coordinates in enumerate(
                        pencil_coordinates
                    )
                },
                **(
                    {"bordered_mod_mu3": bordered_coordinates}
                    if arguments.pivot_scout
                    else {}
                ),
            }
        )
    if candidate_sample_indices:
        pivot_samples = []
        for sample_index in candidate_sample_indices:
            sample = samples[sample_index]
            moments = dict(sample["moments"])
            moments["8"] = serialize_bivariate(
                evaluate_base_to_bivariate(
                    parsed[8],
                    sample["base"],
                    sample["s3"],
                    prime,
                ),
                prime,
            )
            pivot_samples.append({**sample, "moments": moments})
        pivot_values, candidate_matrices, _, _ = run_samples(
            singular,
            prime,
            pivot_samples,
            arguments.timeout,
        )
        assert all(
            all(value == 0 for value in pivot_values[index])
            for index in range(len(pivot_samples))
        )
        assert set(candidate_matrices) == set(range(len(pivot_samples)))
        for index, candidate in enumerate(common_pencil_zero_samples):
            candidate["pivot_signature"] = pivot_signature(
                candidate_matrices[index],
                prime,
            )
    candidate_replays = replay_common_pencil_zeros(
        singular,
        prime,
        samples,
        candidate_sample_indices,
        arguments.timeout,
    )
    assert [
        candidate["sample_index"]
        for candidate in common_pencil_zero_samples
    ] == [
        replay["sample_index"] for replay in candidate_replays
    ]
    for candidate, replay in zip(
        common_pencil_zero_samples,
        candidate_replays,
        strict=True,
    ):
        fiber = unique_fiber_point(replay["standard_basis"], prime)
        if fiber is None:
            replay["higher_moment_evaluation"] = {
                "status": "not evaluated: common quotient is not a unique point"
            }
            continue
        point = full_parameter_point(
            candidate["base"],
            candidate["s3"],
            fiber,
            prime,
        )
        moment_values = {}
        first_nonzero = None
        for order in CORRECTED_ORDERS:
            moment_values[str(order)] = evaluate_moment_at_point(
                order,
                point,
                prime,
            )
            if order >= 8 and moment_values[str(order)]:
                first_nonzero = order
                break
        assert all(
            moment_values[str(order)] == 0
            for order in (2, 3, 4, 5, 6, 7)
        ), moment_values
        replay["higher_moment_evaluation"] = {
            "status": "evaluated directly from the corrected moment formula",
            "full_parameter_point": point,
            "moment_values": moment_values,
            "evaluated_orders": [
                int(order) for order in moment_values
            ],
            "first_nonzero_order": first_nonzero,
            "all_corrected_moments_zero": first_nonzero is None,
        }
    payload = {
        "format": "two-pair-sic-bidegree33-t0-fitting-samples-v2",
        "status": "bounded modular determinant-pencil samples; not a proof",
        "prime": prime,
        "sampling_mode": "random" if generator is not None else "directional",
        "pivot_scout": arguments.pivot_scout,
        "sample_count": len(pairs),
        "attempted_base_count": attempts,
        "common_pencil_zero_sample_count": len(
            common_pencil_zero_samples
        ),
        "common_pencil_zero_samples": common_pencil_zero_samples,
        "selected_pivot_complement_sample_count": len(
            selected_pivot_complement_samples
        ),
        "selected_pivot_complement_samples": (
            selected_pivot_complement_samples
        ),
        "rank_at_most_four_sample_count": len(rank_at_most_four_samples),
        "rank_at_most_four_samples": rank_at_most_four_samples,
        "candidate_replays": candidate_replays,
        "reproduction_command": " ".join(sys.argv),
    }
    if generator is None or arguments.retain_pairs:
        payload["pairs"] = pairs
    if generator is None:
        payload.update(
            {
                "varied_base_variable": arguments.variable,
                "fixed_base_values": {
                    key: value % prime
                    for key, value in base_template.items()
                    if key != arguments.variable
                },
                "candidate_interval": [arguments.start, candidate],
            }
        )
    else:
        payload["random_seed"] = arguments.random_seed
    if arguments.training_count:
        assert generator is None, (
            "rational reconstruction is only available in directional mode"
        )
        assert 1 < arguments.training_count < len(pairs)
        reconstruction = reconstruct_coordinates(
            pairs,
            arguments.training_count,
            prime,
        )
        payload["rational_reconstruction"] = reconstruction
        payload["denominator_model"] = known_denominator_model(
            parsed[3],
            arguments.variable,
            prime,
            reconstruction,
            base_template,
        )
        if arguments.pivot_scout:
            payload["bordered_rational_reconstruction"] = (
                reconstruct_bordered_coordinates(
                    pairs,
                    arguments.training_count,
                    prime,
                )
            )
    if arguments.output is not None:
        output = arguments.output
        if not output.is_absolute():
            output = ROOT / output
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
