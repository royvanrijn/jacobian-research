#!/usr/bin/env python3
"""Inspect specialized t0-open Groebner denominators over finite fields."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from verify_two_pair_sic_bidegree33_boundary_generic_quotient import (
    substitute,
)
from verify_two_pair_sic_bidegree33_corrected_boundary import (
    t0_open_localized_export,
)


ROOT = Path(__file__).resolve().parents[1]
STRATA = (
    "Q",
    "J",
    "K",
    "H",
    "KH",
    "QJH",
    "JH",
    "JK",
    "a2",
)


def square_root_mod(value: int, prime: int) -> int:
    assert prime % 4 == 3
    root = pow(value % prime, (prime + 1) // 4, prime)
    assert root * root % prime == value % prime
    return root


def stratum_data(
    stratum: str,
    prime: int,
) -> tuple[tuple[str, ...], tuple[tuple[str, str], ...]]:
    if stratum == "Q":
        return (
            ("s1", "s3", "t1", "t2", "u"),
            (("s2", "(s1^2*u-(13/3)*u)"),),
        )
    if stratum == "K":
        return (
            ("s1", "s3", "t1", "t2", "u"),
            (("s2", "(s1^2*u-(901/351)*u)"),),
        )
    if stratum == "KH":
        return (
            ("s1", "s3", "t2", "u"),
            (
                ("s2", "(s1^2*u-(901/351)*u)"),
                ("t1", "(s1*u)"),
            ),
        )
    if stratum == "J":
        alpha = square_root_mod(-30420, prime)
        return (
            ("s1", "s3", "ell", "t2", "u"),
            (
                ("s2", f"(s1^2*u-(({alpha}*ell+274*u)/99))"),
                ("t1", "(s1*u-ell)"),
            ),
        )
    if stratum == "QJH":
        alpha = square_root_mod(-30420, prime)
        return (
            ("s1", "s3", "t2", "u"),
            (
                ("s2", "(s1^2*u-(13/3)*u)"),
                ("t1", f"(s1*u-(155*u/{alpha}))"),
            ),
        )
    if stratum == "JH":
        return (
            ("s1", "s3", "t2", "u"),
            (
                ("s2", "(s1^2*u-(274/99)*u)"),
                ("t1", "(s1*u)"),
            ),
        )
    if stratum == "JK":
        alpha = square_root_mod(-30420, prime)
        return (
            ("s1", "s3", "t2", "u"),
            (
                ("s2", "(s1^2*u-(901/351)*u)"),
                (
                    "t1",
                    f"(s1*u-((99*(901/351)*u-274*u)/{alpha}))",
                ),
            ),
        )
    if stratum == "a2":
        return (
            ("s1", "s3", "s2", "t1", "u"),
            (
                (
                    "t2",
                    "(2*s1*t1*u-s2*u-(8/9)*u^2)",
                ),
            ),
        )
    assert stratum == "H"
    return (
        ("s1", "s3", "t2", "u", "r"),
        (
            (
                "s2",
                (
                    "(s1^2*u-((-155*u/33"
                    "+r*(-775*r*u/(1287*r^2+40560))+13*u)/3))"
                ),
            ),
            ("t1", "(s1*u-(-775*r*u/(1287*r^2+40560)))"),
        ),
    )


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=43)
    parser.add_argument("--stratum", choices=STRATA, required=True)
    parser.add_argument("--timeout", type=int, default=20)
    parser.add_argument(
        "--include-basis",
        action="store_true",
        help="include the full lifted Groebner generators in the artifact",
    )
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    assert arguments.prime == 0 or (
        arguments.prime > 7 and arguments.prime % 4 == 3
    )
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    export = t0_open_localized_export(
        singular,
        (2, 3, 4, 5),
        arguments.prime,
        arguments.timeout,
    )
    coefficient_variables, replacements = stratum_data(
        arguments.stratum,
        arguments.prime,
    )
    p4 = substitute(export["polynomials"][1], replacements)
    p5 = substitute(export["polynomials"][2], replacements)
    full_generator_print = (
        """
  print(
    "FULL_GENERATOR "+string(generatorIndex)+" "
    +string(G[generatorIndex])
  );
"""
        if arguments.include_basis
        else ""
    )
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring specialized=(
  {arguments.prime},{",".join(coefficient_variables)}
),(s6,s5),dp;
option(redSB);
poly p4={p4};
poly p5={p5};
ideal inputIdeal=p4,p5;
matrix transformation;
ideal G=liftstd(inputIdeal,transformation);
proc polynomialLcm(poly left, poly right)
{{
  return(left*right/gcd(left,right));
}}
poly commonDenominator=1;
poly transformationDenominator=1;
poly leadingCoefficientLcm=1;
poly cursor;
number coefficient;
int generatorIndex;
int transformationRow;
int transformationColumn;
for(generatorIndex=1;generatorIndex<=size(G);generatorIndex++)
{{
  cursor=G[generatorIndex];
  leadingCoefficientLcm=polynomialLcm(
    leadingCoefficientLcm,
    numerator(leadcoef(cursor))
  );
  print(
    "GENERATOR "+string(generatorIndex)+" "
    +string(leadmonom(cursor))
  );
{full_generator_print}
  while(cursor!=0)
  {{
    coefficient=leadcoef(cursor);
    commonDenominator=polynomialLcm(
      commonDenominator,
      denominator(coefficient)
    );
    cursor=cursor-lead(cursor);
  }}
}}
for(
  transformationRow=1;
  transformationRow<=nrows(transformation);
  transformationRow++
)
{{
  for(
    transformationColumn=1;
    transformationColumn<=ncols(transformation);
    transformationColumn++
  )
  {{
    cursor=transformation[transformationRow,transformationColumn];
    while(cursor!=0)
    {{
      coefficient=leadcoef(cursor);
      transformationDenominator=polynomialLcm(
        transformationDenominator,
        denominator(coefficient)
      );
      cursor=cursor-lead(cursor);
    }}
  }}
}}
print("META "+string(size(G))+" "+string(vdim(G)));
print("DENOMINATOR "+string(commonDenominator));
print("TRANSFORMATION_DENOMINATOR "+string(transformationDenominator));
print("LEADING_COEFFICIENT_LCM "+string(leadingCoefficientLcm));
if(commonDenominator==1)
{{
  print("FACTOR 1 1");
}}
else
{{
  list denominatorFactorization=factorize(commonDenominator);
  ideal denominatorFactors=denominatorFactorization[1];
  intvec denominatorPowers=denominatorFactorization[2];
  int factorIndex;
  for(factorIndex=1;factorIndex<=size(denominatorFactors);factorIndex++)
  {{
    print(
      "FACTOR "+string(denominatorFactors[factorIndex])+" "
      +string(denominatorPowers[factorIndex])
    );
  }}
}}
if(leadingCoefficientLcm==1)
{{
  print("LEADING_FACTOR 1 1");
}}
else
{{
  list leadingFactorization=factorize(leadingCoefficientLcm);
  ideal leadingFactors=leadingFactorization[1];
  intvec leadingPowers=leadingFactorization[2];
  int leadingFactorIndex;
  for(
    leadingFactorIndex=1;
    leadingFactorIndex<=size(leadingFactors);
    leadingFactorIndex++
  )
  {{
    print(
      "LEADING_FACTOR "+string(leadingFactors[leadingFactorIndex])+" "
      +string(leadingPowers[leadingFactorIndex])
    );
  }}
}}
if(transformationDenominator==1)
{{
  print("TRANSFORMATION_FACTOR 1 1");
}}
else
{{
  list transformationFactorization=factorize(transformationDenominator);
  ideal transformationFactors=transformationFactorization[1];
  intvec transformationPowers=transformationFactorization[2];
  int transformationFactorIndex;
  for(
    transformationFactorIndex=1;
    transformationFactorIndex<=size(transformationFactors);
    transformationFactorIndex++
  )
  {{
    print(
      "TRANSFORMATION_FACTOR "
      +string(transformationFactors[transformationFactorIndex])+" "
      +string(transformationPowers[transformationFactorIndex])
    );
  }}
}}
""",
        text=True,
        capture_output=True,
        check=True,
        timeout=arguments.timeout,
    )
    assert "\n   ? " not in completed.stdout, (
        completed.stdout[-4000:],
        completed.stderr[-2000:],
    )
    meta = re.search(r"(?m)^META (\d+) (\d+)$", completed.stdout)
    denominator = re.search(r"(?m)^DENOMINATOR (.*)$", completed.stdout)
    transformation_denominator = re.search(
        r"(?m)^TRANSFORMATION_DENOMINATOR (.*)$",
        completed.stdout,
    )
    leading_coefficient_lcm = re.search(
        r"(?m)^LEADING_COEFFICIENT_LCM (.*)$",
        completed.stdout,
    )
    assert (
        meta is not None
        and denominator is not None
        and transformation_denominator is not None
        and leading_coefficient_lcm is not None
    )
    payload = {
        "format": "two-pair-sic-bidegree33-t0-stratum-leading-v1",
        "status": (
            "exact characteristic-zero rational-function-field "
            "Groebner border"
            if arguments.prime == 0
            else (
                "bounded exact finite-field rational-function-field "
                "Groebner border; not a characteristic-zero certificate"
            )
        ),
        "prime": arguments.prime,
        "stratum": arguments.stratum,
        "coefficient_variables": list(coefficient_variables),
        "groebner_basis_size": int(meta.group(1)),
        "quotient_length": int(meta.group(2)),
        "leading_monomials": [
            monomial
            for _, monomial in re.findall(
                r"(?m)^GENERATOR (\d+) (.*)$",
                completed.stdout,
            )
        ],
        "common_denominator": denominator.group(1),
        "common_denominator_terms": (
            denominator.group(1).count("+")
            + denominator.group(1).count("-")
            + 1
        ),
        "common_denominator_factors": [
            {
                "factor": factor,
                "multiplicity": int(multiplicity),
            }
            for factor, multiplicity in re.findall(
                r"(?m)^FACTOR (.*) (\d+)$",
                completed.stdout,
            )
        ],
        "transformation_denominator": transformation_denominator.group(1),
        "transformation_denominator_factors": [
            {
                "factor": factor,
                "multiplicity": int(multiplicity),
            }
            for factor, multiplicity in re.findall(
                r"(?m)^TRANSFORMATION_FACTOR (.*) (\d+)$",
                completed.stdout,
            )
        ],
        "leading_coefficient_lcm": leading_coefficient_lcm.group(1),
        "leading_coefficient_factors": [
            {
                "factor": factor,
                "multiplicity": int(multiplicity),
            }
            for factor, multiplicity in re.findall(
                r"(?m)^LEADING_FACTOR (.*) (\d+)$",
                completed.stdout,
            )
        ],
        "reproduction_command": " ".join(sys.argv),
    }
    if arguments.include_basis:
        full_generators = [
            generator
            for _, generator in re.findall(
                r"(?m)^FULL_GENERATOR (\d+) (.*)$",
                completed.stdout,
            )
        ]
        assert len(full_generators) == int(meta.group(1))
        payload["basis_polynomials"] = full_generators
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
