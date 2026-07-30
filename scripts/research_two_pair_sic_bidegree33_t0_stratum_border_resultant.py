#!/usr/bin/env python3
"""Project a specialized quotient border against the mu3 base equation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

from research_two_pair_sic_bidegree33_t0_stratum_leading import (
    ROOT,
    stratum_data,
)
from verify_two_pair_sic_bidegree33_boundary_generic_quotient import (
    substitute,
)
from verify_two_pair_sic_bidegree33_corrected_boundary import (
    t0_open_localized_export,
)


STRATA = ("Q", "J", "JK", "QJH", "JH")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=43)
    parser.add_argument("--stratum", choices=STRATA, required=True)
    parser.add_argument("--timeout", type=int, default=20)
    parser.add_argument("--output", type=Path, default=None)
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")
    coefficient_variables, replacements = stratum_data(
        arguments.stratum,
        arguments.prime,
    )
    assert "s3" in coefficient_variables
    base_variables = tuple(
        variable for variable in coefficient_variables if variable != "s3"
    )
    export = t0_open_localized_export(
        singular,
        (2, 3),
        arguments.prime,
        arguments.timeout,
    )
    p3 = substitute(export["polynomials"][0], replacements)
    leading_path = (
        ROOT
        / "artifacts"
        / "generated-results"
        / (
            "two_pair_sic_bidegree33_t0_stratum_"
            f"{arguments.stratum}_leading_"
            + (
                "exact.json"
                if arguments.prime == 0
                else f"mod{arguments.prime}.json"
            )
        )
    )
    leading_payload = json.loads(leading_path.read_text(encoding="utf-8"))
    border = leading_payload["leading_coefficient_lcm"]
    completed = subprocess.run(
        [singular, "-q"],
        input=f"""
ring projected={arguments.prime},(
  {",".join((*base_variables, "s3"))}
),dp;
poly p3={p3};
poly border={border};
poly projectedBorder=resultant(p3,border,s3);
print(
  "META "+string(size(projectedBorder))+" "
  +string(deg(projectedBorder))
);
print("RESULTANT "+string(projectedBorder));
list projectedFactorization=factorize(projectedBorder);
ideal projectedFactors=projectedFactorization[1];
intvec projectedPowers=projectedFactorization[2];
int factorIndex;
int residualTermCount=0;
poly residualFactor=1;
for(factorIndex=1;factorIndex<=size(projectedFactors);factorIndex++)
{{
  print(
    "FACTOR "+string(projectedFactors[factorIndex])+" "
      +string(projectedPowers[factorIndex])
  );
  if(size(projectedFactors[factorIndex])>residualTermCount)
  {{
    residualTermCount=size(projectedFactors[factorIndex]);
    residualFactor=projectedFactors[factorIndex];
  }}
}}
poly p3Quadratic=subst(diff(diff(p3,s3),s3),s3,0)/2;
poly p3Linear=subst(diff(p3,s3),s3,0);
poly p3Constant=subst(p3,s3,0);
poly borderQuadratic=subst(diff(diff(border,s3),s3),s3,0)/2;
poly borderLinear=subst(diff(border,s3),s3,0);
poly borderConstant=subst(border,s3,0);
poly pivotA=borderQuadratic*p3Linear-p3Quadratic*borderLinear;
poly pivotB=borderQuadratic*p3Constant-p3Quadratic*borderConstant;
print(
  "PIVOT_META "+string(size(pivotA))+" "+string(deg(pivotA))+" "
  +string(size(pivotB))+" "+string(deg(pivotB))
);
print("PIVOT_A "+string(pivotA));
print("PIVOT_B "+string(pivotB));
print("PIVOT_GCD_A "+string(gcd(residualFactor,pivotA)));
print("PIVOT_GCD_B "+string(gcd(residualFactor,pivotB)));
for(factorIndex=1;factorIndex<=size(projectedFactors);factorIndex++)
{{
  print(
    "FACTOR_META "+string(factorIndex)+" "
    +string(size(projectedFactors[factorIndex]))+" "
    +string(deg(projectedFactors[factorIndex]))+" "
    +string(gcd(projectedFactors[factorIndex],pivotA))+" "
    +string(gcd(projectedFactors[factorIndex],pivotB))
  );
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
    meta = re.search(r"(?m)^META (\d+) (-?\d+)$", completed.stdout)
    resultant = re.search(r"(?m)^RESULTANT (.*)$", completed.stdout)
    assert meta is not None and resultant is not None
    factors = [
        {"factor": factor, "multiplicity": int(multiplicity)}
        for factor, multiplicity in re.findall(
            r"(?m)^FACTOR (.*) (\d+)$",
            completed.stdout,
        )
    ]
    pivot_meta = re.search(
        r"(?m)^PIVOT_META (\d+) (-?\d+) (\d+) (-?\d+)$",
        completed.stdout,
    )
    pivot_a = re.search(r"(?m)^PIVOT_A (.*)$", completed.stdout)
    pivot_b = re.search(r"(?m)^PIVOT_B (.*)$", completed.stdout)
    pivot_gcd_a = re.search(r"(?m)^PIVOT_GCD_A (.*)$", completed.stdout)
    pivot_gcd_b = re.search(r"(?m)^PIVOT_GCD_B (.*)$", completed.stdout)
    assert all(
        value is not None
        for value in (
            pivot_meta,
            pivot_a,
            pivot_b,
            pivot_gcd_a,
            pivot_gcd_b,
        )
    )
    factor_metadata = {
        int(index): {
            "term_count": int(term_count),
            "total_degree": int(total_degree),
            "gcd_A": gcd_a,
            "gcd_B": gcd_b,
        }
        for index, term_count, total_degree, gcd_a, gcd_b in re.findall(
            r"(?m)^FACTOR_META (\d+) (\d+) (-?\d+) ([^ ]+) ([^ ]+)$",
            completed.stdout,
        )
    }
    assert set(factor_metadata) == set(range(1, len(factors) + 1))
    for index, factor in enumerate(factors, start=1):
        factor.update(factor_metadata[index])
    payload = {
        "format": (
            "two-pair-sic-bidegree33-t0-stratum-border-resultant-v1"
        ),
        "status": (
            "exact characteristic-zero factored border resultant"
            if arguments.prime == 0
            else (
                "bounded exact finite-field resultant; "
                "not a characteristic-zero certificate"
            )
        ),
        "prime": arguments.prime,
        "stratum": arguments.stratum,
        "base_variables": list(base_variables),
        "resultant_term_count": int(meta.group(1)),
        "resultant_total_degree": int(meta.group(2)),
        "resultant": resultant.group(1),
        "factors": factors,
        "linear_subresultant": {
            "A": pivot_a.group(1),
            "A_term_count": int(pivot_meta.group(1)),
            "A_total_degree": int(pivot_meta.group(2)),
            "B": pivot_b.group(1),
            "B_term_count": int(pivot_meta.group(3)),
            "B_total_degree": int(pivot_meta.group(4)),
            "gcd_residual_factor_A": pivot_gcd_a.group(1),
            "gcd_residual_factor_B": pivot_gcd_b.group(1),
            "dense_pivot": "s3=-B/A",
        },
        "reproduction_command": " ".join(sys.argv),
    }
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
