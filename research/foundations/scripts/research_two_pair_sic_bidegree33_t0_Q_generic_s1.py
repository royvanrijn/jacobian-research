#!/usr/bin/env python3
"""Test the residual Q system over GF(p)(s1).

The eight-variable direct solve treats ``s1`` as an ordinary variable and was
too large for the available memory.  This driver instead makes ``s1`` a
transcendental coefficient parameter.  Thus Singular computes in

    GF(p)(s1)[s6,s5,s3,T,v,vinv,ell].

The default ``std`` mode only decides whether the generic fibre is empty.  The
more expensive ``liftstd`` mode also computes a transformation matrix; when
the standard basis is the unit ideal, clearing its coefficient denominators
produces a univariate polynomial whose roots contain every exceptional
``s1``-fibre.  Both modes are modular research calculations, not
characteristic-zero certificates.
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
import time

from research_two_pair_sic_bidegree33_t0_Q_flint import (
    ROOT,
    configure_coefficient_ring,
    direct_residual_system,
    flint_input_data,
    parse_moments,
)


DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_t0_stratum_Q_generic_s1_mod1000003.json"
)


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=1000003)
    parser.add_argument("--timeout", type=int, default=300)
    parser.add_argument(
        "--solver-seconds",
        type=int,
        default=0,
        help="Singular deadline in seconds; 0 means no deadline",
    )
    parser.add_argument(
        "--mode",
        choices=("std", "liftstd"),
        default="std",
    )
    parser.add_argument(
        "--algorithm",
        choices=("std", "slimgb"),
        default="slimgb",
        help="standard-basis engine (std mode only)",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--singular-output",
        type=Path,
        default=None,
        help="optional path for the complete Singular transcript",
    )
    return parser.parse_args()


def denominator_program() -> str:
    return r"""
proc polynomialLcm(poly left, poly right)
{
  return(left*right/gcd(left,right));
}
poly transformationDenominator=1;
poly cursor;
number coefficient;
int transformationRow;
int transformationColumn;
for(
  transformationRow=1;
  transformationRow<=nrows(transformation);
  transformationRow++
)
{
  for(
    transformationColumn=1;
    transformationColumn<=ncols(transformation);
    transformationColumn++
  )
  {
    cursor=transformation[transformationRow,transformationColumn];
    while(cursor!=0)
    {
      coefficient=leadcoef(cursor);
      transformationDenominator=polynomialLcm(
        transformationDenominator,
        denominator(coefficient)
      );
      cursor=cursor-lead(cursor);
    }
  }
}
print(
  "QGENERIC_TRANSFORMATION_SHAPE "
  +string(nrows(transformation))+" "
  +string(ncols(transformation))
);
print(
  "QGENERIC_TRANSFORMATION_DENOMINATOR "
  +string(transformationDenominator)
);
if(transformationDenominator==1)
{
  print("QGENERIC_FACTOR 1 1");
}
else
{
  list denominatorFactorization=factorize(transformationDenominator);
  ideal denominatorFactors=denominatorFactorization[1];
  intvec denominatorPowers=denominatorFactorization[2];
  int factorIndex;
  for(factorIndex=1;factorIndex<=size(denominatorFactors);factorIndex++)
  {
    print(
      "QGENERIC_FACTOR "+string(denominatorFactors[factorIndex])+" "
      +string(denominatorPowers[factorIndex])
    );
  }
}
"""


def singular_program(
    system: dict[str, object],
    prime: int,
    mode: str,
    algorithm: str,
) -> str:
    variables = [
        variable
        for variable in system["ordinary_variables"]
        if variable != "s1"
    ]
    ideal = ",\n".join(system["polynomials"])
    if mode == "liftstd":
        computation = """
matrix transformation;
ideal G=liftstd(I,transformation);
"""
        after = denominator_program()
    else:
        computation = f"ideal G={algorithm}(I);\n"
        after = ""
    return f"""
ring generic=({prime},s1),({",".join(variables)}),dp;
option(redSB);
ideal I=
{ideal};
timer=1;
{computation}
int elapsed=timer;
int isUnit=(reduce(1,G)==0);
print(
  "QGENERIC_META "+string(size(G))+" "+string(dim(G))+" "
  +string(isUnit)+" "+string(elapsed)
);
{after}
"""


def parse_result(
    stdout: str,
    stderr: str,
    returncode: int | None,
    elapsed: float,
    timed_out: bool,
) -> dict[str, object]:
    marker = re.search(
        r"(?m)^QGENERIC_META (-?\d+) (-?\d+) ([01]) (\d+)$",
        stdout,
    )
    if timed_out:
        status = "timeout"
    elif returncode != 0:
        status = f"solver-error-{returncode}"
    elif marker is None:
        status = "missing-result-marker"
    elif marker.group(3) == "1":
        status = "unit"
    else:
        status = "nonunit"
    result: dict[str, object] = {
        "status": status,
        "returncode": returncode,
        "seconds": round(elapsed, 3),
        "stdout_tail": stdout[-8000:],
        "stderr_tail": stderr[-8000:],
    }
    if marker is not None:
        result["basis_size"] = int(marker.group(1))
        result["dimension"] = int(marker.group(2))
        result["unit_ideal"] = marker.group(3) == "1"
        result["singular_timer_ticks"] = int(marker.group(4))
    denominator = re.search(
        r"(?m)^QGENERIC_TRANSFORMATION_DENOMINATOR (.*)$",
        stdout,
    )
    if denominator is not None:
        value = denominator.group(1)
        result["transformation_denominator"] = value
        result["transformation_denominator_sha256"] = hashlib.sha256(
            value.encode()
        ).hexdigest()
        result["transformation_denominator_degree"] = (
            max(
                (
                    int(exponent or "1")
                    for exponent in re.findall(
                        r"s1(?:\^(\d+))?",
                        value,
                    )
                ),
                default=0,
            )
        )
        result["factors"] = [
            {
                "polynomial": factor,
                "multiplicity": int(multiplicity),
            }
            for factor, multiplicity in re.findall(
                r"(?m)^QGENERIC_FACTOR (.*) (\d+)$",
                stdout,
            )
        ]
    return result


def main() -> None:
    arguments = parse_arguments()
    if arguments.prime in (0, 2, 3, 5, 7, 13):
        raise ValueError("choose a prime avoiding the displayed denominators")
    singular = shutil.which("Singular")
    if singular is None:
        raise RuntimeError("Singular is required")

    started = time.monotonic()
    configure_coefficient_ring(arguments.prime)
    arithmetic, a_value, b_value, input_metadata = flint_input_data()
    moments, moment_profiles = parse_moments(singular, arguments.timeout)
    system = direct_residual_system(
        arithmetic,
        a_value,
        b_value,
        moments,
    )
    program = singular_program(
        system,
        arguments.prime,
        arguments.mode,
        arguments.algorithm,
    )
    print(
        "QGENERIC_START "
        f"prime={arguments.prime} mode={arguments.mode} "
        f"algorithm={arguments.algorithm} input_bytes={len(program)}",
        file=sys.stderr,
        flush=True,
    )

    solver_started = time.monotonic()
    timed_out = False
    try:
        completed = subprocess.run(
            [singular, "-q"],
            input=program,
            text=True,
            capture_output=True,
            check=False,
            timeout=(arguments.solver_seconds or None),
        )
        stdout = completed.stdout
        stderr = completed.stderr
        returncode = completed.returncode
    except subprocess.TimeoutExpired as error:
        timed_out = True
        stdout = (
            error.stdout.decode()
            if isinstance(error.stdout, bytes)
            else (error.stdout or "")
        )
        stderr = (
            error.stderr.decode()
            if isinstance(error.stderr, bytes)
            else (error.stderr or "")
        )
        returncode = None
    solver_elapsed = time.monotonic() - solver_started

    transcript = stdout + ("\nSTDERR\n" + stderr if stderr else "")
    if arguments.singular_output is not None:
        transcript_path = arguments.singular_output
        if not transcript_path.is_absolute():
            transcript_path = ROOT / transcript_path
        transcript_path.parent.mkdir(parents=True, exist_ok=True)
        transcript_path.write_text(transcript, encoding="utf-8")
    else:
        transcript_path = None

    payload = {
        "format": "two-pair-sic-bidegree33-t0-Q-generic-s1-v1",
        "status": (
            f"exact finite-field calculation over GF({arguments.prime})(s1); "
            "not a characteristic-zero certificate"
        ),
        "prime": arguments.prime,
        "coefficient_field": f"GF({arguments.prime})(s1)",
        "ordinary_variables": [
            variable
            for variable in system["ordinary_variables"]
            if variable != "s1"
        ],
        "mode": arguments.mode,
        "algorithm": (
            arguments.algorithm
            if arguments.mode == "std"
            else "Singular liftstd"
        ),
        "scope": system["scope"],
        "input": input_metadata,
        "moments": moment_profiles,
        "equation_profiles": system["equation_profiles"],
        "singular_program_sha256": hashlib.sha256(
            program.encode()
        ).hexdigest(),
        "singular_program_bytes": len(program),
        "solve": parse_result(
            stdout,
            stderr,
            returncode,
            solver_elapsed,
            timed_out,
        ),
        "singular_output": (
            None if transcript_path is None else str(transcript_path)
        ),
        "seconds": round(time.monotonic() - started, 3),
        "reproduction_command": " ".join(sys.argv),
        "interpretation": (
            "If solve.status is unit, the generic s1-fibre is empty modulo "
            "the displayed prime. In liftstd mode, every exceptional fibre "
            "is contained among the roots of transformation_denominator. "
            "Neither conclusion alone lifts to characteristic zero."
        ),
    }
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
