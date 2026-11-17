#!/usr/bin/env python3
"""Compute a bounded certificate for the frozen q-period Fitting stratum.

The expensive Weyl/function-field construction is frozen in

    artifacts/generated-results/degree_five_qper_15x16_presentation.sing.

The default run is a deliberately bounded audit at three fixed good primes.
It computes the sixteen maximal minors, removes their visible powers of ``a``
and ``a+1``, and saturates once by the product ``a(a+1)H``.  Saturation by
this product equals the former three sequential saturations.  The explicit
``--prime 0`` path uses Singular's modular saturation over Q.  Both ``sat``
and ``modSat`` already return a standard basis, so no second rational
standard-basis computation is needed.
"""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
PRESENTATION = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_five_qper_15x16_presentation.sing"
)
KERNEL_CHART = (
    "4*a^3*tau^2-24*a^3*tau-72*a^3"
    "+8*a^2*tau^2-54*a^2*tau-216*a^2"
    "+4*a*tau^2-30*a*tau-246*a-105"
)


def fitting_program(
    characteristic: int,
    method: str,
    stage: str,
    print_basis: bool = False,
) -> str:
    matrix_line = next(
        line
        for line in PRESENTATION.read_text().splitlines()
        if line.startswith("matrix M[15][16]=")
    )
    lines = [
        'LIB "elim.lib";',
    ]
    if method in {
        "modular",
        "modular-fast",
        "modular-rebuild",
    } and stage != "minors":
        lines.append('LIB "moddiq.lib";')
    lines.extend(
        [
            f"ring r={characteristic},(a,tau),dp;",
            "short=0;",
        ]
    )
    if method == "modular-fast" and stage != "minors":
        lines.extend(
            [
                "proc fastFinalStd("
                "string command,alias list args,def result)",
                "{",
                "  attrib(result,\"isSB\",1);",
                "  int i;",
                "  for(i=ncols(args[1]);i>0;i--)",
                "  {",
                "    if(reduce(args[1][i],result,1)!=0)",
                "    {",
                "      return(0);",
                "    }",
                "  }",
                "  return(system(\"verifyGB\",result));",
                "}",
                "proc fastModStd(def I)",
                "{",
                "  intvec saved_options=option(get);",
                "  option(redSB);",
                "  I=modular("
                "\"std\",list(I),"
                "Modstd::primeTest_std,"
                "Modstd::deleteUnluckyPrimes_std,"
                "Modstd::pTest_std,"
                "fastFinalStd);",
                "  attrib(I,\"isSB\",1);",
                "  option(set,saved_options);",
                "  return(I);",
                "}",
                "proc fastFinalSat("
                "string command,alias list args,def result)",
                "{",
                "  ideal H=result[1];",
                "  ideal boundary_power=args[2]^result[2];",
                "  ideal membership_tests=H*boundary_power;",
                "  int i;",
                "  for(i=ncols(membership_tests);i>0;i--)",
                "  {",
                "    if(reduce(membership_tests[i],args[1],5)!=0)",
                "    {",
                "      return(0);",
                "    }",
                "  }",
                "  return(system(\"verifyGB\",H));",
                "}",
                "proc fastModSat(def I,def J)",
                "{",
                "  I=fastModStd(I);",
                "  def result=modular("
                "\"sat_with_exp\",list(I,J),"
                "Moddiq::primeTest_sat,"
                "Moddiq::deleteUnluckyPrimes_sat,"
                "Moddiq::pTest_sat,"
                "fastFinalSat);",
                "  return(result);",
                "}",
            ]
        )
    if method == "modular-rebuild" and stage != "minors":
        lines.extend(
            [
                "proc rebuildModStd(def I)",
                "{",
                "  intvec saved_options=option(get);",
                "  option(redSB);",
                "  I=modular("
                "\"std\",list(I),"
                "Modstd::primeTest_std,"
                "Modstd::deleteUnluckyPrimes_std,"
                "Modstd::pTest_std);",
                "  attrib(I,\"isSB\",1);",
                "  option(set,saved_options);",
                "  return(I);",
                "}",
                "proc rebuildModSat(def I,def J)",
                "{",
                "  I=rebuildModStd(I);",
                "  def result=modular("
                "\"sat_with_exp\",list(I,J),"
                "Moddiq::primeTest_sat,"
                "Moddiq::deleteUnluckyPrimes_sat,"
                "Moddiq::pTest_sat);",
                "  return(result);",
                "}",
            ]
        )
    lines.extend(
        [
            "int stage_timer=timer;",
            matrix_line,
            "ideal maximal_minors=minor(M,15);",
            'print("MAXIMAL_MINORS="+string(size(maximal_minors)));',
            'print("MINORS_SECONDS="+string(timer-stage_timer));',
        ]
    )
    if stage == "minors":
        lines.extend(["quit;", ""])
        return "\n".join(lines)

    lines.extend(
        [
            "proc strip_boundary(poly f)",
            "{",
            "  while(subst(f,a,0)==0)",
            "  {",
            "    f=f/a;",
            "  }",
            "  while(subst(f,a,-1)==0)",
            "  {",
            "    f=f/(a+1);",
            "  }",
            "  return(f);",
            "}",
            "ideal primitive_minors;",
            "int minor_index;",
            "for(minor_index=1;"
            "minor_index<=size(maximal_minors);"
            "minor_index++)",
            "{",
            "  primitive_minors[minor_index]="
            "strip_boundary(maximal_minors[minor_index]);",
            "}",
            f"poly H={KERNEL_CHART};",
            "ideal total_boundary=a*(a+1)*H;",
            "stage_timer=timer;",
        ]
    )
    if method == "modular":
        lines.append(
            "list saturation_result="
            "modSat(primitive_minors,total_boundary);"
        )
    elif method == "modular-fast":
        lines.append(
            "list saturation_result="
            "fastModSat(primitive_minors,total_boundary);"
        )
    elif method == "modular-rebuild":
        lines.append(
            "list saturation_result="
            "rebuildModSat(primitive_minors,total_boundary);"
        )
    else:
        lines.append(
            "list saturation_result="
            "sat_with_exp(primitive_minors,total_boundary);"
        )
    lines.extend(
        [
            "ideal G=saturation_result[1];",
        ]
    )
    if method == "modular-rebuild":
        lines.append('attrib(G,"isSB",1);')
    lines.extend(
        [
            'print("SATURATION_SECONDS="+string(timer-stage_timer));',
            'print("STANDARD_BASIS_SIZE="+string(size(G)));',
        ]
    )
    lines.insert(
        -1,
        (
            'print("SATURATION_EXPONENT="'
            "+string(saturation_result[2]));"
        ),
    )
    if stage == "saturation":
        lines.extend(["quit;", ""])
        return "\n".join(lines)

    lines.append("stage_timer=timer;")
    if method == "modular-rebuild":
        lines.extend(
            [
                'print("FITTING_DIMENSION="+string(dim(G)));',
                'print("FITTING_LENGTH="+string(vdim(G)));',
                'print("INPUT_REMAINDER_TERMS=SKIPPED");',
                'print("CERTIFICATE_SECONDS=SKIPPED");',
            ]
        )
    else:
        lines.extend(
            [
                "ideal input_remainders=reduce(primitive_minors,G);",
                'print("FITTING_DIMENSION="+string(dim(G)));',
                'print("FITTING_LENGTH="+string(vdim(G)));',
                (
                    'print("INPUT_REMAINDER_TERMS="'
                    "+string(size(input_remainders)));"
                ),
                'print("CERTIFICATE_SECONDS="+string(timer-stage_timer));',
            ]
        )
    lines.extend(
        [
            "int i;",
            "for(i=1;i<=size(G);i++)",
            "{",
            '  print("BASIS_ELEMENT="+string(i)'
            '    +", TERMS="+string(size(G[i]))'
            '    +", DEG="+string(deg(G[i])));',
            '  print("LEADING_MONOMIAL="+string(leadmonom(G[i])));',
        ]
    )
    if print_basis:
        lines.append(
            '  print("BASIS_POLY="+string(G[i]));'
        )
    lines.extend(
        [
            "}",
        ]
    )
    if method == "modular-rebuild":
        lines.append('print("GROEBNER_VERIFIED=SKIPPED");')
    else:
        lines.append(
            'print("GROEBNER_VERIFIED="+string(system("verifyGB",G)));'
        )
    lines.extend(["quit;", ""])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--prime",
        type=int,
        help=(
            "one coefficient characteristic; zero requests Q. If omitted, "
            "run the bounded three-prime stability audit"
        ),
    )
    parser.add_argument(
        "--method",
        choices=(
            "auto",
            "native",
            "modular",
            "modular-fast",
            "modular-rebuild",
        ),
        default="auto",
        help=(
            "saturation method; auto uses modular saturation over "
            "Q and native saturation over finite fields"
        ),
    )
    parser.add_argument(
        "--basis-output",
        type=Path,
        help=(
            "write the reconstructed characteristic-zero basis as a "
            "generated Singular artifact"
        ),
    )
    parser.add_argument(
        "--stage",
        choices=("minors", "saturation", "certificate"),
        default="certificate",
        help="stop after the named bounded stage",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=min(os.cpu_count() or 1, 8),
        help="maximum Singular CPUs/threads for modular reconstruction",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=600,
        help="hard wall-clock limit in seconds; zero disables it",
    )
    args = parser.parse_args()
    if args.prime is not None and args.prime < 0:
        parser.error("--prime must be zero or positive")
    if args.cores < 1:
        parser.error("--cores must be positive")
    if args.timeout < 0:
        parser.error("--timeout must be nonnegative")
    if args.prime is None and args.method == "modular":
        parser.error("--method modular requires the explicit option --prime 0")
    if args.prime is None and args.method == "modular-fast":
        parser.error(
            "--method modular-fast requires the explicit option --prime 0"
        )
    if args.prime is None and args.method == "modular-rebuild":
        parser.error(
            "--method modular-rebuild requires the explicit option --prime 0"
        )
    if args.basis_output is not None:
        if args.prime != 0:
            parser.error("--basis-output requires --prime 0")
        if args.stage != "certificate":
            parser.error("--basis-output requires --stage certificate")

    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required")

    characteristics = (
        (31991, 32003, 65521)
        if args.prime is None
        else (args.prime,)
    )
    profiles = []
    for characteristic in characteristics:
        method = args.method
        if method == "auto":
            method = "modular" if characteristic == 0 else "native"
        if method == "modular" and characteristic != 0:
            parser.error("--method modular is only available with --prime 0")
        if method == "modular-fast" and characteristic != 0:
            parser.error(
                "--method modular-fast is only available with --prime 0"
            )
        if method == "modular-rebuild" and characteristic != 0:
            parser.error(
                "--method modular-rebuild is only available with --prime 0"
            )

        with tempfile.TemporaryDirectory(
            prefix="degree-five-qper-fitting-",
        ) as directory:
            path = Path(directory) / "fitting.sing"
            path.write_text(
                fitting_program(
                    characteristic=characteristic,
                    method=method,
                    stage=args.stage,
                    print_basis=args.basis_output is not None,
                )
            )
            command = [
                singular,
                f"--cpus={args.cores}",
                f"--threads={args.cores}",
                f"--flint-threads={args.cores}",
                "-q",
                str(path),
            ]
            try:
                result = subprocess.run(
                    command,
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=args.timeout or None,
                )
            except subprocess.TimeoutExpired as error:
                output = error.stdout or ""
                if isinstance(output, bytes):
                    output = output.decode(errors="replace")
                if output:
                    print(output, end="")
                raise SystemExit(
                    f"TIMEOUT: characteristic {characteristic}, stage "
                    f"{args.stage!r} exceeded {args.timeout} seconds"
                ) from error
        combined_output = result.stdout + result.stderr
        if any(
            line.lstrip().startswith("?")
            for line in combined_output.splitlines()
        ):
            raise SystemExit(combined_output)
        print(f"CHARACTERISTIC={characteristic}")
        print(f"METHOD={method}")
        print(f"CORES={args.cores}")
        displayed_stdout = result.stdout
        if args.basis_output is not None:
            displayed_stdout = "\n".join(
                line
                for line in result.stdout.splitlines()
                if not line.startswith("BASIS_POLY=")
            )
            if displayed_stdout:
                displayed_stdout += "\n"
        print(displayed_stdout, end="")
        if result.stderr:
            print(result.stderr, end="")

        values = {}
        leading_monomials = []
        for line in result.stdout.splitlines():
            if line.startswith("LEADING_MONOMIAL="):
                leading_monomials.append(line.split("=", 1)[1])
            elif "=" in line:
                key, value = line.split("=", 1)
                values[key] = value
        profiles.append((values, tuple(leading_monomials)))
        if args.basis_output is not None:
            basis_polynomials = [
                line.split("=", 1)[1]
                for line in result.stdout.splitlines()
                if line.startswith("BASIS_POLY=")
            ]
            if len(basis_polynomials) != 21:
                raise SystemExit(
                    "REBUILD=FAIL\n"
                    f"Expected 21 basis polynomials, got "
                    f"{len(basis_polynomials)}"
                )
            output_path = args.basis_output
            if not output_path.is_absolute():
                output_path = ROOT / output_path
            output_text = "\n".join(
                [
                    "// Reconstructed reduced q-period Fitting basis over Q.",
                    "// Reproduce with:",
                    "// .venv/bin/python "
                    "scripts/compute_degree_five_qper_fitting.py "
                    f"--prime 0 --method {args.method} "
                    f"--basis-output {args.basis_output}",
                    "ring r=0,(a,tau),dp;",
                    "ideal G;",
                    *(
                        f"G[{index}]={polynomial};"
                        for index, polynomial in enumerate(
                            basis_polynomials,
                            start=1,
                        )
                    ),
                    "",
                ]
            )
            output_path.parent.mkdir(parents=True, exist_ok=True)
            output_path.write_text(output_text)
            digest = hashlib.sha256(output_text.encode()).hexdigest()
            print(f"BASIS_OUTPUT={output_path}")
            print(f"BASIS_SHA256={digest}")

    if args.prime is None and args.stage == "certificate":
        expected = {
            "MAXIMAL_MINORS": "16",
            "SATURATION_EXPONENT": "12",
            "STANDARD_BASIS_SIZE": "21",
            "FITTING_DIMENSION": "0",
            "FITTING_LENGTH": "218",
            "INPUT_REMAINDER_TERMS": "0",
        }
        for values, _ in profiles:
            observed = {key: values.get(key) for key in expected}
            if observed != expected:
                raise SystemExit(
                    f"STABLE_AUDIT=FAIL\nEXPECTED={expected}\n"
                    f"OBSERVED={observed}"
                )
        fingerprints = {profile[1] for profile in profiles}
        if len(fingerprints) != 1:
            raise SystemExit(
                "STABLE_AUDIT=FAIL\n"
                "The leading-monomial fingerprints differ by prime"
            )
        print("STABLE_AUDIT=PASS")


if __name__ == "__main__":
    main()
