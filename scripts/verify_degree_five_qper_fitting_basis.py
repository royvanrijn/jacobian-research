#!/usr/bin/env python3
"""Certify a rebuilt characteristic-zero q-period Fitting basis.

For the primitive maximal-minor ideal I and boundary product B=a(a+1)H,
the four checks certify

    I : B^infinity = G.

Indeed, input containment gives I subset G, the boundary-unit check makes G
B-saturated, and reverse membership gives B^12 G subset I.  The Groebner
check makes the reported 218-monomial staircase independently trustworthy.
Each expensive check runs in a separate bounded Singular process.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
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
DEFAULT_BASIS = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_five_qper_fitting_basis_Q.sing"
)
KERNEL_CHART = (
    "4*a^3*tau^2-24*a^3*tau-72*a^3"
    "+8*a^2*tau^2-54*a^2*tau-216*a^2"
    "+4*a*tau^2-30*a*tau-246*a-105"
)
CHECKS = ("shape", "groebner", "input", "boundary-unit", "reverse")
STAIRCASE_ORDER = (*range(14, 22), *range(1, 14))
CRITICAL_PAIRS = tuple(zip(STAIRCASE_ORDER, STAIRCASE_ORDER[1:]))


def frozen_basis_lines(path: Path) -> list[str]:
    lines = [
        line
        for line in path.read_text().splitlines()
        if line.startswith("G[")
    ]
    if len(lines) != 21:
        raise ValueError(f"expected 21 basis assignments in {path}, got {len(lines)}")
    return lines


def primitive_minor_lines() -> list[str]:
    matrix_line = next(
        line
        for line in PRESENTATION.read_text().splitlines()
        if line.startswith("matrix M[15][16]=")
    )
    return [
        matrix_line,
        "ideal maximal_minors=minor(M,15);",
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
    ]


def certificate_program(
    check: str,
    basis_path: Path,
    index: int | None = None,
    pair: tuple[int, int] | None = None,
    jobs: int = 1,
) -> str:
    lines = [
        f'execute(read("{basis_path.as_posix()}"));',
        "short=0;",
        f"poly H={KERNEL_CHART};",
        "poly B=a*(a+1)*H;",
        "int certificate_timer=timer;",
    ]
    if check == "groebner" and pair is not None:
        lines.insert(1, 'LIB "teachstd.lib";')
    if check in {"input", "reverse"}:
        lines.extend(primitive_minor_lines())
    if check == "shape":
        lines.extend(
            [
                'attrib(G,"isSB",1);',
                'print("BASIS_SIZE="+string(size(G)));',
                'print("DIMENSION="+string(dim(G)));',
                'print("LENGTH="+string(vdim(G)));',
                "int i;",
                "for(i=1;i<=size(G);i++)",
                "{",
                '  print("LEADING_MONOMIAL="+string(leadmonom(G[i])));',
                "}",
            ]
        )
    elif check == "groebner":
        if pair is None:
            lines.append(
                'print("GROEBNER_VERIFIED="+string(system("verifyGB",G)));'
            )
        else:
            first, second = pair
            lines.extend(
                [
                    'attrib(G,"isSB",1);',
                    (
                        f"poly spair=spoly(G[{first}],G[{second}]);"
                    ),
                    "poly spair_remainder=reduce(spair,G);",
                    (
                        f'print("SPAIR_{first}_{second}_TERMS="'
                        "+string(size(spair_remainder)));"
                    ),
                ]
            )
    elif check == "input":
        lines.extend(
            [
                'attrib(G,"isSB",1);',
            ]
        )
        if index is None:
            lines.extend(
                [
                    "ideal input_remainders=reduce(primitive_minors,G);",
                    (
                        'print("INPUT_REMAINDER_TERMS="'
                        "+string(size(input_remainders)));"
                    ),
                ]
            )
        else:
            lines.extend(
                [
                    (
                        "poly input_remainder="
                        f"reduce(primitive_minors[{index}],G);"
                    ),
                    (
                        f'print("INPUT_{index}_REMAINDER_TERMS="'
                        "+string(size(input_remainder)));"
                    ),
                ]
            )
    elif check == "boundary-unit":
        lines.extend(
            [
                "ideal boundary_test=std(G+B);",
                "poly unit_remainder=reduce(1,boundary_test);",
                (
                    'print("BOUNDARY_UNIT_REMAINDER_TERMS="'
                    "+string(size(unit_remainder)));"
                ),
            ]
        )
    elif check == "reverse":
        lines.extend(
            [
                "ideal reverse_targets=B^12*G;",
                "matrix reverse_lift=lift(primitive_minors,reverse_targets);",
                (
                    "matrix reverse_residual="
                    "matrix(reverse_targets)"
                    "-matrix(primitive_minors)*reverse_lift;"
                ),
                "ideal reverse_residual_ideal=ideal(reverse_residual);",
                (
                    'print("REVERSE_RESIDUAL_TERMS="'
                    "+string(size(reverse_residual_ideal)));"
                ),
                (
                    'print("REVERSE_LIFT_ROWS="'
                    "+string(nrows(reverse_lift)));"
                ),
                (
                    'print("REVERSE_LIFT_COLS="'
                    "+string(ncols(reverse_lift)));"
                ),
            ]
        )
    else:
        raise ValueError(check)
    lines.extend(
        [
            (
                'print("CERTIFICATE_SECONDS="'
                "+string(timer-certificate_timer));"
            ),
            "quit;",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--basis",
        type=Path,
        default=DEFAULT_BASIS,
        help="rebuilt rational basis artifact",
    )
    parser.add_argument(
        "--check",
        choices=(*CHECKS, "all"),
        default="all",
    )
    parser.add_argument(
        "--index",
        type=int,
        help="check only this input generator; requires --check input",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(os.cpu_count() or 1, 8),
        help="parallel Singular workers for decomposable exact checks",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=900,
        help="hard wall-clock limit in seconds for each check",
    )
    args = parser.parse_args()
    if args.timeout < 1:
        parser.error("--timeout must be positive")
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.index is not None:
        if args.check != "input":
            parser.error("--index requires --check input")
        if not 1 <= args.index <= 16:
            parser.error("--index must lie between 1 and 16")
    basis_path = args.basis
    if not basis_path.is_absolute():
        basis_path = ROOT / basis_path
    if not basis_path.exists():
        raise SystemExit(f"missing rebuilt basis: {basis_path}")
    frozen_basis_lines(basis_path)
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required")

    def run_one(
        check: str,
        label: str,
        index: int | None = None,
        pair: tuple[int, int] | None = None,
    ) -> str:
        with tempfile.TemporaryDirectory(
            prefix=f"degree-five-qper-{label}-",
        ) as directory:
            program = Path(directory) / "certificate.sing"
            program.write_text(
                certificate_program(
                    check,
                    basis_path,
                    index=index,
                    pair=pair,
                    jobs=args.jobs,
                )
            )
            try:
                result = subprocess.run(
                    [
                        singular,
                        f"--cpus={args.jobs}",
                        f"--threads={args.jobs}",
                        "-q",
                        str(program),
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                    timeout=args.timeout,
                )
            except subprocess.TimeoutExpired as error:
                raise RuntimeError(
                    f"TIMEOUT: check {label!r} exceeded "
                    f"{args.timeout} seconds"
                ) from error
        output = result.stdout + result.stderr
        if any(
            line.lstrip().startswith("?")
            for line in output.splitlines()
        ):
            raise RuntimeError(f"CHECK={label}\n{output}")
        return output

    checks = CHECKS if args.check == "all" else (args.check,)
    for check in checks:
        if check == "input" and args.index is not None:
            output = run_one(
                check,
                f"input-{args.index}",
                index=args.index,
            )
            print(f"CHECK={check}")
            print(output, end="")
            continue
        if check == "input" and args.jobs > 1:
            tasks = [
                (f"input-{index}", index, None)
                for index in range(1, 17)
            ]
        elif check == "groebner" and args.jobs > 1:
            tasks = [
                (
                    f"groebner-{first}-{second}",
                    None,
                    (first, second),
                )
                for first, second in CRITICAL_PAIRS
            ]
        else:
            output = run_one(check, check)
            print(f"CHECK={check}")
            print(output, end="")
            continue

        outputs = {}
        failures = {}
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {
                executor.submit(
                    run_one,
                    check,
                    label,
                    index,
                    pair,
                ): label
                for label, index, pair in tasks
            }
            for future in as_completed(futures):
                label = futures[future]
                try:
                    outputs[label] = future.result()
                    print(f"COMPLETED={label}", flush=True)
                except Exception as error:
                    failures[label] = str(error)
                    print(f"FAILED={label}: {error}", flush=True)
        print(f"CHECK={check}")
        for label, _, _ in tasks:
            if label in outputs:
                print(outputs[label], end="")
        print(f"PARALLEL_{check.upper()}_WORKERS={args.jobs}")
        if check == "groebner":
            print(f"CRITICAL_PAIRS={len(CRITICAL_PAIRS)}")
            print("CHAIN_CRITERION=two-variable-adjacent-staircase")
        if failures:
            raise SystemExit(
                "FAILED_CHECKS=" + ",".join(sorted(failures))
            )


if __name__ == "__main__":
    main()
