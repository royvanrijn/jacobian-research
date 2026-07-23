#!/usr/bin/env python3
"""Analyze exported generic bounded-period factors.

The inputs are polynomial expressions in Singular syntax, produced by
``compute_degree_five_fifth_order_function_field.py --factor-output``.
The calculation is performed over a chosen finite field.  It reports the
resultant in ``tau``, its factor degrees, and every intersection fiber whose
``a``-coordinate is rational over that finite field.

The ``--verify-cubic-component`` mode instead works over Q and verifies the
reconstructed common component

    94*a^3+335*a^2+400*a+160 = 0,
    8*tau+658*a^2+1593*a+976 = 0

by exact reduction of every supplied factor.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from fractions import Fraction
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("factors", nargs="+")
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument(
        "--ideal-only",
        action="store_true",
        help="skip the verbose first-pair resultant analysis",
    )
    parser.add_argument(
        "--verify-cubic-component",
        action="store_true",
        help="verify the reconstructed common cubic component over Q",
    )
    parser.add_argument(
        "--local-cubic-multiplicity",
        action="store_true",
        help="compute the local intersection algebra at the cubic component",
    )
    parser.add_argument(
        "--local-rational-point",
        nargs=2,
        metavar=("A", "TAU"),
        help="compute the local algebra at a rational (a,tau) point",
    )
    parser.add_argument(
        "--local-tauzero-quadratic",
        action="store_true",
        help="compute locally at tau=0, 11*a^2+28*a+18=0",
    )
    parser.add_argument(
        "--local-a0-quadratic",
        action="store_true",
        help="compute locally at a=0, 6*tau^2+4*tau+9=0",
    )
    parser.add_argument(
        "--print-local-basis",
        action="store_true",
        help="print the local standard basis as well as its multiplicity",
    )
    parser.add_argument(
        "--exact-decomposition",
        action="store_true",
        help="certify exhaustive characteristic-zero support and total length",
    )
    args = parser.parse_args()

    if len(args.factors) < 2:
        parser.error("provide at least two factor files")
    if args.exact_decomposition and (
        args.verify_cubic_component
        or args.local_cubic_multiplicity
        or args.local_rational_point is not None
        or args.local_tauzero_quadratic
        or args.local_a0_quadratic
    ):
        parser.error("--exact-decomposition is a standalone mode")
    factor_texts = [
        Path(path).read_text().strip()
        for path in args.factors
    ]
    if args.exact_decomposition:
        script_path = str(Path(__file__).resolve())
        factor_paths = [str(Path(path).resolve()) for path in args.factors]
        components = (
            ("a0_tau0", 1, ("--local-rational-point", "0", "0")),
            (
                "a0_tau_quadratic",
                2,
                ("--local-a0-quadratic",),
            ),
            (
                "a_minus_three_halves_tau0",
                1,
                ("--local-rational-point", "-1.5", "0"),
            ),
            (
                "tau0_a_quadratic",
                2,
                ("--local-tauzero-quadratic",),
            ),
            (
                "a_minus_one_half_tau_minus_three",
                1,
                ("--local-rational-point", "-0.5", "-3"),
            ),
            ("cubic", 3, ("--local-cubic-multiplicity",)),
        )
        exact_length = 0
        for label, residue_degree, mode_arguments in components:
            result = subprocess.run(
                [
                    sys.executable,
                    script_path,
                    *mode_arguments,
                    *factor_paths,
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            match = re.search(r"_MULTIPLICITY=(\d+)", result.stdout)
            if match is None:
                raise AssertionError(
                    f"no local multiplicity for {label}:\n{result.stdout}"
                )
            multiplicity = int(match.group(1))
            contribution = residue_degree * multiplicity
            exact_length += contribution
            print(
                f"COMPONENT={label},"
                f"RESIDUE_DEGREE={residue_degree},"
                f"LOCAL_MULTIPLICITY={multiplicity},"
                f"CONTRIBUTION={contribution}"
            )
        modular = subprocess.run(
            [
                sys.executable,
                script_path,
                "--prime",
                "32003",
                "--ideal-only",
                *factor_paths,
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        modular_match = re.search(
            r"ALL_VECTOR_DIMENSION=(\d+)",
            modular.stdout,
        )
        if modular_match is None:
            raise AssertionError(
                "no modular vector dimension:\n" + modular.stdout
            )
        modular_upper_bound = int(modular_match.group(1))
        print(f"EXACT_COMPONENT_LENGTH_SUM={exact_length}")
        print(f"MODULAR_UPPER_BOUND={modular_upper_bound}")
        if exact_length != modular_upper_bound:
            raise AssertionError(
                (exact_length, modular_upper_bound)
            )
        print("EXACT_SUPPORT_EXHAUSTED=1")
        return

    local_mode_count = sum(
        (
            args.local_cubic_multiplicity,
            args.local_rational_point is not None,
            args.local_tauzero_quadratic,
            args.local_a0_quadratic,
        )
    )
    if args.verify_cubic_component and local_mode_count:
        parser.error("choose one characteristic-zero analysis mode")
    if local_mode_count > 1:
        parser.error("choose only one local multiplicity mode")
    if local_mode_count:
        if args.local_cubic_multiplicity:
            ring_declaration = (
                "ring r=(0,a),(x,y),ds;\n"
                "minpoly=94*a^3+335*a^2+400*a+160;"
            )
            center_a = "a"
            center_tau = "((-329/4)*a^2+(-1593/8)*a-122)"
            label = "CUBIC"
        elif args.local_tauzero_quadratic:
            ring_declaration = (
                "ring r=(0,a),(x,y),ds;\n"
                "minpoly=11*a^2+28*a+18;"
            )
            center_a = "a"
            center_tau = "0"
            label = "TAUZERO_QUADRATIC"
        elif args.local_a0_quadratic:
            ring_declaration = (
                "ring r=(0,t),(x,y),ds;\n"
                "minpoly=6*t^2+4*t+9;"
            )
            center_a = "0"
            center_tau = "t"
            label = "A0_QUADRATIC"
        else:
            a_center = Fraction(args.local_rational_point[0])
            tau_center = Fraction(args.local_rational_point[1])

            def rational_text(value: Fraction) -> str:
                if value.denominator == 1:
                    return str(value.numerator)
                return f"({value.numerator}/{value.denominator})"

            ring_declaration = "ring r=0,(x,y),ds;"
            center_a = rational_text(a_center)
            center_tau = rational_text(tau_center)
            label = "RATIONAL"

        local_factors = []
        for factor in factor_texts:
            translated = re.sub(
                r"\ba\b",
                f"({center_a}+x)",
                factor,
            )
            translated = re.sub(
                r"\btau\b",
                f"({center_tau}+y)",
                translated,
            )
            local_factors.append(translated)
        local_definitions = "\n".join(
            f"poly F{index}={factor};"
            for index, factor in enumerate(local_factors, start=1)
        )
        local_names = ",".join(
            f"F{index}" for index in range(1, len(local_factors) + 1)
        )
        basis_lines = (
            f'print("LOCAL_{label}_STANDARD_BASIS=");\nprint(G);'
            if args.print_local_basis
            else ""
        )
        program = f"""
{ring_declaration}
short=0;
{local_definitions}
ideal I={local_names};
ideal G=std(I);
print("LOCAL_{label}_GB_SIZE="+string(size(G)));
print("LOCAL_{label}_DIMENSION="+string(dim(G)));
if(dim(G)==0)
{{
  print("LOCAL_{label}_MULTIPLICITY="+string(vdim(G)));
}}
{basis_lines}
quit;
"""
        with tempfile.TemporaryDirectory(
            prefix="degree-five-period-cubic-local-",
        ) as directory:
            script = Path(directory) / "cubic-local.sing"
            script.write_text(program)
            result = subprocess.run(
                ["Singular", "-q", str(script)],
                check=True,
                capture_output=True,
                text=True,
            )
        print(result.stdout, end="")
        if result.stderr:
            print(result.stderr, end="")
        return

    first, second = factor_texts[:2]
    extra_definitions = "\n".join(
        f"poly F{index}={factor};"
        for index, factor in enumerate(factor_texts[2:], start=3)
    )
    all_names = ",".join(
        ["P", "Q"]
        + [
            f"F{index}"
            for index in range(3, len(factor_texts) + 1)
        ]
    )
    all_fiber_gcd_lines = "\n".join(
        f"      all_fiber=gcd(all_fiber,subst(F{index},a,all_root));"
        for index in range(3, len(factor_texts) + 1)
    )
    if args.verify_cubic_component:
        exact_definitions = "\n".join(
            f"poly F{index}={factor};"
            for index, factor in enumerate(factor_texts, start=1)
        )
        reduction_lines = "\n".join(
            f'print("F{index}_CUBIC_REMAINDER="+string(reduce(F{index},JG)));'
            f'\nprint("F{index}_KAPPA0_TAUM3="'
            f'+string(subst(subst(F{index},a,-1/2),tau,-3)));'
            for index in range(1, len(factor_texts) + 1)
        )
        fiber_gcd_updates = "\n".join(
            f"tau_zero_gcd=gcd(tau_zero_gcd,subst(F{index},tau,0));"
            f"\na_zero_gcd=gcd(a_zero_gcd,subst(F{index},a,0));"
            f"\na_m32_gcd=gcd(a_m32_gcd,subst(F{index},a,-3/2));"
            f"\na_m12_gcd=gcd(a_m12_gcd,subst(F{index},a,-1/2));"
            for index in range(2, len(factor_texts) + 1)
        )
        program = f"""
ring r=0,(a,tau),dp;
short=0;
{exact_definitions}
poly C=94*a^3+335*a^2+400*a+160;
poly L=8*tau+658*a^2+1593*a+976;
ideal J=C,L;
ideal JG=std(J);
print("CUBIC_IDEAL_DIMENSION="+string(dim(JG)));
{reduction_lines}
poly tau_zero_gcd=subst(F1,tau,0);
poly a_zero_gcd=subst(F1,a,0);
poly a_m32_gcd=subst(F1,a,-3/2);
poly a_m12_gcd=subst(F1,a,-1/2);
{fiber_gcd_updates}
print("TAU_ZERO_GCD="+string(tau_zero_gcd));
print("A_ZERO_GCD="+string(a_zero_gcd));
print("A_MINUS_THREE_HALVES_GCD="+string(a_m32_gcd));
print("A_MINUS_ONE_HALF_GCD="+string(a_m12_gcd));
quit;
"""
        with tempfile.TemporaryDirectory(
            prefix="degree-five-period-cubic-check-",
        ) as directory:
            script = Path(directory) / "cubic-check.sing"
            script.write_text(program)
            result = subprocess.run(
                ["Singular", "-q", str(script)],
                check=True,
                capture_output=True,
                text=True,
            )
        print(result.stdout, end="")
        if result.stderr:
            print(result.stderr, end="")
        return

    resultant_program = "" if args.ideal_only else """
poly R=resultant(P,Q,tau);
list RF=factorize(R);
print("RESULTANT_TERMS="+string(size(R)));
print("RESULTANT_DEGREE="+string(deg(R)));
print("RESULTANT_FACTOR_COUNT="+string(size(RF[1])));
int i;
for(i=1;i<=size(RF[1]);i++)
{
  if(deg(RF[1][i])<=12)
  {
    print("RFACTOR="+string(RF[1][i])
      +", EXP="+string(RF[2][i])
      +", TERMS="+string(size(RF[1][i]))
      +", DEG="+string(deg(RF[1][i])));
  }
  else
  {
    print("RFACTOR=<omitted>"
      +", EXP="+string(RF[2][i])
      +", TERMS="+string(size(RF[1][i]))
      +", DEG="+string(deg(RF[1][i])));
  }
  if(deg(RF[1][i])==1)
  {
    number root=-leadcoef(subst(RF[1][i],a,0))/leadcoef(RF[1][i]);
    poly Ps=subst(P,a,root);
    poly Qs=subst(Q,a,root);
    poly fiber_gcd=gcd(Ps,Qs);
    print("LINEAR_A_ROOT="+string(root));
    print("FIBER_GCD="+string(fiber_gcd));
  }
}
"""
    program = f"""
ring r={args.prime},(a,tau),dp;
short=0;
poly P={first};
poly Q={second};
{extra_definitions}
poly common=gcd(P,Q);
print("BIVARIATE_GCD="+string(common));
ideal ALL={all_names};
ideal ALL_GB=std(ALL);
poly unit_test=reduce(1,ALL_GB);
print("ALL_GENERATORS="+string(size(ALL)));
print("ALL_GB_SIZE="+string(size(ALL_GB)));
print("ALL_UNIT_REMAINDER="+string(unit_test));
print("ALL_DIMENSION="+string(dim(ALL_GB)));
if(dim(ALL_GB)==0)
{{
  print("ALL_VECTOR_DIMENSION="+string(vdim(ALL_GB)));
  ideal ALL_ELIM=eliminate(ALL,tau);
  print("ALL_ELIM_GENERATORS="+string(size(ALL_ELIM)));
  if(size(ALL_ELIM)>0)
  {{
    poly ALL_A=ALL_ELIM[1];
    list ALL_AF=factorize(ALL_A);
    print("ALL_ELIM_DEGREE="+string(deg(ALL_A)));
    print("ALL_ELIM_FACTOR_COUNT="+string(size(ALL_AF[1])));
    int all_factor_index;
    number all_root;
    poly all_fiber;
    for(all_factor_index=1;
        all_factor_index<=size(ALL_AF[1]);
        all_factor_index++)
    {{
      print("ALL_EFACTOR="+string(ALL_AF[1][all_factor_index])
        +", DEG="
        +string(deg(ALL_AF[1][all_factor_index]))
        +", EXP="+string(ALL_AF[2][all_factor_index]));
      if(deg(ALL_AF[1][all_factor_index])==1)
      {{
        all_root=-leadcoef(
          subst(ALL_AF[1][all_factor_index],a,0)
        )/leadcoef(ALL_AF[1][all_factor_index]);
        all_fiber=subst(P,a,all_root);
        all_fiber=gcd(all_fiber,subst(Q,a,all_root));
{all_fiber_gcd_lines}
        print("ALL_LINEAR_A_ROOT="+string(all_root));
        print("ALL_FIBER_GCD="+string(all_fiber));
      }}
    }}
  }}
}}
{resultant_program}
quit;
"""
    with tempfile.TemporaryDirectory(
        prefix="degree-five-period-intersection-",
    ) as directory:
        script = Path(directory) / "intersection.sing"
        script.write_text(program)
        result = subprocess.run(
            ["Singular", "-q", str(script)],
            check=True,
            capture_output=True,
            text=True,
        )
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")


if __name__ == "__main__":
    main()
