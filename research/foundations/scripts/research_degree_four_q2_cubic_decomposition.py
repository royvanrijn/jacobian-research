#!/usr/bin/env python3
"""Decompose the cubic normal jet in the quartic q2 synchronization test.

This is an exact finite-field calculation on the same formal normal slice
used by ``research_degree_four_q2_augmented_nullcone.py``.  It does not
prove a characteristic-zero or global nullcone statement.

At the retained prime 32003 the script:

* proves that x7 is nilpotent in the cubic moment-jet algebra;
* decomposes the reduced cubic support into its x6=0 part and the
  saturation away from x6;
* proves the radical of the x6=0 part and factors its binary cubic;
* computes the degree-nine generic fiber of the off-x6 saturation;
* restricts the quartic moment jet to the linear and quadratic factors
  of the dominant cubic component and proves that both radicals collapse
  to the same three-plane.

The off-x6 quartic restriction is deliberately not claimed: the direct
parameter-field standard-basis calculation currently times out.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

from research_degree_four_q2_augmented_nullcone import (
    Polynomial,
    build_moment_jets,
    formal_pivot_elimination,
    polynomial_string,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree_four_q2_cubic_decomposition.json"
)


def run_singular(code: str, timeout: int) -> str:
    executable = shutil.which("Singular")
    if executable is None:
        raise SystemExit("Singular is required on PATH")
    completed = subprocess.run(
        [executable, "-q"],
        input=code,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=True,
    )
    if completed.stderr.strip():
        raise AssertionError(completed.stderr)
    return completed.stdout


def marker(output: str, name: str) -> list[str]:
    prefix = f"{name} "
    for line in output.splitlines():
        stripped = line.strip()
        if stripped.startswith(prefix):
            return stripped.split()[1:]
    raise AssertionError(f"missing Singular marker {name!r}\n{output}")


def specialize_linear(
    polynomial: Polynomial,
    variable_count: int,
    assignments: list[tuple[int, int] | None],
    prime: int,
) -> Polynomial:
    """Apply a homogeneous monomial linear specialization."""
    result: defaultdict[tuple[int, ...], int] = defaultdict(int)
    for monomial, coefficient in polynomial.items():
        target = [0] * variable_count
        value = coefficient
        for source, exponent in enumerate(monomial):
            if not exponent:
                continue
            assignment = assignments[source]
            if assignment is None:
                value = 0
                break
            variable, scalar = assignment
            target[variable] += exponent
            value = value * pow(scalar, exponent, prime) % prime
        if value:
            exponent = tuple(target)
            result[exponent] = (
                result[exponent] + value
            ) % prime
    return {
        exponent: coefficient
        for exponent, coefficient in result.items()
        if coefficient
    }


def cubic_analysis(
    moment_jets: dict[int, Polynomial],
    substitutions: list[Polynomial],
    prime: int,
    timeout: int,
) -> dict[str, object]:
    source_variables = [f"y{index}" for index in range(12)]
    target_variables = [f"x{index}" for index in range(8)]
    source_ideal = ",".join(
        polynomial_string(
            {
                monomial: coefficient
                for monomial, coefficient in moment_jets[order].items()
                if sum(monomial) <= 3
            },
            source_variables,
        )
        for order in range(6, 22)
    )
    map_entries = ",".join(
        polynomial_string(
            {
                monomial: coefficient
                for monomial, coefficient in polynomial.items()
                if sum(monomial) <= 3
            },
            target_variables,
        )
        for polynomial in substitutions
    )
    binary_cubic = (
        "x2^3+6306*x2^2*x5-11062*x2*x5^2+1735*x5^3"
    )
    linear_factor = "x2-7378*x5"
    quadratic_factor = (
        "x2^2+13684*x2*x5+12028*x5^2"
    )
    code = f"""
LIB "elim.lib";
LIB "solve.lib";
proc contained(ideal A, ideal B)
{{
  for (int k=1; k<=size(A); k++)
  {{
    if (reduce(A[k],B)!=0) {{ return(0); }}
  }}
  return(1);
}}
ring ra={prime},({','.join(source_variables)}),dp;
ideal IA={source_ideal};
ring rb={prime},({','.join(target_variables)}),dp;
map phi=ra,{map_entries};
ideal I=slimgb(jet(phi(IA),3));
ideal J=slimgb(I+ideal(x7));
ideal Jx=slimgb(J+ideal(x6));
poly P={binary_cubic};
poly LF={linear_factor};
poly QF={quadratic_factor};
ideal C=std(ideal(x7,x6,x3,P));
list PF=factorize(P,1);
list FS=facstd(J);
ideal Embedded=slimgb(FS[1]);
ideal EmbeddedSupport=std(ideal(x7,x6,x5,x3,x2));
ideal K=sat(J,ideal(x6));
K=slimgb(K);
list IS=indepSet(K,0);
"CUBIC",size(I),dim(I),mult(I);
"X7_CUBE",reduce(x7^3,I)==0;
"X7_REDUCED",size(J),dim(J),mult(J);
"FACSTD_COUNT",size(FS);
for (int k=1; k<=size(FS); k++)
{{
  ideal B=slimgb(FS[k]);
  "FACSTD_BRANCH",k,dim(B),mult(B),size(B);
}}
"X6_ZERO_RADICAL",contained(Jx,C),
  reduce(x3^3,Jx)==0,reduce(P^2,Jx)==0,
  P-LF*QF==0,size(PF[1]),dim(C),mult(C),size(C);
"EMBEDDED_BRANCH",contained(Embedded,EmbeddedSupport),
  reduce(x2^6,Embedded)==0,reduce(x3^3,Embedded)==0;
"OFF_AXIS",dim(K),mult(K),size(K);
"OFF_AXIS_INDEP",IS[1];
ideal K0=subst(K,x7,0);
ring rc=({prime},x1,x4,x6),(x0,x2,x3,x5),dp;
ideal T=imap(rb,K0);
option(redSB);
ideal G=std(T);
int genericBasisSize=size(G);
int genericDegree=vdim(G);
ring rd=({prime},x1,x4,x6),(x0,x2,x3,x5),lp;
ideal L=fglm(rc,G);
list FF=factorize(L[1],1);
"OFF_AXIS_GENERIC",genericBasisSize,genericDegree,
  size(L),deg(L[1]),
  size(FF[1]);
"""
    output = run_singular(code, timeout)

    cubic = [int(value) for value in marker(output, "CUBIC")]
    assert cubic == [85, 4, 9], cubic
    assert marker(output, "X7_CUBE") == ["1"]
    reduced = [
        int(value) for value in marker(output, "X7_REDUCED")
    ]
    assert reduced == [63, 4, 6], reduced
    assert marker(output, "FACSTD_COUNT") == ["3"]
    branches = [
        [int(value) for value in line.split()[1:]]
        for line in output.splitlines()
        if line.strip().startswith("FACSTD_BRANCH ")
    ]
    assert branches == [
        [1, 3, 6, 7],
        [2, 4, 3, 4],
        [3, 3, 9, 38],
    ], branches
    radical = [
        int(value)
        for value in marker(output, "X6_ZERO_RADICAL")
    ]
    assert radical == [1, 1, 1, 1, 2, 4, 3, 4], radical
    embedded = [
        int(value)
        for value in marker(output, "EMBEDDED_BRANCH")
    ]
    assert embedded == [1, 1, 1], embedded
    off_axis = [
        int(value) for value in marker(output, "OFF_AXIS")
    ]
    assert off_axis == [3, 9, 38], off_axis
    independent = [
        int(value)
        for value in ",".join(
            marker(output, "OFF_AXIS_INDEP")
        ).split(",")
    ]
    assert independent == [0, 1, 0, 0, 1, 0, 1, 0]
    generic = [
        int(value)
        for value in marker(output, "OFF_AXIS_GENERIC")
    ]
    assert generic == [10, 9, 4, 9, 1], generic

    return {
        "cubic_ideal": {
            "standard_basis_size": cubic[0],
            "dimension": cubic[1],
            "multiplicity": cubic[2],
            "x7_cubed_is_zero": True,
        },
        "after_adjoining_x7": {
            "standard_basis_size": reduced[0],
            "dimension": reduced[1],
            "multiplicity": reduced[2],
        },
        "factorizing_standard_basis_branches": [
            {
                "branch": branch[0],
                "dimension": branch[1],
                "multiplicity": branch[2],
                "standard_basis_size": branch[3],
                **(
                    {
                        "support_ideal": [
                            "x7",
                            "x6",
                            "x5",
                            "x3",
                            "x2",
                        ],
                        "support_contained_in_dominant_component": (
                            True
                        ),
                        "x2_nilpotence_exponent_certified": 6,
                        "x3_nilpotence_exponent_certified": 3,
                    }
                    if branch[0] == 1
                    else {}
                ),
            }
            for branch in branches
        ],
        "x6_zero_radical": {
            "ideal": [
                "x7",
                "x6",
                "x3",
                binary_cubic,
            ],
            "containment_of_cubic_ideal": True,
            "x3_nilpotence_exponent_certified": 3,
            "binary_cubic_nilpotence_exponent_certified": 2,
            "binary_cubic_factorization": [
                linear_factor,
                quadratic_factor,
            ],
            "factor_count_over_finite_field": radical[4],
            "dimension": radical[5],
            "multiplicity": radical[6],
            "standard_basis_size": radical[7],
        },
        "off_x6_saturation": {
            "dimension": off_axis[0],
            "multiplicity": off_axis[1],
            "standard_basis_size": off_axis[2],
            "independent_variables": ["x1", "x4", "x6"],
            "generic_fiber": {
                "coefficient_field": (
                    "F_32003(x1,x4,x6)"
                ),
                "dependent_variables": [
                    "x0",
                    "x2",
                    "x3",
                    "x5",
                ],
                "degree": generic[1],
                "degree_order_basis_size": generic[0],
                "lex_basis_size": generic[2],
                "primitive_polynomial_degree": generic[3],
                "primitive_polynomial_factor_count": generic[4],
            },
        },
    }


def quartic_sheet_analysis(
    moment_jets: dict[int, Polynomial],
    substitutions: list[Polynomial],
    prime: int,
    timeout: int,
) -> dict[str, object]:
    source_variables = [f"y{index}" for index in range(12)]
    source_ideal = ",".join(
        polynomial_string(moment_jets[order], source_variables)
        for order in range(6, 22)
    )

    linear_variables = [f"z{index}" for index in range(4)]
    # x0,x1,x2,x3,x4,x5,x6,x7 ->
    # z0,z1,7378*z3,0,z2,z3,0,0.
    linear_assignments = [
        (0, 1),
        (1, 1),
        (3, 7378),
        None,
        (2, 1),
        (3, 1),
        None,
        None,
    ]
    linear_entries = ",".join(
        polynomial_string(
            specialize_linear(
                polynomial,
                4,
                linear_assignments,
                prime,
            ),
            linear_variables,
        )
        for polynomial in substitutions
    )
    linear_code = f"""
proc contained(ideal A, ideal B)
{{
  for (int k=1; k<=size(A); k++)
  {{
    if (reduce(A[k],B)!=0) {{ return(0); }}
  }}
  return(1);
}}
ring ra={prime},({','.join(source_variables)}),dp;
ideal IA={source_ideal};
ring rb={prime},({','.join(linear_variables)}),dp;
map phi=ra,{linear_entries};
ideal I=jet(phi(IA),4);
ideal G=slimgb(I);
ideal C=std(ideal(z3));
"LINEAR_SHEET",dim(G),mult(G),size(G),
  contained(I,C),reduce(z3^4,G)==0;
"""
    linear_output = run_singular(linear_code, timeout)
    linear = [
        int(value)
        for value in marker(linear_output, "LINEAR_SHEET")
    ]
    assert linear == [3, 2, 4, 1, 1], linear

    quadratic_variables = [f"z{index}" for index in range(5)]
    # x0,x1,x2,x3,x4,x5,x6,x7 ->
    # z0,z1,z2,0,z3,z4,0,0.
    quadratic_assignments = [
        (0, 1),
        (1, 1),
        (2, 1),
        None,
        (3, 1),
        (4, 1),
        None,
        None,
    ]
    quadratic_entries = ",".join(
        polynomial_string(
            specialize_linear(
                polynomial,
                5,
                quadratic_assignments,
                prime,
            ),
            quadratic_variables,
        )
        for polynomial in substitutions
    )
    quadratic_factor = (
        "z2^2+13684*z2*z4+12028*z4^2"
    )
    quadratic_code = f"""
proc contained(ideal A, ideal B)
{{
  for (int k=1; k<=size(A); k++)
  {{
    if (reduce(A[k],B)!=0) {{ return(0); }}
  }}
  return(1);
}}
ring ra={prime},({','.join(source_variables)}),dp;
ideal IA={source_ideal};
ring rb={prime},({','.join(quadratic_variables)}),dp;
map phi=ra,{quadratic_entries};
poly Q={quadratic_factor};
ideal I=jet(phi(IA),4)+ideal(Q);
ideal G=slimgb(I);
ideal C=std(ideal(z2,z4));
list QF=factorize(Q,1);
"QUADRATIC_SHEET",dim(G),mult(G),size(G),
  contained(I,C),reduce(z4^5,G)==0,size(QF[1]);
"""
    quadratic_output = run_singular(quadratic_code, timeout)
    quadratic = [
        int(value)
        for value in marker(
            quadratic_output, "QUADRATIC_SHEET"
        )
    ]
    assert quadratic == [3, 4, 15, 1, 1, 1], quadratic

    return {
        "common_quartic_radical_in_full_coordinates": [
            "x2",
            "x3",
            "x5",
            "x6",
            "x7",
        ],
        "free_coordinates": ["x0", "x1", "x4"],
        "linear_sheet": {
            "relation": "x2-7378*x5",
            "dimension": linear[0],
            "multiplicity": linear[1],
            "standard_basis_size": linear[2],
            "ideal_contained_in_candidate_radical": True,
            "x5_nilpotence_exponent_certified": 4,
        },
        "quadratic_sheet": {
            "relation": (
                "x2^2+13684*x2*x5+12028*x5^2"
            ),
            "relation_factor_count_over_finite_field": (
                quadratic[5]
            ),
            "dimension": quadratic[0],
            "multiplicity": quadratic[1],
            "standard_basis_size": quadratic[2],
            "ideal_contained_in_candidate_radical": True,
            "x5_nilpotence_exponent_certified": 5,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, default=32003)
    parser.add_argument("--timeout", type=int, default=180)
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    if arguments.prime != 32003:
        raise SystemExit(
            "the retained factor coefficients are specific to "
            "prime 32003"
        )
    moment_jets, allowed, forbidden = build_moment_jets(
        arguments.prime,
        4,
    )
    _, pivots, free, substitutions = formal_pivot_elimination(
        moment_jets,
        arguments.prime,
        4,
        False,
    )
    cubic = cubic_analysis(
        moment_jets,
        substitutions,
        arguments.prime,
        arguments.timeout,
    )
    print("completed cubic support decomposition", flush=True)
    quartic = quartic_sheet_analysis(
        moment_jets,
        substitutions,
        arguments.prime,
        arguments.timeout,
    )
    print("completed dominant-sheet quartic restrictions", flush=True)

    payload = {
        "format": "degree-four-q2-cubic-decomposition-v1",
        "scope": (
            "exact finite-field decomposition of the cubic normal "
            "moment jet and quartic restrictions on its dominant "
            "x6=0 sheets at one deterministic synchronized point"
        ),
        "prime": arguments.prime,
        "normalized_sym2_component": (
            "highest-weight raising matrix E"
        ),
        "allowed_point_coefficients": list(range(2, 11)),
        "allowed_weight_coordinates": [
            list(value) for value in allowed
        ],
        "forbidden_weight_coordinates": [
            list(value) for value in forbidden
        ],
        "linear_pivot_indices": pivots,
        "free_normal_indices": free,
        "cubic_decomposition": cubic,
        "dominant_component_quartic_restrictions": quartic,
        "off_axis_quartic_status": {
            "status": "timeout",
            "timeout_seconds": 240,
            "calculation": (
                "standard basis of the quartic moment ideal plus "
                "the cubic x6-saturation over "
                "F_32003(x1,x4,x6)"
            ),
            "inference": "none",
        },
        "scope_warning": (
            "This is a finite-field normal-jet computation.  It "
            "does not prove that the full moment-zero germ is the "
            "synchronized branch, does not treat the F2=0 boundary, "
            "and does not establish global integrality."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
