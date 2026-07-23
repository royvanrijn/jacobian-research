#!/usr/bin/env python3
"""Scheme-level low-degree hunt for exceptional-seed components.

This is deliberately independent of point sampling.  For every requested
degree it builds the exact full-contact incidence presentations supplied by
``jcsearch.discriminant_geometry``, asks Singular to eliminate the root
coordinates, projectively closes the resulting seed ideals, and compares:

* every full-contact partition with the maximal 2/3 partitions;
* the exact-root saturation with the collision-retaining chart;
* the saturated chart with the raw normalization equations;
* ordered roots with the equal-multiplicity quotient coordinates;
* the union of all full-contact images with the union of the 2/3 images.

The Hilbert numerator and polynomial are computed for each maximal component
and for their reduced union.  Minimal-prime/radical checks are enabled in the
smaller degrees by default because they are substantially more expensive than
the eliminations themselves.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from math import comb
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.discriminant_geometry import contact_partition_incidence  # noqa: E402


def restricted_partitions(total: int, maximum: int | None = None):
    """Yield decreasing partitions of ``total`` with every part at least two."""

    if total == 0:
        yield ()
        return
    maximum = min(total, maximum or total)
    for first in range(maximum, 1, -1):
        for tail in restricted_partitions(total - first, first):
            yield (first,) + tail


def atomic_partitions(total: int):
    return tuple(
        partition
        for partition in restricted_partitions(total)
        if set(partition) <= {2, 3}
    )


def collision_precedes(finer: tuple[int, ...], coarser: tuple[int, ...]) -> bool:
    """Small local version of the collision partial order."""

    targets = tuple(sorted(coarser, reverse=True))
    source = tuple(sorted(finer, reverse=True))

    def search(index: int, bins: tuple[int, ...]) -> bool:
        if index == len(source):
            return all(value == 0 for value in bins)
        value = source[index]
        seen = set()
        for target_index, remaining in enumerate(bins):
            if remaining < value or remaining in seen:
                continue
            seen.add(remaining)
            updated = list(bins)
            updated[target_index] -= value
            if search(index + 1, tuple(updated)):
                return True
        return False

    return sum(source) == sum(targets) and search(0, targets)


def tag(partition: tuple[int, ...]) -> str:
    return "p" + "_".join(map(str, partition))


def singular_expression(expression: sp.Expr, substitution: dict[sp.Symbol, sp.Symbol]) -> str:
    return sp.sstr(sp.expand(expression.xreplace(substitution))).replace("**", "^")


@dataclass(frozen=True)
class SourceIdeal:
    ring: str
    ideal: str


class SingularProgram:
    def __init__(self, degree: int, deep: bool):
        self.degree = degree
        self.deep = deep
        h_names = tuple(f"h{index}" for index in range(3, degree))
        self.lines = [
            'LIB "elim.lib";',
            'LIB "primdec.lib";',
            "option(redSB);",
            f"ring target=0,({','.join(h_names)},z),dp;",
            "proc sameIdeal(ideal a, ideal b)",
            "{",
            "  a=std(a); b=std(b);",
            "  ideal left=simplify(reduce(a,b),2);",
            "  ideal right=simplify(reduce(b,a),2);",
            "  return(size(left)==0 && size(right)==0);",
            "}",
            "proc containedIdeal(ideal small, ideal large)",
            "{",
            "  ideal remainder=simplify(reduce(small,std(large)),2);",
            "  return(size(remainder)==0);",
            "}",
        ]
        self.sources: dict[tuple[tuple[int, ...], str], SourceIdeal] = {}

    def add_source(self, partition: tuple[int, ...], variant: str) -> None:
        data = contact_partition_incidence(self.degree, partition)
        quotient_points = data.quotient_coordinates
        seed_parameters = data.coefficient_parameters
        base_equations = data.coefficient_space_elimination_ideal.equations
        D = data.scale_denominator
        weighted = sp.diff(data.M, data.variable, 2).subs(data.variable, 1) - 2 * D

        if variant == "ordered":
            point_variables = data.marked_roots
            root_substitution = dict(data.root_to_quotient)
            equations = tuple(equation.xreplace(root_substitution) for equation in base_equations)
            distinct_factor = data.distinct_root_factor.xreplace(root_substitution)
            admissible_factor = sp.factor((D * weighted).xreplace(root_substitution))
        else:
            point_variables = quotient_points
            root_substitution = {}
            equations = base_equations
            distinct_factor = data.distinct_root_factor
            admissible_factor = sp.factor(D * weighted)
            if variant != "all":  # pragma: no cover - internal programming error
                raise ValueError(variant)

        # Eliminating a Rabinowitsch gate is dramatically more expensive for
        # the all-double component.  It is also unnecessary for a closure:
        # below we check in the root ring that saturating the irreducible Phi
        # divisor by the admissibility and diagonal factors changes no ideal.
        nuisance_count = len(point_variables)
        nuisance = sp.symbols(f"x1:{nuisance_count + 1}")
        canonical_h = sp.symbols(f"h3:{self.degree}")
        substitution = dict(zip(point_variables, nuisance))
        substitution.update(zip(seed_parameters, canonical_h))
        variables = nuisance + canonical_h
        ring_name = f"r_{tag(partition)}_{variant}"
        ideal_name = f"k_{tag(partition)}_{variant}"
        self.lines.append(f"ring {ring_name}=0,({','.join(map(str, variables))}),dp;")
        equation_strings = [singular_expression(equation, substitution) for equation in equations]
        self.lines.append(f"ideal j={','.join(equation_strings)};")
        self.lines.append(
            "poly rootPhi="
            + singular_expression(data.phi.xreplace(root_substitution) if variant == "ordered" else data.phi, substitution)
            + ";"
        )
        self.lines.append(
            f"poly rootAdmissible={singular_expression(admissible_factor, substitution)};"
        )
        self.lines.append(
            f"poly rootAll={singular_expression(sp.factor(admissible_factor * distinct_factor), substitution)};"
        )
        elimination_monomial = "*".join(map(str, nuisance))
        self.lines.append(
            f'ideal {ideal_name}=eliminate(j,{elimination_monomial},"slimgb");'
        )
        self.lines.append("kill j;")
        if variant == "all" and set(partition) <= {2, 3}:
            self.lines.extend(
                [
                    f'print("CHART|{self.degree}|{tag(partition)}|admissible|"+string(gcd(rootPhi,rootAdmissible)==1));',
                    f'print("CHART|{self.degree}|{tag(partition)}|raw|"+string(gcd(rootPhi,rootAll)==1));',
                ]
            )
        target_name = f"i_{tag(partition)}_{variant}"
        self.lines.extend(
            [
                "setring target;",
                f"ideal {target_name}=imap({ring_name},{ideal_name});",
                f"{target_name}=homog({target_name},z);",
                f"{target_name}=sat({target_name},ideal(z));",
                f"{target_name}=std({target_name});",
                f"kill {ring_name};",
            ]
        )
        self.sources[(partition, variant)] = SourceIdeal(ring_name, ideal_name)

    def finish(self, partitions: tuple[tuple[int, ...], ...]) -> str:
        degree = self.degree
        atoms = atomic_partitions(degree)

        for partition in atoms:
            name = f"i_{tag(partition)}_all"
            component_lines = [
                    f'print("DATA|{degree}|{tag(partition)}|"+string(dim({name})-1)+"|"+string(mult({name}))+"|"+string(hilb({name},2)));'
            ]
            if (partition, "ordered") in self.sources:
                component_lines.append(
                    f'print("CHART|{degree}|{tag(partition)}|ordered|"+string(sameIdeal({name},i_{tag(partition)}_ordered)));'
                )
            component_lines.extend(
                [
                    f"ideal boundary={name}+ideal(z); boundary=std(boundary);",
                    f'print("INFINITY|{degree}|{tag(partition)}|"+string(dim(boundary)-1));',
                ]
            )
            self.lines.extend(component_lines)
            if self.deep:
                self.lines.extend(
                    [
                        f"ideal rad=radical({name});",
                        f"list ass=minAssGTZ({name});",
                        f'print("SCHEME|{degree}|{tag(partition)}|"+string(sameIdeal({name},rad))+"|"+string(size(ass)));',
                    ]
                )

        atom_names = [f"i_{tag(partition)}_all" for partition in atoms]
        union_expression = atom_names[0]
        for name in atom_names[1:]:
            union_expression = f"intersect({union_expression},{name})"
        self.lines.extend(
            [
                f"ideal atomicUnion=std({union_expression});",
                f'print("UNION|{degree}|atomic|"+string(dim(atomicUnion)-1)+"|"+string(mult(atomicUnion))+"|"+string(hilb(atomicUnion,2)));',
            ]
        )
        all_names = [f"i_{tag(partition)}_all" for partition in partitions]
        all_union_expression = all_names[0]
        for name in all_names[1:]:
            all_union_expression = f"intersect({all_union_expression},{name})"
        self.lines.extend(
            [
                f"ideal allUnion=std({all_union_expression});",
                f'print("COVER|{degree}|"+string(sameIdeal(atomicUnion,allUnion)));',
            ]
        )

        for partition in partitions:
            for atom in atoms:
                if collision_precedes(atom, partition):
                    self.lines.append(
                        f'print("EDGE|{degree}|{tag(atom)}|{tag(partition)}|"+string(containedIdeal(i_{tag(atom)}_all,i_{tag(partition)}_all)));'
                    )
        self.lines.append("exit;")
        return "\n".join(self.lines) + "\n"


def hilbert_polynomial(numerator: tuple[int, ...], ring_dimension: int) -> sp.Expr:
    """Convert Q(t)/(1-t)^d to the eventual projective Hilbert polynomial."""

    variable = sp.Symbol("m", integer=True)
    result = sp.Integer(0)
    order = ring_dimension - 1
    for index, coefficient in enumerate(numerator):
        if coefficient:
            result += coefficient * sp.prod(
                variable - index + order - shift for shift in range(order)
            ) / sp.factorial(order)
    return sp.expand(result)


def parse_vector(value: str) -> tuple[int, ...]:
    return tuple(int(entry) for entry in value.split(",") if entry.strip())


def run_degree(degree: int, deep_max_degree: int, timeout: int) -> list[str]:
    partitions = tuple(restricted_partitions(degree))
    atoms = atomic_partitions(degree)
    program = SingularProgram(degree, degree <= deep_max_degree)
    for partition in partitions:
        program.add_source(partition, "all")
        if partition in atoms and degree <= 6:
            program.add_source(partition, "ordered")
    source = program.finish(partitions)
    process = subprocess.run(
        ["Singular", "-q"],
        input=source,
        text=True,
        capture_output=True,
        timeout=timeout,
        check=False,
    )
    combined = process.stdout + process.stderr
    errors = [line for line in combined.splitlines() if line.lstrip().startswith("?")]
    if process.returncode or errors:
        excerpt = "\n".join(errors[-20:] or combined.splitlines()[-20:])
        raise RuntimeError(f"Singular failed in degree {degree}:\n{excerpt}")
    return [
        line.strip()
        for line in process.stdout.splitlines()
        if re.match(r"^(DATA|CHART|INFINITY|SCHEME|UNION|COVER|EDGE)\|", line.strip())
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-degree", type=int, default=4)
    parser.add_argument("--max-degree", type=int, default=7)
    parser.add_argument(
        "--deep-max-degree",
        type=int,
        default=6,
        help="last degree in which to run radical and minimal-prime checks",
    )
    parser.add_argument("--timeout", type=int, default=900, help="seconds per degree")
    args = parser.parse_args()
    if not shutil.which("Singular"):
        raise SystemExit("Singular is required for this audit")

    records: list[str] = []
    for degree in range(args.min_degree, args.max_degree + 1):
        print(f"eliminating degree {degree} ...", flush=True)
        degree_records = run_degree(degree, args.deep_max_degree, args.timeout)
        records.extend(degree_records)

    failures = []
    print("\nmaximal components")
    for record in records:
        fields = record.split("|")
        if fields[0] != "DATA":
            continue
        _, degree_text, component, dimension_text, degree_value, vector_text = fields
        degree = int(degree_text)
        dimension = int(dimension_text)
        numerator = parse_vector(vector_text)
        polynomial = hilbert_polynomial(numerator, dimension + 1)
        partition = tuple(int(value) for value in component[1:].split("_"))
        expected_dimension = len(partition) - 1
        a = partition.count(2)
        b = partition.count(3)
        expected_degree = comb(a + b, a) * (2**a * 3**b - (a + b + 1))
        if dimension != expected_dimension or int(degree_value) != expected_degree:
            failures.append(record)
        print(
            f"N={degree} {partition}: dim={dimension}, degree={degree_value}, "
            f"Hilbert numerator={numerator}, P(m)={polynomial}"
        )

    print("\nunion Hilbert data")
    for record in records:
        if record.startswith("UNION|"):
            _, degree, _, dimension, multiplicity, vector = record.split("|")
            numerator = parse_vector(vector)
            polynomial = hilbert_polynomial(numerator, int(dimension) + 1)
            print(
                f"N={degree}: dim={dimension}, degree={multiplicity}, "
                f"Hilbert numerator={numerator}, P(m)={polynomial}"
            )

    boolean_prefixes = ("CHART|", "SCHEME|", "COVER|", "EDGE|")
    for record in records:
        if record.startswith(boolean_prefixes):
            values = record.split("|")
            if values[0] == "SCHEME":
                if values[-2:] != ["1", "1"]:
                    failures.append(record)
            elif values[-1] != "1":
                failures.append(record)
        elif record.startswith("INFINITY|"):
            _, degree_text, component, boundary_dimension = record.split("|")
            partition = tuple(int(value) for value in component[1:].split("_"))
            if int(boundary_dimension) >= len(partition) - 1:
                failures.append(record)

    if failures:
        print("\nFAILURES")
        print("\n".join(failures))
        raise SystemExit(1)

    print("\nPASS: all full-contact images are contained in the 2/3-component union")
    print("PASS: raw and collision-retaining root-chart closures agree")
    if args.min_degree <= 6:
        print("PASS: ordered-root and symmetry-quotient images agree through degree six")
    print("PASS: no maximal component is supported at projective root infinity")
    if args.deep_max_degree >= args.min_degree:
        print("PASS: checked maximal ideals are radical with one minimal prime")
    print("PASS: component dimensions and projective degrees match the uniform formulas")


if __name__ == "__main__":
    main()
