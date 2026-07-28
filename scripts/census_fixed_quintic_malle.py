#!/usr/bin/env python3
"""Exact bounded census for Malle-style statistics in one Keller family.

The census enumerates primitive projective targets [W:P:B:C] with W > 0,
P != 0, and max(|W|,|P|,|B|,|C|) <= H.  It quotients the exact symmetry
(B,C) ~ (-B,-C), classifies every irreducible quintic with PARI/GP,
computes field and integral-generator discriminants, and deduplicates fields
with exact ``nfisisom`` tests inside discriminant/signature/Galois buckets.

This is a bounded experiment.  It is not evidence for an asymptotic unless
the height and discriminant cutoffs, presentation multiplicities, and cusp
effects are analyzed separately.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import reduce
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_TEMPLATE = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "fixed_quintic_malle_census_h{bound}.json"
)


@dataclass(frozen=True, order=True)
class Target:
    w: int
    pi: int
    b: int
    c: int

    @property
    def height(self) -> int:
        return max(self.w, abs(self.pi), abs(self.b), abs(self.c))

    @property
    def primitive_inverse_coefficients(self) -> tuple[int, ...]:
        coefficients = (
            self.pi**5,
            0,
            -5 * self.pi * self.w**4,
            -2 * self.b * self.w**4,
            4 * self.w**5,
            -2 * self.c * self.w**4,
        )
        content = reduce(math.gcd, (abs(value) for value in coefficients))
        primitive = tuple(value // content for value in coefficients)
        if primitive[0] < 0:
            primitive = tuple(-value for value in primitive)
        return primitive

    @property
    def monic_integral_coefficients(self) -> tuple[int, ...]:
        """Polynomial of a*theta for primitive inverse leading coefficient a."""

        primitive = self.primitive_inverse_coefficients
        leading = primitive[0]
        assert leading > 0
        return (1,) + tuple(
            primitive[index] * leading ** (index - 1)
            for index in range(1, 6)
        )


@dataclass(frozen=True)
class PresentationRecord:
    target: Target
    group: str
    field_discriminant: int
    real_roots: int
    generator_discriminant: int
    generator_index: int
    reduced_polynomial: str


class UnionFind:
    def __init__(self, values: list[str]) -> None:
        self.parent = {value: value for value in values}

    def find(self, value: str) -> str:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[right_root] = left_root


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bound", type=int, default=10)
    parser.add_argument("--gp", default="gp")
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--exceptional-only",
        action="store_true",
        help=(
            "scan the full box but initialize number fields only for the "
            "exact square-discriminant and finite-prime F20 screens"
        ),
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=100,
        help="number of targets per isolated PARI/GP process",
    )
    parser.add_argument(
        "--batch-timeout",
        type=int,
        default=15,
        help="seconds before a PARI batch is killed and recursively bisected",
    )
    parser.add_argument(
        "--gp-stack",
        default="256M",
        help="PARI stack allocated to each isolated process",
    )
    parser.add_argument(
        "--no-write",
        action="store_true",
        help="run and print the summary without writing the JSON artifact",
    )
    return parser.parse_args()


def projective_targets(bound: int):
    """Enumerate primitive targets modulo (B,C) -> (-B,-C)."""

    for w in range(1, bound + 1):
        for pi, b, c in itertools.product(
            range(-bound, bound + 1),
            range(0, bound + 1),
            range(-bound, bound + 1),
        ):
            if not pi or (b == 0 and c < 0):
                continue
            if math.gcd(w, abs(pi), b, abs(c)) != 1:
                continue
            yield Target(w, pi, b, c)


def discriminant_square_class(target: Target) -> int:
    """Discriminant after removing the nonzero square (4*P^4*W^8)^2."""

    w, p, b, c = target.w, target.pi, target.b, target.c
    return (
        432 * b**5 * c * p**2 * w**8
        - 432 * b**4 * p**2 * w**10
        + 12600 * b**3 * c * p**3 * w**9
        - 2000 * b**3 * c * w**12
        + 9000 * b**2 * c**2 * p**7 * w**5
        + 20625 * b**2 * c**2 * p**4 * w**8
        - 11520 * b**2 * p**3 * w**11
        + 2000 * b**2 * w**14
        + 18750 * b * c**3 * p**8 * w**4
        - 25600 * b * c * p**7 * w**7
        + 56000 * b * c * p**4 * w**10
        - 45000 * b * c * p * w**13
        + 3125 * c**4 * p**12
        - 40000 * c**2 * p**8 * w**6
        + 112500 * c**2 * p**5 * w**9
        - 84375 * c**2 * p**2 * w**12
        + 16384 * p**7 * w**9
        - 51200 * p**4 * w**12
        + 40000 * p * w**15
    )


def is_square(value: int) -> bool:
    if value < 0:
        return False
    root = math.isqrt(value)
    return root * root == value


def exceptional_screen(target: Target, square_class: int) -> bool:
    """Necessary global/local screen for A5, C5, D5, and F20."""

    if square_class <= 0:
        return False
    if is_square(square_class):
        return True
    from search_universal_quintic_calculator import (  # local script module
        SCREEN_PRIMES,
        local_f20_allowed,
    )

    return all(
        local_f20_allowed(target.primitive_inverse_coefficients, prime)
        for prime in SCREEN_PRIMES
    )


def polynomial_text(coefficients: tuple[int, ...], variable: str = "x") -> str:
    degree = len(coefficients) - 1
    terms = [
        f"({coefficient})*{variable}^{degree - index}"
        for index, coefficient in enumerate(coefficients)
        if coefficient
    ]
    return "+".join(terms) if terms else "0"


def normalize_group(pari_name: str) -> str:
    if pari_name == "S5":
        return "S5"
    if pari_name == "A5":
        return "A5"
    if pari_name.startswith("C(5)"):
        return "C5"
    if pari_name.startswith("D(5)"):
        return "D5"
    if pari_name.startswith("F(5)"):
        return "F20"
    raise ValueError(f"unrecognized PARI quintic group: {pari_name!r}")


def gp_command(gp: str, gp_stack: str) -> list[str]:
    return [gp, "-s", gp_stack, "-q"]


def gp_version(gp: str, gp_stack: str) -> str:
    completed = subprocess.run(
        gp_command(gp, gp_stack),
        input='print(version());quit\n',
        text=True,
        capture_output=True,
        check=True,
    )
    return completed.stdout.strip()


def parse_classification_output(
    output: str, targets: list[Target]
) -> list[PresentationRecord]:
    records: list[PresentationRecord] = []
    for line in output.splitlines():
        fields = line.split("|", 5)
        if len(fields) != 6:
            raise ValueError(f"unexpected PARI output: {line!r}")
        index_text, pari_group, disc_text, roots_text, poly_disc_text, reduced = (
            fields
        )
        target = targets[int(index_text)]
        field_discriminant = int(disc_text)
        generator_discriminant = int(poly_disc_text)
        quotient, remainder = divmod(
            abs(generator_discriminant), abs(field_discriminant)
        )
        if remainder:
            raise AssertionError(
                f"nonintegral discriminant quotient at target {target}"
            )
        generator_index = math.isqrt(quotient)
        if generator_index * generator_index != quotient:
            raise AssertionError(
                f"nonsquare discriminant quotient at target {target}"
            )
        records.append(
            PresentationRecord(
                target=target,
                group=normalize_group(pari_group),
                field_discriminant=field_discriminant,
                real_roots=int(roots_text),
                generator_discriminant=generator_discriminant,
                generator_index=generator_index,
                reduced_polynomial=reduced.strip(),
            )
        )
    return records


def classify_targets(
    targets: list[Target],
    gp: str,
    gp_stack: str,
    batch_size: int,
    batch_timeout: int,
) -> tuple[list[PresentationRecord], str]:
    """Classify every irreducible target in isolated PARI processes."""

    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    if batch_timeout < 1:
        raise ValueError("batch_timeout must be positive")
    records: list[PresentationRecord] = []
    stderr_chunks: list[str] = []
    batch_count = math.ceil(len(targets) / batch_size)

    def run_range(start: int, stop: int) -> list[PresentationRecord]:
        lines = []
        for index in range(start, stop):
            target = targets[index]
            polynomial = polynomial_text(target.monic_integral_coefficients)
            lines.append(
                f"f={polynomial};"
                "if(polisirreducible(f),"
                "nf=nfinit(f);r=polredabs(f);"
                f'print("{index}|",polgalois(f)[4],"|",nf.disc,"|",'
                'polsturm(f),"|",poldisc(f),"|",r))'
            )

        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".gp", encoding="utf-8"
        ) as source:
            source.write("\n".join(lines))
            source.write("\nquit\n")
            source.flush()
            try:
                completed = subprocess.run(
                    [*gp_command(gp, gp_stack), source.name],
                    text=True,
                    capture_output=True,
                    check=True,
                    timeout=batch_timeout if stop - start > 1 else 300,
                )
            except subprocess.TimeoutExpired:
                if stop - start == 1:
                    raise RuntimeError(
                        f"PARI/GP exceeded 300 seconds at target index {start}"
                    )
                midpoint = start + (stop - start) // 2
                print(
                    f"bisecting slow PARI range [{start},{stop}) into "
                    f"[{start},{midpoint}) and [{midpoint},{stop})",
                    flush=True,
                )
                return run_range(start, midpoint) + run_range(midpoint, stop)
        stderr = completed.stderr.strip()
        if stderr:
            stderr_chunks.append(stderr)
        if "***" in stderr:
            if stop - start > 1:
                midpoint = start + (stop - start) // 2
                print(
                    f"bisecting failed PARI range [{start},{stop})",
                    flush=True,
                )
                return run_range(start, midpoint) + run_range(midpoint, stop)
            raise RuntimeError(
                f"PARI/GP failed at target index {start}:\n{stderr}"
            )
        return parse_classification_output(completed.stdout, targets)

    for batch_index, start in enumerate(range(0, len(targets), batch_size), 1):
        stop = min(start + batch_size, len(targets))
        records.extend(run_range(start, stop))
        if batch_index % 10 == 0 or batch_index == batch_count:
            print(
                f"completed PARI batch {batch_index}/{batch_count}",
                flush=True,
            )
    return records, "\n".join(stderr_chunks)


def exact_isomorphism_edges(
    records: list[PresentationRecord], gp: str, gp_stack: str
) -> set[tuple[str, str]]:
    """Compare distinct reduced-polynomial keys in every viable bucket."""

    buckets: dict[tuple[str, int, int], set[str]] = defaultdict(set)
    for record in records:
        buckets[
            (record.group, record.field_discriminant, record.real_roots)
        ].add(record.reduced_polynomial)

    pairs: list[tuple[str, str]] = []
    for polynomials in buckets.values():
        ordered = sorted(polynomials)
        pairs.extend(itertools.combinations(ordered, 2))
    if not pairs:
        return set()

    lines = []
    for index, (left, right) in enumerate(pairs):
        lines.append(
            f'print("{index}|",nfisisom({left},{right})!=0)'
        )
    completed = subprocess.run(
        gp_command(gp, gp_stack),
        input="\n".join(lines) + "\nquit\n",
        text=True,
        capture_output=True,
        check=True,
    )
    edges: set[tuple[str, str]] = set()
    for line in completed.stdout.splitlines():
        index_text, result = line.split("|", 1)
        if int(result):
            edges.add(pairs[int(index_text)])
    return edges


def field_classes(
    records: list[PresentationRecord], isomorphism_edges: set[tuple[str, str]]
) -> list[dict[str, object]]:
    keys = sorted({record.reduced_polynomial for record in records})
    union_find = UnionFind(keys)
    for left, right in isomorphism_edges:
        union_find.union(left, right)

    grouped: dict[
        tuple[str, int, int, str], list[PresentationRecord]
    ] = defaultdict(list)
    for record in records:
        grouped[
            (
                record.group,
                record.field_discriminant,
                record.real_roots,
                union_find.find(record.reduced_polynomial),
            )
        ].append(record)

    classes = []
    for (
        group,
        field_discriminant,
        real_roots,
        _,
    ), presentations in sorted(grouped.items()):
        presentations.sort(key=lambda item: item.target)
        reduced_polynomials = sorted(
            {item.reduced_polynomial for item in presentations}
        )
        classes.append(
            {
                "group": group,
                "field_discriminant": field_discriminant,
                "absolute_field_discriminant": abs(field_discriminant),
                "real_roots": real_roots,
                "reduced_polynomial": reduced_polynomials[0],
                "equivalent_reduced_polynomials": reduced_polynomials,
                "first_target_height": min(
                    item.target.height for item in presentations
                ),
                "presentation_count": len(presentations),
                "targets": [
                    {
                        "projective": [
                            item.target.w,
                            item.target.pi,
                            item.target.b,
                            item.target.c,
                        ],
                        "height": item.target.height,
                        "generator_discriminant": item.generator_discriminant,
                        "generator_index": item.generator_index,
                    }
                    for item in presentations
                ],
            }
        )
    return classes


def quantiles(values: list[int]) -> dict[str, int]:
    if not values:
        return {}
    ordered = sorted(values)

    def at(numerator: int, denominator: int) -> int:
        index = (len(ordered) - 1) * numerator // denominator
        return ordered[index]

    return {
        "min": ordered[0],
        "q25": at(1, 4),
        "median": at(1, 2),
        "q75": at(3, 4),
        "max": ordered[-1],
    }


def build_summary(
    *,
    bound: int,
    targets: list[Target],
    zero_discriminant_count: int,
    records: list[PresentationRecord],
    classes: list[dict[str, object]],
    isomorphism_edges: set[tuple[str, str]],
) -> dict[str, object]:
    groups = ("S5", "A5", "F20", "D5", "C5")
    summary_by_group = {}
    height_cumulative = []
    for height in range(1, bound + 1):
        height_cumulative.append(
            {
                "height": height,
                "presentations": {
                    group: sum(
                        record.group == group
                        and record.target.height <= height
                        for record in records
                    )
                    for group in groups
                },
                "field_classes": {
                    group: sum(
                        field["group"] == group
                        and field["first_target_height"] <= height
                        for field in classes
                    )
                    for group in groups
                },
            }
        )

    for group in groups:
        group_records = [record for record in records if record.group == group]
        group_classes = [field for field in classes if field["group"] == group]
        multiplicities = Counter(
            int(field["presentation_count"]) for field in group_classes
        )
        signatures_presentations = Counter(
            record.real_roots for record in group_records
        )
        signatures_fields = Counter(
            int(field["real_roots"]) for field in group_classes
        )
        summary_by_group[group] = {
            "presentations": len(group_records),
            "field_classes": len(group_classes),
            "presentation_to_field_ratio": (
                len(group_records) / len(group_classes)
                if group_classes
                else None
            ),
            "multiplicity_histogram": {
                str(key): value for key, value in sorted(multiplicities.items())
            },
            "maximum_presentation_multiplicity": max(multiplicities, default=0),
            "signatures_by_presentations": {
                str(key): value
                for key, value in sorted(signatures_presentations.items())
            },
            "signatures_by_fields": {
                str(key): value
                for key, value in sorted(signatures_fields.items())
            },
            "absolute_field_discriminant_quantiles": quantiles(
                [
                    int(field["absolute_field_discriminant"])
                    for field in group_classes
                ]
            ),
            "minimum_target_height": min(
                (
                    int(field["first_target_height"])
                    for field in group_classes
                ),
                default=None,
            ),
        }

    return {
        "target_bound": bound,
        "symmetry_convention": (
            "W>0 and one representative of (B,C)~(-B,-C): "
            "B>0 or B=0,C>=0"
        ),
        "primitive_projective_targets": len(targets),
        "zero_discriminant_targets": zero_discriminant_count,
        "connected_presentations": len(records),
        "reducible_squarefree_presentations": (
            len(targets) - zero_discriminant_count - len(records)
        ),
        "exact_cross_key_isomorphisms": len(isomorphism_edges),
        "by_group": summary_by_group,
        "height_cumulative": height_cumulative,
    }


def main() -> None:
    args = parse_args()
    if args.bound < 1:
        raise SystemExit("--bound must be positive")
    if args.output:
        output = args.output
    elif args.exceptional_only:
        output = (
            ROOT
            / "artifacts"
            / "generated-results"
            / f"fixed_quintic_malle_exceptional_census_h{args.bound}.json"
        )
    else:
        output = Path(str(DEFAULT_OUTPUT_TEMPLATE).format(bound=args.bound))

    scanned_target_count = 0
    zero_discriminant_count = 0
    targets: list[Target] = []
    for target in projective_targets(args.bound):
        scanned_target_count += 1
        square_class = discriminant_square_class(target)
        if square_class == 0:
            zero_discriminant_count += 1
        if not args.exceptional_only or exceptional_screen(
            target, square_class
        ):
            targets.append(target)
    print(
        f"scanned {scanned_target_count} primitive targets through height "
        f"{args.bound}; classifying {len(targets)} with PARI/GP",
        flush=True,
    )
    records, pari_stderr = classify_targets(
        targets,
        args.gp,
        args.gp_stack,
        args.batch_size,
        args.batch_timeout,
    )
    print(
        f"found {len(records)} connected presentations; "
        "running exact field-isomorphism audit",
        flush=True,
    )
    isomorphism_edges = exact_isomorphism_edges(
        records, args.gp, args.gp_stack
    )
    classes = field_classes(records, isomorphism_edges)
    summary = build_summary(
        bound=args.bound,
        targets=targets,
        zero_discriminant_count=zero_discriminant_count,
        records=records,
        classes=classes,
        isomorphism_edges=isomorphism_edges,
    )
    if args.exceptional_only:
        summary["primitive_projective_targets"] = scanned_target_count
        summary["exceptional_screen_candidates"] = len(targets)
        summary["exceptional_screen_nonconnected"] = (
            len(targets) - len(records)
        )
        summary["reducible_squarefree_presentations"] = None
        summary["classification_scope"] = (
            "complete scan of the target box; PARI initialization only for "
            "square discriminants or targets passing all declared F20 "
            "finite-prime screens"
        )
    artifact = {
        "format": "fixed-quintic-malle-census-v1",
        "status": (
            "bounded exact computation; not an asymptotic theorem or "
            "evidence independent of PARI/GP"
        ),
        "family": (
            "E_(Pi,B,C)(S)=Pi^5*S^5-5*Pi*S^3-2*B*S^2+4*S-2*C"
        ),
        "command": " ".join(
            [
                ".venv/bin/python",
                "scripts/census_fixed_quintic_malle.py",
                "--bound",
                str(args.bound),
                *(["--exceptional-only"] if args.exceptional_only else []),
            ]
        ),
        "software": {
            "python": sys.version.split()[0],
            "pari_gp": gp_version(args.gp, args.gp_stack),
            "pari_stack": args.gp_stack,
        },
        "pari_stderr": pari_stderr,
        "summary": summary,
        "field_classes": classes,
    }

    encoded = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    digest = hashlib.sha256(encoded.encode("utf-8")).hexdigest()
    if not args.no_write:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded, encoding="utf-8")
        try:
            hash_target = output.relative_to(ROOT)
        except ValueError:
            hash_target = output
        output.with_suffix(output.suffix + ".sha256").write_text(
            f"{digest}  {hash_target}\n", encoding="utf-8"
        )
        print(f"wrote {output}")
        print(f"sha256 {digest}")

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
