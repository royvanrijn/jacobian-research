#!/usr/bin/env python3
"""Sharded exact coefficient-torus screen for sparse supports.

After normalizing two nonzero coefficients of distinct diagonal weights,
a mixed-sign support of size ``n`` has ``n-2`` residual coordinates.  One
saturation equation restricts to its dense coefficient torus.  Exact
characteristic-zero msolve computations test the moments there.  For
``n=6``, every boundary is covered by the support-at-most-five certificate.
"""

from __future__ import annotations

import argparse
from itertools import combinations
import json
from math import factorial, gcd
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree33_sparse_six_support_screen.json"
)
POSITIONS = tuple((i, j) for i in range(4) for j in range(4))
H = sp.symbols("h")
RESIDUAL_SYMBOL_POOL = sp.symbols("z0:7")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--support-size", type=int, default=6, choices=(6, 7, 8, 9))
    parser.add_argument("--through", type=int, default=12)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--linear-algebra",
        type=int,
        choices=(1, 2, 42, 44),
        default=2,
    )
    parser.add_argument(
        "--rational-parametrization",
        action="store_true",
        help=(
            "ask msolve for an RUR of every complex solution; without "
            "this flag rational input reports only isolated real boxes"
        ),
    )
    parser.add_argument(
        "--no-saturation",
        action="store_true",
        help=(
            "diagnostic closed-coordinate-space solve; the default "
            "saturates every residual coefficient and is required for "
            "coefficient-torus census claims"
        ),
    )
    parser.add_argument(
        "--combine-input",
        type=Path,
        action="append",
        default=[],
        help=(
            "combine already computed contiguous shard artifacts; repeat "
            "once per input and do not use with --start/--limit"
        ),
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def normalize_support(
    support: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, int], ...]:
    for first, second in combinations(range(len(support)), 2):
        if (
            support[first][0] - support[first][1]
            != support[second][0] - support[second][1]
        ):
            return (
                support[first],
                support[second],
                *(
                    support[index]
                    for index in range(len(support))
                    if index not in (first, second)
                ),
            )
    raise ValueError("support has only one weight")


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first, *tail)


def restricted_moment(
    order: int,
    normalized: tuple[tuple[int, int], ...],
    residual_symbols: tuple[sp.Symbol, ...],
) -> sp.Poly:
    answer: dict[tuple[int, ...], int] = {}
    order_factorial = factorial(order)
    for counts in compositions(order, len(normalized)):
        x_degree = sum(
            count * position[0]
            for count, position in zip(counts, normalized, strict=True)
        )
        y_degree = sum(
            count * position[1]
            for count, position in zip(counts, normalized, strict=True)
        )
        if x_degree != y_degree:
            continue
        denominator = 1
        for count in counts:
            denominator *= factorial(count)
        coefficient = (
            order_factorial
            // denominator
            * factorial(3 * order - x_degree)
            * factorial(x_degree)
        )
        exponent = counts[2:]
        answer[exponent] = answer.get(exponent, 0) + coefficient
    content = 0
    for coefficient in answer.values():
        content = gcd(content, abs(coefficient))
    if content:
        answer = {
            exponent: coefficient // content
            for exponent, coefficient in answer.items()
        }
    return sp.Poly.from_dict(answer, residual_symbols, domain=sp.ZZ)


def verify_restricted_formula() -> dict[str, object]:
    """Compare the count recursion with direct bivariate expansion."""

    support = (
        (0, 1),
        (0, 2),
        (1, 2),
        (2, 0),
        (2, 3),
        (3, 1),
    )
    normalized = normalize_support(support)
    residual_symbols = RESIDUAL_SYMBOL_POOL[: len(support) - 2]
    x_symbol, y_symbol = sp.symbols("x y")
    source = sum(
        coefficient
        * x_symbol ** position[0]
        * y_symbol ** position[1]
        for coefficient, position in zip(
            (1, 1, *residual_symbols),
            normalized,
            strict=True,
        )
    )
    for order in range(1, 7):
        expanded = sp.expand(source**order)
        direct = 0
        for diagonal_degree in range(3 * order + 1):
            direct += (
                factorial(3 * order - diagonal_degree)
                * factorial(diagonal_degree)
                * expanded.coeff(x_symbol, diagonal_degree).coeff(
                    y_symbol,
                    diagonal_degree,
                )
            )
        direct_primitive = sp.Poly(
            direct,
            *residual_symbols,
            domain=sp.ZZ,
        ).primitive()[1]
        if direct_primitive != restricted_moment(
            order,
            normalized,
            residual_symbols,
        ):
            raise AssertionError(
                f"restricted moment mismatch on {support}, mu{order}"
            )
    return {
        "support": [list(position) for position in support],
        "orders": [1, 6],
        "passed": True,
    }


def msolve_expression(polynomial: sp.Poly) -> str:
    return str(polynomial.as_expr()).replace("**", "^")


def screen_support(
    support: tuple[tuple[int, int], ...],
    through: int,
    timeout: int,
    msolve: str,
    saturate: bool = True,
    threads: int = 1,
    rational_parametrization: bool = False,
    linear_algebra: int = 2,
) -> dict[str, object]:
    normalized = normalize_support(support)
    residual_symbols = RESIDUAL_SYMBOL_POOL[: len(support) - 2]
    moments = [
        restricted_moment(order, normalized, residual_symbols)
        for order in range(1, through + 1)
    ]
    saturation_monomial = "*".join(str(symbol) for symbol in residual_symbols)
    equations = [
        *(
            [f"h*{saturation_monomial}-1"]
            if saturate
            else []
        ),
        *(
            msolve_expression(moment)
            for moment in moments
            if not moment.is_zero
        ),
    ]
    started = time.monotonic()
    with tempfile.TemporaryDirectory(
        prefix="sic33-sparse-six-msolve-",
    ) as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        input_path.write_text(
            ",".join(
                (
                    *(("h",) if saturate else ()),
                    *(str(symbol) for symbol in residual_symbols),
                )
            )
            + "\n0\n"
            + ",\n".join(equations)
            + "\n",
            encoding="utf-8",
        )
        try:
            completed = subprocess.run(
                [
                    msolve,
                    "-f",
                    str(input_path),
                    "-o",
                    str(output_path),
                    "-t",
                    str(threads),
                    "-l",
                    str(linear_algebra),
                    "-v",
                    "0",
                    *(
                        ["-P", "1"]
                        if rational_parametrization
                        else []
                    ),
                ],
                text=True,
                capture_output=True,
                timeout=timeout if timeout else None,
                check=False,
            )
        except subprocess.TimeoutExpired as error:
            return {
                "support": [list(position) for position in support],
                "weights": [i - j for i, j in support],
                "seconds": round(time.monotonic() - started, 6),
                "status": "timeout",
                "stdout_tail": (error.stdout or "")[-2000:],
                "stderr_tail": (error.stderr or "")[-2000:],
            }
        result = (
            output_path.read_text(encoding="utf-8").strip()
            if output_path.exists()
            else ""
        )
    unit = result in ("[-1]:", "[-1]")
    record: dict[str, object] = {
        "support": [list(position) for position in support],
        "weights": [i - j for i, j in support],
        "seconds": round(time.monotonic() - started, 6),
        "status": (
            "excluded_on_coefficient_torus"
            if unit
            else (
                "msolve_nonempty"
                if completed.returncode == 0 and result
                else "msolve_failed"
            )
        ),
        "msolve": {
            "returncode": completed.returncode,
            "result": result if unit else None,
        },
    }
    if not unit:
        record["msolve"].update({
            "result_head": result[:2000],
            "result_tail": result[-2000:],
            "stdout_tail": completed.stdout[-2000:],
            "stderr_tail": completed.stderr[-2000:],
        })
        record["profiles"] = [
            {
                "order": order,
                "degree": (
                    -1 if moment.is_zero else int(moment.total_degree())
                ),
                "terms": len(moment.terms()),
            }
            for order, moment in enumerate(moments, start=1)
        ]
    return record


def main() -> None:
    arguments = parse_arguments()
    if (
        arguments.support_size != 6
        and arguments.output == DEFAULT_OUTPUT
    ):
        arguments.output = DEFAULT_OUTPUT.with_name(
            "two_pair_sic_bidegree33_sparse_"
            f"support{arguments.support_size}_screen.json"
        )
    if arguments.combine_input:
        if arguments.start or arguments.limit:
            raise ValueError("--combine-input cannot be used with --start/--limit")
        shards = [
            json.loads(path.read_text(encoding="utf-8"))
            for path in arguments.combine_input
        ]
        shards.sort(key=lambda shard: int(shard["start"]))
        if not shards:
            raise ValueError("no shards supplied")
        through_values = {int(shard["through"]) for shard in shards}
        support_sizes = {
            int(shard.get("support_size", 6))
            for shard in shards
        }
        mixed_counts = {
            int(shard["mixed_support_count"])
            for shard in shards
        }
        if (
            len(through_values) != 1
            or len(support_sizes) != 1
            or len(mixed_counts) != 1
        ):
            raise ValueError(
                "shards disagree on order, support size, or support count"
            )
        support_size = next(iter(support_sizes))
        if support_size != arguments.support_size:
            raise ValueError(
                "combined shard support size does not match --support-size"
            )
        expected_start = 0
        records = []
        source_shards = []
        for path, shard in sorted(
            zip(arguments.combine_input, [
                json.loads(path.read_text(encoding="utf-8"))
                for path in arguments.combine_input
            ], strict=True),
            key=lambda item: int(item[1]["start"]),
        ):
            start = int(shard["start"])
            stop = int(shard["stop"])
            if start != expected_start:
                raise ValueError(
                    f"noncontiguous shards: expected {expected_start}, got {start}"
                )
            if len(shard["records"]) != stop - start:
                raise ValueError(f"record count mismatch in {path}")
            if not shard["independent_formula_check"]["passed"]:
                raise ValueError(f"formula check failed in {path}")
            records.extend(shard["records"])
            source_shards.append({
                "path": str(path),
                "start": start,
                "stop": stop,
            })
            expected_start = stop
        mixed_support_count = next(iter(mixed_counts))
        complete = expected_start == mixed_support_count
        counts: dict[str, int] = {}
        for record in records:
            status = str(record["status"])
            counts[status] = counts.get(status, 0) + 1
        survivors = [
            record
            for record in records
            if record["status"] != "excluded_on_coefficient_torus"
        ]
        payload = {
            "calculation": (
                "two_pair_sic_bidegree33_sparse_six_support_screen"
            ),
            "scope": (
                "combined exact characteristic-zero coefficient-torus "
                "msolve screen"
                + (
                    "; all support boundaries are handled by the "
                    "support-at-most-five certificate"
                    if support_size == 6
                    else
                    "; boundary supports require the separately recorded "
                    "smaller-support results"
                )
            ),
            "through": next(iter(through_values)),
            "support_size": support_size,
            "mixed_support_count": mixed_support_count,
            "start": 0,
            "stop": expected_start,
            "counts": counts,
            "complete_mixed_enumeration": complete,
            "all_supports_excluded": complete and not survivors,
            "survivor_count": len(survivors),
            "survivors": survivors,
            "independent_formula_check": (
                shards[0]["independent_formula_check"]
            ),
            "source_shards": source_shards,
            "records": records,
        }
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        print(json.dumps({
            key: payload[key]
            for key in (
                "through",
                "support_size",
                "mixed_support_count",
                "start",
                "stop",
                "counts",
                "complete_mixed_enumeration",
                "all_supports_excluded",
                "survivor_count",
            )
        }, indent=2, sort_keys=True))
        return

    msolve = shutil.which("msolve")
    if msolve is None:
        raise RuntimeError("msolve is required")
    formula_check = verify_restricted_formula()
    mixed = [
        support
        for support in combinations(POSITIONS, arguments.support_size)
        if min(i - j for i, j in support) < 0
        and max(i - j for i, j in support) > 0
    ]
    stop = (
        len(mixed)
        if not arguments.limit
        else min(len(mixed), arguments.start + arguments.limit)
    )
    selected = mixed[arguments.start:stop]
    records = []
    for offset, support in enumerate(selected, start=arguments.start):
        record = screen_support(
            support,
            arguments.through,
            arguments.timeout,
            msolve,
            saturate=not arguments.no_saturation,
            threads=arguments.threads,
            rational_parametrization=arguments.rational_parametrization,
            linear_algebra=arguments.linear_algebra,
        )
        records.append(record)
        if record["status"] != "excluded_on_coefficient_torus":
            print(
                f"SUPPORT {offset} status={record['status']} "
                f"support={support}",
                flush=True,
            )
    counts: dict[str, int] = {}
    for record in records:
        status = str(record["status"])
        counts[status] = counts.get(status, 0) + 1
    payload = {
        "calculation": "two_pair_sic_bidegree33_sparse_six_support_screen",
        "scope": (
            "exact characteristic-zero "
            + (
                "coefficient-torus"
                if not arguments.no_saturation
                else "closed coefficient-space diagnostic"
            )
            + " msolve screen; sharded enumeration"
            + (
                "; full complex rational parametrization requested"
                if arguments.rational_parametrization
                else "; default rational real-root output"
            )
            + (
                ", with all support boundaries handled by the "
                "support-at-most-five certificate"
                if arguments.support_size == 6
                else
                "; boundary supports require the separately recorded "
                "smaller-support results"
            )
        ),
        "through": arguments.through,
        "support_size": arguments.support_size,
        "mixed_support_count": len(mixed),
        "start": arguments.start,
        "stop": stop,
        "counts": counts,
        "complete_mixed_enumeration": (
            arguments.start == 0 and stop == len(mixed)
        ),
        "all_supports_excluded": (
            arguments.start == 0
            and stop == len(mixed)
            and all(
                record["status"] == "excluded_on_coefficient_torus"
                for record in records
            )
        ),
        "survivor_count": sum(
            record["status"] != "excluded_on_coefficient_torus"
            for record in records
        ),
        "independent_formula_check": formula_check,
        "records": records,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        key: payload[key]
        for key in (
            "through",
            "support_size",
            "mixed_support_count",
            "start",
            "stop",
            "counts",
            "complete_mixed_enumeration",
            "all_supports_excluded",
            "survivor_count",
        )
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
