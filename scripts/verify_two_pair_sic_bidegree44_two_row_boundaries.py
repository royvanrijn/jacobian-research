#!/usr/bin/env python3
"""Classify coordinate boundaries of the separated-row quartic SIC chart.

The eight possible nonzero coefficients are

    {0} x {1,2,3,4}  union  {4} x {0,1,2,3}.

Their contraction-preserving diagonal-torus weights ``i-j`` are precisely
``-4,-3,-2,-1,1,2,3,4``.  Transposition and simultaneous coordinate
reversal both negate these weights, so support types are subsets modulo
``W -> -W``.  On every support of size at least two, overall scaling and
the diagonal torus normalize the first two coefficients to one (over the
algebraic closure); all remaining coefficients are localized away from
zero.  The pure moments are then formed directly over QQ and msolve decides
the localized ideals exactly.

For every mixed proper support this checker proves a sharp finite cutoff:
the moment ideal through the preceding order is nonunit, while adjoining
the moment at the recorded cutoff gives the unit ideal.  One-sided supports
have no zero-weight monomial in any positive power, so every pure moment
vanishes identically.  Supports whose two nonzero rows have identical
column support are additionally split into their rank-one closed locus and
exact-rank-two minor opens.

This is a standalone boundary certificate.  It does not edit the canonical
notes or status ledger.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations
import json
from math import factorial, gcd
from pathlib import Path
import shutil
import subprocess
import tempfile

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "two_pair_sic_bidegree44_two_row_boundaries.json"
)
WEIGHTS = (-4, -3, -2, -1, 1, 2, 3, 4)
MAXIMUM_BOUNDARY_CUTOFF = 7


@dataclass(frozen=True)
class IdealResult:
    unit: bool
    result_sha256: str
    result_bytes: int


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in compositions(total - first, parts - 1):
            yield (first, *tail)


def position(weight: int) -> tuple[int, int]:
    """Return the unique separated-row position having torus weight."""
    if weight < 0:
        return (0, -weight)
    return (4, 4 - weight)


def canonical_support(support: frozenset[int]) -> tuple[int, ...]:
    reflected = frozenset(-weight for weight in support)
    return min(tuple(sorted(support)), tuple(sorted(reflected)))


def support_orbits(*, proper: bool) -> list[tuple[int, ...]]:
    representatives: set[tuple[int, ...]] = set()
    maximum_size = len(WEIGHTS) - int(proper)
    for size in range(maximum_size + 1):
        for subset in combinations(WEIGHTS, size):
            representatives.add(canonical_support(frozenset(subset)))
    return sorted(representatives, key=lambda value: (len(value), value))


def primitive_moment(
    support: tuple[int, ...],
    order: int,
) -> tuple[tuple[sp.Symbol, ...], sp.Poly | int]:
    """Return the primitive normalized moment on one support torus."""
    variable_count = max(0, len(support) - 2)
    variables = sp.symbols(f"x0:{variable_count}")
    terms: dict[tuple[int, ...], int] = {}
    for counts in compositions(order, len(support)):
        if sum(
            weight * count
            for weight, count in zip(support, counts, strict=True)
        ):
            continue
        positive_count = sum(
            count
            for weight, count in zip(support, counts, strict=True)
            if weight > 0
        )
        denominator = 1
        for count in counts:
            denominator *= factorial(count)
        scalar = (
            factorial(order)
            // denominator
            * factorial(4 * positive_count)
            * factorial(4 * (order - positive_count))
        )
        exponent = counts[2:]
        terms[exponent] = terms.get(exponent, 0) + scalar

    content = 0
    for scalar in terms.values():
        content = gcd(content, abs(scalar))
    if content:
        terms = {
            exponent: scalar // content
            for exponent, scalar in terms.items()
        }
    if not variables:
        return variables, terms.get((), 0)
    return variables, sp.Poly.from_dict(terms, variables, domain=sp.ZZ)


def normalized_coefficients(
    support: tuple[int, ...],
) -> tuple[tuple[sp.Symbol, ...], dict[int, sp.Expr]]:
    variables = sp.symbols(f"x0:{max(0, len(support) - 2)}")
    coefficients = (sp.Integer(1), sp.Integer(1), *variables)
    return variables, dict(zip(support, coefficients, strict=True))


def row_minors(support: tuple[int, ...]) -> tuple[sp.Expr, ...]:
    """Return the two-row minors when the row supports agree."""
    variables, coefficient = normalized_coefficients(support)
    del variables
    top_columns = sorted(-weight for weight in support if weight < 0)
    bottom_columns = sorted(4 - weight for weight in support if weight > 0)
    if top_columns != bottom_columns:
        return ()
    answer = []
    for left_index, left in enumerate(top_columns):
        for right in top_columns[left_index + 1 :]:
            answer.append(
                sp.expand(
                    coefficient[-left] * coefficient[4 - right]
                    - coefficient[-right] * coefficient[4 - left]
                )
            )
    return tuple(answer)


def matrix_stratum(support: tuple[int, ...]) -> str:
    if not support:
        return "rank_zero"
    negative_columns = {-weight for weight in support if weight < 0}
    positive_columns = {4 - weight for weight in support if weight > 0}
    if not negative_columns or not positive_columns:
        return "rank_one_one_sided"
    if negative_columns == positive_columns:
        if len(negative_columns) == 1:
            return "rank_one_mixed"
        return "rank_one_and_exact_rank_two"
    return "exact_rank_two"


def run_msolve(
    support: tuple[int, ...],
    through: int,
    *,
    invert: tuple[sp.Expr, ...] = (),
    equations: tuple[sp.Expr, ...] = (),
) -> IdealResult:
    """Decide one localized characteristic-zero ideal exactly."""
    msolve = shutil.which("msolve")
    if msolve is None:
        raise RuntimeError("msolve is required")
    variables, coefficient = normalized_coefficients(support)
    del coefficient
    moments = [primitive_moment(support, order)[1] for order in range(1, through + 1)]

    if not variables:
        assert not invert and not equations
        is_unit = any(int(value) != 0 for value in moments)
        encoded = "unit" if is_unit else "nonunit"
        return IdealResult(
            unit=is_unit,
            result_sha256=sha256(encoded.encode()).hexdigest(),
            result_bytes=len(encoded),
        )

    localizer = sp.prod(variables) * sp.prod(invert)
    h = sp.Symbol("h")
    system = [
        sp.expand(h * localizer - 1),
        *equations,
        *(
            value.as_expr()
            for value in moments
            if isinstance(value, sp.Poly) and not value.is_zero
        ),
    ]
    with tempfile.TemporaryDirectory(prefix="sic44-two-row-boundary-") as directory:
        input_path = Path(directory) / "system.ms"
        output_path = Path(directory) / "result.ms"
        input_path.write_text(
            ",".join(str(value) for value in (h, *variables))
            + "\n0\n"
            + ",\n".join(
                str(value).replace("**", "^") for value in system
            )
            + "\n",
            encoding="utf-8",
        )
        completed = subprocess.run(
            [
                msolve,
                "-f",
                str(input_path),
                "-o",
                str(output_path),
                "-t",
                "2",
                "-l",
                "2",
                "-v",
                "0",
            ],
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        if completed.returncode:
            raise RuntimeError(completed.stdout + completed.stderr)
        result = output_path.read_text(encoding="utf-8").strip()
    if not result:
        raise RuntimeError("msolve returned an empty result")
    compact = " ".join(result.split())
    return IdealResult(
        unit=result in {"[-1]", "[-1]:"},
        result_sha256=sha256(compact.encode()).hexdigest(),
        result_bytes=len(result.encode()),
    )


def sharp_cutoff(
    support: tuple[int, ...],
    *,
    invert: tuple[sp.Expr, ...] = (),
    equations: tuple[sp.Expr, ...] = (),
) -> tuple[int, IdealResult, IdealResult]:
    previous = run_msolve(
        support,
        1,
        invert=invert,
        equations=equations,
    )
    assert not previous.unit
    for order in range(2, MAXIMUM_BOUNDARY_CUTOFF + 1):
        current = run_msolve(
            support,
            order,
            invert=invert,
            equations=equations,
        )
        if current.unit:
            return order, previous, current
        previous = current
    raise AssertionError(f"no cutoff through order seven on {support}")


def compact(result: IdealResult) -> dict[str, object]:
    return {
        "unit": result.unit,
        "result_sha256": result.result_sha256,
        "result_bytes": result.result_bytes,
    }


def require_stored_unit(result: object, *, context: str) -> None:
    if not isinstance(result, dict) or result.get("unit") is not True:
        raise AssertionError(f"stored {context} is not a unit certificate")


def require_stored_nonunit(result: object, *, context: str) -> None:
    if not isinstance(result, dict) or result.get("unit") is not False:
        raise AssertionError(f"stored {context} is not a nonunit certificate")


def validate_existing_artifact(path: Path) -> None:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("format") != "two-pair-sic-bidegree44-two-row-boundaries-v1":
        raise AssertionError("unexpected separated-row boundary artifact format")

    representatives = support_orbits(proper=True)
    records = payload.get("support_orbits")
    if not isinstance(records, list):
        raise AssertionError("boundary artifact has no support_orbits list")
    stored_supports = [tuple(record.get("weights", ())) for record in records]
    if stored_supports != representatives:
        raise AssertionError(
            "stored boundary supports are not the exact ordered canonical census"
        )
    if len(stored_supports) != len(set(stored_supports)):
        raise AssertionError("stored boundary census contains duplicate supports")
    if (
        payload.get("proper_support_orbit_count") != len(representatives)
        or len(representatives) != 135
    ):
        raise AssertionError("unexpected stored boundary orbit count")

    stratum_counts: Counter[str] = Counter()
    mixed_cutoffs: Counter[int] = Counter()
    stratified_supports: set[tuple[int, ...]] = set()
    for record, support in zip(records, representatives, strict=True):
        if record.get("support_size") != len(support):
            raise AssertionError(f"wrong support size for {support}")
        if record.get("positions") != [list(position(weight)) for weight in support]:
            raise AssertionError(f"wrong coefficient positions for {support}")
        support_stratum = matrix_stratum(support)
        if record.get("matrix_stratum") != support_stratum:
            raise AssertionError(f"wrong matrix stratum for {support}")
        stratum_counts[support_stratum] += 1

        if support_stratum in {"rank_zero", "rank_one_one_sided"}:
            if "sharp_cutoff" in record:
                raise AssertionError(f"spurious moment cutoff for {support}")
        else:
            cutoff = record.get("sharp_cutoff")
            if not isinstance(cutoff, int) or not 2 <= cutoff <= MAXIMUM_BOUNDARY_CUTOFF:
                raise AssertionError(f"invalid sharp cutoff for {support}")
            require_stored_nonunit(
                record.get("through_previous"),
                context=f"preceding moment ideal on {support}",
            )
            require_stored_unit(
                record.get("through_cutoff"),
                context=f"terminal moment ideal on {support}",
            )
            mixed_cutoffs[cutoff] += 1

        rank_data = record.get("rank_stratification")
        if support_stratum == "rank_one_and_exact_rank_two":
            stratified_supports.add(support)
            expected_minors = tuple(str(value) for value in row_minors(support))
            if not isinstance(rank_data, dict):
                raise AssertionError(f"missing rank stratification for {support}")
            if tuple(rank_data.get("minors", ())) != expected_minors:
                raise AssertionError(f"stored minors do not cover {support}")
            charts = rank_data.get("exact_rank_two_cover")
            if not isinstance(charts, list):
                raise AssertionError(f"missing exact-rank-two cover for {support}")
            if tuple(chart.get("minor") for chart in charts) != expected_minors:
                raise AssertionError(f"minor-open cover is incomplete for {support}")
            for index, chart in enumerate(charts):
                require_stored_nonunit(
                    chart.get("through_previous"),
                    context=f"preceding minor-open ideal {index} on {support}",
                )
                require_stored_unit(
                    chart.get("through_cutoff"),
                    context=f"terminal minor-open ideal {index} on {support}",
                )
            closed = rank_data.get("rank_one_closed_locus")
            if not isinstance(closed, dict):
                raise AssertionError(f"missing rank-one closed locus for {support}")
            require_stored_nonunit(
                closed.get("through_previous"),
                context=f"preceding rank-one ideal on {support}",
            )
            require_stored_unit(
                closed.get("through_cutoff"),
                context=f"terminal rank-one ideal on {support}",
            )
        elif rank_data is not None:
            raise AssertionError(f"spurious rank stratification for {support}")

    expected_strata = Counter(
        {
            "rank_zero": 1,
            "rank_one_one_sided": 15,
            "rank_one_mixed": 2,
            "rank_one_and_exact_rank_two": 3,
            "exact_rank_two": 114,
        }
    )
    expected_cutoffs = Counter({2: 56, 3: 17, 4: 15, 5: 8, 6: 18, 7: 5})
    if stratum_counts != expected_strata:
        raise AssertionError(f"unexpected reconstructed strata: {stratum_counts}")
    if Counter(payload.get("stratum_counts", {})) != expected_strata:
        raise AssertionError("stored stratum totals do not match the exact records")
    if mixed_cutoffs != expected_cutoffs:
        raise AssertionError(f"unexpected reconstructed cutoffs: {mixed_cutoffs}")
    stored_cutoffs = Counter(
        {
            int(order): count
            for order, count in payload.get(
                "mixed_sharp_cutoff_distribution", {}
            ).items()
        }
    )
    if stored_cutoffs != expected_cutoffs:
        raise AssertionError("stored cutoff totals do not match the exact records")
    if len(stratified_supports) != 3:
        raise AssertionError("not every mixed-rank support was stratified")

    print("PASS exact canonical keys cover all 135 proper support orbits once")
    print("PASS stored strata, cutoffs, and all three rank covers are fail-closed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help=(
            "validate exact stored support keys, strata, sharp cutoffs, and "
            "rank covers without rerunning msolve"
        ),
    )
    arguments = parser.parse_args()
    if arguments.audit_existing_only:
        validate_existing_artifact(arguments.output)
        return

    representatives = support_orbits(proper=True)
    assert len(representatives) == 135
    assert representatives[0] == ()

    records: list[dict[str, object]] = []
    mixed_cutoffs: Counter[int] = Counter()
    stratum_counts: Counter[str] = Counter()
    for support in representatives:
        stratum = matrix_stratum(support)
        stratum_counts[stratum] += 1
        record: dict[str, object] = {
            "weights": list(support),
            "positions": [list(position(weight)) for weight in support],
            "support_size": len(support),
            "matrix_stratum": stratum,
        }
        if stratum == "rank_zero":
            record["moment_status"] = "zero polynomial; all moments vanish"
        elif stratum == "rank_one_one_sided":
            for order in range(1, MAXIMUM_BOUNDARY_CUTOFF + 1):
                assert primitive_moment(support, order)[1] == 0
            record["moment_status"] = (
                "one-sided weight support; every positive pure moment "
                "vanishes identically"
            )
        else:
            cutoff, previous, terminal = sharp_cutoff(support)
            mixed_cutoffs[cutoff] += 1
            record.update(
                {
                    "sharp_cutoff": cutoff,
                    "through_previous": compact(previous),
                    "through_cutoff": compact(terminal),
                }
            )

        if stratum == "rank_one_and_exact_rank_two":
            minors = row_minors(support)
            assert minors
            rank_two_charts = []
            rank_two_cutoff = 0
            for minor in minors:
                cutoff, previous, terminal = sharp_cutoff(
                    support,
                    invert=(minor,),
                )
                rank_two_cutoff = max(rank_two_cutoff, cutoff)
                rank_two_charts.append(
                    {
                        "minor": str(minor),
                        "sharp_cutoff": cutoff,
                        "through_previous": compact(previous),
                        "through_cutoff": compact(terminal),
                    }
                )
            rank_one_cutoff, rank_one_previous, rank_one_terminal = sharp_cutoff(
                support,
                equations=minors,
            )
            record["rank_stratification"] = {
                "minors": [str(value) for value in minors],
                "exact_rank_two_cover": rank_two_charts,
                "exact_rank_two_sharp_cutoff": rank_two_cutoff,
                "rank_one_closed_locus": {
                    "sharp_cutoff": rank_one_cutoff,
                    "through_previous": compact(rank_one_previous),
                    "through_cutoff": compact(rank_one_terminal),
                },
            }
        records.append(record)

    assert stratum_counts == Counter(
        {
            "rank_zero": 1,
            "rank_one_one_sided": 15,
            "rank_one_mixed": 2,
            "rank_one_and_exact_rank_two": 3,
            "exact_rank_two": 114,
        }
    )
    assert mixed_cutoffs == Counter({2: 56, 3: 17, 4: 15, 5: 8, 6: 18, 7: 5})

    stratified = {
        tuple(record["weights"]): record["rank_stratification"]
        for record in records
        if "rank_stratification" in record
    }
    assert {
        support: (
            data["exact_rank_two_sharp_cutoff"],
            data["rank_one_closed_locus"]["sharp_cutoff"],
        )
        for support, data in stratified.items()
    } == {
        (-3, -2, 1, 2): (2, 2),
        (-3, -1, 1, 3): (6, 2),
        (-3, -2, -1, 1, 2, 3): (6, 4),
    }

    payload = {
        "format": "two-pair-sic-bidegree44-two-row-boundaries-v1",
        "field": "QQ (geometric support normalization over its algebraic closure)",
        "ambient_support": [list(position(weight)) for weight in WEIGHTS],
        "weight_model": {
            "weight": "w=i-j",
            "weights": list(WEIGHTS),
            "transpose": "w -> -w",
            "simultaneous_reversal": "w -> -w",
            "torus_normalization": (
                "on supports of size at least two, normalize the first two "
                "coefficients to one and localize every residual coefficient"
            ),
        },
        "proper_support_orbit_count": len(representatives),
        "stratum_counts": dict(sorted(stratum_counts.items())),
        "mixed_sharp_cutoff_distribution": {
            str(order): mixed_cutoffs[order]
            for order in sorted(mixed_cutoffs)
        },
        "maximum_boundary_cutoff": MAXIMUM_BOUNDARY_CUTOFF,
        "dense_chart_comparison": (
            "the omitted full eight-coordinate support is the separately "
            "certified chart whose sharp cutoff is eight"
        ),
        "support_orbits": records,
        "conclusion": (
            "every mixed proper coordinate boundary is removed by a sharp "
            "moment cutoff at most seven; the only all-moment coordinate "
            "boundaries are rank-zero or rank-one one-sided supports"
        ),
        "scope": (
            "complete coordinate-boundary classification for the declared "
            "two-separated-row support, including rank-one loci inside the "
            "three equal-row-support tori; no other rank-two factor chart"
        ),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )

    print("PASS classified 135 proper support orbits under weight negation")
    print("PASS 15 nonzero one-sided orbits have every pure moment zero")
    print("PASS all 119 mixed boundary orbits have sharp cutoff at most seven")
    print("PASS split three coefficient tori into rank-one and rank-two strata")
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()
