#!/usr/bin/env python3
"""Certify quotient escape and exact conductors for control-pair fibres.

Candidates are consumed in increasing projective height of the resulting
``t``.  The generic seventeen and the two defining split-bisection points
must first receive an exact finite-quotient independence certificate.  Only
then does PARI globally minimize the fibre and replay its complete local
Tate/conductor data.  Timeouts and failed escapes remain fail-closed.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import prod
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
from time import perf_counter


ROOT = Path(__file__).resolve().parents[2]
ELLIPTIC = ROOT / "elliptic-curves"
sys.path[:0] = [str(ELLIPTIC), str(ELLIPTIC / "cas")]

from analyze_icarm_7fff_zip_sequence import kodaira_symbol  # noqa: E402
from elliptic_candidate_record import (  # noqa: E402
    build_finite_quotient_certificate,
    is_on_weierstrass_curve,
    verify_finite_quotient_certificate,
)


Q = Fraction
POINT_SEARCH = ROOT / "artifacts/generated-results/elkies-2026-control-pair-base-point-search.json"
MODEL = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_model.json"
SECTIONS = ROOT / "elkies-k3/data/fibrations/elkies_2026_published_r17_sections.json"
BISECTIONS = ROOT / "artifacts/generated-results/elkies-2026-equation-bisections-full.json"
FINITE_HELPER = ROOT / "elliptic-curves/cas/elliptic_candidate_record.py"
OUTPUT = ROOT / "artifacts/generated-results/elkies-2026-control-pair-conductor-gate.json"
ICARM394_LOG_CONDUCTOR = "166.2520985277272016652232895273070674463"
ICARM245_LOG_CONDUCTOR = "150.668907152237"
ICARM394_CONDUCTOR = int(
    "1593562111507190066539814084004447718921281851572777685020200143306222910"
)
ICARM245_CONDUCTOR = int(
    "272066437942638823321634957004224153562929337633497250319389959310"
)
LOCAL_PATTERN = re.compile(
    r"^LOCAL\|([0-9]+)\|(-?[0-9]+)\|(-?[0-9]+)\|(-?[0-9]+)\|(-?[0-9]+)$",
    re.MULTILINE,
)


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def rational(value: object) -> Q:
    return value if isinstance(value, Q) else Q(str(value))


def rational_text(value: Q) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def projective_height(value: Q) -> int:
    return max(abs(value.numerator), value.denominator)


def evaluate(coefficients: list[object] | tuple[object, ...], value: Q) -> Q:
    answer = Q(0)
    for coefficient in reversed(coefficients):
        answer = answer * value + rational(coefficient)
    return answer


def check_short_point(coefficient_a: Q, coefficient_b: Q, point: tuple[Q, Q]) -> None:
    x_value, y_value = point
    if y_value**2 != x_value**3 + coefficient_a * x_value + coefficient_b:
        raise ArithmeticError("a materialized bisection point missed its specialized fibre")


def reconstruct_sections(
    section_document: dict[str, object], t_value: Q, coefficient_a: Q, coefficient_b: Q
) -> tuple[tuple[Q, Q], ...]:
    points: list[tuple[Q, Q]] = []
    for expected, record in enumerate(section_document["sections"]):
        if int(record["basis_index"]) != expected:
            raise ArithmeticError("the published section order changed")
        x_value = evaluate(record["x_coefficients_low_to_high"], t_value)
        if expected == 0:
            y_value = evaluate(record["y_coefficients_low_to_high"], t_value)
        else:
            chord = record["chord"]
            reference = points[int(chord["reference_basis_index"])]
            slope = evaluate(chord["slope_coefficients_low_to_high"], t_value)
            y_value = reference[1] + slope * (x_value - reference[0])
        point = (x_value, y_value)
        check_short_point(coefficient_a, coefficient_b, point)
        points.append(point)
    if len(points) != 17:
        raise ArithmeticError("expected seventeen published sections")
    return tuple(points)


def lifted_point(
    record: dict[str, object], t_value: Q, square_root: Q, coefficient_a: Q, coefficient_b: Q
) -> tuple[Q, Q]:
    lifted = record["lifted_section"]
    x_value = evaluate(lifted["x0_coefficients"], t_value) + evaluate(
        lifted["x1_coefficients"], t_value
    ) * square_root
    y_value = evaluate(lifted["y0_coefficients"], t_value) + evaluate(
        lifted["y1_coefficients"], t_value
    ) * square_root
    point = (x_value, y_value)
    check_short_point(coefficient_a, coefficient_b, point)
    return point


def independence_certificate(
    model: tuple[Q, Q, Q, Q, Q], points: tuple[tuple[Q, Q], ...], prime_bound: int
) -> dict[str, object]:
    if any(not is_on_weierstrass_curve(model, point) for point in points):
        raise ArithmeticError("an exact rank-19 point missed its fibre")
    attempts = []
    success = None
    for relation_prime in (2, 3):
        certificate = build_finite_quotient_certificate(
            model,
            points,
            relation_prime=relation_prime,
            prime_bound=prime_bound,
        )
        attempts.append(
            {
                "relation_prime": relation_prime,
                "combined_rank_over_relation_field": int(
                    certificate["combined_rank_over_relation_field"]
                ),
            }
        )
        if certificate["certified_independent"]:
            verify_finite_quotient_certificate(model, points, certificate)
            success = certificate
            break
    return {
        "certified_rank_lower_bound": len(points) if success is not None else 0,
        "independent_point_indices_zero_based": list(range(len(points))) if success else [],
        "attempts": attempts,
        "successful_certificate": success,
    }


def gp_rational(value: Q) -> str:
    return f"({value.numerator}/{value.denominator})"


def valuation(value: int, prime: int) -> int:
    answer = 0
    while value % prime == 0:
        value //= prime
        answer += 1
    return answer


def exact_tate(
    model: tuple[Q, Q, Q, Q, Q], timeout_seconds: float, stack_bytes: int
) -> dict[str, object]:
    executable = shutil.which("gp")
    if executable is None:
        raise FileNotFoundError("PARI/GP executable 'gp' was not found")
    vector = ",".join(gp_rational(value) for value in model)
    program = f'''\
default(realprecision,80);
E=ellinit([{vector}]);v=0;M=ellminimalmodel(E,&v);G=ellglobalred(M);
print("VERSION|",version());
print("MODEL|",M.a1,"|",M.a2,"|",M.a3,"|",M.a4,"|",M.a6);
print("CHANGE|",v[1],"|",v[2],"|",v[3],"|",v[4]);
print("DISC|",M.disc);
print("CONDUCTOR|",G[1]);
print("LOGN|",log(G[1]));
print("TAMAGAWA|",G[3]);
print("ROOT|",ellrootno(M));
for(i=1,matsize(G[4])[1],p=G[4][i,1];L=elllocalred(M,p);print("LOCAL|",p,"|",L[1],"|",L[2],"|",L[4],"|",ellrootno(M,p)));
quit
'''
    completed = subprocess.run(
        [executable, "-q", "-s", str(stack_bytes), "-f"],
        input=program,
        text=True,
        capture_output=True,
        timeout=timeout_seconds,
    )
    combined = completed.stdout + completed.stderr
    if completed.returncode != 0 or "***" in combined:
        raise RuntimeError(" ".join(combined.split())[:1000])
    scalar: dict[str, list[str]] = {}
    for line in completed.stdout.splitlines():
        fields = line.strip().split("|")
        if len(fields) > 1 and fields[0] != "LOCAL":
            scalar[fields[0]] = fields[1:]
    required = {"VERSION", "MODEL", "CHANGE", "DISC", "CONDUCTOR", "LOGN", "TAMAGAWA", "ROOT"}
    if not required <= scalar.keys():
        raise ArithmeticError(f"PARI omitted exact Tate fields: {sorted(required - scalar.keys())}")

    minimal_discriminant = int(scalar["DISC"][0])
    conductor = int(scalar["CONDUCTOR"][0])
    local_rows = []
    for match in LOCAL_PATTERN.finditer(completed.stdout):
        prime = int(match.group(1))
        local_rows.append(
            {
                "prime": str(prime),
                "minimal_discriminant_valuation": valuation(abs(minimal_discriminant), prime),
                "conductor_exponent": int(match.group(2)),
                "kodaira_code": int(match.group(3)),
                "kodaira_symbol": kodaira_symbol(int(match.group(3))),
                "tamagawa_number": int(match.group(4)),
                "local_root_number": int(match.group(5)),
            }
        )
    reconstructed_conductor = prod(
        int(row["prime"]) ** int(row["conductor_exponent"]) for row in local_rows
    )
    reconstructed_tamagawa = prod(int(row["tamagawa_number"]) for row in local_rows)
    if reconstructed_conductor != conductor:
        raise ArithmeticError("local conductor exponents did not reconstruct the exact conductor")
    if reconstructed_tamagawa != int(scalar["TAMAGAWA"][0]):
        raise ArithmeticError("local Tamagawa numbers did not reconstruct their product")
    return {
        "status": "COMPLETE_EXACT_TATE_REPLAY",
        "pari_version": scalar["VERSION"][0].strip("[]").replace(", ", "."),
        "global_minimal_model": scalar["MODEL"],
        "source_to_minimal_change": scalar["CHANGE"],
        "minimal_discriminant": str(minimal_discriminant),
        "conductor": str(conductor),
        "log_conductor": scalar["LOGN"][0],
        "root_number": int(scalar["ROOT"][0]),
        "tamagawa_product": reconstructed_tamagawa,
        "local_reductions": local_rows,
        "conductor_reconstructed_from_local_exponents": str(reconstructed_conductor),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--point-search", type=Path, default=POINT_SEARCH)
    parser.add_argument("--model", type=Path, default=MODEL)
    parser.add_argument("--sections", type=Path, default=SECTIONS)
    parser.add_argument("--bisections", type=Path, default=BISECTIONS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--limit", type=int, default=0, help="zero means all candidates")
    parser.add_argument("--certificate-prime-bound", type=int, default=300)
    parser.add_argument("--tate-timeout", type=float, default=300.0)
    parser.add_argument("--gp-stack-bytes", type=int, default=512_000_000)
    parser.add_argument("--checkpoint-every", type=int, default=10)
    args = parser.parse_args()
    if args.limit < 0:
        parser.error("limit must be nonnegative")
    started = perf_counter()
    sys.set_int_max_str_digits(0)

    point_document = json.loads(args.point_search.read_text())
    model_document = json.loads(args.model.read_text())
    section_document = json.loads(args.sections.read_text())
    bisection_document = json.loads(args.bisections.read_text())
    if point_document.get("schema") != "elkies-k3.elkies-2026-control-pair-base-point-search.v1":
        raise ValueError("unexpected control-pair point-search schema")
    if len(bisection_document.get("bisections", ())) != 39120:
        raise ArithmeticError("the complete bisection batch changed")
    by_mask = {
        int(record["lattice_orbit_mask"]): record for record in bisection_document["bisections"]
    }
    a_coefficients = tuple(rational(value) for value in model_document["A_coefficients_low_to_high"])
    b_coefficients = tuple(rational(value) for value in model_document["B_coefficients_low_to_high"])
    candidates = sorted(
        point_document["fibres"],
        key=lambda row: (projective_height(rational(row["t"])), rational(row["t"])),
    )
    if args.limit:
        candidates = candidates[: args.limit]

    rows = []
    complete = timeouts = failed_escape = 0
    for ordinal, candidate in enumerate(candidates, start=1):
        t_value = rational(candidate["t"])
        coefficient_a = evaluate(a_coefficients, t_value)
        coefficient_b = evaluate(b_coefficients, t_value)
        generic_points = reconstruct_sections(
            section_document, t_value, coefficient_a, coefficient_b
        )
        source = candidate["sources"][0]
        masks = tuple(int(value) for value in source["orbit_masks"])
        roots = (rational(source["u"]), rational(source["v"]))
        split_points = tuple(
            lifted_point(by_mask[mask], t_value, root, coefficient_a, coefficient_b)
            for mask, root in zip(masks, roots)
        )
        points = generic_points + split_points
        short_model = (Q(0), Q(0), Q(0), coefficient_a, coefficient_b)
        independence = independence_certificate(
            short_model, points, args.certificate_prime_bound
        )
        row = {
            "ordinal_by_resulting_t_height": ordinal,
            "t": rational_text(t_value),
            "t_projective_height": str(projective_height(t_value)),
            "t_projective_bits": int(candidate["t_projective_bits"]),
            "pair_key": source["pair_key"],
            "orbit_masks": list(masks),
            "quotient_escape_certificate": independence,
        }
        if int(independence["certified_rank_lower_bound"]) != 19:
            failed_escape += 1
            row.update(
                status="REJECT_NO_CERTIFIED_RANK19_ESCAPE",
                wider_point_search_authorized=False,
            )
            rows.append(row)
            continue
        try:
            tate = exact_tate(short_model, args.tate_timeout, args.gp_stack_bytes)
        except subprocess.TimeoutExpired:
            timeouts += 1
            row.update(
                status="INCOMPLETE_TATE_TIMEOUT",
                tate_timeout_seconds=args.tate_timeout,
                wider_point_search_authorized=False,
            )
        except Exception as error:  # fail closed and preserve the candidate
            row.update(
                status="INCOMPLETE_TATE_ERROR",
                error_type=type(error).__name__,
                error=str(error),
                wider_point_search_authorized=False,
            )
        else:
            complete += 1
            conductor = int(tate["conductor"])
            beats_394 = conductor < ICARM394_CONDUCTOR
            beats_245 = conductor < ICARM245_CONDUCTOR
            row.update(
                status="PASS_EXACT_RANK19_AND_TATE_GATE",
                tate=tate,
                targets={
                    "beats_icarm394_log_conductor": beats_394,
                    "beats_icarm245_log_conductor": beats_245,
                },
                wider_point_search_authorized=beats_394,
            )
        rows.append(row)
        if ordinal % args.checkpoint_every == 0:
            print(
                "ELKIES2026CONTROLPAIRCONDUCTOR|"
                f"progress={ordinal}/{len(candidates)}|complete={complete}|timeouts={timeouts}|"
                f"authorized={sum(item['wider_point_search_authorized'] for item in rows)}",
                flush=True,
            )

    result = {
        "schema": "elkies-k3.elkies-2026-control-pair-conductor-gate.v1",
        "status": (
            "PASS_COMPLETE_DECLARED_HEIGHT_PREFIX"
            if len(rows) == len(candidates)
            and all(not row["status"].startswith("INCOMPLETE") for row in rows)
            else "PARTIAL_FAIL_CLOSED"
        ),
        "inputs": {
            display_path(Path(__file__).resolve()): digest(Path(__file__).resolve()),
            display_path(args.point_search): digest(args.point_search),
            display_path(args.model): digest(args.model),
            display_path(args.sections): digest(args.sections),
            display_path(args.bisections): digest(args.bisections),
            display_path(FINITE_HELPER): digest(FINITE_HELPER),
        },
        "bounds": {
            "height_order": "max(abs(numerator(t)), denominator(t)), then exact t",
            "candidate_limit_zero_means_all": args.limit,
            "certificate_prime_bound": args.certificate_prime_bound,
            "tate_timeout_seconds_per_candidate": args.tate_timeout,
            "gp_stack_bytes": args.gp_stack_bytes,
        },
        "targets": {
            "icarm394_log_conductor": ICARM394_LOG_CONDUCTOR,
            "icarm394_exact_conductor": str(ICARM394_CONDUCTOR),
            "icarm245_log_conductor": ICARM245_LOG_CONDUCTOR,
            "icarm245_exact_conductor": str(ICARM245_CONDUCTOR),
        },
        "summary": {
            "candidate_count": len(candidates),
            "certified_rank19_count": sum(
                int(row["quotient_escape_certificate"]["certified_rank_lower_bound"]) == 19
                for row in rows
            ),
            "complete_exact_tate_count": complete,
            "tate_timeout_count": timeouts,
            "failed_quotient_escape_count": failed_escape,
            "wider_point_search_authorized_count": sum(
                bool(row["wider_point_search_authorized"]) for row in rows
            ),
        },
        "candidates": rows,
        "runtime_seconds": perf_counter() - started,
        "reproducing_command": shlex.join(sys.argv),
        "proof_boundary": (
            "Each rank lower bound has an exact finite-quotient certificate. Each COMPLETE Tate record "
            "is globally minimal and reconstructs its exact conductor from all local exponents. Timeouts "
            "and errors authorize no point search. A certified rank-19 fibre still needs two additional "
            "independent rational points before it can challenge the rank-at-least-21 record."
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        "ELKIES2026CONTROLPAIRCONDUCTOR|"
        f"candidates={len(rows)}|rank19={result['summary']['certified_rank19_count']}|"
        f"complete_tate={complete}|authorized={result['summary']['wider_point_search_authorized_count']}|"
        f"status={result['status']}|output={display_path(args.output)}"
    )


if __name__ == "__main__":
    main()
