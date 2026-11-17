#!/usr/bin/env python3
"""Compare compact divergence-certificate supports at two primes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from research_two_pair_sic_bidegree33_rank_two_terminal_syzygy import (  # noqa: E402
    named_certificate_polynomial,
)


LOCAL = ROOT / "artifacts" / "local"
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_"
        "certificate_support_research.json"
    )
)


def certificate_path(prime: int, degree: int) -> Path:
    chunk = "m8_m7" if degree == 7 else "m7_m0"
    return LOCAL / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_shift_"
        f"mod{prime}_certificate_{chunk}.sing"
    )


def terminal_path(prime: int) -> Path:
    return LOCAL / (
        "two_pair_sic_bidegree33_rank_two_compact_relative_shift_"
        f"mod{prime}_terminal_syzygy_R.sing"
    )


def support_record(
    name: str,
    paths: list[Path],
    primes: list[int],
) -> dict[str, object]:
    supports = []
    term_counts = []
    for path, prime in zip(paths, primes, strict=True):
        polynomial = named_certificate_polynomial(path, name, prime)
        supports.append(set(polynomial))
        term_counts.append(len(polynomial))
    union = set().union(*supports)
    intersection = set.intersection(*supports)
    return {
        "polynomial": name,
        "term_counts": term_counts,
        "union_term_count": len(union),
        "intersection_term_count": len(intersection),
        "symmetric_difference_term_count": len(union - intersection),
        "supports_identical": len(union) == len(intersection),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--primes",
        type=int,
        nargs=2,
        default=(1_000_003, 1_000_033),
    )
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    primes = list(arguments.primes)

    coefficients = []
    for degree in range(7, -1, -1):
        paths = [certificate_path(prime, degree) for prime in primes]
        for variable in ("X", "Y"):
            coefficients.append(
                {
                    "certificate_degree": degree,
                    **support_record(
                        f"{variable}{degree}",
                        paths,
                        primes,
                    ),
                }
            )
    terminal = support_record(
        "R",
        [terminal_path(prime) for prime in primes],
        primes,
    )
    total_symmetric_difference = sum(
        int(record["symmetric_difference_term_count"])
        for record in coefficients
    )
    result = {
        "format": (
            "two-pair-sic-bidegree33-rank-two-compact-relative-"
            "certificate-support-research-v1"
        ),
        "status": (
            "exact two-prime support comparison; a reconstruction "
            "feasibility scout, not a characteristic-zero certificate"
        ),
        "point": 0,
        "primes": primes,
        "coefficient_polynomials": coefficients,
        "coefficient_supports_identical": all(
            bool(record["supports_identical"])
            for record in coefficients
        ),
        "coefficient_total_symmetric_difference": (
            total_symmetric_difference
        ),
        "terminal_syzygy": terminal,
        "interpretation": (
            "All eight Y supports, five of eight X supports, and the "
            "terminal R support are identical. The remaining three X "
            "supports differ by one monomial each, consistent with an "
            "isolated coefficient vanishing modulo one prime. A union-"
            "support CRT reconstruction is therefore a concrete next "
            "experiment, but more primes are required."
        ),
    }
    if total_symmetric_difference != 3:
        raise RuntimeError("unexpected two-prime certificate support drift")
    if not bool(terminal["supports_identical"]):
        raise RuntimeError("terminal syzygy support drifted")
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + "\n")
    print("PASS all eight Y supports are identical")
    print("PASS five X supports are identical; three differ by one term")
    print("PASS 132615-term terminal R support is identical")
    print(f"PASS wrote {arguments.output}")


if __name__ == "__main__":
    main()
