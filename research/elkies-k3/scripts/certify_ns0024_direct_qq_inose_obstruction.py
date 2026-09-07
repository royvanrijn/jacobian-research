#!/usr/bin/env python3
"""Certify the arithmetic reduction in the direct NS0024 Inose source route.

The mathematical input is the complete Mazur--Kenku classification of cyclic
isogeny degrees over QQ.  This checker does not reprove that external theorem;
it verifies its exact application to the height-950 Inose source recorded by
the foundry certificate and states the descent boundary explicitly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
SOURCE = GENERATED / "elkies-k3-ns0024-new-rootless-source-route-v1.json"
OUTPUT = GENERATED / "elkies-k3-ns0024-direct-qq-inose-obstruction-v1.json"

# The complete degree list in Mazur--Kenku, written without interval notation.
RATIONAL_CYCLIC_ISOGENY_DEGREES = tuple(range(1, 20)) + (
    21,
    25,
    27,
    37,
    43,
    67,
    163,
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_payload():
    source = json.loads(SOURCE.read_text())
    assert source["status"] == (
        "PASS_EXACT_FRAME_AND_SOURCE_ROUTE_WITH_OPEN_EQUATION_TRANSPORT"
    )
    equation_source = source["equation_source"]
    assert equation_source["source_id"] == "NS0024-S001"
    assert equation_source["frame_type"] == "2E8/MW1"
    assert equation_source["mw_height"] == 950
    assert equation_source["inose_equation"]["condition"] == (
        "E1 and E2 are non-isomorphic and joined by a cyclic 475-isogeny"
    )

    isogeny_degree = equation_source["mw_height"] // 2
    assert 2 * isogeny_degree == equation_source["mw_height"]
    assert isogeny_degree == 475 == 5**2 * 19
    assert isogeny_degree not in RATIONAL_CYCLIC_ISOGENY_DEGREES

    return {
        "schema": "elkies-k3.ns0024-direct-qq-inose-obstruction.v1",
        "status": "PASS_DIRECT_QQ_INOSE_SOURCE_OBSTRUCTION",
        "source_id": "NS0024-S001",
        "source_frame": "2E8/MW1",
        "mw_height": 950,
        "required_cyclic_isogeny_degree": isogeny_degree,
        "required_degree_factorization": {"5": 2, "19": 1},
        "mazur_kenku_rational_cyclic_isogeny_degrees": list(
            RATIONAL_CYCLIC_ISOGENY_DEGREES
        ),
        "degree_is_allowed_over_QQ": False,
        "conclusion": (
            "There is no noncuspidal QQ-rational point on X0(475), hence no "
            "pair of elliptic curves over QQ joined by a QQ-rational cyclic "
            "475-isogeny. The direct Utsumi specialization advertised for "
            "NS0024-S001 therefore cannot supply its height-950 rational "
            "Mordell--Weil generator over QQ."
        ),
        "external_theorem": {
            "name": "Mazur--Kenku classification of rational cyclic isogenies",
            "statement_used": (
                "The displayed list is the complete set of degrees of cyclic "
                "isogenies of elliptic curves over QQ."
            ),
            "modern_reference": (
                "Banwait--Najman--Padurariu, Cyclic isogenies of elliptic "
                "curves over fixed quadratic fields, Theorem 1.1 and Table 1.1"
            ),
            "url": "https://arxiv.org/abs/2206.08891",
        },
        "input_hashes": {relative(SOURCE): digest(SOURCE)},
        "proof_boundary": {
            "proved": (
                "Conditional only on the cited complete Mazur--Kenku theorem, "
                "the direct E1,E2,phi over QQ realization of the recorded "
                "height-950 2E8/MW1 Inose source is impossible."
            ),
            "not_proved": (
                "Nonexistence of rational points on Atkin--Lehner quotients, "
                "nonexistence of quadratic Q-curve descents, or nonexistence "
                "of any QQ model of NS0024. Any quotient descent still has to "
                "prove that all nineteen divisor classes descend individually."
            ),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    rendered = json.dumps(build_payload(), indent=2, sort_keys=True) + "\n"
    if args.check:
        assert output.read_text() == rendered
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)
    print("PASS ns0024 direct QQ Inose obstruction")


if __name__ == "__main__":
    main()
