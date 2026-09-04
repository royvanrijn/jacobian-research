#!/usr/bin/env python3
"""Replay the affine-support/Newton bridge obstruction and Kummer gate.

This checker has two roles.

1. Dense triangular automorphisms show that Newton vertex count,
   geometric degree, and reduced nonproperness data cannot upper-bound
   affine-normalized support.
2. A monomial-Jacobian Laurent block ``[P,Q]=x^r`` descends through
   ``u=x^(r+1)/(r+1)`` only when both components lie in ``k[u,y]``.
   The live F2 terminal block is checked exactly and fails this character
   gate, so the constant-Jacobian support-six theorem does not apply to it.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = ROOT / "artifacts/generated-results/jc2_affine_support_newton_bridge.json"
ARTIFACT_SHA256 = (
    "d77c1418a332ae0ecb74e8fc7a3f4ae28e2da0ea137b85e814541c826ca5087b"
)


def polynomial_support(
    polynomial: sp.Expr,
    variables: tuple[sp.Symbol, sp.Symbol],
) -> list[list[int]]:
    return [
        [int(x_degree), int(y_degree)]
        for (x_degree, y_degree), _ in sp.Poly(
            sp.expand(polynomial), *variables
        ).terms()
    ]


def triangular_regression() -> dict[str, object]:
    """Check explicit dense-shear witnesses through degree twelve."""

    t = sp.symbols("t")
    rows: list[dict[str, object]] = []
    digest = hashlib.sha256()
    for degree in range(4, 13):
        polynomial = sum(
            t**power / sp.factorial(power)
            for power in range(2, degree + 1)
        )
        derivatives = {
            order: sp.Poly(sp.diff(polynomial, t, order), t, domain=sp.QQ)
            for order in range(2, degree)
        }
        resultant_count = 0
        for left_order, left in derivatives.items():
            for right_order in range(left_order + 1, degree):
                right = derivatives[right_order]
                gcd = sp.gcd(left, right)
                if gcd.degree() != 0:
                    raise AssertionError(
                        f"degree {degree} derivatives {left_order},{right_order} "
                        "acquired a common root"
                    )
                resultant = sp.resultant(left.as_expr(), right.as_expr(), t)
                if resultant == 0:
                    raise AssertionError("a derivative resultant vanished")
                numerator, denominator = sp.fraction(resultant)
                record = (
                    f"{degree}:{left_order}:{right_order}:"
                    f"{numerator}:{denominator}\n"
                )
                digest.update(record.encode())
                resultant_count += 1
        rows.append(
            {
                "degree": degree,
                "polynomial": f"sum(t^k/k!, k=2..{degree})",
                "pairwise_nonzero_derivative_resultants": resultant_count,
                "affine_support_lower_bound": degree - 2,
            }
        )
    return {
        "theorem": (
            "a Zariski-open set of degree-d triangular automorphisms has "
            "sigma_aff at least d-2 while its first-coordinate Newton polygon "
            "has three vertices and its geometric degree is one"
        ),
        "proof_mechanism": (
            "pairwise root-disjoint derivatives p^(2),...,p^(d-1) force at "
            "least d-2 nonzero Taylor degrees at every translation center"
        ),
        "explicit_regression_rows": rows,
        "resultant_digest_sha256": digest.hexdigest(),
        "claim_boundary": (
            "the generic theorem is proved dimensionally in the note; these "
            "explicit rows are deterministic witnesses, not its proof"
        ),
    }


def kummer_gate_regression() -> dict[str, object]:
    """Check a descending example and the non-descending F2 terminal block."""

    x, y, u = sp.symbols("x y u")

    # A support-four Keller automorphism in (u,y), pulled back by u=x^5/5.
    first = u + y**2
    second = y + first**2
    if sp.expand(
        sp.diff(first, u) * sp.diff(second, y)
        - sp.diff(first, y) * sp.diff(second, u)
    ) != 1:
        raise AssertionError("the descended representative is not Keller")
    pulled_first = sp.expand(first.subs(u, x**5 / 5))
    pulled_second = sp.expand(second.subs(u, x**5 / 5))
    pulled_bracket = sp.expand(
        sp.diff(pulled_first, x) * sp.diff(pulled_second, y)
        - sp.diff(pulled_first, y) * sp.diff(pulled_second, x)
    )
    if pulled_bracket != x**4:
        raise AssertionError("the Kummer pullback chain rule failed")

    # The exact normalized terminal F2 block from F2_75_125_DERIVATION.md.
    s = x**17 * y**5
    f2_first = sp.expand(x**4 * y * (1 + s))
    f2_second = sp.expand(-x * (1 + 3 * s + sp.Rational(9, 5) * s**2))
    f2_bracket = sp.expand(
        sp.diff(f2_first, x) * sp.diff(f2_second, y)
        - sp.diff(f2_first, y) * sp.diff(f2_second, x)
    )
    if f2_bracket != x**4:
        raise AssertionError("the F2 terminal bracket changed")
    f2_first_support = polynomial_support(f2_first, (x, y))
    f2_second_support = polynomial_support(f2_second, (x, y))
    first_characters = sorted({exponent[0] % 5 for exponent in f2_first_support})
    second_characters = sorted({exponent[0] % 5 for exponent in f2_second_support})
    if first_characters != [1, 4] or second_characters != [0, 1, 3]:
        raise AssertionError("the F2 Kummer-character profile changed")
    first_sectors: dict[int, sp.Expr] = {}
    second_sectors: dict[int, sp.Expr] = {}
    for (x_degree, y_degree), coefficient in sp.Poly(
        f2_first, x, y
    ).terms():
        character = x_degree % 5
        first_sectors[character] = (
            first_sectors.get(character, 0)
            + coefficient * x**x_degree * y**y_degree
        )
    for (x_degree, y_degree), coefficient in sp.Poly(
        f2_second, x, y
    ).terms():
        character = x_degree % 5
        second_sectors[character] = (
            second_sectors.get(character, 0)
            + coefficient * x**x_degree * y**y_degree
        )
    bracket_rows: list[dict[str, object]] = []
    sector_sums: dict[int, sp.Expr] = {}
    for p_character, p_sector in sorted(first_sectors.items()):
        for q_character, q_sector in sorted(second_sectors.items()):
            bracket = sp.expand(
                sp.diff(p_sector, x) * sp.diff(q_sector, y)
                - sp.diff(p_sector, y) * sp.diff(q_sector, x)
            )
            bracket_character = (p_character + q_character - 1) % 5
            sector_sums[bracket_character] = sp.expand(
                sector_sums.get(bracket_character, 0) + bracket
            )
            bracket_rows.append(
                {
                    "P_character": p_character,
                    "Q_character": q_character,
                    "bracket_character": bracket_character,
                    "bracket": str(sp.factor(bracket)),
                }
            )
    expected_sector_sums = {0: 0, 1: 0, 3: 0, 4: x**4}
    if sector_sums != expected_sector_sums:
        raise AssertionError("the F2 character-sector bracket cancellation changed")

    return {
        "gate": (
            "if [P,Q]_(x,y)=c*x^r and P,Q lie in "
            "k[x^(r+1)/(r+1),y], then their descended bracket is c"
        ),
        "descending_representative": {
            "r": 4,
            "in_u_y": {
                "P": "u+y^2",
                "Q": "y+(u+y^2)^2",
                "bracket": "1",
            },
            "pulled_back_support": {
                "P": polynomial_support(pulled_first, (x, y)),
                "Q": polynomial_support(pulled_second, (x, y)),
            },
            "pulled_back_bracket": "x^4",
        },
        "f2_75_125_terminal_block": {
            "P_support": f2_first_support,
            "Q_support": f2_second_support,
            "bracket": "x^4",
            "kummer_modulus": 5,
            "P_x_characters": first_characters,
            "Q_x_characters": second_characters,
            "character_bracket_rule": "character([P_a,Q_b])=a+b-1 mod 5",
            "character_bracket_rows": bracket_rows,
            "sector_sums": {
                str(character): str(value)
                for character, value in sorted(sector_sums.items())
            },
            "descends_to_k[u,y]": False,
            "conclusion": (
                "the support-six constant-Jacobian theorem cannot be applied "
                "to the terminal block without the missing lower bands"
            ),
        },
    }


def build_payload() -> dict[str, object]:
    return {
        "schema_version": 1,
        "claim_boundary": {
            "proved": [
                "Newton vertex count, geometric degree, and reduced "
                "nonproperness data do not upper-bound sigma_aff",
                "the Kummer descent chain-rule gate",
                "the exact F2 terminal character obstruction",
            ],
            "not_proved": [
                "a support bound for minimal counterexamples",
                "an exclusion of the (75,125) F2 family",
                "a new degree or geometric-degree frontier",
            ],
        },
        "triangular_affine_support_obstruction": triangular_regression(),
        "kummer_character_gate": kummer_gate_regression(),
        "software": {
            "python": "standard library",
            "sympy": sp.__version__,
        },
    }


def audit_existing_only(artifact: Path) -> None:
    """Audit the committed bridge ledger without symbolic recomputation."""

    actual_hash = hashlib.sha256(artifact.read_bytes()).hexdigest()
    if actual_hash != ARTIFACT_SHA256:
        raise AssertionError("the pinned affine-support bridge artifact bytes changed")
    payload = json.loads(artifact.read_text())
    assert payload["schema_version"] == 1
    assert payload["claim_boundary"] == {
        "proved": [
            "Newton vertex count, geometric degree, and reduced "
            "nonproperness data do not upper-bound sigma_aff",
            "the Kummer descent chain-rule gate",
            "the exact F2 terminal character obstruction",
        ],
        "not_proved": [
            "a support bound for minimal counterexamples",
            "an exclusion of the (75,125) F2 family",
            "a new degree or geometric-degree frontier",
        ],
    }

    triangular = payload["triangular_affine_support_obstruction"]
    rows = triangular["explicit_regression_rows"]
    assert [row["degree"] for row in rows] == list(range(4, 13))
    for row in rows:
        degree = row["degree"]
        assert row["affine_support_lower_bound"] == degree - 2
        assert row["pairwise_nonzero_derivative_resultants"] == (
            (degree - 2) * (degree - 3) // 2
        )
    assert triangular["resultant_digest_sha256"] == (
        "0b3f061bc9d344eae722e8aee30751d61c854bebe7176aa8c1b65e323a98e499"
    )
    assert "not its proof" in triangular["claim_boundary"]

    gate = payload["kummer_character_gate"]
    descending = gate["descending_representative"]
    assert descending["r"] == 4
    assert descending["in_u_y"]["bracket"] == "1"
    assert descending["pulled_back_bracket"] == "x^4"
    assert all(
        x_degree % 5 == 0
        for component in descending["pulled_back_support"].values()
        for x_degree, _y_degree in component
    )

    terminal = gate["f2_75_125_terminal_block"]
    assert terminal["kummer_modulus"] == 5
    assert terminal["P_x_characters"] == sorted(
        {x_degree % 5 for x_degree, _y_degree in terminal["P_support"]}
    ) == [1, 4]
    assert terminal["Q_x_characters"] == sorted(
        {x_degree % 5 for x_degree, _y_degree in terminal["Q_support"]}
    ) == [0, 1, 3]
    assert terminal["descends_to_k[u,y]"] is False
    assert terminal["bracket"] == "x^4"
    assert terminal["sector_sums"] == {
        "0": "0",
        "1": "0",
        "3": "0",
        "4": "x**4",
    }
    bracket_rows = terminal["character_bracket_rows"]
    assert len(bracket_rows) == 6
    assert {
        (row["P_character"], row["Q_character"])
        for row in bracket_rows
    } == {(p_character, q_character) for p_character in (1, 4) for q_character in (0, 1, 3)}
    assert "cannot be applied" in terminal["conclusion"]

    print(
        "PASS committed affine-support bridge audit: pinned artifact, regression "
        "boundary, and complete F2 character-pair ledger; no symbolic replay"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    parser.add_argument(
        "--audit-existing-only",
        action="store_true",
        help="validate the committed artifact without symbolic recomputation",
    )
    args = parser.parse_args()
    artifact = args.artifact.resolve()
    if args.audit_existing_only:
        audit_existing_only(artifact)
        return

    payload = build_payload()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        try:
            display_path = artifact.relative_to(ROOT)
        except ValueError:
            display_path = artifact
        print(f"WROTE {display_path}")
    else:
        expected = json.loads(artifact.read_text())
        current_claim = {key: value for key, value in payload.items() if key != "software"}
        pinned_claim = {key: value for key, value in expected.items() if key != "software"}
        if current_claim != pinned_claim:
            raise AssertionError(
                "pinned bridge artifact is stale; inspect before --refresh"
            )
    print(
        "PASS affine-support bridge obstruction:",
        "Newton vertex count and coarse boundary data cannot bound sigma_aff",
    )
    print(
        "PASS Kummer-character gate:",
        "the F2 terminal block has characters P={1,4}, Q={0,1,3} mod 5",
    )


if __name__ == "__main__":
    main()
