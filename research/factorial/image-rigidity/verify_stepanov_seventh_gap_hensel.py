#!/usr/bin/env python3
"""Replay the seven exact Hensel witnesses for the lower-three seventh gap.

The accompanying written proof supplies the uniform truncation and period.
This program uses only the Python standard library and direct integer products.
It does not generate or assume an all-order coprimality result from a scan.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from math import comb, prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
P = 7
PRECISION = 4
MODULUS = P ** PRECISION
PERIOD = P ** (PRECISION + 1)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def valuation(value: int) -> int:
    require(value != 0, "valuation requested at zero")
    answer = 0
    while value % P == 0:
        value //= P
        answer += 1
    return answer


def normalized_coefficient(order: int, index: int) -> int:
    """Direct exact product, independent of the block-recursion generator."""
    if index > order:
        return 0
    gamma = 3 * order + 1
    numerator = prod(range(gamma, gamma + P * index))
    denominator = P ** index
    require(numerator % denominator == 0, "nonintegral normalized product")
    return comb(order, index) * (numerator // denominator)


def add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return tuple((a + b) % MODULUS for a, b in zip(left, right))


def multiply(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    a, b = left
    c, d = right
    return ((a * c - 3 * b * d) % MODULUS,
            (a * d + b * c + b * d) % MODULUS)


def evaluate(coefficients: list[int], argument: tuple[int, int]) -> tuple[int, int]:
    value = (0, 0)
    for coefficient in reversed(coefficients):
        value = add(multiply(value, argument), (coefficient, 0))
    return value


def norm_mod_seven(value: tuple[int, int]) -> int:
    a, b = value
    return (a * a + a * b + 3 * b * b) % P


def verify(witness: dict, extra_controls: bool = True) -> dict:
    require(witness["prime"] == P, "wrong prime")
    require(witness["precision"] == PRECISION, "wrong precision")
    require(witness["modulus"] == MODULUS, "wrong modulus")
    require(witness["period"] == PERIOD, "wrong period")
    require(witness["extension_polynomial_low_to_high"] == [3, -1, 1],
            "wrong quadratic extension")
    require(witness["coefficient_cutoff"] == 28, "wrong coefficient cutoff")
    require(all((x * x - x + 3) % P for x in range(P)),
            "quadratic is reducible modulo seven")
    records = witness["records"]
    require([row["terminal"] for row in records]
            == [446 + MODULUS * i for i in range(P)], "incomplete residue coverage")
    require(len({row["terminal"] % PERIOD for row in records}) == P,
            "duplicate terminal residues")

    coefficient_count = 0
    tail_count = 0
    period_count = 0
    checked = []
    for row in records:
        terminal = row["terminal"]
        lower = row["lower_coefficients"]
        upper = row["upper_coefficients"]
        for order, coefficients in [(terminal - 1, lower), (terminal, upper)]:
            require(len(coefficients) == 28, "missing coefficient")
            expected = [normalized_coefficient(order, k) % MODULUS
                        for k in range(28)]
            require(coefficients == expected, f"coefficient mismatch at order {order}")
            coefficient_count += 28
            if extra_controls:
                for k in range(28, 35):
                    require(normalized_coefficient(order, k) % MODULUS == 0,
                            f"tail valuation control failed at ({order},{k})")
                    tail_count += 1
                shifted = [normalized_coefficient(order + PERIOD, k) % MODULUS
                           for k in range(28)]
                require(shifted == expected, f"period control failed at order {order}")
                period_count += 28

        lift = tuple(row["lower_root_lift"])
        require(len(lift) == 2 and all(0 <= c < MODULUS for c in lift),
                "invalid root representation")
        require(tuple(c % P for c in lift) == (5, 2), "wrong residue disk")
        require(evaluate(lower, lift) == (0, 0), "lower lift is not a zero")
        derivative = [(k * lower[k]) % MODULUS for k in range(1, len(lower))]
        derivative_value = evaluate(derivative, lift)
        require(norm_mod_seven(derivative_value) != 0, "lower derivative is not a unit")
        upper_value = evaluate(upper, lift)
        require(list(upper_value) == row["upper_value"], "wrong upper value")
        require(all(c % (P ** 3) == 0 for c in upper_value),
                "upper value lacks the stated divisibility")
        divided = tuple((c // (P ** 3)) % P for c in upper_value)
        require(list(divided) == row["upper_value_divided_by_343_mod_7"],
                "wrong divided upper value")
        require(norm_mod_seven(divided) != 0, "upper value vanishes modulo 7^4")
        checked.append({"terminal": terminal, "lower_root_lift": list(lift),
                        "lower_derivative_norm_mod_7": norm_mod_seven(derivative_value),
                        "upper_value_divided_by_343_mod_7": list(divided),
                        "divided_upper_norm_mod_7": norm_mod_seven(divided)})

    binomial_valuations = [valuation(comb(PERIOD, j)) for j in range(1, 28)]
    require(all(v == 5 - valuation(j) and v >= 4
                for j, v in enumerate(binomial_valuations, start=1)),
            "binomial period valuation mismatch")
    return {"residue_classes": len(records), "coefficient_equalities": coefficient_count,
            "tail_controls": tail_count, "shifted_coefficient_controls": period_count,
            "binomial_period_controls": len(binomial_valuations), "records": checked}


def negative_controls(witness: dict) -> list[str]:
    corruptions = []
    altered = copy.deepcopy(witness)
    altered["records"][0]["lower_coefficients"][0] += 1
    corruptions.append(("changed lower coefficient", altered))
    altered = copy.deepcopy(witness)
    altered["records"][1]["lower_root_lift"][0] += 1
    corruptions.append(("changed root lift", altered))
    altered = copy.deepcopy(witness)
    altered["records"][2]["upper_value"][0] += 343
    corruptions.append(("changed upper evaluation", altered))
    altered = copy.deepcopy(witness)
    altered["records"].pop()
    corruptions.append(("missing terminal residue", altered))
    rejected = []
    for name, altered in corruptions:
        try:
            verify(altered, extra_controls=False)
        except ValueError:
            rejected.append(name)
        else:
            raise ValueError(f"corrupt witness accepted: {name}")
    return rejected


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--witness", type=Path,
                        default=HERE / "stepanov_seventh_gap_hensel_witness.json")
    parser.add_argument("--receipt", type=Path,
                        default=HERE / "stepanov_seventh_gap_hensel_verification.json")
    args = parser.parse_args()
    witness = json.loads(args.witness.read_text())
    result = verify(witness)
    rejected = negative_controls(witness)
    proof = HERE / "STEPANOV_SEVENTH_GAP_HENSEL_CERTIFICATE.md"
    receipt = {"status": "PASS", "scope": "seven exact residue certificates plus written uniform truncation and period",
               "base_commit": "19eb9912fb60142f8b620dac24896931d4ff3e46",
               "prime": P, "modulus": MODULUS, "period": PERIOD,
               **result, "negative_controls_rejected": rejected,
               "sha256": {args.witness.name: sha256(args.witness),
                          Path(__file__).name: sha256(Path(__file__)),
                          proof.name: sha256(proof)}}
    args.receipt.write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "residue_classes": result["residue_classes"],
                      "coefficient_equalities": result["coefficient_equalities"],
                      "negative_controls_rejected": len(rejected),
                      "receipt": str(args.receipt)}))


if __name__ == "__main__":
    main()
