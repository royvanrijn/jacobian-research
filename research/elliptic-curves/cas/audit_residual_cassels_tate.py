#!/usr/bin/env python3
"""Audit a residual 2-Selmer Cassels--Tate pairing certificate.

This is deliberately the final handoff after a global K(S,2) computation has
produced a *certified* residual Selmer basis.  It does not manufacture a
pairing: Magma/Sage/Magma-cover code must supply its entries and identify the
algorithm.  The audit verifies the alternating GF(2) linear algebra, derives
the radical, and refuses any cover search recorded before the pairing.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


PROTOCOL = "BNFFREECT"


class PairingError(ValueError):
    pass


def f2_rank(rows: list[list[int]]) -> int:
    pivots: dict[int, int] = {}
    for row in rows:
        packed = sum(int(bit) << index for index, bit in enumerate(row))
        while packed:
            pivot = packed.bit_length() - 1
            if pivot in pivots:
                packed ^= pivots[pivot]
            else:
                pivots[pivot] = packed
                break
    return len(pivots)


def audit(record: dict) -> dict:
    basis = record.get("residual_basis")
    matrix = record.get("cassels_tate_matrix")
    if not isinstance(basis, list) or not isinstance(matrix, list):
        raise PairingError("residual_basis and cassels_tate_matrix are required lists")
    dimension = len(basis)
    if len(matrix) != dimension or any(not isinstance(row, list) or len(row) != dimension for row in matrix):
        raise PairingError("Cassels--Tate matrix must be square on the residual basis")
    if any(entry not in (0, 1) for row in matrix for entry in row):
        raise PairingError("Cassels--Tate entries must be GF(2) bits")
    for i in range(dimension):
        if matrix[i][i] != 0:
            raise PairingError("Cassels--Tate pairing must be alternating")
        for j in range(dimension):
            if matrix[i][j] != matrix[j][i]:
                raise PairingError("Cassels--Tate matrix must be symmetric in characteristic two")

    if record.get("cover_searches") and not record.get("pairing_computed_before_cover_search"):
        raise PairingError("cover search is recorded before the Cassels--Tate pairing")
    pairing_rank = f2_rank(matrix)
    if pairing_rank % 2:
        raise PairingError("an alternating pairing has even rank")
    radical_dimension = dimension - pairing_rank
    known_rank = record["known_mw_rank"]
    if type(known_rank) is not int or known_rank < 0:
        raise PairingError("known_mw_rank must be nonnegative")
    # Alternating linear algebra does not establish that this is the actual
    # arithmetic pairing. Require separate, unconditional input certificates.
    certified = all(record.get(flag) is True for flag in (
        "residual_selmer_basis_certified", "pairing_entries_certified",
        "known_mw_rank_certified", "unconditional",
    )) and all(isinstance(record.get(field), str) and record[field].strip()
               for field in ("pairing_algorithm", "pairing_evidence"))
    classification = (
        "CERTIFIED_RESIDUAL_PAIRING"
        if certified
        else "UNCERTIFIED_BASIS_OR_PAIRING_AUDIT"
    )
    output = {
        "protocol": "BNFFREECT-v1",
        "classification": classification,
        "residual_basis_dimension": dimension,
        "pairing_rank": pairing_rank,
        "radical_dimension": radical_dimension,
        "rank_upper_after_pairing": known_rank + radical_dimension if certified else None,
        "conditional_rank_upper_after_pairing": known_rank + radical_dimension,
        "basis_labels": [str(item.get("label", index)) if isinstance(item, dict) else str(index) for index, item in enumerate(basis)],
        "pairing_algorithm": record.get("pairing_algorithm"),
        "pairing_evidence": record.get("pairing_evidence"),
        "interpretation": (
            "The rank bound requires the supplied residual Selmer basis, actual "
            "pairing entries, and known point rank to be certified unconditionally "
            "on the same curve. This audit checks linear algebra and input "
            "attestations; it does not compute or verify the arithmetic pairing. "
            "Non-radical directions are accounted for by the Cassels--Tate "
            "obstruction, not by a cover-point search."
        ),
    }
    if certified and radical_dimension == 0:
        output["classification"] = "CERTIFIED_EXACT_KNOWN_RANK_AFTER_PAIRING"
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    record = json.loads(args.input.read_text())
    if not isinstance(record, dict):
        raise PairingError("top-level JSON value must be an object")
    output = audit(record)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"{PROTOCOL}|stage=complete|classification={output['classification']}"
        f"|radical_dimension={output['radical_dimension']}"
        f"|rank_upper_after_pairing={output['rank_upper_after_pairing']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
