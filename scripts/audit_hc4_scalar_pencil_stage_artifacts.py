#!/usr/bin/env python3
"""Audit committed HC4RSD29--40 stage ledgers without symbolic replay.

These artifacts preserve the pure-sextic and septic stages of the scalar
constant-Hessian-pencil programme.  The generating scripts import SymPy and
rewrite their outputs, so retrospective maintenance uses this stdlib-only
auditor instead.  It verifies the immutable bytes and the status/format
mapping; it does not re-run any determinant, elimination, or search.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACTS = {
    "HC4RSD29": (
        "hc4_pure_sextic_collision.json",
        "hc4-pure-sextic-collision-v1",
        "d339fc7f9d023d7d5fca9f62df7a02cc5afa0974e49dd278dfe06712c881846f",
    ),
    "HC4RSD30": (
        "hc4_pure_sextic_lower_flag.json",
        "hc4-pure-sextic-lower-flag-v1",
        "244b026e8d2d918e76d049e3ee12ba30ebf86cbb6aff79210ba13c18115866d9",
    ),
    "HC4RSD31": (
        "hc4_pure_sextic_affine_quartic.json",
        "hc4-pure-sextic-affine-quartic-v1",
        "87926fc0d17094f070cc2f0a5fe10770b02291f8da741402077d0d575ea8a5c4",
    ),
    "HC4RSD32": (
        "hc4_pure_sextic_two_linear_tower.json",
        "hc4-pure-sextic-two-linear-tower-v1",
        "27afb422629fb62c1f0e9213053540ff244d5fbec0a133ee98d954123c3913ea",
    ),
    "HC4RSD33": (
        "hc4_nonpure_septic.json",
        "hc4-nonpure-septic-v1",
        "ae4cdf70e00972c3e0996c2210254c10b2dfaa915193b59fecf97b0c7e59c3c3",
    ),
    "HC4RSD34": (
        "hc4_pure_septic_opening.json",
        "hc4-pure-septic-opening-v1",
        "c559bd3a599a7bfa248b5c0b398c5053341d222f8705896b273e1a56f2729d48",
    ),
    "HC4RSD35": (
        "hc4_pure_septic_degree18.json",
        "hc4-pure-septic-degree18-v1",
        "00f15e915c0ffa4583d04db8e950714f5435a4efcce2c60e1f5380a98ce99acd",
    ),
    "HC4RSD36": (
        "hc4_pure_septic_moving_closure.json",
        "hc4-pure-septic-moving-closure-v1",
        "c2130bfafafe98e16c8e6ad379be24edbde997b2bc1e10435851abb1456ddda4",
    ),
    "HC4RSD37": (
        "hc4_pure_septic_kzero.json",
        "hc4-pure-septic-kzero-v1",
        "6701db165bc075bcd93f5b6f0decba57d13bcc5a777ef5662a37f44a39ff0841",
    ),
    "HC4RSD38": (
        "hc4_pure_septic_kzero_wronskian.json",
        "hc4-pure-septic-kzero-wronskian-v1",
        "42cfeb805a3c5085cade6fbea5e406351476669a6446f99c0a4bd39e411c25c9",
    ),
    "HC4RSD39": (
        "hc4_pure_septic_passive_affine.json",
        "hc4-pure-septic-passive-affine-v1",
        "b363c08fdda331043c67720bd56689f3755f6a91dfc326683e045709c385797b",
    ),
    "HC4RSD40": (
        "hc4_pure_septic_quartic_packets.json",
        "hc4-pure-septic-quartic-packets-v1",
        "2fbbdc76c98577158438e24129a33ba634601805ab2ac8d67dd9becda4c3ba06",
    ),
}


def main() -> None:
    artifact_root = ROOT / "artifacts" / "generated-results"
    for entry_id, (filename, expected_format, expected_hash) in ARTIFACTS.items():
        path = artifact_root / filename
        actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual_hash == expected_hash, (entry_id, actual_hash, expected_hash)
        payload = json.loads(path.read_text(encoding="utf-8"))
        assert payload["format"] == expected_format, entry_id
        assert payload["status"]["id"] == entry_id, entry_id
        assert payload["status"]["scope"], entry_id

    assert list(ARTIFACTS) == [f"HC4RSD{index}" for index in range(29, 41)]
    print(
        "PASS: 12 committed HC4RSD29--40 pure-sextic/septic stage ledgers "
        "match their exact hashes and status mappings; no SymPy import, "
        "symbolic replay, search, or artifact rewrite"
    )
    print(
        "SCOPE: HC4RSD32 closes degree six, HC4RSD40 closes degree seven, "
        "and HC4MR1 later closes the complete auxiliary relative-nilpotent "
        "pencil branch; unrestricted HC4 is not claimed"
    )


if __name__ == "__main__":
    main()
