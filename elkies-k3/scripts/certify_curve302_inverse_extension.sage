#!/usr/bin/env sage-python
"""Pack/replay every witness of the complete two-chart 302 extension.

The compact artifact needs no local campaign checkpoints for replay: all trace
words are recovered from their hashed complete tables and each exclusion is
reconstructed from the rational source model modulo its witness prime.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from importlib.machinery import SourceFileLoader
import json
from pathlib import Path

from sage.all import EllipticCurve, GF, PolynomialRing, QQ, ZZ, matrix, vector

ROOT = Path(__file__).resolve().parents[2]
GEN = ROOT / "artifacts/generated-results"
SEARCH = ROOT / "elkies-k3/scripts/search_curve302_inverse_fibrations.sage"
OUTPUT = GEN / "elkies-k3-curve302-inverse-fibration-extension-full-v1.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pack", type=Path)
    parser.add_argument("--artifact", type=Path, default=OUTPUT)
    args = parser.parse_args()
    search = SourceFileLoader("certify_inverse302_search", str(SEARCH)).load_module()
    if args.pack:
        full = json.loads(args.pack.read_text())
        if full["status"] != "BOUNDED_EXACT_NO_SPECIALIZATION" or full["selected_count"] != 127842:
            raise ArithmeticError("campaign is not the declared complete two-table miss")
        sources = []
        for source in full["sources"]:
            records = source["records"]
            count = source["available_classes"]
            if [r["priority_rank"] for r in records] != list(range(1, count+1)):
                raise ArithmeticError("campaign does not contain the full table in order")
            if any(r["status"] != "EXCLUDED_NO_RATIONAL_PARAMETER" for r in records):
                raise ArithmeticError("unresolved record in claimed complete exclusion")
            sources.append({"source_chart": source["source_chart"], "class_count": count,
                            "witness_primes_in_priority_order": [r["witness_prime"] for r in records],
                            "minimum_multiplicity_counts": dict(Counter(r["minimum_unoriented_count"] for r in records)),
                            "witness_prime_counts": dict(Counter(r["witness_prime"] for r in records))})
        artifact = {"schema": "elkies-k3.curve302-inverse-complete-extension.v1",
                    "status": "PASS_REPLAYED_COMPLETE_TWO_NORM8_TABLE_EXCLUSION",
                    "target_curve": 302, "class_count": 127842, "sources": sources,
                    "generic_rank_filter": None, "ranked_specialization_hits": [],
                    "full_campaign_sha256": search.digest(args.pack),
                    "inputs": {**full["inputs"], str(Path(__file__).relative_to(ROOT)): search.digest(Path(__file__))},
                    "proof_boundary": "Complete only for the 63925 minimum-norm-eight classes on 103b2 and the 63917 on 08f72, with no generic-rank filter. No rational target parameter exists in any of these pencils. The classes need not be distinct fibrations across charts. This is not inverse recovery of all NS possibilities and does not exclude a parent of 302. No overlap score zero is assigned to a nonexistent specialization.",
                    "reproducing_command": "sage -python elkies-k3/scripts/certify_curve302_inverse_extension.sage"}
    else:
        artifact = json.loads(args.artifact.read_text())
    for path, expected in artifact["inputs"].items():
        if search.digest(ROOT / path) != expected:
            raise ArithmeticError(f"input hash mismatch: {path}")
    helpers = search.load("certify_inverse302_helpers", search.SCRIPTS / "screen_icarm_curve398_norm8_a1_fibrations.sage")
    chord = search.load("certify_inverse302_chord", search.SCRIPTS / "construct_elkies_2026_bisections.sage")
    public = search.load("certify_inverse302_public", ROOT / "elliptic-curves/cas/icarm_curve302.py")
    curve = EllipticCurve(QQ, [QQ(str(c)) for c in public.GENERAL_WEIERSTRASS_COEFFICIENTS])
    c4, c6 = curve.c_invariants()
    replayed = 0
    for source in artifact["sources"]:
        label = source["source_chart"]
        model = json.loads((GEN / f"elkies-k3-r17-norm12-orbit{label}-direct-fibration-v1.json").read_text())
        table = GEN / f"elkies-k3-r17-norm12-{search.SOURCES[label]}-pencil-priority-v1.tsv"
        rows = list(csv.DictReader(table.open(), delimiter="\t"))
        witnesses = source["witness_primes_in_priority_order"]
        if len(rows) != len(witnesses) or len(rows) != source["class_count"]:
            raise ArithmeticError("incomplete witness vector")
        if [int(row["priority_rank"]) for row in rows] != list(range(1, len(rows)+1)):
            raise ArithmeticError("priority table ordering changed")
        if dict(Counter(int(r["minimal_unoriented_count"]) for r in rows)) != {
            int(k): v for k, v in source["minimum_multiplicity_counts"].items()}:
            raise ArithmeticError("multiplicity census differs")
        gram = matrix(ZZ, model["sections"]["height_gram"])
        contexts = {}
        for position, (row, prime) in enumerate(zip(rows, witnesses), 1):
            if prime not in contexts:
                if not ZZ(prime).is_prime() or curve.discriminant().numerator() % prime == 0:
                    raise ArithmeticError("invalid witness prime")
                contexts[prime] = search.modular_context(model, prime, helpers)
            word = tuple(map(int, row["section_basis_w"].split()))
            v = vector(ZZ, word)
            if len(word) != 17 or v*gram*v != 8:
                raise ArithmeticError("trace norm changed")
            polynomial = search.comparison_for_word(contexts[prime], word, c4, c6, helpers, chord)
            if not search.projective_no_root(polynomial, 24, prime):
                raise ArithmeticError(f"invalid exclusion witness: {label}/{position}")
            replayed += 1
            if position % 8192 == 0:
                print(f"INVERSE302VERIFY|source={label}|replayed={position}/{len(rows)}", flush=True)
    # Regression for the fixed homogeneous degree: a dropped leading
    # coefficient must retain a possible rational point reducing to infinity.
    ring = PolynomialRing(GF(1009), "t")
    t = ring.gen()
    if search.projective_no_root(t**23+1, 24, 1009):
        raise ArithmeticError("projective infinity regression failed")
    if search.projective_no_root(t**24-t, 24, 1009):
        raise ArithmeticError("affine root regression failed")
    if replayed != artifact["class_count"]:
        raise ArithmeticError("global coverage mismatch")
    if args.pack:
        search.write(args.artifact, artifact)
    print(f"INVERSE302VERIFY|reconstructed_and_verified={replayed}|projective_regressions=PASS|artifact={args.artifact}", flush=True)


if __name__ == "__main__":
    main()
