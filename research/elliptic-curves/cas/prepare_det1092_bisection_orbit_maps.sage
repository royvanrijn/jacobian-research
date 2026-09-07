#!/usr/bin/env sage-python
"""Rank all generic norm-ten rational-bisection orbits on one held-out fibre."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from decimal import Decimal, localcontext
from importlib.machinery import SourceFileLoader
from pathlib import Path

from sage.all import ZZ, matrix


CAS = Path(__file__).resolve().parent
ROOT = CAS.parents[1]
sys.path.insert(0, str(CAS))
import certify_compact_r17_candidates as cert
import det1092_bisection_orbit_holdout as campaign
from research_runtime.store import checkpoint, digest


mapper = SourceFileLoader("det1092_bisection_factor_free_mapper", str(CAS / "factor_free_pari_mapping.sage")).load_module()
geometry = SourceFileLoader("det1092_bisection_height_geometry", str(CAS / "prospective_half_lattice_v3.sage")).load_module()


def orbit_rows(parent_gram):
    lattice = cert.read(campaign.LATTICE)
    if lattice["status"] != "PASS_COMPLETE_DEGREE2_TRANSLATION_QUOTIENT" or campaign.sha(campaign.ORBITS) != lattice["orbits_tsv_sha256"]:
        raise ArithmeticError("complete generic bisection table changed")
    rows = []
    with campaign.ORBITS.open() as stream:
        for item in csv.DictReader(stream, delimiter="\t"):
            if item["category"] != "rational" or int(item["minimum_norm"]) != campaign.GENERIC_NORM:
                continue
            word = tuple(map(int, item["parent_MW17_w"].split()))
            if len(word) != campaign.INITIAL_RANK:
                raise ArithmeticError("generic bisection orbit lost MW17 coordinates")
            norm = sum(word[i] * int(parent_gram[i][j]) * word[j] for i in range(len(word)) for j in range(len(word)))
            if norm != campaign.GENERIC_NORM:
                raise ArithmeticError("generic bisection norm changed")
            rows.append({"orbit_mask": int(item["orbit_mask"]), "generic_norm": norm, "representative": list(word),
                         "generic_l1": sum(abs(value) for value in word), "generic_linf": max(abs(value) for value in word)})
    if len(rows) != lattice["rational_bisections"]["translation_orbits"]:
        raise ArithmeticError("generic rational norm-ten orbit enumeration is incomplete")
    return sorted(rows, key=lambda row: row["orbit_mask"])


def specialized_norms(words, gram):
    """Exact all-orbit quadratic values; use a bounded int64 fast path if safe."""
    maximum_word = max(abs(value) for row in words for value in row["representative"])
    maximum_gram = max(abs(int(value)) for row in gram for value in row)
    safe = campaign.INITIAL_RANK**2 * maximum_word**2 * maximum_gram < 2**62
    if safe:
        import numpy as np
        vectors = np.asarray([row["representative"] for row in words], dtype=np.int64)
        form = np.asarray([[int(value) for value in row] for row in gram], dtype=np.int64)
        values = np.einsum("ij,jk,ik->i", vectors, form, vectors, optimize=True)
        return [int(value) for value in values], "numpy-int64-exact-bound-checked"
    values = []
    for row in words:
        word = row["representative"]
        values.append(sum(word[i] * int(gram[i][j]) * word[j] for i in range(len(word)) for j in range(len(word))))
    return values, "python-integer-exact"


def main(index: int) -> None:
    campaign.install_execution_guard()
    row, directory = campaign.configure(index)
    protocol, seed = campaign.campaign(), cert.read(directory / "seed.json")
    output = directory / "maps.json"
    if output.exists():
        raise FileExistsError("preserve frozen full-orbit map roster")
    parent = cert.read(campaign.LOCAL / "parent-sections.json")
    generic_gram = parent["generic_height_gram"]
    all_orbits = orbit_rows(generic_gram)
    model = tuple(map(cert.F, seed["curve"]))
    points = tuple(tuple(map(cert.F, point)) for point in seed["points"])
    gram, asymmetry = geometry.canonical_height_gram(model, points)
    rounded = geometry.rounded_gram(gram, 1_000_000)
    g = matrix(ZZ, rounded)
    if not g.is_symmetric() or not g.is_positive_definite():
        raise ArithmeticError("specialized rounded height form is not positive definite")
    values, evaluation_engine = specialized_norms(all_orbits, g.rows())
    ranked = []
    for orbit, norm in zip(all_orbits, values):
        if norm <= 0:
            raise ArithmeticError("positive specialized height norm expected")
        ranked.append({**orbit, "specialized_metric_norm": norm})
    ranked.sort(key=lambda item: (-item["specialized_metric_norm"], item["generic_l1"], item["generic_linf"], item["orbit_mask"]))
    selected = ranked[:campaign.CHARTS]
    if len(selected) != campaign.CHARTS or len({item["orbit_mask"] for item in selected}) != campaign.CHARTS:
        raise ArithmeticError("fixed equal-budget bisection roster is malformed")
    ranking_payload = json_bytes([{key: item[key] for key in ("orbit_mask", "specialized_metric_norm", "generic_l1", "generic_linf")} for item in ranked])
    data = {
        "schema": "elliptic-curves.det1092-bisection-orbit-maps.v1", "status": "RUNNING_MAPS",
        "protocol_sha256": campaign.sha(campaign.LOCAL / "protocol.json"), "seed_sha256": campaign.sha(directory / "seed.json"),
        "all_norm10_orbit_count": len(all_orbits), "ranking_sha256": hashlib.sha256(ranking_payload).hexdigest(),
        "ranking_evaluation_engine": evaluation_engine, "metric_gram": [[str(value) for value in line] for line in gram],
        "maximum_gram_asymmetry": str(asymmetry), "rounded_gram": [list(map(int, line)) for line in g.rows()],
        "ranking_rule": protocol["centre_policy"]["ranking"], "centres": selected, "rows": [],
    }
    checkpoint(output, data)
    mapper.pari.allocatemem(256000000, silent=True)
    for position, centre in enumerate(selected):
        mapped = mapper.mapping(model, points, centre)
        if mapped["centre"] != centre:
            raise ArithmeticError("factor-free chart lost its full-orbit centre")
        data["rows"].append({"index": position, **mapped})
        checkpoint(output, data)
        print("BISECTION MAP", row["id"], position + 1, "/", campaign.CHARTS, flush=True)
    data["status"] = "COMPLETE_DECLARED_BISECTION_MAPS"
    checkpoint(output, data)
    print("FROZEN ALL-NORM10 TOP49 BISECTION MAPS", row["id"], len(all_orbits), flush=True)


def json_bytes(value):
    import json
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--index", required=True, type=int)
    main(parser.parse_args().index)
