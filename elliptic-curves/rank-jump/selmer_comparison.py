#!/usr/bin/env python3
"""Certify a complete global boundary using local duality and reciprocity."""
import argparse
from pathlib import Path
import retrospective as r
import local_collision as lc
import affine_selmer as af

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "SELMER_COMPARISON_PROTOCOL.json"
OUTPUT = r.OUT / "rank_jump_selmer_comparison_v1.json"


def quotient_rows(rows, subspace):
    """Represent a quotient injectively by its canonical remainder."""
    pivots = r.basis(subspace)
    return [r.reduce(v, pivots) for v in rows]


def analyze(case):
    # The dimension is determined by E[2] and the local field, so it is
    # identical for E0 and Eu under the specified 2-torsion identification.
    boundary = [0] * 21
    old_conditions, new_conditions = [], []
    local = []
    offset = 0
    total_half_dimension = 0
    for entry in case["local"]:
        signatures = list(map(r.pack, entry["class_signature_rows"]))
        width = len(entry["class_signature_rows"][0])
        L0 = lc.canonical(signatures[:20])
        Lu = lc.canonical(map(r.pack, entry["point_signature_rows"]))
        assert len(L0) == len(Lu) == entry["point_dimension"]
        C = lc.intersection(L0, Lu, n=width)
        D = lc.canonical(L0 + Lu)
        half = len(L0) - len(C)
        assert len(D) - len(C) == 2 * half
        assert not any(quotient_rows(signatures, D))
        reduced = quotient_rows(signatures, C)
        for i, value in enumerate(reduced):
            boundary[i] |= value << offset
        old_conditions.extend(r.transpose(
            [[(v >> j) & 1 for j in range(width)]
             for v in quotient_rows(signatures, L0)]))
        new_conditions.extend(r.transpose(
            [[(v >> j) & 1 for j in range(width)]
             for v in quotient_rows(signatures, Lu)]))
        if entry.get("new_relative_to_anchor"):
            # Anchor L0 is unramified here, while eta is ramified.
            assert all(not any(row[::2])
                       for row in entry["class_signature_rows"][:20])
            assert any(entry["class_signature_rows"][20][::2])
        local.append({
            "place": entry["place"], "ambient_signature_width": width,
            "boundary_offset": offset, "L0_basis": L0, "Lu_basis": Lu,
            "intersection_basis": C, "sum_basis": D,
            "local_half_dimension": half,
            "all_21_classes_in_local_sum": True,
            "class_boundary_rows": reduced,
        })
        offset += width
        total_half_dimension += half
    boundary_rank = r.rank(boundary)
    assert boundary_rank == total_half_dimension
    old_kernel = lc.orthogonal(old_conditions, 21)
    new_kernel = lc.orthogonal(new_conditions, 21)
    strict_kernel = lc.orthogonal(old_conditions + new_conditions, 21)
    assert old_kernel == [1 << i for i in range(20)]
    assert strict_kernel == lc.canonical(case["old_inherited_basis"])
    assert len(strict_kernel) == 21 - boundary_rank
    old_image_rank = len(old_kernel) - len(strict_kernel)
    new_image_rank = len(new_kernel) - len(strict_kernel)
    assert old_image_rank == boundary_rank - 1
    assert new_image_rank in (0, 1)
    return {
        "u": case["u"], "local": local,
        "global_class_boundary_rows": boundary,
        "local_boundary_dimension": 2 * total_half_dimension,
        "complete_global_boundary_dimension": boundary_rank,
        "boundary_completeness_status": "PROVED_BY_RECIPROCITY_AND_LOCAL_DUALITY",
        "strict_kernel_anchor_masks": strict_kernel,
        "old_global_space_kernel": old_kernel,
        "new_global_space_kernel": new_kernel,
        "old_selmer_boundary_dimension": old_image_rank,
        "new_selmer_boundary_dimension": new_image_rank,
        "full_selmer_dimension_minus_anchor_dimension": new_image_rank - old_image_rank,
        "known_selmer_intersection_dimension": len(new_kernel),
        "full_selmer_new_ramification_dimension_exact": new_image_rank,
        "unknown_excess_dimension": "dim Sel_2(E0/Q) - 20",
        "unramified_selmer_equals_full_common_selmer": True,
    }


def build(check=False):
    raw = r.read(af.INPUT)
    affine = r.read(af.OUTPUT)
    cases = []
    for u in r.read(PROTOCOL)["parameters"]:
        case = analyze(next(c for c in raw["cases"] if c["u"] == u))
        previous = next(c for c in affine["cases"] if c["u"] == u)
        assert case["new_selmer_boundary_dimension"] == int(
            previous["affine_solution"]["consistent"])
        assert case["known_selmer_intersection_dimension"] == len(
            case["strict_kernel_anchor_masks"]) + case["new_selmer_boundary_dimension"]
        cases.append(case)
    paths = (PROTOCOL, Path(__file__), af.INPUT, af.OUTPUT,
             HERE / "retrospective.py", HERE / "local_collision.py")
    out = {
        "schema": "rank-jump.selmer-comparison.v1",
        "bindings": {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in paths},
        "cases": cases,
        "scope": "Exact relative dimensions of full Selmer groups. Their absolute excess above the retained subspaces is the same unknown anchor excess. No absolute Selmer dimension or unconditional numerical curve-rank upper bound is claimed.",
    }
    if check:
        assert r.read(OUTPUT) == out
        print("PASS complete global boundary and relative full Selmer dimensions")
    else:
        r.write_new(OUTPUT, out)
    for case in cases:
        print(case["u"], "boundary", case["complete_global_boundary_dimension"],
              "Selmer difference", case["full_selmer_dimension_minus_anchor_dimension"],
              "ramification dimension", case["full_selmer_new_ramification_dimension_exact"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("build", "check"))
    build(parser.parse_args().mode == "check")
