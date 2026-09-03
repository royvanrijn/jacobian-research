#!/usr/bin/env sage-python
"""Certify prospective mask-generated H3/NS0024 cores and the Q80 near miss."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
BRIDGES = GENERATED / "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
THETA = GENERATED / "elkies-k3-integral-rank-transfer-theta-convolution-v1.json"
BASE_SCRIPT = ROOT / "elkies-k3/scripts/generate_integral_rank_transfer_masked_core_neighbors.sage"
SEARCH_SCRIPT = ROOT / "elkies-k3/scripts/search_integral_rank_transfer_masked_core_controls.sage"
CORE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_core_generation.sage"
REVERSE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_reverse_theta_masks.sage"
OUTPUT = GENERATED / "elkies-k3-integral-rank-transfer-masked-core-controls-v1.json"

# Each witness is in the basis of the preceding quadratic form.  H3 was
# extracted from the root-descent beam; NS0024 from the bridge-gated capped
# support-diversity beam; Q80 is the best exact near miss from the analogous
# eight-generation support-diversity control.
PATHS = {
    "H3": (
        (7, (6, 6, 2, 5, 1, 6, 2, 3, 0, 3, 4, 4, 4, 1, 1)),
        (7, (5, 1, 1, 1, 2, 0, 0, 0, 4, 5, 4, 3, 3, 2, 5)),
        (13, (5, 0, 4, 5, 6, 3, 8, 8, 5, 6, 1, 10, 5, 7, 12)),
        (19, (4, 6, 1, 14, 2, 18, 5, 1, 17, 10, 1, 9, 17, 6, 16)),
        (19, (6, 10, 17, 4, 7, 10, 14, 5, 0, 14, 11, 16, 7, 8, 7)),
        (17, (7, 7, 7, 13, 14, 2, 3, 0, 5, 1, 0, 2, 13, 8, 13)),
        (7, (6, 2, 4, 2, 5, 4, 1, 5, 5, 4, 0, 4, 0, 6, 0)),
        (11, (10, 2, 0, 1, 8, 0, 4, 8, 2, 7, 6, 1, 5, 3, 6)),
    ),
    "NS0024": (
        (17, (14, 6, 3, 1, 7, 9, 3, 15, 2, 0, 6, 1, 12, 12, 14)),
        (13, (0, 4, 5, 1, 4, 6, 9, 8, 7, 6, 1, 5, 6, 0, 12)),
        (7, (0, 3, 4, 3, 3, 0, 6, 5, 6, 3, 5, 4, 5, 3, 3)),
    ),
    "Q80": (
        (19, (9, 8, 4, 12, 9, 6, 14, 10, 9, 7, 10, 18, 18, 9, 13)),
        (17, (6, 10, 9, 0, 6, 8, 14, 5, 4, 5, 1, 15, 3, 9, 13)),
        (17, (16, 3, 13, 10, 9, 8, 0, 6, 4, 16, 0, 3, 6, 4, 7)),
        (13, (7, 9, 10, 3, 12, 11, 9, 12, 3, 0, 5, 7, 7, 2, 12)),
        (11, (5, 8, 1, 10, 8, 4, 3, 1, 1, 3, 10, 6, 9, 6, 8)),
        (13, (11, 5, 5, 8, 10, 4, 8, 11, 2, 9, 10, 3, 10, 5, 0)),
        (19, (13, 6, 5, 0, 17, 9, 16, 1, 5, 3, 6, 15, 7, 18, 14)),
        (17, (4, 15, 6, 5, 16, 11, 4, 3, 3, 12, 9, 14, 5, 2, 3)),
    ),
}


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def replay(seed, path, base):
    current = base["quadratic_form"](seed)
    states = []
    for step, (prime, raw_witness) in enumerate(path, start=1):
        assert seed.det() % prime
        witness = vector(ZZ, raw_witness)
        assert any(value % prime for value in witness)
        assert current(witness) % prime == 0
        current = current.find_p_neighbor_from_vec(prime, witness)
        gram = current.Hessian_matrix()
        assert gram.det() == seed.det()
        assert all(value % 2 == 0 for value in gram.diagonal())
        states.append(
            {
                "step": step,
                "prime": prime,
                "witness": list(raw_witness),
                "signed_root_count": int(pari(gram).qfminim(2)[0]),
            }
        )
    return base["lll_reduce"](current.Hessian_matrix()), states


def completion(core_gram, prepared, accepted, base, core):
    bridge = next(
        row
        for row in prepared["viable_bridges"]
        if row["bridge_class_index"] == accepted["bridge_class_index"]
    )
    generator = base["primary_generator"](core_gram, prepared["order"])
    multiplier = accepted["isotropic_multipliers"][0]
    glue = vector(QQ, list(multiplier * generator) + list(bridge["generator"]))
    child = core["glued_frame"](core_gram, bridge["gram"], glue)
    assert int(pari(child).qfminim(2)[0]) == 0
    assert core["minimum_norm"](child) == 4
    assert child.det() == prepared["target_frame"].det()
    discriminant_matches = (
        core["discriminant_form_key"](child)
        == core["discriminant_form_key"](prepared["target_frame"])
    )
    assert discriminant_matches
    isometric = bool(
        pari(base["lll_reduce"](child)).qfisom(
            pari(base["lll_reduce"](prepared["target_frame"]))
        )
    )
    return {
        "bridge_class_index": accepted["bridge_class_index"],
        "glue_multiplier": multiplier,
        "rank": child.nrows(),
        "determinant": int(child.det()),
        "minimum": 4,
        "signed_root_count": 0,
        "discriminant_form_matches_target": discriminant_matches,
        "isometric_to_declared_target_frame": isometric,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    bridge_artifact = json.loads(BRIDGES.read_text())
    theta_artifact = json.loads(THETA.read_text())
    base = runpy.run_path(str(BASE_SCRIPT))
    search = runpy.run_path(str(SEARCH_SCRIPT))
    core = runpy.run_path(str(CORE_SCRIPT))
    reverse = runpy.run_path(str(REVERSE_SCRIPT))
    rows = []
    for corridor in ("H3", "NS0024", "Q80"):
        prepared = search["prepare_corridor"](
            corridor,
            bridge_artifact,
            theta_artifact,
            base,
            core,
            reverse,
        )
        search["configure_order"](base, prepared["order"])
        generated_core, states = replay(prepared["seed"], PATHS[corridor], base)
        assert int(pari(generated_core).qfminim(2)[0]) == 0
        assert core["minimum_norm"](generated_core) == 4
        masks, lazy_queries, _ = base["mask_profile"](
            generated_core, prepared["viable_bridges"], reverse
        )
        accepted = [row for row in masks if row["zero_mask_accepts"]]
        if corridor in ("H3", "NS0024"):
            assert len(accepted) == 1
            child = completion(generated_core, prepared, accepted[0], base, core)
        else:
            assert not accepted
            assert [row["occupied_forbidden_cells"] for row in masks] == [2]
            child = None
        historical_isometric = bool(
            pari(generated_core).qfisom(
                pari(base["lll_reduce"](prepared["historical_core"]))
            )
        )
        if corridor in ("H3", "NS0024"):
            assert not historical_isometric
        rows.append(
            {
                "corridor": corridor,
                "cyclic_bridge_order": prepared["order"],
                "canonical_seed_signed_roots": int(
                    pari(prepared["seed"]).qfminim(2)[0]
                ),
                "path": states,
                "generated_core": {
                    "rank": generated_core.nrows(),
                    "determinant": int(generated_core.det()),
                    "minimum": 4,
                    "signed_root_count": 0,
                    "automorphism_group_order": int(
                        pari(generated_core).qfauto()[0]
                    ),
                    "historical_core_isometric": historical_isometric,
                    "gram": [
                        list(map(int, row)) for row in generated_core.rows()
                    ],
                    "viable_bridge_mask_results": masks,
                    "lazy_core_cells_queried": lazy_queries,
                },
                "completion": child,
                "classification": (
                    "rootless_completion"
                    if child is not None
                    else "rootless_two-cell_near_miss"
                ),
            }
        )

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-masked-core-controls.v1",
        "status": "PASS_TWO_PROSPECTIVE_COMPLETIONS_AND_ONE_NEAR_MISS",
        "inputs": {
            relative(BRIDGES): digest(BRIDGES),
            relative(THETA): digest(THETA),
            relative(BASE_SCRIPT): digest(BASE_SCRIPT),
            relative(SEARCH_SCRIPT): digest(SEARCH_SCRIPT),
            relative(CORE_SCRIPT): digest(CORE_SCRIPT),
            relative(REVERSE_SCRIPT): digest(REVERSE_SCRIPT),
        },
        "corridors": rows,
        "proof_boundary": {
            "proved": (
                "The stored good-prime paths start at canonical representatives "
                "of the forced genera. H3 and NS0024 end at new rootless cores "
                "with exact zero-mask completions. The Q80 path ends at a "
                "rootless core with exactly two occupied cells on its sole "
                "viable bridge mask."
            ),
            "not_proved": (
                "The replay does not certify completeness of any beam or genus, "
                "nor nonexistence of a Q80 completion. The 42300-neighbour Q80 "
                "miss remains a bounded search observation."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_integral_rank_transfer_masked_core_controls.sage --check"
        ),
    }
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print("PASS prospective masked-core controls")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
