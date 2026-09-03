#!/usr/bin/env sage-python
"""Certify the defect-directed Q80 zero-mask completion path."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import ZZ, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
DIRECTED_SCRIPT = ROOT / "elkies-k3/scripts/search_integral_rank_transfer_q80_defect_neighbors.sage"
BEAM_SCRIPT = ROOT / "elkies-k3/scripts/search_integral_rank_transfer_q80_defect_beam.sage"
OUTPUT = GENERATED / "elkies-k3-integral-rank-transfer-q80-defect-completion-v1.json"

DIRECTED_PATH = (
    (13, (10, 7, 6, 11, 8, 11, 4, 11, 7, 10, 1, 3, 6, 3, 10)),
    (13, (5, 2, 3, 1, 1, 0, 5, 2, 1, 8, 9, 3, 4, 3, 11)),
    (13, (8, 4, 1, 8, 8, 0, 4, 5, 6, 8, 6, 7, 0, 11, 7)),
    (29, (20, 6, 3, 16, 15, 8, 1, 28, 26, 10, 6, 19, 13, 23, 24)),
)


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    directed = runpy.run_path(str(DIRECTED_SCRIPT))
    base = runpy.run_path(str(directed["BASE_SCRIPT"]))
    search = runpy.run_path(str(directed["SEARCH_SCRIPT"]))
    control = runpy.run_path(str(directed["CONTROL_SCRIPT"]))
    core = runpy.run_path(str(directed["CORE_SCRIPT"]))
    reverse = runpy.run_path(str(directed["REVERSE_SCRIPT"]))
    prepared, bridge, initial_gram, prefix_states, _, _ = directed["initial_q80"](
        base, search, control, core, reverse
    )
    current = base["quadratic_form"](initial_gram)
    transitions = []

    for step, (prime, raw_line) in enumerate(DIRECTED_PATH, start=1):
        parent_gram = current.Hessian_matrix()
        parent_mask, parent_cells, witnesses = directed["masked_witness_data"](
            parent_gram,
            bridge,
            prepared["order"],
            base,
            reverse,
        )
        assert parent_mask["occupied_forbidden_cells"] == 2
        line = vector(ZZ, raw_line)
        assert current(line) % prime == 0
        witness_pairings = [int(witness * line) % prime for witness in witnesses]
        assert all(witness_pairings)

        transform = current.find_p_neighbor_from_vec(
            prime, line, return_matrix=True
        )
        ambient_basis = transform.transpose()
        inverse = parent_gram.inverse()
        for witness in witnesses:
            dual_vector = witness * inverse
            neighbor_pairings = ambient_basis * parent_gram * dual_vector.column()
            assert not all(value in ZZ for value in neighbor_pairings)

        current = current.find_p_neighbor_from_vec(prime, line)
        child_gram = current.Hessian_matrix()
        assert child_gram.det() == initial_gram.det()
        assert int(pari(child_gram).qfminim(2)[0]) == 0
        child_mask, child_cells, child_witnesses = directed["masked_witness_data"](
            child_gram,
            bridge,
            prepared["order"],
            base,
            reverse,
        )
        expected_defect = 0 if step == len(DIRECTED_PATH) else 2
        assert child_mask["occupied_forbidden_cells"] == expected_defect
        transitions.append(
            {
                "step_after_near_miss": step,
                "prime": prime,
                "line": list(raw_line),
                "parent_occupied_cells": len(parent_cells),
                "parent_physical_witnesses": len(witnesses),
                "nonzero_pairings_with_line": witness_pairings,
                "all_parent_witnesses_removed": True,
                "child_occupied_cells": len(child_cells),
                "child_physical_witnesses": len(child_witnesses),
                "replacement_occurred": bool(child_cells),
            }
        )

    final_core = base["lll_reduce"](current.Hessian_matrix())
    final_masks, _, _ = base["mask_profile"](final_core, [bridge], reverse)
    assert len(final_masks) == 1 and final_masks[0]["zero_mask_accepts"]
    completion = control["completion"](
        final_core, prepared, final_masks[0], base, core
    )
    assert not completion["isometric_to_declared_target_frame"]
    historical_core_isometric = bool(
        pari(final_core).qfisom(
            pari(base["lll_reduce"](prepared["historical_core"]))
        )
    )
    assert not historical_core_isometric

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-q80-defect-completion.v1",
        "status": "PASS_EXACT_DEFECT_DIRECTED_Q80_COMPLETION",
        "inputs": {
            relative(DIRECTED_SCRIPT): digest(DIRECTED_SCRIPT),
            relative(BEAM_SCRIPT): digest(BEAM_SCRIPT),
            relative(directed["BRIDGES"]): digest(directed["BRIDGES"]),
            relative(directed["THETA"]): digest(directed["THETA"]),
            relative(directed["BASE_SCRIPT"]): digest(directed["BASE_SCRIPT"]),
            relative(directed["SEARCH_SCRIPT"]): digest(directed["SEARCH_SCRIPT"]),
            relative(directed["CONTROL_SCRIPT"]): digest(directed["CONTROL_SCRIPT"]),
            relative(directed["CORE_SCRIPT"]): digest(directed["CORE_SCRIPT"]),
            relative(directed["REVERSE_SCRIPT"]): digest(directed["REVERSE_SCRIPT"]),
        },
        "survival_theorem": {
            "statement": (
                "If N=M_p(K;l)+Z*y/p and x lies in K^dual, then x lies in "
                "N^dual if and only if <x,y> is zero modulo p."
            ),
            "proof": (
                "The pairing with M_p(K;l), a sublattice of K, is already "
                "integral. Pairing with the remaining generator y/p is "
                "integral exactly when p divides <x,y>."
            ),
        },
        "canonical_seed_prefix": prefix_states,
        "initial_near_miss": {
            "rank": initial_gram.nrows(),
            "determinant": int(initial_gram.det()),
            "minimum": 4,
            "signed_root_count": 0,
            "occupied_mask_cells": 2,
        },
        "directed_transitions": transitions,
        "final_core": {
            "rank": final_core.nrows(),
            "determinant": int(final_core.det()),
            "minimum": 4,
            "signed_root_count": 0,
            "automorphism_group_order": int(pari(final_core).qfauto()[0]),
            "historical_core_isometric": historical_core_isometric,
            "viable_bridge_mask_result": final_masks[0],
            "gram": [list(map(int, row)) for row in final_core.rows()],
        },
        "completion": completion,
        "proof_boundary": {
            "proved": (
                "Every directed edge removes every physical witness of the "
                "parent's two-cell defect. The first three neighbors acquire "
                "replacement witnesses; the fourth acquires none. The final "
                "graph completion is rootless and has the target discriminant form."
            ),
            "not_proved": (
                "There is no monotone scalar defect law, guarantee that a "
                "directed path exists in every forced genus, neighbor-graph "
                "completeness statement, or uniform running-time bound."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_integral_rank_transfer_q80_defect_completion.sage --check"
        ),
    }
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print("PASS exact defect-directed Q80 completion")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
