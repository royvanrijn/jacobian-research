#!/usr/bin/env sage-python
"""Search Q80 neighbors whose defining line kills the exact masked defect."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import random as python_random
import runpy

from sage.all import QQ, ZZ, lcm, matrix, pari, set_random_seed, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
BRIDGES = GENERATED / "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
THETA = GENERATED / "elkies-k3-integral-rank-transfer-theta-convolution-v1.json"
BASE_SCRIPT = ROOT / "elkies-k3/scripts/generate_integral_rank_transfer_masked_core_neighbors.sage"
SEARCH_SCRIPT = ROOT / "elkies-k3/scripts/search_integral_rank_transfer_masked_core_controls.sage"
CONTROL_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_masked_core_controls.sage"
CORE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_core_generation.sage"
REVERSE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_reverse_theta_masks.sage"
OUTPUT = GENERATED / "elkies-k3-integral-rank-transfer-q80-defect-neighbors-v1.json"

RANDOM_SEED = 271828
PRIMES = (7, 11, 13, 17, 19, 29, 31, 37, 41, 43)
DEFAULT_SAMPLES = 10000


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dual_witnesses(gram, residue, norm, reverse):
    """Enumerate exact dual pairings for one occupied theta cell."""

    inverse = gram.inverse()
    denominator = lcm(value.denominator() for value in inverse.list())
    scaled_inverse = (denominator * inverse).change_ring(ZZ)
    enumeration = pari(scaled_inverse).qfminim(2 * denominator)
    representatives = matrix(ZZ, enumeration[2].sage()).columns()
    answer = []
    target = tuple(map(QQ, residue))
    for representative in representatives:
        for pairing in (vector(ZZ, representative), -vector(ZZ, representative)):
            dual_vector = pairing * inverse
            if reverse["discriminant_class"](dual_vector) != target:
                continue
            if pairing * inverse * pairing != norm:
                continue
            answer.append(pairing)
    return answer


def masked_witness_data(gram, bridge, order, base, reverse):
    """Return one bridge mask, its occupied cells, and all physical witnesses."""

    masks, _, generator = base["mask_profile"](gram, [bridge], reverse)
    assert len(masks) == 1
    multiplier = masks[0]["isotropic_multipliers"][0]
    oracle = reverse["CoreCellOracle"](gram)
    mask = reverse["reverse_mask"](
        oracle,
        bridge["theta_profile"],
        generator,
        bridge["generator"],
        multiplier,
        order,
    )
    violations = [
        row for row in mask["requirements"] if row["observed_core_occupied"]
    ]
    witnesses = []
    for violation in violations:
        residue = tuple(QQ(value) for value in violation["core_discriminant_class"])
        norm = QQ(violation["required_core_norm"])
        cell_witnesses = dual_witnesses(gram, residue, norm, reverse)
        assert cell_witnesses
        witnesses.extend(cell_witnesses)
    unique = []
    seen = set()
    for witness in witnesses:
        key = tuple(witness)
        if key not in seen:
            seen.add(key)
            unique.append(witness)
    keys = {tuple(witness) for witness in unique}
    assert all(tuple(-witness) in keys for witness in unique)
    return masks[0], violations, unique


def initial_q80(base, search, control, core, reverse):
    bridge_artifact = json.loads(BRIDGES.read_text())
    theta_artifact = json.loads(THETA.read_text())
    prepared = search["prepare_corridor"](
        "Q80", bridge_artifact, theta_artifact, base, core, reverse
    )
    search["configure_order"](base, prepared["order"])
    gram, states = control["replay"](
        prepared["seed"], control["PATHS"]["Q80"], base
    )
    bridge = prepared["viable_bridges"][0]
    mask, violations, unique = masked_witness_data(
        gram, bridge, prepared["order"], base, reverse
    )
    assert mask["occupied_forbidden_cells"] == 2
    assert len(violations) == 2
    assert unique
    return prepared, bridge, gram, states, violations, unique


def line_key(witness, prime):
    values = [int(value % prime) for value in witness]
    pivot = next(value for value in values if value)
    inverse = pow(pivot, -1, prime)
    return tuple((value * inverse) % prime for value in values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=DEFAULT_SAMPLES)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    if arguments.samples <= 0:
        parser.error("--samples must be positive")

    base = runpy.run_path(str(BASE_SCRIPT))
    search = runpy.run_path(str(SEARCH_SCRIPT))
    control = runpy.run_path(str(CONTROL_SCRIPT))
    core = runpy.run_path(str(CORE_SCRIPT))
    reverse = runpy.run_path(str(REVERSE_SCRIPT))
    prepared, bridge, gram, prefix_states, violations, witnesses = initial_q80(
        base, search, control, core, reverse
    )
    form = base["quadratic_form"](gram)
    set_random_seed(RANDOM_SEED)
    python_random.seed(RANDOM_SEED)

    line_keys = set()
    prime_counts = Counter()
    root_counts = Counter()
    mask_counts = Counter()
    nonorthogonal_lines = 0
    duplicate_lines = 0
    rootless_neighbors = 0
    replacement_defects = 0
    hit = None

    for sample_index in range(1, arguments.samples + 1):
        if sample_index % 1000 == 0:
            print(
                f"sample={sample_index} directed={nonorthogonal_lines} "
                f"rootless={rootless_neighbors} masks={dict(mask_counts)}",
                flush=True,
            )
        prime = python_random.choice(PRIMES)
        raw_witness = form.find_primitive_p_divisible_vector__random(prime)
        key = (prime, line_key(raw_witness, prime))
        if key in line_keys:
            duplicate_lines += 1
            continue
        line_keys.add(key)
        pairings = [int(value * raw_witness) % prime for value in witnesses]
        assert sorted(pairings) == sorted((-value) % prime for value in pairings)
        if any(value == 0 for value in pairings):
            continue
        nonorthogonal_lines += 1
        prime_counts[prime] += 1

        transform = form.find_p_neighbor_from_vec(
            prime, raw_witness, return_matrix=True
        )
        neighbor = form.find_p_neighbor_from_vec(prime, raw_witness)
        neighbor_gram = neighbor.Hessian_matrix()

        # Exact survival lemma: an old dual vector belongs to N^dual exactly
        # when its pairing with the defining isotropic line vanishes mod p.
        ambient_basis = transform.transpose()
        for pairing in witnesses:
            dual_vector = pairing * gram.inverse()
            neighbor_pairings = ambient_basis * gram * dual_vector.column()
            actually_survives = all(value in ZZ for value in neighbor_pairings)
            predicted_survival = int(pairing * raw_witness) % prime == 0
            assert actually_survives == predicted_survival
            assert not actually_survives

        roots = int(pari(neighbor_gram).qfminim(2)[0])
        root_counts[roots] += 1
        if roots:
            continue
        rootless_neighbors += 1
        masks, _, _ = base["mask_profile"](
            neighbor_gram, [bridge], reverse
        )
        defect = masks[0]["occupied_forbidden_cells"]
        mask_counts[defect] += 1
        if defect:
            replacement_defects += 1
            continue
        hit = {
            "prime": prime,
            "witness": list(map(int, raw_witness)),
            "reduced_core": base["lll_reduce"](neighbor_gram),
            "mask": masks[0],
        }
        break

    completion = None
    if hit is not None:
        reduced = hit["reduced_core"]
        final_masks, _, _ = base["mask_profile"](reduced, [bridge], reverse)
        assert final_masks[0]["zero_mask_accepts"]
        completion = control["completion"](
            reduced, prepared, final_masks[0], base, core
        )
        historical_isometric = bool(
            pari(reduced).qfisom(
                pari(base["lll_reduce"](prepared["historical_core"]))
            )
        )
        hit_record = {
            "prime": hit["prime"],
            "witness": hit["witness"],
            "reduced_core_gram": [
                list(map(int, row)) for row in reduced.rows()
            ],
            "historical_core_isometric": historical_isometric,
            "completion": completion,
        }
    else:
        hit_record = None

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-q80-defect-neighbors.v1",
        "status": (
            "PASS_DEFECT_DIRECTED_Q80_COMPLETION"
            if hit is not None
            else "PASS_BOUNDED_DEFECT_DIRECTED_Q80_MISS"
        ),
        "inputs": {
            relative(BRIDGES): digest(BRIDGES),
            relative(THETA): digest(THETA),
            relative(BASE_SCRIPT): digest(BASE_SCRIPT),
            relative(SEARCH_SCRIPT): digest(SEARCH_SCRIPT),
            relative(CONTROL_SCRIPT): digest(CONTROL_SCRIPT),
            relative(CORE_SCRIPT): digest(CORE_SCRIPT),
            relative(REVERSE_SCRIPT): digest(REVERSE_SCRIPT),
        },
        "initial_core": {
            "corridor": "Q80",
            "determinant": int(gram.det()),
            "minimum": 4,
            "signed_root_count": 0,
            "canonical_seed_prefix": prefix_states,
            "viable_bridge_class_index": bridge["bridge_class_index"],
            "occupied_mask_cells": violations,
            "offending_dual_pairings": [
                list(map(int, witness)) for witness in witnesses
            ],
        },
        "survival_theorem": {
            "statement": (
                "For N=M_p(K;l)+Z*y/p and x in K^dual, x lies in N^dual "
                "if and only if <x,y> is zero modulo p."
            ),
            "proof": (
                "Pairing x with M_p(K;l) is integral because M_p is contained "
                "in K. The only additional generator is y/p, whose pairing "
                "with x is integral exactly when <x,y> is divisible by p."
            ),
            "checked_on_every_constructed_neighbor": True,
        },
        "search": {
            "random_seed": RANDOM_SEED,
            "requested_samples": arguments.samples,
            "primes": list(PRIMES),
            "unique_isotropic_lines": len(line_keys),
            "duplicate_lines": duplicate_lines,
            "old_defect_killing_lines": nonorthogonal_lines,
            "old_defect_killing_lines_by_prime": {
                str(key): prime_counts[key] for key in sorted(prime_counts)
            },
            "root_count_histogram": {
                str(key): root_counts[key] for key in sorted(root_counts)
            },
            "rootless_neighbors": rootless_neighbors,
            "new_mask_defect_histogram": {
                str(key): mask_counts[key] for key in sorted(mask_counts)
            },
            "rootless_neighbors_with_replacement_defects": replacement_defects,
            "zero_mask_hit": hit is not None,
        },
        "hit": hit_record,
        "proof_boundary": (
            "The survival equivalence is general and exact. The search samples "
            "only the declared isotropic lines. Killing every old witness is "
            "necessary but not sufficient because a neighbor may acquire new "
            "vectors in the same masked discriminant cells. A miss is not a "
            "genus-wide obstruction."
        ),
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "search_integral_rank_transfer_q80_defect_neighbors.sage --check"
        ),
    }

    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print(payload["status"])
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
