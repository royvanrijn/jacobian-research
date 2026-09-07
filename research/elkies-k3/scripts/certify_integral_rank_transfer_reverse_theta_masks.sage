#!/usr/bin/env sage-python
"""Certify bridge-induced reverse theta masks for rank-15 cores.

For a fixed graph H in A_K+A_C, nonnegativity of theta convolution turns
rootlessness into coefficientwise vanishings on K.  This checker derives the
required core cells from C and H alone, then consults K only to test those
cells.  It replays every oriented terminal graph without constructing a
rank-17 child or evaluating a full convolution.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from time import perf_counter

from sage.all import QQ, ZZ, lcm, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
BRIDGES = GENERATED / "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
THETA = GENERATED / "elkies-k3-integral-rank-transfer-theta-convolution-v1.json"
OUTPUT = GENERATED / "elkies-k3-integral-rank-transfer-reverse-theta-masks-v1.json"


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fractional_part(value):
    value = QQ(value)
    return value - value.floor()


def discriminant_class(coordinates):
    return tuple(fractional_part(value) for value in coordinates)


def canonical_signed_class(residue):
    residue = tuple(map(QQ, residue))
    negative = tuple(fractional_part(-value) for value in residue)
    return min(residue, negative)


def theta_profile(gram):
    """Return every occupied dual-coset theta cell through norm two."""

    inverse = gram.inverse()
    denominator = lcm(value.denominator() for value in inverse.list())
    scaled_inverse = (denominator * inverse).change_ring(ZZ)
    enumeration = pari(scaled_inverse).qfminim(2 * denominator)
    representatives = matrix(ZZ, enumeration[2].sage()).columns()
    profile = defaultdict(Counter)
    zero = vector(QQ, [0] * gram.nrows())
    profile[discriminant_class(zero)][QQ(0)] = 1
    for representative in representatives:
        for dual_pairing in (
            vector(ZZ, representative),
            -vector(ZZ, representative),
        ):
            dual_vector = dual_pairing * inverse
            norm = dual_pairing * inverse * dual_pairing
            assert 0 < norm <= 2
            profile[discriminant_class(dual_vector)][norm] += 1
    assert sum(sum(row.values()) for row in profile.values()) == int(enumeration[0]) + 1
    return profile


def bridge_generator(gram):
    generator = vector(ZZ, [1, 0]) * gram.inverse()
    if generator in ZZ**2:
        generator = vector(ZZ, [0, 1]) * gram.inverse()
    return generator


class CoreCellOracle:
    """Exact lazy occupancy oracle for selected theta cells of one core."""

    def __init__(self, gram):
        self.gram = gram
        self.lll = matrix(ZZ, pari(gram).qflllgram()).transpose()
        self.reduced = self.lll * gram * self.lll.transpose()
        self.inverse_lll = self.lll.inverse().change_ring(ZZ)
        self.cache = {}
        self.rootful = None

    def occupied(self, residue, norm):
        residue = canonical_signed_class(residue)
        key = (residue, QQ(norm))
        if key in self.cache:
            return self.cache[key]
        residue = vector(QQ, residue)
        norm = QQ(norm)
        is_zero = residue in ZZ**self.gram.nrows()
        if is_zero:
            if norm == 0:
                answer = True
            elif norm == 2:
                if self.rootful is None:
                    self.rootful = bool(int(pari(self.gram).qfminim(2)[0]))
                answer = self.rootful
            else:
                answer = False
        else:
            # Norms in one discriminant coset agree modulo 2.  Therefore the
            # requested norm in [0,2] is occupied iff it is the exact coset
            # minimum; no higher shell must be enumerated.
            target = -residue * self.inverse_lll
            target -= vector(ZZ, [value.floor() for value in target])
            # PARI's bound is evaluated in floating point internally.  A tiny
            # rational guard avoids losing vectors exactly on the boundary;
            # exact QQ reevaluation below remains the authority.
            bounded = pari(self.reduced).qfcvp(
                pari(target), B=norm + QQ(1) / 10**8, flag=0
            )
            candidates = matrix(ZZ, bounded[2].sage()).columns()
            exact_norms = [
                (vector(QQ, item) - target)
                * self.reduced
                * (vector(QQ, item) - target)
                for item in candidates
            ]
            assert all(value <= norm for value in exact_norms)
            answer = norm in exact_norms
        self.cache[key] = bool(answer)
        return bool(answer)


def reverse_mask(
    core_oracle,
    bridge_profile,
    core_generator,
    bridge_generator_value,
    multiplier,
    order,
):
    """Derive required zero cells from C,H, then query only those K cells."""

    requirements = []
    for label in range(order):
        core_class = discriminant_class(label * multiplier * core_generator)
        bridge_class = discriminant_class(label * bridge_generator_value)
        for bridge_norm, bridge_multiplicity in sorted(
            bridge_profile.get(bridge_class, {}).items()
        ):
            core_norm = QQ(2) - bridge_norm
            if not 0 <= core_norm <= 2:
                continue
            core_occupied = core_oracle.occupied(core_class, core_norm)
            requirements.append(
                {
                    "graph_label": label,
                    "core_discriminant_class": [str(value) for value in core_class],
                    "required_core_norm": str(core_norm),
                    "bridge_discriminant_class": [
                        str(value) for value in bridge_class
                    ],
                    "bridge_norm": str(bridge_norm),
                    "bridge_multiplicity": int(bridge_multiplicity),
                    "observed_core_occupied": core_occupied,
                }
            )
    queried_cells = {
        (tuple(row["core_discriminant_class"]), row["required_core_norm"])
        for row in requirements
    }
    signed_cells = {
        (
            canonical_signed_class(
                QQ(value) for value in row["core_discriminant_class"]
            ),
            QQ(row["required_core_norm"]),
        )
        for row in requirements
    }
    violations = [row for row in requirements if row["observed_core_occupied"]]
    return {
        "reverse_zero_mask_cell_count": len(queried_cells),
        "sign_symmetry_reduced_mask_cell_count": len(signed_cells),
        "bridge_theta_terms_inducing_mask": len(requirements),
        "mask_fraction_of_graph_classes": str(QQ(len(queried_cells)) / order),
        "occupied_forbidden_cell_count": len(violations),
        "zero_mask_accepts": not violations,
        "requirements": requirements,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    bridges = json.loads(BRIDGES.read_text())
    theta = json.loads(THETA.read_text())
    assert bridges["status"] == "PASS_EXACT_BRIDGE_REGLUE_CERTIFICATES"
    assert theta["status"] == "PASS_EXACT_THETA_CONVOLUTION_ZERO_SUPPORT_ENUMERATOR"
    terminal_edges = {
        row["corridor"]: row
        for row in bridges["edges"]
        if int(row["target_root_rank"]) == 0
    }
    theta_by_corridor = {row["corridor"]: row for row in theta["corridors"]}
    assert set(terminal_edges) == set(theta_by_corridor)

    corridors = []
    all_graphs = []
    for corridor in sorted(terminal_edges):
        edge = terminal_edges[corridor]
        stored = theta_by_corridor[corridor]
        order = int(stored["cyclic_glue_order"])
        core = matrix(ZZ, edge["core"]["gram"])
        split_generator = vector(
            QQ,
            [
                QQ(value)
                for value in edge["new_frame"]["glue_generators"][0][
                    "K_plus_C_dual_coordinates"
                ]
            ],
        )
        core_generator = split_generator[:-2]
        core_oracle = CoreCellOracle(core)
        targeted_started = perf_counter()

        classes = []
        for bridge_class in stored["classes"]:
            bridge = matrix(ZZ, bridge_class["bridge_gram"])
            generator = bridge_generator(bridge)
            bridge_profile = theta_profile(bridge)
            graphs = []
            for graph in bridge_class["admissible_oriented_graphs"]:
                record = reverse_mask(
                    core_oracle,
                    bridge_profile,
                    core_generator,
                    generator,
                    int(graph["core_glue_multiplier"]),
                    order,
                )
                graphs.append(
                    {
                        "core_glue_multiplier": int(graph["core_glue_multiplier"]),
                        "_stored_signed_root_count": int(
                            graph["predicted_signed_root_count"]
                        ),
                        "_stored_zero_support_accepts": bool(
                            graph["zero_support_accepts"]
                        ),
                        **record,
                    }
                )
            classes.append(
                {
                    "bridge_class_index": int(bridge_class["bridge_class_index"]),
                    "bridge_gram": bridge_class["bridge_gram"],
                    "graphs": graphs,
                }
            )
        targeted_seconds = perf_counter() - targeted_started

        # Independent truth phase: only after every zero-mask decision is
        # fixed, enumerate the complete core table to recover multiplicities
        # and compare exact root counts with the earlier convolution record.
        full_started = perf_counter()
        core_profile = theta_profile(core)
        full_seconds = perf_counter() - full_started
        occupied_core_cells = sum(len(row) for row in core_profile.values())
        for bridge_class in classes:
            for graph in bridge_class["graphs"]:
                predicted_root_count = 0
                for requirement in graph["requirements"]:
                    core_class = tuple(
                        QQ(value)
                        for value in requirement["core_discriminant_class"]
                    )
                    core_norm = QQ(requirement["required_core_norm"])
                    core_multiplicity = core_profile.get(core_class, {}).get(
                        core_norm, 0
                    )
                    if bool(core_multiplicity) != requirement[
                        "observed_core_occupied"
                    ]:
                        raise ArithmeticError(
                            f"targeted theta mismatch {corridor=} "
                            f"{core_class=} {core_norm=} "
                            f"oracle={requirement['observed_core_occupied']} "
                            f"truth={core_multiplicity}"
                        )
                    requirement["truth_core_multiplicity"] = int(core_multiplicity)
                    predicted_root_count += (
                        int(requirement["bridge_multiplicity"])
                        * core_multiplicity
                    )
                assert predicted_root_count == graph.pop(
                    "_stored_signed_root_count"
                )
                assert graph["zero_mask_accepts"] == graph.pop(
                    "_stored_zero_support_accepts"
                )
                graph["truth_signed_root_count"] = int(predicted_root_count)
                graph["fraction_of_occupied_core_theta_cells_queried"] = str(
                    QQ(graph["reverse_zero_mask_cell_count"])
                    / occupied_core_cells
                )
                all_graphs.append(graph)

        signed_masks = {
            frozenset(
                (
                    canonical_signed_class(
                        QQ(value)
                        for value in requirement["core_discriminant_class"]
                    ),
                    QQ(requirement["required_core_norm"]),
                )
                for requirement in graph["requirements"]
            )
            for bridge_class in classes
            for graph in bridge_class["graphs"]
        }
        minimal_signed_masks = [
            mask
            for mask in signed_masks
            if not any(other < mask for other in signed_masks)
        ]
        assert len(signed_masks) == len(classes)
        assert len(minimal_signed_masks) == len(signed_masks)
        corridors.append(
            {
                "corridor": corridor,
                "cyclic_order": order,
                "occupied_core_theta_cells_through_norm_two": occupied_core_cells,
                "distinct_core_theta_cells_queried_lazily": len(core_oracle.cache),
                "oriented_graph_count": sum(
                    len(row["graphs"]) for row in classes
                ),
                "unoriented_sign_reduced_mask_count": len(signed_masks),
                "subsumption_minimal_mask_count": len(minimal_signed_masks),
                "sign_reduced_mask_size_distribution": {
                    str(key): value
                    for key, value in sorted(
                        Counter(len(mask) for mask in signed_masks).items()
                    )
                },
                "timing_seconds": {
                    "lazy_targeted_core_cell_queries": targeted_seconds,
                    "independent_full_core_theta_truth_table": full_seconds,
                    "boundary": "Informative workstation timings, not theorem fields.",
                },
                "classes": classes,
            }
        )

    mask_sizes = Counter(row["reverse_zero_mask_cell_count"] for row in all_graphs)
    signed_mask_sizes = Counter(
        row["sign_symmetry_reduced_mask_cell_count"] for row in all_graphs
    )
    accepted = [row for row in all_graphs if row["zero_mask_accepts"]]
    rejected = [row for row in all_graphs if not row["zero_mask_accepts"]]
    assert len(all_graphs) == 28 and len(accepted) == 10
    assert all(row["occupied_forbidden_cell_count"] == 0 for row in accepted)
    assert all(row["occupied_forbidden_cell_count"] > 0 for row in rejected)

    payload = {
        "schema": "elkies-k3.integral-rank-transfer-reverse-theta-masks.v1",
        "status": "PASS_EXACT_REVERSE_THETA_ZERO_MASKS",
        "inputs": {
            relative(BRIDGES): digest(BRIDGES),
            relative(THETA): digest(THETA),
        },
        "theorem": {
            "reverse_mask": (
                "F_(C,H)={(a,2-mu):(a,b) in H and theta_C(b,mu)>0}."
            ),
            "criterion": (
                "W_H is rootless iff theta_K(a,nu)=0 for every (a,nu) "
                "in F_(C,H)."
            ),
            "reason": (
                "Every convolution summand is a nonnegative integer, so a zero "
                "sum is equivalent to coefficientwise vanishing."
            ),
            "enumeration_consequence": (
                "A bridge and graph prescribe a finite zero mask on the core "
                "before K is constructed; only those core theta cells need be "
                "queried for the rootless decision."
            ),
            "support_hypergraph": (
                "Rootlessness depends only on Boolean theta support. Modulo "
                "theta_K(a,nu)=theta_K(-a,nu), the bridge universe is a finite "
                "antichain of forbidden-support masks; a core is compatible "
                "iff its support avoids at least one mask."
            ),
        },
        "aggregate": {
            "oriented_graphs": len(all_graphs),
            "zero_mask_accepting_orientations": len(accepted),
            "zero_mask_rejecting_orientations": len(rejected),
            "root_count_matches": len(all_graphs),
            "mask_size_distribution": {
                str(key): value for key, value in sorted(mask_sizes.items())
            },
            "minimum_mask_cells": min(mask_sizes),
            "maximum_mask_cells": max(mask_sizes),
            "sign_reduced_mask_size_distribution": {
                str(key): value for key, value in sorted(signed_mask_sizes.items())
            },
            "minimum_sign_reduced_mask_cells": min(signed_mask_sizes),
            "maximum_sign_reduced_mask_cells": max(signed_mask_sizes),
            "unoriented_sign_reduced_masks": sum(
                row["unoriented_sign_reduced_mask_count"] for row in corridors
            ),
            "subsumption_minimal_unoriented_masks": sum(
                row["subsumption_minimal_mask_count"] for row in corridors
            ),
            "rank17_children_constructed": 0,
        },
        "corridors": corridors,
        "proof_boundary": {
            "proved": (
                "The reverse-mask equivalence and its exact agreement with all "
                "28 terminal graph orientations."
            ),
            "not_proved": (
                "No direct lattice realization algorithm for a prescribed "
                "allowed theta signature, constrained-genus mass formula, "
                "unbounded determinant cutoff, or speedup theorem is claimed."
            ),
        },
        "reproduce": (
            "sage -python elkies-k3/scripts/"
            "certify_integral_rank_transfer_reverse_theta_masks.sage --check"
        ),
    }
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    if arguments.check:
        if not output.exists():
            raise SystemExit(f"missing artifact: {output}")
        stored = json.loads(output.read_text())
        stored_corridors = {
            row["corridor"]: row for row in stored["corridors"]
        }
        for row in payload["corridors"]:
            row["timing_seconds"] = stored_corridors[row["corridor"]][
                "timing_seconds"
            ]
        encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        if output.read_text() != encoded:
            raise SystemExit(f"stale artifact: {output}")
        print("PASS reverse theta zero masks")
        return
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
