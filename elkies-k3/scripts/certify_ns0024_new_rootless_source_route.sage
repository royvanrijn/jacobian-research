#!/usr/bin/env sage-python
"""Certify the source-level route to the new determinant-950 rootless frame.

This separates three facts which must not be conflated: the exact new frame,
the Kneser path of completed positive-definite frames, and the equation-level
Inose source.  It does not assert that the Kneser steps are elliptic-neighbour
maps or supply rational maps from the Inose pencil to the new fibration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import runpy

from sage.all import Gamma0, ModularSymbols, QQ, ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
MASKED = GENERATED / "elkies-k3-integral-rank-transfer-masked-core-controls-v1.json"
BRIDGES = GENERATED / "elkies-k3-integral-rank-transfer-bridge-reglue-v1.json"
THETA = GENERATED / "elkies-k3-integral-rank-transfer-theta-convolution-v1.json"
FOUNDRY = GENERATED / "elkies-k3-lattice-foundry-v1.json"
SOURCES = GENERATED / "elkies-k3-lattice-foundry-prescribed-root-sources-mw0-mw1-group-b-v1.json"
POLES = GENERATED / "elkies-k3-lattice-foundry-rank1-section-poles-v1.json"
SOURCE_FRAME = GENERATED / "elkies-k3-ns0024-2e8-source-root-adapted.txt"
DEGREE2 = GENERATED / "elkies-k3-ns0024-2e8-zero-mw-degree2-q40-v1.json"
DEGREE3 = GENERATED / "elkies-k3-ns0024-2e8-zero-mw-degree3-q30-v1.json"
BASE_SCRIPT = ROOT / "elkies-k3/scripts/generate_integral_rank_transfer_masked_core_neighbors.sage"
SEARCH_SCRIPT = ROOT / "elkies-k3/scripts/search_integral_rank_transfer_masked_core_controls.sage"
CORE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_core_generation.sage"
REVERSE_SCRIPT = ROOT / "elkies-k3/scripts/certify_integral_rank_transfer_reverse_theta_masks.sage"
OUTPUT = GENERATED / "elkies-k3-ns0024-new-rootless-source-route-v1.json"
PATH = (
    (17, (14, 6, 3, 1, 7, 9, 3, 15, 2, 0, 6, 1, 12, 12, 14)),
    (13, (0, 4, 5, 1, 4, 6, 9, 8, 7, 6, 1, 5, 6, 0, 12)),
    (7, (0, 3, 4, 3, 3, 0, 6, 5, 6, 3, 5, 4, 5, 3, 3)),
)


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def root_data(gram):
    result = pari(gram).qfminim(2)
    count = int(result[0])
    if count == 0:
        return 0, 0, 1
    half = matrix(ZZ, result[2])
    basis = half.transpose().row_module().basis_matrix()
    root_gram = basis * gram * basis.transpose()
    return int(basis.rank()), count, int(abs(root_gram.det()))


def root_name(data):
    return {
        (13, 280, 4): "D5+E8",
        (5, 12, 24): "3A1+A2",
        (0, 0, 1): "rootless",
    }.get(tuple(data), "unclassified")


def completed_frames(core_gram, bridge, order, base, core):
    generator = base["primary_generator"](core_gram, order)
    result = []
    for multiplier in range(1, order):
        glue = vector(QQ, list(multiplier * generator) + list(bridge["generator"]))
        try:
            child = core["glued_frame"](core_gram, bridge["gram"], glue)
        except AssertionError:
            continue
        assert child.det() == 950
        assert core["minimum_norm"](child) in (2, 4)
        result.append((multiplier, child))
    assert result
    return result


def summarize_shell(payload, expected_degree, expected_q):
    summaries = payload["summaries"]
    assert [row["q"] for row in summaries] == list(expected_q)
    assert all(row["factor_order"][1] == expected_degree for row in summaries)
    assert all(row["mw_enumeration_complete"] for row in summaries)
    assert all(row["dominant_orbits_complete"] for row in summaries)
    assert all(row["mw_projection_representatives"] == 1 for row in summaries)
    assert all(row["mw_pari_vector_count"] == 0 for row in summaries)
    assert all(
        entry["root_rank"] == 16
        for row in summaries
        for entry in row["root_histogram"]
    )
    return {
        "degree": expected_degree,
        "q_values": list(expected_q),
        "dominant_orbits": sum(row["dominant_orbits"] for row in summaries),
        "primitive_neighbors": sum(row["primitive_neighbors"] for row in summaries),
        "rank_growing_neighbors": 0,
        "mw_projection": "zero only",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    arguments = parser.parse_args()

    masked = json.loads(MASKED.read_text())
    bridges = json.loads(BRIDGES.read_text())
    theta = json.loads(THETA.read_text())
    foundry = json.loads(FOUNDRY.read_text())
    sources = json.loads(SOURCES.read_text())
    poles = json.loads(POLES.read_text())
    degree2 = json.loads(DEGREE2.read_text())
    degree3 = json.loads(DEGREE3.read_text())
    base = runpy.run_path(str(BASE_SCRIPT))
    search = runpy.run_path(str(SEARCH_SCRIPT))
    core = runpy.run_path(str(CORE_SCRIPT))
    reverse = runpy.run_path(str(REVERSE_SCRIPT))

    prepared = search["prepare_corridor"](
        "NS0024", bridges, theta, base, core, reverse
    )
    search["configure_order"](base, prepared["order"])
    masked_ns0024 = next(row for row in masked["corridors"] if row["corridor"] == "NS0024")
    bridge_index = masked_ns0024["completion"]["bridge_class_index"]
    multiplier = masked_ns0024["completion"]["glue_multiplier"]
    bridge = next(
        row for row in prepared["viable_bridges"]
        if row["bridge_class_index"] == bridge_index
    )

    quadratic = base["quadratic_form"](prepared["seed"])
    cores = [base["lll_reduce"](quadratic.Hessian_matrix())]
    for prime, raw_witness in PATH:
        witness = vector(ZZ, raw_witness)
        assert quadratic(witness) % prime == 0
        quadratic = quadratic.find_p_neighbor_from_vec(prime, witness)
        cores.append(base["lll_reduce"](quadratic.Hessian_matrix()))

    completion_choices = [
        completed_frames(value, bridge, prepared["order"], base, core)
        for value in cores
    ]
    expected_roots = ((13, 280, 4), (5, 12, 24), (5, 12, 24), (0, 0, 1))
    completed = []
    selected_multipliers = []
    valid_multipliers = []
    for index, choices in enumerate(completion_choices):
        valid_multipliers.append([value for value, _ in choices])
        matching = [(value, child) for value, child in choices if root_data(child) == expected_roots[index]]
        assert matching
        if index == len(completion_choices) - 1:
            selected = next(row for row in matching if row[0] == multiplier)
        else:
            selected = min(matching, key=lambda row: row[0])
        selected_multipliers.append(selected[0])
        completed.append(selected[1])
    completed_rows = []
    for index, child in enumerate(completed):
        assert core["discriminant_form_key"](child) == core["discriminant_form_key"](
            prepared["target_frame"]
        )
        data = root_data(child)
        completed_rows.append(
            {
                "stage": index,
                "incoming_core_neighbor_prime": None if index == 0 else PATH[index - 1][0],
                "root_type": root_name(data),
                "root_rank": data[0],
                "signed_root_count": data[1],
                "root_determinant": data[2],
                "mw_rank_for_rho_19": 17 - data[0],
                "determinant": int(child.det()),
                "discriminant_form_matches_ns0024": True,
                "glue_multiplier": selected_multipliers[index],
                "all_integral_glue_multipliers": valid_multipliers[index],
            }
        )
    assert [row["root_type"] for row in completed_rows] == [
        "D5+E8", "3A1+A2", "3A1+A2", "rootless"
    ]

    new_frame = base["lll_reduce"](completed[-1])
    assert new_frame.det() == 950
    assert int(pari(new_frame).qfminim(2)[0]) == 0
    norm4 = int(pari(new_frame).qfminim(4)[0])
    assert norm4 == 2634
    automorphisms = int(pari(new_frame).qfauto()[0])
    assert automorphisms == 2
    assert core["discriminant_form_key"](new_frame) == core["discriminant_form_key"](
        prepared["target_frame"]
    )

    catalog = {
        row["frame_id"]: row["invariants"]["norm_four_vectors"]
        for row in foundry["rootless_targets"]
        if row["ns_id"] == "NS0024"
    }
    assert catalog == {"NS0024-F001": 2640, "NS0024-F002": 2630, "NS0024-F005": 2632}
    assert norm4 not in catalog.values()

    source_entry = next(row for row in sources["sources"] if row["source_id"] == "NS0024-S001")
    source = source_entry["source"]
    source_frame = load_gram(SOURCE_FRAME)
    assert source_frame == matrix(ZZ, source["root_adapted_gram"])
    assert source["root_type"] == "2E8"
    assert source["root_rank"] == 16 and source["root_lattice_primitive"]
    assert source["torsion"] == 1 and source["mw_height_gram"] == [["950"]]
    pole = next(
        row for row in poles["sources"]
        if row.get("source_id") == "NS0024-S001" and row.get("root_type") == "2E8"
    )
    assert pole["minimum_section_frame_norm"] == 950
    assert pole["minimum_section_pole_order"] == 473

    group = Gamma0(475)
    symbols = ModularSymbols(475, 2, sign=0).cuspidal_subspace()
    traces = {str(divisor): int(symbols.atkin_lehner_operator(divisor).trace()) for divisor in (19, 25, 475)}
    assert traces == {"19": -6, "25": 2, "475": -14}
    curve = {
        "level": 475,
        "index": int(group.index()),
        "genus": int(group.genus()),
        "cusps": int(group.ncusps()),
        "elliptic_points_order_2": int(group.nu2()),
        "elliptic_points_order_3": int(group.nu3()),
        "weight_2_cuspidal_modular_symbol_dimension": int(symbols.dimension()),
        "modular_symbol_atkin_lehner_traces": traces,
        "differential_atkin_lehner_traces": {key: value // 2 for key, value in traces.items()},
        "fricke_quotient_genus": 19,
        "full_atkin_lehner_quotient_genus": 9,
    }
    assert curve["index"] == 600 and curve["genus"] == 45 and curve["cusps"] == 12
    assert (curve["genus"] + curve["differential_atkin_lehner_traces"]["475"]) // 2 == 19
    assert (
        curve["genus"]
        + curve["differential_atkin_lehner_traces"]["19"]
        + curve["differential_atkin_lehner_traces"]["25"]
        + curve["differential_atkin_lehner_traces"]["475"]
    ) // 4 == 9

    input_paths = (
        MASKED, BRIDGES, THETA, FOUNDRY, SOURCES, POLES, SOURCE_FRAME,
        DEGREE2, DEGREE3, BASE_SCRIPT, SEARCH_SCRIPT, CORE_SCRIPT, REVERSE_SCRIPT,
    )
    payload = {
        "schema": "elkies-k3.ns0024-new-rootless-source-route.v1",
        "status": "PASS_EXACT_FRAME_AND_SOURCE_ROUTE_WITH_OPEN_EQUATION_TRANSPORT",
        "new_rootless_frame": {
            "rank": 17,
            "determinant": 950,
            "minimum": 4,
            "signed_root_count": 0,
            "norm_four_vectors": norm4,
            "automorphism_group_order": automorphisms,
            "gram": rows(new_frame),
            "distinct_from_catalog_by_norm_four_count": catalog,
        },
        "completed_core_path": {
            "cyclic_bridge_order": int(prepared["order"]),
            "bridge_class_index": bridge_index,
            "terminal_glue_multiplier": multiplier,
            "stages": completed_rows,
            "classification": "exact Kneser path of completed frames; not an elliptic-neighbour route",
        },
        "equation_source": {
            "source_id": "NS0024-S001",
            "frame_type": "2E8/MW1",
            "mw_height": 950,
            "minimum_section_pole_order": 473,
            "modular_base": curve,
            "inose_equation": {
                "elliptic_curves": [
                    "E1: y1^2=x1^3+a2*x1^2+a4*x1+a6",
                    "E2: y2^2=x2^3+a2p*x2^2+a4p*x2+a6p",
                ],
                "A": "(a2^2-3*a4)*(a2p^2-3*a4p)",
                "B": "(32/27)*(2*a2^3-9*a2*a4+27*a6)*(2*a2p^3-9*a2p*a4p+27*a6p)",
                "surface": "Y^2=X^3-(A/3)*X+(Delta1*s+B+Delta2/s)/64",
                "condition": "E1 and E2 are non-isomorphic and joined by a cyclic 475-isogeny",
                "mw_interpretation": "Hom(E1,E2)<2>; the degree-475 isogeny has height 950",
                "reference": "https://arxiv.org/abs/2209.02463",
            },
        },
        "bounded_source_neighbor_audit": [
            summarize_shell(degree2, 2, range(4, 41, 2)),
            summarize_shell(degree3, 3, range(3, 31, 3)),
        ],
        "inputs": {relative(path): digest(path) for path in input_paths},
        "proof_boundary": {
            "proved": (
                "The determinant-950 completion is rootless with minimum 4 and 2634 norm-four vectors; "
                "it is distinct from the three catalogued NS0024 rootless frames by that invariant. "
                "Completing every core on the stored 17,13,7 path with the same order-191 bridge gives "
                "D5+E8/MW4, 3A1+A2/MW12, 3A1+A2/MW12, rootless/MW17. The 2E8/MW1 "
                "source has height 950, and exhaustive zero-MW Weyl-orbit audits find no rank growth "
                "in the stated degree-two and degree-three q ranges."
            ),
            "literature_route": (
                "Utsumi's Inose equation supplies a symbolic equation over the level-475 isogeny locus."
            ),
            "not_proved": (
                "No explicit pair of 475-isogenous elliptic curves over QQ, rational parameterization of "
                "X0(475), equation for the D5+E8 fibration, elliptic-neighbour pencils along the Kneser "
                "path, rational maps, or equation for the new rootless fibration is supplied."
            ),
        },
        "reproduce": [
            "sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage --frame artifacts/generated-results/elkies-k3-ns0024-2e8-source-root-adapted.txt --root-rank 16 --degree 2 --include-zero-mw "
            + " ".join(f"--q {q}" for q in range(4, 41, 2))
            + " --output artifacts/generated-results/elkies-k3-ns0024-2e8-zero-mw-degree2-q40-v1.json",
            "sage -python elkies-k3/scripts/search_root_adapted_weyl_neighbors.sage --frame artifacts/generated-results/elkies-k3-ns0024-2e8-source-root-adapted.txt --root-rank 16 --degree 3 --include-zero-mw "
            + " ".join(f"--q {q}" for q in range(3, 31, 3))
            + " --output artifacts/generated-results/elkies-k3-ns0024-2e8-zero-mw-degree3-q30-v1.json",
            "sage -python elkies-k3/scripts/certify_ns0024_new_rootless_source_route.sage --check",
        ],
    }
    output = arguments.output if arguments.output.is_absolute() else ROOT / arguments.output
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit(f"missing or stale artifact: {output}")
        print("PASS NS0024 new rootless source route")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
    print(relative(output))


if __name__ == "__main__":
    main()
