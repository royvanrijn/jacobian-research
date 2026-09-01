#!/usr/bin/env sage-python
"""Rank exact foundry source certificates by equation-first complexity.

The original foundry shortlist rewarded root rank and target short vectors
before an equation source had been found.  This scorer starts later in the
pipeline: it consumes exact rootful-source certificates, derives only the
source invariants justified by those certificates, and leaves arithmetic
marking and elliptic-neighbour cost explicitly unknown unless separate
artifacts provide them.

The ordering is lexicographic.  Sources of Mordell--Weil rank at most two
form a preferred band.  Inside and after that band the coordinates are MW
rank, number of reducible-fibre supports, semistable compatibility, expected
fibre-stratum dimension, minimum nonzero-section pole order, marking Galois
orbit, expected coefficient conditions, and finally certified neighbour cost.
An audited low-degree multisection vector is used only as a final heuristic
tie-break after corridor cost.

This is a discovery ranking, not a theorem that the coordinates are
statistically predictive or that an A-type root system has a semistable
equation over the desired ground field.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
from pathlib import Path

from sage.all import ZZ, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_PATTERN = (
    "artifacts/generated-results/elkies-k3-lattice-foundry-*-source-*.json"
)
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-source-ranking-v2.json"
)
DATABASE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
MULTISECTION_SPECTRUM = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-multisection-spectrum-v1.json"
)


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def minimum_nonzero_section_pole(source: dict, maximum_norm: int = 40) -> dict:
    """Find the first frame vector outside the root span.

    In root-adapted frame coordinates the first ``root_rank`` coordinates
    span the primitive root lattice.  A vector with a nonzero tail represents
    a nonzero MW class, and its integral frame norm is ``4+2(P.O)`` for the
    associated section representative.  Enumerating complete norm shells
    therefore gives the exact minimum pole order once the first such vector
    appears.
    """

    gram = matrix(ZZ, source["root_adapted_gram"])
    root_rank = int(source["root_rank"])
    for bound in range(2, maximum_norm + 1, 2):
        result = pari(gram).qfminim(bound)
        columns = matrix(ZZ, result[2].sage()).columns()
        candidates = []
        for column in columns:
            value = vector(ZZ, column)
            if any(value[index] for index in range(root_rank, gram.nrows())):
                norm = int(value * gram * value)
                candidates.append((norm, tuple(map(int, value))))
        if candidates:
            norm, witness = min(candidates)
            assert norm >= 4 and (norm - 4) % 2 == 0
            return {
                "status": "PASS_EXACT_COMPLETE_NORM_SHELL",
                "pole_order": (norm - 4) // 2,
                "frame_norm": norm,
                "root_adapted_witness": list(witness),
                "enumerated_through_norm": bound,
            }
    return {
        "status": "OPEN_ABOVE_BOUND",
        "pole_order": None,
        "enumerated_through_norm": maximum_norm,
    }


def route_artifacts() -> dict[str, dict]:
    result = {}
    for manifest_path in sorted(
        (ROOT / "elkies-k3/data/lattice-foundry").glob("*-route-v1.json")
    ):
        manifest = json.loads(manifest_path.read_text())
        if manifest.get("schema") != "elkies-k3.lattice-foundry-route-specification.v1":
            continue
        source_path = (ROOT / manifest["source_hunt"]).resolve()
        output_path = (ROOT / manifest["output"]).resolve()
        if not output_path.is_file():
            continue
        certificate = json.loads(output_path.read_text())
        if "equation_cost_vector" not in certificate:
            continue
        result[str(source_path)] = {
            "status": "PASS_EXACT_CERTIFIED_ROUTE",
            "manifest": relative(manifest_path),
            "certificate": relative(output_path),
            "cost": certificate["equation_cost_vector"],
            "target_frame_ids": [
                row["frame_id"]
                for row in certificate["terminal"].get("catalogue_isometry_matches", [])
            ],
        }
    return result


def route_key(route: dict | None) -> tuple:
    if route is None:
        return (1, 10**9, 10**9, 10**9, 10**9, 10**9, 10**9)
    cost = route["cost"]
    return (
        0,
        int(cost["maximum_old_fibre_degree"]),
        int(cost["maximum_q"]),
        int(cost["physical_weyl_reflections"]),
        int(cost["root_rank_area"]),
        int(cost["edge_count"]),
        int(cost["sum_q"]),
    )


def marking_record(payload: dict) -> dict:
    marking = payload.get("arithmetic_marking", {})
    orbit = marking.get("characteristic_zero_galois_orbit_size")
    if orbit is None:
        return {
            "status": "UNKNOWN_NOT_INFERRED_FROM_LATTICE",
            "characteristic_zero_galois_orbit_size": None,
            "rational_source_marking": None,
            "rational_parameterization": marking.get("rational_parameterization"),
        }
    orbit = int(orbit)
    assert orbit >= 1
    return {
        "status": "PASS_EXPLICIT_ARITHMETIC_MARKING_EVIDENCE",
        "characteristic_zero_galois_orbit_size": orbit,
        "rational_source_marking": orbit == 1,
        "rational_parameterization": marking.get("rational_parameterization"),
        "evidence": marking.get("evidence"),
    }


def multisection_key(spectrum: dict | None) -> tuple:
    if spectrum is None:
        return (1, 0, 0, 0, 0, 0)
    return (
        0,
        -int(spectrum["rational_bisection_orbits_minimum_height"]),
        -int(spectrum["genus_one_bisection_candidate_orbits_minimum_height"]),
        -int(spectrum["sampled_rational_trisection_candidates"]),
        -int(spectrum["sampled_genus_one_trisection_candidates"]),
        -int(spectrum["sampled_low_genus_quadrisection_candidates"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pattern", default=DEFAULT_PATTERN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--database", type=Path, default=DATABASE)
    parser.add_argument(
        "--multisection-spectrum", type=Path, default=MULTISECTION_SPECTRUM
    )
    parser.add_argument("--target-mw-min", type=int, default=15)
    parser.add_argument("--target-mw-max", type=int, default=17)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    paths = [Path(value).resolve() for value in sorted(glob.glob(str(ROOT / arguments.pattern)))]
    if not paths:
        raise SystemExit(f"no source artifacts matched {arguments.pattern!r}")
    routes = route_artifacts()
    database_path = arguments.database.resolve()
    database = json.loads(database_path.read_text())
    multisection_path = arguments.multisection_spectrum.resolve()
    multisection_payload = json.loads(multisection_path.read_text())
    multisection_by_frame = {
        row["frame_id"]: row["richness_coordinates"]
        for row in multisection_payload["targets"]
    }
    high_rank_targets_by_ns = {}
    for ns in database["ns_classes"]:
        targets = []
        for frame in ns["frames"]:
            target_mw = int(frame["mw_rank_for_rho_19"])
            if arguments.target_mw_min <= target_mw <= arguments.target_mw_max:
                targets.append({
                    "frame_id": frame["frame_id"],
                    "mw_rank": target_mw,
                    "root_rank": int(frame["root_rank"]),
                    "root_type": frame["root_type"],
                    "audited_multisection_richness": multisection_by_frame.get(
                        frame["frame_id"]
                    ),
                })
        targets.sort(key=lambda row: (-row["mw_rank"], row["root_rank"], row["frame_id"]))
        high_rank_targets_by_ns[ns["ns_id"]] = targets
    rows = []
    inputs = {}
    for path in paths:
        payload = json.loads(path.read_text())
        if payload.get("schema") != "elkies-k3.lattice-foundry-rootful-source.v1":
            continue
        source = payload["source"]
        components = source["root_components"]
        mw_rank = int(source["mw_rank_for_rho_19"])
        root_rank = int(source["root_rank"])
        assert mw_rank + root_rank == 17
        semistable_compatible = all(row["type"].startswith("A") for row in components)
        pole = minimum_nonzero_section_pole(source)
        marking = marking_record(payload)
        route = routes.get(str(path))
        galois_orbit = marking["characteristic_zero_galois_orbit_size"]
        high_rank_targets = high_rank_targets_by_ns[payload["target"]["ns_id"]]
        assert high_rank_targets
        audited_targets = [
            row for row in high_rank_targets
            if row["audited_multisection_richness"] is not None
        ]
        best_audited_target = min(
            audited_targets,
            key=lambda row: multisection_key(row["audited_multisection_richness"]),
            default=None,
        )
        best_multisection = (
            best_audited_target["audited_multisection_richness"]
            if best_audited_target is not None
            else None
        )
        score = (
            0 if mw_rank <= 2 else 1,
            mw_rank,
            len(components),
            0 if semistable_compatible else 1,
            18 - root_rank,
            pole["pole_order"] if pole["pole_order"] is not None else 10**9,
            0 if galois_orbit is not None else 1,
            galois_orbit if galois_orbit is not None else 10**9,
            mw_rank,
            *route_key(route),
            *multisection_key(best_multisection),
        )
        rows.append(
            {
                "source_artifact": relative(path),
                "ns_id": payload["target"]["ns_id"],
                "target_frame_id": payload["target"]["frame_id"],
                "determinant": int(source["determinant"]),
                "source_root_type": source["root_type"],
                "source_root_rank": root_rank,
                "source_mw_rank": mw_rank,
                "mw_at_most_two_preferred_band": mw_rank <= 2,
                "reducible_fibre_support_count": len(components),
                "semistable_configuration_compatible": semistable_compatible,
                "semistable_status_boundary": (
                    "A-type roots admit a multiplicative Kodaira interpretation, but "
                    "the lattice certificate does not construct or descend that equation."
                ),
                "expected_fibre_stratum_dimension": 18 - root_rank,
                "expected_ns_locus_dimension": 1,
                "one_dimensional_lattice_polarized_family": True,
                "ground_field_and_rational_parameterization_status": (
                    "UNKNOWN_UNLESS_ARITHMETIC_MARKING_EVIDENCE_IS_ATTACHED"
                ),
                "expected_additional_coefficient_conditions": mw_rank,
                "dimension_status": (
                    "EXPECTED_FROM_LATTICE_CODIMENSION_NOT_AN_INDEPENDENCE_PROOF"
                ),
                "minimum_nonzero_section_pole": pole,
                "arithmetic_marking": marking,
                "admissible_high_rank_targets": high_rank_targets,
                "best_audited_multisection_target": best_audited_target,
                "multisection_score_boundary": (
                    "Degree-two coordinates are complete low-height orbit counts. "
                    "Degree-three/four coordinates are exact only in the declared "
                    "sample and are a final discovery tie-break, not a rank theorem."
                ),
                "certified_neighbor_route": route or {
                    "status": "UNKNOWN_NOT_YET_ENUMERATED"
                },
                "score_tuple": list(score),
                "_score": score,
            }
        )
        inputs[relative(path)] = digest(path)

    rows.sort(key=lambda row: (row["_score"], row["ns_id"], row["target_frame_id"]))
    for rank, row in enumerate(rows, 1):
        row["rank"] = rank
        del row["_score"]

    inputs[relative(database_path)] = digest(database_path)
    inputs[relative(multisection_path)] = digest(multisection_path)
    output = {
        "schema": "elkies-k3.lattice-foundry-source-ranking.v2",
        "status": "PASS_EXACT_SOURCE_METRICS_WITH_TYPED_OPEN_ARITHMETIC_AND_ROUTE_GATES",
        "objective_order": [
            "MW<=2 preferred band",
            "MW rank",
            "reducible-fibre support count",
            "semistable compatibility",
            "expected fibre-stratum dimension",
            "minimum nonzero-section pole order",
            "known marking before unknown marking, then Galois orbit size",
            "expected coefficient conditions",
            "certified neighbour cost, unknown last",
            "audited low-degree multisection richness as a final heuristic tie-break",
        ],
        "proof_boundary": {
            "proved": (
                "Root/MW data, support counts, A-type compatibility, and displayed "
                "minimum pole orders are exact consequences of the source lattices."
            ),
            "not_proved": (
                "Expected dimensions and coefficient conditions are deformation-count "
                "heuristics until an equation ansatz is checked. Unknown rational "
                "markings and routes are not imputed from lattice data. The ranking is "
                "not a specialization-rank prediction; sampled degree-three/four "
                "coordinates are not complete censuses."
            ),
        },
        "inputs": inputs,
        "target_mw_range": [arguments.target_mw_min, arguments.target_mw_max],
        "candidates": rows,
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/score_lattice_foundry_sources.sage"
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("source-ranking artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        f"FOUNDRYSOURCESCORE|candidates={len(rows)}|"
        f"mw_le_2={sum(row['source_mw_rank'] <= 2 for row in rows)}|"
        f"best={rows[0]['ns_id']}:{rows[0]['source_root_type']}|status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
