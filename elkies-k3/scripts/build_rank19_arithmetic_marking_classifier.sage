#!/usr/bin/env sage-python
"""Build the fail-closed arithmetic-marking classifier for rank-19 NS lattices.

The candidate set is the exact subset of the rank-seven catalogue carrying an
even rootless rank-17 frame.  For every candidate this script checks the
``NS = U + W(-1)`` lattice model against the catalogue's rank-three primitive
K3 complement, imports the replayed even-Clifford algebra/order, and keeps the
coarse norm-one curve separate from the full discriminant-kernel marking
curve.  Exact positive or negative arithmetic decisions are accepted only
from the small hash-pinned decision registry.

Unknown curve identifications and rational-point problems stay UNKNOWN.  The
equation-agent handoff contains only exact ARITHMETICALLY_POSSIBLE rows.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

from jsonschema import Draft202012Validator
from sage.all import QuadraticForm, ZZ, block_diagonal_matrix, factor, matrix


ROOT = Path(__file__).resolve().parents[2]
GENERATED = ROOT / "artifacts/generated-results"
CATALOGUE = GENERATED / "elkies-k3-rank7-auxiliary-catalogue-v1.json"
T_ARITHMETIC = GENERATED / "elkies-k3-rank7-t-arithmetic-v1.json"
DECISIONS = (
    ROOT
    / "elkies-k3/data/arithmetic/rank19-arithmetic-marking-decisions-v1.json"
)
SCHEMA = ROOT / "schemas/elkies_k3_rank19_arithmetic_marking_classifier.schema.json"
OUTPUT = GENERATED / "elkies-k3-rank19-arithmetic-marking-classifier-v1.json"
DISPATCH = GENERATED / "elkies-k3-rank19-arithmetic-marking-equation-survivors-v1.json"

CLASSIFICATIONS = {
    "ARITHMETICALLY_EXCLUDED",
    "ARITHMETICALLY_POSSIBLE",
    "UNKNOWN",
}
H3_SURFACE_ID = "K3-8188cdcda8c57b2d"


def relative(path):
    return str(path.resolve().relative_to(ROOT))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def nontrivial_elementary_divisors(value):
    return [abs(int(entry)) for entry in value.elementary_divisors() if abs(entry) > 1]


def factorization(value):
    return [[int(prime), int(exponent)] for prime, exponent in factor(ZZ(value))]


def resolve_path(payload, dotted_path):
    value = payload
    for component in dotted_path.split("."):
        if not isinstance(value, dict) or component not in value:
            raise KeyError(f"missing assertion path {dotted_path}")
        value = value[component]
    return value


def exact_rootless_rank17_frames(surface):
    determinant = int(surface["determinant"])
    witnesses = []
    for frame in surface["frames"]:
        if (
            int(frame["mw_rank_for_rho_19"]) != 17
            or frame["root_type"] != "0"
            or int(frame["root_rank"]) != 0
            or int(frame.get("signed_root_count", -1)) != 0
        ):
            continue
        gram = matrix(ZZ, frame["gram"])
        intrinsics = frame.get("rootless_intrinsics") or {}
        if not (
            gram.nrows() == gram.ncols() == 17
            and gram.is_symmetric()
            and all(gram[index, index] % 2 == 0 for index in range(17))
            and gram.det() == determinant
            and int(intrinsics.get("minimum_squared_norm", 0)) >= 4
        ):
            continue
        witnesses.append((frame, gram))
    return witnesses


def validate_certificate(certificate):
    path = ROOT / certificate["path"]
    if not path.is_file():
        raise FileNotFoundError(certificate["path"])
    payload = json.loads(path.read_text())
    for assertion in certificate["assertions"]:
        actual = resolve_path(payload, assertion["path"])
        if actual != assertion["equals"]:
            raise AssertionError(
                f"{certificate['path']}:{assertion['path']} changed: "
                f"{actual!r} != {assertion['equals']!r}"
            )
    return {
        "path": certificate["path"],
        "sha256": digest(path),
        "assertions_replayed": len(certificate["assertions"]),
    }


def default_full_marking_curve(t_row, discriminant_invariants):
    normalization = t_row["similarity_normalization"]
    literal = matrix(ZZ, t_row["literal_transcendental_gram"])
    hessian = matrix(ZZ, normalization["primitive_integral_quadratic_hessian"])
    scale = Fraction(
        int(normalization["literal_content"]),
        int(normalization["quadratic_integrality_scale"]),
    )
    return {
        "status": "UNKNOWN_FULL_DISCRIMINANT_KERNEL_CURVE_NOT_IDENTIFIED",
        "abstract_group": "O^+(T)^* = kernel(O^+(T) -> O(A_T))",
        "explicit_model": None,
        "coarse_norm_one_curve": t_row["arithmetic_source"]["base_curve"],
        "marking_gap": {
            "status": "OPEN_STABLE_ORTHOGONAL_DISCRIMINANT_KERNEL_NOT_COMPUTED",
            "discriminant_group_invariants": discriminant_invariants,
            "coarse_clifford_hessian_equals_literal_T": hessian == literal,
            "literal_T_as_rational_multiple_of_coarse_hessian": str(scale),
            "warning": (
                "The projective norm-one curve of a similar primitive ternary "
                "form is only a coarse orthogonal target. Full rational NS "
                "marking requires the discriminant-kernel subgroup for the "
                "literal integral T lattice."
            ),
        },
    }


def default_quotient_maps(t_row):
    base = t_row["arithmetic_source"]["base_curve"]
    if base["status"].startswith("PASS_EXACT"):
        label = base.get("label") or base.get("group")
        return [
            {
                "status": "CANDIDATE_FORGETFUL_MAP_MARKING_SUBGROUP_NOT_COMPUTED",
                "source": "full discriminant-marking curve",
                "target": label,
                "target_genus": base.get("genus"),
                "boundary": (
                    "The target curve is exact, but the map from the stable "
                    "marked subgroup has not been certified for this row."
                ),
            }
        ]
    return []


def default_arithmetic_tests(t_row):
    tests = []
    isotropy = t_row["rational_isotropy"]
    if not isotropy["isotropic"]:
        tests.append(
            {
                "status": "DIAGNOSTIC_ONLY_NOT_A_MARKING_OBSTRUCTION",
                "method": "ternary rational-isotropy local test",
                "obstruction_prime": isotropy["pari_qfsolve_obstruction_prime"],
                "boundary": (
                    "Anisotropy of T identifies a division quaternion algebra; "
                    "it does not exclude rational points on the marked Shimura curve."
                ),
            }
        )
    base = t_row["arithmetic_source"]["base_curve"]
    if base.get("genus") == 0:
        tests.append(
            {
                "status": "DIAGNOSTIC_ONLY_COARSE_GENUS_ZERO",
                "method": "coarse norm-one curve genus",
                "boundary": (
                    "A rational coarse curve does not prove that its full "
                    "discriminant-marking cover has a rational noncuspidal point."
                ),
            }
        )
    return tests


def unknown_next_gate(t_row):
    source = t_row["arithmetic_source"]
    if source["status"] == "PASS_EXACT_SPLIT_EICHLER_MODULAR_CURVE":
        return (
            "Compute the stable discriminant-kernel subgroup inside the exact "
            "Gamma_0(N) norm-one group, then identify its easiest rational-point quotient."
        )
    if source["status"] == "PARTIAL_EXACT_GENERAL_SPLIT_CLIFFORD_ORDER":
        return (
            "Embed the displayed split order in M_2(QQ), derive its exact congruence "
            "conditions and signature, and then compute the stable marking subgroup."
        )
    if source["status"] == "PARTIAL_EXACT_ANISOTROPIC_CLIFFORD_ORDER":
        return (
            "Certify the local order type and normalizer/Atkin-Lehner quotient before "
            "attempting rational-point or CM tests."
        )
    return "Identify the exact full discriminant-marking period curve."


def build_candidate(surface, t_row, decision):
    witnesses = exact_rootless_rank17_frames(surface)
    if not witnesses:
        raise AssertionError("build_candidate called without a rootless MW17 witness")

    determinant = int(surface["determinant"])
    transcendental = matrix(ZZ, surface["surface_key"]["transcendental_gram"])
    if rows(transcendental) != t_row["literal_transcendental_gram"]:
        raise AssertionError(f"T ledger mismatch for {surface['surface_id']}")
    if transcendental.nrows() != 3 or transcendental.det() != -determinant:
        raise AssertionError(f"bad ternary determinant for {surface['surface_id']}")
    if QuadraticForm(ZZ, transcendental).signature_vector() != (2, 1, 0):
        raise AssertionError(f"bad T signature for {surface['surface_id']}")

    discriminant_invariants = [
        int(value)
        for value in surface["surface_key"]["ns_discriminant_form_key"]["invariants"]
        if int(value) > 1
    ]
    frame_records = []
    for frame, frame_gram in witnesses:
        ns_gram = block_diagonal_matrix(
            matrix(ZZ, [[0, 1], [1, 0]]), -frame_gram
        )
        if QuadraticForm(ZZ, ns_gram).signature_vector() != (1, 18, 0):
            raise AssertionError(f"bad NS signature for {frame['frame_id']}")
        if ns_gram.det() != determinant:
            raise AssertionError(f"bad NS determinant for {frame['frame_id']}")
        if nontrivial_elementary_divisors(ns_gram) != nontrivial_elementary_divisors(
            transcendental
        ):
            raise AssertionError(f"NS/T Smith mismatch for {frame['frame_id']}")
        frame_records.append(
            {
                "frame_id": frame["frame_id"],
                "frame_gram_sha256": frame["gram_sha256"],
                "construction": "NS = U orthogonal_sum W(-1)",
                "rank": 19,
                "signature": [1, 18],
                "determinant": determinant,
                "rootless_mw_rank": 17,
                "minimum_squared_norm": int(
                    frame["rootless_intrinsics"]["minimum_squared_norm"]
                ),
            }
        )

    classification = "UNKNOWN" if decision is None else decision["classification"]
    if classification not in CLASSIFICATIONS:
        raise ValueError(f"bad classification for {surface['surface_id']}")
    full_curve = default_full_marking_curve(t_row, discriminant_invariants)
    quotient_maps = default_quotient_maps(t_row)
    arithmetic_tests = default_arithmetic_tests(t_row)
    certificate_replay = []
    theorem_inputs = []
    decision_text = (
        "No exact positive witness or global exclusion certificate is registered."
    )
    if decision is not None:
        certificate_replay = [
            validate_certificate(certificate)
            for certificate in decision["certificates"]
        ]
        full_curve.update(decision["full_marking_curve"])
        # Retain the automatically exposed coarse curve and marking-gap warning.
        full_curve["coarse_norm_one_curve"] = t_row["arithmetic_source"]["base_curve"]
        quotient_maps = decision["easy_quotient_maps"]
        arithmetic_tests = decision["arithmetic_tests"]
        theorem_inputs = decision["theorem_inputs"]
        decision_text = decision["decision"]

    normalization = t_row["similarity_normalization"]
    row = {
        "surface_id": surface["surface_id"],
        "legacy_ns_ids": surface["legacy_ns_ids"],
        "determinant": determinant,
        "ns_lattice": {
            "status": "PASS_EXACT_CATALOGUE_NS_FROM_ROOTLESS_FRAME",
            "rank": 19,
            "signature": [1, 18],
            "discriminant_group_invariants": discriminant_invariants,
            "frame_witnesses": frame_records,
            "catalogue_boundary": (
                "The catalogue constructs the primitive K3 complement by the "
                "opposite discriminant form; this replay checks rank, signature, "
                "determinant, and Smith invariants."
            ),
        },
        "transcendental_lattice": {
            "status": "PASS_EXACT_PRIMITIVE_K3_COMPLEMENT_FROM_CATALOGUE",
            "construction": "T = NS^perp in the K3 lattice",
            "gram": rows(transcendental),
            "rank": 3,
            "signature": [2, 1],
            "determinant": int(transcendental.det()),
            "smith_invariants_greater_than_one": nontrivial_elementary_divisors(
                transcendental
            ),
            "rational_isotropy": t_row["rational_isotropy"],
        },
        "even_clifford_algebra_and_order": {
            "status": "PASS_EXACT_IMPORTED_T_ARITHMETIC_REPLAY",
            **t_row["clifford"],
            "similarity_normalization": normalization,
        },
        "full_discriminant_marking_curve": full_curve,
        "easy_quotient_maps": quotient_maps,
        "arithmetic_tests": arithmetic_tests,
        "classification": classification,
        "classification_decision": decision_text,
        "theorem_inputs": theorem_inputs,
        "certificate_replay": certificate_replay,
        "equation_agent_eligible": classification == "ARITHMETICALLY_POSSIBLE",
        "different_ns_foundry_equation_eligible": (
            classification == "ARITHMETICALLY_POSSIBLE"
            and surface["surface_id"] != H3_SURFACE_ID
        ),
        "next_arithmetic_gate": (
            None if classification != "UNKNOWN" else unknown_next_gate(t_row)
        ),
    }
    return row


def hospitality_comparison(by_id):
    selected = [
        H3_SURFACE_ID,
        "K3-ebaf00b3723751ba",
        "K3-d1b1381f87d69f1c",
        "K3-f43753fb154e3406",
    ]
    rows_out = []
    for surface_id in selected:
        row = by_id[surface_id]
        clifford = row["even_clifford_algebra_and_order"]
        order = clifford["integral_even_clifford_order"]
        rows_out.append(
            {
                "surface_id": surface_id,
                "determinant": row["determinant"],
                "classification": row["classification"],
                "discriminant_group_invariants": row["ns_lattice"][
                    "discriminant_group_invariants"
                ],
                "quaternion_discriminant": clifford["quaternion_discriminant"],
                "coarse_order_reduced_discriminant": order["reduced_discriminant"],
                "coarse_local_level_index": order["local_level_index"],
                "determinant_factorization": factorization(row["determinant"]),
                "arithmetic_reason": row["classification_decision"],
            }
        )
    return {
        "status": "EVIDENCE_BACKED_STRUCTURAL_COMPARISON_NOT_A_CLASSIFICATION_THEOREM",
        "rows": rows_out,
        "determinant_948_features": [
            "T has a division even-Clifford algebra with the clean coprime squarefree Eichler pair (D,N)=(6,79).",
            "The Atkin-Lehner quotient by w_474 has genus 2 and an explicit non-CM QQ-point.",
            "The rootless QQ equation supplies seventeen rational sections, rho=19, saturation, and trivial torsion, so the point really lifts to the full marking.",
        ],
        "contrast": (
            "Determinant 720 has exact stable curve X_0(60), whose rational points "
            "are all cusps; determinant 950 is forced onto the rigid Fricke quotient "
            "X_0^+(475); and determinant 1184 combines non-split Cartan level 4 with "
            "X_0(37), where both rational 37-isogeny points fail the Frobenius lift. The "
            "observed hospitality of 948 is therefore explained by its unusually "
            "low-genus Atkin-Lehner quotient plus an actual non-CM rational lift, "
            "not by determinant size or Clifford splitting alone."
        ),
    }


def build(catalogue, t_arithmetic, decisions, paths):
    if catalogue["schema"] != "elkies-k3.rank7-auxiliary-catalogue.v1":
        raise ValueError("unexpected catalogue schema")
    if t_arithmetic["schema"] != "elkies-k3.rank7-t-arithmetic.v1":
        raise ValueError("unexpected T-arithmetic schema")
    if decisions["schema"] != "elkies-k3.rank19-arithmetic-marking-decisions.v1":
        raise ValueError("unexpected decision-registry schema")
    if t_arithmetic["input"]["catalogue_sha256"] != digest(paths["catalogue"]):
        raise ValueError("T-arithmetic ledger does not match the catalogue")

    t_by_id = {row["surface_id"]: row for row in t_arithmetic["surfaces"]}
    if set(t_by_id) != {row["surface_id"] for row in catalogue["surfaces"]}:
        raise ValueError("catalogue/T-arithmetic surface sets differ")
    decision_by_id = {row["surface_id"]: row for row in decisions["records"]}
    if len(decision_by_id) != len(decisions["records"]):
        raise ValueError("duplicate decision-registry surface")

    candidate_surfaces = [
        surface
        for surface in catalogue["surfaces"]
        if exact_rootless_rank17_frames(surface)
    ]
    candidate_ids = {surface["surface_id"] for surface in candidate_surfaces}
    if not set(decision_by_id) <= candidate_ids:
        raise ValueError("decision registry contains a noncandidate surface")
    if len(candidate_surfaces) != 66:
        raise AssertionError(f"expected 66 candidate NS lattices, got {len(candidate_surfaces)}")

    candidates = [
        build_candidate(
            surface,
            t_by_id[surface["surface_id"]],
            decision_by_id.get(surface["surface_id"]),
        )
        for surface in candidate_surfaces
    ]
    candidates.sort(key=lambda row: (row["determinant"], row["surface_id"]))
    counts = Counter(row["classification"] for row in candidates)
    if counts != Counter(
        {
            "ARITHMETICALLY_EXCLUDED": 3,
            "ARITHMETICALLY_POSSIBLE": 1,
            "UNKNOWN": 62,
        }
    ):
        raise AssertionError(f"classification count changed: {counts}")

    possible = [
        row["surface_id"]
        for row in candidates
        if row["equation_agent_eligible"]
    ]
    different_ns = [
        row["surface_id"]
        for row in candidates
        if row["different_ns_foundry_equation_eligible"]
    ]
    by_id = {row["surface_id"]: row for row in candidates}
    return {
        "schema": "elkies-k3.rank19-arithmetic-marking-classifier.v1",
        "status": "PASS_FAIL_CLOSED_1_POSSIBLE_3_EXCLUDED_62_UNKNOWN",
        "policy": {
            **decisions["policy"],
            "equation_agent": (
                "Dispatch only ARITHMETICALLY_POSSIBLE rows. UNKNOWN is an "
                "arithmetic-research queue, not an equation queue."
            ),
            "coarse_curve_boundary": (
                "A Clifford norm-one curve, Gamma_0(N) label, genus-zero model, "
                "formal local branch, or bounded point search does not identify "
                "the full discriminant marking and cannot promote UNKNOWN."
            ),
        },
        "proof_scope": {
            "proved": (
                "All 66 exact rootless-MW17 candidate NS lattices are paired with "
                "their catalogue primitive ternary complement and replayed even "
                "Clifford data. The four non-UNKNOWN decisions are backed by exact "
                "registered certificates."
            ),
            "not_proved": (
                "The 62 UNKNOWN rows are not asserted to exist or not exist over QQ. "
                "Their stable discriminant-kernel curves and rational points remain open."
            ),
        },
        "inputs": {relative(path): digest(path) for path in paths.values()},
        "accounting": {
            "catalogue_surfaces": len(catalogue["surfaces"]),
            "candidate_ns_lattices": len(candidates),
            **{classification: counts[classification] for classification in sorted(CLASSIFICATIONS)},
            "equation_agent_survivors": len(possible),
            "different_ns_equation_agent_survivors": len(different_ns),
            "coarse_curve_statuses": dict(
                sorted(
                    Counter(
                        t_by_id[row["surface_id"]]["arithmetic_source"]["status"]
                        for row in candidates
                    ).items()
                )
            ),
        },
        "equation_agent_dispatch": {
            "policy": "ARITHMETICALLY_POSSIBLE only",
            "all_arithmetically_possible_surface_ids": possible,
            "different_ns_objective_surface_ids": different_ns,
            "already_realized_controls_not_requeued": [H3_SURFACE_ID],
            "new_equation_work_queue": different_ns,
        },
        "hospitality_comparison": hospitality_comparison(by_id),
        "candidates": candidates,
        "reproduce": (
            "sage -python elkies-k3/scripts/build_rank19_arithmetic_marking_classifier.sage"
        ),
    }


def dispatch_payload(classifier, classifier_path):
    eligible = set(classifier["equation_agent_dispatch"]["new_equation_work_queue"])
    rows_out = [
        {
            "surface_id": row["surface_id"],
            "legacy_ns_ids": row["legacy_ns_ids"],
            "determinant": row["determinant"],
            "classification": row["classification"],
            "certificate_replay": row["certificate_replay"],
        }
        for row in classifier["candidates"]
        if row["surface_id"] in eligible
    ]
    return {
        "schema": "elkies-k3.rank19-arithmetic-marking-equation-survivors.v1",
        "status": "PASS_EMPTY_NEW_DIFFERENT_NS_EQUATION_QUEUE",
        "input": {
            "path": relative(classifier_path),
            "sha256": hashlib.sha256(
                (json.dumps(classifier, indent=2, sort_keys=True) + "\n").encode()
            ).hexdigest(),
        },
        "policy": (
            "This handoff contains only new different-NS rows with an exact "
            "ARITHMETICALLY_POSSIBLE classification. UNKNOWN rows are intentionally absent."
        ),
        "survivors": rows_out,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalogue", type=Path, default=CATALOGUE)
    parser.add_argument("--t-arithmetic", type=Path, default=T_ARITHMETIC)
    parser.add_argument("--decisions", type=Path, default=DECISIONS)
    parser.add_argument("--schema", type=Path, default=SCHEMA)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--dispatch", type=Path, default=DISPATCH)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    paths = {
        "catalogue": arguments.catalogue.resolve(),
        "t_arithmetic": arguments.t_arithmetic.resolve(),
        "decisions": arguments.decisions.resolve(),
        "schema": arguments.schema.resolve(),
    }
    payloads = {
        name: json.loads(path.read_text())
        for name, path in paths.items()
        if name != "schema"
    }
    result = build(
        payloads["catalogue"],
        payloads["t_arithmetic"],
        payloads["decisions"],
        paths,
    )
    Draft202012Validator(json.loads(paths["schema"].read_text())).validate(result)
    encoded = json.dumps(result, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    dispatch = dispatch_payload(result, output)
    dispatch_encoded = json.dumps(dispatch, indent=2, sort_keys=True) + "\n"
    dispatch_output = arguments.dispatch.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit(f"stale classifier artifact: {output}")
        if not dispatch_output.exists() or dispatch_output.read_text() != dispatch_encoded:
            raise SystemExit(f"stale equation-survivor artifact: {dispatch_output}")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
        dispatch_output.parent.mkdir(parents=True, exist_ok=True)
        dispatch_output.write_text(dispatch_encoded)
    accounting = result["accounting"]
    print(
        "ARITHMARK|candidates={}|possible={}|excluded={}|unknown={}|"
        "new_equation_survivors={}|status=PASS".format(
            accounting["candidate_ns_lattices"],
            accounting["ARITHMETICALLY_POSSIBLE"],
            accounting["ARITHMETICALLY_EXCLUDED"],
            accounting["UNKNOWN"],
            accounting["different_ns_equation_agent_survivors"],
        )
    )


if __name__ == "__main__":
    main()
