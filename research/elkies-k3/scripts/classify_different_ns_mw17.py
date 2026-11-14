#!/usr/bin/env python3
"""Join retained K3 evidence into a bounded MW17 classification experiment.

No lattice enumeration, curve search, proof replay or equation compilation.
Input bytes and supporting claim scopes are pinned in the reviewed manifest.
Missing inputs, stale claims and contradictory evidence stop the build.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "elkies-k3/data/arithmetic/mw17-classification-inputs-v1.json"
OUT = ROOT / "artifacts/generated-results/elkies-k3-mw17-classification-v1"
NOTE = ROOT / "elkies-k3/DIFFERENT_NS_CLASSIFICATION_2026-09-15.md"
GATES = ("geometric_admissibility", "galois_compatible_marking", "rational_NS",
         "rational_elliptic_fibration", "arithmetic_MW17")
OBSTRUCTION_ORDER = ("REAL_CONTENT", "REAL_PERIOD", "LOCAL_SQUARE_TWO",
                     "SCALED_SPLIT_MARKED_CURVE", "SPECIFIC_MARKED_CURVE")
CLAIMS = {
    "catalogue": "EC-K3-RANK19-ARITHMETIC-MARKING-CLASSIFIER",
    "content": "EC-K3-REAL-MARKING-CONTENT-OBSTRUCTION-20260914",
    "period": "EC-K3-REAL-RANK19-PERIOD-OBSTRUCTION-20260914",
    "structural": "EC-K3-RATIONAL-MARKING-STRUCTURAL-GATES-20260914",
    "det800": "EC-K3-DET800-FULL-MARKING-OBSTRUCTION-20260914",
    "948": "EC-K3-R17-ALTERNATE-Q80-ARITHMETIC-RANK17",
    "388": "EC-K3-DET388-RATIONAL-MARKING-SOURCE-20260914",
    "852": "EC-K3-DET852-RATIONAL-MARKING-SOURCE-20260915",
    "622": "EC-K3-DET622-INOSE-RATIONAL-MARKING-SOURCE-20260915",
    "1020": "EC-K3-DET1020-RATIONAL-MARKING-SOURCE-20260915",
    "root388": "EC-K3-DET388-GLOBAL-ROOT-OBSTRUCTION-20260914",
    "root852": "EC-K3-DET852-GLOBAL-ROOT-OBSTRUCTION-20260915",
    "root622": "EC-K3-DET622-GLOBAL-ROOT-OBSTRUCTION-20260915",
    "root1020": "EC-K3-DET1020-ARITHMETIC-MW17-EXISTENCE-20260915",
    "physical1020": "EC-K3-DET1020-PHYSICAL-NEF-MW17-20260915",
    "1236": "EC-K3-DET1236-CANDIDATE-DOUBLE-COVER",
}
SPECIFIC = {
    "K3-04b86146cc6b284b": "EC-K3-DET500-DET750-QQ-MARKING-OBSTRUCTIONS",
    "K3-10a14a46c14b3150": "EC-K3-DET500-DET750-QQ-MARKING-OBSTRUCTIONS",
    "K3-7b71a1cc00b0c6e2": "EC-K3-DET378-QQ-MARKING-OBSTRUCTION",
    "K3-d1b1381f87d69f1c": "EC-K3-NS0031-QQ-MARKING-OBSTRUCTION",
    "K3-ebaf00b3723751ba": "EC-K3-NS0024-QQ-MARKING-OBSTRUCTION",
    "K3-f43753fb154e3406": "EC-K3-GOLAY-DET720-QQ-MARKING-OBSTRUCTION",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fingerprint(entry):
    return hashlib.sha256("\0".join(entry[k] for k in ("state", "title", "scope"))
                          .encode()).hexdigest()[:16]


def unique(rows, key):
    result = {row[key]: row for row in rows}
    require(len(result) == len(rows), f"duplicate {key}")
    return result


def integral_gram(rows):
    values = [[Fraction(x) for x in row] for row in rows]
    require(all(x.denominator == 1 for row in values for x in row), "nonintegral Gram")
    return [[int(x) for x in row] for row in values]


def read_inputs(manifest_path=MANIFEST):
    manifest = json.loads(manifest_path.read_text())
    require(manifest["schema"] == "elkies-k3.mw17-classification-inputs.v1", "manifest schema")
    data = {}
    for name, receipt in manifest["inputs"].items():
        path = ROOT / receipt["path"]
        require(path.is_file(), f"missing retained input: {receipt['path']}; no reconstruction")
        require(sha(path) == receipt["sha256"], f"input drift: {receipt['path']}")
        if path.suffix == ".json":
            data[name] = json.loads(path.read_text())
    entries = unique(json.loads((ROOT / "MATH_STATUS.json").read_text())["entries"], "id")
    needed = set(CLAIMS.values()) | set(SPECIFIC.values())
    require(set(manifest["claims"]) == needed, "claim coverage differs")
    for identifier, receipt in manifest["claims"].items():
        entry = entries[identifier]
        require(fingerprint(entry) == receipt["scope_fingerprint"], f"claim needs review: {identifier}")
        require(entry["state"] == receipt["state"], f"claim state: {identifier}")
        require(entry["canonical_source"] == receipt["canonical_source"], f"claim source: {identifier}")
        require((ROOT / entry["canonical_source"]).is_file(), f"missing canonical source: {identifier}")
        require(sha(ROOT / entry["canonical_source"]) == receipt["canonical_source_sha256"],
                f"canonical source needs review: {identifier}")
        if identifier not in {CLAIMS["catalogue"], CLAIMS["1236"]}:
            require(entry["state"] == "proved", f"unproved decision: {identifier}")
    return manifest, data


def cell(state, evidence=(), reason=None):
    result = {"state": state, "evidence": list(evidence)}
    if reason:
        result["reason"] = reason
    return result


def finish(row):
    failures = [g for g in GATES if row["gates"][g]["state"] == "FAIL"]
    unknowns = [g for g in GATES if row["gates"][g]["state"] == "UNKNOWN"]
    row["first_certified_failure"] = failures[0] if failures else None
    row["first_unresolved_gate"] = unknowns[0] if unknowns and not failures else None
    row["classification"] = "EXCLUDED" if failures else "UNKNOWN" if unknowns else "PASS"
    if row["classification"] == "PASS":
        require(all(row["gates"][g]["state"] == "PASS" for g in GATES), "incomplete pass")
    return row


def counts(rows, field="classification"):
    c = Counter(r[field] for r in rows)
    return {state: c[state] for state in ("PASS", "EXCLUDED", "UNKNOWN")}


def build(manifest, data):
    catalogue = data["catalogue"]["surfaces"]
    surfaces = unique(catalogue, "surface_id")
    arithmetic = unique(data["arithmetic"]["surfaces"], "surface_id")
    require(set(surfaces) == set(arithmetic), "T/NS population mismatch")
    require(len(surfaces) == manifest["population"]["NS_classes"], "NS population changed")
    planner = data["planner"]
    specific = unique(planner["excluded_before_NS_or_equation_work"], "surface_id")
    require(set(specific) == set(SPECIFIC), "specific exclusions changed")
    require(planner["positive_controls"] == [{"already_realized": True, "determinant": 948,
             "surface_id": "K3-8188cdcda8c57b2d"}], "positive control identity changed")
    exclusions = {sid: [] for sid in surfaces}

    def exclude(record, code, evidence, detail):
        sid = record["surface_id"]
        require(sid in surfaces, f"unknown excluded row {sid}")
        require(record["determinant"] == surfaces[sid]["determinant"], f"determinant drift {sid}")
        if "gram" in record:
            require(record["gram"] == arithmetic[sid]["literal_transcendental_gram"], f"Gram drift {sid}")
        exclusions[sid].append({"code": code, "evidence": evidence, "detail": detail})

    content = unique(data["content"]["rows"], "surface_id")
    require(set(content) == set(surfaces), "content coverage")
    for r in content.values():
        if r["full_real_marking_excluded"]:
            require(r["content"] > 2, "invalid content exclusion")
            exclude(r, "REAL_CONTENT", CLAIMS["content"], {"content": r["content"]})
    for r in unique(data["period"]["rows"], "surface_id").values():
        require(r["hilbert_symbol"] == -1 and r["all_integral_squares_divisible_by_four"], "period evidence")
        exclude(r, "REAL_PERIOD", CLAIMS["period"], {"anisotropy_prime": r["obstruction_prime"]})
    for r in unique(data["structural"]["local_square2_rows"], "surface_id").values():
        require(r["hilbert_symbol"] == -1 and r["residues_exhausted"] == r["obstructing_modulus"] ** 3,
                "local square-two evidence")
        exclude(r, "LOCAL_SQUARE_TWO", CLAIMS["structural"],
                {"modulus": r["obstructing_modulus"], "anisotropy_prime": r["anisotropy_prime"]})
    for r in unique(data["structural"]["scaled_split_rows"], "surface_id").values():
        require(r["N"] >= 5, "scaled-split theorem range")
        exclude(r, "SCALED_SPLIT_MARKED_CURVE", CLAIMS["structural"], {"level": 4 * r["N"]})
    for sid, r in specific.items():
        require(r["full_marking_curve"]["rational_non_CM_point"] is False, "specific exclusion evidence")
        exclude(r, "SPECIFIC_MARKED_CURVE", SPECIFIC[sid], {"reason": r["decision"]})
    r = data["det800"]
    require(r["rational_non_CM_marking_exists"] is False, "det800 evidence")
    # Keep its independent proof, while assigning one primary cause per row.
    exclude({"surface_id": r["surface_id"], "determinant": 800},
            "SPECIFIC_MARKED_CURVE", CLAIMS["det800"], {"curve": "X_0(200)"})

    admitted = {"K3-8188cdcda8c57b2d": CLAIMS["948"]}
    globally_rootful = {}
    for det in ("388", "852"):
        admission, root = data[f"admit{det}"], data[f"root{det}"]
        require(admission["full_rational_NS_existence"] ==
                ("PROVED_BY_WRITTEN_DESCENT_ARGUMENT" if det == "388" else
                 "PROVED_BY_WRITTEN_MODULI_AND_DIVISOR_DESCENT"), "admission evidence")
        sid = admission["surface_id"]
        require(root.get("surface_id", root.get("source_surface_id")) == sid, "global root identity")
        require(root["status"] == "PASS_ALL_FRAMES_ROOTFUL", "global root evidence")
        require(integral_gram(admission["T_gram"]) == arithmetic[sid]["literal_transcendental_gram"], "admission T identity")
        admitted[sid] = CLAIMS[det]
        globally_rootful[sid] = CLAIMS[f"root{det}"]
    require(not set(admitted) & {sid for sid, es in exclusions.items() if es}, "admitted/excluded contradiction")

    rows, frames = [], []
    for sid, s in sorted(surfaces.items()):
        t = arithmetic[sid]
        require(s["surface_key"]["transcendental_gram"] == t["literal_transcendental_gram"], "literal T mismatch")
        require(s["determinant"] == t["determinant"], "T determinant mismatch")
        rootless = []
        for f in s["frames"]:
            require(f["root_rank"] + f["mw_rank_for_rho_19"] == 17, "rank balance")
            require((f["signed_root_count"] == 0) == (f["root_rank"] == 0), "root data contradiction")
            require(f["determinant"] == s["determinant"] and len(f["gram"]) == 17, "frame identity")
            if f["root_rank"] == 0:
                rootless.append(f["frame_id"])
        require(not (rootless and sid in globally_rootful), "rootless/global root contradiction")
        obs = sorted(exclusions[sid], key=lambda e: (OBSTRUCTION_ORDER.index(e["code"]), e["evidence"]))
        g = {gate: cell("UNKNOWN") for gate in GATES}
        if rootless:
            g[GATES[0]] = cell("PASS", ["input:catalogue", CLAIMS["catalogue"]])
        elif sid in globally_rootful:
            g[GATES[0]] = cell("FAIL", [globally_rootful[sid]])
        else:
            g[GATES[0]] = cell("UNKNOWN", reason="No retained rootless witness or all-frames exclusion.")
        if obs:
            g[GATES[1]] = cell("FAIL", sorted({e["evidence"] for e in obs}))
            g[GATES[2]] = cell("NOT_REACHED", reason="Full rational marking is excluded.")
        elif sid in admitted:
            g[GATES[1]] = cell("PASS", [admitted[sid]])
            g[GATES[2]] = cell("PASS", [admitted[sid]])
        else:
            reason = ("Identify the characteristic-zero cubic CM branch orbit; the candidate cover is not the marked curve."
                      if sid == "K3-6d288cfad55e0d15" else
                      "Full canonical marked curve and rational non-CM period are not certified by these inputs.")
            g[GATES[1]] = cell("UNKNOWN", [CLAIMS["1236"]] if sid == "K3-6d288cfad55e0d15" else (), reason)
            g[GATES[2]] = cell("UNKNOWN", reason="Effective K3 descent and full saturated rational divisor descent are unproved.")
        if obs or sid in globally_rootful:
            for gate in GATES[3:]:
                g[gate] = cell("NOT_REACHED", reason="An exact necessary gate failed; see first_certified_failure.")
        elif sid in admitted and rootless:
            for gate in GATES[3:]:
                g[gate] = cell("PASS", [admitted[sid]], "Rootless U, full rational divisors and the nef/fibration theorem.")
        row = finish({"surface_id": sid, "determinant": s["determinant"],
            "determinant_band": s["determinant_band"], "legacy_ns_ids": s["legacy_ns_ids"],
            "transcendental_gram": t["literal_transcendental_gram"],
            "identity_reference": "input:catalogue surfaces[surface_id].surface_key (T and full NS discriminant form)",
            "literal_content": content[sid]["content"],
            "rationally_isotropic_T": t["rational_isotropy"]["isotropic"],
            "quaternion_discriminant": t["clifford"]["quaternion_discriminant"],
            "retained_frame_ids": [f["frame_id"] for f in s["frames"]],
            "retained_rootless_frame_ids": rootless, "gates": g,
            "obstructions": obs, "primary_arithmetic_obstruction": obs[0]["code"] if obs else None,
            "arithmetic_admission": "PASS" if sid in admitted else "EXCLUDED" if obs else "UNKNOWN"})
        rows.append(row)
        for f in s["frames"]:
            frame_gates = dict(g)
            frame_gates[GATES[0]] = cell("PASS" if f["root_rank"] == 0 else "FAIL", ["input:catalogue"],
                                       "This fixed frame only; a rootful frame does not exclude its NS genus.")
            if f["root_rank"]:
                for gate in GATES[3:]:
                    frame_gates[gate] = cell("NOT_REACHED", reason="This frame has reducible-fibre roots.")
            frames.append(finish({"frame_id": f["frame_id"], "surface_id": sid,
                "gram_sha256": f["gram_sha256"], "root_rank": f["root_rank"],
                "signed_root_count": f["signed_root_count"], "geometric_MW_rank": f["mw_rank_for_rho_19"],
                "gates": frame_gates}))
    unique(frames, "frame_id")
    geometric = [r for r in rows if r["retained_rootless_frame_ids"]]
    geometric_frames = [f for f in frames if f["root_rank"] == 0]
    require((len(frames), len(geometric), len(geometric_frames)) ==
            tuple(manifest["population"][k] for k in ("frames", "rootless_NS_classes", "rootless_frames")),
            "frame/cohort coverage changed")
    # Reconcile the independent historical frontier without using its counts as decisions.
    unknown = {r["surface_id"] for r in rows if r["arithmetic_admission"] == "UNKNOWN"}
    require(unknown == set(data["structural"]["surviving_queue_ids"]) - set(admitted), "frontier disagreement")
    controls = []
    for det in ("622", "1020"):
        require(data[f"admit{det}"]["status"] ==
                ("PASS_RATIONAL_INOSE_MODEL_AND_EXISTENTIAL_FULL_MARKING_TWIST" if det == "622" else
                 "PASS_FULL_RATIONAL_SATURATED_NS_BY_INTRINSIC_NONCM_AND_DESCENT"), "external admission evidence")
        require(not any(s["surface_key"]["transcendental_gram"] == integral_gram(data[f"admit{det}"]["T_gram"])
                        for s in catalogue), "external literal control now in catalogue; review identity")
        root = data[f"root{det}"]
        passed = det == "1020"
        require(root["status"] == ("PASS_ROOTLESS_FRAME_AND_DISCRIMINANT_MATCH" if passed else
                                  "PASS_ALL_FRAMES_ROOTFUL"), "external root evidence")
        if passed:
            require(root["rank"] == 17 and root["norm2_vectors"] == 0, "external rootless evidence")
            require(data["physical1020"]["status"] == "PASS_PHYSICAL_NEF_MW17_CLASSES", "physical pencil evidence")
        controls.append({"control_id": f"OUTSIDE-CATALOGUE-DET{det}", "cohort": "external_control_only",
            "determinant": int(det), "T_gram": data[f"admit{det}"]["T_gram"],
            "arithmetic_admission": "PASS", "geometric_MW17": "PASS" if passed else "EXCLUDED",
            "arithmetic_MW17": "PASS" if passed else "EXCLUDED",
            "evidence": [CLAIMS[det], CLAIMS[f"root{det}"]],
            "equation_compilation": "MW17 equation and section coordinates uncomputed" if passed else
                                    "No MW17 fibration; rationalizing twist squareclass uncomputed"})
    primary = {}
    for code in OBSTRUCTION_ORDER:
        chosen = [r for r in rows if r["primary_arithmetic_obstruction"] == code]
        primary[code] = {"NS_classes": len(chosen),
                         "rootless_NS_classes": sum(bool(r["retained_rootless_frame_ids"]) for r in chosen),
                         "rootless_frames": sum(len(r["retained_rootless_frame_ids"]) for r in chosen)}
    summary = {"NS_classes": len(rows), "frames": len(frames),
        "arithmetic_admission": counts(rows, "arithmetic_admission"), "MW17_NS_outcomes": counts(rows),
        "rootless_NS_classes": len(geometric), "rootless_NS_outcomes": counts(geometric),
        "rootless_frames": len(geometric_frames), "rootless_frame_outcomes": counts(geometric_frames),
        "primary_arithmetic_obstructions": primary,
        "rootless_NS_gate_states": {g: dict(sorted(Counter(r["gates"][g]["state"] for r in geometric).items())) for g in GATES},
        "by_content": {str(c): counts([r for r in rows if r["literal_content"] == c], "arithmetic_admission")
                       for c in sorted({r["literal_content"] for r in rows})},
        "by_determinant_band": {b: counts([r for r in rows if r["determinant_band"] == b], "arithmetic_admission")
                                for b in sorted({r["determinant_band"] for r in rows})}}
    return {"schema": "elkies-k3.mw17-classification.v1", "authority": "MATH_STATUS.json",
        "boundary": "Bounded retained-evidence join. No exhaustive lattice census, new obstruction theorem or independent proof replay. UNKNOWN is not a failure.",
        "implementation_sha256": sha(Path(__file__)), "inputs": manifest["inputs"], "claims": manifest["claims"],
        "gate_order": list(GATES), "primary_obstruction_order": list(OBSTRUCTION_ORDER),
        "first_failure_semantics": "Earliest certified failure in this order; an earlier UNKNOWN may have been bypassed. Independent evidence is retained at later gates.",
        "summary": summary, "surfaces": rows, "frames": frames, "external_controls": controls}


def tsv(rows, identity):
    output = io.StringIO()
    fields = [identity, "surface_id"] if identity == "frame_id" else [identity, "determinant", "arithmetic_admission", "primary_arithmetic_obstruction"]
    fields += list(GATES) + ["classification", "first_certified_failure", "first_unresolved_gate"]
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        values = {k: row["gates"][k]["state"] if k in GATES else row.get(k) for k in fields}
        writer.writerow({k: "NONE" if value is None else value for k, value in values.items()})
    return output.getvalue()


def table(summary):
    lines = ["| Cohort | Total | Pass | Excluded | Unknown |", "|---|---:|---:|---:|---:|"]
    for label, total, key in (("Retained NS: full rational marking", "NS_classes", "arithmetic_admission"),
        ("Geometric MW17 NS classes: arithmetic MW17", "rootless_NS_classes", "rootless_NS_outcomes"),
        ("Geometric MW17 frames: arithmetic MW17", "rootless_frames", "rootless_frame_outcomes")):
        c = summary[key]
        lines.append(f"| {label} | {summary[total]} | {c['PASS']} | {c['EXCLUDED']} | {c['UNKNOWN']} |")
    lines += ["", "Primary arithmetic causes below are disjoint; all overlapping proofs remain in the JSON.", "",
              "| First arithmetic obstruction | All retained NS | Geometric MW17 NS | Geometric MW17 frames |",
              "|---|---:|---:|---:|"]
    for code, c in summary["primary_arithmetic_obstructions"].items():
        lines.append(f"| {code} | {c['NS_classes']} | {c['rootless_NS_classes']} | {c['rootless_frames']} |")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    manifest, data = read_inputs()
    result = build(manifest, data)
    outputs = {OUT / "classification.json": json.dumps(result, indent=2, sort_keys=True) + "\n",
               OUT / "surfaces.tsv": tsv(result["surfaces"], "surface_id"),
               OUT / "frames.tsv": tsv(result["frames"], "frame_id")}
    note = NOTE.read_text()
    start, stop = "<!-- BEGIN CLASSIFICATION COUNTS -->", "<!-- END CLASSIFICATION COUNTS -->"
    require(note.count(start) == note.count(stop) == 1, "report markers")
    before, rest = note.split(start)
    _, after = rest.split(stop)
    outputs[NOTE] = before + start + "\n\n" + table(result["summary"]) + "\n\n" + stop + after
    for path, expected in outputs.items():
        if args.check:
            require(path.is_file() and path.read_text() == expected, f"classification drift: {path.relative_to(ROOT)}")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(expected)
    print(json.dumps(result["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
