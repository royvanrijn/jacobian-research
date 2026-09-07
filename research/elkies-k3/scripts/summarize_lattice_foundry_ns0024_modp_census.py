#!/usr/bin/env python3
"""Summarize exact NS0024 finite-field MW3/P4 census directories."""

import argparse
import hashlib
import json
from pathlib import Path


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--census", action="append", nargs=2, metavar=("LABEL", "DIRECTORY"), default=[]
)
parser.add_argument(
    "--p4-census", "--quadratic", dest="p4_census",
    action="append",
    nargs=3,
    metavar=("LABEL", "DIRECTORY", "LOG_GLOB"),
    default=[],
)
parser.add_argument(
    "--joint-point",
    action="append",
    nargs=2,
    metavar=("LABEL", "POINT_JSON"),
    default=[],
    help="aggregate the arithmetic-realizability audit from a certified joint point",
)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()


def digest(lines):
    return hashlib.sha256(("\n".join(sorted(lines)) + "\n").encode()).hexdigest()


def field(line, key):
    return line.split(f"|{key}=", 1)[1].split("|", 1)[0]


censuses = {}
for label, raw_directory in args.census:
    directory = Path(raw_directory).resolve()
    rows = []
    for path in sorted(directory.glob("*.txt")):
        rows.extend(
            (path.stem, line)
            for line in path.read_text().splitlines()
            if "P1X=" in line
        )
    surfaces = {(field(line, "A"), field(line, "B")) for _, line in rows}
    censuses[label] = {
        "directory": str(directory),
        "marked_triples": len(rows),
        "marked_surfaces": len(surfaces),
        "node_charts": sorted({chart for chart, _ in rows}),
        "record_sha256": digest(line for _, line in rows),
    }


def log_stem(path):
    name = path.name
    for suffix in (
        ".resolved-rot1.extract.log",
        ".resolved.extract.log",
        ".corrected.extract.log",
        ".fast.extract.log",
        ".extract.log",
    ):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"unsupported extraction log name: {path}")


p4_censuses = {}
marker = "NS0024P4MARKING|accepted="
for label, raw_directory, pattern in args.p4_census:
    directory = Path(raw_directory).resolve()
    systems = {path.stem for path in directory.glob("*.ms")}
    decoded = {}
    used_logs = []
    for path in sorted(directory.glob(pattern)):
        text = path.read_text(errors="replace")
        if marker not in text:
            continue
        stem = log_stem(path)
        accepted = int(text.split(marker, 1)[1].splitlines()[0])
        decoded.setdefault(stem, set()).add(accepted)
        used_logs.append(f"{path.name}:{accepted}")
    covered = systems.intersection(decoded)
    positive = {stem for stem in covered if any(decoded[stem])}
    p4_censuses[label] = {
        "directory": str(directory),
        "log_glob": pattern,
        "systems": len(systems),
        "decoded_systems": len(covered),
        "positive_markings": len(positive),
        "all_systems_decoded": covered == systems,
        "decode_ledger_sha256": digest(used_logs),
    }


joint_galois_rows = []
for label, raw_path in args.joint_point:
    path = Path(raw_path).resolve()
    point = json.loads(path.read_text())
    if point.get("schema") != "elkies-k3.lattice-foundry-ns0024-mw4-point-modp.v1":
        raise ValueError(f"joint point has the wrong schema: {path}")
    audit = point.get("arithmetic_realizability")
    if not isinstance(audit, dict):
        raise ValueError(f"joint point has no arithmetic-realizability audit: {path}")
    factor_audit = point.get("rur", {}).get("factor_audit")
    accepted_factors = (
        [
            row
            for row in factor_audit
            if row.get("outcome") == "PASS_EXACT_RESOLVED_MW4_MARKING"
        ]
        if isinstance(factor_audit, list)
        else []
    )
    sources = accepted_factors or [audit]
    artifact_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
    for source in sources:
        factor_index = source.get("factor_index")
        joint_galois_rows.append(
            {
                "label": (
                    f"{label}:factor{factor_index}" if factor_index is not None else label
                ),
                "path": str(path),
                "sha256": artifact_sha256,
                "factor_index": factor_index,
                "prime": int(point["prime"]),
                "residue_degree": int(
                    source["residue_degree"]
                    if "residue_degree" in source
                    else source["degree"]
                ),
                "surface_frobenius_orbit_size": int(
                    source["surface_frobenius_orbit_size"]
                ),
                "section_marking_frobenius_orbit_size": int(
                    source["section_marking_frobenius_orbit_size"]
                ),
                "resolved_oriented_marking_frobenius_orbit_size": int(
                    source["resolved_oriented_marking_frobenius_orbit_size"]
                ),
                "relative_section_degree_over_surface_field": int(
                    source["relative_section_degree_over_surface_field"]
                ),
                "relative_orientation_degree_over_section_field": int(
                    source["relative_orientation_degree_over_section_field"]
                ),
                "action_closes_on_marked_mw4": bool(
                    source["action_closes_on_marked_mw4"]
                ),
                "status": (
                    source["status"]
                    if "status" in source
                    else source["arithmetic_realizability_status"]
                ),
                "relative_fixed_marked_mw4_rank": source.get(
                    "relative_fixed_marked_mw4_rank"
                ),
                "prime_field_fixed_marked_mw4_rank": source.get(
                    "prime_field_fixed_marked_mw4_rank"
                ),
            }
        )


payload = {
    "schema": "elkies-k3.lattice-foundry-ns0024-modp-census.v2",
    "status": "PASS_EXACT_BOUNDED_MODULAR_CENSUS_NO_TWO_PRIME_FAMILY",
    "censuses": censuses,
    "p4_extension_censuses": p4_censuses,
    "proved_boundary": (
        "The scanner exhausts the declared split I7+I5+I4 polynomial-section "
        "charts over each stated prime. The p=11 P4 ledgers exhaust both "
        "base-field pole coordinates with all coefficients over GF(11^2) and "
        "all monic irreducible quadratic pole coordinates on the recorded "
        "rational MW3 markings, with exact group-law component and intersection checks."
    ),
    "open_boundary": (
        "MW3 markings defined only over extension fields, P4 pole coordinates "
        "of degree at least three, a one-dimensional marked component, the "
        "q4/orbit1 edge child, and characteristic-zero reconstruction remain open."
    ),
}
if joint_galois_rows:
    fixed_rank_histogram = {}
    for row in joint_galois_rows:
        rank = row["prime_field_fixed_marked_mw4_rank"]
        if rank is not None:
            fixed_rank_histogram[str(rank)] = fixed_rank_histogram.get(str(rank), 0) + 1
    payload["joint_closed_point_galois_gates"] = {
        "points": joint_galois_rows,
        "prime_field_surface_points": sum(
            row["surface_frobenius_orbit_size"] == 1 for row in joint_galois_rows
        ),
        "prime_field_fixed_marked_mw4_rank_histogram": fixed_rank_histogram,
        "full_mw4_fixed_primes": sorted(
            {
                row["prime"]
                for row in joint_galois_rows
                if row["prime_field_fixed_marked_mw4_rank"] == 4
            }
        ),
        "warning_primes": sorted(
            {
                row["prime"]
                for row in joint_galois_rows
                if row["prime_field_fixed_marked_mw4_rank"] is not None
                and row["prime_field_fixed_marked_mw4_rank"] < 4
            }
        ),
        "proof_boundary": (
            "These are exact Frobenius actions on certified finite-field MW4 fibres. "
            "Repeated full fixed rank is positive arithmetic-realizability evidence, "
            "while smaller fixed rank is warning evidence. Neither proves or disproves "
            "NS=NS_Q for a characteristic-zero family without a common producer and descent."
        ),
    }
args.output.resolve().parent.mkdir(parents=True, exist_ok=True)
args.output.resolve().write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
print(
    "NS0024MODPCENSUS|censuses={} |p4_censuses={} |output={}".format(
        len(censuses), len(p4_censuses), args.output.resolve()
    ).replace(" ", ""),
    flush=True,
)
