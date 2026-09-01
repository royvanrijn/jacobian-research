#!/usr/bin/env sage
"""Build a certificate-first Picard-19 elliptic-K3 lattice foundry.

The first checked search shell deliberately stays small: it mutates either
non-root generator of the two complete H3 J2 control embeddings by one signed
ambient root in N(2A7+2D5).  The auxiliary generator is heuristic outside
that declared shell.  Everything after generation -- saturation, complement,
root data, integral-isometry deduplication, discriminant forms, and ternary
realizability -- is exact.

This script consumes, rather than reimplements, the hash-pinned Niemeier
catalogue and the completed D5/residual-Weyl/primitive-closure control census.
It also loads ``exact_neighbor_engine.sage`` for the shared root/frame kernel.

status: ACTIVE_SEARCH
claim: exact census inside the JSON-declared one-root mutation shell; no claim
  of completeness for all rank-seven auxiliaries below the determinant bound.
inputs: elkies-k3/data/lattice-foundry/one-root-control-shell-v1.json and the
  existing rooted-Niemeier/J2 control artifacts.
outputs: artifacts/generated-results/elkies-k3-lattice-foundry-v1.json; the
  cumulative human report is maintained separately after later route phases.
"""

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from sage.all import Genus, QQ, ZZ, block_diagonal_matrix, matrix, pari, vector
from sage.quadratic_forms.genera.genus import genera


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_CONFIG = (
    ROOT / "elkies-k3/data/lattice-foundry/one-root-control-shell-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
)
DEFAULT_REPORT = ROOT / "elkies-k3/LATTICE_FOUNDRY_REPORT_2026-09-01.md"
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
CONTROLS = ROOT / "artifacts/generated-results/elkies-k3-rootless-j2-niemeier-first.json"
CONTROL_GLUE = ROOT / "artifacts/generated-results/elkies-k3-rootless-j2-niemeier-controls.json"
PUBLISHED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
ALTERNATE = ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"

# Reuse the exact root/frame implementation used by the neighbour engine.  It
# is ordinary Python-with-Sage imports, so executing the source avoids Sage's
# preparser when this driver is launched with ``sage -python``.
_engine_path = HERE / "exact_neighbor_engine.sage"
exec(compile(_engine_path.read_text(), str(_engine_path), "exec"), globals())


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def gram_digest(value):
    payload = "\n".join(" ".join(map(str, row)) for row in value.rows()) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def normal_form_key(discriminant_form):
    normal = discriminant_form.normal_form()
    return {
        "invariants": list(map(int, normal.invariants())),
        "quadratic_gram": rows_as_strings(normal.gram_matrix_quadratic()),
        "value_module": str(normal.value_module_qf()),
    }


def rows_as_strings(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def discriminant_length(gram):
    smith = gram.smith_form()[0]
    invariants = [abs(int(smith[index, index])) for index in range(gram.nrows())]
    return sum(value > 1 for value in invariants), invariants


def signed_roots(gram):
    minimum = pari(gram).qfminim(2)
    half = [vector(ZZ, column) for column in matrix(ZZ, minimum[2].sage()).columns()]
    result = half + [-root for root in half]
    assert len(result) == int(minimum[0])
    return result


def root_type(frame):
    roots, unused_basis, data = roots_and_data(frame)
    if not roots:
        return "0", 0, 0, 1, []
    unseen = set(range(len(roots)))
    components = []
    while unseen:
        component = {min(unseen)}
        pending = list(component)
        unseen.difference_update(component)
        while pending:
            current = pending.pop()
            adjacent = {
                index
                for index in unseen
                if roots[current] * frame * roots[index] != 0
            }
            component.update(adjacent)
            unseen.difference_update(adjacent)
            pending.extend(sorted(adjacent))
        component_roots = matrix(ZZ, [roots[index] for index in sorted(component)])
        rank = int(component_roots.rank())
        count = len(component)
        if count == rank * (rank + 1):
            label = f"A{rank}"
            determinant = rank + 1
        elif rank >= 4 and count == 2 * rank * (rank - 1):
            label = f"D{rank}"
            determinant = 4
        elif (rank, count) in {(6, 72), (7, 126), (8, 240)}:
            label = f"E{rank}"
            determinant = {6: 3, 7: 2, 8: 1}[rank]
        else:
            raise AssertionError(f"unrecognized ADE component {(rank, count)}")
        components.append((label, rank, count, determinant))
    components.sort()
    multiplicities = Counter(label for label, unused_rank, unused_count, unused_det in components)
    label = "+".join(
        (f"{multiplicity}{name}" if multiplicity > 1 else name)
        for name, multiplicity in sorted(multiplicities.items())
    )
    root_determinant = math.prod(item[3] for item in components)
    return label, int(data[0]), int(data[1]), int(root_determinant), components


def exact_isometric(left, right):
    if left.nrows() != right.nrows() or left.det() != right.det():
        return False
    return pari(left).qfisom(pari(right)) != 0


def canonical_reduced_gram(gram):
    change = gram.LLL_gram()
    reduced = change.transpose() * gram * change
    assert abs(change.det()) == 1 and reduced.det() == gram.det()
    return reduced, change


def primitive_closure(basis):
    original = matrix(ZZ, basis)
    saturated = original.row_module(ZZ).saturation().basis_matrix()
    assert saturated.nrows() == original.nrows()
    coordinates = original * saturated.pseudoinverse()
    assert all(value.denominator() == 1 for value in coordinates.list())
    closure_index = abs(int(matrix(ZZ, coordinates).det()))
    return saturated, closure_index


def complement_record(ambient, auxiliary_basis):
    complement_basis = (auxiliary_basis * ambient).right_kernel_matrix()
    frame = complement_basis * ambient * complement_basis.transpose()
    assert complement_basis.nrows() == 17
    assert frame.is_positive_definite()
    assert frame.det() == (auxiliary_basis * ambient * auxiliary_basis.transpose()).det()
    return complement_basis, frame


def find_ternary_realizations(frame):
    determinant = abs(ZZ(frame.det()))
    target = normal_form_key(Genus(frame).discriminant_form())
    matches = []
    all_genera = genera((2, 1), determinant, even=True)
    for genus_index, genus in enumerate(all_genera):
        if normal_form_key(genus.discriminant_form()) != target:
            continue
        representative = matrix(ZZ, genus.representative())
        assert representative.det() == -determinant
        matches.append(
            {
                "genus_index": genus_index,
                "gram": rows(representative),
                "determinant": int(representative.det()),
                "local_symbols": [str(symbol) for symbol in genus.local_symbols()],
                "discriminant_form_normal_key": target,
            }
        )
    return len(all_genera), matches


def frame_intrinsics(frame, norm_bound):
    minimum = pari(frame).qfminim(norm_bound)
    columns = matrix(ZZ, minimum[2].sage()).columns()
    by_norm = Counter(int(vector(ZZ, column) * frame * vector(ZZ, column)) for column in columns)
    signed_theta = {"0": 1}
    signed_theta.update(
        {str(norm): 2 * count for norm, count in sorted(by_norm.items())}
    )
    minimum_norm = min(by_norm) if by_norm else None
    norm_four = 2 * by_norm.get(4, 0)
    short_cosets = {
        tuple(int(entry % 2) for entry in vector(ZZ, column))
        for column in columns
    }
    automorphism_order = int(pari(frame).qfauto()[0])
    return {
        "minimum_squared_norm": minimum_norm,
        "theta_coefficients_by_squared_norm_through_bound": signed_theta,
        "theta_squared_norm_bound": norm_bound,
        "norm_four_vectors": norm_four,
        "norm_four_unoriented_pairs": norm_four // 2,
        "automorphism_group_order": automorphism_order,
        "short_cosets_mod_2_hit_through_bound": len(short_cosets),
        "short_coset_squared_norm_bound": norm_bound,
        "hermite_invariant": float(QQ(minimum_norm) / (QQ(frame.det()) ** (QQ(1) / 17))),
    }


def report_text(payload, output_path, config_path):
    accounting = payload["accounting"]
    shortlist = payload["equation_shortlist"]
    lines = [
        "# Lattice foundry report — 2026-09-01",
        "",
        "## Outcome",
        "",
        (
            "The first deterministic foundry shell found "
            f"**{accounting['new_rootless_frame_classes']} rootless rank-17 frame "
            "classes on determinant-changing Picard-19 lattice classes** after exact "
            "saturation and ternary-realizability gates."
        ),
        "",
        (
            f"The run contains {accounting['auxiliary_isometry_classes']} auxiliary "
            f"classes, {accounting['frame_isometry_classes']} frame classes, and "
            f"{accounting['rootless_frame_classes']} rootless frames including the "
            "two locked H3 controls."
        ),
        "",
        "No candidate is promoted to equation work: the best companion found in this "
        "one-root shell has root rank only "
        f"{accounting['maximum_companion_root_rank_for_new_target']}, and no certified "
        "low-cost lattice-neighbour corridor has yet been run. This is a strong lattice "
        "discovery catalogue, not success level A.",
        "",
        "## Exact gates passed",
        "",
        "- Every generated rank-seven auxiliary is saturated in its Niemeier ambient.",
        "- Determinant and discriminant length are filtered at the declared bounds.",
        "- Every rank-17 complement is the exact integral orthogonal complement.",
        "- Root systems and Mordell--Weil ranks use exact norm-two enumeration.",
        "- Auxiliary and frame classes are deduplicated with PARI integral isometry.",
        "- Every retained NS class has an explicit even ternary signature-(2,1) genus representative with the required opposite discriminant form.",
        "- The existing complete H3 control branch is replayed: published R17 and alternate Q80 are recovered, and there is no third rootless H3 J2 frame.",
        "",
        "## Search boundary",
        "",
        payload["completeness_boundary"]["complete"],
        "",
        payload["completeness_boundary"]["heuristic"],
        "",
        "In particular, this report does **not** classify every primitive rank-seven "
        "Niemeier auxiliary of determinant at most 5,000. A negative result from this "
        "shell would not have been a theorem. The exact H3 J2 result remains separate "
        "and complete only for the pinned determinant-948 auxiliary.",
        "",
        "## Ranked discovery shortlist",
        "",
        "| rank | target | disc | N4 pairs | Aut | best companion | source MW | route |",
        "| ---: | --- | ---: | ---: | ---: | --- | ---: | --- |",
    ]
    for row in shortlist:
        lines.append(
            "| {rank} | `{target}` | {disc} | {n4} | {aut} | {root} | {mw} | {route} |".format(
                rank=row["shortlist_rank"],
                target=row["target_frame_id"],
                disc=row["determinant"],
                n4=row["target_invariants"]["norm_four_unoriented_pairs"],
                aut=row["target_invariants"]["automorphism_group_order"],
                root=row["source_root_type"],
                mw=row["source_mw_rank"],
                route=row["route_status"],
            )
        )
    lines.extend(
        [
            "",
            "The score is exploratory and combines determinant, norm-four sections, "
            "automorphism-reduced shell size, and companion root rank. It is not a "
            "specialization-rank prediction or theorem status.",
            "",
            "## Route and equation frontier",
            "",
            "The route ledger intentionally distinguishes a same-auxiliary Niemeier "
            "companion from a certified elliptic-neighbour path. No edge is emitted "
            "without a primitive isotropic class, physical Weyl replay, horizontal-wall "
            "test, and integral marking transport. The next calculation is therefore a "
            "multi-objective neighbour search from the highest-root companion/target "
            "pairs, using q in 4, 6, 8, 12, 24 and old-fibre degree at most four.",
            "",
            "No characteristic-zero ansatz has been started. The current sources have "
            "too much MW rank for the preferred equation gate, so Phase 7 correctly "
            "remains closed.",
            "",
            "## Reproduce",
            "",
            "```bash",
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \\",
            "  elkies-k3/scripts/build_lattice_foundry.sage",
            "",
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \\",
            "  elkies-k3/scripts/build_lattice_foundry.sage --check",
            "```",
            "",
            f"Configuration: [`{config_path.relative_to(report_path.parent)}`]({config_path.relative_to(report_path.parent)})",
            "",
            f"Database: [`{output_path.relative_to(ROOT)}`](../{output_path.relative_to(ROOT)})",
            "",
        ]
    )
    return "\n".join(lines)


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

config_path = arguments.config.resolve()
output_path = arguments.output.resolve()
report_path = arguments.report.resolve()
config = json.loads(config_path.read_text())
assert config["schema"] == "elkies-k3.lattice-foundry-search.v1"

catalog_payload = json.loads(CATALOG.read_text())
control_payload = json.loads(CONTROLS.read_text())
control_glue_payload = json.loads(CONTROL_GLUE.read_text())
alternate_payload = json.loads(ALTERNATE.read_text())
published = load_matrix(PUBLISHED)
alternate = matrix(ZZ, alternate_payload["rootless_frame"])

ambient_by_label = {
    entry["label"]: matrix(ZZ, entry["gram"])
    for entry in catalog_payload["rooted_niemeier_lattices"]
}
assert config["ambient_labels"] == ["2A7_2D5"]
ambient_label = config["ambient_labels"][0]
ambient = ambient_by_label[ambient_label]
ambient_roots = signed_roots(ambient)
assert len(ambient_roots) == 192

seed_rows = control_payload["rootless_classes"]
assert len(seed_rows) == 2
assert sum(row["matches_published_R17"] for row in seed_rows) == 1
assert sum(row["matches_alternate_Q80"] for row in seed_rows) == 1

generated = []
raw_attempts = 0
rank_rejected = 0
determinant_rejected = 0
length_rejected = 0
for seed in seed_rows:
    seed_name = (
        "published_R17" if seed["matches_published_R17"] else "alternate_Q80"
    )
    basis = matrix(ZZ, seed["representative_embedding"]["auxiliary_basis_in_ambient"])
    for mutable_row in config["extension"]["mutable_rows_zero_based"]:
        for root_index, root in enumerate(ambient_roots):
            raw_attempts += 1
            candidate = matrix(ZZ, basis)
            candidate[mutable_row] += root
            if candidate.rank() != 7:
                rank_rejected += 1
                continue
            saturated, closure_index = primitive_closure(candidate)
            auxiliary = saturated * ambient * saturated.transpose()
            determinant = int(auxiliary.det())
            if not (0 < determinant <= config["determinant_bound"]):
                determinant_rejected += 1
                continue
            length, smith = discriminant_length(auxiliary)
            if length > config["discriminant_length_bound"]:
                length_rejected += 1
                continue
            complement_basis, frame = complement_record(ambient, saturated)
            root_label, root_rank, root_count, root_det, root_components = root_type(frame)
            generated.append(
                {
                    "seed": seed_name,
                    "mutable_row_zero_based": mutable_row,
                    "ambient_root_index": root_index,
                    "ambient_root": list(map(int, root)),
                    "closure_index": closure_index,
                    "auxiliary_basis": saturated,
                    "auxiliary_gram": auxiliary,
                    "auxiliary_smith": smith,
                    "discriminant_length": length,
                    "complement_basis": complement_basis,
                    "frame": frame,
                    "root_type": root_label,
                    "root_rank": root_rank,
                    "root_count": root_count,
                    "root_determinant": root_det,
                    "root_components": root_components,
                }
            )

print(
    "FOUNDRY|stage=generated|attempts={}|retained={}".format(
        raw_attempts, len(generated)
    ),
    flush=True,
)

# Group by exact auxiliary isometry.  This is the NS/T class key in the
# Kneser--Nishiyama construction; every embedding remains attached below.
auxiliary_classes = []
for record in generated:
    matched = None
    for class_row in auxiliary_classes:
        if exact_isometric(record["auxiliary_gram"], class_row["auxiliary_gram"]):
            matched = class_row
            break
    if matched is None:
        matched = {
            "auxiliary_gram": record["auxiliary_gram"],
            "records": [],
            "frame_classes": [],
        }
        auxiliary_classes.append(matched)
    matched["records"].append(record)
    frame_match = None
    for frame_class in matched["frame_classes"]:
        if exact_isometric(record["frame"], frame_class["frame"]):
            frame_match = frame_class
            break
    if frame_match is None:
        frame_match = {"frame": record["frame"], "records": []}
        matched["frame_classes"].append(frame_match)
    frame_match["records"].append(record)

print(
    "FOUNDRY|stage=isometry|auxiliary_classes={}|frame_classes={}".format(
        len(auxiliary_classes),
        sum(len(row["frame_classes"]) for row in auxiliary_classes),
    ),
    flush=True,
)

# Locked H3 regression.  This is not inferred from the mutation shell: it is
# cross-checked against the complete control census and its exact accounting.
pinned_auxiliary = matrix(ZZ, control_payload["standard_auxiliary"]["pinned_gram"])
pinned_class = next(
    row for row in auxiliary_classes if exact_isometric(row["auxiliary_gram"], pinned_auxiliary)
)
pinned_rootless = [
    row for row in pinned_class["frame_classes"] if root_type(row["frame"])[1] == 0
]
assert len(pinned_rootless) == 2
assert sum(exact_isometric(row["frame"], published) for row in pinned_rootless) == 1
assert sum(exact_isometric(row["frame"], alternate) for row in pinned_rootless) == 1
assert control_payload["accounting"]["rootless_complement_isometry_classes"] == 2
assert control_payload["accounting"]["primitive_rootless_embeddings_in_cover"] == 12
assert sorted(
    row["embedding_count_in_enumerated_cover"] for row in control_payload["rootless_classes"]
) == [4, 8]
assert control_glue_payload["niemeier_accounting"]["control_ambient_root_system"] == "2A7+2D5"

ns_rows = []
rootless_rows = []
companion_rows = []
frame_count = 0
new_rootless = 0
for ns_index, class_row in enumerate(auxiliary_classes, start=1):
    auxiliary = class_row["auxiliary_gram"]
    determinant = int(auxiliary.det())
    reduced_auxiliary, auxiliary_change = canonical_reduced_gram(auxiliary)
    representative_frame = class_row["frame_classes"][0]["frame"]
    genus_count, ternary = find_ternary_realizations(representative_frame)
    if not ternary:
        # Keep the rejected class visible rather than silently dropping it.
        realizability = "REJECT_NO_TERNARY_REALIZATION"
    else:
        realizability = "PASS_EXACT_TERNARY_DISCRIMINANT_FORM"
    ns_id = f"NS{ns_index:04d}"
    frames = []
    for local_index, frame_class in enumerate(class_row["frame_classes"], start=1):
        frame_count += 1
        frame = frame_class["frame"]
        reduced_frame, frame_change = canonical_reduced_gram(frame)
        label, root_rank, root_count, root_det, components = root_type(frame)
        minimized = minimize_child_frame(frame)
        mw_rank = 17 - root_rank
        frame_id = f"{ns_id}-F{local_index:03d}"
        intrinsics = None
        if root_rank == 0:
            intrinsics = frame_intrinsics(frame, config["theta_squared_norm_bound"])
        mw_height = minimized["mw_height"]
        predicted_regulator = None
        if mw_height is not None:
            predicted_regulator = str(abs(mw_height.det()))
        embedding_records = []
        for embedding in frame_class["records"]:
            embedding_records.append(
                {
                    "seed": embedding["seed"],
                    "mutable_row_zero_based": embedding["mutable_row_zero_based"],
                    "ambient_root_index": embedding["ambient_root_index"],
                    "ambient_root": embedding["ambient_root"],
                    "closure_index": embedding["closure_index"],
                    "ambient": ambient_label,
                    "auxiliary_basis_in_ambient": rows(embedding["auxiliary_basis"]),
                    "complement_basis_in_ambient": rows(embedding["complement_basis"]),
                }
            )
        frame_payload = {
            "frame_id": frame_id,
            "gram": rows(frame),
            "gram_sha256": gram_digest(frame),
            "reduced_gram": rows(reduced_frame),
            "reduced_basis_columns_in_frame_basis": rows(frame_change),
            "determinant": determinant,
            "root_type": label,
            "root_rank": root_rank,
            "signed_root_count": root_count,
            "root_determinant": root_det,
            "mw_rank_for_rho_19": mw_rank,
            "root_lattice_primitive": minimized["root_lattice_primitive"],
            "root_smith_invariants": list(map(int, minimized["root_smith_invariants"])),
            "torsion_possibilities": (
                [1] if minimized["root_lattice_primitive"] else "requires_glue_analysis"
            ),
            "determinant_predicted_mw_regulator": predicted_regulator,
            "embedding_count_in_declared_shell": len(embedding_records),
            "embeddings": embedding_records,
            "rootless_intrinsics": intrinsics,
        }
        frames.append(frame_payload)
        companion_rows.append(
            {
                "ns_id": ns_id,
                "frame_id": frame_id,
                "root_type": label,
                "root_rank": root_rank,
                "mw_rank_for_rho_19": mw_rank,
                "torsion_possibilities": frame_payload["torsion_possibilities"],
                "determinant_predicted_mw_regulator": predicted_regulator,
            }
        )
        if root_rank == 0 and ternary:
            is_h3 = determinant == 948 and (
                exact_isometric(frame, published) or exact_isometric(frame, alternate)
            )
            if not is_h3:
                new_rootless += 1
            rootless_rows.append(
                {
                    "ns_id": ns_id,
                    "frame_id": frame_id,
                    "determinant": determinant,
                    "is_existing_H3_control": is_h3,
                    "invariants": intrinsics,
                }
            )
    ns_gram = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -representative_frame)
    ns_rows.append(
        {
            "ns_id": ns_id,
            "auxiliary_gram": rows(auxiliary),
            "auxiliary_gram_sha256": gram_digest(auxiliary),
            "auxiliary_reduced_gram": rows(reduced_auxiliary),
            "auxiliary_reduced_basis_columns": rows(auxiliary_change),
            "determinant": determinant,
            "discriminant_length": discriminant_length(auxiliary)[0],
            "discriminant_form_normal_key": normal_form_key(
                Genus(representative_frame).discriminant_form()
            ),
            "ns_gram_representative": rows(ns_gram),
            "ns_signature": [1, 18],
            "ternary_genus_count_at_determinant": genus_count,
            "ternary_realizations": ternary,
            "k3_primitive_embedding_certificate": {
                "status": realizability,
                "method": (
                    "The displayed T has signature (2,1), determinant -disc(NS), "
                    "and discriminant form opposite to NS. Gluing along the certified "
                    "finite quadratic-module anti-isometry gives an even unimodular "
                    "lattice of signature (3,19), hence the K3 lattice by uniqueness."
                ),
                "explicit_22_dimensional_glue_basis": "not_exported_v1",
            },
            "embedding_count_in_declared_shell": len(class_row["records"]),
            "frames": frames,
        }
    )

# Pair every new rootless target with every rootful companion on the same exact
# auxiliary class.  This is an inventory, not a neighbour-path claim.
route_ledger = []
pairs = []
ns_by_id = {row["ns_id"]: row for row in ns_rows}
for target in rootless_rows:
    if target["is_existing_H3_control"]:
        continue
    sources = [
        row
        for row in companion_rows
        if row["ns_id"] == target["ns_id"] and row["root_rank"] > 0
    ]
    for source in sources:
        route_id = f"R-{source['frame_id']}-{target['frame_id']}"
        route_ledger.append(
            {
                "route_id": route_id,
                "ns_id": target["ns_id"],
                "source_frame_id": source["frame_id"],
                "target_frame_id": target["frame_id"],
                "status": "NOT_YET_ENUMERATED",
                "same_auxiliary_companion_proved": True,
                "certified_neighbor_edges": [],
                "search_bounds_planned": {
                    "q": [4, 6, 8, 12, 24],
                    "old_fibre_degree_at_most": 4,
                    "objective": [
                        "resolved_RR_dimension",
                        "old_fibre_degree",
                        "physical_Weyl_reflections",
                        "horizontal_P_dot_O",
                        "edge_count"
                    ]
                },
                "proof_boundary": (
                    "Sharing an exact auxiliary/NS class does not itself exhibit an "
                    "elliptic-neighbour route or transport the full U marking."
                )
            }
        )
        pairs.append((target, source, route_id))

def discovery_score(target, source):
    invariants = target["invariants"]
    return (
        100.0 * source["root_rank"]
        + 5000.0 / target["determinant"]
        + invariants["norm_four_unoriented_pairs"] / 100.0
        - math.log2(invariants["automorphism_group_order"])
    )


best_pair_by_target = {}
for target, source, route_id in pairs:
    score = discovery_score(target, source)
    old = best_pair_by_target.get(target["frame_id"])
    if old is None or score > old[0]:
        best_pair_by_target[target["frame_id"]] = (score, target, source, route_id)

ranked = sorted(best_pair_by_target.values(), key=lambda row: (-row[0], row[1]["determinant"]))
shortlist = []
for shortlist_rank, (score, target, source, route_id) in enumerate(
    ranked[: config["shortlist_size"]], start=1
):
    shortlist.append(
        {
            "shortlist_rank": shortlist_rank,
            "discovery_score": score,
            "score_status": "EXPLORATORY_NOT_THEOREM",
            "ns_id": target["ns_id"],
            "target_frame_id": target["frame_id"],
            "determinant": target["determinant"],
            "target_invariants": target["invariants"],
            "source_frame_id": source["frame_id"],
            "source_root_type": source["root_type"],
            "source_root_rank": source["root_rank"],
            "source_mw_rank": source["mw_rank_for_rho_19"],
            "route_id": route_id,
            "route_status": "OPEN",
            "equation_gate": "CLOSED_SOURCE_MW_TOO_HIGH_AND_ROUTE_UNCERTIFIED",
            "weierstrass_ansatz": None,
            "moduli_dimension": 1,
            "parameter_space_plausibility": "not_assessed_before_route_gate",
        }
    )

max_companion_root = max(
    (row["root_rank"] for row in companion_rows if row["root_rank"] > 0),
    default=0,
)
payload = {
    "schema": "elkies-k3.lattice-foundry-database.v1",
    "status": "PASS_EXACT_DECLARED_SHELL_NEW_K3_TARGETS_ROUTE_GATE_OPEN",
    "search_specification": config,
    "positive_controls": {
        "status": "PASS_COMPLETE_H3_J2_REGRESSION",
        "ambient": "N(2A7+2D5)",
        "primitive_rootless_embeddings_in_complete_cover": 12,
        "cover_counts_by_frame": {"published_R17": 8, "alternate_Q80": 4},
        "rootless_frame_classes": 2,
        "third_rootless_frame": False,
        "published_R17_recovered": True,
        "alternate_Q80_recovered": True,
        "source_artifact": str(CONTROLS.relative_to(ROOT)),
    },
    "accounting": {
        "raw_mutations_attempted": raw_attempts,
        "rank_rejected": rank_rejected,
        "determinant_rejected": determinant_rejected,
        "discriminant_length_rejected": length_rejected,
        "retained_primitive_embeddings": len(generated),
        "auxiliary_isometry_classes": len(auxiliary_classes),
        "frame_isometry_classes": frame_count,
        "rootless_frame_classes": len(rootless_rows),
        "new_rootless_frame_classes": new_rootless,
        "ns_classes_without_ternary_realization": sum(
            not row["ternary_realizations"] for row in ns_rows
        ),
        "target_source_pairs": len(pairs),
        "maximum_companion_root_rank_for_new_target": max_companion_root,
    },
    "ns_classes": ns_rows,
    "rootless_targets": rootless_rows,
    "companion_fibrations": companion_rows,
    "route_ledger": route_ledger,
    "equation_shortlist": shortlist,
    "completeness_boundary": config["completeness_boundary"],
    "theorem_status": {
        "proved": (
            "The census and every arithmetic/lattice gate are exact inside the "
            "declared one-root mutation shell. The H3 positive control is independently "
            "complete at J2 and has exactly two rootless frame classes."
        ),
        "not_proved": (
            "The auxiliary generator is not complete under determinant 5000, the "
            "companion inventory is not the full fibration inventory, and no new "
            "source-to-target neighbour route or equation is yet certified."
        ),
    },
    "inputs": {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in (config_path, CATALOG, CONTROLS, CONTROL_GLUE, PUBLISHED, ALTERNATE)
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/build_lattice_foundry.sage"
    ),
}

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("lattice-foundry database is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)

# The report began as a Phase-1 rendering, but it now synthesizes the later
# exact source and route certificates too.  Do not overwrite that cumulative
# mathematical note when replaying the narrower mutation-shell database.

print(
    "FOUNDRY|stage=done|auxiliary_classes={}|frames={}|rootless={}|new={}|"
    "pairs={}|shortlist={}|status=PASS".format(
        len(auxiliary_classes), frame_count, len(rootless_rows), new_rootless,
        len(pairs), len(shortlist)
    )
)
