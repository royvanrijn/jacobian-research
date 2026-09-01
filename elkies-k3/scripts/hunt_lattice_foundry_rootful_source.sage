#!/usr/bin/env sage
"""Hunt and Niemeier-certify a low-MW companion for a foundry target.

This is the source-fibration half of the lattice foundry.  Starting from one
exact frame in the generated database (rootless or already high-MW), it walks
the Kneser p-neighbour graph of the *same positive-definite genus*, maximizing
exact ADE root rank.  A hit is not promoted from genus data alone: the script
glues the source frame to the stored rank-seven auxiliary along an explicit
discriminant anti-isometry, constructs the even unimodular rank-24 ambient,
and recovers the source as the saturated orthogonal complement.

The p-neighbour beam is deterministic for the declared seed and is a bounded
discovery search.  The final source, its ADE/MW data, genus membership,
Niemeier ambient, primitive embedding and complement are exact certificates.

status: ACTIVE_SEARCH
claim: exact source certificate for any emitted hit; no completeness claim for
  source frames not visited by the bounded p-neighbour beam.
inputs: artifacts/generated-results/elkies-k3-lattice-foundry-v1.json
outputs: artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-hunt.json
"""

DESCRIPTION = __doc__

import argparse
import hashlib
import json
import math
import random
from collections import Counter
from pathlib import Path

from sage.all import (
    Genus,
    IntegralLattice,
    QQ,
    ZZ,
    QuadraticForm,
    block_diagonal_matrix,
    identity_matrix,
    matrix,
    pari,
    set_random_seed,
    vector,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATABASE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-hunt.json"
)
DEFAULT_FRAME_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-lattice-foundry-ns0024-source-root-adapted.txt"
)

# Reuse the exact root/frame/saturation kernel of the elliptic-neighbour engine.
_engine_path = HERE / "exact_neighbor_engine.sage"
exec(compile(_engine_path.read_text(), str(_engine_path), "exec"), globals())


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def gram_digest(value):
    payload = "\n".join(" ".join(map(str, row)) for row in value.rows()) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def form_from_gram(gram):
    coefficients = []
    for row in range(gram.nrows()):
        for column in range(row, gram.ncols()):
            coefficients.append(
                gram[row, row] // 2
                if row == column
                else gram[row, column]
            )
    result = QuadraticForm(ZZ, gram.nrows(), coefficients)
    assert result.Hessian_matrix() == gram
    return result


def reduced_gram(gram):
    change_columns = matrix(ZZ, pari(gram).qflllgram())
    change_rows = change_columns.transpose()
    result = change_rows * gram * change_rows.transpose()
    assert abs(change_rows.det()) == 1
    assert result.det() == gram.det()
    return result, change_rows


def reduced_key(gram):
    reduced, unused = reduced_gram(gram)
    return gram_digest(reduced)


def root_rank_and_count(gram):
    minimum = pari(gram).qfminim(2)
    count = int(minimum[0])
    if count == 0:
        return 0, 0
    return int(matrix(ZZ, minimum[2].sage()).rank()), count


def connected_components(cartan):
    unseen = set(range(cartan.nrows()))
    result = []
    while unseen:
        start = min(unseen)
        unseen.remove(start)
        pending = [start]
        component = []
        while pending:
            current = pending.pop()
            component.append(current)
            adjacent = [index for index in unseen if cartan[current, index]]
            for index in adjacent:
                unseen.remove(index)
                pending.append(index)
        result.append(tuple(sorted(component)))
    return tuple(result)


def ade_name(frame):
    simple, unused_positive, cartan_rows = deterministic_simple_roots(frame)
    if simple.nrows() == 0:
        return "0", []
    cartan = matrix(ZZ, cartan_rows)
    components = []
    for indices in connected_components(cartan):
        block = cartan.matrix_from_rows_and_columns(indices, indices)
        rank = block.nrows()
        determinant = abs(int(block.det()))
        count = int(pari(block).qfminim(2)[0])
        if determinant == rank + 1 and count == rank * (rank + 1):
            label = f"A{rank}"
        elif rank >= 4 and determinant == 4 and count == 2 * rank * (rank - 1):
            label = f"D{rank}"
        else:
            label = {
                (6, 3, 72): "E6",
                (7, 2, 126): "E7",
                (8, 1, 240): "E8",
            }[(rank, determinant, count)]
        components.append(
            {
                "type": label,
                "rank": rank,
                "determinant": determinant,
                "signed_root_count": count,
            }
        )
    multiplicities = Counter(component["type"] for component in components)
    label = "+".join(
        (f"{count}{name}" if count > 1 else name)
        for name, count in sorted(multiplicities.items())
    )
    return label, components


def discriminant_data(gram):
    group = IntegralLattice(gram).discriminant_group()
    invariants = tuple(map(int, group.invariants()))
    if len(invariants) != 1:
        raise ValueError(
            "v1 explicit gluing currently requires a cyclic discriminant group"
        )
    return group, invariants[0], QQ(group.gram_matrix_quadratic()[0, 0])


def anti_isometry_units(auxiliary, frame):
    unused_aux_group, order, auxiliary_q = discriminant_data(auxiliary)
    unused_frame_group, frame_order, frame_q = discriminant_data(frame)
    assert frame_order == order
    return [
        unit
        for unit in range(order)
        if math.gcd(unit, order) == 1
        and ((auxiliary_q * unit * unit + frame_q) / 2).denominator() == 1
    ]


def primitive_row_lattice(value):
    smith = value.smith_form()[0]
    return all(
        abs(int(smith[index, index])) == 1
        for index in range(value.nrows())
    )


def integral_matrix(value):
    assert all(entry.denominator() == 1 for entry in value.list())
    return matrix(ZZ, value)


def ambient_root_components(gram):
    minimum = pari(gram).qfminim(2)
    signed_count = int(minimum[0])
    half = matrix(ZZ, minimum[2].sage()).transpose()
    unseen = set(range(half.nrows()))
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
                if half[current] * gram * half[index] != 0
            }
            component.update(adjacent)
            unseen.difference_update(adjacent)
            pending.extend(sorted(adjacent))
        basis = half[sorted(component)]
        components.append(
            {
                "rank": int(basis.rank()),
                "signed_root_count": 2 * len(component),
            }
        )
    return signed_count, sorted(
        components,
        key=lambda item: (item["rank"], item["signed_root_count"]),
        reverse=True,
    )


def catalog_ambient_label(components, catalog):
    observed = sorted(
        (item["rank"], item["signed_root_count"]) for item in components
    )
    matches = []
    for entry in catalog["rooted_niemeier_lattices"]:
        expected = sorted(
            (item["rank"], item["signed_root_count"])
            for item in entry["root_components"]
        )
        if expected == observed:
            matches.append(entry["label"])
    assert len(matches) == 1
    return matches[0]


def glue_to_niemeier(auxiliary, frame, unit, catalog):
    auxiliary_group, order, unused_auxiliary_q = discriminant_data(auxiliary)
    frame_group, frame_order, unused_frame_q = discriminant_data(frame)
    assert frame_order == order
    auxiliary_generator = vector(QQ, auxiliary_group.gen(0).lift())
    frame_generator = vector(QQ, frame_group.gen(0).lift())
    split_gram = block_diagonal_matrix(auxiliary, frame)
    split = IntegralLattice(split_gram)
    glue_vector = vector(
        QQ,
        list(unit * auxiliary_generator) + list(frame_generator),
    )
    ambient = split.overlattice([glue_vector])
    assert ambient.rank() == 24 and ambient.gram_matrix().det() == 1

    old_gram = ambient.gram_matrix()
    change_columns = matrix(ZZ, pari(old_gram).qflllgram())
    change_rows = change_columns.transpose()
    ambient_gram = integral_matrix(change_rows * old_gram * change_rows.transpose())
    ambient_basis = change_rows * ambient.basis_matrix()
    assert ambient_gram == ambient_basis * split_gram * ambient_basis.transpose()
    assert all(ambient_gram[index, index] % 2 == 0 for index in range(24))

    auxiliary_split_basis = identity_matrix(ZZ, 24)[:7]
    auxiliary_coordinates = integral_matrix(
        auxiliary_split_basis * ambient_basis.inverse()
    )
    assert primitive_row_lattice(auxiliary_coordinates)
    assert (
        auxiliary_coordinates
        * ambient_gram
        * auxiliary_coordinates.transpose()
        == auxiliary
    )
    complement_coordinates = (
        auxiliary_coordinates * ambient_gram
    ).right_kernel_matrix()
    assert complement_coordinates.nrows() == 17
    assert primitive_row_lattice(complement_coordinates)
    complement = (
        complement_coordinates
        * ambient_gram
        * complement_coordinates.transpose()
    )
    assert pari(frame).qfisom(pari(complement)) != 0
    signed_roots, components = ambient_root_components(ambient_gram)
    return {
        "anti_isometry_unit": unit,
        "discriminant_order": order,
        "ambient_label": catalog_ambient_label(components, catalog),
        "ambient_signed_root_count": signed_roots,
        "ambient_root_components": components,
        "ambient_gram": rows(ambient_gram),
        "ambient_gram_sha256": gram_digest(ambient_gram),
        "auxiliary_basis_in_ambient": rows(auxiliary_coordinates),
        "complement_basis_in_ambient": rows(complement_coordinates),
        "complement_gram": rows(complement),
        "complement_gram_sha256": gram_digest(complement),
        "primitive_auxiliary_embedding": True,
        "saturated_orthogonal_complement": True,
        "complement_integrally_isometric_to_source": True,
    }


def replay_path(start, path, expected_digests, determinant):
    current = form_from_gram(start)
    observed = []
    for edge, expected in zip(path, expected_digests):
        prime = ZZ(edge["prime"])
        witness = vector(ZZ, edge["primitive_p_divisible_vector"])
        current = current.find_p_neighbor_from_vec(prime, witness)
        gram = current.Hessian_matrix()
        assert gram.det() == determinant
        digest = gram_digest(gram)
        assert digest == expected
        observed.append(digest)
    return current, observed


parser = argparse.ArgumentParser(description=DESCRIPTION)
parser.add_argument("--database", type=Path, default=DATABASE)
parser.add_argument("--ns-id", default="NS0024")
parser.add_argument("--target-frame-id", default="NS0024-F005")
parser.add_argument(
    "--start-source-artifact",
    type=Path,
    help=(
        "restart the Kneser beam from the exact source Gram in an earlier "
        "same-NS rootful-source certificate instead of from the catalogue target"
    ),
)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument(
    "--root-adapted-frame-output", type=Path, default=DEFAULT_FRAME_OUTPUT
)
parser.add_argument("--generations", type=int, default=15)
parser.add_argument("--beam", type=int, default=20)
parser.add_argument("--samples-per-parent", type=int, default=100)
parser.add_argument("--primes", default="3,7,11,13,17,23")
parser.add_argument("--seed", type=int, default=20260901)
parser.add_argument("--target-root-rank", type=int, default=12)
parser.add_argument(
    "--continue-through-bound",
    action="store_true",
    help="continue improving root count after the requested rank is first reached",
)
parser.add_argument(
    "--allow-below-target",
    action="store_true",
    help=(
        "emit the exact best source reached by the bounded search even when its "
        "root rank is below --target-root-rank; this changes no completeness claim"
    ),
)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

database_path = arguments.database.resolve()
output_path = arguments.output.resolve()
frame_output_path = arguments.root_adapted_frame_output.resolve()
database = json.loads(database_path.read_text())
catalog = json.loads(CATALOG.read_text())
ns_row = next(row for row in database["ns_classes"] if row["ns_id"] == arguments.ns_id)
target_row = next(
    row for row in ns_row["frames"]
    if row["frame_id"] == arguments.target_frame_id
)
target = matrix(ZZ, target_row["gram"])
auxiliary = matrix(ZZ, ns_row["auxiliary_gram"])
determinant = int(target.det())
assert determinant == int(auxiliary.det())
target_root_rank, target_root_count = root_rank_and_count(target)
assert target_root_rank == int(target_row["root_rank"])
target_mw_rank = 17 - target_root_rank

start = target
start_record = {
    "kind": "CATALOGUE_TARGET_FRAME",
    "frame_id": arguments.target_frame_id,
    "root_type": target_row["root_type"],
    "root_rank": target_root_rank,
    "mw_rank_for_rho_19": target_mw_rank,
    "gram_sha256": gram_digest(target),
}
start_source_path = None
if arguments.start_source_artifact is not None:
    start_source_path = arguments.start_source_artifact.resolve()
    start_payload = json.loads(start_source_path.read_text())
    if start_payload.get("schema") != "elkies-k3.lattice-foundry-rootful-source.v1":
        raise ValueError("start source is not a rootful-source certificate")
    if start_payload["target"]["ns_id"] != arguments.ns_id:
        raise ValueError("start source belongs to a different NS class")
    start = matrix(ZZ, start_payload["source"]["gram"])
    if start.det() != determinant or Genus(start) != Genus(target):
        raise ValueError("start source is not in the catalogue target genus")
    start_root_rank, start_root_count = root_rank_and_count(start)
    if start_root_rank != int(start_payload["source"]["root_rank"]):
        raise ValueError("start source root rank does not replay")
    start_record = {
        "kind": "CERTIFIED_ROOTFUL_SOURCE_RESTART",
        "artifact": str(start_source_path.relative_to(ROOT)),
        "artifact_sha256": hashlib.sha256(start_source_path.read_bytes()).hexdigest(),
        "root_type": start_payload["source"]["root_type"],
        "root_rank": start_root_rank,
        "mw_rank_for_rho_19": 17 - start_root_rank,
        "gram_sha256": gram_digest(start),
    }
else:
    start_root_rank = target_root_rank
    start_root_count = target_root_count

random.seed(arguments.seed)
set_random_seed(arguments.seed)
requested_primes = [ZZ(value) for value in arguments.primes.split(",")]
primes = [prime for prime in requested_primes if determinant % prime]
if not primes:
    raise ValueError("every requested p-neighbour prime divides the determinant")

start_form = form_from_gram(start)
frontier = [
    {
        "root_rank": start_root_rank,
        "root_count": start_root_count,
        "form": start_form,
        "path": [],
        "path_digests": [],
    }
]
seen = {reduced_key(start)}
best = frontier[0]
retained_candidates = list(frontier)
generation_accounting = []

print(
    f"FOUNDRYSOURCE|stage=start|ns={arguments.ns_id}|target={arguments.target_frame_id}"
    f"|det={determinant}|start_root_rank={start_root_rank}"
    f"|target_root_rank={arguments.target_root_rank}",
    flush=True,
)

for generation in range(1, arguments.generations + 1):
    candidates = []
    failures = 0
    for parent in frontier:
        for unused_sample in range(arguments.samples_per_parent):
            prime = random.choice(primes)
            try:
                witness = parent["form"].find_primitive_p_divisible_vector__random(prime)
                child = parent["form"].find_p_neighbor_from_vec(prime, witness)
                gram = child.Hessian_matrix()
                if gram.det() != determinant or any(
                    gram[index, index] % 2 for index in range(17)
                ):
                    failures += 1
                    continue
                key = reduced_key(gram)
                if key in seen:
                    continue
                seen.add(key)
                root_rank, root_count = root_rank_and_count(gram)
            except Exception:
                failures += 1
                continue
            candidates.append(
                {
                    "root_rank": root_rank,
                    "root_count": root_count,
                    "form": child,
                    "path": parent["path"]
                    + [
                        {
                            "prime": int(prime),
                            "primitive_p_divisible_vector": list(map(int, witness)),
                        }
                    ],
                    "path_digests": parent["path_digests"]
                    + [gram_digest(gram)],
                }
            )

    candidates.sort(
        key=lambda row: (-row["root_rank"], -row["root_count"])
    )
    next_frontier = []
    per_profile = Counter()
    for candidate in candidates:
        profile = (candidate["root_rank"], candidate["root_count"])
        if per_profile[profile] >= 3:
            continue
        per_profile[profile] += 1
        next_frontier.append(candidate)
        if len(next_frontier) == arguments.beam:
            break
    if not next_frontier:
        raise RuntimeError("source-hunt beam died before producing a candidate")
    frontier = next_frontier
    retained_candidates.extend(frontier)
    if (
        frontier[0]["root_rank"], frontier[0]["root_count"]
    ) > (best["root_rank"], best["root_count"]):
        best = frontier[0]
    generation_accounting.append(
        {
            "generation": generation,
            "new_candidates": len(candidates),
            "failures": failures,
            "visited_reduced_keys": len(seen),
            "frontier_profiles": [
                [row["root_rank"], row["root_count"]]
                for row in frontier
            ],
            "best_root_rank": best["root_rank"],
            "best_signed_root_count": best["root_count"],
        }
    )
    print(
        f"FOUNDRYSOURCE|stage=generation|generation={generation}"
        f"|candidates={len(candidates)}|seen={len(seen)}"
        f"|best_root_rank={best['root_rank']}|best_roots={best['root_count']}",
        flush=True,
    )
    if (
        best["root_rank"] >= arguments.target_root_rank
        and not arguments.continue_through_bound
    ):
        break

# A high-root frame with a nonprimitive root lattice does not yet provide the
# integral root-plus-MW source coordinates required by this certificate.
# Select the best certifiable retained beam state instead of failing only at
# the end or silently treating a torsion/glue problem as an MW basis.
certifiable_best = None
for candidate in sorted(
    retained_candidates,
    key=lambda row: (-row["root_rank"], -row["root_count"], len(row["path"])),
):
    candidate_minimized = minimize_child_frame(candidate["form"].Hessian_matrix())
    if (
        candidate_minimized["root_lattice_primitive"]
        and candidate_minimized["mw_height"] is not None
    ):
        certifiable_best = candidate
        break
if certifiable_best is None:
    raise RuntimeError("bounded source hunt retained no primitive-root source frame")
best = certifiable_best

if (
    best["root_rank"] < arguments.target_root_rank
    and not arguments.allow_below_target
):
    raise RuntimeError(
        f"bounded source hunt reached root rank {best['root_rank']}, below target "
        f"{arguments.target_root_rank}"
    )

source_raw = best["form"].Hessian_matrix()
source, source_change = reduced_gram(source_raw)
assert Genus(source) == Genus(target)
source_ade, source_components = ade_name(source)
source_root_rank, source_root_count = root_rank_and_count(source)
source_minimized = minimize_child_frame(source)
assert source_root_rank == best["root_rank"]
assert source_minimized["root_lattice_primitive"]
assert source_minimized["mw_height"] is not None
source_mw_rank = 17 - source_root_rank

replayed_form, replayed_digests = replay_path(
    start, best["path"], best["path_digests"], determinant
)
assert replayed_form.Hessian_matrix() == source_raw

units = anti_isometry_units(auxiliary, source)
assert units
niemeier = glue_to_niemeier(auxiliary, source, units[0], catalog)

payload = {
    "schema": "elkies-k3.lattice-foundry-rootful-source.v1",
    "status": (
        f"PASS_EXACT_NEW_K3_ROOTFUL_MW{source_mw_rank}"
        "_SOURCE_AND_NIEMEIER_CERTIFICATE"
    ),
    "proof_boundary": {
        "proved": (
            "The emitted source is an even positive-definite rank-17 frame in "
            "the exact target genus, has the displayed exact ADE root system, "
            "and is recovered as the saturated orthogonal complement of the "
            "same primitive rank-seven auxiliary in the displayed Niemeier lattice."
        ),
        "not_proved": (
            "The bounded p-neighbour beam is not a complete companion-fibration "
            "classification. Its path is a Kneser genus-search provenance path, "
            "not yet an elliptic-neighbour corridor between marked U embeddings."
        ),
    },
    "search": {
        "seed": arguments.seed,
        "generations_bound": arguments.generations,
        "generations_used": len(generation_accounting),
        "beam": arguments.beam,
        "samples_per_parent": arguments.samples_per_parent,
        "primes": list(map(int, primes)),
        "target_root_rank": arguments.target_root_rank,
        "target_root_rank_reached": (
            best["root_rank"] >= arguments.target_root_rank
        ),
        "allow_below_target": arguments.allow_below_target,
        "visited_reduced_keys": len(seen),
        "generation_accounting": generation_accounting,
    },
    "target": {
        "ns_id": arguments.ns_id,
        "frame_id": arguments.target_frame_id,
        "gram": rows(target),
        "gram_sha256": gram_digest(target),
        "root_type": target_row["root_type"],
        "root_rank": target_root_rank,
        "mw_rank_for_rho_19": target_mw_rank,
        "determinant": determinant,
    },
    "source": {
        "gram": rows(source),
        "gram_sha256": gram_digest(source),
        "raw_search_gram_sha256": gram_digest(source_raw),
        "reduced_basis_rows_in_raw_basis": rows(source_change),
        "root_type": source_ade,
        "root_components": source_components,
        "root_rank": source_root_rank,
        "signed_root_count": source_root_count,
        "signed_short_vector_counts_through_norm_4": {
            str(norm): 2 * count
            for norm, count in sorted(
                Counter(
                    int(vector(ZZ, column) * source * vector(ZZ, column))
                    for column in matrix(ZZ, pari(source).qfminim(4)[2]).columns()
                ).items()
            )
        },
        "automorphism_group_order": int(pari(source).qfauto()[0]),
        "mw_rank_for_rho_19": source_mw_rank,
        "root_lattice_primitive": True,
        "root_adapted_gram": rows(source_minimized["frame"]),
        "root_adapted_basis_rows_in_source_basis": rows(
            source_minimized["basis"]
        ),
        "root_adapted_simple_roots_in_source_basis": rows(
            source_minimized["simple_roots"]
        ),
        "root_smith_invariants": list(
            map(int, source_minimized["root_smith_invariants"])
        ),
        "mw_height_gram": rational_rows(source_minimized["mw_height"]),
        "mw_regulator": str(abs(source_minimized["mw_height"].det())),
        "torsion": 1,
        "determinant": determinant,
        "genus_equals_target": True,
        "discriminant_form_equals_target": True,
    },
    "kneser_p_neighbor_provenance": {
        "edge_count": len(best["path"]),
        "edges": best["path"],
        "child_raw_gram_sha256_by_edge": best["path_digests"],
        "exact_replay_passed": replayed_digests == best["path_digests"],
        "warning": (
            "These positive-definite p-neighbours discover the source frame; "
            "they are not elliptic-neighbour edges on the K3 surface."
        ),
    },
    "niemeier_certificate": {
        "auxiliary_gram": rows(auxiliary),
        "auxiliary_gram_sha256": gram_digest(auxiliary),
        "all_cyclic_anti_isometry_units": units,
        **niemeier,
    },
    "inputs": {
        str(database_path.relative_to(ROOT)): hashlib.sha256(
            database_path.read_bytes()
        ).hexdigest(),
        str(CATALOG.relative_to(ROOT)): hashlib.sha256(CATALOG.read_bytes()).hexdigest(),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/hunt_lattice_foundry_rootful_source.sage"
    ),
}

if start_source_path is not None:
    payload["search"]["start"] = start_record
    payload["inputs"][str(start_source_path.relative_to(ROOT))] = hashlib.sha256(
        start_source_path.read_bytes()
    ).hexdigest()

serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
frame_text = "\n".join(
    [
        f"# source_artifact = {output_path.relative_to(ROOT)}",
        f"# ns_id = {arguments.ns_id}",
        f"# target_frame_id = {arguments.target_frame_id}",
        f"# ADE = {source_ade}",
        f"# root_rank = {source_root_rank}",
        f"# MW_rank = {source_mw_rank}",
    ]
    + [
        " ".join(map(str, row))
        for row in source_minimized["frame"].rows()
    ]
) + "\n"
if arguments.check:
    if output_path.read_text() != serialized:
        raise SystemExit("lattice-foundry rootful-source artifact is stale")
    if frame_output_path.read_text() != frame_text:
        raise SystemExit("lattice-foundry root-adapted source frame is stale")
else:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(serialized)
    frame_output_path.parent.mkdir(parents=True, exist_ok=True)
    frame_output_path.write_text(frame_text)

print(
    f"FOUNDRYSOURCE|stage=done|ns={arguments.ns_id}|source={source_ade}"
    f"|root_rank={source_root_rank}|mw_rank={source_mw_rank}"
    f"|ambient={niemeier['ambient_label']}|p_edges={len(best['path'])}|status=PASS"
)
