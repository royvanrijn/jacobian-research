#!/usr/bin/env sage
"""Canonicalize completed N(24A1) shards under W(24A1) semidirect M24.

Each input completion shard is already canonicalized under residual M24 on
positive octad supports.  This script merges contiguous prefix ranges and
applies the missing Weyl sign action without enumerating 2^24 signs.

For a primitive auxiliary K, write its doubled physical coordinate matrix as
C (seven rows and 24 columns), so C*C^t=2*Gram(K).  An intrinsic isometry of
K and an ambient signed coordinate permutation identify two embeddings
exactly when the 24 columns agree up to independent signs and an M24
permutation.  Thus the Weyl quotient is the M24 orbit problem for the ordered
partition of coordinates by signed column covectors.  GAP
RepresentativeAction(...,OnTuplesSets) supplies an exact transporter.

Every retained orbit includes explicit row-isometry, M24-permutation, and
coordinate-sign witnesses for all input members.  Full stabilizer orders in
2^24 semidirect M24 are computed independently from intrinsic automorphisms,
ordered-partition stabilizers, and zero-column sign freedom.

status: EXACT_FULL_WEYL_M24_CANONICALIZATION_OF_DECLARED_INPUT_SHARDS
"""

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from sage.all import (
    GF,
    QQ,
    ZZ,
    identity_matrix,
    libgap,
    matrix,
    pari,
    vector,
)
from sage.coding.golay_code import GolayCode


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_MANIFEST = (
    ROOT
    / "artifacts/generated-results/elkies-k3-24a1-octad-completion-manifest-v1.json"
)


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def construction_a_physical_basis():
    code = GolayCode(GF(2), extended=True)
    generator = code.generator_matrix().echelon_form()
    assert generator[:, :12] == identity_matrix(GF(2), 12)
    basis = matrix(QQ, 24, 24)
    for index in range(12):
        basis[index] = vector(
            QQ, [QQ(int(entry)) / 2 for entry in generator[index]]
        )
    for index in range(12):
        basis[12 + index, 12 + index] = 1
    return code, basis


def doubled_physical_coordinates(record, physical_basis):
    value = 2 * matrix(ZZ, record["auxiliary_ambient_basis"]) * physical_basis
    assert all(entry in ZZ for entry in value.list())
    value = matrix(ZZ, value)
    gram = matrix(ZZ, record["auxiliary_gram"])
    assert value * value.transpose() == 2 * gram
    return value


def row_isometry(source, target):
    """Return U with U*source*U^t=target, or None."""
    raw = pari(source).qfisom(pari(target))
    if raw == 0:
        return None
    result = matrix(ZZ, raw).inverse().transpose().change_ring(ZZ)
    assert abs(result.det()) == 1
    assert result * source * result.transpose() == target
    return result


def row_automorphisms(gram):
    data = pari(gram).qfauto()
    generators = [matrix(ZZ, item).transpose() for item in data[1]]
    identity = identity_matrix(ZZ, gram.nrows())
    found = {tuple(identity.list()): identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            candidate = current * generator
            key = tuple(candidate.list())
            if key not in found:
                assert candidate * gram * candidate.transpose() == gram
                found[key] = candidate
                frontier.append(candidate)
    assert len(found) == int(data[0])
    ordered = sorted(found.values(), key=lambda item: tuple(item.list()))
    ordered.remove(identity)
    return [identity] + ordered


def signed_column_key(column):
    result = tuple(map(int, column))
    for entry in result:
        if entry:
            return tuple(-value for value in result) if entry < 0 else result
    return result


def signed_column_partition(value):
    result = defaultdict(list)
    for index, column in enumerate(value.columns(), start=1):
        result[signed_column_key(column)].append(index)
    return result


def partition_transporter(group_gap, source, target):
    source_partition = signed_column_partition(source)
    target_partition = signed_column_partition(target)
    keys = sorted(source_partition)
    if keys != sorted(target_partition):
        return None
    if any(
        len(source_partition[key]) != len(target_partition[key]) for key in keys
    ):
        return None
    transporter = libgap.RepresentativeAction(
        group_gap,
        libgap([libgap.Set(source_partition[key]) for key in keys]),
        libgap([libgap.Set(target_partition[key]) for key in keys]),
        libgap.OnTuplesSets,
    )
    return None if transporter == libgap.fail else transporter


def signed_permutation_witness(source, target, transporter):
    images = [
        int(libgap.OnPoints(index, transporter)) for index in range(1, 25)
    ]
    signs = []
    transformed = matrix(ZZ, source.nrows(), 24)
    for source_index, target_index_one_based in enumerate(images):
        source_column = source.column(source_index)
        target_column = target.column(target_index_one_based - 1)
        if source_column == target_column:
            sign = 1
        else:
            assert source_column == -target_column
            sign = -1
        signs.append(sign)
        transformed[:, target_index_one_based - 1] = sign * source_column
    assert transformed == target
    return images, signs


def embedding_equivalence(
    candidate_standard,
    representative_standard,
    automorphisms,
    group_gap,
):
    for automorphism_index, automorphism in enumerate(automorphisms):
        moved = automorphism * candidate_standard
        transporter = partition_transporter(
            group_gap, moved, representative_standard
        )
        if transporter is None:
            continue
        images, signs = signed_permutation_witness(
            moved, representative_standard, transporter
        )
        return {
            "standard_automorphism_index_zero_based": automorphism_index,
            "standard_automorphism": automorphism,
            "coordinate_permutation_images_one_based": images,
            "coordinate_signs_at_source_coordinates": signs,
        }
    return None


def full_stabilizer_data(standard, automorphisms, group, group_gap):
    partition = signed_column_partition(standard)
    keys = sorted(partition)
    gap_partition = libgap(
        [libgap.Set(partition[key]) for key in keys]
    )
    partition_stabilizer_order = int(
        libgap.Size(
            libgap.Stabilizer(
                group_gap, gap_partition, libgap.OnTuplesSets
            )
        )
    )
    compatible_intrinsic_automorphisms = sum(
        partition_transporter(group_gap, automorphism * standard, standard)
        is not None
        for automorphism in automorphisms
    )
    zero_columns = sum(column == 0 for column in standard.columns())
    stabilizer_order = (
        compatible_intrinsic_automorphisms
        * partition_stabilizer_order
        * (2 ** zero_columns)
    )
    full_group_order = (2 ** 24) * int(group.order())
    assert full_group_order % stabilizer_order == 0
    return {
        "compatible_intrinsic_automorphisms": (
            compatible_intrinsic_automorphisms
        ),
        "ordered_partition_stabilizer_order_in_M24": (
            partition_stabilizer_order
        ),
        "zero_coordinate_covectors": zero_columns,
        "free_weyl_sign_factor": 2 ** zero_columns,
        "full_embedding_stabilizer_order": stabilizer_order,
        "full_embedding_orbit_size": full_group_order // stabilizer_order,
    }


def merge_completion_shards(paths, payloads):
    assert len(paths) == len(payloads) and paths
    ordered_inputs = sorted(
        zip(paths, payloads),
        key=lambda item: item[1]["parameters"][
            "prefix_start_zero_based_inclusive"
        ],
    )
    shards = []
    records = []
    determinant_bound = None
    prefix_artifact_hash = None
    for shard_index, (path, payload) in enumerate(ordered_inputs):
        assert payload["schema"] == (
            "elkies-k3.24a1-octad-rank7-completion-shard.v1"
        )
        assert payload["status"] == (
            "PASS_EXACT_DECLARED_24A1_OCTAD_COMPLETION_SHARD"
        )
        parameters = payload["parameters"]
        start = parameters["prefix_start_zero_based_inclusive"]
        stop = parameters["prefix_stop_zero_based_exclusive"]
        assert start < stop
        if determinant_bound is None:
            determinant_bound = parameters["determinant_bound"]
            prefix_artifact_hash = payload["input"]["prefix_artifact_sha256"]
        else:
            assert parameters["determinant_bound"] == determinant_bound
            assert payload["input"]["prefix_artifact_sha256"] == prefix_artifact_hash
        shards.append(
            {
                "shard_index_zero_based": shard_index,
                "artifact": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "prefix_start_zero_based_inclusive": start,
                "prefix_stop_zero_based_exclusive": stop,
                "residual_m24_orbits": len(payload["orbits"]),
            }
        )
        for local_index, source in enumerate(payload["orbits"]):
            record = dict(source)
            record["_origin"] = {
                "shard_index_zero_based": shard_index,
                "local_residual_m24_orbit_index_zero_based": local_index,
                "artifact": str(path.relative_to(ROOT)),
            }
            records.append(record)
    assert all(
        left["prefix_stop_zero_based_exclusive"]
        == right["prefix_start_zero_based_inclusive"]
        for left, right in zip(shards, shards[1:])
    )
    return {
        "parameters": {
            "prefix_start_zero_based_inclusive": shards[0][
                "prefix_start_zero_based_inclusive"
            ],
            "prefix_stop_zero_based_exclusive": shards[-1][
                "prefix_stop_zero_based_exclusive"
            ],
            "determinant_bound": determinant_bound,
        },
        "shards": shards,
        "orbits": records,
    }


def canonicalize(payload):
    records = payload["orbits"]
    parameters = payload["parameters"]
    prefix_start = parameters["prefix_start_zero_based_inclusive"]
    prefix_stop = parameters["prefix_stop_zero_based_exclusive"]

    code, physical_basis = construction_a_physical_basis()
    group = code.permutation_automorphism_group()
    assert int(group.order()) == 244823040
    group_gap = libgap(group)
    full_group_order = (2 ** 24) * int(group.order())

    intrinsic_classes = []
    standardized = {}
    for input_index, record in enumerate(records):
        gram = matrix(ZZ, record["auxiliary_gram"])
        matched = None
        for class_index, item in enumerate(intrinsic_classes):
            if gram.det() != item["gram"].det():
                continue
            isometry = row_isometry(gram, item["gram"])
            if isometry is not None:
                matched = class_index, item, isometry
                break
        if matched is None:
            item = {
                "gram": gram,
                "automorphisms": row_automorphisms(gram),
                "members": [],
                "embedding_orbits": [],
            }
            intrinsic_classes.append(item)
            class_index = len(intrinsic_classes) - 1
            isometry = identity_matrix(ZZ, 7)
        else:
            class_index, item, isometry = matched
        physical = doubled_physical_coordinates(record, physical_basis)
        standard = isometry * physical
        assert standard * standard.transpose() == 2 * item["gram"]
        entry = {
            "input_index": input_index,
            "intrinsic_class_index": class_index,
            "standardizing_isometry": isometry,
            "standard_coordinates": standard,
        }
        item["members"].append(entry)
        standardized[input_index] = entry

    for class_index, item in enumerate(intrinsic_classes):
        for member in item["members"]:
            input_index = member["input_index"]
            found = None
            for orbit_index, orbit in enumerate(item["embedding_orbits"]):
                witness = embedding_equivalence(
                    member["standard_coordinates"],
                    orbit["representative_standard_coordinates"],
                    item["automorphisms"],
                    group_gap,
                )
                if witness is not None:
                    found = orbit_index, orbit, witness
                    break
            if found is None:
                item["embedding_orbits"].append(
                    {
                        "representative_input_index": input_index,
                        "representative_standard_coordinates": member[
                            "standard_coordinates"
                        ],
                        "members": [input_index],
                        "standard_witnesses": {
                            input_index: {
                                "standard_automorphism_index_zero_based": 0,
                                "standard_automorphism": identity_matrix(ZZ, 7),
                                "coordinate_permutation_images_one_based": list(
                                    range(1, 25)
                                ),
                                "coordinate_signs_at_source_coordinates": [1] * 24,
                            }
                        },
                    }
                )
            else:
                unused_orbit_index, orbit, witness = found
                orbit["members"].append(input_index)
                orbit["standard_witnesses"][input_index] = witness
        print(
            f"WEYLM24|intrinsic_class={class_index + 1}/{len(intrinsic_classes)}"
            f"|det={item['gram'].det()}|input={len(item['members'])}"
            f"|orbits={len(item['embedding_orbits'])}",
            flush=True,
        )

    output_orbits = []
    output_classes = []
    for class_index, item in enumerate(intrinsic_classes):
        class_orbit_ids = []
        for orbit in item["embedding_orbits"]:
            representative_index = orbit["representative_input_index"]
            representative_record = records[representative_index]
            representative_entry = standardized[representative_index]
            witnesses = []
            for input_index in orbit["members"]:
                source_entry = standardized[input_index]
                standard_witness = orbit["standard_witnesses"][input_index]
                standard_automorphism = standard_witness[
                    "standard_automorphism"
                ]
                row_map = (
                    representative_entry["standardizing_isometry"].inverse()
                    * standard_automorphism
                    * source_entry["standardizing_isometry"]
                ).change_ring(ZZ)
                source_gram = matrix(ZZ, records[input_index]["auxiliary_gram"])
                representative_gram = matrix(
                    ZZ, representative_record["auxiliary_gram"]
                )
                assert abs(row_map.det()) == 1
                assert (
                    row_map * source_gram * row_map.transpose()
                    == representative_gram
                )
                source_physical = doubled_physical_coordinates(
                    records[input_index], physical_basis
                )
                representative_physical = doubled_physical_coordinates(
                    representative_record, physical_basis
                )
                moved = row_map * source_physical
                images = standard_witness[
                    "coordinate_permutation_images_one_based"
                ]
                signs = standard_witness[
                    "coordinate_signs_at_source_coordinates"
                ]
                transformed = matrix(ZZ, 7, 24)
                for source_coordinate, target_coordinate in enumerate(images):
                    transformed[:, target_coordinate - 1] = (
                        signs[source_coordinate]
                        * moved.column(source_coordinate)
                    )
                assert transformed == representative_physical
                # The displayed signed ambient coordinate identity is stronger
                # than a separate frame qfisom call: an ambient isometry that
                # carries K_source to K_representative necessarily carries
                # their saturated orthogonal complements isometrically.
                assert (
                    records[input_index]["mordell_weil_rank"]
                    == representative_record["mordell_weil_rank"]
                )
                assert (
                    records[input_index]["k3_discriminant_gate"]
                    == representative_record["k3_discriminant_gate"]
                )
                witnesses.append(
                    {
                        "input_index_zero_based": input_index,
                        "input_origin": records[input_index]["_origin"],
                        "row_isometry_source_to_representative": rows(row_map),
                        "coordinate_permutation_images_one_based": images,
                        "coordinate_signs_at_source_coordinates": signs,
                    }
                )

            stabilizer = full_stabilizer_data(
                orbit["representative_standard_coordinates"],
                item["automorphisms"],
                group,
                group_gap,
            )
            orbit_id = f"W24A1-{len(output_orbits) + 1:04d}"
            class_orbit_ids.append(orbit_id)
            output_orbits.append(
                {
                    "orbit_id": orbit_id,
                    "intrinsic_auxiliary_class_id": (
                        f"K24A1-{class_index + 1:04d}"
                    ),
                    "representative_input_index_zero_based": representative_index,
                    "representative_input_origin": representative_record["_origin"],
                    "representative_octad_indices_one_based": representative_record[
                        "octad_indices_one_based"
                    ],
                    "representative_auxiliary_basis_in_ambient": representative_record[
                        "auxiliary_ambient_basis"
                    ],
                    "representative_auxiliary_gram": representative_record[
                        "auxiliary_gram"
                    ],
                    "representative_complement_basis_in_ambient": representative_record[
                        "frame_ambient_basis"
                    ],
                    "representative_frame_gram": representative_record[
                        "frame_gram"
                    ],
                    "frame_discriminant_form_normal_key": representative_record[
                        "k3_discriminant_gate"
                    ]["frame_discriminant_form_normal_key"],
                    "ternary_genus_representatives": representative_record[
                        "k3_discriminant_gate"
                    ]["ternary_genus_representatives"],
                    "input_shard_local_record_indices_zero_based": orbit[
                        "members"
                    ],
                    "input_shard_local_records_collapsed": len(orbit["members"]),
                    "determinant": representative_record["determinant"],
                    "frame_root_system": representative_record[
                        "frame_root_system"
                    ],
                    "mordell_weil_rank": representative_record[
                        "mordell_weil_rank"
                    ],
                    "matching_even_ternary_genera": representative_record[
                        "k3_discriminant_gate"
                    ]["matching_even_ternary_genera"],
                    "full_group_stabilizer": stabilizer,
                    "member_witnesses": witnesses,
                }
            )
        output_classes.append(
            {
                "intrinsic_auxiliary_class_id": f"K24A1-{class_index + 1:04d}",
                "determinant": int(item["gram"].det()),
                "representative_gram": rows(item["gram"]),
                "automorphism_group_order": len(item["automorphisms"]),
                "input_shard_local_records": len(item["members"]),
                "full_weyl_m24_embedding_orbits": len(
                    item["embedding_orbits"]
                ),
                "embedding_orbit_ids": class_orbit_ids,
            }
        )

    mw_distribution = Counter(
        record["mordell_weil_rank"] for record in output_orbits
    )
    determinant_distribution = Counter(
        record["determinant"] for record in output_orbits
    )
    compatible = sum(
        record["matching_even_ternary_genera"] > 0
        for record in output_orbits
    )
    if prefix_start == 0 and prefix_stop == 500 and len(payload["shards"]) == 2:
        assert len(records) == 291
        assert len(intrinsic_classes) == 5
        assert len(output_orbits) == 16
        assert compatible == 13
        assert determinant_distribution == Counter(
            {384: 1, 448: 1, 480: 10, 486: 1, 500: 3}
        )
        assert mw_distribution == Counter({12: 4, 13: 10, 14: 2})
    if prefix_start == 0 and prefix_stop == 1000 and len(payload["shards"]) == 4:
        assert len(records) == 675
        assert len(intrinsic_classes) == 5
        assert len(output_orbits) == 21
        assert compatible == 16
        assert determinant_distribution == Counter(
            {384: 3, 448: 1, 480: 13, 486: 1, 500: 3}
        )
        assert mw_distribution == Counter({12: 5, 13: 13, 14: 3})
    if prefix_start == 0 and prefix_stop == 2000 and len(payload["shards"]) == 8:
        assert len(records) == 1267
        assert len(intrinsic_classes) == 5
        assert len(output_orbits) == 23
        assert compatible == 18
        assert determinant_distribution == Counter(
            {384: 3, 448: 1, 480: 13, 486: 1, 500: 5}
        )
        assert mw_distribution == Counter({12: 5, 13: 13, 14: 5})
    if (
        prefix_start == 0
        and prefix_stop == 10547
        and len(payload["shards"]) == 43
    ):
        assert len(records) == 3051
        assert len(intrinsic_classes) == 5
        assert len(output_orbits) == 24
        assert compatible == 18
        assert determinant_distribution == Counter(
            {384: 3, 448: 2, 480: 13, 486: 1, 500: 5}
        )
        assert mw_distribution == Counter({12: 5, 13: 13, 14: 6})
    return {
        "schema": "elkies-k3.24a1-weyl-m24-canonicalization.v2",
        "status": "PASS_EXACT_FULL_WEYL_M24_CANONICALIZATION_OF_DECLARED_INPUT_SHARDS",
        "proof_scope": {
            "proved": (
                f"The {len(records)} residual-M24 embedding records in the declared "
                "contiguous input shards are partitioned exactly into full W(24A1) semidirect "
                "M24 embedding orbits. Every membership has an explicit "
                "intrinsic row isometry, M24 coordinate permutation, and Weyl "
                "sign witness, and every full stabilizer order is certified."
            ),
            "not_proved": (
                "The result inherits the input boundary: only positive seven-"
                f"octad generators and five-prefix indices {prefix_start}:{prefix_stop} are covered. "
                "It is not a complete determinant-500 24A1 auxiliary census."
            ),
        },
        "ambient_group": {
            "name": "W(24A1) semidirect M24",
            "weyl_order": 2 ** 24,
            "m24_order": int(group.order()),
            "full_order": full_group_order,
        },
        "parameters": parameters,
        "completion_shards": payload["shards"],
        "method": {
            "intrinsic_isometry": "PARI qfisom and full qfauto enumeration",
            "embedding_invariant": (
                "24 doubled physical coordinate covectors modulo independent signs"
            ),
            "exact_residual_transporter": (
                "GAP RepresentativeAction(M24,ordered signed-column partitions,OnTuplesSets)"
            ),
            "stabilizer_formula": (
                "compatible intrinsic automorphisms times M24 ordered-partition "
                "stabilizer times 2^(zero coordinate covectors)"
            ),
        },
        "accounting": {
            "input_shard_local_residual_m24_records": len(records),
            "intrinsic_auxiliary_isometry_classes": len(intrinsic_classes),
            "full_weyl_m24_embedding_orbits": len(output_orbits),
            "k3_compatible_full_embedding_orbits_by_ternary_genus_gate": compatible,
            "determinant_distribution": {
                str(key): value for key, value in sorted(determinant_distribution.items())
            },
            "mordell_weil_rank_distribution": {
                str(key): value for key, value in sorted(mw_distribution.items())
            },
        },
        "intrinsic_auxiliary_classes": output_classes,
        "embedding_orbits": output_orbits,
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--input", type=Path, action="append")
parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
parser.add_argument("--output", type=Path)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
manifest_metadata = None
if arguments.input:
    input_paths = arguments.input
else:
    manifest_bytes = arguments.manifest.read_bytes()
    manifest = json.loads(manifest_bytes)
    assert manifest["schema"] == (
        "elkies-k3.24a1-octad-completion-manifest.v1"
    )
    assert manifest["status"] == (
        "PASS_EXACT_CONTIGUOUS_24A1_OCTAD_COMPLETION_SHARD_MANIFEST"
    )
    input_paths = [ROOT / row["artifact"] for row in manifest["shards"]]
    assert all(
        hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
        for path, row in zip(input_paths, manifest["shards"])
    )
    manifest_metadata = {
        "artifact": str(arguments.manifest.relative_to(ROOT)),
        "sha256": hashlib.sha256(manifest_bytes).hexdigest(),
    }
input_payloads = [json.loads(path.read_text()) for path in input_paths]
merged = merge_completion_shards(input_paths, input_payloads)
if manifest_metadata is not None:
    assert merged["parameters"] == {
        key: manifest["parameters"][key]
        for key in (
            "prefix_start_zero_based_inclusive",
            "prefix_stop_zero_based_exclusive",
            "determinant_bound",
        )
    }
    assert len(merged["orbits"]) == manifest["accounting"][
        "shard_local_residual_m24_records"
    ]
payload = canonicalize(merged)
if manifest_metadata is not None:
    payload["completion_manifest"] = manifest_metadata
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
output = arguments.output or (
    ROOT
    / "artifacts/generated-results"
    / (
        "elkies-k3-24a1-weyl-m24-canonicalization-"
        f"{merged['parameters']['prefix_start_zero_based_inclusive']:05d}-"
        f"{merged['parameters']['prefix_stop_zero_based_exclusive']:05d}-v2.json"
    )
)
if arguments.check:
    if not output.exists() or output.read_text() != encoded:
        raise SystemExit("24A1 Weyl-M24 canonicalization artifact is stale")
else:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded)
print(
    "WEYLM24|input={}|intrinsic={}|full_orbits={}|k3_compatible={}"
    "|status=PASS_EXACT".format(
        payload["accounting"]["input_shard_local_residual_m24_records"],
        payload["accounting"]["intrinsic_auxiliary_isometry_classes"],
        payload["accounting"]["full_weyl_m24_embedding_orbits"],
        payload["accounting"][
            "k3_compatible_full_embedding_orbits_by_ternary_genus_gate"
        ],
    ),
    flush=True,
)
