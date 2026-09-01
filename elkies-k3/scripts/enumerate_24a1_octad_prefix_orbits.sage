#!/usr/bin/env sage
"""Enumerate M24-orbits of Golay-octad subsets through size five.

This is the canonical-prefix engine for the 24A1 rank-seven auxiliary
backend.  It induces the exact M24 action from the extended Golay code on all
759 octads and uses set-stabilizer augmentation.  Candidates are bucketed by
exact M24-invariants and deduplicated only by an exact group transporter.

For each size k, every (k+1)-set has a deletion in one of the retained
k-orbits.  Extending a retained representative by the setwise-stabilizer
orbits therefore covers every next orbit.  The independent orbit-stabilizer
mass identity proves both coverage and absence of duplicate retained orbits.

The size-five frontier is intended for determinant-aware two-octad completion
to rank seven.  This script does not claim that every rank-seven auxiliary is
octad-generated and does not close a determinant band.

status: EXACT_FINITE_ORBIT_ENUMERATION
output: artifacts/generated-results/elkies-k3-24a1-octad-prefix-orbits-v1.json
"""

import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sage.all import (
    GF,
    Graph,
    Permutation,
    PermutationGroup,
    Set,
    ZZ,
    libgap,
    matrix,
)
from sage.coding.golay_code import GolayCode


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-24a1-octad-prefix-orbits-v1.json"
)


def mask(word):
    return sum(1 << index for index, entry in enumerate(word) if entry)


def image_mask(value, permutation):
    result = 0
    for coordinate in range(24):
        if (value >> coordinate) & 1:
            result |= 1 << (permutation(coordinate + 1) - 1)
    return result


def intersection_gram(subset, octads):
    selected = [octads[index - 1] for index in subset]
    return matrix(
        ZZ,
        len(selected),
        len(selected),
        lambda row, column: (
            4
            if row == column
            else (selected[row] & selected[column]).bit_count() // 2
        ),
    )


def canonical_edge_key(gram):
    graph = Graph(gram.nrows())
    for row in range(gram.nrows()):
        for column in range(row + 1, gram.ncols()):
            graph.add_edge(row, column, int(gram[row, column]))
    canonical = graph.canonical_label(edge_labels=True)
    return tuple(
        int(label)
        for unused_left, unused_right, label in canonical.edges(
            sort=True, labels=True
        )
    )


def build_action():
    code = GolayCode(GF(2), extended=True)
    octads = sorted(mask(word) for word in code if word.hamming_weight() == 8)
    assert len(octads) == 759 and len(set(octads)) == 759
    index = {value: position + 1 for position, value in enumerate(octads)}
    coordinate_group = code.permutation_automorphism_group()
    assert coordinate_group.order() == 244823040
    induced_generators = []
    for generator in coordinate_group.gens():
        images = [index[image_mask(value, generator)] for value in octads]
        induced_generators.append(Permutation(images))
    octad_group = PermutationGroup(induced_generators, domain=range(1, 760))
    assert octad_group.order() == coordinate_group.order()
    assert len(octad_group.orbits()) == 1
    relation = np.array(
        [
            [(left & right).bit_count() // 2 for right in octads]
            for left in octads
        ],
        dtype=np.uint8,
    )
    relation_counts = Counter(map(int, relation[0]))
    assert relation_counts == Counter({0: 30, 1: 448, 2: 280, 4: 1})
    return code, octads, coordinate_group, octad_group, relation


def invariant_signature(subset, octads, relation):
    gram = intersection_gram(subset, octads)
    union = 0
    for index in subset:
        union |= octads[index - 1]

    # For each one of the 759 ambient octads, sort its intersection numbers
    # with the selected subset.  The histogram is unchanged by M24 and by a
    # reordering of the selected octads.  It is only a bucket key: final orbit
    # equality is always decided by RepresentativeAction(...,OnSets).
    values = np.sort(relation[np.array(subset) - 1, :], axis=0)
    codes = np.zeros(759, dtype=np.int64)
    for row in values:
        codes = 5 * codes + row
    histogram = tuple(
        map(int, np.bincount(codes, minlength=5 ** len(subset)))
    )
    return (
        canonical_edge_key(gram),
        union.bit_count(),
        histogram,
    )


def equivalent(group_gap, left, right):
    return (
        libgap.RepresentativeAction(
            group_gap,
            libgap.Set(list(left)),
            libgap.Set(list(right)),
            libgap.OnSets,
        )
        != libgap.fail
    )


def extend_orbits(representatives, group, group_gap, octads, relation, size):
    candidates = set()
    for parent_index, subset in enumerate(representatives, start=1):
        stabilizer = group.stabilizer(Set(subset), action="OnSets")
        selected = set(subset)
        for orbit in stabilizer.orbits():
            representative = orbit[0]
            if representative in selected:
                continue
            candidates.add(tuple(sorted(subset + (representative,))))
        if parent_index % 1000 == 0:
            print(
                f"OCTADPREFIX|stage=augment|target_size={size}"
                f"|parents={parent_index}/{len(representatives)}"
                f"|candidates={len(candidates)}",
                flush=True,
            )

    buckets = defaultdict(list)
    retained = []
    duplicates = 0
    for candidate_index, subset in enumerate(sorted(candidates), start=1):
        signature = invariant_signature(subset, octads, relation)
        bucket = buckets[signature]
        if any(equivalent(group_gap, subset, prior) for prior in bucket):
            duplicates += 1
        else:
            bucket.append(subset)
            retained.append(subset)
        if candidate_index % 5000 == 0:
            print(
                f"OCTADPREFIX|stage=dedup|size={size}"
                f"|candidates={candidate_index}/{len(candidates)}"
                f"|orbits={len(retained)}|duplicates={duplicates}",
                flush=True,
            )
    return retained, {
        "size": size,
        "parent_orbits": len(representatives),
        "augmentation_candidates": len(candidates),
        "invariant_buckets": len(buckets),
        "exact_duplicate_candidates": duplicates,
        "retained_orbits": len(retained),
    }


def orbit_records(representatives, group, octads, size):
    records = []
    mass = 0
    determinant_histogram = Counter()
    union_histogram = Counter()
    group_order = int(group.order())
    for index, subset in enumerate(representatives, start=1):
        gram = intersection_gram(subset, octads)
        determinant = int(gram.det())
        union = 0
        for octad_index in subset:
            union |= octads[octad_index - 1]
        stabilizer_order = int(
            group.stabilizer(Set(subset), action="OnSets").order()
        )
        assert group_order % stabilizer_order == 0
        orbit_size = group_order // stabilizer_order
        mass += orbit_size
        determinant_histogram[determinant] += 1
        union_histogram[union.bit_count()] += 1
        records.append(
            {
                "orbit_id": f"O{size}-{index:05d}",
                "octad_indices_one_based": list(map(int, subset)),
                "octad_masks_hex": [
                    format(octads[octad_index - 1], "06x")
                    for octad_index in subset
                ],
                "intersection_gram_upper_triangle": [
                    int(gram[row, column])
                    for row in range(size)
                    for column in range(row, size)
                ],
                "raw_span_rank": int(gram.rank()),
                "raw_span_determinant": determinant,
                "coordinate_union_size": union.bit_count(),
                "coordinate_root_rank_in_orthogonal_complement": (
                    24 - union.bit_count()
                ),
                "setwise_stabilizer_order_in_M24": stabilizer_order,
                "orbit_size": orbit_size,
            }
        )
        if index % 1000 == 0:
            print(
                f"OCTADPREFIX|stage=mass|size={size}"
                f"|orbits={index}/{len(representatives)}|mass={mass}",
                flush=True,
            )
    expected_mass = math.comb(759, size)
    assert mass == expected_mass
    return records, {
        "orbit_count": len(records),
        "orbit_stabilizer_mass": mass,
        "expected_subset_count": expected_mass,
        "mass_identity_passed": True,
        "raw_span_determinant_histogram": {
            str(key): determinant_histogram[key]
            for key in sorted(determinant_histogram)
        },
        "coordinate_union_size_histogram": {
            str(key): union_histogram[key] for key in sorted(union_histogram)
        },
    }


def build(maximum_size):
    assert 1 <= maximum_size <= 5
    code, octads, coordinate_group, group, relation = build_action()
    group_gap = group._libgap_()
    representatives = [(1,)]
    layers = []
    accounting = []
    records, layer_accounting = orbit_records(
        representatives, group, octads, 1
    )
    layers.append({"size": 1, "orbits": records})
    accounting.append(
        {
            "size": 1,
            "parent_orbits": 0,
            "augmentation_candidates": 1,
            "invariant_buckets": 1,
            "exact_duplicate_candidates": 0,
            "retained_orbits": 1,
            **layer_accounting,
        }
    )
    for size in range(2, maximum_size + 1):
        representatives, augmentation = extend_orbits(
            representatives, group, group_gap, octads, relation, size
        )
        records, layer_accounting = orbit_records(
            representatives, group, octads, size
        )
        layers.append({"size": size, "orbits": records})
        accounting.append({**augmentation, **layer_accounting})
        print(
            f"OCTADPREFIX|stage=layer|size={size}"
            f"|orbits={len(representatives)}"
            f"|mass={layer_accounting['orbit_stabilizer_mass']}|status=PASS",
            flush=True,
        )

    generator_images = [
        [int(generator(point)) for point in range(1, 760)]
        for generator in group.gens()
    ]
    action_digest = hashlib.sha256(
        json.dumps(generator_images, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "schema": "elkies-k3.24a1-octad-prefix-orbits.v1",
        "status": f"PASS_EXACT_M24_OCTAD_SUBSET_ORBITS_THROUGH_SIZE_{maximum_size}",
        "proof_scope": {
            "proved": (
                f"All M24 orbits on unordered subsets of the 759 Golay octads "
                f"through size {maximum_size} are enumerated exactly. Setwise-"
                "stabilizer augmentation proves coverage, exact M24 transporters "
                "deduplicate candidates, and every layer passes the independent "
                "orbit-stabilizer mass identity."
            ),
            "not_proved": (
                "The enumeration stops before rank seven. No determinant band, "
                "primitive rank-seven saturation, orthogonal-complement frame, "
                "or ternary transcendental realization is classified. Auxiliaries "
                "not generated by octad glue vectors are outside this prefix engine."
            ),
        },
        "golay_design": {
            "code_parameters": [24, 12, 8],
            "octad_count": len(octads),
            "m24_order": int(coordinate_group.order()),
            "induced_octad_action_order": int(group.order()),
            "induced_action_transitive": True,
            "induced_action_generator_images_sha256": action_digest,
            "fixed_octad_intersection_distribution": {
                "0": 30,
                "2": 448,
                "4": 280,
                "8": 1,
            },
        },
        "method": {
            "augmentation": "setwise stabilizer orbits on the remaining octads",
            "bucket_invariants": [
                "canonical edge-labelled octad-intersection Gram graph",
                "coordinate union size",
                "histogram over all 759 ambient octads of sorted intersection profiles",
            ],
            "final_equivalence": "GAP RepresentativeAction(M24,left,right,OnSets)",
            "independent_completeness_check": (
                "sum over retained representatives of |M24|/|Stab(S)| = binomial(759,k)"
            ),
        },
        "accounting": accounting,
        "layers": layers,
        "next_gate": {
            "description": (
                "Complete every retained five-octad prefix by two octads, using "
                "the rank-seven determinant formula before exact M24 canonicalization."
            ),
            "determinant_band": [1, 500],
            "target_coordinate_union_minimum_for_mw_12_through_17": 19,
            "state": "OPEN",
        },
        "reproduction": {
            "command": (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elkies-k3/scripts/enumerate_24a1_octad_prefix_orbits.sage"
            ),
            "check_command": (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elkies-k3/scripts/enumerate_24a1_octad_prefix_orbits.sage --check"
            ),
        },
    }


parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--maximum-size", type=int, default=5)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
payload = build(arguments.maximum_size)
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"

if arguments.check:
    if not arguments.output.exists() or arguments.output.read_text() != encoded:
        raise SystemExit("24A1 octad-prefix orbit artifact is stale")
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded)

print(
    "OCTADPREFIX|maximum_size={}|orbits={}|mass={}|status=PASS_EXACT".format(
        arguments.maximum_size,
        payload["accounting"][-1]["orbit_count"],
        payload["accounting"][-1]["orbit_stabilizer_mass"],
    ),
    flush=True,
)
