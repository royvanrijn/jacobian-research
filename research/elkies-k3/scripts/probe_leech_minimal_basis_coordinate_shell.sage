#!/usr/bin/env sage-python
"""Exhaust a declared minimal-basis rank-seven shell in the Leech lattice.

Recover the 98,280 antipodal norm-four pairs from the exact Leech backend and
certify that 24 pinned representatives form an ambient unimodular basis.  Test
all ``binomial(24,7)`` coordinate direct summands of that basis.  Quotient the
declared generator language exactly by independent generator signs and
permutations, apply the determinant/discriminant-length gates, construct the
saturated rootless complements, and run the even ternary-genus gate.

The resulting T/NS keys are complete for this finite coordinate language.
This is deliberately a *pre-Co1* probe: it does not identify full Conway
embedding orbits and does not claim completeness for all primitive rank-seven
sublattices of determinant at most 5,000.

status: EXACT_DECLARED_LEECH_COORDINATE_LANGUAGE_PRE_CO1
output: artifacts/generated-results/elkies-k3-leech-minimal-basis-coordinate-shell-v1.json
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import runpy
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from sage.all import Graph, ZZ, matrix, pari


ROOT = Path(__file__).resolve().parents[2]
BACKEND = ROOT / "artifacts/generated-results/elkies-k3-leech-co0-backend-v1.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-leech-minimal-basis-coordinate-shell-v1.json"
)
COMMON_SOURCE = (
    Path(__file__).resolve().parent
    / "enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage"
)
COMMON = runpy.run_path(str(COMMON_SOURCE), run_name="_rank7_leech_probe_common")

rows = COMMON["rows"]
digest = COMMON["digest"]
compact_digest = COMMON["compact_digest"]
discriminant_invariants = COMMON["discriminant_invariants"]
negative_form_key = COMMON["negative_form_key"]
ternary_gate = COMMON["ternary_gate"]


# These are zero-based row indices in the deterministic PARI qfminim(4)
# representative matrix of the pinned backend.  The exact determinant-one
# check below, rather than the indices alone, is the certificate that they
# form an ambient Z-basis.
MINIMAL_VECTOR_BASIS_INDICES_ZERO_BASED = [
    44,
    893,
    78,
    76,
    80,
    79,
    6,
    49,
    48,
    50,
    10,
    45,
    46,
    13,
    14,
    15,
    16,
    17,
    18,
    19,
    20,
    21,
    4601,
    23,
]


def encoded_digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def modular_power(values, exponent, prime):
    result = np.ones(values.shape, dtype=np.int64)
    base = values.astype(np.int64, copy=True) % prime
    while exponent:
        if exponent & 1:
            result = result * base % prime
        base = base * base % prime
        exponent >>= 1
    return result


def exact_small_positive_determinants(grams, prime=65537, hadamard_bound=4**7):
    """Return exact determinants using one prime larger than Hadamard's bound."""

    assert grams.ndim == 3 and grams.shape[1:] == (7, 7)
    values = grams.astype(np.int64, copy=True) % prime
    count = values.shape[0]
    indices = np.arange(count)
    determinants = np.ones(count, dtype=np.int64)
    for column in range(7):
        nonzero = values[:, column:, column] != 0
        assert np.all(np.any(nonzero, axis=1))
        pivot_rows = column + np.argmax(nonzero, axis=1)
        swapping = np.flatnonzero(pivot_rows != column)
        if len(swapping):
            temporary = values[swapping, column, :].copy()
            values[swapping, column, :] = values[
                swapping, pivot_rows[swapping], :
            ]
            values[swapping, pivot_rows[swapping], :] = temporary
            determinants[swapping] = -determinants[swapping]
        pivots = values[:, column, column]
        assert np.all(pivots != 0)
        determinants = determinants * pivots % prime
        if column + 1 < 7:
            inverses = modular_power(pivots, prime - 2, prime)
            factors = values[:, column + 1 :, column] * inverses[:, None] % prime
            values[:, column + 1 :, column:] = (
                values[:, column + 1 :, column:]
                - factors[:, :, None] * values[:, None, column, column:]
            ) % prime
    determinants %= prime
    # Every Gram matrix is positive definite with diagonal entries four, so
    # Hadamard gives 0 < det <= 4^7 < prime.  The residues are therefore the
    # exact integer determinants, not merely modular fingerprints.
    assert np.all((0 < determinants) & (determinants <= hadamard_bound))
    return determinants


def canonical_signed_permutation_key(gram):
    """Exact key modulo generator permutations and independent sign changes."""

    assert gram.nrows() == gram.ncols() == 7
    graph = Graph(14)
    signs = (1, -1)
    for left in range(7):
        # Label 3 is outside the possible inner-product range [-2,2] and
        # intrinsically marks the two signs of one generator.
        graph.add_edge(2 * left, 2 * left + 1, 3)
        for right in range(left + 1, 7):
            pairing = int(gram[left, right])
            assert -2 <= pairing <= 2
            for left_sign_index, left_sign in enumerate(signs):
                for right_sign_index, right_sign in enumerate(signs):
                    graph.add_edge(
                        2 * left + left_sign_index,
                        2 * right + right_sign_index,
                        left_sign * right_sign * pairing,
                    )
    assert graph.size() == 91
    canonical = graph.canonical_label(edge_labels=True)
    return tuple(
        int(label)
        for unused_left, unused_right, label in canonical.edges(
            sort=True, labels=True
        )
    )


def counter_rows(counter):
    return {str(key): value for key, value in sorted(counter.items())}


def build(backend, determinant_bound, discriminant_length_bound):
    assert backend["schema"] == "elkies-k3.leech-co0-backend.v1"
    assert backend["status"] == (
        "PASS_EXACT_LEECH_GRAM_AND_CO0_ATLAS_ACTION_BACKEND_FOUNDATION"
    )
    assert backend["leech_lattice"]["minimum_squared_norm"] == 4
    assert backend["leech_lattice"]["roots"] == 0
    ambient = matrix(ZZ, backend["leech_lattice"]["gram"])
    minimal_representatives = matrix(
        ZZ, pari(ambient).qfminim(4)[2].sage()
    ).transpose()
    assert minimal_representatives.nrows() == 98280
    assert all(
        int(vector * ambient * vector) == 4
        for vector in minimal_representatives.rows()
    )
    basis = matrix(
        ZZ,
        [
            minimal_representatives[index]
            for index in MINIMAL_VECTOR_BASIS_INDICES_ZERO_BASED
        ],
    )
    assert abs(int(basis.det())) == 1
    basis_gram = basis * ambient * basis.transpose()
    assert basis_gram.det() == 1
    assert set(map(int, basis_gram.diagonal())) == {4}

    combinations = np.array(
        list(itertools.combinations(range(24), 7)), dtype=np.uint8
    )
    assert combinations.shape == (346104, 7)
    basis_gram_numpy = np.array(basis_gram, dtype=np.int16)
    all_grams = basis_gram_numpy[
        combinations[:, :, None], combinations[:, None, :]
    ]
    flattened = all_grams.reshape(len(combinations), 49)
    unique_flattened, first_indices, literal_inverse, literal_counts = np.unique(
        flattened,
        axis=0,
        return_index=True,
        return_inverse=True,
        return_counts=True,
    )
    assert len(unique_flattened) == 129064
    unique_grams = unique_flattened.reshape(-1, 7, 7)
    determinants = exact_small_positive_determinants(unique_grams)
    assert int(determinants.min()) == 486
    assert int(determinants.max()) == 1984
    assert np.all(determinants <= determinant_bound)

    signed_types = {}
    literal_signed_keys = []
    for literal_index, (gram_numpy, first_index, multiplicity, determinant) in enumerate(
        zip(unique_grams, first_indices, literal_counts, determinants)
    ):
        gram = matrix(ZZ, gram_numpy.tolist())
        assert int(gram.det()) == int(determinant)
        key = canonical_signed_permutation_key(gram)
        literal_signed_keys.append(key)
        combination = tuple(map(int, combinations[int(first_index)]))
        if key not in signed_types:
            signed_types[key] = {
                "representative_subset": combination,
                "representative_gram": gram,
                "determinant": int(determinant),
                "coordinate_embedding_multiplicity": 0,
                "literal_ordered_gram_patterns": 0,
            }
        row = signed_types[key]
        assert row["determinant"] == int(determinant)
        row["coordinate_embedding_multiplicity"] += int(multiplicity)
        row["literal_ordered_gram_patterns"] += 1
    assert len(signed_types) == 221
    assert sum(
        row["coordinate_embedding_multiplicity"] for row in signed_types.values()
    ) == len(combinations)
    assert sum(
        row["literal_ordered_gram_patterns"] for row in signed_types.values()
    ) == len(unique_grams)

    # Compute the complete norm-four spectrum of every complement without
    # invoking 346,104 separate 17-dimensional short-vector enumerations.
    # For each ambient minimal line, record the 24 basis vectors to which it
    # is orthogonal. A superset zeta transform then answers every seven-subset
    # query exactly.
    minimal_numpy = np.array(minimal_representatives, dtype=np.int64)
    ambient_numpy = np.array(ambient, dtype=np.int64)
    basis_numpy = np.array(basis, dtype=np.int64)
    pairings = minimal_numpy @ ambient_numpy @ basis_numpy.transpose()
    assert pairings.shape == (98280, 24)
    bit_weights = np.left_shift(
        np.uint32(1), np.arange(24, dtype=np.uint32)
    )
    zero_masks = np.sum(
        (pairings == 0).astype(np.uint32) * bit_weights[None, :],
        axis=1,
        dtype=np.uint32,
    )
    orthogonality_frequency = np.bincount(
        zero_masks.astype(np.int64), minlength=1 << 24
    ).astype(np.int32)
    assert int(orthogonality_frequency.sum()) == 98280
    superset_counts = orthogonality_frequency.copy()
    for bit in range(24):
        step = 1 << bit
        blocks = superset_counts.reshape(-1, 2 * step)
        blocks[:, :step] += blocks[:, step:]
    assert int(superset_counts[0]) == 98280
    subset_masks = np.sum(
        bit_weights[combinations], axis=1, dtype=np.uint32
    )
    complement_norm_four_pairs = superset_counts[subset_masks]
    assert np.all(complement_norm_four_pairs > 0)

    sorted_signed_keys = sorted(signed_types)
    signed_index = {
        key: index for index, key in enumerate(sorted_signed_keys)
    }
    literal_type_indices = np.array(
        [signed_index[key] for key in literal_signed_keys], dtype=np.uint16
    )
    coordinate_type_indices = literal_type_indices[literal_inverse]
    complement_spectra = [Counter() for unused in sorted_signed_keys]
    overall_complement_spectrum = Counter()
    for type_index, norm_four_pairs in zip(
        coordinate_type_indices, complement_norm_four_pairs
    ):
        value = int(norm_four_pairs)
        complement_spectra[int(type_index)][value] += 1
        overall_complement_spectrum[value] += 1
    assert sum(map(sum, (row.values() for row in complement_spectra))) == len(
        combinations
    )
    assert sum(overall_complement_spectrum.values()) == len(combinations)

    genera_cache = {}
    form_cache = {}
    length_distribution = Counter()
    weighted_length_distribution = Counter()
    determinant_type_distribution = Counter()
    determinant_embedding_distribution = Counter()
    accepted_determinant_distribution = Counter()
    type_records = []
    surface_groups = defaultdict(list)
    surface_key_payloads = {}
    sorted_types = [
        (key, signed_types[key]) for key in sorted_signed_keys
    ]
    for type_index, (canonical_key, source) in enumerate(sorted_types, start=1):
        type_id = f"LEECH-MB-C{type_index:03d}"
        subset = source["representative_subset"]
        auxiliary_basis = matrix(ZZ, [basis[index] for index in subset])
        auxiliary = auxiliary_basis * ambient * auxiliary_basis.transpose()
        assert auxiliary == source["representative_gram"]
        determinant = int(auxiliary.det())
        invariants = discriminant_invariants(auxiliary)
        length = len(invariants)
        length_distribution[length] += 1
        weighted_length_distribution[length] += source[
            "coordinate_embedding_multiplicity"
        ]
        determinant_type_distribution[determinant] += 1
        determinant_embedding_distribution[determinant] += source[
            "coordinate_embedding_multiplicity"
        ]
        record = {
            "type_id": type_id,
            "canonical_signed_permutation_key": list(canonical_key),
            "canonical_signed_permutation_key_sha256": encoded_digest(
                list(canonical_key)
            ),
            "coordinate_embedding_multiplicity": source[
                "coordinate_embedding_multiplicity"
            ],
            "literal_ordered_gram_patterns": source[
                "literal_ordered_gram_patterns"
            ],
            "representative_coordinate_subset_zero_based": list(subset),
            "representative_auxiliary_basis_in_ambient": rows(auxiliary_basis),
            "auxiliary_gram": rows(auxiliary),
            "determinant": determinant,
            "discriminant_invariants_greater_than_one": invariants,
            "discriminant_length": length,
            "mw_rank_for_rho_19": 17,
            "rootless_complement_norm_four_pair_distribution": counter_rows(
                complement_spectra[type_index - 1]
            ),
        }
        if length > discriminant_length_bound:
            record["gate"] = {
                "status": "REJECT_DISCRIMINANT_LENGTH",
                "matching_even_ternary_genera": 0,
            }
            type_records.append(record)
            continue

        # A subset of a Z-basis is a primitive direct summand.  The right
        # kernel is consequently the saturated integral complement.
        complement_basis = (auxiliary_basis * ambient).right_kernel_matrix()
        assert complement_basis.nrows() == 17
        frame = complement_basis * ambient * complement_basis.transpose()
        assert int(frame.det()) == determinant
        gate = ternary_gate(frame, genera_cache, form_cache)
        matching = gate["matching_even_ternary_genera"]
        assert matching in (0, 1)
        record["representative_complement_basis_in_ambient"] = rows(
            complement_basis
        )
        record["representative_rootless_frame_gram"] = rows(frame)
        record["gate"] = {
            "status": (
                "PASS_EXACT_ONE_TERNARY_GENUS"
                if matching == 1
                else "REJECT_NO_TERNARY_GENUS"
            ),
            **gate,
        }
        if matching == 1:
            accepted_determinant_distribution[determinant] += 1
            ternary = gate["ternary_genus_representatives"][0]
            # Use the literal global-catalogue key schema so cross-backend
            # comparisons are byte-exact rather than a later relabelling.
            surface_key_data = {
                "ns_discriminant_form_key": negative_form_key(
                    gate["frame_discriminant_form_normal_key"]
                ),
                "transcendental_gram": ternary["gram"],
            }
            surface_key = compact_digest(surface_key_data)
            if surface_key in surface_key_payloads:
                assert surface_key_payloads[surface_key] == surface_key_data
            else:
                surface_key_payloads[surface_key] = surface_key_data
            record["preliminary_surface_id"] = f"K3-{surface_key}"
            surface_groups[surface_key].append(type_id)
        type_records.append(record)

    assert counter_rows(length_distribution) == {
        "1": 94,
        "2": 30,
        "3": 83,
        "4": 1,
        "5": 12,
        "7": 1,
    }
    length_admissible = sum(
        count
        for length, count in length_distribution.items()
        if length <= discriminant_length_bound
    )
    ternary_accepted = sum(
        record["gate"]["status"] == "PASS_EXACT_ONE_TERNARY_GENUS"
        for record in type_records
    )
    assert length_admissible == 207
    assert ternary_accepted == 194
    assert all(
        record["gate"]["matching_even_ternary_genera"] == 1
        for record in type_records
        if record["gate"]["status"] == "PASS_EXACT_ONE_TERNARY_GENUS"
    )

    surface_records = []
    by_type_id = {record["type_id"]: record for record in type_records}
    for surface_key, type_ids in sorted(surface_groups.items()):
        representatives = [by_type_id[type_id] for type_id in type_ids]
        determinant_values = {row["determinant"] for row in representatives}
        assert len(determinant_values) == 1
        first = representatives[0]
        surface_spectrum = Counter()
        for row in representatives:
            surface_spectrum.update(
                {
                    int(norm_four_pairs): multiplicity
                    for norm_four_pairs, multiplicity in row[
                        "rootless_complement_norm_four_pair_distribution"
                    ].items()
                }
            )
        surface_records.append(
            {
                "preliminary_surface_id": f"LEECH-PRE-{surface_key}",
                "catalogue_surface_id_if_imported": f"K3-{surface_key}",
                "surface_key": surface_key_payloads[surface_key],
                "surface_key_sha256_prefix": surface_key,
                "determinant": first["determinant"],
                "ternary_gram": first["gate"][
                    "ternary_genus_representatives"
                ][0]["gram"],
                "ns_discriminant_form": negative_form_key(
                    first["gate"]["frame_discriminant_form_normal_key"]
                ),
                "signed_basis_type_ids": type_ids,
                "coordinate_embedding_multiplicity": sum(
                    row["coordinate_embedding_multiplicity"]
                    for row in representatives
                ),
                "mw_rank_for_rho_19": 17,
                "rootless_complement_norm_four_pair_distribution": counter_rows(
                    surface_spectrum
                ),
                "rootless_complement_norm_four_pairs_minimum": min(
                    surface_spectrum
                ),
                "rootless_complement_norm_four_pairs_maximum": max(
                    surface_spectrum
                ),
                "orbit_status": "PRE_CO1_EMBEDDING_QUOTIENT",
            }
        )

    return {
        "schema": "elkies-k3.leech-minimal-basis-coordinate-shell.v1",
        "status": "PASS_EXACT_DECLARED_LEECH_COORDINATE_LANGUAGE_PRE_CO1",
        "proof_scope": {
            "proved": (
                "The pinned 24 norm-four representatives form a determinant-one "
                "ambient basis. All 346,104 rank-seven coordinate direct summands "
                "are therefore primitive and are exhausted. Their determinants "
                "are computed exactly, their generator bases are quotiented exactly "
                "by signed permutations, every complement is saturated and rootless, "
                "and every length-admissible type passes through the exact even "
                "ternary-genus and T/NS-first gates."
            ),
            "not_proved": (
                "The 221 signed-basis types are not full Co1 embedding orbits. This "
                "coordinate language is not all primitive rank-seven Leech "
                "sublattices through determinant 5,000, and its preliminary T/NS "
                "records are not imported into the global catalogue before the "
                "Conway quotient and cross-backend provenance gates are completed."
            ),
        },
        "parameters": {
            "ambient_label": "Leech",
            "determinant_bound": determinant_bound,
            "discriminant_length_bound": discriminant_length_bound,
            "minimum_mw_rank": 12,
            "maximum_mw_rank": 17,
            "seed_language": (
                "all 7-of-24 coordinate direct summands of one certified "
                "norm-four determinant-one ambient basis"
            ),
            "pre_ambient_orbit_equivalence": (
                "independent generator sign changes and generator permutations"
            ),
            "ambient_orbit_group_remaining": "Co1 on antipodal lines (Co0 on vectors)",
        },
        "minimal_vector_basis": {
            "qfminim_representative_indices_zero_based": (
                MINIMAL_VECTOR_BASIS_INDICES_ZERO_BASED
            ),
            "basis_in_pinned_ambient_coordinates": rows(basis),
            "gram": rows(basis_gram),
            "basis_matrix_determinant": int(basis.det()),
            "gram_determinant": int(basis_gram.det()),
            "all_basis_vector_norms": [
                int(basis_gram[index, index]) for index in range(24)
            ],
            "primitive_coordinate_summand_reason": (
                "every subset of a Z-basis spans a primitive direct summand"
            ),
        },
        "accounting": {
            "coordinate_subsets_tested": len(combinations),
            "determinant_minimum": int(determinants.min()),
            "determinant_maximum": int(determinants.max()),
            "distinct_determinants": len(determinant_embedding_distribution),
            "determinant_rejected": 0,
            "literal_ordered_gram_patterns": len(unique_grams),
            "signed_permutation_basis_types": len(type_records),
            "discriminant_length_type_distribution": counter_rows(
                length_distribution
            ),
            "discriminant_length_coordinate_embedding_distribution": counter_rows(
                weighted_length_distribution
            ),
            "length_admissible_signed_basis_types": length_admissible,
            "ternary_compatible_signed_basis_types": ternary_accepted,
            "preliminary_T_NS_surface_keys": len(surface_records),
            "mw_rank_distribution_for_all_coordinate_embeddings": {
                "17": len(combinations)
            },
            "rootless_complement_norm_four_pair_distribution": counter_rows(
                overall_complement_spectrum
            ),
            "rootless_complement_norm_four_pairs_minimum": min(
                overall_complement_spectrum
            ),
            "rootless_complement_norm_four_pairs_maximum": max(
                overall_complement_spectrum
            ),
            "determinant_signed_basis_type_distribution": counter_rows(
                determinant_type_distribution
            ),
            "determinant_coordinate_embedding_distribution": counter_rows(
                determinant_embedding_distribution
            ),
            "ternary_accepted_signed_basis_type_determinant_distribution": counter_rows(
                accepted_determinant_distribution
            ),
        },
        "signed_basis_types": type_records,
        "preliminary_T_NS_surfaces": surface_records,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", type=Path, default=BACKEND)
    parser.add_argument("--determinant-bound", type=int, default=5000)
    parser.add_argument("--discriminant-length-bound", type=int, default=3)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    assert arguments.determinant_bound == 5000
    assert arguments.discriminant_length_bound == 3
    payload = build(
        json.loads(arguments.backend.read_text()),
        arguments.determinant_bound,
        arguments.discriminant_length_bound,
    )
    payload["inputs"] = {
        str(arguments.backend.resolve().relative_to(ROOT)): digest(
            arguments.backend
        ),
        str(COMMON_SOURCE.resolve().relative_to(ROOT)): digest(COMMON_SOURCE),
    }
    payload["reproduce"] = (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/probe_leech_minimal_basis_coordinate_shell.sage"
    )
    payload["check_command"] = payload["reproduce"] + " --check"
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output = arguments.output.resolve()
    if arguments.check:
        if not output.exists() or output.read_text() != encoded:
            raise SystemExit("Leech minimal-basis coordinate-shell artifact is stale")
    else:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded)
    print(
        "LEECHMBSHELL|subsets={}|signed_types={}|ternary={}|T_NS={}|"
        "status=PASS_EXACT_PRE_CO1".format(
            payload["accounting"]["coordinate_subsets_tested"],
            payload["accounting"]["signed_permutation_basis_types"],
            payload["accounting"]["ternary_compatible_signed_basis_types"],
            payload["accounting"]["preliminary_T_NS_surface_keys"],
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
