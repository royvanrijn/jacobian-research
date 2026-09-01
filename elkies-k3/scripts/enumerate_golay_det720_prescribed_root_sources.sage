#!/usr/bin/env sage-python
"""Search the Golay determinant-720 auxiliary directly for low-MW sources.

The older prescribed-root foundry enumerator starts from a D5 root subsystem
inside each auxiliary.  The Golay-octad auxiliary is rootless, so that anchor
does not exist.  This script instead embeds an ordered norm-four basis of the
fixed determinant-720 auxiliary one vector at a time in every selected rooted
Niemeier lattice.  After each vector, the residual Weyl group is used to put
the next vector in its dominant chamber.  At the seventh vector, zero Dynkin
labels prescribe a rank-15 or rank-16 complement root system before the small
fixed-space ellipsoid is solved.

Every retained row is checked in the full integral Niemeier lattice: the
auxiliary embedding is primitive, its saturated orthogonal complement has the
declared complete root system, and its determinant and genus agree with the
pinned Golay target.  Equal deterministic reduced Grams are merged.  Distinct
reduced Grams are not asserted to be distinct integral-isometry classes.

For primitive-root MW1/MW2 sources, complete frame shells through norm eight
also audit whether a Mordell--Weil basis can be represented by physical
sections of pole order at most two.  This is a lattice/source certificate; it
does not construct a rational marking, an equation, or a neighbour corridor.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from pathlib import Path

from fpylll import Enumeration, FPLLL, GSO, IntegerMatrix
from sage.all import Genus, QQ, ZZ, matrix, pari, vector


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TARGET = ROOT / "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json"
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json"
)


_shared_path = HERE / "enumerate_lattice_foundry_prescribed_root_sources.sage"
_shared = {"__file__": str(_shared_path), "__name__": "golay_prescribed_shared"}
exec(compile(_shared_path.read_text(), str(_shared_path), "exec"), _shared)

rows = _shared["rows"]
rational_rows = _shared["rational_rows"]
gram_digest = _shared["gram_digest"]
reduced_gram = _shared["reduced_gram"]
primitive_rows = _shared["primitive_rows"]
signed_roots = _shared["signed_roots"]
residual_simple_data = _shared["residual_simple_data"]
dominant_labels_up_to = _shared["dominant_labels_up_to"]
prescribed_label_rows = _shared["prescribed_label_rows"]
coordinate_model = _shared["coordinate_model"]
shifted_ellipsoid_shell = _shared["shifted_ellipsoid_shell"]
frame_root_record = _shared["frame_root_record"]
minimize_child_frame = _shared["minimize_child_frame"]
cartan_components = _shared["cartan_components"]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def label_combinations(simple_components, gram, bound):
    """Enumerate componentwise dominant labels without a large Cartesian product."""
    choices = [
        dominant_labels_up_to(
            (component * gram * component.transpose()).inverse(), bound
        )
        for component in simple_components
    ]
    result = []
    labels = []

    def extend(index, norm):
        if index == len(choices):
            result.append(
                {
                    "labels": tuple(entry for block in labels for entry in block),
                    "label_norm": norm,
                }
            )
            return
        for block, block_norm in choices[index]:
            if norm + block_norm > bound:
                continue
            labels.append(block)
            extend(index + 1, norm + block_norm)
            labels.pop()

    extend(0, QQ(0))
    return result


def enumerate_extensions(
    ambient,
    fixed,
    ambient_roots,
    pairings,
    norm,
    arguments,
    final=False,
):
    unused_roots, simple_components, simple, cartan = residual_simple_data(
        ambient, fixed, ambient_roots
    )
    residual_rank = simple.nrows()
    if residual_rank < arguments.source_root_rank_min:
        return [], {
            "residual_root_rank": residual_rank,
            "label_rows": 0,
            "ellipsoid_solutions": 0,
        }
    support_capacity = sum(
        sorted(
            (component.nrows() for component in simple_components),
            reverse=True,
        )[: arguments.source_support_max]
    )
    if support_capacity < arguments.source_root_rank_min:
        return [], {
            "residual_root_rank": residual_rank,
            "support_capacity": support_capacity,
            "label_rows": 0,
            "support_capacity_pruned": True,
            "ellipsoid_solutions": 0,
        }

    fixed_rank = fixed.nrows()
    pairings = vector(QQ, pairings)
    if fixed_rank:
        fixed_gram = fixed * ambient * fixed.transpose()
        base_norm = pairings * fixed_gram.inverse() * pairings
    else:
        base_norm = QQ(0)
    label_budget = QQ(norm) - base_norm
    if label_budget < 0:
        return [], {
            "residual_root_rank": residual_rank,
            "label_rows": 0,
            "ellipsoid_solutions": 0,
        }

    if final:
        label_rows, total_label_rows = prescribed_label_rows(
            cartan,
            label_budget,
            arguments.source_root_rank_min,
            arguments.source_root_rank_max,
            arguments.source_support_min,
            arguments.source_support_max,
            arguments.all_a_only,
            arguments.final_label_limit,
        )
    else:
        label_rows = label_combinations(
            simple_components, ambient, label_budget
        )
        total_label_rows = len(label_rows)

    inverse, coordinate_norm, fixed_quadratic, extra_coordinates = coordinate_model(
        ambient, fixed, simple
    )
    fixed_dimension = 24 - fixed_rank - residual_rank
    assert fixed_quadratic.nrows() == fixed_dimension
    if fixed_dimension > arguments.fixed_dimension_max:
        return [], {
            "residual_root_rank": residual_rank,
            "fixed_dimension": fixed_dimension,
            "label_rows": total_label_rows,
            "skipped_fixed_dimension": True,
            "ellipsoid_solutions": 0,
        }

    fixed_inverse = (
        fixed_quadratic.inverse() if fixed_dimension else matrix(QQ, 0, 0)
    )
    prefix_size = fixed_rank + residual_rank
    candidates = {}
    ellipsoid_count = 0
    for label_row in label_rows:
        labels = vector(QQ, label_row["labels"])
        prefix = vector(QQ, list(pairings) + list(labels))
        prefix_norm = (
            prefix * coordinate_norm[:prefix_size, :prefix_size] * prefix
        )
        if fixed_dimension:
            linear = prefix * coordinate_norm[:prefix_size, prefix_size:]
            centre = -linear * fixed_inverse
            minimum_norm = prefix_norm - linear * fixed_inverse * linear
        else:
            centre = vector(QQ, [])
            minimum_norm = prefix_norm
        assert minimum_norm == base_norm + label_row["label_norm"]
        fixed_solutions = shifted_ellipsoid_shell(
            fixed_quadratic, centre, QQ(norm) - minimum_norm
        )
        ellipsoid_count += len(fixed_solutions)
        for fixed_coordinates in fixed_solutions:
            coordinates = vector(
                QQ, list(pairings) + list(labels) + list(fixed_coordinates)
            )
            extension = coordinates * inverse
            if not all(entry.denominator() == 1 for entry in extension):
                continue
            extension = vector(ZZ, extension)
            enlarged = fixed.stack(matrix(ZZ, [extension]))
            expected = (fixed * ambient * fixed.transpose()).augment(
                matrix(ZZ, fixed_rank, 1, list(pairings))
            ) if fixed_rank else matrix(ZZ, 0, 1)
            if extension * ambient * extension != norm:
                continue
            if fixed_rank and extension * ambient * fixed.transpose() != pairings:
                continue
            key = tuple(map(int, extension))
            candidates[key] = {
                "vector": extension,
                "labels": tuple(map(int, labels)),
                "fixed_coordinates": tuple(map(int, fixed_coordinates)),
                "fixed_dimension": fixed_dimension,
                "extra_coordinates": extra_coordinates,
                "prescribed_root_type": label_row.get("prescribed_root_type"),
                "prescribed_root_rank": label_row.get("prescribed_root_rank"),
                "prescribed_support_count": label_row.get(
                    "prescribed_support_count"
                ),
            }
    ordered = [candidates[key] for key in sorted(candidates)]
    if arguments.extension_limit:
        ordered = ordered[: arguments.extension_limit]
    return ordered, {
        "residual_root_rank": residual_rank,
        "fixed_dimension": fixed_dimension,
        "label_rows": total_label_rows,
        "retained_label_rows": len(label_rows),
        "ellipsoid_solutions": ellipsoid_count,
        "integral_extensions": len(candidates),
        "retained_extensions": len(ordered),
    }


def make_gso(root_gram, float_type, precision):
    if float_type == "mpfr":
        FPLLL.set_precision(precision)
    integer_matrix = IntegerMatrix.from_matrix(
        [[int(entry) for entry in row] for row in root_gram.rows()]
    )
    gso = GSO.Mat(
        integer_matrix, gram=True, float_type=float_type, update=True
    )
    rank = root_gram.nrows()
    mu = [
        [gso.get_mu(i, j) if i > j else 0.0 for j in range(rank)]
        for i in range(rank)
    ]
    return gso, mu


def affine_closest(frame, root_gram, target, tail, gso, mu):
    rank = root_gram.nrows()
    target_gso = [
        float(target[i])
        + sum(float(target[j]) * mu[j][i] for j in range(i + 1, rank))
        for i in range(rank)
    ]
    zero_distance = target * root_gram * target
    solutions = Enumeration(gso).enumerate(
        0,
        rank,
        float(zero_distance) + 1.0,
        0,
        target=target_gso,
    )
    if not solutions:
        raise RuntimeError("Golay source affine root CVP returned no solution")
    reported_distance, coordinates = solutions[0]
    root_coordinates = [int(round(value)) for value in coordinates]
    if any(
        abs(value - integer) > 1e-7
        for value, integer in zip(coordinates, root_coordinates)
    ):
        raise RuntimeError("Golay source affine CVP coordinates are not integral")
    section = vector(ZZ, root_coordinates + list(map(int, tail)))
    norm = int(section * frame * section)
    exact_distance = (
        (vector(QQ, root_coordinates) - target)
        * root_gram
        * (vector(QQ, root_coordinates) - target)
    )
    return section, norm, abs(float(reported_distance) - float(exact_distance))


def dominant_section(section, frame, root_rank):
    section = vector(ZZ, section)
    for unused in range(10000):
        pairings = section * frame
        negative = next(
            (index for index in range(root_rank) if pairings[index] < 0), None
        )
        if negative is None:
            return section
        reflected = vector(ZZ, section)
        reflected[negative] -= pairings[negative]
        assert reflected * frame * reflected == section * frame * section
        section = reflected
    raise RuntimeError("Golay source root-Weyl dominance did not terminate")


def physical_sections_through_pole_two(
    root_adapted, root_rank, components, height
):
    """Return exact affine-CVP section representatives through frame norm eight."""
    mw_rank = root_adapted.nrows() - root_rank
    denominator = math.lcm(
        *(entry.denominator() for entry in height.list())
    )
    integral_height = (denominator * height).change_ring(ZZ)
    quotient_shell = pari(integral_height).qfminim(8 * denominator)
    half_tails = matrix(ZZ, quotient_shell[2].sage()).columns()
    tails = {
        tuple(map(int, sign * tail))
        for tail in half_tails
        for sign in (1, -1)
        if any(tail)
    }

    root_gram = matrix(QQ, root_adapted[:root_rank, :root_rank])
    cross = matrix(QQ, root_adapted[:root_rank, root_rank:])
    root_inverse = root_gram.inverse()
    primary_gso, primary_mu = make_gso(root_gram, "dd", 0)
    audit_gso, audit_mu = make_gso(root_gram, "mpfr", 256)
    candidates = {}
    component_indices = []
    cursor = 0
    for component in components:
        indices = tuple(range(cursor, cursor + component["rank"]))
        component_indices.append(indices)
        cursor += component["rank"]
    assert cursor == root_rank
    maximum_error = 0.0
    for tail_tuple in sorted(tails):
        tail = vector(QQ, tail_tuple)
        target = -root_inverse * cross * tail
        primary, primary_norm, primary_error = affine_closest(
            root_adapted,
            root_gram,
            target,
            tail_tuple,
            primary_gso,
            primary_mu,
        )
        audit, audit_norm, audit_error = affine_closest(
            root_adapted,
            root_gram,
            target,
            tail_tuple,
            audit_gso,
            audit_mu,
        )
        if primary_norm != audit_norm:
            raise ValueError("cross-precision Golay source affine CVP mismatch")
        maximum_error = max(maximum_error, primary_error, audit_error)
        if primary_norm not in (4, 6, 8):
            continue
        section = dominant_section(primary, root_adapted, root_rank)
        pairings = section * root_adapted
        root_pairings = tuple(map(int, pairings[:root_rank]))
        if any(value not in (0, 1) for value in root_pairings):
            continue
        if any(
            sum(root_pairings[index] for index in indices) > 1
            for indices in component_indices
        ):
            continue
        quotient = tuple(map(int, section[root_rank:]))
        key = (primary_norm, quotient, tuple(map(int, section)))
        candidates[key] = {
            "frame_vector": list(map(int, section)),
            "mw_quotient_coordinates": list(quotient),
            "frame_norm": primary_norm,
            "pole_order": (primary_norm - 4) // 2,
            "simple_root_pairings": list(root_pairings),
            "maximum_cross_precision_cvp_distance_error": maximum_error,
        }
    return [candidates[key] for key in sorted(candidates)]


def pole_audit(minimized, components, mw_rank):
    if mw_rank not in (1, 2):
        return {
            "status": "NOT_APPLICABLE_OUTSIDE_MW1_MW2_SEARCH_BAND",
            "minimum_nonzero_section_pole_order": None,
            "basis_with_all_poles_at_most_two": None,
        }
    if not minimized["root_lattice_primitive"]:
        return {
            "status": "OPEN_NONPRIMITIVE_ROOT_GLUE_ANALYSIS_REQUIRED",
            "minimum_nonzero_section_pole_order": None,
            "basis_with_all_poles_at_most_two": None,
        }
    root_rank = 17 - mw_rank
    frame = minimized["frame"]
    sections = physical_sections_through_pole_two(
        frame, root_rank, components, minimized["mw_height"]
    )
    if not sections:
        return {
            "status": "OPEN_NO_PHYSICAL_SECTION_THROUGH_POLE_TWO",
            "minimum_nonzero_section_pole_order": None,
            "basis_with_all_poles_at_most_two": False,
            "physical_sections_through_pole_two": 0,
        }
    basis = None
    if mw_rank == 1:
        basis = next(
            (
                [section]
                for section in sections
                if abs(section["mw_quotient_coordinates"][0]) == 1
            ),
            None,
        )
    elif mw_rank == 2:
        for left, right in itertools.combinations(sections, 2):
            quotient = matrix(
                ZZ,
                [
                    left["mw_quotient_coordinates"],
                    right["mw_quotient_coordinates"],
                ],
            )
            if abs(quotient.det()) == 1:
                basis = [left, right]
                break
    minimum_pole = min(section["pole_order"] for section in sections)
    return {
        "status": (
            "PASS_EXACT_PHYSICAL_MW_BASIS_THROUGH_POLE_TWO"
            if basis is not None
            else "OPEN_NO_UNIMODULAR_PHYSICAL_BASIS_THROUGH_POLE_TWO"
        ),
        "minimum_nonzero_section_pole_order": minimum_pole,
        "basis_with_all_poles_at_most_two": basis is not None,
        "physical_sections_through_pole_two": len(sections),
        "basis": basis,
        "complete_shell_through_frame_norm": 8,
    }


def source_record(frame):
    label, root_rank, root_count, root_det, components = frame_root_record(frame)
    minimized = minimize_child_frame(frame)
    mw_rank = 17 - root_rank
    height = minimized["mw_height"]
    return {
        "gram": rows(frame),
        "gram_sha256": gram_digest(frame),
        "determinant": int(frame.det()),
        "root_type": label,
        "root_components": components,
        "root_rank": root_rank,
        "signed_root_count": root_count,
        "root_determinant": root_det,
        "support_count": len(components),
        "mw_rank_for_rho_19": mw_rank,
        "root_lattice_primitive": minimized["root_lattice_primitive"],
        "root_smith_invariants": list(
            map(int, minimized["root_smith_invariants"])
        ),
        "root_adapted_gram": rows(minimized["frame"]),
        "root_adapted_basis_rows_in_source_basis": rows(minimized["basis"]),
        "mw_height_gram": None if height is None else rational_rows(height),
        "mw_regulator": None if height is None else str(abs(height.det())),
        "torsion": (
            1
            if minimized["root_lattice_primitive"]
            else "REQUIRES_GLUE_ANALYSIS"
        ),
        "pole_audit": pole_audit(minimized, components, mw_rank),
    }


def search_ambient(label, ambient, auxiliary, ambient_roots, target_genus, arguments):
    counters = Counter()
    depth_records = {str(depth): Counter() for depth in range(1, 8)}
    sources = {}

    def descend(fixed, depth):
        if depth == 7:
            counters["complete_auxiliary_embeddings"] += 1
            if not primitive_rows(fixed):
                counters["rejected_nonprimitive_auxiliary"] += 1
                return
            complement_basis = (fixed * ambient).right_kernel_matrix()
            frame = complement_basis * ambient * complement_basis.transpose()
            if frame.det() != auxiliary.det() or Genus(frame) != target_genus:
                counters["rejected_target_genus"] += 1
                return
            reduced, change = reduced_gram(frame)
            key = gram_digest(reduced)
            if key in sources:
                sources[key]["embedding_count_merged"] += 1
                counters["retained_embeddings"] += 1
                return
            source = source_record(frame)
            if not (
                arguments.source_root_rank_min
                <= source["root_rank"]
                <= arguments.source_root_rank_max
            ):
                counters["rejected_final_root_rank"] += 1
                return
            if not (
                arguments.source_support_min
                <= source["support_count"]
                <= arguments.source_support_max
            ):
                counters["rejected_final_support_count"] += 1
                return
            if arguments.all_a_only and any(
                not component["type"].startswith("A")
                for component in source["root_components"]
            ):
                counters["rejected_non_a_component"] += 1
                return
            sources[key] = {
                "source": source,
                "reduced_gram": rows(reduced),
                "reduced_gram_sha256": key,
                "source_to_reduced_basis": rows(change),
                "auxiliary_basis_in_ambient": rows(fixed),
                "complement_basis_in_ambient": rows(complement_basis),
                "ambient_label": label,
                "embedding_count_merged": 0,
            }
            sources[key]["embedding_count_merged"] += 1
            counters["retained_embeddings"] += 1
            return

        pairings = list(auxiliary.row(depth)[:depth])
        norm = auxiliary[depth, depth]
        extensions, accounting = enumerate_extensions(
            ambient,
            fixed,
            ambient_roots,
            pairings,
            norm,
            arguments,
            final=(depth == 6),
        )
        record = depth_records[str(depth + 1)]
        record["prefixes"] += 1
        record["label_rows"] += accounting.get("label_rows", 0)
        record["ellipsoid_solutions"] += accounting.get(
            "ellipsoid_solutions", 0
        )
        record["integral_extensions"] += accounting.get(
            "integral_extensions", 0
        )
        record["retained_extensions"] += len(extensions)
        record["maximum_fixed_dimension"] = max(
            record["maximum_fixed_dimension"],
            accounting.get("fixed_dimension", 0),
        )
        if accounting.get("skipped_fixed_dimension"):
            record["skipped_fixed_dimension_prefixes"] += 1
        if accounting.get("support_capacity_pruned"):
            record["support_capacity_pruned_prefixes"] += 1
        for extension in extensions:
            enlarged = fixed.stack(matrix(ZZ, [extension["vector"]]))
            if enlarged * ambient * enlarged.transpose() != auxiliary[: depth + 1, : depth + 1]:
                raise AssertionError("Golay auxiliary prefix Gram mismatch")
            descend(enlarged, depth + 1)

    empty = matrix(ZZ, 0, 24)
    descend(empty, 0)
    ordered = [sources[key] for key in sorted(sources)]
    return ordered, {
        "ambient_label": label,
        "ambient_signed_roots": len(ambient_roots),
        "depth_accounting": {
            depth: dict(record) for depth, record in depth_records.items()
        },
        "totals": dict(counters),
        "distinct_reduced_gram_sources": len(ordered),
    }


def build(arguments):
    target_payload = json.loads(arguments.target.read_text())
    catalog_payload = json.loads(arguments.catalog.read_text())
    auxiliary = matrix(
        ZZ, target_payload["support_design"]["raw_octad_intersection_gram"]
    )
    assert auxiliary.nrows() == 7 and auxiliary.det() == 720
    target_frame = matrix(ZZ, target_payload["frame"]["gram"])
    target_genus = Genus(target_frame)

    requested = set(arguments.ambient_label)
    ambient_rows = [
        row
        for row in catalog_payload["rooted_niemeier_lattices"]
        if not requested or row["label"] in requested
    ]
    known = {row["label"] for row in catalog_payload["rooted_niemeier_lattices"]}
    unknown = requested - known
    if unknown:
        raise ValueError(f"unknown rooted Niemeier labels: {sorted(unknown)}")

    all_sources = []
    ambient_accounting = []
    for row in ambient_rows:
        label = row["label"]
        ambient = matrix(ZZ, row["gram"])
        roots = signed_roots(ambient)
        sources, accounting = search_ambient(
            label, ambient, auxiliary, roots, target_genus, arguments
        )
        all_sources.extend(sources)
        ambient_accounting.append(accounting)
        print(
            "GOLAY720SOURCE|"
            f"ambient={label}|sources={len(sources)}|"
            f"embeddings={accounting['totals'].get('complete_auxiliary_embeddings', 0)}|"
            "status=PASS_AMBIENT_SEARCH",
            flush=True,
        )

    merged = {}
    for row in all_sources:
        key = row["reduced_gram_sha256"]
        if key not in merged:
            merged[key] = row
            merged[key]["ambient_provenance"] = []
            merged[key]["embedding_count_across_ambients"] = 0
        merged[key]["ambient_provenance"].append(row["ambient_label"])
        merged[key]["embedding_count_across_ambients"] += row[
            "embedding_count_merged"
        ]
    sources = [merged[key] for key in sorted(merged)]
    for index, row in enumerate(sources, 1):
        row["source_id"] = f"G720-S{index:04d}"
        row["ambient_provenance"] = sorted(set(row["ambient_provenance"]))

    success = [
        row
        for row in sources
        if row["source"]["mw_rank_for_rho_19"] in (1, 2)
        and row["source"]["support_count"] <= arguments.source_support_max
        and row["source"]["pole_audit"][
            "minimum_nonzero_section_pole_order"
        ]
        is not None
        and row["source"]["pole_audit"][
            "minimum_nonzero_section_pole_order"
        ]
        <= 2
    ]
    return {
        "schema": "elkies-k3.golay-octad-det720-prescribed-root-sources.v1",
        "status": (
            "PASS_SUCCESS_CONDITION_HIT" if success else "PASS_EXACT_SEARCH_NO_SUCCESS_HIT"
        ),
        "inputs": {
            relative(arguments.target): digest(arguments.target),
            relative(arguments.catalog): digest(arguments.catalog),
        },
        "fixed_auxiliary": {
            "name": "Golay-octad determinant-720 rank-seven auxiliary",
            "gram": rows(auxiliary),
            "gram_sha256": gram_digest(auxiliary),
            "determinant": 720,
            "basis_choice": "the seven original octad vectors, all of norm four",
        },
        "search_scope": {
            "rooted_niemeier_ambients": [row["label"] for row in ambient_rows],
            "source_root_rank": [
                arguments.source_root_rank_min,
                arguments.source_root_rank_max,
            ],
            "source_support_count": [
                arguments.source_support_min,
                arguments.source_support_max,
            ],
            "all_a_only": arguments.all_a_only,
            "fixed_dimension_max": arguments.fixed_dimension_max,
            "extension_limit": arguments.extension_limit,
            "final_label_limit": arguments.final_label_limit,
            "complete": not arguments.extension_limit and not arguments.final_label_limit,
        },
        "method": {
            "embedding": (
                "ordered norm-four Golay basis, sequential residual-Weyl dominant "
                "extension, prescribed final zero Dynkin labels, exact shifted ellipsoids"
            ),
            "acceptance": (
                "primitive full auxiliary, saturated orthogonal complement, exact "
                "complete roots, determinant and target-genus equality"
            ),
            "pole_audit": (
                "complete physical section shells through frame norm eight; "
                "unimodular MW quotient basis test for primitive-root MW1/MW2 rows"
            ),
        },
        "accounting": {
            "ambient_searches": ambient_accounting,
            "distinct_reduced_gram_sources": len(sources),
            "success_condition_hits": len(success),
        },
        "success_condition": {
            "definition": "MW rank 1 or 2, minimum nonzero-section P.O. <= 2, and support count within the declared few-support bound",
            "source_ids": [row["source_id"] for row in success],
        },
        "proof_boundary": {
            "proved": (
                "Every retained row is an exact primitive embedding of the fixed "
                "Golay-720 auxiliary in a hash-pinned full rooted Niemeier lattice, "
                "with exact saturated complement, roots, MW rank, and pole shell."
            ),
            "search_completeness": (
                "With zero explicit limits, sequential dominance is complete modulo "
                "the residual Weyl groups for the fixed ordered auxiliary basis and "
                "declared root-rank/support window. Diagram/umbral automorphisms are "
                "not quotiented, so duplicate embeddings may remain."
            ),
            "not_proved": (
                "Distinct reduced Grams need not be distinct integral-isometry or "
                "J2 classes. No rational marking, Weierstrass family, arithmetic "
                "descent, or elliptic-neighbour corridor is constructed."
            ),
        },
        "sources": sources,
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/enumerate_golay_det720_prescribed_root_sources.sage"
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, default=TARGET)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ambient-label", action="append", default=[])
    parser.add_argument("--source-root-rank-min", type=int, default=15)
    parser.add_argument("--source-root-rank-max", type=int, default=16)
    parser.add_argument("--source-support-min", type=int, default=1)
    parser.add_argument("--source-support-max", type=int, default=3)
    parser.add_argument("--all-a-only", action="store_true")
    parser.add_argument("--fixed-dimension-max", type=int, default=8)
    parser.add_argument("--extension-limit", type=int, default=0)
    parser.add_argument("--final-label-limit", type=int, default=0)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    arguments.target = arguments.target.resolve()
    arguments.catalog = arguments.catalog.resolve()
    arguments.output = arguments.output.resolve()

    result = build(arguments)
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if arguments.output.read_text() != serialized:
            raise SystemExit("Golay-720 prescribed-root source artifact is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(serialized)
    print(
        "GOLAY720SOURCE|"
        f"sources={len(result['sources'])}|"
        f"success={result['accounting']['success_condition_hits']}|"
        f"status={result['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
