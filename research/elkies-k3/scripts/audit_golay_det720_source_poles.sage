#!/usr/bin/env sage-python
"""Audit Golay-720 MW1/MW2 sources by complete quotient-class shells.

The prescribed-root hunt can produce a nonprimitive root lattice, hence
torsion in the frame/root quotient.  A root-plus-free-tail CVP is then not an
integral coordinate model.  This auditor instead uses the Smith quotient by
the simple roots, enumerates every torsion/free class whose exact height lower
bound can permit frame norm at most eight, and solves the corresponding
affine root-lattice closest-vector problem.  Each closest lift is moved to the
deterministic ADE chamber and checked against the physical-section component
gate.

Thus the audit covers primitive and nonprimitive root lattices uniformly.  A
rank-one or rank-two low-pole basis is certified only when the free quotient
coordinates of physical sections of norm 4, 6, or 8 contain a unimodular
basis.  No equation or rational marking is constructed.
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
from sage.all import QQ, ZZ, matrix, pari, vector


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_INPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-golay-octad-det720-prescribed-root-sources-v1.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-golay-octad-det720-source-poles-v1.json"
)

_engine_path = HERE / "exact_neighbor_engine.sage"
exec(compile(_engine_path.read_text(), str(_engine_path), "exec"), globals())

_shared_path = HERE / "enumerate_lattice_foundry_prescribed_root_sources.sage"
_shared = {"__file__": str(_shared_path), "__name__": "golay_pole_shared"}
exec(compile(_shared_path.read_text(), str(_shared_path), "exec"), _shared)
cartan_components = _shared["cartan_components"]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(ROOT))
    except ValueError:
        return str(path.resolve())


def dominate(value, frame, simple):
    value = vector(ZZ, value)
    for unused in range(10000):
        pairings = value * frame * simple.transpose()
        negative = next(
            (index for index, pairing in enumerate(pairings) if pairing < 0),
            None,
        )
        if negative is None:
            return value, tuple(map(int, pairings))
        root = simple.row(negative)
        value -= pairings[negative] * root
    raise RuntimeError("Golay-720 section Weyl dominance did not terminate")


def quotient_data(simple):
    smith, left, right = simple.smith_form()
    if smith != left * simple * right:
        raise AssertionError("unexpected Smith normal form convention")
    rank = simple.nrows()
    invariants = tuple(abs(int(smith[index, index])) for index in range(rank))
    return right, invariants


def make_gso(gram, float_type, precision):
    if float_type == "mpfr":
        FPLLL.set_precision(precision)
    integer = IntegerMatrix.from_matrix(
        [[int(entry) for entry in row] for row in gram.rows()]
    )
    gso = GSO.Mat(integer, gram=True, float_type=float_type, update=True)
    rank = gram.nrows()
    mu = [
        [gso.get_mu(i, j) if i > j else 0.0 for j in range(rank)]
        for i in range(rank)
    ]
    return gso, mu


def affine_cvp(gram, target, gso, mu):
    rank = gram.nrows()
    target_gso = [
        float(target[i])
        + sum(float(target[j]) * mu[j][i] for j in range(i + 1, rank))
        for i in range(rank)
    ]
    zero_distance = target * gram * target
    solutions = Enumeration(gso).enumerate(
        0,
        rank,
        float(zero_distance) + 1.0,
        0,
        target=target_gso,
    )
    if not solutions:
        raise RuntimeError("Golay-720 torsion affine CVP returned no solution")
    reported_distance, coordinates = solutions[0]
    integral = tuple(int(round(value)) for value in coordinates)
    if any(
        abs(value - integer) > 1e-7
        for value, integer in zip(coordinates, integral)
    ):
        raise RuntimeError(
            "Golay-720 torsion affine CVP returned nonintegral coordinates"
        )
    displacement = vector(QQ, integral) - target
    exact_distance = displacement * gram * displacement
    return integral, abs(float(reported_distance) - float(exact_distance))


def section_shell(source):
    frame = matrix(ZZ, source["gram"])
    simple, unused_positive, cartan_rows = deterministic_simple_roots(frame)
    root_rank = simple.nrows()
    mw_rank = frame.nrows() - root_rank
    if mw_rank not in (1, 2):
        raise ValueError("Golay-720 pole audit expects MW rank one or two")
    cartan = matrix(ZZ, cartan_rows)
    components = cartan_components(cartan)
    smith_right, invariants = quotient_data(simple)
    smith_inverse = smith_right.inverse()
    quotient_frame = smith_inverse * frame * smith_inverse.transpose()
    diagonal = matrix.diagonal(ZZ, invariants)
    root_block = matrix(QQ, quotient_frame[:root_rank, :root_rank])
    cross = matrix(QQ, quotient_frame[:root_rank, root_rank:])
    tail = matrix(QQ, quotient_frame[root_rank:, root_rank:])
    height = tail - cross.transpose() * root_block.inverse() * cross
    scaled_root_gram = diagonal * root_block * diagonal
    primary_gso, primary_mu = make_gso(scaled_root_gram, "dd", 0)
    audit_gso, audit_mu = make_gso(scaled_root_gram, "mpfr", 256)

    denominator = math.lcm(*(entry.denominator() for entry in height.list()))
    integral_height = (denominator * height).change_ring(ZZ)
    tail_result = pari(integral_height).qfminim(8 * denominator)
    half_tails = [
        tuple(map(int, column))
        for column in matrix(ZZ, tail_result[2].sage()).columns()
    ]
    tails = {(0,) * mw_rank}
    for tail_value in half_tails:
        tails.add(tail_value)
        tails.add(tuple(-entry for entry in tail_value))
    torsion_ranges = [range(invariant) for invariant in invariants]
    sections = {}
    norm_histogram = Counter()
    maximum_error = 0.0
    tested_classes = 0
    for torsion_full in itertools.product(*torsion_ranges):
        for free in sorted(tails):
            if not any(torsion_full) and not any(free):
                continue
            linear = diagonal * (
                root_block * vector(QQ, torsion_full) + cross * vector(QQ, free)
            )
            target = -(scaled_root_gram.inverse() * linear)
            primary, primary_error = affine_cvp(
                scaled_root_gram, target, primary_gso, primary_mu
            )
            audited, audit_error = affine_cvp(
                scaled_root_gram, target, audit_gso, audit_mu
            )
            if primary != audited:
                primary_vector = vector(QQ, primary) - target
                audit_vector = vector(QQ, audited) - target
                primary_norm = primary_vector * scaled_root_gram * primary_vector
                audit_norm = audit_vector * scaled_root_gram * audit_vector
                if primary_norm != audit_norm:
                    raise ValueError(
                        "cross-precision Golay-720 torsion affine CVP mismatch"
                    )
            maximum_error = max(maximum_error, primary_error, audit_error)
            root_coordinates = vector(ZZ, torsion_full) + diagonal * vector(
                ZZ, primary
            )
            quotient = vector(ZZ, tuple(root_coordinates) + tuple(free))
            raw = quotient * smith_inverse
            norm = int(raw * frame * raw)
            norm_histogram[norm] += 1
            tested_classes += 1
            if norm not in (4, 6, 8):
                continue
            dominant, pairings = dominate(raw, frame, simple)
            if any(pairing not in (0, 1) for pairing in pairings):
                continue
            if any(
                sum(pairings[index] for index in component) > 1
                for component in components
            ):
                continue
            checked_quotient = dominant * smith_right
            checked_torsion = tuple(
                int(checked_quotient[index] % invariant)
                for index, invariant in enumerate(invariants)
            )
            checked_free = tuple(map(int, checked_quotient[root_rank:]))
            expected_torsion = tuple(
                int(value % invariant)
                for value, invariant in zip(torsion_full, invariants)
            )
            if checked_torsion != expected_torsion or checked_free != tuple(
                free
            ):
                raise AssertionError(
                    "Weyl reduction changed the Smith quotient class"
                )
            torsion = tuple(
                value
                for value, invariant in zip(checked_torsion, invariants)
                if invariant > 1
            )
            key = (norm, checked_free, torsion, tuple(map(int, dominant)))
            sections[key] = {
                "frame_vector": list(map(int, dominant)),
                "frame_norm": norm,
                "pole_order": (norm - 4) // 2,
                "simple_root_pairings": list(pairings),
                "torsion_coordinates": list(torsion),
                "free_mw_coordinates": list(checked_free),
            }
    ordered = [sections[key] for key in sorted(sections)]
    basis = None
    if mw_rank == 1:
        basis = next(
            (
                [section]
                for section in ordered
                if abs(section["free_mw_coordinates"][0]) == 1
            ),
            None,
        )
    else:
        for left_section, right_section in itertools.combinations(ordered, 2):
            free = matrix(
                ZZ,
                [
                    left_section["free_mw_coordinates"],
                    right_section["free_mw_coordinates"],
                ],
            )
            if abs(free.det()) == 1:
                basis = [left_section, right_section]
                break
    non_torsion = [
        section for section in ordered if any(section["free_mw_coordinates"])
    ]
    return {
        "status": (
            "PASS_EXACT_PHYSICAL_MW_BASIS_THROUGH_POLE_TWO"
            if basis is not None
            else "PASS_EXACT_NO_PHYSICAL_MW_BASIS_THROUGH_POLE_TWO"
        ),
        "root_smith_invariants": list(invariants),
        "torsion_order": math.prod(invariants),
        "complete_shell_through_frame_norm": 8,
        "height_gram": [[str(entry) for entry in row] for row in height.rows()],
        "height_eligible_free_classes": len(tails),
        "torsion_classes": math.prod(invariants),
        "affine_classes_tested": tested_classes,
        "minimum_norm_histogram_by_affine_class": {
            str(norm): count for norm, count in sorted(norm_histogram.items())
        },
        "maximum_cross_precision_cvp_distance_error": maximum_error,
        "physical_sections_through_pole_two": len(ordered),
        "minimum_nonidentity_section_pole_order": (
            min(section["pole_order"] for section in ordered) if ordered else None
        ),
        "minimum_non_torsion_section_pole_order": (
            min(section["pole_order"] for section in non_torsion)
            if non_torsion
            else None
        ),
        "basis_with_all_poles_at_most_two": basis is not None,
        "basis_sorted_pole_profile": (
            sorted(section["pole_order"] for section in basis)
            if basis is not None
            else None
        ),
        "basis": basis,
        "sections": ordered,
    }


def build(arguments):
    payload = json.loads(arguments.input.read_text())
    expected = "elkies-k3.golay-octad-det720-prescribed-root-sources.v1"
    if payload.get("schema") != expected:
        raise ValueError(
            f"unexpected Golay-720 source schema: {payload.get('schema')}"
        )
    audits = []
    for index, entry in enumerate(payload["sources"], 1):
        source = entry["source"]
        audit = section_shell(source)
        audits.append(
            {
                "source_id": entry["source_id"],
                "source_gram_sha256": source["gram_sha256"],
                "root_type": source["root_type"],
                "mw_rank_for_rho_19": source["mw_rank_for_rho_19"],
                "support_count": source["support_count"],
                "ambient_provenance": entry["ambient_provenance"],
                "audit": audit,
            }
        )
        if arguments.progress_every and index % arguments.progress_every == 0:
            print(
                f"GOLAY720POLE|audited={index}|total={len(payload['sources'])}",
                flush=True,
            )
    hits = [
        row
        for row in audits
        if row["mw_rank_for_rho_19"] in (1, 2)
        and row["support_count"] <= 3
        and row["audit"]["basis_with_all_poles_at_most_two"]
    ]
    return {
        "schema": "elkies-k3.golay-octad-det720-source-poles.v1",
        "status": (
            "PASS_SUCCESS_CONDITION_HIT"
            if hits
            else "PASS_EXACT_NO_LOW_POLE_MW_BASIS"
        ),
        "inputs": {relative(arguments.input): digest(arguments.input)},
        "method": (
            "complete exact-height torsion/free quotient-class shell through "
            "norm eight, cross-precision affine root CVP, deterministic ADE "
            "Weyl dominance, physical component gate, and exact Smith coordinates"
        ),
        "proof_boundary": {
            "proved": (
                "For every supplied MW1/MW2 source, every torsion/free quotient "
                "class capable by exact height of frame norm at most eight is "
                "tested, and every retained norm, physical pairing, Smith class, "
                "and unimodular free MW basis is checked exactly."
            ),
            "not_proved": (
                "Closest-vector branch decisions have the same independent "
                "double-double/MPFR-256 numerical boundary as the repository's "
                "rank-one and rank-two pole auditors. The input source enumeration "
                "has its own declared boundary. No rational marking, equation, or "
                "neighbour route is built."
            ),
        },
        "accounting": {
            "source_rows": len(audits),
            "success_condition_hits": len(hits),
        },
        "success_condition": {
            "definition": (
                "MW rank one or two with a unimodular physical section basis whose "
                "every member has pole order at most two"
            ),
            "source_ids": [row["source_id"] for row in hits],
        },
        "audits": audits,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--progress-every", type=int, default=25)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    arguments.input = arguments.input.resolve()
    arguments.output = arguments.output.resolve()
    result = build(arguments)
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if arguments.output.read_text() != serialized:
            raise SystemExit("Golay-720 pole artifact is stale")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(serialized)
    print(
        "GOLAY720POLE|"
        f"sources={result['accounting']['source_rows']}|"
        f"success={result['accounting']['success_condition_hits']}|"
        f"status={result['status']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
