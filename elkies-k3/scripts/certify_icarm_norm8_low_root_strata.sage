#!/usr/bin/env sage-python
"""Certify the root-rank strata inside the 11952 norm-eight screen.

The historical ICARM screen called all 63,917 minimum-norm-eight classes an
``A1/MW16 layer``.  That is correct only when the parity class has one
minimum representative up to sign.  If it has ``m`` such representatives,
the associated residual-chord pencil has ``m`` distinct split members.  The
complete singular-pencil certificate proves that these account for the full
even discriminant part, so its semistable Jacobian has

    m I2 + (24 - 2m) I1,

root lattice ``m A1``, and geometric Mordell--Weil rank ``17-m`` at Picard
rank 19.  This checker joins the exact shell enumeration, the priority table,
and the complete discriminant certificate.  It also proves that every
norm-eight parity class contains norm-twelve vectors, hence supplies
degree-one sections, and replays one primitive-U/frame control per stratum.
"""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path
import runpy

from sage.all import QQ, ZZ, block_diagonal_matrix, matrix, pari, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
ENUMERATOR = ROOT / "elkies-k3/scripts/enumerate_rootless_bisection_orbits.sage"
HISTORICAL = ROOT / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
TABLE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-norm8-pencil-priority-v1.tsv"
SINGULAR = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-singular-bisection-search-complete-v1.json"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-icarm-11952-norm8-low-root-strata-v1.json"


def digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def parity_mask(value) -> int:
    return sum((int(entry) & 1) << index for index, entry in enumerate(value))


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def section_class(value, gram):
    value = vector(ZZ, value)
    first = (value * gram * value - 2) / 2
    if first not in ZZ:
        raise ArithmeticError("section class is not integral")
    return vector(ZZ, [ZZ(first), 1] + list(value))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    historical_payload = json.loads(HISTORICAL.read_text())
    direct = json.loads(DIRECT.read_text())
    singular = json.loads(SINGULAR.read_text())
    historical_gram = matrix(ZZ, historical_payload["rootless_frame"])
    direct_gram = matrix(ZZ, direct["frame_certificate"]["frame_gram"])
    height_gram = matrix(ZZ, direct["sections"]["height_gram"])
    historical_to_direct = matrix(
        ZZ, direct["frame_certificate"]["integral_isometry_to_alternate_Q80"]
    )
    section_coordinates = matrix(
        ZZ, direct["sections"]["coordinate_matrix_in_compiled_frame"]
    )
    if historical_to_direct * historical_gram * historical_to_direct.transpose() != direct_gram:
        raise ArithmeticError("historical/direct rootless isometry changed")
    if abs(section_coordinates.det()) != 1:
        raise ArithmeticError("equation section marking is not unimodular")
    if height_gram.det() != 948 or int(pari(height_gram).qfminim(2)[0]):
        raise ArithmeticError("direct equation basis is no longer rootless R17")

    short_change = historical_gram.LLL_gram().transpose()
    short_gram = short_change * historical_gram * short_change.transpose()
    short_to_section_qq = (
        short_change * historical_to_direct.inverse() * section_coordinates.inverse()
    )
    if any(entry.denominator() != 1 for entry in short_to_section_qq):
        raise ArithmeticError("short-to-section transport is not integral")
    short_to_section = matrix(ZZ, short_to_section_qq)

    enumerator = runpy.run_path(str(ENUMERATOR))
    shell = enumerator["streaming_short_vectors"](short_gram, bound=12)
    signed_counts = {int(key): int(value) for key, value in shell["signed_counts"].items()}
    mask_counts = {int(key): len(value) for key, value in shell["masks_by_norm"].items()}
    expected_signed = {4: 2626, 6: 53290, 8: 460360, 10: 2472050, 12: 9618310}
    expected_masks = {2: 0, 4: 1313, 6: 26645, 8: 63917, 10: 65792, 12: 65279}
    if signed_counts != expected_signed or mask_counts != expected_masks:
        raise ArithmeticError("alternate R17 shell census changed")

    masks_by_norm = shell["masks_by_norm"]
    norm_eight_masks = masks_by_norm[8] - masks_by_norm[4] - masks_by_norm[6]
    if len(norm_eight_masks) != 63917:
        raise ArithmeticError("norm-eight minimum class count changed")
    if not norm_eight_masks <= masks_by_norm[12]:
        raise ArithmeticError("a norm-eight pencil class has no degree-one section")
    deep_norm_twelve = masks_by_norm[12] - masks_by_norm[4] - masks_by_norm[8]
    if len(deep_norm_twelve) != 49:
        raise ArithmeticError("alternate norm-twelve deep-hole count changed")

    table_rows = []
    with TABLE.open(newline="") as stream:
        for expected_rank, row in enumerate(csv.DictReader(stream, delimiter="\t"), start=1):
            if int(row["priority_rank"]) != expected_rank:
                raise ArithmeticError("priority table is not rank-contiguous")
            short_word = tuple(map(int, row["short_basis_w"].split()))
            mask = parity_mask(short_word)
            if mask not in norm_eight_masks:
                raise ArithmeticError("priority row is not a norm-eight minimum class")
            table_rows.append(
                {
                    "priority_rank": expected_rank,
                    "mask": mask,
                    "minimum_unoriented_count": int(row["minimal_unoriented_count"]),
                    "section_basis_w": tuple(map(int, row["section_basis_w"].split())),
                }
            )
    if len(table_rows) != 63917 or {row["mask"] for row in table_rows} != norm_eight_masks:
        raise ArithmeticError("priority table does not cover the norm-eight minimum classes")

    multiplicity_histogram = Counter(
        row["minimum_unoriented_count"] for row in table_rows
    )
    expected_multiplicity = {
        1: 1266,
        2: 8410,
        3: 20348,
        4: 21405,
        5: 9861,
        6: 2280,
        7: 331,
        8: 16,
    }
    if dict(sorted(multiplicity_histogram.items())) != expected_multiplicity:
        raise ArithmeticError("minimum-representative histogram changed")
    if singular.get("status") != "PASS_COMPLETE_NO_NONSPLIT_RATIONAL_SINGULAR_MEMBER":
        raise ArithmeticError("complete singular-pencil certificate is not passing")
    singular_histogram = {
        int(key): int(value)
        for key, value in singular["minimal_unoriented_split_member_count_histogram"].items()
    }
    even_histogram = {
        int(key): int(value)
        for key, value in singular["finite_even_discriminant_degree_histogram"].items()
    }
    if singular_histogram != expected_multiplicity:
        raise ArithmeticError("singular-pencil split-member histogram changed")
    if even_histogram != {
        2 * (multiplicity - 1): count
        for multiplicity, count in expected_multiplicity.items()
    }:
        raise ArithmeticError("split members no longer exhaust the finite even discriminant")

    norm_twelve_multiplicity = shell["unoriented_multiplicities"]
    degree_one_counts = {
        row["priority_rank"]: 2 * int(norm_twelve_multiplicity[row["mask"]])
        for row in table_rows
    }
    if not degree_one_counts or min(degree_one_counts.values()) <= 0:
        raise ArithmeticError("a pencil class has no degree-one section")
    degree_one_histogram = Counter(degree_one_counts.values())

    ns_gram = block_diagonal_matrix(matrix(ZZ, [[0, 1], [1, 0]]), -height_gram)
    controls = []
    for multiplicity in sorted(expected_multiplicity):
        row = next(
            row for row in table_rows
            if row["minimum_unoriented_count"] == multiplicity
        )
        norm_twelve_short = vector(ZZ, shell["representatives"][row["mask"]])
        norm_twelve_section = norm_twelve_short * short_to_section
        trace = vector(ZZ, row["section_basis_w"])
        difference = trace - norm_twelve_section
        if any(entry % 2 for entry in difference):
            raise ArithmeticError("norm-eight/norm-twelve representatives differ modulo 2M")
        new_zero_word = vector(ZZ, [entry // 2 for entry in difference])
        fibre = vector(ZZ, [2, 2] + list(trace))
        new_zero = section_class(new_zero_word, height_gram)
        if fibre * ns_gram * fibre or fibre * ns_gram * new_zero != 1:
            raise ArithmeticError("primitive-U control has wrong intersections")
        mate = fibre + new_zero
        complement = matrix(
            ZZ, [list(fibre * ns_gram), list(mate * ns_gram)]
        ).right_kernel_matrix()
        transport = matrix(
            ZZ,
            [list(fibre), list(mate)] + [list(value) for value in complement.rows()],
        )
        child_frame = -(complement * ns_gram * complement.transpose())
        minimum = pari(matrix(ZZ, child_frame)).qfminim(2).sage()
        root_count = int(minimum[0])
        positive_roots = matrix(ZZ, minimum[2]).transpose()
        root_gram = positive_roots * child_frame * positive_roots.transpose()
        if (
            abs(transport.det()) != 1
            or child_frame.det() != 948
            or root_count != 2 * multiplicity
            or positive_roots.rank() != multiplicity
            or root_gram != 2 * matrix.identity(ZZ, multiplicity)
        ):
            raise ArithmeticError(f"{multiplicity}A1 primitive-frame control failed")
        controls.append(
            {
                "minimum_unoriented_split_member_count": multiplicity,
                "priority_rank": row["priority_rank"],
                "degree_one_section_count": degree_one_counts[row["priority_rank"]],
                "primitive_u_transport_determinant": int(transport.det()),
                "child_frame_determinant": int(child_frame.det()),
                "signed_norm_two_root_count": root_count,
                "root_span_rank": int(positive_roots.rank()),
                "positive_root_gram": rows(root_gram),
            }
        )

    strata = []
    for multiplicity, count in sorted(expected_multiplicity.items()):
        strata.append(
            {
                "minimum_unoriented_split_member_count": multiplicity,
                "class_count": count,
                "finite_even_discriminant_degree": 2 * (multiplicity - 1),
                "semistable_fibre_configuration": (
                    f"{multiplicity}I2+{24 - 2 * multiplicity}I1"
                ),
                "root_lattice": "A1" if multiplicity == 1 else f"{multiplicity}A1",
                "root_rank": multiplicity,
                "geometric_mw_rank_at_rho_19": 17 - multiplicity,
            }
        )

    payload = {
        "schema": "elkies-k3.icarm-11952-norm8-low-root-strata.v1",
        "status": "PASS_EXACT_COMPLETE_NORM8_LOW_ROOT_STRATIFICATION",
        "source_chart": "norm12-orbit-11952 alternate-Q80 rootless/MW17",
        "scope": {
            "minimum_norm_eight_old_degree_two_class_count": len(table_rows),
            "stratum_count": len(strata),
            "minimum_root_rank": 1,
            "maximum_root_rank": 8,
            "degree_one_section_count_minimum": min(degree_one_counts.values()),
            "degree_one_section_count_maximum": max(degree_one_counts.values()),
            "all_classes_have_degree_one_sections": True,
        },
        "shell_census": {
            "method": "LLL-reduced Fincke--Pohst traversal with exact leaf norms",
            "bound": 12,
            "signed_shell_counts": {
                str(key): value for key, value in sorted(signed_counts.items())
            },
            "parity_masks_hit_by_shell": {
                str(key): value for key, value in sorted(mask_counts.items())
            },
            "norm_twelve_deep_hole_count": len(deep_norm_twelve),
            "degree_one_section_count_histogram": {
                str(key): value for key, value in sorted(degree_one_histogram.items())
            },
        },
        "strata": strata,
        "primitive_frame_controls": controls,
        "argument": {
            "split_members": (
                "A parity class with m minimum norm-eight representatives up to sign "
                "has m distinct split members in its residual-chord pencil."
            ),
            "even_discriminant_exhaustion": singular["even_multiplicity_exhaustion"],
            "fibre_conclusion": (
                "The regular chord gauge places one split I2 at infinity. The other "
                "m-1 split members account for the entire finite even discriminant "
                "degree 2(m-1); all remaining discriminant roots are simple. Hence the "
                "semistable configuration is mI2+(24-2m)I1 and the root lattice is mA1."
            ),
            "rank_conclusion": (
                "Every class has a norm-twelve representative in the same parity coset, "
                "so it has a degree-one section and a primitive U splitting. Picard rank "
                "19 and Shioda--Tate give geometric MW rank 17-m."
            ),
        },
        "proof_boundary": (
            "This is complete for the 63,917 minimum-norm-eight, old-degree-two "
            "residual-chord classes on source chart 11952. It corrects the earlier "
            "uniform A1/MW16 label: only the m=1 stratum has that type. It is not a "
            "classification of all elliptic fibrations on X948, and it does not include "
            "A2/MW15 or classes of another old degree, trace norm, or source chart."
        ),
        "inputs": {
            relative(path): digest(path)
            for path in (Path(__file__).resolve(), ENUMERATOR, HISTORICAL, DIRECT, TABLE, SINGULAR)
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "pari_version": ".".join(map(str, pari.version())),
            "required_features": ["exact Fincke--Pohst traversal", "PARI qfminim"],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/certify_icarm_norm8_low_root_strata.sage"
        ),
    }
    serialized = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    output_path = args.output.resolve()
    if args.check:
        if not output_path.is_file() or output_path.read_text() != serialized:
            raise ArithmeticError("stored low-root stratum certificate differs from replay")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "ICARMLOWROOTSTRATA|classes={}|strata={}|mw16={}|mw15={}|mw9={}|status={}|output={}".format(
            len(table_rows),
            len(strata),
            expected_multiplicity[1],
            expected_multiplicity[2],
            expected_multiplicity[8],
            payload["status"],
            relative(output_path),
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
