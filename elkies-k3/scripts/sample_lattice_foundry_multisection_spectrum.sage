#!/usr/bin/env sage-python
"""Measure low-degree, low-genus multisection richness of rootless targets.

For a rootless marked fibration ``NS = U + M(-1)``, write a divisor class as
``D=(a,d,w)`` with ``D.F=d`` and arithmetic genus ``g``.  Then

    a = (w.M.w + 2*g - 2)/(2*d),

and section translation sends ``w`` to ``w+d*x``.  The minimum intersection
with an old section in that translation orbit is nonnegative exactly when

    min_{x in M} (w-d*x).M.(w-d*x) >= 2*d^2 - 2*g + 2.

Degree two is enumerated exactly through the first admissible low-height
shells.  Degree three and four use deterministic residue-class samples and an
exact rational LDL branch-and-bound inside a declared norm cap.  Those rows
are exact statements about the sampled lattice cosets, not a complete orbit
census.  Except for section-nonnegative rational bisections on a rootless
fibration, irreducibility and arithmetic descent remain separate gates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from collections import Counter
from pathlib import Path

from sage.all import QQ, ZZ, ceil, matrix, pari, vector


ROOT = Path(__file__).resolve().parents[2]
DATABASE = ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-v1.json"
PINNED_R17 = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-lattice-foundry-multisection-spectrum-v1.json"
)


def load_matrix(path: Path):
    return matrix(
        ZZ,
        [
            [ZZ(entry) for entry in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def parity_mask(value) -> int:
    return sum((int(entry) % 2) << index for index, entry in enumerate(value))


def ldl_data(gram):
    """Return exact unit-lower LDL data in the convention used by the shell walker."""

    dimension = gram.nrows()
    lower = matrix(QQ, dimension, dimension)
    diagonal = []
    for index in range(dimension):
        lower[index, index] = 1
        value = QQ(gram[index, index]) - sum(
            lower[index, prior] ** 2 * diagonal[prior]
            for prior in range(index)
        )
        assert value > 0
        diagonal.append(value)
        for row in range(index + 1, dimension):
            lower[row, index] = (
                QQ(gram[row, index])
                - sum(
                    lower[row, prior]
                    * lower[index, prior]
                    * diagonal[prior]
                    for prior in range(index)
                )
            ) / value
    return lower, diagonal


def minimum_coset_norm_through_bound(gram, lower, diagonal, degree, residue, bound):
    """Return the exact coset minimum if it is at most ``bound``.

    This enumerates all integral ``x`` satisfying
    ``(residue+degree*x).G.(residue+degree*x) <= bound``.  Coordinate ranges
    use an intentionally broad exact bound and every branch is then tested
    with rational arithmetic; no floating-point pruning enters the result.
    """

    dimension = gram.nrows()
    shift = [QQ(entry) / degree for entry in residue]
    coordinates = [ZZ.zero()] * dimension
    best = None

    def visit(index, used):
        nonlocal best
        if best is not None and used > QQ(best) / (degree * degree):
            return
        if index < 0:
            scaled = vector(ZZ, [
                int(residue[i] + degree * coordinates[i])
                for i in range(dimension)
            ])
            exact = int(scaled * gram * scaled)
            assert exact <= bound
            if best is None or exact < best:
                best = exact
            return

        remaining = QQ(bound) / (degree * degree) - used
        if remaining < 0:
            return
        center = shift[index] + sum(
            lower[row, index] * (coordinates[row] + shift[row])
            for row in range(index + 1, dimension)
        )
        radius_squared = remaining / diagonal[index]
        # If |x+center| <= sqrt(R), then |x| is at most
        # ceil(sqrt(R))+ceil(|center|).  This wider integer interval avoids
        # any algebraic/floating endpoint decision; the exact square test
        # below removes its harmless excess.
        floor_radius = ZZ(radius_squared.floor()).isqrt()
        radius_ceiling = floor_radius
        if QQ(radius_ceiling * radius_ceiling) < radius_squared:
            radius_ceiling += 1
        coordinate_bound = int(radius_ceiling + ceil(abs(center)))
        for entry in range(-coordinate_bound, coordinate_bound + 1):
            contribution = diagonal[index] * (QQ(entry) + center) ** 2
            if contribution > remaining:
                continue
            coordinates[index] = ZZ(entry)
            visit(index - 1, used + contribution)

    visit(dimension - 1, QQ.zero())
    return best


def exact_degree_two_spectrum(gram):
    """Enumerate the norm shells through ten and deduplicate modulo ``2M``."""

    reduced_change = gram.LLL_gram().transpose()
    reduced = reduced_change * gram * reduced_change.transpose()
    assert abs(reduced_change.det()) == 1
    result = pari(reduced).qfminim(10)
    columns = matrix(ZZ, result[2].sage()).columns()
    masks_by_norm = {norm: set() for norm in range(2, 11, 2)}
    signed_counts = Counter()
    for column in columns:
        value = vector(ZZ, column)
        norm = int(value * reduced * value)
        signed_counts[norm] += 2
        masks_by_norm[norm].add(parity_mask(value))
    assert sum(signed_counts.values()) == int(result[0])
    assert not masks_by_norm[2]

    rational_minimum_ten = masks_by_norm[10] - masks_by_norm[6] - masks_by_norm[2]
    genus_one_minimum_eight = masks_by_norm[8] - masks_by_norm[4] - {0}
    return {
        "method": "complete PARI short-vector shells and exact M/2M masks",
        "enumerated_through_norm": 10,
        "signed_shell_counts": {
            str(norm): int(signed_counts.get(norm, 0))
            for norm in range(2, 11, 2)
        },
        "rational_bisections": {
            "status": "PASS_EXACT_LOW_HEIGHT_TRANSLATION_ORBITS_AND_GEOMETRY",
            "arithmetic_genus": 0,
            "section_nonnegative_threshold_norm": 10,
            "translation_orbits_with_minimum_norm_ten": len(rational_minimum_ten),
            "geometric_consequence": (
                "On any K3 realization of this rootless fibration, each class is an "
                "irreducible smooth rational bisection. Arithmetic field of definition "
                "is not inferred from the geometric lattice."
            ),
        },
        "genus_one_bisections": {
            "status": "PASS_EXACT_LOW_HEIGHT_LATTICE_CANDIDATES",
            "arithmetic_genus": 1,
            "section_nonnegative_threshold_norm": 8,
            "translation_orbits_with_minimum_norm_eight": len(genus_one_minimum_eight),
            "geometric_boundary": (
                "Section nonnegativity is proved, but global nefness, irreducibility, "
                "and arithmetic descent are not asserted."
            ),
        },
    }


def sample_degree(gram, frame_id, degree, genera, sample_count, slack, seed):
    change = gram.LLL_gram().transpose()
    reduced = change * gram * change.transpose()
    lower, diagonal = ldl_data(reduced)
    dimension = reduced.nrows()
    frame_seed = int(hashlib.sha256(frame_id.encode()).hexdigest()[:16], 16)
    rng = random.Random(seed ^ frame_seed ^ (degree << 24))
    residues = set()
    zero = (0,) * dimension
    residues.add(zero)
    while len(residues) < sample_count:
        residues.add(tuple(rng.randrange(degree) for unused in range(dimension)))

    thresholds = {genus: 2 * degree * degree - 2 * genus + 2 for genus in genera}
    bound = max(thresholds.values()) + slack
    categories = {
        genus: {
            "arithmetic_genus": genus,
            "section_nonnegative_threshold_norm": thresholds[genus],
            "low_height_norm_cap": thresholds[genus] + slack,
            "sampled_qualifying_translation_cosets": 0,
            "minimum_norm_histogram": Counter(),
        }
        for genus in genera
    }
    minima_above_bound = 0
    minima = Counter()
    for residue_tuple in sorted(residues):
        residue = vector(ZZ, residue_tuple)
        minimum = minimum_coset_norm_through_bound(
            reduced, lower, diagonal, degree, residue, bound
        )
        if minimum is None:
            minima_above_bound += 1
            continue
        minima[minimum] += 1
        residue_norm = int(residue * reduced * residue)
        for genus in genera:
            if (residue_norm + 2 * genus - 2) % (2 * degree):
                continue
            threshold = thresholds[genus]
            if threshold <= minimum <= threshold + slack:
                row = categories[genus]
                row["sampled_qualifying_translation_cosets"] += 1
                row["minimum_norm_histogram"][minimum] += 1

    for row in categories.values():
        row["minimum_norm_histogram"] = {
            str(key): value for key, value in sorted(row["minimum_norm_histogram"].items())
        }
        row["sample_fraction_of_all_cosets"] = (
            row["sampled_qualifying_translation_cosets"] / len(residues)
        )
        row["geometric_boundary"] = (
            "Exact for the sampled translation cosets. Section nonnegativity follows "
            "from the exact coset minimum, but effectivity, irreducibility, arithmetic "
            "field of definition, and rank gain are not asserted."
        )
    return {
        "method": "deterministic uniform residue sample with exact rational-LDL CVP bound",
        "degree": degree,
        "total_translation_cosets": degree ** dimension,
        "sample_size": len(residues),
        "seed": seed,
        "norm_bound": bound,
        "sampled_cosets_with_minimum_above_bound": minima_above_bound,
        "sampled_minimum_norm_histogram": {
            str(key): value for key, value in sorted(minima.items())
        },
        "categories": [categories[genus] for genus in genera],
    }


def target_rows(database):
    frame_by_id = {
        frame["frame_id"]: frame
        for ns in database["ns_classes"]
        for frame in ns["frames"]
    }
    return [
        (target, frame_by_id[target["frame_id"]])
        for target in database["rootless_targets"]
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DATABASE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--frame-id", action="append", default=[])
    parser.add_argument("--limit", type=int)
    parser.add_argument("--sample-count", type=int, default=256)
    parser.add_argument("--height-slack", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--pari-stack-gb", type=int, default=4)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.sample_count <= 0 or arguments.height_slack < 0:
        parser.error("sample count must be positive and height slack nonnegative")
    if arguments.pari_stack_gb <= 0:
        parser.error("PARI stack size must be positive")
    pari.allocatemem(arguments.pari_stack_gb * 1024**3)

    database_path = arguments.database.resolve()
    database = json.loads(database_path.read_text())
    selected = target_rows(database)
    if arguments.frame_id:
        wanted = set(arguments.frame_id)
        selected = [row for row in selected if row[0]["frame_id"] in wanted]
        missing = wanted - {row[0]["frame_id"] for row in selected}
        if missing:
            parser.error(f"unknown or non-rootless frame ids: {sorted(missing)}")
    else:
        selected.sort(
            key=lambda row: (
                row[0]["is_existing_H3_control"],
                row[0]["determinant"],
                row[0]["frame_id"],
            )
        )
    if arguments.limit is not None:
        selected = selected[: arguments.limit]
    if not selected:
        parser.error("no targets selected")

    pinned = load_matrix(PINNED_R17)
    rows = []
    for index, (target, frame) in enumerate(selected, 1):
        gram = matrix(ZZ, frame["gram"])
        assert gram.nrows() == 17 and pari(gram).qfminim(2)[0] == 0
        is_published = pari(gram).qfisom(pari(pinned)) != 0
        degree_two = exact_degree_two_spectrum(gram)
        degree_three = sample_degree(
            gram, target["frame_id"], 3, (0, 1),
            arguments.sample_count, arguments.height_slack, arguments.seed,
        )
        degree_four = sample_degree(
            gram, target["frame_id"], 4, (0, 1, 2),
            arguments.sample_count, arguments.height_slack, arguments.seed,
        )
        spectrum = {
            "rational_bisection_orbits_minimum_height": degree_two[
                "rational_bisections"
            ]["translation_orbits_with_minimum_norm_ten"],
            "genus_one_bisection_candidate_orbits_minimum_height": degree_two[
                "genus_one_bisections"
            ]["translation_orbits_with_minimum_norm_eight"],
            "sampled_rational_trisection_candidates": degree_three[
                "categories"
            ][0]["sampled_qualifying_translation_cosets"],
            "sampled_genus_one_trisection_candidates": degree_three[
                "categories"
            ][1]["sampled_qualifying_translation_cosets"],
            "sampled_low_genus_quadrisection_candidates": sum(
                row["sampled_qualifying_translation_cosets"]
                for row in degree_four["categories"]
            ),
        }
        rows.append(
            {
                "frame_id": target["frame_id"],
                "ns_id": target["ns_id"],
                "determinant": target["determinant"],
                "published_R17_isometry_control": bool(is_published),
                "degree_two": degree_two,
                "degree_three_sample": degree_three,
                "degree_four_sample": degree_four,
                "richness_coordinates": spectrum,
            }
        )
        print(
            f"FOUNDRYMULTI|index={index}/{len(selected)}|frame={target['frame_id']}|"
            f"bisections={spectrum['rational_bisection_orbits_minimum_height']}|"
            f"g1_bisections={spectrum['genus_one_bisection_candidate_orbits_minimum_height']}|"
            f"status=PASS",
            flush=True,
        )

    output = {
        "schema": "elkies-k3.lattice-foundry-multisection-spectrum.v1",
        "status": "PASS_EXACT_D2_AND_EXACT_BOUNDED_SAMPLED_D3_D4_LATTICE_SPECTRA",
        "proof_boundary": {
            "proved": (
                "Degree-two low-height translation-orbit counts are complete through "
                "norm ten. Every reported degree-three/four coset minimum is exact "
                "inside a deterministic, explicitly sized sample."
            ),
            "not_proved": (
                "The degree-three/four samples are not complete censuses. Lattice "
                "candidates are not called irreducible or defined over QQ without "
                "separate geometric and arithmetic evidence. Richness is a search "
                "heuristic, not a proved predictor of exceptional specialization rank."
            ),
        },
        "formulae": {
            "class": "D=((w.M.w+2*g-2)/(2*d),d,w)",
            "translation": "w -> w+d*x",
            "section_nonnegative_threshold": "min_norm(w+dM) >= 2*d^2-2*g+2",
        },
        "search": {
            "selected_frame_ids": [row[0]["frame_id"] for row in selected],
            "limit": arguments.limit,
            "degree_three_four_sample_count": arguments.sample_count,
            "height_slack": arguments.height_slack,
            "seed": arguments.seed,
        },
        "inputs": {
            relative(database_path): digest(database_path),
            relative(PINNED_R17): digest(PINNED_R17),
        },
        "targets": rows,
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/sample_lattice_foundry_multisection_spectrum.sage "
            f"--sample-count {arguments.sample_count} --height-slack {arguments.height_slack}"
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("multisection-spectrum artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        f"FOUNDRYMULTI|targets={len(rows)}|"
        f"published_controls={sum(row['published_R17_isometry_control'] for row in rows)}|"
        "status=PASS",
        flush=True,
    )


if __name__ == "__main__":
    main()
