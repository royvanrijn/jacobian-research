#!/usr/bin/env sage-python
"""Measure low-degree, low-genus multisection richness of marked targets.

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
census. For a primitive A1 target, the finite quotient is divided further by
its exact root reflection. Except for section-nonnegative rational bisections
on a rootless fibration, irreducibility and arithmetic descent remain separate
gates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
from collections import Counter
from itertools import product
from pathlib import Path

from sage.all import QQ, ZZ, block_diagonal_matrix, ceil, gcd, identity_matrix, matrix, pari, vector


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
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


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


def roots_and_data(gram):
    result = pari(gram).qfminim(2)
    count = ZZ(result[0])
    if count == 0:
        return (), matrix(ZZ, 0, gram.nrows()), (0, 0, 1)
    half = [vector(ZZ, column) for column in matrix(ZZ, result[2]).columns()]
    roots = tuple(half + [-root for root in half])
    root_basis = matrix(ZZ, [list(root) for root in roots]).row_module().basis_matrix()
    root_gram = root_basis * gram * root_basis.transpose()
    return roots, root_basis, (root_basis.rank(), count, abs(ZZ(root_gram.det())))


def root_adapt_a1(gram):
    """Put the primitive A1 root first and LLL-reduce its MW quotient."""

    roots, root_basis, root_data = roots_and_data(gram)
    assert root_data == (1, 2, 2) and len(roots) == 2
    smith, smith_left, smith_right = root_basis.smith_form()
    assert smith == smith_left * root_basis * smith_right
    assert abs(ZZ(smith[0, 0])) == 1
    completion = smith_right.inverse()
    basis = root_basis.stack(completion[1:])
    assert abs(basis.det()) == 1
    adapted = basis * gram * basis.transpose()
    assert adapted[0, 0] == 2

    coupling = adapted[:1, 1:]
    tail = adapted[1:, 1:]
    height = tail - coupling.transpose() * coupling / 2
    scale = ZZ(2)
    quotient_lll = matrix(ZZ, pari((scale * height).change_ring(ZZ)).qflllgram())
    change = block_diagonal_matrix(identity_matrix(ZZ, 1), quotient_lll.transpose())
    basis = change * basis
    adapted = basis * gram * basis.transpose()
    assert abs(basis.det()) == 1 and adapted[0, 0] == 2
    return adapted, basis


def a1_reflected_residue(residue, gram, degree):
    """Reflect a coordinate residue in the first (A1) root."""

    tail = residue[1:]
    reflected_zero = -residue[0] - sum(
        gram[0, index] * tail[index - 1]
        for index in range(1, gram.nrows())
    )
    return tuple(
        [int(reflected_zero) % degree]
        + [int(entry) % degree for entry in tail]
    )


def a1_orbit_key(residue, gram, degree):
    literal = tuple(int(entry) % degree for entry in residue)
    return min(literal, a1_reflected_residue(literal, gram, degree))


def a1_quotient_orbit_count(gram, degree):
    """Burnside count for (M/dM)/W(A1)."""

    coefficients = [ZZ(2)] + [ZZ(gram[0, index]) for index in range(1, 17)]
    common = ZZ(degree)
    for coefficient in coefficients:
        common = gcd(common, coefficient)
    fixed = degree**16 * int(common)
    return (degree**17 + fixed) // 2, fixed


def exact_degree_two_spectrum_a1(gram, orbit_cap=None):
    """Complete A1-Weyl quotient minima through norm ten."""

    lower, diagonal = ldl_data(gram)
    minima = {}
    total_orbits, fixed_cosets = a1_quotient_orbit_count(gram, 2)
    tested = 0
    for residue_tuple in product(range(2), repeat=gram.nrows()):
        key = a1_orbit_key(residue_tuple, gram, 2)
        if residue_tuple != key:
            continue
        minimum = minimum_coset_norm_through_bound(
            gram, lower, diagonal, 2, vector(ZZ, residue_tuple), 10
        )
        if minimum is not None:
            minima[key] = minimum
        tested += 1
        if tested % 8192 == 0:
            print(
                f"FOUNDRYMULTI_A1_D2|orbits={tested}/{total_orbits}|"
                f"covered={len(minima)}|status=RUNNING",
                flush=True,
            )
        if orbit_cap is not None and tested >= orbit_cap:
            break
    complete = tested == total_orbits

    def category(genus, threshold):
        histogram = Counter(
            minimum
            for key, minimum in minima.items()
            if minimum >= threshold
            and (minimum + 2 * genus - 2) % 4 == 0
        )
        return {
            "arithmetic_genus": genus,
            "section_nonnegative_threshold_norm": threshold,
            "low_height_norm_cap": 10,
            "weyl_translation_orbits_in_cap": sum(histogram.values()),
            "minimum_norm_histogram": {
                str(norm): count for norm, count in sorted(histogram.items())
            },
            "status": (
                "PASS_EXACT_A1_WEYL_TRANSLATION_ORBITS_THROUGH_NORM_TEN"
                if complete
                else "PASS_BOUNDED_A1_WEYL_TRANSLATION_ORBIT_PREFIX"
            ),
        }

    return {
        "method": "complete exact rational-LDL CVP over every (M/2M)/W(A1) orbit",
        "enumerated_through_norm": 10,
        "total_weyl_translation_orbits": total_orbits,
        "tested_weyl_translation_orbits": tested,
        "complete_orbit_census": complete,
        "reflection_fixed_cosets": fixed_cosets,
        "orbits_with_minimum_at_most_cap": len(minima),
        "orbits_with_minimum_above_cap": (
            total_orbits - len(minima) if complete else None
        ),
        "rational_bisections": category(0, 10),
        "genus_one_bisections": category(1, 8),
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


def sample_degree_a1(gram, frame_id, degree, genera, sample_count, slack, seed):
    """Sample (M/dM)/W(A1), using an exact CVP for every sampled orbit."""

    lower, diagonal = ldl_data(gram)
    dimension = gram.nrows()
    frame_seed = int(hashlib.sha256(frame_id.encode()).hexdigest()[:16], 16)
    rng = random.Random(seed ^ frame_seed ^ (degree << 24) ^ 0xA1)
    residues = {(0,) * dimension}
    while len(residues) < sample_count:
        raw = tuple(rng.randrange(degree) for unused in range(dimension))
        residues.add(a1_orbit_key(raw, gram, degree))

    thresholds = {genus: 2 * degree * degree - 2 * genus + 2 for genus in genera}
    bound = max(thresholds.values()) + slack
    categories = {
        genus: {
            "arithmetic_genus": genus,
            "section_nonnegative_threshold_norm": thresholds[genus],
            "low_height_norm_cap": thresholds[genus] + slack,
            "sampled_qualifying_weyl_translation_orbits": 0,
            "minimum_norm_histogram": Counter(),
        }
        for genus in genera
    }
    minima_above_bound = 0
    minima = Counter()
    for residue_tuple in sorted(residues):
        residue = vector(ZZ, residue_tuple)
        minimum = minimum_coset_norm_through_bound(
            gram, lower, diagonal, degree, residue, bound
        )
        if minimum is None:
            minima_above_bound += 1
            continue
        minima[minimum] += 1
        residue_norm = int(residue * gram * residue)
        for genus in genera:
            if (residue_norm + 2 * genus - 2) % (2 * degree):
                continue
            threshold = thresholds[genus]
            if threshold <= minimum <= threshold + slack:
                row = categories[genus]
                row["sampled_qualifying_weyl_translation_orbits"] += 1
                row["minimum_norm_histogram"][minimum] += 1

    for row in categories.values():
        row["minimum_norm_histogram"] = {
            str(key): value for key, value in sorted(row["minimum_norm_histogram"].items())
        }
        row["sample_fraction_of_sampled_orbits"] = (
            row["sampled_qualifying_weyl_translation_orbits"] / len(residues)
        )
        row["geometric_boundary"] = (
            "Exact for the sampled A1 Weyl/section-translation orbits. The lattice "
            "minimum and adjunction congruence are certified; effectivity, "
            "irreducibility, arithmetic descent, and rank gain are not asserted."
        )
    total_orbits, fixed_cosets = a1_quotient_orbit_count(gram, degree)
    return {
        "method": (
            "deterministic uniform (M/dM)/W(A1) orbit sample with exact "
            "rational-LDL CVP bound"
        ),
        "degree": degree,
        "total_weyl_translation_orbits": total_orbits,
        "reflection_fixed_cosets": fixed_cosets,
        "sample_size": len(residues),
        "seed": seed,
        "norm_bound": bound,
        "sampled_orbits_with_minimum_above_bound": minima_above_bound,
        "sampled_minimum_norm_histogram": {
            str(key): value for key, value in sorted(minima.items())
        },
        "categories": [categories[genus] for genus in genera],
    }


def target_rows(database):
    if "ns_classes" in database:
        frame_by_id = {
            frame["frame_id"]: frame
            for ns in database["ns_classes"]
            for frame in ns["frames"]
        }
        return [
            (target, frame_by_id[target["frame_id"]])
            for target in database["rootless_targets"]
        ]
    if "surfaces" in database:
        return [
            (
                {
                    "frame_id": frame["frame_id"],
                    "ns_id": surface["surface_id"],
                    "determinant": int(frame["determinant"]),
                    "is_existing_H3_control": False,
                },
                frame,
            )
            for surface in database["surfaces"]
            for frame in surface["frames"]
            if int(frame["root_rank"]) == 0
            and int(frame["mw_rank_for_rho_19"]) == 17
        ]
    raise ValueError("database has neither ns_classes nor surfaces")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DATABASE)
    parser.add_argument(
        "--target-artifact",
        type=Path,
        help="an extracted single-frame artifact; supports rootless or primitive A1 targets",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--frame-id", action="append", default=[])
    parser.add_argument("--limit", type=int)
    parser.add_argument("--sample-count", type=int, default=256)
    parser.add_argument("--height-slack", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--pari-stack-gb", type=int, default=4)
    parser.add_argument(
        "--a1-degree-two-orbit-cap",
        type=int,
        help="diagnostic prefix cap for the otherwise complete A1 degree-two census",
    )
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.sample_count <= 0 or arguments.height_slack < 0:
        parser.error("sample count must be positive and height slack nonnegative")
    if (
        arguments.a1_degree_two_orbit_cap is not None
        and arguments.a1_degree_two_orbit_cap <= 0
    ):
        parser.error("A1 degree-two orbit cap must be positive")
    if arguments.pari_stack_gb <= 0:
        parser.error("PARI stack size must be positive")
    pari.allocatemem(arguments.pari_stack_gb * 1024**3)

    database_path = arguments.database.resolve()
    if arguments.target_artifact is None:
        database = json.loads(database_path.read_text())
        selected = target_rows(database)
    else:
        target_path = arguments.target_artifact.resolve()
        target_payload = json.loads(target_path.read_text())
        frame = target_payload["frame"]
        selected = [
            (
                {
                    "frame_id": frame["frame_id"],
                    "ns_id": target_payload.get("surface_id", frame["frame_id"]),
                    "determinant": int(frame["determinant"]),
                    "is_existing_H3_control": False,
                },
                frame,
            )
        ]
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
        root_rank = int(frame.get("root_rank", roots_and_data(gram)[2][0]))
        assert gram.nrows() == 17 and root_rank in (0, 1)
        is_published = root_rank == 0 and pari(gram).qfisom(pari(pinned)) != 0
        if root_rank == 0:
            spectrum_gram = gram
            degree_two = exact_degree_two_spectrum(spectrum_gram)
            degree_three = sample_degree(
                spectrum_gram, target["frame_id"], 3, (0, 1),
                arguments.sample_count, arguments.height_slack, arguments.seed,
            )
            degree_four = sample_degree(
                spectrum_gram, target["frame_id"], 4, (0, 1, 2),
                arguments.sample_count, arguments.height_slack, arguments.seed,
            )
            d3_count_key = "sampled_qualifying_translation_cosets"
            d4_count_key = d3_count_key
            quotient_kind = "section_translation"
            root_adapted_gram = None
        else:
            spectrum_gram, root_adapted_basis = root_adapt_a1(gram)
            degree_two = exact_degree_two_spectrum_a1(
                spectrum_gram, arguments.a1_degree_two_orbit_cap
            )
            degree_three = sample_degree_a1(
                spectrum_gram, target["frame_id"], 3, (0, 1),
                arguments.sample_count, arguments.height_slack, arguments.seed,
            )
            degree_four = sample_degree_a1(
                spectrum_gram, target["frame_id"], 4, (0, 1, 2),
                arguments.sample_count, arguments.height_slack, arguments.seed,
            )
            d3_count_key = "sampled_qualifying_weyl_translation_orbits"
            d4_count_key = d3_count_key
            quotient_kind = "section_translation_and_A1_Weyl"
            root_adapted_gram = [
                [int(entry) for entry in row] for row in spectrum_gram.rows()
            ]
        spectrum = {
            "rational_bisection_orbits_minimum_height": degree_two[
                "rational_bisections"
            ].get(
                "translation_orbits_with_minimum_norm_ten",
                degree_two["rational_bisections"].get("weyl_translation_orbits_in_cap"),
            ),
            "genus_one_bisection_candidate_orbits_minimum_height": degree_two[
                "genus_one_bisections"
            ].get(
                "translation_orbits_with_minimum_norm_eight",
                degree_two["genus_one_bisections"].get("weyl_translation_orbits_in_cap"),
            ),
            "sampled_rational_trisection_candidates": degree_three[
                "categories"
            ][0][d3_count_key],
            "sampled_genus_one_trisection_candidates": degree_three[
                "categories"
            ][1][d3_count_key],
            "sampled_low_genus_quadrisection_candidates": sum(
                row[d4_count_key]
                for row in degree_four["categories"]
            ),
        }
        rows.append(
            {
                "frame_id": target["frame_id"],
                "ns_id": target["ns_id"],
                "determinant": target["determinant"],
                "root_rank": root_rank,
                "quotient_kind": quotient_kind,
                "root_adapted_gram": root_adapted_gram,
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

    output_path = arguments.output.resolve()
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
            "rootful_A1_quotient": "(M/dM)/W(A1)",
        },
        "search": {
            "selected_frame_ids": [row[0]["frame_id"] for row in selected],
            "limit": arguments.limit,
            "degree_three_four_sample_count": arguments.sample_count,
            "height_slack": arguments.height_slack,
            "seed": arguments.seed,
        },
        "inputs": {
            relative(
                arguments.target_artifact.resolve()
                if arguments.target_artifact is not None
                else database_path
            ): digest(
                arguments.target_artifact.resolve()
                if arguments.target_artifact is not None
                else database_path
            ),
            relative(PINNED_R17): digest(PINNED_R17),
        },
        "targets": rows,
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/sample_lattice_foundry_multisection_spectrum.sage "
            + (
                f"--target-artifact {relative(arguments.target_artifact)} "
                if arguments.target_artifact is not None
                else ""
            )
            + (
                f"--database {relative(database_path)} "
                if arguments.target_artifact is None and database_path != DATABASE.resolve()
                else ""
            )
            + (
                f"--output {relative(output_path)} "
                if output_path != DEFAULT_OUTPUT.resolve()
                else ""
            )
            + (
                f"--a1-degree-two-orbit-cap {arguments.a1_degree_two_orbit_cap} "
                if arguments.a1_degree_two_orbit_cap is not None
                else ""
            )
            +
            f"--sample-count {arguments.sample_count} --height-slack {arguments.height_slack} "
            + " ".join(
                f"--frame-id {row[0]['frame_id']}" for row in selected
            )
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
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
