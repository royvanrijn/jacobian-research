#!/usr/bin/env sage-python
"""Classify the bounded alternate-Q80 bisection arithmetic laboratory.

The laboratory combines the exact 121 inherited covers with a priority prefix
of native alternate-Q80 bisections.  Source artifacts use different parity
bases, so deduplication is performed only after transporting every record to
the direct compiled alternate frame.  The script then checks conic solubility,
geometric branch incidence, squareclass collisions, pair-product characters,
and independent triples exactly.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from functools import reduce
from math import gcd
from pathlib import Path

from sage.all import Conic, PolynomialRing, QQ, ZZ, lcm, matrix, vector
from sage.env import SAGE_VERSION


ROOT = Path(__file__).resolve().parents[2]
INHERITED = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-inherited-bisection-covers-v1.json"
INHERITED_PRODUCTS = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-inherited-product-characters-v1.json"
NATIVE = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisections-cheapest-1024-v1.json"
NATIVE_SQUARECLASSES = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-bisection-collisions-cheapest-1024-v1.json"
DIRECT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json"
OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-r17-norm12-11952-alternate-arithmetic-laboratory-cheapest-1024-v1.json"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def rational_text(value) -> str:
    value = QQ(value)
    return str(value.numerator()) if value.denominator() == 1 else f"{value.numerator()}/{value.denominator()}"


def primitive_quadratic(coefficients) -> tuple[ZZ, ZZ, ZZ]:
    values = tuple(QQ(value) for value in coefficients)
    if len(values) != 3 or not values[2]:
        raise ArithmeticError("branch equation is not a genuine quadratic")
    denominator = lcm(value.denominator() for value in values)
    integers = [ZZ(value * denominator) for value in values]
    content = reduce(gcd, (abs(int(value)) for value in integers))
    integers = tuple(value // content for value in integers)
    if integers[2] < 0:
        integers = tuple(-value for value in integers)
    return integers


def rational_square(value) -> bool:
    value = QQ(value)
    return value >= 0 and value.numerator().is_square() and value.denominator().is_square()


def parity_mask(entries) -> int:
    return sum((int(value) % 2) << index for index, value in enumerate(entries))


def choose(count: int, size: int) -> int:
    if size == 2:
        return count * (count - 1) // 2
    if size == 3:
        return count * (count - 1) * (count - 2) // 6
    raise ValueError("only pair and triple counts are used")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inherited", type=Path, default=INHERITED)
    parser.add_argument("--inherited-products", type=Path, default=INHERITED_PRODUCTS)
    parser.add_argument("--native", type=Path, default=NATIVE)
    parser.add_argument("--native-squareclasses", type=Path, default=NATIVE_SQUARECLASSES)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    inherited = json.loads(arguments.inherited.read_text())
    inherited_products = json.loads(arguments.inherited_products.read_text())
    native = json.loads(arguments.native.read_text())
    native_squareclasses = json.loads(arguments.native_squareclasses.read_text())
    direct = json.loads(DIRECT.read_text())

    if inherited.get("status") != "PASS_EXACT_121_INHERITED_ALTERNATE_Q80_BISECTION_COVERS":
        raise ArithmeticError("the inherited input is not the exact 121-cover certificate")
    if len(inherited.get("bisections", [])) != 121:
        raise ArithmeticError("the inherited certificate does not contain 121 covers")
    inherited_group = inherited_products.get("inherited_character_group", {})
    if (
        inherited_products.get("status") != "PASS_EXACT_NO_INHERITED_PRODUCT_CHARACTER_CLOSURE"
        or inherited_group.get("pair_product_count") != 7260
        or inherited_group.get("distinct_pair_product_count") != 7260
        or inherited_group.get("matches_another_inherited_character_count") != 0
    ):
        raise ArithmeticError("the exact 7,260 inherited-product certificate changed")
    interval = native.get("interval", {})
    if (
        native.get("status") != "PASS_EXACT_ALTERNATE_BISECTION_EQUATION_CHUNK"
        or interval != {"start_zero_based": 0, "stop_exclusive": 1024}
        or len(native.get("bisections", [])) != 1024
        or not native["construction"]["all_branch_fibres_smooth"]
        or not native["construction"]["all_lifted_sections_verified"]
    ):
        raise ArithmeticError("the native input is not the exact cheapest-1024 chunk")
    if (
        native_squareclasses.get("status") != "PASS_EXTENSION_CANONICALIZATION"
        or native_squareclasses.get("input_bisection_count") != 1024
        or native_squareclasses.get("distinct_quadratic_extensions") != 1024
        or native_squareclasses.get("collision_count") != 0
        or native_squareclasses.get("input", {}).get("sha256") != digest(arguments.native)
    ):
        raise ArithmeticError("the native 1,024-squareclass certificate changed")

    frame = matrix(ZZ, direct["frame_certificate"]["frame_gram"])
    if frame.det() != 948:
        raise ArithmeticError("the direct alternate frame determinant changed")

    entries = []
    for source, payload in (("inherited", inherited), ("native", native)):
        for record in payload["bisections"]:
            if source == "inherited":
                coefficients = record["canonical_squareclass"]["q_coefficients_low_to_high"]
                direct_vector = vector(ZZ, record["alternate_rank17_w"])
                if not record["lifted_section"]["equation_verified"]:
                    raise ArithmeticError("an inherited lifted section is not verified")
            else:
                coefficients = record["branch"]["numerator_coefficients"]
                direct_vector = vector(ZZ, record["direct_alternate_w"])
                if (
                    not record["lifted_section"]["constant_and_linear_identities_verified"]
                    or not record["lifted_section"]["two_branches_verified"]
                    or not record["residual_chord_certificate"]["branch_fibres_smooth"]
                ):
                    raise ArithmeticError("a native lifted-section certificate is incomplete")
            if direct_vector * frame * direct_vector != 10:
                raise ArithmeticError("a laboratory bisection does not have norm ten")
            q = tuple(QQ(value) for value in coefficients)
            primitive = primitive_quadratic(q)
            discriminant = primitive[1] ** 2 - 4 * primitive[2] * primitive[0]
            if not discriminant or rational_square(discriminant):
                raise ArithmeticError("a laboratory branch quadratic is reducible or singular")
            if int(record["lifted_section"]["anti_invariant_height"]) != 12:
                raise ArithmeticError("a laboratory anti-invariant height changed")
            entries.append({
                "source": source,
                "label": str(record["label"]),
                "priority_rank": record.get("priority_rank"),
                "direct_vector": direct_vector,
                "direct_mask": parity_mask(direct_vector),
                "q": q,
                "primitive": primitive,
            })

    by_direct_mask = {}
    for entry in entries:
        by_direct_mask.setdefault(entry["direct_mask"], []).append(entry)
    overlaps = []
    representatives = []
    for direct_mask, group in sorted(by_direct_mask.items()):
        inherited_group_records = [entry for entry in group if entry["source"] == "inherited"]
        representative = (inherited_group_records or group)[0]
        representatives.append(representative)
        if len(group) == 1:
            continue
        if len(group) != 2 or {entry["source"] for entry in group} != {"inherited", "native"}:
            raise ArithmeticError("unexpected within-source or higher-multiplicity orbit overlap")
        left, right = group
        if left["primitive"] != right["primitive"]:
            raise ArithmeticError("one direct translation class has two branch divisors")
        ratio = left["q"][2] / right["q"][2]
        if not rational_square(ratio):
            raise ArithmeticError("one direct translation class has two squareclasses")
        difference = left["direct_vector"] - right["direct_vector"]
        if any(value % 2 for value in difference):
            raise ArithmeticError("purported source overlap is not a section translation")
        overlaps.append({
            "direct_alternate_parity_mask": f"0x{direct_mask:05x}",
            "records": [
                {
                    "source": entry["source"],
                    "label": entry["label"],
                    "source_stored_lattice_orbit_mask": next(
                        record["lattice_orbit_mask"]
                        for payload in (inherited, native)
                        for record in payload["bisections"]
                        if record["label"] == entry["label"]
                    ),
                }
                for entry in group
            ],
            "same_primitive_branch_divisor": True,
            "squareclass_ratio": rational_text(ratio),
            "squareclass_ratio_is_rational_square": True,
            "direct_vector_difference_is_in_2M": True,
        })

    branch_buckets = {}
    for entry in representatives:
        branch_buckets.setdefault(entry["primitive"], []).append(entry)
    shared_distinct_class_branches = [
        group for group in branch_buckets.values() if len(group) > 1
    ]
    if shared_distinct_class_branches:
        raise ArithmeticError("distinct direct translation classes share a branch divisor")

    # Each branch polynomial is irreducible over QQ.  Therefore two distinct
    # primitive quadratics cannot share even one geometric root: a common root
    # would have the same degree-two minimal polynomial.  Unique factorization
    # also makes all unordered pair products distinct and prevents a quartic
    # pair product from matching any degree-two catalogued character.
    unique_count = len(representatives)
    pair_count = choose(unique_count, 2)
    triple_count = choose(unique_count, 3)

    P = PolynomialRing(QQ, names=("U", "S", "Z"))
    U, S, Z = P.gens()
    conic_points = []
    for position, entry in enumerate(representatives, start=1):
        q0, q1, q2 = entry["q"]
        conic = Conic(S**2 - (q2 * U**2 + q1 * U * Z + q0 * Z**2))
        soluble, point = conic.has_rational_point(point=True)
        point_record = None
        if soluble:
            coordinates = tuple(QQ(value) for value in point)
            if coordinates[1] ** 2 != (
                q2 * coordinates[0] ** 2
                + q1 * coordinates[0] * coordinates[2]
                + q0 * coordinates[2] ** 2
            ):
                raise ArithmeticError("Hasse--Minkowski conic point failed substitution")
            point_record = [rational_text(value) for value in coordinates]
        conic_points.append({
            "source": entry["source"],
            "label": entry["label"],
            "direct_alternate_parity_mask": f"0x{entry['direct_mask']:05x}",
            "has_Q_point": bool(soluble),
            "point_U_S_Z": point_record,
        })
        if position == 1 or position % 100 == 0 or position == unique_count:
            print(
                f"ALTLABCONIC|completed={position}/{unique_count}|label={entry['label']}|Qpoint={int(soluble)}",
                flush=True,
            )
    rational_conics = sum(record["has_Q_point"] for record in conic_points)

    native_ranked = sorted(
        (entry for entry in entries if entry["source"] == "native"),
        key=lambda entry: int(entry["priority_rank"]),
    )
    cheapest_independent_triple = [
        {
            "label": entry["label"],
            "priority_rank": int(entry["priority_rank"]),
            "direct_alternate_parity_mask": f"0x{entry['direct_mask']:05x}",
        }
        for entry in native_ranked[:3]
    ]

    result = {
        "schema": "elkies-k3.r17-norm12-11952-alternate-arithmetic-laboratory.v1",
        "status": "PASS_EXACT_BOUNDED_ALTERNATE_ARITHMETIC_LABORATORY",
        "scope": {
            "inherited_records": 121,
            "native_priority_interval": [0, 1024],
            "raw_record_count": len(entries),
            "unique_direct_translation_class_count": unique_count,
        },
        "coordinate_aware_deduplication": {
            "common_coordinate_system": "direct compiled alternate-Q80 frame modulo 2M",
            "source_overlap_count": len(overlaps),
            "source_overlaps": overlaps,
            "deduplicated_record_count": len(entries) - unique_count,
            "consequence": (
                "These are repeated realizations of the same translation classes, not "
                "independent sections on colliding quadratic covers."
            ),
        },
        "branch_and_squareclass_classification": {
            "irreducible_squarefree_quadratic_count": unique_count,
            "distinct_geometric_branch_divisor_count": len(branch_buckets),
            "distinct_quadratic_squareclass_count": unique_count,
            "shared_branch_fibres_between_distinct_classes": 0,
            "shared_branch_divisors_between_distinct_classes": 0,
            "squareclass_collisions_between_distinct_classes": 0,
            "proof": (
                "Every primitive branch quadratic is irreducible over QQ and the primitive "
                "quadratics are pairwise distinct after direct-frame orbit deduplication. "
                "A common geometric root would force equality of their degree-two minimal polynomials."
            ),
        },
        "individual_cover_arithmetic": {
            "cover_count": unique_count,
            "Q_rational_conic_count": rational_conics,
            "anisotropic_conic_count": unique_count - rational_conics,
            "all_conics_Q_rational": rational_conics == unique_count,
            "proof_method": "exact Sage Hasse--Minkowski conic solver plus point substitution",
            "points": conic_points,
            "generic_rank_lower_bound_on_each_Q_rational_cover": 18,
        },
        "pair_and_product_search": {
            "unordered_distinct_pair_count": pair_count,
            "distinct_pair_product_squareclass_count": pair_count,
            "pair_product_matches_catalogued_section_character_count": 0,
            "v4_nontrivial_character_triple_count": 0,
            "shared_branch_genus_zero_pair_count": 0,
            "disjoint_branch_genus_one_pair_count": pair_count,
            "minimum_connected_pair_fibre_product_genus": 1,
            "pair_cover_anti_invariant_height_matrix": [[24, 0], [0, 24]],
            "pair_cover_generic_rank_lower_bound": 19,
            "product_twist_boundary": (
                "The exact search tests whether a pair product equals any compiled "
                "degree-two section character. It does not compute Mordell--Weil groups of "
                "the uncatalogued quartic product twists, so undiscovered sections there remain UNKNOWN."
            ),
        },
        "independent_character_triples": {
            "count": triple_count,
            "all_catalogue_triples_F2_independent": True,
            "common_cover_galois_group": "(Z/2Z)^3",
            "geometric_branch_point_count": 6,
            "common_cover_genus": 5,
            "riemann_hurwitz": "2g-2=8*(-2)+6*4=8",
            "pulled_anti_invariant_height_matrix": [[48, 0, 0], [0, 48, 0], [0, 0, 48]],
            "generic_rank_lower_bound": 20,
            "cheapest_native_example": cheapest_independent_triple,
            "target_boundary": (
                "These are three independent characters on an eightfold genus-five cover. "
                "They are not the three dependent nontrivial characters of one V4 cover."
            ),
        },
        "foundry_branch_incidence_score": {
            "unique_classes": unique_count,
            "shared_branch_fibre_pairs": 0,
            "genus_zero_pair_bases": 0,
            "catalogued_product_character_sections": 0,
            "v4_three_character_closures": 0,
            "interpretation": (
                "The cheapest-1024 alternate prefix adds no controlled branch incidence "
                "after the inherited/native source overlap is removed."
            ),
        },
        "inputs": {
            relative(path): digest(path)
            for path in (
                arguments.inherited,
                arguments.inherited_products,
                arguments.native,
                arguments.native_squareclasses,
                DIRECT,
            )
        },
        "software_assumptions": {
            "sage_version": SAGE_VERSION,
            "required_features": [
                "exact rational arithmetic",
                "exact Hasse--Minkowski conic solver",
            ],
        },
        "reproducing_command": (
            "sage -python elkies-k3/scripts/"
            "analyze_r17_norm12_11952_alternate_laboratory.sage"
        ),
        "proof_boundary": (
            "This is an exact classification of the 121 inherited records together with "
            "only the cheapest 1,024 native alternate classes. It proves the displayed "
            "deduplication, branch-incidence, pair-product, conic, and character-decomposition "
            "statements inside that bounded laboratory. It makes no claim about the remaining "
            "38,123 native classes or about sections on uncatalogued quartic product twists."
        ),
    }
    serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if arguments.check:
        if not arguments.output.exists() or arguments.output.read_text() != serialized:
            raise ArithmeticError("stored alternate laboratory artifact differs from exact replay")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(serialized)
    print(
        "ALTLAB|raw={}|classes={}|overlaps={}|Qconics={}|pairs={}|"
        "shared_branch=0|product_sections=0|v4_closures=0|triples={}|status={}".format(
            len(entries), unique_count, len(overlaps), rational_conics, pair_count,
            triple_count, result["status"],
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
