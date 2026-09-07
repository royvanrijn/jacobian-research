#!/usr/bin/env sage
"""Certify the rank-7 auxiliary lattice and Niemeier control embeddings.

status: ACTIVE_PROOF
claim: The pinned rank-7 lattice has discriminant form opposite to the two
  rootless J2 controls.  Every discriminant anti-isometry glues either control
  to the Niemeier lattice N(2A7+2D5), and the primitive complement recovers the
  requested rank-17 frame exactly up to integral isometry.
inputs: elkies-k3/data/lattice/rootless_j2_auxiliary_rank7_gram.txt,
  elkies-k3/data/lattice/rank17_gram.txt,
  artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json,
  elkies-k3/seeds/**/*.txt
outputs: artifacts/generated-results/elkies-k3-rootless-j2-niemeier-controls.json
supersedes/superseded-by: none
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

from sage.all import (
    Genus,
    IntegralLattice,
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    matrix,
    pari,
    vector,
)
from sage.quadratic_forms.genera.genus import genera


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
AUXILIARY = ROOT / "elkies-k3/data/lattice/rootless_j2_auxiliary_rank7_gram.txt"
PUBLISHED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
ALTERNATE = (
    ROOT
    / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
)
SEEDS = ROOT / "elkies-k3/seeds"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rootless-j2-niemeier-controls.json"
)


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def load_seed_gram(path):
    lines = Path(path).read_text().splitlines()
    try:
        start = next(
            index for index, line in enumerate(lines) if line.strip() == "gram ="
        ) + 1
    except StopIteration:
        return None
    result = []
    for line in lines[start:]:
        stripped = line.strip()
        if not stripped:
            break
        values = [
            ZZ(value)
            for value in stripped.replace("[", " ").replace("]", " ").split()
        ]
        if values:
            result.append(values)
    try:
        return matrix(ZZ, result)
    except (TypeError, ValueError):
        return None


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def gram_sha256(value):
    payload = "\n".join(" ".join(map(str, row)) for row in value.rows()) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def discriminant_data(gram):
    group = IntegralLattice(gram).discriminant_group()
    assert tuple(map(int, group.invariants())) == (948,)
    return group, group.gram_matrix_quadratic()[0, 0]


def anti_isometry_units(left_q, right_q):
    return [
        unit
        for unit in range(948)
        if math.gcd(unit, 948) == 1
        and ((left_q * unit * unit + right_q) / 2).denominator() == 1
    ]


def automorphism_discriminant_image(gram):
    group, unused = discriminant_data(gram)
    generator = vector(QQ, group.gen(0).lift())
    pari_automorphisms = pari(gram).qfauto()
    image_generators = []
    for raw in pari_automorphisms[1].sage():
        automorphism = matrix(ZZ, raw)
        assert automorphism.transpose() * gram * automorphism == gram
        transformed = generator * automorphism.transpose()
        units = [
            unit
            for unit in range(948)
            if all(
                (transformed[index] - unit * generator[index]).denominator() == 1
                for index in range(gram.nrows())
            )
        ]
        assert len(units) == 1
        image_generators.append(units[0])

    image = {1}
    while True:
        enlarged = set(image)
        for value in image:
            enlarged.update((value * unit) % 948 for unit in image_generators)
        if enlarged == image:
            break
        image = enlarged
    return int(pari_automorphisms[0]), sorted(image)


def double_cosets(units, left_image, right_image):
    remaining = set(units)
    result = []
    while remaining:
        representative = min(remaining)
        orbit = {
            (representative * left * right) % 948
            for left in left_image
            for right in right_image
        } & set(units)
        assert orbit
        result.append(sorted(orbit))
        remaining.difference_update(orbit)
    return result


def root_components(gram):
    minimum_data = pari(gram).qfminim(2)
    signed_count = int(minimum_data[0])
    representatives = matrix(ZZ, minimum_data[2].sage()).transpose()
    assert 2 * representatives.nrows() == signed_count

    unseen = set(range(representatives.nrows()))
    components = []
    while unseen:
        component = {min(unseen)}
        frontier = list(component)
        unseen.difference_update(component)
        while frontier:
            current = frontier.pop()
            neighbours = {
                index
                for index in unseen
                if representatives[current] * gram * representatives[index] != 0
            }
            component.update(neighbours)
            unseen.difference_update(neighbours)
            frontier.extend(neighbours)
        basis = representatives[sorted(component)]
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


def integral_matrix(value):
    assert all(entry.denominator() == 1 for entry in value.list())
    return matrix(ZZ, value)


def primitive_row_lattice(value):
    diagonal, unused_left, unused_right = value.smith_form()
    return all(abs(diagonal[index, index]) == 1 for index in range(value.nrows()))


def glue(auxiliary, frame, unit):
    auxiliary_group, unused = discriminant_data(auxiliary)
    frame_group, unused = discriminant_data(frame)
    auxiliary_generator = vector(QQ, auxiliary_group.gen(0).lift())
    frame_generator = vector(QQ, frame_group.gen(0).lift())

    split_gram = block_diagonal_matrix(auxiliary, frame)
    split_lattice = IntegralLattice(split_gram)
    glue_vector = vector(
        QQ,
        list(unit * auxiliary_generator) + list(frame_generator),
    )
    ambient = split_lattice.overlattice([glue_vector])
    assert ambient.rank() == 24
    assert ambient.gram_matrix().det() == 1

    old_gram = ambient.gram_matrix()
    change_columns = matrix(ZZ, pari(old_gram).qflllgram())
    change_rows = change_columns.transpose()
    assert abs(change_rows.det()) == 1
    ambient_gram = integral_matrix(
        change_rows * old_gram * change_rows.transpose()
    )
    ambient_basis = change_rows * ambient.basis_matrix()
    assert ambient_gram == ambient_basis * split_gram * ambient_basis.transpose()
    assert all(ambient_gram[index, index] % 2 == 0 for index in range(24))

    split_auxiliary_basis = identity_matrix(ZZ, 24)[:7]
    auxiliary_coordinates = integral_matrix(
        split_auxiliary_basis * ambient_basis.inverse()
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
    complement_gram = (
        complement_coordinates
        * ambient_gram
        * complement_coordinates.transpose()
    )
    assert complement_gram.det() == 948
    assert int(pari(complement_gram).qfminim(2)[0]) == 0
    assert Genus(complement_gram) == Genus(frame)
    assert pari(frame).qfisom(pari(complement_gram)) != 0

    signed_roots, components = root_components(ambient_gram)
    assert signed_roots == 192
    assert components == [
        {"rank": 7, "signed_root_count": 56},
        {"rank": 7, "signed_root_count": 56},
        {"rank": 5, "signed_root_count": 40},
        {"rank": 5, "signed_root_count": 40},
    ]

    return {
        "anti_isometry_unit": unit,
        "ambient_root_system": "2A7+2D5",
        "ambient_signed_root_count": signed_roots,
        "ambient_root_components": components,
        "ambient_gram_sha256": gram_sha256(ambient_gram),
        "complement_gram_sha256": gram_sha256(complement_gram),
        "ambient_gram": rows(ambient_gram),
        "auxiliary_basis_in_ambient": rows(auxiliary_coordinates),
        "complement_basis_in_ambient": rows(complement_coordinates),
        "complement_gram": rows(complement_gram),
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

auxiliary = load_matrix(AUXILIARY)
published = load_matrix(PUBLISHED)
alternate_payload = json.loads(ALTERNATE.read_text())
alternate = matrix(ZZ, alternate_payload["rootless_frame"])
unused_published_group, published_q = discriminant_data(published)

assert auxiliary.nrows() == auxiliary.ncols() == 7
assert auxiliary.is_positive_definite()
assert auxiliary.det() == 948
assert all(auxiliary[index, index] % 2 == 0 for index in range(7))

auxiliary_group, auxiliary_q = discriminant_data(auxiliary)
assert auxiliary_q == QQ(1267) / 948
auxiliary_roots, auxiliary_components = root_components(auxiliary)
assert auxiliary_roots == 40
assert auxiliary_components == [{"rank": 5, "signed_root_count": 40}]
auxiliary_mass = Genus(auxiliary).mass()
assert auxiliary_mass == QQ(80119) / 20736

rank_seven_genera = genera((7, 0), ZZ(948), even=True)
assert len(rank_seven_genera) == 4
assert sum(Genus(genus.representative()) == Genus(auxiliary) for genus in rank_seven_genera) == 1

auxiliary_automorphism_order, auxiliary_discriminant_image = (
    automorphism_discriminant_image(auxiliary)
)
assert auxiliary_automorphism_order == 7680
assert auxiliary_discriminant_image == [1, 157, 791, 947]

old_seed_rank_seven = []
old_seed_anti_discriminant_matches = []
for seed_path in sorted(SEEDS.rglob("*.txt")):
    seed_gram = load_seed_gram(seed_path)
    if (
        seed_gram is None
        or seed_gram.nrows() != 7
        or seed_gram.ncols() != 7
        or seed_gram.det() != 948
    ):
        continue
    old_seed_rank_seven.append(seed_path)
    unused_group, seed_q = discriminant_data(seed_gram)
    if anti_isometry_units(seed_q, published_q):
        old_seed_anti_discriminant_matches.append(seed_path)
assert len(old_seed_rank_seven) == 228
assert not old_seed_anti_discriminant_matches

controls = []
for label, path, frame in (
    ("published_R17", PUBLISHED, published),
    ("alternate_Q80", ALTERNATE, alternate),
):
    assert frame.nrows() == frame.ncols() == 17
    assert frame.det() == 948
    assert int(pari(frame).qfminim(2)[0]) == 0
    frame_group, frame_q = discriminant_data(frame)
    units = anti_isometry_units(auxiliary_q, frame_q)
    assert len(units) == 8
    frame_automorphism_order, frame_discriminant_image = (
        automorphism_discriminant_image(frame)
    )
    orbits = double_cosets(
        units,
        auxiliary_discriminant_image,
        frame_discriminant_image,
    )

    all_gluings = []
    representative_gluings = []
    for unit in units:
        certificate = glue(auxiliary, frame, unit)
        all_gluings.append(
            {
                "anti_isometry_unit": unit,
                "ambient_root_system": certificate["ambient_root_system"],
                "ambient_signed_root_count": certificate[
                    "ambient_signed_root_count"
                ],
                "ambient_gram_sha256": certificate["ambient_gram_sha256"],
                "complement_gram_sha256": certificate[
                    "complement_gram_sha256"
                ],
            }
        )
        if unit in {orbit[0] for orbit in orbits}:
            representative_gluings.append(certificate)

    controls.append(
        {
            "label": label,
            "source": str(path.relative_to(ROOT)),
            "source_sha256": file_sha256(path),
            "frame_discriminant_generator_q": str(frame_q),
            "frame_automorphism_group_order": frame_automorphism_order,
            "frame_automorphism_discriminant_image": frame_discriminant_image,
            "anti_isometry_units": units,
            "glue_double_cosets_under_AutK_times_AutR": orbits,
            "primitive_embedding_orbit_count": len(orbits),
            "all_gluings": all_gluings,
            "representative_gluings": representative_gluings,
        }
    )

assert controls[0]["primitive_embedding_orbit_count"] == 2
assert controls[1]["primitive_embedding_orbit_count"] == 1

cost_ambient = matrix(
    ZZ,
    controls[0]["representative_gluings"][0]["ambient_gram"],
)
cost_auxiliary = matrix(
    ZZ,
    controls[0]["representative_gluings"][0]["auxiliary_basis_in_ambient"],
)
embedded_d5 = cost_auxiliary[:5]
d5_orthogonal_coordinates = (
    embedded_d5 * cost_ambient
).right_kernel_matrix()
d5_orthogonal_gram = (
    d5_orthogonal_coordinates
    * cost_ambient
    * d5_orthogonal_coordinates.transpose()
)
assert d5_orthogonal_gram.nrows() == 19
assert d5_orthogonal_gram.det() == 4
norm_pair_counts = list(map(int, pari(d5_orthogonal_gram).qfrep(6, 1)))
assert norm_pair_counts == [76, 28899, 913368, 10487734, 70163352, 329206692]

result = {
    "schema": "elkies-k3.rootless-j2-niemeier-controls.v1",
    "status": "PASS_EXACT_AUXILIARY_AND_NIEMEIER_CONTROL_PROVENANCE_NOT_COMPLETE",
    "classification_scope": {
        "proved": (
            "The pinned rank-7 auxiliary has the required anti-discriminant "
            "form. For each mandatory rank-17 control, all eight cyclic "
            "discriminant anti-isometries give primitive gluings in "
            "N(2A7+2D5), with saturated rootless complement integrally "
            "isometric to that control."
        ),
        "not_proved": (
            "Primitive embeddings of the auxiliary into all 24 Niemeier "
            "lattices have not yet been enumerated, so the rootless J2 "
            "classification is not complete."
        ),
    },
    "auxiliary": {
        "source": str(AUXILIARY.relative_to(ROOT)),
        "source_sha256": file_sha256(AUXILIARY),
        "rank": 7,
        "determinant": 948,
        "gram": rows(auxiliary),
        "discriminant_group_invariants": [948],
        "discriminant_generator_q": str(auxiliary_q),
        "root_system": "D5",
        "signed_root_count": auxiliary_roots,
        "genus_count_among_even_rank7_det948": len(rank_seven_genera),
        "genus_mass": {
            "numerator": int(auxiliary_mass.numerator()),
            "denominator": int(auxiliary_mass.denominator()),
        },
        "automorphism_group_order": auxiliary_automorphism_order,
        "automorphism_discriminant_image": auxiliary_discriminant_image,
        "old_seed_rank7_det948_grams_checked": len(old_seed_rank_seven),
        "old_seed_anti_discriminant_matches": len(
            old_seed_anti_discriminant_matches
        ),
    },
    "niemeier_accounting": {
        "number_of_even_unimodular_positive_rank24_classes": 24,
        "rooted_classes": 23,
        "leech_excluded_because_auxiliary_has_norm_two_vectors": True,
        "control_ambient_root_system": "2A7+2D5",
        "root_system_identification": (
            "The 192 signed roots split into irreducible components with "
            "(rank,root-count) (7,56),(7,56),(5,40),(5,40)."
        ),
    },
    "naive_enumeration_cost_probe": {
        "source_embedding": "published_R17 first glue double-coset representative",
        "d5_orthogonal_rank": 19,
        "d5_orthogonal_determinant": 4,
        "sign_pair_counts_by_norm": {
            str(2 * (index + 1)): count
            for index, count in enumerate(norm_pair_counts)
        },
        "consequence": (
            "The 329206692 sign-pairs of norm 12 rule out materializing all "
            "choices for the sixth auxiliary basis vector. A complete "
            "Niemeier traversal must enumerate stabilizer or Weyl orbits."
        ),
    },
    "controls": controls,
}

payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
if arguments.check:
    assert arguments.output.read_text() == payload
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(payload)

print(
    "ROOTLESSJ2NIEMEIER|aux_q={}|aux_roots={}|published_orbits={}|"
    "alternate_orbits={}|ambient=2A7+2D5|complete=0|status=PASS".format(
        auxiliary_q,
        auxiliary_roots,
        controls[0]["primitive_embedding_orbit_count"],
        controls[1]["primitive_embedding_orbit_count"],
    )
)
