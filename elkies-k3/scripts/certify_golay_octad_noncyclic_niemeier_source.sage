#!/usr/bin/env sage
"""Glue the determinant-720 Golay auxiliary to its best source genus mate.

The determinant-720 discriminant module has invariants (2,6,60), so the
cyclic-unit shortcut used by the first foundry source certifier does not
apply.  This script enumerates exact anti-isometries on Smith generators,
builds the graph overlattice from three glue vectors, identifies the resulting
Niemeier lattice by its complete root system, and recovers the source as the
saturated orthogonal complement of the pinned Golay auxiliary.

status: ACTIVE_PROOF
claim: exact noncyclic discriminant gluing and named-Niemeier source embedding.
inputs: Golay octad target certificate and bounded source-hunt certificate.
outputs: artifacts/generated-results/elkies-k3-golay-octad-det720-source-niemeier.json
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

from sage.all import (
    IntegralLattice,
    QQ,
    ZZ,
    block_diagonal_matrix,
    identity_matrix,
    matrix,
    pari,
    vector,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_TARGET = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-rank17-det720.json"
)
DEFAULT_SOURCE = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-source-hunt.json"
)
CATALOG = ROOT / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-golay-octad-det720-source-niemeier.json"
)

_shared_path = HERE / "hunt_lattice_foundry_rootful_source.sage"
_shared_source = _shared_path.read_text().split(
    "parser = argparse.ArgumentParser", 1
)[0]
_shared = {"__file__": str(_shared_path)}
exec(compile(_shared_source, str(_shared_path), "exec"), _shared)

rows = _shared["rows"]
gram_digest = _shared["gram_digest"]
primitive_row_lattice = _shared["primitive_row_lattice"]
integral_matrix = _shared["integral_matrix"]
ambient_root_components = _shared["ambient_root_components"]
catalog_ambient_label = _shared["catalog_ambient_label"]


def module_order(module):
    return math.prod(map(int, module.invariants()))


def anti_isometries(left, right):
    left_generators = left.gens()
    right_elements = list(right)
    candidates = [
        [
            element
            for element in right_elements
            if element.additive_order() == generator.additive_order()
            and element.q() == -generator.q()
        ]
        for generator in left_generators
    ]
    assert len(left_generators) == 3
    result = []
    for image0 in candidates[0]:
        for image1 in candidates[1]:
            if image0.inner_product(image1) != -left_generators[0].inner_product(
                left_generators[1]
            ):
                continue
            for image2 in candidates[2]:
                if image0.inner_product(image2) != -left_generators[0].inner_product(
                    left_generators[2]
                ):
                    continue
                if image1.inner_product(image2) != -left_generators[1].inner_product(
                    left_generators[2]
                ):
                    continue
                images = (image0, image1, image2)
                if module_order(right.submodule_with_gens(images)) != module_order(right):
                    continue
                result.append(images)
    return result


def glue(auxiliary, source, images, catalog):
    auxiliary_module = IntegralLattice(auxiliary).discriminant_group()
    source_module = IntegralLattice(source).discriminant_group()
    assert auxiliary_module.invariants() == source_module.invariants() == (2, 6, 60)

    split_gram = block_diagonal_matrix(auxiliary, source)
    split_lattice = IntegralLattice(split_gram)
    glue_vectors = [
        vector(QQ, list(generator.lift()) + list(image.lift()))
        for generator, image in zip(auxiliary_module.gens(), images)
    ]
    ambient = split_lattice.overlattice(glue_vectors)
    assert ambient.rank() == 24
    old_gram = ambient.gram_matrix()
    assert old_gram.det() == 1

    change_columns = matrix(ZZ, pari(old_gram).qflllgram())
    change_rows = change_columns.transpose()
    ambient_gram = integral_matrix(
        change_rows * old_gram * change_rows.transpose()
    )
    ambient_basis = change_rows * ambient.basis_matrix()
    assert ambient_gram == ambient_basis * split_gram * ambient_basis.transpose()
    assert all(ambient_gram[index, index] % 2 == 0 for index in range(24))

    auxiliary_split_basis = identity_matrix(ZZ, 24)[:7]
    auxiliary_coordinates = integral_matrix(
        auxiliary_split_basis * ambient_basis.inverse()
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
    complement = (
        complement_coordinates
        * ambient_gram
        * complement_coordinates.transpose()
    )
    isometry = pari(source).qfisom(pari(complement))
    assert isometry != 0

    signed_roots, components = ambient_root_components(ambient_gram)
    label = catalog_ambient_label(components, catalog)
    return {
        "ambient_label": label,
        "ambient_signed_root_count": signed_roots,
        "ambient_root_components": components,
        "ambient_gram": rows(ambient_gram),
        "ambient_gram_sha256": gram_digest(ambient_gram),
        "auxiliary_basis_in_ambient": rows(auxiliary_coordinates),
        "complement_basis_in_ambient": rows(complement_coordinates),
        "complement_gram": rows(complement),
        "complement_gram_sha256": gram_digest(complement),
        "source_to_complement_isometry": rows(matrix(ZZ, isometry.sage()).transpose()),
        "primitive_auxiliary_embedding": True,
        "saturated_orthogonal_complement": True,
        "complement_integrally_isometric_to_source": True,
    }


def certify(arguments):
    target_payload = json.loads(arguments.target.read_text())
    source_payload = json.loads(arguments.source.read_text())
    catalog = json.loads(CATALOG.read_text())
    auxiliary = matrix(ZZ, target_payload["auxiliary"]["gram"])
    source = matrix(ZZ, source_payload["source"]["gram"])
    assert auxiliary.det() == source.det() == 720
    assert source_payload["source"]["genus_equals_target"]

    auxiliary_module = IntegralLattice(auxiliary).discriminant_group()
    source_module = IntegralLattice(source).discriminant_group()
    assert auxiliary_module.invariants() == source_module.invariants() == (2, 6, 60)
    maps = anti_isometries(auxiliary_module, source_module)
    assert maps
    selected = maps[0]
    niemeier = glue(auxiliary, source, selected, catalog)

    return {
        "schema": "elkies-k3.golay-octad-noncyclic-niemeier-source.v1",
        "status": "PASS_EXACT_DET720_NONCYCLIC_GLUE_AND_NIEMEIER_SOURCE",
        "proof_scope": {
            "proved": (
                "An explicit anti-isometry of the noncyclic discriminant "
                "modules glues the pinned determinant-720 Golay auxiliary and "
                "the bounded-search source frame to an even unimodular rank-24 "
                "lattice. Its complete root system identifies the named "
                "Niemeier ambient. The auxiliary embedding is primitive and "
                "its saturated orthogonal complement is integrally isometric "
                "to the source frame."
            ),
            "not_proved": (
                "The bounded Kneser beam is not a source classification and "
                "found no MW0--2 source in its declared window. No marked "
                "elliptic-neighbour corridor, rational source equation, "
                "arithmetic descent, or target section basis is constructed."
            ),
        },
        "inputs": {
            str(arguments.target.resolve().relative_to(ROOT)): hashlib.sha256(
                arguments.target.read_bytes()
            ).hexdigest(),
            str(arguments.source.resolve().relative_to(ROOT)): hashlib.sha256(
                arguments.source.read_bytes()
            ).hexdigest(),
            str(CATALOG.relative_to(ROOT)): hashlib.sha256(CATALOG.read_bytes()).hexdigest(),
        },
        "discriminant_glue": {
            "invariants": [2, 6, 60],
            "anti_isometry_count_in_smith_generator_enumeration": len(maps),
            "selected_images_in_source_smith_coordinates": [
                list(map(int, image)) for image in selected
            ],
            "selected_image_orders": [
                int(image.additive_order()) for image in selected
            ],
            "graph_glue_generators": 3,
            "graph_glue_order": 720,
        },
        "auxiliary": {
            "gram": rows(auxiliary),
            "gram_sha256": gram_digest(auxiliary),
            "rank": 7,
            "determinant": 720,
        },
        "source": {
            "gram": rows(source),
            "gram_sha256": gram_digest(source),
            "root_type": source_payload["source"]["root_type"],
            "root_rank": source_payload["source"]["root_rank"],
            "mw_rank_for_rho_19": source_payload["source"]["mw_rank_for_rho_19"],
            "mw_height_gram": source_payload["source"]["mw_height_gram"],
            "mw_regulator": source_payload["source"]["mw_regulator"],
        },
        "niemeier_certificate": niemeier,
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--target", type=Path, default=DEFAULT_TARGET)
parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

result = certify(arguments)
rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
if arguments.check:
    assert arguments.output.read_text() == rendered
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(rendered)

print(
    f"GOLAYGLUE|ambient={result['niemeier_certificate']['ambient_label']}"
    f"|source={result['source']['root_type']}"
    f"|root_rank={result['source']['root_rank']}"
    f"|mw_rank={result['source']['mw_rank_for_rho_19']}"
    f"|anti_isometries={result['discriminant_glue']['anti_isometry_count_in_smith_generator_enumeration']}"
    "|status=PASS_EXACT_NONCYCLIC_NIEMEIER_SOURCE"
)
