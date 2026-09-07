#!/usr/bin/env sage
"""Audit the stored rootless rank-17 corpus and the full genus mass.

status: ACTIVE_SEARCH
claim: Exact lower bound for rootless J2 frame classes in the stored
  good-prime neighbour corpus, together with the exact genus mass obstruction
  to an unfiltered exhaustive traversal.  This does not certify completeness.
inputs: elkies-k3/data/lattice/rank17_gram.txt,
  artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json,
  elkies-k3/seeds/target-genus-rootless-pneighbor/*.txt
outputs: artifacts/generated-results/elkies-k3-rootless-j2-completeness-track.json
supersedes/superseded-by: none
"""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from sage.all import ZZ, matrix, pari
from sage.quadratic_forms.genera.genus import Genus


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PINNED = ROOT / "elkies-k3/data/lattice/rank17_gram.txt"
ALTERNATE = (
    ROOT
    / "artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json"
)
CORPUS = ROOT / "elkies-k3/seeds/target-genus-rootless-pneighbor"
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rootless-j2-completeness-track.json"
)


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_matrix(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def gram_sha256(value):
    payload = "\n".join(" ".join(map(str, row)) for row in value.rows()) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def lll_gram(value):
    transform = matrix(ZZ, pari(value).qflllgram())
    reduced = transform.transpose() * value * transform
    assert abs(transform.det()) == 1
    assert reduced.det() == value.det()
    return reduced


def source_metadata(path):
    result = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith("# ") or " = " not in line:
            continue
        key, value = line[2:].split(" = ", 1)
        try:
            result[key] = int(value)
        except ValueError:
            result[key] = value
    return result


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()

pinned = load_matrix(PINNED)
alternate_payload = json.loads(ALTERNATE.read_text())
alternate = matrix(ZZ, alternate_payload["rootless_frame"])

sources = [
    {
        "label": "published_R17",
        "path": PINNED,
        "gram": pinned,
        "control": True,
        "metadata": {"provenance": "published R17 frame"},
    },
    {
        "label": "alternate_Q80",
        "path": ALTERNATE,
        "gram": alternate,
        "control": True,
        "metadata": {"provenance": "alternate Q80-derived frame"},
    },
]

for path in sorted(CORPUS.glob("*.txt")):
    sources.append(
        {
            "label": path.name,
            "path": path,
            "gram": load_matrix(path),
            "control": False,
            "metadata": source_metadata(path),
        }
    )

target_genus = Genus(pinned)
mass = target_genus.mass()
local_symbols = {
    str(prime): str(target_genus.local_symbol(prime)).split(":", 1)[1].strip()
    for prime in (2, 3, 79)
}

representatives = []
excluded_representatives = []
source_rows = []

for source in sources:
    gram = lll_gram(source["gram"])
    assert gram.nrows() == gram.ncols() == 17
    assert gram.is_positive_definite()
    assert gram.det() == 948
    assert all(gram[index, index] % 2 == 0 for index in range(17))
    source_genus = Genus(gram)
    genus_match = source_genus == target_genus

    minimum_data = pari(gram).qfminim()
    minimum = int(minimum_data[1])
    minimal_pairs = int(ZZ(minimum_data[0]) // 2)
    norm_two_vectors = int(pari(gram).qfminim(2)[0])
    assert minimum == 4
    assert norm_two_vectors == 0

    class_bucket = representatives if genus_match else excluded_representatives
    class_index = None
    for index, representative in enumerate(class_bucket, start=1):
        if pari(representative["gram"]).qfisom(pari(gram)) != 0:
            class_index = index
            representative["members"].append(source["label"])
            break

    if class_index is None:
        class_index = len(class_bucket) + 1
        automorphism_order = int(pari(gram).qfauto()[0])
        class_bucket.append(
            {
                "class_index": class_index,
                "representative": source["label"],
                "members": [source["label"]],
                "gram": gram,
                "reduced_gram_sha256": gram_sha256(gram),
                "norm_four_pairs": minimal_pairs,
                "automorphism_group_order": automorphism_order,
            }
        )

    source_rows.append(
        {
            "label": source["label"],
            "class_index": class_index,
            "class_namespace": "target_genus" if genus_match else "excluded_genus",
            "target_genus_match": genus_match,
            "control": source["control"],
            "minimum": minimum,
            "norm_two_vectors": norm_two_vectors,
            "norm_four_pairs": minimal_pairs,
            "source": str(source["path"].relative_to(ROOT)),
            "source_sha256": file_sha256(source["path"]),
            "metadata": source["metadata"],
        }
    )

assert source_rows[0]["class_index"] == 1
assert source_rows[1]["class_index"] == 2
assert source_rows[0]["target_genus_match"]
assert source_rows[1]["target_genus_match"]
assert len(representatives) == 2
assert len(excluded_representatives) == 19
assert all(not row["target_genus_match"] for row in source_rows[2:])

excluded_genus = Genus(lll_gram(sources[2]["gram"]))
assert all(Genus(lll_gram(source["gram"])) == excluded_genus for source in sources[2:])
excluded_local_symbols = {
    str(prime): str(excluded_genus.local_symbol(prime)).split(":", 1)[1].strip()
    for prime in (2, 3, 79)
}

observed_rootless_mass = sum(
    (ZZ(1) / representative["automorphism_group_order"])
    for representative in representatives
)

for representative in representatives:
    representative["reduced_gram"] = rows(representative.pop("gram"))
    representative["mass_contribution"] = {
        "numerator": 1,
        "denominator": representative["automorphism_group_order"],
    }

for representative in excluded_representatives:
    representative.pop("gram")
    representative["mass_contribution"] = {
        "numerator": 1,
        "denominator": representative["automorphism_group_order"],
    }

twice_mass = 2 * mass
genus_class_lower_bound = int(twice_mass.ceil())
prime_histogram = Counter(
    row["metadata"].get("prime")
    for row in source_rows
    if not row["control"] and "prime" in row["metadata"]
)

result = {
    "schema": "elkies-k3.rootless-j2-completeness-track.v1",
    "status": "PASS_EXACT_ROOTLESS_J2_CONTROLS_AND_OFF_GENUS_REJECTION_NOT_COMPLETE",
    "classification_scope": {
        "proved": (
            "The published R17 and alternate Q80 controls lie in the same even "
            "positive-definite rank-17 genus of determinant 948 and are not "
            "integrally isometric. All 65 old files named target-genus-rootless-"
            "pneighbor instead lie in one different 2-adic and 79-adic genus; "
            "they are excluded from the J2 corpus."
        ),
        "not_proved": (
            "The two valid controls are not an exhaustive genus or Niemeier-"
            "embedding enumeration, so no J2 completeness claim is made."
        ),
    },
    "genus": {
        "rank": 17,
        "signature": [17, 0],
        "even": True,
        "determinant": 948,
        "local_symbols": local_symbols,
        "mass": {
            "numerator": int(mass.numerator()),
            "denominator": int(mass.denominator()),
            "decimal": float(mass),
        },
        "class_number_lower_bound_from_minus_identity": genus_class_lower_bound,
        "class_number_lower_bound_reason": (
            "Every lattice automorphism group contains plus/minus identity, "
            "so every mass summand is at most 1/2 and h is at least ceil(2*mass)."
        ),
    },
    "corpus": {
        "source_count": len(sources),
        "mandatory_control_count": 2,
        "stored_neighbour_output_count": len(sources) - 2,
        "target_genus_source_count": sum(
            row["target_genus_match"] for row in source_rows
        ),
        "target_genus_rootless_isometry_class_count": len(representatives),
        "excluded_off_genus_source_count": sum(
            not row["target_genus_match"] for row in source_rows
        ),
        "excluded_off_genus_isometry_class_count": len(excluded_representatives),
        "observed_rootless_mass": {
            "numerator": int(observed_rootless_mass.numerator()),
            "denominator": int(observed_rootless_mass.denominator()),
            "decimal": float(observed_rootless_mass),
        },
        "neighbour_prime_histogram": {
            str(key): value for key, value in sorted(prime_histogram.items())
        },
        "sources": source_rows,
        "target_genus_classes": representatives,
        "excluded_genus": {
            "local_symbols": excluded_local_symbols,
            "reason": (
                "The 2-adic and 79-adic symbols differ from the R17/Q80 "
                "target; determinant, rank, parity, and rootlessness alone "
                "do not determine the genus."
            ),
            "classes": excluded_representatives,
        },
    },
    "control_fingerprints": {
        "published_R17": {
            "class_index": source_rows[0]["class_index"],
            "norm_four_pairs": source_rows[0]["norm_four_pairs"],
            "automorphism_group_order": representatives[0][
                "automorphism_group_order"
            ],
        },
        "alternate_Q80": {
            "class_index": source_rows[1]["class_index"],
            "norm_four_pairs": source_rows[1]["norm_four_pairs"],
            "automorphism_group_order": representatives[1][
                "automorphism_group_order"
            ],
        },
        "integrally_isometric": False,
        "observed_ordering": (
            "The controls have 1311 and 1313 norm-four pairs. The old files "
            "beginning at 1320 cannot be used to rank target-genus frames "
            "because the local-symbol audit excludes them."
        ),
    },
    "route_decision": {
        "unfiltered_full_genus_traversal": (
            "Rejected as the cheap first formulation: the exact mass forces "
            "at least 167967 total isometry classes before rootless filtering."
        ),
        "next_completeness_gate": (
            "Enumerate primitive embeddings of the pinned rank-seven "
            "auxiliary lattice into all 24 Niemeier lattices modulo Weyl and "
            "glue automorphisms; compute primitive-closure complements, filter "
            "rootless determinant-948 frames, verify the full local genus, and "
            "deduplicate against the two mandatory controls."
        ),
        "construction_independence": (
            "This classification track does not alter or block the generic "
            "characteristic-zero lift of the alternate Q80 frame."
        ),
    },
    "inputs": {
        str(PINNED.relative_to(ROOT)): file_sha256(PINNED),
        str(ALTERNATE.relative_to(ROOT)): file_sha256(ALTERNATE),
        "rootless_neighbour_corpus": str(CORPUS.relative_to(ROOT)),
    },
    "reproduce": (
        "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
        "elkies-k3/scripts/audit_rootless_j2_completeness_track.sage"
    ),
}

serialized = json.dumps(result, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not arguments.output.exists() or arguments.output.read_text() != serialized:
        raise SystemExit("rootless-J2 completeness-track artifact is stale")
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(serialized)

print(
    "ROOTLESSJ2|genus_mass={}/{}|genus_class_lower_bound={}|"
    "corpus_sources={}|target_rootless_classes={}|off_genus_sources={}|"
    "off_genus_classes={}|observed_rootless_mass={}/{}|"
    "complete=0|status=PASS".format(
        mass.numerator(),
        mass.denominator(),
        genus_class_lower_bound,
        len(sources),
        len(representatives),
        len(sources) - 2,
        len(excluded_representatives),
        observed_rootless_mass.numerator(),
        observed_rootless_mass.denominator(),
    ),
    flush=True,
)
