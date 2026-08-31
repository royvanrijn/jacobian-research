#!/usr/bin/env sage
"""Build and certify the 23 rooted Niemeier Gram models.

status: ACTIVE_PROOF
claim: Hash-pinned catalogue Gram matrices give all 23 rooted Niemeier
  lattices.  Their exact root decompositions exclude the ten rooted classes
  having no D5 root subsystem; the Leech class is excluded separately because
  the pinned auxiliary has roots.
inputs: hash-pinned pages in the Nebe--Sloane Catalogue of Lattices
outputs: artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json
supersedes/superseded-by: none
"""

import argparse
import hashlib
import html
import json
import re
import urllib.request
from pathlib import Path

from sage.all import ZZ, matrix, pari


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUTPUT = (
    ROOT
    / "artifacts/generated-results/elkies-k3-rooted-niemeier-catalog.json"
)
BASE_URL = (
    "https://www.math.rwth-aachen.de/~Gabriele.Nebe/"
    "LATTICES/NIEM_{}.html"
)


CATALOG = [
    ("D24", [("D", 24)], "b6d584efc6524d69ca8911f69df617cd9c717c1b58cecefe96a8fc9e8c8a4696"),
    ("D16_E8", [("D", 16), ("E", 8)], "2b784245ff5503301aa5acfd5bdfa5673869d61b9adc2aeeb24c9aafd4f1b6ff"),
    ("3E8", [("E", 8)] * 3, "addc568b9c57c90717ab084ef8753e3f0d848a41b72f0636b8287925ce3e4d2e"),
    ("A24", [("A", 24)], "9b2a83855182ff134ffe138c1203e36b952cd4d0e6006108380886b6b5577b78"),
    ("2D12", [("D", 12)] * 2, "90eaf3857efd531badbf88190312cacce2df4bd488a4e4fe901995f4c02372f2"),
    ("A17_E7", [("A", 17), ("E", 7)], "7c0a4f32049e26f8a4bf645e8123152c9e5b8e66cdc6041b76ef1e654b04188b"),
    ("D10_2E7", [("D", 10), ("E", 7), ("E", 7)], "d4876c9e759431ee870d21817d24469cd41f27307be80ef765b825d531e60593"),
    ("A15_D9", [("A", 15), ("D", 9)], "598c4c4358ee63d5f85899f74e741042839669290b3cfb63169bb7746f4ba982"),
    ("3D8", [("D", 8)] * 3, "85894834a3a2fed3749c2467f683e7613682f49542643383d2dcab58aaabe369"),
    ("2A12", [("A", 12)] * 2, "6e7a2670a14d6731e3a7e5d8f0557d732ff5494c3b9d3b1f7cab670a01c1bd54"),
    ("A11_D7_E6", [("A", 11), ("D", 7), ("E", 6)], "ef982659d7703e83d32a0e920509a52081f353e450c32adea9d156692182c44d"),
    ("4E6", [("E", 6)] * 4, "4f30cf8a20b8dfeb9cc4e2cf68954d8873703b700fe455a88c3b31fcebd04b15"),
    ("2A9_D6", [("A", 9), ("A", 9), ("D", 6)], "7456c90a410bcd644e8dcd3c546322dca12f921b1f5beb2530cd3f0989cbd2bd"),
    ("4D6", [("D", 6)] * 4, "af4d5ce08e87cb34778dab58be548ad8dbd549a5fc13a59fee93fa66fa80a24e"),
    ("3A8", [("A", 8)] * 3, "dbc586c9fb81edfa444cbc9315ce65923d3ba93ecf2ec748b7efd744ad2e83fd"),
    ("2A7_2D5", [("A", 7), ("A", 7), ("D", 5), ("D", 5)], "691dc36df2c8704116e9e14fdfdffe2fc7386ed471218b90ee966b68470411f1"),
    ("4A6", [("A", 6)] * 4, "bba89dd64c5610268f22465027eb3a4cde1b2081d44d58d028b3a0f360301039"),
    ("4A5_D4", [("A", 5)] * 4 + [("D", 4)], "052548bb6917df4a8b922b72c7f642efc19d66da14c99cf8b8179c4c68144637"),
    ("6D4", [("D", 4)] * 6, "d53ba07b29e82b0c1b31ca68af0d062b91964a06de3a3b00a1a0e90c7464d226"),
    ("6A4", [("A", 4)] * 6, "844d9da6ac22a33f8dbc4309c8c500dbae5b0e8bd99f2a73058699139205e69f"),
    ("8A3", [("A", 3)] * 8, "1e4f5875c80c986cd0976b4ec2aa2fbd45086cfaa1aeeb0d4a6b4096c7af1d9a"),
    ("12A2", [("A", 2)] * 12, "cd356ada4f840c68ba5092e80beba00a9d6ab9541f9ca8c7733093d3a614980b"),
    ("24A1", [("A", 1)] * 24, "2299ffce39af7704b28df4e88c832f2deeb916deb8af10e6fea8b174ec361efc"),
]


def rows(value):
    return [list(map(int, row)) for row in value.rows()]


def component_invariant(component):
    family, rank = component
    if family == "A":
        roots = rank * (rank + 1)
    elif family == "D":
        roots = 2 * rank * (rank - 1)
    else:
        roots = {6: 72, 7: 126, 8: 240}[rank]
    return (rank, roots)


def root_components(gram):
    minimum_data = pari(gram).qfminim(2)
    assert int(minimum_data[1]) == 2
    representatives = matrix(ZZ, minimum_data[2].sage()).transpose()
    assert 2 * representatives.nrows() == int(minimum_data[0])
    unseen = set(range(representatives.nrows()))
    result = []
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
        vectors = representatives[sorted(component)]
        result.append((int(vectors.rank()), 2 * len(component)))
    return sorted(result, reverse=True)


def parse_catalog_page(raw):
    source = raw.decode("latin1")
    match = re.search(
        r'<a NAME="GRAM"><STRONG>GRAM</STRONG></a><br>(.*?)<p><li>',
        source,
        re.IGNORECASE | re.DOTALL,
    )
    assert match is not None
    block = re.sub(r"<br\s*/?>", "\n", match.group(1), flags=re.IGNORECASE)
    block = re.sub(r"<[^>]+>", "", block)
    lines = [
        html.unescape(line).strip()
        for line in block.splitlines()
        if line.strip()
    ]
    assert lines[0] == "24 0"
    triangular = [[ZZ(value) for value in line.split()] for line in lines[1:]]
    assert len(triangular) == 24
    assert all(len(row) == index + 1 for index, row in enumerate(triangular))
    gram = matrix(ZZ, 24)
    for row, values in enumerate(triangular):
        for column, value in enumerate(values):
            gram[row, column] = value
            gram[column, row] = value
    return gram


def validate_entry(label, expected_components, expected_hash, entry):
    assert entry["label"] == label
    assert entry["url"] == BASE_URL.format(label)
    assert entry["catalog_page_sha256"] == expected_hash
    gram = matrix(ZZ, entry["gram"])
    assert gram.nrows() == gram.ncols() == 24
    assert gram.det() == 1
    assert gram.is_positive_definite()
    assert all(gram[index, index] % 2 == 0 for index in range(24))
    observed = root_components(gram)
    expected = sorted(map(component_invariant, expected_components), reverse=True)
    assert observed == expected
    assert entry["root_components"] == [
        {
            "family": family,
            "rank": rank,
            "signed_root_count": component_invariant((family, rank))[1],
        }
        for family, rank in expected_components
    ]
    admits_d5 = any(
        (family == "D" and rank >= 5) or family == "E"
        for family, rank in expected_components
    )
    assert entry["admits_D5_root_subsystem"] is admits_d5
    return gram, admits_d5


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=OUTPUT)
parser.add_argument(
    "--refresh",
    action="store_true",
    help="Download the hash-pinned catalogue pages and rebuild the artifact.",
)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
assert not (arguments.refresh and arguments.check)

if arguments.refresh:
    entries = []
    for label, components, expected_hash in CATALOG:
        url = BASE_URL.format(label)
        with urllib.request.urlopen(url) as response:
            raw = response.read()
        assert hashlib.sha256(raw).hexdigest() == expected_hash
        gram = parse_catalog_page(raw)
        entries.append(
            {
                "label": label,
                "url": url,
                "catalog_page_sha256": expected_hash,
                "gram": rows(gram),
                "root_components": [
                    {
                        "family": family,
                        "rank": rank,
                        "signed_root_count": component_invariant((family, rank))[1],
                    }
                    for family, rank in components
                ],
                "admits_D5_root_subsystem": any(
                    (family == "D" and rank >= 5) or family == "E"
                    for family, rank in components
                ),
            }
        )
else:
    entries = json.loads(arguments.output.read_text())["rooted_niemeier_lattices"]

assert len(entries) == len(CATALOG) == 23
admissible = []
excluded = []
for specification, entry in zip(CATALOG, entries):
    unused_gram, admits_d5 = validate_entry(*specification, entry)
    (admissible if admits_d5 else excluded).append(entry["label"])

assert len(admissible) == 13
assert len(excluded) == 10
assert "2A7_2D5" in admissible

result = {
    "schema": "elkies-k3.rooted-niemeier-catalog.v1",
    "status": "PASS_EXACT_23_ROOTED_NIEMEIER_MODELS_AND_D5_FILTER_NOT_EMBEDDING_COMPLETE",
    "source": {
        "name": "Nebe--Sloane Catalogue of Lattices",
        "base_url": BASE_URL,
        "page_hash_algorithm": "sha256",
    },
    "accounting": {
        "rooted_niemeier_classes": 23,
        "leech_classes": 1,
        "all_niemeier_classes": 24,
        "D5_admissible_rooted_classes": len(admissible),
        "D5_excluded_rooted_classes": len(excluded),
        "leech_excluded_because_auxiliary_contains_D5_roots": True,
        "classes_requiring_embedding_enumeration": admissible,
        "classes_excluded_before_embedding_enumeration": excluded + ["Leech"],
    },
    "classification_scope": {
        "proved": (
            "All 23 rooted Niemeier Gram models are even, positive-definite, "
            "unimodular of rank 24 and have the stated exact ADE root "
            "decomposition. Ten rooted classes and the Leech class cannot "
            "contain the pinned auxiliary because its connected D5 root "
            "subsystem has nowhere to embed."
        ),
        "not_proved": (
            "The primitive embeddings in the thirteen surviving rooted "
            "classes have not yet been enumerated modulo symmetry."
        ),
    },
    "rooted_niemeier_lattices": entries,
}

payload = json.dumps(result, indent=2, sort_keys=True) + "\n"
if arguments.check:
    assert arguments.output.read_text() == payload
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(payload)

print(
    "NIEMEIERCATALOG|rooted=23|D5_admissible={}|excluded={}|"
    "leech_excluded=1|remaining={}|status=PASS".format(
        len(admissible),
        len(excluded),
        len(admissible),
    )
)
