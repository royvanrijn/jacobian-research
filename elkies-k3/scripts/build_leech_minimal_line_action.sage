#!/usr/bin/env sage
"""Build the exact Co1 action on antipodal Leech minimal-vector lines.

The pinned Atlas matrices generate Co0 on the Leech lattice.  Independent
sign changes of generators do not change a generated sublattice, so the
natural prefix domain is the 98,280 antipodal pairs of minimal vectors.  The
central -1 acts trivially and the induced permutation group is Co1.

status: EXACT_LEECH_MINIMAL_LINE_ACTION_FOUNDATION
output: artifacts/generated-results/elkies-k3-leech-minimal-line-action-v1.json
"""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from sage.all import ZZ, libgap, matrix, pari


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BACKEND = ROOT / "artifacts/generated-results/elkies-k3-leech-co0-backend-v1.json"
DEFAULT_OUTPUT = (
    ROOT / "artifacts/generated-results/elkies-k3-leech-minimal-line-action-v1.json"
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_lines(values):
    nonzero = np.argmax(values != 0, axis=1)
    signs = np.where(
        values[np.arange(len(values)), nonzero] < 0,
        -1,
        1,
    )
    return (values * signs[:, None]).astype(np.int16)


def permutation_digest(images):
    encoded = "\n".join(map(str, images)) + "\n"
    return hashlib.sha256(encoded.encode()).hexdigest()


def build(backend):
    assert backend["schema"] == "elkies-k3.leech-co0-backend.v1"
    assert backend["status"] == (
        "PASS_EXACT_LEECH_GRAM_AND_CO0_ATLAS_ACTION_BACKEND_FOUNDATION"
    )
    gram = matrix(ZZ, backend["leech_lattice"]["gram"])
    generators = [
        matrix(ZZ, value)
        for value in backend["atlas_representation"]["generators"]
    ]
    raw = matrix(ZZ, pari(gram).qfminim(4)[2].sage()).transpose()
    assert raw.nrows() == 98280 and raw.ncols() == 24
    lines = canonical_lines(
        np.array([[int(entry) for entry in row] for row in raw.rows()], dtype=np.int16)
    )
    assert int(lines.min()) == -4 and int(lines.max()) == 4
    keys = [row.tobytes() for row in lines]
    assert len(set(keys)) == 98280
    index = {key: position for position, key in enumerate(keys)}

    gap_generators = []
    generator_rows = []
    for generator in generators:
        moved = canonical_lines(
            lines.astype(np.int64)
            @ np.array(generator, dtype=np.int64).transpose()
        )
        images_zero_based = [index[row.tobytes()] for row in moved]
        assert len(set(images_zero_based)) == len(lines)
        gap_generators.append(
            libgap.PermList([value + 1 for value in images_zero_based])
        )
        generator_rows.append(
            {
                "permutation_sha256": permutation_digest(images_zero_based),
                "fixed_lines": sum(
                    value == position
                    for position, value in enumerate(images_zero_based)
                ),
            }
        )

    group = libgap.Group(gap_generators)
    co0_order = backend["atlas_representation"]["group_order"]
    assert int(group.Size()) == co0_order // 2
    assert bool(libgap.IsTransitive(group))
    stabilizer = libgap.Stabilizer(group, 1)
    assert int(stabilizer.Size()) == 42305421312000
    domain = libgap(list(range(1, len(lines) + 1)))
    suborbits = libgap.Orbits(stabilizer, domain)

    fixed = lines[0].astype(np.int64)
    pairing = (
        lines.astype(np.int64)
        @ np.array(gram, dtype=np.int64)
        @ fixed
    )
    pairing_distribution = Counter(map(abs, map(int, pairing)))
    assert pairing_distribution == Counter({0: 46575, 1: 47104, 2: 4600, 4: 1})
    suborbit_rows = []
    seen_pairings = set()
    for orbit in suborbits:
        representative = int(orbit[0]) - 1
        absolute_pairing = abs(int(pairing[representative]))
        assert absolute_pairing not in seen_pairings
        seen_pairings.add(absolute_pairing)
        assert int(orbit.Length()) == pairing_distribution[absolute_pairing]
        suborbit_rows.append(
            {
                "absolute_inner_product": absolute_pairing,
                "size": int(orbit.Length()),
                "representative_line_index_zero_based": representative,
                "representative_coordinates": list(
                    map(int, lines[representative])
                ),
            }
        )
    suborbit_rows.sort(key=lambda row: row["absolute_inner_product"])
    assert [row["size"] for row in suborbit_rows] == [46575, 47104, 4600, 1]

    return {
        "schema": "elkies-k3.leech-minimal-line-action.v1",
        "status": "PASS_EXACT_CO1_ACTION_ON_98280_ANTIPODAL_MINIMAL_LINES",
        "proof_scope": {
            "proved": (
                "All 98,280 antipodal pairs in the norm-four shell are recovered "
                "in the pinned Leech basis. The two Atlas Co0 generators induce "
                "permutations whose group is transitive of order |Co0|/2=|Co1|; "
                "a line stabilizer has the exact Co2 order and four suborbits, "
                "identified intrinsically by absolute inner product."
            ),
            "not_proved": (
                "This action foundation does not enumerate rank-seven subsets or "
                "claim a determinant-band-complete Leech auxiliary census."
            ),
        },
        "input": {
            "backend_artifact": str(BACKEND.relative_to(ROOT)),
            "backend_sha256": digest(BACKEND),
        },
        "minimal_lines": {
            "oriented_minimal_vectors": 196560,
            "antipodal_lines": 98280,
            "coordinate_minimum": int(lines.min()),
            "coordinate_maximum": int(lines.max()),
            "canonical_orientation": "first nonzero ambient coordinate positive",
        },
        "induced_group": {
            "identified_group": "Co1=Co0/{+1,-1}",
            "order": int(group.Size()),
            "transitive": True,
            "generator_permutations": generator_rows,
            "line_stabilizer_identified_group": "Co2",
            "line_stabilizer_order": int(stabilizer.Size()),
            "line_stabilizer_suborbits": suborbit_rows,
        },
        "reproduction": {
            "command": (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elkies-k3/scripts/build_leech_minimal_line_action.sage"
            ),
            "check_command": (
                "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
                "elkies-k3/scripts/build_leech_minimal_line_action.sage --check"
            ),
        },
    }


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
parser.add_argument("--check", action="store_true")
arguments = parser.parse_args()
payload = build(json.loads(BACKEND.read_text()))
encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
if arguments.check:
    if not arguments.output.exists() or arguments.output.read_text() != encoded:
        raise SystemExit("Leech minimal-line action artifact is stale")
else:
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(encoded)
print(
    "LEECHLINES|lines=98280|group=Co1|suborbits=1,4600,46575,47104|"
    "status=PASS_EXACT"
)
