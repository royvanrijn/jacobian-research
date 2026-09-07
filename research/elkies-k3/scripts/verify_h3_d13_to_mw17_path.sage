#!/usr/bin/env sage -python
"""Replay the selected low-degree H3/D13-to-rootless-MW17 lattice path.

The H3 q=6 and q=8 entrance to the pinned D13/MW4 frame is certified by
``analyze_h3_first_q6_chamber.sage``.  Starting at that D13 frame, this
checker independently reconstructs the selected eleven primitive U-neighbors

    D13/MW4 --q24--> D12/MW5 --q6--> A11/MW6
              --q8--> 2A5/MW7 --q4--> 3A3/MW8
              --q4--> A3+2A2/MW10 --q4--> 5A1/MW12
              --q4--> 4A1/MW13 --q4--> 3A1/MW14
              --q4--> 2A1/MW15 --q4--> A1/MW16
              --q6--> rootless/MW17.

All eleven stored factor presentations have old-fiber degree two.  This is an
exact integral lattice certificate.  Except where a separate chamber checker
is cited, it does not assert that the raw degree-two class is already nef or
that the corresponding characteristic-zero pencil has been executed.
"""

from sage.all import *

import argparse
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
ENGINE = ROOT / "elkies-k3/scripts/exact_neighbor_engine.sage"
load(str(ENGINE))
U = matrix(ZZ, ((0, 1), (1, 0)))
SOURCE = ROOT / "elkies-k3/data/fibrations/h3_q6_q8_d13_mw4_root_adapted_frame.txt"
DEFAULT_OUTPUT = ROOT / "artifacts/generated-results/elkies-k3-h3-d13-to-mw17-path.json"

STEPS = (
    {
        "artifact": "artifacts/generated-results/elkies-k3-h3-q6-q8-d13-q24-degree2.json",
        "sha256": "66d5a7ff6ec26f8aa8344cdbd779a6c96707b041ba4f89d7dbfe460c95485a93",
        "orbit": 85,
        "q": 24,
        "root_data": (12, 264, 4),
        "ade": "D12",
        "mw": 5,
    },
    {
        "artifact": "artifacts/generated-results/elkies-k3-h3-d12-o85-q6-degree2.json",
        "sha256": "1b8d7f37794bcf49c48949cec3bffe7baba69cb55e369f7af9f5002908a75b7f",
        "orbit": 42,
        "q": 6,
        "root_data": (11, 132, 12),
        "ade": "A11",
        "mw": 6,
    },
    {
        "artifact": "artifacts/generated-results/elkies-k3-h3-a11-middle-q8-degree2.json",
        "sha256": "d336b0d32a07907ec61464c7d5ace4c76f257ddcf16b904a96a3d9064f408323",
        "orbit": 922,
        "q": 8,
        "root_data": (10, 60, 36),
        "ade": "A5+A5",
        "mw": 7,
    },
    {
        "artifact": "artifacts/generated-results/elkies-k3-h3-a5a5-c2-q4-degree2.json",
        "sha256": "98fdd553768b27d5800f247b41a6e2a28f0ee2787ad5959d96ef420a2eb09185",
        "orbit": 472,
        "q": 4,
        "root_data": (9, 36, 64),
        "ade": "A3+A3+A3",
        "mw": 8,
    },
    {
        "artifact": "artifacts/generated-results/elkies-k3-h3-a3x3-q4-degree2.json",
        "sha256": "8f1d5105831cc3356bc4598932380295bcc6e91629e0d05a6a9d64c0d840d29d",
        "orbit": 323,
        "q": 4,
        "root_data": (7, 24, 36),
        "ade": "A2+A2+A3",
        "mw": 10,
    },
    {
        "artifact": "artifacts/generated-results/elkies-k3-h3-mw10-a3a2a2-q4-degree2.json",
        "sha256": "e204b81e6fae699f89e5abb97e632318809249b057309d2e31db5911348c2254",
        "orbit": 207,
        "q": 4,
        "root_data": (5, 10, 32),
        "ade": "A1+A1+A1+A1+A1",
        "mw": 12,
    },
    {
        "artifact": "artifacts/generated-results/elkies-k3-h3-mw12-5a1-q4-degree2-first-hit.json",
        "sha256": "4b2b182b55df980c730860d29d3cc7c4e3c6735f36c533283d9747456742ce73",
        "orbit": 52,
        "q": 4,
        "root_data": (4, 8, 16),
        "ade": "A1+A1+A1+A1",
        "mw": 13,
    },
    {
        "artifact": "artifacts/generated-results/elkies-k3-h3-mw13-4a1-q4-degree2-first-hit.json",
        "sha256": "271bf7c1ec2ed8634d503b9aba6fd0564ba46dab00bb3b5a9b361225b7ceeb89",
        "orbit": 114,
        "q": 4,
        "root_data": (3, 6, 8),
        "ade": "A1+A1+A1",
        "mw": 14,
    },
    {
        "artifact": "artifacts/generated-results/elkies-k3-h3-mw14-3a1-q4-degree2-first-hit.json",
        "sha256": "ed9e110aa7cd430097561b1bc033a381656e4e983064c923a859fd5a1fdb9cc3",
        "orbit": 498,
        "q": 4,
        "root_data": (2, 4, 4),
        "ade": "A1+A1",
        "mw": 15,
    },
    {
        "artifact": "artifacts/generated-results/elkies-k3-h3-mw15-2a1-q4-degree2-first-hit.json",
        "sha256": "7ef18e1c855eb62c6460bbd5bd0c861113e2a2c289060cbd85b7360965d5633f",
        "orbit": 981,
        "q": 4,
        "root_data": (1, 2, 2),
        "ade": "A1",
        "mw": 16,
    },
    {
        "artifact": "artifacts/generated-results/elkies-k3-h3-mw16-a1-q6-degree2-cap10000-stream-chunk001.json",
        "sha256": "1d4ca7cfc8ebee80d7364d87ac9ca19e18742130925fa1b54cf57abbe774de64",
        "orbit": 2247,
        "q": 6,
        "root_data": (0, 0, 1),
        "ade": "rootless",
        "mw": 17,
    },
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display_path(path):
    """Use a repository-relative path when possible, including CLI relatives."""
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(path)


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ],
    )


def matrix_rows(value):
    return [[int(entry) for entry in row] for row in value.rows()]


def rational_rows(value):
    return [[str(entry) for entry in row] for row in value.rows()]


def matrix_digest(value):
    payload = ";".join(
        ",".join(str(entry) for entry in row) for row in value.rows()
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def child_frame(ns, fiber):
    old_fiber = vector(ZZ, [1, 0] + [0] * (ns.nrows() - 2))
    result = degree_two_neighbor(ns, fiber, old_fiber, curves=())
    return result["child_frame"], result["transport"]


def root_data(frame):
    return tuple(int(value) for value in roots_and_data(frame)[2])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()

    source = load_gram(SOURCE)
    assert source.nrows() == 17 and source.det() == 948
    assert root_data(source) == (13, 312, 4)
    current = source
    composite = identity_matrix(ZZ, 19)
    replay = []

    for index, specification in enumerate(STEPS, 1):
        artifact_path = ROOT / specification["artifact"]
        assert digest(artifact_path) == specification["sha256"]
        data = json.loads(artifact_path.read_text())
        assert data["status"] == "PASS_ROOT_ADAPTED_WEYL_NEIGHBORS"
        declared_source = ROOT / data["frame"]
        assert load_gram(declared_source) == current

        matches = [
            record
            for record in data["neighbors"]
            if record["orbit_index"] == specification["orbit"]
        ]
        assert len(matches) == 1
        record = matches[0]
        assert record["q"] == specification["q"]
        assert record["factor_order"] == [specification["q"] // 2, 2]
        assert record["old_fiber_degree"] == 2
        assert tuple(record["child_root_data"]) == specification["root_data"]
        assert record["child_ade"] == specification["ade"]
        assert record["child_mw_rank"] == specification["mw"]

        witness = vector(ZZ, record["witness"])
        assert witness * current * witness == 2 * specification["q"]
        fiber = vector(ZZ, record["fiber"])
        assert fiber == vector(
            ZZ,
            [specification["q"] // 2, 2] + list(witness),
        )
        assert gcd(tuple(fiber)) == 1

        ns = block_diagonal_matrix(U, -current)
        child, neighbor_basis = child_frame(ns, fiber)
        assert child == matrix(ZZ, record["child_frame"])
        adapted_basis = matrix(ZZ, record["child_root_adapted_basis"])
        assert abs(adapted_basis.det()) == 1
        adapted = adapted_basis * child * adapted_basis.transpose()
        assert adapted == matrix(ZZ, record["child_root_adapted_frame"])
        assert root_data(adapted) == specification["root_data"]

        root_rank = specification["root_data"][0]
        cartan = adapted[:root_rank, :root_rank]
        coupling = adapted[:root_rank, root_rank:]
        tail = adapted[root_rank:, root_rank:]
        height = tail - coupling.transpose() * cartan.inverse() * coupling
        assert height == matrix(QQ, record["child_mw_height"])

        full_adaptation = block_diagonal_matrix(
            identity_matrix(ZZ, 2), adapted_basis
        )
        step_basis = full_adaptation * neighbor_basis
        assert abs(step_basis.det()) == 1
        composite = step_basis * composite
        current = adapted
        assert (
            composite
            * block_diagonal_matrix(U, -source)
            * composite.transpose()
            == block_diagonal_matrix(U, -current)
        )

        replay.append(
            {
                "step": index,
                "q": specification["q"],
                "factor_order": [specification["q"] // 2, 2],
                "orbit_index": specification["orbit"],
                "witness": list(map(int, witness)),
                "root_data": list(specification["root_data"]),
                "ade": specification["ade"],
                "mw_rank": specification["mw"],
                "mw_height": rational_rows(height),
                "artifact": specification["artifact"],
                "artifact_sha256": specification["sha256"],
            }
        )
        print(
            "H3D13MW17|step={}|q={}|ab={},2|orbit={}|ADE={}|MW={}|status=PASS".format(
                index,
                specification["q"],
                specification["q"] // 2,
                specification["orbit"],
                specification["ade"],
                specification["mw"],
            ),
            flush=True,
        )

    assert root_data(current) == (0, 0, 1)
    assert abs(composite.det()) == 1
    payload = {
        "schema": "elkies-k3.h3-d13-to-mw17-path.v1",
        "status": "PASS_H3_D13_TO_MW17_LATTICE_PATH",
        "proof_boundary": (
            "This is an exact integral U-neighbor and composite-transport "
            "certificate.  Companion exact chamber checkers certify nefness "
            "and old-fiber degree two, but neither layer executes the "
            "corresponding characteristic-zero pencils."
        ),
        "companion_chamber_certificates": [
            {
                "checker": "elkies-k3/scripts/analyze_h3_d13_q4_chamber.sage",
                "scope": "D13/MW4 to D12/MW5",
            },
            {
                "checker": "elkies-k3/scripts/analyze_h3_rank_growing_degree2_chain.sage",
                "scope": "D12/MW5 to A3+2A2/MW10",
                "status": "PASS_H3_RANK_GROWING_DEGREE2_CHAIN",
            },
            {
                "checker": "elkies-k3/scripts/analyze_h3_mw10_to_rootless_chambers.sage",
                "scope": "A3+2A2/MW10 to rootless/MW17",
                "status": "PASS_H3_MW10_TO_ROOTLESS_NEF",
            },
        ],
        "source": str(SOURCE.relative_to(ROOT)),
        "source_root_data": [13, 312, 4],
        "source_ade": "D13",
        "source_mw_rank": 4,
        "steps": replay,
        "final_root_data": [0, 0, 1],
        "final_ade": "rootless",
        "final_mw_rank": 17,
        "final_frame": matrix_rows(current),
        "composite_transport": matrix_rows(composite),
        "composite_transport_sha256": matrix_digest(composite),
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=int) + "\n"
    )
    print(
        "H3D13MW17|steps=11|final=rootless|MW=17|transport_det={}|"
        "transport_sha256={}|artifact={}|status={}".format(
            composite.det(),
            payload["composite_transport_sha256"],
            display_path(arguments.output),
            payload["status"],
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
