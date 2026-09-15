"""Replay only binary compatibility in the retained RR pair, not descent."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "artifacts/generated-results/elliptic-curves"


def rank(rows):
    basis = {}
    for row in rows:
        value = sum(bit << i for i, bit in enumerate(row))
        while value:
            pivot = value.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = value
                break
            value ^= basis[pivot]
    return len(basis)


def audit(path):
    data = json.loads(path.read_text())
    rows = data["matrix_rows"]
    assert rows and all(len(r) == 18 and set(r) <= {0, 1} for r in rows)
    blocks = [t for t in data["trials"] if "block_rows" in t]
    assert rows == [r for t in blocks for r in t["block_rows"]]
    assert all(rank(t["block_rows"]) == rank([r[:-1] for r in t["block_rows"]])
               for t in blocks)
    separator = data["separator"]
    assert len(separator) == len(rows) and set(separator) <= {0, 1}
    product = [sum(a * r[j] for a, r in zip(separator, rows)) % 2
               for j in range(18)]
    assert product == [0] * 17 + [1]
    assert rank([r[:-1] for r in rows]) == 16 and rank(rows) == 17
    offset = 0
    support = []
    for trial in blocks:
        count = len(trial["block_rows"])
        if any(separator[offset:offset + count]):
            support.append(trial["p"])
        offset += count
    return {
        "input": str(path.relative_to(ROOT)),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "rows": len(rows), "prime_blocks": len(blocks),
        "every_individual_block_consistent": True,
        "stacked_inherited_rank": 16, "stacked_augmented_rank": 17,
        "saved_separator_prime_support": support,
        "separator_product": product,
    }


if __name__ == "__main__":
    print(json.dumps({
        "status": "PASS",
        "boundary": "Binary replay of retained matrices only; no re-evaluation of divisors, global Selmer completeness, or new independence theorem.",
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "cases": {name: audit(ART / path) for name, path in [
            ("historical09", "det1092_rr_full_inherited_jacobian_v1.json"),
            ("control08", "det1092_rr_generic_point_controls_v2/case-08.json"),
        ]},
    }, indent=2))
