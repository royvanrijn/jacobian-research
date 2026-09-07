#!/usr/bin/env sage-python
"""Compare compact and expanded NS0007 inputs as Sage polynomials.

The compact file prints the last I4 and exact-I2 equations in factored form.
This auditor parses both inputs over their declared finite field and compares
all polynomial equations exactly after expansion.  This does not certify that
msolve's input parser or preprocessing treats the two syntaxes identically.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from sage.all import GF, PolynomialRing


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_EXPANDED = (
    ROOT / "artifacts/local/elkies-k3/ns0007-pole0-reduced-modp/p7-lambda2.json"
)
DEFAULT_COMPACT = (
    ROOT
    / "artifacts/local/elkies-k3/ns0007-pole0-reduced-modp/"
    "p7-lambda2-compact.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / "artifacts/generated-results/"
    "elkies-k3-lattice-foundry-ns0007-compact-msolve-encoding-audit-mod7.json"
)
SYSTEM_SCHEMA = "elkies-k3.lattice-foundry-ns0007-pole0-reduced-modp-system.v1"


def relative(path: Path) -> str:
    return str(path.resolve().relative_to(ROOT))


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(metadata_path: Path) -> tuple[dict, Path, list[str], int, list[str]]:
    metadata_path = metadata_path.resolve()
    metadata = json.loads(metadata_path.read_text())
    if metadata.get("schema") != SYSTEM_SCHEMA:
        raise ValueError(f"unexpected system schema: {metadata_path}")
    system_path = (ROOT / metadata["system"]["msolve_input"]).resolve()
    if digest(system_path) != metadata["system"]["msolve_sha256"]:
        raise ArithmeticError(f"msolve input digest mismatch: {system_path}")
    lines = system_path.read_text().strip().splitlines()
    names = lines[0].split(",")
    prime = int(lines[1])
    equation_texts = "\n".join(lines[2:]).split(",\n")
    return metadata, system_path, names, prime, equation_texts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expanded", type=Path, default=DEFAULT_EXPANDED)
    parser.add_argument("--compact", type=Path, default=DEFAULT_COMPACT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()

    expanded = load(arguments.expanded)
    compact = load(arguments.compact)
    expanded_metadata, expanded_path, names, prime, expanded_texts = expanded
    compact_metadata, compact_path, compact_names, compact_prime, compact_texts = compact
    if compact_names != names or compact_prime != prime:
        raise ArithmeticError("compact and expanded rings differ")
    if expanded_metadata["input"] != compact_metadata["input"]:
        raise ArithmeticError("compact and expanded source inputs differ")
    if expanded_metadata["lambda"] != compact_metadata["lambda"]:
        raise ArithmeticError("compact and expanded lambda values differ")
    if len(expanded_texts) != len(compact_texts):
        raise ArithmeticError("compact and expanded equation counts differ")

    ring = PolynomialRing(GF(prime), names=names, order="degrevlex")
    expanded_polynomials = [ring(text.replace("^", "**")) for text in expanded_texts]
    compact_polynomials = [ring(text.replace("^", "**")) for text in compact_texts]
    mismatches = [
        index
        for index, (left, right) in enumerate(
            zip(expanded_polynomials, compact_polynomials), start=1
        )
        if left != right
    ]
    if mismatches:
        raise ArithmeticError(f"polynomial encoding mismatch: {mismatches}")

    output = {
        "schema": "elkies-k3.lattice-foundry-ns0007-compact-msolve-encoding-audit.v1",
        "status": "PASS_EXACT_POLYNOMIAL_ENCODING_IDENTITY",
        "prime": prime,
        "lambda": expanded_metadata["lambda"],
        "variables": names,
        "equation_count": len(expanded_polynomials),
        "expanded": {
            "metadata": relative(arguments.expanded.resolve()),
            "metadata_sha256": digest(arguments.expanded.resolve()),
            "msolve_input": relative(expanded_path),
            "msolve_sha256": digest(expanded_path),
        },
        "compact": {
            "metadata": relative(arguments.compact.resolve()),
            "metadata_sha256": digest(arguments.compact.resolve()),
            "msolve_input": relative(compact_path),
            "msolve_sha256": digest(compact_path),
        },
        "audit": {
            "equal_equations": len(expanded_polynomials),
            "mismatching_one_based_indices": mismatches,
            "expanded_term_counts": [
                len(polynomial.monomials()) for polynomial in expanded_polynomials
            ],
            "total_degrees": [
                int(polynomial.total_degree()) for polynomial in expanded_polynomials
            ],
        },
        "proved": (
            "Parsing both declared msolve inputs over the same GF(p) polynomial "
            "ring gives equal polynomials equation by equation. The compact final "
            "equations are therefore only an exact factored printing of the expanded "
            "system in Sage. This does not certify msolve parser or preprocessing "
            "equivalence; only fully expanded msolve inputs are used for the census."
        ),
        "reproduce": (
            "/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python "
            "elkies-k3/scripts/"
            "audit_lattice_foundry_ns0007_compact_msolve_encoding.sage"
        ),
    }
    serialized = json.dumps(output, indent=2, sort_keys=True) + "\n"
    output_path = arguments.output.resolve()
    if arguments.check:
        if output_path.read_text() != serialized:
            raise SystemExit("NS0007 compact-encoding audit artifact is stale")
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(serialized)
    print(
        "FOUNDRYNS0007ENCODING|"
        f"equations={len(expanded_polynomials)}|status="
        "PASS_EXACT_POLYNOMIAL_ENCODING_IDENTITY"
    )


if __name__ == "__main__":
    main()
