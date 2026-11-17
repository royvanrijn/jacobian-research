#!/usr/bin/env python3
"""Fast local adapter for the archived (72,108) hard-certificate replay.

The archived source remains byte-for-byte pinned.  This adapter reuses its
checker, changes only the decimal-to-integer parser to GMP's native parser,
and then performs the archive's branch-transport checks without replaying the
same 89 MB certificate a second time.
"""

from __future__ import annotations

import importlib.util
import os
import pickle
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REPLAY = (
    REPO_ROOT
    / "plane-jc"
    / "external"
    / "zenodo-21479814"
    / "bilLkarkariy-jc2-72-108-exact-certificates-d9ea4fd"
    / "release_bundle"
    / "exact_replay"
)
REPLAY = Path(os.environ.get("JC2_72_108_REPLAY_ROOT", DEFAULT_REPLAY)).resolve()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load archived checker {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_hard_replay() -> None:
    hard = load_module(
        "_jc2_archived_hard_verifier",
        REPLAY / "hard" / "verify_certificate_gmpy2.py",
    )

    def parse_field(expression: str):
        if expression == "0":
            return tuple(hard.zero() for _ in range(5))
        out = [hard.zero() for _ in range(5)]
        for item in expression.split("+"):
            match = hard.TERM.match(item)
            if not match:
                raise RuntimeError(
                    f"cannot parse certificate coefficient: {item[:120]}"
                )
            numerator = hard.mpz(match.group(1))
            denominator = hard.mpz(match.group(2) or 1)
            degree = int(match.group(3) or (1 if "*w" in item else 0))
            if not 0 <= degree < 5:
                raise RuntimeError(
                    "certificate coefficient is not reduced modulo the minpoly"
                )
            out[degree] += hard.mpq(numerator, denominator)
        return tuple(out)

    hard.parse_field = parse_field
    hard.main()


def verify_branch_transport() -> None:
    sys.path.insert(0, str(REPLAY))
    symmetry = load_module(
        "_jc2_archived_branch_symmetry",
        REPLAY / "verify_hne0_branch_symmetry.py",
    )
    branch1 = [
        symmetry.decode_poly(q)
        for q in pickle.loads(symmetry.BRANCH1.read_bytes())
    ]
    branch2 = [
        symmetry.decode_poly(q)
        for q in pickle.loads(symmetry.BRANCH2.read_bytes())
    ]
    for index, (left, right, scale) in enumerate(
        zip(branch1, branch2, symmetry.ROW_SCALES)
    ):
        transported = symmetry.sign_substitution(left, (1, -1, -1)) * scale
        if transported.terms != right.terms:
            raise RuntimeError(
                f"degree-five branch symmetry fails on equation {index}"
            )
    print("SYSTEM_SYMMETRY_PASS")

    supports = []
    with symmetry.CERTIFICATE.open() as handle:
        for line in handle:
            if line.startswith("C|"):
                parts = line.split("|", 5)
                supports.append(tuple(map(int, parts[1:5])))
    if len(supports) != 385 or any(
        not 1 <= item[0] <= 4 for item in supports
    ):
        raise RuntimeError("unexpected support in branch-1 certificate")
    print("BRANCH1_EXACT_IDENTITY_PASS")
    print("BRANCH2_EXACT_IDENTITY_PASS")


def main() -> None:
    run_hard_replay()
    verify_branch_transport()


if __name__ == "__main__":
    main()
