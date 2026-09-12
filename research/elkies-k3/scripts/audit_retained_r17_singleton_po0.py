#!/usr/bin/env python3
"""Audit the retained top-200 singleton records without regenerating any shell.

The strict campaign auditor expects the generation-time exporter and process
wrapper. Both were archived during the shared-runtime migration. Resolve only
those two source pins to their verified archived bytes while the unchanged
auditors check all recorded data hashes, coverage, aggregates and boundaries.
No archived program, finite-field enumeration or Hensel calculation is run.
This is a retained-record audit, not independent mathematical replay.
"""

from __future__ import annotations

from contextlib import contextmanager
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "elliptic-curves/cas"))
from research_runtime.witnesses import retained_source
import audit_r17_norm12_11952_singleton_po0_two_prime_top200 as audit
import run_r17_norm12_11952_singleton_po0_top150 as campaign


SOURCE_NAMES = (
    "elkies-k3/scripts/export_elkies_2026_twist_polynomial_sections_modp.sage",
    "elkies-k3/scripts/run_twist_polynomial_sections_bruteforce.py",
)


@contextmanager
def source_digest_view(module, sources: dict[Path, Path]):
    """Resolve source pins and hash each unchanged file once during this audit."""
    original = module.digest
    cache = {}

    def signature(path):
        stat = path.stat()
        return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns

    def recorded_digest(path):
        resolved = sources.get(path.resolve(), path)
        before = signature(resolved)
        key = resolved, before
        if key not in cache:
            value = original(resolved)
            if signature(resolved) != before:
                raise ArithmeticError(f"input changed while hashing: {resolved}")
            cache[key] = value
        return cache[key]

    module.digest = recorded_digest
    try:
        yield
    finally:
        module.digest = original


def main() -> None:
    if len(sys.argv) != 1:
        raise SystemExit("This audit takes no arguments and never writes output.")
    stored = json.loads(audit.OUTPUT.read_text())
    # The programs interpreting the retained records must still be exactly the
    # audited generation. Only the two nonexecuted producers may be relocated.
    for module in (audit, campaign):
        path = Path(module.__file__).resolve()
        name = str(path.relative_to(ROOT))
        if sha256(path.read_bytes()).hexdigest() != stored["inputs"][name]:
            raise ArithmeticError(f"retained-record auditor changed: {name}")
    primary = json.loads(audit.PRIMARY.read_text())
    secondary = json.loads(audit.SECONDARY.read_text())
    sources = {}
    for name in SOURCE_NAMES:
        expected = primary["inputs"][name]
        if secondary["inputs"][name] != expected:
            raise ArithmeticError(f"campaign source versions disagree: {name}")
        sources[(ROOT / name).resolve()] = retained_source(ROOT, name, expected)
    argv = sys.argv
    try:
        sys.argv = [str(audit.__file__), "--check"]
        with source_digest_view(campaign, sources):
            audit.main()
    finally:
        sys.argv = argv
    print("PASS retained records; two generation-time sources resolved; no shell or lift rerun")


if __name__ == "__main__":
    main()
