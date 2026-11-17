#!/usr/bin/env python3
"""Compile local polynomial specifications to JSON and generated Lean."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shlex
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from jcsearch.arithmetic_certificate import (  # noqa: E402
    build_certificate,
    write_certificate,
    write_lean_instantiation,
)


DEFAULT_SPEC = ROOT / "arithmetic/specifications/ramified_quintic.json"
DEFAULT_CERTIFICATE = (
    ROOT / "artifacts/generated-results/arithmetic_keller_quintic.json"
)
DEFAULT_LEAN = (
    ROOT
    / "formal/finite-etale-keller/FiniteEtaleKeller/"
    "GeneratedArithmeticQuintic.lean"
)


def display_path(path: Path) -> str:
    """Prefer a repository-relative path in recorded regeneration commands."""
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spec", type=Path, default=DEFAULT_SPEC)
    parser.add_argument("--certificate", type=Path, default=DEFAULT_CERTIFICATE)
    parser.add_argument("--lean", type=Path, default=DEFAULT_LEAN)
    parser.add_argument(
        "--lean-module",
        help="override the generated Lean namespace recorded in the certificate",
    )
    parser.add_argument(
        "--stable-parameter",
        type=int,
        help="override the specification with a nonnegative stable-family parameter",
    )
    arguments = parser.parse_args()
    if arguments.stable_parameter is not None and arguments.stable_parameter < 0:
        parser.error("--stable-parameter must be nonnegative")
    if arguments.stable_parameter is not None and arguments.certificate == DEFAULT_CERTIFICATE:
        parser.error(
            "--certificate is required when --stable-parameter is supplied"
        )

    specification = json.loads(arguments.spec.read_text(encoding="utf-8"))
    if arguments.lean_module is not None:
        specification["lean_module"] = arguments.lean_module
    if arguments.stable_parameter is not None:
        specification["stable_parameter"] = arguments.stable_parameter
        specification["name"] = (
            f"{specification['name']}-stable-k{arguments.stable_parameter}"
        )
    if specification.get("stable_parameter") is not None:
        if "lean_module" not in specification:
            parser.error(
                "a stable compilation requires lean_module in the specification "
                "or --lean-module"
            )
        if arguments.lean == DEFAULT_LEAN:
            parser.error(
                "--lean is required for a stable compilation to avoid "
                "overwriting the minimal-gauge specialization"
            )
    certificate, polynomial = build_certificate(specification)
    if (
        arguments.spec != DEFAULT_SPEC
        or arguments.certificate != DEFAULT_CERTIFICATE
        or arguments.stable_parameter is not None
    ):
        command = [
            ".venv/bin/python",
            "scripts/compile_arithmetic_keller_certificate.py",
            "--spec",
            display_path(arguments.spec),
            "--certificate",
            display_path(arguments.certificate),
        ]
        if arguments.stable_parameter is not None:
            command.extend(
                ["--stable-parameter", str(arguments.stable_parameter)]
            )
        if arguments.lean_module is not None:
            command.extend(["--lean-module", arguments.lean_module])
        if "lean_instantiation" in certificate:
            command.extend(["--lean", display_path(arguments.lean)])
        certificate["provenance"]["regeneration_command"] = " ".join(
            shlex.quote(part) for part in command
        )
    write_certificate(certificate, arguments.certificate)
    if "lean_instantiation" in certificate:
        write_lean_instantiation(certificate, polynomial, arguments.lean)

    print(f"WROTE: {arguments.certificate}")
    if "lean_instantiation" in certificate:
        print(f"WROTE: {arguments.lean}")
    print(
        "MAP_SHA256:",
        certificate["keller_map"]["expanded_map_hash"]["digest"],
    )


if __name__ == "__main__":
    main()
