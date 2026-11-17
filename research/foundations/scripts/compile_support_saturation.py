#!/usr/bin/env python3
"""Compile one support-saturation problem from the shared JSON schema."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from jcsearch.support_saturation import (  # noqa: E402
    CompilerOptions,
    SupportSaturationCompiler,
    SupportSaturationProblem,
    certificate_json,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        type=Path,
        help="support-saturation-input.v1 JSON file",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="write the certificate here instead of standard output",
    )
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    problem = SupportSaturationProblem.from_mapping(payload)
    options = CompilerOptions(**payload.get("compiler_options", {}))
    certificate = SupportSaturationCompiler(options).compile_problem(problem)
    serialized = certificate_json(certificate)
    if args.output is None:
        sys.stdout.write(serialized)
        return
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialized)
    print(
        json.dumps(
            {
                "input_sha256": certificate["problem_sha256"],
                "output": str(args.output),
                "schema": certificate["schema"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
