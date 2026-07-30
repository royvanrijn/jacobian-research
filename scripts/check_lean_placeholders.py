#!/usr/bin/env python3
"""Audit placeholders and explicit axioms in every checked-in Lean project."""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEAN_ROOT = ROOT / "formal"
PLACEHOLDER = re.compile(r"\b(sorry|admit)\b")
AXIOM = re.compile(r"\baxiom\s+([A-Za-z_][A-Za-z0-9_']*)")
ALLOWED_AXIOMS = {
    (
        Path("formal/gmc2/GMC2/DuistermaatVanDerKallen.lean"),
        "duistermaat_van_der_kallen",
    ),
    (
        Path("formal/gmc2/GMC2/Specialization.lean"),
        "exists_goodReduction",
    ),
}


def strip_comments_and_strings(source: str) -> str:
    """Replace nested Lean comments and strings while preserving line numbers."""
    out: list[str] = []
    i = 0
    block_depth = 0
    in_line_comment = False
    in_string = False
    escaped = False

    while i < len(source):
        ch = source[i]
        nxt = source[i + 1] if i + 1 < len(source) else ""

        if in_line_comment:
            if ch == "\n":
                in_line_comment = False
                out.append("\n")
            else:
                out.append(" ")
            i += 1
            continue

        if block_depth:
            if ch == "/" and nxt == "-":
                block_depth += 1
                out.extend((" ", " "))
                i += 2
            elif ch == "-" and nxt == "/":
                block_depth -= 1
                out.extend((" ", " "))
                i += 2
            else:
                out.append("\n" if ch == "\n" else " ")
                i += 1
            continue

        if in_string:
            if ch == "\n" and not escaped:
                # Keep line numbering stable even for malformed literals.
                out.append("\n")
                in_string = False
            else:
                out.append(" ")
                if escaped:
                    escaped = False
                elif ch == "\\":
                    escaped = True
                elif ch == '"':
                    in_string = False
            i += 1
            continue

        if ch == "-" and nxt == "-":
            in_line_comment = True
            out.extend((" ", " "))
            i += 2
        elif ch == "/" and nxt == "-":
            block_depth = 1
            out.extend((" ", " "))
            i += 2
        elif ch == '"':
            in_string = True
            out.append(" ")
            i += 1
        else:
            out.append(ch)
            i += 1

    return "".join(out)


def main() -> int:
    failures: list[str] = []
    found_axioms: set[tuple[Path, str]] = set()
    files = sorted(
        path
        for path in LEAN_ROOT.rglob("*.lean")
        if ".lake" not in path.relative_to(LEAN_ROOT).parts
    )
    if not files:
        print(f"ERROR: no Lean files found under {LEAN_ROOT}", file=sys.stderr)
        return 2

    for path in files:
        cleaned = strip_comments_and_strings(path.read_text(encoding="utf-8"))
        relative_path = path.relative_to(ROOT)
        for match in PLACEHOLDER.finditer(cleaned):
            line = cleaned.count("\n", 0, match.start()) + 1
            failures.append(f"{path.relative_to(ROOT)}:{line}: forbidden `{match.group(1)}`")
        for match in AXIOM.finditer(cleaned):
            line = cleaned.count("\n", 0, match.start()) + 1
            declaration = (relative_path, match.group(1))
            found_axioms.add(declaration)
            if declaration not in ALLOWED_AXIOMS:
                failures.append(
                    f"{relative_path}:{line}: unexpected explicit axiom `{match.group(1)}`"
                )

    for path, name in sorted(ALLOWED_AXIOMS - found_axioms):
        failures.append(f"{path}: missing documented explicit axiom `{name}`")

    if failures:
        print("Lean source audit failed:", file=sys.stderr)
        print("\n".join(failures), file=sys.stderr)
        return 1

    print(
        f"PASS: {len(files)} Lean files contain no sorry/admit placeholders; "
        f"the {len(found_axioms)} explicit axioms match the documented GMC(2) allowlist"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
