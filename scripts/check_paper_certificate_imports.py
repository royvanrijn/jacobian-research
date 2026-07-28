#!/usr/bin/env python3
"""Keep both publication Lean certificates at their intended boundaries."""

from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_DIR = (
    ROOT
    / "formal"
    / "finite-etale-keller"
    / "FiniteEtaleKeller"
)
EXPECTED = {
    "PaperCertificate.lean": (
        "FiniteEtaleKeller.GeneralGaugeMap",
        "FiniteEtaleKeller.GeneralGaugeFunctionFieldComparison",
        "FiniteEtaleKeller.PageOneTheorem",
        "FiniteEtaleKeller.AbstractFiniteEtale",
    ),
    "FixedHassePaperCertificate.lean": (
        "FiniteEtaleKeller.FixedHasseArithmetic",
    ),
}


def main() -> int:
    for filename, expected in EXPECTED.items():
        certificate = CERTIFICATE_DIR / filename
        imports = tuple(
            line.removeprefix("import ").strip()
            for line in certificate.read_text(encoding="utf-8").splitlines()
            if line.startswith("import ")
        )
        if imports != expected:
            print(f"{filename} import boundary changed.")
            print(f"Expected: {expected}")
            print(f"Actual:   {imports}")
            return 1
    print("PASS: both paper certificates have their publication imports only")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
