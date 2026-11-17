#!/usr/bin/env python3
"""Verify completed splitting of the cut-14 and cut-21 conormal towers."""

from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from research_degree42_ritt_completed_presentations import (  # noqa: E402
    completed_presentation_output,
    presentation_cache,
)
from research_degree42_tensor_extension import parse_sections  # noqa: E402


ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "degree42_ritt_completed_splits.json"
)
SECTORS = {
    (2, 3, 7): {
        "omission": 14,
        "lambda": (
            "6*u^15",
            "6*u^8",
            "3*u",
            "1",
            "0",
            "0",
            "0",
        ),
        "correction_component": 6,
        "correction": "-3*u^2+2*zeta",
        "relation_identity": (
            "zeta*s(e4)=-r4-u*zeta*r2"
        ),
    },
    (3, 2, 7): {
        "omission": 21,
        "lambda": (
            "12*u^22",
            "6*u^8",
            "2*u",
            "1",
            "0",
            "0",
            "0",
        ),
        "correction_component": 7,
        "correction": "-4*u^3+8*u*zeta",
        "relation_identity": (
            "3*zeta*s(e4)="
            "-r4-u*zeta*r1-(3*u^2*zeta-4*zeta^2)*r2"
        ),
    },
}


def completed_split_audit(word: tuple[int, int, int]) -> dict[str, object]:
    """Reduce the explicit section identities against both presentations."""

    output = completed_presentation_output(word)
    sections = parse_sections(
        output, ("PRESENTATION_TOTAL", "PRESENTATION_SPECTATOR")
    )
    total = sections["PRESENTATION_TOTAL"]
    spectator = sections["PRESENTATION_SPECTATOR"]
    assert total and spectator
    data = SECTORS[word]
    assignments = []
    correction_row = data["correction_component"]
    for column, coefficient in enumerate(data["lambda"], 1):
        assignments.append(f"A[4,{column}]={coefficient};")
        assignments.append(
            f"A[{correction_row},{column}]="
            f"({data['correction']})*({coefficient});"
        )
    program = f"""
ring b=0,(tau,zeta),dp;
poly u=1+tau;
module PT={",".join(total)};
module PS={",".join(spectator)};
matrix A[7][7];
{"".join(assignments)}
matrix identity[7][7];
identity=identity+1;
module sectionRelations=A*PS;
module retractionResidual=A-identity;
module wellDefinedRemainder=simplify(reduce(sectionRelations,std(PT)),2);
module retractionRemainder=simplify(reduce(retractionResidual,std(PS)),2);
print("COMPLETED_SPLIT_AUDIT");
print(size(PT));
print(size(PS));
print(size(wellDefinedRemainder));
print(size(retractionRemainder));
"""
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=300,
    )
    if "? error" in result.stdout:
        raise RuntimeError(result.stdout)
    compact = " ".join(
        line.strip()
        for line in result.stdout.splitlines()
        if line.strip() and not line.strip().startswith("//")
    )
    match = re.search(
        r"COMPLETED_SPLIT_AUDIT ([0-9]+) ([0-9]+) ([0-9]+) ([0-9]+)",
        compact,
    )
    assert match is not None, result.stdout
    values = tuple(map(int, match.groups()))
    assert values[2:] == (0, 0)
    return {
        "word": word,
        "thick_composite_omission": data["omission"],
        "presentation_generators": {
            "total": values[0],
            "spectator": values[1],
        },
        "spectator_normal_form": {
            "module": "Q[[tau]]*e4",
            "zeta_action": "0",
            "free_generator_coordinates": list(data["lambda"]),
        },
        "completed_section": {
            "generator_image": (
                f"e4+({data['correction']})*e{correction_row}"
            ),
            "relation_identity": data["relation_identity"],
            "well_defined_remainder_generators": values[2],
            "retraction_remainder_generators": values[3],
        },
        "extension_obstruction": {
            "class": (
                "[id_S] in coker(Hom_B(S,T) -> End_B(S))"
            ),
            "value": "0",
        },
    }


def main() -> None:
    sectors = [
        completed_split_audit(word)
        for word in SECTORS
    ]
    caches = [presentation_cache(word) for word in SECTORS]
    finite_tower = (
        ROOT
        / "artifacts"
        / "generated-results"
        / "degree42_ritt_inverse_limit_sections_q5_q7.json"
    )
    finite = json.loads(finite_tower.read_text())
    assert finite["schema"] == (
        "degree42-ritt-inverse-limit-sections-q5-q7.v1"
    )
    output = {
        "schema": "degree42-ritt-completed-splits.v1",
        "status": "exact completed module splitting theorem",
        "completed_base_ring": "B=Q[[tau,zeta]]",
        "presentation_caches": [
            {
                "path": str(cache.relative_to(ROOT)),
                "sha256": hashlib.sha256(cache.read_bytes()).hexdigest(),
            }
            for cache in caches
        ],
        "finite_tower_artifact": str(finite_tower.relative_to(ROOT)),
        "finite_tower_sha256": hashlib.sha256(
            finite_tower.read_bytes()
        ).hexdigest(),
        "sectors": sectors,
        "inverse_limit_consequence": (
            "The polynomial completed sections reduce to one compatible "
            "section at every base order. Hence the affine section-torsor "
            "inverse-limit obstruction and its lim^1 class vanish."
        ),
        "sector_asymmetry": (
            "Before presentation-first tensoring, cut 21 has one extra "
            "fourth-jet quotient-ring sector dimension. The completed "
            "first-conormal tensor modules instead have the common Hilbert "
            "functions K_q=3q-3, T_q=4q-3, S_q=q through q=2,...,7. Thus "
            "the extra dimension is non-flat base-change Tor. A genuine "
            "labelled asymmetry remains in the splitting formulas: cut 14 "
            "uses the correction (-3*(1+tau)^2+2*zeta)*e6, whereas cut 21 "
            "uses (-4*(1+tau)^3+8*(1+tau)*zeta)*e7."
        ),
        "theorem_boundary": (
            "This proves completed splitting and truncation compatibility "
            "for the two rotated first-conormal extensions. It does not "
            "construct the restriction maps between different factor "
            "charts, so full braid-cell coefficient coherence remains a "
            "separate calculation."
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_degree42_ritt_completed_splits.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: cut-14 completed extension splits")
    print("PASS: cut-21 completed extension splits")
    print("PASS: inverse-limit section obstructions vanish")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
