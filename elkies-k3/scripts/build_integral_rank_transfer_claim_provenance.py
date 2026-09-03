#!/usr/bin/env python3
"""Build and check the integral-rank-transfer claim provenance artifact.

The human-readable table is authoritative for classification.  This builder
also reads the canonical theorem headings, so a renamed, added, or omitted
label fails closed instead of silently leaving the provenance ledger stale.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
THEOREMS = ROOT / "elkies-k3" / "RANK_MUTATION_AND_LIFT_THEOREMS.md"
MAP = ROOT / "elkies-k3" / "LITERATURE_AND_NOVELTY_MAP_2026-09-03.md"
BIB = ROOT / "elkies-k3" / "references" / "integral-rank-transfer.bib"
OUTPUT = (
    ROOT
    / "artifacts/generated-results"
    / "elkies-k3-integral-rank-transfer-claim-provenance-v1.json"
)

CLASSES = {
    "ESTABLISHED",
    "TAILORED_COROLLARY",
    "LIKELY_NEW_ALGORITHM",
    "NEW_COMPUTATION",
    "OPEN_CONJECTURE",
}

SOURCE_LOCI = {
    "SHIODA-MW": {
        "citation_keys": ["Shioda1990MWLattices", "SchuettShioda2010EllipticSurfaces"],
        "locus": "Shioda (1990), Mordell--Weil quotient, height and discriminant formula; Schuett--Shioda (2010), Sections 6 and 11",
    },
    "SHIODA-GAL": {
        "citation_keys": ["Shioda1989GaloisI", "Shioda1989GaloisIII", "SchuettShioda2019MWLattices"],
        "locus": "Shioda (1989) III, Section 5; Schuett--Shioda (2019), chapter Galois Representations and Algebraic Equations",
    },
    "NIKULIN": {
        "citation_keys": ["Nikulin1980IntegralForms"],
        "locus": "Nikulin (1980), Propositions 1.4.1 and 1.5.1; indefinite uniqueness theorem where explicitly invoked",
    },
    "KN": {
        "citation_keys": ["Nishiyama1996JacobFibrations", "BraunKimuraWatari2015Classification"],
        "locus": "Nishiyama (1996), frame-complement method; Braun--Kimura--Watari (2015), Sections 2--3 and Proposition C'",
    },
    "HOP": {
        "citation_keys": ["BrandhorstElkies2023Lehmer", "Kumar2014KummerFibrations", "ElkiesKumar2014Hilbert"],
        "locus": "Brandhorst--Elkies, Section 2, Lemmas 2.5--2.6, Remark 2.10, Sections 2.3--2.4; Kumar, Appendix; Elkies--Kumar, Section 5",
    },
    "KNEIGH": {
        "citation_keys": ["Chenevier2022Statistics", "Voight2023Neighbors"],
        "locus": "Chenevier, equations (1.2)--(1.3), Examples 5.2--5.4, Theorem 5.9, and Remarks 5.10--5.11; Voight, Section 3 and Theorem 3.18",
    },
    "THETA": {
        "citation_keys": ["BruinierStein2009Weil", "KaneKim2026Cosets", "Mueller2024WeilBasis"],
        "locus": "Weil representations and Hecke operators; lattice-coset theta series and neighbour algorithms; Weil-representation basis problem",
    },
    "PASTEN-SALGADO": {
        "citation_keys": ["PastenSalgado2024NonThin"],
        "locus": "Pasten--Salgado (2024), Theorem 1.1",
    },
    "TATE": {
        "citation_keys": ["Kubert1976UniversalTorsion"],
        "locus": "Kubert (1976), classical marked-point Tate normal form; the displayed two-point Bézout packaging is proved directly",
    },
    "MASS": {
        "citation_keys": ["ConwaySloane1982LowDimensional"],
        "locus": "Minkowski--Siegel mass formula; Conway--Sloane (1982), neighbour enumeration and exact mass closure",
    },
    "ROOT-MASS": {
        "citation_keys": ["King2003RootlessMass"],
        "locus": "King (2003), Proposition 1 and Sections 7--9: prescribed-root-system mass inversion from Siegel representation averages",
    },
    "CT-CLASS": {
        "citation_keys": ["ChenevierTaibi2026Rank29"],
        "locus": "Chenevier--Taibi (2026), inductive classification, isometry invariants, automorphism data, and complete-list verification",
    },
    "ELKIES-2026": {
        "citation_keys": ["Elkies2026R17I"],
        "locus": "Elkies, arXiv:2608.25406, abstract and announced construction sequel",
    },
}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_bibliography_keys() -> None:
    keys = re.findall(r"^@\w+\{([^,]+),", BIB.read_text(), flags=re.MULTILINE)
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate bibliography key")
    required = {
        key for source in SOURCE_LOCI.values() for key in source["citation_keys"]
    }
    missing = sorted(required - set(keys))
    if missing:
        raise ValueError(f"missing bibliography keys: {missing}")


def headings() -> dict[str, dict[str, object]]:
    pattern = re.compile(
        r"^(?P<marks>#{3,4}) "
        r"(?P<kind>Theorem|Corollary|Proposition|Lemma|Exact finite control) "
        r"(?P<id>[A-Z][A-Za-z0-9.-]*): (?P<title>.+)$"
    )
    found: dict[str, dict[str, object]] = {}
    for line_number, line in enumerate(THEOREMS.read_text().splitlines(), 1):
        match = pattern.match(line)
        if not match:
            continue
        claim_id = match.group("id")
        if claim_id in found:
            raise ValueError(f"duplicate theorem label {claim_id}")
        found[claim_id] = {
            "kind": match.group("kind"),
            "current_wording": match.group("title"),
            "line": line_number,
        }
    return found


def ledger_rows() -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    in_table = False
    for line in MAP.read_text().splitlines():
        if line.startswith("| ID | Canonical term |"):
            in_table = True
            continue
        if not in_table or not line.startswith("|"):
            if in_table and rows:
                break
            continue
        if line.startswith("| ---"):
            continue
        cells = [cell.strip() for cell in line.strip().strip("|").split("|")]
        if len(cells) != 6:
            raise ValueError(f"malformed ledger row: {line}")
        claim_id = cells[0].strip("`")
        novelty_class = cells[2].strip("`")
        if novelty_class not in CLASSES:
            raise ValueError(f"unknown novelty class for {claim_id}: {novelty_class}")
        if claim_id in rows:
            raise ValueError(f"duplicate ledger row {claim_id}")
        rows[claim_id] = {
            "canonical_literature_term": cells[1],
            "novelty_class": novelty_class,
            "strongest_source": cells[3],
            "project_adaptation": cells[4],
            "new_axes": cells[5],
        }
    return rows


def statement_preambles(theorem_rows: dict[str, dict[str, object]]) -> dict[str, str]:
    """Return the exact text from each label through its proof boundary."""
    lines = THEOREMS.read_text().splitlines()
    ordered = sorted(
        ((int(row["line"]), claim_id) for claim_id, row in theorem_rows.items())
    )
    result: dict[str, str] = {}
    for index, (line_number, claim_id) in enumerate(ordered):
        end = ordered[index + 1][0] - 1 if index + 1 < len(ordered) else len(lines)
        block = lines[line_number:end]
        for offset, line in enumerate(block):
            if re.match(r"^#{3,5} Proof(?:$|\s)", line) or line.startswith("## "):
                block = block[:offset]
                break
        result[claim_id] = "\n".join(block).strip()
    return result


def build() -> dict[str, object]:
    check_bibliography_keys()
    theorem_rows = headings()
    preambles = statement_preambles(theorem_rows)
    map_rows = ledger_rows()
    missing = sorted(set(theorem_rows) - set(map_rows))
    extra = sorted(set(map_rows) - set(theorem_rows))
    if missing or extra:
        raise ValueError(f"claim coverage mismatch: missing={missing}, extra={extra}")

    claims = []
    for claim_id, heading in theorem_rows.items():
        row = map_rows[claim_id]
        cited_loci = sorted(key for key in SOURCE_LOCI if f"`{key}`" in row["strongest_source"])
        has_no_antecedent_note = "antecedent" in row["strongest_source"].lower()
        if not cited_loci and not has_no_antecedent_note:
            raise ValueError(
                f"{claim_id} needs a primary-source locus or an explicit antecedent note"
            )
        axes = {key: False for key in ("mathematical", "algorithmic", "computational", "certification")}
        for short, key in {
            "math": "mathematical",
            "alg": "algorithmic",
            "comp": "computational",
            "cert": "certification",
        }.items():
            axes[key] = short in {part.strip() for part in row["new_axes"].split(",")}

        likely_new = row["novelty_class"] == "LIKELY_NEW_ALGORITHM"
        claims.append(
            {
                "repository_id": claim_id,
                "kind": heading["kind"],
                "current_wording": heading["current_wording"],
                "canonical_literature_term": row["canonical_literature_term"],
                "canonical_statement": {
                    "path": str(THEOREMS.relative_to(ROOT)),
                    "line": heading["line"],
                    "hypotheses_and_statement": preambles[claim_id],
                    "statement_sha256": hashlib.sha256(
                        preambles[claim_id].encode()
                    ).hexdigest(),
                    "scope_rule": "No broader scope is imported from the cited source.",
                },
                "strongest_source": row["strongest_source"],
                "resolved_source_loci": cited_loci,
                "project_adaptation": row["project_adaptation"],
                "novelty_class": row["novelty_class"],
                "novelty_axes": axes,
                "confidence": "medium" if likely_new else "high",
                "unresolved_prior_art_query": (
                    "No explicit antecedent was located in the sources checked; search MathSciNet/zbMATH, theses, software, and the announced Elkies construction sequel before any priority claim."
                    if likely_new
                    else None
                ),
            }
        )

    return {
        "schema": "elkies-k3-integral-rank-transfer-claim-provenance-v1",
        "repository_snapshot_audited": "fb056eaa153f70da49e5268b6cefaca039b9b819",
        "generated_from": {
            "theorem_note": str(THEOREMS.relative_to(ROOT)),
            "theorem_note_sha256": digest(THEOREMS),
            "literature_map": str(MAP.relative_to(ROOT)),
            "literature_map_sha256": digest(MAP),
            "bibliography": str(BIB.relative_to(ROOT)),
            "bibliography_sha256": digest(BIB),
        },
        "claim_classes": sorted(CLASSES),
        "source_loci": SOURCE_LOCI,
        "publication_risk": "The ledger records search outcomes, not proof of priority. Recheck arXiv:2608.25406 and its announced sequel and obtain specialist review before publication.",
        "claims": claims,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text() != rendered:
            raise SystemExit(f"stale provenance artifact: run {Path(__file__).relative_to(ROOT)}")
        print(f"OK: {OUTPUT.relative_to(ROOT)}")
        return
    OUTPUT.write_text(rendered)
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
