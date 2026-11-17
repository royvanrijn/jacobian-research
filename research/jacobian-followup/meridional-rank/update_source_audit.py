#!/usr/bin/env python3
"""Reconcile this programme's sourced family and 2026 dataset audit entries."""

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    path = HERE / "sources.json"
    catalogue = json.loads(path.read_text())
    rows = {r["id"]: r for r in catalogue["sources"]}
    rows["BZ85"].update({"url": "https://doi.org/10.1007/BF02567413",
                        "audit": "Original publication metadata checked; Theorem1.1 explicitly restated as Musick2012 Theorem1 and applied to the independently identified sixteen five-tangle Montesinos knots. Original proof not fully audited.",
                        "scope": "Montesinos links with r>=3 noninteger rational tangles have meridional rank=bridge=r; integer twist parameter allowed. Exact later generalized-Montesinos variants still need their own hypotheses."})
    rows["BHAT26"].update({"audit": "All243 algebraic/four-seed certificates pass separate standard-library replay and exact named-diagram matching to an independent HT census. All243 PDs match literally in KnotInfo snapshots2026.10.5 and2026.11.2; both uncertainty sets equal the author243+293 partition. Full paper and publication priority not independently audited.",
                          "scope": "Author-deposited243 cases with mr=bridge=4 via transvections/reflections are independently replayed as named knot types;293 cases remain in that study. Both archived KnotInfo snapshots still list all536 as bridge[3,4], so database lag is distinguished from the certificate proof.",
                          "status": "Unrefereed author release; credit retained with author, named finite certificates independently checked, full paper not certified."})
    added = [
        {"id": "MUSICK12", "title": "Minimal bridge projections for 11-crossing prime knots", "authors": "Chad Musick", "year": 2012,
         "url": "https://arxiv.org/abs/1208.4233v3", "location": "Page1 Theorem1 and page2 Montesinos definition",
         "audit": "Relevant first two pages read; explicit Boileau-Zieschang theorem statement extracted. Full86-page diagram appendix not audited.",
         "scope": "Prime11crossing bridge classification and explicit Montesinos r>=3 noninteger-tangle bridge/meridional theorem statement; used as primary later source for the known-family audit."},
        {"id": "KL04", "title": "On the classification of rational tangles", "authors": "Kauffman Lambropoulou", "year": 2004,
         "url": "https://arxiv.org/html/math/0311499v2", "doi": "10.1016/j.aam.2003.06.002", "location": "Theorem1, Lemma2(iii), continued-fraction constructions",
         "audit": "Definitions, classification statement and flipping lemma read; full article not independently audited.",
         "scope": "Rational-tangle fraction classification and invariance under a planar half-turn; justifies the fixed half-turn variants in the exact five-Montesinos identity audit."},
        {"id": "HT-CENSUS", "title": "SnapPy independently packaged Hoste-Thistlethwaite named link census", "authors": "SnapPy and Spherogram contributors; Hoste Thistlethwaite census", "year": 2026,
         "url": "https://snappy.computop.org/censuses.html", "versions": {"snappy": "3.3.2", "spherogram": "2.4.1", "snappy_15_knots": "1.2.1", "snappy_manifolds": "1.4"},
         "audit": "Selected named PD/DT records and original installed census file hashes retained. All16 five-bridge certificates transferred directly; all243 October diagrams matched by exact dart bijections. No numerical isometry acceptance.",
         "scope": "Second-source name/diagram verification for the retained finite cases, not a whole-census or software correctness audit."},
        {"id": "KNOTINFO-SNAPSHOTS", "title": "Database KnotInfo snapshots 2026.10.5 and 2026.11.2", "authors": "KnotInfo curators; Sebastian Oehms snapshot packaging", "year": 2026,
         "url": "https://pypi.org/project/database-knotinfo/2026.11.2/", "versions": ["2026.10.5", "2026.11.2"],
         "audit": "Pinned wheel SHA256 and CSV SHA256 checked; all536 selected name/bridge/PD rows retained read-only. In both versions exactly536 thirteen-crossing entries are[3,4], exactly the pinned author partition; all243 certificate PDs agree literally.",
         "scope": "Archived source status, not a claim about every later live KnotInfo update or an independent proof of lower bounds."}
    ]
    for row in added:
        if row["id"] in rows:
            assert rows[row["id"]] == row
        else:
            catalogue["sources"].append(row)
    catalogue["unresolved_scope_extractions"] = [
        "Exact generalized Montesinos parameter restrictions in LM93", "Exact many-twig hypotheses in BBKM23",
        "Original low-rank orbifold proof audit", "Reported Dutra subgroup gap and any subsequent repair",
        "Full proof audits of November and October2026 releases; October243 algebra, named-diagram and archived KnotInfo checks completed",
        "Priority scope of fixed-fibre meridional bounds in the quandle literature"]
    path.write_text(json.dumps(catalogue, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": "UPDATED", "sources": len(catalogue["sources"]), "sharp16_family": "all already-known five-tangle Montesinos"}))


if __name__ == "__main__":
    main()
