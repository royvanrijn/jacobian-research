# Completion audit for Zenodo 21479814

This ledger maps the original mission to concrete repository artifacts.  It
also records limitations that must not be hidden by the phrase “locally
reproduced.”

| Objective | Status | Evidence / limitation |
| --- | --- | --- |
| Exact source provenance | complete | [PROVENANCE.md](PROVENANCE.md) pins the sole Zenodo attachment, PDF, source snapshot, versions, relations, licenses, and hashes. |
| Exact theorem and scope | complete | [DEGREE_FRONTIER_125.md](DEGREE_FRONTIER_125.md) quotes the manuscript theorem and conditional corollary and states the characteristic-zero, larger-coordinate interpretation. |
| Arbitrary-to-standard reduction chain | complete as an external reduction audit | Every implication is tied to the labelled GGV/GGHV theorem used; this repository does not reprove those long normal-form theorems. |
| Historical frontier | complete with a disclosed source-access qualification | The 2014--2022 primary sources were read in full.  Heitmann's primary text was inspected; Moh's publisher exposes only the opening scan without authentication, so later section claims are cross-checked through exact citations and caveats in the later primary papers. |
| \((72,108)\) elimination | complete | [PAIR_72_108_REPRODUCTION.md](PAIR_72_108_REPRODUCTION.md) reconstructs supports, bracket layers, branches, fields, ideals, and final identities.  Case 1 is correctly described as a complete necessary truncation, not a full infinite-band system. |
| Primary exact CAS replay | complete | The pinned `exact_replay/verify_all.sh` run ended `JC2_72_108_EXACT_REPLAY_PASS`; commands and versions are in [cas/README.md](cas/README.md). |
| Independent certificate | complete | [cas/verify_h_certificate_independent.py](cas/verify_h_certificate_independent.py) independently checks the hard identity and now rejects any system or certificate input whose SHA-256 is not pinned. |
| Division and saturation audit | complete | Section 7 of the reproduction note accounts for every inverse, factor split, discarded equation, and the complementary \(h=0\) strata. |
| Status and source integration | complete | `STATUS.md`, the root index, reproduction guide, marked-point note, source index, and search policy are updated with dimensionally scoped language. |
| Newton--boundary dictionary | complete as a qualified comparison | [NEWTON_BOUNDARY_DICTIONARY.md](NEWTON_BOUNDARY_DICTIONARY.md) marks exact, standard, interpretive, and speculative correspondences. |
| Conceptual replacement | investigated; negative result | No intersection, conductor, ramification, determinant, or Smith-normal-form proof is established.  The order-seven grading/value-semigroup question remains repository-original conjectural work. |
| Frontier 125--150 | complete as a deterministic published-table regression | [NEXT_DEGREE_FRONTIER.md](NEXT_DEGREE_FRONTIER.md) records all 13 pairs and all 24 chain realizations.  It does not claim to independently reimplement the complete-chain generator or to have derived the new Laurent systems. |
| Family extension | tested; no theorem | Ratio \(2:3\) and gcd data do not suffice.  The \((96,144)\) repeated-tail case is the first proposed reuse experiment. |
| Search-policy consequence | complete | [SEARCH_POLICY.md](SEARCH_POLICY.md) requires admissibility and boundary filters before coefficient solving and reserves sub-frontier searches for tests and automorphisms. |

Accordingly, the degree-125 consequence is locally reproduced in the precise
conditional sense used by this directory.  The next research objectives are
deriving the six \((96,144)\) Laurent systems and testing whether the repeated
\((8,28)\to(11/4,7)\) tail carries the old valuation obstruction.
