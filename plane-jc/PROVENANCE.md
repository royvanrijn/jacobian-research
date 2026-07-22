# Plane-JC external provenance

## Zenodo record 21479814 — audit target

Audit status on 2026-07-22: **externally reduced and locally reproduced.**
This status means that the archived exact identities were replayed, the
Newton-polygon transcription and all localization strata were checked, and
the degree conclusion was traced through the published normal-form reduction.
It is not external peer review or author review of this repository audit.

Pinned record metadata, obtained from the Zenodo Records API endpoint
`https://zenodo.org/api/records/21479814` on 2026-07-22:

| Field | Pinned value |
| --- | --- |
| Title | *Exact Computer-Assisted Exclusion of the (72,108) Frontier in the Two-Dimensional Jacobian Problem* |
| Creator | Billel Helali |
| Publication date | 2026-07-21 |
| Record creation time | 2026-07-21T20:01:47.665915+00:00 |
| Version | 1.0.1 |
| Version DOI | `10.5281/zenodo.21479814` |
| Concept DOI | `10.5281/zenodo.21479813` |
| Resource type | Software |
| Language | English |
| Zenodo metadata license | CC BY 4.0 |
| Internal archive licenses | CC BY 4.0 for manuscript, documentation, logs, and certificate data; MIT for Python, shell, and computer-algebra source |
| Related manuscript | Guccione–Guccione–Horruitiner–Valqui, arXiv:2204.14178 (`isSupplementTo`) |
| Related source repository | `https://github.com/bilLkarkariy/jc2-72-108-exact-certificates` (`isIdenticalTo`) |
| Journal submission | None identified in the record, manuscript, `CITATION.cff`, or linked repository |
| Zenodo archive | `bilLkarkariy/jc2-72-108-exact-certificates-v1.0.1.zip` |
| Archive size | 90,732,544 bytes |
| Zenodo archive checksum | MD5 `91255150c8c689b26dc6fb61f9d80aec` |

The record exposes one ZIP archive rather than separate top-level manuscript
and data files.  The exact attachment and its extracted snapshot are pinned
under `plane-jc/external/zenodo-21479814/`.  The attachment matched Zenodo's
MD5 and has SHA-256
`f7f0876de12d35badbed2be6a773d4a9dada50aff778c080126aab541deefcde`.
The archive is the GitHub source snapshot at commit
`d9ea4fd088eb7b0f3fe4e881a93ead88c9b4e4a7` and contains:

| Role | Archive member | Bytes | SHA-256 |
| --- | --- | ---: | --- |
| Audited manuscript | `paper/jc2_72_108_exact_exclusion.pdf` | 85,717 | `86e942f2f85b4dc88d49dbc34af695435dd2ee53ea03160b96dbe5187107d300` |
| Manuscript source | `paper/jc2_72_108_exact_exclusion.tex` | 17,218 | `449f2cce45414b64c873e613a74168719d5f72981ec9fae1e9101a1a5ee88f0d` |
| Complete replay bundle | `release_bundle/jc2_72_108_exact_replay_v1.0.1.zip` | 90,647,957 | `232204bdb598cc2ea0368e154c8573e18bbfdc69fa631c8878de4b884b38bb18` |
| Saved replay output | `verification/verify_all_final.out` | 3,512 | recorded inside the pinned outer archive |
| Certificate manifest | `verification/RECONSTRUCTED_CERTIFICATES.sha256` | 755 | recorded inside the pinned outer archive |

The complete replay bundle contains the generated Singular ideals, serialized
Nullstellensatz certificates, the 89,105,967-byte hard certificate, generation
scripts, and the independent `gmpy2` checker.  The principal certificate has
SHA-256
`0e48ffab32469ef8405a6945b16cf1521ddeb3c592ae4e5051968110a4dc656a`.

### Version relationship

Zenodo lists only this one deposited version (`relations.version.index=0`).
GitHub has tags `v1.0.0` (published 2026-07-21T19:41:30Z) and `v1.0.1`
(published 2026-07-21T20:01:38Z).  The PDF, TeX source, selected-certificate
archive, and complete replay assets are byte-identical across those tags; even
the replay ZIP retains the same SHA-256 under both filenames.  Version 1.0.1
adds the archival/Zenodo packaging and DOI metadata.  It is therefore a
packaging revision of the same independent 2026 manuscript, not a mathematical
revision of the 2022 Guccione--Guccione--Horruitiner--Valqui paper.

The latter is a separate prior manuscript, arXiv:2204.14178, titled
*Increasing the degree of a possible counterexample to the Jacobian Conjecture
from 100 to 108*, by Jorge Alberto Guccione, Juan José Guccione, Rodrigo
Horruitiner, and Christian Valqui.  The audited arXiv PDF has SHA-256
`ac18e80cc2391f204f73b908a6a6557eb1141d9fbb5ebdb9f6e0a22121db80bd`.

The deposit description itself limits the claimed computation to two
Laurent/Newton coefficient systems transcribed from Proposition 4.3 of
arXiv:2204.14178 and says that the degree conclusion remains conditional on
the correctness, exhaustiveness, and faithful transcription of that published
reduction.  This description is metadata, not proof, and will not be used as a
substitute for the manuscript, reduction theorems, or certificate replay.

## Local replay record

The complete archived command was replayed on 2026-07-22 and ended with
`JC2_72_108_EXACT_REPLAY_PASS`.  The environment was:

```text
Python 3.14.2
gmpy2 2.3.1
numpy 2.5.1
python-flint 0.9.0
sympy 1.14.0
Singular 4.4.1p5 / GMP 6.3.0 / NTL 11.6.0 / FLINT 3.6.0
```

The first-block Singular Gröbner/FGLM input was also rerun from scratch and
returned `DP_SIZE=56`, `LEX_SIZE=6`, and the same degree-35 eliminant.  On
Singular 4.4.1 the later `nfmodStd` scripts return an ideal whose first
generator is literally `1` but which also contains a redundant field
polynomial, so their version-sensitive `size==1` diagnostic prints `NONUNIT`.
This display discrepancy is not used: the proof artifacts are the explicit
identities checked by the serialized-certificate verifier and the independent
`gmpy2` checker.
