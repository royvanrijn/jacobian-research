# From Euler's restricted-prime sums to Eisenstein norms

**Unpublished · draft.**

This directory contains the manuscript, its current PDF, paper metadata and
review record, and the generated standalone source archive. The branch's
[research notes, checkers, receipts and supplied inputs](../../research/arithmetic/eisenstein-restricted-prime-sums/README.md)
live under `research/`.

The paper begins with the Gaussian restricted-prime theorem, explores its
Eisenstein variation and the mixed case forced by the prime 3, then derives
primitive-norm sumsets and sharp additive-basis bounds. It proves the exact
binary identities for `A+A` and `A+(3·A)`, with a written large-target proof
and complete finite certificates. The current PDF retains the balanced finite
experiment; the separate [angular and finite-polynomial proofs](../../research/arithmetic/EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md)
now establish eventual balance without supplying an explicit onset.

[PDF](paper.pdf) · [LaTeX](main.tex) · [Metadata](metadata.json) ·
[Source archive](source.zip) · [Verification and review](verification.json)

## Build the manuscript

From this directory, run:

```sh
latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
```

The verifier appendix reads its canonical source from the research directory.
In the standalone source archive, the same unchanged verifier is included
locally. Both layouts build the same paper. Numerical replay instructions,
software requirements, output conventions and the archived Java-replay limits
are in the research README linked above.

## Build the standalone source archive

Run the following from this directory. It includes the manuscript sources,
paper records, verification programs and their finite evidence. Exploratory
research notes and the original supplied input archive remain in `research/`.
A generated archive is a paper distribution snapshot; the canonical research
files have one current home.

```sh
python3 - <<'PY'
from pathlib import Path
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED
paper = Path.cwd()
research = (paper / '../../research/arithmetic/eisenstein-restricted-prime-sums').resolve()
code_home = research if research.is_dir() else paper
members = {p.name: p for p in (paper / 'main.tex', paper / 'references.tex',
           paper / 'README.md', paper / 'metadata.json', paper / 'verification.json')}
for p in (paper / 'sections').glob('*.tex'):
    members[p.relative_to(paper).as_posix()] = p
for p in (*code_home.glob('*.py'), *code_home.glob('*.java')):
    members[p.name] = p
for p in (code_home / 'receipts').rglob('*.json'):
    members[p.relative_to(code_home).as_posix()] = p
with ZipFile(paper / 'source.zip', 'w') as archive:
    for name, path in sorted(members.items()):
        info = ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
        info.compress_type = ZIP_DEFLATED
        info.external_attr = 0o100644 << 16
        archive.writestr(info, path.read_bytes())
PY
```

## Review and remaining work

The review record distinguishes written arguments, bounded diagnostics, full
finite replays and internal PDF-only reading. No external specialist review,
external replay, formal verification, public deposit or established novelty is
claimed. The research relocation changes file locations and packaging only.

Next: specialist review of the Hecke/Jacquet--Langlands normalization and the
ramified calculation, a complete priority audit, external replay, and a
persistent source-archive identifier before submission. The newer qualitative
balanced result awaits integration into a later manuscript revision. Proving
that `8,44` are the complete ordinary balanced exception list still requires
an effective onset and finite verification covering the remaining range.
