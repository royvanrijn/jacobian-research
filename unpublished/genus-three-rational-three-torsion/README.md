# A genus-three hyperelliptic curve with rational Jacobian three-torsion of rank three

This standalone paper extracts the genus-three existence result from
[HG1 and HG5](../../research/arithmetic/GENUS_THREE_THREE_MARKED_ELLIPTIC_QUOTIENTS.md).
Its main theorem gives an explicit hyperelliptic genus-three curve over
\(\mathbb Q\) with \(J(\mathbb Q)[3]\simeq(\mathbb Z/3)^3\), answering the
precise existence question following Levin–Yan–Wiljanen Theorem 9.
Its second theorem transfers multiplication-by-\(3^n\) Selmer groups and
three-primary Tate–Shafarevich groups for every quadratic twist through
dual degree-\(2^9\) isogenies. No finiteness assumption on Sha is used.

[Manuscript PDF](genus-three-rational-three-torsion.pdf) ·
[LaTeX source](main.tex) · [Source and checker bundle](genus-three-rational-three-torsion-source.zip) ·
[Metadata](metadata.json) · [Verification record](verification.json) ·
[Checksums](SHA256SUMS.txt)

The manuscript is prepared for specialist review. The mathematical claims
retain the repository's written-proof status and internal checking scope.
External review, independent replay, formal verification and historical
priority are unestablished. The PDF-only blind-reader test is unperformed.
There is no submission or public deposit.

This extraction is the current presentation priority for HG1/HG5.
Further class-group counts, conductor refinements, cubic covers or
higher-rank constructions are not prerequisites for reviewing it.
The current research proof remains the source of the broader HG1–HG5
results; this manuscript presents the existence result and the transfer
for this one curve. It makes no claim about the full rational torsion
group or geometric simplicity.

## Build and replay

From this directory, or the root of the extracted source ZIP:

```sh
SOURCE_DATE_EPOCH=1795737600 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=genus-three-rational-three-torsion main.tex
sage -python verify_construction.py
sha256sum -c SHA256SUMS.txt
```

The source requires standard AMS, Geometry, Microtype, Xcolor and Hyperref
packages. The exact checker was run with SageMath 10.9. It has no
external inputs, seed or search bound: it works in rational function
fields over \(\mathbb Q\), checks squarefreeness over \(\mathbb F_7\),
and verifies the marked points using exact elliptic-curve arithmetic.
It includes a repeated-root negative control. No new number-field
search or extension calculation is part of this extraction.

The existing broader replay also passed:

```sh
sage -python ../../research/arithmetic/verify_genus_three_marked_elliptic_quotients.py
```

That replay checks additional local, conductor and unit fixtures which
are outside this paper. Neither checker computes all Selmer groups;
the universal transfer is a written Kummer argument in Section 4.

## Specialist review brief

The review should first decide the existence theorem on its own merits,
then review the transfer theorem. The current dependency chain is
squarefree octic → conic model → three degree-two quartic quotients →
three degree-four elliptic maps → nonzero independent order-three
pullbacks → the Weil-pairing upper bound. The transfer uses character
orthogonality, duality and local Kummer functoriality.

1. Reconstruct the function fields of all three quotients in Section 2,
   especially the simultaneous-sign quotient \((x/y,-w/(3y^2))\).
   Check extension over poles on smooth projective models and the
   degree-four argument in Lemma 3; each composite must have degree eight.
2. Verify the exact divisor of \(Y+h_iX+h_i+1\), nonvanishing of its
   divided pullback, and all signs in (9). Check that the three character
   projectors separate the classes over \(\mathbb F_3\), and that the
   principally polarized Weil pairing gives the upper bound three.
3. Check the off-diagonal push–pull vanishing, the passage from finite
   kernel to an isogeny, duality under the canonical polarizations and
   the degree \(8^3\). No pairwise nonisogeny calculation is needed here.
4. For Theorem 2, check descent of both homomorphisms on every quadratic
   twist, inverses on torsion and on local point quotients, all local
   Kummer conditions including three and infinity, and the passage to
   the full three-primary Sha without a finiteness hypothesis.
5. Compare this exact model and its isomorphism/isogeny classes with
   later literature before any historical-priority claim. Distinguish
   rational three-torsion from full torsion and geometrically simple
   Jacobians. The classical gluing mechanism is credited to HLP.

The inspected primary comparison is narrow and explicit:

| Source | Inspected statement | Relevance |
|---|---|---|
| [Levin–Yan–Wiljanen, arXiv:1910.12276v1](../../research/arithmetic/sources/levin-yan-wiljanen-large-three-rank.pdf) | Theorem 9, printed pp. 12–14; question on printed p. 13 | Gives rank three in genus four and asks for hyperelliptic genus three with rational three-rank three. |
| [Howe–Leprévost–Poonen](../../research/arithmetic/sources/howe-leprevost-poonen-split-jacobians.pdf) | Proposition 14 and quotient calculation in its proof, PDF pp. 21–24 | Classical hyperelliptic genus-three gluing at vanishing twisting factor. The manuscript verifies its particular model directly. |
| [Müller–Reitsma, arXiv:2211.03372v4](../../research/arithmetic/sources/muller-reitsma-genus-three-torsion-2211.03372v4.pdf) | Theorem 1.1 and Remark 1.2, pp. 1–2 | Exact full torsion and geometric simplicity are stronger, separate questions. |

These sources and their digests are also indexed in the existing
[research reference catalogue](../../research/REFERENCES.md).
The bounded source comparison establishes the cited question and the
scope distinction. It does not settle historical priority.

For the PDF-only reader test, supply only the PDF to a reader without
workspace or conversation context. Record the exact PDF hash, supplied
context, contribution as understood, conceptual mechanism, hardest step,
limitations and verification obstacles in `verification.json`. The
current author-side agent has workspace context and cannot serve as that
blind reader. No specialist has been contacted or attributed a review.
