# A genus-three hyperelliptic curve with rational Jacobian three-torsion of rank three

**Hexagon submission under review.** Reserved `hexagon:2611.01278v1`. Updated November 28, 2026.

This standalone paper extracts the genus-three existence result from
[HG1 and HG5](../../research/arithmetic/GENUS_THREE_THREE_MARKED_ELLIPTIC_QUOTIENTS.md).
Its main theorem gives an explicit hyperelliptic genus-three curve over
\(\mathbb Q\) with \(J(\mathbb Q)[3]\simeq(\mathbb Z/3)^3\), answering the
precise existence question following Levin–Yan–Wiljanen Theorem 9.
Its second theorem transfers multiplication-by-\(3^n\) Selmer groups and
three-primary Tate–Shafarevich groups for every quadratic twist through
dual degree-\(2^9\) isogenies. No finiteness assumption on Sha is used.
Theorem 3 determines all eight rational points and the full finite group
\(J(\mathbb Q)=\mathbb Z/3\oplus\mathbb Z/6\oplus\mathbb Z/12\),
of order 216, together with conductor 243000. It is the manuscript
version of HG6.

[Manuscript PDF](genus-three-rational-three-torsion.pdf) ·
[LaTeX source](main.tex) · [Source and checker bundle](genus-three-rational-three-torsion-source.zip) ·
[Metadata](metadata.json) · [Verification record](verification.json) ·
[Checksums](SHA256SUMS.txt)

The mathematical claims retain the repository's written-proof status and
internal checking scope. External specialist review, independent replay,
formal verification and historical priority remain open. The PDF-only
blind-reader test is unperformed.

This extraction is the current presentation priority for HG1/HG5.
Further class-group counts, conductor refinements, cubic covers or
higher-rank constructions are not prerequisites for reviewing it.
The current research proof remains the source of the broader HG1–HG5
results; this manuscript presents the existence result, the transfer
and the full rational arithmetic theorem for this one curve. It makes
no geometric simplicity claim. The broader Sha group and quadratic-twist
rational points are not computed.

## Build and replay

From this directory, or the root of the extracted source ZIP:

```sh
SOURCE_DATE_EPOCH=1795737600 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=genus-three-rational-three-torsion main.tex
sage -python verify_construction.py
sage -python verify_arithmetic.py
sha256sum -c SHA256SUMS.txt
```

The source and checker ZIP preserves the bundled preparation records.
The submission status and updated checksums are recorded beside the ZIP.

The source requires standard AMS, Geometry, Microtype, Xcolor and Hyperref
packages. The exact checker was run with SageMath 10.9. It has no
external inputs, seed or search bound: it works in rational function
fields over \(\mathbb Q\), checks squarefreeness over \(\mathbb F_7\),
and verifies the marked points using exact elliptic-curve arithmetic.
It includes a repeated-root negative control. No new number-field
search or extension calculation is part of this extraction.

The arithmetic checker obtains proved elliptic rank bounds with the rank
database disabled, enumerates the first elliptic target and its fibres,
checks the principal-divisor certificates, exhausts the sixteen-element
Sylow group at seven, and checks the square-class obstruction at 17.
Its final block also checks the marked-norm coefficient and involution
tangent matrices in the separate
[deformation note](../../research/arithmetic/GENUS_THREE_MARKED_TORSION_DEFORMATIONS.md).
Those local deformation calculations are outside the manuscript's theorems.
They produce no rational geometrically simple example or infinite family.

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
6. Audit Theorem 3 separately: proved elliptic rank bounds, the complete
   rational fibres including the pole exclusions, the finite Sylow group
   at seven, the rational Pell certificate, and Poonen–Schaefer's descent
   map at 17. In particular, its nontrivial image is an obstruction to
   halving; no injectivity of fake descent is assumed. The order bound
   432 is strict: the full rational order is 216.

The inspected primary comparison is narrow and explicit:

| Source | Inspected statement | Relevance |
|---|---|---|
| [Levin–Yan–Wiljanen, arXiv:1910.12276v1](../../research/arithmetic/sources/levin-yan-wiljanen-large-three-rank.pdf) | Theorem 9, printed pp. 12–14; question on printed p. 13 | Gives rank three in genus four and asks for hyperelliptic genus three with rational three-rank three. |
| [Howe–Leprévost–Poonen](../../research/arithmetic/sources/howe-leprevost-poonen-split-jacobians.pdf) | Proposition 14 and quotient calculation in its proof, PDF pp. 21–24 | Classical hyperelliptic genus-three gluing at vanishing twisting factor. The manuscript verifies its particular model directly. |
| [Müller–Reitsma, arXiv:2211.03372v4](../../research/arithmetic/sources/muller-reitsma-genus-three-torsion-2211.03372v4.pdf) | Theorem 1.1 and Remark 1.2, pp. 1–2 | Exact full torsion and geometric simplicity are stronger, separate questions. |
| [Poonen–Schaefer](../../research/arithmetic/sources/poonen-schaefer-explicit-descent.pdf) | Section 5, Proposition 5.1 and (8), printed pp. 8–9 | The even-degree descent homomorphism kills doubles and rational infinities; its image obstructs halving at 17. |

These sources and their digests are also indexed in the existing
[research reference catalogue](../../research/REFERENCES.md).
The bounded source comparison establishes the cited question and the
scope distinction. It does not settle historical priority.

## Construction priorities from the primary literature

A focused literature check on 27 November 2026 inspected the actual statements
below. The search located no later resolution of the Müller–Reitsma exceptions;
this bounded observation does not establish their present open status.

| Precise target | What is available here | Remaining step |
|---|---|---|
| Levin–Yan–Wiljanen, question after Theorem 9: a hyperelliptic genus-three curve over $\mathbb Q$ with rational Jacobian three-rank three | HG1 supplies such a curve. HG6 also determines its full torsion as [3,6,12]. | Independent proof review and a priority comparison; no additional existence construction is needed. |
| Müller–Reitsma, Theorem 1.1: exact full torsion [3,3,3] on a geometrically simple hyperelliptic genus-three Jacobian | HG7 gives marked norm equations and rational tangent directions. | An actual rational smooth solution, a geometric-simplicity certificate, and exclusion of all additional torsion. Formal deformations alone do not suffice. |
| The same theorem: exact full torsion [6,6] | HG11 already supplies a geometrically simple family with two independent rational three-torsion classes. | Preserve those classes while producing exactly two independent rational two-torsion classes, then bound the entire torsion group. The existing HG11 congruence subfamily has exact torsion [3,3], so a different specialization or family is required. |
| A genus-three rational-Weierstrass construction with three independent virtual directions | HG11–HG12 supply two directions and simultaneous local virtuality. | A third independent class and compatible local conditions. This is a construction target motivated by the class-group application, rather than a separately verified published open question. |

The Sha comparison needs a further distinction. Theorem 1.5 of
[Flynn–Shnidman, with an appendix by Fisher](https://people.maths.ox.ac.uk/flynn/arts/art42.pdf)
proves arbitrarily large elementary $p$-torsion in Sha for absolutely simple
abelian varieties over $\mathbb Q$, for every prime $p$. Its introduction
explicitly credits the earlier $p=3$ case to Cassels. Consequently, existence
of a subgroup $(\mathbb Z/3)^{2r}$ for arbitrary $r$, without further
restrictions on the variety, is not a new open question.
The paper also poses a different question in the paragraph following
Theorem 1.5: order-$p$ Sha elements for varieties with Mumford–Tate group
$\mathrm{GSp}_{2g}$ and irreducible mod-$p$ representation. A split gluing
does not meet those hypotheses. Neither does HG11 at $p=3$, because its
rational three-torsion makes that representation reducible.

MCT16's prescribed alternating blocks, the exact entire three-primary
group in its later refinements, and HG14's rank-zero genus-three Jacobians
are narrower constructions. They should be compared under their exact rank,
dimension, polarization, splitting and group-structure hypotheses before
any novelty claim. A subgroup assertion does not determine the entire
three-primary Sha. The present literature check imports no new theorem
into HG1–HG6 and upgrades no review assurance.

For the PDF-only reader test, supply only the PDF to a reader without
workspace or conversation context. Record the exact PDF hash, supplied
context, contribution as understood, conceptual mechanism, hardest step,
limitations and verification obstacles in `verification.json`. The
current author-side agent has workspace context and cannot serve as that
blind reader. No specialist has been contacted or attributed a review.
