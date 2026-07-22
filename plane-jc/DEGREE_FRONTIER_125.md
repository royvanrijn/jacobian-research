# The plane degree frontier at 125

> **Status.**  The external coefficient identities have been replayed locally,
> and the reduction chain below has been checked against the cited primary
> manuscripts.  “Locally reproduced” here means a reproduction of the final
> algebraic exclusion conditional on the published normal-form and candidate
> enumeration theorems.  It is not a proof of JC(2), a claim that degree 125 is
> attainable, or a claim that both coordinate degrees are at least 125.

Exact artifact and version provenance precedes all claims in
[PROVENANCE.md](PROVENANCE.md).

## What the new manuscript actually proves

Helali's theorem is, verbatim (notation restored from the TeX source):

> **Theorem (Exact exclusion of the transcribed systems).** The two
> characteristic-zero coefficient systems transcribed from Cases 1 and 2 of
> the \((72,108)\) Laurent/Newton reduction have no solutions.

Its following corollary is explicitly conditional: if Proposition 4.3 of
Guccione--Guccione--Horruitiner--Valqui is exhaustive for the \((72,108)\)
configuration and the local normalization and transcription are faithful,
then the pair is impossible and

\[
  \max\{\deg P,\deg Q\}\ge 125.
\]

Thus the manuscript itself proves a statement about two transcribed systems,
not an unconditional new reduction theorem.  The local work in
[PAIR_72_108_REPRODUCTION.md](PAIR_72_108_REPRODUCTION.md) checks the missing
interface from the two published polygons to those systems.

## Repository statement

Let \(k\) be an algebraically closed field of characteristic zero and let
\((P,Q)\in k[x,y]^2\) satisfy \([P,Q]\in k^\times\).  If \((P,Q)\) is not a
polynomial automorphism, then the reproduced external reduction and exact
certificate imply

\[
  \boxed{\max\{\deg P,\deg Q\}\ge 125.}
\]

The same degree statement holds over any characteristic-zero field after
base change to an algebraic closure: a polynomial inverse over the closure,
if it existed, is unique and descends coefficientwise.  The Newton-polygon
papers normally work over an algebraically closed characteristic-zero field
after their initial reductions.  No positive-characteristic claim is made.

This is a lower bound for the **larger** coordinate degree.  It neither says
\(\deg P\ge125\) and \(\deg Q\ge125\), nor asserts that 125 occurs.

## Exact logical chain

1. **Arbitrary counterexample to a minimal standard pair.**  In *On the shape
   of possible counterexamples to the Jacobian Conjecture*, Guccione, Guccione,
   and Valqui define
   \(B=\min\gcd(\deg P,\deg Q)\) over counterexamples.  Their Corollary 5.21
   (label `B finito` in arXiv:1401.1784) says that if JC(2) is false, there is
   a minimal standard \((m,n)\)-pair with \(m,n>1\) coprime, preserving the
   relevant total degrees.  Proposition 4.5 supplies the standardizing
   automorphism; the definition of minimal pair is in Section 4.

2. **Minimal standard pair to an admissible complete chain.**  In *Some
   algorithms related to the Jacobian Conjecture*, Guccione, Guccione,
   Horruitiner, and Valqui, the theorem labelled `standard pair generates
   complete chain` in Section 2 and the divisibility conditions summarized
   immediately before the `Main algorithm` imply that the first corner of a minimal standard pair occurs
   in the output of the admissible-complete-chain algorithm for any bound
   \(M\ge B\).  Section 5 then converts each final corner into the allowed
   coprime \((m,n)\)-families.

3. **Admissible chain to the list below 125.**  The same paper's Section 7
   lists all 34 oriented possibilities with maximum degree at most 150.  The
   2022 paper's Theorem 2.1 combines that enumeration with the cited
   exclusions and proves:
   \[
   \text{counterexample}\Longrightarrow
   \max\{\deg P,\deg Q\}\ge125
   \quad\text{or}\quad
   (\deg P,\deg Q)\in\{(72,108),(108,72)\}.
   \]

4. **The remaining pair to two Laurent polygons.**  In the 2022 paper,
   Proposition 4.3 (the case \((8,28)\)) transforms the remaining standard
   pair into \(P,Q\in k[x,x^{-1},y]\) with \([P,Q]=x^2\) and exactly one of
   two displayed Newton polygons.  In the paper's orientation
   \(A_0=(8,28)\), \((m,n)=(3,2)\), so the original degrees are
   \((108,72)\); swapping the two coordinates gives \((72,108)\).

5. **Two Laurent polygons to contradiction.**  The coordinate change
   \(t=xy^2,z=y^{-1}\), the full lattice supports, the five bracket equations,
   all compatibility equations used, both sign branches, and the exact ideal
   certificates are reconstructed in the companion note.  Case 2 has an
   exact unit certificate.  Case 1 splits exhaustively as \(s=\pm c\); the
   pre-division \(h=0\) strata have unit certificates, while on \(h\ne0\) an
   exact identity forces \(h=0\).  Hence neither polygon can occur.

Every arrow is implication-only.  No assertion is made that every arbitrary
counterexample already has the displayed Laurent form without passing through
minimality and the standard-pair/chain theorems.

## Prior frontier and credit

The immediate (108)-frontier claim was checked in the complete 2022 GGHV
source, and its enumerative input in the complete 2017 algorithm source.
Heitmann's 1990 primary text states both
“\(\gcd(\deg f,\deg g)\ge16\)” and that his computation reaches the same four
cases as Moh but does **not** reproduce Moh's reduction-of-degree argument.
The publisher exposes Moh's introduction and bibliographic scan but restricts
the remaining article text; the claims attributed to its later sections below
are therefore cross-checked against the exact section/proposition citations in
the 2014, 2017, and 2022 primary papers rather than silently presented as a
fresh full rereading of the paywalled article.

| Date | Authors | Result used here | Hypotheses / normalization | Computational status |
| --- | --- | --- | --- | --- |
| 1983 | T. T. Moh | Proved \(\max(\deg P,\deg Q)>100\); only the maximum-64 case has a complete proof in print according to the later GGHV audit | Algebraically closed characteristic-zero field; approximate-root/minimal-counterexample analysis | Computer search followed by four reduction-of-degree cases; later authors flag that only one case is fully printed |
| 1990 | Raymond C. Heitmann | Proved \(\gcd(\deg P,\deg Q)\ge16\), recovered the same four computational cases as Moh, and developed structural corner restrictions; explicitly did not reprove Moh's \(>100\) theorem | Complex minimal counterexample / monomial valuations | Theorem 2.23 is by hand; later restrictions use a computer search |
| 2013--2017 | Jorge A. Guccione, Juan J. Guccione, Christian Valqui | Shape theorem, standard minimal pair, differential-equation exclusions, and lower-side restrictions | Characteristic zero; algebraic closure after reduction; standard \((m,n)\)-pairs | Exact algebra plus bounded corner computations |
| 2017 | Jorge A. Guccione, Juan J. Guccione, Rodrigo Horruitiner, Christian Valqui | Admissible-chain algorithm and all 34 oriented cases through maximum degree 150 | Minimal standard pair and complete-chain conditions | Published pseudocode and tables; no separate source-code archive located |
| 2019 | the same four authors | Approximate-root/intersection inequality used to remove the degree-84 configuration | Their standard-pair and \(\pi\)-root conventions | Theoretical |
| 2022 | the same four authors | Theorem 2.1: only \((72,108)\) and its reversal survive below 125; Proposition 4.3: two Laurent polygons for its remaining \((8,28)\) case | Plane counterexample reduced to their standard notation | Mix of theoretical reductions and exact/hand coefficient elimination |
| 2026 | Billel Helali | Exact certificates excluding the two transcribed Proposition 4.3 systems | Exact characteristic-zero number fields; faithful transcription required | Archived FLINT/SymPy/Singular generation plus explicit certificates |
| 2026, this repository | repository audit | Replayed the archived identities, reran the first-block Singular FGLM calculation, checked supports/normalization/branches/divisions, and supplied an independent checker and bounded frontier regression | Same published normal form; no author review claimed | Local reproduction, not external peer review |

The correct pre-2026 historical summary is therefore

\[
  \max\{\deg P,\deg Q\}\ge108,
\]

with \((72,108)\) (up to reversal) the sole unresolved degree pair whose
maximum is below 125.  The bound was not \(>108\): the unresolved pair itself
has maximum 108.  The 125 frontier depends on the whole earlier enumeration
and exclusion program and is not attributable to the 2026 certificate alone.

## The ten cases used for the below-125 theorem

The 2022 Theorem 2.1 starts from the following table (orientation suppressed
only in this repository summary):

| First corner \(A_0\) | \((m,n)\) | maximum degree | Disposition before/within 2022 |
| --- | --- | ---: | --- |
| \((4,12)\) | \((3,4)\) | 64 | Moh; Heitmann; GGV differential equation |
| \((4,12)\) | \((5,7)\) | 112 | GGV differential equation |
| \((5,20)\) | \((2,3)\) | 75 | GGV polynomial system |
| \((5,20)\) | \((3,2)\) | 75 | GGV polynomial system |
| \((7,21)\) | \((2,3)\) | 84 | 2022 via intersection inequality; second coefficient proof |
| \((8,24)\) | \((2,3)\) | 96 | 2017 algorithm paper, Proposition 6.1 |
| \((8,28)\) | \((3,2)\) | 108 | one subcase excluded in 2022; Proposition 4.3 subcase audited here |
| \((8,32)\) | \((3,2)\) | 120 | 2022 via lower-side Proposition 3.29 |
| \((9,24)\) | \((2,3)\) | 99 | 2022 coefficient/approximate-root argument |
| \((9,27)\) | \((2,3)\) | 108 | 2022 coefficient/approximate-root argument |

The two occurrences of maximum 108 correspond to different corner chains;
the \((9,27)\) chain was already excluded, leaving the \((8,28)\) chain.

## Dimension separation

The repository status remains:

\[
 n=1:\text{ true},\qquad
 n=2:\text{ open with degree and architecture restrictions},\qquad
 n\ge3:\text{ false via the repository's announced explicit construction}.
\]

The present note belongs only to the JC(2) constraint program.
