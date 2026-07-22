# Quadratic ladder and rank-two Poisson audit

This note audits the two supplied cards separately.  Their evidentiary status
is different: the quadratic marked-root ladder is exactly identifiable from
the displayed formulas, while the rank-two Poisson card omits most of the
polynomials needed to check its claim.

## The ladder is the `m=1` cancellation column

To avoid conflicting indices, write `ell` for the rung on the card.  Set

\[
 D(s)=1-Bs+As^2,
 \qquad
 P_\ell(\lambda)=-J_\ell\int_0^\lambda D(s)^{\ell-1}\,ds-C.
\]

Then

\[
 P_\ell'(\lambda)=-J_\ell D(\lambda)^{\ell-1},
 \qquad \deg P_\ell=2\ell-1.
\]

This is precisely the generalized cancellation inverse polynomial with

\[
 (m_{\rm cancel},r)=(1,\ell-1),
 \qquad N=r(m_{\rm cancel}+1)+1=2\ell-1.
\]

Thus the card supplies a useful geometric interpretation, not a disjoint new
family: a target is a degree-`2 ell - 1` polynomial with one simple root
marked, and its two stationary points have multiplicity `ell - 1`.  For
`ell=2`, substituting `A=3a`, `B=b`, and `C=-3c`, then dividing by three,
gives exactly

\[
 c-2s+bs^2-2as^3,
\]

the repository's foundational marked-root cubic.

This identification has two useful consequences.

1. The generic degree and inequivalence of different rungs are already
   covered by the cancellation degree theorem.
2. The full source-boundary intersection diagram is certified for every
   rung, not merely in a bounded experiment: the `m=1` branch polynomial is
   irreducible, and the critical-contact coefficient is not an associate, so
   their resultant is nonzero.

The exact regression is
[`verify_counterexample_ladder.py`](../scripts/verify_counterexample_ladder.py).

## What is new in the Poisson card

The first card claims a standard rank-two Poisson endomorphism on
`Q[x,q,p,z]`, with four output polynomials `(R,T,D,S)`, and displays only

\[
 R=x(2-3xq).
\]

If correct, this would improve the repository's current cotangent-lift
consequence from three canonical pairs to two.  It would also require a
terminology correction: with two canonical pairs it is a rank-two Weyl
algebra in the usual convention (four polynomial generators in the associated
graded algebra), not a fourth Weyl algebra unless “fourth” is explicitly
being used to count generators.

The screenshot alone is not a certificate.  To admit the claim, the audit
needs the exact formulas for `T,D,S` and must independently verify

\[
 \{D,R\}=1,\quad \{S,T\}=1,
 \quad \{R,S\}=\{R,T\}=\{D,S\}=\{D,T\}=0,
\]

plus a displayed collision or another proof of nonsurjectivity.  The proposed
Weyl lift must then be checked for ordering corrections: commuting principal
symbols do not by themselves prove that the corresponding Weyl operators
commute.

No public source matching the exact text was located by exact-phrase search
on 22 July 2026.  Until the missing formulas or source are supplied, this is
a high-priority candidate audit, not an active theorem or a dependency of the
boundary work.
