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
consequence from three canonical pairs to two.  It also needs a terminology
explanation.  In the usual convention `A_n` has `2n` associated-graded
polynomial generators.  The four displayed Poisson coordinates have the
symbol size of `A_2`.  The abstract does not say that their direct
quantization automatically gives `A_4`: it advertises a separate construction
using the same four polynomials as commuting base-coordinate multiplication
operators, together with four Hamiltonian duals.  Read as an
inverse-Jacobian/cotangent construction, that
produces eight Weyl generators and therefore lands naturally in `A_4`.  The
full appendix is still needed to check that this is the intended construction
and that its operator identities and non-surjectivity proof are correct.

The screenshot alone is not a certificate.  To admit the claim, the audit
needs the exact formulas for `T,D,S` and must independently verify

\[
 \{D,R\}=1,\quad \{S,T\}=1,
 \quad \{R,S\}=\{R,T\}=\{D,S\}=\{D,T\}=0,
\]

plus a displayed collision or another proof of nonsurjectivity.  The proposed
Weyl lift must then be checked at the operator level; principal-symbol
identities alone do not prove all Weyl commutators.

## An exact foundational fingerprint

One part of the provenance question can be answered without guessing the
missing outputs.  Put

\[
 \begin{aligned}
 X&=x,\\
 Z&=3x^2p+(2-6xq)z,\\
 Y&=q-\frac{xZ}{3},\\
 E&=\frac{1+3xq}{2}p-3q^2z.
 \end{aligned}                                           \tag{1}
\]

This is a polynomial automorphism of `Q[x,q,p,z]`, with determinant `-1`.
Indeed, writing `q=Y+XZ/3`, its inverse is

\[
 p=-6EXq+2E+3Zq^2,
 \qquad
 z=-3EX^2+\frac32ZXq+\frac12Z.                         \tag{2}
\]

For the foundational map `F=(F_1,F_2,F_3)`, substitution of (1) gives

\[
 F_3(X,Y,Z)=2X-3X^2Y-X^3Z=x(2-3xq)=R.                 \tag{3}
\]

Thus the displayed `R` is an exact algebraic fingerprint of the foundational
counterexample, not merely a coincident degree or fiber count.  It is strong
evidence that the announced construction is derived from that map.  It is not
yet a proof that the unavailable `T,D,S` are the manuscript's particular
lift, nor does it identify the author or source.

There is also a useful negative result.  The obvious invariant choices

\[
 S_0=\frac12F_1(X,Y,Z),\qquad T_0=F_2(X,Y,Z)
\]

satisfy exactly

\[
 \{S_0,T_0\}=1,
 \qquad \{R,S_0\}=\{R,T_0\}=0,
 \qquad \{E,R\}=1.                                    \tag{4}
\]

Nevertheless there is **no** polynomial correction
`D=E+f(X,Y,Z)` for which both `{D,S_0}` and `{D,T_0}` vanish.  Since (1) is a
polynomial coordinate system and the kernel of `{-,R}` is `Q[X,Y,Z]`, this
covers every possible polynomial partner `D` for this particular pair
`S_0,T_0`.

Here is the exact obstruction.  With `Q=Y+XZ/3`, the first completion
equation can be integrated explicitly; the two remaining equations reduce to

\[
 \left(-3X^2\partial_X+(6XQ-2)\partial_Q\right)h
 =\frac{Q^3}{2}
  \left(54Q^3X^3+189Q^2X^2+222QX+89\right).           \tag{5}
\]

After localizing at `X`, set `v=1/X` and
`rho=2X-3X^2Q`.  The derivation on the left sends `v` to `3` and `rho` to
zero.  Its forced antiderivative, modulo a polynomial in `rho`, becomes

\[
 \frac{54Q^6X^6+234Q^5X^5+375Q^4X^4+270Q^3X^3
       +90Q^2X^2+24QX+4}{60X^4}.                      \tag{6}
\]

A polynomial in `rho` has no negative power of `X`, so it cannot cancel the
term `1/(15X^4)`.  Equation (5) has no polynomial solution.  This proves that
the manuscript's missing formulas, if correct, contain a genuinely different
choice or construction; they cannot responsibly be filled in by the naive
foundational substitution.  The exact certificate is
[`verify_rank_two_poisson_preaudit.py`](../scripts/verify_rank_two_poisson_preaudit.py).

## What remains unavailable

The complete audit must identify the paper and compare its actual formulas
with (1)--(6).  Generic degree three and a three-point fiber alone still do
not establish polynomial left--right equivalence or derivation.

No public source matching the exact text was located on 22 July 2026 after
exact-phrase and formula searches, arXiv metadata searches, and inspection of
the public Omniscience Project papers index.  Its separately published `A_3`
paper is not this manuscript.  Until the missing formulas or source are
supplied, the rank-two claim remains a high-priority **external announced
manuscript under audit**, not an active theorem, not an attributed Long
result, and not a dependency of the boundary work.  The locally proved
fingerprint and obstruction are repository results, not substitutes for the
external audit.  The broader consequence comparison is in
[External consequences and provenance](EXTERNAL_CONSEQUENCES_AND_PROVENANCE.md).
