# Support-saturation principle

Let \(S\) be a Noetherian ring, let \(\mathfrak a\subset S\) be an ideal,
and let \(M\) be a finite \(S\)-module.  Write
\[
 H^0_{\mathfrak a}(M)
 =\{m\in M:\mathfrak a^n m=0\text{ for some }n\}.
\]
This module is the universal obstruction to extending a vanishing statement
from \(\operatorname{Spec}S\setminus V(\mathfrak a)\) across
\(V(\mathfrak a)\).

## The theorem

> **Theorem 1 (support-saturation principle).**  The following conditions
> are equivalent:
>
> 1. \(H^0_{\mathfrak a}(M)=0\);
> 2. \(0:_M\mathfrak a^\infty=0\);
> 3. no associated prime of \(M\) contains \(\mathfrak a\);
> 4. \(\operatorname{grade}(\mathfrak a,M)\ge1\);
> 5. \(\mathfrak a\) contains an \(M\)-regular element.
>
> If \(F\) is finite free and \(M=F/N\), these are also equivalent to
> \[
>  N:_F\mathfrak a^\infty=N.                              \tag{1}
> \]

**Proof.**  Noetherianity makes the ascending chain
\[
 0:_M\mathfrak a\subseteq0:_M\mathfrak a^2\subseteq\cdots
\]
stationary, and its union is \(H^0_{\mathfrak a}(M)\).  This proves the
equivalence of the first two conditions.  The associated primes of this
submodule are exactly
\[
 \operatorname{Ass}(M)\cap V(\mathfrak a).
\]
Thus it vanishes exactly when no associated prime of \(M\) contains
\(\mathfrak a\).  Prime avoidance for the finite set
\(\operatorname{Ass}(M)\) then identifies this with the existence of an
\(M\)-regular element in \(\mathfrak a\), equivalently positive grade.
Finally, an element \(f+N\in F/N\) is killed by a power of
\(\mathfrak a\) exactly when
\(\mathfrak a^n f\subseteq N\) for some \(n\), which gives
\[
 H^0_{\mathfrak a}(F/N)
 \simeq (N:_F\mathfrak a^\infty)/N
\]
and proves (1).  \(\square\)

## Defects and completion

> **Corollary 2 (extension across a support).**  If
> \(d\in M\) vanishes after localization at every prime outside
> \(V(\mathfrak a)\), then
> \[
>  d\in H^0_{\mathfrak a}(M).
> \]
> Consequently any condition in Theorem 1 forces \(d=0\).

**Proof.**  The cyclic module \(Sd\) has support in
\(V(\mathfrak a)\).  Hence
\(\mathfrak a\subseteq\sqrt{\operatorname{Ann}(d)}\), and finite
generation of \(\mathfrak a\) gives
\(\mathfrak a^n d=0\) for some \(n\).  \(\square\)

Let \(\widehat S\) be a flat adic completion of \(S\).  Since the
annihilator chain stabilizes, flat base change gives
\[
 H^0_{\mathfrak a}(M)\otimes_S\widehat S
 \simeq
 H^0_{\mathfrak a\widehat S}(M\otimes_S\widehat S).         \tag{2}
\]
In particular, the precompletion saturation equality (1) is sufficient
for the completed module to have no \(\mathfrak a\)-torsion.

## Structural shortcuts

Theorem 1 is deliberately weaker than flatness.  Useful sufficient
conditions are:

- \(M\) is torsion-free over an integral base and
  \(\mathfrak a\) contains a nonzero base element;
- \(M\) has no embedded associated primes and every minimal component
  avoids \(V(\mathfrak a)\);
- \(M\) is a finite maximal Cohen--Macaulay module over a regular base:
  Auslander--Buchsbaum makes it locally free, hence torsion-free.

The last statement requires module-finiteness and full depth over the
regular base.  Cohen--Macaulayness of an algebra without these hypotheses
does not by itself imply support saturation.

## Two applications

1. In the degree-forty-two Hessian residual problem, the synchronization
   defect vanishes away from an explicitly determined closed support.
   The remaining global target is the saturation of the residual ideal by
   the ideal defining that support.
2. In the cubic-normalization frontend, the closed-point cotangent torsion
   is already presented as
   \[
   (N:_F\mathfrak a^\infty)/N.
   \]

They are the algebra and module versions of the same theorem: a defect
known to vanish off a closed set extends across it exactly when the
ambient module has no torsion supported there.
