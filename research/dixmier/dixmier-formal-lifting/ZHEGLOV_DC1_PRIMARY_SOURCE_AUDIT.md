# Zheglov's DC1 claim: a scoped primary-source audit

The claimed DC1 proof in Alexander Zheglov's
[*The Conjecture of Dixmier for the first Weyl algebra is true*,
arXiv:2410.06959v5](https://arxiv.org/html/2410.06959v5) is **not adopted
as a proved dependency** here. This audit independently checks a local
polynomial differential-equation lemma and a conditional Laurent
coefficient-invariance step. The earlier constructions supplying the
latter's hypotheses still require a complete proof review. This establishes
neither a proof nor a refutation of the paper's main theorem.

## Exact sources and the scope of inspection

The inspected Zheglov version is v5, dated 19 January 2026, 78 pages.
The [cached exact PDF](sources/zheglov-2410.06959v5.pdf) has SHA256
`897913ef95e011332e86e348803c6c9da5dbd196db17a04f8f951110e9b095e3`.
Both the primary HTML and the corresponding PDF passages were read.
Printed PDF page numbers are used below.

A later primary assessment is Kanel-Belov--Fedorov,
[*The generalized Dixmier conjecture and the two-dimensional Jacobian
conjecture*, arXiv:2610.33919v1](https://arxiv.org/html/2610.33919v1),
24 October 2026, Section 1.1. They include DC1 among their remaining
questions and write: “We have not been able to fully understand the proof.”
This reports their assessment of Zheglov's article; it is not a
mathematical disproof or an external review of the present audit.

The exact version, locations and adoption boundary are pinned in the
[source audit record](zheglov_dc1_source_audit.json).

## What the claimed proof needs

Theorems 4.1--4.2 concern subrectangular counterexample pairs and a finite
head-chopping sequence. They assert termination, alignment of two orders,
an exact total distance, and preservation of the rectangular Newton shape.
Lemma 4.3 supplies a polynomial-ODE restriction inside that argument.

The final comparison passes through the tail-cutting construction of
Theorem 5.1 and its correspondence with rank-one commuting differential
operators. In Section 6.2, Lemma 6.3 asserts independence of selected
Puiseux coefficients from an auxiliary integer \(N\); Steps 4--7 calculate
one of those coefficients as nonconstant in \(N\). The claimed contradiction
depends on both assertions referring to the same normalized coefficient.

Reading those statements and the final formulas is not a complete audit
of their dependencies. In particular, this audit has not independently
proved the all-stage Newton-shape assertions or the passage from the
generalized normal form to ordinary differential operators after cutting
the tail. The scalar Laurent invariance step is now independently proved
under explicit head-order, final-order, trace and normalization hypotheses;
verifying that the preceding constructions supply every hypothesis
remains outstanding. The paper's main theorem is not imported.

## An independently proved strengthening of Lemma 4.3

Let \(k\) have characteristic zero, let \(g,H\in k[x]\), and let
\(c\in k^\times\), \(A,d\geq1\) be integers. Suppose \(g\) has degree
\(\ell\) and \(s\geq2\) distinct roots over an algebraic closure, and
\[
           c g^A=gH'-\beta g'H,
                 \qquad \beta=A-1+\frac1d.                    \tag{1}
\]
Then
\[
 \begin{split}
  &g^{A-1}\operatorname{rad}(g)\mid H,\\
  &\deg H=\left(A-1+\frac1d\right)\ell,\\
  &\frac\ell d\in\mathbb Z,\qquad \frac\ell d\geq s.
 \end{split}                                                   \tag{2}
\]
Divisibility is understood up to nonzero constants. In particular
\(\ell/d>1\), the necessary conclusion of Zheglov's Lemma 4.3.
No Schur operator, DC-pair, subrectangularity assumption or formal
completion is needed for (2).

Extend scalars to an algebraic closure; polynomial identities and degrees
are unchanged. At a root \(\alpha\) of \(g\), write
\(j=\operatorname{ord}_\alpha g\) and
\(q=\operatorname{ord}_\alpha H\). The polynomial \(H\) is nonzero by
(1). If \(q\leq(A-1)j\), then \(q<\beta j\), so the leading terms in
the right side of (1) do not cancel. Its valuation is \(j+q-1<Aj\),
whereas the left side has valuation \(Aj\), a contradiction. Hence
\[
                           q\geq(A-1)j+1.
\]
Summing over the distinct roots proves the divisibility and
\(\deg H\geq(A-1)\ell+s\).

Put \(m=\deg H\). If \(m\ne\beta\ell\), the top terms at infinity
do not cancel either, and the right side of (1) has degree
\(m+\ell-1\). Equality with the left side would give
\(m=(A-1)\ell+1\), contradicting \(s\geq2\). Thus
\(m=\beta\ell\). Its integrality implies \(d\mid\ell\), and the
root lower bound gives \(\ell/d\geq s\). This proves (2), including
the possible resonant cancellations at roots and at infinity.

The root-count bound is sharp as a necessary condition. For any
\(d\geq2\), take
\[
 g=x^{d+1}(x-1)^{d-1},\qquad
 H=d\,g^{A-1}x(x-1),\qquad c=1.
\]
Substitution in (1) gives equality. Here \(\ell/d=s=2\).
For instance \(d=A=2\) gives
\(g=x^3(x-1)\) and \(H=2x^4(x-1)^2\).
The lemma is only a necessary condition for other \(g\); it does
not assert existence whenever its degree conditions hold.

## The one-root branch and its remaining normalization input

The [complete one-root ODE calculation](ONE_ROOT_CHOPPING_ODE_HOMOGENEOUS_BRANCH.md)
shows an extra homogeneous polynomial when \(l\mid d\). The displayed
ODE inference in Theorem 4.2, Case 2b on printed page 46 and induction
case b on page 49, retains only its particular power. A higher-degree
branch can instead be retained in the corner argument, so its existence
alone does not refute subrectangularity.

The exact two-stage local control also satisfies the next coefficient
chop, both ODEs, and the displayed numerical order, gcd and distance
conditions. Its terminal polynomial has several powers. In the case
\(d_2=2,l_2=1\), a common translation leaves a normalized subleading
coefficient difference \(3/(92C)\ne0\). Thus the simultaneous centering
used in Lemma 6.1 on printed page 64 does not follow from those local
conditions alone. This is precisely relevant to the exceptional
missing-head normalization required by the conditional Laurent lemma.

The control deliberately does not realize the initial simultaneous
Weyl pair and first chop. The concrete remaining obligation is to show
that the actual pair or its normal forms force the extra coefficient
\(C\) to vanish, or to repair the transport and terminal-edge argument
while retaining it. This is a missing local justification to audit,
not a counterexample to Theorem 4.2 under all its global hypotheses or
a refutation of the main DC1 claim.

## Repaired indexing and a verified conditional invariance step

The proof of Lemma 6.3, printed pages 70--71, first uses the order scale
\[
 \operatorname{ord}Q_{N,i}
       =\frac{d_2N+l_2}{d_2}\operatorname{ord}Q_i
                          \qquad(i<I),                        \tag{3}
\]
and then displays the coefficient identification
\[
                         c_{n,(d_2N+l_2)i}=c_i.
                                                               \tag{4}
\]
The preceding expansions use raw Laurent exponents. Put
\(D=d_2N+l_2\). The consistent lattice formula is
\[
                         c_{n,Dj}=c_{d_2j},               \tag{5}
\]
or \(c_{n,(D/d_2)i}=c_i\) on the relevant raw indices divisible by
\(d_2\). The thresholds must use the same convention. Thus (4)
mixes compressed and raw indexing as printed. This admits a notational
repair and does not establish a failure of the invariance conclusion.

The independent
[conditional Laurent-coefficient lemma](CONDITIONAL_LAURENT_COEFFICIENT_INVARIANCE.md)
bypasses that displayed formula. For a fixed chopping chain and its
fixed intermediate heads, a compressed formal head root \(W\) gives
\(W_N(t)=W(t^D)\). A derivative calculation identifies the distinguished
exponent and a nonlinear-order bound proves \(N\)-invariance when the
last head gap exceeds \(N+1\). Integrality proves that the only other
possibility is the exceptional gap \(N\); that case is handled using
the required missing coefficient one below the final head. The proof
applies to every prefix, and states all head, final-order, trace and
normalization assumptions explicitly. No separate gap assumption or
unresolved \(N+1\) boundary remains.

The accompanying written cyclic-vector trace, filtered tail and
Puiseux trace calculations verify the stated normalization implications
under their own explicit assumptions. They do not prove that the
source's earlier normal-form and commuting-curve constructions exist
with all those properties. The stronger wording about integer-polynomial
dependence is not adopted by this conditional argument.

The [exact finite controls](check_conditional_laurent_invariance.py)
include a nonzero invariant coefficient \(-231/512\) at three values
of \(N\), and resonant controls showing precisely why the missing-head
condition is needed. They are commutative Laurent controls, not
DC-pairs, and the all-\(N\) result follows from the written proof.

The retained local results are (2), the complete one-root ODE
classification and countercontrol, and the conditional invariance lemma. The earlier application inputs and the main DC1 claim remain
unverified. The final rational expression's dependence on \(N\), by
itself, cannot supply those inputs. No metadata is upgraded, and no
formal verification, independent whole-proof review, external review
or novelty priority is claimed.
