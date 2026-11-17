# Paired operator powers and polynomial phase recomposition

16 November 2026. Internal all-degree preservation and reflection theorem.
No independent replay, external review, formal verification or novelty is
claimed. See the [checker](verify_gvc_paired_recomposition.py),
[receipt](gvc_paired_recomposition_verification.json) and
[source audit](gvc_paired_recomposition_source_audit.json).

For monic H of degree e>=1, qualitative all-marker safety is EQUIVALENT
for (T,P) and (T^e,H(P)), for every formal constant-coefficient operator T
with ZERO constant symbol coefficient, over every commutative Q-algebra.
On the original pure locus, linear cutoffs transfer both ways with explicit
coefficients and ordinary-ideal containments. Over fields the asymptotic
nilpotence-height slope is exactly unchanged. The general formal-linearity
question remains OPEN. This theorem differs from fixed-functional
Gaussian/Haar recomposition: the operator is simultaneously raised to e.

## 1. Assumptions and the finite power basis

Let A be any commutative Q-algebra, P in A[x] of degree<=D with D>=1,
and T=lambda(partial), lambda in A[[Z]], lambda(0)=0. T lowers polynomial
degree by at least1. Zero operators/phases, zero divisors and nilpotent
coefficients are allowed. Retain polynomial outputs without specializing x.
For monic H(t) in A[t] of degree e>=1, put

    R=H(P), S=T^e.

The polynomials t^i H(t)^k, 0<=i<e,k>=0, form an A-module basis of A[t]:
their leading degrees ek+i are distinct and their leading coefficients1.
Triangular reduction expresses t^N as a finite combination with ek+i<=N.
Substitution t=P remains valid even when P has relations; independence
after substitution is not used. A unit leading coefficient can be normalized
to1, since unit phase scaling changes no zero predicate or coefficient
ideal. Reflection is not asserted for a nonunit leading coefficient.

## 2. Qualitative preservation, reflection and bad markers

**GVCRECOMPQUAL20261116.** Without assuming purity,

    every Q has T^m(QP^m)=0 eventually
       iff every Q has S^m(QR^m)=0 eventually.              (1)

Forward: expand H(P)^m into P^j, 0<=j<=em. Fix Q. For all j above its
original eventuality threshold, T^j(QP^j)=0 and hence T^(em)(QP^j)=0.
The finitely many smaller j give fixed polynomials killed by all sufficiently
large powers T^(em). Thus every term vanishes.

Reverse: fix Q and take a common eventuality threshold for the FINITE list
of new markers QP^i, 0<=i<e. Expand P^N using the basis of section1.
For k above that threshold, with ek+i<=N,

    T^N(QP^i R^k)=T^(N-ek)[S^k(QP^i R^k)]=0.               (2)

For the finitely many smaller (i,k), local nilpotence kills those fixed
polynomials when N is large. This proves (1) over the stated rings.

If one original marker Q is nonzero at infinitely many indices, at least
one of Q,QP,...,QP^(e-1) is a bad new marker, of degree<=degQ+D(e-1).
Otherwise the reverse proof would make Q eventual. The same Q alone need
not remain bad, and no every-late-index assertion follows.

Original purity implies new purity: for 1<=j<=em,
T^(em)P^j=T^(em-j)(T^jP^j)=0; also T^(em)1=0. Expand H(P)^m.
No converse from new PURE tests alone is asserted; the reverse proof of
(1) uses markers. Original purity will be kept explicitly below.

## 3. Explicit all-ring cutoffs and ordinary ideals

Assume the ORIGINAL pair is pure. A cutoff C means T^m(QP^m)=0 when
m>CdegQ, including constant markers.

**GVCRECOMPCUT20261116.** If C>=0 is an integer original cutoff, then

    C_forward=ceil((CD+1)/e)                               (3)

is a sufficient new cutoff. For a marker of degree delta, terms QP^j
with j>Cdelta vanish as above. The remaining terms have input degree
<=delta+Dj<=(CD+1)delta and are killed when em>(CD+1)delta. Constant
markers use the pure identity. No optimality is claimed.

Conversely, if C>=0 is an integer new cutoff, then, under ORIGINAL purity,

    C_backward=(1+eDC)(1+D(e-1))                            (4)

is sufficient for the original pair. For delta>=1 use the phase-power
basis. If k>C(delta+Di), the bracket in (2) vanishes. Otherwise

    deg(QP^i R^k)<=delta+Di+eDk
                 <=(1+eDC)(delta+Di)
                 <=C_backward delta,

so T^N kills the term for N>C_backward delta. Constants use original
purity. Thus SOME linear cutoff exists for one pair iff it exists for
the other, on the original pure locus, over EVERY Q-algebra.

In the universal coefficient ring with lambda(0)=0 and freely adjoined
monic-H coefficients, let I_C denote the full ORDINARY cutoff ideal and
J the full ORDINARY pure ideal. Apply these ring theorems in the actual
quotient coefficient rings to obtain literal containments

    J(S,R) subset J(T,P),
    I_(C_forward)(S,R) subset I_C(T,P),
    I_(C_backward)(T,P) subset J(T,P)+I_C(S,R).              (5)

These retain nilpotents and embedded structure. The reverse sum with
J(T,P) is compulsory. Each coefficient uses finitely many jets; no
infinite-ring ACC, radical replacement or same-index mixed-ideal equality
is used. General ordinary-ideal producers remain unimplemented.

For e=1,H(t)=t+c, EVERY finite-prefix ordinary PURE ideal is exactly
unchanged. In T^m(P+c)^m, the j=m term is T^mP^m and each 1<=j<m term
is a scalar multiple of T^(m-j)(T^jP^j); j=0 is zero. Triangular induction
gives both inclusions prefix by prefix. No least mixed-cutoff equality is
claimed. Unit affine scaling preserves the same conclusion.

The constant-symbol hypothesis is essential: T=Id,P=0 is safe, whereas
H(t)=t+1 gives R=1,S=Id and fails safety. The unit-leading hypothesis
also matters. Over A=K x K take H(t)=a t with a=(1,0) and put the retained
[DVG1 pure failure](../../foundations/extended-geometry/DVORSKY_GVC5_COUNTEREXAMPLE.md)
in the second coefficient factor (the first phase can be0).
Recomposition erases it, so reflection fails for such a nonunit leading
coefficient. This is a coefficient-ring boundary, not a new frontier.

## 4. Exact asymptotic height invariance

Now K is a characteristic-zero FIELD and P is NONCONSTANT. Continue
lambda(0)=0; purity is not needed until the last consequence. Define

    h_m=max{k>=0:T^kP^m!=0}, h_0=0,
    h'_m=max{k>=0:T^kR^m!=0}, h'_0=0,
    g_m=max{k>=0:S^kR^m!=0}=floor(h'_m/e),
    M_N=max_(0<=j<=N)h_j.

All maxima exist by degree lowering, and their phase powers are nonzero.
There are exact inequalities

    h'_m<=M_(em),
    M_N<=max_(0<=k<=floor(N/e))h'_k+D(e-1).                (6)

The first uses expansion into P^j. For the second, multiplication by
P^i increases T-height by at most Di, by the retained
[formal commutator lemma](FORMAL_SHIFTED_POWERS_AND_APOLAR_NILPOTENCE_GAPS.md):
T^(h'_k+1)R^k=0 implies T^(h'_k+1+Di)(P^iR^k)=0.
Apply this to the basis expansion of each P^j, j<=N.

For any nonnegative at-most-linear sequence,

    limsup M_N/N=limsup h_N/N.                              (7)

The lower inequality is immediate. For the upper, bound h_j by
(rho+epsilon)j for all sufficiently large j, and absorb the finitely many
earlier values into a constant. This proves (7) directly, without fitting
or a rational-slope hypothesis. The same argument applies to h'_k.

**GVCRECOMPHEIGHT20261116.** Equations (6)-(7) imply

    rho(T,P):=limsup h_m/m
              =limsup g_m/m=:rho(T^e,H(P)).                 (8)

One direction follows by dividing the first inequality in (6) by em.
For the other, divide the second by N and bound its running maximum by
the eventual linear bound for h'_k/(ek), plus a fixed early constant.
Finally |g_m-h'_m/e|<1, so the paired-power floor does not change limsup.

On the ORIGINAL pure locus the new pair is pure too. Hence

    liminf (m-h_m)/m=1-rho(T,P)                             (9)

is exactly preserved. Qualitative eventuality is preserved by (1), and
linear-cutoff existence by (3)-(4). Thus a pure eventual-but-nonlinear
pair stays such after EVERY paired recomposition. Conversely, new linear
safety would force original linear safety under the retained original
purity hypothesis. This does not reflect original purity from new tests.

For e=1 the spans of 1,P,...,P^N and of 1,H(P),...,H(P)^N are exactly
equal, so their running maximum heights agree for EVERY N. Since
sup_N M_N/N=sup_m h_m/m, even inf_m(m-h_m)/m is exactly affine-phase
invariant. Individual h_m need not be identified by this argument.

## 5. Scope, verification and the remaining gap

This closes all-degree qualitative preservation/reflection, finite bad-marker
block transfer, linear-cutoff transfer, ordinary-ideal comparisons and exact
asymptotic slope invariance for paired recompositions. The ring and field
scopes are separated above. Neither finite symbol order, algebraic symbol
divisors, balancing weights nor factorial moment laws are required.

The general gap question remains OPEN. Recomposition preserves its actual
obstruction and cannot supply the missing implication. No arbitrary abstract
height sequence is a realized counterexample. Source discovery this turn
did not yield a general gap theorem. Formal-local Noetherianity, bounded
height changes and polynomial dimension estimates are not substituted for
one; unrelated search snippets are not dependencies.

Exact controls check the finite phase-power basis, height comparisons,
paired-power rounding, cutoff arithmetic, forward annihilation and ordinary
pure-prefix identities, including a square-zero coefficient model. The
universal arguments are written proofs. No general source, ideal or gap
solver is run. The already closed binary/quadratic/algebraic-divisor/balanced
classes and all protected frontiers stay retained. Both full exceptional
Haar MZ cases remain open at actual complement logarithms. Closure-first
priority and deferred external review/outreach remain unchanged.
