# One eventual cutoff for every ordinary binomial: finite-regime synthesis

SCRATCH theorem candidate. The three new boundary packets cited below are
whole-proof-reviewed internal drafts; their source guards are being built.
This note does not change the canonical registry or claim external review.
All existence constants are unevaluated.

For integers a>=0,d>=1,m>=1 put

    A_(m;a,d)(x)=sum_(j=0)^m binom(m,j)(am+1)_(dj)x^(m-j).

Call (a,d,N) bad when N>=2 and the actual rational polynomials
A_(N;a,d), A_(N-1;a,d) have a nonconstant common divisor. Equivalently
they have a common complex zero. The polynomials are monic and have
nonzero constant terms, so all common roots are nonzero and finite.

## 1. Precise inputs and their scopes

(G0) Adjacent exponent d=1: the retained raw theorem thm:main in
UNIFORM_PADE_PARABOLA_SFC.tex excludes bad pairs for every a>=0,N>=2.

(G1) Constant lower exponent a=0: the raw prime-power and uniform
analytic/large-exponent arguments give one uniform cutoff N=10000
for ALL d. This is the alternative proof below, not a citation to the
large completed finite-prefix computation.

(G2) Uniform quadratic high-lower-exponent tail, retained401:
ALL a>=10^24*d^2,d>=2,N>=2 are good.

(G3) UNRESTRICTED_SMALL_RATIO_GAMMA_TWO_WELL_EXCLUSION.md:
there exist epsilon_s>0,N_s such that ALL a>=1,d>=2,c<=epsilon_s,
N>=N_s are good, where c=a/(d-1)^2. This packet imports the old raw
sparse theorem as a logical premise. It does not use the relative-gap
all-order theorem or the lower-one finite boundary certificate.

(G4) COMPACT_SUBCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md, canonical414:
for each fixed compact I inside (0,1/4), one N_-(I) works for all
supports with c in I. Its finite-gap fallback is the retained strict
subcritical theorem, not general407.

(G5) UNIFORM_BOUNDED_THERMAL_TRANSITION_EVENTUAL_SFC.md, canonical413:
for each fixed finite U, one N_b(U) excludes all supports with
sqrt(N)*|r|<=U, where

    s=(d-1)/2, e=a-s^2, C=a+s, r=e/C.

This includes exact r=0 via canonical412.

(G6) LARGE_THERMAL_QUARTIC_FOLD_MATCHING.md:
there is no bad sequence with N->infinity, 0<|r|->0 and
sqrt(N)*|r|->infinity. Its exact critical actions and original cycles
are retained. The analytic model is not substituted at finite support.

(G7) COMPACT_SUPERCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md:
for each fixed compact I inside (1/4,infinity), one N_+(I) works for
all supports with c in I. The large-gap core has a narrow classical
Kabluchko2.3 inverse-existence input. The full statement uses407 only
for its finitely many bounded-gap supports.

The new main mathematical work is in G3,G6,G7. The synthesis below is
an elementary finite covering and has no independent analytic estimate.

## 2. Uniform critical strip from the sequential matching theorem

G6 implies there exist delta>0, U>=1 and N_l such that every support
with N>=N_l,

    0<|r|<delta, sqrt(N)*|r|>U

is good. If no such triple of constants existed, for each integer j
choose a bad triple with N>=j,0<|r|<1/j and sqrt(N)*|r|>j.
That sequence contradicts G6. Increasing U to at least1 is harmless.
Together with G5 at THIS ONE fixed U, it follows that

    N>=N_c:=max(N_l,N_b(U)), |r|<delta

implies good, including r=0. No assertion uniform in unbounded U is
inferred from G5 alone.

The exact ratio relation is

    r=(c-1/4)/(c+1/[2(d-1)]).

Choose h=min(1/8,delta/16)>0. If |c-1/4|<=h then c>=1/8 and
|r|<=8h<=delta/2<delta. The critical-strip cutoff N_c therefore
covers the entire CLOSED c interval [1/4-h,1/4+h], uniformly even
when d is bounded. There is no unproved interchange of limiting
support and limiting order in this step.

## 3. A finite covering of every positive-lower-exponent support

Set epsilon=min(epsilon_s,1/16) and K=4*10^24. Then epsilon<1/4-h,
and K>1/4+h. Five fixed closed regions cover all c>0:

    (0,epsilon],
    [epsilon,1/4-h],
    [1/4-h,1/4+h],
    [1/4+h,K],
    [K,infinity).

Their respective inputs are G3, G4, Section2, G7 and G2. The last
use is exact: for every integer d>=2, (d-1)^2>=d^2/4, so

    c>=K => a=c(d-1)^2>=10^24*d^2.

Let N_* be the maximum of2,10000 and the four finite eventual cutoffs
of the first four regions. This is ONE integer, independent of a,d.
G0 and G1 cover d=1 and a=0. Thus:

**Theorem SupportUniformEventualOrdinaryBinomialSFC.** There exists
one unevaluated integer N_* such that for ALL integers a>=0,d>=1
and ALL N>=N_*, the ACTUAL pair A_(N;a,d),A_(N-1;a,d) is coprime over Q.

A rational Bezout identity remains a Bezout identity after extension to
any characteristic-zero field. For B!=0, exactly

    Lambda((X^a(A+B X^d))^m)=B^m(am)! A_(m;a,d)(A/B).

A zero coefficient leaves a nonzero monomial with nonzero factorial
moments over a characteristic-zero field. Consequently this SAME
N_* works for every nonzero ordinary binomial coefficient pair over
every characteristic-zero field. The converse from the complex-field
formulation is literal: a common root z gives the phase X^a(z+X^d).

## 4. Old raw fixed-terminal support finiteness

Fix N>=2. For a>=1 let e_N be the least integer with
2^e_N>=4N(N-1), D_N=2(N-1)e_N. The raw retained thm:gapbound in
lower_exponent_initial_and_boundary.tex makes all d>=D_N good.
The raw Gamma soft criterion makes a>=d(2N-3) good. Thus every bad
positive support lies in the FINITE containing set

    1<=d<D_N, 1<=a<d(2N-3).

This also covers d=1 without requiring its all-order proof.
For a=0, the raw large-exponent simplicity theorem makes all gaps
with 2^d>=4N^2 good: the exact Appell identity
A_(N;0,d)'=N*A_(N-1;0,d) converts simplicity to coprimality.
Only finitely many smaller d remain. Therefore each fixed-terminal
bad-support fiber is finite using just these old raw criteria.
No huge boundary prefix or terminal-ten certificate is a premise.

## 5. Genuine reverse implications and exact limits

Let B be the FULL actual bad triple set, a>=0,d>=1,N>=2.
The following statements are equivalent:

(i) One integer N_* excludes all bad triples at every N>=N_*.
(ii) The set B is finite.
(iii) One terminal cutoff works for all nonzero ordinary-binomial
      coefficient pairs over every characteristic-zero field.
(iv) One terminal cutoff works for all such coefficient pairs over C.

(i)=>(ii) uses Section4 for the finitely many N<N_*; without that
fixed-terminal support input the implication is not formal.
(ii)=>(i) takes a larger integer than the maximum of the finite
N-projection (and takes2 if B is empty).
(i)=>(iii) is rational Bezout, (iii)=>(iv) is specialization, and
(iv)=>(i) uses a complex common root as in Section3.
All four predicates are meaningful independently of the analytic
route and can be represented as a genuinely equivalent DAG component.
The component must have an actual seeded proof route; no cycle can
seed itself. No analytic tail estimate follows in reverse from (i).

**Corollary FiniteProjectiveExceptionalPhases.** Up to nonzero scalar
multiplication of the whole phase, only finitely many ordinary complex
binomial phases can have ANY consecutive zero factorial moments.
Indeed B is finite, and each bad triple contributes at most N-1
projective coefficient ratios, the finitely many roots of its gcd.
Every such ratio is algebraic over Q. Supports and coefficient ratios
are both part of the finite data. This does not say there are finitely
many unscaled coefficient pairs: multiplying a phase by an arbitrary
nonzero scalar preserves its vanishing moments.

There is also a genuine converse, with its extra retained premise explicit.
Let E be the set of projective complex phases that have at least one
consecutive zero pair, and let S be the support-pair projection of B.
Then B finite implies E finite implies S finite: every bad support
contributes a phase by a complex gcd root. Conversely, if S is finite,
retained407 supplies for each support in S a terminal cutoff independent
of coefficients. Taking the maximum of these finitely many support
cutoffs bounds all bad terminal orders. There are then only finitely
many triples in B. Thus, USING407, the four equivalent statements above
are also equivalent to each of

(v) Only finitely many projective ordinary complex binomial phases
    have any consecutive zero pair.
(vi) Only finitely many ordinary support pairs (a,d) have any bad order.

The implication from (v) or (vi) to B finite uses fixed-support eventual
coprimality; finiteness of the phase/support projection alone does not
bound orders for an arbitrary abstract bad set. The proof of (i)=>(ii)
in Section4 does NOT require407. These two reverse mechanisms must
remain distinct in the DAG and in source-withdrawal checks.

The finite global exceptional set has NOT been shown empty. N_* is
NOT evaluated. The first-two-moment assertion at every N, the full
all-order ordinary binomial Strong Factorial statement, and every
larger polynomial-support conjecture remain open. No replication of
a bad triple into infinitely many ordinary-binomial bad triples is
supplied by the theorem. A finite computation without an explicit
proved containing box cannot settle the full conjecture from this
non-effective finiteness result alone.

## 6. Minimal-input constant-family uniform tail

For completeness, this is G1 without the completed finite prefix.
At a=0 the polynomial is A_(N;0,d)=sum_j binom(N,j)(dj)! x^(N-j).
For d=1 use G0; for d=2,3,4,5 use the raw prime-power simplicity
proof PRIME_POWER_TWO_TERM_SFC.tex, thm:sfc at shape1. For 6<=d<=sqrt(N),N>=10000, raw thm:analytic in
ALL_UNIVARIATE_TWO_TERM_SFC.tex gives squarefreeness. For d>sqrt(N),
raw thm:large in ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex gives
squarefreeness because 2^d>2^sqrt(N)>=4N^2. To prove the final
universal inequality put t=sqrt(N)>=100. At t=100 the integer
comparison 2^100>4*100^4 holds. The derivative of
log(2^t/(4t^4)) is log2-4/t>0 since log2=integral_1^2 du/u>=1/2.
The exact Appell identity then gives consecutive coprimality.
The analytic source reproduces the needed exponent-six entropy bound
and endpoint calculation, with a credit link to the RAW
EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex proof. That raw proof can be
authenticated for the analytic lemma; no finite-prefix computation is
imported. This uses just the named raw analytic/large-gap/prime-power
proofs, not the92935 or9998 finite-prefix computations.

## 7. Provenance and publication scope

This is a scratch internal synthesis, conditional on retaining and
accepting the exact whole analytic packets G3,G6,G7 and old inputs
G0--G2,G4,G5. It is not yet a canonical theorem declaration. All
source guards, portable replays, explicit source-role dependencies,
whole mathematical review and scope-withdrawal DAG tests must be
completed before promotion. A source hash confirms identity, not a
correct analytic estimate; a finite literal gcd does not establish
any displayed universal theorem. No independent external or formal
verification is claimed. No historical-priority claim is made here.
