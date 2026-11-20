# Fiberwise-uniform terminal bounds for ordinary binomials

The fixed-support eventual theorem combines with retained all-order support tails to give stronger quantifiers. For every fixed gap there is one eventual terminal cutoff for ALL lower exponents. For every fixed lower exponent there is one eventual cutoff for ALL gaps. These cutoffs are uniform over coefficients and characteristic-zero fields, but are not evaluated and need not be uniform as both support parameters vary.

For integers a>=0,d>=1,m>=1 put

    A_(m;a,d)(x)
       =sum_(j=0)^m binom(m,j)*(am+1)_(dj)*x^(m-j).

Define the set of actual bad triples

    B={(a,d,N): a>=0,d>=1,N>=2 integers,
          gcd_(Q[x])(A_(N;a,d),A_(N-1;a,d))!=1}.          (F1)

For a fixed gap or lower exponent write

    B_d={(a,N):(a,d,N) in B},
    B^a={(d,N):(a,d,N) in B},
    T_(a,d)={N:(a,d,N) in B}.

The projections are BadLower_d={a:B_d has an element with lower exponent a} and BadGap_a={d:B^a has an element with gap d}. These are actual rational gcd statements, not prescribed-prime modular tests.

## 1. Exact retained inputs

The following quantified statements are the premises.

1. ALL_ORDINARY_BINOMIAL_EVENTUAL_STRONG_FACTORIAL.md, Theorem EveryFixedOrdinaryBinomialEventualStrongFactorial, proves T_(a,d) finite for EVERY fixed integer a>=0,d>=1. Its equivalent moment cutoff is uniform over all coefficients and all characteristic-zero fields at that support.

2. UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md, Theorem UniformQuadraticStrongFactorialTail, proves actual coprimality at EVERY terminal N>=2 whenever

       d>=2, a>=C*d^2, C=10^24.                       (F2)

   In particular BadLower_d is contained in {0,...,C*d^2-1}. The older qualitative FixedGapLargeLowerExponent theorem in FIXED_GAP_LARGE_LOWER_EXPONENT_EXTENSION.md also implies that this projection is finite. Either all-order large-a statement suffices for the finite-union argument; (F2) additionally retains an explicit bound on the possible lower exponents.

3. RELATIVE_GAP_THREE_STRONG_FACTORIAL.md, Theorem RelativeGapThreeStrongFactorial, proves actual coprimality at EVERY terminal whenever

       a>=2, d>=3a.                                  (F3)

   Thus BadGap_a is contained in {1,...,3a-1} for each a>=2.

4. The boundary families are known at every terminal. UNIFORM_PADE_PARABOLA_SFC.tex covers d=1 for every a>=0. ALL_UNIVARIATE_TWO_TERM_SFC.tex, Theorem thm:all and its exact Appell derivative identity, covers a=0 for every d. ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md covers a=1 for every d. These are separate already retained all-order inputs; the use of their full conclusions here is declared explicitly.

5. For the reverse interpretation of a terminal-only cutoff at fixed d, the retained ImprovedLinearLowerExponentCriterion in GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md gives

       a>=d*(N-2) => the actual pair at terminal N is coprime. (F4)

   Equality is included. The older sufficient criterion a>=d*(2N-3) in GAMMA_PRODUCT_ROOT_MODULUS.tex would also suffice to prove finiteness in that reverse direction, with a weaker finite bound. Neither soft criterion is needed for the forward finite-union proof in Sections2--3.

6. Theorem thm:gapbound in lower_exponent_initial_and_boundary.tex already gives an explicit bound at EACH fixed terminal. Let e_N be the least integer with2^(e_N)>=4N(N-1), and D_N=2(N-1)e_N. Then d>=D_N and a>=1 imply that all first N actual moment polynomials are pairwise coprime. This is the old exact high-gap theorem, also used by POLYNOMIAL_TRIANGLE_AND_TERMINAL_SIX.md; no new factorial asymptotic or rational-cutoff existence argument is needed.

All parameters a,d,N in this note are integers. The rising-factorial normalizations have nonzero leading and constant terms, so reciprocal polynomial conventions in the input sources preserve these actual gcd assertions.

## 2. One terminal cutoff at each fixed gap, uniform over every lower exponent

**Theorem FixedGapUniformTerminalCutoff.**
For each fixed d>=1 the entire set B_d is finite. Consequently there exists an integer L(d)>=2 such that

    for every a>=0 and every N>=L(d),
    gcd_Q(A_(N;a,d),A_(N-1;a,d))=1.                     (F5)

The same L(d) works simultaneously over every characteristic-zero field and every nonzero coefficient pair of every support with this gap.

Proof. For d=1, Input4 makes B_1 empty, so L(1)=2 is sufficient.

Fix d>=2. By Input2 a bad pair must have0<=a<C*d^2. There are only finitely many such lower exponents. For each one, Input1 says that T_(a,d) is finite. Therefore

    B_d=union_(0<=a<C*d^2) ({a} times T_(a,d))           (F6)

is a FINITE union of finite sets, hence finite. The equality uses (F2) to exclude every larger a. If B_d is empty take L(d)=2; otherwise take L(d)=1+max{N:(a,N) in B_d}. This proves(F5).

Equivalently, choose one support-uniform terminal bound from Input1 for each of the finitely many a<C*d^2 and take their maximum with2. Input2 covers every remaining a at every terminal. This explains the quantifier order: a SINGLE bound is chosen for the whole a-fiber before varying coefficients or fields. No effective procedure for this maximum is provided because the individual eventual cutoffs are unevaluated.

The exact equivalent statement over fields follows by the rational Bezout argument in Section5.

## 3. One terminal cutoff at each fixed lower exponent, uniform over every gap

**Theorem FixedLowerExponentUniformTerminalCutoff.**
For each fixed a>=0 the entire set B^a is finite. Consequently there is an integer U(a)>=2 such that

    for every d>=1 and every N>=U(a),
    gcd_Q(A_(N;a,d),A_(N-1;a,d))=1.                     (F7)

The same U(a) works simultaneously over every characteristic-zero field and all nonzero coefficients. One may take U(0)=U(1)=2 using Input4.

Proof. The boundary cases a=0,1 are the all-order assertions in Input4.

Fix a>=2. Input3 excludes EVERY gap d>=3a at EVERY terminal. Only the finitely many gaps1<=d<3a remain. Input1 makes T_(a,d) finite at each of them, so

    B^a=union_(1<=d<3a) ({d} times T_(a,d))              (F8)

is finite. Choose U(a)=2 if this set is empty and otherwise1+max of its terminal coordinates. Alternatively take the maximum of the finitely many support cutoffs from Input1 and2. This proves(F7). The choice depends only on a, and Input3 covers all other gaps uniformly at every order.

The core statement for a>=2 uses only Inputs1 and3. The extensions at a=0,1 use the separate all-order boundary theorems stated in Input4. They are not being inferred from fixed-support eventuality alone.

## 4. Genuine reverse equivalences and the necessary finiteness controls

**Lemma FiniteProjectionFiberEquivalence.**
For each fixed d,

    B_d finite
      iff [BadLower_d finite and every T_(a,d) finite]. (F9)

For each fixed a,

    B^a finite
      iff [BadGap_a finite and every T_(a,d) finite].   (F10)

Proof. A finite set has a finite projection and finite fibers. Conversely, a finite projection leaves only finitely many nonempty fibers, and the assumed finiteness of each of them makes their union finite. These are set-theoretic equivalences, not converses to one-way analytic or modular certificates.

**Corollary FixedGapFiniteBadPairsTerminalEquivalence.**
For each fixed gap d, the actual set B_d is finite if and only if there is a terminal bound L(d) satisfying(F5).

The forward direction takes the maximum terminal coordinate as in Section2. For the reverse, suppose(F5) holds. Any bad terminal has2<=N<L(d). By(F4), its lower exponent satisfies

    0<=a<d*(N-2).                                     (F11)

There are finitely many such integers for each of those finitely many N. Thus B_d is finite. For L(d)>=2 the resulting numerical-in-L(d) bound is

    |B_d|<=d*(L(d)-2)*(L(d)-3)/2,                       (F12)

with the value0 at L(d)=2,3; it is only a bound expressed in an UNEVALUATED cutoff, not a computed size of the bad set. Intersecting(F11) with a<C*d^2 is also valid. The older sufficient bound a<d*(2N-3) gives an independent weaker reverse estimate |B_d|<=d*(L(d)-2)^2, without using the sharper annulus input.

It would be false for an arbitrary abstract set that a uniform terminal cutoff alone implies the whole set finite: infinitely many lower exponents could have failed at one small terminal. The actual soft bounds in(F4), or their older version, exclude precisely that escape.

**Corollary FixedLowerFiniteBadPairsTerminalEquivalence.**
For each fixed a, the actual set B^a is finite if and only if there is a bound U(a) satisfying(F7).

Again the forward direction uses the maximum terminal coordinate. Conversely, at a>=2 Input3 confines every bad pair to d<3a, and(F7) confines it to N<U(a). This is a finite rectangle, giving

    |B^a|<=(3a-1)*(U(a)-2).                            (F13)

At a=0,1 the bad set is empty by Input4. The large-gap all-order premise is essential to this reverse argument; an abstract infinite set {(d,2):d>=1} would have a terminal-only cutoff3 while remaining infinite. This is a logical countermodel, not actual moment data.

## 5. Exact translation to field-valued moment windows

For a characteristic-zero field K and a nonzero phase f=X^a*(A+B X^d), the actual normalization is

    Lambda_K(f^m)=B^m*(am)!*A_(m;a,d)(A/B), B!=0.        (F14)

When B=0 the phase is a nonzero monomial and every positive moment is nonzero. When A=0 the same holds.

At any actual rational coprime pair, a Bezout identity over Q remains valid over K and can be evaluated at A/B. The scalar factors in(F14) are units. Therefore(F5) is equivalent to one large-terminal bound at the fixed gap working for every lower exponent, every characteristic-zero field and every coefficient pair. The analogous statement holds for(F7) at a fixed lower exponent.

Conversely a universal coefficient assertion over C already excludes a rational common divisor: choose a complex root x of that divisor and use the nonzero polynomial X^a*(x+X^d). Equivalently choose a root in an algebraic closure of Q and use that characteristic-zero field. A statement only over one specified non-algebraically-closed field is not substituted. No embedding of an arbitrary characteristic-zero field into C is assumed.

Thus the finite-bad-pair equivalences also have genuine uniform field-and-coefficient moment formulations, with the same L(d) or U(a).

## 6. The third coordinate: each fixed terminal has finitely many bad supports

**Theorem FixedTerminalFiniteBadSupports.**
For each fixed N>=2, the set

    B_N={(a,d):(a,d,N) in B}

is finite. More precisely define the exact integers

    e_N=min{e>=0:2^e>=4N(N-1)}, D_N=2(N-1)e_N.         (F15)

Input4 excludes a=0 and d=1. Input6 excludes d>=D_N whenever a>=1. Input5 excludes a>=d(N-2). Hence every bad support lies in the explicit finite triangle

    2<=d<D_N, 1<=a<d(N-2).                            (F16)

For N=2 this triangle is EMPTY, as the lower-exponent inequality is impossible. For N>=3 its number of support pairs is

    (N-2)*[(D_N-1)*D_N/2-1]-(D_N-2).                  (F17)

Indeed for each gap d there are exactly d(N-2)-1 positive lower exponents in(F16); sum over d=2,...,D_N-1. Thus the maximal gap is O(N log N), the maximal lower exponent is O(N^2 log N), and the NUMBER of candidate supports is O(N^3 log^2 N). The triangle is a containing set; its members are not asserted to be bad.

Fixed-terminal finiteness itself is an antecedent already supplied by the retained polynomial-triangle theorem, even with its older weaker lower-exponent bound a<d(2N-3). Formula(F16) combines its old high-gap input with the newer sharp soft bound. No new terminal-coprimality scan is part of this note.

**Corollary AllCoordinateFiberEscape.**
If there are infinitely many DISTINCT actual bad triples, then in every infinite sequence of distinct bad triples all three coordinates a,d,N eventually exceed each prescribed bound.

Proof. For a fixed upper bound on a, the finitely many lower-exponent fibers B^0,...,B^A are finite by Section3; their union is finite. The same argument with the finitely many fixed-gap fibers applies to bounded d. For bounded N, Section6 gives a finite union of finite terminal-support fibers. A sequence of distinct elements can visit each finite union only finitely many times. This proves simultaneous coordinatewise escape.

Conversely, for a countable triple set, that escape property for EVERY infinite sequence of distinct elements implies all fixed-coordinate fibers are finite: an infinite fiber could be enumerated while its fixed coordinate remained bounded. Thus this is a genuine set-theoretic equivalence with finiteness of all three types of coordinate fibers. It is not equivalent to finiteness of the whole triple set.

## 7. What these stronger quantifiers do not imply

We have proved BOTH

    for each d>=1 exists L(d), for ALL a>=0 and N>=L(d), coprimality;
    for each a>=0 exists U(a), for ALL d>=1 and N>=U(a), coprimality.

These statements do NOT imply one terminal bound working for all supports, nor finiteness of the global set B of triples. The abstract diagonal model

    B_model={(n,n,n):n>=M}, with any fixed M>=4              (F18)

has one bad pair in each fixed gap fiber and each fixed lower-exponent fiber, and at most one bad terminal for each fixed support. It also avoids the all-order regions a>=C*d^2 and d>=3a, and satisfies the sharp necessary condition a<d*(N-2). Its gap d=n is less than D_n because D_n>=2(n-1)>n for n>=4, so it lies inside the fixed-terminal containing triangle as well. Every fixed-terminal fiber has one element. Nevertheless it is infinite and its terminal coordinates are unbounded. The start M may be taken beyond any prescribed finite terminal prefix. This is a countermodel to a logical inference, not a claimed counterexample to Strong Factorial.

In particular the note supplies no effective L(d) or U(a), no finite executable global support box, and no proof that any finite bad set is empty. Full all-order ordinary binomial Strong Factorial remains a separate problem. The new content is a finite-union connection between independently retained quantified theorems, not a new saddle estimate or a newly proved all-order tail.


## Exact evidence and reuse boundary

    python3 verify_fiberwise_uniform_eventual_binomial.py --repo-root /path/to/repository --output review.json

The guarded replay authenticates the retained401,407 and sharp-soft core
packets and freshly runs their existing guards. It reuses the relative-gap
all-terminal theorem's already completed full finite source and receipt,
checking every recorded support and terminal-stream digest without
rerunning its millions of gcd computations. The existing all-order
boundary theorems are explicitly reused as mathematical inputs; their old
finite computations are not rerun by this note. The source guard directly
authenticates the two necessary a=0 prefix certificate files (92935
general-prefix records and9998 exponent-six records), their defining code
and selected current receipts, and the a=1 completed171812-case primary
bridge with its input, ledger, receipt and analytical source packet.
Historical top-level receipt entries superseded by current_replay are not
used. Authenticated older density context is not a mathematical premise.

Finite set models, exact triangle/rectangle counts, binary cutoffs and
actual moment normalizations support the written finite-union proof and
its scope. They are not a formal verification of the imported analytic
theorems or an inference from finite samples. All review remains internal;
no external review, historical priority or effective new cutoff is claimed.
