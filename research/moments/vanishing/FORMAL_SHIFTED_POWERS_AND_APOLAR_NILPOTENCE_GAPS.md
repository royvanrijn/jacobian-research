# Formal shifted powers, apolar nilpotence gaps and ordinary cutoff ideals

16 November 2026. Internal written reductions covering the ENTIRE formal
symbol question, in every dimension and polynomial phase degree. The general
eventuality-to-linearity implication remains OPEN. No independent replay,
formal verification, external review or novelty is claimed. See the
[checker](verify_formal_apolar_gaps.py),
[receipt](formal_apolar_gaps_verification.json) and
[source audit](formal_apolar_gaps_source_audit.json).

This route does not assume a holonomic source or an algebraic principal
symbol divisor. It replaces all polynomial markers by shifted phase powers,
then expresses the full remaining question as one asymptotic nilpotence-gap
problem in the ACTUAL apolar algebras of P^m. Ordinary cutoff ideals are
compared before taking radicals, including over nonreduced coefficient rings.

## 1. A formal-symbol commutator lemma over every Q-algebra

Let A be any commutative Q-algebra. A formal constant-coefficient operator
T=lambda(partial) acts locally finitely on A[x]. For a polynomial Q let M_Q
denote multiplication by Q and put delta=degQ. For any polynomial a,

    [T,M_a]=sum_(gamma!=0) M_(partial^gamma a)
                              (partial_Z^gamma lambda)(partial)/gamma!. (1)

The sum is finite because gamma differentiates a. More generally an operator
sum_i M_(a_i) f_i(partial), with polynomial a_i of degrees<=e and formal
f_i, has commutator with T of the same form with polynomial degrees<=e-1.
Thus

    (ad T)^(delta+1)(M_Q)=0.                                (2)

All these are identities of endomorphisms of A[x]; they require neither
finite differential order nor convergence. The binomial commutator identity
therefore gives, when m>=delta,

    T^m M_Q=sum_(j=0)^delta binomial(m,j)
                            (ad T)^j(M_Q) T^(m-j).           (3)

If T^(m-delta) f=0, each right-side term kills f. Consequently

    T^(m-delta) f=0  implies  T^m(Qf)=0.                    (4)

The finite-operator form of (4) is de Bondt's Theorem1 in
[1206.2836v2](https://arxiv.org/html/1206.2836v2), with the shifted-power
reduction in Corollary2. Its full proof and current version history were
rechecked. Equations (1)-(3) give the extension to arbitrary formal symbols
and arbitrary Q-algebras directly. No truncation order uniform in m is
assumed. This extension is not a new claim that the paper itself treated
formal symbols or nonreduced coefficient families.

## 2. All markers reduce to shifted powers, with their quantifiers retained

**GVCFORMALCOMM20261116.** For any polynomial P over A and any formal T,
the following qualitative conditions are equivalent:

* For every polynomial Q, T^m(QP^m)=0 for all sufficiently large m.
* For every fixed integer d>=0, T^k P^(k+d)=0 for all sufficiently large k.

The first implies the second by taking Q=P^d. Conversely, for a marker of
degree delta, take k=m-delta and d=delta. The second condition gives
T^(m-delta)P^m=0 for sufficiently large m; apply (4). On the pure locus the
d=0 condition already holds for every positive index, so it may be omitted.
The eventuality threshold still depends on d. This equivalence does NOT
replace all d by one shift or supply a linear bound by itself.

## 3. A single nilpotence-height sequence captures the full field question

Now work over a characteristic-zero FIELD K, with NONZERO polynomial P,
and assume T^mP^m=0 for every m>=1. If lambda(0)!=0, T is invertible on
each bounded-degree polynomial space (degree-triangular with nonzero scalar
diagonal), contradicting TP=0. Hence lambda(0)=0. T lowers degree by at
least1, unless T=0, which also causes no difficulty.

Define

    h_m=max{k>=0:T^k P^m!=0},    b_m=m-h_m,      m>=1.        (5)

The maximum exists and includes k=0, since P^m!=0 over a field.
Purity is exactly h_m<=m-1, so 1<=b_m<=m. The nilpotence INDEX is h_m+1;
the distinction of1 matters in the strict inequalities below.

Equivalently set

    A_m=K[[Z]]/Ann(P^m),
    Ann(P^m)={f:f(partial)P^m=0}.                            (6)

This is the actual finite local apolar algebra of P^m, with maximal-ideal
power (Z)^(Dm+1)=0 for D=degP. Then h_m is the largest k with the class of
lambda^k nonzero in A_m. Each h_m is determined by the finite symbol jet
through order Dm; no joint-source presentation is used. Formal units preserve
the entire sequence, since (U lambda)(partial)^k=U(partial)^k T^k and the
unit action on polynomials is invertible.

**GVCAPOLARGAP20261116.** Under these field and purity hypotheses:

    individual-marker eventuality  iff  b_m tends to infinity;            (7)
    one common linear cutoff       iff  inf_(m>=1) b_m/m >0
                                   iff  liminf_(m->infty) b_m/m >0.       (8)

For (7), shifted eventuality at shift d says T^(m-d)P^m=0 for all large m,
which says precisely b_m>d. Apply section2. Conversely, whenever b_m>delta,
equation (4) kills EVERY marker of degree<=delta at that index m.

For a quantitative forward implication in (8), put

    epsilon=inf_(m>=1) b_m/m.

If epsilon>0, then m>degQ/epsilon implies b_m>degQ; (4) gives a sufficient
integer cutoff ceil(1/epsilon). Conversely, if a common cutoff C>=0 is valid
and D>=1 bounds degP, then T^(h_m)P^m!=0. For h_m>0 use the marker
Q=P^(m-h_m) at moment index h_m. Its degree is at most D b_m. Therefore

    h_m<=CD b_m,    b_m/m>=1/(CD+1)  for EVERY m>=1.          (9)

For h_m=0 these inequalities hold directly. This proves necessity.
Finally b_m>=1 makes a positive asymptotic lower bound equivalent to a
positive infimum over all m: the finitely many early ratios are positive.

Thus the ORIGINAL OP-GVC-FORMAL-LINEAR is exactly the following question,
without any restriction to a smaller family:

    For the actual power-apolar algebras (6) of a fixed polynomial P,
    does b_m -> infinity force liminf b_m/m >0?              (10)

A counterexample would require an actual formal lambda and polynomial P
with unbounded b_m but liminf b_m/m=0. An arbitrary integer sequence such as
ceil(sqrt(m)) is NOT a realization and does not disprove (10).
No assertion that all possible height sequences have rational limits,
eventual linear formulas or sub/superadditivity is made.

## 4. Ordinary cutoff ideals have an excess-power presentation up to cofinality

Fix n and a phase degree bound D>=1. In the universal coefficient ring

    R=Q[p_beta (|beta|<=D), l_gamma (gamma in N^n)]

take the universal polynomial P and formal symbol lambda. This polynomial
coefficient ring has infinitely many symbol variables. Every individual
output still uses only finitely many of them. No infinite-ring ACC is used.

For a real c>=0 define I_c to be the ordinary ideal generated by all output
coefficients of T^m(x^alpha P^m) for m>c|alpha|, m>=1. For integer N>=0 let

    E_N=(coefficients of T^k P^(k+d): d>=0,k>=1,k>Nd).       (11)

Both include the entire ordinary pure ideal through d=0 or alpha=0.

**GVCEXCESSIDEALS20261116.** The literal ordinary containments are

    I_(N+1) subset E_N;       E_(CD) subset I_C              (12)

for integers N,C>=0. In particular

    I_(CD+1) subset E_(CD) subset I_C.                      (13)

For the first containment work in R/E_N. If delta=|alpha| and
m>(N+1)delta, let k=m-delta. Then k>Ndelta and T^kP^(k+delta)=0 in that
quotient. Formula (4), valid over this arbitrary Q-algebra, kills the
corresponding marked output. The delta=0 case is already a pure generator.
This proves ordinary membership, with all nilpotents retained.

For the second, expand P^d as a polynomial marker of degree<=Dd. If k>CDd,
every monomial term of that marker satisfies k>C|alpha|, so all its marked
output coefficients belong to I_C. This proves (12)-(13) before radicalizing.
The same argument works in any specialized finite or infinite coefficient
ring. It does not assert equality at the same cutoff or finite generation.

Consequently the two decreasing ordinary ideal filtrations are cofinal.
Their quotient inverse limits are canonically isomorphic, and on arbitrary
Q-algebras their unions of prescribed-cutoff functors coincide. This is a
description of the LINEAR-cutoff class. Section3's gap remains necessary to
identify it with separate individual eventuality; cofinality alone does not
close that gap, force stabilization or provide finite universal constants.

## 5. What the recent symmetric-powers paper does and does not supply

Flavi--Jelisiejew--Michalek,
[2408.02754v2](https://arxiv.org/html/2408.02754v2), revised29August2025
and published in IMRN2025 as rnaf277, Section3, compares power-apolar
algebras with symmetric tensor constructions. The relevant diagram has a
subalgebra of the symmetric tensor algebra mapping onto Ap(P^m); it is not
an unconditional identification. Their maximum-dimension corollary makes
the maps isomorphisms under its stated dimension equality. Their example
P=x^2 versus P=x^2+y has isomorphic first apolar algebras but square-apolar
dimensions5 versus6. Thus the abstract algebra Ap(P) alone does not
determine the family needed in (10). No smoothability, rank or dimension
theorem from this paper is imported to infer the missing nilpotence gap.

The checker independently confirms those two small dimensions as a source
transcription control. For future use, preserve the embedding, distinguished
generators and the actual quotient in (6); replacing the family by symmetric
powers of one abstract algebra without checking the extra hypothesis could
change the nilpotence sequence. No conclusion about (10) follows just from
a polynomial dimension bound.

The positive-weight balanced filtered reduction considered during this turn
was already proved in
[GVCFILTEREDPAIR/GVCFILTEREDCUT](FILTERED_PRINCIPAL_PAIR_GVC_TRANSFER.md).
It is retained rather than registered again. The present reduction (7)-(13)
has no balancing-weight hypothesis and covers the entire unresolved class.

## 6. Verification and next decisive gate

Exact controls check formal-jet commutators, binomial normal ordering,
local nilpotence heights, the strict shift/index conversions, a square-zero
coefficient specialization, and the stated small apolar dimensions. These
finite controls do not prove (10), realize arbitrary abstract sequences,
or certify a general ordinary-ideal solver. The universal formal and ring
arguments are written above; no quantitative asymptotic fitting is used.

The decisive gate is now the gap alternative (10) for the actual family (6).
The already closed binary, quadratic, algebraic-divisor and balanced classes
must not be mistaken for a general proof. The raw holonomic-source route
has the exact algebraicity boundary proved in
[FORMAL_SYMBOL_SOURCE_RIGIDITY.md](FORMAL_SYMBOL_SOURCE_RIGIDITY.md).
Both full exceptional Haar MZ cases remain open at their actual complement
zero-value logarithmic extension. Closure-first priority, truthful assurance
flags and deferred external review/outreach remain unchanged.
