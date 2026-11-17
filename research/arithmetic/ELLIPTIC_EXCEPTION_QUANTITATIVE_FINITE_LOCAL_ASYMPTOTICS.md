# Quantitative finite-place counts and rational constants on the elliptic slice

For the existing elliptic exception, every fixed finite-place coefficient
is rational, and its logarithmic target count has a square-root
logarithmic error multiplied by a fixed power of the iterated logarithm.
The proof uses the explicit local tests and a lattice bound uniform in
the congruence index. It does not require a displayed Mordell--Weil basis.

## 1. Scope and notation

Use `E`, `Gamma=E(Q)`, the four-class squareclass map `alpha:Gamma->A`,
the target map, actual projective target height and coefficient `C_E`
from [EFA1--EFA2](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md).
For a fixed finite set `S` of places, write

\[
 H_S=\overline{\{(\alpha(P),(P)_p):P\in\Gamma\}}
       \subseteq A\times\prod_{p\in S_f}E(\mathbb Q_p),
 \qquad \mu_{S,a}=\operatorname{Haar}_{H_S}(U_{S,a}).             \tag{1.1}
\]

Haar measure is normalized on the whole joint closure. The count
`N_{S,a}(B)` concerns admissible full field-fibre targets of class `a`
soluble at every place in `S`, without requiring other places.

The real place and primes `p=1 mod3` impose no condition. Put

\[
 T=S_f\cap\bigl(\{2,3,5\}
                  \cup\{p:p\ge11,\ p=2\pmod3\}\bigr),
 \qquad r=|T|.                                      \tag{1.2}
\]

The local tests below are on the actual projected rational closures.
At the finitely many inadmissible points, extend their truth values in
any convenient way. This changes neither Haar measure nor the target
count beyond finitely many global exceptions. The same convention is
used for the finite reducible-cubic exceptions.

## 2. The joint closure is open

**Lemma EQA1 (an open product inside the joint closure).** For every
fixed finite `S`, there are integers `d_p>=1` such that

\[
 W=\{1_A\}\times\prod_{p\in S_f}E_{d_p}(p)
                           \subseteq H_S.           \tag{2.1}
\]

In particular `H_S` has finite index in the ambient compact product,
and all correlations between the factors are recorded in the finite
quotient `H_S/W`.

**Proof.** Fix the known rational nontorsion point `P_1`. The group
`E(Q_p)/E_1(p)` is finite at each of the finitely many primes in `S`.
Multiply `P_1` by a nonzero even integer `M` annihilating all these
finite quotients. Then `Q=MP_1` has `alpha(Q)=1` and lies in `E_1(p)`
at every such prime. It is nonzero at each prime, because a rational
nontorsion point remains of infinite order locally. Let its exact
formal depth at `p` be `d_p`.

Milne's filtration theorem gives `E_n(p)/E_{n+1}(p)=F_p`. His
short-model parameter estimate gives, for `t=-X/y`,

\[
 t(P+Q)=t(P)+t(Q)\pmod {p^{5n}}
                   \quad(P,Q\in E_n(p)).           \tag{2.2}
\]

Consequently multiplication by `p` raises every exact depth `n>=1`
by one: the linear term has valuation `n+1`, strictly below `5n`.
Intermediate prime-to-`p` multiples preserve depth. A point of depth
`d_p` therefore generates each cyclic quotient
`E_{d_p}(p)/E_N(p)`, of order `p^{N-d_p}`. Its closed integer span
is `E_{d_p}(p)`.

The integers are dense in the finite product of the groups `Z_p`, by
the Chinese remainder theorem. Applying this to multiples of the
single point `Q` proves that their closure contains the entire product
in (2.1). This is an open subgroup of the ambient product, so its
index there, and in `H_S`, is finite. For an empty finite-prime set,
the assertion just says that the trivial subgroup is open in `A`.
\(\square\)

The precise primary inputs are Milne, *Elliptic Curves* (2006),
Chapter II, Theorem 4.1 and Lemma 5.9, printed pp.62--63 and 67,
with Corollary 4.2 for prime-to-`p` multiplication:
[primary book](https://www.jmilne.org/math/Books/ectext6.pdf),
[retained copy](sources/milne-elliptic-curves-2006.pdf).
These statements apply to the integral short model at bad reduction
as well. No strong-approximation or generator-saturation assertion is
used in Lemma EQA1.

## 3. Explicit local annuli and rational Haar weights

Here a ball means a translate of a formal subgroup `E_n(p)`. An annulus
of exact depth `n` about a point is a translate of
`E_n(p) minus E_{n+1}(p)`, a disjoint union of `p-1` balls.

The previously proved local tests have the following descriptions.

- At two, the projected rational closure is
  `<T_0> times E_4(2)`. The translated component `T_0+E_4(2)` is
  soluble in full. In `E_4(2)`, precisely depths divisible by three
  are soluble, by ETF1--ETF2.
- At three, the closure is `E(Q_3)[2]+E_2(3)`. The `X=2` disc is
  soluble in full, the `X=0` disc is insoluble, and the soluble
  part of `X=1` is the local two-torsion translate of `E_3(3)`.
  In the identity disc, one requires depth divisible by three and
  normalized parameter equal to `1 or -1 mod9`, by ETC1 and ETC3a.
  At each such depth this is two of the six unit classes in
  `E_n(3)/E_{n+2}(3)`.
- At five, the complement of `E_1(5)` is soluble in full, and within
  `E_1(5)` precisely depths divisible by three are soluble, by ETF3.

The cited results are
[ETF1--ETF3](ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md) and
[ETC1--ETC3a](ELLIPTIC_EXCEPTION_COMPLETE_THREE_ADIC_CLASSIFICATION.md).
They are exact local polynomial tests; the finite global irreducibility
exceptions do not enter their local calculations.

For every remaining active prime `p>=11`, use
[ELV1](ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md): solubility is
equivalent to `v_p(F(x)) in3 Z`, where

\[
 F(x)=128x-25865,\qquad x_0=25865/128.
\]

There are at most two local zeros `B_+,B_-` above `x=x_0`. Whenever
they exist, their `y` coordinates are units. Indeed

\[
 f_E(x_0)=\frac{10585532408369625}{2097152},\qquad
 10585532408369625=3^2 5^3\cdot7\cdot739\cdot42649^2.
                                                               \tag{3.1}
\]

The primes `7,739,42649` are `1 mod3`, so none is an active prime
at least eleven. Existence of the two zeros is equivalent to `51730`
being a square over `Q_p`. Their reductions, if present, are distinct
smooth points because `y` and `2y` are units. In their residue discs
the coordinate `x-x_0` is an isometric local parameter. Equivalently
translation by `B_+` or `B_-` sends formal depth `n` to
`v_p(F)=n`: its leading coefficient is `128*2y(B_+)` or
`128*2y(B_-)`, a unit, and its higher coefficients are integral
because the center is a smooth integral point. At the identity,
formal depth `n` gives
`v_p(F)=-2n`.

Thus the only possibly insoluble residue discs are `E_1(p)` and
`B_++E_1(p),B_-+E_1(p)` when the latter exist. In each, precisely
depths divisible by three are soluble. All other residue classes
are soluble. If a local zero does not exist, an integral point cannot
have `F=0 modp`: its unit `y` coordinate would lift such a zero by
Hensel's lemma. This proves the description also at the bad-reduction
prime 1889. Smoothness of these particular residue points suffices;
good reduction of the whole elliptic curve is unnecessary.

It follows that, for each active prime and every sufficiently large
integer `N`, the local soluble set has inner and outer approximations
by disjoint unions of `O_p(N+1)` balls at varying subgroup depths.
Their difference is contained in `O_p(1)` balls of depth at least
`N-O_p(1)` about the relevant centers. At three a fixed two-level
angular refinement handles the unit condition. No refinement of every
piece to the common deepest subgroup is made.

The projected rational closure contains `E_{d_p}(p)` by EQA1. Hence
the normalized measure of any nonempty sufficiently deep ball in it
is `O_{E,p}(p^{-N})`. The inner and outer local measures therefore
differ by `O_{E,p}(p^{-N})`.

**Theorem EQA2 (rational finite-place coefficients).** Every coefficient
`mu_{S,a}` in EFA2 is a rational number. More precisely, choose `W`
as in EQA1 and write `J=[H_S:W]`. With representatives
`(a_j,(P_{j,p}))`,

\[
 \boxed{\mu_{S,a}=\frac1J
       \sum_{j:a_j=a}\ \prod_{p\in S_f} r_{j,p},
       \qquad r_{j,p}\in\mathbb Q\cap[0,1].}        \tag{3.2}
\]

Here `r_{j,p}` is the conditional soluble Haar measure in
`P_{j,p}+E_{d_p}(p)`. A certified description of this one finite
joint quotient determines (3.2) exactly by finite local arithmetic
and rational geometric series.

**Proof.** On a fixed ball, the descriptions above involve finitely
many initial subballs and eventual annuli whose allowed depths are
periodic modulo three. The only angular restriction is a fixed
finite residue condition. Haar masses of balls are rational, and
each infinite tail is a finite sum of geometric series with ratio
`p^{-3}`. Thus every `r_{j,p}` is rational. This remains true if a
branch center is outside the given ball; its contribution is then
constant or determined by a finite angular refinement.

Each coset of `W` in `H_S` has measure `1/J`. Conditional Haar measure
on that coset is the translated product Haar measure on `W`, which
is an actual product of formal groups. Multiplying the conditional
local measures and summing the allowed finite cosets proves (3.2).
The finite-quotient sum retains the correlations, including those
with the global squareclass; it is not an unconditional Euler product.
No bound on the work needed to certify the finite rational image for
an arbitrary `S` is supplied. \(\square\)

**Corollary EQA2a (finite-local dichotomy).** Fix `S` and `a`.
Then `mu_(S,a)=0` if and only if no admissible full field-fibre target
of class `a` passes the places in `S`. Equivalently, the count
`N_(S,a)(B)` is either identically zero or has a positive rational
leading coefficient and the quantitative asymptotic of EQA3.

**Proof.** Any such target has a rational elliptic preimage whose
image lies in the open set `U_(S,a)` of EFA1. That set is therefore
nonempty and has positive normalized Haar measure. Conversely a
positive coefficient gives infinitely many full field targets by
EFA2, or by EQA3 below. This proves both implications. Rationality
is EQA2. This dichotomy does not decide which finite-place sets are
feasible in the three classes without a certified all-places seed.
\(\square\)

## 4. A lattice estimate uniform in the congruence index

Fix any integral basis of the free rank-two group `Gamma/Gamma_tors`,
only for the proof, and its positive definite canonical quadratic
height form. For the ellipse of radius `R` in these fixed coordinates,
every coset of every sublattice `Lambda subset Z^2` satisfies

\[
 \#((u+\Lambda)\cap\mathcal E_R)
       =\frac{\operatorname{area}(\mathcal E_R)}{[\mathbb Z^2:\Lambda]}
                              +O_E(R+1),           \tag{4.1}
\]

with one constant independent of the sublattice, its index and the
coset.

Here is an elementary proof of the uniformity. Write a Hermite basis
`(a,0),(b,c)`, where `a,c>=1`. On a row `y=y_0+jc`, the allowed
`x` coordinates form one residue class modulo `a`. If `ell_R(y)`
is the horizontal length of the ellipse, the row contributes
`ell_R(y)/a+O(1)`. There are `O_E(R/c+1)` nonempty rows. The
compactly supported row-length function has total variation `O_E(R)`.
The elementary bounded-variation quadrature estimate gives

\[
 \sum_j\ell_R(y_0+jc)
       =\frac1c\int\ell_R(y)\,dy+O_E(R).
\]

Consequently the total error is
`O_E(R/a+R/c+1)`, bounded by `O_E(R+1)` as asserted. The horizontal
skew `b` affects only the residue on each row, not the bound.

For a coset of a finite-index subgroup of `Gamma`, the number of
torsion lifts of each free-lattice point is one or two. Its index
is the free-lattice index times the torsion quotient index. Thus
(4.1) gives the correct main fraction and still the same uniform
error. Finally

\[
 \log H(y(P))=3q_x(P)+O_E(1)
\]

from EEC1 brackets the actual height set between ellipses with
squared radii differing by `O_E(1)`. Their main areas differ by
`O_E(1)`, while both boundary errors remain `O_E(R+1)`.
We have therefore proved, uniformly in every subgroup coset `C`,

\[
 \#\{P\in C:H(y(P))\le B\}
   =\frac{2C_E}{[\Gamma:\Lambda_C]}\log B
                                  +O_E(\sqrt{\log B}+1).
                                                               \tag{4.2}
\]

The finitely many inadmissible points can be removed uniformly as
well. Existence of a basis for this argument does not assert that
the previously displayed points form one.

## 5. The quantitative finite-place asymptotic

**Theorem EQA3 (a quantitative form of EFA2).** For every fixed finite
`S` and every `a in A`, with `r` as in (1.2),

\[
 \boxed{N_{S,a}(B)=\mu_{S,a}C_E\log B
       +O_{E,S}\!\left(\sqrt{\log B}\,
                            (1+\log\log B)^r\right).}          \tag{5.1}
\]

This holds for sufficiently large `B`. For `r=0` the iterated
logarithm factor is one. The constant `mu_{S,a}` is the rational
number of EQA2; it may be zero for some classes.

**Proof.** Ignore the places with automatic local points, and use the
joint closure for `T` together with `A`. Put `L=log B` and `R=sqrt(L)`.
For every active prime choose

\[
 N_p=\lceil\log_p R\rceil+O_{E,p}(1),
\]

with fixed extra depths making the local descriptions of section 3
valid. Their inner and outer approximations have `O_p(N_p+1)`
pieces and measure difference `O_{E,p}(p^{-N_p})`.

The product approximations have
`O_{E,S}(prod_{p in T}(N_p+1))` disjoint rectangular pieces. The
intersection of any such rectangle with the fixed class in `H_T`
is empty or one coset of an open subgroup: all its defining balls
are subgroup translates. Its pullback to `Gamma` is therefore one
finite-index coset. The subgroup depth can differ from piece to
piece. Apply the uniform estimate (4.2) to each one, without
replacing them by cosets of a common deepest subgroup.

The total lattice error for either approximation is

\[
 O_{E,S}\left(R\prod_{p\in T}(N_p+1)\right)
            =O_{E,S}(R(1+\log R)^r).                \tag{5.2}
\]

The measure gap is at most the sum of the local gaps, using Haar
pushforward to each projected local closure. No independence is
needed for this union bound. The gap in main terms is

\[
 O_{E,S}\left(L\sum_{p\in T}p^{-N_p}\right)=O_{E,S}(R).
                                                               \tag{5.3}
\]

The two estimates squeeze the point count to its main term with
error (5.2). The finite global exceptions change only `O_E(1)`.
The original soluble set is invariant under negation, and targets
generically have exactly the two preimages `P,-P`. Dividing the
point count by two proves (5.1). \(\square\)

For a fixed class the exponent can be reduced further. Let
`P_a={p>=11:p=2 mod3, a is a square in Q_p}` as in
[EPE2](ELLIPTIC_EXCEPTION_PRIME_ELIGIBILITY.md), and put

\[
 r_a=\#\bigl(S_f\cap(\{2,3,5\}\cup\mathcal P_a)\bigr).
\]

EPE2 proves exact equality of the counts after erasing the other
places, together with equality of their Haar coefficients. Applying
EQA3 to that smaller prime set therefore replaces `r` in (5.1) by
`r_a`, with the same `mu_{S,a}`. This uses a class-specific redundant
prime criterion, not independence of the remaining prime tests.

The constants in (5.1) are for each fixed `S`. No estimate uniform in
growing prime sets, no effective numerical bound for the error
constant, and no all-primes tail estimate is asserted. In particular
the additional all-primes obligation of EFA3 remains separate.

## 6. Proof support

This deduction uses the existing exact local classifications, the
retained Milne local-group results, the canonical-height comparison,
and the uniform row-count proof above. Its unrestricted statements
do not rest on a new finite experiment. The finite quotient in (3.2)
must be certified when an explicit constant is evaluated; only its
existence is needed for rationality and (5.1). No formal verification
or external review is claimed.
