# Full finite-place asymptotics on the elliptic exceptional slice

Every fixed finite set of local conditions on the existing elliptic
slice has a leading logarithmic target count. Its constant is the Haar
measure of the soluble subset in the actual joint closure of the
rational elliptic group. The certified Hasse class has a positive
constant for every finite set of places. This result keeps the
correlations between places and does not impose the remaining places.

## 1. The joint closure and the counted targets

Retain the curve, parameter, actual target height and coefficient from
[the elliptic height count](ELLIPTIC_EXCEPTION_COUNT.md):

\[
 E:T^2=X^3-10758382X+6600759744,\quad
 v=\frac{3(X+3984)}{6400},\quad
 y=y(v),\quad C_E=\frac{\pi}{3\sqrt{R_x}}>0.          \tag{1.1}
\]

The normalization of `R_x` is that of EEC1. The target depends on `v`,
so `y(P)=y(-P)`. Put `Gamma=E(Q)` and retain the global two-isogeny
squareclass homomorphism of
[the elliptic arithmetic note](ELLIPTIC_EXCEPTION_ARITHMETIC.md):

\[
 \alpha:\Gamma\longrightarrow
 A=\{1,2,5p_0,10p_0\},\qquad p_0=2709173.           \tag{1.2}
\]

The group `A` has order four and exponent two. At an admissible point,
the quadratic field is `K_a=Q(sqrt(-3a))` when `alpha(P)=a`.

Let `S` be any fixed finite set of places of `Q`, and let `S_f` be its
finite primes. The real place needs no condition: a real cubic always
has a real root. Define the compact group

\[
 H_S=\overline{\{(\alpha(P),(P)_p):P\in\Gamma\}}
       \ \subseteq\ A\times\prod_{p\in S_f}E(\mathbb Q_p),       \tag{1.3}
\]

with normalized Haar measure `mu_S` on the **whole** group `H_S`.
The closure is in the displayed product topology. If `S_f` is empty,
take `H_S=A`. All these groups are profinite; the local elliptic groups
are compact and have the open formal subgroups used in ETS1.

For each finite prime let `B_p` be the finite subset of `E(Q_p)` where
the affine parameter or factors have a pole, or where

\[
 D(v)\Delta(v)\operatorname{Res}(q,h)=0.             \tag{1.4}
\]

On its complement the monic factors `q,h` are separable and give the
complete admissible fibre. Let `V_p` be the points in that complement
where at least one of `q,h` has a `Q_p` root. Define

\[
 U_{S,a}=\{(a',(P_p))\in H_S:
              a'=a,\ P_p\in V_p\text{ for every }p\in S_f\},
 \qquad \mu_{S,a}=\mu_S(U_{S,a}).                   \tag{1.5}
\]

For an empty product this just means `U_{S,a}={a}` in `A`.
The finite global inadmissible set is removed separately in the target
count. This convention avoids putting a global admissibility condition
on a tuple with no local coordinates.

Let `N_{S,a}(B)` count distinct targets on this slice of actual
projective height at most `B`, with full quadratic-times-cubic field
fibre, quadratic class `a`, and a local source point at every place of
`S`. Both field factors have degree greater than one, so these targets
have no rational source point. They are required to pass the specified
finite set, rather than to be everywhere locally soluble.

## 2. Local continuity and a null boundary

**Lemma EFA1 (the soluble set has Haar-null boundary).** Each `V_p` is
open in `E(Q_p)` and is both open and closed relative to
`E(Q_p) minus B_p`. Its boundary is contained in `B_p`. The set
`U_{S,a}` is open in `H_S`, and

\[
 \mu_S(\partial U_{S,a})=0.                         \tag{2.1}
\]

**Proof of the local assertion.** At a point of `V_p`, a separable
factor has a root with nonzero derivative. The coefficients are
continuous functions of the elliptic point. Strict Hensel lifting
therefore preserves the root in a sufficiently small neighbourhood.
This is the same lifting argument used in ELS1 of
[the finite-place sieve](ELLIPTIC_EXCEPTION_LOCAL_SIEVE.md).

Nonexistence of a root is also open where the coefficients are finite.
Here is the compactness step explicitly. Near any one point, all
coefficients of each monic polynomial are bounded by a fixed `p`-adic
constant `R>=1`. Every possible root has absolute value at most `R`:
outside that ball the leading term strictly dominates all other terms.
The ball is compact. If points without a local root could be approached
by soluble points, choose roots in that fixed ball. After passing to
one factor and a convergent subsequence, continuity would give a root
at the limiting point. This is a contradiction. Thus the complement
of the root-existence set is open as well. Including the finite
excluded set leaves only the boundary asserted above.

For the measure assertion, projection of `H_S` to any `E(Q_p)` is the
closure of the rational image there. This image is infinite: a rational
nontorsion point remains of infinite order in `E(Q_p)`, since equality
of rational points after embedding in `Q_p` is equality over `Q`.
The pushforward of Haar measure is normalized Haar measure on that
infinite compact image. Its singletons have measure zero; all translates
have equal measure, and arbitrarily many distinct singletons force
that common measure to vanish. Consequently each inverse image of
the finite set `B_p` has measure zero in `H_S`.

The fixed `A`-class condition is clopen. Away from those finitely many
inverse images, every local root condition is locally constant. Thus
the boundary of (1.5) lies in their union and has measure zero. The
empty-product case is clopen from the outset. \(\square\)

The only local lifting input is Hensel's lemma, retained in
[Milne, *Algebraic Number Theory*, version 3.00](sources/milne-ant300.pdf),
Proposition 7.31 and Theorem 7.32, printed pp.113--114. Compactness, local elliptic
subgroups and their finite quotients are covered by Milne,
[*Elliptic Curves* (2006)](sources/milne-elliptic-curves-2006.pdf),
Chapter II, Theorem 4.1, printed pp.62--63, as used in ETS1.
The null-boundary argument above is part of this proof and is not
inferred from finite computation.

## 3. The leading asymptotic for the entire finite-place set

**Theorem EFA2 (full finite-place target asymptotics).** For every
fixed finite `S` and every `a in A`,

\[
 \boxed{N_{S,a}(B)=\mu_{S,a}C_E\log B
                            +o_{E,S}(\log B).}     \tag{3.1}
\]

The measure is normalized on `H_S`, not separately on the class `a`.
In particular

\[
 0\le\mu_{S,a}\le\frac14.                          \tag{3.2}
\]

For the certified class `a_0=10p_0`, corresponding to
`K_0=Q(sqrt(-81275190))`, one has

\[
 0<\mu_{S,a_0}\le\frac14                           \tag{3.3}
\]

for every fixed finite `S`. The separate
[class-one theorem](ELLIPTIC_EXCEPTION_FINITE_LOCAL_CLASS_ONE.md)
also proves positivity for `a=1`; positivity for all finite sets in
the two remaining classes is unresolved.
For `S_f` empty, every class has measure `1/4`. The total count over
the four fields has leading coefficient

\[
 \mu_S\left(\bigcup_{a\in A}U_{S,a}\right)C_E,
 \qquad 0<\sum_{a\in A}\mu_{S,a}\le1.              \tag{3.4}
\]

**Proof of equidistribution in finite quotients.** Let `W` be an open
subgroup of `H_S`. Its quotient is finite, and the dense rational image
surjects onto it. The pullback `Lambda` in `Gamma` therefore has

\[
 [\Gamma:\Lambda]=[H_S:W].
\]

The canonical-height lattice argument of EEC1 counts the points in
each coset of `Lambda` with leading term

\[
 \frac{2C_E}{[\Gamma:\Lambda]}\log B
                         +O_{E,\Lambda}(\sqrt{\log B}).         \tag{3.5}
\]

The factor two counts elliptic points rather than target images.
This coset formula includes the rational torsion subgroup: its index
is the product of the projected free-lattice index and the torsion
quotient index. It does not require a displayed saturated generator
basis. Thus any fixed finite union `C` of cosets in `H_S` has point
main term `2 mu_S(C) C_E log B`.

**Proof of passage to the soluble set.** A compact profinite group has
a basis of clopen cosets of open subgroups. Any clopen subset is a
finite union of cosets of one common open subgroup. By EFA1 and
regularity of Haar measure, for every `epsilon>0` there are clopen
sets

\[
 C_-\subseteq U_{S,a}\subseteq C_+,
 \qquad \mu_S(C_+)-\mu_S(C_-)<\epsilon.             \tag{3.6}
\]

One can construct these approximations from progressively finer finite
quotients. Cosets wholly inside or wholly outside the soluble set are
decided; the undecided cosets shrink towards its null boundary. Their
total measure tends to zero. Equivalently the finite-quotient inner
and outer soluble proportions converge to the same `mu_{S,a}`.

Apply (3.5) to the two fixed sets in (3.6), take the lower and upper
limits as `B` grows, and then let `epsilon` tend to zero. This proves
the point-count limit with coefficient `2 mu_{S,a} C_E`.

Remove the finite global inadmissible set and the finite reducible-cubic
exceptions of the genus-four/Faltings theorem in the elliptic arithmetic
note. These changes contribute only `O_E(1)`. Negation preserves
`U_{S,a}`: the class group `A` has exponent two and every factor depends
on `v(P)=v(-P)`. Away from finitely many identifications, each target
has exactly the two preimages `P,-P`, by EEC1. Dividing the point count
by two proves (3.1).

Projection of Haar measure to `A` is uniform because that projection
is surjective. This proves (3.2) and the empty-product assertion.
The exact [elliptic Hasse seed](ELLIPTIC_EXCEPTION_HASSE_EXAMPLE.md)
has class `a_0` and is soluble at every place. Its image is in
`U_{S,a_0}`, which is open by EFA1. Every nonempty open subset of a
compact group has positive Haar measure, proving (3.3). Equivalently,
the finite-index cosets of ELS1 and its explicit even-multiplication
corollary give a positive logarithmic subcount contained in this full
finite-place set. Summing the disjoint class sets proves (3.4).
\(\square\)

[ESQ3](ELLIPTIC_EXCEPTION_SHARP_FINITE_LOCAL_ERROR.md) sharpens the
error in (3.1) to `O_(E,S)(sqrt(log B))`, using a primary bound for
two p-adic elliptic logarithms. The rational coefficient and
finite-local dichotomy are proved in EQA2--EQA2a of
[the quantitative note](ELLIPTIC_EXCEPTION_QUANTITATIVE_FINITE_LOCAL_ASYMPTOTICS.md).
These stronger finite-place statements leave the all-primes limit
and arithmetic tail in the following section separate.

[EOP1--EOP3](ELLIPTIC_EXCEPTION_ONE_PRIME_POSITIVITY.md) establish
positivity in all four classes at every individual place. They also
give `0<=mu_(S,a)-mu_(S union {p},a)<=K_E/log p` uniformly in the
finite prefix and class for large good `p`. Thus each fixed feasible
prefix retains proportion tending to one at a new prime, while a
growing prefix can have a shrinking denominator. This absolute
coefficient bound does not settle the conditional loss series or
the arithmetic tail below.

## 4. The all-primes limit and the exact tail condition

Fix \(a_0=10p_0\), so the quadratic field throughout this section is
\(K_0=\mathbb Q(\sqrt{-81275190})\). All counts use the same universe
of admissible full quadratic-times-cubic field targets on this
elliptic slice. Put
\[
 S_Y=\{p:p\text{ prime},\ p\le Y\},\qquad Y\in\mathbb Z,\ Y\ge5,
 \qquad N_Y(B)=N_{S_Y,a_0}(B),\quad \mu_Y=\mu_{S_Y,a_0}.
                                                               \tag{4.1}
\]
The real place is automatic. Write \(A_{K_0}(B)\) for the targets in
this universe that are soluble at every place, hence Hasse-failing.
Set \(L(B)=C_E\log B\), for \(B>1\).

The constants can be realized in one compact group. Let
\[
 H_\infty=
 \overline{\{(\alpha(P),(P)_p):P\in\Gamma\}}
       \subseteq A\times\prod_{p\text{ finite}}E(\mathbb Q_p)
                                                               \tag{4.2}
\]
with normalized Haar measure \(\nu_\infty\). Its projection
\(\rho_Y\) onto the finite coordinates is surjective onto
\(H_{S_Y}\): the compact image is closed and contains the dense
rational image. Thus its pushforward measure is \(\mu_{S_Y}\).
The sets
\[
 \mathscr U_Y=\rho_Y^{-1}(U_{S_Y,a_0}),\qquad
 \mathscr U_\infty=\bigcap_{Y\ge5}\mathscr U_Y
\]
are measurable and decrease. Continuity of probability measure
from above gives
\[
 \boxed{\mu_\infty:=\lim_{Y\to\infty}\mu_Y
        =\nu_\infty(\mathscr U_\infty),\qquad
                     0\le\mu_\infty\le\tfrac18.}      \tag{4.3}
\]
The sharper upper bound uses
[J2174](ELLIPTIC_EXCEPTION_JOINT_2_3_17_SIEVE.md): the fixed-class
coefficient at the finite set `2,3,17` is `1/8`, and all larger
prime sets impose additional conditions.
Every finite \(\mu_Y\) is positive by EFA2. This does not show that
\(\mu_\infty\) is positive.

**Theorem EFA3 (all-primes upper limit and tail equivalence).**
Define the nonnegative tail count
\[
 T_Y(B)=N_Y(B)-A_{K_0}(B).
\]
Then
\[
 0\le\liminf_{B\to\infty}\frac{A_{K_0}(B)}{L(B)}
 \le\limsup_{B\to\infty}\frac{A_{K_0}(B)}{L(B)}
 \le\mu_\infty.                                      \tag{4.4}
\]
Moreover the iterated tail limit exists and satisfies the exact identity
\[
 \boxed{\tau_\infty:=
   \lim_{Y\to\infty}\limsup_{B\to\infty}\frac{T_Y(B)}{L(B)}
   =\mu_\infty-
       \liminf_{B\to\infty}\frac{A_{K_0}(B)}{L(B)}.}    \tag{4.5}
\]
Consequently
\[
 \boxed{\tau_\infty=0
 \quad\Longleftrightarrow\quad
 \frac{A_{K_0}(B)}{C_E\log B}\longrightarrow\mu_\infty.} \tag{4.6}
\]

If \(\mu_\infty=0\), then
\(A_{K_0}(B)=o(\log B)\), and the tail condition in (4.6) holds
automatically. This leaves open whether the count is finite or
infinite with sublogarithmic growth. If \(\mu_\infty>0\) and
\(\tau_\infty=0\), then
\[
 A_{K_0}(B)\sim\mu_\infty C_E\log B,
\]
which proves infinitude. No value, vanishing, or positivity of
\(\mu_\infty\), and no vanishing of \(\tau_\infty\) in the positive
case, is established here.

**Proof.** EFA2 gives \(N_Y(B)/L(B)\to\mu_Y\) for every fixed \(Y\).
The counts \(N_Y(B)\) decrease in \(Y\), and
\(0\le A_{K_0}(B)\le N_Y(B)\). Taking upper limits in \(B\), and
then letting \(Y\) grow, proves (4.4).

Put \(a_-=\liminf_{B\to\infty}A_{K_0}(B)/L(B)\). It is finite by
(4.4). For fixed \(Y\), subtracting a sequence with a limit gives
\[
 \limsup_{B\to\infty}\frac{T_Y(B)}{L(B)}=\mu_Y-a_-.
\]
These quantities are nonnegative and decrease with \(Y\).
Their limit is \(\mu_\infty-a_-\), proving (4.5).
If that limit is zero, (4.4) squeezes the lower and upper limits
of the normalized count to \(\mu_\infty\). Conversely convergence
to \(\mu_\infty\) makes the right side of (4.5) zero. This proves
(4.6) and both consequences. \(\square\)

For each fixed height bound \(B\), there are only finitely many
targets. Each target failing an all-places test fails at some
finite prime. Hence \(T_Y(B)\) eventually vanishes for that fixed
\(B\). This pointwise fact does not interchange the two limits in
(4.5); the required tail condition first takes the upper limit in
height.

The arithmetic meaning of this tail is explicit. By
[ELV1](ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md), once the
tests at \(2,3,5\) have passed, an additional failure occurs exactly
at a prime \(p\equiv2\pmod3\), \(p\ge11\), with
\[
 v_p(128x-25865)\notin3\mathbb Z.
\]
Thus \(T_Y(B)\) counts the finite-place soluble targets in (4.1)
with at least one such bad prime \(p>Y\). No independence between
these prime conditions, or Euler product for (4.3), is asserted.

[EPE2](ELLIPTIC_EXCEPTION_PRIME_ELIGIBILITY.md) narrows this tail
in each fixed quadratic class. Write `x=M/C^2` and
`N=128M-25865C^2` in reduced coordinates. A denominator obstruction
can occur only at primes where the global class `a` is locally square;
a numerator obstruction additionally requires `51730` locally square.
Every other prime test is automatic, including a prime ramifying in
the class `a`. For the certified nontrivial class these two eligible
prime sets have respective densities `1/4` and `1/8` among all primes.
These are prime eligibility densities, not target densities or
independence statements. The uniform tail obligation remains.

The finiteness for each fixed rational cubeclass of \(128x-25865\)
in ELV4 does not prove \(\mu_\infty=0\). That theorem concerns
rational points with one rational cubeclass. The adelic tuples in
\(H_\infty\) need not give one rational value of this function.
In fact the entire embedded rational group is already a countable,
Haar-null, dense subset of the infinite compact group \(H_\infty\).
Its fixed-cubeclass finiteness therefore supplies no measure bound
for the full adelic set \(\mathscr U_\infty\).

There is a separate all-class version. Sum \(N_{S_Y,a}(B)\) over
the four classes, replace \(A_{K_0}(B)\) by the total Hasse count
\(A_E(B)\), and replace \(\mu_Y\) by
\(\mu_Y^{\rm tot}=\sum_{a\in A}\mu_{S_Y,a}\).
The same proof gives (4.3)--(4.6) with
\(0\le\mu_\infty^{\rm tot}\le2300089/10392564\), using
[J235173](ELLIPTIC_EXCEPTION_JOINT_2_3_5_17_SIEVE.md) at the finite
set `2,3,5,17`. This uses the total normalization
throughout; it does not mix a fixed-class numerator with a
total-class coefficient.

## 5. The meaning of the constant

The measure in (3.1) is an intrinsic limit of actual finite-quotient
proportions for this fixed elliptic group. It is not asserted to be an
Euler product of separate local densities. The diagonal rational image
can correlate the squareclass, component groups and local tests; the
two/three-adic correlation of
[ETF4--ETF5](ELLIPTIC_EXCEPTION_TWO_FIVE_ADIC_CLASSIFICATION.md) is an
explicit example. There is no error estimate uniform in growing `S`
or in the chosen finite-quotient resolution.

As `S` grows, the positive measure for `a_0` can shrink. Finite-place
asymptotics alone do not permit interchanging that limit with the
target-height limit. EFA3 gives the exact additional tail condition
for passage to the limiting coefficient. Its fulfillment in the
positive-measure case remains open, as does Hasse infinitude on this
elliptic slice.

## 6. Conditional prime losses and the limiting measure

Fix the certified class `a_0` throughout. Start with
`S_0={2,3,5}` and enumerate the remaining primes as
`p_1<p_2<...`. Put `S_j=S_0 union {p_1,...,p_j}` and
`m_j=mu_(S_j,a_0)`. Every `m_j` is positive by EFA2. Define the
actual conditional loss at the next prime by

\[
 \lambda_j=1-\frac{m_j}{m_{j-1}},\qquad 0\le\lambda_j<1.
                                                               \tag{6.1}
\]

**Proposition EFA4 (conditional-loss criterion).** With the same
limiting measure as EFA3,

\[
 \boxed{\mu_\infty=m_0\prod_{j\ge1}(1-\lambda_j),\qquad
 \mu_\infty>0\ \Longleftrightarrow\
                         \sum_{j\ge1}\lambda_j<\infty.}        \tag{6.2}
\]

Every prime erased by EPE2 has `lambda_j=0`, even after all the
preceding finite conditions have been imposed. If the actual
conditional losses satisfy `lambda_j>=c/p_j`, with fixed `c>0`,
on a subset of eligible primes having positive natural density
among all primes, then `mu_infinity=0` and
`A_(K0)(B)=o(log B)`.

**Proof.** The finite product in (6.2) telescopes to `m_j/m_0`;
letting `j` increase proves the first equality. The inequality
`-log(1-u)>=u` for `0<=u<1` proves vanishing of the product when
the series diverges. If the series converges, eventually `u<=1/2`
and `-log(1-u)<=2u`, so the product is positive. EPE2 gives exact
equality of the before-and-after finite-place counts at each erased
prime, hence equality of their leading coefficients. This proves
the zero-loss assertion. Positive prime density gives divergence of
the reciprocal-prime sum by partial summation. The asserted lower
bound then makes the loss series diverge, and EFA3 gives the
sublogarithmic Hasse bound. \(\square\)

Each factor in (6.2) is conditioned on all preceding tests in the
actual joint rational closure. A lower bound for an unconditioned
local bad proportion does not establish the hypothesis of EFA4.
The [ambient Haar diagnostic](ELLIPTIC_EXCEPTION_AMBIENT_HAAR_DIAGNOSTIC.md)
makes this distinction explicit. Summability of the actual losses
would prove positive limiting measure; a positive Hasse asymptotic
would still require EFA3's uniform arithmetic tail estimate.

The [triangular sieve](ELLIPTIC_EXCEPTION_TRIANGULAR_PRIME_SIEVE.md)
proves more than an ambient comparison: at sufficiently large primes
the new formal component is independent of the entire past, while
the finite reduction is measurable in lower-prime components.
When `e_p=1`, the actual loss is
`(1-1/(p^2+p+1))` times the conditional probability of reduction
at the pole or a zero of `128x-25865`. Those return probabilities
remain correlated and undetermined. Its predictable-return theorem
gives a precise recurrence criterion for a zero limiting measure.

ETP4--ETP5 remove the depth-one restriction: their exact hazard
kernel retains both the rational formal depth and each branch offset.
The intrinsic all-primes soluble set has finite accumulated hazard
almost surely. At every fixed cutoff they also give a count bound
from the finite low-hazard region, with an explicit exponentially
small remainder in the hazard threshold. Almost-sure weighted
divergence remains an unproved arithmetic input. The conditional
Wieferich-equidistribution comparison ETP6 is context, not a premise
of these unconditional hazard statements.

## 7. Proof support

The unrestricted proof is the local compactness/Hensel argument,
the null-boundary lemma, the torsion-aware lattice approximation,
and the elementary limiting identities in EFA3--EFA4. Its rank, height, finite-exception and local-seed inputs are
the existing canonical sources linked in sections 1--3. No new finite
experiment is needed for this topological and counting deduction.
This is an internal written proof, with no formal verification or
external review claimed.
